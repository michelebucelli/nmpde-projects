#include "Heat.hpp"

// Constructor definition
Heat::Heat(const std::string                              &mesh_file_name_, 
           const unsigned int                              &r_,
           const double                                    &T_,
           const double                                    &theta_,
           const double                                    &delta_t_,
           const std::function<double(const Point<dim> &)> &mu_,
           const std::function<double(const Point<dim> &, const double &)> &f_)
  : mesh_file_name(mesh_file_name_)
  , r(r_)
  , T(T_)
  , theta(theta_)
  , delta_t(delta_t_)
  , mu(mu_)
  , f(f_)
  , mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
  , mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
  , mesh(MPI_COMM_WORLD)
  , pcout(std::cout, mpi_rank == 0)
  , computing_timer(MPI_COMM_WORLD, pcout, TimerOutput::never, TimerOutput::wall_times)
{}

// Reads the meh from file fopr the iniziaitazion of the domain
void Heat::init_mesh() {
  TimerOutput::Scope s(computing_timer, "init_mesh");
  pcout << "Initializing the mesh from file" << std::endl;
  Triangulation<dim> mesh_serial;
  GridIn<dim> grid_in;
  grid_in.attach_triangulation(mesh_serial);
  std::ifstream mesh_file(mesh_file_name);
  if (!mesh_file.is_open()) {
    throw std::runtime_error("Could not open mesh file: " + mesh_file_name);
  }
  grid_in.read_msh(mesh_file);

  GridTools::partition_triangulation(mpi_size, mesh_serial);
  const auto construction_data = TriangulationDescription::Utilities::
  create_description_from_triangulation(mesh_serial, MPI_COMM_WORLD);
  mesh.create_triangulation(construction_data);
}

void Heat::setup_system() {
  TimerOutput::Scope s(computing_timer, "setup_system");

  if (!fe) {
    fe = std::make_unique<FE_SimplexP<dim>>(r);
    quadrature = std::make_unique<QGaussSimplex<dim>>(r + 1);
    dof_handler.reinit(mesh); 
  }

  dof_handler.distribute_dofs(*fe);

  const IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
  const IndexSet locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);

  // Re-construction of the sparsity pattern
  TrilinosWrappers::SparsityPattern sparsity(locally_owned_dofs, MPI_COMM_WORLD);
  DoFTools::make_sparsity_pattern(dof_handler, sparsity);
  sparsity.compress();

  // Re-initialization of matrices and vectors
  mass_matrix.reinit(sparsity);
  stiffness_matrix.reinit(sparsity);
  system_rhs.reinit(locally_owned_dofs, MPI_COMM_WORLD);
  solution_owned.reinit(locally_owned_dofs, MPI_COMM_WORLD);
  solution_owned_old.reinit(locally_owned_dofs, MPI_COMM_WORLD);
  solution.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
}

// Computes the mass and stifness matrices to be re-used
void Heat::assemble_matrices() {
  TimerOutput::Scope s(computing_timer, "assemble_matrices");
  
  // Puts at zero the global matrices
  mass_matrix = 0.0;
  stiffness_matrix = 0.0;

  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q = quadrature->size();

  FEValues<dim> fe_values(*fe, *quadrature,
                          update_values | update_gradients | update_quadrature_points | update_JxW_values);

  FullMatrix<double> cell_mass(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_stiffness(dofs_per_cell, dofs_per_cell);
  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators()) {
    if (!cell->is_locally_owned()) continue;

    fe_values.reinit(cell);
    cell_mass = 0.0;
    cell_stiffness = 0.0;

    for (unsigned int q = 0; q < n_q; ++q) {
      const double mu_loc = mu(fe_values.quadrature_point(q));

      for (unsigned int i = 0; i < dofs_per_cell; ++i) {
        for (unsigned int j = 0; j < dofs_per_cell; ++j) {
          // Matrice di Massa M
          cell_mass(i, j) += fe_values.shape_value(i, q) * fe_values.shape_value(j, q) * fe_values.JxW(q);

          // Matrice di Rigidezza A (Diffusione)
          cell_stiffness(i, j) += mu_loc * scalar_product(fe_values.shape_grad(i, q), 
                                                 fe_values.shape_grad(j, q)) * fe_values.JxW(q);
        }
      }
    }
    
    cell->get_dof_indices(dof_indices);
    mass_matrix.add(dof_indices, cell_mass);
    stiffness_matrix.add(dof_indices, cell_stiffness);
  }

  mass_matrix.compress(VectorOperation::add);
  stiffness_matrix.compress(VectorOperation::add);
}

// System assembled from the computed matrices and assembling the rhs
void Heat::assemble_system_and_rhs() {
  TimerOutput::Scope s(computing_timer, "assemble_system_and_rhs");

  // 1. Compute the system matrix: S = (1/dt) * M + theta * A
  // Copiamo M in S per efficienza, la scaliamo e aggiungiamo A
  system_matrix.copy_from(mass_matrix);
  system_matrix *= (1.0 / delta_t);
  system_matrix.add(theta, stiffness_matrix);

  // 2. Compute the contribution of the old solution at the RHS by mvm
  // RHS = ( (1/dt)*M - (1 - theta)*A ) * u_old
  system_rhs = 0.0;
  TrilinosWrappers::MPI::Vector tmp(solution_owned); // Vettore temporaneo con stesso layout

  // + (1/dt) * M * u_old
  mass_matrix.vmult(tmp, solution_owned_old);
  system_rhs.add(1.0 / delta_t, tmp);

  // - (1 - theta) * A * u_old
  stiffness_matrix.vmult(tmp, solution_owned_old);
  system_rhs.add(-(1.0 - theta), tmp);

  // 3. Compute f
  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q = quadrature->size();
  FEValues<dim> fe_values(*fe, *quadrature, update_values | update_quadrature_points | update_JxW_values);

  Vector<double> cell_rhs(dofs_per_cell);
  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators()) {
    if (!cell->is_locally_owned()) continue;

    fe_values.reinit(cell);
    cell_rhs = 0.0;

    for (unsigned int q = 0; q < n_q; ++q) {
      const double f_old_loc = f(fe_values.quadrature_point(q), time);
      const double f_new_loc = f(fe_values.quadrature_point(q), time + delta_t);

      for (unsigned int i = 0; i < dofs_per_cell; ++i) {
        cell_rhs(i) += (theta * f_new_loc + (1.0 - theta) * f_old_loc) * fe_values.shape_value(i, q) * fe_values.JxW(q);
      }
    }
    cell->get_dof_indices(dof_indices);
    system_rhs.add(dof_indices, cell_rhs);
  }
  
  system_rhs.compress(VectorOperation::add);
}

// Solving the system with the AMG precoditioner
void Heat::solve_linear_system()
{
  TimerOutput::Scope s(computing_timer, "solve_linear_system");
  
  
  TrilinosWrappers::PreconditionAMG preconditioner;
  
  // Configuriamo i parametri per un problema di tipo parabolico
  TrilinosWrappers::PreconditionAMG::AdditionalData amg_data;
  amg_data.elliptic = true;
  amg_data.smoother_sweeps = 2;
  amg_data.smoother_type = "ML symmetric Gauss-Seidel";
  
  preconditioner.initialize(system_matrix, amg_data);

  ReductionControl solver_control(/* maxiter = */ 10000,
                                  /* tolerance = */ 1.0e-12,
                                  /* reduce = */ 1.0e-6);

  // Using CG for system with SPD matrix (b = 0)
  SolverCG<TrilinosWrappers::MPI::Vector> solver(solver_control);

  solver.solve(system_matrix, solution_owned, system_rhs, preconditioner);
  
  pcout << "  -> " << solver_control.last_step() << " CG iterations (AMG Preconditioned)" << std::endl;
}

// Function called to dynamically refine the mesh based on the spatial error (Kelly estimator)
void Heat::refine_mesh() {
  TimerOutput::Scope s(computing_timer, "refine_mesh");
  pcout << "  Estimating error and refining..." << std::endl;

  Vector<float> estimated_error_per_cell(mesh.n_active_cells());
  KellyErrorEstimator<dim>::estimate(dof_handler,
                                     QGauss<dim - 1>(r + 1),
                                     {},
                                     solution, 
                                     estimated_error_per_cell);

  GridRefinement::refine_and_coarsen_fixed_number(mesh, 
                                                  estimated_error_per_cell, 
                                                  0.3, 0.1);

  // 1. Saving a copy of the current solution
  TrilinosWrappers::MPI::Vector old_solution = solution_owned;

  // 2. SolutionTransfer catches the vectors on the old mesh
  SolutionTransfer<dim, TrilinosWrappers::MPI::Vector> sol_trans(dof_handler);
  sol_trans.prepare_for_coarsening_and_refinement(old_solution);

  // 3. change the mesh
  mesh.execute_coarsening_and_refinement();

  // 4. DoF re-allocated
  //    viene riallocato con le nuove dimensioni (pieno di zeri)
  setup_system(); 
  assemble_matrices(); // Compute A and M in the new mesh

  // 5. SolutionTransfer interpolatesthe values
  sol_trans.interpolate(old_solution, solution_owned);
  
  // 6. Upload the solution with ghost values
  solution = solution_owned;
}

// output gives the solution at equally spaced time step
void Heat::output() const
{
  TimerOutput::Scope s(computing_timer, "output");
  DataOut<dim> data_out;

  data_out.add_data_vector(dof_handler, solution, "solution");

  // Add vector for parallel partition.
  std::vector<unsigned int> partition_int(mesh.n_active_cells());
  GridTools::get_subdomain_association(mesh, partition_int);
  const Vector<double> partitioning(partition_int.begin(), partition_int.end());
  data_out.add_data_vector(partitioning, "partitioning");

  data_out.build_patches();

  const std::filesystem::path mesh_path(mesh_file_name);
  const std::string output_file_name = "output-" + mesh_path.stem().string();

  data_out.write_vtu_with_pvtu_record(/* folder = */ "./",
                                      /* basename = */ output_file_name,
                                      /* index = */ timestep_number,
                                      MPI_COMM_WORLD);
}

// Function run where you copute the L2 norm error and compute the new delta_t if needed
void Heat::run() {
  TimerOutput::Scope s(computing_timer, "run");

  init_mesh();
  setup_system();
  assemble_matrices();
  VectorTools::interpolate(dof_handler, FunctionU0(), solution_owned);
  solution = solution_owned;
  solution_owned_old = solution_owned;

  const double output_interval = 0.05; // intervallo in cui vengono prese due soluzioni successive
  double next_output_time = 0.0;       // Forza il salvataggio del dato iniziale a t=0

  while (time < T) {
    bool step_accepted = false;
    double error_t = 0.0;
    
    while (!step_accepted) {
      assemble_system_and_rhs();
      solve_linear_system();

      // Calcoliamo u_new - u_old
      TrilinosWrappers::MPI::Vector diff = solution_owned;
      diff.add(-1.0, solution_owned_old);

      // ||e||_{L^2}^2 = e^T * M * e
      TrilinosWrappers::MPI::Vector M_diff(solution_owned);
      mass_matrix.vmult(M_diff, diff); // M_diff = M * e
      
      
      double error_l2_sq = diff * M_diff; 
      error_t = std::sqrt(error_l2_sq);


      if (error_t > tol_time_max && delta_t > dt_min) {
        // Error too big, reduce delta_t and go ahead 
        delta_t /= 2.0;
        pcout << "  Step rejected! New delta_t = " << delta_t << std::endl;
        solution_owned = solution_owned_old;
      } 
      else {
        // Acceptable error
        step_accepted = true;
        time += delta_t;
        timestep_number++;
        
        if (error_t < tol_time_min && delta_t < dt_max)
          delta_t *= 1.2; 
      }
    }

    // after a accepted step refine the mesh
    if (timestep_number % 5 == 0) {
      refine_mesh();
    }

    solution = solution_owned;
    solution_owned_old = solution_owned;
    
    if (time >= next_output_time || time >= T) {
      pcout << "  --> Salvataggio output a t = " << time << std::endl;
      output();
      
      // Ulpoad the solution with a while for safety reason
      while (next_output_time <= time && next_output_time < T) {
        next_output_time += output_interval;
      }
    }
  }
}


double
Heat::compute_error(const VectorTools::NormType &norm_type,
                                 const Function<dim> &exact_solution) const
{
  const QGaussSimplex<dim> quadrature_error(r + 2);

  // For triangular meshes, we need to explicitly provide a Mapping (the mapping
  // between the reference element and the individual elements in the mesh) to
  // integrate_difference. These two commands construct the most basic linear
  // mapping.
  FE_SimplexP<dim> fe_linear(1);
  MappingFE        mapping(fe_linear);

  Vector<double> error_per_cell(mesh.n_active_cells());
  VectorTools::integrate_difference(mapping,
                                    dof_handler,
                                    solution,
                                    exact_solution,
                                    error_per_cell,
                                    quadrature_error,
                                    norm_type);

  const double error =
    VectorTools::compute_global_error(mesh, error_per_cell, norm_type);

  return error;
}

