#include "Heat.hpp"

void Heat::init_mesh() {
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

  if (!fe) {
    fe = std::make_unique<FE_SimplexP<dim>>(r);
    quadrature = std::make_unique<QGaussSimplex<dim>>(r + 1);
    dof_handler.reinit(mesh); 
  }

  
  dof_handler.distribute_dofs(*fe);

  const IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
  const IndexSet locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);

  // Ricostruzione dello sparsity pattern
  TrilinosWrappers::SparsityPattern sparsity(locally_owned_dofs, MPI_COMM_WORLD);
  DoFTools::make_sparsity_pattern(dof_handler, sparsity);
  sparsity.compress();

  // Re-inizializzazione di matrice e vettori
  system_matrix.reinit(sparsity);
  system_rhs.reinit(locally_owned_dofs, MPI_COMM_WORLD);
  solution_owned.reinit(locally_owned_dofs, MPI_COMM_WORLD);
  solution_owned_old.reinit(locally_owned_dofs, MPI_COMM_WORLD);
  solution.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
}

void
Heat::assemble()
{
  // Number of local DoFs for each element.
  const unsigned int dofs_per_cell = fe->dofs_per_cell;

  // Number of quadrature points for each element.
  const unsigned int n_q = quadrature->size();

  FEValues<dim> fe_values(*fe,
                          *quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  // Local matrix and vector.
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  // Reset the global matrix and vector, just in case.
  system_matrix = 0.0;
  system_rhs    = 0.0;

  // Evaluation of the old solution on quadrature nodes of current cell.
  std::vector<double> solution_old_values(n_q);

  // Evaluation of the gradient of the old solution on quadrature nodes of
  // current cell.
  std::vector<Tensor<1, dim>> solution_old_grads(n_q);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);

      cell_matrix = 0.0;
      cell_rhs    = 0.0;

      // Evaluate the old solution and its gradient on quadrature nodes.
      fe_values.get_function_values(solution, solution_old_values);
      fe_values.get_function_gradients(solution, solution_old_grads);

      for (unsigned int q = 0; q < n_q; ++q)
        {
          const double mu_loc = mu(fe_values.quadrature_point(q));

          // time è il tempo "vecchio" (t_n), time + delta_t è il tempo "nuovo" (t_{n+1})
          const double f_old_loc = f(fe_values.quadrature_point(q), time);
          const double f_new_loc = f(fe_values.quadrature_point(q), time + delta_t);


          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                  // Time derivative.
                  cell_matrix(i, j) += (1.0 / delta_t) *             //
                                       fe_values.shape_value(i, q) * //
                                       fe_values.shape_value(j, q) * //
                                       fe_values.JxW(q);

                  // Diffusion.
                  cell_matrix(i, j) += theta * mu_loc *                             //
                    scalar_product(fe_values.shape_grad(i, q),   //
                                   fe_values.shape_grad(j, q)) * //
                    fe_values.JxW(q);

                }

              // Time derivative.
              cell_rhs(i) += (1.0 / delta_t) *             //
                             fe_values.shape_value(i, q) * //
                             solution_old_values[q] *      //
                             fe_values.JxW(q);

                             

              // Diffusion.
              cell_rhs(i) -= (1.0 - theta) * mu_loc *                   //
                             scalar_product(fe_values.shape_grad(i, q), //
                                            solution_old_grads[q]) *    //
                             fe_values.JxW(q);

              // Forcing term.
              cell_rhs(i) +=
                (theta * f_new_loc + (1.0 - theta) * f_old_loc) * //
                fe_values.shape_value(i, q) *                     //
                fe_values.JxW(q);
            }
        }


        
      cell->get_dof_indices(dof_indices);

      system_matrix.add(dof_indices, cell_matrix);
      system_rhs.add(dof_indices, cell_rhs);
    }

  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);

 
}

void
Heat::solve_linear_system()
{
  
  TrilinosWrappers::PreconditionSSOR preconditioner;
  preconditioner.initialize(
    system_matrix, TrilinosWrappers::PreconditionSSOR::AdditionalData(1.0));

  ReductionControl solver_control(/* maxiter = */ 10000,
                                  /* tolerance = */ 1.0e-12,
                                  /* reduce = */ 1.0e-6);

  SolverCG<TrilinosWrappers::MPI::Vector> solver(solver_control);

  solver.solve(system_matrix, solution_owned, system_rhs, preconditioner);
  
  pcout << solver_control.last_step() << " CG iterations" << std::endl;
}


void Heat::refine_mesh() {
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

  SolutionTransfer<dim, TrilinosWrappers::MPI::Vector> sol_trans(dof_handler);
  sol_trans.prepare_for_coarsening_and_refinement(solution_owned);

  mesh.execute_coarsening_and_refinement();

  // FIX CRITICO: Salviamo il vecchio vettore prima che setup_system lo distrugga!
  TrilinosWrappers::MPI::Vector old_solution_owned = solution_owned;

  // Ricostruzione sistema (questo riempirà solution_owned di zeri)
  setup_system(); 

  TrilinosWrappers::MPI::Vector interpolated_sol;
  interpolated_sol.reinit(dof_handler.locally_owned_dofs(), MPI_COMM_WORLD);
  
  // Interpoliamo usando la vecchia soluzione salvata (old_solution_owned)
  sol_trans.interpolate(old_solution_owned, interpolated_sol);
  
  solution_owned = interpolated_sol;
  solution = solution_owned;
}


void
Heat::output() const
{
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

void Heat::run() {
  init_mesh();
  setup_system();
  VectorTools::interpolate(dof_handler, FunctionU0(), solution_owned);
  solution = solution_owned;

  while (time < T) {
    bool step_accepted = false;
    double error_t = 0.0;
    
    while (!step_accepted) {
      assemble();
      solve_linear_system();

      // Calcoliamo la differenza tra u_new e u_old (errore temporale stimato)
      TrilinosWrappers::MPI::Vector diff = solution_owned;
      diff.add(-1.0, solution_owned_old);
      error_t = diff.l2_norm() / std::sqrt(dof_handler.n_dofs()); //L2 norm

      if (error_t > tol_time_max && delta_t > dt_min) {
        // Errore troppo grande: riduciamo dt e non avanziamo nel tempo
        delta_t /= 2.0;
        pcout << "  Step rejected! New delta_t = " << delta_t << std::endl;
        solution_owned = solution_owned_old;
      } 
      else {
        // Errore accettabile o dt minimo raggiunto: procediamo
        step_accepted = true;
        time += delta_t;
        timestep_number++;
        
        if (error_t < tol_time_min && delta_t < dt_max)
          delta_t *= 1.2; 
      }
    }

    // Dopo uno step accettato, gestiamo l'adattività spaziale
    if (timestep_number % 5 == 0) {
      refine_mesh();
    }

    solution = solution_owned;
    solution_owned_old = solution_owned;
    output();
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
