#include "NavierStokes3D.hpp"

constexpr unsigned int PRECONDITIONER = 4;
// 1 = Identity, 2 = BlockTriangular, 3 = BlockDiagonal, 4 = aSIMPLE, 5 = aYosida

void
NavierStokes3D::setup()
{
  // Create the mesh.
  {
    pcout << "\nInitializing the mesh" << std::endl;

    Triangulation<dim> mesh_serial;

    GridIn<dim> grid_in;
    grid_in.attach_triangulation(mesh_serial);

    std::ifstream grid_in_file(mesh_file_name);
    grid_in.read_msh(grid_in_file);

    GridTools::partition_triangulation(mpi_size, mesh_serial);
    const auto construction_data = TriangulationDescription::Utilities::
      create_description_from_triangulation(mesh_serial, MPI_COMM_WORLD);
    mesh.create_triangulation(construction_data);

    pcout << "  Number of elements = " << mesh.n_global_active_cells()
          << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the finite element space.
  {
    pcout << "Initializing the finite element space" << std::endl;

    const FE_SimplexP<dim> fe_scalar_velocity(degree_velocity);
    const FE_SimplexP<dim> fe_scalar_pressure(degree_pressure);
    fe = std::make_unique<FESystem<dim>>(fe_scalar_velocity,
                                         dim,
                                         fe_scalar_pressure,
                                         1);

    pcout << "  Velocity degree:           = " << fe_scalar_velocity.degree
          << std::endl;
    pcout << "  Pressure degree:           = " << fe_scalar_pressure.degree
          << std::endl;
    pcout << "  DoFs per cell              = " << fe->dofs_per_cell
          << std::endl;

    quadrature = std::make_unique<QGaussSimplex<dim>>(fe->degree + 1);

    pcout << "  Quadrature points per cell = " << quadrature->size()
          << std::endl;

    quadrature_face = std::make_unique<QGaussSimplex<dim - 1>>(fe->degree + 1);

    pcout << "  Quadrature points per face = " << quadrature_face->size()
          << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the DoF handler.
  {
    pcout << "Initializing the DoF handler" << std::endl;

    dof_handler.reinit(mesh);
    dof_handler.distribute_dofs(*fe);

    // We want to reorder DoFs so that all velocity DoFs come first, and then
    // all pressure DoFs.
    std::vector<unsigned int> block_component(dim + 1, 0);
    block_component[dim] = 1;
    DoFRenumbering::component_wise(dof_handler, block_component);

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);

    // Besides the locally owned and locally relevant indices for the whole
    // system (velocity and pressure), we will also need those for the
    // individual velocity and pressure blocks.
    std::vector<types::global_dof_index> dofs_per_block =
      DoFTools::count_dofs_per_fe_block(dof_handler, block_component);
    const unsigned int n_u = dofs_per_block[0];
    const unsigned int n_p = dofs_per_block[1];

    block_owned_dofs.resize(2);
    block_relevant_dofs.resize(2);
    block_owned_dofs[0]    = locally_owned_dofs.get_view(0, n_u);
    block_owned_dofs[1]    = locally_owned_dofs.get_view(n_u, n_u + n_p);
    block_relevant_dofs[0] = locally_relevant_dofs.get_view(0, n_u);
    block_relevant_dofs[1] = locally_relevant_dofs.get_view(n_u, n_u + n_p);

    pcout << "  Number of DoFs: " << std::endl;
    pcout << "    velocity = " << n_u << std::endl;
    pcout << "    pressure = " << n_p << std::endl;
    pcout << "    total    = " << n_u + n_p << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the linear system.
  {
    pcout << "Initializing the linear system" << std::endl;

    pcout << "  Initializing the sparsity pattern" << std::endl;

    // Velocity DoFs interact with other velocity DoFs (the weak formulation has
    // terms involving u times v), and pressure DoFs interact with velocity DoFs
    // (there are terms involving p times v or u times q). However, pressure
    // DoFs do not interact with other pressure DoFs (there are no terms
    // involving p times q). We build a table to store this information, so that
    // the sparsity pattern can be built accordingly.
    Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1);
    for (unsigned int c = 0; c < dim + 1; ++c)
      {
        for (unsigned int d = 0; d < dim + 1; ++d)
          {
            if (c == dim && d == dim) // pressure-pressure term
              coupling[c][d] = DoFTools::none;
            else // other combinations
              coupling[c][d] = DoFTools::always;
          }
      }

    TrilinosWrappers::BlockSparsityPattern sparsity(block_owned_dofs,
                                                    MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler, coupling, sparsity);
    sparsity.compress();

    // We also build a sparsity pattern for the pressure mass matrix.
    for (unsigned int c = 0; c < dim + 1; ++c)
      {
        for (unsigned int d = 0; d < dim + 1; ++d)
          {
            if (c == dim && d == dim) // pressure-pressure term
              coupling[c][d] = DoFTools::always;
            else // other combinations
              coupling[c][d] = DoFTools::none;
          }
      }
    TrilinosWrappers::BlockSparsityPattern sparsity_pressure_mass(
      block_owned_dofs, MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler,
                                    coupling,
                                    sparsity_pressure_mass);
    sparsity_pressure_mass.compress();

    pcout << "  Initializing the matrices" << std::endl;
    system_matrix.reinit(sparsity);
    pressure_mass.reinit(sparsity_pressure_mass);
    velocity_mass.reinit(sparsity);

    pcout << "  Initializing the system right-hand side" << std::endl;
    system_rhs.reinit(block_owned_dofs, MPI_COMM_WORLD);
    pcout << "  Initializing the solution vector" << std::endl;
    solution_owned.reinit(block_owned_dofs, MPI_COMM_WORLD);
    solution.reinit(block_owned_dofs, block_relevant_dofs, MPI_COMM_WORLD);
  }
}

void
NavierStokes3D::assemble()
{
  // pcout << "===============================================" << std::endl;
  pcout << "  Assembling the system" << std::endl;

  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q           = quadrature->size();
  const unsigned int n_q_face      = quadrature_face->size();

  FEValues<dim>     fe_values(*fe,
                          *quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);
  FEFaceValues<dim> fe_face_values(*fe,
                                   *quadrature_face,
                                   update_values | update_normal_vectors |
                                     update_JxW_values);

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_pressure_mass_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_velocity_mass(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  system_matrix = 0.0;
  system_rhs    = 0.0;
  pressure_mass = 0.0;
  velocity_mass = 0.0;

  FEValuesExtractors::Vector velocity(0);
  FEValuesExtractors::Scalar pressure(dim);

  // Evaluation of the old solution on quadrature nodes of current cell.
  std::vector<Tensor<1, dim>> solution_old_values(n_q);

  // Evaluation of the gradient of the old solution on quadrature nodes of
  // current cell.
  std::vector<Tensor<dim-1, dim>> solution_old_grads(n_q);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);

      cell_matrix               = 0.0;
      cell_rhs                  = 0.0;
      cell_pressure_mass_matrix = 0.0;
      cell_velocity_mass        = 0.0;

      // Evaluate the old solution and its gradient on quadrature nodes.
      fe_values[velocity].get_function_values(solution, solution_old_values);
      fe_values[velocity].get_function_gradients(solution, solution_old_grads);

      for (unsigned int q = 0; q < n_q; ++q)
        {

           const Tensor<1, dim> f_old_loc = data.f(fe_values.quadrature_point(q), time - delta_t);
           const Tensor<1, dim> f_new_loc = data.f(fe_values.quadrature_point(q), time);

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                  // Time derivative: (1/dt) (u^{n+1}, v)
                  cell_matrix(i, j) += (1.0 / delta_t) *             //
                                       fe_values[velocity].value(i, q) * //
                                       fe_values[velocity].value(j, q) * //
                                       fe_values.JxW(q);

                  // Viscosity term: nu * theta * (grad u^{n+1}, grad v)
                  cell_matrix(i, j) += data.nu * theta *
                                       scalar_product(fe_values[velocity].gradient(i, q),
                                                      fe_values[velocity].gradient(j, q)) *
                                       fe_values.JxW(q);


                  // Convection term: ((u^n · grad) u^{n+1}, v)
                  cell_matrix(i, j) += theta *
                                      (fe_values[velocity].gradient(j, q) * solution_old_values[q]) 
                                      * fe_values[velocity].value(i, q) 
                                      * fe_values.JxW(q);

                  // Pressure term in the momentum equation: -(p^{n+1}, div v)
                  cell_matrix(i, j) -= fe_values[velocity].divergence(i, q) *
                                       fe_values[pressure].value(j, q) *
                                       fe_values.JxW(q);

                  // Pressure term in the continuity equation: -(div u^{n+1}, q)
                  cell_matrix(i, j) -= fe_values[velocity].divergence(j, q) *
                                       fe_values[pressure].value(i, q) *
                                       fe_values.JxW(q);

                  // Pressure mass matrix
                   cell_pressure_mass_matrix(i, j) += // (data.nu + (1.0 / delta_t)) *
                                                      1.0 / delta_t *
                                                      fe_values[pressure].value(i, q) *
                                                      fe_values[pressure].value(j, q) *
                                                      fe_values.JxW(q);
                
                  // Only for Yosida preconditioner
                  if (PRECONDITIONER == 5)
                    {
                      // Velocity mass matrix: (u, v)
                      cell_velocity_mass(i, j) += fe_values[velocity].value(i, q) *
                                                  fe_values[velocity].value(j, q) *
                                                  fe_values.JxW(q);
                    }

                  }

                  // Time derivative: (1/dt)(u^n, v)
                  cell_rhs(i) += (1.0 / delta_t) *             //
                             fe_values[velocity].value(i, q) * //
                             solution_old_values[q] *      //
                             fe_values.JxW(q);

                  // - nu (1-theta) (grad u^n, grad v)
                  cell_rhs(i) -= (1.0 - theta) * data.nu *                   //
                             scalar_product(fe_values[velocity].gradient(i, q), //
                                            solution_old_grads[q]) *    //
                             fe_values.JxW(q);
                  
                  // Forcing term: (f^{n+theta}, v)
                  cell_rhs(i) +=  (theta * f_new_loc + (1.0 - theta) * f_old_loc) *                      //
                                fe_values[velocity].value(i, q) * //
                                fe_values.JxW(q);
              
                  // RHS convection explicit part
                  const Tensor<1, dim> conv_old = solution_old_grads[q] * solution_old_values[q];

                  cell_rhs(i) -= (1.0 - theta) *
                                (conv_old * fe_values[velocity].value(i, q)) *
                                fe_values.JxW(q);
                      
            }
        }

      // Boundary integral for Neumann BCs.
      if (cell->at_boundary())
        {
          for (unsigned int f = 0; f < cell->n_faces(); ++f)
            {
              if (cell->face(f)->at_boundary() &&
                  (cell->face(f)->boundary_id() == 2 ))
                {
                  fe_face_values.reinit(cell, f);

                  for (unsigned int q = 0; q < n_q_face; ++q)
                    {
                      const Tensor<1, dim> h_val = data.h(fe_face_values.quadrature_point(q), time);

                    for (unsigned int i = 0; i < dofs_per_cell; ++i)
                        {
                             cell_rhs(i) += scalar_product(h_val,
                              fe_face_values[velocity].value(i, q)) *
                              fe_face_values.JxW(q);
                        }
                    }
                }
            }
        }

      cell->get_dof_indices(dof_indices);

      system_matrix.add(dof_indices, cell_matrix);
      system_rhs.add(dof_indices, cell_rhs);
      pressure_mass.add(dof_indices, cell_pressure_mass_matrix);
      if (PRECONDITIONER == 5)
        velocity_mass.add(dof_indices, cell_velocity_mass);
    }

  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
  pressure_mass.compress(VectorOperation::add);
  if (PRECONDITIONER == 5)
    velocity_mass.compress(VectorOperation::add);

  // Dirichlet boundary conditions.
  {
    std::map<types::global_dof_index, double>           boundary_values;
    std::map<types::boundary_id, const Function<dim> *> boundary_functions;
    
    InletVelocity inlet(data.u_in);
    inlet.set_time(time);
    Functions::ZeroFunction<dim> zero_function(dim + 1);
 

    boundary_functions[1] = &inlet; // inlet
    boundary_functions[3] = &zero_function; // top/bottom walls
    boundary_functions[4] = &zero_function; // cylinder

    ComponentMask mask_velocity(dim + 1, true);
    mask_velocity.set(dim, false);

    VectorTools::interpolate_boundary_values(dof_handler,
                                             boundary_functions,
                                             boundary_values,
                                             mask_velocity);
    

    MatrixTools::apply_boundary_values(
      boundary_values, system_matrix, solution_owned, system_rhs, false);
  }
}

void
NavierStokes3D::solve()
{
  //pcout << "===============================================" << std::endl;

  SolverControl solver_control(50000, 1e-4 * system_rhs.l2_norm());
  SolverFGMRES<TrilinosWrappers::MPI::BlockVector> solver(solver_control);
 

  if (PRECONDITIONER == 1) 
  {
    pcout << "  Using: Identity preconditioner" << std::endl;
    advanced_preconditioner = std::make_unique<NSPreconditioners::PreconditionIdentity>();
  }
  else if (PRECONDITIONER == 2) 
  {
    pcout << "  Using: Block Triangular preconditioner" << std::endl;
    auto concrete_prec = std::make_unique<NSPreconditioners::PreconditionBlockTriangular>();
    concrete_prec->initialize(system_matrix.block(0, 0),
                              pressure_mass.block(1, 1),
                              system_matrix.block(1, 0));
    advanced_preconditioner = std::move(concrete_prec);
  }
  else if (PRECONDITIONER == 3) 
  {
    pcout << "  Using: Block Diagonal preconditioner" << std::endl;
    auto concrete_prec = std::make_unique<NSPreconditioners::PreconditionBlockDiagonal>();
    concrete_prec->initialize(system_matrix.block(0, 0), 
                              pressure_mass.block(1, 1));
    advanced_preconditioner = std::move(concrete_prec);
  }
  else if (PRECONDITIONER == 4) 
  {
    pcout << "  Using: aSIMPLE preconditioner" << std::endl;
    auto concrete_prec = std::make_unique<NSPreconditioners::PreconditionaSIMPLE>();
    concrete_prec->initialize(system_matrix.block(0, 0),
                              system_matrix.block(0, 1),
                              system_matrix.block(1, 0),
                              solution_owned,
                              0.6 /* alpha */);
    advanced_preconditioner = std::move(concrete_prec);
  }
  else if (PRECONDITIONER == 5) 
  {
    pcout << "  Using: aYosida preconditioner" << std::endl;
    auto concrete_prec = std::make_unique<NSPreconditioners::PreconditionaYosida>();
    concrete_prec->initialize(system_matrix.block(0, 0),
                               system_matrix.block(0, 1),
                               system_matrix.block(1, 0),
                               velocity_mass.block(0, 0), 
                               solution_owned, 
                               delta_t);                    
  advanced_preconditioner = std::move(concrete_prec);
  }

  pcout << "  Solving the linear system" << std::endl;
  solver.solve(system_matrix, solution_owned, system_rhs, *advanced_preconditioner);
  pcout << "   " << (solver_control.last_step()+1) << "  FGMRES iteration(s)"
        << std::endl;

  solution = solution_owned;
}

void
NavierStokes3D::output()
{
  pcout << "\n===============================================" << std::endl;

  DataOut<dim> data_out;

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    interpretation(dim,
                   DataComponentInterpretation::component_is_part_of_vector);
  interpretation.push_back(DataComponentInterpretation::component_is_scalar);

  std::vector<std::string> names(dim, "velocity");
  names.push_back("pressure");

  data_out.add_data_vector(dof_handler, solution, names, interpretation);

  std::vector<unsigned int> partition_int(mesh.n_active_cells());
  GridTools::get_subdomain_association(mesh, partition_int);
  const Vector<double> partitioning(partition_int.begin(), partition_int.end());
  data_out.add_data_vector(partitioning, "partitioning");

  data_out.build_patches();

  const std::string output_file_name = "output-navier-stokes-3D";
  data_out.write_vtu_with_pvtu_record(/* folder = */ "./",
                                      /* basename = */ output_file_name,
                                      /* index = */ timestep_number,
                                      MPI_COMM_WORLD);

  pcout << "Output written to " << output_file_name << std::endl;
  pcout << "===============================================\n" << std::endl;
}

void
NavierStokes3D::run()
{
  // Setup initial conditions.
  {
    setup();

   Functions::ZeroFunction<dim> zero_initial(dim + 1);
   VectorTools::interpolate(dof_handler, zero_initial, solution_owned);
   solution = solution_owned;

    time            = 0.0;
    timestep_number = 0;

    // Output initial condition.
    output();
  }

  pcout << "\n===============================================" << std::endl;

  // Vector to store the wall-clock time taken by each timestep.
  std::vector<double> timestep_wall_times;

  // Time-stepping loop.
  while (time < T - 0.5 * delta_t)  // Avoid round-off errors.
    {
      Timer timestep_timer;

      time += delta_t;
      ++timestep_number;

      pcout << "Timestep " << std::setw(3) << timestep_number
            << ", time = " << std::setw(6) << std::fixed << std::setprecision(3)
            << time << ": " << std::endl;

      assemble();
      solve();

      // Perform parallel communication to update the ghost values of the
      // solution vector.
      solution = solution_owned;

        if (data.choice == 2 || data.choice == 3)
      {
        const auto [drag, lift] = compute_drag_lift_forces();
        const double U   = reference_velocity();
        const double fac = 2.0 / (data.rho * U * U * D_cylinder * H_channel);

        const double cD = fac * drag;
        const double cL = fac * lift;
        const double dP = compute_pressure_difference();

        cD_history.push_back(cD);
        cL_history.push_back(cL);
        dP_history.push_back(dP);
        time_history.push_back(time);

        if (data.choice == 2)
        {
            if (time > 4.0) // skip transient phase (first 4 seconds)
            {
                cD_max = std::max(cD_max,std::abs(cD));
                cL_max = std::max(cL_max,std::abs(cL));
            }
        }

        if (data.choice == 3)
        {
            cD_max = std::max(cD_max,std::abs(cD));
            cL_max = std::max(cL_max,std::abs(cL));
        }
      }

      if (timestep_number % 10 == 0)
      { 
        output();
      }

      // Stop the timer and take the maximum wall time across all MPI ranks
      timestep_timer.stop();
      const double timestep_wall_time =
        Utilities::MPI::max(timestep_timer.wall_time(), MPI_COMM_WORLD);
  
      timestep_wall_times.push_back(timestep_wall_time);
  
      pcout << "  Timestep wall time: " << std::fixed << std::setprecision(4)
            << timestep_wall_time << " s\n" << std::endl;

    }
  
  // Export the per-timestep wall time to a CSV file, so it can be
  // plotted or compared across different preconditioners/mesh sizes.
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      std::ofstream timing_out("timestep_timings.csv");
      timing_out << "timestep,wall_time_s\n";
      for (unsigned int i = 0; i < timestep_wall_times.size(); ++i)
        timing_out << (i + 1) << "," << timestep_wall_times[i] << "\n";
    }
 
  // Computation
  const BenchmarkResult res = compute_benchmark_result();

  // Write csv
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
  {
    {
      std::ofstream out("benchmark_quantities.csv");
      if (data.choice == 1)
      {
        out << "Re,cD,cL,DeltaP\n"
            << res.Re << "," << res.cD << "," << res.cL << ","
            << res.dP << "\n";
      }
      else if (data.choice == 2)
      {
        out << "Re,cDmax,cLmax,St\n"
            << res.Re << "," << res.cD_max << "," << res.cL_max << ","
            << res.St << "\n";
      }
      else if (data.choice == 3)
      {
        out << "cDmax,cLmax,DeltaP_t8\n"
            << res.cD_max << "," << res.cL_max << "," << res.dP << "\n";
      }
    }

    // Time series 
    if (data.choice == 2)
    {
      std::ofstream out("time_series_3D2.csv");
      out << "time,cD,cL\n";
      for (unsigned int i = 0; i < time_history.size(); ++i)
        out << time_history[i] << "," << cD_history[i] << ","
            << cL_history[i] << "\n";
    }
    if (data.choice == 3)
    {
      std::ofstream out("time_series_3D3.csv");
      out << "time,cD,cL,deltaP\n";
      for (unsigned int i = 0; i < time_history.size(); ++i)
        out << time_history[i] << "," << cD_history[i] << ","
            << cL_history[i] << "," << dP_history[i] << "\n";
    }
  }

  print_benchmark_quantities(res);

}


std::pair<double, double>
NavierStokes3D::compute_drag_lift_forces() const
{

  // FEFaceValues object for integration over boundary faces
  // (updates values, gradients, outward normals, and quadrature weights JxW)
   FEFaceValues<dim> fe_face_values(*fe,
                                   *quadrature_face,
                                   update_values | update_gradients |
                                     update_normal_vectors |
                                     update_JxW_values);
  // Extractors to access velocity (dim components) and pressure from the FE solution
  FEValuesExtractors::Vector velocity(0);
  FEValuesExtractors::Scalar pressure(dim);

  const unsigned int n_q_face = quadrature_face->size();

  // Buffers for solution values at face quadrature points
  //std::vector<Tensor<1, dim>> velocity_values(n_q_face); // u at quadrature point q
  std::vector<Tensor<2, dim>> velocity_grads(n_q_face); // ∇u at quadrature point q
  std::vector<double> pressure_values(n_q_face); // p at quadrature point q

  // Local accumulator for the drag and lift force
  double drag = 0.0;
  double lift = 0.0;
  double side = 0.0;  // z component to check symmetry

  for (const auto &cell : dof_handler.active_cell_iterators())
    {

      // skip cells not owned by this MPI process
      if (!cell->is_locally_owned()) 
        continue;  

      // skip cells that do not touch the boundary
      if (!cell->at_boundary())
        continue; 

      for (unsigned int f = 0; f < cell->n_faces(); ++f)
        {

          // Only consider faces on the boundary with ID=4 (cylinder surface)
          if (!(cell->face(f)->at_boundary() && cell->face(f)->boundary_id() == 4))
            continue;

          // initialize FEFaceValues on the current face
          fe_face_values.reinit(cell, f);

          // Extract velocity values, gradients, and pressure from the solution vector
          //fe_face_values[velocity].get_function_values(solution, velocity_values);
          fe_face_values[velocity].get_function_gradients(solution, velocity_grads);
          fe_face_values[pressure].get_function_values(solution, pressure_values);

               
         for (unsigned int q = 0; q < n_q_face; ++q)
            {

              // The outward normal from the fluid domain points into the cylinder,
              // so we flip it to get the normal pointing outward from the cylinder into the fluid.
              const Tensor<1, dim> n = -fe_face_values.normal_vector(q);

              // Pressure at the quadrature point
              const double p = pressure_values[q];

              // Gradient of velocity at the quadrature point
              const Tensor<2, dim> &grad_u = velocity_grads[q];

              // Traction: t = sigma . n = -p n + rho*nu*(grad_u + grad_u^T) n
              Tensor<1, dim> traction;
              for (unsigned int i = 0; i < dim; ++i)
                {
                  traction[i] = -p * n[i];
                  for (unsigned int j = 0; j < dim; ++j)
                    traction[i] += data.rho * data.nu *
                                   (grad_u[i][j] + grad_u[j][i]) * n[j];
                }

              drag += traction[0] * fe_face_values.JxW(q);
              lift += traction[1] * fe_face_values.JxW(q);
              side += traction[2] * fe_face_values.JxW(q);
            }

        }

    }

  // Sum contributions across MPI tasks.
  const double drag_global = Utilities::MPI::sum(drag, MPI_COMM_WORLD);
  const double lift_global = Utilities::MPI::sum(lift, MPI_COMM_WORLD);
  const double side_global = Utilities::MPI::sum(side, MPI_COMM_WORLD);

  // We can ignore the side force, but we can check that it is small (ideally zero) to verify symmetry.
  // pcout << "  [check] side force (z) = " << side_global << std::endl;
  (void) side_global;

  return {drag_global, lift_global};
}


double
NavierStokes3D::compute_pressure_difference() const
{
  const FE_SimplexP<dim> fe_map(1);
  const MappingFE<dim>   mapping(fe_map);

  // Extractor for the pressure component
  const FEValuesExtractors::Scalar pressure(dim);

  // Lambda: returns the pressure value at a given geometric point p
  auto get_pressure_at_point = [&](const Point<dim> &p) -> double
  {
    double local_value = 0.0;

    try
      {
        Vector<double> values(fe->n_components());
        // Evaluates the solution at point p; throws if the point is not on this MPI process
        VectorTools::point_value(mapping, dof_handler, solution, p, values);
        // extract only the pressure component (index dim)
        local_value = values(dim);
      }
    catch (const VectorTools::ExcPointNotAvailableHere &)
      {
          // The point belongs to another MPI process: local_value stays 0
      }

    // MPI reduction: sum across processes (only one will have local_value ≠ 0)
    return Utilities::MPI::sum(local_value, MPI_COMM_WORLD);

  };

  // Sampling point at the front of the cylinder 
  const Point<dim> p_front(x_front_cylinder, y_probe, z_probe);
  // Sampling point at the rear of the cylinder 
  const Point<dim> p_back(x_back_cylinder, y_probe, z_probe);

  // pressure upstream and downstream of the cylinder
  const double p1 = get_pressure_at_point(p_front);
  const double p2 = get_pressure_at_point(p_back);

  // ΔP = p_front − p_back (positive when there is a pressure drop)
  return p1 - p2;
}


/*
double NavierStokes3D::compute_recirculation_length() const 
{ 
  const FE_SimplexP<dim> fe_map(1);  
  const MappingFE<dim>   mapping(fe_map);
  const unsigned int n_samples = 400; 
  // fixed y-coordinate (cylinder centerline)
  const double y = y_probe; 
  const double z = z_probe;

  // Lambda: returns the x-component of velocity at point (x, y, z)
  auto u_x_at = [&](const double x) -> double { 
    double local_value = 0.0;
    try
      {
        Vector<double> values(fe->n_components()); 
        const Point<dim> p(x, y, z); 
        VectorTools::point_value(mapping, dof_handler, solution, p, values); 
        local_value = values[0]; // x-component u_x
      }
    catch (const VectorTools::ExcPointNotAvailableHere &)
      {
        // Point not on this MPI process
      }
    return Utilities::MPI::sum(local_value, MPI_COMM_WORLD);
  }; 

  // First sampling point: just past the rear edge of the cylinder (small offset to avoid boundary)
  double x_prev = x_wake_start + 1e-6; 
  double u_prev = u_x_at(x_prev); 

  for (unsigned int k = 1; k <= n_samples; ++k) { 
    // Sample uniformly from x_back_cylinder+0.005 to x_wake_end
    const double x =
        (x_back_cylinder + 0.005) +
        (x_wake_end - (x_back_cylinder + 0.005)) * static_cast<double>(k) / n_samples;
    const double u = u_x_at(x);

    // Detect sign change of u_x: u_prev < 0 and u >= 0
    // -> the recirculation zone (u_x < 0) ends here
    if (u_prev < 0.0 && u >= 0.0) 
    { 
      // Linear interpolation to find the exact point xr where u_x = 0
      const double xr = x_prev - u_prev * (x - x_prev) / (u - u_prev); 
      return xr - x_back_cylinder; 
    } 

    // advance to the next sampling point
    x_prev = x; 
    u_prev = u; 
  } 
  
  // No sign change found: La is undefined (steady flow with no recirculation bubble)
  return std::numeric_limits<double>::quiet_NaN();
}
*/


double
NavierStokes3D::compute_reynolds_number() const
{
  const double U = reference_velocity();

  // Re = U · D / ν   (Reynolds number based on cylinder diameter)
  return U * D_cylinder / data.nu;
}

double
NavierStokes3D::compute_strouhal_number() const
{
  // need at least a time history to work with
  if (cL_history.size() < 2)
    return std::numeric_limits<double>::quiet_NaN();

  // Determine the raw max and min of cL after the transient phase, and
  // figure out which extremum has the larger magnitude. In 3D
  // (circular cylinder), cL is nearly symmetric around zero and may be
  // negatively biased.
  double raw_max = -std::numeric_limits<double>::max();
  double raw_min =  std::numeric_limits<double>::max();

  for (unsigned int i = 0; i < cL_history.size(); ++i)
    {
      if (time_history[i] < 4.0) continue; // skip transient phase (first 4 seconds)
      raw_max = std::max(raw_max, cL_history[i]);
      raw_min = std::min(raw_min, cL_history[i]);
    }

  const bool use_maxima = std::abs(raw_max) >= std::abs(raw_min);
  const double extreme_val = use_maxima ? raw_max : -raw_min;

  if (extreme_val <= 0.0)
    return std::numeric_limits<double>::quiet_NaN();

  // Threshold for peak detection: 50% of the dominant extremum's magnitude.
  const double threshold = 0.5 * extreme_val;

  // Minimum time gap between consecutive peaks to avoid detecting small oscillations
  const double min_time_gap = 0.15; // seconds

  std::vector<double> peak_times;
  double last_peak_time = -min_time_gap; // initialize to a time before the start

  for (unsigned int i = 1; i + 1 < cL_history.size(); ++i)
    {
      if (time_history[i] < 4.0) continue; // skip transient phase (first 4 seconds)
      
      bool is_peak = false;
      if (use_maxima)
        {
          // Look for local maxima above +threshold.
          is_peak = (cL_history[i] >= threshold) &&
                    (cL_history[i] > cL_history[i - 1]) &&
                    (cL_history[i] > cL_history[i + 1]);
        }
      else
        {
          // Look for local minima below -threshold.
          is_peak = (cL_history[i] <= -threshold) &&
                    (cL_history[i] < cL_history[i - 1]) &&
                    (cL_history[i] < cL_history[i + 1]);
        }

      if (!is_peak) continue;
      if (time_history[i] - last_peak_time < min_time_gap) continue;
      
      peak_times.push_back(time_history[i]);
      last_peak_time = time_history[i];
    }

  // need at least 2 peaks to estimate the period
  if (peak_times.size() < 2)
    return std::numeric_limits<double>::quiet_NaN();

  // Estimate the mean period T as the average gap between consecutive peaks
  double period = 0.0;
  for (unsigned int i = 1; i < peak_times.size(); ++i)
    period += peak_times[i] - peak_times[i - 1];
  period /= static_cast<double>(peak_times.size() - 1);

  // vortex shedding frequency: f = 1/T
  const double f  = 1.0 / period;
  const double U  = reference_velocity();

  // St = f · D / U   (Strouhal number)
  return D_cylinder * f / U;
}


NavierStokes3D::BenchmarkResult
NavierStokes3D::compute_benchmark_result() const
{
  BenchmarkResult r;
  r.Re = compute_reynolds_number();

  if (data.choice == 1)
  {
    auto [drag, lift] = compute_drag_lift_forces();
    const double U    = reference_velocity();
    const double fac  = 2.0 / (data.rho * U * U * D_cylinder * H_channel);  

    r.cD = fac * drag;
    r.cL = fac * lift;
    r.dP = compute_pressure_difference();
  }
  else if (data.choice == 2)
  {
    r.cD_max  = cD_max;
    r.cL_max  = cL_max;
    r.St      = compute_strouhal_number();
  }
  else if (data.choice == 3)
  {
    r.cD_max = cD_max;
    r.cL_max = cL_max;
    r.dP     = compute_pressure_difference(); // at t=8
  }
  return r;
}


void
NavierStokes3D::print_benchmark_quantities(const BenchmarkResult &res) const
{
  pcout << std::fixed << std::setprecision(4);
  pcout << "Benchmark quantities at t = " << time << '\n';
  pcout << "  Re = " << res.Re << '\n';

  if (data.choice == 1)
  {
    pcout << "  cD = " << res.cD << '\n';
    pcout << "  cL = " << res.cL << '\n';
    pcout << "  ΔP = " << res.dP << '\n';
  }
  else if (data.choice == 2)
  {
    pcout << "  cD_max     = " << res.cD_max  << '\n';
    pcout << "  cL_max     = " << res.cL_max  << '\n';
    pcout << "  St         = " << res.St      << '\n';
  }
  else if (data.choice == 3)
  {
    pcout << "  cD_max    = " << res.cD_max << '\n';
    pcout << "  cL_max    = " << res.cL_max << '\n';
    pcout << "  ΔP(t=8)   = " << res.dP    << '\n';
  }
} 

