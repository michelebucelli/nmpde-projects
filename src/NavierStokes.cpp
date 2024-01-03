#include "NavierStokes.hpp"

// Static members have to be defined in the .cpp file.
// Forcing term.
NavierStokes::ForcingTerm NavierStokes::forcing_term;

// Inlet velocity.
NavierStokes::InletVelocity NavierStokes::inlet_velocity;

// Initial conditions.
NavierStokes::InitialConditions NavierStokes::initial_conditions;

// Reynolds number.
NavierStokes::ReynoldsNumber NavierStokes::reynolds_number;


//Current issues:
//-nothing was tested
//-MPI communication has to be checked
//-the mass term in the right hand side is missing (refer to assemble_time_dependent)
//-the convection term in the matrix is missing (refer to assemble_time_dependent)
//-the initial condition is not applied


void
NavierStokes::setup()
{
  // Create the mesh.
  {
    pcout << "Initializing the mesh" << std::endl;

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
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

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

    // We also build a sparsity pattern for the mass matrix.
    for (unsigned int c = 0; c < dim + 1; ++c)
      {
        for (unsigned int d = 0; d < dim + 1; ++d)
          {
            if (c != dim && d != dim) // velocity-velocity term
              coupling[c][d] = DoFTools::always;
            else // other combinations
              coupling[c][d] = DoFTools::none;
          }
      }
    TrilinosWrappers::BlockSparsityPattern sparsity_mass(
      block_owned_dofs, MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler,
                                    coupling,
                                    sparsity_mass);
    sparsity_mass.compress();

    pcout << "  Initializing the matrices" << std::endl;
    constant_lhs_matrix.reinit(sparsity);
    lhs_matrix.reinit(sparsity);
    mass_matrix.reinit(sparsity_mass);
    pressure_mass.reinit(sparsity_pressure_mass);

    pcout << "  Initializing the system right-hand side" << std::endl;
    system_rhs.reinit(block_owned_dofs, MPI_COMM_WORLD);
    pcout << "  Initializing the solution vector" << std::endl;
    solution_owned.reinit(block_owned_dofs, MPI_COMM_WORLD);
    solution.reinit(block_owned_dofs, block_relevant_dofs, MPI_COMM_WORLD);
  }
}

void
NavierStokes::assemble_constant_matrices()
{
  pcout << "===============================================" << std::endl;
  pcout << "Assembling the system" << std::endl;

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

  FullMatrix<double> cell_lhs_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_mass_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_pressure_matrix(dofs_per_cell, dofs_per_cell);

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  constant_lhs_matrix = 0.0;
  mass_matrix = 0.0;
  pressure_mass = 0.0;

  FEValuesExtractors::Vector velocity(0);
  FEValuesExtractors::Scalar pressure(dim);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);

      cell_lhs_matrix = 0.0;
      cell_mass_matrix = 0.0;
      cell_pressure_matrix = 0.0;

      for (unsigned int q = 0; q < n_q; ++q)
        {
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                  // Mass term.
                  cell_mass_matrix(i, j) +=
                    scalar_product(fe_values[velocity].value(i, q),
                                   fe_values[velocity].value(j, q)) *
                    fe_values.JxW(q) / deltat;
                  cell_lhs_matrix(i, j) += cell_mass_matrix(i, j);

                  // Viscosity term.
                  cell_lhs_matrix(i, j) +=
                    nu *
                    scalar_product(fe_values[velocity].gradient(i, q),
                                   fe_values[velocity].gradient(j, q)) *
                    fe_values.JxW(q);

                  // Pressure term in the momentum equation.
                  cell_lhs_matrix(i, j) -= fe_values[velocity].divergence(i, q) *
                                       fe_values[pressure].value(j, q) *
                                       fe_values.JxW(q);

                  // Pressure term in the continuity equation.
                  cell_lhs_matrix(i, j) -= fe_values[velocity].divergence(j, q) *
                                       fe_values[pressure].value(i, q) *
                                       fe_values.JxW(q);

                  // Pressure mass matrix.
                  cell_pressure_matrix(i, j) +=
                    fe_values[pressure].value(i, q) *
                    fe_values[pressure].value(j, q) / nu * fe_values.JxW(q);
                }
            }
        }

      cell->get_dof_indices(dof_indices);

      constant_lhs_matrix.add(dof_indices, cell_lhs_matrix);
      mass_matrix.add(dof_indices, cell_mass_matrix);
      pressure_mass.add(dof_indices, cell_pressure_matrix);
    }

  constant_lhs_matrix.compress(VectorOperation::add);
  mass_matrix.compress(VectorOperation::add);
  pressure_mass.compress(VectorOperation::add);

  // Dirichlet boundary conditions.
  // If time dependent boundary conditions are used instead, this has to be moved to the assemble_rhs() method.
  // Assuming inlet has boundary tag 0, outlet has tag 1 and other boundaries have tags 2 to 9
  {
    std::map<types::global_dof_index, double>           boundary_values;
    std::map<types::boundary_id, const Function<dim> *> boundary_functions;

    // We interpolate first the inlet velocity condition alone, then the wall
    // condition alone, so that the latter "win" over the former where the two
    // boundaries touch.
    boundary_functions[0] = &inlet_velocity;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             boundary_functions,
                                             boundary_values,
                                             ComponentMask(
                                               {true, true, true, false}));

    boundary_functions.clear();
    Functions::ZeroFunction<dim> zero_function(dim + 1);
    for (unsigned int i = 2; i < 10; i++)
      boundary_functions[i] = &zero_function;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             boundary_functions,
                                             boundary_values,
                                             ComponentMask(
                                               {true, true, true, false}));
  }
}

void
NavierStokes::assemble_time_dependent(const double &time_step)
{
  pcout << "===============================================" << std::endl;
  pcout << "Assembling the system" << std::endl;

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

  Vector<double>     cell_rhs(dofs_per_cell);
  FullMatrix<double> cell_lhs_matrix(dofs_per_cell, dofs_per_cell);

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  system_rhs    = 0.0;
  lhs_matrix = constant_lhs_matrix;

  FEValuesExtractors::Vector velocity(0);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);

      cell_rhs = 0.0;
      cell_lhs_matrix = 0.0;

      for (unsigned int q = 0; q < n_q; ++q)
        {
          Vector<double> forcing_term_loc(dim);
          forcing_term.vector_value(fe_values.quadrature_point(q),
                                    forcing_term_loc);
          Tensor<1, dim> forcing_term_tensor;
          for (unsigned int d = 0; d < dim; ++d)
            forcing_term_tensor[d] = forcing_term_loc[d];

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                  // Convection term. Other types of disrectization are possible. 
                  Tensor<1, dim> convection_term;
                  // convection_term needs the value of the latest solution to be computable.
                  cell_lhs_matrix(i, j) += scalar_product(
                                       convection_term,
                                       fe_values[velocity].value(i, q)) *
                                       fe_values.JxW(q);
                }

              // Forcing term.
              cell_rhs(i) += scalar_product(forcing_term_tensor,
                                            fe_values[velocity].value(i, q)) *
                             fe_values.JxW(q);

              // Mass matrix.
              Tensor<1, dim> mass_term;
              // mass_term should be the product of the mass matrix and the latest velocity solution. Each device/CPU needs to have the entire solution vector.
              //mass_matrix.vmult(mass_term, solution);
              cell_rhs(i) += scalar_product(mass_term,
                             fe_values[velocity].value(i, q)) *
                             fe_values.JxW(q);
            }
        }

      // Boundary integral for Neumann BCs.
      // Assuming tag for outlet boundary is 1.
      if (cell->at_boundary())
        {
          for (unsigned int f = 0; f < cell->n_faces(); ++f)
            {
              if (cell->face(f)->at_boundary() &&
                  cell->face(f)->boundary_id() == 1)
                {
                  fe_face_values.reinit(cell, f);

                  for (unsigned int q = 0; q < n_q_face; ++q)
                    {
                      for (unsigned int i = 0; i < dofs_per_cell; ++i)
                        {
                          cell_rhs(i) +=
                            -p_out *
                            scalar_product(fe_face_values.normal_vector(q),
                                           fe_face_values[velocity].value(i,
                                                                          q)) *
                            fe_face_values.JxW(q);
                        }
                    }
                }
            }
        }

      cell->get_dof_indices(dof_indices);

      system_rhs.add(dof_indices, cell_rhs);
      lhs_matrix.add(dof_indices, cell_lhs_matrix);
    }

  system_rhs.compress(VectorOperation::add);
  lhs_matrix.compress(VectorOperation::add);
}

void
NavierStokes::solve_time_step()
{
  pcout << "===============================================" << std::endl;

  SolverControl solver_control(2000, 1e-6 * system_rhs.l2_norm());

  SolverGMRES<TrilinosWrappers::MPI::BlockVector> solver(solver_control);

  PreconditionIdentity preconditioner;

  pcout << "Solving the linear system" << std::endl;
  solver.solve(lhs_matrix, solution_owned, system_rhs, preconditioner);
  pcout << "  " << solver_control.last_step() << " GMRES iterations"
        << std::endl;

  solution = solution_owned;
}

void
NavierStokes::solve()
{
  pcout << "===============================================" << std::endl;

  time = 0.0;

  // Apply the initial condition.
  {
    // TODO: apply the real initial condition.
    pcout << "Applying the initial condition" << std::endl;

    solution = solution_owned;

    // Output the initial solution.
    output(0);
    pcout << "-----------------------------------------------" << std::endl;
  }

  unsigned int time_step = 0;

  while (time < T - 0.5 * deltat)
    {
      time += deltat;
      ++time_step;

      pcout << "n = " << std::setw(3) << time_step << ", t = " << std::setw(5)
            << time << ":" << std::flush;

      assemble_time_dependent(time);
      solve_time_step();
      output(time_step);
    }
}

void
NavierStokes::output(const unsigned int &time_step) const
{
  pcout << "===============================================" << std::endl;

  DataOut<dim> data_out;

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(
      dim, DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(
    DataComponentInterpretation::component_is_scalar);
  std::vector<std::string> names = {"velocity",
                                    "velocity",
                                    "velocity",
                                    "pressure"};

  data_out.add_data_vector(dof_handler,
                           solution,
                           names,
                           data_component_interpretation);

  std::vector<unsigned int> partition_int(mesh.n_active_cells());
  GridTools::get_subdomain_association(mesh, partition_int);
  const Vector<double> partitioning(partition_int.begin(), partition_int.end());
  data_out.add_data_vector(partitioning, "partitioning");

  data_out.build_patches();

  const std::string output_file_name = "output-navier-stokes";
  data_out.write_vtu_with_pvtu_record("./",
                                      output_file_name,
                                      0,
                                      MPI_COMM_WORLD);

  pcout << "Output written to " << output_file_name << std::endl;
  pcout << "===============================================" << std::endl;
}
