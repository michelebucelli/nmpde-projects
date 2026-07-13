#include "H_Heat.hpp"

void
Heat::setup()
{

  // Create the mesh.
  {

    GridGenerator::hyper_cube(mesh);
    //refining level
    mesh.refine_global(5);

  }


  // Initialize the finite element space.
  {

    fe = std::make_unique<FE_Q<dim>>(r);

    quadrature = std::make_unique<QGauss<dim>>(r + 1);

  }


  // Initialize the DoF handler.
  {

    dof_handler.reinit(mesh);
    dof_handler.distribute_dofs(*fe);

  }

  // Initialize the linear system.
  {
    const IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
    const IndexSet locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);

    TrilinosWrappers::SparsityPattern sparsity(locally_owned_dofs,
                                               MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler, sparsity);
    sparsity.compress();

    system_matrix.reinit(sparsity);

    system_rhs.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    solution_owned.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    solution.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
  }
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

          const double f_old_loc =
            f(fe_values.quadrature_point(q), time - delta_t);
          const double f_new_loc = f(fe_values.quadrature_point(q), time);

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
                  cell_matrix(i, j) +=
                    theta * mu_loc *                             //
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

  // Homogeneous Neumann boundary conditions: we do nothing.
}

void
Heat::solve_linear_system()
{
  TrilinosWrappers::PreconditionSSOR preconditioner;
  preconditioner.initialize(
    system_matrix, TrilinosWrappers::PreconditionSSOR::AdditionalData(1.0));

  ReductionControl solver_control(/* maxiter = */ 10000,
                                  /* tolerance = */ 1.0e-16,
                                  /* reduce = */ 1.0e-6);

  SolverCG<TrilinosWrappers::MPI::Vector> solver(solver_control);

  solver.solve(system_matrix, solution_owned, system_rhs, preconditioner);
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

  const std::filesystem::path mesh_path(output_file_name);
  const std::string output_file_name = "output-" + mesh_path.stem().string();

  data_out.write_vtu_with_pvtu_record(/* folder = */ "./",
                                      /* basename = */ output_file_name,
                                      /* index = */ timestep_number,
                                      MPI_COMM_WORLD);
}

void
Heat::run()
{
  // Setup initial conditions.
  {
    profiler.tic("setup");
    setup();
    profiler.toc("setup");

    VectorTools::interpolate(dof_handler, FunctionU0(), solution_owned);
    solution = solution_owned;

    time            = 0.0;
    timestep_number = 0;

    // Output initial condition.
    output();
  }

  // Time-stepping loop.
  while (time < T - 0.5 * delta_t)
    {
      num_of_steps++;

      time += delta_t;
      ++timestep_number;


      profiler.tic("assemble");
      assemble();
      profiler.toc("assemble");

      profiler.tic("solve");
      solve_linear_system();
      profiler.toc("solve");

      // Perform parallel communication to update the ghost values of the
      // solution vector.
      profiler.tic("update-ghost");
      solution = solution_owned;
      profiler.toc("update-ghost");

      profiler.tic("output");
      output();
      profiler.toc("output");
    }
}

const DoFHandler<Heat::dim> & 
Heat::get_dof_handler() const
{
  return dof_handler;
}

Vector<double> 
Heat::get_serial_solution() const
{
  Vector<double> serial;
  serial.reinit(solution.size());
  for(unsigned int i = 0; i < solution.size(); ++i)
    {
      serial(i) = solution(i);
    }
  return serial;
}

void Heat::print_results(){
  pcout << "Baseline Results" << std::endl;
  pcout << "Number of steps: " << num_of_steps << std::endl 
        << "Elapsed Time: " << std::endl;
        profiler.report(MPI_COMM_WORLD);
  pcout << "Number of DoFs: " << dof_handler.n_dofs() << std::endl;
}