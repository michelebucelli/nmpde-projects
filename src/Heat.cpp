#include "Heat.hpp"

#include <array>

void
Heat::setup_system()
{
  if (!fe)
    {
      fe         = std::make_unique<FE_Q<dim>>(r);
      quadrature = std::make_unique<QGauss<dim>>(r + 1);
      dof_handler.reinit(mesh);
    }

  dof_handler.distribute_dofs(*fe);

  const IndexSet locally_owned_dofs    = dof_handler.locally_owned_dofs();
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
  solution_old.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
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
      fe_values.get_function_values(solution_old, solution_old_values);
      fe_values.get_function_gradients(solution_old, solution_old_grads);

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
  preconditioner.initialize(system_matrix,
                             TrilinosWrappers::PreconditionSSOR::AdditionalData(1.0));

  ReductionControl solver_control(/* maxiter = */ 10000,
                                   /* tolerance = */ 1.0e-12,
                                   /* reduce = */ 1.0e-6);

  SolverCG<TrilinosWrappers::MPI::Vector> solver(solver_control);

  solver.solve(system_matrix, solution_owned, system_rhs, preconditioner);

  pcout << solver_control.last_step() << " CG iterations" << std::endl;
}

void
Heat::refine_mesh()
{
  pcout << "  Estimating error and refining..." << std::endl;

  Vector<float> estimated_error_per_cell(mesh.n_active_cells());
  KellyErrorEstimator<dim>::estimate(dof_handler,
                                      QGauss<dim - 1>(r + 1),
                                      {},
                                      solution,
                                      estimated_error_per_cell);

  parallel::distributed::GridRefinement::refine_and_coarsen_fixed_fraction(
    mesh, estimated_error_per_cell, 0.6, 0.4);

  if (mesh.n_levels() > max_grid_level)
    for (const auto &cell : mesh.active_cell_iterators_on_level(max_grid_level))
      cell->clear_refine_flag();

  mesh.prepare_coarsening_and_refinement();

  parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::Vector> sol_trans(
    dof_handler);

  std::vector<const TrilinosWrappers::MPI::Vector *> transfer_in = {&solution, &solution_old};
  sol_trans.prepare_for_coarsening_and_refinement(transfer_in);

  mesh.execute_coarsening_and_refinement();

  // Ricostruzione sistema (questo riempirà solution_owned di zeri)
  setup_system();

  TrilinosWrappers::MPI::Vector interpolated_new(solution_owned);
  TrilinosWrappers::MPI::Vector interpolated_old(solution_owned_old);

  std::vector<TrilinosWrappers::MPI::Vector *> transfer_out = {&interpolated_new,
                                                                &interpolated_old};
  sol_trans.interpolate(transfer_out);

  solution_owned     = interpolated_new;
  solution_owned_old = interpolated_old;
  solution           = solution_owned;
  solution_old       = solution_owned_old;
}

void
Heat::output() const
{
  DataOut<dim> data_out;

  data_out.set_flags(DataOutBase::VtkFlags(time, timestep_number));
  data_out.add_data_vector(dof_handler, solution, "solution");

  // Add vector for parallel partition.
  std::vector<unsigned int> partition_int(mesh.n_active_cells());
  GridTools::get_subdomain_association(mesh, partition_int);
  const Vector<double> partitioning(partition_int.begin(), partition_int.end());
  data_out.add_data_vector(partitioning, "partitioning");

  data_out.build_patches();

  const std::string output_file_name = "output-" + run_basename();

  data_out.write_vtu_with_pvtu_record(/* folder = */ "./",
                                      /* basename = */ output_file_name,
                                      /* index = */ timestep_number,
                                      MPI_COMM_WORLD);
}

double
Heat::compute_time_error_indicator() const
{
  const unsigned int n_q = quadrature->size();

  FEValues<dim> fe_values(*fe,
                          *quadrature,
                          update_gradients | update_quadrature_points | update_JxW_values);

  std::vector<Tensor<1, dim>> grad_new(n_q);
  std::vector<Tensor<1, dim>> grad_old(n_q);

  const std::array<double, 3> gauss_pts = {{0.1127016653792583, 0.5, 0.8872983346207417}};
  const std::array<double, 3> gauss_wts = {{5.0 / 18.0, 8.0 / 18.0, 5.0 / 18.0}};

  double eps_sq   = 0.0;
  double gamma_sq = 0.0;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);
      fe_values.get_function_gradients(solution, grad_new);
      fe_values.get_function_gradients(solution_old, grad_old);

      for (unsigned int q = 0; q < n_q; ++q)
        {
          const Tensor<1, dim> diff_grad = grad_new[q] - grad_old[q];
          eps_sq += delta_t * diff_grad.norm_square() * fe_values.JxW(q);

          const Point<dim> &x_q         = fe_values.quadrature_point(q);
          const double       f_new_here = f(x_q, time + delta_t);

          double gamma_q = 0.0;
          for (unsigned int g = 0; g < gauss_pts.size(); ++g)
            {
              const double t_g  = time + gauss_pts[g] * delta_t;
              const double diff = f(x_q, t_g) - f_new_here;
              gamma_q += gauss_wts[g] * diff * diff;
            }
          gamma_sq += gamma_q * delta_t * fe_values.JxW(q);
        }
    }

  eps_sq   = Utilities::MPI::sum(eps_sq, MPI_COMM_WORLD);
  gamma_sq = Utilities::MPI::sum(gamma_sq, MPI_COMM_WORLD);

  return std::sqrt(eps_sq + gamma_sq);
}

double
Heat::compute_forcing_integral(const double t) const
{
  FEValues<dim> fe_values(*fe, *quadrature, update_quadrature_points | update_JxW_values);

  double integral = 0.0;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);
      for (unsigned int q = 0; q < quadrature->size(); ++q)
        integral += f(fe_values.quadrature_point(q), t) * fe_values.JxW(q);
    }

  return Utilities::MPI::sum(integral, MPI_COMM_WORLD);
}

double
Heat::compute_mass() const
{
  return VectorTools::compute_mean_value(dof_handler, *quadrature, solution, 0);
}

void
Heat::run()
{
  GridGenerator::subdivided_hyper_cube(mesh, n_elements, 0.0, 1.0);
  setup_system();

  VectorTools::interpolate(dof_handler, FunctionU0(), solution_owned);
  solution           = solution_owned;
  solution_owned_old = solution_owned;
  solution_old       = solution;

  std::ofstream        mass_balance_file;
  if (mpi_rank == 0)
    {
      mass_balance_file.open("mass_balance-" + run_basename() + ".csv");
      mass_balance_file << "computed_mass,expected_mass\n";
    }

  double expected_mass = 0.0;

  while (time < T - 1e-12)
    {
      delta_t = std::min(delta_t, T - time);

      if (adaptive)
        {
          bool   step_accepted = false;
          double error_t       = 0.0;

          while (!step_accepted)
            {
              assemble();
              solve_linear_system();
              solution = solution_owned;

              error_t = compute_time_error_indicator();

              if (error_t > tol_time_max && delta_t > dt_min)
                {
                  // Errore troppo grande: riduciamo dt e non avanziamo nel tempo
                  delta_t /= 2.0;
                  pcout << "  Step rejected! New delta_t = " << delta_t << std::endl;
                }
              else
                {
                  // Errore accettabile o dt minimo raggiunto: procediamo
                  step_accepted = true;

                  const double f_old_integral = compute_forcing_integral(time);
                  const double f_new_integral = compute_forcing_integral(time + delta_t);
                  expected_mass +=
                    delta_t * (theta * f_new_integral + (1.0 - theta) * f_old_integral);

                  time += delta_t;
                  timestep_number++;

                  if (error_t < tol_time_min && delta_t < dt_max)
                    delta_t *= 1.2;
                }
            }
        }
      else
        {
          // Non-adaptive baseline: always accept, delta_t stays fixed.
          assemble();
          solve_linear_system();
          solution = solution_owned;

          const double f_old_integral = compute_forcing_integral(time);
          const double f_new_integral = compute_forcing_integral(time + delta_t);
          expected_mass += delta_t * (theta * f_new_integral + (1.0 - theta) * f_old_integral);

          time += delta_t;
          timestep_number++;
        }

      // Dopo uno step accettato, gestiamo l'adattività spaziale
      if (adaptive && timestep_number % 5 == 0)
        refine_mesh();

      solution_owned_old = solution_owned;
      solution_old       = solution;
      output();
    }

  if (mpi_rank == 0)
    mass_balance_file << compute_mass() << "," << expected_mass << "\n";
}

double
Heat::compute_error(const VectorTools::NormType &norm_type,
                     const Function<dim>         &exact_solution) const
{
  const QGauss<dim> quadrature_error(r + 2);

  Vector<double> error_per_cell(mesh.n_active_cells());
  VectorTools::integrate_difference(
    dof_handler, solution, exact_solution, error_per_cell, quadrature_error, norm_type);

  return VectorTools::compute_global_error(mesh, error_per_cell, norm_type);
}
