#include "TimeADR.hpp"

#include <algorithm>
#include <chrono>
#include <limits>

namespace
{
  using Clock = std::chrono::steady_clock;

  double
  elapsed_seconds(const Clock::time_point &start, const Clock::time_point &end)
  {
    return std::chrono::duration_cast<std::chrono::duration<double>>(end -
                                                                     start)
      .count();
  }

  // Wraps the (space + time) Dirichlet datum g(x, t) into a deal.II
  // Function<dim> with a settable time, so it can feed
  // VectorTools::interpolate_boundary_values when building the lifting.
  class GFunction : public Function<TimeADR::dim>
  {
  public:
    explicit GFunction(
      const std::function<double(const Point<TimeADR::dim> &, double)> &g)
      : g_(g)
    {}

    double
    value(const Point<TimeADR::dim> &p,
          const unsigned int /*component*/ = 0) const override
    {
      return g_(p, this->get_time());
    }

  private:
    std::function<double(const Point<TimeADR::dim> &, double)> g_;
  };
} // namespace

void
TimeADR::setup()
{
  const auto t0 = Clock::now();

  pcout << "===============================================" << std::endl;

  // -------------------------------------------------------------------------
  // Mesh.
  // -------------------------------------------------------------------------
  {
    pcout << "Initializing the mesh" << std::endl;

    Triangulation<dim> mesh_serial;

    std::visit(
      [&](auto &&arg) {
        using T = std::decay_t<decltype(arg)>;

        if constexpr (std::is_same_v<T, unsigned int>)
          {
            GridGenerator::subdivided_hyper_cube_with_simplices(mesh_serial,
                                                                N_el,
                                                                0.0,
                                                                1.0);
          }
        else if constexpr (std::is_same_v<T, std::string>)
          {
            const std::string &mesh_file_name = arg;
            GridIn<dim>        grid_in;
            grid_in.attach_triangulation(mesh_serial);
            std::ifstream mesh_file(mesh_file_name);
            grid_in.read_msh(mesh_file);
          }
      },
      mesh_spec);

    GridTools::partition_triangulation(mpi_size, mesh_serial);
    const auto construction_data = TriangulationDescription::Utilities::
      create_description_from_triangulation(mesh_serial, MPI_COMM_WORLD);
    mesh.create_triangulation(construction_data);

    pcout << "  Number of elements = " << mesh.n_global_active_cells()
          << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // -------------------------------------------------------------------------
  // Finite element space and quadrature.
  // -------------------------------------------------------------------------
  {
    pcout << "Initializing the finite element space" << std::endl;

    fe = std::make_unique<FE_SimplexP<dim>>(r);

    pcout << "  Degree                     = " << fe->degree << std::endl;
    pcout << "  DoFs per cell              = " << fe->dofs_per_cell
          << std::endl;

    quadrature      = std::make_unique<QGaussSimplex<dim>>(r + 1);
    quadrature_face = std::make_unique<QGaussSimplex<dim - 1>>(fe->degree + 1);

    pcout << "  Quadrature points per cell = " << quadrature->size()
          << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // -------------------------------------------------------------------------
  // DoF handler.
  // -------------------------------------------------------------------------
  {
    pcout << "Initializing the DoF handler" << std::endl;

    dof_handler.reinit(mesh);
    dof_handler.distribute_dofs(*fe);

    pcout << "  Number of DoFs = " << dof_handler.n_dofs() << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // -------------------------------------------------------------------------
  // Constraints: homogeneous Dirichlet on the whole boundary.
  // (Constant in time, so we build them once here.)
  // -------------------------------------------------------------------------
  const IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
  const IndexSet locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);

  {
    pcout << "Initializing the constraints" << std::endl;

    constraints.clear();
    // NOTE: on very recent deal.II versions use the two-argument form
    //   constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
    constraints.reinit(locally_relevant_dofs);

    DoFTools::make_hanging_node_constraints(dof_handler, constraints);

    // Homogeneous Dirichlet on boundary ids 0..3 (covers both the
    // non-colorized case, where only id 0 exists, and the colorized one).
    for (const types::boundary_id id : {types::boundary_id(0),
                                        types::boundary_id(1),
                                        types::boundary_id(2),
                                        types::boundary_id(3)})
      VectorTools::interpolate_boundary_values(dof_handler,
                                               id,
                                               Functions::ZeroFunction<dim>(),
                                               constraints);

    constraints.close();
  }

  // Non-homogeneous case only: hanging-node-only constraints, used to assemble
  // the "full" matrices while KEEPING the boundary coupling.
  if (g_fun)
    {
      constraints_hn.clear();
      constraints_hn.reinit(locally_relevant_dofs);
      DoFTools::make_hanging_node_constraints(dof_handler, constraints_hn);
      constraints_hn.close();
    }

  pcout << "-----------------------------------------------" << std::endl;

  // -------------------------------------------------------------------------
  // Sparsity pattern, matrices and vectors.
  // -------------------------------------------------------------------------
  {
    pcout << "Initializing the linear system" << std::endl;

    TrilinosWrappers::SparsityPattern sparsity(locally_owned_dofs,
                                               MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler,
                                    sparsity,
                                    constraints,
                                    /* keep_constrained_dofs = */ false);
    sparsity.compress();

    mass_matrix.reinit(sparsity);
    stiffness_matrix.reinit(sparsity);
    lhs_matrix.reinit(sparsity);

    system_rhs.reinit(locally_owned_dofs, MPI_COMM_WORLD);

    u_old.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    u.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    u_new.reinit(locally_owned_dofs, MPI_COMM_WORLD);

    tmp.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    tmp2.reinit(locally_owned_dofs, MPI_COMM_WORLD);

    solution.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
  }

  // Non-homogeneous case only: allocate the "full" matrices with a sparsity
  // pattern that keeps the boundary coupling (constraints_hn, and we keep the
  // constrained DoFs).
  if (g_fun)
    {
      TrilinosWrappers::SparsityPattern sparsity_full(locally_owned_dofs,
                                                      MPI_COMM_WORLD);
      DoFTools::make_sparsity_pattern(dof_handler,
                                      sparsity_full,
                                      constraints_hn,
                                      /* keep_constrained_dofs = */ true);
      sparsity_full.compress();

      mass_full.reinit(sparsity_full);
      stiffness_full.reinit(sparsity_full);
    }

  pcout << "-----------------------------------------------" << std::endl;

  // -------------------------------------------------------------------------
  // Assemble the constant matrices and build the preconditioner.
  // -------------------------------------------------------------------------
  {
    pcout << "Assembling the constant matrices" << std::endl;

    assemble_matrices();
    if (g_fun)
      assemble_full_matrices();

    pcout << "Initializing the preconditioner" << std::endl;
    preconditioner.initialize(
      lhs_matrix, TrilinosWrappers::PreconditionSSOR::AdditionalData(1.0));
  }

  perf.setup_seconds += elapsed_seconds(t0, Clock::now());
}

void
TimeADR::assemble_matrices()
{
  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q           = quadrature->size();

  FEValues<dim> fe_values(*fe,
                          *quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  FullMatrix<double> cell_mass(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_stiffness(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_lhs(dofs_per_cell, dofs_per_cell);

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  mass_matrix      = 0.0;
  stiffness_matrix = 0.0;
  lhs_matrix       = 0.0;

  const double inv_dt2 = 1.0 / (delta_t * delta_t);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);

      cell_mass      = 0.0;
      cell_stiffness = 0.0;
      cell_lhs       = 0.0;

      for (unsigned int q = 0; q < n_q; ++q)
        {
          const double mu_loc = mu(fe_values.quadrature_point(q));
          const double JxW    = fe_values.JxW(q);

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              {
                const double m_ij = fe_values.shape_value(i, q) *
                                    fe_values.shape_value(j, q) * JxW;

                const double a_ij = mu_loc *
                                    scalar_product(fe_values.shape_grad(i, q),
                                                   fe_values.shape_grad(j, q)) *
                                    JxW;

                cell_mass(i, j) += m_ij;
                cell_stiffness(i, j) += a_ij;
                cell_lhs(i, j) += inv_dt2 * m_ij + theta * a_ij;
              }
        }

      cell->get_dof_indices(dof_indices);

      constraints.distribute_local_to_global(cell_mass,
                                             dof_indices,
                                             mass_matrix);
      constraints.distribute_local_to_global(cell_stiffness,
                                             dof_indices,
                                             stiffness_matrix);
      constraints.distribute_local_to_global(cell_lhs, dof_indices, lhs_matrix);
    }

  mass_matrix.compress(VectorOperation::add);
  stiffness_matrix.compress(VectorOperation::add);
  lhs_matrix.compress(VectorOperation::add);
}

// Same cell integrals as assemble_matrices(), but distributed with the
// hanging-node-only constraints so that the boundary coupling is retained.
// Used only in the non-homogeneous case, to compute the lifting correction.
void
TimeADR::assemble_full_matrices()
{
  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q           = quadrature->size();

  FEValues<dim> fe_values(*fe,
                          *quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  FullMatrix<double> cell_mass(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_stiffness(dofs_per_cell, dofs_per_cell);

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  mass_full      = 0.0;
  stiffness_full = 0.0;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);

      cell_mass      = 0.0;
      cell_stiffness = 0.0;

      for (unsigned int q = 0; q < n_q; ++q)
        {
          const double mu_loc = mu(fe_values.quadrature_point(q));
          const double JxW    = fe_values.JxW(q);

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              {
                cell_mass(i, j) += fe_values.shape_value(i, q) *
                                   fe_values.shape_value(j, q) * JxW;

                cell_stiffness(i, j) +=
                  mu_loc *
                  scalar_product(fe_values.shape_grad(i, q),
                                 fe_values.shape_grad(j, q)) *
                  JxW;
              }
        }

      cell->get_dof_indices(dof_indices);

      constraints_hn.distribute_local_to_global(cell_mass,
                                                dof_indices,
                                                mass_full);
      constraints_hn.distribute_local_to_global(cell_stiffness,
                                                dof_indices,
                                                stiffness_full);
    }

  mass_full.compress(VectorOperation::add);
  stiffness_full.compress(VectorOperation::add);
}

// Build R(t): the Dirichlet datum g on the boundary, 0 in the interior.
void
TimeADR::build_lifting(const double t, TrilinosWrappers::MPI::Vector &Rg) const
{
  GFunction g(g_fun);
  g.set_time(t);

  const IndexSet relevant =
    DoFTools::extract_locally_relevant_dofs(dof_handler);

  AffineConstraints<double> bc;
  bc.reinit(relevant);
  for (const types::boundary_id id : {types::boundary_id(0),
                                      types::boundary_id(1),
                                      types::boundary_id(2),
                                      types::boundary_id(3)})
    VectorTools::interpolate_boundary_values(dof_handler, id, g, bc);
  bc.close();

  // Start from 0 (interior stays 0) and let distribute() set the boundary DoFs
  // to their prescribed value g.
  Rg = 0.0;
  bc.distribute(Rg);
}

void
TimeADR::assemble_rhs()
{
  const auto t0 = Clock::now();

  // -------------------------------------------------------------------------
  // Forcing contribution: theta f^{n+1} + (1 - 2 theta) f^n + theta f^{n-1}.
  // (At this point "time" already equals t_{n+1}.)
  // -------------------------------------------------------------------------
  system_rhs = 0.0;

  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q           = quadrature->size();

  FEValues<dim> fe_values(*fe,
                          *quadrature,
                          update_values | update_quadrature_points |
                            update_JxW_values);

  Vector<double>                       cell_rhs(dofs_per_cell);
  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);
      cell_rhs = 0.0;

      for (unsigned int q = 0; q < n_q; ++q)
        {
          const Point<dim> &x_q = fe_values.quadrature_point(q);

          const double f_np1 = f(x_q, time);
          const double f_n   = f(x_q, time - delta_t);
          const double f_nm1 = f(x_q, time - 2.0 * delta_t);

          const double f_theta =
            theta * f_np1 + (1.0 - 2.0 * theta) * f_n + theta * f_nm1;

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            cell_rhs(i) +=
              f_theta * fe_values.shape_value(i, q) * fe_values.JxW(q);
        }

      cell->get_dof_indices(dof_indices);
      constraints.distribute_local_to_global(cell_rhs, dof_indices, system_rhs);
    }

  system_rhs.compress(VectorOperation::add);

  // -------------------------------------------------------------------------
  // Mass and stiffness contributions via matrix-vector products:
  //   (1/dt^2) M (2 u^n - u^{n-1})  -  A ( (1 - 2 theta) u^n + theta u^{n-1} ).
  // -------------------------------------------------------------------------
  const double inv_dt2 = 1.0 / (delta_t * delta_t);

  // Mass part.
  tmp = 0.0;
  tmp.add(2.0, u);
  tmp.add(-1.0, u_old);
  mass_matrix.vmult(tmp2, tmp);
  system_rhs.add(inv_dt2, tmp2);

  // Stiffness part.
  tmp = 0.0;
  tmp.add(1.0 - 2.0 * theta, u);
  tmp.add(theta, u_old);
  stiffness_matrix.vmult(tmp2, tmp);
  system_rhs.add(-1.0, tmp2);

  // -------------------------------------------------------------------------
  // Non-homogeneous Dirichlet: subtract the lifting correction
  //   (1/dt^2) M_full (R^{n+1} - 2 R^n + R^{n-1})
  //   +        A_full (theta R^{n+1} + (1-2theta) R^n + theta R^{n-1}),
  // with R^k the discrete lifting at t_k. Here "time" == t_{n+1}.
  // -------------------------------------------------------------------------
  if (g_fun)
    {
      TrilinosWrappers::MPI::Vector R_np1(u), R_n(u), R_nm1(u);
      build_lifting(time, R_np1);
      build_lifting(time - delta_t, R_n);
      build_lifting(time - 2.0 * delta_t, R_nm1);

      // Mass part.
      tmp = 0.0;
      tmp.add(1.0, R_np1);
      tmp.add(-2.0, R_n);
      tmp.add(1.0, R_nm1);
      mass_full.vmult(tmp2, tmp);
      system_rhs.add(-inv_dt2, tmp2);

      // Stiffness part.
      tmp = 0.0;
      tmp.add(theta, R_np1);
      tmp.add(1.0 - 2.0 * theta, R_n);
      tmp.add(theta, R_nm1);
      stiffness_full.vmult(tmp2, tmp);
      system_rhs.add(-1.0, tmp2);
    }

  // Enforce homogeneous Dirichlet on the right-hand side.
  constraints.set_zero(system_rhs);

  perf.assembly_seconds += elapsed_seconds(t0, Clock::now());
}

void
TimeADR::solve_linear_system()
{
  const auto t0 = Clock::now();

  // The system matrix (1/dt^2) M + theta A is symmetric positive definite:
  // use CG. The preconditioner was built once in setup().
  ReductionControl solver_control(/* maxiter   = */ 10000,
                                  /* tolerance = */ 1.0e-16,
                                  /* reduce    = */ 1.0e-8);

  SolverCG<TrilinosWrappers::MPI::Vector> solver(solver_control);
  solver.solve(lhs_matrix, u_new, system_rhs, preconditioner);

  constraints.distribute(u_new);

  perf.solve_seconds += elapsed_seconds(t0, Clock::now());
  perf.total_solver_iterations += solver_control.last_step();

  if (verbose)
    pcout << "  " << solver_control.last_step() << " CG iterations"
          << std::endl;
}

void
TimeADR::output() const
{
  const auto t0 = Clock::now();

  DataOut<dim> data_out;
  data_out.add_data_vector(dof_handler, solution, "solution");

  std::vector<unsigned int> partition_int(mesh.n_active_cells());
  GridTools::get_subdomain_association(mesh, partition_int);
  const Vector<double> partitioning(partition_int.begin(), partition_int.end());
  data_out.add_data_vector(partitioning, "partitioning");

  data_out.build_patches();

  const std::string output_file_name = std::visit(
    [](auto &&arg) -> std::string {
      using T = std::decay_t<decltype(arg)>;
      if constexpr (std::is_same_v<T, unsigned int>)
        return "output-Nel-" + std::to_string(arg);
      else if constexpr (std::is_same_v<T, std::string>)
        {
          const std::filesystem::path mesh_path(arg);
          return "output-" + mesh_path.stem().string();
        }
    },
    mesh_spec);

  data_out.write_vtu_with_pvtu_record("./",
                                      output_file_name,
                                      timestep_number,
                                      MPI_COMM_WORLD);

  perf.output_seconds += elapsed_seconds(t0, Clock::now());
}

void
TimeADR::run()
{
  const auto run_t0 = Clock::now();

  setup();

  // -------------------------------------------------------------------------
  // Initial data.
  // -------------------------------------------------------------------------
  // u^0 = u0.  ScalarFunctionFromFunctionObject wraps the plain std::function
  // (u0_fun) into a deal.II Function<dim> that interpolate() can consume.
  VectorTools::interpolate(dof_handler,
                           ScalarFunctionFromFunctionObject<dim>(u0_fun),
                           u);
  constraints.distribute(u);

  // u1 (initial velocity).
  TrilinosWrappers::MPI::Vector u1_vec(u);
  VectorTools::interpolate(dof_handler,
                           ScalarFunctionFromFunctionObject<dim>(u1_fun),
                           u1_vec);

  // -------------------------------------------------------------------------
  // Second-order consistent startup.
  //
  // From the PDE, u_tt(0) = mu Delta u0 + f(.,0), i.e. the initial acceleration
  // a0 solves the FE system   M a0 = F^0 - A u0.
  // Then a Taylor expansion gives the ghost step
  //   u^{-1} = u0 - dt u1 + (dt^2 / 2) a0,
  // which is O(dt^2) accurate (a first-order startup would spoil the global
  // second-order convergence in time).
  // -------------------------------------------------------------------------
  {
    // Right-hand side F^0 - A u0.
    tmp = 0.0; // will hold F^0
    {
      const unsigned int dofs_per_cell = fe->dofs_per_cell;
      const unsigned int n_q           = quadrature->size();

      FEValues<dim> fe_values(*fe,
                              *quadrature,
                              update_values | update_quadrature_points |
                                update_JxW_values);

      Vector<double>                       cell_rhs(dofs_per_cell);
      std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          if (!cell->is_locally_owned())
            continue;

          fe_values.reinit(cell);
          cell_rhs = 0.0;

          for (unsigned int q = 0; q < n_q; ++q)
            {
              const double f0 = f(fe_values.quadrature_point(q), 0.0);
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                cell_rhs(i) +=
                  f0 * fe_values.shape_value(i, q) * fe_values.JxW(q);
            }

          cell->get_dof_indices(dof_indices);
          constraints.distribute_local_to_global(cell_rhs, dof_indices, tmp);
        }
      tmp.compress(VectorOperation::add);
    }

    stiffness_matrix.vmult(tmp2, u); // tmp2 = A u0
    tmp.add(-1.0, tmp2);             // tmp  = F^0 - A u0

    // Non-homogeneous case: the physical initial acceleration also feels the
    // boundary lifting at t = 0, i.e. M a0 = F^0 - A (u0 + R^0). Subtract the
    // extra A_full R^0 coupling here.
    if (g_fun)
      {
        TrilinosWrappers::MPI::Vector R0(u);
        build_lifting(0.0, R0);
        stiffness_full.vmult(tmp2, R0);
        tmp.add(-1.0, tmp2);
      }

    constraints.set_zero(tmp);

    // Solve M a0 = tmp.
    TrilinosWrappers::MPI::Vector a0(u);
    a0 = 0.0;

    TrilinosWrappers::PreconditionSSOR mass_prec;
    mass_prec.initialize(
      mass_matrix, TrilinosWrappers::PreconditionSSOR::AdditionalData(1.0));

    ReductionControl solver_control(10000, 1.0e-16, 1.0e-10);
    SolverCG<TrilinosWrappers::MPI::Vector> solver(solver_control);
    solver.solve(mass_matrix, a0, tmp, mass_prec);
    constraints.distribute(a0);

    // u^{-1} = u0 - dt u1 + (dt^2 / 2) a0.
    u_old = u;
    u_old.add(-delta_t, u1_vec);
    u_old.add(0.5 * delta_t * delta_t, a0);
    constraints.distribute(u_old);
  }

  // Output the initial condition.
  time            = 0.0;
  timestep_number = 0;
  solution        = u;
  if (write_output)
    output();

  if (verbose)
    pcout << "===============================================" << std::endl;

  // -------------------------------------------------------------------------
  // Energy monitoring!!!
  // The discrete energy  E^n = 1/2 (v^n)^T M v^n + 1/2 (u^n)^T A u^n , with the
  // velocity reconstructed by central differences v^n =
  // (u^{n+1}-u^{n-1})/(2dt), is conserved for gamma = 1/2 (no numerical
  // dissipation).
  // -------------------------------------------------------------------------
  std::ofstream energy_csv;
  if (write_output && mpi_rank == 0)
    {
      energy_csv.open("energy.csv");
      energy_csv << "time,energy,kinetic,potential\n";
    }
  double energy0    = 0.0;
  double energy_min = std::numeric_limits<double>::max();
  double energy_max = -std::numeric_limits<double>::max();

  while (time < T - 0.5 * delta_t)
    {
      time += delta_t;
      ++timestep_number;
      perf.n_steps += 1;

      if (verbose)
        pcout << "Timestep " << std::setw(4) << timestep_number
              << ", time = " << std::setw(6) << std::fixed
              << std::setprecision(4) << time << std::endl;

      assemble_rhs();
      solve_linear_system(); // fills u_new

      // Discrete energy at time level n (t = time - dt)
      // At this point: u_old = u^{n-1}, u = u^n, u_new = u^{n+1}.
      {
        tmp = 0.0;
        tmp.add(1.0 / (2.0 * delta_t), u_new);
        tmp.add(-1.0 / (2.0 * delta_t), u_old);
        mass_matrix.vmult(tmp2, tmp);
        const double kinetic = 0.5 * (tmp * tmp2);

        // Potential: 1/2 u^T A u.
        stiffness_matrix.vmult(tmp2, u);
        const double potential = 0.5 * (u * tmp2);

        const double energy = kinetic + potential;

        if (timestep_number == 1)
          energy0 = energy;
        energy_min = std::min(energy_min, energy);
        energy_max = std::max(energy_max, energy);

        if (write_output && mpi_rank == 0)
          energy_csv << (time - delta_t) << "," << energy << "," << kinetic
                     << "," << potential << "\n";
      }

      // Shift the time levels: u^{n-1} <- u^n, u^n <- u^{n+1}.
      u_old = u;
      u     = u_new;

      // Update the ghosted copy and write output every "output_every" steps.
      if (write_output && timestep_number % output_every == 0)
        {
          solution = u;
          output();
        }
    }

  // Make sure the ghosted copy holds the final solution (for compute_error).
  // Non-homogeneous case: the loop solved for the boundary-zero part w, so the
  // physical solution is u = w + R(T).
  if (g_fun)
    {
      build_lifting(time, tmp); // R(T)
      tmp.add(1.0, u);          // w + R(T)
      solution = tmp;
    }
  else
    solution = u;

  if (write_output && mpi_rank == 0)
    energy_csv.close();

  // Report the energy conservation: the relative drift should be O(dt^2),
  // i.e. tiny and NOT growing (no numerical dissipation for gamma = 1/2).
  if (verbose && energy0 != 0.0)
    {
      const double drift = (energy_max - energy_min) / energy0;
      pcout << "\nEnergy: E0 = " << energy0 << ", relative drift = " << drift
            << std::endl;
    }

  perf.total_run_seconds = elapsed_seconds(run_t0, Clock::now());
  if (verbose)
    print_performance_report();
}

double
TimeADR::compute_error(const VectorTools::NormType &cell_norm,
                       Function<dim>               &exact_solution) const
{
  exact_solution.set_time(time);

  const QGaussSimplex<dim> quadrature_error(r + 2);

  FE_SimplexP<dim> fe_linear(1);
  MappingFE<dim>   mapping(fe_linear);

  Vector<float> error_per_cell(mesh.n_active_cells());

  VectorTools::integrate_difference(mapping,
                                    dof_handler,
                                    solution,
                                    exact_solution,
                                    error_per_cell,
                                    quadrature_error,
                                    cell_norm);

  const double global_error =
    VectorTools::compute_global_error(mesh, error_per_cell, cell_norm);

  return global_error;
}

void
TimeADR::print_performance_report() const
{
  const double total_run_max =
    Utilities::MPI::max(perf.total_run_seconds, MPI_COMM_WORLD);
  const double setup_max =
    Utilities::MPI::max(perf.setup_seconds, MPI_COMM_WORLD);
  const double assembly_max =
    Utilities::MPI::max(perf.assembly_seconds, MPI_COMM_WORLD);
  const double solve_max =
    Utilities::MPI::max(perf.solve_seconds, MPI_COMM_WORLD);
  const double output_max =
    Utilities::MPI::max(perf.output_seconds, MPI_COMM_WORLD);
  const double solver_iters_max =
    Utilities::MPI::max(static_cast<double>(perf.total_solver_iterations),
                        MPI_COMM_WORLD);

  if (mpi_rank == 0)
    {
      std::cout << "\n===============================================\n";
      std::cout << "Performance report\n";
      std::cout << "-----------------------------------------------\n";
      std::cout << "Total run time         : " << total_run_max << " s\n";
      std::cout << "Setup time             : " << setup_max << " s\n";
      std::cout << "RHS assembly time      : " << assembly_max << " s\n";
      std::cout << "Solve time             : " << solve_max << " s\n";
      std::cout << "Output time            : " << output_max << " s\n";
      std::cout << "Number of timesteps    : " << perf.n_steps << "\n";
      std::cout << "Total CG iterations    : " << solver_iters_max << "\n";
      std::cout << "===============================================\n";
    }
}