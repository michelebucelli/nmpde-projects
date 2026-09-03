#pragma once

#include <cstdlib>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_gmres.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <memory>

#include "cdr_assembled.hpp"
#include "cdr_operator.hpp"
#include "cdr_parameters.hpp"
#include "cdr_problem.hpp"
#include "mg_preconditioner.hpp"

using namespace dealii;

// Beta is stored dimension independently in the input file, so only the
// leading `dim` components are used here.
template <int dim>
auto make_coefficients(const Parameters &parameters) -> Coefficients<dim> {
  Tensor<1, dim> beta;
  for (unsigned int d = 0; d < dim; ++d) {
    beta[d] = parameters.beta[d];
  }
  return Coefficients<dim>(parameters.mu, parameters.gamma, beta);
}

// Distributed matrix-free CdrProblem.
template <int dim, int fe_degree> class CdrProblem {
public:
  using Number = double;
  using LevelNumber = float; // mixed precision on the multigrid levels.
  using SystemMatrixType = Operator<dim, fe_degree, Number>;
  using VectorType = LinearAlgebra::distributed::Vector<Number>;
  using Multigrid = MultigridPreconditioner<dim, fe_degree, LevelNumber>;

  CdrProblem(const Parameters &parameters);

  auto run() -> void;

private:
  auto setup_system() -> void;
  auto assemble_rhs() -> void;
  auto local_assemble_rhs(
      const MatrixFree<dim, Number> &data, VectorType &destination,
      const VectorType &source,
      const std::pair<unsigned int, unsigned int> &cell_range) const -> void;
  auto solve() -> SolveResult;
  auto benchmark_operator() const -> void;
  auto memory_consumption() const -> std::size_t;
  auto is_matrix_free() const -> bool;
  auto compute_errors() const -> std::pair<double, double>;
  auto output_results(const unsigned int cycle) const -> void;

  auto dirichlet_boundary_ids() const -> std::set<types::boundary_id>;
  auto is_root() const -> bool;
  auto run_cycle(const unsigned int cycle) -> void;
  auto print_header() const -> void;
  auto print_convergence_table() -> void;

  const Parameters parameters;
  const Coefficients<dim> coefficients;
  const MPI_Comm communicator;
  ConditionalOStream pcout;
  TimerOutput timer;
  ConvergenceTable convergence_table;

  parallel::distributed::Triangulation<dim> triangulation;
  const MappingQ1<dim> mapping;
  const FE_Q<dim> fe;
  DoFHandler<dim> dof_handler;

  AffineConstraints<Number> constraints;
  SystemMatrixType system_matrix;
  VectorType solution;
  VectorType system_rhs;

  Multigrid multigrid;
  AssembledSystem<dim, fe_degree> assembled;
  CycleTimings timings;
};

template <int dim, int fe_degree>
CdrProblem<dim, fe_degree>::CdrProblem(const Parameters &parameters)
    : parameters(parameters), coefficients(make_coefficients<dim>(parameters)),
      communicator(MPI_COMM_WORLD),
      pcout(std::cout, parameters.verbose && Utilities::MPI::this_mpi_process(
                                                 MPI_COMM_WORLD) == 0),
      timer(communicator, pcout, TimerOutput::never, TimerOutput::wall_times),
      triangulation(communicator,
                    Triangulation<dim>::limit_level_difference_at_vertices,
                    parallel::distributed::Triangulation<
                        dim>::construct_multigrid_hierarchy),
      fe(fe_degree), dof_handler(triangulation) {}

template <int dim, int fe_degree>
auto CdrProblem<dim, fe_degree>::is_matrix_free() const -> bool {
  return parameters.backend == Backend::MatrixFree;
}

template <int dim, int fe_degree>
auto CdrProblem<dim, fe_degree>::is_root() const -> bool {
  return Utilities::MPI::this_mpi_process(communicator) == 0;
}

template <int dim, int fe_degree>
auto CdrProblem<dim, fe_degree>::dirichlet_boundary_ids() const
    -> std::set<types::boundary_id> {
  std::set<types::boundary_id> ids;
  for (const types::boundary_id id : triangulation.get_boundary_ids()) {
    if (id != neumann_boundary_id) {
      ids.insert(id);
    }
  }
  return ids;
}

template <int dim, int fe_degree>
auto CdrProblem<dim, fe_degree>::setup_system() -> void {
  TimerOutput::Scope scope(timer, "setup");
  Timer wall;

  const bool use_multigrid =
      parameters.preconditioner == PreconditionerType::Multigrid &&
      is_matrix_free();

  dof_handler.distribute_dofs(fe);
  if (use_multigrid) {
    dof_handler.distribute_mg_dofs();
  }

  constraints.clear();
  constraints.reinit(DoFTools::extract_locally_relevant_dofs(dof_handler));

  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  for (const types::boundary_id id : dirichlet_boundary_ids()) {
    VectorTools::interpolate_boundary_values(mapping, dof_handler, id,
                                             ExactSolution<dim>(), constraints);
  }
  constraints.close();

  if (!is_matrix_free()) {
    assembled.setup(dof_handler, constraints, communicator);
    solution.reinit(dof_handler.locally_owned_dofs(),
                    DoFTools::extract_locally_relevant_dofs(dof_handler),
                    communicator);
    timings.setup = wall.wall_time();
    return;
  }

  auto additional_data = matrix_free_data<dim, Number>(
      update_gradients | update_JxW_values | update_quadrature_points);

  additional_data.mapping_update_flags_boundary_faces =
      update_values | update_JxW_values | update_quadrature_points |
      update_normal_vectors;

  auto matrix_free = std::make_shared<MatrixFree<dim, Number>>();
  matrix_free->reinit(mapping, dof_handler, constraints,
                      QGauss<1>(fe_degree + 1), additional_data);

  system_matrix.initialize(matrix_free);
  system_matrix.set_coefficients(coefficients);
  system_matrix.initialize_dof_vector(solution);
  system_matrix.initialize_dof_vector(system_rhs);

  if (parameters.preconditioner == PreconditionerType::Jacobi) {
    system_matrix.compute_diagonal();
  }

  if (use_multigrid) {
    multigrid.build(dof_handler, mapping, coefficients,
                    dirichlet_boundary_ids());
  }

  timings.setup = wall.wall_time();
}

template <int dim, int fe_degree>
auto CdrProblem<dim, fe_degree>::assemble_rhs() -> void {
  TimerOutput::Scope scope(timer, "assemble rhs");
  Timer wall;

  if (!is_matrix_free()) {
    assembled.assemble(mapping, dof_handler, constraints, coefficients);
    timings.rhs = wall.wall_time();
    return;
  }

  const MatrixFree<dim, Number> &matrix_free = *system_matrix.get_matrix_free();
  const ExactSolution<dim> exact_solution;

  solution = 0.0;
  constraints.distribute(solution);
  solution.update_ghost_values();

  system_rhs = 0.0;
  matrix_free.cell_loop(&CdrProblem::local_assemble_rhs, this, system_rhs,
                        solution);

  FEFaceEvaluation<dim, fe_degree, fe_degree + 1, 1, Number> phi_face(
      matrix_free, true);
  const unsigned int first_face = matrix_free.n_inner_face_batches();
  const unsigned int last_face =
      first_face + matrix_free.n_boundary_face_batches();

  for (unsigned int face = first_face; face < last_face; ++face) {
    if (matrix_free.get_boundary_id(face) != neumann_boundary_id) {
      continue;
    }

    phi_face.reinit(face);
    for (unsigned int q = 0; q < phi_face.n_q_points; ++q) {
      const auto normal = phi_face.normal_vector(q);
      const auto point = phi_face.quadrature_point(q);

      VectorizedArray<Number> h;
      for (unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v) {
        Point<dim> p;
        Tensor<1, dim> n;
        for (unsigned int d = 0; d < dim; ++d) {
          p[d] = point[d][v];
          n[d] = normal[d][v];
        }
        h[v] = exact_solution.gradient(p) * n;
      }

      phi_face.submit_value(coefficients.mu * h, q);
    }
    phi_face.integrate_scatter(EvaluationFlags::values, system_rhs);
  }
  system_rhs.compress(VectorOperation::add);

  solution.zero_out_ghost_values();
  timings.rhs = wall.wall_time();
}

template <int dim, int fe_degree>
auto CdrProblem<dim, fe_degree>::local_assemble_rhs(
    const MatrixFree<dim, Number> &data, VectorType &destination,
    const VectorType &source,
    const std::pair<unsigned int, unsigned int> &cell_range) const -> void {
  const RightHandSide<dim> right_hand_side(coefficients);
  FEEvaluation<dim, fe_degree, fe_degree + 1, 1, Number> phi(data);

  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell) {
    phi.reinit(cell);
    phi.read_dof_values_plain(source);
    phi.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);

    system_matrix.do_quadrature_points(
        phi, -1.0, [&](const auto &phi_, const unsigned int q) {
          return evaluate_lanewise<dim, Number>(
              [&](const Point<dim> &p) { return right_hand_side.value(p); },
              phi_.quadrature_point(q));
        });

    phi.integrate_scatter(EvaluationFlags::values | EvaluationFlags::gradients,
                          destination);
  }
}

template <int dim, int fe_degree>
auto CdrProblem<dim, fe_degree>::solve() -> SolveResult {
  TimerOutput::Scope scope(timer, "solve");
  Timer wall;

  if (!is_matrix_free()) {
    const SolveResult result = assembled.solve(parameters, constraints);
    timings.solve = wall.wall_time();

    // The error and output paths are shared, so bring the solution over into
    // the deal.II vector they expect.
    const auto &assembled_solution = assembled.get_solution();
    for (const types::global_dof_index i : dof_handler.locally_owned_dofs()) {
      solution(i) = assembled_solution(i);
    }
    solution.compress(VectorOperation::insert);

    return result;
  }

  SolverControl solver_control(5000,
                               parameters.tolerance * system_rhs.l2_norm());

  typename SolverFGMRES<VectorType>::AdditionalData fgmres_data;
  fgmres_data.max_basis_size = 50;
  SolverFGMRES<VectorType> solver(solver_control, fgmres_data);

  VectorType correction;
  system_matrix.initialize_dof_vector(correction);

  const auto solve_with = [&](const auto &preconditioner) {
    solver.solve(system_matrix, correction, system_rhs, preconditioner);
  };

  switch (parameters.preconditioner) {
  case PreconditionerType::None:
    solve_with(PreconditionIdentity());
    break;
  case PreconditionerType::Jacobi:
    solve_with(*system_matrix.get_matrix_diagonal_inverse());
    break;
  case PreconditionerType::Multigrid:
    solve_with(multigrid.get());
    break;
  }

  timings.solve = wall.wall_time();

  VectorType residual;
  system_matrix.initialize_dof_vector(residual);
  system_matrix.vmult(residual, correction);
  residual -= system_rhs;

  solution += correction;
  constraints.distribute(solution);

  return {solver_control.last_step(),
          residual.l2_norm() / system_rhs.l2_norm()};
}

// Times bare operator applications, which isolates the matrix-free cell loop
// from the Krylov vector algebra and the multigrid.
template <int dim, int fe_degree>
auto CdrProblem<dim, fe_degree>::benchmark_operator() const -> void {
  if (!is_matrix_free()) {
    const double time =
        assembled.benchmark(parameters.n_benchmarks_applications, communicator);
    if (is_root()) {
      std::cout << "VMULT backend=matrix_based"
                << " ranks=" << Utilities::MPI::n_mpi_processes(communicator)
                << " threads=" << MultithreadInfo::n_threads()
                << " dofs=" << dof_handler.n_dofs() << " time=" << time
                << " mdofs_per_s=" << dof_handler.n_dofs() / time / 1e6
                << std::endl;
    }
    return;
  }

  VectorType src, dst;
  system_matrix.initialize_dof_vector(src);
  system_matrix.initialize_dof_vector(dst);
  src = 1.0;

  system_matrix.vmult(dst, src); // warmup the caches
  MPI_Barrier(communicator);

  Timer wall;
  for (unsigned int i = 0; i < parameters.n_benchmarks_applications; ++i) {
    system_matrix.vmult(dst, src);
  }

  const double time = Utilities::MPI::max(wall.wall_time(), communicator) /
                      parameters.n_benchmarks_applications;

  if (is_root()) {
    std::cout << "VMULT backend=matrix_free"
              << " ranks=" << Utilities::MPI::n_mpi_processes(communicator)
              << " threads=" << MultithreadInfo::n_threads()
              << " dofs=" << dof_handler.n_dofs() << " time=" << time
              << " mdofs_per_s=" << dof_handler.n_dofs() / time / 1e6
              << std::endl;
  }
}

template <int dim, int fe_degree>
auto CdrProblem<dim, fe_degree>::compute_errors() const
    -> std::pair<double, double> {
  const QGauss<dim> quadrature(fe_degree + 2);
  Vector<double> cellwise_error(triangulation.n_active_cells());

  const auto global_error = [&](const VectorTools::NormType norm) {
    VectorTools::integrate_difference(mapping, dof_handler, solution,
                                      ExactSolution<dim>(), cellwise_error,
                                      quadrature, norm);
    return VectorTools::compute_global_error(triangulation, cellwise_error,
                                             norm);
  };

  solution.update_ghost_values();
  const std::pair<double, double> errors = {
      global_error(VectorTools::L2_norm),
      global_error(VectorTools::H1_seminorm),
  };
  solution.zero_out_ghost_values();

  return errors;
}

template <int dim, int fe_degree>
auto CdrProblem<dim, fe_degree>::output_results(const unsigned int cycle) const
    -> void {
  if (dof_handler.n_dofs() > 5e6) {
    return;
  }

  VectorType ghosted(solution);
  ghosted.update_ghost_values();

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(ghosted, "u");

  Vector<double> subdomain(triangulation.n_active_cells());
  subdomain = Utilities::MPI::this_mpi_process(communicator);
  data_out.add_data_vector(subdomain, "rank");

  data_out.build_patches(mapping, fe_degree);
  data_out.write_vtu_with_pvtu_record(
      "./", "solution-" + std::to_string(dim) + "d", cycle, communicator, 2);
}

// Bytes the operator itself holds, summed over the ranks. For matrix-free
// this is the MatrixFree bookkeeping plus the multigrid levels; for the
// matrix-based backend it is dominated by the sparse matrix entries.
template <int dim, int fe_degree>
auto CdrProblem<dim, fe_degree>::memory_consumption() const -> std::size_t {
  if (!is_matrix_free()) {
    return assembled.memory_consumption(communicator);
  }

  const std::size_t local =
      system_matrix.get_matrix_free()->memory_consumption() +
      multigrid.memory_consumption();
  return Utilities::MPI::sum(local, communicator);
}

template <int dim, int fe_degree>
auto CdrProblem<dim, fe_degree>::print_header() const -> void {
  pcout << "Convection-diffusion-reaction, " << to_string(parameters.backend)
        << ", " << dim << "D, Q" << fe_degree << ": mu = " << coefficients.mu
        << ", gamma = " << coefficients.gamma
        << ", beta = " << coefficients.beta << std::endl
        << Utilities::MPI::n_mpi_processes(communicator) << " MPI rank(s), "
        << MultithreadInfo::n_threads() << " TBB thread(s), "
        << VectorizedArray<Number>::size() << " SIMD lanes" << std::endl;
}

template <int dim, int fe_degree>
auto CdrProblem<dim, fe_degree>::run_cycle(const unsigned int cycle) -> void {
  const unsigned int n_refinements = parameters.n_refinements + cycle;

  // Drop MatrixFree objects that still refer to the old mesh before it is
  // destroyed.
  multigrid.clear();
  triangulation.clear();
  GridGenerator::hyper_cube(triangulation, 0.0, 1.0, true);
  triangulation.refine_global(n_refinements);

  setup_system();
  assemble_rhs();

  if (parameters.n_benchmarks_applications > 0)
    benchmark_operator();

  const SolveResult result = solve();
  const auto [l2_error, h1_error] = compute_errors();

  if (parameters.write_output)
    output_results(cycle);

  // The timings are wall times of this rank; a scaling table needs the
  // slowest rank, not whichever one happens to be the root.
  timings.setup = Utilities::MPI::max(timings.setup, communicator);
  timings.rhs = Utilities::MPI::max(timings.rhs, communicator);
  timings.solve = Utilities::MPI::max(timings.solve, communicator);

  const std::size_t memory = memory_consumption();

  convergence_table.add_value("cells", triangulation.n_global_active_cells());
  convergence_table.add_value("dofs", dof_handler.n_dofs());
  convergence_table.add_value("L2", l2_error);
  convergence_table.add_value("H1", h1_error);
  convergence_table.add_value("its", result.n_iterations);
  convergence_table.add_value("res", result.relative_residual);
  convergence_table.add_value("solve[s]", timings.solve);

  pcout << "cycle " << cycle << ": " << triangulation.n_global_active_cells()
        << " cells, " << dof_handler.n_dofs() << " dofs, "
        << result.n_iterations << " iterations, " << timings.setup
        << " s setup, " << timings.rhs << " s rhs, " << timings.solve
        << " s solve" << std::endl;

  if (is_root()) {
    std::cout << "SUMMARY"
              << " backend=" << to_string(parameters.backend)
              << " precond=" << to_string(parameters.preconditioner)
              << " dim=" << dim << " degree=" << fe_degree
              << " refine=" << n_refinements
              << " ranks=" << Utilities::MPI::n_mpi_processes(communicator)
              << " threads=" << MultithreadInfo::n_threads()
              << " dofs=" << dof_handler.n_dofs()
              << " iters=" << result.n_iterations << " setup=" << timings.setup
              << " rhs=" << timings.rhs << " solve=" << timings.solve
              << " memory_mb=" << memory / 1024.0 / 1024.0 << " l2=" << l2_error
              << " h1=" << h1_error << std::endl;
  }
}

template <int dim, int fe_degree>
auto CdrProblem<dim, fe_degree>::print_convergence_table() -> void {
  for (const std::string column : {"L2", "H1", "res"}) {
    convergence_table.set_precision(column, 3);
    convergence_table.set_scientific(column, true);
  }
  convergence_table.set_precision("solve[s]", 3);

  if (parameters.n_cycles > 1) {
    for (const std::string column : {"L2", "H1"})
      convergence_table.evaluate_convergence_rates(
          column, "dofs", ConvergenceTable::reduction_rate_log2, dim);
  }

  if (parameters.verbose && is_root()) {
    std::cout << std::endl;
    convergence_table.write_text(std::cout);
    std::cout << std::endl;
  }
}

template <int dim, int fe_degree>
auto CdrProblem<dim, fe_degree>::run() -> void {
  print_header();

  for (unsigned int cycle = 0; cycle < parameters.n_cycles; ++cycle)
    run_cycle(cycle);

  print_convergence_table();
  timer.print_wall_time_statistics(communicator);
}
