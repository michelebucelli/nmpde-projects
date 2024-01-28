#include <deal.II/lac/solver_gmres.h>
#include <deal.II/numerics/vector_tools.h>

#include <chrono>

#include "NavierStokes.hpp"
#include "Precondition.hpp"

template class NavierStokes<2U>;
template class NavierStokes<3U>;

template <unsigned int dim>
void NavierStokes<dim>::apply_initial_conditions() {
  pcout << "  Applying the initial conditions" << std::endl;

  // Since only velocity has an initial condition, create a mask.
  ComponentMask mask;
  if constexpr (dim == 2) {
    mask = ComponentMask({true, true, false});
  } else {
    mask = ComponentMask({true, true, true, false});
  }

  VectorTools::interpolate(dof_handler, *initial_conditions, solution_owned,
                           mask);
  solution = solution_owned;
}

template <unsigned int dim>
void NavierStokes<dim>::solve_time_step() {
  auto t0 = std::chrono::high_resolution_clock::now();

  pcout << "  Initializing the preconditioner" << std::endl;

  // Choose the preconditioner and initialize it.
  std::shared_ptr<BlockPrecondition> precondition;
  switch (solver_options.id) {
    case BLOCK_DIAGONAL: {
      std::shared_ptr<PreconditionBlockDiagonal> actual_precondition =
          std::make_shared<PreconditionBlockDiagonal>();
      actual_precondition->initialize(
          system_matrix.block(0, 0), pressure_mass.block(1, 1),
          solver_options.maxiter_inner, solver_options.tol_inner);
      precondition = actual_precondition;
      break;
    }
    case SIMPLE: {
      std::shared_ptr<PreconditionSIMPLE> actual_precondition =
          std::make_shared<PreconditionSIMPLE>();
      actual_precondition->initialize(
          system_matrix.block(0, 0), system_matrix.block(1, 0),
          system_matrix.block(0, 1), solution_owned, solver_options.alpha,
          solver_options.maxiter_inner, solver_options.tol_inner);
      precondition = actual_precondition;
      break;
    }
    case ASIMPLE: {
      std::shared_ptr<PreconditionaSIMPLE> actual_precondition =
          std::make_shared<PreconditionaSIMPLE>();
      actual_precondition->initialize(
          system_matrix.block(0, 0), system_matrix.block(1, 0),
          system_matrix.block(0, 1), solution_owned, solver_options.alpha,
          solver_options.use_inner_solver, solver_options.maxiter_inner,
          solver_options.tol_inner);
      precondition = actual_precondition;
      break;
    }
    case YOSHIDA: {
      std::shared_ptr<PreconditionYoshida> actual_precondition =
          std::make_shared<PreconditionYoshida>();
      actual_precondition->initialize(
          system_matrix.block(0, 0), system_matrix.block(1, 0),
          system_matrix.block(0, 1), velocity_mass.block(0, 0), solution_owned,
          solver_options.maxiter_inner, solver_options.tol_inner);
      precondition = actual_precondition;
      break;
    }
  }

  // Solve the system.
  pcout << "  Solving the linear system" << std::endl;

  SolverControl solver_control(solver_options.maxiter,
                               solver_options.tol * system_rhs.l2_norm());
  SolverGMRES<TrilinosWrappers::MPI::BlockVector> solver(solver_control);
  solver.solve(system_matrix, solution_owned, system_rhs, *precondition);
  pcout << "  " << solver_control.last_step() << " GMRES iterations"
        << std::endl;

  solution = solution_owned;

  // Print the execution time.
  auto t1 = std::chrono::high_resolution_clock::now();
  auto dt =
      std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
  pcout << "  Total linear system resolution time: " << dt << "μs" << std::endl;
}

template <unsigned int dim>
void NavierStokes<dim>::solve() {
  assemble_constant();

  pcout << "===============================================" << std::endl;

  // Calculate and output the initial solution.
  time_step = 0;
  pcout << "n = " << std::setw(3) << time_step << ", t = " << std::setw(5) << 0
        << ": " << std::endl;
  apply_initial_conditions();
  output();

  time_step = 1;
  double time = deltat;

  // Solve the problem at each time step.
  while (time < T + deltat) {
    pcout << "n = " << std::setw(3) << time_step << ", t = " << std::setw(5)
          << time << ": " << std::endl;

    assemble_time_dependent();
    solve_time_step();
    output();

    time_step += 1;
    time = time_step * deltat;
  }
}
