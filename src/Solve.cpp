#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/numerics/vector_tools.h>

#include "NavierStokes.hpp"

template class NavierStokes<2U>;
template class NavierStokes<3U>;

template <unsigned int dim>
void NavierStokes<dim>::apply_initial_conditions() {
  pcout << "Applying the initial conditions" << std::endl;

  initial_conditions->set_time(time_step * deltat);
  VectorTools::interpolate(dof_handler, *initial_conditions, solution_owned);
  solution = solution_owned;
}

template <unsigned int dim>
void NavierStokes<dim>::solve_time_step() {
  pcout << "  Solving the linear system" << std::endl;

  SolverControl solver_control(100000, 1e-3 * system_rhs.l2_norm());

  SolverGMRES<TrilinosWrappers::MPI::BlockVector> solver(solver_control);

  solver.solve(system_matrix, solution_owned, system_rhs,
               PreconditionIdentity());
  pcout << "  " << solver_control.last_step() << " GMRES iterations"
        << std::endl;

  solution = solution_owned;
}

template <unsigned int dim>
void NavierStokes<dim>::solve() {
  assemble_constant();

  pcout << "===============================================" << std::endl;

  // Calculate and output the initial solution.
  time_step = 0;
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

    time += deltat;
    ++time_step;
  }
}
