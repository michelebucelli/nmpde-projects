#include "NavierStokes.hpp"

#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/numerics/vector_tools.h>

template class NavierStokes<2U>;
template class NavierStokes<3U>;

template <unsigned int dim>
void NavierStokes<dim>::solve_time_step()
{
  SolverControl solver_control(10000, 1e-6 * system_rhs.l2_norm());

  SolverGMRES<TrilinosWrappers::MPI::BlockVector> solver(solver_control);

  solver.solve(system_matrix, solution_owned, system_rhs, PreconditionIdentity());
  pcout << "  " << solver_control.last_step() << " GMRES iterations" << std::endl;

  solution = solution_owned;
}

template <unsigned int dim>
void NavierStokes<dim>::solve()
{
  assemble_constant();

  pcout << "===============================================" << std::endl;

  // Apply the initial condition.
  {
    pcout << "Applying the initial condition" << std::endl;

    VectorTools::interpolate(dof_handler, *initial_conditions, solution_owned);
    solution = solution_owned;

    // Output the initial solution.
    output(0);
    pcout << "-----------------------------------------------" << std::endl;
  }

  unsigned int time_step = 0;
  double time            = 0;

  // Solvethe problem at each time step.
  while (time < T)
    {
      time += deltat;
      ++time_step;

      pcout << "n = " << std::setw(3) << time_step << ", t = " << std::setw(5)
            << time << ":" << std::flush;

      assemble_time_dependent();
      solve_time_step();
      output(time_step);
    }
}
