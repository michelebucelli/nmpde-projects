#include "FK_solver.hpp"
#include <deal.II/base/convergence_table.h>

// Main function.
int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const unsigned int degree = 1;

  double T = 5;
  double deltat = 0.1;
  FK_solver problem("/home/samuele/Primo_semestre/Numerical_methods_for_partial_differential_equations/temp/FK_solver/mesh/mesh-cube-10.msh", degree, T, deltat);
  problem.setup();
  problem.solve();

  return 0;
}
