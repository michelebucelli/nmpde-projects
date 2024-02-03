#include "FK_solver.hpp"
#include <deal.II/base/convergence_table.h>

// Main function.
int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const unsigned int degree = 1;

  double T = 2.0;
  double deltat = 0.01;
  FK_solver problem("../mesh/brain-h3.0.msh", degree, T, deltat);
  problem.setup();
  problem.solve();

  return 0;
}
