#include "FK_solver.hpp"
#include <deal.II/base/convergence_table.h>


// Main function.
int
main(int argc, char *argv[])
{
    Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
   
  const unsigned int degree = 2;

   double T      = 3.0;
   double deltat = 0.1;
   FK_solver problem("../mesh/brain-h3.0.msh", degree, T, deltat);
   problem.setup();
   problem.solve();




  return 0;
}
