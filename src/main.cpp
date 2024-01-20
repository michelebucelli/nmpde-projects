#include "Problem.hpp"

// Main function.
int main(int argc, char *argv[]) {
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const std::string mesh_file_name = "../mesh/mesh-square-h0.100000.msh";
  constexpr unsigned int degree_velocity = 2;
  constexpr unsigned int degree_pressure = 1;
  constexpr double T = 10;
  constexpr double deltat = 0.1;

  Cylinder2D problem(mesh_file_name, degree_velocity, degree_pressure, T,
                     deltat);

  problem.setup();
  problem.solve();

  return 0;
}