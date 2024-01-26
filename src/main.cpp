#include "Problem.hpp"

// Main function.
int main(int argc, char *argv[]) {
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const std::string mesh_file_name = "../mesh/3d-flow-factor-1.msh";
  constexpr unsigned int degree_velocity = 2;
  constexpr unsigned int degree_pressure = 1;
  constexpr double T = 1;
  constexpr double deltat = 0.1;
  // constexpr double nu = 0.01;

  Cylinder3D problem(mesh_file_name, degree_velocity, degree_pressure, T,
                     deltat);

  /*
  EthierSteinman problem(mesh_file_name, degree_velocity, degree_pressure, T,
                         deltat, nu);
                         */

  problem.setup();
  problem.solve();

  return 0;
}