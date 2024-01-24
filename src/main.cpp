#include "Problem.hpp"

// Main function.
int main(int argc, char *argv[]) {
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const std::string mesh_file_name = "../mesh/cube-0.25.msh";
  constexpr unsigned int degree_velocity = 2;
  constexpr unsigned int degree_pressure = 1;
  constexpr double T = 1e-2;
  constexpr double deltat = 1e-3;
  constexpr double nu = 0.01;

  EthierSteinman problem(mesh_file_name, degree_velocity, degree_pressure, T,
                         deltat, nu);

  problem.setup();
  problem.solve();

  return 0;
}