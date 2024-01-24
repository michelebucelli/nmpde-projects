#include "Problem.hpp"

// Main function.
int main(int argc, char *argv[]) {
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const std::string mesh_file_name =
      "../mesh/2d-flow-factor-1.msh";  //"../mesh/cube-0.1.msh";
  constexpr unsigned int degree_velocity = 2;
  constexpr unsigned int degree_pressure = 1;
  constexpr double T = 10;      // 1e-2;
  constexpr double deltat = 1;  // 1e-3;
  // constexpr double nu = 0.01;

  // EthierSteinman problem(mesh_file_name, degree_velocity, degree_pressure, T,
  //                        deltat, nu);

  Cylinder2D problem(mesh_file_name, degree_velocity, degree_pressure, T,
                     deltat);

  problem.setup();
  problem.solve();

  return 0;
}