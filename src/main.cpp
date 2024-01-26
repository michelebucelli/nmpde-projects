#include <getopt.h>

#include "Problem.hpp"
// Main function.
int main(int argc, char* argv[]) {
  int problem_id = 0;
  double deltat = 0.0;
  std::string mesh_file;
  std::string err_msg =
      "Required Arguments: \n"
      "-p problem_id(1: Cylinder2D, 2: Cylinder3D, 3: "
      "EthierSteinman, 4: Step) \n"
      "-t deltat \n"
      "-m path_to_mesh_file";

  if (argc < 4) {
    std::cerr << err_msg << std::endl;
    return 1;
  }
  const char* const short_opts = "p:t:m:";
  const option long_opts[] = {{"problem", required_argument, nullptr, 'p'},
                              {"deltat", required_argument, nullptr, 't'},
                              {"mesh_file", required_argument, nullptr, 'm'},
                              {nullptr, 0, nullptr, 0}};

  int opt;
  while ((opt = getopt_long(argc, argv, short_opts, long_opts, nullptr)) !=
         -1) {
    switch (opt) {
      case 'p':
        problem_id = std::stoi(optarg);
        break;
      case 't':
        deltat = std::stod(optarg);
        break;
      case 'm':
        mesh_file = optarg;
        break;
      default:
        std::cerr << err_msg << std::endl;
        return 1;
    }
  }
  // Check if all required parameters are provided
  if (problem_id == 0 || deltat == 0.0 || mesh_file.empty()) {
    std::cerr << "Error: All parameters (-p, -d, -m) must be provided."
              << std::endl;
    return 1;
  }

  // Print the parsed arguments
  std::cout << "arg2: " << problem_id << std::endl;
  std::cout << "arg3: " << deltat << std::endl;
  std::cout << "arg4: " << mesh_file << std::endl;

  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const std::string mesh_file_name = "../mesh/cube-0.25.msh";
  constexpr unsigned int degree_velocity = 2;
  constexpr unsigned int degree_pressure = 1;
  constexpr double T = 1e-2;
  // constexpr double deltat = 1e-4;
  constexpr double nu = 0.01;

  EthierSteinman problem(mesh_file_name, degree_velocity, degree_pressure, T,
                         deltat, nu);

  problem.setup();
  problem.solve();

  return 0;
}