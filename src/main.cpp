#include <getopt.h>

#include "EthierSteinman.hpp"
#include "Cylinder.hpp"
#include "Step.hpp"

// Main function.
int main(int argc, char* argv[]) {
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  int problem_id = 0;
  double deltat = 0.0;
  std::string mesh_file_name;
  std::string err_msg =
      "Required Arguments: \n"
      "-p problem_id(1: Cylinder2D, 2: Cylinder3D, 3: "
      "EthierSteinman, 4: Step) \n"
      "-t deltat \n"
      "-m path_to_mesh_file";

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
        mesh_file_name = optarg;
        break;
      default:
        std::cerr << err_msg << std::endl;
        return 1;
    }
  }

  // Check if all required parameters are provided.
  if (problem_id == 0 || deltat == 0.0 || mesh_file_name.empty()) {
    std::cerr << "Error: All parameters (-p, -d, -m) must be provided."
              << std::endl;
    return 1;
  }

  constexpr unsigned int degree_velocity = 2;
  constexpr unsigned int degree_pressure = 1;
  constexpr double T = 1e-2;

  // Run the chosen problem.
  if (problem_id == 1) {
    Cylinder2D problem(mesh_file_name, degree_velocity, degree_pressure, T,
                       deltat);
    problem.setup();
    problem.solve();

  } else if (problem_id == 2) {
    Cylinder3D problem(mesh_file_name, degree_velocity, degree_pressure, T,
                       deltat);
    problem.setup();
    problem.solve();

  } else if (problem_id == 3) {
    constexpr double nu = 0.01;
    EthierSteinman problem(mesh_file_name, degree_velocity, degree_pressure, T,
                           deltat, nu);
    problem.setup();
    problem.solve();

    problem.setup();
    problem.solve();

  } else {
    std::cerr << "Error: Problem ID must be 1, 2, 3 or 4." << std::endl;
    return 1;
  }

  return 0;
}