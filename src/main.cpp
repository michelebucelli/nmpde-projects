#include <getopt.h>

#include "Cylinder.hpp"
#include "EthierSteinman.hpp"
#include "NavierStokes.hpp"
#include "Step.hpp"

// Main function.
int main(int argc, char* argv[]) {
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  int problem_id = 0;
  double deltat = 0.0;
  std::string mesh_file_name;

  ConditionalOStream pcout(
      std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0);

  // User flags are --problem-id -p, --deltat -t and --mesh-file -m.
  // There's also an optional flag --help -h and --convergence-check -c.
  std::string err_msg =
      "Usage: " + std::string(argv[0]) +
      " -problem-id <id> -deltat <deltat> -mesh-file <file> \n" +
      "  -p, --problem-id <id>    Problem ID (1, 2, 3 or 4)\n" +
      "                           1: 2D cylinder\n"
      "                           2: 3D cylinder\n"
      "                           3: Ethier-Steinman\n"
      "                           4: Step\n" +
      "  -t, --deltat <deltat>    Time step size\n" +
      "  -m, --mesh-file <file>   Mesh file name\n" +
      "  -h, --help               Display this message\n" +
      "  -c, --convergence-check  Check convergence (performs only one "
      "step)\n ";

  const char* const short_opts = "p:t:m:hc";
  const option long_opts[] = {{"problem-id", required_argument, nullptr, 'p'},
                              {"deltat", required_argument, nullptr, 't'},
                              {"mesh-file", required_argument, nullptr, 'm'},
                              {"help", no_argument, nullptr, 'h'},
                              {"convergence-check", no_argument, nullptr, 'c'},
                              {nullptr, no_argument, nullptr, 0}};

  // Parse the command line arguments.
  while (true) {
    const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

    if (opt == -1) {
      break;
    }

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

      case 'h':
        pcout << err_msg << std::endl;
        return 0;

      case 'c':
        pcout << "Convergence check is not implemented yet." << std::endl;
        return 0;

      case '?':
        std::cerr << err_msg << std::endl;
        return 1;

      default:
        std::cerr << err_msg << std::endl;
        return 1;
    }
  }

  // Check if all required parameters are provided.
  if (problem_id == 0 || deltat == 0.0 || mesh_file_name.empty()) {
    std::cerr << err_msg << std::endl;
    return 1;
  }

  constexpr unsigned int degree_velocity = 2;
  constexpr unsigned int degree_pressure = 1;
  constexpr double T = 1e-2;
  const PreconditionerType precondition(ASIMPLE, 1.0);

  // Run the chosen problem.
  switch (problem_id) {
    case 1: {
      Cylinder2D problem(mesh_file_name, degree_velocity, degree_pressure, T,
                         deltat, precondition);
      problem.setup();
      problem.solve();
      break;
    }

    case 2: {
      Cylinder3D problem(mesh_file_name, degree_velocity, degree_pressure, T,
                         deltat, precondition);
      problem.setup();
      problem.solve();
      break;
    }

    case 3: {
      constexpr double nu = 0.01;
      EthierSteinman problem(mesh_file_name, degree_velocity, degree_pressure,
                             T, deltat, precondition, nu);
      problem.setup();
      problem.solve();
      break;
    }

    case 4: {
      constexpr double alpha = 1.0;
      Step problem(mesh_file_name, degree_velocity, degree_pressure, T, deltat,
                   precondition, alpha);
      problem.setup();
      problem.solve();
      break;
    }

    default:
      std::cerr << err_msg << std::endl;
      return 1;
  }

  return 0;
}