#include <deal.II/base/convergence_table.h>
#include <getopt.h>

#include "Cylinder.hpp"
#include "EthierSteinman.hpp"
#include "NavierStokes.hpp"
#include "Step.hpp"

// Main function.
int main(int argc, char* argv[]) {
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  int problem_id = 0;
  int precondition_id = 0;
  double deltat = 0.0;
  double T = 0;
  std::string mesh_file_name;
  preconditioner_id preconditioner = ASIMPLE;

  ConditionalOStream pcout(
      std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0);

  // User flags are --problem-id -p, --deltat -t and --mesh-file -m.
  // There's also an optional flag --help -h and --convergence-check -c.
  std::string err_msg =
      "Usage: " + std::string(argv[0]) +
      " -problem-id <id> -deltat <deltat> -mesh-file <file> \n" +
      "  -P, --problem-id <id>      Problem ID (1, 2, 3 or 4)\n" +
      "                             1: 2D cylinder\n"
      "                             2: 3D cylinder\n"
      "                             3: Ethier-Steinman\n"
      "                             4: Step\n" +
      "  -p, --precondition-id <id> Problem ID (1, 2, 3, 4 or 5)\n" +
      "                             1: Block diagonal\n"
      "                             2: SIMPLE\n"
      "                             3: aSIMPLE\n"
      "                             4: Yoshida\n"
      "                             5: aYoshida\n" +
      "  -T, --end-time <T>         End of the resolution time range\n" +
      "  -t, --deltat <deltat>      Length of a time step\n" +
      "  -m, --mesh-file <file>     Mesh file name\n" +
      "  -h, --help                 Display this message\n" +
      "  -c, --convergence-check    Check convergence (performs only one "
      "step of the Ethier-Steinman problem)\n ";

  bool convergence_check = false;
  const char* const short_opts = "P:p:T:t:m:h:c";
  const option long_opts[] = {
      {"problem-id", required_argument, nullptr, 'P'},
      {"precondition-id", required_argument, nullptr, 'p'},
      {"end-time", required_argument, nullptr, 'T'},
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
      case 'P':
        problem_id = std::stoi(optarg);
        break;

      case 'p':
        precondition_id = std::stoi(optarg);
        break;

      case 'T':
        T = std::stod(optarg);
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
        convergence_check = true;
        break;

      case '?':
        pcout << err_msg << std::endl;
        return 1;

      default:
        pcout << err_msg << std::endl;
        return 1;
    }
  }

  // Check if all required parameters are provided.
  if (problem_id == 0 || precondition_id == 0 || T == 0.0 || deltat == 0.0 ||
      mesh_file_name.empty()) {
    pcout << err_msg << std::endl;
    return 1;
  }

  constexpr unsigned int degree_velocity = 2;
  constexpr unsigned int degree_pressure = 1;
  constexpr unsigned int maxit = 10000;
  constexpr double tol = 1e-6;  // Relative tolerance.

  // Get the correct preconditioner.
  switch (precondition_id) {
    case 1:
      preconditioner = BLOCK_DIAGONAL;
      break;
    case 2:
      preconditioner = SIMPLE;
      break;
    case 3:
      preconditioner = ASIMPLE;
      break;
    case 4:
      preconditioner = YOSHIDA;
      break;
    case 5:
      preconditioner = AYOSHIDA;
      break;
    default:
      pcout << err_msg << std::endl;
      return 1;
  }

  const SolverOptions solver_options(maxit, tol, preconditioner, 1.0);

  // Run the chosen problem.
  switch (problem_id) {
    case 1: {
      Cylinder2D problem(mesh_file_name, degree_velocity, degree_pressure, T,
                         deltat, solver_options);
      problem.setup();
      problem.solve();

      problem.update_lift_drag();
      break;
    }

    case 2: {
      Cylinder3D problem(mesh_file_name, degree_velocity, degree_pressure, T,
                         deltat, solver_options);
      problem.setup();
      problem.solve();

      problem.update_lift_drag();
      pcout << "Lift coefficient: " << problem.get_lift() << std::endl;
      pcout << "Drag coefficient: " << problem.get_drag() << std::endl;
      break;
    }

    case 3: {
      if (convergence_check) {
        pcout << "Convergence check is being performed" << std::endl;
        pcout << "===================================" << std::endl;
        pcout << "Please note that the provided problem ID is ignored"
              << std::endl;
        pcout << "and we're defaulting to problem 3 (Ethier-Steinman)"
              << std::endl;
        pcout << "===================================" << std::endl;

        std::vector<std::string> mesh_factors = {"0.8", "0.4", "0.2", "0.1"};
        // We're setting T to deltat so that only one time step is performed.
        constexpr double nu = 0.01;

        std::vector<double> h_values;
        std::vector<double> errors_L2;
        std::vector<double> errors_H1;
        for (auto& mesh_factor : mesh_factors) {
          std::string mesh_full_name =
              mesh_file_name + "-" + mesh_factor + ".msh";

          EthierSteinman problem(mesh_full_name, degree_velocity,
                                 degree_pressure, deltat, deltat,
                                 solver_options, nu);

          problem.setup();
          problem.solve();

          h_values.emplace_back(std::stod(mesh_factor));

          errors_L2.emplace_back(
              problem.compute_error(VectorTools::L2_norm, false));

          errors_H1.emplace_back(
              problem.compute_error(VectorTools::H1_norm, true));
        }

        const unsigned int mpi_rank =
            Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
        if (mpi_rank == 0) {
          ConvergenceTable table;

          std::ofstream convergence_file("convergence.csv");
          convergence_file << "h,eL2,eH1" << std::endl;

          for (size_t i = 0; i < mesh_factors.size(); i++) {
            table.add_value("h", h_values[i]);
            table.add_value("L2", errors_L2[i]);
            table.add_value("H1", errors_H1[i]);

            convergence_file << h_values[i] << "," << errors_L2[i] << ","
                             << errors_H1[i] << std::endl;
          }

          table.evaluate_all_convergence_rates(
              ConvergenceTable::reduction_rate_log2);

          table.set_scientific("L2", true);
          table.set_scientific("H1", true);

          table.write_text(std::cout);
        }
      }

      else {
        constexpr double nu = 0.01;
        EthierSteinman problem(mesh_file_name, degree_velocity, degree_pressure,
                               T, deltat, solver_options, nu);
        problem.setup();
        problem.solve();
      }
      break;
    }

    case 4: {
      constexpr double alpha = 1.0;
      Step problem(mesh_file_name, degree_velocity, degree_pressure, T, deltat,
                   solver_options, alpha);
      problem.setup();
      problem.solve();
      break;
    }

    default:
      pcout << err_msg << std::endl;
      return 1;
  }

  return 0;
}