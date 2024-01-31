#include <deal.II/base/convergence_table.h>
#include <getopt.h>

#include "Cylinder.hpp"
#include "EthierSteinman.hpp"
#include "NavierStokes.hpp"
#include "Step.hpp"

// Main function.
int main(int argc, char* argv[]) {
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  // Parameters with no command line option assigned.
  constexpr unsigned int degree_velocity = 2;
  constexpr unsigned int degree_pressure = 1;
  constexpr unsigned int maxit = 10000;
  constexpr unsigned int maxit_inner = 10000;

  // Default arguments.
  int problem_id = 1;
  int precondition_id = 2;
  double deltat = 0.01;
  double T = 1.0;
  double U_m = 0.0;
  bool varying_inlet = false;
  std::string mesh_file_name;
  bool convergence_check = false;
  bool use_inner_solver = true;
  double tol = 1e-7;
  double tol_inner = 1e-5;
  double alpha = 1.0;
  bool use_ilu = false;
  bool compute_lift_drag = false;

  ConditionalOStream pcout(
      std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0);

  std::string err_msg =
      "Usage: " + std::string(argv[0]) +
      " -problem-id <id> -deltat <deltat> -mesh-file <file> \n" +
      "  -P, --problem-id <id>       Problem ID (1, 2, 3 or 4)\n" +
      "                              1: 2D cylinder (default)\n"
      "                              2: 3D cylinder\n"
      "                              3: Ethier-Steinman\n"
      "                              4: Step\n" +
      "  -p, --precondition-id <id>  Preconditioner ID (1, 2, 3, 4)\n" +
      "                              1: Block diagonal\n"
      "                              2: SIMPLE (default)\n"
      "                              3: aSIMPLE\n"
      "                              4: Yoshida\n" +
      "  -T, --end-time <T>          End of the time range\n" +
      "  -t, --deltat <deltat>       Length of a time step\n" +
      "  -m, --mesh-file <file>      Mesh file name\n" +
      "  -h, --help                  Display this message\n" +
      "  -c, --convergence-check     Check convergence\n" +
      "  -d  --lift-drag             Whether to compute the lift and drag "
      "coefficients for flow past a cylinder\n" +
      "  -u, --inlet-velocity <U>    Reference inlet velocity for flow past a "
      "cylinder [m/s]\n" +
      "  -v, --varying-inlet         Use a non-constant inlet velocity in flow "
      "past a cylinder\n" +
      "  -i  --no-inner-solver       Use a single preconditioner sweep instead "
      "of an inner solver for aSIMPLE\n" +
      "  -l  --ilu-preconditioner    Use ILU as inner preconditioner instead "
      "of AMG\n" +
      "  -a  --alpha <alpha>         Value of alpha in (0, 1] for the SIMPLE "
      "or aSIMPLE preconditioners\n" +
      "  -x  --tol <tol>             Relative tolerance for the main solver\n" +
      "  -y  --tol-inner <tol-inner> Relative tolerance for the inner "
      "solvers\n.";

  const char* const short_opts = "P:p:T:t:m:hcdu:vila:x:y:";
  const option long_opts[] = {
      {"problem-id", required_argument, nullptr, 'P'},
      {"precondition-id", required_argument, nullptr, 'p'},
      {"end-time", required_argument, nullptr, 'T'},
      {"deltat", required_argument, nullptr, 't'},
      {"mesh-file", required_argument, nullptr, 'm'},
      {"help", no_argument, nullptr, 'h'},
      {"convergence-check", no_argument, nullptr, 'c'},
      {"lift-drag", no_argument, nullptr, 'd'},
      {"inlet-velocity", required_argument, nullptr, 'u'},
      {"ilu-preconditioner", no_argument, nullptr, 'l'},
      {"varying-inlet", no_argument, nullptr, 'v'},
      {"no-inner-solver", no_argument, nullptr, 'i'},
      {"alpha", required_argument, nullptr, 'a'},
      {"tol", required_argument, nullptr, 'x'},
      {"tol-inner", required_argument, nullptr, 'y'},
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

      case 'd':
        compute_lift_drag = true;
        break;

      case 'u':
        U_m = std::stod(optarg);
        break;

      case 'v':
        varying_inlet = true;
        break;

      case 'i':
        use_inner_solver = false;
        break;

      case 'l':
        use_ilu = true;
        break;

      case 'a':
        alpha = std::stod(optarg);
        break;

      case 'x':
        tol = std::stod(optarg);
        break;

      case 'y':
        tol_inner = std::stod(optarg);
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
  if (mesh_file_name.empty() || alpha <= 0 || alpha > 1) {
    pcout << err_msg << std::endl;
    return 1;
  }

  // Get the correct preconditioner and set solver options.
  preconditioner_id preconditioner = SIMPLE;
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
    default:
      pcout << err_msg << std::endl;
      return 1;
  }
  const SolverOptions solver_options(preconditioner, maxit, tol,
                                     use_inner_solver, maxit_inner, tol_inner,
                                     use_ilu, alpha);

  // Run the chosen problem.
  switch (problem_id) {
    case 1: {
      if (U_m == 0.0) U_m = 1.5;
      Cylinder2D problem(mesh_file_name, degree_velocity, degree_pressure, T,
                         deltat, U_m, varying_inlet, solver_options,
                         compute_lift_drag);
      problem.setup();
      problem.solve();
      break;
    }

    case 2: {
      if (U_m == 0.0) U_m = 2.25;
      Cylinder3D problem(mesh_file_name, degree_velocity, degree_pressure, T,
                         deltat, U_m, varying_inlet, solver_options,
                         compute_lift_drag);
      problem.setup();
      problem.solve();
      break;
    }

    case 3: {
      if (convergence_check) {
        pcout << "Convergence check is being performed" << std::endl;
        pcout << "===============================================" << std::endl;

        constexpr unsigned int file_precision = 5;
        constexpr unsigned int num_files = 4;

        // We're setting T to 4*deltat so that only 4 time steps are performed,
        // this is experimentally enough to get a precise enough result.
        // For low values of T the stopping criterion on the residual can often
        // lead to imprecise pressure solutions, especially for large meshes.
        // Over multiple time steps, the phenomenon is mitigated.
        const double T = 4.0 * deltat;
        constexpr double nu = 0.01;

        std::vector<double> h_values;
        std::vector<double> errors_L2;
        std::vector<double> errors_H1;

        for (unsigned int i = 0; i < num_files; i++) {
          // Extract the factor from the file name
          size_t pos = mesh_file_name.find_last_of('-');
          if (pos == std::string::npos) {
            std::cerr << "Invalid file format. Exiting." << std::endl;
            return 1;
          }

          std::string factorString = mesh_file_name.substr(pos + 1);

          double mesh_factor = std::stod(factorString);

          // Halve the factor
          mesh_factor /= pow(2, i);

          // Generate the new file name
          std::ostringstream new_mesh_name;
          new_mesh_name << mesh_file_name.substr(0, pos + 1) << std::fixed
                        << std::setprecision(file_precision) << mesh_factor
                        << ".msh";

          // Print the new file name
          pcout << "New file name: " << new_mesh_name.str() << std::endl;

          EthierSteinman problem(new_mesh_name.str(), degree_velocity,
                                 degree_pressure, T, deltat, solver_options,
                                 nu);

          problem.setup();
          problem.solve();

          h_values.emplace_back(mesh_factor);

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
          convergence_file << "h,p-L2,u-H1" << std::endl;

          for (unsigned int i = 0; i < num_files; i++) {
            table.add_value("h", h_values[i]);
            table.add_value("p-L2", errors_L2[i]);
            table.add_value("u-H1", errors_H1[i]);

            convergence_file << h_values[i] << "," << errors_L2[i] << ","
                             << errors_H1[i] << std::endl;
          }

          table.evaluate_all_convergence_rates(
              ConvergenceTable::reduction_rate_log2);

          table.set_scientific("p-L2", true);
          table.set_scientific("u-H1", true);

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
      constexpr double alpha_step = 1.0;
      Step problem(mesh_file_name, degree_velocity, degree_pressure, T, deltat,
                   solver_options, alpha_step);
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