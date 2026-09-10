#include "../include/NavierStokes.hpp"
#include "../include/ConfigHandler.hpp"

#include <vector>
#include <string>
#include <filesystem>

namespace fs = std::filesystem;

// Main function.
int main(int argc, char *argv[])
{
    Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

    MPI_Comm mpi_comm = MPI_COMM_WORLD;

    int size;
    MPI_Comm_size(mpi_comm, &size);

    int rank;
    MPI_Comm_rank(mpi_comm, &rank);

    // Check command line arguments
    if (argc != 3)
    {
        if (rank == 0)
        {
            std::cerr << "ERROR! Usage: " << argv[0] << " <config_file_path> <output_dir>" << std::endl;
        }
        return -1;
    }

    std::string configFile = argv[1];
    std::string outputDir = argv[2];

    // Variables to store configuration parameters
    std::string mesh_file_name;
    unsigned int degree_velocity;
    unsigned int degree_pressure;
    InflowRegime regime;
    double peak_velocity;
    Preconditioner preconditioner;
    double T;
    double dt;

    try
    {
        parseConfigFile(configFile, mesh_file_name, degree_velocity, degree_pressure,
                        regime, peak_velocity, preconditioner, T, dt);
    }
    catch (const std::exception &e)
    {
        if (rank == 0)
        {
            std::cerr << "Error parsing config file: " << e.what() << std::endl;
        }
        return -1;
    }

    dealii::Timer timer;
    timer.restart();

    if (mesh_file_name == "../mesh/Navier_Stokes_2D_fine.msh" || mesh_file_name == "../mesh/Navier_Stokes_2D_coarse.msh")
    {
        // 2D
        NavierStokes<2> problem(mesh_file_name, degree_velocity, degree_pressure,
                                T, dt, peak_velocity, regime, preconditioner, outputDir);

        problem.setup();
        problem.run();

        timer.stop();

        if (rank == 0)
        {
            std::cout << "Time taken to solve Navier Stokes problem: "
                      << timer.wall_time() << " seconds" << std::endl;

            const std::string output_filename = (fs::path(outputDir) / "results.csv").string();

            // --- Validate vector sizes ---
            const size_t num_iterations = problem.vec_drag_force.size();
            if (num_iterations == 0 ||
                num_iterations != problem.vec_drag_coeff.size() ||
                num_iterations != problem.vec_lift_coeff.size() ||
                num_iterations != problem.time_preconditioning.size() ||
                num_iterations != problem.time_solve.size())
            {
                std::cerr << "Error: Vector sizes are inconsistent or empty." << std::endl;
                return -1;
            }

            // --- Open file ---
            std::ofstream outputFile(output_filename);
            if (!outputFile.is_open())
            {
                std::cerr << "Error: Failed to open " << output_filename << std::endl;
                return -1;
            }

            // --- Set floating-point precision (e.g., 10 digits) ---
            outputFile << std::scientific << std::setprecision(10);

            // --- Write header ---
            outputFile << "Iteration,Coeff Drag,CoeffLift,time prec,time solve\n";

            // --- Write data ---
            for (size_t i = 0; i < num_iterations; ++i)
            {
                outputFile << i << ","
                           << problem.vec_drag_coeff[i] << ","
                           << problem.vec_lift_coeff[i] << ","
                           << problem.time_preconditioning[i] << ","
                           << problem.time_solve[i] << "\n";
            }

            // --- Check for write errors ---
            if (!outputFile)
            {
                std::cerr << "Error: Failed to write to " << output_filename << std::endl;
                return -1;
            }

            std::cout << "Results written to " << output_filename << std::endl;
        }
    }
    else if (mesh_file_name == "../mesh/Navier_Stokes_3D_fine.msh" || mesh_file_name == "../mesh/Navier_Stokes_3D_coarse.msh")
    {
        // 3D
        NavierStokes<3> problem(mesh_file_name, degree_velocity, degree_pressure,
                                T, dt, peak_velocity, regime, preconditioner, outputDir);

        problem.setup();
        problem.run();

        timer.stop();

        if (rank == 0)
        {
            std::cout << "Time taken to solve Navier Stokes problem: "
                      << timer.wall_time() << " seconds" << std::endl;

            const std::string output_filename = (fs::path(outputDir) / "results.csv").string();

            // --- Validate vector sizes ---
            const size_t num_iterations = problem.vec_drag_force.size();
            if (num_iterations == 0 ||
                num_iterations != problem.vec_drag_coeff.size() ||
                num_iterations != problem.vec_lift_coeff.size() ||
                num_iterations != problem.time_preconditioning.size() ||
                num_iterations != problem.time_solve.size())
            {
                std::cerr << "Error: Vector sizes are inconsistent or empty." << std::endl;
                return -1;
            }

            // --- Open file ---
            std::ofstream outputFile(output_filename);
            if (!outputFile.is_open())
            {
                std::cerr << "Error: Failed to open " << output_filename << std::endl;
                return -1;
            }

            // --- Set floating-point precision (e.g., 10 digits) ---
            outputFile << std::scientific << std::setprecision(10);

            // --- Write header ---
            outputFile << "Iteration,Coeff Drag,CoeffLift,time prec,time solve\n";

            // --- Write data ---
            for (size_t i = 0; i < num_iterations; ++i)
            {
                outputFile << i << ","
                           << problem.vec_drag_coeff[i] << ","
                           << problem.vec_lift_coeff[i] << ","
                           << problem.time_preconditioning[i] << ","
                           << problem.time_solve[i] << "\n";
            }

            // --- Check for write errors ---
            if (!outputFile)
            {
                std::cerr << "Error: Failed to write to " << output_filename << std::endl;
                return -1;
            }

            std::cout << "Results written to " << output_filename << std::endl;
        }
    }
    else
    {
        if (rank == 0)
        {
            std::cerr << "Error: Unsupported mesh file " << std::endl;
            std::cerr << "Supported files: Navier_Stokes_2D_fine.msh, Navier_Stokes_2D_coarse.msh, Navier_Stokes_3D_fine.msh, Navier_Stokes_3D_coarse.msh.msh" << std::endl;
        }
        return -1;
    }

    return 0;
}