#include "../include/SteadyNavierStokes.hpp"
#include "../include/MonolithicNavierStokes.hpp"
#include "../include/UncoupledNavierStokes.hpp"
#include "../include/ConfigReader.hpp"

int main(int argc, char *argv[])
{
    Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
    const unsigned int mpi_rank =
        Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

    int mpi_size; 

    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    if (mpi_rank == 0)
    {
        std::cout << "Welcome to the Navier-Stokes solver" << std::endl;
    }

    ConfigReader configReader;

    std::filesystem::path mesh2DPath = configReader.getMesh2DPath();
    std::filesystem::path mesh3DPath = configReader.getMesh3DPath();
    int degreeVelocity = configReader.getDegreeVelocity();
    int degreePressure = configReader.getDegreePressure();
    double simulationPeriod = configReader.getSimulationPeriod();
    double timeStep = configReader.getTimeStep();
    double Re = configReader.getRe();

    int choice = 0;

    if (mpi_rank == 0)
    {                // std::cout << "Div u_star: " << div_u_star[q] << std::endl;

        std::cout << "Please choose the problem to solve:" << std::endl;
        std::cout << "(1) Steady Navier-Stokes Solver 2D" << std::endl;
        std::cout << "(2) Steady Navier-Stokes Solver 3D" << std::endl;
        std::cout << "(3) Monolithic Navier-Stokes Solver 2D" << std::endl;
        std::cout << "(4) Monolithic Navier-Stokes Solver 3D" << std::endl;
        std::cout << "(5) Uncoupled Navier-Stokesm Solver 2D" << std::endl;
        std::cout << "(6) Uncoupled Navier-Stokesm Solver 3D" << std::endl;
        std::cout << std::endl;
        std::cout << "Enter your choice: ";

        while (choice < 1 || choice > 6)
        {
            std::cin >> choice;
            if (choice < 1 || choice > 6)
            {
                std::cout << "Invalid choice. Please enter a valid choice: ";
            }
        }
    }

    MPI_Bcast(&choice, 1, MPI_INT, 0, MPI_COMM_WORLD);

    auto start = std::chrono::high_resolution_clock::now();

    switch (choice)
    {
    case 1:
    {
        if (mpi_rank == 0) std::cout << "Solving the Steady Navier-Stokesm Problem 2D" << std::endl;
        SteadyNavierStokes<2> steadyNavierStokes2D(mesh2DPath, degreeVelocity, degreePressure , Re);
        steadyNavierStokes2D.run_full_problem_pipeline();
        break;
    }
    case 2:
    {
        if (mpi_rank == 0) std::cout << "Solving the Steady Navier-Stokesm Problem 3D" << std::endl;
        SteadyNavierStokes<3> steadyNavierStokes3D(mesh3DPath, degreeVelocity, degreePressure , Re);
        steadyNavierStokes3D.run_full_problem_pipeline();
        break;
    }
    case 3:
    {
        MonolithicNavierStokes<2> monolithicNavierStokes(mesh2DPath, degreeVelocity, degreePressure, simulationPeriod, timeStep, Re);
        monolithicNavierStokes.run();
        break;
    }
    case 4:
    {
        MonolithicNavierStokes<3> monolithicNavierStokes(mesh3DPath, degreeVelocity, degreePressure, simulationPeriod, timeStep, Re);
        monolithicNavierStokes.run();
        break;
    }
    case 5:
    {
        UncoupledNavierStokes<2> uncoupledNavierStokes(mesh2DPath, degreeVelocity, degreePressure, simulationPeriod, timeStep, Re);
        uncoupledNavierStokes.run();
        break;
    }
    case 6:
    {
        UncoupledNavierStokes<3> uncoupledNavierStokes(mesh3DPath, degreeVelocity, degreePressure, simulationPeriod, timeStep, Re);
        uncoupledNavierStokes.run();
        break; 
    }

        return 0;
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    if (mpi_rank == 0)
    {
        std::cout << "Elapsed time: " << elapsed.count() << " s\n";
        std::cout << std::endl
                  << "THE END" << std::endl;

        std::cout << "Number of Processors: " << mpi_size << std::endl;
    }
}
