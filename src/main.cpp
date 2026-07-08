#include "../include/CardiacProblem.hpp"
#include <iostream>

int main(int argc, char *argv[])
{
    try
    {
        using namespace dealii;
        using namespace CardiacProject; 

        Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

        // Impostiamo il default: Modalità Benchmark Interna
        MeshMode mode = MeshMode::InternalGeneration;
        std::string mesh_file = "";

        // Se da terminale hai scritto qualcosa dopo il nome del programma
        if (argc > 1) {
            mode = MeshMode::ExternalLoad;
            mesh_file = argv[1]; // Prendi il nome del file dal terminale
        }

        // Passiamo tutto al costruttore aggiornato
        CardiacProblem<3> cardiac_simulator(mode, mesh_file);
        cardiac_simulator.run();
    }
    catch (std::exception &exc)
    {
        std::cerr << std::endl
                  << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        std::cerr << "Exception on processing: " << std::endl
                  << exc.what() << std::endl
                  << "Aborting!" << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        return 1;
    }
    catch (...)
    {
        std::cerr << std::endl
                  << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        std::cerr << "Unknown exception!" << std::endl
                  << "Aborting!" << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        return 1;
    }
    return 0;
}
