#include "Timing.hpp"

void Profiler::tic(const std::string& name)
{
    timers[name].tic();
}

void Profiler::toc(const std::string& name)
{
    timers[name].toc();
}

void Profiler::report(MPI_Comm comm) const
{
    int rank;
    MPI_Comm_rank(comm, &rank);

    // Only rank 0 prints final table
    struct AggTimer {
        double total = 0.0;
        double min = 0.0;
        double max = 0.0;
        int calls = 0;
    };

    if (rank == 0)
    {
        std::cout << "---------------------------------------------------------------\n";
        std::cout << std::left
                  << std::setw(18) << "Timer"
                  << std::setw(10) << "Calls"
                  << std::setw(10) << "Min(s)"
                  << std::setw(10) << "Avg(s)"
                  << std::setw(10) << "Max(s)"
                  << std::setw(10) << "Total(s)"
                  << "\n";
        std::cout << "---------------------------------------------------------------\n";
    }

    for (const auto& [name, local] : timers)
    {
        AggTimer global;

        // Reduce TOTAL
        MPI_Reduce(&local.total, &global.total, 1,
                   MPI_DOUBLE, MPI_MAX, 0, comm);

        // Reduce CALLS
        MPI_Reduce(&local.calls, &global.calls, 1,
                   MPI_INT, MPI_SUM, 0, comm);

        // Reduce MIN (best case across ranks)
        MPI_Reduce(&local.min, &global.min, 1,
                   MPI_DOUBLE, MPI_MIN, 0, comm);

        // Reduce MAX (worst case across ranks)
        MPI_Reduce(&local.max, &global.max, 1,
                   MPI_DOUBLE, MPI_MAX, 0, comm);

        if (rank == 0)
        {
            double avg = global.total / global.calls;

            std::cout << std::left
                      << std::setw(18) << name
                      << std::setw(10) << global.calls
                      << std::setw(10) << global.min
                      << std::setw(10) << avg
                      << std::setw(10) << global.max
                      << std::setw(10) << global.total
                      << "\n";
        }
    }

    if (rank == 0)
    {
        std::cout << "---------------------------------------------------------------\n";
    }
}