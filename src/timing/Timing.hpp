#pragma once

#include <mpi.h>
#include <limits>
#include <string>
#include <unordered_map>
#include <iostream>
#include <iomanip>


struct Timer {
    double total = 0.0;
    double start = 0.0;
    
    int calls = 0;

    double min = std::numeric_limits<double>::max();
    double max = 0.0;

    void tic() {
        start = MPI_Wtime(); // Start timer
    }

    void toc() {
        double end = MPI_Wtime();
        total += end - start; // Update total time
        ++calls;
        min = std::min(min, end - start); // update min
        max = std::max(max, end - start); // update max
    }
};

class Profiler {
public:
    void tic(const std::string& name);
    void toc(const std::string& name);
    void report(MPI_Comm) const;

private:
    std::unordered_map<std::string, Timer> timers;
};