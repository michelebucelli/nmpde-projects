# Project Source

This folder contains the WaveFEM solver implementation, split by execution
model while keeping the simple course-laboratory structure.

- `sequential/`: the baseline solver. `WaveEquation.hpp` declares the solver,
  `WaveEquation.cpp` implements assembly, Newmark time stepping, diagnostics and
  output, and `main.cpp` parses parameters and starts a run.
- `mpi/`: the lightweight distributed-memory version. `WaveEquationMPI.hpp` and
  `WaveEquationMPI.cpp` contain the distributed solver, while `main_mpi.cpp`
  provides its command-line entry point.

The executables are `wave_fem` for the sequential solver and `wave_fem_mpi` for
the MPI solver. Both use the same manufactured free/forced wave problem.
