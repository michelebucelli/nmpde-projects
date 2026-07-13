# Space-Time Adaptivity

## Project Information

- **Course:** Numerical Methods for Partial Differential Equations
- **Authors:** Giacomo Maglio, Francesca Morciano, Giuseppe Hares
- **Advisor:** Alfio Maria Quarteroni
- **Co-Advisor:** Michele Bucelli
- **Academic year:** 2025/2026
- **Institution:** Politecnico di Milano

This project contains two deal.II-based solvers for the heat equation:

- `H_Heat`: homogeneous reference solver
- `STA_Heat`: adaptive space-time solver

The code is organized under `src/` and builds with CMake.

## Repository layout

- `src/homogeneous/`: reference solver implementation
- `src/adaptive/`: adaptive solver implementation
- `src/timing/`: profiling utilities

## Requirements

You need a working C++ toolchain and deal.II with MPI support.

Environment setup:

```bash
$ apptainer shell *.sif
$ source /u/sw/etc/bash.bashrc
$ module load gcc-glibc dealii
```

## Build

From the repository root:

```bash
$ mkdir build
$ cd build
$ cmake ..
$ make
```

This creates the executables inside `build/`.

## Run

Available executables:

- `build/H_Heat`
- `build/STA_Heat`

Run them from the build directory, for example:

```bash
$ cd build
$ mpirun -n <number of threads> ./H_Heat  # to run the homogenous solver
$ mpirun -n <number of threads> ./STA_Heat # to run the adaptive solver
```

Both executables currently use hardcoded parameters in their respective `main.cpp` files.
To change time step, final time, or forcing term, edit:

- `src/homogeneous/main.cpp`
- `src/adaptive/main.cpp`

## Output

The solvers write output files in the build directory.

## Notes

- The adaptive executable contains an optional baseline-comparison path behind the `COMPARE_WITH_BASE` macro.
- If you want to compare adaptive and homogeneous results, enable that macro in `src/adaptive/main.cpp` and ensure the baseline mesh is available.
