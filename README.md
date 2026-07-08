# Cardiac Electrophysiology Solver

Parallel finite-element solver for cardiac electrical wave propagation, built on deal.II and MPI.

## Overview

This project solves a monodomain-type reaction-diffusion system coupled with a minimal ventricular ionic model. The code is organized around two main pieces:

- `CardiacProblem`, which manages mesh setup, DoF distribution, assembly, linear solves, time stepping, and VTU/PVTU output.
- `IonicModel`, which implements the Bueno-Orovio ventricular action potential model with SIMD-friendly evaluation.

The executable supports two mesh modes:

- internal mesh generation for a simple benchmark geometry
- external loading from a `.msh` file, such as the included `cardiac_mesh_optimized.msh`

## Repository Layout

```
├── CMakeLists.txt
├── README.md
├── cardiac_mesh_optimized.msh
├── include/
│   ├── CardiacProblem.hpp
│   └── IonicModel.hpp
└── src/
	├── CardiacProblem.cpp
	└── main.cpp
```

## Requirements

- C++17 compiler
- CMake 3.13.4 or newer
- deal.II 9.0 or newer built with MPI support
- MPI runtime such as OpenMPI or MPICH
- Trilinos, as provided through the deal.II build

The project’s `CMakeLists.txt` currently expects deal.II to be available and checks that MPI support is enabled.

## Build

If you are using the same module environment as the original setup, load the required modules first:

```bash
module load gcc-glibc dealii
```

Then configure and build:

```bash
mkdir build
cd build
cmake ..
make
```

## Run

The executable is named `cardiac_solver`.

Run with the internal benchmark mesh:

```bash
mpirun -n 4 ./cardiac_solver
```

Run with an external mesh file:

```bash
mpirun -n 4 ./cardiac_solver ../cardiac_mesh_optimized.msh
```

The program treats the first command-line argument as a mesh filename and switches to external-mesh mode.

## Output

The solver writes VTU/PVTU output for ParaView. Files are emitted periodically during the time loop and include:

- `V`, the transmembrane potential
- `ActivationTime`, the first time each DoF crosses the activation threshold

Open the `.pvtu` files in ParaView to inspect the distributed output.