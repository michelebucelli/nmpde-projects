# Navier-Stokes Solver 
This project implements a high-performance simulator for solving steady and unsteady Navier-Stokes equations using the finite element method (FEM) with parallelization for computational efficiency. The main branch provides solvers capable of handling both time-independent and time-dependent fluid dynamics problems, specifically focusing on the classical flow past a cylinder scenario in both 2D and 3D. The implementation includes a monolithic solver with various preconditioners to enhance performance and an uncoupled solver using a second-order incremental projection method. The code structure consists of three main classes: `MonolithicNavierStokes.cpp`, `SteadyNavierStokes.cpp`, and `UncoupledNavierStokes.cpp`, each addressing different numerical strategies for solving the equations. Additionally, the simulator provides outputs compatible with ParaView and computes essential flow characteristics such as lift and drag coefficients. This project aims to offer a scalable and flexible approach to solving Navier-Stokes equations, making it suitable for benchmarking and performance analysis.
## Authors
- __Maria Aurora Bertasini__  Master's Degree student in High-Performance Computing Engineering at Politecnico di Milano\
GitHub: [aurorabertasini](https://github.com/aurorabertasini)\
Email: [mariaaurora.bertasini@mail.polimi.it](mailto:mariaaurora.bertasini@mail.polimi.it)

- __Marco Cioci__ Master's Degree student in High-Performance Computing Engineering at Politecnico di Milano\
GitHub: [MarcoCioci](https://github.com/MarcoCioci)\
Email: [marco.cioci@mail.polimi.it](mailto:marco.cioci@mail.polimi.it)

- __Leonardo Ignazio Pagliochini__ Master's Degree student in High-Performance Computing Engineering at Politecnico di Milano\
GitHub: [leonardopagliochini](https://github.com/leonardopagliochini)\
Mail: [leonardoignazio.pagliochini@mail.polimi.it](mailto:leonardoignazio.pagliochini@mail.polimi.it)

This project was developed for the course Numerical Method for Partial Differential Equations\
Professor: Alfio Maria Quarteroni\
Assistant Professor: Michele Bucelli\
Politecnico di Milano

## Branch structure
The project is divided into two main branches: `main` and `Numerical-Tests`.

In the first branch, the main focus is on the implementation of the Navier-Stokes solver to solve the problem of the fluid past a cylinder.
In the second branch, the focus is on the implementation of the numerical tests to verify the correctness of the implemented solver and the theorical results converning time and space convergence.

## Prerequisites
### Deal.II
The project is based on the finite element library [Deal.II](https://www.dealii.org/). To install it, follow the instructions on the [official website](https://www.dealii.org/current_release/download/).\

### Gmsh
The mesh is generated using the software [Gmsh](https://gmsh.info/). To install it, follow the instructions on the [official website](https://gmsh.info/).

### OpenMPI
The project uses the Message Passing Interface (MPI) for parallelization. To install it, follow the instructions on the [official website](https://www.open-mpi.org/).

### Python
The project uses Python for the generation of the mesh. To install it, follow the instructions on the [official website](https://www.python.org/).

## Getting Started
### Generate the mesh
The mesh is generated using the `gmsh` software. To generate all the meshes with the same command run the python script `generate_mesh.py` in the folder scripts.
For instance, if you are located in the root folder of the project, you can run the following command:
```bash
python scripts/generate_mesh.py
``` 
### Modify the parameters of the simulation
The parameters of the simulation can be modified in the file `parameters.config`. The file presents a custom format. Lines that start with `#` are considered comments and are ignored. By this file the following parameters can be modified:
- `mesh_2d_path`: path to the 2D mesh file
- `mesh_3d_path`: path to the 3D mesh file
- `degree_velocity`: degree of the velocity space
- `degree_pressure`: degree of the pressure space
- `T`: final time of the simulation
- `deltat`: time step
- `Re`: Reynolds number

If one of these parameters is not present in the file, it will be asked to the user at the beginning of the simulation.

### Compiling
To build the executable, make sure you have loaded the needed modules with
```bash
$ module load gcc-glibc dealii
```
Then run the following commands:
```bash
$ mkdir build
$ cd build
$ cmake ..
$ make
```
The executable will be created into `build`, and can be executed through
```bash
$ ./executable-name
```

### Executing
After having generated the executable the user can start the computation by execute from the build directory the command:
```
./main
```
or
```
mpirun ./main
```
to execute it in parallel using MPI.

### Outputs
During the execution a directory `outputs` is generated inside the build folder containing the result of the computation.
The  `outputs` is automatically divided in subfolder depending on the type of simulation executed and the reynolds number