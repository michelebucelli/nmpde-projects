# Navier-Stokes Solver 
## Authors
- __Lorenzo Terzi__  Master's Degree student in High-Performance Computing Engineering at Politecnico di Milano\
Email: [lorenzo1.terzi@mail.polimi.it](mailto:lorenzo1.terzi@mail.polimi.it)

- __Daniele Cursano__ Master's Degree student in High-Performance Computing Engineering at Politecnico di Milano\
Email: [daniele.cursano@mail.polimi.it](mailto:daniele.cursano@mail.polimi.it)

- __Marco Melzi__ Master's Degree student in High-Performance Computing Engineering at Politecnico di Milano\
Mail: [marco2.melzi@mail.polimi.it](mailto:marco2.melzi@mail.polimi.it)

This project was developed for the course Numerical Method for Partial Differential Equations\
Professor: Alfio Maria Quarteroni\
Assistant Professor: Michele Bucelli\
Politecnico di Milano

## Prerequisites
### Deal.II
The project is based on the finite element library [Deal.II](https://www.dealii.org/). To install it, follow the instructions on the [official website](https://www.dealii.org/current_release/download/).

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
python mesh/mesh_gen_2d.py
python mesh/mesh_gen_3d.py
``` 
### Elect the test case to run
There are in total 6 configuration files related to the tests presented in the benchmark paper.
The filename is: config_(dim)D(test_case_number).txt
Select the test you want to perform by command line (see later).
Lines that start with `#` are considered comments and are ignored. By this file the following parameters can be modified:
- `mesh_file_name`: path to the mesh file (can be changed)
- `degree_velocity`: degree of the velocity space
- `degree_pressure`: degree of the pressure space
- `regime`: regime of the inflow
- `peak_velocity`: U_m as peak veloocity in the test
- `preconditioner`: preconditioner adopted (can be changed)
- `T`: final time of the simulation
- `dt`: time step

### Compiling
To build the executable, make sure you have loaded the needed modules with
```bash
$ module load gcc-glibc dealii
```
Then run the following commands:
```bash
$ mkdir build
$ cd build
$ mkdir output
$ cmake ..
$ make
```
The executable will be created into `build`, and can be executed through
```bash
$ ./navier-stokes
```

### Executing
After having generated the executable the user can start the computation by execute from the build directory the command:
```
./navier-stokes <config_file> <output_dir>
```


or
```
mpirun ./navier-stokes <config_file> <output_dir>
```
to execute it in parallel using MPI.

**Arguments:**
- `<config_file>`: Path to the selected configuration file containing the simulation parameters (e.g., `../test/config_2D1.txt`).
- `<output_dir>`: Path to the directory where the resulting output files will be written (e.g., `./output`).

### Outputs
During the execution the output is generated inside the output folder containing the result of the computation.
