**Numerical Methods for Partial Differential Equations**  
Politecnico di Milano, A.Y. 2025/2026

Authors:
- Micaela Perlini, 10860443
- Francesca Marina Pozzi, 10837189
- Martina Rusconi, 10857811

<br>

# Navier-Stokes equations

This project solves, by means of the *finite element method*, the unsteady, incompressible Navier–Stokes equations to simulate the benchmark problem "flow past a cylinder" by M. Schäfer and S. Turek for different values of the Reynolds number $Re \le 200$, together with the computation of the drag and lift coefficients acting on the obstacle over
time.

More information about the problem, the discretization, the methods employed, and the results can be found in the project [report](NMPDE_NavierStokes_Project_Report_PerliniPozziRusconi.pdf).

## Repository structure

```text
.
├── CMakeLists.txt
├── NMPDE_NavierStokes_Project_Report_PerliniPozziRusconi.pdf
├── README.md
├── common
│   └── cmake-common.cmake
├── mesh
│   ├── Makefile
│   ├── README.md
│   ├── mesh-2D-cylinder-circular-cs-0.0205.geo
│   ├── mesh-2D-cylinder-circular-cs-0.0410.geo
│   ├── mesh-2D-cylinder-circular-cs-0.0820.geo
│   ├── mesh-3D-cylinder-circular-cs-0.0205.geo
│   ├── mesh-3D-cylinder-circular-cs-0.0410.geo
│   ├── mesh-3D-cylinder-circular-cs-0.0820.geo
│   ├── mesh-3D-cylinder-rectangular-cs-0.0205.geo
│   ├── mesh-3D-cylinder-rectangular-cs-0.0410.geo
│   └── mesh-3D-cylinder-rectangular-cs-0.0820.geo
└── src
    ├── NavierStokes2D.cpp
    ├── NavierStokes2D.hpp
    ├── NavierStokes3D.cpp
    ├── NavierStokes3D.hpp
    ├── main2D.cpp
    ├── main3D.cpp
    └── preconditioners.hpp
```

## Mesh generation

Instructions about mesh generation can be retrieved [here](mesh/README.md). The `.geo` scripts inside `mesh/` folder relate to the geometry described by M. Schäfer and S. Turek benchmark.


## Compiling and running
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
The executables will be created into `build/`. Based on the  problem desired to solve:

 * Use
    ```bash
    $ ./main2D
    ```
    or
    ```bash
    $ mpirun -n <NPROC> ./main2D
    ```
    if the test cases to run refer to 2D M. Schäfer and S. Turek benchmark.

 * Use
    ```bash
    $ ./main3D
    ```
    or
    ```bash
    $ mpirun -n <NPROC> ./main3D
    ```
    if the test cases to run refer to 3D M. Schäfer and S. Turek benchmark.





