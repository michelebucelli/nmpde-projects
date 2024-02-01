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

### Running
The executable will be created into `build`, and can be executed through
```bash
$ ./navier_stokes <arguments>
```
or through
```bash
$ mpirun -n <N> navier_stokes <arguments>
```
which will run the executable using MPI with `N` processes.
The executable supports multiple command line arguments, of which some are mandatory. Here is a list of all the arguments:
- **`-P, --problem-id <id>`**: specifies the problem to solve. The possible values are:
    - `1`: 2D problem with a circular hole. This is a flat plate with a circular hole in the left half of it. The inflow is on the left side.
    - `2`: 3D problem with a cylindrical hole. This is a cuboid with a cylindrical hole in the left half of it. The inflow is on the left side.
    - `3`: The Ethier-Steinman problem. This problem has a cube domain, and its exact solution is known. Therefore, we can use it to test the correctness of our implementation.
    - `4`: This is the Step problem (very similar to one of the problems tackled in the laboratory sessions).
- **`-p, --precondition-id <id>`**: specifies the preconditioner to use. The preconditioner is used to solve the linear system of equation that arises from the linearization step of the algorithm. The possible <id> values are:
    - `1`: Block diagonal preconditioner.
    - `2`: SIMPLE preconditioner.
    - `3`: aSIMPLE preconditioner.
    - `4`: Yoshida preconditioner. 
- **`-T, --end-time <T>`**: specifies the end time of the simulation, or the overall length of the simulation. The default value is 1. If it is not a multiple of `<deltat>`, the last simulation time will exceed `<T>`.
- **`-t, --delta-t <deltat>`**: specifies the time step of the simulation. The lower the value, the more accurate the simulation will be. The default value is 0.01. Some preconditioners have better convergence behaviour for lower values of this parameter, e.g. 1e-3 or 1e-4.
- **`-m --mesh-file <file>`**: specifies the mesh file to use. This argument is mandatory. Be careful to use the correct mesh file for the problem you're trying to solve, as the mesh files are specific to each problem. Mesh file generators for the Step problem are not provided, as mesh files from the laboratory sessions can be used.
- **`-h --help`**: prints the help message and exits.
- **`-c --convergence-check`**: enables the convergence check. If this argument is provided, the program will compute for each mesh factor h, the H1 norm of the error on the velocity and the L2 norm of the error on the pressure after 4 time steps, and will print them on screen. This is useful to check if the solution is converging with the expected rate, to help verify corrected of the solution. This option is only compatible with the Ethier-Steinman problem. This option requires a specific format for the mesh file name passed to the option `-m`, which should be in the format `*-h.msh`, where `*` indicates any string, `.msh` is the file extension and `h` should be a decimal number with 5 digits after the point, indicating the maximum diameter of a finite element. For example, `cube-0.06250.msh` and `cube-107.28372.msh` are both valid file names. The program will then automatically find and use three other files, with paths in the same format but using `h/2`, `h/4` and `h/8` in place of `h`. Make sure the precision of the linear solvers is high enough when using this, as the pressure is very sensitive.
- **`-d --lift-drag`**: compute the lift and drag coefficients. This option is only relevant for the 2D and 3D flow past a cylinder problems.
- **`-u, --inlet-velocity` <u>**: specifies the reference velocity of the inflow. This parameter is only relevant for the 2D and 3D flow past a cylinder problems.
- **`-v, --varying-inlet`**: this is a boolean flag. If provided, the inflow velocity will change in time. This parameter is only relevant for the 2D and 3D flow past a cylinder problems.
- **`-i, --no-inner-solver`**: use a single sweep of the preconditioner instead of an inner solver for the aSIMPLE preconditioner. This parameter is only relevant if aSIMPLE is selected.
- **`-l, --ilu-preconditioner`**: use an ILU preconditioner instead of AMG for the inner solvers.
- **`-a, --alpha` <a>**: value of alpha for the SIMPLE or aSIMPLE preconditioners, it must be greater than 0 and at most 1, which is the default value. This parameter is only relevant if SIMPLE or aSIMPLE are selected.
- **`-x --tol` <tol>**: relative tolerance for the main solver, the default value is `1e-7`.
- **`-y --tol-inner` <tol-inner>**: relative tolerance for the inner solvers, the default value is `1e-5`. Some preconditioners might need this value to be below a certain amount in some problems.

Note that the only strictly necessary argument is `--mesh-file`. If no other argument is provided, the program will run with the default values for the other arguments.

### Results
Execution results will be generated in the `results` folder.
In particular:
    - `output*` files will contain the computed solutions to the problem.
    - `lift_drag.csv` will contain the computed lift and drag coefficients. This file is only updated if the 2D or 3D flow past a cylinder problems are run.

### Mesh Generation
To generate meshes for the simulations, we use the [Gmsh](http://gmsh.info/) software. The meshes are generated from the `.geo` files in the `gmsh` folder. To generate meshes, make sure you're in the `scripts` folder. All generated meshes will be saved in the `mesh` folder. To generate a mesh, run the following commands:

- **Generate a mesh for a the 2D problem**: run `bash generate_2d.sh <fe_size>`, where `<fe_size>` is the size of the finite elements. If no argument is provided, the default value is 0.02. The resulting mesh will consist in a rectangle with a circular hole positioned in the left half of it. The boundaries are tagged as follows:
    - 1: bottom
    - 2: right
    - 3: top
    - 4: left
    - 5: circle
- **Generate a mesh for a the 3D problem**: run `bash generate_3d.sh <fe_size>`, where `<fe_size>` is the size of the finite elements. If no argument is provided, the default value is 0.02. The resulting mesh will consist in a cuboid with a cylindrical hole positioned in the left half of it. The boundaries are tagged as follows:
    - 1: bottom
    - 2: front
    - 3: right
    - 4: back
    - 5: left
    - 6: top
    - 7: cylinder
- **Generate a simple cube**: run `bash generate_cube.sh <fe_size>`, where `<fe_size>` is the size of the finite elements. If no argument is provided, the default value is 0.1. A cube centered on the origin with a side size of 2 will be generated. The boundaries are the same as the 3D problem, minus the cylinder.
- **Generate a variety of meshes**: run `bash generate_all.sh <factor>`, where `<factor>` is a multiplicative factor for the size of the finite elements. For each problem type (two dimensional, three dimensional and simple cube), four meshes are generated with different levels of detail. The default value for `<factor>` is 0.1, which is divided by two each time the level of detail is increases. This script is a fast way to get some meaningful meshes for testing.

Please note that running `generate_all.sh` will preventively delete all the meshes with the word `flow` or `cube` in their name. This is to avoid accumulating too many meshes in the `mesh` folder, which would make it difficult to find the ones you need.

### Data Analysis
The `scripts` folder contains a Python script to plot the lift and drag coefficients over time. To use it, make sure you have the `matplotlib` and `numpy` packages installed, and that the mk modules are not loaded. Then, run the following command from the `scripts` folder:
```bash
python3 plot_lift_drag.py <left-cut> <right-cut>
```

This will plot the lift and drag coefficients, provided that the simulation has been run. After running the simulation, a file named `lift_drag.csv` will be created in the `results` folder. This file contains the lift and drag coefficients over time, and is used by the Python script to plot the coefficients.

The two arguments are optional, and specify the amount of values to ignore at the beginning and at the end of the simulation. This is useful because those samples are often not representative of the actual lift and drag coefficients. The default values are 0 for both arguments.

### Generating the report
The report is generated using LaTeX. To generate the report, make sure you have the `pdflatex` command installed, and run the following command from the `report` folder:
```bash
pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex
```