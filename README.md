### Current issues
- Non negligible error for Ethier Steinman problem after the first step
- 2D flow past a cylinder requires 0 GMRES iterations for convergence: learn why
- Solving the linear system for F in the preconditioner can sometimes result in lack of convergence, find a better preconditioner for it
- Missing many ad hoc preconditioners
- Missing computation of the lift and drag coefficients
- Check if serial computations (e.g. initialization of the triangulation) or dense matrices (e.g. in assembly) cause problems for performance/memory usage, remove them is necessary

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
$ ./navier_stokes
```