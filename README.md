### Current issues
- Comments need improvement
- Non negligible L2 error for Ethier Steinman problem
- SIMPLE preconditioner not working (runtime error)
- Missing many ad hoc preconditioners
- 2D and 3D flow past a cylinder do not converge (missing preconditioner)
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