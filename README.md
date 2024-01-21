### Current issues
- Segmentation fault in Ethier Steinman problem (calculation of the gradient of u^n)
- Missing most mesh files
- Flow past a cyclinder was not tested (waiting for mesh files)
- Boundary tags for flow past a cyclinder are placeholders (waiting for mesh files)
- No preconditioner is used
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