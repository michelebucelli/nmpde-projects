### Current issues
- CMakeLists should be improved so that .cpp and .hpp files are automatically added as sources and includes
- Missing mesh files
- Nothing was tested (waiting for mesh files)
- Boundary tags are placeholders (waiting for mesh files)
- No preconditioner is used
- Missing computation of the requested quantities

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
$ ./main
```