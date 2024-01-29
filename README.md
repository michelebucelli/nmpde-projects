### OInstruction to build and run the code


The brain mesh can be downloaded from this link and placed inside the foder "FK_solver/mesh":
https://drive.google.com/file/d/1PJTaHAU-kgxId5_C6Zd4nj27HRcyHZS2/view

### Compiling
To build the executable, make sure you have loaded the needed modules with
```bash
$  cd FK_solver/src/
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
