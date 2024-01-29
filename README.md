### Organizing the source code
Please place all your sources into the `src` folder.

Binary files must not be uploaded to the repository (including executables).

Mesh files should not be uploaded to the repository. If applicable, upload `gmsh` scripts with suitable instructions to generate the meshes (and ideally a Makefile that runs those instructions). If not applicable, consider uploading the meshes to a different file sharing service, and providing a download link as part of the building and running instructions.

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
