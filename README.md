### Current issues
- Non negligible error for Ethier Steinman problem after the first step
- Building the aYoshida preconditioner is extremely slow
- Missing PCD and aPCD preconditioners
- Missing inner preconditioners
- Preconditioner code could be improved in readability and possibly efficiency
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
- **Generate two simple cubes**: run `bash generate_cubes.sh <fe_size>`, where `<fe_size>` is the size of the finite elements. If no argument is provided, the default value is 0.02. Two cubes are generated. The first will be centered on the origin, with a side size of 2. The second will have positive vertices, and one of the vertices will be at the origin. The side length will be 2 pi. The boundaries are the same as the 3D problem, minus the cylinder.
- **Generate a variety of meshes**: run `bash generate_all.sh <factor>`, where `<factor>` is a multiplicative factor for the size of the finite elements. For each problem type (two dimensional, three dimensional and simple cubes), four meshes are generated with different levels of detail. The default value for `<factor>` is 0.1, which is multiplied by all the standard values for the generated meshes. This script is a fast way to get some meaningful meshes for testing.

Please note that running `generate_all.sh` will preventively delete all the meshes with the word `flow` or `cube` in their name. This is to avoid accumulating too many meshes in the `mesh` folder, which would make it difficult to find the ones you need.