### PROJECT 3: NAVIER-STOKES EQUATIONS
Group members:

* Cenzato Matteo
* Pisante Giuseppe
* Procaccio Arianna

For an history of the commits refer to this [repository](https://github.com/mattecenz/3-NAVIER-STOKES-Cenzato-Pisante-Procaccio).

It was made public only recently (on the 31/01/2024), if it needs to stay private I will change it as soon as possible.

### SOME INFORMATIONS FOR THE TESTS

As instructed no mesh files were uploaded but only the .geo files used to generate them. For the 3D mesh in order to put the customized boundary it is a .geo_unrolled but it should not make any difference.

To run the construction of the mesh the command is 

```
./gmsh cilinder_2D.geo
// OR
./gmsh cilinder_3D.geo_unrolled
```

These outputs the files called "cilinder_XD_fine.msh".

When compiling (standard methodology) two executables are generated, in order to run them the command is 

```
./navier_stokesXD mesh_name.msh
```

Where X can be either 2 or 3, depending on the dimension of the mesh. The mesh name is optional if the file created above ( in `../mesh/`) is present.
