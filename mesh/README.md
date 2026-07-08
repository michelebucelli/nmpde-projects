# Instructions for mesh generation

This folder contains Gmsh scripts (`.geo`) needed to generate the computational domain and a `Makefile` to automate the process.

## Prerequisites

To generate the meshes, it is **necessary** to have **Gmsh** installed. 

This can be done by running
``` bash
sudo apt-get install gmsh
```
on Ubuntu terminal.

## How to manage meshes

1. Open terminal (do not enter the container yet).
2. Navigate to `nmpde-NavierStokes/mesh`.
3. Execute the following command:

   ```bash
   make
   ```

<br>

If you need to clean up all generated mesh files, run

```bash
make clean
```