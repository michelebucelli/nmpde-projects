# This is a bash script that generates the meshes for the
# 2d and 3d stokes problems. Meshes are generated using
# gmsh. The meshes are saved in the mesh folder.

# We clear the mesh folder from files that contain the
# string "flow", so that we don't accidentally delete
# other meshes that we might have in the folder.
rm ../mesh/*flow*

# 2D mesh.
gmsh -2 -clmax 0.1  ../gmsh/2d-flow.geo -o ../mesh/2d-flow.0.1.msh
gmsh -2 -clmax 0.05 ../gmsh/2d-flow.geo -o ../mesh/2d-flow.0.05.msh
gmsh -2 -clmax 0.02 ../gmsh/2d-flow.geo -o ../mesh/2d-flow.0.02.msh
gmsh -2 -clmax 0.01 ../gmsh/2d-flow.geo -o ../mesh/2d-flow.0.01.msh

# 3D mesh.
gmsh -3 -clmax 0.1  ../gmsh/3d-flow.geo -o ../mesh/3d-flow.0.1.msh
gmsh -3 -clmax 0.05 ../gmsh/3d-flow.geo -o ../mesh/3d-flow.0.05.msh
gmsh -3 -clmax 0.02 ../gmsh/3d-flow.geo -o ../mesh/3d-flow.0.02.msh
gmsh -3 -clmax 0.01 ../gmsh/3d-flow.geo -o ../mesh/3d-flow.0.01.msh

# We're using the -clmax flag to specify the maximum element size,
# as a simple way to control the mesh size. Inside the .geo files,
# we specify the minimum element size as 1, which is pretty big,
# so the -clmax flag is necessary to get a reasonable mesh size.