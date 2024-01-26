# This is a bash script that generates the meshes for the
# 2d and 3d stokes problems, as well as the "simple cube"
# mesh. Meshes are generated using gmsh. The meshes are 
# saved in the mesh folder.

# We clear the mesh folder from files that contain the
# string "flow", or "cube", so that we don't accidentally 
# delete other meshes that we might have in the folder.
rm ../mesh/*flow*
rm ../mesh/*cube*

# This is the standard factor that we're using to begin
# with. We can change this to get different mesh sizes.
factor=0.1

# Use the first (optional) argument to change the factor.
if [ $# -eq 1 ]
  then
    factor=$1
fi

# 2D mesh.
gmsh -2 -clmax $factor                            ../gmsh/2d-flow.geo -o ../mesh/2d-flow-factor-1.msh
gmsh -2 -clmax $(echo "$factor * 0.5" | bc -l)    ../gmsh/2d-flow.geo -o ../mesh/2d-flow-factor-0.5.msh
gmsh -2 -clmax $(echo "$factor * 0.25" | bc -l)   ../gmsh/2d-flow.geo -o ../mesh/2d-flow-factor-0.25.msh
gmsh -2 -clmax $(echo "$factor * 0.125" | bc -l)  ../gmsh/2d-flow.geo -o ../mesh/2d-flow-factor-0.125.msh

# 3D mesh.
gmsh -3 -clmax $factor                            ../gmsh/3d-flow.geo -o ../mesh/3d-flow-factor-1.msh
gmsh -3 -clmax $(echo "$factor * 0.5" | bc -l)    ../gmsh/3d-flow.geo -o ../mesh/3d-flow-factor-0.5.msh
gmsh -3 -clmax $(echo "$factor * 0.25" | bc -l)   ../gmsh/3d-flow.geo -o ../mesh/3d-flow-factor-0.25.msh
gmsh -3 -clmax $(echo "$factor * 0.125" | bc -l)  ../gmsh/3d-flow.geo -o ../mesh/3d-flow-factor-0.125.msh

# Cube mesh.
gmsh -3 -clmax $factor                            ../gmsh/cube.geo -o ../mesh/cube-factor-1.msh
gmsh -3 -clmax $(echo "$factor * 0.5" | bc -l)    ../gmsh/cube.geo -o ../mesh/cube-factor-0.5.msh
#gmsh -3 -clmax $(echo "$factor * 0.25" | bc -l)  ../gmsh/cube.geo -o ../mesh/cube-factor-0.25.msh
#gmsh -3 -clmax $(echo "$factor * 0.125" | bc -l) ../gmsh/cube.geo -o ../mesh/cube-factor-0.125.msh

# Feel free to uncomment the above two lines if you need a finer mesh for the cube.
# Lines are commented out because the cube mesh is pretty big, and it takes a while
# to generate the mesh.

# We're using the -clmax flag to specify the maximum element size,
# as a simple way to control the mesh size. Inside the .geo files,
# we specify the standard element size as 1, which is pretty big,
# so the -clmax flag is necessary to get a reasonable mesh size.