# This is a bash script that generates just one mesh,
# for the 2d stokes problem. Meshes are generated using
# gmsh. The meshes are saved in the mesh folder.

# This is the standard factor that we're using to begin
# with. We can change this to get different mesh sizes.
factor=0.02

# Use the first (optional) argument to change the factor.
if [ $# -eq 1 ]
  then
    factor=$1
fi

# Format the factor with two digits after the decimal point
formatted_factor=$(printf "%.3f" $factor)

# 2D mesh generation.
gmsh -2 -clmax $factor ../gmsh/2d-flow.geo -o ../mesh/2d-flow-$formatted_factor.msh

# We're using the -clmax flag to specify the maximum element size,
# as a simple way to control the mesh size. Inside the .geo files,
# we specify the standard element size as 1, which is pretty big,
# so the -clmax flag is necessary to get a reasonable mesh size.