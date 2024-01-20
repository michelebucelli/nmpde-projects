// Using this alternative geometry kernel helped reduce the clutter in mesh generation.
SetFactory("OpenCASCADE");

// This is an extension of "square-inside-square.geo" in the third dimension, in order to get
// a deeper understanding of how to define 3D shapes in gmsh. This is my last exercise before
// I begin tackling the actual problem scenario, for the navier-stokes equations.

// This value is used to define the characteristic length of the mesh elements.
// The characteristic length is used to define the size of the mesh elements, so
// if we want a finer mesh, we can decrease this value.
characteristic_length = 1 / 6;

// The cube is defined by four points, that will be extruded into a 3D shape.
// The points are defined by their x, y, and z coordinates, and the characteristic length.
// Below, we define the edge size of the cube.
cube_size = 2.0;

// The size of the hole that represents the obstacle, with respect to the size of the domain.
hole_ratio = 0.5;

// These instructions require the OpenCASCADE kernel.
hole_size = cube_size * hole_ratio;
Box(1) = {-cube_size / 2, -cube_size / 2, -cube_size / 2, cube_size, cube_size, cube_size};
Box(2) = {-hole_size / 2, -hole_size / 2, -hole_size / 2, hole_size, hole_size, hole_size};

// Boolean operations with OpenCASCADE create new entities, that's why we're adding "Delete"
// in the arguments. Source: https://gmsh.info/doc/texinfo/gmsh.html#t16
BooleanDifference(3) = { Volume{1}; Delete; }{ Volume{2}; Delete; };

// I found this workaround to set the characteristic length value, but I don't think
// we'll be using the OpenCASCADE factory in the actual navier-stokes mesh generation,
// as boundary definition appears to be more complex.
MeshSize{ PointsOf{ Volume{:}; } } = characteristic_length;
Physical Volume(1) = {3};

// This script should be used with the command
// gmsh cube-inside-cube.geo -3 -format vtk && gmsh cube-inside-cube.vtk