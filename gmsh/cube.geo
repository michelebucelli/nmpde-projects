//I've been trying to learn the gmsh scripting language. This is a small playground for me to test out some of the features.
// This script in particular defines the geometry of a 3D cube, making sure that the origin is at the center of the cube.

// This value is used to define the characteristic length of the mesh elements.
// The characteristic length is used to define the size of the mesh elements, so
// if we want a finer mesh, we can decrease this value.
characteristic_length = 0.25;

// The cube is defined by four points, that will be extruded into a 3D shape.
// The points are defined by their x, y, and z coordinates, and the characteristic length.
// Below, we define the edge size of the cube.
cube_size = 2.0;

// The points are defined here. You might be expecting eight points, but we only need four,
// since we will be extruding them into a 3D shape.
Point(1) = {-cube_size / 2, -cube_size / 2, -cube_size / 2, characteristic_length};
Point(2) = {+cube_size / 2, -cube_size / 2, -cube_size / 2, characteristic_length};
Point(3) = {+cube_size / 2, +cube_size / 2, -cube_size / 2, characteristic_length};
Point(4) = {-cube_size / 2, +cube_size / 2, -cube_size / 2, characteristic_length};

// Connect the points to form a square.
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Define the surface of the square, by connecting the lines with a loop.
// Loops are used to define the surface of a shape.
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

// Extrusion is used to define the 3D shape of the cube, by adding four more points.
tmp[] = Extrude {0.0, 0.0, cube_size} {
  Surface{6};
};

// The term "Physical" is used to define the physical properties of the shape.
// In this case, we are defining the physical volume of the cube.
Physical Volume(1) = tmp[1];

// If you're wondering why how we can choose the format of the output file, or
// how to tell gmsh that we also want vertices inside the cube, keep in mind
// that these are not defined in the .geo file. These are part of the options
// that are passed to gmsh when we run the script.

// In this case, I'd use the following command to generate and see the mesh:
// gmsh cube.geo -3 -format vtk && gmsh cube.vtk

// The -3 flag tells gmsh to generate a 3D mesh, with tetrahedral elements, instead of a 2D surface
// (the shell). The -format vtk flag tells gmsh to output the mesh in the vtk format, which is
// what we've been using in the labs. The && is used to run two commands in sequence, so we can
// generate the mesh and then open it in gmsh's built-in viewer.