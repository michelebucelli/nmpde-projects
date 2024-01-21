// I've been trying to learn the gmsh scripting language. This is a small playground for me to test out some of the features.
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

// To add the six faces as physical surfaces, we need to define the physical
// surface of each face. This is done by defining the physical surface of each
// line that makes up the face.

// TAG #1 is BOTTOM
Physical Surface(1) = {6};

// Side face 1 (1, 2, 6, 5)
Line(101) = {1, 2};
Line(102) = {2, 6};
Line(103) = {6, 5};
Line(104) = {5, 1};

Line Loop(2) = {101, 102, 103, 104};
Plane Surface(3) = {2};

// TAG #2 is FRONT
Physical Surface(2) = {3};

// Side face 2 (2, 3, 10, 6)
Line(105) = {2, 3};
Line(106) = {3, 10};
Line(107) = {10, 6};
Line(108) = {6, 2};

Line Loop(3) = {105, 106, 107, 108};
Plane Surface(4) = {3};

// TAG #3 is RIGHT
Physical Surface(3) = {4};

// Side face 3 (3, 4, 14, 10)
Line(109) = {3, 4};
Line(110) = {4, 14};
Line(111) = {14, 10};
Line(112) = {10, 3};

Line Loop(4) = {109, 110, 111, 112};
Plane Surface(5) = {4};

// TAG #4 is BACK
Physical Surface(4) = {5};

// Side face 4 (4, 1, 5, 14)
Line(113) = {4, 1};
Line(114) = {1, 5};
Line(115) = {5, 14};
Line(116) = {14, 4};

Line Loop(105) = {113, 114, 115, 116};
Plane Surface(106) = {105};

// TAG #5 is LEFT
Physical Surface(5) = {106};

// Top face (5, 6, 10, 14)
Line(117) = {5, 6};
Line(118) = {6, 10};
Line(119) = {10, 14};
Line(120) = {14, 5};

Line Loop(107) = {117, 118, 119, 120};
Plane Surface(108) = {107};

// TAG #6 is TOP
Physical Surface(6) = {108};

// If you're wondering why how we can choose the format of the output file, or
// how to tell gmsh that we also want vertices inside the cube, keep in mind
// that these are not defined in the .geo file. These are part of the options
// that are passed to gmsh when we run the script.

// In this case, I'd use the following command to generate and see the mesh:
// gmsh cube.geo -3 && gmsh cube.vtk

// The -3 flag tells gmsh to generate a 3D mesh, with tetrahedral elements, instead of a 2D surface
// (the shell). The && is used to run two commands in sequence, so we can
// generate the mesh and then open it in gmsh's built-in viewer.

// The tag order is:
// 1. BOTTOM
// 2. FRONT
// 3. RIGHT
// 4. BACK
// 5. LEFT
// 6. TOP

// Where by "bottom" we mean the face with the lowest z coordinate, and by
// "front" we mean the face with the lowest y coordinate.