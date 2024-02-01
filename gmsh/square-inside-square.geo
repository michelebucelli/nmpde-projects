// To better understand how to create a good mesh file for a stokes flow simulation,
// I'm going to create a "square inside a square". The exterior square will be the  
// domain boundary, and the interior square will be the obstacle. I will later
// extend this to a 3D simulation with a cube inside a cube.

// I think I understand that these properties can be achieved by creating a mesh
// that represents the domain boundary, and then adding a "hole" to the mesh that
// represents the obstacle.
square_size = 2.0;

// The size of the hole that represents the obstacle, with respect to the size of the domain.
square_ratio = 0.25;

// This value is used to define the characteristic length of the mesh elements.
// The characteristic length is used to define the size of the mesh elements, so
// if we want a finer mesh, we can decrease this value.
characteristic_length = 1 / 8;

// This section creates the domain boundary.
Point(1) = {-square_size / 2, -square_size / 2, 0, characteristic_length};
Point(2) = {+square_size / 2, -square_size / 2, 0, characteristic_length};
Point(3) = {+square_size / 2, +square_size / 2, 0, characteristic_length};
Point(4) = {-square_size / 2, +square_size / 2, 0, characteristic_length};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(6) = {1, 2, 3, 4};
Plane Surface(6) = {6};

// This section creates a smaller "obstacle" square.
hole_size = square_size * square_ratio;

Point(5) = {-hole_size / 2, -hole_size / 2, 0, characteristic_length};
Point(6) = {+hole_size / 2, -hole_size / 2, 0, characteristic_length};
Point(7) = {+hole_size / 2, +hole_size / 2, 0, characteristic_length};
Point(8) = {-hole_size / 2, +hole_size / 2, 0, characteristic_length};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

Line Loop(7) = {5, 6, 7, 8};
Plane Surface(7) = {7};

// The first curve loop defines the exterior boundary of the surface; all other curve loops define holes in the surface.
// Source: https://gmsh.info/doc/texinfo/gmsh.html#t5
Plane Surface(8) = {6, 7};
Physical Surface(1) = {8};

// Don't forget to use dimension=2 when creating the mesh.
// gmsh square-inside-square.geo -2 && gmsh square-inside-square.msh