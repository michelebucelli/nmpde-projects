// To better understand how to create a good mesh file for a stokes flow simulation,
// I'm going to create a "square inside a square". The exterior square will be the  
// domain boundary, and the interior square will be the obstacle.

// I think I understand that these properties can be achieved by creating a mesh
// that represents the domain boundary, and then adding a "hole" to the mesh that
// represents the obstacle.

// This section creates the domain boundary.
Point(1) = {-1, -1, 0, 0.125};
Point(2) = {1, -1, 0, 0.125};
Point(3) = {1, 1, 0, 0.125};
Point(4) = {-1, 1, 0, 0.125};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(6) = {1, 2, 3, 4};
Plane Surface(6) = {6};

// This section creates a smaller "obstacle" square.
Point(5) = {-0.5, -0.5, 0.0, 0.125};
Point(6) = {0.5, -0.5, 0.0, 0.125};
Point(7) = {0.5, 0.5, 0.0, 0.125};
Point(8) = {-0.5, 0.5, 0.0, 0.125};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

Line Loop(7) = {8, 5, 6, 7};
Plane Surface(7) = {7};

// The first curve loop defines the exterior boundary of the surface; all other curve loops define holes in the surface.
// Source: https://gmsh.info/doc/texinfo/gmsh.html#t5
Plane Surface(8) = {6, 7};
Physical Surface(1) = {8};