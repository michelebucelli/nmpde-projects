// These are the parameters for the mesh. In particular, we're setting the
// width and height of the domain.
domain_width = 2.2;
domain_height = 0.15 + 0.1 + 0.16;

// These are the parameters for the circle. In particular, we're setting the
// radius of the circle, and the x and y coordinates of the center of the
// circle.
circle_radius = 0.1 / 2;
circle_x = 0.15 + circle_radius;
circle_y = 0.15 + circle_radius;

// We define the height of the domain above the circle.
H = 0.4;

// This is the characteristic length of the mesh. It's used to control the
// size of the mesh elements. The smaller the characteristic length, the
// smaller the mesh elements.
characteristic_length = 1 / 40;

// Oddly enough, to define a circle, gmsh requires you to list the center,
// then a point on the circle, then another point on the circle. The
// following lines define these points.
Point(1) = {circle_x, circle_y, 0, characteristic_length};
Point(2) = {circle_x + circle_radius, circle_y, 0, characteristic_length};
Point(3) = {circle_x - circle_radius, circle_y, 0, characteristic_length};

// The following lines define the vertices of the domain.
Point(4) = {0, 0, 0, characteristic_length};
Point(5) = {domain_width, 0, 0, characteristic_length};
Point(6) = {domain_width, domain_height, 0, characteristic_length};
Point(7) = {0, domain_height, 0, characteristic_length};

// Here's where I got stuck. I couldn't figure out how to define the
// circle, but I did figure out that gmsh actually generates an arc, not a
// circle. Moreoever, the arc's maximum angle is 180 degrees. So, I
// defined two arcs, each with a maximum angle of 180 degrees.
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 2};

// The following lines define the lines that make up the domain.
Line(3) = {4, 5};
Line(4) = {5, 6};
Line(5) = {6, 7};
Line(6) = {7, 4};

// The two arcs that make up the circle are combined into a single loop.
// As well as the lines that make up the domain.
Line Loop(11) = {2, 1};
Line Loop(12) = {3, 4, 5, 6};

// We're putting loop 11 after loop 12 because we want the circle to be
// carved out of the domain.
Plane Surface(14) = {12, 11};

Translate {0.0, 0.0, H} {
	Duplicata { Surface{14}; }
}

// Half-cylinder loop
Circle(1001) = {2, 1, 3};
Line(1002) = {3, 24};
Circle(1003) = {26, 25, 24};
Line(1004) = {26, 2};

Line Loop(1005) = {1001, 1002, -1003, 1004};
Surface(1006) = {1005};

// Half-cylinder loop
Circle(2001) = {3, 1, 2};
Circle(2002) = {24, 25, 26};

Line Loop(2003) = {2001, -1004, -2002, -1002};
Surface(2004) = {2003};

// First face 13 6 7 17
Line(3001) = {13, 6};
Line(3002) = {6, 7};
Line(3003) = {7, 17};
Line(3004) = {17, 13};

Line Loop(3005) = {3001, 3002, 3003, 3004};
Surface(3006) = {3005};

// Second face 4 7 17 8 
Line(4001) = {4, 7};
Line(4002) = {7, 17};
Line(4003) = {17, 8};
Line(4004) = {8, 4};

Line Loop(4005) = {4001, 4002, 4003, 4004};
Surface(4006) = {4005};

// Third face 5 4 8 9
Line(5001) = {5, 4};
Line(5002) = {4, 8};
Line(5003) = {8, 9};
Line(5004) = {9, 5};

Line Loop(5005) = {5001, 5002, 5003, 5004};
Surface(5006) = {5005};

// Fourth face 6 5 9 13
Line(6001) = {6, 5};
Line(6002) = {5, 9};
Line(6003) = {9, 13};
Line(6004) = {13, 6};

Line Loop(6005) = {6001, 6002, 6003, 6004};
Surface(6006) = {6005};

// Use all surfaces to create a volume
Surface Loop(7001) = {14, 15, 5006, 3006, 4006, 6006, -1006, -2004};
Physical Volume(7002) = {7001};