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

// This is the characteristic length of the mesh. It's used to control the
// size of the mesh elements. The smaller the characteristic length, the
// smaller the mesh elements. We're setting to 1, because we want the mesh
// elements to be controlled by the flag -clmax, which is set in the
// command line. See the script generate_mesh.sh for more details.
characteristic_length = 1;

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
Physical Surface(1) = {14};

// This code should hopefully add tags to the four lines that make up the 
// domain, and the circle.
Physical Line(1) = {3}; // (1) Bottom.
Physical Line(2) = {4}; // (2) Right.
Physical Line(3) = {5}; // (3) Top.
Physical Line(4) = {6}; // (4) Left.

// (5) Circle.
Physical Line(5) = {1, 2};