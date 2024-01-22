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

// This time, I'm defining the circle points. Below, I'll define the
// center of the circle, and then the points on the circle.
Point(1) = {circle_x, circle_y, 0, characteristic_length};

Point(2) = {circle_x, circle_y + circle_radius, 0, characteristic_length};
Point(3) = {circle_x + circle_radius, circle_y, 0, characteristic_length};
Point(4) = {circle_x, circle_y - circle_radius, 0, characteristic_length};
Point(5) = {circle_x - circle_radius, circle_y, 0, characteristic_length};

// The following lines define the vertices of the domain.
Point(6) = {0, 0, 0, characteristic_length};
Point(7) = {domain_width, 0, 0, characteristic_length};
Point(8) = {domain_width, domain_height, 0, characteristic_length};
Point(9) = {0, domain_height, 0, characteristic_length};

// Note that, for some reason, using four arcs to define the circle instead
// of two arcs is the best way to make sure that the circle is meshed
// properly after extrusion. I haven't quite figured out why this is the
// case, but it took me a while to figure out, so I'm leaving this comment
// here to help anyone else who might be trying to do the same thing.
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

// The following lines define the lines that make up the domain.
Line(5) = {6, 7};
Line(6) = {7, 8};
Line(7) = {8, 9};
Line(8) = {9, 6};

// The two arcs that make up the circle are combined into a single loop.
// As well as the lines that make up the domain.
Line Loop(11) = {1, 2, 3, 4};
Line Loop(12) = {5, 6, 7, 8};

// We're putting loop 11 after loop 12 because we want the circle to be
// carved out of the domain.
Plane Surface(14) = {12, 11};

out[] = Extrude {0, 0, H} { Surface{14}; };
Physical Volume(1) = {out[0]};
