
/* domain_width = 2.2;
domain_height = 0.15 + 0.1 + 0.16;

circle_radius = 0.1 / 2;
circle_x = 0.15 + circle_radius;
circle_y = 0.15 + circle_radius; */

h = 1 / 4;

Point(1) = {0, 0, 0, h};
Point(2) = {5, 0, 0, h};
Point(3) = {-5, 0, 0, h};

Point(4) = {-10, -10, 0, h};
Point(5) = {10, -10, 0, h};
Point(6) = {-10, 10, 0, h};
Point(7) = {10, 10, 0, h};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 2};

Line(3) = {5, 4};
Line(4) = {4, 6};
Line(5) = {6, 7};
Line(6) = {7, 5};

Line Loop(11) = {2, 1};
Line Loop(12) = {4, 5, 6, 3};

Plane Surface(14) = {12, 11};
Physical Surface(1) = {14};