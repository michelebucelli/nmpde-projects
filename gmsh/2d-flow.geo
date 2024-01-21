
domain_width = 2.2;
domain_height = 0.15 + 0.1 + 0.16;

circle_radius = 0.1 / 2;
circle_x = 0.15 + circle_radius;
circle_y = 0.15 + circle_radius;

h = 1 / 80;

Point(1) = {circle_x, circle_y, 0, h};
Point(2) = {circle_x + circle_radius, circle_y, 0, h};
Point(3) = {circle_x - circle_radius, circle_y, 0, h};

Point(4) = {0, 0, 0, h};
Point(5) = {domain_width, 0, 0, h};
Point(6) = {domain_width, domain_height, 0, h};
Point(7) = {0, domain_height, 0, h};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 2};

Line(3) = {4, 5};
Line(4) = {5, 6};
Line(5) = {6, 7};
Line(6) = {7, 4};

Line Loop(11) = {2, 1};
Line Loop(12) = {3, 4, 5, 6};

Plane Surface(14) = {12, 11};
Physical Surface(1) = {14};