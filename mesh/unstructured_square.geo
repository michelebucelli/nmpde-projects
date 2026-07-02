Point(1) = {0, 0, 0, 0.02};
Point(2) = {1, 0, 0, 0.04};
Point(3) = {1, 1, 0, 0.02};
Point(4) = {0, 1, 0, 0.03};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Physical Surface(1) = {1};

Mesh.RecombineAll = 1;
