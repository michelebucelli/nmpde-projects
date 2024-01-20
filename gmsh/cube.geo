characteristic_length = 0.1;

Point(1) = {0.0, 0.0, 0.0, characteristic_length};
Point(2) = {1.0, 0.0, 0.0, characteristic_length};
Point(3) = {1.0, 1.0, 0.0, characteristic_length};
Point(4) = {0.0, 1.0, 0.0, characteristic_length};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

tmp[] = Extrude {0,0.0,1} {
  Surface{6};
};

Physical Volume(1) = tmp[1];