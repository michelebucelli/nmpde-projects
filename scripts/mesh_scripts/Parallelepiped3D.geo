

// Gmsh script to define a 3D flow domain with an inner parallelepiped obstacle
// without using the OpenCASCADE kernel

// Mesh size parameter
lc = 0.05; // Adjust this value to refine or coarsen the mesh
Printf("mesh3D-Paralleliped: lc = %g", lc);

// Domain dimensions
H = 0.41;
L = 2.5;

// Define corner points of the outer domain
Point(1) = {0, 0, 0, lc};
Point(2) = {0, 0, H, lc};
Point(3) = {0, H, H, lc};
Point(4) = {0, H, 0, lc};
Point(5) = {L, 0, 0, lc};
Point(6) = {L, 0, H, lc};
Point(7) = {L, H, H, lc};
Point(8) = {L, H, 0, lc};


// Define corner points of the inner obstacle
x_start = 0.45;
x_end = 0.55;
y_start = 0.15;
y_end = 0.25;

Point(9) = {x_start, y_start, 0, lc};
Point(10) = {x_start, y_start, H, lc};
Point(11) = {x_start, y_end, H, lc};
Point(12) = {x_start, y_end, 0, lc};
Point(13) = {x_end, y_start, 0, lc};
Point(14) = {x_end, y_start, H, lc};
Point(15) = {x_end, y_end, H, lc};
Point(16) = {x_end, y_end, 0, lc};


// Define lines for the outer domain in anticlockwise ordering // BLUE
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,5};

Line(9) = {1,5};
Line(10) = {2,6};
Line(11) = {7,3};
Line(12) = {8,4};

// Define lines for the inner obstacle
Line(101) = {9,10};
Line(102) = {10,11};
Line(103) = {11,12};
Line(104) = {12,9};

Line(105) = {13,14};
Line(106) = {14,15};
Line(107) = {15,16};
Line(108) = {16,13};

Line(109) = {10,14};
Line(110) = {15,11};
Line(111) = {9,13};
Line(112) = {16,12};


// Define surfaces for the fluid domain, taking into account the obstacle

// Inlet face: x=0 (no hole)
Line Loop(13) = {1,2,3,4};
Plane Surface(1) = {13};

// Outflow face: x=L (no hole)
Line Loop(14) = {5,6,7,8};
Plane Surface(2) = {14};

// Left face: z=0, with hole for obstacle
Line Loop(15) = {4,9,-8,12};              // Outer boundary of bottom face
Line Loop(16) = {104,111,-108,112};        // Inner boundary (obstacle base), reversed order
Plane Surface(3) = {15,16};               // Bottom face with hole

// Right face: z=H, with hole for obstacle
Line Loop(17) = {10,6,11,-2};             // Outer boundary of top face
Line Loop(18) = {109,106,110,-102};      // Inner boundary (obstacle top), normal order
Plane Surface(4) = {17,18};              // Top face with hole (inner loop reversed)

// Bottom face y=0, (no hole)
Line Loop(19) = {1,10,-5,-9};            // Outer boundary of side face y=0
Plane Surface(5) = {19};                 // Side face y=0 with no hole

// Side face y=H, (no hole)
Line Loop(20) = {-3,-11,7,12};           // Outer boundary of side face y=H
Plane Surface(6) = {20};                 // Side face y=H with no hole

// Obstacle front face at x = x_start
Line Loop(21) = {101,102,103,104};       // Corrected line sequence
Plane Surface(7) = {21};

// Obstacle back at x = x_end
Line Loop(22) = {105,106,107,108};     // Corrected line sequence with reversed lines
Plane Surface(8) = {22};

// Obstacle top face at y = y_end
Line Loop(23) = {-103,-110,107,112};     // Corrected line sequence
Plane Surface(9) = {23};

// Obstacle bottom face at y = y_start
Line Loop(24) = {101,109,-105,-111};     // Corrected line sequence
Plane Surface(10) = {24};

// Define the Surface Loop for the fluid domain, including obstacle surfaces with negative signs
Surface Loop(1) = {1,2,3,4,5,6,-7,-8,-9,-10};

// Define the fluid volume
Volume(1) = {1};

// Define Physical Groups for boundary conditions
// Define Physical Groups for boundary conditions
// Surfaces:
Physical Surface(0) = {1};          // Inlet
Physical Surface(1) = {2};          // Outlet
Physical Surface(2) = {3, 4, 5, 6}; // Walls
Physical Surface(3) = {7, 8, 9, 10}; // Obstacle
// Volume:
Physical Volume(4) = {1};            // Fluid domain
