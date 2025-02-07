// Gmsh script to define a 3D flow domain with an inner cylindrical obstacle
// without using the OpenCASCADE kernel

// Mesh size parameter
lc = 0.1; // Adjust this value to refine or coarsen the mesh
Printf("mesh3D-Cylinder: lc = %g", lc);

// Domain dimensions
H = 0.41;
L = 2.5;

// Cylinder parameters
x_c = 0.5;    // Cylinder center x-coordinate
y_c = 0.2;    // Cylinder center y-coordinate
R   = 0.05;    // Cylinder radius

// Define corner points of the outer domain
Point(1) = {0, 0, 0, 1.75*lc};
Point(2) = {0, 0, H, 1.75*lc};
Point(3) = {0, H, H, 1.75*lc};
Point(4) = {0, H, 0, 1.75*lc};
Point(5) = {L, 0, 0, 1.75*lc};
Point(6) = {L, 0, H, 1.75*lc};
Point(7) = {L, H, H, 1.75*lc};
Point(8) = {L, H, 0, 1.75*lc};

// Define points around the cylinder at z=0 (base)
Point(9)  = {x_c + R, y_c,      0, 0.5*lc}; // Right point
Point(10) = {x_c,      y_c + R, 0, 0.5*lc}; // Top point
Point(11) = {x_c - R, y_c,      0, 0.5*lc}; // Left point
Point(12) = {x_c,      y_c - R, 0, 0.5*lc}; // Bottom point

// Define points around the cylinder at z=H (top)
Point(13) = {x_c + R, y_c,      H, 0.5*lc};
Point(14) = {x_c,      y_c + R, H, 0.5*lc};
Point(15) = {x_c - R, y_c,      H, 0.5*lc};
Point(16) = {x_c,      y_c - R, H, 0.5*lc};

// Center points for circle arcs
Point(17) = {x_c, y_c, 0, 0.5*lc}; // Center at z=0
Point(18) = {x_c, y_c, H, 0.5*lc}; // Center at z=H

// Define circle arcs at z=0 (base of cylinder)
Circle(101) = {9, 17, 10};
Circle(102) = {10, 17, 11};
Circle(103) = {11, 17, 12};
Circle(104) = {12, 17, 9};

// Define circle arcs at z=H (top of cylinder)
Circle(105) = {13, 18, 14};
Circle(106) = {14, 18, 15};
Circle(107) = {15, 18, 16};
Circle(108) = {16, 18, 13};

// Define vertical lines connecting the base and top of the cylinder
Line(109) = {9, 13};
Line(110) = {10, 14};
Line(111) = {11, 15};
Line(112) = {12, 16};

// Define lines for the outer domain
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

Line(9)  = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {8, 4};

// Define surfaces for the fluid domain, taking into account the obstacle

// Inlet face at x=0 (no hole)
Line Loop(13) = {1, 2, 3, 4};
Plane Surface(1) = {13};

// Outlet face at x=L (no hole)
Line Loop(14) = {5, 6, 7, 8};
Plane Surface(2) = {14};

// Bottom face z=0, with hole for the cylinder
Line Loop(15) = {4, 9, -8, 12};                 // Outer boundary of the bottom face
Line Loop(16) = {-104, -103, -102, -101};       // Inner boundary (hole), reversed order
Plane Surface(3) = {15, 16};                    // Bottom face with hole

// Top face z=H, with hole for the cylinder
Line Loop(17) = {2, 11, -6, -10};               // Outer boundary of the top face
Line Loop(18) = {-105, -106, -107, -108};       // Inner boundary (hole), reversed order
Plane Surface(4) = {17, 18};                    // Top face with hole

// Side face y=0 (no hole)
Line Loop(19) = {1, 10, -5, -9};
Plane Surface(5) = {19};

// Side face y=H (no hole)
Line Loop(20) = {3, -11, -7, -12};
Plane Surface(6) = {20};

// Define the side surfaces of the cylinder
Line Loop(21) = {101, 110, -105, -109};
Ruled Surface(7) = {21};

Line Loop(22) = {102, 111, -106, -110};
Ruled Surface(8) = {22};

Line Loop(23) = {103, 112, -107, -111};
Ruled Surface(9) = {23};

Line Loop(24) = {104, 109, -108, -112};
Ruled Surface(10) = {24};

// Define the Surface Loop for the fluid domain, including the cylinder surfaces (negative sign)
Surface Loop(1) = {1, 2, 3, 4, 5, 6, -7, -8, -9, -10};

// Define the fluid volume
Volume(1) = {1};

// Define Physical Groups for boundary conditions
// Surfaces:
Physical Surface(0) = {1};          // Inlet
Physical Surface(1) = {2};          // Outlet
Physical Surface(2) = {3, 4, 5, 6}; // Walls
Physical Surface(3) = {7, 8, 9, 10}; // Obstacle
// Volume:
Physical Volume(4) = {1};            // Fluid domain


