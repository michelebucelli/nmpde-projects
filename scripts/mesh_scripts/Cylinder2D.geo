// Mesh size parameter
lc = 0.05; // Adjust this value to refine or coarsen the mesh
Printf("mesh2D: lc = %g", lc);

// Define the outer rectangle points
Point(1) = {0, 0, 0, 1.5 * lc};
Point(2) = {2.2, 0, 0, 1.5 * lc};
Point(3) = {2.2, 0.41, 0, 1.5 * lc};
Point(4) = {0, 0.41, 0, 1.5 * lc};

// Define the lines of the rectangle
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

// Define the circle (cylinder) points
Point(5) = {0.25, 0.20, 0, 0.65*lc};
Point(6) = {0.20, 0.25, 0, 0.65*lc};
Point(7) = {0.15, 0.20, 0, 0.65*lc};
Point(8) = {0.20, 0.15, 0, 0.65*lc};
Point(9) = {0.20, 0.20, 0, 0.65*lc};  // Center point of the circle
Point(10) = {1.0, 0.20, 0, 1.5*lc}; // Coarser mesh far from the cylinder
Point(11) = {2.0, 0.20, 0, 1.5*lc}; // Coarser mesh far from the cylinder

// Define circle arcs to form the cylinder boundary
Circle(5) = {5,9,6};
Circle(6) = {6,9,7};
Circle(7) = {7,9,8};
Circle(8) = {8,9,5};

// Create line loops for the outer boundary and the hole (cylinder)
Line Loop(1) = {1,2,3,4};       // Outer rectangle
Line Loop(2) = {5,6,7,8};       // Circle (cylinder)

// Define the plane surface with the cylinder as a hole
Plane Surface(1) = {1,2};

// Define physical groups for boundary conditions
Physical Line(0) = {4};         // Inlet
Physical Line(1) = {2};         // Outlet
Physical Line(2) = {1,3};       // Walls
Physical Line(3) = {5,6,7,8};   // Cylinder
Physical Surface(4) = {1};      // Fluid domain


// You can refine the mesh near the cylinder if needed
