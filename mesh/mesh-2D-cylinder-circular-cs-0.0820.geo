// ============================================================
// 2D channel flow around a circle - Gmsh .geo template
// Domain: [0, 2.2] x [0, 0.41]
// Circle: center (0.2, 0.2), radius 0.05
// Physical groups:
//   1 "Inlet"    - left wall  (x = 0)
//   2 "Outlet"   - right wall (x = 2.2)
//   3 "Boundary" - top + bottom walls
//   4 "Circle"   - cylinder boundary
//   5 "Domain"   - 2D fluid domain
// ============================================================

// --- Mesh size parameters ---
h_far    = 0.082;  // element size far from circle (background)
h_circle = 0.041;  // element size on circle boundary

// --- Domain corners ---
Point(1) = {0.0,  0.0,  0, h_far};
Point(2) = {2.2,  0.0,  0, h_far};
Point(3) = {2.2,  0.41, 0, h_far};
Point(4) = {0.0,  0.41, 0, h_far};

// --- Rectangle edges ---
Line(1) = {1, 2};  // bottom
Line(2) = {2, 3};  // right  (Outlet)
Line(3) = {3, 4};  // top
Line(4) = {4, 1};  // left   (Inlet)

// --- Circle (cylinder cross-section) ---
// Center point
Point(5) = {0.2, 0.2, 0, h_circle};
// Cardinal points on circumference
Point(6) = {0.25, 0.2,  0, h_circle};  // right
Point(7) = {0.2,  0.25, 0, h_circle};  // top
Point(8) = {0.15, 0.2,  0, h_circle};  // left
Point(9) = {0.2,  0.15, 0, h_circle};  // bottom

// Four circular arcs (counter-clockwise)
Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};

// --- Curve loops ---
Curve Loop(1) = {1, 2, 3, 4};  // outer rectangle 
Curve Loop(2) = {5, 6, 7, 8};  // circle 

// --- Surface (domain minus circle) ---
Plane Surface(1) = {1, -2};  // -2 means circle is a hole

// --- Physical groups ---
Physical Curve("Inlet",    1) = {4};
Physical Curve("Outlet",   2) = {2};
Physical Curve("Boundary", 3) = {1, 3};
Physical Curve("Circle",   4) = {5, 6, 7, 8};
Physical Surface("Domain", 5) = {1};

// --- Meshing options ---
Mesh.Algorithm = 6;     // Frontal-Delaunay algorithm
Mesh.RecombineAll = 0;  // keep triangles
Mesh.CharacteristicLengthMin = 0;
Mesh.CharacteristicLengthMax = h_far;
