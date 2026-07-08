// ============================================================
// 3D channel flow around a cylinder - Gmsh .geo template
// Domain: [0, 2.5] x [0, 0.41] x [0, 0.41]
// Circle: center (0.5, 0.2), radius 0.05
// Physical groups:
//   1 "Inflow"   
//   2 "Outflow"   
//   3 "Wall"
//   4 "Cylinder"   
//   5 "Domain"   
// ============================================================

// --- Mesh size parameters ---
h_far    = 0.082;  // element size far from circle (background)
h_cyl    = 0.041;  // element size on circle boundary

// --- Domain parameters ---
L = 2.5;
H = 0.41;
W = 0.41;

// --- Cylinder parameters ---
cx = 0.5; 
cy = 0.2;  
r  = 0.05; 

// --- Domain corners ---
Point(1) = {0.0, 0.0, 0.0, h_far};
Point(2) = {L,   0.0, 0.0, h_far};
Point(3) = {L,   H,   0.0, h_far};
Point(4) = {0.0, H,   0.0, h_far};

// --- Rectangle edges ---
Line(1) = {1, 2};  // bottom
Line(2) = {2, 3};  // right (Outlet)
Line(3) = {3, 4};  // top
Line(4) = {4, 1};  // left  (Inlet)

// --- Cylinder basis ---
Point(5) = {cx,   cy,   0.0, h_cyl};  // central point
Point(6) = {cx+r, cy,   0.0, h_cyl};  // right
Point(7) = {cx,   cy+r, 0.0, h_cyl};  // top
Point(8) = {cx-r, cy,   0.0, h_cyl};  // left
Point(9) = {cx,   cy-r, 0.0, h_cyl};  // bottom

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

// --- 3D extrusion ---
// Extrude surface 1 along Z axis by a quantity W = 0.41
ext[] = Extrude {0, 0, W} {
  Surface{1};
};

// Mapping of ext[] returned by Extrude:
//   -  ext[0] = new upper surface (Z = 0.41)
//   -  ext[1] = new volume (3D fluid domain)
//   -  ext[2] = surface extruded from Line 1 (lower wall, Y = 0)
//   -  ext[3] = surface extruded from Line 2 (outflow, X = 2.5)
//   -  ext[4] = surface extruded from Line 3 (upper wall, Y = 0.41)
//   -  ext[5] = surface extruded from Line 4 (inflow, X = 0)
//   - ext[6], ext[7], ext[8], ext[9] = surfaces extruded from circular arcs (cylinder walls)

// --- Physical groups ---
Physical Surface("Inflow",   1) = {ext[5]};
Physical Surface("Outflow",  2) = {ext[3]};
Physical Surface("Wall",     3) = {1, ext[0], ext[2], ext[4]};  // Z=0, Z=0.41, Y=0, Y=0.41
Physical Surface("Cylinder", 4) = {ext[6], ext[7], ext[8], ext[9]};
Physical Volume("Domain",    5) = {ext[1]};

// --- Meshing options ---
Mesh.Algorithm3D = 10;  // HXT algorithm
Mesh.RecombineAll = 0;  // keep tetrahedra
Mesh.CharacteristicLengthMin = 0;
Mesh.CharacteristicLengthMax = h_far;