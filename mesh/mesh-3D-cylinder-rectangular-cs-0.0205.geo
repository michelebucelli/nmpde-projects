// ============================================================
// 3D channel flow around a cylinder - Gmsh .geo template
// Domain: [0, 2.5] x [0, 0.41] x [0, 0.41]
// Rectangle vertices: (0.45, 0.15), (0.55, 0.15), (0.45, 0.25), (0.55, 0.25)
// Physical groups:
//   1 "Inflow"   
//   2 "Outflow"   
//   3 "Wall"
//   4 "Cylinder"   
//   5 "Domain"   
// ============================================================

// --- Mesh size parameters ---
h_far    = 0.0205;  // element size far from rectangle (background)
h_cyl    = 0.0102;  // element size on circle boundary

// --- Domain parameters ---
L = 2.5;
H = 0.41;
W = 0.41;

// --- Cylinder parameters ---
x0    = 0.45; 
y0    = 0.15;  
edge  = 0.1; 

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
Point(5) = {x0,        y0,        0.0, h_cyl};  // bottom-left
Point(6) = {x0 + edge, y0,        0.0, h_cyl};  // bottom-right
Point(7) = {x0 + edge, y0 + edge, 0.0, h_cyl};  // top-right
Point(8) = {x0,        y0 + edge, 0.0, h_cyl};  // top-left

// Four edges
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

// --- Curve loops ---
Curve Loop(1) = {1, 2, 3, 4};  // outer rectangle 
Curve Loop(2) = {5, 6, 7, 8};  // inner rectangle 

// --- Surface (domain minus inner rectangle) ---
Plane Surface(1) = {1, -2};  // -2 means inner rectangle is a hole

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
//   - ext[6], ext[7], ext[8], ext[9] = surfaces extruded from inner ractangle edges (cylinder walls)

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