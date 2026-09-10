import gmsh
import sys
import os

# ==============================================================================
# INITIALIZE THE 2D MESH 
# ==============================================================================
gmsh.initialize()
gmsh.model.add("DFG_Navier_Stokes_2D")

# ==============================================================================
# DEFINE THE GEOMETRY PARAMETERS 
# ==============================================================================
L = 2.2 # [m] Length of the channel
H = 0.41 # [m] Height of the channel 
x_circle = 0.2 # [m] x-coordinate of the circle center
y_circle = 0.2 # [m] y-coordinate of the circle center
D = 0.1 # [m] Diameter of the circle
R = D / 2 # [m] Radius of the circle

# ==============================================================================
# DEFINE THE POINTS
# ==============================================================================
''' Rectangle '''
p1 = gmsh.model.geo.addPoint(0, 0, 0)
p2 = gmsh.model.geo.addPoint(L, 0, 0)
p3 = gmsh.model.geo.addPoint(L, H, 0)
p4 = gmsh.model.geo.addPoint(0, H, 0)

''' Circle '''
p_center_circle = gmsh.model.geo.addPoint(x_circle, y_circle, 0)
p_right_circle = gmsh.model.geo.addPoint(x_circle + R, y_circle, 0)
p_left_circle = gmsh.model.geo.addPoint(x_circle - R, y_circle, 0)

# ==============================================================================
# DEFINE THE LINES
# ==============================================================================
''' Rectangle '''
l_bottom = gmsh.model.geo.addLine(p1,p2)
l_outlet = gmsh.model.geo.addLine(p2,p3)
l_top = gmsh.model.geo.addLine(p3,p4)
l_inlet = gmsh.model.geo.addLine(p4,p1)

''' Circle '''
arc_top = gmsh.model.geo.addCircleArc(p_right_circle,p_center_circle,p_left_circle)
arc_bottom = gmsh.model.geo.addCircleArc(p_left_circle,p_center_circle,p_right_circle)

# ==============================================================================
# BUILD THE SURFACE
# ==============================================================================
external_loop = gmsh.model.geo.addCurveLoop([l_bottom, l_outlet, l_top, l_inlet])
circle_loop = gmsh.model.geo.addCurveLoop([arc_top, arc_bottom])

surface = gmsh.model.geo.addPlaneSurface([external_loop, circle_loop])

gmsh.model.geo.synchronize()

# ==============================================================================
# DEFINE THE PHYSICAL GROUPS
# ==============================================================================
ID_INLET = 1
ID_OUTLET = 2
ID_WALLS = 3
ID_CIRCLE = 4
ID_DOMAIN = 10

# a) BCs (1 dim (lines))
gmsh.model.addPhysicalGroup(1, [l_inlet], ID_INLET, name="Inlet")
gmsh.model.addPhysicalGroup(1, [l_outlet], ID_OUTLET, name="Outlet")
gmsh.model.addPhysicalGroup(1, [l_top, l_bottom], ID_WALLS, name="Walls")
gmsh.model.addPhysicalGroup(1, [arc_top, arc_bottom], ID_CIRCLE, name="Circle")

# b) Surface (2 dim)
gmsh.model.addPhysicalGroup(2, [surface], ID_DOMAIN, name="Domain")

# ==============================================================================
# ADVANCED MESH REFINING
# ==============================================================================
# a) Field Distance from the Circle
circle_distance = gmsh.model.mesh.field.add("Distance")
gmsh.model.mesh.field.setNumbers(circle_distance, "CurvesList", [arc_top, arc_bottom])
gmsh.model.mesh.field.setNumber(circle_distance, "Sampling", 100)

# b) Applying Threshold Function
threshold = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(threshold, "IField", circle_distance)
gmsh.model.mesh.field.setNumber(threshold, "DistMin", 0.001)       # Stay ultra-finitesimal up to 1mm
gmsh.model.mesh.field.setNumber(threshold, "DistMax", 0.25)        # Smooth transition up to 25cm

# c) Set the Field as Background Mesh
gmsh.model.mesh.field.setAsBackgroundMesh(threshold)

gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

# ==============================================================================
# MULTI-MESH GENERATION LOOP
# ==============================================================================
script_dir = os.path.dirname(os.path.abspath(__file__))

mesh_configs = [
    {"name": "Coarse", "cl_channel": 0.05,  "cl_circle": 0.01,  "filename": "Navier_Stokes_2D_coarse.msh"},
    {"name": "Fine",   "cl_channel": 0.025, "cl_circle": 0.01, "filename": "Navier_Stokes_2D_fine.msh"},
]

for cfg in mesh_configs:
    print(f"\n---> Generating [{cfg['name']}] Mesh...")
    
    # Update field parameters dynamically
    gmsh.model.mesh.field.setNumber(threshold, "LcMin", cfg["cl_circle"])
    gmsh.model.mesh.field.setNumber(threshold, "LcMax", cfg["cl_channel"])

    # Generate 2D mesh
    gmsh.model.mesh.generate(2)

    # Write file to the script directory
    output_path = os.path.join(script_dir, cfg["filename"])
    gmsh.write(output_path)
    print(f"Saved: {output_path}")

# ==============================================================================
# FINALIZE
# ==============================================================================
gmsh.finalize()