import gmsh
import sys
import os

# ==============================================================================
# INITIALIZE THE 3D MESH 
# ==============================================================================
gmsh.initialize()
gmsh.model.add("DFG_Navier_Stokes_3D")

# ==============================================================================
# DEFINE THE GEOMETRY PARAMETERS
# ==============================================================================
L = 2.5  # [m] Length of the channel (along X-axis)
H = 0.41 # [m] Height of the channel (along Y-axis)
W = 0.41 # [m] Width of the channel (along Z-axis)

x_cylinder = 0.50 # [m] x-coordinate of the cylinder center
y_cylinder = 0.20 # [m] y-coordinate of the cylinder center
D = 0.1 # [m] Diameter of the cylinder
R = D / 2 # [m] Radius of the cylinder

# ==============================================================================
# BUILD THE 3D GEOMETRY
# ==============================================================================
# a) Create the main fluid domain box
channel = gmsh.model.occ.addBox(0, 0, 0, L, H, W)

# b) Create the cylinder
cylinder = gmsh.model.occ.addCylinder(x_cylinder, y_cylinder, 0, 0, 0, W, R)

# c) Subtract the cylinder from the box to get the final fluid domain
domain_t, _ = gmsh.model.occ.cut([(3, channel)], [(3, cylinder)])
fluid_domain_tag = domain_t[0][1]

gmsh.model.occ.synchronize()

# ==============================================================================
# DEFINE THE PHYSICAL GROUPS
# ==============================================================================
# BCs are 2D Surfaces (dim=2), and the Domain volume is 3D (dim=3)

# Define the various IDs
ID_INLET = 1
ID_OUTLET = 2
ID_WALLS = 3
ID_CYLINDER = 4

ID_DOMAIN = 10


surfaces = gmsh.model.getBoundary([(3, fluid_domain_tag)], oriented=False)

inlet_surfaces = []
outlet_surfaces = []
wall_surfaces = []
cylinder_surfaces = []

# Loop through all boundary surfaces to isolate them based on their center position
for dim, tag in surfaces:
    x_c, y_c, z_c = gmsh.model.occ.getCenterOfMass(dim, tag)
    
    if abs(x_c - 0.0) < 1e-4:          # Surface at z = 0 is the Inflow
        inlet_surfaces.append(tag)
    elif abs(x_c - L) < 1e-4:          # Surface at z = L is the Outflow
        outlet_surfaces.append(tag)
    elif abs(y_c - 0.0) < 1e-4 or abs(y_c - H) < 1e-4 or abs(z_c - 0.0) < 1e-4 or abs(z_c - W) < 1e-4:
        wall_surfaces.append(tag)      # External walls of the channel (Top, Bottom, Front, Back)
    else:
        cylinder_surfaces.append(tag)  # The remaining internal surface is the cylinder

# Assign the Physical Groups
gmsh.model.addPhysicalGroup(2, inlet_surfaces, ID_INLET, name="Inlet")
gmsh.model.addPhysicalGroup(2, outlet_surfaces, ID_OUTLET, name="Outlet")
gmsh.model.addPhysicalGroup(2, wall_surfaces, ID_WALLS, name="Walls")
gmsh.model.addPhysicalGroup(2, cylinder_surfaces, ID_CYLINDER, name="Cylinder")


# Volume Group
gmsh.model.addPhysicalGroup(3, [fluid_domain_tag], ID_DOMAIN, name="Domain")

# ==============================================================================
# ADVANCED MESH REFINING
# ==============================================================================
# a) Field Distance from the Cylinder surface
cylinder_distance = gmsh.model.mesh.field.add("Distance")
gmsh.model.mesh.field.setNumbers(cylinder_distance, "FacesList", cylinder_surfaces)
gmsh.model.mesh.field.setNumber(cylinder_distance, "Sampling", 100)

# b) Applying Threshold Function around the cylinder
threshold = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(threshold, "IField", cylinder_distance)
gmsh.model.mesh.field.setNumber(threshold, "DistMin", 0.01)        # Up to 1cm
gmsh.model.mesh.field.setNumber(threshold, "DistMax", 0.30)        # Transition up to 30cm

# c) Set the Field as Background Mesh
gmsh.model.mesh.field.setAsBackgroundMesh(threshold)

gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

# ==============================================================================
# MULTI-MESH GENERATION LOOP (3D)
# ==============================================================================
script_dir = os.path.dirname(os.path.abspath(__file__))

mesh_configs = [
    {"name": "Coarse 3D", "cl_channel": 0.2, "cl_cylinder": 0.025,  "filename": "Navier_Stokes_3D_coarse.msh"},
    {"name": "Fine 3D",   "cl_channel": 0.1, "cl_cylinder": 0.0125, "filename": "Navier_Stokes_3D_fine.msh"},
]

for cfg in mesh_configs:
    print(f"\n---> Generating [{cfg['name']}] Mesh...")
    
    # Dynamic field update
    gmsh.model.mesh.field.setNumber(threshold, "LcMin", cfg["cl_cylinder"])
    gmsh.model.mesh.field.setNumber(threshold, "LcMax", cfg["cl_channel"])

    # Generate 3D volume mesh
    gmsh.model.mesh.generate(3)

    # Save to script folder
    output_path = os.path.join(script_dir, cfg["filename"])
    gmsh.write(output_path)
    print(f"Saved: {output_path}")

# ==============================================================================
# FINALIZE
# ==============================================================================
gmsh.finalize()