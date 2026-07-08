# Configuration
GMSH = gmsh

# .geo files for 2D and 3D meshes
MESH_2D = mesh-2D-cylinder-circular-cs-0.0205.geo mesh-2D-cylinder-circular-cs-0.0410.geo mesh-2D-cylinder-circular-cs-0.0820.geo
MESH_3D = mesh-3D-cylinder-circular-cs-0.0205.geo mesh-3D-cylinder-circular-cs-0.0410.geo mesh-3D-cylinder-circular-cs-0.0820.geo \
		  mesh-3D-cylinder-rectangular-cs-0.0205.geo mesh-3D-cylinder-rectangular-cs-0.0410.geo mesh-3D-cylinder-rectangular-cs-0.0820.geo

# Target names
TARGETS_2D = $(MESH_2D:.geo=.msh)
TARGETS_3D = $(MESH_3D:.geo=.msh)

# Generate all meshes using make
all: $(TARGETS_2D) $(TARGETS_3D)

# Rule for 2D meshes
$(TARGETS_2D): %.msh: %.geo
	@echo "Generating 2D mesh... $<"
	$(GMSH) -2 $< -o $@

# Rule for 3D meshes
$(TARGETS_3D): %.msh: %.geo
	@echo "Generating 3D mesh... $<"
	$(GMSH) -3 $< -o $@

# Clean up generated .msh files
clean:
	@echo "Removing .msh files..."
	rm -f *.msh

.PHONY: all clean