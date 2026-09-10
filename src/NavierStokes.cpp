#include "../include/NavierStokes.hpp"
/**
 * @file NavierStokes.cpp
 * @brief Implementation of the NavierStokes class template for solving incompressible flow problems.
 *
 * @details
 * This file contains the implementation of all methods for the NavierStokes class template.
 * It solves the incompressible Navier-Stokes equations using the Finite Element Method (FEM)
 * with deal.II. The solver supports:
 * - 2D and 3D spatial dimensions
 * - Steady and unsteady flows
 * - Taylor-Hood elements (P2 for velocity, P1 for pressure)
 * - MPI parallel computation
 * - Multiple preconditioners
 * - Force and pressure difference calculations
 */

// -----------------------------------------------------------------------
// SETUP METHOD
// -----------------------------------------------------------------------

/**
 * @brief Initializes the entire finite element system for the simulation.
 *
 * @details
 * This method performs all one-time initialization required before the time-stepping loop:
 * 1. **Mesh initialization**: Loads the mesh from file, partitions it across MPI processes
 * 2. **FE space setup**: Creates finite element spaces for velocity and pressure
 * 3. **DoF distribution**: Distributes degrees of freedom and reorders them (velocity first, then pressure)
 * 4. **Quadrature rules**: Initializes quadrature for volume and face integrals
 * 5. **Linear system setup**: Creates sparsity patterns and initializes all matrices/vectors
 *
 * @tparam dim Spatial dimension (2 or 3)
 */
template <unsigned int dim>
void NavierStokes<dim>::setup()
{
    // ======================================================================
    // MESH INITIALIZATION
    // ======================================================================
    {
        pcout << "Initializing the mesh" << std::endl;

        // Create a serial triangulation to read the mesh file
        Triangulation<dim> mesh_serial;

        // GridIn is deal.II's class for reading mesh files
        GridIn<dim> grid_in;
        grid_in.attach_triangulation(mesh_serial);

        // Open and read the mesh file in Gmsh format
        std::ifstream grid_in_file(mesh_file_name);
        grid_in.read_msh(grid_in_file);

        // Partition the serial mesh for parallel computation
        // This distributes cells across MPI processes
        GridTools::partition_triangulation(mpi_size, mesh_serial);

        // Create a description of the partitioned mesh
        // This allows creating a fully distributed triangulation
        const auto construction_data = TriangulationDescription::Utilities::
            create_description_from_triangulation(mesh_serial, MPI_COMM_WORLD);

        // Create the fully distributed triangulation
        mesh.create_triangulation(construction_data);

        pcout << "  Number of elements = " << mesh.n_global_active_cells() << std::endl;

        // Smallest cell diameter across all locally owned cells, reduced
        // globally over MPI. Used later for the CFL number diagnostic.
        double local_min_diameter = std::numeric_limits<double>::max();
        for (const auto &cell : mesh.active_cell_iterators())
        {
            if (cell->is_locally_owned())
                local_min_diameter = std::min(local_min_diameter, cell->diameter());
        }
        mesh_h_min = Utilities::MPI::min(local_min_diameter, MPI_COMM_WORLD);
        pcout << "  Minimum cell diameter      = " << mesh_h_min << std::endl;
    }

    pcout << "-----------------------------------------------" << std::endl;

    // ======================================================================
    // FINITE ELEMENT SPACE INITIALIZATION
    // ======================================================================
    {
        pcout << "Initializing the finite element space" << std::endl;

        // Create scalar FE spaces for velocity and pressure components
        // Using SimplexP (Lagrange on simplices) for better geometry handling
        const FE_SimplexP<dim> fe_scalar_velocity(degree_velocity);
        const FE_SimplexP<dim> fe_scalar_pressure(degree_pressure);

        // Create a mixed FE space: [velocity_1, velocity_2, ..., velocity_dim, pressure]
        // This uses Taylor-Hood elements (P2 for velocity, P1 for pressure)
        // which are LBB-stable for incompressible flows
        fe = std::make_unique<FESystem<dim>>(
            fe_scalar_velocity, // Velocity FE space
            dim,                // Number of velocity components (2 for 2D, 3 for 3D)
            fe_scalar_pressure, // Pressure FE space
            1                   // Number of pressure components
        );

        pcout << "  Velocity degree:           = " << fe_scalar_velocity.degree << std::endl;
        pcout << "  Pressure degree:           = " << fe_scalar_pressure.degree << std::endl;
        pcout << "  DoFs per cell              = " << fe->dofs_per_cell << std::endl;

        // Initialize quadrature rule for volume integrals
        // Uses (fe_degree + 1) points per direction for accurate integration
        quadrature = std::make_unique<QGaussSimplex<dim>>(fe->degree + 1);
        pcout << "  Quadrature points per cell = " << quadrature->size() << std::endl;

        // Initialize quadrature rule for face integrals
        quadrature_face = std::make_unique<QGaussSimplex<dim - 1>>(fe->degree + 1);
        pcout << "  Quadrature points per face = " << quadrature_face->size() << std::endl;
    }

    pcout << "-----------------------------------------------" << std::endl;

    // ======================================================================
    // DoF HANDLER INITIALIZATION
    // ======================================================================
    {
        pcout << "Initializing the DoF handler" << std::endl;

        // Initialize the DoF handler with the mesh and FE system
        dof_handler.reinit(mesh);
        dof_handler.distribute_dofs(*fe);

        // Reorder DoFs: all velocity DoFs first (block 0), then pressure DoFs (block 1)
        // This is crucial for efficient block matrix operations
        std::vector<unsigned int> block_component(dim + 1, 0); // All in block 0 by default
        block_component[dim] = 1;                              // Pressure component goes to block 1

        DoFRenumbering::component_wise(dof_handler, block_component);

        // Get DoF sets for parallel computation
        locally_owned_dofs = dof_handler.locally_owned_dofs();                        // DoFs owned by this process
        locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler); // + ghost DoFs

        // Count DoFs per block (velocity and pressure)
        std::vector<types::global_dof_index> dofs_per_block =
            DoFTools::count_dofs_per_fe_block(dof_handler, block_component);
        const unsigned int n_u = dofs_per_block[0]; // Total velocity DoFs
        const unsigned int n_p = dofs_per_block[1]; // Total pressure DoFs

        // Create block-wise DoF sets for efficient block operations
        block_owned_dofs.resize(2);
        block_relevant_dofs.resize(2);

        // Block 0: Velocity DoFs
        block_owned_dofs[0] = locally_owned_dofs.get_view(0, n_u);
        block_relevant_dofs[0] = locally_relevant_dofs.get_view(0, n_u);

        // Block 1: Pressure DoFs
        block_owned_dofs[1] = locally_owned_dofs.get_view(n_u, n_u + n_p);
        block_relevant_dofs[1] = locally_relevant_dofs.get_view(n_u, n_u + n_p);

        pcout << "  Number of DoFs: " << std::endl;
        pcout << "    velocity = " << n_u << std::endl;
        pcout << "    pressure = " << n_p << std::endl;
        pcout << "    total    = " << n_u + n_p << std::endl;
    }

    pcout << "-----------------------------------------------" << std::endl;

    // ======================================================================
    // LINEAR SYSTEM INITIALIZATION
    // ======================================================================
    {
        pcout << "Initializing the linear system" << std::endl;
        pcout << "  Initializing the sparsity pattern" << std::endl;

        // Define coupling between DoF components for the Navier-Stokes system
        // The system has a saddle-point structure:
        // [ A  B ] [u]   = [f]
        // [B^T 0 ] [p]     [g]
        // where:
        // - A: velocity-velocity block (convection + diffusion + mass)
        // - B: pressure-velocity block (gradient)
        // - B^T: velocity-pressure block (divergence)
        Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1);

        // Initialize coupling table
        for (unsigned int c = 0; c < dim + 1; ++c)
        {
            for (unsigned int d = 0; d < dim + 1; ++d)
            {
                if (c == dim && d == dim)              // pressure-pressure coupling
                    coupling[c][d] = DoFTools::none;   // No pressure-pressure terms
                else                                   // All other couplings
                    coupling[c][d] = DoFTools::always; // Full coupling
            }
        }

        // Create sparsity pattern for the main system matrix
        TrilinosWrappers::BlockSparsityPattern sparsity(block_owned_dofs, MPI_COMM_WORLD);
        DoFTools::make_sparsity_pattern(dof_handler, coupling, sparsity);
        sparsity.compress();

        // Create a separate coupling table for the pressure mass matrix
        // Only pressure-pressure terms are non-zero
        for (unsigned int c = 0; c < dim + 1; ++c)
        {
            for (unsigned int d = 0; d < dim + 1; ++d)
            {
                if (c == dim && d == dim) // pressure-pressure term
                    coupling[c][d] = DoFTools::always;
                else // All other combinations
                    coupling[c][d] = DoFTools::none;
            }
        }

        // Create sparsity pattern for pressure mass matrix
        TrilinosWrappers::BlockSparsityPattern sparsity_pressure_mass(
            block_owned_dofs, MPI_COMM_WORLD);
        DoFTools::make_sparsity_pattern(dof_handler, coupling, sparsity_pressure_mass);
        sparsity_pressure_mass.compress();

        pcout << "  Initializing the matrices" << std::endl;
        // Initialize all matrices with their sparsity patterns
        system_matrix.reinit(sparsity);               // Main system matrix
        mass_matrix.reinit(sparsity);                 // Mass matrix (time derivative)
        convection_matrix.reinit(sparsity);           // Convection matrix
        stiffness_matrix.reinit(sparsity);            // Stiffness matrix (viscosity)
        pressure_mass.reinit(sparsity_pressure_mass); // Pressure mass matrix

        pcout << "  Initializing the system right-hand side" << std::endl;
        system_rhs.reinit(block_owned_dofs, MPI_COMM_WORLD); // RHS vector

        pcout << "  Initializing the solution vector" << std::endl;
        solution_owned.reinit(block_owned_dofs, MPI_COMM_WORLD);                // Solution (owned DoFs)
        solution.reinit(block_owned_dofs, block_relevant_dofs, MPI_COMM_WORLD); // Solution (+ ghosts)
    }
}

// -----------------------------------------------------------------------
// ASSEMBLE METHOD
// -----------------------------------------------------------------------

/**
 * @brief Assembles the full Navier-Stokes system for the first time step.
 *
 * @details
 * This method computes all contributions to the linear system matrix and RHS vector:
 * - **Viscous terms**: ν∫∇u:∇v dΩ (diffusion)
 * - **Convective terms**: ∫(u·∇)u·v dΩ + 0.5∫(∇·u)u·v dΩ (Temam stabilization)
 * - **Mass matrix**: ∫u·v/Δt dΩ (time derivative)
 * - **Pressure terms**: -∫p∇·v dΩ (momentum) + ∫q∇·u dΩ (continuity)
 * - **Time derivative RHS**: ∫u_n·v/Δt dΩ
 *
 * The assembly uses standard FEM approach with numerical integration over cells.
 *
 * @tparam dim Spatial dimension (2 or 3)
 * @param time Current simulation time
 *
 * @note This is called only once at the beginning
 * @note Uses the current solution for nonlinear terms (Picard iteration)
 * @note Applies Dirichlet BCs: inlet (parabolic profile), walls/obstacle (zero)
 * @note The Temam stabilization term ensures kinetic energy conservation in the discrete system
 */
template <unsigned int dim>
void NavierStokes<dim>::assemble(const double &time)
{
    pcout << "===============================================" << std::endl;
    pcout << "Assembling the system" << std::endl;

    const unsigned int dofs_per_cell = fe->dofs_per_cell;
    const unsigned int n_q = quadrature->size();

    // FEValues object for evaluating FE functions at quadrature points
    FEValues<dim> fe_values(
        *fe, *quadrature,
        update_values | update_gradients | update_quadrature_points | update_JxW_values);

    // Local matrices and vectors for cell contributions
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);               // Full system matrix
    FullMatrix<double> cell_mass_matrix(dofs_per_cell, dofs_per_cell);          // Mass matrix
    FullMatrix<double> cell_stiffness_matrix(dofs_per_cell, dofs_per_cell);     // Stiffness
    FullMatrix<double> cell_convection_matrix(dofs_per_cell, dofs_per_cell);    // Convection
    FullMatrix<double> cell_pressure_mass_matrix(dofs_per_cell, dofs_per_cell); // Pressure mass
    Vector<double> cell_rhs(dofs_per_cell);                                     // RHS vector

    std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

    // Reset global matrices and vectors
    system_matrix = 0.0;
    mass_matrix = 0.0;
    stiffness_matrix = 0.0;
    convection_matrix = 0.0;
    system_rhs = 0.0;
    pressure_mass = 0.0;

    // Extractors for velocity and pressure components
    FEValuesExtractors::Vector velocity(0);   // First 'dim' components are velocity
    FEValuesExtractors::Scalar pressure(dim); // Last component is pressure

    // Vectors to store current solution at quadrature points
    std::vector<Tensor<1, dim>> current_velocity_values(n_q);    // Velocity vector u
    std::vector<Tensor<2, dim>> current_velocity_gradients(n_q); // Velocity gradient ∇u
    std::vector<double> current_velocity_divergence(n_q);        // Velocity divergence ∇·u

    // Loop over all locally owned cells
    for (const auto &cell : dof_handler.active_cell_iterators())
    {
        if (!cell->is_locally_owned())
            continue; // Skip cells not owned by this MPI process

        fe_values.reinit(cell); // Reinitialize FEValues for current cell

        // Reset local matrices and vectors
        cell_matrix = 0.0;
        cell_mass_matrix = 0.0;
        cell_stiffness_matrix = 0.0;
        cell_convection_matrix = 0.0;
        cell_rhs = 0.0;
        cell_pressure_mass_matrix = 0.0;

        // Evaluate current solution at all quadrature points in this cell
        fe_values[velocity].get_function_values(solution, current_velocity_values);
        fe_values[velocity].get_function_gradients(solution, current_velocity_gradients);
        fe_values[velocity].get_function_divergences(solution, current_velocity_divergence);

        // Loop over all quadrature points in the cell
        for (unsigned int q = 0; q < n_q; ++q)
        {
            // Get forcing term at this quadrature point
            Vector<double> forcing_term_loc(dim);
            forcing_term.vector_value(fe_values.quadrature_point(q), forcing_term_loc);

            // Convert to tensor for easier manipulation
            Tensor<1, dim> forcing_term_tensor;
            for (unsigned int i = 0; i < dim; ++i)
                forcing_term_tensor[i] = forcing_term_loc[i];

            // Loop over all DoFs in the cell (i = test function index)
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
                // Loop over all DoFs in the cell (j = trial function index)
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                    // =========================================================
                    // VISCOSITY TERM (Diffusion): ν ∫ ∇u : ∇v dΩ
                    // =========================================================
                    cell_stiffness_matrix(i, j) +=
                        kinematic_viscosity *
                        scalar_product(
                            fe_values[velocity].gradient(i, q), // ∇φ_i (test function)
                            fe_values[velocity].gradient(j, q)  // ∇φ_j (trial function)
                            ) *
                        fe_values.JxW(q);

                    // =========================================================
                    // TIME DERIVATIVE (Backward Euler): ∫ u^{n+1}·v / Δt dΩ
                    // =========================================================
                    cell_mass_matrix(i, j) +=
                        scalar_product(
                            fe_values[velocity].value(i, q), // φ_i
                            fe_values[velocity].value(j, q)  // φ_j
                            ) /
                        time_step_size * fe_values.JxW(q);

                    // =========================================================
                    // CONVECTIVE TERM: ∫ (u·∇)u · v dΩ
                    // Uses current solution for u (nonlinear term)
                    // =========================================================
                    cell_convection_matrix(i, j) +=
                        scalar_product(
                            fe_values[velocity].gradient(j, q) * current_velocity_values[q], // (∇φ_j)·u
                            fe_values[velocity].value(i, q)                                  // φ_i
                            ) *
                        fe_values.JxW(q);

                    // =========================================================
                    // TEMAM STABILIZATION: 0.5 ∫ (∇·u) u · v dΩ
                    // Ensures kinetic energy conservation for the discrete system
                    // =========================================================
                    cell_convection_matrix(i, j) +=
                        0.5 * current_velocity_divergence[q] *
                        scalar_product(
                            fe_values[velocity].value(i, q), // φ_i
                            fe_values[velocity].value(j, q)  // φ_j
                            ) *
                        fe_values.JxW(q);

                    // =========================================================
                    // PRESSURE TERM IN MOMENTUM EQUATION: -∫ p ∇·v dΩ
                    // =========================================================
                    cell_matrix(i, j) -=
                        fe_values[pressure].value(j, q) *      // ψ_j (pressure trial)
                        fe_values[velocity].divergence(i, q) * // ∇·φ_i (velocity test)
                        fe_values.JxW(q);

                    // =========================================================
                    // PRESSURE TERM IN CONTINUITY EQUATION: ∫ q ∇·u dΩ
                    // =========================================================
                    cell_matrix(i, j) +=
                        fe_values[pressure].value(i, q) *      // ψ_i (pressure test)
                        fe_values[velocity].divergence(j, q) * // ∇·φ_j (velocity trial)
                        fe_values.JxW(q);

                    // =========================================================
                    // PRESSURE MASS MATRIX: ∫ p q / ν dΩ (for preconditioning)
                    // =========================================================
                    cell_pressure_mass_matrix(i, j) +=
                        fe_values[pressure].value(i, q) * // ψ_i
                        fe_values[pressure].value(j, q) * // ψ_j
                        (1.0 / kinematic_viscosity) *
                        fe_values.JxW(q);
                }

                // =========================================================
                // TIME DERIVATIVE ON RHS: ∫ u_n · v / Δt dΩ
                // u_n is the solution from the previous time step
                // =========================================================
                cell_rhs(i) +=
                    scalar_product(
                        current_velocity_values[q],     // u_n
                        fe_values[velocity].value(i, q) // φ_i
                        ) *
                    fe_values.JxW(q) / time_step_size;
            }
        }

        // No Neumann boundary conditions are applied (all boundaries use Dirichlet)

        // Get global DoF indices for this cell
        cell->get_dof_indices(dof_indices);

        // Add cell contributions to global matrices/vectors
        system_matrix.add(dof_indices, cell_matrix);
        mass_matrix.add(dof_indices, cell_mass_matrix);
        convection_matrix.add(dof_indices, cell_convection_matrix);
        stiffness_matrix.add(dof_indices, cell_stiffness_matrix);
        system_rhs.add(dof_indices, cell_rhs);
        pressure_mass.add(dof_indices, cell_pressure_mass_matrix);
    }

    // Compress matrices and vectors for parallel computation
    // This performs the actual MPI communication to assemble the global matrices
    system_matrix.compress(VectorOperation::add);
    mass_matrix.compress(VectorOperation::add);
    convection_matrix.compress(VectorOperation::add);
    stiffness_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);
    pressure_mass.compress(VectorOperation::add);

    // ======================================================================
    // BUILD FINAL SYSTEM MATRIX: M + A + C(u_n) + B
    // Where:
    // - M: Mass matrix (time derivative)
    // - A: Stiffness matrix (viscosity)
    // - C(u_n): Convection matrix (nonlinear, depends on current solution)
    // - B: Pressure terms (divergence constraints)
    // ======================================================================
    system_matrix.add(1., mass_matrix);
    system_matrix.add(1., convection_matrix);
    system_matrix.add(1., stiffness_matrix);

    // ======================================================================
    // APPLY DIRICHLET BOUNDARY CONDITIONS
    // ======================================================================
    {
        std::map<types::global_dof_index, double> boundary_values;
        std::map<types::boundary_id, const Function<dim> *> boundary_functions;

        // Create mask for velocity components (exclude pressure from boundary conditions)
        std::vector<bool> velocity_mask_vec(dim + 1, true);
        velocity_mask_vec[dim] = false; // Exclude pressure component
        const ComponentMask mask_velocity(velocity_mask_vec);

        // Apply inlet boundary condition first
        // Uses parabolic velocity profile from InletVelocity class
        inlet_velocity.set_time(time);
        boundary_functions[id_inlet] = &inlet_velocity;
        VectorTools::interpolate_boundary_values(
            dof_handler, boundary_functions, boundary_values, mask_velocity);

        // Clear and apply wall/obstacle conditions
        // These will override inlet conditions at shared nodes (corners)
        boundary_functions.clear();
        Functions::ZeroFunction<dim> zero_function(dim + 1);

        boundary_functions[id_walls] = &zero_function;
        boundary_functions[id_obstacle] = &zero_function;
        VectorTools::interpolate_boundary_values(
            dof_handler, boundary_functions, boundary_values, mask_velocity);

        // Apply boundary values to the linear system
        // false = do not eliminate constrained DoFs from the matrix
        MatrixTools::apply_boundary_values(
            boundary_values, system_matrix, solution_owned, system_rhs, false);
    }
}

// -----------------------------------------------------------------------
// ASSEMBLE TIME STEP METHOD
// -----------------------------------------------------------------------

/**
 * @brief Assembles time-dependent terms for subsequent time steps.
 *
 * @details
 * This method updates the system for time steps after the first one.
 * It only recomputes the time-dependent parts:
 * - Convective term (nonlinear, depends on current solution)
 * - Time derivative term on RHS
 *
 * This is more efficient than full assembly as it reuses the constant parts
 * (viscosity, pressure terms) from the initial assembly.
 *
 * @tparam dim Spatial dimension (2 or 3)
 * @param time Current simulation time
 *
 * @note Assumes 'assemble()' has been called at least once before
 * @note Uses Backward Euler time discretization
 * @note Applies Dirichlet BCs for velocity at inlet, walls, and obstacle
 */
template <unsigned int dim>
void NavierStokes<dim>::assemble_time_step(const double &time)
{
    pcout << "===============================================" << std::endl;
    pcout << "Assembling the system" << std::endl;

    const unsigned int dofs_per_cell = fe->dofs_per_cell;
    const unsigned int n_q = quadrature->size();

    FEValues<dim> fe_values(
        *fe, *quadrature,
        update_values | update_gradients | update_quadrature_points | update_JxW_values);

    // Local matrices and vectors for cell contributions
    FullMatrix<double> cell_convection_matrix(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_mass_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

    // Remove previous convection contribution from system matrix
    system_matrix.add(-1., convection_matrix);

    // Reset matrices and vectors
    convection_matrix = 0.0;
    system_rhs = 0.0;

    FEValuesExtractors::Vector velocity(0);
    FEValuesExtractors::Scalar pressure(dim);

    // Vectors to store current and previous solution values
    std::vector<Tensor<1, dim>> current_velocity_values(n_q);
    std::vector<Tensor<2, dim>> current_velocity_gradients(n_q);
    std::vector<double> current_velocity_divergence(n_q);

    // Loop over all locally owned cells
    for (const auto &cell : dof_handler.active_cell_iterators())
    {
        if (!cell->is_locally_owned())
            continue;

        fe_values.reinit(cell);

        cell_mass_matrix = 0.0;
        cell_convection_matrix = 0.0;
        cell_rhs = 0.0;

        // Get current solution values
        fe_values[velocity].get_function_values(solution, current_velocity_values);
        fe_values[velocity].get_function_gradients(solution, current_velocity_gradients);
        fe_values[velocity].get_function_divergences(solution, current_velocity_divergence);

        // Loop over all quadrature points in the cell
        for (unsigned int q = 0; q < n_q; ++q)
        {
            // Loop over all DoFs in the cell
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                    // =========================================================
                    // CONVECTIVE TERM: ∫ (u·∇)u · v dΩ
                    // Uses current solution for u (nonlinear term)
                    // =========================================================
                    cell_convection_matrix(i, j) +=
                        scalar_product(
                            fe_values[velocity].gradient(j, q) * current_velocity_values[q],
                            fe_values[velocity].value(i, q)) *
                        fe_values.JxW(q);
                }

                // =========================================================
                // TIME DERIVATIVE ON RHS (Backward Euler):
                // ∫ u^{n+1} · v / Δt dΩ = ∫ u^n · v / Δt dΩ + ∫ (u^{n+1} - u^n) · v / Δt dΩ
                // The second term is handled by the mass matrix on LHS,
                // so we only need u^n · v / Δt on RHS
                // =========================================================
                cell_rhs(i) +=
                    scalar_product(
                        current_velocity_values[q],
                        fe_values[velocity].value(i, q)) *
                    fe_values.JxW(q) / time_step_size;
            }
        }

        cell->get_dof_indices(dof_indices);
        convection_matrix.add(dof_indices, cell_convection_matrix);
        system_rhs.add(dof_indices, cell_rhs);
    }

    // Compress matrices and vectors for parallel computation
    convection_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);
    pressure_mass.compress(VectorOperation::add);

    // Add the new convection matrix to the system matrix
    system_matrix.add(1., convection_matrix);

    // ======================================================================
    // APPLY DIRICHLET BOUNDARY CONDITIONS
    // ======================================================================
    {
        std::map<types::global_dof_index, double> boundary_values;
        std::map<types::boundary_id, const Function<dim> *> boundary_functions;

        std::vector<bool> velocity_mask_vec(dim + 1, true);
        velocity_mask_vec[dim] = false;
        const ComponentMask velocity_mask(velocity_mask_vec);

        // Apply inlet boundary condition
        inlet_velocity.set_time(time);
        boundary_functions[id_inlet] = &inlet_velocity;
        VectorTools::interpolate_boundary_values(
            dof_handler, boundary_functions, boundary_values, velocity_mask);

        // Clear and apply wall/obstacle conditions
        boundary_functions.clear();
        Functions::ZeroFunction<dim> zero_function(dim + 1);

        boundary_functions[id_walls] = &zero_function;
        boundary_functions[id_obstacle] = &zero_function;
        VectorTools::interpolate_boundary_values(
            dof_handler, boundary_functions, boundary_values, velocity_mask);

        // Apply boundary values to the linear system
        MatrixTools::apply_boundary_values(
            boundary_values, system_matrix, solution_owned, system_rhs, false);
    }
}

// -----------------------------------------------------------------------
// SOLVE TIME STEP METHOD
// -----------------------------------------------------------------------

/**
 * @brief Solves the linear system for the current time step.
 *
 * @details
 * This method solves the block linear system using GMRES with a specified preconditioner.
 * It handles the saddle-point structure of the Navier-Stokes equations and applies
 * the selected preconditioner.
 *
 * The system matrix has the block structure:
 *
 * [ A    B  ]
 * [ B^T  0  ]
 *
 * where:
 * - A = velocity-velocity block (convection + diffusion + mass)
 * - B = pressure-velocity block
 * - B^T = velocity-pressure block
 *
 * @tparam dim Spatial dimension (2 or 3)
 * @param preconditioner The preconditioner to use (YOSIDA or SIMPLE)
 *
 * @note Uses TrilinosWrappers for parallel linear algebra
 * @note Timing information is collected for performance profiling
 */
template <unsigned int dim>
void NavierStokes<dim>::solve_time_step(const Preconditioner &preconditioner)
{
    pcout << "===============================================" << std::endl;

    const unsigned int maxiter = 100000;
    const double tol = 1e-4 * system_rhs.l2_norm(); // Relative tolerance based on RHS norm
    SolverControl solver_control(maxiter, tol, true);
    SolverGMRES<TrilinosWrappers::MPI::BlockVector> solver(solver_control);

    previous_solution = solution;

    {
        dealii::Timer preconditioner_timer;
        preconditioner_timer.restart();
        dealii::Timer solve_timer;

        auto solve_with_preconditioner = [&](const auto &prec)
        {
            preconditioner_timer.stop();
            pcout << "Time taken to initialize preconditioner: "
                  << preconditioner_timer.wall_time() << " seconds" << std::endl;
            time_preconditioning.push_back(preconditioner_timer.wall_time());

            solve_timer.restart();
            solver.solve(system_matrix, solution_owned, system_rhs, prec);
            solve_timer.stop();

            pcout << "Time taken to solve Navier Stokes problem: "
                  << solve_timer.wall_time() << " seconds" << std::endl;
            time_solve.push_back(solve_timer.wall_time());
        };

        TrilinosWrappers::SparseMatrix negB;
        negB.copy_from(system_matrix.block(1, 0));
        negB *= -1.0; // negB = -B

        TrilinosWrappers::SparseMatrix B_transpose;
        B_transpose.copy_from(system_matrix.block(0, 1));
        B_transpose *= -1.0; // B_transpose = +B^T

        switch (preconditioner)
        {
        // ======================================================================
        // IDENTITY — scalar identity, not usable on this block system
        // ======================================================================
        case (Preconditioner::IDENTITY):
        {
            throw std::runtime_error(
                "Preconditioner::IDENTITY acts on a scalar TrilinosWrappers::MPI::Vector, "
                "not the BlockVector system solved here. Use BLOCK_IDENTITY instead.");
        }

        // ======================================================================
        // BLOCK IDENTITY — no-op preconditioner, useful as a GMRES baseline
        // ======================================================================
        case (Preconditioner::BLOCK_IDENTITY):
        {
            PreconditionBlockIdentity prec;
            solve_with_preconditioner(prec);
            break;
        }

        // ======================================================================
        // SIMPLE PRECONDITIONER
        // ======================================================================
        case (Preconditioner::SIMPLE):
        {
            PreconditionSIMPLE prec;
            const double alpha = 0.5; // SIMPLE relaxation factor, in (0,1]

            prec.initialize(
                system_matrix.block(0, 0), // F: velocity-velocity block
                negB,                      // -B: pressure-velocity block
                B_transpose,               // +B^T: velocity-pressure block
                solution_owned,
                alpha);

            solve_with_preconditioner(prec);
            break;
        }

        // ======================================================================
        // APPROXIMATE SIMPLE — single ILU/AMG applications, no
        // inner GMRES solves; cheaper per iteration, more GMRES iterations
        // ======================================================================
        case (Preconditioner::APPROX_SIMPLE):
        {
            PreconditionApproxSIMPLE prec;
            const double alpha = 0.5;

            prec.initialize(
                system_matrix.block(0, 0),
                negB,
                B_transpose,
                solution_owned,
                alpha);

            solve_with_preconditioner(prec);
            break;
        }

        // ======================================================================
        // YOSIDA PRECONDITIONER
        // ======================================================================
        case (Preconditioner::YOSIDA):
        {
            PreconditionYosida prec;
            prec.initialize(
                system_matrix.block(0, 0), // F: velocity-velocity block
                negB,                      // -B: pressure-velocity block
                B_transpose,               // +B^T: velocity-pressure block
                mass_matrix.block(0, 0),   // Velocity mass matrix
                solution_owned,
                time_step_size);

            solve_with_preconditioner(prec);
            break;
        }

        // ======================================================================
        // APPROXIMATE YOSIDA
        // ======================================================================
        case (Preconditioner::APPROX_YOSIDA):
        {
            PreconditionApproxYosida prec;
            prec.initialize(
                system_matrix.block(0, 0),
                negB,
                B_transpose,
                mass_matrix.block(0, 0),
                solution_owned,
                time_step_size);

            solve_with_preconditioner(prec);
            break;
        }

        // ======================================================================
        // PCD / APPROX_PCD
        // ======================================================================
        case (Preconditioner::PCD):
        {
            throw std::runtime_error(
                "Preconditioner::PCD not yet implemented");
        }
        case (Preconditioner::APPROX_PCD):
        {
            throw std::runtime_error(
                "Preconditioner::APPROX_PCD not yet implemented");
        }

        default:
            throw std::runtime_error("Preconditioner not yet implemented in solve_time_step().");
        }
    }

    pcout << "Result:  " << solver_control.last_step() << " GMRES iterations" << std::endl;

    solution = solution_owned;
}

// -----------------------------------------------------------------------
// OUTPUT METHOD
// -----------------------------------------------------------------------

/**
 * @brief Outputs the solution fields to VTU files for visualization.
 *
 * @details
 * This method writes the solution (velocity and pressure) to disk in VTU format
 * (compatible with ParaView) at the specified time step. It includes:
 * - Velocity as a vector field
 * - Pressure as a scalar field
 * - Mesh partitioning information for parallel visualization
 *
 * The output uses a PVTU record for multi-piece visualization in parallel.
 *
 * @tparam dim Spatial dimension (2 or 3)
 * @param time The current time step index
 *
 * @note Files are written to the "./output/" directory
 * @note The filename includes the dimension and time step for uniqueness
 * @note Uses DataOut for efficient parallel I/O
 */
template <unsigned int dim>
void NavierStokes<dim>::output(const unsigned int &time)
{
    pcout << "===============================================" << std::endl;

    DataOut<dim> data_out;

    // Define how to interpret the data components
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
        interpretation(dim,
                       DataComponentInterpretation::component_is_part_of_vector);
    interpretation.push_back(DataComponentInterpretation::component_is_scalar);

    // Name the data components
    std::vector<std::string> names(dim, "velocity");
    names.push_back("pressure");

    // Add solution data to the output
    data_out.add_data_vector(dof_handler, solution, names, interpretation);

    // Add mesh partitioning information for visualization
    std::vector<unsigned int> partition_int(mesh.n_active_cells());
    GridTools::get_subdomain_association(mesh, partition_int);
    const Vector<double> partitioning(partition_int.begin(), partition_int.end());
    data_out.add_data_vector(partitioning, "partitioning");

    // Build patches for output
    data_out.build_patches();

    // Write output to file
    const std::string output_file_name = "output-Navier-Stokes-" + std::to_string(dim) + "d";
    data_out.write_vtu_with_pvtu_record(
        output_dir,                    // Output directory
        output_file_name,              // Base filename
        time,                          // Time step index
        MPI_COMM_WORLD,                // MPI communicator
        numbers::invalid_unsigned_int, // No group files
        1                              // Write one file per process
    );

    pcout << "Output written to " << output_file_name << std::endl;
    pcout << "===============================================" << std::endl;
}

// -----------------------------------------------------------------------
// RUN METHOD
// -----------------------------------------------------------------------

/**
 * @brief Executes the main simulation loop from t=0 to final_time.
 *
 * @details
 * This is the main driver method that:
 * 1. Applies initial conditions (zero velocity and pressure)
 * 2. Performs the time-stepping loop:
 *    - Advances time
 *    - Assembles the system (full assembly for first step, partial for others)
 *    - Solves the linear system
 *    - Computes forces (drag and lift) on the obstacle
 *    - Outputs solution at regular intervals
 * 3. Tracks maximum drag and minimum lift coefficients
 * 4. Computes pressure difference at the final time step
 *
 * @tparam dim Spatial dimension (2 or 3)
 *
 * @note Uses Backward Euler time discretization
 * @note Outputs solution every 10 time steps
 * @note Avoids initial transient (time > 0.1) for force coefficient tracking
 */
template <unsigned int dim>
void NavierStokes<dim>::run()
{
    pcout << "===============================================" << std::endl;

    // ======================================================================
    // APPLY INITIAL CONDITION
    // ======================================================================
    {
        pcout << "Applying the initial condition" << std::endl;

        // Interpolate zero initial condition (defined in the class)
        VectorTools::interpolate(dof_handler, initial_solution_function, solution_owned);
        solution = solution_owned;

        // Output initial solution
        output(0);
        pcout << "===============================================" << std::endl;
    }

    // Variables for tracking forces and coefficients
    std::vector<double> coefficients;
    double drag_coefficient_max = std::numeric_limits<double>::lowest(); // Smallest possible double
    double lift_coefficient_min = std::numeric_limits<double>::max();    // Largest possible double
    unsigned int time_step = 0;
    double time = 0;

    // Calculate total number of time steps
    const unsigned int total_steps = static_cast<unsigned int>(
        std::round(final_time / time_step_size));

    // ======================================================================
    // TIME STEPPING LOOP
    // ======================================================================
    while (time < final_time)
    {
        // Advance time
        time += time_step_size;
        ++time_step;
        inlet_velocity.set_time(time); // Update time for inlet velocity function

        pcout << "n = " << std::setw(3) << time_step << ", t = " << std::setw(5) << time << ":" << std::flush;

        // First time step: full assembly
        if (time_step == 1)
            assemble(time);
        // Subsequent time steps: only update time-dependent terms
        else
            assemble_time_step(time);

        // Solve the linear system
        solve_time_step(this->preconditioner);

        // ==================================================================
        // CFL DIAGNOSTIC
        // ==================================================================
        // Approximate advective Courant number using the max nodal velocity
        // component and the smallest cell diameter in the mesh. This mixes
        // velocity components rather than the true |u| magnitude and ignores
        // local mesh size variation, so treat it as an order-of-magnitude
        // indicator, not an exact per-cell CFL.

        const double velocity_linfty = solution.block(0).linfty_norm();
        const double cfl = velocity_linfty * time_step_size / mesh_h_min;
        vec_cfl.push_back(cfl);

        if (cfl > 1.0)
            pcout << "  [WARNING] CFL = " << cfl << " > 1 at t = " << time << std::endl;

        // Compute pressure difference at the final time step
        if (time_step == total_steps - 1)
            compute_pressure_difference();

        // Compute forces (drag and lift) on the obstacle
        coefficients = compute_forces();

        // Track maximum drag and minimum lift (after initial transient)
        // The initial condition is zero, so we skip the first part of the simulation
        // where the solution is not yet physical
        if (time > 0.1)
        {
            drag_coefficient_max = std::max(coefficients[0], drag_coefficient_max);
            lift_coefficient_min = std::min(coefficients[1], lift_coefficient_min);
        }

        // Output solution every 10 time steps
        if (time_step % 10 == 0)
            output(time_step);
    }
}

// -----------------------------------------------------------------------
// COMPUTE FORCES METHOD
// -----------------------------------------------------------------------

/**
 * @brief Computes drag and lift forces on the obstacle boundary.
 *
 * @details
 * This method integrates the stress tensor over the obstacle boundary
 * (boundary_id == id_obstacle) to compute:
 * - **Drag force**: F_D = ∫ (ν ∂u_t/∂n * n_y - p * n_x) dS
 * - **Lift force**: F_L = -∫ (ν ∂u_t/∂n * n_x + p * n_y) dS
 *
 * where:
 * - u_t is the tangential velocity
 * - n is the outward normal vector
 * - p is the pressure
 *
 * The forces are then normalized to compute coefficients:
 * - c_D = 2 * F_D / (ρ * Ū² * D)
 * - c_L = 2 * F_L / (ρ * Ū² * D)
 *
 * For 3D, the cylinder spans the full channel height, and the integration
 * accounts for the z-direction as well.
 *
 * @tparam dim Spatial dimension (2 or 3)
 * @return std::vector<double> Containing {drag_coefficient, lift_coefficient}
 *
 * @note Uses FEFaceValues for boundary integration
 * @note The tangent vector is obtained by rotating the normal vector by 90 degrees
 * @note For 3D, assumes the cylinder is aligned with the z-axis (n_z = 0 on lateral surface)
 */
template <unsigned int dim>
std::vector<double> NavierStokes<dim>::compute_forces()
{
    pcout << "===============================================" << std::endl;
    pcout << "Computing forces: " << std::endl;

    // FEFaceValues for evaluating FE functions on boundary faces
    FEFaceValues<dim> fe_face_values(
        *fe, *quadrature_face,
        update_values | update_quadrature_points | update_gradients |
            update_normal_vectors | update_JxW_values);

    const unsigned int n_q_face = quadrature_face->size();

    FEValuesExtractors::Vector velocity(0);
    FEValuesExtractors::Scalar pressure(dim);

    // Vectors to store solution values at face quadrature points
    std::vector<double> current_pressure_values(n_q_face);
    std::vector<Tensor<2, dim>> current_velocity_gradients(n_q_face);

    // Local force values (before MPI sum)
    double local_lift = 0.0;
    double local_drag = 0.0;

    // Loop over all locally owned cells
    for (const auto &cell : dof_handler.active_cell_iterators())
    {
        if (!cell->is_locally_owned())
            continue;

        // Check if cell has any boundary faces
        if (cell->at_boundary())
        {
            // Loop over all faces of the cell
            for (unsigned int f = 0; f < cell->n_faces(); ++f)
            {
                // Check if this face is on the boundary and is the obstacle
                if (cell->face(f)->at_boundary() &&
                    (cell->face(f)->boundary_id() == id_obstacle))
                {
                    fe_face_values.reinit(cell, f);

                    // Get pressure and velocity gradient values at face quadrature points
                    fe_face_values[pressure].get_function_values(solution, current_pressure_values);
                    fe_face_values[velocity].get_function_gradients(solution, current_velocity_gradients);

                    // Loop over all quadrature points on this face
                    for (unsigned int q = 0; q < n_q_face; ++q)
                    {
                        // Get the outward normal vector (negated for inward normal)
                        Tensor<1, dim> n = -fe_face_values.normal_vector(q);
                        const double nx = n[0];
                        const double ny = n[1];

                        // Compute tangent vector by rotating normal by 90 degrees
                        // For 2D: t = (ny, -nx)
                        // For 3D: t = (ny, -nx, 0) - assumes cylinder aligned with z-axis
                        Tensor<1, dim> tangent;
                        tangent[0] = ny;
                        tangent[1] = -nx;
                        if constexpr (dim == 3)
                            tangent[2] = 0.;

                        // Compute shear stress: ν * (∇u · n) · t / |t|²
                        // This is the tangential component of the viscous stress
                        const double shear = (n * current_velocity_gradients[q]) *
                                             (tangent / tangent.norm_square());

                        // =========================================================
                        // DRAG FORCE: ∫ (ν * shear * n_y - p * n_x) dS
                        // =========================================================
                        local_drag += (density * kinematic_viscosity * shear * ny -
                                       current_pressure_values[q] * nx) *
                                      fe_face_values.JxW(q);

                        // =========================================================
                        // LIFT FORCE: -∫ (ν * shear * n_x + p * n_y) dS
                        // =========================================================
                        local_lift -= (density * kinematic_viscosity * shear * nx +
                                       current_pressure_values[q] * ny) *
                                      fe_face_values.JxW(q);
                    }
                }
            }
        }
    }

    // Sum local contributions across all MPI processes
    const double drag = Utilities::MPI::sum(local_drag, MPI_COMM_WORLD);
    const double lift = Utilities::MPI::sum(local_lift, MPI_COMM_WORLD);

    pcout << "Drag :\t " << drag << " Lift :\t " << lift << std::endl;

    // Compute mean velocity (used for coefficient normalization)
    // Note: getMeanVelocity() uses dimension-dependent coefficients:
    // - 2D: (2/3) * peak_velocity * temporal_envelope
    // - 3D: (4/9) * peak_velocity * temporal_envelope
    const double mean_velocity = inlet_velocity.getMeanVelocity();

    // Cylinder geometry parameters (from DFG benchmark)
    const double cylinder_diameter = 0.1;
    const double channel_height = 0.41;

    // Denominator for coefficient normalization
    double denom;
    if constexpr (dim == 3)
        // 3D: force is integrated over cylinder height, so divide by (D * H)
        denom = density * mean_velocity * mean_velocity *
                cylinder_diameter * channel_height;
    else
        // 2D: force is per unit depth, so divide by D only
        denom = density * mean_velocity * mean_velocity * cylinder_diameter;

    // Compute drag and lift coefficients
    const double drag_coefficient = (2. * drag) / denom;
    const double lift_coefficient = (2. * lift) / denom;
    std::vector<double> coefficients = {drag_coefficient, lift_coefficient};
    pcout << "Coeff:\t " << drag_coefficient << " Coeff:\t " << lift_coefficient << std::endl;

    vec_drag_force.push_back(drag);
    vec_lift_force.push_back(lift);
    vec_drag_coeff.push_back(drag_coefficient);
    vec_lift_coeff.push_back(lift_coefficient);

    pcout << "===============================================" << std::endl;
    return coefficients;
}

// -----------------------------------------------------------------------
// COMPUTE PRESSURE DIFFERENCE METHOD
// -----------------------------------------------------------------------

/**
 * @brief Computes the pressure difference between two points across the obstacle.
 *
 * @details
 * This method evaluates the pressure at two specific points:
 * - **Upstream**: (0.15, 0.2) for 2D or (0.45, 0.2, 0.205) for 3D
 * - **Downstream**: (0.25, 0.2) for 2D or (0.55, 0.2, 0.205) for 3D
 *
 * The pressure difference ΔP = P(upstream) - P(downstream) is computed and
 * printed by the root process.
 *
 * @tparam dim Spatial dimension (2 or 3)
 *
 * @note Uses VectorTools::point_value to evaluate pressure at specific points
 * @note Handles cases where the point is not on the current process's domain
 * @note Uses MPI_Reduce to collect pressure values from all processes
 * @note Points are chosen according to DFG benchmark specifications
 */
template <unsigned int dim>
void NavierStokes<dim>::compute_pressure_difference()
{
    // Define upstream and downstream points based on dimension
    Point<dim> pressure_point_upstream;
    Point<dim> pressure_point_downstream;

    if constexpr (dim == 3)
    {
        // 3D points (from DFG benchmark)
        pressure_point_upstream = Point<dim>(0.45, 0.2, 0.205);
        pressure_point_downstream = Point<dim>(0.55, 0.2, 0.205);
    }
    else
    {
        // 2D points (from DFG benchmark)
        pressure_point_upstream = Point<dim>(0.15, 0.2);
        pressure_point_downstream = Point<dim>(0.25, 0.2);
    }

    // Vectors to store solution values at the points
    Vector<double> solution_values_upstream(dim + 1);
    Vector<double> solution_values_downstream(dim + 1);

    bool upstream_point_available = true;
    bool downstream_point_available = true;

    // Try to evaluate solution at upstream point
    try
    {
        VectorTools::point_value(
            this->dof_handler, this->solution,
            pressure_point_upstream, solution_values_upstream);
    }
    catch (const dealii::VectorTools::ExcPointNotAvailableHere &)
    {
        upstream_point_available = false;
    }

    // Try to evaluate solution at downstream point
    try
    {
        VectorTools::point_value(
            this->dof_handler, this->solution,
            pressure_point_downstream, solution_values_downstream);
    }
    catch (const dealii::VectorTools::ExcPointNotAvailableHere &)
    {
        downstream_point_available = false;
    }

    // Initialize pressure variables
    double pressure_upstream = 0.0;
    double pressure_downstream = 0.0;

    // Assign pressure values if available (pressure is the last component)
    if (upstream_point_available)
        pressure_upstream = solution_values_upstream(dim);
    if (downstream_point_available)
        pressure_downstream = solution_values_downstream(dim);

    // Variables for global pressure values (on rank 0)
    double global_pressure_upstream = 0.0;
    double global_pressure_downstream = 0.0;

    // Sum pressure contributions from all processes
    MPI_Reduce(
        &pressure_upstream, &global_pressure_upstream,
        1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(
        &pressure_downstream, &global_pressure_downstream,
        1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // Compute and print pressure difference on root process
    if (this->mpi_rank == 0)
    {
        double pressure_difference = global_pressure_upstream - global_pressure_downstream;
        pcout << "Pressure difference (P(A) - P(B)) = "
              << pressure_difference << std::endl;
    }

    // Ensure all processes have completed the reductions
    MPI_Barrier(MPI_COMM_WORLD);
}

// -----------------------------------------------------------------------
// EXPLICIT TEMPLATE INSTANTIATION
// -----------------------------------------------------------------------

/**
 * @brief Explicit template instantiation for 2D and 3D.
 *
 * @details
 * This forces the compiler to generate code for both dim=2 and dim=3.
 * Without this, the compiler would only generate code for the dimensions
 * actually used in the code.
 */
template class NavierStokes<2>;
template class NavierStokes<3>;