/**
 * @file NavierStokes.hpp
 * @brief Header file for a Navier-Stokes solver implemented using the deal.II library.
 *
 * This file defines a C++ class template for solving the incompressible Navier-Stokes equations
 * in 2D or 3D using the Finite Element Method (FEM). It supports both steady and unsteady inflow
 * conditions and includes functionality for computing drag and lift forces, as well as
 * pressure differences. The implementation leverages MPI for parallel computation.
 */

#ifndef NAVIER_STOKES_HPP
#define NAVIER_STOKES_HPP

#include <fstream>
#include <filesystem>
#include <iostream>
#include <functional>
#include <limits>
#include <algorithm>
#include <stdexcept>

/**
 * @brief Include for custom preconditioner definitions used in the solver.
 */
#include "./Preconditioners.hpp"
#include "./IncludeFiles.hpp"

/**
 * @brief Namespace alias for deal.II library components.
 */
using namespace dealii;

/**
 * @enum InflowRegime
 * @brief Enumeration defining the type of inflow boundary condition for the fluid flow.
 *
 * This enum class specifies whether the inflow velocity profile is steady (constant in time)
 * or unsteady (time-dependent, e.g., sinusoidal).
 */
enum class InflowRegime
{
    Steady,
    Unsteady
};

/**
 * @class NavierStokes
 * @brief A template class for solving the incompressible Navier-Stokes equations.
 *
 * @tparam dim The spatial dimension (2 or 3). The class is templated to support both 2D and 3D problems.
 *
 * This class implements a finite element solver for the incompressible Navier-Stokes equations.
 * It supports parallel computation via MPI and includes methods for assembling the system,
 * solving the time-stepping scheme, and computing aerodynamic forces (drag and lift).
 * The solver uses a block matrix structure to handle the coupled velocity-pressure system.
 *
 * @note The class enforces a static assertion to ensure 'dim' is either 2 or 3.
 */
template <unsigned int dim>
class NavierStokes
{
    /**
     * @brief Static assertion to ensure the solver is only instantiated for 2D or 3D.
     *
     * This prevents accidental instantiation of the class for unsupported dimensions.
     */
    static_assert(dim == 2 || dim == 3, "Navier Stokes has been implemented only for 2D and 3D case");

public:
    // -----------------------------------------------------------------------
    // Boundary ID constants
    // -----------------------------------------------------------------------

    /**
     * @brief Boundary ID for the inlet boundary.
     *
     * Used to identify the inlet boundary in the mesh for applying inlet velocity conditions.
     */
    static constexpr types::boundary_id id_inlet = 1;

    /**
     * @brief Boundary ID for the outlet boundary.
     *
     * Used to identify the outlet boundary in the mesh for applying outlet pressure conditions.
     */
    static constexpr types::boundary_id id_outlet = 2;

    /**
     * @brief Boundary ID for the wall boundaries.
     *
     * Used to identify the wall boundaries in the mesh for applying no-slip conditions.
     */
    static constexpr types::boundary_id id_walls = 3;

    /**
     * @brief Boundary ID for the obstacle boundary.
     *
     * Used to identify the obstacle boundary in the mesh for applying no-slip conditions
     * or for computing forces (drag and lift) on the obstacle.
     */
    static constexpr types::boundary_id id_obstacle = 4;

    // -----------------------------------------------------------------------
    // Nested classes
    // -----------------------------------------------------------------------

    /**
     * @class InletVelocity
     * @brief A nested class defining the inlet velocity profile for the Navier-Stokes problem.
     *
     * This class inherits from 'Function<dim>' (a deal.II base class for analytical functions)
     * and provides a custom velocity profile at the inlet boundary. The profile is parabolic
     * in 2D and 3D, and can be either steady or unsteady (time-dependent).
     */
    class InletVelocity : public Function<dim>
    {
    public:
        /**
         * @brief Constructor for the InletVelocity class.
         *
         * @param regime_ The inflow regime (steady or unsteady). Default is 'InflowRegime::Steady'.
         * @param peak_velocity_ The peak velocity of the flow at the inlet. Default is 1.5.
         *
         * Initializes the inlet velocity profile with the specified regime and peak velocity.
         * The constructor of the base class 'Function<dim>' is called with 'dim + 1' to account
         * for both velocity components (dim) and pressure (1).
         */
        InletVelocity(const InflowRegime regime_ = InflowRegime::Steady, const double peak_velocity_ = 1.5)
            : Function<dim>(dim + 1), regime(regime_), peak_velocity(peak_velocity_)
        {
        }

        /**
         * @brief Computes the vector-valued function at a given point.
         *
         * @param p The point in physical space where the function is evaluated.
         * @param values The output vector where the function values are stored.
         *               'values[0]' corresponds to the x-velocity component,
         *               'values[1]' to the y-velocity component (and z in 3D),
         *               and 'values[dim]' to the pressure (always zero here).
         *
         * This method sets the x-velocity component using the 'inflow_profile' method
         * and sets all other components (y, z, and pressure) to zero.
         */
        virtual void
        vector_value(const Point<dim> &p, Vector<double> &values) const override
        {
            values[0] = inflow_profile(p);
            for (unsigned int i = 1; i < dim + 1; ++i)
            {
                values[i] = 0.0;
            }
        }

        /**
         * @brief Computes the scalar value of a specific component of the function at a given point.
         *
         * @param p The point in physical space where the function is evaluated.
         * @param component The component of the function to evaluate (0 for x-velocity, 1 for y-velocity, etc.).
         * @return The value of the specified component at point 'p'.
         *
         * For 'component == 0', it returns the inflow profile value (x-velocity).
         * For all other components, it returns 0.0.
         */
        virtual double
        value(const Point<dim> &p, const unsigned int component = 0) const override
        {
            if (component == 0)
            {
                return inflow_profile(p);
            }
            else
            {
                return 0.0;
            }
        }

        /**
         * @brief Computes the mean velocity of the inlet flow.
         *
         * @return The mean velocity across the inlet boundary.
         *
         * The mean velocity is computed using a dimension-dependent coefficient:
         * - For 2D: '2.0 / 3.0 * peak_velocity * temporal_envelope(time)'
         * - For 3D: '4.0 / 9.0 * peak_velocity * temporal_envelope(time)'
         */
        double getMeanVelocity() const
        {
            const double mean_coefficient = (dim == 2) ? (2.0 / 3.0) : (4.0 / 9.0);
            return mean_coefficient * peak_velocity * temporal_envelope(this->get_time());
        }

    protected:
        /**
         * @brief Computes the temporal envelope for the inlet velocity.
         *
         * @param time The current simulation time.
         * @return The temporal envelope value at the given time.
         *
         * The temporal envelope modifies the velocity profile over time:
         * - For 'InflowRegime::Steady': returns 1.0 (constant in time).
         * - For 'InflowRegime::Unsteady': returns 'sin(π * time / 8.0)' (sinusoidal variation).
         */
        double temporal_envelope(double time) const
        {
            switch (regime)
            {
            case InflowRegime::Steady:
                return 1.;
            case InflowRegime::Unsteady:
                return std::sin(M_PI * time / 8.0);
            default:
                throw std::runtime_error("Unknown InflowRegime.");
            }
        }

        /**
         * @brief Computes the spatial profile of the inlet velocity.
         *
         * @param p The point in physical space where the profile is evaluated.
         * @return The velocity value at point 'p' based on the spatial profile.
         *
         * The spatial profile is parabolic:
         * - In 2D: '4 * peak_velocity * p[1] * (channel_height - p[1]) / (channel_height^2)'
         * - In 3D: '16 * peak_velocity * p[1] * p[2] * (channel_height - p[1]) * (channel_height - p[2]) / (channel_height^4)'
         *
         * The profile is scaled by the temporal envelope to account for time dependence.
         */
        double inflow_profile(const Point<dim> &p) const
        {
            const double time_component = temporal_envelope(this->get_time());

            if constexpr (dim == 2)
            {
                return 4. * peak_velocity * p[1] * (channel_height - p[1]) * time_component / (channel_height * channel_height);
            }
            else
            {
                return 16. * peak_velocity * p[1] * p[2] * (channel_height - p[1]) * (channel_height - p[2]) * time_component /
                       (channel_height * channel_height * channel_height * channel_height);
            }
        }
        /**
         * @brief Inflow regime (steady or unsteady).
         */
        const InflowRegime regime;

        /**
         * @brief Peak velocity of the inflow profile.
         */
        const double peak_velocity;

        /**
         * @brief Height of the channel for the inlet velocity profile.
         *
         * This is a constant value used to define the parabolic profile.
         */
        const double channel_height = 0.41;
    };

    // -----------------------------------------------------------------------
    // Constructor and public methods
    // -----------------------------------------------------------------------

    /**
     * @brief Constructor for the NavierStokes class.
     *
     * @param mesh_file_name_ The path to the mesh file to be loaded.
     * @param degree_velocity_ The polynomial degree for the velocity finite element space.
     * @param degree_pressure_ The polynomial degree for the pressure finite element space.
     * @param final_time_ The final simulation time.
     * @param time_step_size_ The size of each time step.
     * @param peak_velocity_ The peak velocity for the inlet boundary condition.
     * @param regime_ The inflow regime (steady or unsteady). Default is 'InflowRegime::Steady'.
     * @param preconditioner_ The preconditioner to use for solving the linear system. Default is 'Preconditioner::YOSIDA'.
     *
     * Initializes the solver with the provided parameters and sets up MPI-related variables.
     * The mesh is initialized using 'MPI_COMM_WORLD' for parallel computation.
     */
    NavierStokes(const std::string &mesh_file_name_,
                 const unsigned int &degree_velocity_,
                 const unsigned int &degree_pressure_,
                 const double &final_time_,
                 const double &time_step_size_,
                 const double &peak_velocity_,
                 const InflowRegime regime_ = InflowRegime::Steady,
                 const Preconditioner preconditioner_ = Preconditioner::YOSIDA,
                 const std::string &output_dir_ = "output")
        : mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)),
          mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)),
          pcout(std::cout, mpi_rank == 0),
          inlet_velocity(regime_, peak_velocity_),
          final_time(final_time_),
          mesh_file_name(mesh_file_name_),
          degree_velocity(degree_velocity_),
          degree_pressure(degree_pressure_),
          preconditioner(preconditioner_),
          time_step_size(time_step_size_),
          mesh(MPI_COMM_WORLD),
          output_dir(output_dir_)

    {
    }

    /**
     * @brief Sets up the finite element system, including the mesh, DoF handlers, and matrices.
     *
     * This method initializes the finite element space, quadrature rules, and distributed
     * triangulation. It also sets up the block structure for the velocity and pressure fields.
     */
    void
    setup();

    /**
     * @brief Runs the simulation from the initial time to the final time.
     *
     * This method performs the time-stepping loop, assembling the system at each time step,
     * solving the linear system, and outputting results. It also computes and stores drag and lift forces.
     */
    void
    run();

    /**
     * @brief Vector storing the drag force values computed at each time step.
     *
     * Drag force is typically computed on the obstacle boundary ('id_obstacle').
     */
    std::vector<double> vec_drag_force;

    /**
     * @brief Vector storing the lift force values computed at each time step.
     *
     * Lift force is typically computed on the obstacle boundary ('id_obstacle').
     */
    std::vector<double> vec_lift_force;

    /**
     * @brief Vector storing the drag coefficient values computed at each time step.
     *
     * The drag coefficient is a dimensionless quantity representing the drag force
     * normalized by dynamic pressure and reference area.
     */
    std::vector<double> vec_drag_coeff;

    /**
     * @brief Vector storing the lift coefficient values computed at each time step.
     *
     * The lift coefficient is a dimensionless quantity representing the lift force
     * normalized by dynamic pressure ad reference area.
     */
    std::vector<double> vec_lift_coeff;

    /**
     * @brief Vector storing the time taken for preconditioning at each time step.
     *
     * This is used for performance profiling of the preconditioner.
     */
    std::vector<double> time_preconditioning;

    /**
     * @brief Vector storing the time taken for solving the linear system at each time step.
     *
     * This is used for performance profiling of the solver.
     */
    std::vector<double> time_solve;

    // -----------------------------------------------------------------------
    // Stability / accuracy diagnostics
    // -----------------------------------------------------------------------

    /**
     * @brief Minimum cell diameter of the mesh, computed once in setup().
     *
     * Used together with time_step_size and the velocity field to compute the
     * advective CFL number at each time step.
     */
    double mesh_h_min = 0.0;

    /**
     * @brief Vector storing the advective CFL number at each time step,
     * CFL = ||u||_inf * time_step_size / mesh_h_min.
     *
     * This is a diagnostic only: the scheme is implicit in the diffusion and
     * time-derivative terms, so it does not need CFL <= 1 to remain stable in
     * the linear-algebra sense, but a large CFL still degrades the temporal
     * accuracy of the semi-implicit (Picard-lagged) convection term.
     */
    std::vector<double> vec_cfl;

    /**
     * @brief Pressure difference P(upstream) - P(downstream), computed by
     * compute_pressure_difference() near the final time step, per the DFG
     * benchmark point definitions.
     */
    double pressure_difference = 0.0;

protected:
    // -----------------------------------------------------------------------
    // Protected methods
    // -----------------------------------------------------------------------

    /**
     * @brief Assembles the system matrices and right-hand side vector at a given time.
     *
     * @param time The current simulation time.
     *
     * This method computes the contributions to the system matrix and right-hand side vector
     * from the Navier-Stokes equations (convection, diffusion, and pressure terms).
     */
    void assemble(const double &time);

    /**
     * @brief Assembles the time-dependent terms for the current time step.
     *
     * @param time The current simulation time.
     *
     * This method updates the system matrices and right-hand side vector to account for
     * the time discretization (backward Euler).
     */
    void assemble_time_step(const double &time);

    /**
     * @brief Solves the linear system for the current time step.
     *
     * @param preconditioner The preconditioner to use for solving the linear system.
     *
     * This method solves the block linear system for the velocity and pressure fields
     * using the specified preconditioner.
     */
    void solve_time_step(const Preconditioner &preconditioner);

    /**
     * @brief Outputs the solution (velocity and pressure) at a given time step.
     *
     * @param time The current time step index.
     *
     * This method writes the solution fields to disk in VTU format for visualization.
     */
    void output(const unsigned int &time);

    /**
     * @brief Computes the drag and lift forces acting on the obstacle.
     *
     * @return A vector containing the drag and lift forces (and possibly other force-related quantities).
     *
     * This method integrates the stress tensor over the obstacle boundary to compute
     * the aerodynamic forces.
     */
    std::vector<double> compute_forces();

    /**
     * @brief Computes the pressure difference between the inlet and outlet.
     *
     * This method calculates the average pressure at the inlet and outlet boundaries
     * and stores the difference, which is useful for analyzing pressure drops in the flow.
     */
    void compute_pressure_difference();

    // -----------------------------------------------------------------------
    // MPI parallel member variables
    // -----------------------------------------------------------------------

    /**
     * @brief Total number of MPI processes in the parallel computation.
     */
    const unsigned int mpi_size;

    /**
     * @brief Rank (ID) of the current MPI process.
     */
    const unsigned int mpi_rank;

    /**
     * @brief Parallel output stream for MPI.
     *
     * This stream ensures that only the root process (rank 0) outputs to 'std::cout',
     * avoiding duplicated output in parallel runs.
     */
    ConditionalOStream pcout;

    // -----------------------------------------------------------------------
    // Problem definition member variables
    // -----------------------------------------------------------------------

    /**
     * @brief Kinematic viscosity of the fluid [m²/s].
     *
     * This is a constant property of the fluid, set to 1e-3.
     */
    const double kinematic_viscosity = 1e-3;

    /**
     * @brief Density of the fluid [kg/m³].
     *
     * This is a constant property of the fluid, set to 1.0 (dimensionless or normalized).
     */
    const double density = 1.;

    /**
     * @brief Inlet velocity function.
     *
     * This object defines the velocity profile at the inlet boundary.
     */
    InletVelocity inlet_velocity;

    /**
     * @brief Forcing term for the Navier-Stokes equations.
     *
     * This is a zero function by default, meaning no external body forces are applied.
     */
    Functions::ZeroFunction<dim> forcing_term;

    /**
     * @brief Final simulation time.
     *
     * The solver will run until this time is reached.
     */
    const double final_time;

    // -----------------------------------------------------------------------
    // Discretization member variables
    // -----------------------------------------------------------------------

    /**
     * @brief Path to the mesh file.
     *
     * The mesh is loaded from this file at the beginning of the simulation.
     */
    const std::string mesh_file_name;

    /**
     * @brief Polynomial degree for the velocity finite element space.
     */
    const unsigned int degree_velocity;

    /**
     * @brief Polynomial degree for the pressure finite element space.
     */
    const unsigned int degree_pressure;

    /**
     * @brief Preconditioner type to use for solving the linear system.
     *
     * This is set during construction and can be, e.g., 'Preconditioner::YOSIDA'.
     */
    Preconditioner preconditioner;

    /**
     * @brief Time step size for the simulation.
     *
     * This is a constant time step used for the time discretization.
     */
    const double time_step_size;

    /**
     * @brief Initial condition for the solution (velocity and pressure).
     *
     * By default, this is a zero function, meaning the flow starts from rest.
     */
    Functions::ZeroFunction<dim> initial_solution_function;

    /**
     * @brief Dirichlet lifting function (g(x)) for enforcing Dirichlet boundary conditions.
     *
     * This is a zero function by default, meaning no lifting is applied.
     * It can be used to enforce non-homogeneous Dirichlet conditions.
     */
    Functions::ZeroFunction<dim> dirichlet_lifting_function;

    // -----------------------------------------------------------------------
    // Mesh and finite element member variables
    // -----------------------------------------------------------------------

    /**
     * @brief Parallel distributed triangulation (mesh).
     *
     * This object represents the computational mesh, distributed across MPI processes.
     */
    parallel::fullydistributed::Triangulation<dim> mesh;

    const std::string output_dir;

    /**
     * @brief Finite element space.
     *
     * This object defines the finite element type (Taylor-Hood elements for velocity and pressure).
     */
    std::unique_ptr<FiniteElement<dim>> fe;

    /**
     * @brief Quadrature formula for volume integrals.
     *
     * This defines the integration points and weights for computing integrals over cells.
     */
    std::unique_ptr<Quadrature<dim>> quadrature;

    /**
     * @brief Quadrature formula for face integrals.
     *
     * This defines the integration points and weights for computing integrals over faces (boundaries).
     */
    std::unique_ptr<Quadrature<dim - 1>> quadrature_face;

    // -----------------------------------------------------------------------
    // Degree of Freedom (DoF) member variables
    // -----------------------------------------------------------------------

    /**
     * @brief Degree of Freedom (DoF) handler for the finite element system.
     */
    DoFHandler<dim> dof_handler;

    /**
     * @brief IndexSet of DoFs owned by the current MPI process.
     *
     * These are the degrees of freedom that the current process is responsible for.
     */
    IndexSet locally_owned_dofs;

    /**
     * @brief Vector of IndexSets for DoFs owned by the current process in each block (velocity and pressure).
     *
     * This is used for block-wise operations on the velocity and pressure fields.
     */
    std::vector<IndexSet> block_owned_dofs;

    /**
     * @brief IndexSet of DoFs relevant to the current MPI process (including ghost DoFs).
     *
     * These include both owned DoFs and ghost DoFs (DoFs owned by other processes but needed for local computations).
     */
    IndexSet locally_relevant_dofs;

    /**
     * @brief Vector of IndexSets for DoFs relevant to the current process in each block (velocity and pressure).
     *
     * This is used for block-wise operations, including ghost DoFs.
     */
    std::vector<IndexSet> block_relevant_dofs;

    // -----------------------------------------------------------------------
    // System matrices and vectors
    // -----------------------------------------------------------------------

    /**
     * @brief Block sparse matrix representing the system matrix for the Navier-Stokes equations.
     *
     * This matrix is structured in blocks to handle the coupled velocity-pressure system:
     * ```
     * [ A  B ]
     * [ B^T C ]
     * ```
     * where `A` is the velocity-velocity block, `B` is the velocity-pressure block,
     * and `C` is the pressure-pressure block (often zero for incompressible flow).
     */
    TrilinosWrappers::BlockSparseMatrix system_matrix;

    /**
     * @brief Block sparse matrix representing the mass matrix.
     *
     * The mass matrix is used for time discretization.
     */
    TrilinosWrappers::BlockSparseMatrix mass_matrix;

    /**
     * @brief Block sparse matrix representing the stiffness matrix (diffusion term).
     *
     * This matrix corresponds to the viscous (diffusion) term in the Navier-Stokes equations.
     */
    TrilinosWrappers::BlockSparseMatrix stiffness_matrix;

    /**
     * @brief Block sparse matrix representing the convection matrix.
     *
     * This matrix corresponds to the nonlinear convection term '(u · ∇) u' in the Navier-Stokes equations.
     */
    TrilinosWrappers::BlockSparseMatrix convection_matrix;

    /**
     * @brief Block sparse matrix representing the pressure mass matrix.
     *
     * This matrix is used for preconditioning and corresponds to the pressure-pressure block.
     * Although stored as a block matrix, only the pressure-pressure block is used in practice.
     */
    TrilinosWrappers::BlockSparseMatrix pressure_mass;

    // -----------------------------------------------------------------------
    // Right-hand side and solution vectors
    // -----------------------------------------------------------------------

    /**
     * @brief Block vector representing the right-hand side of the linear system.
     *
     * This vector stores the residual or source terms for the Navier-Stokes equations.
     */
    TrilinosWrappers::MPI::BlockVector system_rhs;

    /**
     * @brief Block vector storing the solution (velocity and pressure) without ghost elements.
     *
     * This vector contains the solution for the degrees of freedom owned by the current process.
     */
    TrilinosWrappers::MPI::BlockVector solution_owned;

    /**
     * @brief Block vector storing the solution (velocity and pressure) including ghost elements.
     *
     * This vector contains the solution for all degrees of freedom relevant to the current process,
     * including ghost DoFs for parallel computations.
     */
    TrilinosWrappers::MPI::BlockVector solution;

    /**
     * @brief Block vector storing the solution from the previous time step.
     *
     * This is used for time discretization schemes that require the solution from the previous step.
     */
    TrilinosWrappers::MPI::BlockVector previous_solution;
};

#endif