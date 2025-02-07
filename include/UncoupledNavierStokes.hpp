#ifndef UNCOUPLED_NAVIER_STOKES_HPP
#define UNCOUPLED_NAVIER_STOKES_HPP

#include "includes_file.hpp"

using namespace dealii;

// ==================================================================
// Class: UncoupledNavierStokes
//
// Description:
//   This class solves the incompressible Navier-Stokes equations
//   using an uncoupled approach. The class is templated on the
//   dimensionality of the problem in order to handle 2D and 3D
//   problems.
//
//  =================================================================

template <unsigned int dim>
class UncoupledNavierStokes
{

public:
    // ---------------------------------------------------------------
    // Class: InletVelocity
    //
    // Description:
    //   This class defines an inlet velocity function.
    //   It sets the maximum velocity value (uM) based on the
    //   dimensionality of the problem: 1.5 for 2D problems, 2.25 for 3D problems.
    //
    // Parameters:
    //   H - characteristic length of the domain.
    //   uM - maximum velocity value.
    // ---------------------------------------------------------------
    class InletVelocity : public Function<dim>
    {
    public:
        InletVelocity(const double H)
            : Function<dim>(dim)
        {
            this->H = H;
            if constexpr (dim == 2)
                this->uM = 1.5;
            else
                this->uM = 2.25;
        }

        virtual void
        vector_value(const Point<dim> &p, Vector<double> &values) const override
        {
            if constexpr (dim == 2)
                values[0] = 4.0 * uM * p[1] * (H - p[1]) / (H * H);
            else
                values[0] = 16.0 * uM * p[1] * (H - p[1]) * p[2] * (H - p[2]) / (H * H * H * H);
            for (unsigned int i = 1; i < dim; ++i)
                values[i] = 0.0;
        }

        virtual double
        value(const Point<dim> &p, const unsigned int component = 0) const override
        {
            if (component == 0)
                if constexpr (dim == 2)
                    return 4.0 * uM * p[1] * (H - p[1]) / (H * H);
                else
                    return 16.0 * uM * p[1] * (H - p[1]) * p[2] * (H - p[2]) / (H * H * H * H);
            else
                return 0.0;
        }

        double get_u_max() const
        {
            return uM;
        }

    protected:
        double uM;
        double H;
    };

    // ============================== PUBLIC FUNCTIONS ===============================
    // ............................................................
    // Constructor
    // ............................................................
    // Parameters:
    //   mesh_file_name_ - name of the mesh file.
    //   degree_velocity_ - velocity polynomial degree.
    //   degree_pressure_ - pressure polynomial degree.
    //   T_ - final time.
    //   deltat_ - time step size.
    //   reynolds_number_ - Reynolds number.
    // ............................................................

UncoupledNavierStokes(
        const std::string &mesh_file_name_,
        const unsigned int &degree_velocity_,
        const unsigned int &degree_pressure_,
        const double &T_,
        const double &deltat_,
        const double &reynolds_number_)
        : reynolds_number(reynolds_number_),
        T(T_),
        deltat(deltat_),
        mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)),
        mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)),
        pcout(std::cout, mpi_rank == 0),
        mesh(MPI_COMM_WORLD),
        mesh_file_name(mesh_file_name_),
        triangulation(),
        fe_velocity(FE_SimplexP<dim>(2), dim),
        dof_handler_velocity(triangulation),
        fe_pressure(1),
        dof_handler_pressure(triangulation),
        degree_velocity(degree_velocity_),
        degree_pressure(degree_pressure_),
        inlet_velocity(H),
        computing_timer(MPI_COMM_WORLD, pcout,
                        TimerOutput::summary,
                        TimerOutput::wall_times)
    {
        this->nu = (2. / 3.) * inlet_velocity.get_u_max() * cylinder_radius / reynolds_number;
    }

    auto run() -> void;

    // ============================== PRIVATE FUNCTIONS ==============================
private:

    auto setup() -> void; // Setup the problem by initializing the mesh, DoF handler, and finite element spaces.

    auto assemble_system_velocity() -> void; // Assemble the system matrix and right-hand side for the velocity problem.

    auto solve_velocity_system() -> void; // Solve the velocity system.

    auto assemble_system_pressure() -> void; // Assemble the system matrix and right-hand side for the pressure problem.

    auto solve_pressure_system() -> void; // Solve the pressure system.

    auto update_velocity() -> void; // Assemble the system matrix and right-hand side for the velocity update problem.

    auto solve_update_velocity_system() -> void; // Solve the velocity update system.

    auto pressure_update(bool rotational) -> void; // Update the pressure field. If therotational flag is true, the rotational term is included.

    auto output_results() -> void; // Save the output of the computation in a pvtk format.

    auto compute_lift_drag() -> void; // Compute lift and drag coefficients

    auto get_output_directory() -> std::string; // Defines the path of the directory where the outputs will be stored

    // ================================ PRIVATE VARIABLES ===============================

    // ================================
    // Geometrical and Physical Parameters

    const double H = 0.41;                                      // Height of the channel
    const double cylinder_radius = 0.1;                         // Cylinder radius of the obstacle
    const double rho = 1.0;                                     // Fluid density
    const double reynolds_number;                               // Reynolds number, governing flow characteristics
    double nu;                                                  // Kinematic viscosity of the fluid

    // ================================
    // Time and Simulation Control

    const double T;                                             // Final simulation time
    double deltat;                                              // Time step size
    double time;                                                // Current simulation time
    unsigned int time_step = 0;                                 // Time step counter
    unsigned int timestep_number;                               // Number of time steps performed

    // ================================
    // MPI and Parallelization

    const unsigned int mpi_size;                                // Number of MPI processes
    const unsigned int mpi_rank;                                // Rank of the current MPI process
    ConditionalOStream pcout;                                   // Parallel output stream for controlled logging
    parallel::fullydistributed::Triangulation<dim> mesh;        // Fully distributed parallel mesh

    // ================================
    // Mesh and Finite Element Setup

    const std::string mesh_file_name;                           // Name of the mesh file
    Triangulation<dim> triangulation;                           // Main triangulation object for the domain

    FESystem<dim> fe_velocity;                                  // Velocity finite element: Uses a Q2 (quadratic) vector-valued basis
    DoFHandler<dim> dof_handler_velocity;                       // DoF handler for velocity field

    FE_SimplexP<dim> fe_pressure;                               // Pressure finite element: Uses a Q1 (linear) scalar-valued basis
    DoFHandler<dim> dof_handler_pressure;                       // DoF handler for pressure field

    const unsigned int degree_velocity;                         // Polynomial degree for velocity field
    const unsigned int degree_pressure;                         // Polynomial degree for pressure field

    // ================================
    // Boundary and Initial Conditions

    Functions::ZeroFunction<dim> forcing_term;                  // External forcing term (zero in this case)
    Functions::ZeroFunction<dim> neumann_function;              // Neumann boundary condition function
    Functions::ZeroFunction<dim> initial_condition;             // Initial velocity and pressure condition
    InletVelocity inlet_velocity;                               // Prescribed velocity profile at the inlet

    // Note that forcing_term and neumann_function are not applied in this case 
    // To see an example of implementation please check the Numerical-Test Branch

    // ================================
    // Timing and Computational Monitoring

    TimerOutput computing_timer;                                // Timer for performance monitoring

    // ================================
    // Constraints and Degrees of Freedom

    AffineConstraints<double> constraints_velocity;             // Affine constraints for velocity field
    AffineConstraints<double> constraints_pressure;             // Affine constraints for pressure field

    IndexSet locally_owned_velocity;                            // Velocity DoFs owned by the current process
    IndexSet locally_relevant_velocity;                         // Velocity DoFs relevant for the current process
    IndexSet locally_owned_pressure;                            // Pressure DoFs owned by the current process
    IndexSet locally_relevant_pressure;                         // Pressure DoFs relevant for the current process

    // ================================
    // System Matrices

    TrilinosWrappers::SparseMatrix velocity_matrix;             // System matrix for velocity field
    TrilinosWrappers::SparseMatrix pressure_matrix;             // System matrix for pressure field
    TrilinosWrappers::SparseMatrix velocity_update_matrix;      // Matrix used for velocity updates

    // ================================
    // System Vectors

    TrilinosWrappers::MPI::Vector old_velocity;                 // Velocity field at previous time step
    TrilinosWrappers::MPI::Vector old_old_velocity;             // Velocity field two time steps ago
    TrilinosWrappers::MPI::Vector u_star;                       // Intermediate velocity field in fractional step method
    TrilinosWrappers::MPI::Vector u_star_divergence;            // Divergence of u_star field
    TrilinosWrappers::MPI::Vector velocity_solution;            // Solution vector for velocity field
    TrilinosWrappers::MPI::Vector update_velocity_solution;     // Update for velocity field
    TrilinosWrappers::MPI::Vector velocity_system_rhs;          // Right-hand side of the velocity system
    TrilinosWrappers::MPI::Vector velocity_update_rhs;          // Right-hand side of the velocity update system

    TrilinosWrappers::MPI::Vector old_pressure;                 // Pressure field at previous time step
    TrilinosWrappers::MPI::Vector deltap;                       // Change in pressure between iterations
    TrilinosWrappers::MPI::Vector pressure_solution;            // Solution vector for pressure field
    TrilinosWrappers::MPI::Vector pressure_system_rhs;          // Right-hand side of the pressure system

    // ================================
    // Post-Processing Data

    bool rotational = false;                                    // Flag to indicate whether rotation effects are included

    std::vector<double> vec_drag;                               // History of drag force values over time
    std::vector<double> vec_lift;                               // History of lift force values over time
    std::vector<double> vec_drag_coeff;                         // History of drag coefficient values
    std::vector<double> vec_lift_coeff;                         // History of lift coefficient values

    double lift;                                                // Current lift force value
    double drag;                                                // Current drag force value

};

#endif // UNCOUPLED_NAVIER_STOKES_HPP

