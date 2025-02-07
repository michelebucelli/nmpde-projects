#ifndef MONOLITHICNAVIERSTOKES_HPP
#define MONOLITHICNAVIERSTOKES_HPP

#include "includes_file.hpp"
using namespace dealii;

// ==================================================================
// Class: MonolithicNavierStokes
//
// Description:
//   This class solves the incompressible Navier-Stokes equations
//   using a monolithic approach. The class is templated on the
//   dimensionality of the problem in order to handle 2D and 3D
//   problems.
//
//  =================================================================

template <unsigned int dim>
class MonolithicNavierStokes
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
        InletVelocity(const double H_)
            : Function<dim>(dim), H(H_)
        {
            static_assert(dim == 2 || dim == 3,
                          "Dimensions other than 2 or 3 are not supported");

            if constexpr (dim == 2)
                this->uM = 1.5;
            else if constexpr (dim == 3)
                this->uM = 2.25;
        }

        virtual void
        vector_value(const Point<dim> &p, Vector<double> &values) const override
        {
            static_assert(dim == 2 || dim == 3,
                          "Dimensions other than 2 or 3 are not supported");
            if constexpr (dim == 2)
                values[0] = 4.0 * uM * p[1] * (H - p[1]) / (H * H);
            else if constexpr (dim == 3)
                values[0] = 16.0 * uM * p[1] * (H - p[1]) * p[2] * (H - p[2]) / (H * H * H * H);
        }

        virtual double
        value(const Point<dim> &p, const unsigned int component = 0) const override
        {
            if (component == 0)
            {
                if constexpr (dim == 2)
                    return 4.0 * uM * p[1] * (H - p[1]) / (H * H);
                else if constexpr (dim == 3)
                    return 16.0 * uM * p[1] * (H - p[1]) * p[2] * (H - p[2]) / (H * H * H * H);
            }
            else
                return 0.0;
        }

        double get_u_max() const
        {
            return uM;
        }

    protected:
        double uM;
        const double H;
    };

    // ============================== PUBLIC FUNCTIONS ===============================

    // ............................................................
	// Constructor 
	// ............................................................
	// Parameters:
	//   mesh_file_name_ - name of the mesh file.
    //   degree_velocity_ - polynomial degree for velocity.
    //   degree_pressure_ - polynomial degree for pressure.
    //   T_ - final time.
    //   deltat_ - time step.
    //   re_ - Reynolds number.
    // ............................................................

    MonolithicNavierStokes(
        const std::string &mesh_file_name_,
        const unsigned int &degree_velocity_,
        const unsigned int &degree_pressure_,
        const double &T_,
        const double &deltat_,
        const double &re_)
        : mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
        , mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
        , pcout(std::cout, mpi_rank == 0)
        , mesh_file_name(mesh_file_name_)
        , mesh(MPI_COMM_WORLD)
        , reynolds_number(re_)
        , T(T_)
        , deltat(deltat_)
        , degree_velocity(degree_velocity_)
        , degree_pressure(degree_pressure_)
        , inlet_velocity(H)
        , velocity(0) 
        , pressure(dim)
    {
        this->nu = (2. / 3.) * inlet_velocity.get_u_max() * cylinder_radius / reynolds_number;
    }

    auto run() -> void; // Function to run the full problem pipeline.
    

    // ============================== PRIVATE FUNCTIONS ==============================
private:
    auto setup() -> void; // Setup the problem.

    auto assemble_base_matrix() -> void; // Assemble the time independent part of the system matrix.

    auto add_convective_term() -> void; // Assemble the convective part of the system matrix.

    auto assemble_rhs() -> void; // Assemble the right-hand side of the problem.

    auto solve_time_step() -> void; // Group all the instructions that need to be executed at each time step.

    auto solve() -> void; // Solve the entire problem by looping over time steps.

    auto output(const unsigned int &time_step) -> void; // Save the output of the computation in a pvtk format.

    auto get_output_directory() const -> std::string; // Defines the path of the directory where the outputs will be stored


    // ================================ PRIVATE VARIABLES ===============================

    // ================================
    // MPI and Parallelization
    
    const unsigned int mpi_size;                            // Number of MPI processes.

    const unsigned int mpi_rank;                            // This MPI process.

    ConditionalOStream pcout;                               // Parallel output stream.

    // ================================
    // Mesh and Geometry Parameters

    const std::string mesh_file_name;                       // Mesh file name.

    static constexpr double H = 0.41;                       // Height of the channel.

    static constexpr double cylinder_radius = 0.1;          // Cylinder radius.

    parallel::fullydistributed::Triangulation<dim> mesh;    // Mesh.

    // ================================
    // Physical Parameters

    const double reynolds_number;                           // Reynolds number.

    double nu;                                              // Diffusion coefficient.

    // ================================
    // Time and Simulation Control

    const double T;                                         // Final time.

    const double deltat;                                    // Time step.

    double time;                                            // Current time.

    // ================================
    // Finite Element and Discretization

    const unsigned int degree_velocity;                     // Polynomial degree for velocity.

    const unsigned int degree_pressure;                     // Polynomial degree for pressure.

    std::unique_ptr<FiniteElement<dim>> fe;                 // Finite element space.

    std::unique_ptr<Quadrature<dim>> quadrature;            // Quadrature formula.

    std::unique_ptr<Quadrature<dim - 1>> quadrature_face;   // Quadrature formula for boundary.

    // ================================
    // Degrees of Freedom (DoFs) and Constraints

    DoFHandler<dim> dof_handler;                            // DoF handler.

    IndexSet locally_owned_dofs;                            // DoFs owned by the current process.

    std::vector<IndexSet> block_owned_dofs;                 // Block-wise partitioning of owned DoFs.

    IndexSet locally_relevant_dofs;                         // DoFs relevant to the current process (including ghosts).

    std::vector<IndexSet> block_relevant_dofs;              // Block-wise partitioning of relevant DoFs.

    AffineConstraints<double> constraints;                  // Affine constraints.

    // ================================
    // Boundary and Initial Conditions

    InletVelocity inlet_velocity;                           // Inlet velocity.

    Functions::ZeroFunction<dim> forcing_term;              // Forcing term.

    const bool zero_forcing = 0;                            // Indicator for zero forcing term.

    Functions::ZeroFunction<dim> neumann_function;          // Neumann boundary condition.

    const bool zero_neumann = 0;                            // Indicator for zero Neumann condition.

    Functions::ZeroFunction<dim> initial_condition;         // Initial condition.

    // ================================
    // Extractors for Field Variables

    FEValuesExtractors::Vector velocity;                    // Extractor for velocity.

    FEValuesExtractors::Scalar pressure;                    // Extractor for pressure.

    // ================================
    // System Matrices and Vectors

    TrilinosWrappers::BlockSparseMatrix system_matrix;      // Time-independent matrix.

    TrilinosWrappers::BlockSparseMatrix velocity_mass;      // Velocity mass matrix.

    TrilinosWrappers::BlockSparseMatrix pressure_mass;      // Pressure mass matrix.

    TrilinosWrappers::BlockSparseMatrix lhs_matrix;         // Complete system matrix.

    TrilinosWrappers::MPI::BlockVector system_rhs;          // Right-hand side vector.

    TrilinosWrappers::MPI::BlockVector solution_owned;      // System solution without ghosts.

    TrilinosWrappers::MPI::BlockVector solution;            // System solution with ghosts.
};

#endif