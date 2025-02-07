#ifndef STEADYNAVIERSTOKES_HPP
#define STEADYNAVIERSTOKES_HPP

#include "includes_file.hpp"

using namespace dealii;

// ==================================================================
// Base Class: SteadyNavierStokes
//
// Description:
//   This class solves the steady incompressible Navier-Stokes equations
//   using a monolithic approach. The class is templated on the
//   dimensionality of the problem in order to handle 2D and 3D
//   problems.
//
//  =================================================================

template <int dim>
class SteadyNavierStokes
{
public:
	// ---------------------------------------------------------------
	// Class: PreconditionIdentity
	//
	// Description:
	//   This class defines a trivial identity preconditioner.
	// ---------------------------------------------------------------
	class PreconditionIdentity
	{
	public:
		void vmult(TrilinosWrappers::MPI::BlockVector &dst,
				   const TrilinosWrappers::MPI::BlockVector &src) const
		{
			dst = src;
		}
	};

	// ---------------------------------------------------------------
	// Class: PreconditionBlockTriangularStokes
	//
	// Description:
	//   This class defines a block triangular preconditioner for the
	//   Stokes system. The preconditioner is based on the block structure:
	//
	//       [ A  B^T ]
	//       [ B  0   ]
	//
	//   where A is the velocity stiffness matrix, B is the pressure
	//   mass matrix, and 0 is the zero matrix.
	// ---------------------------------------------------------------
	class PreconditionBlockTriangularStokes
	{
	public:
		void initialize(const TrilinosWrappers::SparseMatrix &velocity_stiffness_,
						const TrilinosWrappers::SparseMatrix &pressure_mass_,
						const TrilinosWrappers::SparseMatrix &B_)
		{
			velocity_stiffness = &velocity_stiffness_;
			pressure_mass = &pressure_mass_;
			B = &B_;

			preconditioner_velocity.initialize(velocity_stiffness_);
			preconditioner_pressure.initialize(pressure_mass_);
		}

		void vmult(TrilinosWrappers::MPI::BlockVector &dst,
				   const TrilinosWrappers::MPI::BlockVector &src) const
		{
			// 1) Solve velocity block
			{
				SolverControl solver_control_velocity(10000,
													  1e-2 * src.block(0).l2_norm());
				SolverCG<TrilinosWrappers::MPI::Vector> solver_cg_velocity(solver_control_velocity);
				solver_cg_velocity.solve(*velocity_stiffness,
										 dst.block(0),
										 src.block(0),
										 preconditioner_velocity);
			}
			// 2) Solve pressure block
			{
				tmp.reinit(src.block(1));
				B->vmult(tmp, dst.block(0));
				tmp.sadd(-1.0, src.block(1));

				SolverControl solver_control_pressure(10000,
													  1e-2 * src.block(1).l2_norm());
				SolverCG<TrilinosWrappers::MPI::Vector> solver_cg_pressure(solver_control_pressure);
				solver_cg_pressure.solve(*pressure_mass,
										 dst.block(1),
										 tmp,
										 preconditioner_pressure);
			}
		}

	protected:
		const TrilinosWrappers::SparseMatrix *velocity_stiffness = nullptr;
		const TrilinosWrappers::SparseMatrix *pressure_mass = nullptr;
		const TrilinosWrappers::SparseMatrix *B = nullptr;

		TrilinosWrappers::PreconditionILU preconditioner_velocity;
		TrilinosWrappers::PreconditionILU preconditioner_pressure;

		mutable TrilinosWrappers::MPI::Vector tmp;
	};

	// ---------------------------------------------------------------
	// Class: InletVelocity
	//
	// Description:
	//   This class defines an inlet velocity function.
	//   It sets the maximum velocity value (uM) based on the
	//   dimensionality of the problem: 0.3 for 2D problems, 0.45 for 3D problems.
	//
	// Parameters:
	//   H - characteristic length of the domain.
	//   uM - maximum velocity value.
	// ---------------------------------------------------------------
	class InletVelocity : public Function<dim>
	{
	public:
		InletVelocity(const double H_in, const double uMax_in)
			: Function<dim>(dim + 1), H(H_in), uM(uMax_in)
		{
		}

		void vector_value(const Point<dim> &p, Vector<double> &values) const override
		{
			if constexpr (dim == 2)
			{
				values[0] = 4.0 * uM * p[1] * (H - p[1]) / (H * H);
				for (unsigned int i = 1; i < dim + 1; ++i)
					values[i] = 0.0;
			}
			else if constexpr (dim == 3)
			{
				values[0] = 16.0 * uM * p[1] * p[2] * (H - p[1]) * (H - p[2]) / (H * H * H * H);
				for (unsigned int i = 1; i < dim + 1; ++i)
					values[i] = 0.0;
			}
		}

		double value(const Point<dim> &p, const unsigned int comp = 0) const override
		{
			if (comp == 0)
			{
				if constexpr (dim == 2)
					return 4.0 * uM * p[1] * (H - p[1]) / (H * H);
				else if constexpr (dim == 3)
					return 16.0 * uM * p[1] * p[2] * (H - p[1]) * (H - p[2]) / (H * H * H * H);
			}
			return 0.0;
		}

	private:
		const double H;
		const double uM;
	};



	// ============================== PUBLIC FUNCTIONS ===============================

	// ............................................................
	// Constructor
	// ............................................................
	// Parameters:
	//   mesh_file_name_in - name of the mesh file.
	//   degree_velocity_in - polynomial degree for velocity.
	//   degree_pressure_in - polynomial degree for pressure.
	//   Re_in - Reynolds number.
	// ............................................................

	SteadyNavierStokes(
		const std::string &mesh_file_name_in,
		const unsigned int degree_velocity_in,
		const unsigned int degree_pressure_in,
		const double Re_in)
		: mesh_file_name(mesh_file_name_in), degree_velocity(degree_velocity_in), degree_pressure(degree_pressure_in), Re(Re_in), H(0.41), D(0.1), uMax(dim == 2 ? 0.3 : 0.45), uMean(dim == 2 ? 2. / 3. * uMax : 4. / 9. * uMax), nu(uMean * D / Re), p_out(0.0), forcing_term(), inlet_velocity(H, uMax), mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)), mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)), pcout(std::cout, mpi_rank == 0), mesh(MPI_COMM_WORLD), dof_handler(mesh)
		{}

	// Disallow copy construction for clarity.
	SteadyNavierStokes(const SteadyNavierStokes<dim> &) = delete;
	virtual ~SteadyNavierStokes() = default;

	auto run_full_problem_pipeline() -> void; // Instantiates and runs the full problem pipeline and computes lift/drag.

	auto get_mesh_file_name() const -> const std::string& // returns the mesh file name
	{return mesh_file_name;}

	auto get_degree_velocity() const -> unsigned int // returns the velocity polynomial degree
	{return degree_velocity;}

	auto get_degree_pressure() const -> unsigned int // returns the pressure polynomial degree
	{return degree_pressure;}

	auto get_uMax() const -> double // returns the maximum inflow velocity
	{return uMax;}

	auto get_H() const -> double // returns the height of the channel
	{return H;}

	auto get_p_out() const -> double // returns the outlet Neumann BC value
	{return p_out;}

	auto get_Re() const -> double // returns the Reynolds number
	{return Re;}

	// Provide access to the mesh
	const parallel::fullydistributed::Triangulation<dim> &get_mesh() const
	{
		return mesh;
	}



	//=============================== PROTECTED FUNCTIONS ===============================
protected:
	virtual auto setup() -> void; // Setup the problem by initializing the mesh, DoF handler, and finite element spaces.

	virtual auto assemble() -> void; // Assemble the system matrix and right-hand side.

	virtual auto solve() -> void; // Solve the linear system.	

	virtual auto output() -> void; // Save the output of the computation in a pvtk format.
	
	virtual auto get_output_directory() const -> std::string; // Defines the path of the directory where the outputs will be stored



	// ================================ PRIVATE VARIABLES ===============================

	// ================================
	// Problem Parameters

	std::string mesh_file_name;  							// Mesh file name
	const unsigned int degree_velocity;  					// Velocity polynomial degree
	const unsigned int degree_pressure;  					// Pressure polynomial degree
	const double Re;  										// Reynolds number

	// ================================
	// Geometrical & Physical Values

	const double H;      									// Height of the channel
	const double D;      									// Diameter of the obstacle
	const double uMax;   									// Maximum inflow velocity
	const double uMean;  									// Mean inflow velocity
	const double nu;     									// Viscosity
	const double p_out;  									// Outlet Neumann BC value
	const bool   zero_neumann = true;  							// Homogeneous Neumann BC flag

	// ================================
	// Post-processing Data

	double lift; 											// Lift coefficient
	double drag; 											// Drag coefficient
	double deltaP; 											// Pressure drop

	// ================================
	// Problem-specific Objects

	Functions::ZeroFunction<dim> forcing_term;  			// Zero forcing term - fixed by problem
	InletVelocity inlet_velocity;  							// Inlet velocity

	// ================================
	// MPI Tools

	const unsigned int mpi_size;  							// Number of MPI processes
	const unsigned int mpi_rank;  							// Rank of the process
	ConditionalOStream pcout;  								// Conditional output stream

	// ================================
	// FEM Objects

	parallel::fullydistributed::Triangulation<dim> mesh;  	// Parallel mesh
	DoFHandler<dim> dof_handler;  							// DoF handler
	std::unique_ptr<FiniteElement<dim>> fe;  				// Finite element handler
	std::unique_ptr<Quadrature<dim>> quadrature;  			// Domain quadrature
	std::unique_ptr<Quadrature<dim - 1>> quadrature_face;  	// Face quadrature
	IndexSet locally_owned_dofs;  							// Locally owned DoFs
	std::vector<IndexSet> block_owned_dofs;  				// Block-wise owned DoFs
	IndexSet locally_relevant_dofs;  						// Locally relevant DoFs
	std::vector<IndexSet> block_relevant_dofs;  			// Block-wise relevant DoFs

	// ================================
	// Matrices & Vectors

	TrilinosWrappers::BlockSparseMatrix system_matrix;  	// System LHS matrix
	TrilinosWrappers::BlockSparseMatrix pressure_mass;  	// Pressure mass matrix
	TrilinosWrappers::MPI::BlockVector system_rhs;  		// System RHS vector
	TrilinosWrappers::MPI::BlockVector solution_owned;  	// Solution vector (owned)
	TrilinosWrappers::MPI::BlockVector solution;  			// Solution vector (complete)

};





// ==================================================================
// Derived class: Stokes
//
// Description:
//   This class solves the stokes problem under the same 
// 	 physical conditions without taking into account the
//   convective term.
// ==================================================================
template <int dim>
class Stokes : public SteadyNavierStokes<dim>
{
public:

	// ============================== PUBLIC FUNCTIONS ===============================

	// ............................................................
	// Constructor
	// ............................................................
	// Parameters:
	//   mesh_file_name_in - name of the mesh file.
	//   degree_velocity_in - polynomial degree for velocity.
	//   degree_pressure_in - polynomial degree for pressure.
	//   Re_in - Reynolds number.
	// ............................................................
	
	Stokes(const std::string &mesh_file_name_in,
		   unsigned int degree_velocity_in,
		   unsigned int degree_pressure_in,
		   double Re_in)
		: SteadyNavierStokes<dim>(mesh_file_name_in,
								  degree_velocity_in,
								  degree_pressure_in,
								  Re_in)
	{
	}

	auto get_solution() const -> TrilinosWrappers::MPI::BlockVector
	{return this->solution;};

	auto setup() -> void override; // Setup the problem by initializing the mesh, DoF handler, and finite element spaces.

	auto assemble() -> void override; // Assemble the system matrix and right-hand side.

	auto solve() -> void override; // Solve the linear system.

	auto output() -> void override; // Save the output of the computation in a pvtk format.

	auto get_output_directory() const -> std::string override; // Defines the path of the directory where the outputs will be stored
	
};

// ==================================================================
// Derived class: NonLinearCorrection
//
// Description:
//   This class solves the non-linear correction problem by
//   taking into account the convective term.
// ==================================================================

template <int dim>
class NonLinearCorrection : public SteadyNavierStokes<dim>
{
public:
	// ============================== PUBLIC FUNCTIONS ===============================

	// ............................................................
	// Constructor
	// ............................................................
	// Parameters:
	//   stokes_obj - Stokes object to copy the mesh and parameters.
	// ............................................................
	
	NonLinearCorrection(const Stokes<dim> &stokes_obj)
		: SteadyNavierStokes<dim>(stokes_obj.get_mesh_file_name(),
								  stokes_obj.get_degree_velocity(),
								  stokes_obj.get_degree_pressure(),
								  stokes_obj.get_Re()),
		  u_k(0), p_k(dim)
	{
		// Copy the fully-distributed Triangulation
		this->mesh.copy_triangulation(stokes_obj.get_mesh());

		// Compute scaling_factor for lift and drag computation, based on dimension
		/*
		 *  Scaling factor for flow past cylinder test case
		 *  2D: scaling_factor = 2 / (Umean^2 * D)
		 *  3D: scaling_factor = 2 / (Umean^2 * D * H)
		 */
		double uMean = this->uMean;
		double D = this->D;
		if constexpr (dim == 2)
		{
			scaling_factor = 2.0 / (uMean * uMean * D);
		}
		else if constexpr (dim == 3)
		{
			double H = this->H;
			scaling_factor = 2.0 / (uMean * uMean * D * H);
		}
		else
		{
			static_assert(dim == 2 || dim == 3, "Dimension must be 2 or 3.");
		}
	}

	auto setup() -> void override; // Setup the problem by initializing the mesh, DoF handler, and finite element spaces.

	auto assemble() -> void override; // Assemble the system matrix and right-hand side.

	auto solve() -> void override; // Solve the linear system.

	auto output() -> void override; // Save the output of the computation in a pvtk format.

	auto get_output_directory() const -> std::string override; // Defines the path of the directory where the outputs will be stored

	// Set the initial conditions for the iterative scheme
	auto set_initial_conditions(const TrilinosWrappers::MPI::BlockVector solution_stokes_) -> void
	{
		solution_old.reinit(solution_stokes_);
		solution_old = solution_stokes_;
	}

	auto compute_lift_drag() -> void; // Compute lift and drag coefficients

protected:
	// ================================ PROTECTED FUNCTIONS ===============================

	// ================================
	// Newton Iteration Parameters
	unsigned int iter = 0;  							// Newton iterations counter
	static constexpr unsigned int maxIter = 20;  		// Maximum iterations
	static constexpr double tolerance = 1e-7;  			// Update tolerance

	// ================================
	// Extractors and Constraints
	FEValuesExtractors::Vector u_k;  					// Velocity extractor
	FEValuesExtractors::Scalar p_k;  					// Pressure extractor
	AffineConstraints<double> constraints;  			// Constraints

	// ================================
	// Iterative Scheme Data
	TrilinosWrappers::MPI::BlockVector solution_old;  	// Old solution
	TrilinosWrappers::MPI::BlockVector new_res;  		// Residual at current iteration

	// ================================
	// Post-Processing Data
	double scaling_factor;  				// Scaling factor for lift and drag
	std::vector<Point<dim>> points_of_interest;  // Points of interest for pressure difference

};

#endif // STEADYNAVIERSTOKES_HPP