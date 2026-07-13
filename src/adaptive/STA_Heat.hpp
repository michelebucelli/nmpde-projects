#ifndef STA_HEAT
#define STA_HEAT

#include <deal.II/base/utilities.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/conditional_ostream.h> //

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/affine_constraints.h> //

#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h> //
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_refinement.h> //
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/distributed/solution_transfer.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/error_estimator.h> //
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/matrix_creator.h>
#include <deal.II/numerics/matrix_tools.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "../timing/Timing.hpp"

using namespace dealii;

/**
 * Class managing the Heat problem.
 */ 
class AdaptiveHeat
{
public:
	// Physical dimension (1D, 2D, 3D)
	static constexpr unsigned int dim = 3;

	// Class Constructor
	AdaptiveHeat(const unsigned int 																						 &r_,
							 const double 		  																						 &T_,
							 const double 		  																						 &theta_,
							 const double 		  																						 &delta_t_,
							 const std::function<double(const Point<dim> &)> 								 &mu_, // in this project 1.0
							 const std::function<double(const Point<dim> &, const double &)> &f_)
		: r(r_)
		, T(T_)
		, theta(theta_)
		, delta_t(delta_t_)
		, mu(mu_)
		, f(f_)
		, mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
		, mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
		, mesh(MPI_COMM_WORLD)
		, pcout(std::cout, mpi_rank == 0)
	{}

	// Run the simulation for the Heat problem
	void run();

	// Computer L2 error against a baseline solution
	double l2_against_base(const Function<dim> &);

	void print_results();

protected:
	// Initialize the mesh and the FE space
	void setup();

	// Initialize the DoFs handler and the system
	void setup_system();

	// Assemble the system
	void assemble();

	// Solve the time step linear system 
	void solve_time_step();

	// Apply adaptive mesh refinement
	bool refine_grid(const unsigned int, const unsigned int);

	// Output
	void output() const;

	// Polynomial degree
	const unsigned int r;

	// Final time
	const double T;

	// Theta parameter for theta method
	const double theta;

	// Timestep
	double delta_t;

	// Current time
	double time = 0.0;

	// Current timestep number
	unsigned int timestep_number = 0;

	// For spatial refinement
	double spatial_tol = 5.0e-4; // start point, we can refine 
    unsigned int min_steps_between_refine = 5; // to not refine too often
    unsigned int last_refine_step = 0; // when was the last refinement

	// Diffusion coefficient
	std::function<double(const Point<dim> &)> mu;

	// Forcing term
	std::function<double(const Point<dim> &, const double &)> f;

	// Number of MPI processes
	const unsigned int mpi_size;

	// Rank of the current MPI process
	const unsigned int mpi_rank;

	// Number of steps
	unsigned int num_of_steps = 0;
	unsigned int max_dofs = 0;

	// Collect data about execution
	Profiler profiler;

	// Note: deal.II tutorial uses parallel::distributed::Triangulation.
	// fullydistributed is newer and more scalable (does not replicate coarse mesh
	// in all mpi processes) however it is more difficult to implement Adaptive Mesh Refinement
	// since parallel::distributed::Triangulation already supports many useful features.
	// Triangulation
	parallel::distributed::Triangulation<dim> mesh;

	// Finite element space
	std::unique_ptr<FE_Q<dim>> fe;

	// Quadrature formula
	std::unique_ptr<QGauss<dim>> quadrature;

	// DoF handler
	DoFHandler<dim> dof_handler;

	// Holds a list of constraints to hold the hanging nodes and the boundary conditions.
	AffineConstraints<double> constraints;

	// System matrix
	TrilinosWrappers::SparseMatrix system_matrix;

	// System right-hand side
	TrilinosWrappers::MPI::Vector system_rhs;

	// System solution, no ghost elements
	TrilinosWrappers::MPI::Vector solution_owned;

	// System solution, with ghost elements
	TrilinosWrappers::MPI::Vector solution;

	// Output stream, only for process 0
	ConditionalOStream pcout;
};

#endif