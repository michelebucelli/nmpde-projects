#ifndef NAVIER_STOKES_HPP
#define NAVIER_STOKES_HPP

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>

using namespace dealii;

// Class implementing a solver for the Navier Stokes problem.
template <unsigned int dim>
class NavierStokes {
 public:
  // Constructor.
  NavierStokes(const std::string &mesh_file_name_,
               const unsigned int &degree_velocity_,
               const unsigned int &degree_pressure_, const double &T_,
               const double &deltat_)
      : mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)),
        mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)),
        pcout(std::cout, mpi_rank == 0),
        mesh_file_name(mesh_file_name_),
        degree_velocity(degree_velocity_),
        degree_pressure(degree_pressure_),
        T(T_),
        deltat(deltat_),
        time_step(0),
        mesh(MPI_COMM_WORLD) {}

  // Destructor.
  virtual ~NavierStokes() = default;

  // Initialization.
  void setup();

  // Solve the problem.
  void solve();

 protected:
  // Assemble the constant part of the matrix.
  void assemble_constant();

  // Assemble the right-hand side and nonlinear term.
  void assemble_time_dependent();

  // Apply the initial condition.
  virtual void apply_initial_conditions();

  // Solve the problem for one time step.
  virtual void solve_time_step();

  // Output.
  void output() const;

  // MPI parallel. /////////////////////////////////////////////////////////////

  // Number of MPI processes.
  const unsigned int mpi_size;

  // This MPI process.
  const unsigned int mpi_rank;

  // Parallel output stream.
  ConditionalOStream pcout;

  // Problem. //////////////////////////////////////////////////////////////////

  // Fluid density [kg/m^3].
  double ro;

  // Kinematic viscosity [m^2/s].
  double nu;

  // Initial conditions.
  std::shared_ptr<Function<dim>> initial_conditions;

  // Boundary conditions.
  std::map<types::boundary_id, const Function<dim> *>
      dirichlet_boundary_functions;
  std::map<types::boundary_id, double> neumann_boundary_functions;

  // Discretization. ///////////////////////////////////////////////////////////

  // Mesh file name.
  const std::string mesh_file_name;

  // Polynomial degree used for velocity.
  const unsigned int degree_velocity;

  // Polynomial degree used for pressure.
  const unsigned int degree_pressure;

  // Final time.
  const double T;

  // Time step.
  const double deltat;

  // Current time step.
  unsigned int time_step;

  // Mesh.
  parallel::fullydistributed::Triangulation<dim> mesh;

  // Finite element space.
  std::unique_ptr<FiniteElement<dim>> fe;

  // Quadrature formula.
  std::unique_ptr<Quadrature<dim>> quadrature;

  // Quadrature formula for Neumann BC.
  std::unique_ptr<Quadrature<dim - 1>> quadrature_face;

  // DoF handler.
  DoFHandler<dim> dof_handler;

  // DoFs owned by current process.
  IndexSet locally_owned_dofs;

  // DoFs owned by current process in the velocity and pressure blocks.
  std::vector<IndexSet> block_owned_dofs;

  // DoFs relevant to the current process (including ghost DoFs).
  IndexSet locally_relevant_dofs;

  // DoFs relevant to current process in the velocity and pressure blocks.
  std::vector<IndexSet> block_relevant_dofs;

  // Velocity mass matrix (M/deltat).
  TrilinosWrappers::BlockSparseMatrix velocity_mass;

  // Constant part of the matrix (M/deltat + A B^T; -B 0).
  TrilinosWrappers::BlockSparseMatrix constant_matrix;

  // System matrix (constant_matrix + [C 0; 0 0])
  TrilinosWrappers::BlockSparseMatrix system_matrix;

  // Right-hand side vector in the linear system at the current time step.
  TrilinosWrappers::MPI::BlockVector system_rhs;

  // System solution (without ghost elements).
  TrilinosWrappers::MPI::BlockVector solution_owned;

  // System solution (including ghost elements).
  TrilinosWrappers::MPI::BlockVector solution;
};

#endif