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

#include "SolverOptions.hpp"

using namespace dealii;

// This template class defines a solver for the Navier-Stokes problem. It is
// templated on the dimension, which can be either 2 or 3. The reason for this
// is that we need the dimension to be known at compile time.
template <unsigned int dim>
class NavierStokes {
 public:
  // Virtual destructor.
  virtual ~NavierStokes() = default;

  // The setup function is used to initialize the mesh, the finite element
  // space, the DoF handler, the DoF indices, and the matrices and vectors
  // for the linear system.
  void setup();

  // This calls the problem solution execution, which mainly consists of
  // calculating the initial solution and then performing the time stepping
  // until the final time is reached.
  void solve();

 protected:
  // Methods. //////////////////////////////////////////////////////////////////

  // Constructor.
  NavierStokes(const std::string &mesh_file_name_,
               const unsigned int &degree_velocity_,
               const unsigned int &degree_pressure_, const double &T_,
               const double &deltat_, const SolverOptions &solver_options_)
      : mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)),
        mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)),
        pcout(std::cout, mpi_rank == 0),
        mesh_file_name(mesh_file_name_),
        degree_velocity(degree_velocity_),
        degree_pressure(degree_pressure_),
        T(T_),
        deltat(deltat_),
        time_step(0),
        mesh(MPI_COMM_WORLD),
        solver_options(solver_options_) {}

  // The matrix assembly is split into two parts: the constant part, which is
  // assembled only once, and the time-dependent part, which is assembled at
  // each time step. This is done to avoid assembling the constant part at each
  // time step, which would be a waste of computational resources. The constant
  // part is assembled in this function.
  void assemble_constant();

  // The time dependent part of the matrix is assembled in this function, which
  // is called at each time step.
  void assemble_time_dependent();

  // Apply the initial conditions to the system.
  virtual void apply_initial_conditions();

  // This function is used to solve the linear system for one single time step.
  // In the `solve` function, this function is called in a loop until the final
  // time is reached.
  virtual void solve_time_step();

  // This is a simple function to output the solution to the /results folder.
  void output() const;

  // MPI parallel. /////////////////////////////////////////////////////////////

  // We use MPI to parallelize the code. This is the number of MPI processes.
  const unsigned int mpi_size;

  // This is the rank of the current MPI process.
  const unsigned int mpi_rank;

  // In order not to print the same output multiple times, we use a conditional
  // output stream, which prints only if the rank is 0.
  ConditionalOStream pcout;

  // Problem. //////////////////////////////////////////////////////////////////

  // Fluid density.
  double ro;

  // Kinematic viscosity.
  double nu;

  // Initial conditions.
  std::shared_ptr<Function<dim>> initial_conditions;

  // Boundary conditions.
  std::map<types::boundary_id, Function<dim> *> dirichlet_boundary_functions;
  std::map<types::boundary_id, Function<dim> *> neumann_boundary_functions;

  // Discretization. ///////////////////////////////////////////////////////////

  // Mesh file name.
  const std::string mesh_file_name;

  // Polynomial degree used for velocity.
  const unsigned int degree_velocity;

  // Polynomial degree used for pressure.
  const unsigned int degree_pressure;

  // This is the overall length of the simulation [s].
  const double T;

  // Duration of a single time step [s].
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

  // Velocity mass matrix divided by deltat (M/deltat). The symbols M and M_u
  // will be used as synonyms.
  TrilinosWrappers::BlockSparseMatrix velocity_mass;

  // Pressure mass matrix divided by nu (M_p/nu).
  TrilinosWrappers::BlockSparseMatrix pressure_mass;

  // Constant part of the system matrix (M/deltat + A  B^T; -B 0).
  TrilinosWrappers::BlockSparseMatrix constant_matrix;

  // System matrix [F  B^T; -B 0], where F = M/deltat + A + C.
  TrilinosWrappers::BlockSparseMatrix system_matrix;

  // Right-hand side vector in the linear system at the current time step.
  TrilinosWrappers::MPI::BlockVector system_rhs;

  // System solution (without ghost elements).
  TrilinosWrappers::MPI::BlockVector solution_owned;

  // System solution (including ghost elements).
  TrilinosWrappers::MPI::BlockVector solution;

  // Solver options.
  const SolverOptions solver_options;
};

#endif