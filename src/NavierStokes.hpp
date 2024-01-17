#ifndef NAVIER_STOKES_HPP
#define NAVIER_STOKES_HPP

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_in.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

using namespace dealii;

// Class implementing a solver for the Navier Stokes problem.
class NavierStokes
{
public:
  // Physical parameters //////////////////////////////////////////////////////
  // Physical dimension (currently only dim=3 is supported).
  static constexpr unsigned int dim = 3;
  // Kinematic viscosity [m^2/s].
  static constexpr double nu = 1e-3;
  // Fluid density [kg/m^3].
  static constexpr double ro = 1.0;
  // Cyclinder diameter [m].
  static constexpr double D = 0.1;
  // Inlet side length [m].
  static constexpr double H = 0.41;
  // Amplitude of inlet velocity.
  static constexpr double U_m = (dim == 2) ? 1.5 : 2.25;
  // Outlet pressure [Pa] (the outflow condition can be changed freely).
  static constexpr double p_out = 10.0;


  // Function for inlet velocity. This actually returns an object with four
  // components (one for each velocity component, and one for the pressure), but
  // then only the first three are really used (see the component mask when
  // applying boundary conditions at the end of assembly). If we only return
  // three components, however, we may get an error message due to this function
  // being incompatible with the finite element space.
  // This is the function for the steady case.
  class InletVelocity : public Function<dim>
  {
  public:
    InletVelocity()
      : Function<dim>(dim + 1)
    {}

    virtual void
    vector_value(const Point<dim> &p, Vector<double> &values) const override
    {
      if constexpr(dim == 2) {
        values[0] = 4.0 * U_m * p[1] * (H - p[1]) / (H * H);
      } else {
        values[0] = 16.0 * U_m * p[1] * p[2] * (H - p[1]) * (H - p[2]) / (H * H * H * H);
      }

      for (unsigned int i = 1; i < dim + 1; ++i)
        values[i] = 0.0;
    }

    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override
    {
      if (component == 0) {
        if constexpr(dim == 2) {
          return 4.0 * U_m * p[1] * (H - p[1]) / (H * H);
        } else {
          return 16.0 * U_m * p[1] * p[2] * (H - p[1]) * (H - p[2]) / (H * H * H * H);
        }
      }
      else {
        return 0.0;
      }
    }
  };

  // Function for the initial conditions (u=p=0).
  class InitialConditions : public Function<dim>
  {
  public:
    InitialConditions()
      : Function<dim>(dim + 1)
    {}

    virtual void
    vector_value(const Point<dim> &/*p*/, Vector<double> &values) const override
    {
      for (unsigned int i = 0; i < dim + 1; ++i)
        values[i] = 0.0;
    }

    virtual double
    value(const Point<dim> &/*p*/, const unsigned int /*component*/ = 0) const override
    {
      return 0.0;
    }
  };

  // A class to calculate the Reynolds number of the system.
  class ReynoldsNumber
  {
  public:
    double
    getValue() const
    {
      if constexpr(dim == 2) {
        return 2.0 * inlet_velocity.value(Point<dim>(0, H/2.0), 0) / 3.0 * D / nu;
      } else {
        return 4.0 * inlet_velocity.value(Point<dim>(0, H/2.0, H/2.0), 0) / 9.0 * D / nu;
      }
    }
  protected:
    InletVelocity inlet_velocity;
  };

  // Constructor.
  NavierStokes(const std::string  &mesh_file_name_,
         const unsigned int &degree_velocity_,
         const unsigned int &degree_pressure_,
         const double &T_,
         const double &deltat_)
    : mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
    , mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    , pcout(std::cout, mpi_rank == 0)
    , mesh_file_name(mesh_file_name_)
    , degree_velocity(degree_velocity_)
    , degree_pressure(degree_pressure_)
    , T(T_)
    , deltat(deltat_)
    , mesh(MPI_COMM_WORLD)
  {}

  // Initialization.
  void
  setup();

  // Solve the problem.
  void
  solve();

protected:
  // Assemble the constant part of the matrix.
  void
  assemble_constant();

  // Assemble the right-hand side and nonlinear term.
  void
  assemble_time_dependent();

  // Solve the problem for one time step.
  void
  solve_time_step();

  // Output.
  void
  output(const unsigned int &time_step) const;

  // MPI parallel. /////////////////////////////////////////////////////////////

  // Number of MPI processes.
  const unsigned int mpi_size;

  // This MPI process.
  const unsigned int mpi_rank;

  // Parallel output stream.
  ConditionalOStream pcout;

  // Problem. //////////////////////////////////////////////////////////////////

  // Inlet velocity.
  InletVelocity inlet_velocity;

  // Initial conditions.
  InitialConditions initial_conditions;

  // Reynolds number.
  ReynoldsNumber reynolds_number;

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