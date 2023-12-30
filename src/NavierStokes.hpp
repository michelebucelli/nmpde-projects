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
  // Physical dimension. (can be 2 or 3)
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
  static constexpr double p_out = 10;

  // Function for the forcing term.
  // While it is always 0, it was implemented in case changes have to be made.
  class ForcingTerm : public Function<dim>
  {
  public:
    virtual void
    vector_value(const Point<dim> & /*p*/,
                 Vector<double> &values) const override
    {
      for (unsigned int i = 0; i < dim; ++i)
        values[i] = 0.0;
    }

    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int component = 0) const override
    {
      return 0.0;
    }
  };

  // Function for inlet velocity. This actually returns an object with four
  // components (one for each velocity component, and one for the pressure), but
  // then only the first three are really used (see the component mask when
  // applying boundary conditions at the end of assembly). If we only return
  // three components, however, we may get an error message due to this function
  // being incompatible with the finite element space.
  class InletVelocity : public Function<dim>
  {
  public:
    InletVelocity()
      : Function<dim>(dim + 1)
    {}

    // For the unsteady case, multiply the first component by sin(pi*t/8.0) and add time dependency (1 more dimension to p).
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
  
  // Function for the initial conditions. Refer to explanation of InletVelocity for details.
  class InitialConditions : public Function<dim>
  {
  public:
    InitialConditions()
      : Function<dim>(dim + 1)
    {}

    virtual void
    vector_value(const Point<dim> &p, Vector<double> &values) const override
    {
      for (unsigned int i = 0; i < dim + 1; ++i)
        values[i] = 0.0;
    }

    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override
    {
      return 0.0;
    }
  };

  // Since we're working with block matrices, we need to make our own
  // preconditioner class. A preconditioner class can be any class that exposes
  // a vmult method that applies the inverse of the preconditioner.

  // Identity preconditioner.
  class PreconditionIdentity
  {
  public:
    // Application of the preconditioner: we just copy the input vector (src)
    // into the output vector (dst).
    void
    vmult(TrilinosWrappers::MPI::BlockVector       &dst,
          const TrilinosWrappers::MPI::BlockVector &src) const
    {
      dst = src;
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
  };

  // Constructor.
  NavierStokes(const std::string  &mesh_file_name_,
         const unsigned int &degree_velocity_,
         const unsigned int &degree_pressure_)
    : mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
    , mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    , pcout(std::cout, mpi_rank == 0)
    , mesh_file_name(mesh_file_name_)
    , degree_velocity(degree_velocity_)
    , degree_pressure(degree_pressure_)
    , mesh(MPI_COMM_WORLD)
  {}

  // Setup system.
  void
  setup();

  // Assemble system. We also assemble the pressure mass matrix (needed for the
  // preconditioner).
  void
  assemble();

  // Solve system.
  void
  solve();

  // Output results.
  void
  output();

protected:
  // MPI parallel. /////////////////////////////////////////////////////////////

  // Number of MPI processes.
  const unsigned int mpi_size;

  // This MPI process.
  const unsigned int mpi_rank;

  // Parallel output stream.
  ConditionalOStream pcout;

  // Problem definition. ///////////////////////////////////////////////////////

  // Forcing term.
  static ForcingTerm forcing_term;

  // Inlet velocity.
  static InletVelocity inlet_velocity;
  
  // Initial conditions.
  static InitialConditions initial_conditions;
  
  // Reynolds number.
  static ReynoldsNumber reynolds_number;

  // Discretization. ///////////////////////////////////////////////////////////

  // Mesh file name.
  const std::string mesh_file_name;

  // Polynomial degree used for velocity.
  const unsigned int degree_velocity;

  // Polynomial degree used for pressure.
  const unsigned int degree_pressure;

  // Mesh.
  parallel::fullydistributed::Triangulation<dim> mesh;

  // Finite element space.
  std::unique_ptr<FiniteElement<dim>> fe;

  // Quadrature formula.
  std::unique_ptr<Quadrature<dim>> quadrature;

  // Quadrature formula for face integrals.
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

  // System matrix.
  TrilinosWrappers::BlockSparseMatrix system_matrix;

  // Pressure mass matrix, needed for preconditioning. We use a block matrix for
  // convenience, but in practice we only look at the pressure-pressure block.
  TrilinosWrappers::BlockSparseMatrix pressure_mass;

  // Right-hand side vector in the linear system.
  TrilinosWrappers::MPI::BlockVector system_rhs;

  // System solution (without ghost elements).
  TrilinosWrappers::MPI::BlockVector solution_owned;

  // System solution (including ghost elements).
  TrilinosWrappers::MPI::BlockVector solution;
};

#endif
