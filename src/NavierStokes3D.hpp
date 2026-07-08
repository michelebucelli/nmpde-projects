#ifndef NAVIER_STOKES_3D_HPP
#define NAVIER_STOKES_3D_HPP

#include "preconditioners.hpp"

#include <limits>
#include <utility>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <fstream>

#include <deal.II/base/timer.h>

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
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

using namespace dealii;

// Class implementing a solver for the Stokes problem.
class NavierStokes3D
{
public:
  // Physical dimension (1D, 2D, 3D)
  static constexpr unsigned int dim = 3;

    struct ProblemData
  {
    double nu  = 0.001;
    double rho = 1.0;

    double Um = 0.0;

    unsigned int choice = 1;

    // Forcing term
    std::function<Tensor<1, dim>(const Point<dim> &, const double &)> f;

    // Neumann BC on outlet
    std::function<Tensor<1, dim>(const Point<dim> &, const double &)> h;

    // Inlet profile
    std::function<Tensor<1, dim>(const Point<dim> &, const double &)> u_in;

    // Initial condition
    std::function<Tensor<1, dim>(const Point<dim> &)> u0;
  };


  // Function for inlet velocity.
   class InletVelocity : public Function<dim>
  {
  public:
    InletVelocity(const std::function<Tensor<1, dim>(const Point<dim> &, const double &)> &fun)
      : Function<dim>(dim + 1)
      , fun(fun)
    {}

    void set_time(const double t) override
    { 
      time = t; 
    }

    void vector_value(const Point<dim> &p, Vector<double> &values) const override
    {
      const Tensor<1, dim> u = fun(p, time);

      values = 0.0;
      for (unsigned int d = 0; d < dim; ++d)
        values[d] = u[d];
    }

  private:
    std::function<Tensor<1, dim>(const Point<dim> &, const double &)> fun;
    double time = 0.0;
  };


  // Constructor.
  NavierStokes3D(const std::string                               &mesh_file_name_,
                 const unsigned int                              &degree_velocity_,
                 const unsigned int                              &degree_pressure_,
                 const double                                    &T_,
                 const double                                    &theta_,
                 const double                                    &delta_t_,
                 const ProblemData                               &problem_data_)
    : mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
    , mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    , pcout(std::cout, mpi_rank == 0)
    , mesh_file_name(mesh_file_name_)
    , degree_velocity(degree_velocity_)
    , degree_pressure(degree_pressure_)
    , T(T_)
    , theta(theta_)
    , delta_t(delta_t_)
    , data(problem_data_)
    , mesh(MPI_COMM_WORLD)
  {}

  // Run the time-dependent simulation.
  void
  run();

  // Setup system.
  void
  setup();

  // Assemble system.
  void
  assemble();

  // Solve system.
  void
  solve();

  // Output results.
  void
  output();

  // Benchmark quantities.
  double compute_reynolds_number() const;
  std::pair<double, double> compute_drag_lift_forces() const;
  double compute_pressure_difference() const;
  double compute_recirculation_length() const;
  double compute_strouhal_number() const;
  double compute_delta_p_at_half_period() const;

  struct BenchmarkResult
  {
    double Re     = 0.0;
    double cD     = 0.0;
    double cL     = 0.0;
    double dP     = 0.0;
    double La     = std::numeric_limits<double>::quiet_NaN();
    double cD_max = 0.0;
    double cL_max = 0.0;
    double St     = std::numeric_limits<double>::quiet_NaN();
    double dP_half = std::numeric_limits<double>::quiet_NaN();
  };

  BenchmarkResult compute_benchmark_result() const;
  void print_benchmark_quantities(const BenchmarkResult &res) const;

protected:

  // MPI parallel. /////////////////////////////////////////////////////////////

  // Number of MPI processes.
  const unsigned int mpi_size;

  // This MPI process.
  const unsigned int mpi_rank;

  // Parallel output stream.
  ConditionalOStream pcout;

  double cD_max = 0.0;
  double cL_max = 0.0;

  // Time series of drag and lift coefficients and pressure difference.
  std::vector<double> time_history;
  std::vector<double> cL_history;
  std::vector<double> cD_history;
  std::vector<double> dP_history;


  // Discretization. ///////////////////////////////////////////////////////////

  // Mesh file name.
  const std::string mesh_file_name;

  // Polynomial degree used for velocity.
  const unsigned int degree_velocity;

  // Polynomial degree used for pressure.
  const unsigned int degree_pressure;
  // Final time.
  const double T;

  // Theta parameter for the theta method.
  const double theta;

  // Time step.
  const double delta_t;

  // Data for the problem.
  ProblemData data;

  // Current time.
  double time = 0.0;

  // Current timestep number.
  unsigned int timestep_number = 0;

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

  // Pressure mass matrix, needed for preconditioning. 
  TrilinosWrappers::BlockSparseMatrix pressure_mass;

  // Velocity mass matrix for Yosida preconditioner.
  TrilinosWrappers::BlockSparseMatrix velocity_mass;

  // Base pointer for managing any block preconditioner. 
  std::unique_ptr<NSPreconditioners::NavierStokesPreconditionerBase> advanced_preconditioner;

  // Right-hand side vector in the linear system.
  TrilinosWrappers::MPI::BlockVector system_rhs;

  // System solution (without ghost elements).
  TrilinosWrappers::MPI::BlockVector solution_owned;

  // System solution (including ghost elements).
  TrilinosWrappers::MPI::BlockVector solution;

  // Physical constants for the benchmark.
  static constexpr double D_cylinder = 0.1;  // 2*R_cylinder
  static constexpr double H_channel  = 0.41; // channel height

  // 3D benchmark points.
  static constexpr double x_front_cylinder = 0.45;   // cx - R_cylinder
  static constexpr double y_probe          = 0.20;   // cy
  static constexpr double z_probe          = 0.205;  // W\2
  static constexpr double x_back_cylinder  = 0.55;   // cx + R_cylinder

  // Search interval for the recirculation zone.
  static constexpr double x_wake_start = x_back_cylinder;
  static constexpr double x_wake_end   = 2.20;

  double reference_velocity() const
  {
    // From the benchmark: U = 4 U(0,H/2,W/2,t) / 9.
    return 4.0 * data.Um / 9.0;
  }

};

#endif