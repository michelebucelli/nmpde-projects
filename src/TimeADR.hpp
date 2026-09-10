#ifndef TIMEADR_HPP
#define TIMEADR_HPP

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <variant>

using namespace dealii;

/**
 Finite element solver for the 2D wave equation

 Key implementation choices:
  - M, A and the system matrix (1/dt^2) M + theta A are assembled ONCE, since
    mesh, mu, dt and theta are constant in time. Each time step only needs a
    few matrix-vector products to build the right-hand side
  - The SPD system solved with CG
 */
class TimeADR
{
public:
  static constexpr unsigned int dim = 2;

  class FunctionDirichlet : public Function<dim>
  {
  public:
    FunctionDirichlet() = default;

    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int /*component*/ = 0) const override
    {
      return 0.0;
    }
  };

  // Constructor
  TimeADR(const unsigned int                              &N_el_,
          const unsigned int                              &r_,
          const double                                    &T_,
          const double                                    &theta_,
          const double                                    &delta_t_,
          const std::function<double(const Point<dim> &)> &mu_,
          const std::function<double(const Point<dim> &, const double &)> &f_,
          const std::function<double(const Point<dim> &)>                 &u0_,
          const std::function<double(const Point<dim> &)>                 &u1_,
          const std::function<double(const Point<dim> &)>                 &psi_)
    : N_el(N_el_)
    , mesh_spec(N_el_)
    , r(r_)
    , T(T_)
    , theta(theta_)
    , delta_t(delta_t_)
    , u0_fun(u0_)
    , u1_fun(u1_)
    , mu(mu_)
    , f(f_)
    , psi(psi_)
    , mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
    , mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    , mesh(MPI_COMM_WORLD)
    , pcout(std::cout, mpi_rank == 0)
  {}

  void
  run();

  double
  compute_error(const VectorTools::NormType &norm_type,
                Function<dim>               &exact_solution) const;

  void
  set_verbose(const bool v)
  {
    verbose = v;
  }
  void
  set_write_output(const bool w)
  {
    write_output = w;
  }
  void
  set_output_every(const unsigned int k)
  {
    output_every = k;
  }

  void
  set_dirichlet_datum(std::function<double(const Point<dim> &, double)> g)
  {
    g_fun = std::move(g);
  }

protected:
  // Initialization: mesh, FE space, DoFs, constraints, the constant matrices
  // M, A and (1/dt^2) M + theta A, and the preconditioner.
  void
  setup();

  // Assemble the constant matrices (called once from setup()).
  void
  assemble_matrices();

  void
  assemble_full_matrices();

  // Build the discrete lifting vector R(t): equal to the Dirichlet datum g on
  // the boundary and 0 in the interior, at time t. Non-homogeneous case only.
  void
  build_lifting(const double t, TrilinosWrappers::MPI::Vector &Rg) const;

  void
  assemble_rhs();

  void
  solve_linear_system();

  void
  output() const;

  const unsigned int                            N_el;
  const std::variant<unsigned int, std::string> mesh_spec;

  const unsigned int r; // grado polinomiale

  const double T;

  const double theta; // parametro di Newmark

  const double delta_t;

  std::function<double(const Point<dim> &)> u0_fun;
  std::function<double(const Point<dim> &)> u1_fun;

  unsigned int output_every = 10; // ogni quanto scrive l'output

  bool verbose      = true;
  bool write_output = true;

  double       time            = 0.0;
  unsigned int timestep_number = 0;

  std::function<double(const Point<dim> &)>
    mu; // COEFF DI DIFFUSIONE (quadrato della wave speed, c^2)

  std::function<double(const Point<dim> &, const double &)> f;

  std::function<double(const Point<dim> &)> psi; // Neumann (non usato)

  std::function<double(const Point<dim> &, double)> g_fun;

  const unsigned int mpi_size;
  const unsigned int mpi_rank;

  parallel::fullydistributed::Triangulation<dim> mesh;

  std::unique_ptr<FiniteElement<dim>>  fe;
  std::unique_ptr<Quadrature<dim>>     quadrature;
  std::unique_ptr<Quadrature<dim - 1>> quadrature_face;

  DoFHandler<dim> dof_handler;

  AffineConstraints<double> constraints;

  AffineConstraints<double> constraints_hn;

  // Matrici costanti
  TrilinosWrappers::SparseMatrix mass_matrix;      // M
  TrilinosWrappers::SparseMatrix stiffness_matrix; // A
  TrilinosWrappers::SparseMatrix lhs_matrix;       // (1/dt^2) M + theta A

  TrilinosWrappers::SparseMatrix mass_full;      // per g di Diriclet
  TrilinosWrappers::SparseMatrix stiffness_full; // per g di Diriclet

  TrilinosWrappers::PreconditionSSOR preconditioner;

  TrilinosWrappers::MPI::Vector system_rhs;

  TrilinosWrappers::MPI::Vector u_old;
  TrilinosWrappers::MPI::Vector u;
  TrilinosWrappers::MPI::Vector u_new;

  TrilinosWrappers::MPI::Vector solution;

  TrilinosWrappers::MPI::Vector tmp;
  TrilinosWrappers::MPI::Vector tmp2;

  ConditionalOStream pcout;

  void
  print_performance_report() const;

  struct PerformanceData
  {
    double setup_seconds     = 0.0;
    double assembly_seconds  = 0.0;
    double solve_seconds     = 0.0;
    double output_seconds    = 0.0;
    double total_run_seconds = 0.0;

    unsigned int       n_steps                 = 0;
    unsigned long long total_solver_iterations = 0;
  };

  mutable PerformanceData perf;
};

#endif