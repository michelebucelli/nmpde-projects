#ifndef HEAT_HPP
#define HEAT_HPP

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/data_out_base.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>

using namespace dealii;

/**
 * Class managing the differential problem.
 */
class Heat
{
public:
  // Physical dimension (1D, 2D, 3D)
  static constexpr unsigned int dim = 3;

  // Initial condition.
  class FunctionU0 : public Function<dim>
  {
  public:
    // Constructor.
    FunctionU0() = default;

    // Evaluation of the function. u(t = 0) = 0, per Project 5.
    virtual double
    value(const Point<dim> &, const unsigned int /*component*/ = 0) const override
    {
      return 0.0;
    }
  };

  class g_function : public Function<dim>
  {
  public:
    g_function() : Function<dim>() {}

    virtual double
    value(const Point<dim> &, const unsigned int /*component*/ = 0) const override
    {
      const double t = this->get_time();
      const double a = 5.0; // Parametro 'a' dal progetto
      const double N = 5.0; // Numero di impulsi

      // g(t) = exp(-a * cos(2 * N * pi * t)) / exp(a)
      return std::exp(-a * std::cos(2.0 * N * numbers::PI * t)) / std::exp(a);
    }
  };

  class h_function : public Function<dim>
  {
  public:
    h_function() : Function<dim>() {}

    virtual double
    value(const Point<dim> &p, const unsigned int /*component*/ = 0) const override
    {
      Point<dim> x0_point;
      for (unsigned int d = 0; d < dim; ++d)
        x0_point[d] = 0.5;

      const double sigma_val = 0.1;
      // h(x) = exp(-(x-x0)^2 / sigma^2)
      return std::exp(-std::pow(p.distance(x0_point), 2) / std::pow(sigma_val, 2));
    }
  };

  // Constructor.
  Heat(const unsigned int                                              &n_elements_,
       const unsigned int                                              &r_,
       const double                                                    &T_,
       const double                                                    &theta_,
       const double                                                    &delta_t_,
       const std::function<double(const Point<dim> &)>                 &mu_,
       const std::function<double(const Point<dim> &, const double &)> &f_,
       const bool                                                      &adaptive_      = true,
       const double                                                    &tol_time_max_  = 1e-4,
       const double                                                    &tol_time_min_  = 2e-6,
       const double                                                    &dt_max_        = 0.01)
    : n_elements(n_elements_)
    , r(r_)
    , T(T_)
    , theta(theta_)
    , delta_t(delta_t_)
    , adaptive(adaptive_)
    , tol_time_max(tol_time_max_)
    , tol_time_min(tol_time_min_)
    , dt_max(dt_max_)
    , mu(mu_)
    , f(f_)
    , mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
    , mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    , mesh(MPI_COMM_WORLD)
    , pcout(std::cout, mpi_rank == 0)
  {}

  // Run the time-dependent simulation.
  void
  run();

  // Compute the error against a given exact solution.
  double
  compute_error(const VectorTools::NormType &norm_type,
                const Function<dim>         &exact_solution) const;

  const DoFHandler<dim> &
  get_dof_handler() const
  {
    return dof_handler;
  }

  const TrilinosWrappers::MPI::Vector &
  get_solution() const
  {
    return solution;
  }

protected:
  // Initialization.
  void
  setup_system();

  // System assembly.
  void
  assemble();

  // System solution.
  void
  solve_linear_system();

  // La funzione core per l'adattività
  void
  refine_mesh();

  // Output.
  void
  output() const;

  double
  compute_time_error_indicator() const;
  double
  compute_forcing_integral(const double t) const;
  double
  compute_mass() const;

  // Basename shared by output() and run(), so adaptive/non-adaptive runs at
  // the same n_elements never collide.
  std::string
  run_basename() const
  {
    return "n" + std::to_string(n_elements) + (adaptive ? "-adaptive" : "-nonadaptive");
  }

  const unsigned int n_elements;

  // Polynomial degree.
  const unsigned int r;

  // Final time.
  const double T;

  // Theta parameter for the theta method.
  const double theta;

  // Time step.
  double delta_t;

  // If false, the mesh is never refined and delta_t stays fixed at its
  // constructor value for the whole run (the non-adaptive baseline).
  const bool adaptive;

  // Per-instance, since the pulsed problem needs tight tolerances to avoid
  // missing a narrow pulse, but a smooth problem (e.g. manufactured-solution-test)
  // has no sharp feature to miss and would just pay the extra rejection cost.
  const double tol_time_max; // Soglia massima errore temporale
  const double tol_time_min; // Soglia minima per aumentare delta_t
  const double dt_max;       // Limite superiore per delta_t
  const double dt_min = 1e-5; // Limite inferiore per delta_t

  const unsigned int max_grid_level = 2;

  // Current time.
  double time = 0.0;

  // Current timestep number.
  unsigned int timestep_number = 0;

  // Diffusion coefficient.
  std::function<double(const Point<dim> &)> mu;

  // Forcing term.
  std::function<double(const Point<dim> &, const double &)> f;

  // Number of MPI processes.
  const unsigned int mpi_size;

  // Rank of the current MPI process.
  const unsigned int mpi_rank;

  // Triangulation.
  parallel::distributed::Triangulation<dim> mesh;

  // Finite element space.
  std::unique_ptr<FiniteElement<dim>> fe;

  // Quadrature formula.
  std::unique_ptr<Quadrature<dim>> quadrature;

  // DoF handler.
  DoFHandler<dim> dof_handler;

  // System matrix.
  TrilinosWrappers::SparseMatrix system_matrix;

  // System right-hand side.
  TrilinosWrappers::MPI::Vector system_rhs;

  // System solution, without ghost elements.
  TrilinosWrappers::MPI::Vector solution_owned;
  TrilinosWrappers::MPI::Vector solution_owned_old;

  // System solution, with ghost elements.
  TrilinosWrappers::MPI::Vector solution;
  TrilinosWrappers::MPI::Vector solution_old;

  // Output stream for process 0.
  ConditionalOStream pcout;
};

#endif
