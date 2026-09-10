#ifndef HETEROGENEOUS_DIFFUSION_HPP
#define HETEROGENEOUS_DIFFUSION_HPP

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_description.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <string>
#include <vector>

using namespace dealii;

/**
 * 3D Poisson with heterogeneous diffusion coefficient, discretized on a
 * tetrahedral (P1 simplex) mesh.
 *
 *   -div(mu grad u) = f    in Omega = (0, 1)^3,
 *              u = 0       on d Omega,
 *
 * with mu(x) = 10^p if x lies in any of a fixed set of spheres, 1 otherwise,
 * and f = 1.
 */
class HeterogeneousDiffusion
{
public:
  static constexpr unsigned int dim = 3;

  // Ball struct for the heterogeneous domain. 
  struct Ball
  {
    Point<dim> center;
    double     radius;
  };

  // Diffusion coefficient: 10^p inside any ball, 1 otherwise.
  class DiffusionCoefficient : public Function<dim>
  {
  public:
    DiffusionCoefficient(const double p_, const std::vector<Ball> &balls_)
      : p(p_)
      , balls(balls_)
    {}

    double
    value(const Point<dim>  &point,
          const unsigned int /*component*/ = 0) const override
    {
      for (const auto &ball : balls)
        if (point.distance(ball.center) <= ball.radius)
          return std::pow(10.0, p);
      return 1.0;
    }

  private:
    const double             p;
    const std::vector<Ball> &balls;
  };

  // Forcing term (constant).
  class RightHandSide : public Function<dim>
  {
  public:
    double
    value(const Point<dim> & /*point*/,
          const unsigned int /*component*/ = 0) const override
    {
      return 1.0;
    }
  };

    class ManufacturedSolution : public Function<dim>
    {
    public:
      double
      value(const Point<dim>  &point,
            const unsigned int /*component*/ = 0) const override
      {
        return std::sin(M_PI * point[0]) * std::sin(M_PI * point[1]) *
               std::sin(M_PI * point[2]);
      }

      Tensor<1, dim>
      gradient(const Point<dim>  &point,
               const unsigned int /*component*/ = 0) const override
      {
        Tensor<1, dim> g;
        const double   sx = std::sin(M_PI * point[0]);
        const double   sy = std::sin(M_PI * point[1]);
        const double   sz = std::sin(M_PI * point[2]);
        const double   cx = std::cos(M_PI * point[0]);
        const double   cy = std::cos(M_PI * point[1]);
        const double   cz = std::cos(M_PI * point[2]);
        g[0] = M_PI * cx * sy * sz;
        g[1] = M_PI * sx * cy * sz;
        g[2] = M_PI * sx * sy * cz;
        return g;
      }
    };

    class ManufacturedRightHandSide : public Function<dim>
    {
    public:
      double
      value(const Point<dim>  &point,
            const unsigned int /*component*/ = 0) const override
      {
        return 3.0 * M_PI * M_PI * std::sin(M_PI * point[0]) *
               std::sin(M_PI * point[1]) * std::sin(M_PI * point[2]);
      }
    };

  // Constructor.
    HeterogeneousDiffusion(const unsigned int refinements_,
                           const double       p_,
                           const MPI_Comm     mpi_communicator_,
                           const bool         verify_mode_ = false);

  // Mesh, DoFs, sparsity, and matrix/rhs assembly. Call once per case.
  void
  setup_and_assemble();

  void
  solve(const std::string &preconditioner_name);

  double
  compute_error(const VectorTools::NormType &norm_type) const;

  void
  output(const std::string &tag) const;

  static std::vector<Ball>
  make_balls();

protected:
  template <class Preconditioner>
  unsigned int
  solve_with(const Preconditioner &preconditioner,
             const std::string    &name,
             const double          setup_time);

  // Appends one row of solver results to results.csv (rank 0 only),
  // writing the header the first time the file is created.
  void
  log_result_csv(const std::string &name,
                 const unsigned int iterations,
                 const double       setup_time,
                 const double       solve_time,
                 const double       condition_number) const;

  const unsigned int refinements;
  const double       p;
  MPI_Comm           mpi_communicator;
  const unsigned int mpi_size;
  const unsigned int mpi_rank;
  const bool verify_mode;
  ConditionalOStream pcout;

  // Tracks whether the results-table header has been printed for the
  // current case, so it only appears once even though solve_with() is a
  // template instantiated separately per preconditioner type.
  mutable bool result_header_printed = false;

  std::vector<Ball> balls;

  parallel::fullydistributed::Triangulation<dim> mesh;
  FE_SimplexP<dim>                     fe;
  const MappingFE<dim>                 mapping;
  DoFHandler<dim>                      dof_handler;
  AffineConstraints<double>            constraints;

  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;

  TrilinosWrappers::SparsityPattern sparsity_pattern;
  TrilinosWrappers::SparseMatrix    system_matrix;
  TrilinosWrappers::MPI::Vector     system_rhs;
  TrilinosWrappers::MPI::Vector     solution;
};

#endif