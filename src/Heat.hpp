#ifndef HEAT_HPP
#define HEAT_HPP

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

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/numerics/error_estimator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/numerics/solution_transfer.h>

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

    // Evaluation of the function. u( t = 0); p[0] = x, p[1] = y, p[2] = z.
    virtual double
    value(const Point<dim> &p,
          const unsigned int /*component*/ = 0) const override
    {
      return 0;
    }
  };

  class g_function : public Function<dim>
  {
  public:
    g_function() : Function<dim>() {}

    virtual double value(const Point<dim> & /*p*/, const unsigned int /*component*/ = 0) const override
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

    virtual double value(const Point<dim> &p, const unsigned int /*component*/ = 0) const override
    {
      Point<dim> x0_point; 
      x0_point[0] = 0.5; // x dove viene applicato l'impulso, cambiare per test
      
      const double sigma_val = 0.1; // 
      // h(x) = exp(-(x-x0)^2 / sigma^2) 
      return std::exp(-std::pow(p.distance(x0_point), 2) / std::pow(sigma_val, 2));
    }
  };

  // Constructor.
  Heat(const std::string                              &mesh_file_name_, 
      const unsigned int                              &r_,
      const double                                    &T_,
      const double                                    &theta_,
      const double                                    &delta_t_,
      const std::function<double(const Point<dim> &)> &mu_,
      const std::function<double(const Point<dim> &, const double &)> &f_)
    : mesh_file_name(mesh_file_name_)
    , r(r_)
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

  // Run the time-dependent simulation.
  void
  run();

  // Compute the error against a given exact solution.
  double
  compute_error(const VectorTools::NormType &norm_type,
                const Function<dim>         &exact_solution) const;

protected:
  
  /*
  // Initialization.
  void
  setup();
  */

  // System assembly.
  void
  assemble();

  // System solution.
  void
  solve_linear_system();

  void setup_system(); // Ex setup(), ora gestisce solo matrici/DoF
  void init_mesh();    // Nuova funzione per leggere la mesh solo una volta
  void refine_mesh();  // La funzione core per l'adattività

  // Output.
  void
  output() const;

  // Name of the mesh.
  const std::string mesh_file_name;

  // Polynomial degree.
  const unsigned int r;

  // Final time.
  const double T;

  // Theta parameter for the theta method.
  const double theta;

  // Time step.
  double delta_t;

  
  const double tol_time_max = 5e-3; // Soglia massima errore temporale
  const double tol_time_min = 1e-4; // Soglia minima per aumentare delta_t
  const double dt_min = 1e-5;       // Limite inferiore per delta_t
  const double dt_max = 0.1;        // Limite superiore per delta_t

  // Current time.
  double time = 0.0;

  // Current timestep number.
  unsigned int timestep_number = 0;

  // Diffusion coefficient.
  std::function<double(const Point<dim> &)> mu;

   // Reaction coefficient.
  std::function<double(const Point<dim> &)> sigma;

  // Forcing term.
  std::function<double(const Point<dim> &, const double &)> f;

  // Advection coefficient.
  std::function<Tensor<1, dim>(const Point<dim> &)> b;

  // Number of MPI processes.
  const unsigned int mpi_size;

  // Rank of the current MPI process.
  const unsigned int mpi_rank;

  // Triangulation.
  parallel::fullydistributed::Triangulation<dim> mesh;

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
  
  // Vettore per salvare la soluzione prima della rifinitura <-----
  TrilinosWrappers::MPI::Vector solution_owned_old;

  // System solution, without ghost elements.
  TrilinosWrappers::MPI::Vector solution_owned;

  // System solution, with ghost elements.
  TrilinosWrappers::MPI::Vector solution;

  // Output stream for process 0.
  ConditionalOStream pcout;
};

#endif
