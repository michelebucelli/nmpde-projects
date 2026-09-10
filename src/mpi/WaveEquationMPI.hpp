#ifndef WAVE_EQUATION_MPI_HPP
#define WAVE_EQUATION_MPI_HPP

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <string>
#include <vector>

using namespace dealii;

class WaveEquationMPI
{
public:
  static constexpr unsigned int dim = 2;
  WaveEquationMPI(unsigned int n, unsigned int degree, double final_time,
                  double dt, unsigned int output_frequency,
                  std::string output_directory, std::string case_name);
  void run();

private:
  void setup();
  void assemble_matrices();
  void initialize();
  void assemble_force(double t, TrilinosWrappers::MPI::Vector &) const;
  void solve_initial_acceleration();
  void solve_step(double t);
  void output();
  void diagnostics();
  double energy() const;
  std::pair<double, double> errors() const;
  double omega() const;

  const unsigned int n, degree, output_frequency;
  const double final_time, dt;
  const std::string output_directory, case_name;
  double current_time = 0.;
  unsigned int current_step = 0;
  const unsigned int rank, n_ranks;
  parallel::distributed::Triangulation<dim> mesh;
  FE_Q<dim> fe;
  QGauss<dim> quadrature;
  DoFHandler<dim> dof_handler;
  IndexSet owned, relevant;
  AffineConstraints<double> constraints;
  TrilinosWrappers::SparsityPattern sparsity;
  TrilinosWrappers::SparseMatrix mass, stiffness, effective;
  TrilinosWrappers::MPI::Vector u_owned, v_owned, a_owned, rhs;
  TrilinosWrappers::MPI::Vector u, v, a;
  ConditionalOStream pcout;
};

#endif
