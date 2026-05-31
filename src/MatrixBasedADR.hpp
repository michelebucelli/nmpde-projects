#ifndef MATRIX_BASED_ADR_HPP
#define MATRIX_BASED_ADR_HPP


#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/numerics/data_out.h>


using namespace dealii;


template <int dim>
class MatrixBasedADR
{
public:
    MatrixBasedADR(int n, int fe_degree);

    void make_grid();
    void setup_system();
    void assemble_system();

    void assemble();
    int solve(int nmax, double rtol);
    void output() const;

    double
    compute_error(const VectorTools::NormType &norm_type,
                  const Function<dim>         &exact_solution) const;

protected:
    const int          m_n;
    const int          m_fe_degree;
    Triangulation<dim> m_triangulation;
    const FE_Q<dim>    m_fe;
    DoFHandler<dim>    m_dof_handler;

    SparsityPattern    m_sparsity_pattern;
    SparseMatrix<double> m_system_matrix;

    Vector<double> m_solution;
    Vector<double> m_system_rhs;
};


template <int dim>
MatrixBasedADR<dim>::MatrixBasedADR(int n, int fe_degree)
  : m_n(n)
  , m_fe_degree(fe_degree)
  , m_fe(fe_degree)
  , m_dof_handler(m_triangulation)
{}


template <int dim>
void
MatrixBasedADR<dim>::make_grid() {
    GridGenerator::subdivided_hyper_cube(
        m_triangulation,
        m_n
    );
}


template <int dim>
void
MatrixBasedADR<dim>::setup_system() {
    m_dof_handler.distribute_dofs(m_fe);

    DynamicSparsityPattern dsp(m_dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(m_dof_handler, dsp);
    m_sparsity_pattern.copy_from(dsp);

    m_system_matrix.reinit(m_sparsity_pattern);

    m_solution.reinit(m_dof_handler.n_dofs());
    m_system_rhs.reinit(m_dof_handler.n_dofs());
}


template <int dim>
void
MatrixBasedADR<dim>::assemble_system() {
    const QGauss<dim> quadrature_formula(m_fe.degree + 1);
    FEValues<dim> fe_values(
        m_fe,
        quadrature_formula,
        update_values | update_gradients | update_JxW_values | update_quadrature_points
    );

    const unsigned int dofs_per_cell = m_fe.n_dofs_per_cell();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : m_dof_handler.active_cell_iterators()) {
        fe_values.reinit(cell);

        cell_matrix = 0;
        cell_rhs = 0;

        for (const unsigned int q_index : fe_values.quadrature_point_indices()) {
            for (const unsigned int i : fe_values.dof_indices()) {
                for (const unsigned int j : fe_values.dof_indices()) {
                    // TODO: insert here the weak formulation for the ADR problem
                    cell_matrix(i,j) += 0;
                }


                // TODO: insert here the rhs part of the weak formulation
                // Not the boundary term, though...
                const auto &x_q = fe_values.quadrature_point(q_index);
                cell_rhs(i) += 0;
            }
        }

        cell->get_dof_indices(local_dof_indices);

        for (const unsigned int i : fe_values.dof_indices()) {
            for (const unsigned int j : fe_values.dof_indices()) {
                m_system_matrix.add(local_dof_indices[i],
                                    local_dof_indices[j],
                                    cell_matrix(i, j));
            }

            m_system_rhs(local_dof_indices[i]) += cell_rhs(i);
        }
    }

    // here we have to add the boundary terms to the system matrix
    // we have Dirichlet boundary conditions and Neumann ones
    std::map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values(m_dof_handler,
                                             types::boundary_id(0),
                                             Functions::ZeroFunction<dim>(),
                                             boundary_values);
    MatrixTools::apply_boundary_values(boundary_values,
                                       m_system_matrix,
                                       m_solution,
                                       m_system_rhs);
}


template <int dim>
int
MatrixBasedADR<dim>::solve(int nmax, double rtol) {
    SolverControl               solver_control(nmax, rtol * m_system_rhs.l2_norm());
    SolverGMRES<Vector<double>> solver(solver_control);
    solver.solve(m_system_matrix, m_solution, m_system_rhs, PreconditionIdentity());

    return solver_control.last_step();
}


template <int dim>
double
MatrixBasedADR<dim>::compute_error(const VectorTools::NormType &norm_type,
                                   const Function<dim>         &exact_solution) const {

    Vector<double> error_per_cell(m_triangulation.n_active_cells());
    VectorTools::integrate_difference(
        m_dof_handler,
        m_solution,
        exact_solution,
        error_per_cell,
        QGauss<dim>(m_fe_degree + 2),
        norm_type
    );

    return VectorTools::compute_global_error(m_triangulation, error_per_cell, norm_type);
}


#endif // MATRIX_BASED_ADR_HPP
