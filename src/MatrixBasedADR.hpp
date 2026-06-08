#ifndef MATRIX_BASED_ADR_HPP
#define MATRIX_BASED_ADR_HPP

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include <iostream>

#include "ADR.hpp"
#include "ADRSpec.hpp"
#include "ADRFactory.hpp"

using namespace dealii;

// WARNING: TrilinosWrappers::MPI::Vector and TrilinosWrappers::SparseMatrix
// use TrilinosScalar as value_type which is hardcoded to double.
// We can't use the spec defined ADRSpec::ScalarType.

template<typename ADRSpec>
class MatrixBasedADR :
    public ADR<
        ADRSpec, typename ADRSpec::ScalarType, TrilinosWrappers::MPI::Vector>
{
    using number = typename ADRSpec::ScalarType;
    // We need to explicitly bring into this scope the base template
    // class members due to the two-phase lookup (or use this->member)
    using ADR<ADRSpec, number, TrilinosWrappers::MPI::Vector>::dim;

    friend class ADRFactory;

    MatrixBasedADR(const parallel::distributed::Triangulation<dim> &mesh,
                   unsigned int fe_degree);

public:
    void setup() override;

    void assemble() override;

    unsigned int solve(unsigned int max_iter,
                       number rel_tol,
                       bool apply_preconditioner) override;

protected:
    IndexSet locally_owned_dofs;

    TrilinosWrappers::SparseMatrix system_matrix;

    TrilinosWrappers::MPI::Vector solution_owned;
};

template<typename ADRSpec>
MatrixBasedADR<ADRSpec>::MatrixBasedADR(
    const parallel::distributed::Triangulation<dim> &mesh,
    unsigned int fe_degree
)
    : ADR<ADRSpec, number, TrilinosWrappers::MPI::Vector>(mesh, fe_degree)
{}

template<typename ADRSpec>
void MatrixBasedADR<ADRSpec>::setup()
{
    this->dof_handler.distribute_dofs(this->fe);
    locally_owned_dofs = this->dof_handler.locally_owned_dofs();

    TrilinosWrappers::SparsityPattern sparsity(locally_owned_dofs,
                                               MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(this->dof_handler, sparsity);
    sparsity.compress();

    system_matrix.reinit(sparsity);

    IndexSet locally_relevant_dofs;
    DoFTools::extract_locally_relevant_dofs(this->dof_handler,
                                            locally_relevant_dofs);

    this->solution.reinit(
        locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
    solution_owned.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    this->system_rhs.reinit(locally_owned_dofs, MPI_COMM_WORLD);
}

template<typename ADRSpec>
void MatrixBasedADR<ADRSpec>::assemble()
{
    // TODO: quadratures could be class members
    const QGauss<dim> quadrature(this->fe.degree + 1);
    const QGauss<dim-1> quadrature_face(this->fe.degree + 1);

    FEValues<dim> fe_values(this->fe,
                            quadrature,
                            update_values | update_gradients |
                            update_JxW_values | update_quadrature_points);

    FEFaceValues<dim> fe_face_values(this->fe,
                                     quadrature_face,
                                     update_values | update_JxW_values |
                                     update_quadrature_points);

    const unsigned int dofs_per_cell = this->fe.n_dofs_per_cell();
    const unsigned int n_q = quadrature.size();
    const unsigned int n_q_face = quadrature_face.size();

    FullMatrix<number> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<number> cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    BoundaryFuncMap<NeumannFunction<number, Point<dim, number>>>
    neumann_func_map;
    this->local_spec.fill_neumann_map(neumann_func_map);

    for (const auto &cell : this->dof_handler.active_cell_iterators()) {
        if (!cell->is_locally_owned()) {
            continue;
        }

        fe_values.reinit(cell);

        cell_matrix = 0;
        cell_rhs = 0;

        for (unsigned int q = 0; q < n_q; ++q) {
            const Point<dim> &x_q = fe_values.quadrature_point(q);

            const typename ADRSpec::TensorType beta =
                this->local_spec.advection(x_q);
            const number mu = this->local_spec.diffusion(x_q);
            const number gamma = this->local_spec.reaction(x_q);
            const number f = this->local_spec.forcing(x_q);

            for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                for (unsigned int j = 0; j < dofs_per_cell; ++j) {

                    cell_matrix(i, j) += (
                          mu * fe_values.shape_grad(i, q)
                             * fe_values.shape_grad(j, q)

                        + beta * fe_values.shape_grad(j, q)
                               * fe_values.shape_value(i, q)

                        + gamma * fe_values.shape_value(i, q)
                                * fe_values.shape_value(j, q)
                    ) * fe_values.JxW(q);
                }

                cell_rhs(i) += f * fe_values.shape_value(i, q)
                                 * fe_values.JxW(q);
            }
        }

        if (cell->at_boundary()) {
            for (const auto &face : cell->face_iterators()) {
                const types::boundary_id b_id = face->boundary_id();

                if (face->at_boundary() && neumann_func_map.count(b_id)) {
                    fe_face_values.reinit(cell, face);

                    for (unsigned int q = 0; q < n_q_face; ++q) {
                        const Point<dim> &x_q =
                            fe_face_values.quadrature_point(q);

                        const number mu = this->local_spec.diffusion(x_q);
                        const number h = neumann_func_map[b_id]->value(x_q);

                        for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                            cell_rhs(i) +=
                                mu * h * fe_face_values.shape_value(i, q)
                                       * fe_face_values.JxW(q);
                        }
                    }
                }
            }
        }

        cell->get_dof_indices(local_dof_indices);
        system_matrix.add(local_dof_indices, cell_matrix);
        this->system_rhs.add(local_dof_indices, cell_rhs);
    }

    system_matrix.compress(VectorOperation::add);
    this->system_rhs.compress(VectorOperation::add);

    std::map<types::global_dof_index, number> boundary_values;

    BoundaryFuncMap<Function<dim>> tmp_dirichlet_map;
    std::map<types::boundary_id, const Function<dim> *> dirichlet_map;
    this->local_spec.fill_dirichlet_map(tmp_dirichlet_map);
    for (const auto& [boundary_id, dirichlet_function] : tmp_dirichlet_map) {
        dirichlet_map[boundary_id] = dirichlet_function.get();
    }

    VectorTools::interpolate_boundary_values(this->dof_handler,
                                             dirichlet_map,
                                             boundary_values);

    MatrixTools::apply_boundary_values(boundary_values,
                                       system_matrix,
                                       solution_owned,
                                       this->system_rhs);
}

template<typename ADRSpec>
unsigned int MatrixBasedADR<ADRSpec>::solve(
    unsigned int max_iter, number rel_tol, bool apply_preconditioner)
{
    SolverControl solver_control(max_iter,
                                 rel_tol * this->system_rhs.l2_norm());
    SolverGMRES<TrilinosWrappers::MPI::Vector> solver(solver_control);

    if (not apply_preconditioner) {
        solver.solve(system_matrix,
                     solution_owned,
                     this->system_rhs,
                     TrilinosWrappers::PreconditionIdentity());
    } else {
        TrilinosWrappers::PreconditionJacobi preconditioner;
        preconditioner.initialize(system_matrix, 1.0);

        solver.solve(system_matrix,
                     solution_owned,
                     this->system_rhs,
                     preconditioner);
    }

    this->solution = solution_owned;

    return solver_control.last_step();
}

#endif
