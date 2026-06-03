#ifndef MATRIX_BASED_ADR_HPP
#define MATRIX_BASED_ADR_HPP


#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include <iostream>

#include "AdrSpec.hpp"

using namespace dealii;


template <int dim, typename AdrSpec>
class MatrixBasedAdr {
	public:
		MatrixBasedAdr(const Triangulation<dim>& mesh,
		               int fe_degree);

		void setup_system();
		void assemble_system();
		int  solve(int nmax, double rtol, bool apply_preconditioner);

		double compute_error(const VectorTools::NormType &norm_type,
		                     const Function<dim>         &exact_solution,
		                     const Triangulation<dim>    &mesh) const;

	protected:
		const int          fe_degree;
		const FE_Q<dim>    fe;
		DoFHandler<dim>    dof_handler;

		SparsityPattern      sparsity_pattern;
		SparseMatrix<double> system_matrix;

		Vector<double> solution;
		Vector<double> system_rhs;

		AdrSpec local_spec;
};


template <int dim, typename AdrSpec>
MatrixBasedAdr<dim, AdrSpec>::
MatrixBasedAdr(const Triangulation<dim>& mesh,
               int fe_degree_)
	: fe_degree(fe_degree_)
	, fe(fe_degree)
	, dof_handler(mesh)
{}


template <int dim, typename AdrSpec>
void
MatrixBasedAdr<dim, AdrSpec>::
setup_system() {
	dof_handler.distribute_dofs(fe);

	std::cout << "Number of DoF: "
	          << dof_handler.n_dofs()
	          << std::endl;

	DynamicSparsityPattern dsp(dof_handler.n_dofs());
	DoFTools::make_sparsity_pattern(dof_handler, dsp);
	sparsity_pattern.copy_from(dsp);

	system_matrix.reinit(sparsity_pattern);

	solution.reinit(dof_handler.n_dofs());
	system_rhs.reinit(dof_handler.n_dofs());
}


template <int dim, typename AdrSpec>
void
MatrixBasedAdr<dim, AdrSpec>::
assemble_system()
{
	// TODO: quadrature_formula could be a class member
	const QGauss<dim> quadrature_formula(fe.degree + 1);
	FEValues<dim> fe_values(fe,
	                        quadrature_formula,
	                        ( update_values
	                        | update_gradients
	                        | update_JxW_values
	                        | update_quadrature_points ));

	const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
	const unsigned int n_q = quadrature_formula.size();

	FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
	Vector<double>     cell_rhs   (dofs_per_cell);

	std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

	for (const auto &cell : dof_handler.active_cell_iterators()) {
		fe_values.reinit(cell);

		cell_matrix = 0;
		cell_rhs = 0;

		for (unsigned int q = 0; q < n_q; ++q) {
			const Point<dim> &x_q = fe_values.quadrature_point(q);

			const double         mu    = local_spec.diffusion(x_q);
			const Tensor<1, dim> beta  = local_spec.advection(x_q);
			const double         gamma = local_spec.reaction(x_q);
			const double         f     = local_spec.forcing(x_q);

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

		for (const auto &face : cell->face_iterators()) {
			if (face->at_boundary()) {

			}
		}

		cell->get_dof_indices(local_dof_indices);
		system_matrix.add(local_dof_indices, cell_matrix);
		system_rhs.add(local_dof_indices, cell_rhs);
	}

	std::map<types::global_dof_index, double> boundary_values;


	BoundaryFuncMap<dim> tmp_dirichlet_map;
	std::map<types::boundary_id, const Function<dim> *> dirichlet_map;
	local_spec.fill_dirichlet_map(tmp_dirichlet_map);
	for (const auto& [boundary_id, dirichlet_function] : tmp_dirichlet_map) {
		dirichlet_map[boundary_id] = dirichlet_function.get();
	}

	VectorTools::interpolate_boundary_values(dof_handler,
	                                         dirichlet_map,
	                                         boundary_values);
	MatrixTools::apply_boundary_values(boundary_values,
	                                   system_matrix,
	                                   solution,
	                                   system_rhs);
}


template <int dim, typename AdrSpec>
int
MatrixBasedAdr<dim, AdrSpec>::
solve(int nmax, double rtol, bool apply_preconditioner) {
	SolverControl solver_control(nmax, rtol * system_rhs.l2_norm());
	SolverGMRES<Vector<double>> solver(solver_control);

	if (not apply_preconditioner) {
		solver.solve(system_matrix,
		             solution,
		             system_rhs,
		             PreconditionIdentity());
	} else {
		PreconditionJacobi preconditioner;
		preconditioner.initialize(system_matrix, 1.0);

		solver.solve(system_matrix,
		             solution,
		             system_rhs,
		             preconditioner);
	}

	return solver_control.last_step();
}


template <int dim, typename AdrSpec>
double
MatrixBasedAdr<dim, AdrSpec>::
compute_error(const VectorTools::NormType &norm_type,
              const Function<dim>         &exact_solution,
              const Triangulation<dim>    &mesh) const {

	Vector<double> error_per_cell(mesh.n_active_cells());
	VectorTools::integrate_difference(dof_handler,
	                                  solution,
	                                  exact_solution,
	                                  error_per_cell,
	                                  QGauss<dim>(fe_degree + 2),
	                                  norm_type);

	return VectorTools::compute_global_error(mesh, error_per_cell, norm_type);
}


#endif // MATRIX_BASED_ADR_HPP
