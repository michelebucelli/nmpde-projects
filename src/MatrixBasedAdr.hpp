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
		const FE_Q<dim>    fe;
		const int          fe_degree;
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
	: dof_handler(mesh)
	, fe(fe_degree_)
	, fe_degree(fe_degree_) {
}


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
assemble_system() {
	const QGauss<dim> quadrature_formula(fe.degree + 1);
	FEValues<dim> fe_values(fe,
	                        quadrature_formula,
	                        ( update_values
	                        | update_gradients
	                        | update_JxW_values
	                        | update_quadrature_points ));

	const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

	FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
	Vector<double>     cell_rhs   (dofs_per_cell);

	std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

	for (const auto &cell : dof_handler.active_cell_iterators()) {
		fe_values.reinit(cell);

		cell_matrix = 0;
		cell_rhs = 0;

		for (const auto q_index : fe_values.quadrature_point_indices()) {
			const Point<dim> &x_q = fe_values.quadrature_point(q_index);

			const double         mu    = local_spec.mu(x_q);
			const Tensor<1, dim> beta  = local_spec.beta(x_q);
			const double         gamma = local_spec.gamma(x_q);
			const double         f     = local_spec.forcing(x_q);

			for (const auto i : fe_values.dof_indices()) {
				for (const auto j : fe_values.dof_indices()) {
					cell_matrix(i,j) += (
						  mu * fe_values.shape_grad(i, q_index)
						     * fe_values.shape_grad(j, q_index)

						+ (beta * fe_values.shape_grad(j, q_index))
						        * fe_values.shape_value(i, q_index)

						+ gamma * fe_values.shape_value(i, q_index)
						        * fe_values.shape_value(j, q_index)
					) * fe_values.JxW(q_index);
				}

				cell_rhs(i) += (f * fe_values.shape_value(i, q_index))
				                * fe_values.JxW(q_index);
			}
		}

		for (const auto &face : cell->face_iterators()) {
			if (face->at_boundary()) {

			}
		}

		cell->get_dof_indices(local_dof_indices);

		for (const auto i : fe_values.dof_indices()) {
			for (const auto j : fe_values.dof_indices()) {
				system_matrix.add(local_dof_indices[i],
				                  local_dof_indices[j],
				                  cell_matrix(i,j));
			}

			system_rhs(local_dof_indices[i]) += cell_rhs(i);
		}
	}

	std::map<types::global_dof_index, double> boundary_values;


	std::map<types::boundary_id, std::unique_ptr<Function<dim>>> tmp_dirichlet_map;
	std::map<types::boundary_id, const Function<dim> *> dirichlet_map;
	local_spec->fill_dirichlet_map(tmp_dirichlet_map);
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
	solver.solve(system_matrix,
	             solution,
	             system_rhs,
	             PreconditionIdentity());

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
