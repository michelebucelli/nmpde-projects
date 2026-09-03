#pragma once

#include "cdr_parameters.hpp"
#include "cdr_problem.hpp"

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

using namespace dealii;

// The matrix-based reference the matrix-free operator is measured against.
// Same bilinear form, same quadrature, same FGMRES and tolerance; the only
// difference is that the operator is assembled into a distributed sparse
// matrix once and then applied as a sparse matrix-vector product.
//
// Inhomogeneous Dirichlet data are lifted by the constraints during
// distribute_local_to_global, rather than by the correction solve the
// matrix-free path uses. Both end at the same discrete solution.
template <int dim, int fe_degree> class AssembledSystem {
public:
  using VectorType = TrilinosWrappers::MPI::Vector;

  auto setup(const DoFHandler<dim> &dof_handler,
             const AffineConstraints<double> &constraints,
             const MPI_Comm communicator) -> void {
    locally_owned = dof_handler.locally_owned_dofs();
    locally_relevant = DoFTools::extract_locally_relevant_dofs(dof_handler);

    DynamicSparsityPattern sparsity(locally_relevant);
    DoFTools::make_sparsity_pattern(dof_handler, sparsity, constraints, false);
    SparsityTools::distribute_sparsity_pattern(sparsity, locally_owned,
                                               communicator, locally_relevant);

    matrix.reinit(locally_owned, locally_owned, sparsity, communicator);
    rhs.reinit(locally_owned, communicator);
    solution.reinit(locally_owned, communicator);
  }

  // Fills the matrix and the right hand side in one cell loop, which is how a
  // matrix-based code is normally written and therefore what the matrix-free
  // rhs assembly should be timed against.
  auto assemble(const Mapping<dim> &mapping, const DoFHandler<dim> &dof_handler,
                const AffineConstraints<double> &constraints,
                const Coefficients<dim> &coefficients) -> void {
    matrix = 0.0;
    rhs = 0.0;

    const QGauss<dim> quadrature(fe_degree + 1);
    const QGauss<dim - 1> face_quadrature(fe_degree + 1);
    const FiniteElement<dim> &fe = dof_handler.get_fe();

    FEValues<dim> fe_values(mapping, fe, quadrature,
                            update_values | update_gradients |
                                update_quadrature_points | update_JxW_values);
    FEFaceValues<dim> fe_face_values(mapping, fe, face_quadrature,
                                     update_values | update_quadrature_points |
                                         update_normal_vectors |
                                         update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs(dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    const RightHandSide<dim> right_hand_side(coefficients);
    const ExactSolution<dim> exact_solution;

    for (const auto &cell : dof_handler.active_cell_iterators()) {
      if (!cell->is_locally_owned()) {
        continue;
      }

      fe_values.reinit(cell);
      cell_matrix = 0.0;
      cell_rhs = 0.0;

      for (const unsigned int q : fe_values.quadrature_point_indices()) {
        const double weight = fe_values.JxW(q);
        const double source =
            right_hand_side.value(fe_values.quadrature_point(q));

        for (unsigned int i = 0; i < dofs_per_cell; ++i) {
          const double value_i = fe_values.shape_value(i, q);
          const Tensor<1, dim> gradient_i = fe_values.shape_grad(i, q);

          for (unsigned int j = 0; j < dofs_per_cell; ++j) {
            const double value_j = fe_values.shape_value(j, q);
            const Tensor<1, dim> gradient_j = fe_values.shape_grad(j, q);

            cell_matrix(i, j) += (coefficients.mu * (gradient_i * gradient_j) +
                                  (coefficients.beta * gradient_j) * value_i +
                                  coefficients.gamma * value_j * value_i) *
                                 weight;
          }

          cell_rhs(i) += source * value_i * weight;
        }
      }

      // Neumann data mu * (grad u . n) on the one flagged face.
      for (const auto &face : cell->face_iterators()) {
        if (!face->at_boundary() ||
            face->boundary_id() != neumann_boundary_id) {
          continue;
        }

        fe_face_values.reinit(cell, face);
        for (const unsigned int q : fe_face_values.quadrature_point_indices()) {
          const double flux =
              coefficients.mu *
              (exact_solution.gradient(fe_face_values.quadrature_point(q)) *
               fe_face_values.normal_vector(q));

          for (unsigned int i = 0; i < dofs_per_cell; ++i) {
            cell_rhs(i) +=
                flux * fe_face_values.shape_value(i, q) * fe_face_values.JxW(q);
          }
        }
      }

      cell->get_dof_indices(local_dof_indices);
      constraints.distribute_local_to_global(cell_matrix, cell_rhs,
                                             local_dof_indices, matrix, rhs);
    }

    matrix.compress(VectorOperation::add);
    rhs.compress(VectorOperation::add);
  }

  auto solve(const Parameters &parameters,
             const AffineConstraints<double> &constraints) -> SolveResult {
    SolverControl solver_control(5000, parameters.tolerance * rhs.l2_norm());

    typename SolverFGMRES<VectorType>::AdditionalData fgmres_data;
    fgmres_data.max_basis_size = 50;
    SolverFGMRES<VectorType> solver(solver_control, fgmres_data);

    solution = 0.0;

    const auto solve_with = [&](const auto &preconditioner) {
      solver.solve(matrix, solution, rhs, preconditioner);
    };

    switch (parameters.preconditioner) {
    case PreconditionerType::None:
      solve_with(PreconditionIdentity());
      break;
    case PreconditionerType::Jacobi: {
      TrilinosWrappers::PreconditionJacobi jacobi;
      jacobi.initialize(matrix);
      solve_with(jacobi);
      break;
    }
    case PreconditionerType::Multigrid: {
      // Algebraic multigrid is the matrix-based counterpart of the geometric
      // V-cycle; elliptic = false picks the smoother for a non-symmetric
      // operator, as convection makes this one.
      TrilinosWrappers::PreconditionAMG amg;
      TrilinosWrappers::PreconditionAMG::AdditionalData data;
      data.elliptic = false;
      data.higher_order_elements = fe_degree > 1;
      amg.initialize(matrix, data);
      solve_with(amg);
      break;
    }
    }

    VectorType residual(rhs);
    matrix.vmult(residual, solution);
    residual -= rhs;

    constraints.distribute(solution);

    return {solver_control.last_step(), residual.l2_norm() / rhs.l2_norm()};
  }

  // Bare sparse matrix-vector products, the direct counterpart of the
  // matrix-free VMULT benchmark.
  auto benchmark(const unsigned int n_applications,
                 const MPI_Comm communicator) const -> double {
    VectorType source(rhs);
    VectorType destination(rhs);
    source = 1.0;

    matrix.vmult(destination, source); // warm the caches
    MPI_Barrier(communicator);

    Timer wall;
    for (unsigned int i = 0; i < n_applications; ++i) {
      matrix.vmult(destination, source);
    }

    return Utilities::MPI::max(wall.wall_time(), communicator) / n_applications;
  }

  // Bytes held by the matrix and its sparsity, summed over the ranks. This is
  // the number matrix-free does not pay.
  auto memory_consumption(const MPI_Comm communicator) const -> std::size_t {
    return Utilities::MPI::sum(matrix.memory_consumption(), communicator);
  }

  auto n_nonzeros() const -> types::global_dof_index {
    return matrix.n_nonzero_elements();
  }

  auto get_solution() const -> const VectorType & { return solution; }

private:
  IndexSet locally_owned;
  IndexSet locally_relevant;
  TrilinosWrappers::SparseMatrix matrix;
  VectorType solution;
  VectorType rhs;
};
