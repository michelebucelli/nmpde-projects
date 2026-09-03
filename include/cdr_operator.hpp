#pragma once

#include "cdr_problem.hpp"

#include <deal.II/base/vectorization.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/tools.h>

using namespace dealii;

// Cell-loop parallelism is hybrid: MPI owns mesh partitions, TBB splits the
// local cell batches, and VectorizedArray fills SIMD lanes inside a batch.
// partition_partition is deal.II's default nested TBB scheme.
template <int dim, typename Number>
auto matrix_free_data(const UpdateFlags cell_flags) ->
    typename MatrixFree<dim, Number>::AdditionalData {
  typename MatrixFree<dim, Number>::AdditionalData data;
  data.tasks_parallel_scheme =
      MatrixFree<dim, Number>::AdditionalData::partition_partition;
  data.mapping_update_flags = cell_flags;
  return data;
}

// Matrix-free convention-diffussion-reaction operator
//
// a(u, v) = (mu grad u, grad v) + (beta . grad u, v) + (gamma u, v).
//
template <int dim, int fe_degree, typename Number>
class Operator : public MatrixFreeOperators::Base<
                     dim, LinearAlgebra::distributed::Vector<Number>> {
public:
  using value_type = Number;
  using VectorType = LinearAlgebra::distributed::Vector<Number>;
  using Base = MatrixFreeOperators::Base<dim, VectorType>;
  using FEEval = FEEvaluation<dim, fe_degree, fe_degree + 1, 1, Number>;

  auto set_coefficients(const Coefficients<dim, Number> &c) -> void {
    coefficients = c;
  }

  // Core quadrature kernel evaluating local cell integrals for the CDR
  // operator.
  template <typename Source>
  auto do_quadrature_points(FEEval &phi, const Number sign,
                            const Source &source) const -> void {
    for (unsigned int q = 0; q < phi.n_q_points; ++q) {
      const auto value = phi.get_value(q);
      const auto gradient = phi.get_gradient(q);

      VectorizedArray<Number> convection = VectorizedArray<Number>();
      for (unsigned int d = 0; d < dim; ++d) {
        convection += coefficients.beta[d] * gradient[d];
      }

      phi.submit_gradient((sign * coefficients.mu) * gradient, q);
      phi.submit_value(
          sign * (convection + coefficients.gamma * value) + source(phi, q), q);
    }
  }

  // Convenience overload for standard operator evaluation (A*u) without a
  // source term.
  auto do_quadrature_points(FEEval &phi, const Number sign = 1.0) const
      -> void {
    do_quadrature_points(phi, sign, [](const FEEval &, const unsigned int) {
      return VectorizedArray<Number>();
    });
  }

  // Computes and stores the inverse diagonal entries of the operator.
  //
  // Assembles A_ii for each locally ownd DoF by evaluating the operator agains
  // unit basis functions, then stores 1 / A_ii inside inverse_diagonal_entries.
  auto compute_diagonal() -> void override {
    this->inverse_diagonal_entries =
        std::make_shared<DiagonalMatrix<VectorType>>();
    VectorType &inverse_diagonal = this->inverse_diagonal_entries->get_vector();
    this->data->initialize_dof_vector(inverse_diagonal);

    MatrixFreeTools::compute_diagonal<dim, fe_degree, fe_degree + 1, 1, Number,
                                      VectorizedArray<Number>, VectorType>(
        *this->data, inverse_diagonal, [this](FEEval &phi) {
          phi.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
          do_quadrature_points(phi);
          phi.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
        });

    for (unsigned int i = 0; i < inverse_diagonal.locally_owned_size(); ++i) {
      const Number d = inverse_diagonal.local_element(i);
      inverse_diagonal.local_element(i) =
          (std::abs(d) > 1e-10 ? Number(1.0) / d : Number(1.0));
    }
  }

private:
  // Applies the matrix-vector operation (destination += A * src)
  // Delegates cell-loop scheduling to MatrixFree (MPI + TBB + SIMD).
  auto apply_add(VectorType &destination, const VectorType &source) const
      -> void override {
    this->data->cell_loop(&Operator::local_apply, this, destination, source);
  }

  // Thread/SIMD worker kernel executing matrix-vector multiplication on a range
  // of cells.
  //
  // Performs gather -> evaluate -> PDE quadrature kernel -> integrate ->
  // scatter operations sequentially over the assigned cell batch.
  auto
  local_apply(const MatrixFree<dim, Number> &data, VectorType &destination,
              const VectorType &source,
              const std::pair<unsigned int, unsigned int> &cell_range) const
      -> void {
    FEEval phi(data);

    for (unsigned int cell = cell_range.first; cell < cell_range.second;
         ++cell) {
      phi.reinit(cell);
      phi.gather_evaluate(source,
                          EvaluationFlags::values | EvaluationFlags::gradients);
      do_quadrature_points(phi);
      phi.integrate_scatter(
          EvaluationFlags::values | EvaluationFlags::gradients, destination);
    }
  }

  Coefficients<dim, Number> coefficients;
};
