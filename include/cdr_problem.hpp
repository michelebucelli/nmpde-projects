#pragma once

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/multigrid.h>

using namespace dealii;

constexpr types::boundary_id neumann_boundary_id = 1;

// Coefficients of	-div(mu grad u) + beta . grad u + gamma u = f.
template <int dim, typename Number = double> struct Coefficients {
  Coefficients(const Number mu = 1.0, const Number gamma = 1.0)
      : mu(mu), gamma(gamma) {
    for (unsigned int d = 0; d < dim; ++d) {
      beta[d] = static_cast<Number>(1.0 + 0.5 * d);
    }
  }

  Coefficients(const Number mu, const Number gamma,
               const Tensor<1, dim, Number> &beta)
      : mu(mu), gamma(gamma), beta(beta) {}

  template <typename OtherNumber>
  Coefficients(const Coefficients<dim, OtherNumber> &other)
      : mu(static_cast<Number>(other.mu)),
        gamma(static_cast<Number>(other.gamma)) {
    for (unsigned int d = 0; d < dim; ++d) {
      beta[d] = static_cast<Number>(other.beta[d]);
    }
  }

  Number mu;
  Number gamma;
  Tensor<1, dim, Number> beta;
};

// Manufactured solution u(x) = prod_d cos(pi x_d / 2) on (0,1)^d. It is non
// zero on the faces x_d = 0, so the Dirichlet data are inhomogeneous, and its
// normal derivative on x_0 = 1 is non zero as well.
template <int dim> class ExactSolution : public Function<dim> {
public:
  static constexpr double k = 0.5 * numbers::PI;

  auto value(const Point<dim> &p, const unsigned int = 0) const
      -> double override {
    double v = 1.0;
    for (unsigned int d = 0; d < dim; ++d) {
      v *= std::cos(k * p[d]);
    }
    return v;
  }

  auto gradient(const Point<dim> &p, const unsigned int = 0) const
      -> Tensor<1, dim> override {
    Tensor<1, dim> g;
    for (unsigned int i = 0; i < dim; ++i) {
      double v = -k * std::sin(k * p[i]);
      for (unsigned int j = 0; j < dim; ++j) {
        if (j != i) {
          v *= std::cos(k * p[j]);
        }
        g[i] = v;
      }
    }
    return g;
  }

  // u is an eigenfunction of the Laplacian: -delta u = dim k^2 u.
  auto laplacian(const Point<dim> &p, const unsigned int = 0) const
      -> double override {
    return -dim * k * k * value(p);
  }
};

template <int dim> class RightHandSide : public Function<dim> {
public:
  RightHandSide(const Coefficients<dim> &coefficients)
      : coefficients(coefficients) {}

  auto value(const Point<dim> &p, const unsigned int = 0) const
      -> double override {
    return -coefficients.mu * exact.laplacian(p) +
           coefficients.beta * exact.gradient(p) +
           coefficients.gamma * exact.value(p);
  }

private:
  const Coefficients<dim> coefficients;
  const ExactSolution<dim> exact;
};

template <int dim, typename Number, typename Callable>
auto evaluate_lanewise(const Callable &f,
                       const Point<dim, VectorizedArray<Number>> &p)
    -> VectorizedArray<Number> {
  VectorizedArray<Number> result;
  for (unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v) {
    Point<dim> q;
    for (unsigned int d = 0; d < dim; ++d) {
      q[d] = p[d][v];
    }
    result[v] = f(q);
  }
  return result;
}
