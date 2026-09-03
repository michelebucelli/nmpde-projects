// Unit tests for the pointwise data in cdr_problem.hpp: the coefficient
// struct, the manufactured solution and its derivatives, the right-hand side
// derived from them, and the lanewise evaluation helper.

#include "cdr_problem.hpp"
#include "test_harness.hpp"

#include <deal.II/base/mpi.h>

using namespace dealii;

namespace {

// Central difference of a scalar function of a point, used to check the
// analytic gradient and laplacian against the value function alone.
template <int dim, typename Callable>
auto finite_difference(const Callable &f, const Point<dim> &p,
                       const unsigned int d, const double h) -> double {
  Point<dim> forward = p, backward = p;
  forward[d] += h;
  backward[d] -= h;
  return (f(forward) - f(backward)) / (2.0 * h);
}

template <int dim, typename Callable>
auto second_finite_difference(const Callable &f, const Point<dim> &p,
                              const unsigned int d, const double h) -> double {
  Point<dim> forward = p, backward = p;
  forward[d] += h;
  backward[d] -= h;
  return (f(forward) - 2.0 * f(p) + f(backward)) / (h * h);
}

template <int dim> auto sample_point() -> Point<dim> {
  Point<dim> p;
  for (unsigned int d = 0; d < dim; ++d) {
    p[d] = 0.17 + 0.23 * d;
  }
  return p;
}

auto test_coefficients(TestSuite &suite) -> void {
  const Coefficients<3> defaults;
  suite.close(defaults.mu, 1.0, 1e-15, "Coefficients default mu");
  suite.close(defaults.gamma, 1.0, 1e-15, "Coefficients default gamma");
  for (unsigned int d = 0; d < 3; ++d) {
    suite.close(defaults.beta[d], 1.0 + 0.5 * d, 1e-15,
                "Coefficients default beta[" + std::to_string(d) + "]");
  }

  const Coefficients<2> custom(2.5, -0.75);
  suite.close(custom.mu, 2.5, 1e-15, "Coefficients mu from constructor");
  suite.close(custom.gamma, -0.75, 1e-15,
              "Coefficients gamma from constructor");

  // The multigrid levels run in single precision, so a double-valued
  // Coefficients has to convert.
  const Coefficients<2, float> converted(custom);
  suite.close(converted.mu, 2.5, 1e-6, "converted mu");
  suite.close(converted.gamma, -0.75, 1e-6, "converted gamma");
  for (unsigned int d = 0; d < 2; ++d) {
    suite.close(converted.beta[d], custom.beta[d], 1e-6,
                "converted beta[" + std::to_string(d) + "]");
  }
  suite.check(std::is_same_v<decltype(Coefficients<2, float>::beta),
                             Tensor<1, 2, float>>,
              "converted beta is stored in the target number type");
  suite.check(std::is_same_v<decltype(Coefficients<2, float>::mu), float>,
              "converted mu is stored in the target number type");
  suite.check(std::is_same_v<decltype(Coefficients<2, float>::gamma), float>,
              "converted gamma is stored in the target number type");
}

template <int dim> auto test_exact_solution(TestSuite &suite) -> void {
  const std::string tag = std::to_string(dim) + "D ";
  const ExactSolution<dim> exact;
  const double k = ExactSolution<dim>::k;

  // u = prod_d cos(k x_d) is one at the origin and vanishes on x_d = 1.
  suite.close(exact.value(Point<dim>()), 1.0, 1e-14, tag + "u at the origin");

  Point<dim> on_far_face;
  on_far_face[0] = 1.0;
  suite.close(exact.value(on_far_face), 0.0, 1e-14,
              tag + "u vanishes on x_0 = 1");

  const Point<dim> p = sample_point<dim>();
  double product = 1.0;
  for (unsigned int d = 0; d < dim; ++d) {
    product *= std::cos(k * p[d]);
  }
  suite.close(exact.value(p), product, 1e-14,
              tag + "u is a product of cosines");

  const auto value = [&exact](const Point<dim> &x) { return exact.value(x); };
  const Tensor<1, dim> gradient = exact.gradient(p);
  double laplacian_by_differences = 0.0;
  for (unsigned int d = 0; d < dim; ++d) {
    suite.close(gradient[d], finite_difference(value, p, d, 1e-5), 1e-8,
                tag + "grad u[" + std::to_string(d) + "] matches differences");
    laplacian_by_differences += second_finite_difference(value, p, d, 1e-4);
  }

  suite.close(exact.laplacian(p), laplacian_by_differences, 1e-6,
              tag + "laplacian matches differences");
  // u is an eigenfunction of the Laplacian with eigenvalue -dim k^2.
  suite.close(exact.laplacian(p), -dim * k * k * exact.value(p), 1e-14,
              tag + "laplacian is the eigenvalue relation");
}

template <int dim> auto test_right_hand_side(TestSuite &suite) -> void {
  const std::string tag = std::to_string(dim) + "D ";
  const Coefficients<dim> coefficients(0.7, 1.3);
  const RightHandSide<dim> right_hand_side(coefficients);
  const ExactSolution<dim> exact;

  // f has to be exactly the strong form applied to the manufactured solution.
  for (const double shift : {0.0, 0.31, 0.62}) {
    Point<dim> p = sample_point<dim>();
    for (unsigned int d = 0; d < dim; ++d) {
      p[d] = std::fmod(p[d] + shift, 1.0);
    }

    const double expected = -coefficients.mu * exact.laplacian(p) +
                            coefficients.beta * exact.gradient(p) +
                            coefficients.gamma * exact.value(p);
    suite.close(right_hand_side.value(p), expected, 1e-13,
                tag + "f is the strong form at shift " + std::to_string(shift));
  }

  // With zero convection and reaction the source reduces to -mu delta u.
  Coefficients<dim> diffusion_only(2.0, 0.0);
  diffusion_only.beta = Tensor<1, dim>();
  const RightHandSide<dim> pure_diffusion(diffusion_only);
  const Point<dim> p = sample_point<dim>();
  suite.close(pure_diffusion.value(p), -2.0 * exact.laplacian(p), 1e-13,
              tag + "f reduces to -mu delta u");
}

template <int dim> auto test_evaluate_lanewise(TestSuite &suite) -> void {
  const std::string tag = std::to_string(dim) + "D ";
  constexpr unsigned int n_lanes = VectorizedArray<double>::size();

  // Give every lane a different point so that a broadcast of lane 0 would show
  // up as a mismatch.
  Point<dim, VectorizedArray<double>> vectorized;
  for (unsigned int lane = 0; lane < n_lanes; ++lane) {
    for (unsigned int d = 0; d < dim; ++d) {
      vectorized[d][lane] = 0.1 * (d + 1) + 0.05 * lane;
    }
  }

  const ExactSolution<dim> exact;
  const auto scalar = [&exact](const Point<dim> &p) { return exact.value(p); };
  const VectorizedArray<double> result =
      evaluate_lanewise<dim, double>(scalar, vectorized);

  for (unsigned int lane = 0; lane < n_lanes; ++lane) {
    Point<dim> p;
    for (unsigned int d = 0; d < dim; ++d) {
      p[d] = vectorized[d][lane];
    }
    suite.close(result[lane], scalar(p), 1e-14,
                tag + "lane " + std::to_string(lane) +
                    " is evaluated on its "
                    "own point");
  }
}

} // namespace

int main(int argc, char *argv[]) {
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  TestSuite suite("problem_data");

  test_coefficients(suite);
  test_exact_solution<2>(suite);
  test_exact_solution<3>(suite);
  test_right_hand_side<2>(suite);
  test_right_hand_side<3>(suite);
  test_evaluate_lanewise<2>(suite);
  test_evaluate_lanewise<3>(suite);

  return suite.report();
}
