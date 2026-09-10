#ifndef TESTCASES_HPP
#define TESTCASES_HPP

#include <deal.II/base/function.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <cmath>
#include <functional>
#include <stdexcept>
#include <string>

using namespace dealii;

struct TestCase
{
  static constexpr unsigned int dim = 2;

  std::string name;

  // u_exact and grad_exact both depend on space and time
  std::function<double(const Point<dim> &, double)>         u_exact;
  std::function<Tensor<1, dim>(const Point<dim> &, double)> grad_exact;

  std::function<double(const Point<dim> &, double)> f;
  std::function<double(const Point<dim> &)>         u0;
  std::function<double(const Point<dim> &)>         u1;
};

class ExactSolutionFunction : public Function<TestCase::dim>
{
public:
  static constexpr unsigned int dim = TestCase::dim;

  explicit ExactSolutionFunction(const TestCase &tc)
    : u_exact(tc.u_exact)
    , grad_exact(tc.grad_exact)
  {}

  double
  value(const Point<dim> &p,
        const unsigned int /*component*/ = 0) const override
  {
    return u_exact(p, this->get_time());
  }

  Tensor<1, dim>
  gradient(const Point<dim> &p,
           const unsigned int /*component*/ = 0) const override
  {
    return grad_exact(p, this->get_time());
  }

private:
  std::function<double(const Point<dim> &, double)>         u_exact;
  std::function<Tensor<1, dim>(const Point<dim> &, double)> grad_exact;
};

inline TestCase
make_test_case(const int id)
{
  constexpr unsigned int dim = TestCase::dim;
  const double           pi  = numbers::PI;

  TestCase tc;

  switch (id)
    {
        // =====================================================================
        // 1) Proper mode of the square (free vibration, f = 0).
        //    u = sin(pi x) sin(pi y) cos(sqrt(2) pi t)
        // =====================================================================
        case 1: {
          tc.name            = "1 - proper mode of the square";
          const double omega = std::sqrt(2.0) * pi;

          tc.u_exact = [=](const Point<dim> &p, double t) {
            return std::sin(pi * p[0]) * std::sin(pi * p[1]) *
                   std::cos(omega * t);
          };
          tc.grad_exact = [=](const Point<dim> &p, double t) {
            Tensor<1, dim> g;
            g[0] = pi * std::cos(pi * p[0]) * std::sin(pi * p[1]) *
                   std::cos(omega * t);
            g[1] = pi * std::sin(pi * p[0]) * std::cos(pi * p[1]) *
                   std::cos(omega * t);
            return g;
          };
          tc.f  = [](const Point<dim> &, double) { return 0.0; };
          tc.u0 = [=](const Point<dim> &p) {
            return std::sin(pi * p[0]) * std::sin(pi * p[1]);
          };
          tc.u1 = [](const Point<dim> &) { return 0.0; };
          break;
        }

        // =====================================================================
        // 2) Manufactured solution with forcing (decaying in time).
        //    u = sin(pi x) sin(pi y) exp(-t)
        //    f = (1 + 2 pi^2) sin(pi x) sin(pi y) exp(-t)
        //    u1 = u_t(.,0) = -sin(pi x) sin(pi y) (NON-zero initial velocity)
        // =====================================================================
        case 2: {
          tc.name = "2 - manufactured solution with forcing";

          tc.u_exact = [=](const Point<dim> &p, double t) {
            return std::sin(pi * p[0]) * std::sin(pi * p[1]) * std::exp(-t);
          };
          tc.grad_exact = [=](const Point<dim> &p, double t) {
            Tensor<1, dim> g;
            g[0] =
              pi * std::cos(pi * p[0]) * std::sin(pi * p[1]) * std::exp(-t);
            g[1] =
              pi * std::sin(pi * p[0]) * std::cos(pi * p[1]) * std::exp(-t);
            return g;
          };
          tc.f = [=](const Point<dim> &p, double t) {
            return (1.0 + 2.0 * pi * pi) * std::sin(pi * p[0]) *
                   std::sin(pi * p[1]) * std::exp(-t);
          };
          tc.u0 = [=](const Point<dim> &p) {
            return std::sin(pi * p[0]) * std::sin(pi * p[1]);
          };
          tc.u1 = [=](const Point<dim> &p) {
            return -std::sin(pi * p[0]) * std::sin(pi * p[1]);
          };
          break;
        }

        // =====================================================================
        // 3) Polynomial solution.
        //    u = t^2 x(1-x) y(1-y)
        //    f = 2 x(1-x) y(1-y) + 2 t^2 [ x(1-x) + y(1-y) ]
        //    u0 = 0, u1 = 0.
        // =====================================================================
        case 3: {
          tc.name = "3 - polynomial solution";

          tc.u_exact = [=](const Point<dim> &p, double t) {
            const double x = p[0], y = p[1];
            return t * t * x * (1.0 - x) * y * (1.0 - y);
          };
          tc.grad_exact = [=](const Point<dim> &p, double t) {
            const double   x = p[0], y = p[1];
            Tensor<1, dim> g;
            g[0] = t * t * (1.0 - 2.0 * x) * y * (1.0 - y);
            g[1] = t * t * x * (1.0 - x) * (1.0 - 2.0 * y);
            return g;
          };
          tc.f = [=](const Point<dim> &p, double t) {
            const double x = p[0], y = p[1];
            return 2.0 * x * (1.0 - x) * y * (1.0 - y) +
                   2.0 * t * t * (x * (1.0 - x) + y * (1.0 - y));
          };
          tc.u0 = [](const Point<dim> &) { return 0.0; };
          tc.u1 = [](const Point<dim> &) { return 0.0; };
          break;
        }

        // =====================================================================
        // 5) Multi-frequency mode (free vibration, f = 0).
        //    u = sin(2 pi x) sin(3 pi y) cos(sqrt(13) pi t)
        // =====================================================================
        case 5: {
          tc.name            = "5 - multi-frequency mode";
          const double omega = std::sqrt(13.0) * pi;

          tc.u_exact = [=](const Point<dim> &p, double t) {
            return std::sin(2.0 * pi * p[0]) * std::sin(3.0 * pi * p[1]) *
                   std::cos(omega * t);
          };
          tc.grad_exact = [=](const Point<dim> &p, double t) {
            Tensor<1, dim> g;
            g[0] = 2.0 * pi * std::cos(2.0 * pi * p[0]) *
                   std::sin(3.0 * pi * p[1]) * std::cos(omega * t);
            g[1] = 3.0 * pi * std::sin(2.0 * pi * p[0]) *
                   std::cos(3.0 * pi * p[1]) * std::cos(omega * t);
            return g;
          };
          tc.f  = [](const Point<dim> &, double) { return 0.0; };
          tc.u0 = [=](const Point<dim> &p) {
            return std::sin(2.0 * pi * p[0]) * std::sin(3.0 * pi * p[1]);
          };
          tc.u1 = [](const Point<dim> &) { return 0.0; };
          break;
        }

        // =====================================================================
        // 6) Sum of two modes (free vibration, f = 0).
        //    u = sin(pi x) sin(pi y) cos(sqrt(2) pi t)
        //      + 1/2 sin(4 pi x) sin(2 pi y) cos(sqrt(20) pi t)
        // =====================================================================
        case 6: {
          tc.name             = "6 - sum of two modes";
          const double omega1 = std::sqrt(2.0) * pi;
          const double omega2 = std::sqrt(20.0) * pi;

          tc.u_exact = [=](const Point<dim> &p, double t) {
            return std::sin(pi * p[0]) * std::sin(pi * p[1]) *
                     std::cos(omega1 * t) +
                   0.5 * std::sin(4.0 * pi * p[0]) * std::sin(2.0 * pi * p[1]) *
                     std::cos(omega2 * t);
          };
          tc.grad_exact = [=](const Point<dim> &p, double t) {
            Tensor<1, dim> g;
            g[0] = pi * std::cos(pi * p[0]) * std::sin(pi * p[1]) *
                     std::cos(omega1 * t) +
                   2.0 * pi * std::cos(4.0 * pi * p[0]) *
                     std::sin(2.0 * pi * p[1]) * std::cos(omega2 * t);
            g[1] = pi * std::sin(pi * p[0]) * std::cos(pi * p[1]) *
                     std::cos(omega1 * t) +
                   pi * std::sin(4.0 * pi * p[0]) * std::cos(2.0 * pi * p[1]) *
                     std::cos(omega2 * t);
            return g;
          };
          tc.f  = [](const Point<dim> &, double) { return 0.0; };
          tc.u0 = [=](const Point<dim> &p) {
            return std::sin(pi * p[0]) * std::sin(pi * p[1]) +
                   0.5 * std::sin(4.0 * pi * p[0]) * std::sin(2.0 * pi * p[1]);
          };
          tc.u1 = [](const Point<dim> &) { return 0.0; };
          break;
        }

        // =====================================================================
        // 4) Wave with NON-homogeneous, time-dependent Dirichlet boundary.
        //    u = cos(t) (x + y)
        //    f = u_tt - Delta u = -cos(t)(x+y) - 0 = -cos(t)(x+y)
        //    On the boundary u != 0, so g(x,y,t) = cos(t)(x+y) is imposed via
        //    the exact solution itself (handled by the solver's lifting).
        //    u0 = u(.,0) = x + y ;   u1 = u_t(.,0) = 0.
        // =====================================================================
        case 4: {
          tc.name = "4 - non-homogeneous, time-dependent Dirichlet boundary";

          tc.u_exact = [=](const Point<dim> &p, double t) {
            return std::cos(t) * (p[0] + p[1]);
          };
          tc.grad_exact = [=](const Point<dim> &p, double t) {
            (void)p;
            Tensor<1, dim> g;
            g[0] = std::cos(t);
            g[1] = std::cos(t);
            return g;
          };
          tc.f = [=](const Point<dim> &p, double t) {
            return -std::cos(t) * (p[0] + p[1]);
          };
          tc.u0 = [](const Point<dim> &p) { return p[0] + p[1]; };
          tc.u1 = [](const Point<dim> &) { return 0.0; };
          break;
        }

      default:
        throw std::runtime_error(
          "Unknown test case id (use 1, 2, 3, 4, 5 or 6).");
    }

  return tc;
}

#endif
