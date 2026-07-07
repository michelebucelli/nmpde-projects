#include "Heat.hpp"
#include <deal.II/base/convergence_table.h>
#include <fstream>
#include <cmath>

using namespace dealii;

constexpr unsigned int dim = 3;

// =====================================================
// Truncated Spectral Exact Solution (Neumann, 3D)
// =====================================================
class ExactSolution : public Function<dim>
{
public:
  ExactSolution(Heat::g_function &g_fun,
                Heat::h_function &h_fun,
                const unsigned int M = 4)
    : Function<dim>(),
      g(g_fun),
      h(h_fun),
      Mmax(M)
  {}

  virtual double value(const Point<dim> &p,
                       const unsigned int = 0) const override
  {
    const double x = p[0];
    const double y = p[1];
    const double z = p[2];
    const double t = this->get_time();

    double result = 0.0;

    for (unsigned int m = 0; m <= Mmax; ++m)
      for (unsigned int n = 0; n <= Mmax; ++n)
        for (unsigned int l = 0; l <= Mmax; ++l)
        {
          const double lambda =
            M_PI * M_PI * (m*m + n*n + l*l);

          const double h_mnl = compute_h_coeff(m,n,l);

          const double time_part =
            time_integral(lambda, h_mnl, t);

          result += time_part *
                    std::cos(m * M_PI * x) *
                    std::cos(n * M_PI * y) *
                    std::cos(l * M_PI * z);
        }

    return result;
  }

private:
  Heat::g_function &g;
  Heat::h_function &h;
  unsigned int Mmax;

  double compute_h_coeff(unsigned int m,
                         unsigned int n,
                         unsigned int l) const
  {
    const unsigned int Q = 6;
    double sum = 0.0;

    for (unsigned int i = 0; i < Q; ++i)
      for (unsigned int j = 0; j < Q; ++j)
        for (unsigned int k = 0; k < Q; ++k)
        {
          Point<dim> p((i+0.5)/Q,
                       (j+0.5)/Q,
                       (k+0.5)/Q);

          sum += h.value(p) *
                 std::cos(m*M_PI*p[0]) *
                 std::cos(n*M_PI*p[1]) *
                 std::cos(l*M_PI*p[2]);
        }

    return sum / (Q*Q*Q);
  }

  double time_integral(double lambda,
                       double h_coeff,
                       double t) const
  {
    const unsigned int Q = 20;
    double sum = 0.0;

    for (unsigned int i = 0; i < Q; ++i)
    {
      double s = t * (i+0.5)/Q;
      g.set_time(s);

      sum += std::exp(-lambda*(t-s))
             * g.value(Point<dim>())
             * h_coeff;
    }

    return sum * t / Q;
  }
};


// =====================================================
// Main
// =====================================================
int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  // mu
  const auto mu = [](const Point<dim> & /*p*/) { 
    return 1.0; 
  };

  Heat::g_function g_pulsation;
  Heat::h_function h_spatial;

  // f = g(t) * h(x)
  const auto f = [&g_pulsation, &h_spatial](const Point<dim> &p, const double &t) {
  g_pulsation.set_time(t);

  return g_pulsation.value(p) * h_spatial.value(p);
  };

  ConvergenceTable table;
  std::ofstream convergence_file("convergence.csv");
  convergence_file << "h,eL2,eH1\n";

  const std::vector<unsigned int> N_el_values =
    {5, 10, 20};
    // {5, 10, 20, 40};
    // {10};

  for (const auto &N_el : N_el_values)
  {
    const std::string mesh_file =
      "../mesh/mesh-cube-" +
      std::to_string(N_el) + ".msh";

    Heat problem(/*mesh_filename = */ mesh_file,
               /* degree = */ 2,
               /* T = */ 1.0,
               /* theta = */ 1.0,    
               /* delta_t = */ 0.001,
               mu,    
               f);    

    problem.run();

    const double h = 1.0 / N_el;

    ExactSolution exact(g_pulsation,
                        h_spatial,
                        4);

    exact.set_time(1.0);

    const double error_L2 =
      problem.compute_error(VectorTools::L2_norm,
                            exact);

    const double error_H1 =
      problem.compute_error(VectorTools::H1_norm,
                            exact);

    table.add_value("h", h);
    table.add_value("L2", error_L2);
    table.add_value("H1", error_H1);

    convergence_file << h << ","
                     << error_L2 << ","
                     << error_H1 << "\n";
  }

  table.evaluate_all_convergence_rates(
    ConvergenceTable::reduction_rate_log2);

  table.write_text(std::cout);

  return 0;
}
