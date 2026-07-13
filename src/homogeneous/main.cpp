#include "H_Heat.hpp"

// Main function.
int
main(int argc, char *argv[])
{
  constexpr unsigned int dim = Heat::dim;

  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const auto mu = [](const Point<dim> & /*p*/) { return 1.0; };
  double a = 1.5; 
  int N = 3; 
  Point<dim> x0(0, 0, 0); 
  double sigma = 0.5; 
  const auto g  = [&]( const double  &t) {
    return (std::exp(-a * (std::cos(2*N*M_PI*t) + 1)));
  };
  const auto h = [&]( const Point<dim> &p) {
    return std::exp(-((p-x0)*(p-x0)/std::pow(sigma,2)));
  };

  const auto f  = [&](const Point<dim>  &p, const double  &t) {
    return g(t)*h(p);
  };

  Heat problem(/*ouput_filename = */ "homogeneus.msh",
               /* degree = */ 1,
               /* T = */ 1.0,
               /* theta = */ 0.5,
               /* delta_t = */ 0.05,
               mu,
               f);

  problem.run();

  return 0;
}