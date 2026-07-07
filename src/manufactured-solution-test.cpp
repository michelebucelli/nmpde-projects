#include "Heat.hpp"
#include <deal.II/base/convergence_table.h>
#include <iomanip>
#include <sstream>

constexpr unsigned int dim = Heat::dim;

// This test uses the method of manufactured solutions. We define an exact
// solution that satisfies the boundary and initial conditions, and then
// compute the forcing term `f` that makes it a solution to the PDE.

// Manufactured solution: u(x, y, z, t) = t^2 * cos(pi*x) * cos(pi*y) * cos(pi*z)
// This solution satisfies:
// - Homogeneous Neumann BCs on the unit cube [0,1]^3, because its normal
//   derivative is zero on all faces.
// - Zero initial condition u(x, y, z, 0) = 0.
class ManufacturedSolution : public Function<dim>
{
public:
  virtual double
  value(const Point<dim> &p, const unsigned int = 0) const override
  {
    const double t = get_time();
    return t * t * std::cos(M_PI * p[0]) * std::cos(M_PI * p[1]) *
           std::cos(M_PI * p[2]);
  }
};

// Forcing term f = u_t - div(mu * grad(u))
// For mu=1, f = u_t - laplacian(u)
class ForcingTerm : public Function<dim>
{
public:
  virtual double
  value(const Point<dim> &p, const unsigned int = 0) const override
  {
    const double t = get_time();
    const double c_x = std::cos(M_PI * p[0]);
    const double c_y = std::cos(M_PI * p[1]);
    const double c_z = std::cos(M_PI * p[2]);

    // Time derivative part: du/dt = 2t * cos(pi*x) * cos(pi*y) * cos(pi*z)
    const double time_derivative = 2.0 * t * c_x * c_y * c_z;

    // Laplacian part: laplacian(u) = -3 * pi^2 * u
    const double laplacian =
      -3.0 * M_PI * M_PI * (t * t * c_x * c_y * c_z);

    return time_derivative - laplacian;
  }
};

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  ForcingTerm forcing_term;
  const auto  f = & {
    forcing_term.set_time(t);
    return forcing_term.value(p);
  };

  const auto mu =  { return 1.0; };

  Heat problem(/*mesh_filename = */ "../mesh/mesh-cube-10.msh",
               /* degree = */ 1,
               /* T = */ 1.0,
               /* theta = */ 0.5, // Crank-Nicolson
               /* delta_t = */ 0.01,
               mu,
               f);

  problem.run();

  ManufacturedSolution exact_solution;
  exact_solution.set_time(1.0); // Evaluate error at final time

  const double error_L2 =
    problem.compute_error(VectorTools::L2_norm, exact_solution);

  std::cout << "L2 error at T=1.0: " << error_L2 << std::endl;

  return 0;
}