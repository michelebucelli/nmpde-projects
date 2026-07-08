#include "Heat.hpp"

#include <cmath>

constexpr double mu_value = 0.05;

// This test uses the method of manufactured solutions. We define an exact
// solution that satisfies the boundary and initial conditions, and then
// compute the forcing term `f` that makes it a solution to the PDE.

// Manufactured solution: u(x, y, z, t) = sin(t) * cos(pi*x) * cos(pi*y) * cos(pi*z)
// This solution satisfies:
// - Homogeneous Neumann BCs on the unit cube [0,1]^3, because its normal
//   derivative is zero on all faces.
// - Zero initial condition u(x, y, z, 0) = 0.
class ManufacturedSolution : public Function<Heat::dim>
{
public:
  virtual double
  value(const Point<Heat::dim> &p, const unsigned int = 0) const override
  {
    const double t = get_time();
    return std::sin(t) * std::cos(numbers::PI * p[0]) * std::cos(numbers::PI * p[1]) *
           std::cos(numbers::PI * p[2]);
  }
};

// Forcing term f = u_t - div(mu * grad(u))
class ForcingTerm : public Function<Heat::dim>
{
public:
  virtual double
  value(const Point<Heat::dim> &p, const unsigned int = 0) const override
  {
    const double t   = get_time();
    const double c_x = std::cos(numbers::PI * p[0]);
    const double c_y = std::cos(numbers::PI * p[1]);
    const double c_z = std::cos(numbers::PI * p[2]);

    // Time derivative part: du/dt = cos(t) * cos(pi*x) * cos(pi*y) * cos(pi*z)
    const double time_derivative = std::cos(t) * c_x * c_y * c_z;

    // Laplacian part: laplacian(u) = -3 * pi^2 * u
    const double laplacian = -3.0 * numbers::PI * numbers::PI * (std::sin(t) * c_x * c_y * c_z);

    return time_derivative - mu_value * laplacian;
  }
};

namespace
{
  // Runs the manufactured-solution problem once, with the given adaptivity
  // setting, and returns the L2 error at T=1. Uses looser time-adaptivity
  // tolerances than the main pulsed experiment: this forcing is smooth (no
  // sharp feature to miss), so the tight pulse-detection tolerances would
  // only cost extra rejections/runtime here for no accuracy benefit.
  double
  run_case(bool adaptive, ConditionalOStream &pcout)
  {
    ForcingTerm forcing_term;
    const auto  f = [&forcing_term](const Point<Heat::dim> &p, const double &t) {
      forcing_term.set_time(t);
      return forcing_term.value(p);
    };

    const auto mu = [](const Point<Heat::dim> &) { return mu_value; };

    Heat problem(4,
                 /* degree = */ 1,
                 /* T = */ 1.0,
                 /* theta = */ 0.5, // Crank-Nicolson
                 /* delta_t = */ 0.01,
                 mu,
                 f,
                 adaptive,
                 /* tol_time_max = */ 5e-3,
                 /* tol_time_min = */ 1e-4,
                 /* dt_max = */ 0.1);

    problem.run();

    ManufacturedSolution exact_solution;
    exact_solution.set_time(1.0); // Evaluate error at final time

    const double error_L2 = problem.compute_error(VectorTools::L2_norm, exact_solution);

    pcout << (adaptive ? "adaptive" : "non-adaptive") << ": L2 error at T=1.0: " << error_L2
          << std::endl;

    // Write u_h and u_exact side by side, for a real-vs-constructed comparison.
    const IndexSet locally_owned_dofs = problem.get_dof_handler().locally_owned_dofs();
    const IndexSet locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(problem.get_dof_handler());

    TrilinosWrappers::MPI::Vector exact_owned(locally_owned_dofs, MPI_COMM_WORLD);
    VectorTools::interpolate(problem.get_dof_handler(), exact_solution, exact_owned);

    TrilinosWrappers::MPI::Vector exact_ghosted(locally_owned_dofs,
                                                 locally_relevant_dofs,
                                                 MPI_COMM_WORLD);
    exact_ghosted = exact_owned;

    DataOut<Heat::dim> data_out;
    data_out.add_data_vector(problem.get_dof_handler(), problem.get_solution(), "u_h");
    data_out.add_data_vector(problem.get_dof_handler(), exact_ghosted, "u_exact");
    data_out.build_patches();
    data_out.write_vtu_with_pvtu_record(/* folder = */ "./",
                                        /* basename = */ std::string("manufactured-solution-") +
                                          (adaptive ? "adaptive" : "nonadaptive"),
                                        /* index = */ 0,
                                        MPI_COMM_WORLD);

    return error_L2;
  }
} // namespace

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
  ConditionalOStream pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0);

  const double tolerance = 1.0e-1;
  bool         passed    = true;

  for (const bool adaptive : {true, false})
    {
      const double error_L2 = run_case(adaptive, pcout);
      if (error_L2 > tolerance)
        {
          if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
            std::cerr << "Manufactured-solution test FAILED ("
                      << (adaptive ? "adaptive" : "non-adaptive") << "): L2 error " << error_L2
                      << " exceeds tolerance " << tolerance << std::endl;
          passed = false;
        }
    }

  if (!passed)
    return 1;

  pcout << "Manufactured-solution test PASSED (adaptive and non-adaptive)" << std::endl;
  return 0;
}
