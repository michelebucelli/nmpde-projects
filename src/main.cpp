#include <iomanip>
#include <iostream>

#include "TestCases.hpp"
#include "TimeADR.hpp"

int
main(int argc, char *argv[])
{
  constexpr unsigned int dim = TimeADR::dim;

  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const int      chosen_case = 5;
  const TestCase tc          = make_test_case(chosen_case);

  const unsigned int Nel     = 80;
  const unsigned int degree  = 1;
  const double       T       = 1.0;
  const double       theta   = 0.25;
  const double       delta_t = 0.005;

  const auto mu = [](const Point<dim> & /*p*/) { return 1.0; };

  const auto psi = [](const Point<dim> & /*p*/) { return 0.0; }; // non usato

  TimeADR problem(Nel,
                  degree,
                  T,
                  theta,
                  delta_t,
                  mu,
                  tc.f,
                  tc.u0,
                  tc.u1,
                  psi); // prendiamo i dati da testcase.hpp

  if (chosen_case == 4)
    problem.set_dirichlet_datum(tc.u_exact);

  const bool is_root = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0;
  if (is_root)
    std::cout << "Running test case: " << tc.name << "\n";

  problem.run();

  ExactSolutionFunction exact(tc);

  const double error_L2 = problem.compute_error(VectorTools::L2_norm, exact);
  const double error_H1 =
    problem.compute_error(VectorTools::H1_seminorm, exact);
  const double error_Linfty =
    problem.compute_error(VectorTools::Linfty_norm, exact);

  if (is_root)
    {
      std::cout << std::scientific << std::setprecision(4);
      std::cout << "L2 error       = " << error_L2 << std::endl;
      std::cout << "H1 seminorm    = " << error_H1 << std::endl;
      std::cout << "Linfty error   = " << error_Linfty << std::endl;
    }

  return 0;
}
