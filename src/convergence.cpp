#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "TestCases.hpp"
#include "TimeADR.hpp"

struct Errors
{
  double l2;
  double h1;
};

Errors
run_one(const unsigned int Nel,
        const unsigned int degree,
        const double       dt,
        const int          case_id)
{
  constexpr unsigned int dim = TimeADR::dim;

  const TestCase tc = make_test_case(case_id);

  const auto mu  = [](const Point<dim> &) { return 1.0; };
  const auto psi = [](const Point<dim> &) { return 0.0; };

  TimeADR problem(Nel,
                  degree,
                  /* T = */ 1.0,
                  /* theta = */ 0.25,
                  dt,
                  mu,
                  tc.f,
                  tc.u0,
                  tc.u1,
                  psi);

  problem.set_verbose(false);      // no per-step console spam
  problem.set_write_output(false); // no .vtu files during the sweep

  problem.run();

  ExactSolutionFunction exact(tc);
  Errors                e;
  e.l2 = problem.compute_error(VectorTools::L2_norm, exact);
  e.h1 = problem.compute_error(VectorTools::H1_seminorm, exact);
  return e;
}

double
rate(const double e_prev,
     const double e_curr,
     const double x_prev,
     const double x_curr)
{
  return std::log(e_prev / e_curr) / std::log(x_prev / x_curr);
}


// Space study
void
space_study(const unsigned int               degree,
            const std::vector<unsigned int> &nels,
            const double                     dt,
            const int                        case_id,
            const std::string               &csv_name)
{
  const bool is_root = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0;

  std::ofstream csv;
  if (is_root)
    {
      csv.open(csv_name);
      csv << "Nel,h,L2,rateL2,H1,rateH1\n";
      std::cout << "\n=== Space study, case " << case_id << ", degree "
                << degree << ", dt = " << dt << " ===\n";
      std::cout << std::setw(6) << "Nel" << std::setw(12) << "h"
                << std::setw(14) << "L2" << std::setw(8) << "rate"
                << std::setw(14) << "H1" << std::setw(8) << "rate" << "\n";
    }

  double prev_h = 0.0, prev_l2 = 0.0, prev_h1 = 0.0;

  for (std::size_t i = 0; i < nels.size(); ++i)
    {
      const double h = 1.0 / static_cast<double>(nels[i]);
      const Errors e = run_one(nels[i], degree, dt, case_id);

      double rl2 = 0.0, rh1 = 0.0;
      if (i > 0)
        {
          rl2 = rate(prev_l2, e.l2, prev_h, h);
          rh1 = rate(prev_h1, e.h1, prev_h, h);
        }

      if (is_root)
        {
          std::cout << std::setw(6) << nels[i] << std::setw(12)
                    << std::scientific << std::setprecision(3) << h
                    << std::setw(14) << e.l2;
          if (i > 0)
            std::cout << std::setw(8) << std::fixed << std::setprecision(2)
                      << rl2;
          else
            std::cout << std::setw(8) << "-";
          std::cout << std::setw(14) << std::scientific << std::setprecision(3)
                    << e.h1;
          if (i > 0)
            std::cout << std::setw(8) << std::fixed << std::setprecision(2)
                      << rh1;
          else
            std::cout << std::setw(8) << "-";
          std::cout << "\n";

          csv << nels[i] << "," << h << "," << e.l2 << ",";
          if (i > 0)
            csv << rl2;
          csv << "," << e.h1 << ",";
          if (i > 0)
            csv << rh1;
          csv << "\n";
        }

      prev_h  = h;
      prev_l2 = e.l2;
      prev_h1 = e.h1;
    }

  if (is_root)
    csv.close();
}

// Time convergence study
void
time_study(const unsigned int         Nel,
           const unsigned int         degree,
           const std::vector<double> &dts,
           const int                  case_id,
           const std::string         &csv_name)
{
  const bool is_root = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0;

  std::ofstream csv;
  if (is_root)
    {
      csv.open(csv_name);
      csv << "dt,L2,rateL2,H1,rateH1\n";
      std::cout << "\n=== Time study, case " << case_id << ", Nel = " << Nel
                << ", degree " << degree << " ===\n";
      std::cout << std::setw(12) << "dt" << std::setw(14) << "L2"
                << std::setw(8) << "rate" << std::setw(14) << "H1"
                << std::setw(8) << "rate" << "\n";
    }

  double prev_dt = 0.0, prev_l2 = 0.0, prev_h1 = 0.0;

  for (std::size_t i = 0; i < dts.size(); ++i)
    {
      const Errors e = run_one(Nel, degree, dts[i], case_id);

      double rl2 = 0.0, rh1 = 0.0;
      if (i > 0)
        {
          rl2 = rate(prev_l2, e.l2, prev_dt, dts[i]);
          rh1 = rate(prev_h1, e.h1, prev_dt, dts[i]);
        }

      if (is_root)
        {
          std::cout << std::setw(12) << std::scientific << std::setprecision(3)
                    << dts[i] << std::setw(14) << e.l2;
          if (i > 0)
            std::cout << std::setw(8) << std::fixed << std::setprecision(2)
                      << rl2;
          else
            std::cout << std::setw(8) << "-";
          std::cout << std::setw(14) << std::scientific << std::setprecision(3)
                    << e.h1;
          if (i > 0)
            std::cout << std::setw(8) << std::fixed << std::setprecision(2)
                      << rh1;
          else
            std::cout << std::setw(8) << "-";
          std::cout << "\n";

          csv << dts[i] << "," << e.l2 << ",";
          if (i > 0)
            csv << rl2;
          csv << "," << e.h1 << ",";
          if (i > 0)
            csv << rh1;
          csv << "\n";
        }

      prev_dt = dts[i];
      prev_l2 = e.l2;
      prev_h1 = e.h1;
    }

  if (is_root)
    csv.close();
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const std::vector<int> cases_to_run = {1, 2, 3, 5, 6};

  for (const int c : cases_to_run)
    {
      const std::string tag = "_case" + std::to_string(c);

      space_study(/* degree = */ 1,
                  {10, 20, 40, 80, 160},
                  /* dt = */ 1.0e-3,
                  c,
                  "conv_space_p1" + tag + ".csv");

      space_study(/* degree = */ 2,
                  {10, 20, 40, 80},
                  /* dt = */ 2.5e-4,
                  c,
                  "conv_space_p2" + tag + ".csv");

      time_study(/* Nel = */ 100,
                 /* degree = */ 2,
                 {0.02, 0.01, 0.005, 0.0025, 0.00125},
                 c,
                 "conv_time" + tag + ".csv");
    }

  return 0;
}