#include "HeterogeneousDiffusion.hpp"

#include <string>
#include <vector>
#include <deal.II/base/convergence_table.h>

int
main(int argc, char **argv)
{
  using namespace dealii;

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  const MPI_Comm mpi_communicator = MPI_COMM_WORLD;

  if (argc > 1 && std::string(argv[1]) == "verify")
    {
      const unsigned int ref_min = argc > 2 ? std::stoi(argv[2]) : 3;
      const unsigned int ref_max = argc > 3 ? std::stoi(argv[3]) : 5;

      ConditionalOStream pcout(std::cout,
                               Utilities::MPI::this_mpi_process(
                                 mpi_communicator) == 0);
      ConvergenceTable table;

      for (unsigned int r = ref_min; r <= ref_max; ++r)
        {
          HeterogeneousDiffusion problem(r,
                                         /* p (ignored) = */ 0.0,
                                         mpi_communicator,
                                         /* verify_mode = */ true);
          problem.setup_and_assemble();
          problem.solve("amg");   // fastest preconditioner

          const double e_L2 = problem.compute_error(VectorTools::L2_norm);
          const double e_H1 = problem.compute_error(VectorTools::H1_norm);

          table.add_value("h",  1.0 / std::pow(2.0, r));
          table.add_value("L2", e_L2);
          table.add_value("H1", e_H1);
        }

      table.evaluate_convergence_rates(
        "L2", ConvergenceTable::reduction_rate_log2);
      table.evaluate_convergence_rates(
        "H1", ConvergenceTable::reduction_rate_log2);
      table.set_scientific("L2", true);
      table.set_scientific("H1", true);
      table.set_scientific("h",  true);

      if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
        {
          pcout << "--- convergence verification (mu = 1) ---" << std::endl;
          table.write_text(std::cout);
        }
      return 0;
    }

  const std::vector<double>       study_p_values      = {1.0, 2.0, 3.0, 4.0, 5.0};
  const std::vector<unsigned int> study_refinements   = {3, 4, 5};
  const std::vector<std::string>  study_preconditioners = {
    "none", "jacobi", "ssor", "ilu", "amg"};

  if (argc > 1)
    {
      const double       p           = std::stod(argv[1]);
      const unsigned int refinements =
        argc > 2 ? static_cast<unsigned int>(std::stoi(argv[2])) : 4;
      const std::vector<std::string> preconditioners =
        argc > 3 ? std::vector<std::string>{argv[3]} : study_preconditioners;

      HeterogeneousDiffusion problem(refinements, p, mpi_communicator);
      problem.setup_and_assemble();
      for (const auto &prec : preconditioners)
        problem.solve(prec);
      problem.output(preconditioners.back());
      return 0;
    }

  for (const unsigned int refinements : study_refinements)
    for (const double p : study_p_values)
      {
        HeterogeneousDiffusion problem(refinements, p, mpi_communicator);
        problem.setup_and_assemble();
        for (const auto &prec : study_preconditioners)
          problem.solve(prec);
      }

  return 0;
}