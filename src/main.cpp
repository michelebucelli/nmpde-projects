#include "cdr_parameters.hpp"
#include "cdr_solver.hpp"
#include "cli.hpp"

#include <deal.II/base/multithread_info.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/utilities.h>

#include <exception>
#include <fstream>
#include <iostream>

using namespace dealii;

namespace {

// dim and degree are template arguments of the whole solver, so the input
// file selects an instantiation rather than a runtime value.
template <int dim, int degree> auto run(const Parameters &parameters) -> void {
  CdrProblem<dim, degree> problem(parameters);
  problem.run();
}

auto dispatch(const Parameters &parameters) -> void {
  const auto unsupported = [&]() {
    AssertThrow(false,
                ExcMessage("no instantiation for dimension " +
                           std::to_string(parameters.dim) + " and degree " +
                           std::to_string(parameters.degree)));
  };

  switch (parameters.dim) {
  case 2:
    switch (parameters.degree) {
    case 1:
      return run<2, 1>(parameters);
    case 2:
      return run<2, 2>(parameters);
    case 3:
      return run<2, 3>(parameters);
    default:
      return unsupported();
    }
  case 3:
    switch (parameters.degree) {
    case 1:
      return run<3, 1>(parameters);
    case 2:
      return run<3, 2>(parameters);
    case 3:
      return run<3, 3>(parameters);
    default:
      return unsupported();
    }
  default:
    return unsupported();
  }
}

} // namespace

int main(int argc, char *argv[]) {
  try {
    const CommandLine cli = CommandLine::parse(argc, argv);

    if (cli.action == CommandLine::Action::Help) {
      std::cout << CommandLine::usage();
      return 0;
    }

    // Reading the input file needs no MPI, and its Threads entry has to be
    // known before MPI_InitFinalize sizes the TBB pool.
    Parameters parameters;
    ParameterHandler prm;
    parameters.declare(prm);

    if (cli.action == CommandLine::Action::Generate) {
      if (cli.generate_file.empty()) {
        prm.print_parameters(std::cout, ParameterHandler::OutputStyle::PRM);
      } else {
        std::ofstream file(cli.generate_file);
        AssertThrow(file,
                    ExcMessage("cannot write '" + cli.generate_file + "'"));
        prm.print_parameters(file, ParameterHandler::OutputStyle::PRM);
        std::cout << "wrote " << cli.generate_file << std::endl;
      }
      return 0;
    }

    cli.apply(prm);
    parameters.validate();

    int forwarded_argc = static_cast<int>(cli.forwarded.size());
    char **forwarded_argv = const_cast<char **>(cli.forwarded.data());

    const unsigned int thread_limit = parameters.n_threads == 0
                                          ? numbers::invalid_unsigned_int
                                          : parameters.n_threads;

    // deal.II calls Kokkos::initialize() before it applies the
    // MPI_InitFinalize thread argument. Kokkos then keeps a persistent
    // thread pool. If that pool has the same size as TBB, the two runtimes
    // oversubscribe the cores and MatrixFree's cell_loop does not speed up.
    // Pin Kokkos to one thread first; MPI_InitFinalize then sizes TBB.
    MultithreadInfo::set_thread_limit(1);

    Utilities::MPI::MPI_InitFinalize mpi_initialization(
        forwarded_argc, forwarded_argv, thread_limit);

    dispatch(parameters);
  } catch (const std::exception &exception) {
    std::cerr << "error: " << exception.what() << std::endl;
    return 1;
  }

  return 0;
}
