#pragma once

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/patterns.h>

#include <stdexcept>
#include <string>
#include <vector>

using namespace dealii;

enum class PreconditionerType { None, Jacobi, Multigrid };

enum class Backend { MatrixFree, MatrixBased };

// Wall times of a single refinement cycle, for the scaling study.
struct CycleTimings {
  double setup = 0.0;
  double rhs = 0.0;
  double solve = 0.0;
};

struct SolveResult {
  unsigned int n_iterations = 0;
  double relative_residual = 0.0;
};

// Every knob of a run, populated from a .prm input file. The defaults below
// are the ones written by --generate, so a generated file always reproduces
// the behaviour of running with no input file at all.
struct Parameters {
  // Problem
  unsigned int dim = 2;
  unsigned int degree = 2;
  double mu = 1.0;
  double gamma = 1.0;
  // Convection field. Entries past `dim` are ignored; the default reproduces
  // the beta_d = 1 + d/2 hard-coded in Coefficients.
  std::vector<double> beta = {1.0, 1.5, 2.0};

  // Discretization
  unsigned int n_refinements = 5;
  unsigned int n_cycles = 1;

  // Solver
  Backend backend = Backend::MatrixFree;
  PreconditionerType preconditioner = PreconditionerType::Multigrid;
  double tolerance = 1e-9;

  // Parallelism. TBB threads per MPI rank; zero means the deal.II/TBB
  // default, which is one thread per core.
  unsigned int n_threads = 0;

  // Output
  bool verbose = true;
  bool write_output = false;
  // Number of operator applications timed by the builtin benchmark; zero
  // skips the benchmark.
  unsigned int n_benchmarks_applications = 0;

  auto declare(ParameterHandler &prm) -> void;
  auto validate() const -> void;
};

inline auto to_backend(const std::string &name) -> Backend {
  if (name == "matrix_free")
    return Backend::MatrixFree;
  if (name == "matrix_based")
    return Backend::MatrixBased;
  throw std::runtime_error("unknown backend '" + name + "'");
}

inline auto to_string(const Backend backend) -> std::string {
  return backend == Backend::MatrixFree ? "matrix_free" : "matrix_based";
}

inline auto to_preconditioner(const std::string &name) -> PreconditionerType {
  if (name == "none")
    return PreconditionerType::None;
  if (name == "jacobi")
    return PreconditionerType::Jacobi;
  if (name == "multigrid")
    return PreconditionerType::Multigrid;
  throw std::runtime_error("unknown preconditioner '" + name + "'");
}

inline auto to_string(const PreconditionerType preconditioner) -> std::string {
  switch (preconditioner) {
  case PreconditionerType::None:
    return "none";
  case PreconditionerType::Jacobi:
    return "jacobi";
  case PreconditionerType::Multigrid:
    return "multigrid";
  }
  return "multigrid";
}

// Binds every field to an entry of `prm`. add_parameter() writes straight
// back into the member when the file is parsed or an entry is overridden with
// ParameterHandler::set, so there is no separate read-back step.
inline auto Parameters::declare(ParameterHandler &prm) -> void {
  prm.enter_subsection("Problem");
  {
    prm.add_parameter("Dimension", dim, "Space dimension, 2 or 3.",
                      Patterns::Integer(2, 3));
    prm.add_parameter("Degree", degree, "Lagrange element degree, 1 to 3.",
                      Patterns::Integer(1, 3));
    prm.add_parameter("Mu", mu, "Diffusion coefficient, positive.",
                      Patterns::Double(0.0));
    prm.add_parameter("Gamma", gamma, "Reaction coefficient, non negative.",
                      Patterns::Double(0.0));
    prm.add_parameter("Beta", beta,
                      "Convection field, comma separated. Entries past the "
                      "space dimension are ignored.",
                      Patterns::List(Patterns::Double(), 1, 3));
  }
  prm.leave_subsection();

  prm.enter_subsection("Discretization");
  {
    prm.add_parameter("Refinements", n_refinements,
                      "Global refinements of the coarsest cycle.",
                      Patterns::Integer(1));
    prm.add_parameter("Cycles", n_cycles,
                      "Number of refinement cycles. Convergence rates are "
                      "only tabulated when this is larger than one.",
                      Patterns::Integer(1));
  }
  prm.leave_subsection();

  prm.enter_subsection("Solver");
  {
    prm.declare_entry("Backend", to_string(backend),
                      Patterns::Selection("matrix_free|matrix_based"),
                      "Operator representation. 'matrix_based' assembles a "
                      "sparse matrix and is the reference for the "
                      "matrix-free benchmark.");
    prm.add_action("Backend", [this](const std::string &value) {
      backend = to_backend(value);
    });

    // Held as an enum, so the action converts instead of add_parameter.
    prm.declare_entry("Preconditioner", to_string(preconditioner),
                      Patterns::Selection("none|jacobi|multigrid"),
                      "FGMRES preconditioner. 'multigrid' means the "
                      "geometric V-cycle for the matrix-free backend and "
                      "algebraic multigrid for the matrix-based one.");
    prm.add_action("Preconditioner", [this](const std::string &value) {
      preconditioner = to_preconditioner(value);
    });
    prm.add_parameter("Tolerance", tolerance,
                      "Relative residual stopping tolerance.",
                      Patterns::Double(0.0));
  }
  prm.leave_subsection();

  prm.enter_subsection("Parallelism");
  {
    prm.add_parameter("Threads", n_threads,
                      "TBB threads per MPI rank. Zero means one per core. "
                      "The MPI rank count comes from mpirun, not from here.",
                      Patterns::Integer(0));
  }
  prm.leave_subsection();

  prm.enter_subsection("Output");
  {
    prm.add_parameter("Verbose", verbose,
                      "Print the per cycle log and the convergence table.");
    prm.add_parameter("Write VTU", write_output,
                      "Write the solution of each cycle to a .vtu record.");
    prm.add_parameter("Benchmark applications", n_benchmarks_applications,
                      "Timed bare operator applications per cycle. Zero "
                      "skips the benchmark.",
                      Patterns::Integer(0));
  }
  prm.leave_subsection();
}

// Patterns cover the per entry ranges; this catches what only makes sense
// once several entries are read together.
inline auto Parameters::validate() const -> void {
  if (beta.size() < dim) {
    throw std::runtime_error("Problem/Beta has " + std::to_string(beta.size()) +
                             " components, but Problem/Dimension is " +
                             std::to_string(dim));
  }
  if (mu <= 0.0) {
    throw std::runtime_error("Problem/Mu must be positive");
  }
  if (tolerance <= 0.0) {
    throw std::runtime_error("Solver/Tolerance must be positive");
  }
}
