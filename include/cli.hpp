#pragma once

#include "cdr_parameters.hpp"

#include <deal.II/base/parameter_handler.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Command line of the solver. The run itself is configured by a .prm input
// file; the command line only says which file to read, what to override for a
// single run of a sweep, and which arguments to hand on to MPI.
struct CommandLine {
  // What to do once the command line is understood.
  enum class Action { Run, Help, Generate };

  Action action = Action::Run;

  // Input file. Empty means run with the built in defaults.
  std::string input_file;

  // Destination of --generate. Empty means standard output.
  std::string generate_file;

  // "Section/Entry=value" overrides applied after the file is read, so a
  // scaling sweep can vary one entry without writing a file per data point.
  std::vector<std::string> overrides;

  // argv entries that are none of our business, forwarded to MPI_InitFinalize
  // with argv[0] in front.
  std::vector<char *> forwarded;

  static auto parse(int argc, char *argv[]) -> CommandLine;
  static auto usage() -> std::string;

  // Reads `input_file` into `prm`, then applies `overrides`. Both write
  // straight into the Parameters bound to `prm` by Parameters::declare.
  auto apply(ParameterHandler &prm) const -> void;
};

inline auto CommandLine::usage() -> std::string {
  return R"(Hybrid MPI + TBB matrix-free CDR solver.

Usage:
  hybrid_matrix_free_solver [options] [-- mpi arguments]

Options:
  -i, --input FILE        Read parameters from FILE. Without it the built in
                          defaults are used.
  -s, --set PATH=VALUE    Override one entry after FILE is read, for instance
                          --set Solver/Preconditioner=jacobi. Repeatable.
  -g, --generate [FILE]   Write a fully documented parameter file to FILE, or
                          to standard output, and exit.
      --threads N         Shorthand for --set Parallelism/Threads=N.
      --refine N          Shorthand for --set Discretization/Refinements=N.
  -h, --help              Print this text and exit.

The MPI rank count is not a parameter; it comes from mpirun -n.

Examples:
  hybrid_matrix_free_solver --generate input.prm
  hybrid_matrix_free_solver -i input.prm
  mpirun -n 4 hybrid_matrix_free_solver -i input.prm --threads 2 \
      --set Solver/Preconditioner=jacobi
)";
}

inline auto CommandLine::parse(int argc, char *argv[]) -> CommandLine {
  CommandLine cli;
  cli.forwarded.push_back(argv[0]);

  // Matches "-x value", "--long value", "-x=value" or "--long=value" and
  // returns the value. Returns false when `arg` is not this option at all.
  const auto option = [&](int &i, const std::string &arg,
                          const std::string &shrt, const std::string &lng,
                          std::string &value) -> bool {
    for (const std::string &name : {shrt, lng}) {
      if (arg == name) {
        if (i + 1 >= argc) {
          throw std::runtime_error(name + " needs a value");
        }
        value = argv[++i];
        return true;
      }
      if (arg.rfind(name + "=", 0) == 0) {
        value = arg.substr(name.size() + 1);
        return true;
      }
    }
    return false;
  };

  std::string value;

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];

    if (arg == "-h" || arg == "--help") {
      cli.action = Action::Help;
      return cli;
    }

    if (arg == "--") { // everything after this belongs to MPI
      for (int j = i + 1; j < argc; ++j) {
        cli.forwarded.push_back(argv[j]);
      }
      break;
    }

    if (option(i, arg, "-i", "--input", value)) {
      cli.input_file = value;
      continue;
    }

    if (option(i, arg, "-s", "--set", value)) {
      cli.overrides.push_back(value);
      continue;
    }

    if (arg == "-g" || arg == "--generate" ||
        arg.rfind("--generate=", 0) == 0 || arg.rfind("-g=", 0) == 0) {
      cli.action = Action::Generate;
      if (arg.rfind("--generate=", 0) == 0) {
        cli.generate_file = arg.substr(11);
      } else if (arg.rfind("-g=", 0) == 0) {
        cli.generate_file = arg.substr(3);
      } else if (i + 1 < argc && argv[i + 1][0] != '-') {
        cli.generate_file = argv[++i];
      }
      continue;
    }

    // Shorthands for the two entries a scaling sweep varies most often.
    if (option(i, arg, "--threads", "--threads", value)) {
      cli.overrides.push_back("Parallelism/Threads=" + value);
      continue;
    }

    if (option(i, arg, "--refine", "--refine", value)) {
      cli.overrides.push_back("Discretization/Refinements=" + value);
      continue;
    }

    cli.forwarded.push_back(argv[i]);
  }

  return cli;
}

inline auto CommandLine::apply(ParameterHandler &prm) const -> void {
  if (!input_file.empty()) {
    std::ifstream file(input_file);
    if (!file) {
      throw std::runtime_error("cannot open input file '" + input_file + "'");
    }
    prm.parse_input(file, input_file);
  }

  for (const std::string &override : overrides) {
    const auto equals = override.find('=');
    if (equals == std::string::npos) {
      throw std::runtime_error("--set expects PATH=VALUE, got '" + override +
                               "'");
    }

    const std::string path = override.substr(0, equals);
    const std::string value = override.substr(equals + 1);

    // Split "Section/Subsection/Entry" and walk into the subsections, since
    // ParameterHandler::set only sees the current one.
    std::vector<std::string> parts;
    std::stringstream stream(path);
    std::string part;
    while (std::getline(stream, part, '/')) {
      parts.push_back(part);
    }
    if (parts.size() < 2) {
      throw std::runtime_error("--set path '" + path +
                               "' needs a section, as in Solver/Tolerance");
    }

    for (unsigned int i = 0; i + 1 < parts.size(); ++i) {
      prm.enter_subsection(parts[i]);
    }
    try {
      prm.set(parts.back(), value); // runs the action bound by declare()
    } catch (const std::exception &) {
      // deal.II reports an undeclared entry or a pattern mismatch with a
      // stack trace, which is noise for a typo on the command line.
      throw std::runtime_error("cannot set '" + path + "' to '" + value +
                               "'. Run --generate to list the valid entries");
    }
    for (unsigned int i = 0; i + 1 < parts.size(); ++i) {
      prm.leave_subsection();
    }
  }
}
