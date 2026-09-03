// Unit tests for the parameter file and the command line. They check that the
// defaults round trip through a generated file, that an override reaches the
// bound member, and that malformed input is rejected rather than silently
// ignored.

#include "cdr_parameters.hpp"
#include "cli.hpp"
#include "test_harness.hpp"

#include <deal.II/base/parameter_handler.h>

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace dealii;

namespace {

// Builds a command line the way main() sees one, from argv[0] plus arguments.
auto command_line(const std::vector<std::string> &arguments) -> CommandLine {
  static std::vector<std::string> storage;
  storage = {"hybrid_matrix_free_solver"};
  storage.insert(storage.end(), arguments.begin(), arguments.end());

  std::vector<char *> argv;
  for (std::string &argument : storage) {
    argv.push_back(argument.data());
  }
  return CommandLine::parse(static_cast<int>(argv.size()), argv.data());
}

auto test_defaults(TestSuite &suite) -> void {
  const Parameters parameters;

  suite.check(parameters.dim == 2, "default dimension is 2");
  suite.check(parameters.degree == 2, "default degree is 2");
  suite.check(parameters.preconditioner == PreconditionerType::Multigrid,
              "multigrid is the default preconditioner");
  suite.check(parameters.beta.size() >= 3,
              "the default convection field covers 3D");
  suite.close(parameters.beta[1], 1.5, 1e-12,
              "the default beta matches the Coefficients default");
}

// The file written by --generate has to reproduce the defaults exactly, or a
// sweep started from a generated file would not be the documented baseline.
auto test_generate_round_trip(TestSuite &suite) -> void {
  Parameters defaults;
  ParameterHandler writer;
  defaults.declare(writer);

  std::ostringstream written;
  writer.print_parameters(written, ParameterHandler::OutputStyle::PRM);

  Parameters read;
  ParameterHandler reader;
  read.declare(reader);

  // Perturb every field first, so a silently skipped entry shows up.
  read.dim = 3;
  read.degree = 1;
  read.mu = 42.0;
  read.n_refinements = 99;
  read.preconditioner = PreconditionerType::None;
  read.verbose = false;

  std::istringstream stream(written.str());
  reader.parse_input(stream, "generated");

  suite.check(read.dim == defaults.dim, "round trip restores the dimension");
  suite.check(read.degree == defaults.degree, "round trip restores the degree");
  suite.close(read.mu, defaults.mu, 1e-12, "round trip restores mu");
  suite.check(read.n_refinements == defaults.n_refinements,
              "round trip restores the refinements");
  suite.check(read.preconditioner == defaults.preconditioner,
              "round trip restores the preconditioner");
  suite.check(read.verbose == defaults.verbose, "round trip restores verbose");
}

auto test_overrides(TestSuite &suite) -> void {
  Parameters parameters;
  ParameterHandler prm;
  parameters.declare(prm);

  command_line({"--set", "Solver/Preconditioner=jacobi", "--set",
                "Problem/Dimension=3", "--threads", "4", "--refine=7", "--set",
                "Problem/Beta=2.0, 0.5, 0.25"})
      .apply(prm);

  suite.check(parameters.preconditioner == PreconditionerType::Jacobi,
              "--set reaches the preconditioner");
  suite.check(parameters.dim == 3, "--set reaches the dimension");
  suite.check(parameters.n_threads == 4, "--threads is a Parallelism override");
  suite.check(parameters.n_refinements == 7,
              "--refine is a Discretization override");
  suite.close(parameters.beta[2], 0.25, 1e-12, "--set reaches a list entry");
}

// Overrides are applied after the file, so a sweep can vary one entry of a
// committed input file without editing it.
auto test_file_then_override(TestSuite &suite) -> void {
  const std::string path = "test_parameters_tmp.prm";
  {
    std::ofstream file(path);
    file << "subsection Discretization\n"
         << "  set Refinements = 3\n"
         << "end\n"
         << "subsection Solver\n"
         << "  set Preconditioner = none\n"
         << "end\n";
  }

  Parameters parameters;
  ParameterHandler prm;
  parameters.declare(prm);

  command_line({"-i", path, "--set", "Discretization/Refinements=8"})
      .apply(prm);

  suite.check(parameters.preconditioner == PreconditionerType::None,
              "the file sets the preconditioner");
  suite.check(parameters.n_refinements == 8, "the override wins over the file");

  std::remove(path.c_str());
}

auto test_actions(TestSuite &suite) -> void {
  const CommandLine help = command_line({"--help"});
  suite.check(help.action == CommandLine::Action::Help, "--help asks for help");

  const CommandLine generate = command_line({"--generate", "out.prm"});
  suite.check(generate.action == CommandLine::Action::Generate &&
                  generate.generate_file == "out.prm",
              "--generate takes a file name");

  const CommandLine to_stdout = command_line({"-g"});
  suite.check(to_stdout.action == CommandLine::Action::Generate &&
                  to_stdout.generate_file.empty(),
              "-g without a name writes to standard output");

  // Unknown arguments belong to MPI, not to us.
  const CommandLine forwarded = command_line({"--bind-to", "core"});
  suite.check(forwarded.forwarded.size() == 3,
              "unknown arguments are forwarded to MPI");
}

auto test_rejects_bad_input(TestSuite &suite) -> void {
  const auto throws = [](const auto &callable) {
    try {
      callable();
    } catch (const std::exception &) {
      return true;
    }
    return false;
  };

  Parameters parameters;
  ParameterHandler prm;
  parameters.declare(prm);

  suite.check(
      throws([&] { command_line({"--set", "Bogus/Entry=1"}).apply(prm); }),
      "an undeclared entry is rejected");
  suite.check(throws([&] {
                command_line({"--set", "Solver/Preconditioner=amg"}).apply(prm);
              }),
              "an unknown preconditioner is rejected");
  suite.check(
      throws([&] { command_line({"--set", "Tolerance=1"}).apply(prm); }),
      "a path without a section is rejected");
  suite.check(
      throws([&] { command_line({"-i", "does_not_exist.prm"}).apply(prm); }),
      "a missing input file is rejected");
  suite.check(throws([&] { command_line({"--input"}); }),
              "an option without its value is rejected");

  Parameters inconsistent;
  inconsistent.dim = 3;
  inconsistent.beta = {1.0, 1.0};
  suite.check(throws([&] { inconsistent.validate(); }),
              "a beta shorter than the dimension is rejected");
}

} // namespace

int main() {
  TestSuite suite("parameters");

  test_defaults(suite);
  test_generate_round_trip(suite);
  test_overrides(suite);
  test_file_then_override(suite);
  test_actions(suite);
  test_rejects_bad_input(suite);

  return suite.report();
}
