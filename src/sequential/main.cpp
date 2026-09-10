#include "WaveEquation.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

namespace
{
void
print_usage(const char *program_name)
{
  std::cout << "Usage:\n"
            << "  " << program_name
            << " [n_subdivisions] [degree] [final_time] [time_step] "
               "[output_frequency] [output_directory] [validation_case] "
               "[solver] [preconditioner]\n\n"
            << "Defaults:\n"
            << "  n_subdivisions   = 24\n"
            << "  degree           = 1\n"
            << "  final_time       = 1.0\n"
            << "  time_step        = 0.005\n"
            << "  output_frequency = 10\n"
            << "  output_directory = output\n"
            << "  validation_case  = free\n"
            << "  solver = cg (or gmres), preconditioner = ssor (or jacobi)\n"
            << "  T/dt must be an integer; WAVE_CHECK_MATRIX=1 enables sparse "
               "checks.\n\n"
            << "Validation cases:\n"
            << "  free   : exact solution with f = 0\n"
            << "  forced : same exact solution shape, but with nonzero f\n";
}

unsigned int
positive_integer(const char *argument)
{
  const std::string text(argument);
  if (text.empty() || text.find_first_not_of("0123456789") != std::string::npos)
    throw std::invalid_argument("Expected a positive integer: " + text);
  const auto value = std::stoul(text);
  if (value == 0 || value > std::numeric_limits<unsigned int>::max())
    throw std::invalid_argument("Integer out of range: " + text);
  return static_cast<unsigned int>(value);
}

double
positive_real(const char *argument)
{
  std::size_t  end = 0;
  const double value = std::stod(argument, &end);
  if (end != std::string(argument).size() || !std::isfinite(value)
      || value <= 0)
    throw std::invalid_argument("Expected a finite positive number");
  return value;
}
} // namespace

int
main(int argc, char *argv[])
{
  try
    {
      // Default parameters. They are intentionally modest so the first run is
      // quick in the course Docker environment.
      unsigned int n_subdivisions = 24;
      unsigned int degree = 1;
      double       final_time = 1.0;
      double       time_step = 0.005;
      unsigned int output_frequency = 10;
      std::string  output_directory = "output";
      std::string  validation_case = "free";
      std::string  solver_name = "cg";
      std::string  preconditioner_name = "ssor";

      if (argc == 2 && std::string(argv[1]) == "--help")
        {
          print_usage(argv[0]);
          return 0;
        }

      if (argc > 1)
        n_subdivisions = positive_integer(argv[1]);
      if (argc > 2)
        degree = positive_integer(argv[2]);
      if (argc > 3)
        final_time = positive_real(argv[3]);
      if (argc > 4)
        time_step = positive_real(argv[4]);
      if (argc > 5)
        output_frequency = positive_integer(argv[5]);
      if (argc > 6)
        output_directory = argv[6];
      if (argc > 7)
        validation_case = argv[7];

      if (argc > 8)
        solver_name = argv[8];
      if (argc > 9)
        preconditioner_name = argv[9];

      if (argc > 10)
        {
          print_usage(argv[0]);
          return 1;
        }

      const double steps = final_time / time_step;
      if (!std::isfinite(steps) || steps < 1
          || steps > std::numeric_limits<unsigned int>::max()
          || std::abs(steps - std::round(steps))
                 > 1.e-10 * std::max(1.0, steps))
        throw std::invalid_argument("final_time/time_step must be a positive "
                                    "integer (fixed Newmark matrix)");
      if (output_directory.empty())
        throw std::invalid_argument("Output directory cannot be empty");
      if (solver_name != "cg" && solver_name != "gmres")
        throw std::invalid_argument("Solver must be cg or gmres");
      if (preconditioner_name != "ssor" && preconditioner_name != "jacobi")
        throw std::invalid_argument("Preconditioner must be ssor or jacobi");

      if (validation_case != "free" && validation_case != "forced")
        {
          std::cerr << "Unknown validation_case '" << validation_case
                    << "'. Use 'free' or 'forced'." << std::endl;
          return 1;
        }

      WaveEquation problem(n_subdivisions, degree, final_time, time_step,
                           output_frequency, output_directory, validation_case,
                           solver_name, preconditioner_name);

      problem.run();

      return 0;
    }
  catch (const std::exception &error)
    {
      std::cerr << "Error: " << error.what() << '\n';
      return 1;
    }
}
