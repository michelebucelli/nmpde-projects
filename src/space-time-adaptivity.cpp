#include "Heat.hpp"

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/numerics/fe_field_function.h>

#include <algorithm>
#include <fstream>
#include <vector>

namespace
{
  struct Params
  {
    bool                      adaptive     = true;
    bool                      convergence  = false;
    unsigned int              n_elements   = 12;
    unsigned int              degree       = 1;
    double                    T            = 1.0;
    double                    theta        = 1.0;
    double                    delta_t      = 1.0e-3;
    double                    tol_time_max = 1.0e-4;
    double                    tol_time_min = 2.0e-6;
    double                    dt_max       = 0.01;
    std::vector<unsigned int> n_el_list    = {5, 10, 20, 40};
  };

  Params
  parse_params(int argc, char *argv[])
  {
    ParameterHandler prm;
    prm.declare_entry("mode",
                       "adaptive",
                       Patterns::Selection("adaptive|nonadaptive"),
                       "Whether the run uses AMR + adaptive time-stepping, or a fixed "
                       "uniform mesh + fixed delta_t baseline.");
    prm.declare_entry("run",
                       "default",
                       Patterns::Selection("default|convergence"),
                       "'default' is a single run at n_elements. 'convergence' runs a "
                       "reference (the largest value in n_el_list) plus a sweep over the "
                       "rest of n_el_list, and is slow.");
    prm.declare_entry("n_elements",
                       "12",
                       Patterns::Integer(1),
                       "Base mesh resolution for the 'default' run.");
    prm.declare_entry("degree", "1", Patterns::Integer(1), "Polynomial degree.");
    prm.declare_entry("T", "1.0", Patterns::Double(0.0), "Final time.");
    prm.declare_entry(
      "theta", "1.0", Patterns::Double(0.0, 1.0), "Theta-method parameter (1=Backward Euler).");
    prm.declare_entry(
      "delta_t", "1e-3", Patterns::Double(0.0), "Initial (adaptive) or fixed (non-adaptive) time step.");
    prm.declare_entry("tol_time_max",
                       "1e-4",
                       Patterns::Double(0.0),
                       "Reject/halve delta_t if the time error indicator exceeds this. Keep tight "
                       "for problems with a sharp feature (like the pulsed forcing) that a "
                       "too-large step could step over undetected.");
    prm.declare_entry("tol_time_min",
                       "2e-6",
                       Patterns::Double(0.0),
                       "Grow delta_t by 1.2x once accepted if the error indicator is below this.");
    prm.declare_entry(
      "dt_max", "0.01", Patterns::Double(0.0), "Ceiling on delta_t during adaptive growth.");
    prm.declare_entry("n_el_list",
                       "5,10,20,40",
                       Patterns::List(Patterns::Integer(1)),
                       "Mesh sizes for the 'convergence' run; the largest becomes the "
                       "reference, the rest are the sweep.");

    if (argc > 1)
      {
        prm.parse_input(argv[1]);
      }
    else if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
      {
        std::ofstream default_prm("space-time-adaptivity.prm");
        prm.print_parameters(default_prm, ParameterHandler::Text);
      }

    Params params;
    params.adaptive    = (prm.get("mode") == "adaptive");
    params.convergence = (prm.get("run") == "convergence");
    params.n_elements  = prm.get_integer("n_elements");
    params.degree      = prm.get_integer("degree");
    params.T           = prm.get_double("T");
    params.theta       = prm.get_double("theta");
    params.delta_t     = prm.get_double("delta_t");
    params.tol_time_max = prm.get_double("tol_time_max");
    params.tol_time_min = prm.get_double("tol_time_min");
    params.dt_max        = prm.get_double("dt_max");

    params.n_el_list.clear();
    for (const auto &token : Utilities::split_string_list(prm.get("n_el_list")))
      params.n_el_list.push_back(Utilities::string_to_int(token));

    return params;
  }

  // f = g(t) * h(x); g_pulsation/h_spatial must outlive the returned lambda's use.
  std::function<double(const Point<Heat::dim> &, const double &)>
  make_forcing(Heat::g_function &g_pulsation, Heat::h_function &h_spatial)
  {
    return [&g_pulsation, &h_spatial](const Point<Heat::dim> &x, const double &t) {
      g_pulsation.set_time(t);
      return g_pulsation.value(x) * h_spatial.value(x);
    };
  }

  void
  run_default(const Params &params, ConditionalOStream &pcout)
  {
    const auto       mu = [](const Point<Heat::dim> &) { return 1.0; };
    Heat::g_function g_pulsation;
    Heat::h_function h_spatial;
    const auto       f = make_forcing(g_pulsation, h_spatial);

    pcout << "Running default " << (params.adaptive ? "adaptive" : "non-adaptive")
          << " simulation (n_elements=" << params.n_elements << ", degree=" << params.degree
          << ", T=" << params.T << ", theta=" << params.theta << ", delta_t=" << params.delta_t
          << ")" << std::endl;

    Heat problem(params.n_elements,
                 params.degree,
                 params.T,
                 params.theta,
                 params.delta_t,
                 mu,
                 f,
                 params.adaptive,
                 params.tol_time_max,
                 params.tol_time_min,
                 params.dt_max);
    problem.run();
  }

  void
  run_convergence(const Params &params, ConditionalOStream &pcout)
  {
    const auto       mu = [](const Point<Heat::dim> &) { return 1.0; };
    Heat::g_function g_pulsation;
    Heat::h_function h_spatial;
    const auto       f = make_forcing(g_pulsation, h_spatial);

    std::vector<unsigned int> n_el_list = params.n_el_list;
    std::sort(n_el_list.begin(), n_el_list.end());
    const unsigned int reference_n_el = n_el_list.back();
    n_el_list.pop_back();

    pcout << "Running " << (params.adaptive ? "adaptive" : "non-adaptive")
          << " convergence test: reference n_elements=" << reference_n_el << ", sweep = {";
    for (const auto &n : n_el_list)
      pcout << n << " ";
    pcout << "}" << std::endl;

    Heat reference(reference_n_el,
                    params.degree,
                    params.T,
                    params.theta,
                    params.delta_t,
                    mu,
                    f,
                    params.adaptive,
                    params.tol_time_max,
                    params.tol_time_min,
                    params.dt_max);
    reference.run();

    Functions::FEFieldFunction<Heat::dim, TrilinosWrappers::MPI::Vector> reference_solution(
      reference.get_dof_handler(), reference.get_solution());

    ConvergenceTable  table;
    const std::string csv_name =
      "convergence-" + std::string(params.adaptive ? "adaptive" : "nonadaptive") + ".csv";
    std::ofstream convergence_file(csv_name);
    convergence_file << "h,eL2,eH1,type\n";

    for (const auto &N_el : n_el_list)
      {
        Heat problem(N_el,
                     params.degree,
                     params.T,
                     params.theta,
                     params.delta_t,
                     mu,
                     f,
                     params.adaptive,
                     params.tol_time_max,
                     params.tol_time_min,
                     params.dt_max);

        problem.run();

        const double h = 1.0 / N_el;

        const double error_L2 = problem.compute_error(VectorTools::L2_norm, reference_solution);
        const double error_H1 = problem.compute_error(VectorTools::H1_norm, reference_solution);

        table.add_value("h", h);
        table.add_value("L2", error_L2);
        table.add_value("H1", error_H1);

        convergence_file << h << "," << error_L2 << "," << error_H1 << ",sweep\n";
      }

    table.set_precision("L2", 6);
    table.set_scientific("L2", true);
    table.set_precision("H1", 6);
    table.set_scientific("H1", true);
    // Use the "h"-referenced overload: n_el_list isn't necessarily a clean
    // doubling sequence (e.g. 3,6,9 has a 1.5x step), and the no-reference
    // overload silently assumes h halves every row, which would misreport
    // the rate on any non-doubling step.
    table.evaluate_all_convergence_rates("h", ConvergenceTable::reduction_rate_log2);
    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
      table.write_text(std::cout);

    // The reference mesh is the ground truth the sweep is compared against,
    // so its own error is trivially zero and not a real convergence-rate
    // data point; record it separately (not folded into the rate table)
    // purely so downstream table/plot can mark where it sits on the h axis.
    const double h_reference = 1.0 / reference_n_el;
    convergence_file << h_reference << ",,,reference\n";
    pcout << "(reference: n_elements=" << reference_n_el << ", h=" << h_reference
          << " -- ground truth, not a sweep/rate point)" << std::endl;
  }
} // namespace

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
  ConditionalOStream pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0);

  const Params params = parse_params(argc, argv);

  if (params.convergence)
    run_convergence(params, pcout);
  else
    run_default(params, pcout);

  return 0;
}
