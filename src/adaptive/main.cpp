#include "STA_Heat.hpp"
//#define COMPARE_WITH_BASE

#ifdef COMPARE_WITH_BASE
#include "../homogeneous/H_Heat.hpp"
#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/fe/mapping_q1.h> 
#endif

int
main(int argc, char *argv[])
{
  constexpr unsigned int dim = AdaptiveHeat::dim;

  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const auto mu = [](const Point<dim> & /*p*/) { return 1.0; };
  double a = 1.5; 
  int N = 3; 
  Point<dim> x0(0, 0, 0); 
  double sigma = 0.5; 
  const auto g  = [&]( const double  &t) {
    return (std::exp(-a * (std::cos(2*N*M_PI*t) + 1)));
  };
  const auto h = [&]( const Point<dim> &p) {
    return std::exp(-((p-x0)*(p-x0)/std::pow(sigma,2)));
  };

  const auto f  = [&](const Point<dim>  &p, const double  &t) {
    return g(t)*h(p);
  };

  
#ifdef COMPARE_WITH_BASE
  
  if(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) > 1){
    if(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
      std::cout << "Baseline comparing only runs with one process for simplicity." << std::endl
      << "Either run with one process or disable comparing with baseline" << std::endl;
    return 0;
  }
  std::unique_ptr<dealii::Functions::FEFieldFunction<dim>> baseline_function;
  Heat baseline_heat(/*output_filename = */ "homogeneous.msh",
                     /* degree = */ 1,
                     /* T = */ 1.0,
                     /* theta = */ 0.5,
                     /* delta_t = */ 0.025, // use a small timestep 
                     mu,
                     f);

  Vector<double> baseline_serial_solution;

  baseline_heat.run();
  
  baseline_serial_solution.reinit(baseline_heat.get_dof_handler().n_dofs());
  baseline_serial_solution = baseline_heat.get_serial_solution();
  
  MappingQ1<dim> mapping;
  baseline_function = std::make_unique<dealii::Functions::FEFieldFunction<dim>> (
    baseline_heat.get_dof_handler(), 
    baseline_serial_solution,
    mapping
  );
  baseline_function->set_time(1.0);

#endif

  AdaptiveHeat problem(/* degree = */ 1,
               /* T = */ 1.0,
               /* theta = */ 0.5,
               /* delta_t = */ 0.05,
               mu,
               f);

  problem.run();


#ifdef COMPARE_WITH_BASE

  if(baseline_function){
    baseline_heat.print_results();
    problem.print_results();

    double L2_err = 0;

    if(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0){
      L2_err = problem.l2_against_base(*baseline_function);
    }
    MPI_Bcast(&L2_err, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }

#endif

  return 0;
}