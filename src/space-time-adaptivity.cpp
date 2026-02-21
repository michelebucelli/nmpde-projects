#include "Heat.hpp"

//  du/dt -div(mu * grad(u)) = f
int
main(int argc, char *argv[])
{
  constexpr unsigned int dim = Heat::dim;

  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  // Coefficiente di Diffusione (mu)
  const auto mu = [](const Point<dim> & /*p*/) { 
    return 1.0; 
  };

  Heat::g_function g_pulsation;
  Heat::h_function h_spatial;

  // Forzante (f)
  const auto f = [&g_pulsation, &h_spatial](const Point<dim> &p, const double &t) {
  // Impostiamo il tempo corrente per la funzione g(t)
  g_pulsation.set_time(t);
  
  // Il risultato è il prodotto g(t) * h(x) 
  return g_pulsation.value(p) * h_spatial.value(p);
  };

  // Istanziazione del problema
  
  Heat problem(/*mesh_file_name = */ "../mesh/mesh-cube-10.msh",
               /* degree = */ 2,
               /* T = */ 1.0,
               /* theta = */ 1.0,    
               /* delta_t = */ 0.001,
               mu,
               f);    

  problem.run();

  return 0;
}
