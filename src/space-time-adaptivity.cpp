#include "Heat.hpp"

//  du/dt -div(mu * grad(u)) + div(b * u) + sigma * u = f
int
main(int argc, char *argv[])
{
  constexpr unsigned int dim = Heat::dim;

  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  // 1. Coefficiente di Diffusione (mu)
  const auto mu = [](const Point<dim> & /*p*/) { 
    return 1.0; 
  };

  // 2. Coefficiente di Reazione (sigma) 
  const auto sigma = [](const Point<dim> & /*p*/) { 
    return 0.0; 
  };

  Heat::g_function g_pulsation;
  Heat::h_function h_spatial;

  // 3. Forzante (f)
  const auto f = [&g_pulsation, &h_spatial](const Point<dim> &p, const double &t) {
  // Impostiamo il tempo corrente per la funzione g(t)
  g_pulsation.set_time(t);
  
  // Il risultato è il prodotto g(t) * h(x) 
  return g_pulsation.value(p) * h_spatial.value(p);
  };

  // 4. Campo di Velocità (b) 
  const auto b = [](const Point<dim> & /*p*/) {
    Tensor<1, dim> velocity;
    for (unsigned int d = 0; d < dim; ++d)
      velocity[d] = 0.0; 
    return velocity;
  };

  // Istanziazione del problema
  
  // Heat problem(/*mesh_file_name = */ "../mesh/mesh-cube-10.msh",
  // Heat problem(/*mesh_file_name = */ "../mesh/mesh-cube-20.msh", // 109 iters
  Heat problem(/*mesh_file_name = */ "../mesh/mesh-cube-40.msh", //201 iters
               /* degree = */ 2,
               /* T = */ 1.0,
               /* theta = */ 1.0,    
              //  /* delta_t = */ 0.001,
               /* delta_t = */ 0.05,
               mu,    
               sigma, 
               f,    
               b);    

  problem.run();

  return 0;
}
