#include "Heat.hpp"
#include <deal.II/base/convergence_table.h>

// class ExactSolution : public Function<3>
// {
// public:
//   ExactSolution()
//     : Function<3>()
//   {
//     precompute_projection();
//   }

//   virtual double value(const Point<3> &p,
//                        const unsigned int = 0) const override
//   {
//     const double t = this->get_time();
//     double u = 0.0;

//     for (const auto &mode : modes)
//     {
//       const double lambda =
//         numbers::PI * numbers::PI *
//         (mode.k1*mode.k1 +
//          mode.k2*mode.k2 +
//          mode.k3*mode.k3);

//       const double phi =
//         std::cos(mode.k1*numbers::PI*p[0]) *
//         std::cos(mode.k2*numbers::PI*p[1]) *
//         std::cos(mode.k3*numbers::PI*p[2]);

//       const double a_k =
//         mode.Hk * time_convolution(lambda, t);

//       u += a_k * phi;
//     }

//     return u;
//   }

// private:

//   struct Mode
//   {
//     unsigned int k1,k2,k3;
//     double Hk;
//   };

//   std::vector<Mode> modes;

//   double g(const double t) const
//   {
//     const unsigned int N = 5;
//     return std::sin(N*numbers::PI*t)
//            * std::sin(N*numbers::PI*t);
//   }

//   double h(const Point<3> &p) const
//   {
//     const Point<3> x0(0.5,0.5,0.5);
//     const double sigma = 0.08;
//     return std::exp(-(p-x0).norm_square()
//                      /(sigma*sigma));
//   }

//   void precompute_projection()
//   {
//     const unsigned int max_k = 3;
//     QGauss<3> quad(6);

//     for (unsigned int k1=0;k1<=max_k;++k1)
//       for (unsigned int k2=0;k2<=max_k;++k2)
//         for (unsigned int k3=0;k3<=max_k;++k3)
//         {
//           double Hk = 0.0;

//           for (unsigned int q=0;q<quad.size();++q)
//           {
//             const Point<3> p = quad.point(q);

//             const double phi =
//               std::cos(k1*numbers::PI*p[0]) *
//               std::cos(k2*numbers::PI*p[1]) *
//               std::cos(k3*numbers::PI*p[2]);

//             Hk += h(p) * phi * quad.weight(q);
//           }

//           modes.push_back({k1,k2,k3,Hk});
//         }
//   }

//   double time_convolution(const double lambda,
//                           const double t) const
//   {
//     if (t == 0.0)
//       return 0.0;

//     const unsigned int steps = 100;
//     const double dt = t / steps;

//     double integral = 0.0;

//     for (unsigned int i=0;i<steps;++i)
//     {
//       const double s = i*dt;
//       integral += std::exp(-lambda*(t-s))
//                   * g(s) * dt;
//     }

//     return integral;
//   }
// };

using namespace dealii;

constexpr unsigned int dim = 3;

// =====================================================
// Truncated Spectral Exact Solution (Neumann, 3D)
// =====================================================
class ExactSolution : public Function<dim>
{
public:
  ExactSolution(Heat::g_function &g_fun,
                Heat::h_function &h_fun,
                const unsigned int M = 4)
    : Function<dim>(),
      g(g_fun),
      h(h_fun),
      Mmax(M)
  {}

  virtual double value(const Point<dim> &p,
                       const unsigned int = 0) const override
  {
    const double x = p[0];
    const double y = p[1];
    const double z = p[2];
    const double t = this->get_time();

    double result = 0.0;

    for (unsigned int m = 0; m <= Mmax; ++m)
      for (unsigned int n = 0; n <= Mmax; ++n)
        for (unsigned int l = 0; l <= Mmax; ++l)
        {
          const double lambda =
            M_PI * M_PI * (m*m + n*n + l*l);

          const double h_mnl = compute_h_coeff(m,n,l);

          const double time_part =
            time_integral(lambda, h_mnl, t);

          result += time_part *
                    std::cos(m * M_PI * x) *
                    std::cos(n * M_PI * y) *
                    std::cos(l * M_PI * z);
        }

    return result;
  }

private:
  Heat::g_function &g;
  Heat::h_function &h;
  unsigned int Mmax;

  double compute_h_coeff(unsigned int m,
                         unsigned int n,
                         unsigned int l) const
  {
    const unsigned int Q = 6;
    double sum = 0.0;

    for (unsigned int i = 0; i < Q; ++i)
      for (unsigned int j = 0; j < Q; ++j)
        for (unsigned int k = 0; k < Q; ++k)
        {
          Point<dim> p((i+0.5)/Q,
                       (j+0.5)/Q,
                       (k+0.5)/Q);

          sum += h.value(p) *
                 std::cos(m*M_PI*p[0]) *
                 std::cos(n*M_PI*p[1]) *
                 std::cos(l*M_PI*p[2]);
        }

    return sum / (Q*Q*Q);
  }

  double time_integral(double lambda,
                       double h_coeff,
                       double t) const
  {
    const unsigned int Q = 20;
    double sum = 0.0;

    for (unsigned int i = 0; i < Q; ++i)
    {
      double s = t * (i+0.5)/Q;
      g.set_time(s);

      sum += std::exp(-lambda*(t-s))
             * g.value(Point<dim>())
             * h_coeff;
    }

    return sum * t / Q;
  }
};


//  du/dt -div(mu * grad(u)) = f
int
main(int argc, char *argv[])
{
  constexpr unsigned int dim = Heat::dim;

  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  // Check if the user provided the mesh name argument
  if (argc < 2)
  {
    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      std::cerr << "Error: Missing mesh name argument.\n"
                << "Usage:   " << argv[0] << " <mesh-name>\n"
                << "Example: mpirun -np 4 " << argv[0] << " mesh-cube-10\n";
    }
    return 1;
  }

  // Retrieve the mesh name from the command line
  std::string mesh_name = argv[1];
  
  // Automatically append ".msh" if the user didn't include it
  if (mesh_name.find(".msh") == std::string::npos)
  {
    mesh_name += ".msh";
  }

  // Construct the full path looking directly inside the ../mesh/ folder
  std::string mesh_path = "../mesh/" + mesh_name;

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
  {
    std::cout << "Starting simulation using mesh: " << mesh_path << std::endl;
  }

  // Diffusion coefficient mu
  const auto mu = [](const Point<dim> & /*p*/) { 
    return 1.0; 
  };

  Heat::g_function g_pulsation;
  Heat::h_function h_spatial;

  // Forcing term f
  const auto f = [&g_pulsation, &h_spatial](const Point<dim> &p, const double &t) {
    // Set the current time for g(t)
    g_pulsation.set_time(t);
  
    // Return g(t) * h(x) 
    return g_pulsation.value(p) * h_spatial.value(p);
  };


  ConvergenceTable table;
  std::ofstream convergence_file("convergence.csv");
  convergence_file << "h,eL2,eH1\n";

  const std::vector<unsigned int> N_el_values =
    {5, 10, 20};
    // {5, 10, 20, 40};
    // {10};

  for (const auto &N_el : N_el_values)
  {
    const std::string mesh_file =
      "../mesh/mesh-cube-" +
      std::to_string(N_el) + ".msh";

    Heat problem(/*mesh_filename = */ mesh_file,
               /* degree = */ 2,
               /* T = */ 1.0,
               /* theta = */ 1.0,    
               /* delta_t = */ 0.001,
               mu,    
               f);    

    problem.run();

    const double h = 1.0 / N_el;

    // ExactSolution exact;
		ExactSolution exact(g_pulsation,
                        h_spatial,
                        10);

    exact.set_time(1.0);

    const double error_L2 =
      problem.compute_error(VectorTools::L2_norm,
                            exact);

    const double error_H1 =
      problem.compute_error(VectorTools::H1_norm,
                            exact);

    table.add_value("h", h);
    table.add_value("L2", error_L2);
    table.add_value("H1", error_H1);

    convergence_file << h << ","
                     << error_L2 << ","
                     << error_H1 << "\n";

  // Print final global statistics
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
  {
    std::cout << "\n===============================================" << std::endl;
    std::cout << "GLOBAL EXECUTION STATISTICS" << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << "  MPI Processes:        " << Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) << std::endl;
    std::cout << "  Threads Used:         " << MultithreadInfo::n_threads() << std::endl;
    std::cout << "  Output Files:         " << problem.get_timestep_number() << " VTU+PVTU records" << std::endl;
    std::cout << "===============================================\n" << std::endl;
  }
  }

  table.evaluate_all_convergence_rates(
    ConvergenceTable::reduction_rate_log2);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
	  table.write_text(std::cout);	

  return 0;
}
