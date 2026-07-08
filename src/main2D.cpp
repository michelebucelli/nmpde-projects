#include "NavierStokes2D.hpp"

namespace
{
  constexpr double H  = 0.41;

  enum class Case2D { case1 = 1, case2 = 2, case3 = 3 };

  NavierStokes2D::ProblemData make_problem_data(const unsigned int choice)
  {
    NavierStokes2D::ProblemData data;

    auto zero_vec = [](const Point<NavierStokes2D::dim> &, const double &) {
      Tensor<1, NavierStokes2D::dim> F;
      F = 0.0;
      return F;
    };

    auto zero_u0 = [](const Point<NavierStokes2D::dim> &) {
      Tensor<1, NavierStokes2D::dim> U;
      U = 0.0;
      return U;
    };

    switch (static_cast<Case2D>(choice))
    {
      case Case2D::case1:
      {
        // 2D-1: steady
        data.choice = choice;
        data.nu = 0.001;  
        data.rho = 1.0;
        data.Um = 0.3;

        data.f = zero_vec;
        data.h = zero_vec;
        data.u0 = zero_u0;
        

        data.u_in = [data](const Point<NavierStokes2D::dim> &p, const double &) {
          Tensor<1, NavierStokes2D::dim> U;
          const double y = p[1];
          const double profile = 4.0 * y * (H - y) / (H * H);
          U[0] = data.Um * profile;
          U[1] = 0.0;
          return U;
        };
        break;
      }

      case Case2D::case2:
      {
        // 2D-2: unsteady
        data.choice = choice;
        data.nu = 0.001;
        data.rho = 1.0;
        data.Um = 1.5;
        
        data.f = zero_vec;
        data.h = zero_vec;
        data.u0 = zero_u0;
        

        data.u_in = [data](const Point<NavierStokes2D::dim> &p, const double &) {
          Tensor<1, NavierStokes2D::dim> U;
          const double y = p[1];
          const double profile = 4.0 * y * (H - y) / (H * H);
          U[0] = data.Um * profile;
          U[1] = 0.0;
          return U;
        };
        break;
      }

      case Case2D::case3:
      {
        // 2D-3: unsteady
        data.choice = choice;
        data.nu = 0.001;
        data.rho = 1.0;
        data.Um = 1.5;

        data.f = zero_vec;
        data.h = zero_vec;
        data.u0 = zero_u0;
        

        data.u_in = [data](const Point<NavierStokes2D::dim> &p, const double &t) {
          Tensor<1, NavierStokes2D::dim> U;
          const double y = p[1];
          const double profile = 4.0 * y * (H - y) / (H * H);
          const double amp = std::sin(M_PI * t / 8.0);
          U[0] = data.Um * amp * profile;
          U[1] = 0.0;
          return U;
        };
        break;
      }

      default:
        throw std::runtime_error("Not valid choice: use 1, 2 or 3.");
    }

    return data;
  }
}


// Main function.
int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  unsigned int choice = 1;

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      std::cout << "Which case you want to test:\n";
      std::cout << "  1 = case 2D-1\n";
      std::cout << "  2 = case 2D-2\n";
      std::cout << "  3 = case 2D-3\n";
      std::cout << "Insert 1, 2 or 3: ";
      std::cin >> choice;
    }

  // Broadcast the choice to all MPI processes
  choice = Utilities::MPI::broadcast(MPI_COMM_WORLD, choice, 0);

  const std::string  mesh_file_name  = "../mesh/mesh-2D-cylinder-circular-cs-0.0205.msh";
  const unsigned int degree_velocity = 2;
  const unsigned int degree_pressure = 1;
  
  const double T = 8.0;
  const double theta = 1.0;
  const double delta_t = 0.001;

  auto fdata = make_problem_data(choice);

  NavierStokes2D problem(mesh_file_name, degree_velocity, degree_pressure, T, theta, delta_t, fdata);

  problem.run();

  return 0;
}