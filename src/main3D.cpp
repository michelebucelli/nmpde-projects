#include "NavierStokes3D.hpp"

namespace
{
  constexpr double H  = 0.41;

  enum class Case3D { case1 = 1, case2 = 2, case3 = 3 };

  NavierStokes3D::ProblemData make_problem_data(const unsigned int choice)
  {
    NavierStokes3D::ProblemData data;

    auto zero_vec = [](const Point<NavierStokes3D::dim> &, const double &) {
      Tensor<1, NavierStokes3D::dim> F;
      F = 0.0;
      return F;
    };

    auto zero_u0 = [](const Point<NavierStokes3D::dim> &) {
      Tensor<1, NavierStokes3D::dim> U;
      U = 0.0;
      return U;
    };

    switch (static_cast<Case3D>(choice))
    {
      case Case3D::case1:
      {
        // 3D-1: steady
        data.choice = choice;
        data.nu = 0.001;  
        data.rho = 1.0;
        data.Um = 0.45;

        data.f = zero_vec;
        data.h = zero_vec;
        data.u0 = zero_u0;
        

        data.u_in = [data](const Point<NavierStokes3D::dim> &p, const double &) {
          Tensor<1, NavierStokes3D::dim> U;
          const double y = p[1], z = p[2];
          const double profile = 16.0 * y * z * (H - y) * (H - z) / (H * H * H * H);
          U[0] = data.Um * profile;
          U[1] = 0.0;
          U[2] = 0.0;
          return U;
        };
        break;
      }

      case Case3D::case2:
      {
        // 3D-2: unsteady
        data.choice = choice;
        data.nu = 0.001;
        data.rho = 1.0;
        data.Um = 2.25;
        
        data.f = zero_vec;
        data.h = zero_vec;
        data.u0 = zero_u0;
        

        data.u_in = [data](const Point<NavierStokes3D::dim> &p, const double &) {
          Tensor<1, NavierStokes3D::dim> U;
          const double y = p[1], z = p[2];
          const double profile = 16.0 * y * z * (H - y) * (H - z) / (H * H * H * H);
          U[0] = data.Um * profile;
          U[1] = 0.0;
          U[2] = 0.0;
          return U;
        };
        break;
      }

      case Case3D::case3:
      {
        // 3D-3: unsteady
        data.choice = choice;
        data.nu = 0.001;
        data.rho = 1.0;
        data.Um = 2.25;

        data.f = zero_vec;
        data.h = zero_vec;
        data.u0 = zero_u0;
        

        data.u_in = [data](const Point<NavierStokes3D::dim> &p, const double &t) {
          Tensor<1, NavierStokes3D::dim> U;
          const double y = p[1], z = p[2];
          const double profile = 16.0 * y * z * (H - y) * (H - z) / (H * H * H * H);
          const double amp = std::sin(M_PI * t / 8.0);
          U[0] = data.Um * amp * profile;
          U[1] = 0.0;
          U[2] = 0.0;
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
      std::cout << "  1 = case 3D-1Z\n";
      std::cout << "  2 = case 3D-2Z\n";
      std::cout << "  3 = case 3D-3Z\n";
      std::cout << "Insert 1, 2 or 3: ";
      std::cin >> choice;
    }

  // Broadcast the choice to all MPI processes
  choice = Utilities::MPI::broadcast(MPI_COMM_WORLD, choice, 0);

  const std::string  mesh_file_name  = "../mesh/mesh-3D-cylinder-circular-cs-0.0410.msh";
  const unsigned int degree_velocity = 2;
  const unsigned int degree_pressure = 1;
  
  const double T = 8.0;
  const double theta = 1.0;
  const double delta_t = 0.01;

  auto fdata = make_problem_data(choice);

  NavierStokes3D problem(mesh_file_name, degree_velocity, degree_pressure, T, theta, delta_t, fdata);

  problem.run();

  return 0;
}