#include "NavierStokes.hpp"
#include "NavierStokes.cpp"

// Main function.
int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const std::string mesh_file_name = argc > 1 ? argv[1] : "../mesh/cilinder_3D_fine.msh";
  const unsigned int degree_velocity = 2;
  const unsigned int degree_pressure = 1;

  const double T = 8;
  const double deltat = 0.05;

  dealii::Timer timer;
  // Start the timer
  timer.restart();

  NavierStokes<3> problem(mesh_file_name, degree_velocity, degree_pressure, T, deltat);

  problem.setup();
  problem.solve();

  // Stop the timer
  timer.stop();

  // Output the elapsed time
  if(rank == 0)
    std::cout << "Time taken to solve ENTIRE Navier Stokes problem: " << timer.wall_time() << " seconds" << std::endl;

  if (rank == 0)
  {
    const std::string output_filename = "forces_results_3D_2case.csv";
    std::ofstream outputFile(output_filename);

    if (!outputFile.is_open())
    {
      std::cerr << "Error opening output file" << std::endl;
      return -1;
    }
    outputFile << "Iteration, Drag, Lift, Coeff Drag, CoeffLift, time prec, time solve" << std::endl;

    for (size_t ite = 0; ite < problem.get_result_size(); ite++)
    {
      outputFile << ite * deltat << ", " << problem.get_drag(ite) << ", " << problem.get_lift(ite) << ", " 
                << problem.get_drag_coeff(ite) << ", " << problem.get_lift_coeff(ite) << ", "
                << problem.get_time_prec(ite) << ", " << problem.get_time_solve(ite)
                << std::endl;
    }
    outputFile.close();
  }

  return 0;
}
