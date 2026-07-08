#include "Current.hpp"
#include "args.hpp"

// Main function.
int main(int argc, char *argv[]){
  // constexpr unsigned int dim = Current::dim;

  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
  Args args(argc, argv);
  
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    std::cout << args << std::endl;

  Current problem(/*mesh_filename = */ args.get_mesh_filename(),
               /* output_file_name = */ args.get_output_file_name(),
               /* degree = */ 1,
               /* solvertype */ args.get_solver(),
               /* T = */ args.get_max_time(),
               /* theta = */ args.get_theta(),
               /* delta_t = */ args.get_delta_t());
  
  problem.set_model(args.get_model());

  problem.run();

  return 0;
}
