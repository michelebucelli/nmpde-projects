#include "FK_solver.hpp"
#include <deal.II/base/convergence_table.h>


// Main function.
int
main(int argc, char *argv[])
{
    Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
    const unsigned int               mpi_rank =
    Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  const unsigned int degree = 2;

   double T      = 3.0;
   double deltat = 0.1;
   FK_solver problem("../mesh/brain-h3.0.msh", degree, T, deltat);
   problem.setup();
   problem.solve();


  //     ConvergenceTable table;

  // const std::vector<std::string> meshes = {"../mesh/mesh-square-h0.100000.msh",
  //                                          "../mesh/mesh-square-h0.050000.msh",
  //                                          "../mesh/mesh-square-h0.025000.msh",
  //                                          "../mesh/mesh-square-h0.012500.msh"};
  // const std::vector<double>      h_vals = {0.1,
  //                                          0.05,
  //                                          0.025,
  //                                          0.0125};

  // std::ofstream convergence_file("convergence.csv");
  // convergence_file << "h,eL2,eH1" << std::endl;

  // T = 0.00000001;
  // deltat = 0.000000001;

  // for (unsigned int i = 0; i < meshes.size(); ++i)
  //   {
  
  // FK_solver problem(meshes[i], degree, T, deltat);

  //     problem.setup();
  //     problem.solve();

  //     const double error_L2 = problem.compute_error(VectorTools::L2_norm);
  //     const double error_H1 = problem.compute_error(VectorTools::H1_norm);

  //     table.add_value("h", h_vals[i]);
  //     table.add_value("L2", error_L2);
  //     table.add_value("H1", error_H1);

  //     convergence_file << h_vals[i] << "," << error_L2 << "," << error_H1
  //                      << std::endl;
  //   }
  

  //  if (mpi_rank == 0)
  //   {
  //       table.evaluate_all_convergence_rates(ConvergenceTable::reduction_rate_log2);
  //       table.set_scientific("L2", true);
  //       table.set_scientific("H1", true);
  //       table.write_text(std::cout);
  //   }

  return 0;
}
