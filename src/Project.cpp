#include <fstream>
#include <iostream>
#include <vector>

#include "Poisson3D.hpp"


int
main(int argc, char * argv[])
{

  Poisson3D problem;
  problem.manage_flags(argc, argv);


  problem.setup();
  problem.assemble();
  problem.solve();
  problem.output();

  return 0;
}