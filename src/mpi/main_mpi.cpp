#include "WaveEquationMPI.hpp"
#include <deal.II/base/mpi.h>
#include <cmath>
#include <stdexcept>
#include <string>

int main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  try {
    const unsigned int n=argc>1?std::stoul(argv[1]):24;
    const double T=argc>2?std::stod(argv[2]):.1;
    const double dt=argc>3?std::stod(argv[3]):.005;
    const std::string out=argc>4?argv[4]:"output-mpi";
    const std::string test=argc>5?argv[5]:"free";
    if((test!="free"&&test!="forced")||T<=0||dt<=0||std::abs(T/dt-std::round(T/dt))>1.e-10)
      throw std::invalid_argument("Usage: wave_fem_mpi [N] [T] [dt] [output] [free|forced]");
    WaveEquationMPI problem(n,1,T,dt,10,out,test); problem.run(); return 0;
  } catch(const std::exception &error) { std::cerr<<"Error: "<<error.what()<<'\n'; return 1; }
}
