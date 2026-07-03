#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

using namespace dealii;

template <int dim>
void generate_and_output_grid(unsigned int N)
{
    Triangulation<dim> triangulation;
    GridGenerator::subdivided_hyper_cube(triangulation, N, 0.0, 1.0, true);
    GridTools::distort_random(0.3, triangulation, true);

    std::string filename = "distorted_hyper_cube-" + std::to_string(dim) + "-" + std::to_string(N) + ".msh";
    GridOut grid_out;
    grid_out.write_msh(triangulation, filename);
    std::cout << "Successfully wrote mesh to: " << filename << std::endl;
}


int main(int argc, char* argv[])
{
    unsigned int dim = atoi(argv[1]);
    unsigned int N = atoi(argv[2]);
    

    if(dim == 2) generate_and_output_grid<2>(N);
    if(dim == 3) generate_and_output_grid<3>(N);
}