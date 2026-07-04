#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/base/mpi.h>

#include <iostream>
#include <fstream>
#include <optional>

#include "ADRFactory.hpp"
#include "TubeSpec.hpp"

void cylindric_wall_adapted_mesh(Triangulation<3> &mesh, double length, double radius, unsigned int length_subdivision, unsigned int global_refinement_steps, unsigned int wall_refinement_steps)
{
    GridGenerator::subdivided_cylinder(mesh, length_subdivision, radius, length/2);
    //to have the inlet around the origin
    //Tensor<1, 3, double> shift_right{{length/2, 0.0, 0.0}};
    //GridTools::shift(shift_right, mesh);

    CylindricalManifold<dim> cylindrical_manifold(0);
    mesh.set_manifold(0, cylindrical_manifold);
    mesh.refine_global(global_refinement_steps);

    for (unsigned int step = 0; step < wall_refinement_steps; ++step)
    {
        for (const auto &cell : mesh.active_cell_iterators())
        {
            if (cell->at_boundary())
            {
                cell->set_refine_flag();
                //if we want specifically to refine the curved walls
                for (const auto &face : cell->face_iterators())
                {
                    if (face->at_boundary() && (face->boundary_id() == 0))
                    {
                        

                        break;
                    }
                }
            }
        }
        mesh.execute_coarsening_and_refinement();
    }
}

int main(int argc, char* argv[])
{
    Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

    double radius = 1.0; //TubeSpec<>::tube_radius;
    double length = radius*28;

    parallel::distributed::Triangulation<dim> mesh(MPI_COMM_WORLD);
    cylindric_wall_adapted_mesh(mesh, length, radius, 100, 2, 3);

    auto adr_mf = ADRFactory::matrix_free<dim, 3, double, TubeSpec, 2>(mesh);
    adr_mf.setup();
    adr_mf.assemble();
    adr_mf.solve(5000, 1e-9, true);
    adr_mf.output("./", "tube_experiment");

}