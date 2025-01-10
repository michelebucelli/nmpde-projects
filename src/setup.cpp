#include "setup.h"

template<int dim>
void setup_problem(
    parallel::distributed::Triangulation<dim> &mesh,
    FE_Q<dim>                                 &fe,
    DoFHandler<dim>                           &dof_handler,
    const unsigned int                         level)
{

    GridGenerator::hyper_cube(mesh);
    mesh.refine_global(level);

    for (auto &cell : mesh.active_cell_iterators())
    {
        if (cell->at_boundary())
        {
            for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
            {
                if (cell->face(f)->at_boundary()) 
                {
                    auto center = cell->face(f)->center();

                    if (std::abs(center[0]) < 1e-12)
                        cell->face(f)->set_boundary_id(0);     // Left boundary ---> g
                    else if (std::abs(center[0] - 1.0) < 1e-12)
                        cell->face(f)->set_boundary_id(1);     // Right boundary --> h
                    else if (std::abs(center[1]) < 1e-12)
                        cell->face(f)->set_boundary_id(1);     // Bottom boundary -> h
                    else if (std::abs(center[1] - 1.0) < 1e-12)
                        cell->face(f)->set_boundary_id(1);     // Top boundary ----> h

                    if constexpr (dim == 3)
                    {
                        if (std::abs(center[2]) < 1e-12)
                            cell->face(f)->set_boundary_id(1); // Front boundary --> h
                        else if (std::abs(center[2] - 1.0) < 1e-12)
                            cell->face(f)->set_boundary_id(1); // Back boundary ---> h
                    }
                }
            }
        }
    }

    dof_handler.reinit(mesh);
    dof_handler.distribute_dofs(fe);
}

template void setup_problem<2>(parallel::distributed::Triangulation<2> &, FE_Q<2> &, DoFHandler<2> &, const unsigned int);
template void setup_problem<3>(parallel::distributed::Triangulation<3> &, FE_Q<3> &, DoFHandler<3> &, const unsigned int);