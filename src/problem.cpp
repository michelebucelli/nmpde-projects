#include "problem.h"

template<int dim, int fe_degree>
void Problem<dim, fe_degree>::setup(const unsigned int ref_level)
{
    GridGenerator::hyper_cube(mesh);
    mesh.refine_global(ref_level);

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
                        cell->face(f)->set_boundary_id(0); // Left boundary ---> g
                    else if (std::abs(center[0] - 1.0) < 1e-12)
                        cell->face(f)->set_boundary_id(1); // Right boundary --> h
                    else if (std::abs(center[1]) < 1e-12)
                        cell->face(f)->set_boundary_id(1); // Bottom boundary -> h
                    else if (std::abs(center[1] - 1.0) < 1e-12)
                        cell->face(f)->set_boundary_id(1); // Top boundary ----> h

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

    fe = std::make_unique<FE_Q<dim>>(fe_degree);
    quadrature = std::make_unique<QGauss<dim>>(fe_degree + 1);
    quadrature_boundary = std::make_unique<QGauss<dim - 1>>(fe_degree + 1);

    dof_handler.reinit(mesh);
    dof_handler.distribute_dofs(*fe);
    locally_owned_dofs = dof_handler.locally_owned_dofs();

    TrilinosWrappers::SparsityPattern sparsity(locally_owned_dofs, MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler, sparsity);
    sparsity.compress();

    system_matrix.reinit(sparsity);
    system_rhs.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    solution.reinit(locally_owned_dofs, MPI_COMM_WORLD);
}

template<int dim, int fe_degree>
void Problem<dim, fe_degree>::assemble()
{
    const unsigned int dofs_per_cell = fe->dofs_per_cell;
    const unsigned int n_q           = quadrature->size();

    FEValues<dim> fe_values(*fe, *quadrature, update_values | update_gradients | update_quadrature_points | update_JxW_values);
    FEFaceValues<dim> fe_values_boundary(*fe, *quadrature_boundary, update_values | update_quadrature_points | update_JxW_values);

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

    system_matrix = 0.0;
    system_rhs    = 0.0;

    for (const auto &cell : dof_handler.active_cell_iterators())
    {
        if (!cell->is_locally_owned())
            continue;

        fe_values.reinit(cell);

        cell_matrix = 0.0;
        cell_rhs    = 0.0;

        for (unsigned int q = 0; q < n_q; ++q)
        {
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                    // Diffusion term.
                    cell_matrix(i, j) += mu_coefficient.value(fe_values.quadrature_point(q)) *
                                                              fe_values.shape_grad(i, q) *
                                                              fe_values.shape_grad(j, q) *
                                                              fe_values.JxW(q);

                    // Advection term.
                    for (unsigned int d = 0; d < dim; ++d)
                        cell_matrix(i, j) -= beta_coefficient.value(fe_values.quadrature_point(q), d) *
                                                                    fe_values.shape_value(j, q) *
                                                                    fe_values.shape_grad(i, q)[d] *
                                                                    fe_values.JxW(q);


                    // Reaction term.
                    cell_matrix(i, j) += gamma_coefficient.value(fe_values.quadrature_point(q)) *
                                                                 fe_values.shape_value(i, q) *
                                                                 fe_values.shape_value(j, q) *
                                                                 fe_values.JxW(q);
                }

                cell_rhs(i) += source_function.value(fe_values.quadrature_point(q)) * 
                                                     fe_values.shape_value(i, q) * 
                                                     fe_values.JxW(q);
            }
        }

        if (cell->at_boundary())
        {
            for (unsigned int face_number = 0; face_number < cell->n_faces(); ++face_number)
            {
                if (cell->face(face_number)->at_boundary() && (cell->face(face_number)->boundary_id() == 1))
                {
                    fe_values_boundary.reinit(cell, face_number);

                    for (unsigned int q = 0; q < quadrature_boundary->size(); ++q)
                        for (unsigned int i = 0; i < dofs_per_cell; ++i)
                            cell_rhs(i) += h_function.value(fe_values_boundary.quadrature_point(q)) *
                                                            fe_values_boundary.shape_value(i, q) *
                                                            fe_values_boundary.JxW(q);
                }
            }
        }

        cell->get_dof_indices(dof_indices);

        system_matrix.add(dof_indices, cell_matrix);
        system_rhs.add(dof_indices, cell_rhs);
    }

    system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);

    std::map<types::global_dof_index, double> boundary_values;
    std::map<types::boundary_id, const Function<dim> *> boundary_functions;

    boundary_functions[0] = &g_function;

    VectorTools::interpolate_boundary_values(dof_handler, boundary_functions, boundary_values);
    MatrixTools::apply_boundary_values(boundary_values, system_matrix, solution, system_rhs, true);
}

template<int dim, int fe_degree>
void Problem<dim, fe_degree>::solve()
{
    SolverControl solver_control(10000, 1e-12);
    SolverCG<TrilinosWrappers::MPI::Vector> solver(solver_control);
    solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
    pcout << "N-iterations: " << solver_control.last_step() << std::endl;
}

template<int dim, int fe_degree>
void Problem<dim, fe_degree>::output() const
{
    IndexSet locally_relevant_dofs;
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

    TrilinosWrappers::MPI::Vector solution_ghost(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
    solution_ghost = solution;

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution_ghost, "solution");
    data_out.add_data_vector(system_rhs, "rhs");
    data_out.build_patches();

    data_out.write_vtu_with_pvtu_record("./", "solution_classic" + std::to_string(mpi_size), 0, MPI_COMM_WORLD, 5);
}

template class Problem<2, 1>;
template class Problem<3, 1>;
