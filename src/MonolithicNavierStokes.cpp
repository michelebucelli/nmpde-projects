#include "../include/MonolithicNavierStokes.hpp"
#include "../include/preconditioners.hpp"

template <unsigned int dim>
void MonolithicNavierStokes<dim>::setup()
{
    // Create the mesh.
    {
        pcout << "Initializing the mesh" << std::endl;

        Triangulation<dim> mesh_serial;

        GridIn<dim> grid_in;
        grid_in.attach_triangulation(mesh_serial);

        std::ifstream grid_in_file(mesh_file_name);
        grid_in.read_msh(grid_in_file);

        GridTools::partition_triangulation(mpi_size, mesh_serial);
        const auto construction_data = TriangulationDescription::Utilities::
            create_description_from_triangulation(mesh_serial, MPI_COMM_WORLD);
        mesh.create_triangulation(construction_data);
        pcout << "-----------------------------------------------" << std::endl;
        pcout << "Number of elements = " << mesh.n_global_active_cells()
              << std::endl;
    }

    // Initialize the finite element space.
    {

        const FE_SimplexP<dim> fe_scalar_velocity(degree_velocity);
        const FE_SimplexP<dim> fe_scalar_pressure(degree_pressure);

        fe = std::make_unique<FESystem<dim>>(
            fe_scalar_velocity,
            dim,
            fe_scalar_pressure,
            1);

        quadrature = std::make_unique<QGaussSimplex<dim>>(fe->degree + 1);

        quadrature_face = std::make_unique<QGaussSimplex<dim - 1>>(fe->degree + 1);
    }

    // Initialize the DoF handler.
    {
        dof_handler.reinit(mesh);
        dof_handler.distribute_dofs(*fe);

        std::vector<unsigned int> block_component(dim + 1, 0);
        block_component[dim] = 1;
        DoFRenumbering::component_wise(dof_handler, block_component);

        locally_owned_dofs = dof_handler.locally_owned_dofs();
        DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

        std::vector<types::global_dof_index> dofs_per_block =
            DoFTools::count_dofs_per_fe_block(dof_handler, block_component);
        const unsigned int n_u = dofs_per_block[0];
        const unsigned int n_p = dofs_per_block[1];

        block_owned_dofs.resize(2);
        block_relevant_dofs.resize(2);
        block_owned_dofs[0] = locally_owned_dofs.get_view(0, n_u);
        block_owned_dofs[1] = locally_owned_dofs.get_view(n_u, n_u + n_p);
        block_relevant_dofs[0] = locally_relevant_dofs.get_view(0, n_u);
        block_relevant_dofs[1] = locally_relevant_dofs.get_view(n_u, n_u + n_p);

        pcout << "  Number of DoFs: " << std::endl;
        pcout << "    velocity = " << n_u << std::endl;
        pcout << "    pressure = " << n_p << std::endl;
        pcout << "    total    = " << n_u + n_p << std::endl;
    }

    pcout << "-----------------------------------------------" << std::endl;

    // Initialize the linear system.
    {

        Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1);
        for (unsigned int c = 0; c < dim + 1; ++c)
        {
            for (unsigned int d = 0; d < dim + 1; ++d)
            {
                if (c == dim && d == dim)
                    coupling[c][d] = DoFTools::none;
                else
                    coupling[c][d] = DoFTools::always;
            }
        }

        TrilinosWrappers::BlockSparsityPattern sparsity(block_owned_dofs,
                                                        MPI_COMM_WORLD);
        DoFTools::make_sparsity_pattern(dof_handler, coupling, sparsity);
        sparsity.compress();

        for (unsigned int c = 0; c < dim + 1; ++c)
        {
            for (unsigned int d = 0; d < dim + 1; ++d)
            {
                if (c == dim || d == dim)
                    coupling[c][d] = DoFTools::none;
                else
                    coupling[c][d] = DoFTools::always;
            }
        }

        TrilinosWrappers::BlockSparsityPattern velocity_mass_sparsity(
            block_owned_dofs, MPI_COMM_WORLD);
        DoFTools::make_sparsity_pattern(dof_handler, coupling,
                                        velocity_mass_sparsity);
        velocity_mass_sparsity.compress();

        TrilinosWrappers::BlockSparsityPattern pressure_mass_sparsity(
            block_owned_dofs, MPI_COMM_WORLD);
        if (true)
        {
            for (unsigned int c = 0; c < dim + 1; ++c)
            {
                for (unsigned int d = 0; d < dim + 1; ++d)
                {
                    if (c == dim && d == dim)
                        coupling[c][d] = DoFTools::always;
                    else
                        coupling[c][d] = DoFTools::none;
                }
            }

            DoFTools::make_sparsity_pattern(dof_handler, coupling,
                                            pressure_mass_sparsity);
            pressure_mass_sparsity.compress();
        }
        system_matrix.reinit(sparsity);
        velocity_mass.reinit(velocity_mass_sparsity);
        pressure_mass.reinit(pressure_mass_sparsity);
        lhs_matrix.reinit(sparsity);
        system_rhs.reinit(block_owned_dofs, MPI_COMM_WORLD);
        solution_owned.reinit(block_owned_dofs, MPI_COMM_WORLD);
        solution.reinit(block_owned_dofs, block_relevant_dofs, MPI_COMM_WORLD);
    }
}

template <unsigned int dim>
void MonolithicNavierStokes<dim>::assemble_base_matrix()
{
    const unsigned int dofs_per_cell = fe->dofs_per_cell;
    const unsigned int n_q = quadrature->size();

    FEValues<dim> fe_values(*fe,
                            *quadrature,
                            update_values | update_gradients |
                                update_quadrature_points | update_JxW_values);

    FullMatrix<double> cell_system_matrix(dofs_per_cell, dofs_per_cell);

    std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

    FullMatrix<double> velocity_mass_cell_matrix(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> pressure_mass_cell_matrix(dofs_per_cell, dofs_per_cell);

    velocity_mass = 0.0;
    pressure_mass = 0.0;
    system_matrix = 0.0;

    for (const auto &cell : dof_handler.active_cell_iterators())
    {
        if (!cell->is_locally_owned())
            continue;

        fe_values.reinit(cell);

        std::vector<double> div_phi_u(dofs_per_cell);
        std::vector<Tensor<1, dim>> phi_u(dofs_per_cell);
        std::vector<Tensor<2, dim>> grad_phi_u(dofs_per_cell);
        std::vector<double> phi_p(dofs_per_cell);

        pressure_mass_cell_matrix = 0.0;
        velocity_mass_cell_matrix = 0.0;
        cell_system_matrix = 0.0;

        for (unsigned int q = 0; q < n_q; ++q)
        {
            for (unsigned int k = 0; k < dofs_per_cell; k++)
            {
                div_phi_u[k] = fe_values[velocity].divergence(k, q);
                grad_phi_u[k] = fe_values[velocity].gradient(k, q);
                phi_u[k] = fe_values[velocity].value(k, q);
                phi_p[k] = fe_values[pressure].value(k, q);
            }

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                    // Mass Component
                    // ------
                    // M_ij = ∫ φ_i·φ_j dx
                    // ------
                    cell_system_matrix(i, j) += scalar_product(fe_values[velocity].value(i, q), fe_values[velocity].value(j, q)) / deltat * fe_values.JxW(q);

                    // Stiffness Component
                    // ------
                    // A_ij = ∫ ν ∇φ_i:∇φ_j dx
                    // ------
                    cell_system_matrix(i, j) += nu * scalar_product(fe_values[velocity].gradient(i, q), fe_values[velocity].gradient(j, q)) * fe_values.JxW(q);

                    // Pressure term in the momentum equation
                    // ------
                    // B_ij = -∫ ψ_j ∇·φ_i dx
                    // ------
                    cell_system_matrix(i, j) -= fe_values[pressure].value(j, q) * fe_values[velocity].divergence(i, q) * fe_values.JxW(q);

                    // Pressure term in the continuity equation
                    // ------
                    // B_ij^T = ∫ ψ_i ∇·φ_j dx
                    // ------
                    cell_system_matrix(i, j) += fe_values[pressure].value(i, q) * fe_values[velocity].divergence(j, q) * fe_values.JxW(q);

                    // Mass matrix for the pressure
                    // ------
                    // Mp_ij = ∫ (1/ν)ψ_i ψ_j dx
                    // ------
                    pressure_mass_cell_matrix(i, j) += fe_values[pressure].value(j, q) * fe_values[pressure].value(i, q) / nu * fe_values.JxW(q);

                    // Mass matrix for the velocity
                    // ------
                    // Mv_ij = ∫ (1/Δt) φ_i·φ_j dx
                    // ------
                    velocity_mass_cell_matrix(i, j) += scalar_product(fe_values[velocity].value(i, q), fe_values[velocity].value(j, q)) / deltat * fe_values.JxW(q);
                }
            }
        }
        cell->get_dof_indices(dof_indices);

        system_matrix.add(dof_indices, cell_system_matrix);
        velocity_mass.add(dof_indices, velocity_mass_cell_matrix);
        pressure_mass.add(dof_indices, pressure_mass_cell_matrix);
    }
    system_matrix.compress(VectorOperation::add);
    pressure_mass.compress(VectorOperation::add);
    velocity_mass.compress(VectorOperation::add);
}

template <unsigned int dim>
void MonolithicNavierStokes<dim>::add_convective_term()
{
    const unsigned int dofs_per_cell = fe->dofs_per_cell;
    const unsigned int n_q = quadrature->size();

    FEValues<dim> fe_values(*fe,
                            *quadrature,
                            update_values | update_gradients |
                                update_quadrature_points | update_JxW_values);

    FullMatrix<double> cell_lhs_matrix(dofs_per_cell, dofs_per_cell);

    std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

    std::vector<Tensor<1, dim>> previous_velocity_values(n_q);
    std::vector<double> previous_velocity_divergence(n_q);

    lhs_matrix = 0.0;

    lhs_matrix.copy_from(system_matrix);

    for (const auto &cell : dof_handler.active_cell_iterators())
    {
        if (!cell->is_locally_owned())
            continue;

        fe_values.reinit(cell);

        cell_lhs_matrix = 0.0;

        fe_values[velocity].get_function_values(solution, previous_velocity_values);
        fe_values[velocity].get_function_divergences(solution, previous_velocity_divergence);

        for (unsigned int q = 0; q < n_q; ++q)
        {
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                    // Non-linear term
                    // ------
                    // N_ij(u*) = ∫ (u*·∇u)·φ_j dx
                    // ------
                    cell_lhs_matrix(i, j) += scalar_product(fe_values[velocity].gradient(j, q) * previous_velocity_values[q], fe_values[velocity].value(i, q)) * fe_values.JxW(q);

                    // skew-symmetric term
                    // ------
                    // S_ij = (1/2) ∫ ∇·u* φ_i·φ_j dx
                    // ------
                    cell_lhs_matrix(i, j) += 0.5 * previous_velocity_divergence[q] * scalar_product(fe_values[velocity].value(i, q), fe_values[velocity].value(j, q)) * fe_values.JxW(q);
                }
            }
        }

        cell->get_dof_indices(dof_indices);

        lhs_matrix.add(dof_indices, cell_lhs_matrix);
    }
    lhs_matrix.compress(VectorOperation::add);
}

template <unsigned int dim>
void MonolithicNavierStokes<dim>::assemble_rhs()
{
    const unsigned int dofs_per_cell = fe->dofs_per_cell;
    const unsigned int n_q = quadrature->size();
    const unsigned int n_q_face = quadrature_face->size();

    FEValues<dim> fe_values(*fe,
                            *quadrature,
                            update_values | update_gradients | update_quadrature_points |
                                update_JxW_values);

    FEFaceValues<dim> fe_face_values(*fe,
                                     *quadrature_face,
                                     update_values |
                                         update_normal_vectors |
                                         update_JxW_values);
    Vector<double> cell_rhs(dofs_per_cell);

    std::vector<Tensor<1, dim>> previous_velocity_values(n_q);

    std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

    Vector<double> f_neumann_loc(dim + 1);

    Tensor<1, dim> f_neumann_tensor;

    system_rhs = 0.0;

    for (const auto &cell : dof_handler.active_cell_iterators())
    {
        if (!cell->is_locally_owned())
            continue;

        fe_values.reinit(cell);

        fe_values[velocity].get_function_values(solution, previous_velocity_values);

        cell_rhs = 0.0;

        for (unsigned int q = 0; q < n_q; ++q)
        {
            // Compute f(tn+1)
            Vector<double> f_new_loc(dim);
            forcing_term.set_time(time);
            forcing_term.vector_value(fe_values.quadrature_point(q),
                                      f_new_loc);
            Tensor<1, dim> forcing_term_new_tensor;
            for (unsigned int d = 0; d < dim; ++d)
                forcing_term_new_tensor[d] = f_new_loc[d];

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
                // Time dependent term
                // ------
                // ∫ (1/Δt)(u*·φ_i) dx
                // ------
                cell_rhs(i) += scalar_product(previous_velocity_values[q], fe_values[velocity].value(i, q)) * fe_values.JxW(q) / deltat;

                // Forcing Term
                // ------
                // ∫ f(t+1)·φ_i dx
                // ------
                if (!zero_forcing)
                    cell_rhs(i) += scalar_product(forcing_term_new_tensor, fe_values[velocity].value(i, q)) * fe_values.JxW(q);
            }
        }

        if (cell->at_boundary())
        {
            for (unsigned int f = 0; f < cell->n_faces(); ++f)
            {
                if (cell->face(f)->at_boundary() &&
                    cell->face(f)->boundary_id() == 1)
                {

                    fe_face_values.reinit(cell, f);

                    for (unsigned int q = 0; q < n_q_face; ++q)
                    {
                        neumann_function.vector_value(
                            fe_face_values.quadrature_point(q), f_neumann_loc);

                        for (unsigned int d = 0; d < dim; ++d)
                            f_neumann_tensor[d] = f_neumann_loc(d);

                        for (unsigned int i = 0; i < dofs_per_cell; ++i)
                        {
                            // Neumann boundary term
                            // ------
                            // ∫ f_neumann·φ_i dx
                            // ------
                            if (!zero_neumann)
                                cell_rhs(i) += scalar_product(f_neumann_tensor, fe_face_values[velocity].value(i, q)) * fe_face_values.JxW(q);
                        }
                    }
                }
            }
        }
        cell->get_dof_indices(dof_indices);
        system_rhs.add(dof_indices, cell_rhs);
    }

    system_rhs.compress(VectorOperation::add);

    // We apply boundary conditions to the algebraic system.
    {
        std::map<types::global_dof_index, double> boundary_values;
        std::map<types::boundary_id, const Function<dim> *> boundary_functions;

        ComponentMask mask;

        static_assert(dim == 2 || dim == 3,
                      "Dimensions other than 2 or 3 are not supported");

        if constexpr (dim == 2)
            mask = ComponentMask({true, true, false});
        else if constexpr (dim == 3)
            mask = ComponentMask({true, true, true, false});

        boundary_functions[0] = &inlet_velocity;
        VectorTools::interpolate_boundary_values(dof_handler,
                                                 boundary_functions,
                                                 boundary_values,
                                                 mask);

        boundary_functions.clear();
        Functions::ZeroFunction<dim> zero_function(dim + 1);
        boundary_functions[2] = &zero_function;
        boundary_functions[3] = &zero_function;
        VectorTools::interpolate_boundary_values(dof_handler,
                                                 boundary_functions,
                                                 boundary_values,
                                                 mask);

        MatrixTools::apply_boundary_values(
            boundary_values, lhs_matrix, solution_owned, system_rhs, false);
    }
}

template <unsigned int dim>
void MonolithicNavierStokes<dim>::solve_time_step()
{
    // Choose the preconditioner type:
    // 1 = SIMPLE, 2 = ASIMPLE, 3 = YOSIDA, 4 = AYOSIDA.
    static constexpr int precond_type = 1; // Set this value as needed

    // Local parameters for inner solvers and preconditioner initialization.
    static constexpr double alpha = 1;                   // Damping parameter for SIMPLE-like preconditioners
    static constexpr unsigned int maxiter_inner = 10000; // Maximum iterations for inner solvers
    static constexpr double tol_inner = 1e-5;            // Tolerance for inner solvers
    static constexpr bool use_ilu = false;               // Flag: true for ILU, false for AMG

    // Set up the outer GMRES solver.
    SolverControl solver_control(10000, 1e-7);
    SolverGMRES<TrilinosWrappers::MPI::BlockVector> solver(solver_control);

    std::shared_ptr<BlockPrecondition> block_precondition;

    // Select and initialize the preconditioner based on precond_type.
    switch (precond_type)
    {
    case 1:
    {
        auto simple_precondition = std::make_shared<PreconditionSIMPLE>();
        simple_precondition->initialize(
            lhs_matrix.block(0, 0),
            lhs_matrix.block(1, 0),
            lhs_matrix.block(0, 1),
            solution_owned,
            alpha,
            maxiter_inner,
            tol_inner,
            use_ilu);
        block_precondition = simple_precondition;
        break;
    }
    case 2:
    {
        auto asimple_precondition = std::make_shared<PreconditionaSIMPLE>();
        asimple_precondition->initialize(
            lhs_matrix.block(0, 0),
            lhs_matrix.block(1, 0),
            lhs_matrix.block(0, 1),
            solution_owned,
            alpha,
            maxiter_inner,
            tol_inner,
            use_ilu);
        block_precondition = asimple_precondition;
        break;
    }
    case 3:
    {
        auto yosida_precondition = std::make_shared<PreconditionYosida>();
        yosida_precondition->initialize(
            lhs_matrix.block(0, 0),
            lhs_matrix.block(1, 0),
            lhs_matrix.block(0, 1),
            velocity_mass.block(0, 0),
            solution_owned,
            maxiter_inner,
            tol_inner,
            use_ilu);
        block_precondition = yosida_precondition;
        break;
    }
    default:
        Assert(false, ExcNotImplemented());
    }

    solver.solve(lhs_matrix,
                 solution_owned,
                 system_rhs,
                 *block_precondition);

    pcout << "  " << solver_control.last_step() << " GMRES iterations" << std::endl;

    solution = solution_owned;

    solution.update_ghost_values();
}

template <unsigned int dim>
void MonolithicNavierStokes<dim>::solve()
{
    time = 0.0;

    VectorTools::interpolate(dof_handler, initial_condition, solution_owned);
    solution = solution_owned;

    unsigned int time_step = 0;

    assemble_base_matrix();

    while (time < T - 0.5 * deltat)
    {
        time += deltat;
        ++time_step;

        pcout << "n = " << std::setw(3) << time_step << ", t = " << std::setw(5)
              << time << ":" << std::flush;

        add_convective_term();
        assemble_rhs();
        solve_time_step();
        output(time_step);
    }
}

template <unsigned int dim>
void MonolithicNavierStokes<dim>::output(const unsigned int &time_step)
{
    DataOut<dim> data_out;

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
        data_component_interpretation(
            dim, DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation.push_back(
        DataComponentInterpretation::component_is_scalar);

    std::vector<std::string> names(dim, "velocity");
    names.push_back("pressure");

    data_out.add_data_vector(dof_handler, solution, names, data_component_interpretation);

    std::vector<unsigned int> partition_int(mesh.n_active_cells());
    GridTools::get_subdomain_association(mesh, partition_int);
    const Vector<float> partitioning(partition_int.begin(), partition_int.end());
    data_out.add_data_vector(partitioning, "partitioning");

    data_out.build_patches();

    std::string output_dir = get_output_directory();

    data_out.write_vtu_with_pvtu_record(
        output_dir, "output_", time_step, MPI_COMM_WORLD, 3);
}

template <unsigned int dim>
void MonolithicNavierStokes<dim>::run()
{
    setup();
    solve();
}

template <unsigned int dim>
std::string MonolithicNavierStokes<dim>::get_output_directory() const
{
    namespace fs = std::filesystem;

    if (!fs::exists("./outputs"))
    {
        fs::create_directory("./outputs");
    }
    if constexpr (dim == 2)
    {
        if (!fs::exists("./outputs/monolithicNavierStokes2D"))
        {
            fs::create_directory("./outputs/monolithicNavierStokes2D");
        }
    }
    else if constexpr (dim == 3)
    {
        if (!fs::exists("./outputs/monolithicNavierStokes3D"))
        {
            fs::create_directory("./outputs/monolithicNavierStokes3D");
        }
    }

    std::string sub_dir_name = "outputs_reynolds_" + std::to_string(static_cast<int>(reynolds_number));
    fs::path sub_dir_path = "";
    if constexpr (dim == 2)
        sub_dir_path = "./outputs/monolithicNavierStokes2D/" + sub_dir_name + "/";
    if constexpr (dim == 3)
        sub_dir_path = "./outputs/monolithicNavierStokes3D/" + sub_dir_name + "/";
    if (!fs::exists(sub_dir_path))
    {
        fs::create_directory(sub_dir_path);
    }

    return sub_dir_path.string();
}

template class MonolithicNavierStokes<2>;
template class MonolithicNavierStokes<3>;