#include "../include/UncoupledNavierStokes.hpp"

template <unsigned int dim>
void UncoupledNavierStokes<dim>::setup()
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
    pcout << "  Number of elements = " << mesh.n_global_active_cells() << std::endl;

    dof_handler_velocity.reinit(mesh);
    dof_handler_velocity.distribute_dofs(fe_velocity);
    locally_owned_velocity = dof_handler_velocity.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs(dof_handler_velocity, locally_relevant_velocity);

    constraints_velocity.clear();
    constraints_velocity.reinit(locally_relevant_velocity);
    DoFTools::make_hanging_node_constraints(dof_handler_velocity, constraints_velocity);

    VectorTools::interpolate_boundary_values(dof_handler_velocity,
                                             /*boundary_id=*/0,
                                             InletVelocity(H),
                                             constraints_velocity);

    VectorTools::interpolate_boundary_values(dof_handler_velocity,
                                             /*boundary_id=*/2,
                                             Functions::ZeroFunction<dim>(dim),
                                             constraints_velocity);

    VectorTools::interpolate_boundary_values(dof_handler_velocity,
                                             /*boundary_id=*/3,
                                             Functions::ZeroFunction<dim>(dim),
                                             constraints_velocity);

    constraints_velocity.close();

    //-----------------------------
    // Pressure dofs
    //-----------------------------
    dof_handler_pressure.reinit(mesh);
    dof_handler_pressure.distribute_dofs(fe_pressure);
    locally_owned_pressure = dof_handler_pressure.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs(dof_handler_pressure, locally_relevant_pressure);

    constraints_pressure.clear();
    constraints_pressure.reinit(locally_relevant_pressure);
    DoFTools::make_hanging_node_constraints(dof_handler_pressure, constraints_pressure);

    // Homogeneus dirichlet boundary conditions for the pressure since there is neumann on the velocity
    VectorTools::interpolate_boundary_values(
        dof_handler_pressure,
        1,
        Functions::ZeroFunction<dim>(1), 
        constraints_pressure);
    constraints_pressure.close();

    TrilinosWrappers::SparsityPattern dsp_v(locally_owned_velocity,
                                            MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler_velocity, dsp_v);
    dsp_v.compress();

    velocity_matrix.reinit(dsp_v);
    velocity_update_matrix.reinit(dsp_v);

    TrilinosWrappers::SparsityPattern dsp_p(locally_owned_pressure,
                                            MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler_pressure, dsp_p);
    dsp_p.compress();

    pressure_matrix.reinit(dsp_p);

    old_velocity.reinit(locally_owned_velocity, locally_relevant_velocity, MPI_COMM_WORLD);
    old_old_velocity.reinit(locally_owned_velocity, locally_relevant_velocity, MPI_COMM_WORLD);
    velocity_solution.reinit(locally_owned_velocity, locally_relevant_velocity, MPI_COMM_WORLD);
    update_velocity_solution.reinit(locally_owned_velocity, locally_relevant_velocity, MPI_COMM_WORLD);
    velocity_system_rhs.reinit(locally_owned_velocity, MPI_COMM_WORLD);
    velocity_update_rhs.reinit(locally_owned_velocity, MPI_COMM_WORLD);

    // old_pressure.reinit(locally_owned_pressure, locally_relevant_pressure, MPI_COMM_WORLD);
    deltap.reinit(locally_owned_pressure, locally_relevant_pressure, MPI_COMM_WORLD);
    pressure_solution.reinit(locally_owned_pressure, locally_relevant_pressure, MPI_COMM_WORLD);
    pressure_system_rhs.reinit(locally_owned_pressure, MPI_COMM_WORLD);

    pcout << "  Number of DoFs: " << std::endl;
    pcout << "    velocity = " << dof_handler_velocity.n_dofs() << std::endl;
    pcout << "    pressure = " << dof_handler_pressure.n_dofs() << std::endl;
    pcout << "    total    = " << dof_handler_velocity.n_dofs() + dof_handler_pressure.n_dofs() << std::endl;
    pcout << "-----------------------------------------------" << std::endl;
}

template <unsigned int dim>
void UncoupledNavierStokes<dim>::assemble_system_velocity()
{
    TimerOutput::Scope t(computing_timer, "assemble_velocity");

    velocity_matrix = 0;
    velocity_system_rhs = 0;

    const unsigned int quad_deg = std::max<unsigned int>(2u, fe_velocity.degree + 1u);
    QGaussSimplex<dim> quadrature_formula(quad_deg);

    FEValues<dim> fe_values(fe_velocity, quadrature_formula,
                            update_values | update_gradients |
                                update_quadrature_points | update_JxW_values);

    FEValues<dim> fe_values_pressure(fe_pressure, quadrature_formula,
                                     update_gradients | update_quadrature_points);

    const unsigned int dofs_per_cell = fe_velocity.dofs_per_cell;
    const unsigned int n_q = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs(dofs_per_cell);
    std::vector<types::global_dof_index> local_indices(dofs_per_cell);

    std::vector<Tensor<1, dim>> old_val(n_q);
    std::vector<double> old_div(n_q);
    std::vector<Tensor<1, dim>> old_old_val(n_q);
    std::vector<double> old_old_div(n_q);
    std::vector<Tensor<1, dim>> pressure_grad(n_q);

    auto cell_v = dof_handler_velocity.begin_active();
    auto cell_p = dof_handler_pressure.begin_active();
    const auto end_v = dof_handler_velocity.end();

    for (; cell_v != end_v; ++cell_v, ++cell_p)
    {
        if (!cell_v->is_locally_owned())
            continue;

        fe_values.reinit(cell_v);
        fe_values_pressure.reinit(cell_p);

        cell_matrix = 0;
        cell_rhs = 0;

        const auto &vel_extract = fe_values[FEValuesExtractors::Vector(0)];
        vel_extract.get_function_values(old_velocity, old_val);
        vel_extract.get_function_divergences(old_velocity, old_div);
        vel_extract.get_function_values(old_old_velocity, old_old_val);
        vel_extract.get_function_divergences(old_old_velocity, old_old_div);

        fe_values_pressure.get_function_gradients(pressure_solution, pressure_grad);

        for (unsigned int q = 0; q < n_q; ++q)
        {
            // ------
            // u* = 2*u^n - u^{n-1}
            // ------
            const Tensor<1, dim> &u_star = 2.0 * old_val[q] - old_old_val[q]; 

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {

                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                    double lhs = 0.0;
                    
                    // Mass Term
                    // ------
                    // M_ij = (3/2) * (1/Δt) ∫ φ_i·φ_j dx
                    // ------
                    lhs += (3.0/2.0) * (1./deltat) * scalar_product(vel_extract.value(j, q), vel_extract.value(i, q));
                    
                    // Viscous
                    // ------
                    // A_ij = ∫ ν ∇φ_i:∇φ_j dx
                    // ------
                    lhs += nu * scalar_product(vel_extract.gradient(j, q), vel_extract.gradient(i, q));
                    
                    // Convection
                    // ------
                    // N_ij(u*) = ∫ (u* φ_j)·φ_i dx
                    // ------
                    lhs += (scalar_product(vel_extract.gradient(j, q) * u_star, vel_extract.value(i, q))); 

                    cell_matrix(i, j) += lhs * fe_values.JxW(q);
                }

                // Time dependent term
                // ------
                // ∫ (1/2Δt)(4u^n·φ_i- u^{n-1}·φ_i) dx
                // ------
                cell_rhs(i) += (1.0/(2*deltat)) * scalar_product(4.0 * old_val[q], vel_extract.value(i, q)) * fe_values.JxW(q);
                cell_rhs(i) -= (1.0/(2*deltat)) * scalar_product(old_old_val[q], vel_extract.value(i, q)) * fe_values.JxW(q);
                
                // Pressure term
                // ------
                // - ∫ φ_j·∇p dx
                cell_rhs(i) -= scalar_product(pressure_grad[q], vel_extract.value(i, q)) * fe_values.JxW(q);
            }
        }
        cell_v->get_dof_indices(local_indices);
        constraints_velocity.distribute_local_to_global(cell_matrix, cell_rhs,
                                                        local_indices,
                                                        velocity_matrix,
                                                        velocity_system_rhs);
    }

    velocity_matrix.compress(VectorOperation::add);
    velocity_system_rhs.compress(VectorOperation::add);
}

template <unsigned int dim>
void UncoupledNavierStokes<dim>::solve_velocity_system()
{
    TimerOutput::Scope t(computing_timer, "solve_velocity");

    TrilinosWrappers::MPI::Vector tmp(locally_owned_velocity, MPI_COMM_WORLD);

    SolverControl solver_control(1000000, 1e-7 * velocity_system_rhs.l2_norm());

    // Create and initialize preconditioner:
    TrilinosWrappers::PreconditionSSOR prec;
    prec.initialize(velocity_matrix);

    // Create GMRES solver *without* specifying any restart parameter:
    SolverGMRES<TrilinosWrappers::MPI::Vector> solver_gmres(solver_control);

    // Solve the linear system:
    solver_gmres.solve(velocity_matrix, tmp, velocity_system_rhs, prec);

    if (mpi_rank == 0)
        std::cout << "Velocity GMRES iterations: " << solver_control.last_step() << std::endl;

    // Distribute constraints (apply hanging-node constraints, Dirichlet BC, etc.):
    constraints_velocity.distribute(tmp);

    // Update the global velocity solution:
    velocity_solution = tmp;

    constraints_velocity.distribute(velocity_solution);

    velocity_solution.update_ghost_values();
}
template <unsigned int dim>
void UncoupledNavierStokes<dim>::assemble_system_pressure()
{

    TimerOutput::Scope t(computing_timer, "assemble_pressure");
    pressure_matrix = 0;
    pressure_system_rhs = 0;

    const unsigned int quad_deg = std::max<unsigned int>(2u, fe_pressure.degree + 1u);
    QGaussSimplex<dim> quad(quad_deg);
    QGaussSimplex<dim - 1> quad_face(quad_deg);

    FEValues<dim> fe_values_p(fe_pressure, quad,
                              update_values | update_gradients |
                                  update_quadrature_points | update_JxW_values);

    FEValues<dim> fe_values_v(fe_velocity, quad,
                              update_values | update_gradients | update_quadrature_points);

    FEFaceValues<dim> fe_face_values(fe_pressure, quad_face,
                                     update_values | update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe_pressure.dofs_per_cell;
    const unsigned int n_q = quad.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs(dofs_per_cell);
    std::vector<types::global_dof_index> local_indices(dofs_per_cell);

    std::vector<double> div_u_tilde(n_q);

    auto cell_p = dof_handler_pressure.begin_active();
    auto cell_v = dof_handler_velocity.begin_active();
    const auto end_p = dof_handler_pressure.end();

    for (; cell_p != end_p; ++cell_p, ++cell_v)
    {
        if (!cell_p->is_locally_owned())
            continue;
        fe_values_p.reinit(cell_p);
        fe_values_v.reinit(cell_v);

        cell_matrix = 0;
        cell_rhs = 0;

        const auto &vel_extract = fe_values_v[FEValuesExtractors::Vector(0)];
        vel_extract.get_function_divergences(velocity_solution, div_u_tilde);

        for (unsigned int q = 0; q < n_q; ++q)
        {
            const double JxW = fe_values_p.JxW(q);
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {

                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                    // LSH
                    // ------
                    // L_ij = ∫ ∇ψ_i·∇ψ_j dx
                    // ------
                    cell_matrix(i, j) += scalar_product(fe_values_p.shape_grad(j, q), fe_values_p.shape_grad(i, q)) * JxW;
                }
                
                // RHS
                // ------
                // - (3/2) * (1/Δt) ∫ div(u~)ψ_i dx
                cell_rhs(i) -= 3.0 / (2.0 * deltat) * (div_u_tilde[q] * fe_values_p.shape_value(i, q)) * JxW;
            }
        }

        cell_p->get_dof_indices(local_indices);
        constraints_pressure.distribute_local_to_global(cell_matrix, cell_rhs,
                                                        local_indices,
                                                        pressure_matrix,
                                                        pressure_system_rhs);
    }

    pressure_matrix.compress(VectorOperation::add);
    pressure_system_rhs.compress(VectorOperation::add);
}

template <unsigned int dim>
void UncoupledNavierStokes<dim>::solve_pressure_system()
{
    TimerOutput::Scope t(computing_timer, "solve_pressure");

    TrilinosWrappers::MPI::Vector tmp(locally_owned_pressure, MPI_COMM_WORLD);
    SolverControl solver_control(2000000, 1e-7 * pressure_system_rhs.l2_norm());

    TrilinosWrappers::PreconditionIC prec;
    prec.initialize(pressure_matrix);

    SolverCG<TrilinosWrappers::MPI::Vector> solver_cg(solver_control);
    solver_cg.solve(pressure_matrix, tmp, pressure_system_rhs, prec);

    if (mpi_rank == 0)
        std::cout << "Pressure CG iterations: " << solver_control.last_step() << std::endl;

    constraints_pressure.distribute(tmp);
    deltap = tmp;
    constraints_pressure.distribute(deltap);

    deltap.update_ghost_values();
}

template <unsigned int dim>
void UncoupledNavierStokes<dim>::update_velocity()
{
    TimerOutput::Scope t(computing_timer, "assemble_update");

    velocity_update_matrix = 0;
    velocity_update_rhs = 0;

    const unsigned int quad_deg = std::max<unsigned int>(2u, fe_velocity.degree + 1u);
    QGaussSimplex<dim> quad(quad_deg);

    FEValues<dim> fe_values_vel(fe_velocity, quad,
                                update_values | update_gradients |
                                    update_quadrature_points | update_JxW_values);

    FEValues<dim> fe_values_p(fe_pressure, quad,
                              update_gradients | update_quadrature_points);

    const unsigned int dofs_per_cell = fe_velocity.dofs_per_cell;
    const unsigned int n_q = quad.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs(dofs_per_cell);
    std::vector<types::global_dof_index> local_indices(dofs_per_cell);

    std::vector<Tensor<1, dim>> u_tilde_vals(n_q);
    std::vector<Tensor<1, dim>> grad_delta_p(n_q);

    auto cell_v = dof_handler_velocity.begin_active();
    auto cell_p = dof_handler_pressure.begin_active();
    const auto end_v = dof_handler_velocity.end();

    for (; cell_v != end_v; ++cell_v, ++cell_p)
    {
        if (!cell_v->is_locally_owned())
            continue;

        fe_values_vel.reinit(cell_v);
        fe_values_p.reinit(cell_p);

        cell_matrix = 0;
        cell_rhs = 0;

        const auto &vel_extract = fe_values_vel[FEValuesExtractors::Vector(0)];
        vel_extract.get_function_values(velocity_solution, u_tilde_vals);

        fe_values_p.get_function_gradients(deltap, grad_delta_p);

        for (unsigned int q = 0; q < n_q; ++q)
        {
            const double JxW = fe_values_vel.JxW(q);

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {

                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                    // ------
                    // L_ij = ∫ φ_i·φ_j dx
                    // ------
                    cell_matrix(i, j) += scalar_product(vel_extract.value(i, q), vel_extract.value(j, q)) * JxW;
                }
                // RHS
                // ------
                // ∫ u~·φ_i dx - (2/3) Δt ∫ ∇δp·φ_i dx
                // ------
                cell_rhs(i) += u_tilde_vals[q] * vel_extract.value(i, q) * JxW;
                cell_rhs(i) -= (2.0 / 3.0) * deltat * grad_delta_p[q] * vel_extract.value(i, q) * JxW;
            }
        }

        cell_v->get_dof_indices(local_indices);
        constraints_velocity.distribute_local_to_global(cell_matrix, cell_rhs,
                                                        local_indices,
                                                        velocity_update_matrix,
                                                        velocity_update_rhs);
    }
    velocity_update_matrix.compress(VectorOperation::add);
    velocity_update_rhs.compress(VectorOperation::add);
}

template <unsigned int dim>
void UncoupledNavierStokes<dim>::solve_update_velocity_system()
{

    TimerOutput::Scope t(computing_timer, "solve_update");

    TrilinosWrappers::MPI::Vector tmp(locally_owned_velocity, MPI_COMM_WORLD);
    SolverControl solver_control(2000, 1e-7 * velocity_update_rhs.l2_norm());

    // Jacobi or SSOR
    TrilinosWrappers::PreconditionJacobi::AdditionalData data;
    data.omega = 0.7;
    data.n_sweeps = 5;

    TrilinosWrappers::PreconditionJacobi prec;
    prec.initialize(velocity_update_matrix, data);

    SolverCG<TrilinosWrappers::MPI::Vector> solver_cg(solver_control);

    solver_cg.solve(velocity_update_matrix, tmp, velocity_update_rhs, prec);

    if (mpi_rank == 0)
        std::cout << "Velocity update CG iters: " << solver_control.last_step() << std::endl;

    constraints_velocity.distribute(tmp);
    update_velocity_solution = tmp;

    constraints_velocity.distribute(update_velocity_solution);

    update_velocity_solution.update_ghost_values();
}

template <unsigned int dim>
void UncoupledNavierStokes<dim>::output_results()
{
    DataOut<dim> data_out;

    std::vector<std::string> velocity_names(dim, "velocity");
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
        velocity_interpretation(dim,
                                DataComponentInterpretation::component_is_part_of_vector);

    data_out.attach_dof_handler(dof_handler_velocity);
    data_out.add_data_vector(update_velocity_solution,
                             velocity_names,
                             DataOut<dim>::type_dof_data,
                             velocity_interpretation);

    data_out.add_data_vector(dof_handler_pressure, pressure_solution, "pressure");
    std::vector<unsigned int> partition_int(this->mesh.n_active_cells());
    GridTools::get_subdomain_association(this->mesh, partition_int);
    const Vector<double> partitioning(partition_int.begin(), partition_int.end());
    data_out.add_data_vector(partitioning, "partitioning");

    data_out.build_patches();

    std::string numProcessors = std::to_string(this->mpi_size);
    numProcessors += (this->mpi_size == 1) ? "_processor" : "_processors";
    std::string output_dir = get_output_directory();

    std::ostringstream fname;
    fname << "solution-0_";

    std::ofstream out(fname.str());
    data_out.write_vtu_with_pvtu_record(output_dir,
                                        fname.str(),
                                        time_step,
                                        MPI_COMM_WORLD,
                                        3);
}

template <unsigned int dim>
void UncoupledNavierStokes<dim>::run()
{
    setup();

    {
        TrilinosWrappers::MPI::Vector tmp(locally_owned_velocity, MPI_COMM_WORLD);
        VectorTools::interpolate(dof_handler_velocity,
                                 Functions::ZeroFunction<dim>(dim),
                                 tmp);
        constraints_velocity.distribute(tmp);
        old_velocity = tmp;
        old_old_velocity = tmp;
        velocity_solution = tmp;
    }

    {
        TrilinosWrappers::MPI::Vector tmp(locally_owned_pressure, MPI_COMM_WORLD);
        VectorTools::interpolate(dof_handler_pressure,
                                 Functions::ZeroFunction<dim>(1),
                                 tmp);
        constraints_pressure.distribute(tmp);
        deltap = tmp;
        pressure_solution = tmp;
    }

    // Possibly output initial
    output_results();

    while (time < T - 0.5 * deltat)
    {
        time += deltat;
        ++time_step;

        time = deltat * time_step;
        if (mpi_rank == 0)
            std::cout << "\nTime step " << time_step << " at t=" << time << std::endl;

        // 1) Intermediate velocity
        assemble_system_velocity();
        solve_velocity_system();

        // 2) Pressure
        assemble_system_pressure();
        solve_pressure_system();

        // 3) Velocity update
        update_velocity();
        solve_update_velocity_system();

        // 4) Update pressure
        pressure_update(rotational);

        // Shift old velocities
        old_old_velocity = old_velocity;
        old_velocity = update_velocity_solution;

        compute_lift_drag();

        output_results();
    }
}

template <unsigned int dim>
void UncoupledNavierStokes<dim>::compute_lift_drag()
{
    // -------------------------------------------------
    // 1) Setup for face integration
    // -------------------------------------------------
    const unsigned int face_quad_degree = 3;
    QGaussSimplex<dim - 1> face_quadrature_formula(face_quad_degree);

    // FEFaceValues for velocity & pressure
    FEFaceValues<dim> fe_face_values_velocity(
        fe_velocity,
        face_quadrature_formula,
        update_values | update_gradients | update_JxW_values | update_normal_vectors);

    FEFaceValues<dim> fe_face_values_pressure(
        fe_pressure,
        face_quadrature_formula,
        update_values | update_JxW_values);

    const unsigned int n_face_q_points = face_quadrature_formula.size();

    // We'll store velocity gradients and pressure values at quadrature points
    std::vector<Tensor<2, dim>> velocity_gradients(n_face_q_points);
    std::vector<double> pressure_values(n_face_q_points);

    double local_drag = 0.0;
    double local_lift = 0.0;

    // Mean inflow velocity
    double u_max = inlet_velocity.get_u_max();
    double u_mean = 2.0 * u_max / 3.0;

    // Cylinder diameter
    double D = 0.1;

    // This factor  = 2/(ρ U^2 D)  for computing c_D and c_L
    double coefficient = 2.0 / (rho * u_mean * u_mean * D);

    // -------------------------------------------------
    // 2) Loop over cells and faces to integrate traction
    // -------------------------------------------------
    auto cell_v = dof_handler_velocity.begin_active();
    auto cell_p = dof_handler_pressure.begin_active();
    const auto end_v = dof_handler_velocity.end();

    for (; cell_v != end_v; ++cell_v, ++cell_p)
    {
        // Only integrate on locally owned cells
        if (!cell_v->is_locally_owned())
            continue;

        // Loop over faces
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
        {
            // Check boundary condition ID
            if (cell_v->face(f)->at_boundary() && (cell_v->face(f)->boundary_id() == 3))
            {
                // Reinit face-values for velocity and pressure
                fe_face_values_velocity.reinit(cell_v, f);
                fe_face_values_pressure.reinit(cell_p, f);

                // Extract velocity gradients
                fe_face_values_velocity[FEValuesExtractors::Vector(0)]
                    .get_function_gradients(velocity_solution, velocity_gradients);

                // Extract pressure values
                fe_face_values_pressure[FEValuesExtractors::Scalar(0)]
                    .get_function_values(pressure_solution, pressure_values);

                for (unsigned int q = 0; q < n_face_q_points; ++q)
                {
                    const double p = pressure_values[q];
                    const Tensor<2, dim> grad_u = velocity_gradients[q];
                    const Tensor<1, dim> normal_vec = fe_face_values_velocity.normal_vector(q);
                    const double JxW = fe_face_values_velocity.JxW(q);

                    // Build fluid_stress = -p I + 2 nu e(u)
                    // (if nu is kinematic viscosity, multiply by rho if you want μ=ρν)
                    Tensor<2, dim> fluid_stress;
                    // First, set diagonal entries to -p
                    for (unsigned int d = 0; d < dim; ++d)
                        fluid_stress[d][d] = -p;

                    // Add the viscous part = 2 nu sym_grad_u
                    for (unsigned int i = 0; i < dim; ++i)
                        for (unsigned int j = 0; j < dim; ++j)
                            fluid_stress[i][j] += nu * (grad_u[i][j] + grad_u[j][i]);

                    // Traction = stress * normal
                    const Tensor<1, dim> traction = fluid_stress * normal_vec;

                    // Multiply by area element
                    const Tensor<1, dim> force_contribution = traction * JxW;

                    // x-component => drag, y-component => lift
                    local_drag += -coefficient * force_contribution[0];
                    local_lift += coefficient * force_contribution[1];
                } // q loop
            } // boundary check
        } // face loop
    } // cell loop

    // -------------------------------------------------
    // 3) MPI: sum up partial forces to rank 0
    // -------------------------------------------------
    double global_drag = 0.0, global_lift = 0.0;
    MPI_Allreduce(&local_drag, &global_drag, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_lift, &global_lift, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // -------------------------------------------------
    // 4) Compute pressure difference between points p1 & p2
    // -------------------------------------------------
    Point<dim> p1, p2;
    if constexpr (dim == 2)
    {
        p1[0] = 0.15;
        p1[1] = 0.20;
        p2[0] = 0.25;
        p2[1] = 0.20;
    }
    else if constexpr (dim == 3)
    {
        p1[0] = 0.45;
        p1[1] = 0.20;
        p1[2] = 0.205;
        p2[0] = 0.55;
        p2[1] = 0.20;
        p2[2] = 0.205;
    }

    double local_p1 = 0.0, local_p2 = 0.0;
    bool have_p1 = false, have_p2 = false;

    try
    {
        local_p1 = VectorTools::point_value(dof_handler_pressure, pressure_solution, p1);
        have_p1 = true;
    }
    catch (...)
    {
    }
    try
    {
        local_p2 = VectorTools::point_value(dof_handler_pressure, pressure_solution, p2);
        have_p2 = true;
    }
    catch (...)
    {
    }

    double p1_on_root = 0.0, p2_on_root = 0.0;
    int rank_p1 = -1, rank_p2 = -1;

    if (have_p1)
    {
        MPI_Send(&local_p1, 1, MPI_DOUBLE, 0, 111, MPI_COMM_WORLD);
        MPI_Send(&mpi_rank, 1, MPI_INT, 0, 112, MPI_COMM_WORLD);
    }
    if (have_p2)
    {
        MPI_Send(&local_p2, 1, MPI_DOUBLE, 0, 221, MPI_COMM_WORLD);
        MPI_Send(&mpi_rank, 1, MPI_INT, 0, 222, MPI_COMM_WORLD);
    }

    if (mpi_rank == 0)
    {
        // We expect exactly one sender for p1 and p2
        MPI_Status status;
        MPI_Recv(&p1_on_root, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 111, MPI_COMM_WORLD, &status);
        MPI_Recv(&rank_p1, 1, MPI_INT, status.MPI_SOURCE, 112, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&p2_on_root, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 221, MPI_COMM_WORLD, &status);
        MPI_Recv(&rank_p2, 1, MPI_INT, status.MPI_SOURCE, 222, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        const double p_diff = p1_on_root - p2_on_root;

        // -------------------------------------------------
        // 5) Print/write output
        // -------------------------------------------------
        // e.g. append to CSV file
        std::string output_dir = get_output_directory();
        std::string filename = output_dir + "lift_drag_output.csv";

        std::ofstream out_file(filename.c_str(), std::ios::app);
        out_file << time << ","
                 << global_drag << ","
                 << global_lift << ","
                 << p_diff << "\n";
        out_file.close();
    }
}

template <unsigned int dim>
std::string UncoupledNavierStokes<dim>::get_output_directory()
{

    namespace fs = std::filesystem;

    if (!fs::exists("./outputs"))
        fs::create_directory("./outputs");
    if constexpr (dim == 2)
    {
        if (!fs::exists("./outputs/UncoupledNavierStokes2D"))
        {
            fs::create_directory("./outputs/UncoupledNavierStokes2D");
        }
    }
    else if constexpr (dim == 3)
    {
        if (!fs::exists("./outputs/UncoupledNavierStokes3D"))
            fs::create_directory("./outputs/UncoupledNavierStokes3D");
    }

    fs::path sub_dir_path = "";

    std::string sub_dir_name = "outputs_reynolds_" + std::to_string(static_cast<int>(reynolds_number));
    if constexpr (dim == 2)
        sub_dir_path = "./outputs/UncoupledNavierStokes2D/" + sub_dir_name + "/";
    else if constexpr (dim == 3)
        sub_dir_path = "./outputs/UncoupledNavierStokes3D/" + sub_dir_name + "/";

    if (!fs::exists(sub_dir_path))
    {
        fs::create_directory(sub_dir_path);
    }

    return sub_dir_path.string();
}

template <unsigned int dim>
void UncoupledNavierStokes<dim>::pressure_update(const bool rotational)
{
    // Non-rotational variant: just p^{n+1} = p^n + deltap
    if (!rotational)
    {
        pressure_solution.add(deltap); // p^{n+1} = p^n + Δp
        return;
    }

    // -------------------------
    // Rotational variant:
    // p^{n+1} = p^n + Δp - ν * div(u_tilde)
    // We must project div(u_tilde) onto the same FE space as p.
    // -------------------------

    // 1) Build an L2 mass matrix for the pressure FE space
    TrilinosWrappers::SparseMatrix mass_matrix;
    {
        TrilinosWrappers::SparsityPattern dsp_p(locally_owned_pressure, MPI_COMM_WORLD);
        DoFTools::make_sparsity_pattern(dof_handler_pressure,
                                        dsp_p,
                                        constraints_pressure,
                                        /*keep_constrained_dofs*/ false);
        dsp_p.compress();
        mass_matrix.reinit(dsp_p);
    }
    mass_matrix = 0.0;

    // 2) Build the RHS for the L2-projection of div(u_tilde)
    TrilinosWrappers::MPI::Vector rhs(locally_owned_pressure, MPI_COMM_WORLD);
    rhs = 0.0;

    const unsigned int quad_deg = std::max<unsigned int>(2u, fe_pressure.degree + 1u);
    QGaussSimplex<dim> quadrature_formula(quad_deg);

    FEValues<dim> fe_values_p(fe_pressure,
                              quadrature_formula,
                              update_values | update_quadrature_points | update_JxW_values);

    FEValues<dim> fe_values_v(fe_velocity,
                              quadrature_formula,
                              update_values | update_gradients | update_quadrature_points);

    const unsigned int n_q = quadrature_formula.size();
    const unsigned int dofs_per_cell_p = fe_pressure.dofs_per_cell;

    FullMatrix<double> cell_mass(dofs_per_cell_p, dofs_per_cell_p);
    Vector<double> cell_rhs(dofs_per_cell_p);
    std::vector<types::global_dof_index> local_indices(dofs_per_cell_p);

    // To read the velocity divergence:
    std::vector<double> local_div_u_tilde(n_q);

    // Extractor for the velocity field in the velocity FE:
    const FEValuesExtractors::Vector velocity_extract(0);

    // Loop over cells of the pressure and velocity DoFHandlers in parallel
    auto cell_p = dof_handler_pressure.begin_active();
    auto cell_v = dof_handler_velocity.begin_active();
    const auto end_p = dof_handler_pressure.end();

    for (; cell_p != end_p; ++cell_p, ++cell_v)
    {
        if (!cell_p->is_locally_owned())
            continue;

        fe_values_p.reinit(cell_p);
        fe_values_v.reinit(cell_v);

        // Get div(u_tilde) at quadrature points
        fe_values_v[velocity_extract].get_function_divergences(velocity_solution, local_div_u_tilde);

        cell_mass = 0.0;
        cell_rhs = 0.0;

        cell_p->get_dof_indices(local_indices);

        for (unsigned int q = 0; q < n_q; ++q)
        {
            const double div_val = local_div_u_tilde[q];
            const double JxW = fe_values_p.JxW(q);

            for (unsigned int i = 0; i < dofs_per_cell_p; ++i)
            {
                for (unsigned int j = 0; j < dofs_per_cell_p; ++j)
                {
                    // Mass matrix contribution
                    // ------
                    // ∫ ψ_i·ψ_j dx
                    // ------
                    cell_mass(i, j) += (fe_values_p.shape_value(i, q) * fe_values_p.shape_value(j, q)) * JxW;
                }
                // RHS 
                // ------
                // ∫ div(u~) ψ_i dx
                // ------
                cell_rhs(i) += (div_val * fe_values_p.shape_value(i, q)) * JxW;
            }
        }

        // Add local contributions to global mass matrix & RHS
        constraints_pressure.distribute_local_to_global(cell_mass,
                                                        cell_rhs,
                                                        local_indices,
                                                        mass_matrix,
                                                        rhs);
    } // end cell loop

    mass_matrix.compress(VectorOperation::add);
    rhs.compress(VectorOperation::add);

    // 3) Solve M * (divProj) = rhs for the L2-projection of div(u_tilde)
    TrilinosWrappers::MPI::Vector div_projected(locally_owned_pressure, MPI_COMM_WORLD);

    {
        SolverControl solver_control(2000, 1e-7 * rhs.l2_norm());
        SolverCG<TrilinosWrappers::MPI::Vector> solver_cg(solver_control);

        TrilinosWrappers::PreconditionIC prec;
        prec.initialize(mass_matrix);

        solver_cg.solve(mass_matrix, div_projected, rhs, prec);
        constraints_pressure.distribute(div_projected);

        if (mpi_rank == 0)
            std::cout << "Pressure update CG iterations: " << solver_control.last_step() << std::endl;
    }

    // 4) Now apply the update:
    // p^{n+1} = p^n + Δp - ν * div_projected
    // (The user formula does not have a Δt factor, but add it if your scheme needs it.)
    pressure_solution.add(deltap);             // p^{n+1} = p^n + Δp
    pressure_solution.add(-nu, div_projected); // p^{n+1} -= ν * div(u_tilde)

    pressure_solution.update_ghost_values();
}

template class UncoupledNavierStokes<2>;
template class UncoupledNavierStokes<3>;