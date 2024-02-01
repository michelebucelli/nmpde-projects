#include "Poisson3D_Parallel.hpp"


void
Poisson3DParallel::write_csv(const long int elapsed_time, int iterations, double cond_number) const{

    /*
     * We use this function to write result to the csv file.
     */
    if (mpi_rank == 0)
    {
        std::ofstream csv_file;
        csv_file.open("/u/par1/NMPDE-Project/results.csv", std::ios_base::app);

        if (!csv_file.is_open()) {
            std::cerr<<"Error! Cannot open csv file! ";
            return ;
        }

        csv_file << extractFileName(mesh_file_name) << ","
                 << p_value <<","
                 << r <<","
                 << preconditioner_name << ","
                 << elapsed_time<<","
                 << iterations <<","
                 << p_or_c <<","
                 << mpi_size <<","
                 << cond_number <<","
                 << overlap << ","
                 << sweeps<< ","
                 << omega<< ","
                 << n_spheres
                 << std::endl;
        csv_file.close();

        std::cout<< "CSV write successful."<<std::endl;

    }else{
        pcout << "Cannot use write_csv function in global scope"<<std::endl;
    }

}

std::string
Poisson3DParallel::extractFileName(const std::string& filePath) {

    /*
     * This function is only needed to write the mesh filename into the csv without its path.
     */

    size_t lastSlash = filePath.find_last_of("/\\");


    std::string fileName = filePath.substr(lastSlash + 1);


    size_t lastDot = fileName.find_last_of(".");


    return fileName.substr(0, lastDot);
}

void
Poisson3DParallel::manage_flags(int argc, char **argv) {


    /*
     * In order to efficiently test different preconditioners and their parameters we needed a way to modify their values
     * without recompiling the source files. We decided to define those values at runtime and pass them to the compiled
     * program using command line flags.
     * This function reads the given values at runtime and store its value accordingly in the class.
     * To read command line values we use the c function: getopt()
     */

    //We first default initialize the parameters

    std::string mesh_name_no_path = "mesh-cube-20.msh";
    preconditioner_name = "identity";
    p_value = 2;
    r = 1;
    p_or_c = "cube";
    std::string user_choice_for_coefficient_symmetry;

    // We modify the values

    int opt;
    const char *options = "hp:m:r:P:s:o:e:w:n:";

    while ((opt = getopt(argc, argv, options))!=-1){
        switch (opt) {
            case 'h':
                pcout << "For help see help.txt file "<<std::endl;
                break;

            case 'p':
                p_value = std::stod(optarg);
                pcout<<"Setting p="<<p_value<<std::endl;
                break;

            case 'm':
                mesh_name_no_path = optarg;
                pcout<<"Setting mesh="<<mesh_name_no_path<<std::endl;
                break;

            case 'P':
                preconditioner_name = optarg;
                pcout<<"Setting preconditioner=" <<preconditioner_name<<std::endl;
                break;

            case 'r':
                r = std::stoi(optarg);
                pcout<<"Setting r="<<r<<std::endl;
                break;

            case 's':
                user_choice_for_coefficient_symmetry = optarg;
                if (user_choice_for_coefficient_symmetry == "cube") {
                    p_or_c = "cube";
                    pcout<<"Initializing a symmetric diffusion coefficient for a cube" << std::endl;
                }
                else
                    pcout<<"Initializing a symmetric diffusion coefficient for a paralepiped"<< std::endl;

                break;

            case 'o':
                overlap = std::stoi(optarg);
                pcout<<"Setting overlap = "<<overlap<<std::endl;
                break;

            case 'e':
                sweeps = std::stoi(optarg);
                pcout<<"Setting sweeps = "<<sweeps<<std::endl;
                break;

            case 'w':
                omega = std::stod(optarg);
                pcout<<"Setting omega = " << omega << std::endl;
                break;

            case 'n':
                n_spheres = std::stoi(optarg);
                pcout<<"Setting n_shperes = " << n_spheres  <<"^3" <<std::endl;
                break;

            case '?':
                std::cerr << "Unknown command line option\n";
                return;

            default:
                pcout<<"default"<<std::endl;
                break;
        }

    }

    mesh_file_name = "/u/par1/NMPDE-Project/mesh/input_mesh/" + mesh_name_no_path;

    if(p_or_c=="cube")
        initialize_diffusion_coefficient_symmetric(p_value, n_spheres);
    else
        initialize_diffusion_coefficient_symmetric_paralepiped(p_value, n_spheres);

}


void
Poisson3DParallel::setup()
{
  pcout << "===============================================" << std::endl;

  // Create the mesh.
  {
    pcout << "Initializing the mesh" << std::endl;

    // First we read the mesh from file into a serial (i.e. not parallel)
    // triangulation.
    Triangulation<dim> mesh_serial;

    {
      GridIn<dim> grid_in;
      grid_in.attach_triangulation(mesh_serial);

      std::ifstream grid_in_file(mesh_file_name);
      grid_in.read_msh(grid_in_file);
    }

    // Then, we copy the triangulation into the parallel one.
    {
      GridTools::partition_triangulation(mpi_size, mesh_serial);
      const auto construction_data = TriangulationDescription::Utilities::
        create_description_from_triangulation(mesh_serial, MPI_COMM_WORLD);
      mesh.create_triangulation(construction_data);
    }

    // Notice that we write here the number of *global* active cells (across all
    // processes).
    pcout << "  Number of elements = " << mesh.n_global_active_cells()
          << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the finite element space. This is the same as in serial codes.
  {
    pcout << "Initializing the finite element space" << std::endl;

    fe = std::make_unique<FE_SimplexP<dim>>(r);

    pcout << "  Degree                     = " << fe->degree << std::endl;
    pcout << "  DoFs per cell              = " << fe->dofs_per_cell
          << std::endl;

    quadrature = std::make_unique<QGaussSimplex<dim>>(r + 1);

    pcout << "  Quadrature points per cell = " << quadrature->size()
          << std::endl;
   }
    quadrature_boundary = std::make_unique<QGaussSimplex<dim - 1>>(r + 1);

    pcout << "  Quadrature points per boundary cell = "
              << quadrature_boundary->size() << std::endl;

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the DoF handler.
  {
    pcout << "Initializing the DoF handler" << std::endl;

    dof_handler.reinit(mesh);
    dof_handler.distribute_dofs(*fe);

    // We retrieve the set of locally owned DoFs, which will be useful when
    // initializing linear algebra classes.
    locally_owned_dofs = dof_handler.locally_owned_dofs();

    pcout << "  Number of DoFs = " << dof_handler.n_dofs() << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the linear system.
  {
    pcout << "Initializing the linear system" << std::endl;

    pcout << "  Initializing the sparsity pattern" << std::endl;

    // To initialize the sparsity pattern, we use Trilinos' class, that manages
    // some of the inter-process communication.
    TrilinosWrappers::SparsityPattern sparsity(locally_owned_dofs,
                                               MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler, sparsity);

    // After initialization, we need to call compress, so that all process
    // retrieve the information they need for the rows they own (i.e. the rows
    // corresponding to locally owned DoFs).
    sparsity.compress();

    // Then, we use the sparsity pattern to initialize the system matrix. Since
    // the sparsity pattern is partitioned by row, so will the matrix.
    pcout << "  Initializing the system matrix" << std::endl;
    system_matrix.reinit(sparsity);

    // Finally, we initialize the right-hand side and solution vectors.
    pcout << "  Initializing the system right-hand side" << std::endl;
    system_rhs.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    pcout << "  Initializing the solution vector" << std::endl;
    solution.reinit(locally_owned_dofs, MPI_COMM_WORLD);
  }
}

void
Poisson3DParallel::assemble()
{
  pcout << "===============================================" << std::endl;

  pcout << "  Assembling the linear system" << std::endl;

  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q           = quadrature->size();

  FEValues<dim> fe_values(*fe,
                          *quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);


  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  system_matrix = 0.0;
  system_rhs    = 0.0;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      // If current cell is not owned locally, we skip it.
      if (!cell->is_locally_owned())
        continue;

      // On all other cells (which are owned by current process), we perform the
      // assembly as usual.

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
                  cell_matrix(i, j) += diffusion_coefficient.value(
                                         fe_values.quadrature_point(q)) // mu(x)
                                       * fe_values.shape_grad(i, q)     // (I)
                                       * fe_values.shape_grad(j, q)     // (II)
                                       * fe_values.JxW(q);              // (III)
                }

              cell_rhs(i) += forcing_term.value(fe_values.quadrature_point(q)) *
                             fe_values.shape_value(i, q) * fe_values.JxW(q);
            }
        }
        // If the cell is adjacent to the boundary...
      

      cell->get_dof_indices(dof_indices);

      system_matrix.add(dof_indices, cell_matrix);
      system_rhs.add(dof_indices, cell_rhs);
    }

  // Each process might have written to some rows it does not own (for instance,
  // if it owns elements that are adjacent to elements owned by some other
  // process). Therefore, at the end of the assembly, processes need to exchange
  // information: the compress method allows to do this.
  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);

  // Boundary conditions.
  {
    std::map<types::global_dof_index, double> boundary_values;

    std::map<types::boundary_id, const Function<dim> *> boundary_functions;

    for (unsigned int i = 0; i < 6; ++i)
      boundary_functions[i] = &function_g;

    VectorTools::interpolate_boundary_values(dof_handler,
                                             boundary_functions,
                                             boundary_values);

    MatrixTools::apply_boundary_values(
      boundary_values, system_matrix, solution, system_rhs, true);
  }
}

void
Poisson3DParallel::solve()
{
    /*
     * In order to choose at runtime the preconditoner we added in the solve method
     * a big if-else statement that solves the system matrix as requested by the user.
     * We initialize two instances of solver methods: one for symmetric precondtioners (SolverGC),
     * and one for unsymmetric one (solverGMRES) . Both of this solver have a method called
     * connect_condition_number_slot(). It allows to retrieve the condition number of the system matrix, once
     * the FEM linear system has been solved.
     */



  pcout << "===============================================" << std::endl;

  SolverControl solver_control(10000, 1e-6 * system_rhs.l2_norm());

  // The linear solver is basically the same as in serial, in terms of
  // interface: we only have to use appropriate classes, compatible with
  // Trilinos linear algebra.
  SolverCG<TrilinosWrappers::MPI::Vector> solver(solver_control);
  SolverGMRES<TrilinosWrappers::MPI::Vector> GMRESsolver(solver_control);
  pcout << "  Solving the linear system" << std::endl;




  double condition_number;
    /*
   * In order to compute the condition number of the system matrix we used the solvers' method
   */

    auto conditionNumberSlot = [&condition_number](double conditionNumber) {

        condition_number = conditionNumber;

        return conditionNumber;
    };

    solver.connect_condition_number_slot(conditionNumberSlot, false);
    GMRESsolver.connect_condition_number_slot(conditionNumberSlot, false);



  /*
   * We start sampling the time. The big if-else loop for sure increments execution time due to the branch hazards it
   * introduces. However, in our analysis we want to study not the execution time per-se but the trend of the execution
   * time as the number of processes increases. Therefore, we believe we can safely ignore this overhead.
   */



  const auto t0 = std::chrono::high_resolution_clock::now();


    if (preconditioner_name == "identity")
    {
        pcout<<"Using preconditioner identity"<<std::endl;
        TrilinosWrappers::PreconditionIdentity preconditioner;
        solver.solve(system_matrix, solution, system_rhs, preconditioner);

    }else if (preconditioner_name == "jacobi") {
        pcout << "Using preconditioner Jacobi" << std::endl;

        // Configuring parameters for the SOR preconditioner
        TrilinosWrappers::PreconditionJacobi::AdditionalData jacobi_data;


        jacobi_data.omega = omega;
        jacobi_data.min_diagonal=0;
        jacobi_data.n_sweeps=sweeps;


        TrilinosWrappers::PreconditionJacobi preconditioner;
        preconditioner.initialize(system_matrix, jacobi_data);


        solver.solve(system_matrix, solution, system_rhs, preconditioner);
    }else if (preconditioner_name == "ssor") {
        pcout << "Using preconditioner SSOR" << std::endl;

        // Configuring parameters for the SOR preconditioner
        TrilinosWrappers::PreconditionSSOR::AdditionalData ssor_data;


        ssor_data.omega = omega;
        ssor_data.min_diagonal=0.29;
        ssor_data.overlap=overlap;
        ssor_data.n_sweeps=sweeps;

        // Create and initialize the Successive Overrelaxation (SOR) preconditioner
        TrilinosWrappers::PreconditionSSOR preconditioner;
        preconditioner.initialize(system_matrix, ssor_data);

        // Solve the linear system using GMRES and the Successive Overrelaxation (SOR) preconditioner
        GMRESsolver.solve(system_matrix, solution, system_rhs, preconditioner);
    }else if (preconditioner_name == "amg") {
        pcout << "Using Algebraic Multigrid (AMG) preconditioner" << std::endl;

        // Configuring parameters for the AMG preconditioner
        TrilinosWrappers::PreconditionAMG::AdditionalData amg_data;

        // Set parameters based on your problem characteristics
        amg_data.elliptic = true;  // Adjust based on the nature of your problem
        amg_data.higher_order_elements = false;  // Adjust if using higher-order elements
        amg_data.n_cycles = 1;  // Number of multigrid cycles
        amg_data.w_cycle = false;  // Use w-cycle if needed
        amg_data.aggregation_threshold = 1e-4;  // Threshold for coarsening
        amg_data.constant_modes = std::vector<std::vector<bool>>(0);  // Constant modes for null space
        amg_data.smoother_sweeps = sweeps;  // Number of smoother sweeps
        amg_data.smoother_overlap = 0;  // Smoother overlap in parallel
        amg_data.output_details = false;  // Output internal details
        amg_data.smoother_type = "Chebyshev";  // Smoother type
        amg_data.coarse_type = "Amesos-KLU";  // Coarsest level solver type

        // Create and initialize the Algebraic Multigrid (AMG) preconditioner
        TrilinosWrappers::PreconditionAMG preconditioner;
        preconditioner.initialize(system_matrix, amg_data);

        TrilinosWrappers::PreconditionAMG::size_type memoryUsage = preconditioner.memory_consumption();

        // Solve the linear system using GMRES and the AMG preconditioner
        GMRESsolver.solve(system_matrix, solution, system_rhs, preconditioner);
        pcout << "Memory Consumption: " << memoryUsage << " bytes" << std::endl;
    }else{
        pcout<<"Error! Preconditioner \" " <<preconditioner_name <<" \" not supported!"<<std::endl;
        std::exit(-1);
    }

    const auto t1= std::chrono::high_resolution_clock::now();


    const auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

    pcout<<"Elapsed time for solve phase: "
             << dt << " [ms] " << std::endl;

    pcout << "  " << solver_control.last_step() << "  iterations" << std::endl;
    pcout << "Condition Number: " << condition_number << std::endl;

  if (mpi_rank == 0)
      write_csv(dt, solver_control.last_step(), condition_number);

}

void
Poisson3DParallel::output() const
{
  pcout << "===============================================" << std::endl;

  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

  // To correctly export the solution, each process needs to know the solution
  // DoFs it owns, and the ones corresponding to elements adjacent to the ones
  // it owns (the locally relevant DoFs, or ghosts). We create a vector to store
  // them.
  TrilinosWrappers::MPI::Vector solution_ghost(locally_owned_dofs,
                                               locally_relevant_dofs,
                                               MPI_COMM_WORLD);

  // This performs the necessary communication so that the locally relevant DoFs
  // are received from other processes and stored inside solution_ghost.
  solution_ghost = solution;

  // Then, we build and fill the DataOut class as usual.
  DataOut<dim> data_out;
  data_out.add_data_vector(dof_handler, solution_ghost, "solution");

  // We also add a vector to represent the parallel partitioning of the mesh.
  std::vector<unsigned int> partition_int(mesh.n_active_cells());
  GridTools::get_subdomain_association(mesh, partition_int);
  const Vector<double> partitioning(partition_int.begin(), partition_int.end());
  data_out.add_data_vector(partitioning, "partitioning");

  data_out.build_patches();

  const std::filesystem::path mesh_path(mesh_file_name);
  const std::string output_file_name = "/u/par1/NMPDE-Project/mesh/output_mesh/output-" + mesh_path.stem().string();

  // Finally, we need to write in a format that supports parallel output. This
  // can be achieved in multiple ways (e.g. XDMF/H5). We choose VTU/PVTU files,
  // because the interface is nice and it is quite robust.
  data_out.write_vtu_with_pvtu_record("./",
                                      output_file_name,
                                      0,
                                      MPI_COMM_WORLD);

  pcout << "Output written to " << output_file_name << std::endl;

  pcout << "===============================================" << std::endl;
}
