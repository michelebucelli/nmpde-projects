#include "Poisson3D.hpp"


void
Poisson3D::write_csv(const long int elapsed_time, int iterations) const{
    std::ofstream csv_file;
    csv_file.open("../results/results.csv", std::ios_base::app);

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
             << symmetric <<","
             << 0
             << std::endl;
    csv_file.close();

    std::cout<< "CSV write successful."<<std::endl;
}

std::string
Poisson3D::extractFileName(const std::string& filePath) {

    size_t lastSlash = filePath.find_last_of("/\\");


    std::string fileName = filePath.substr(lastSlash + 1);


    size_t lastDot = fileName.find_last_of(".");


    return fileName.substr(0, lastDot);
}

void
Poisson3D::manage_flags(int argc, char **argv) {


    /*
     * DEFAULT INITIALIZATIONS OF PARAMETERS
     */

    std::string mesh_name_no_path = "mesh-cube-20.msh";
    preconditioner_name = "identity";
    p_value = 2;
    r = 1;
    symmetric = true;
    std::string user_choice_for_coefficient_symmetry;


    /*
     * MANAGING COMMAND LINE FLAGS
     */

    int opt;
    const char *options = "hp:m:r:P:s:";

    while ((opt = getopt(argc, argv, options))!=-1){
        switch (opt) {
            case 'h':
                std::cout << "Help message "<<std::endl;
                break;

            case 'p':
                p_value = std::stod(optarg);
                std::cout<<"Setting p="<<p_value<<std::endl;
                break;

            case 'm':
                mesh_name_no_path = optarg;
                std::cout<<"Setting mesh="<<mesh_name_no_path<<std::endl;
                break;

            case 'P':
                preconditioner_name = optarg;
                std::cout<<"Setting preconditioner=" <<preconditioner_name<<std::endl;
                break;

            case 'r':
                r = std::stoi(optarg);
                std::cout<<"Setting r="<<r<<std::endl;
                break;

            case 's':
                user_choice_for_coefficient_symmetry = optarg;
                if (user_choice_for_coefficient_symmetry == "no") {
                    symmetric = false;
                    std::cout<<"Initializing randomly an unsymmetric diffusion coefficient " << std::endl;
                }
                else
                    std::cout<<"Initializing non-randomly a symmetric diffusion coefficient "<< std::endl;

                break;
            case '?':
                std::cerr << "Unknown command line option\n";
                return;

            default:
                std::cout<<"default"<<std::endl;
                break;
        }

    }

   mesh_file_name = "../mesh/" + mesh_name_no_path;

   if(symmetric)
       initialize_diffusion_coefficient_symmetric(p_value);
   else
       initialize_diffusion_coefficient(p_value);

}

void
Poisson3D::setup()
{
  std::cout << "===============================================" << std::endl;

  // Create the mesh.
  {
    std::cout << "Initializing the mesh from " << mesh_file_name << std::endl;

    // Or read the mesh from file:
    GridIn<dim> grid_in;
    grid_in.attach_triangulation(mesh);

    std::ifstream mesh_file(mesh_file_name);
    grid_in.read_msh(mesh_file);

    std::cout << "  Number of elements = " << mesh.n_active_cells()
              << std::endl;
  }

  std::cout << "-----------------------------------------------" << std::endl;

  // Initialize the finite element space.
  {
    std::cout << "Initializing the finite element space" << std::endl;

    // Construct the finite element object. Notice that we use the FE_SimplexP
    // class here, that is suitable for triangular (or tetrahedral) meshes.
    fe = std::make_unique<FE_SimplexP<dim>>(r);

    std::cout << "  Degree                     = " << fe->degree << std::endl;
    std::cout << "  DoFs per cell              = " << fe->dofs_per_cell
              << std::endl;

    // Construct the quadrature formula of the appopriate degree of exactness.
    quadrature = std::make_unique<QGaussSimplex<dim>>(r + 1);

    std::cout << "  Quadrature points per cell = " << quadrature->size()
              << std::endl;

    quadrature_boundary = std::make_unique<QGaussSimplex<dim - 1>>(r + 1);

    std::cout << "  Quadrature points per boundary cell = "
              << quadrature_boundary->size() << std::endl;
  }

  std::cout << "-----------------------------------------------" << std::endl;

  // Initialize the DoF handler.
  {
    std::cout << "Initializing the DoF handler" << std::endl;

    // Initialize the DoF handler with the mesh we constructed.
    dof_handler.reinit(mesh);

    // "Distribute" the degrees of freedom. For a given finite element space,
    // initializes info on the control variables (how many they are, where
    // they are collocated, their "global indices", ...).
    dof_handler.distribute_dofs(*fe);

    std::cout << "  Number of DoFs = " << dof_handler.n_dofs() << std::endl;
  }

  std::cout << "-----------------------------------------------" << std::endl;

  // Initialize the linear system.
  {
    std::cout << "Initializing the linear system" << std::endl;

    // We first initialize a "sparsity pattern", i.e. a data structure that
    // indicates which entries of the matrix are zero and which are different
    // from zero. To do so, we construct first a DynamicSparsityPattern (a
    // sparsity pattern stored in a memory- and access-inefficient way, but
    // fast to write) and then convert it to a SparsityPattern (which is more
    // efficient, but cannot be modified).
    std::cout << "  Initializing the sparsity pattern" << std::endl;
    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    // Then, we use the sparsity pattern to initialize the system matrix
    std::cout << "  Initializing the system matrix" << std::endl;
    system_matrix.reinit(sparsity_pattern);

    // Finally, we initialize the right-hand side and solution vectors.
    std::cout << "  Initializing the system right-hand side" << std::endl;
    system_rhs.reinit(dof_handler.n_dofs());
    std::cout << "  Initializing the solution vector" << std::endl;
    solution.reinit(dof_handler.n_dofs());
  }
}

void
Poisson3D::assemble()
{
  std::cout << "===============================================" << std::endl;

  std::cout << "  Assembling the linear system" << std::endl;

  // Number of local DoFs for each element.
  const unsigned int dofs_per_cell = fe->dofs_per_cell;

  // Number of quadrature points for each element.
  const unsigned int n_q = quadrature->size();

  // FEValues instance. This object allows to compute basis functions, their
  // derivatives, the reference-to-current element mapping and its
  // derivatives on all quadrature points of all elements.
  FEValues<dim> fe_values(
    *fe,
    *quadrature,
    // Here we specify what quantities we need FEValues to compute on
    // quadrature points. For our test, we need:
    // - the values of shape functions (update_values);
    // - the derivative of shape functions (update_gradients);
    // - the position of quadrature points (update_quadrature_points);
    // - the product J_c(x_q)*w_q (update_JxW_values).
    update_values | update_gradients | update_quadrature_points |
      update_JxW_values);

  // Since we need to compute integrals on the boundary for Neumann conditions,
  // we also need a FEValues object to compute quantities on boundary edges
  // (faces).
  FEFaceValues<dim> fe_values_boundary(*fe,
                                       *quadrature_boundary,
                                       update_values |
                                         update_quadrature_points |
                                         update_JxW_values);

  // Local matrix and right-hand side vector. We will overwrite them for
  // each element within the loop.
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  // We will use this vector to store the global indices of the DoFs of the
  // current element within the loop.
  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  // Reset the global matrix and vector, just in case.
  system_matrix = 0.0;
  system_rhs    = 0.0;

  //Set to get the boundaryIds
  std::set<types::boundary_id> boundary_ids;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      for (const auto &face : cell->face_iterators()) {
        if (face->at_boundary()) {
            boundary_ids.insert(face->boundary_id());
        }
    }
      
      // Reinitialize the FEValues object on current element. This
      // precomputes all the quantities we requested when constructing
      // FEValues (see the update_* flags above) for all quadrature nodes of
      // the current cell.
      fe_values.reinit(cell);

      // We reset the cell matrix and vector (discarding any leftovers from
      // previous element).
      cell_matrix = 0.0;
      cell_rhs    = 0.0;

      for (unsigned int q = 0; q < n_q; ++q)
        {
          // Here we assemble the local contribution for current cell and
          // current quadrature point, filling the local matrix and vector.

          // Here we iterate over *local* DoF indices.
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                  // FEValues::shape_grad(i, q) returns the gradient of the i-th
                  // basis function at the q-th quadrature node, already mapped
                  // on the physical element: we don't have to deal with the
                  // mapping, it's all hidden inside FEValues.
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


      // At this point the local matrix and vector are constructed: we
      // need to sum them into the global matrix and vector. To this end,
      // we need to retrieve the global indices of the DoFs of current
      // cell.
      cell->get_dof_indices(dof_indices);

      // Then, we add the local matrix and vector into the corresponding
      // positions of the global matrix and vector.
      system_matrix.add(dof_indices, cell_matrix);
      system_rhs.add(dof_indices, cell_rhs);
    }

  // Boundary conditions.
  {
    // We construct a map that stores, for each DoF corresponding to a
    // Dirichlet condition, the corresponding value. E.g., if the Dirichlet
    // condition is u_i = b, the map will contain the pair (i, b).
    std::map<types::global_dof_index, double> boundary_values;

    // Then, we build a map that, for each boundary tag, stores the
    // corresponding boundary function.
    std::map<types::boundary_id, const Function<dim> *> boundary_functions;
    for (const auto &boundary_id : boundary_ids) {
      boundary_functions[boundary_id] = &function_g;
    }

    //interpolate_boundary_values fills the boundary_values map.
    VectorTools::interpolate_boundary_values(dof_handler,
                                             boundary_functions,
                                             boundary_values);

    // Finally, we modify the linear system to apply the boundary
    // conditions. This replaces the equations for the boundary DoFs with
    // the corresponding u_i = 0 equations.
    MatrixTools::apply_boundary_values(
      boundary_values, system_matrix, solution, system_rhs, true);
  }
}

void
Poisson3D::solve()
{
  std::cout << "===============================================" << std::endl;

  // Here we specify the maximum number of iterations of the iterative solver,
  // and its tolerance.
  SolverControl solver_control(10000, 1e-6 * system_rhs.l2_norm());
  for (unsigned int i = 0; i < solution.size(); ++i) {
      solution[i] = 0.0;
  }


  const auto t0 = std::chrono::high_resolution_clock::now();

  SolverCG<Vector<double>> solver(solver_control);
  SolverGMRES<Vector<double>> GMRESsolver(solver_control);


  if (preconditioner_name == "identity")
  {
      std::cout<<"Using preconditioner identity"<<std::endl;
      //TrilinosWrappers::PreconditionIdentity preconditioner;
      PreconditionIdentity preconditioner;
      solver.solve(system_matrix, solution, system_rhs, preconditioner);
  }else if (preconditioner_name == "jacobi"){
      std::cout<<"Using preconditioner jacobi"<<std::endl;
      PreconditionJacobi preconditioner;
      preconditioner.initialize(system_matrix);
      solver.solve(system_matrix, solution, system_rhs, preconditioner);
  }else if (preconditioner_name == "ssor"){
      std::cout<<"Using preconditioner ssor"<<std::endl;
      PreconditionSSOR preconditioner;
      preconditioner.initialize(
        system_matrix,
        PreconditionSSOR<SparseMatrix<double>>::AdditionalData(1.0)
        );
      solver.solve(system_matrix, solution, system_rhs, preconditioner);

  }else if (preconditioner_name == "sor"){

      std::cout<<"Using preconditioner sor"<<std::endl;
      std::cout<<"Not symmetric preconditioner, hence solving with GMRES and not GC"<<std::endl;
      PreconditionSOR preconditioner;
      preconditioner.initialize(
        system_matrix, PreconditionSOR<SparseMatrix<double>>::AdditionalData(1.0)
        );
      GMRESsolver.solve(system_matrix, solution, system_rhs, preconditioner);
  }
  else{
      std::cerr<<"Error! Preconditioner not supported!"<<std::endl;
      std::exit(-1);
  }


  const auto t1 = std::chrono::high_resolution_clock::now();
  const auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
  std::cout<<"Elapsed time for solve phase: "
             << dt << " [ms] " << std::endl;



  if (isMatrixSingular(system_matrix)) {
        std::cout << "The matrix is singular." << std::endl;
    } else {
        std::cout << "The matrix isn't singular." << std::endl;
    }
  std::cout << "  Solving the linear system" << std::endl;


  std::cout << "  " << solver_control.last_step() << " iterations"
            << std::endl;

  write_csv(dt, solver_control.last_step());

}

void
Poisson3D::output() const
{
  std::cout << "===============================================" << std::endl;

  // The DataOut class manages writing the results to a file.
  DataOut<dim> data_out;

  // It can write multiple variables (defined on the same mesh) to a single
  // file. Each of them can be added by calling add_data_vector, passing the
  // associated DoFHandler and a name.
  data_out.add_data_vector(dof_handler, solution, "solution");

  // Once all vectors have been inserted, call build_patches to finalize the
  // DataOut object, preparing it for writing to file.
  data_out.build_patches();

  // Then, use one of the many write_* methods to write the file in an
  // appropriate format.
  const std::filesystem::path mesh_path(mesh_file_name);
  const std::string           output_file_name =
    "output-" + mesh_path.stem().string() + ".vtk";
  std::ofstream output_file(output_file_name);
  data_out.write_vtk(output_file);

  std::cout << "Output written to " << output_file_name << std::endl;

  std::cout << "===============================================" << std::endl;
}



