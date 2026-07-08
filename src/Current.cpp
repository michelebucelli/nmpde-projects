#include "Current.hpp"
#include <fstream>
#include <string>
#include <cstdio>

const char* model_file = "models/models.cvs";

void Current::set_model(const std::string& model){
  std::string current_model;
  const char* csv_keys[]{"model","v_o","v_u","theta_w1","theta_w2","theta_w1_m","theta_o","tau_w1_1_m","tau_w1_2_m","tau_w1_p","tau_w2_1_m","tau_w2_2_m","k_w2_m","v_w2_m","tau_w2_p","tau_fi","tau_o1","tau_o2","tau_so1","tau_so2","k_so","v_so","tau_s1","tau_w3_2","k_w3","v_w3","tau_si","tau_w2_inf","w2_inf_star"};
  const void* parameters_list[]{&current_model,&v_o,&v_u,&theta_w1,&theta_w2,&theta_w1_m,&theta_o,&tau_w1_1_m,&tau_w1_2_m,&tau_w1_p,&tau_w2_1_m,&tau_w2_2_m,&k_w2_m,&v_w2_m,&tau_w2_p,&tau_fi,&tau_o1,&tau_o2,&tau_so1,&tau_so2,&k_so,&v_so,&tau_w3_1,&tau_w3_2,&k_w3,&v_w3,&tau_si,&tau_w2_inf,&w2_inf_star};
  constexpr size_t num_parameters = sizeof(csv_keys) / sizeof(char*);
  // std::unordered_map<std::string, size_t> parameter_indices;

  // for (int i = 0; i < sizeof(csv_keys) / sizeof(char*); ++i)
  //   parameter_indices.insert({csv_keys[i], -1});
  
  std::ifstream model_file("models/models.csv");
  std::string line;
  std::string cell;
  std::istringstream line_stream;
  if (std::getline(model_file, line)){
    line_stream.clear();
    line_stream.str(line);
    
    for (size_t i = 0; i < num_parameters; ++i){
      line_stream >> std::ws;
      if(std::getline(line_stream, cell, ',')){
        {
          size_t end = cell.find_first_of(" \t");
          if (end < cell.size())
            cell.erase(end);
        }
        if(cell.compare(csv_keys[i]) != 0){
          std::cerr << "Unexpected key " << cell << std::endl;
          std::abort();
        }
      }
    }
  }

  while (model.compare(current_model) != 0 && std::getline(model_file, line)){
    line_stream.clear();
    line_stream.str(line);
    line_stream >> std::ws;
    std::getline(line_stream, current_model, ',');
    {
      size_t end = current_model.find_first_of(" \t");
      if (end < current_model.size())
        current_model.erase(end);
    }
    if (model.compare(current_model) == 0){
      for(size_t i = 1; i < num_parameters; ++i){
      std::getline(line_stream, cell, ',');
        *((double*) parameters_list[i]) = std::atof(cell.c_str());
        // if (!line_stream.good()){
        //   pcout << "Read leads to an error" << std::endl;
        //   abort();
        // }
      }
    }
  }

  if (model.compare(current_model) != 0){
    std::cerr << "Model " << model << " not found" << std::endl;
    abort();
  }

}

void
Current::setup()
{
  pcout << "===============================================" << std::endl;

  // Create the mesh.
  {
    pcout << "Initializing the mesh" << std::endl;

    // Read serial mesh.
    Triangulation<dim> mesh_serial;

    {
      GridIn<dim> grid_in;
      grid_in.attach_triangulation(mesh_serial);

      std::ifstream mesh_file(mesh_file_name);
      grid_in.read_msh(mesh_file);
    }

    // Copy the serial mesh into the parallel one.
    {
      GridTools::partition_triangulation(mpi_size, mesh_serial);

      const auto construction_data = TriangulationDescription::Utilities::
        create_description_from_triangulation(mesh_serial, MPI_COMM_WORLD);
      mesh.create_triangulation(construction_data);
    }
    pcout << "  Number of elements = " << mesh.n_global_active_cells()
          << std::endl;
      
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize anisotropic diffusion tensor from monodomain conductivities.
  {
    const double sigma_long = sigma_i_long * sigma_e_long /
                              (sigma_i_long + sigma_e_long);
    const double sigma_trans = sigma_i_trans * sigma_e_trans /
                               (sigma_i_trans + sigma_e_trans);

    const double sigma_long_mm = sigma_long * 1e-3; // S/mm
    const double sigma_trans_mm = sigma_trans * 1e-3; // S/mm

    const double D_long = sigma_long_mm / (Cm * beta_mm) / 1000.0;   // mm^2/ms
    const double D_trans = sigma_trans_mm / (Cm * beta_mm) / 1000.0;  // mm^2/ms

    diffusion_tensor.clear();
    diffusion_tensor[0][0] = D_long;
    diffusion_tensor[1][1] = D_trans;
    diffusion_tensor[2][2] = D_trans;

    pcout << "Diffusion tensor initialized:" << std::endl
          << "  D_long = " << D_long << std::endl
          << "  D_trans = " << D_trans << std::endl;
  }

  // Initialize the finite element space.
  {
    pcout << "Initializing the finite element space" << std::endl;

    fe = std::make_unique<FE_Q<dim>>(r);

    pcout << "  Degree                     = " << fe->degree << std::endl;
    pcout << "  DoFs per cell              = " << fe->dofs_per_cell
          << std::endl;

    quadrature = std::make_unique<QGauss<dim>>(r + 1);

    pcout << "  Quadrature points per cell = " << quadrature->size()
          << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the DoF handler.
  {
    pcout << "Initializing the DoF handler" << std::endl;

    dof_handler.reinit(mesh);
    dof_handler.distribute_dofs(*fe);

    pcout << "  Number of DoFs = " << dof_handler.n_dofs() << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the linear system.
  {
    pcout << "Initializing the linear system" << std::endl;

    const IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
    const IndexSet locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);

    pcout << "  Initializing the sparsity pattern" << std::endl;
    TrilinosWrappers::SparsityPattern sparsity(locally_owned_dofs,
                                               MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler, sparsity);
    sparsity.compress();

    pcout << "  Initializing the system matrix" << std::endl;
    system_matrix.reinit(sparsity);

    pcout << "  Initializing vectors" << std::endl;
    system_rhs.reinit(locally_owned_dofs, MPI_COMM_WORLD);

    v.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    v_ghost.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);

    activation_time.reinit(locally_owned_dofs, MPI_COMM_WORLD);

    w1.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    w2.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    w3.reinit(locally_owned_dofs, MPI_COMM_WORLD);

    J_fi.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    J_so.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    J_si.reinit(locally_owned_dofs, MPI_COMM_WORLD);

    J_fi_ghost.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
    J_so_ghost.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
    J_si_ghost.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
  }

  {
   pcout << "  Assembly system matrix" << std::endl;

   // Number of local DoFs for each element.
   const unsigned int dofs_per_cell = fe->dofs_per_cell;

   // Number of quadrature points for each element.
   const unsigned int n_q = quadrature->size();

   FEValues<dim> fe_values(*fe,
                           *quadrature,
                           update_values | update_gradients |
                             update_quadrature_points | update_JxW_values);

   // Matrix
   FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

   std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

   // Reset the global matrix and vector, just in case.
   system_matrix = 0.0;


   for (const auto &cell : dof_handler.active_cell_iterators())
     {
       if (!cell->is_locally_owned())
         continue;

       fe_values.reinit(cell);

       cell_matrix = 0.0;

       for (unsigned int q = 0; q < n_q; ++q)
         {
           for (unsigned int i = 0; i < dofs_per_cell; ++i)
             {
               for (unsigned int j = 0; j < dofs_per_cell; ++j)
                 {
                   // Time derivative term: (u - u_old) / delta_t
                   cell_matrix(i, j) += (1.0 / delta_t) *             //
                                        fe_values.shape_value(i, q) * //
                                        fe_values.shape_value(j, q) * //
                                        fe_values.JxW(q);

                   // Diffusion.
                   cell_matrix(i, j) +=
                     theta * scalar_product(fe_values.shape_grad(i, q),
                                           diffusion_tensor *
                                             fe_values.shape_grad(j, q)) *
                     fe_values.JxW(q);

                 }
             }
         }

       cell->get_dof_indices(dof_indices);
       system_matrix.add(dof_indices, cell_matrix);
     }

   system_matrix.compress(VectorOperation::add);
   solver->initialize(system_matrix);
  }
}

void Current::integrate_auxiliar_variables(){
  double w1_old;
  double w2_old;
  double w3_old;
  for (auto i: v.locally_owned_elements()){
    w1_old = w1[i];
    w2_old = w2[i];
    w3_old = w3[i];
    
    if(v[i]< theta_w1_m)
      w1[i] = w1_old + delta_t * (1-w1_old)/tau_w1_1_m;
    else if(v[i]<theta_w1)
      w1[i] = w1_old + delta_t *(-w1_old/tau_w1_2_m);
    else
      w1[i] = w1_old + delta_t *(-w1_old/tau_w1_p);


    {
      double tau_w2_m
      = tau_w2_1_m + 0.5*(tau_w2_2_m - tau_w2_1_m)*(1+ tanh(k_w2_m*(v[i]-v_w2_m)));

      if(v[i]< theta_o)
        w2[i] = w2_old + delta_t * (1 - v[i] /tau_w2_inf - w2_old)/ tau_w2_m;
      else if(v[i]<theta_w2)    
        w2[i] = w2_old + delta_t *  (w2_inf_star-w2_old) / tau_w2_m;
      else
        w2[i] = w2_old + delta_t * (-w2_old/tau_w2_p);
    }

    {
      double tau_w3 = v[i] < theta_w2 ? tau_w3_1 : tau_w3_2;

      w3[i] = w3_old + delta_t * (1 + tanh(k_w3*(v[i]-v_w3)) - 2*w3_old)/(2 * tau_w3);
    }

  }
}

void Current::compute_ionic_currents(){
  for (auto i: v.locally_owned_elements()){
    if (v[i] >= theta_w1)
      J_fi[i] = w1[i] * (v[i] - theta_w1) * (v[i] - v_u) / tau_fi;
    else
      J_fi[i] = 0.0;

    if(v[i] < theta_o)  
      J_so[i] = (v[i] - v_o)/tau_o1;
    else if(v[i] < theta_w2)  
      J_so[i] = (v[i] - v_o)/tau_o2;
    else
      J_so[i] = 2/(2*tau_so1 + (tau_so2-tau_so1)*( 1 + tanh(k_so*(v[i] - v_so))));
    
    if(v[i]>= theta_w2)
      J_si[i] = -w2[i]*w3[i] / tau_si;
    else
      J_si[i] = 0.0; 

  } 

  J_fi_ghost = J_fi;
  J_so_ghost = J_so;
  J_si_ghost = J_si;
}

 void
 Current::assemble()
 {
   // Number of local DoFs for each element.
   const unsigned int dofs_per_cell = fe->dofs_per_cell;

   // Number of quadrature points for each element.
   const unsigned int n_q = quadrature->size();

   FEValues<dim> fe_values(*fe,
                           *quadrature,
                           update_values | update_gradients |
                             update_quadrature_points | update_JxW_values);

   // Local vector.
   Vector<double>     cell_rhs(dofs_per_cell);

   std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

   // Reset the global matrix and vector, just in case.
   system_rhs    = 0.0;

   // Evaluation of the old solution on quadrature nodes of current cell.
   std::vector<double> solution_old_values(n_q);

   // Evaluation of the gradient of the old solution on quadrature nodes of
   // current cell.
   std::vector<Tensor<1, dim>> solution_old_grads(n_q);

   std::vector<double> J_fi_values(n_q);
   std::vector<double> J_so_values(n_q);
   std::vector<double> J_si_values(n_q);

   for (const auto &cell : dof_handler.active_cell_iterators())
     {
       if (!cell->is_locally_owned())
         continue;

       fe_values.reinit(cell);

       cell_rhs    = 0.0;

       // Evaluate the old solution and its gradient on quadrature nodes.
       fe_values.get_function_values(v_ghost, solution_old_values);
       fe_values.get_function_gradients(v_ghost, solution_old_grads);

       fe_values.get_function_values(J_fi_ghost, J_fi_values);
       fe_values.get_function_values(J_so_ghost, J_so_values);
       fe_values.get_function_values(J_si_ghost, J_si_values);

       for (unsigned int q = 0; q < n_q; ++q)
         {
           
           for (unsigned int i = 0; i < dofs_per_cell; ++i)
             {
               // Time derivative.
               cell_rhs(i) += (1.0 / delta_t) *             //
                              fe_values.shape_value(i, q) * //
                              solution_old_values[q] *      //
                              fe_values.JxW(q);

               // Diffusion.
               cell_rhs(i) -= (1.0 - theta) *
                              scalar_product(fe_values.shape_grad(i, q),
                                             diffusion_tensor *
                                               solution_old_grads[q]) *
                              fe_values.JxW(q);

               // ionic currents
               cell_rhs(i) -=
                 (J_fi_values[q] + J_so_values[q] + J_si_values[q]) * //
                 fe_values.shape_value(i, q) *                     //
                 fe_values.JxW(q);

              // the forcing term is a current applied on a cubic region
              //on the sideof the domain
              Point<dim> q_point = fe_values.quadrature_point(q);
              if(q_point[0] <= 1.5 && q_point[1] <= 1.5 && q_point[2] <= 1.5 && time <= 2.0){
                cell_rhs[i] += 0.416 * fe_values.shape_value(i,q)* fe_values.JxW(q); // microA/cm^3
              }
             }
         }

       cell->get_dof_indices(dof_indices);

       system_rhs.add(dof_indices, cell_rhs);
     }

   system_rhs.compress(VectorOperation::add);

   // Homogeneous Neumann boundary conditions: we do nothing.
 }

void
Current::solve_linear_system()
{
  solver->solve(v, system_rhs);
  pcout << solver->get_iterations() << " " << solver->get_name() << " iterations" << std::endl;
}

void
Current::output() const
{
  DataOut<dim> data_out;

  data_out.add_data_vector(dof_handler, v_ghost, "membrane tension");

  data_out.add_data_vector(dof_handler, w1, "w1");
  data_out.add_data_vector(dof_handler, w2, "w2");
  data_out.add_data_vector(dof_handler, w3, "w3");
  
  data_out.add_data_vector(dof_handler, J_fi, "J_fi");
  data_out.add_data_vector(dof_handler, J_so, "J_so");
  data_out.add_data_vector(dof_handler, J_si, "J_si");

  // Add vector for parallel partition.
  std::vector<unsigned int> partition_int(mesh.n_active_cells());
  GridTools::get_subdomain_association(mesh, partition_int);
  const Vector<double> partitioning(partition_int.begin(), partition_int.end());
  data_out.add_data_vector(partitioning, "partitioning");

  data_out.build_patches();

  const std::filesystem::path mesh_path(mesh_file_name);
  const std::string output_file_name = "output-" + mesh_path.stem().string();

  data_out.write_vtu_with_pvtu_record(/* folder = */ "./",
                                      /* basename = */ output_file_name,
                                      /* index = */ timestep_number,
                                      MPI_COMM_WORLD);
}


void
Current::check_activation_time()
{  
  has_non_activated_cells = false;
  for (auto i: activation_time.locally_owned_elements())
    if (std::isnan(activation_time[i])){
      if (v_ghost[i] > 84/85.7)
        activation_time[i] = time;
      else
        has_non_activated_cells = true;
    }
  
  MPI_Allreduce(MPI_IN_PLACE, &has_non_activated_cells, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
}

void
Current::output_activation_time() const
{
  DataOut<dim> data_out;

  data_out.add_data_vector(dof_handler, activation_time, "Activation Time");

  std::filesystem::path output_file_path(output_file_name);
  std::filesystem::current_path(output_file_path.parent_path());

  // Add vector for parallel partition.
  std::vector<unsigned int> partition_int(mesh.n_active_cells());
  GridTools::get_subdomain_association(mesh, partition_int);
  const Vector<double> partitioning(partition_int.begin(), partition_int.end());
  data_out.add_data_vector(partitioning, "partitioning");

  data_out.build_patches();

  data_out.write_vtu_with_pvtu_record(/* folder = */ "./",
                                      /* basename = */ output_file_path.filename().string(),
                                      /* index = */ 0,
                                      MPI_COMM_WORLD);
}




void
Current::run()
{
  // Setup initial conditions.
  {
    setup();

    VectorTools::interpolate(dof_handler, Functions::ZeroFunction<dim>(), v);
    v_ghost = v;
      
   
    for (auto i : activation_time.locally_owned_elements())
      activation_time[i] = NAN;

    // -------------------------------------------------------
    // Iniziliaze the auxilial values to initial conditions.
    //   w1 = 1 -> canale sodium channel ready to be opened
    //   w2 = 1 -> calcium channel ready to be opened  
    //   w3 = 0 -> 4-th setted to zero
    // -------------------------------------------------------
    VectorTools::interpolate(dof_handler, Functions::ConstantFunction<dim>(1.0), w1);
    VectorTools::interpolate(dof_handler, Functions::ConstantFunction<dim>(1.0), w2);
    VectorTools::interpolate(dof_handler, Functions::ZeroFunction<dim>(), w3);

    pcout << "  Ionic variables initialized (w1=1, w2=1, w3=0)" << std::endl;

    compute_ionic_currents();

    time            = 0.0;
    timestep_number = 0;

    // Output initial condition.
//    output();
  }

  pcout << "===============================================" << std::endl;

  // Time-stepping loop.
  while (has_non_activated_cells && time < T - 0.5 * delta_t)
    {
      time += delta_t;
      ++timestep_number;

      pcout << "Timestep " << std::setw(3) << timestep_number
            << ", time = " << std::setw(4) << std::fixed << std::setprecision(2)
            << time << " : ";

      integrate_auxiliar_variables();
      compute_ionic_currents();

      assemble();
      solve_linear_system();

      // Perform parallel communication to update the ghost values of the
      // v vector.
      v_ghost = v;

      check_activation_time();

//      output();
    }
  output_activation_time();

}
