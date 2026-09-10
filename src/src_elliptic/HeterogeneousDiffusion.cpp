#include "HeterogeneousDiffusion.hpp"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <filesystem>

namespace
{
  // Shared width for the dashed separators printed around the results
  // table, so the footer under "Setup/Assembly time" and the separator
  // under the "Precond | Iters | ..." header always line up.
  constexpr std::size_t results_table_width = 74;
}

std::vector<HeterogeneousDiffusion::Ball>
HeterogeneousDiffusion::make_balls()
{
  return {{{0.18, 0.20, 0.23}, 0.090},
          {{0.34, 0.28, 0.72}, 0.075},
          {{0.52, 0.18, 0.46}, 0.080},
          {{0.76, 0.24, 0.32}, 0.070},
          {{0.24, 0.52, 0.55}, 0.085},
          {{0.48, 0.48, 0.22}, 0.070},
          {{0.68, 0.56, 0.78}, 0.095},
          {{0.86, 0.62, 0.50}, 0.060},
          {{0.30, 0.78, 0.30}, 0.075},
          {{0.55, 0.82, 0.60}, 0.090},
          {{0.78, 0.84, 0.20}, 0.070},
          {{0.14, 0.86, 0.82}, 0.065}};
}

HeterogeneousDiffusion::HeterogeneousDiffusion(
  const unsigned int refinements_,
  const double       p_,
  const MPI_Comm     mpi_communicator_,
  const bool         verify_mode_)
  : refinements(refinements_)
  , p(p_)
  , mpi_communicator(mpi_communicator_)
  , mpi_size(Utilities::MPI::n_mpi_processes(mpi_communicator_))
  , mpi_rank(Utilities::MPI::this_mpi_process(mpi_communicator_))
  , verify_mode(verify_mode_)
  , pcout(std::cout, mpi_rank == 0)
  , balls(make_balls())
  , mesh(mpi_communicator_)
  , fe(1)
  , mapping(FE_SimplexP<dim>(1))
  , dof_handler(mesh)
{}

void
HeterogeneousDiffusion::setup_and_assemble()
{
  {
    std::ostringstream header;
    header << "Case:  p = " << p << "   |  refinements = " << refinements
           << "   |  MPI ranks = " << mpi_size
           << "   |  balls = " << balls.size();

    pcout << std::string(header.str().size(), '=') << std::endl;
    pcout << header.str() << std::endl;
    pcout << std::string(header.str().size(), '=') << std::endl;
  }

  Timer timer;

  // Build a tetrahedral mesh. The parameter `refinements` means
  // h = 2^{-r}; each Cartesian subdivision is split into six tetrahedra.
  {
    Triangulation<dim> serial_tet_mesh;

    const unsigned int subdivisions = 1u << refinements; // "<<" is a bit-wise left shift operator, equivalent to 2^refinements

    GridGenerator::subdivided_hyper_cube_with_simplices(
      serial_tet_mesh, subdivisions, 0.0, 1.0);

    // Set subdomain ids on the serial tetrahedral mesh. These ids are then
    // used when constructing the fully distributed triangulation.
    GridTools::partition_triangulation(mpi_size, serial_tet_mesh);

    const auto construction_data =
      TriangulationDescription::Utilities::create_description_from_triangulation(
        serial_tet_mesh, mpi_communicator);

    mesh.create_triangulation(construction_data);
  }

  dof_handler.distribute_dofs(fe);

  locally_owned_dofs = dof_handler.locally_owned_dofs();
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

  // Set up constraints for Dirichlet boundary conditions.
  {
    constraints.clear();
    constraints.reinit(locally_relevant_dofs);
    VectorTools::interpolate_boundary_values(mapping,
                                            dof_handler,
                                            0, // id of the boundary (0 = all boundaries)
                                            Functions::ZeroFunction<dim>(),
                                            constraints);
    constraints.close();
  }

  sparsity_pattern.reinit(locally_owned_dofs,
                          locally_owned_dofs,
                          locally_relevant_dofs,
                          mpi_communicator);
  DoFTools::make_sparsity_pattern(
    dof_handler, sparsity_pattern, constraints, false, mpi_rank);
  sparsity_pattern.compress();
  system_matrix.reinit(sparsity_pattern);

  system_rhs.reinit(locally_owned_dofs, mpi_communicator);
  solution.reinit(locally_owned_dofs, mpi_communicator);

  pcout << std::left << std::setw(22) << "  Cells"
        << mesh.n_global_active_cells() << std::endl;
  pcout << std::left << std::setw(22) << "  DoFs"
        << dof_handler.n_dofs() << std::endl;

  timer.stop();
  const double setup_time =
    Utilities::MPI::max(timer.wall_time(), mpi_communicator);

  timer.restart();

  const QGaussSimplex<dim> quadrature(fe.degree + 1);
  FEValues<dim> fe_values(mapping,
                          fe,
                          quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
  const unsigned int n_q           = quadrature.size();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  DiffusionCoefficient      heterogeneous_diffusion(p, balls);
  RightHandSide             heterogeneous_rhs;
  ManufacturedRightHandSide manufactured_rhs;

  const Function<dim> *rhs_ptr = verify_mode
    ? static_cast<const Function<dim> *>(&manufactured_rhs)
    : static_cast<const Function<dim> *>(&heterogeneous_rhs);

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
          const double mu = verify_mode
            ? 1.0
            : heterogeneous_diffusion.value(fe_values.quadrature_point(q));

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              cell_rhs(i) += rhs_ptr->value(fe_values.quadrature_point(q)) *
                             fe_values.shape_value(i, q) * fe_values.JxW(q);

              for (unsigned int j = 0; j < dofs_per_cell; ++j)

                cell_matrix(i, j) += mu * fe_values.shape_grad(i, q) *
                                     fe_values.shape_grad(j, q) *
                                     fe_values.JxW(q);
            }
        }

      cell->get_dof_indices(local_dof_indices);
      // Distribute the Dirichlet boundary conditions and assemble the local contributions into the global system.
      constraints.distribute_local_to_global(
        cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
    }

  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);

  timer.stop();
  const double assembly_time =
    Utilities::MPI::max(timer.wall_time(), mpi_communicator);

  static constexpr std::size_t table_width = results_table_width;

  std::ostringstream setup_line, assembly_line;
  setup_line    << std::left << std::setw(22) << "  Setup time [s]" << setup_time;
  assembly_line << std::left << std::setw(22) << "  Assembly time [s]" << assembly_time;

  pcout << setup_line.str() << std::endl;
  pcout << assembly_line.str() << std::endl;
  pcout << std::string(table_width, '-') << std::endl;
}

template <class Preconditioner>
unsigned int
HeterogeneousDiffusion::solve_with(const Preconditioner &preconditioner,
                                   const std::string    &name,
                                   const double          setup_time)
{
  SolverControl solver_control(20000, 1e-10 * system_rhs.l2_norm());
  SolverCG<TrilinosWrappers::MPI::Vector> solver(solver_control);

  double condition_number = -1.0;
  solver.connect_condition_number_slot(
    [&condition_number](const double cn) { condition_number = cn; });

  Timer timer;
  solver.solve(system_matrix, solution, system_rhs, preconditioner);
  timer.stop();

  const double solve_time =
    Utilities::MPI::max(timer.wall_time(), mpi_communicator);
  const double setup_time_global =
    Utilities::MPI::max(setup_time, mpi_communicator);

  if (!result_header_printed)
    {
      std::ostringstream header;
      header << std::left << std::setw(10) << "Precond" << std::right
             << std::setw(8)  << "Iters"
             << std::setw(14) << "SetupT[s]"
             << std::setw(14) << "SolveT[s]"
             << std::setw(14) << "TotalT[s]"
             << std::setw(14) << "CondNumber";

      pcout << header.str() << std::endl;
      pcout << std::string(results_table_width, '-') << std::endl;
      result_header_printed = true;
    }

  pcout << std::left << std::setw(10) << name << std::right
        << std::setw(8)  << solver_control.last_step()
        << std::setw(14) << std::scientific << std::setprecision(3)
        << setup_time_global
        << std::setw(14) << solve_time
        << std::setw(14) << setup_time_global + solve_time
        << std::setw(14) << condition_number
        << std::defaultfloat << std::endl;

  log_result_csv(
    name, solver_control.last_step(), setup_time_global, solve_time,
    condition_number);

  return solver_control.last_step();
}

void
HeterogeneousDiffusion::log_result_csv(const std::string &name,
                                       const unsigned int iterations,
                                       const double       setup_time,
                                       const double       solve_time,
                                       const double condition_number) const
{
  if (mpi_rank != 0)
    return;

  const std::string filename    = "results.csv";
  const bool        file_exists = std::filesystem::exists(filename);

  std::ofstream csv(filename, std::ios::app);
  if (!file_exists)
    csv << "refinements,p,mpi_procs,cells,dofs,preconditioner,iterations,"
           "setup_time,solve_time,total_time,condition_number"
        << std::endl;

  csv << refinements << "," << p << "," << mpi_size << ","
      << mesh.n_global_active_cells() << "," << dof_handler.n_dofs() << ","
      << name << "," << iterations << "," << setup_time << "," << solve_time
      << "," << (setup_time + solve_time) << "," << condition_number
      << std::endl;
}

void
HeterogeneousDiffusion::solve(const std::string &name)
{
  solution = 0.0;

  try
    {
      if (name == "none")
        {
          solve_with(PreconditionIdentity(), name, 0.0);
          return;
        }

      Timer timer;

      if (name == "jacobi")
        {
          TrilinosWrappers::PreconditionJacobi preconditioner;
          preconditioner.initialize(system_matrix);
          timer.stop();
          solve_with(preconditioner, name, timer.wall_time());
          return;
        }

      if (name == "ssor")
        {
          TrilinosWrappers::PreconditionSSOR preconditioner;
          preconditioner.initialize(system_matrix);
          timer.stop();
          solve_with(preconditioner, name, timer.wall_time());
          return;
        }

      if (name == "ilu")
        {
          TrilinosWrappers::PreconditionILU preconditioner;
          preconditioner.initialize(system_matrix);
          timer.stop();
          solve_with(preconditioner, name, timer.wall_time());
          return;
        }

      if (name == "amg")
        {
          TrilinosWrappers::PreconditionAMG                 preconditioner;
          TrilinosWrappers::PreconditionAMG::AdditionalData data;
          data.elliptic              = true;
          data.higher_order_elements = false;
          data.n_cycles              = 1; // default 1 for benchmark
          preconditioner.initialize(system_matrix, data);
          timer.stop();
          solve_with(preconditioner, name, timer.wall_time());
          return;
        }

      pcout << "failed " << name << " unknown preconditioner" << std::endl;
    }
  catch (const std::exception &exception)
    {
      pcout << "failed " << name << " " << exception.what() << std::endl;
    }
}

double
HeterogeneousDiffusion::compute_error(
  const VectorTools::NormType &norm_type) const
{
  TrilinosWrappers::MPI::Vector solution_ghost(locally_owned_dofs,
                                               locally_relevant_dofs,
                                               mpi_communicator);
  solution_ghost = solution;

  const QGaussSimplex<dim> quadrature_error(fe.degree + 2);

  ManufacturedSolution manufactured_solution;

  Vector<double> error_per_cell(mesh.n_active_cells());
  VectorTools::integrate_difference(mapping,
                                    dof_handler,
                                    solution_ghost,
                                    manufactured_solution,
                                    error_per_cell,
                                    quadrature_error,
                                    norm_type);

  return VectorTools::compute_global_error(mesh, error_per_cell, norm_type);
}

void
HeterogeneousDiffusion::output(const std::string &tag) const
{
  TrilinosWrappers::MPI::Vector solution_ghost(locally_owned_dofs,
                                               locally_relevant_dofs,
                                               mpi_communicator);

  const std::string folder = "output/";
  if (mpi_rank == 0)
    std::filesystem::create_directories(folder);
  MPI_Barrier(mpi_communicator);
  
  
  solution_ghost = solution;
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution_ghost, "solution");

  
  Vector<double>       mu_field(mesh.n_active_cells());
  DiffusionCoefficient diffusion(p, balls);
  for (const auto &cell : mesh.active_cell_iterators())
    if (cell->is_locally_owned())
      mu_field[cell->active_cell_index()] = diffusion.value(cell->center());
  data_out.add_data_vector(mu_field, "mu");

  data_out.build_patches(mapping, fe.degree);

  const std::string basename = "heterogeneous-" + tag + "-p" +
                               std::to_string(int(p)) + "-ref" +
                               std::to_string(refinements);
  data_out.write_vtu_with_pvtu_record(folder, basename, 0, mpi_communicator);
  pcout << "vtk " << folder << basename << ".pvtu" << std::endl;
}
