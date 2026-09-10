#include "WaveEquation.hpp"

#include <cmath>
#include <chrono>
#include <cstdlib>
#include <limits>
#include <random>
#include <iomanip>
#include <sstream>

namespace
{
constexpr double pi = numbers::PI;

std::string
make_output_name(const unsigned int timestep_number)
{
  std::ostringstream name;
  name << "wave_solution_" << std::setw(4) << std::setfill('0')
       << timestep_number << ".vtu";
  return name.str();
}
} // namespace

WaveEquation::ExactSolution::ExactSolution(const double &omega_)
  : omega(omega_)
{}

double
WaveEquation::ExactSolution::value(const Point<dim> &p,
                                   const unsigned int /*component*/) const
{
  return std::sin(pi * p[0]) * std::sin(pi * p[1]) *
         std::cos(omega * this->get_time());
}

Tensor<1, WaveEquation::dim>
WaveEquation::ExactSolution::gradient(const Point<dim> &p,
                                      const unsigned int /*component*/) const
{
  Tensor<1, dim> grad;
  grad[0] = pi * std::cos(pi * p[0]) * std::sin(pi * p[1]) *
            std::cos(omega * this->get_time());
  grad[1] = pi * std::sin(pi * p[0]) * std::cos(pi * p[1]) *
            std::cos(omega * this->get_time());

  return grad;
}

double
WaveEquation::InitialVelocity::value(const Point<dim> & /*p*/,
                                     const unsigned int /*component*/) const
{
  return 0.0;
}

WaveEquation::WaveEquation(const unsigned int &n_subdivisions_,
                           const unsigned int &fe_degree_,
                           const double       &final_time_,
                           const double       &time_step_,
                           const unsigned int &output_frequency_,
                           const std::string  &output_directory_,
                           const std::string  &validation_case_,
                           const std::string  &solver_name_,
                           const std::string  &preconditioner_name_)
  : n_subdivisions(n_subdivisions_)
  , fe_degree(fe_degree_)
  , final_time(final_time_)
  , time_step(time_step_)
  , output_frequency(output_frequency_)
  , output_directory(output_directory_)
  , validation_case(validation_case_)
  , solver_name(solver_name_)
  , preconditioner_name(preconditioner_name_)
  , dof_handler(mesh)
{}

void
WaveEquation::setup()
{
  std::cout << "===============================================" << std::endl;

  // Mesh.
  //
  // We start with the unit square because the manufactured solution is defined
  // there and no external .msh file is needed. This is the safest first target
  // before adding more geometries.
  {
    std::cout << "Initializing the mesh" << std::endl;
    GridGenerator::subdivided_hyper_cube(mesh,
                                         n_subdivisions,
                                         0.0,
                                         1.0,
                                         /*colorize = */ false);

    std::cout << "  Number of elements = " << mesh.n_active_cells()
              << std::endl;

    std::filesystem::create_directories(output_directory);

    GridOut       grid_out;
    std::ofstream mesh_file(output_directory + "/mesh.vtk");
    grid_out.write_vtk(mesh, mesh_file);
    std::cout << "  Mesh written to " << output_directory + "/mesh.vtk"
              << std::endl;
  }

  std::cout << "-----------------------------------------------" << std::endl;

  // Finite element space.
  //
  // FE_Q is the standard Lagrange element on quadrilateral meshes. The labs use
  // FE_SimplexP when reading triangular meshes; the assembly pattern is the
  // same, only the element and quadrature classes change.
  {
    std::cout << "Initializing the finite element space" << std::endl;

    fe         = std::make_unique<FE_Q<dim>>(fe_degree);
    quadrature = std::make_unique<QGauss<dim>>(fe_degree + 1);

    std::cout << "  Degree                     = " << fe->degree << std::endl;
    std::cout << "  DoFs per cell              = " << fe->dofs_per_cell
              << std::endl;
    std::cout << "  Quadrature points per cell = " << quadrature->size()
              << std::endl;
  }

  std::cout << "-----------------------------------------------" << std::endl;

  // DoF handler.
  {
    std::cout << "Initializing the DoF handler" << std::endl;

    dof_handler.reinit(mesh);
    dof_handler.distribute_dofs(*fe);

    std::cout << "  Number of DoFs = " << dof_handler.n_dofs() << std::endl;
  }

  std::cout << "-----------------------------------------------" << std::endl;

  // Linear algebra structures.
  {
    std::cout << "Initializing matrices and vectors" << std::endl;

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    mass_matrix.reinit(sparsity_pattern);
    stiffness_matrix.reinit(sparsity_pattern);
    newmark_matrix.reinit(sparsity_pattern);
    mass_matrix_constrained.reinit(sparsity_pattern);

    displacement.reinit(dof_handler.n_dofs());
    velocity.reinit(dof_handler.n_dofs());
    acceleration.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
  }

  // Homogeneous Dirichlet boundary values.
  //
  // All boundary faces have boundary id 0 because the mesh was generated with
  // colorize=false. This map is reused whenever we need to constrain a system.
  {
    Functions::ZeroFunction<dim> zero_boundary;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             /*boundary_id = */ 0,
                                             zero_boundary,
                                             boundary_values);
  }
}

void
WaveEquation::assemble_matrices()
{
  std::cout << "===============================================" << std::endl;
  std::cout << "Assembling mass, stiffness and Newmark matrices"
            << std::endl;

  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q           = quadrature->size();

  FEValues<dim> fe_values(*fe,
                          *quadrature,
                          update_values | update_gradients |
                            update_JxW_values);

  FullMatrix<double> cell_mass(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_stiffness(dofs_per_cell, dofs_per_cell);

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  mass_matrix      = 0.0;
  stiffness_matrix = 0.0;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);

      cell_mass      = 0.0;
      cell_stiffness = 0.0;

      for (unsigned int q = 0; q < n_q; ++q)
        {
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                  // Mass matrix:
                  //   M_ij = int phi_j phi_i dx.
                  // This comes from the u_tt term in the weak form.
                  cell_mass(i, j) += fe_values.shape_value(i, q) *
                                     fe_values.shape_value(j, q) *
                                     fe_values.JxW(q);

                  // Stiffness matrix:
                  //   K_ij = int grad(phi_j) dot grad(phi_i) dx.
                  // This comes from integrating -Delta u by parts.
                  cell_stiffness(i, j) +=
                    scalar_product(fe_values.shape_grad(i, q),
                                   fe_values.shape_grad(j, q)) *
                    fe_values.JxW(q);
                }
            }
        }

      cell->get_dof_indices(dof_indices);
      mass_matrix.add(dof_indices, cell_mass);
      stiffness_matrix.add(dof_indices, cell_stiffness);
    }

  // Matrix for Newmark acceleration solve:
  //
  //   (M + beta dt^2 K) a^{n+1} = rhs.
  //
  // The matrix is constant because M and K are constant and dt is fixed, so we
  // assemble it once.
  newmark_matrix.copy_from(mass_matrix);
  newmark_matrix.add(beta * time_step * time_step, stiffness_matrix);

  // Matrix for the initial acceleration solve:
  //
  //   M a(0) = F(0) - K u(0).
  mass_matrix_constrained.copy_from(mass_matrix);

  // Apply homogeneous Dirichlet rows to the two matrices that are directly
  // solved. The original M and K are left untouched for energy and rhs products.
  Vector<double> dummy_solution(dof_handler.n_dofs());
  Vector<double> dummy_rhs(dof_handler.n_dofs());
  MatrixTools::apply_boundary_values(boundary_values,
                                     newmark_matrix,
                                     dummy_solution,
                                     dummy_rhs,
                                     true);
  MatrixTools::apply_boundary_values(boundary_values,
                                     mass_matrix_constrained,
                                     dummy_solution,
                                     dummy_rhs,
                                     true);

  if (std::getenv("WAVE_CHECK_MATRIX") != nullptr)
    {
      check_matrix(mass_matrix_constrained, "mass");
      check_matrix(newmark_matrix, "newmark");
    }
}

void
WaveEquation::initialize_solution()
{
  std::cout << "===============================================" << std::endl;
  std::cout << "Interpolating initial displacement and velocity"
            << std::endl;

  ExactSolution exact_displacement(exact_omega());
  exact_displacement.set_time(0.0);
  VectorTools::interpolate(dof_handler, exact_displacement, displacement);

  InitialVelocity initial_velocity;
  VectorTools::interpolate(dof_handler, initial_velocity, velocity);

  // Boundary values are zero for both displacement and velocity.
  for (const auto &[dof_index, value] : boundary_values)
    {
      displacement[dof_index] = value;
      velocity[dof_index]     = 0.0;
    }

  solve_initial_acceleration();
}

void
WaveEquation::solve_initial_acceleration()
{
  // Initial acceleration comes from the PDE at t=0:
  //
  //   M a(0) + K u(0) = F(0).
  //
  // Therefore we solve M a(0) = F(0) - K u(0).
  Vector<double> force(dof_handler.n_dofs());
  assemble_force(0.0, force);

  Vector<double> stiffness_times_displacement(dof_handler.n_dofs());
  stiffness_matrix.vmult(stiffness_times_displacement, displacement);

  system_rhs = force;
  system_rhs -= stiffness_times_displacement;

  for (const auto &[dof_index, value] : boundary_values)
    {
      acceleration[dof_index] = value;
      system_rhs[dof_index]   = value;
    }

  solve_linear_system(mass_matrix_constrained,
                      acceleration,
                      system_rhs,
                      "initial acceleration");
}

void
WaveEquation::solve_time_step()
{
  // Newmark predictor:
  //
  //   u_predict = u^n + dt v^n + dt^2 (1/2 - beta) a^n.
  //
  // Then solve:
  //
  //   (M + beta dt^2 K) a^{n+1} = F^{n+1} - K u_predict.
  Vector<double> force(dof_handler.n_dofs());
  assemble_force(time + time_step, force);

  Vector<double> displacement_predictor(displacement);
  displacement_predictor.add(time_step, velocity);
  displacement_predictor.add(time_step * time_step * (0.5 - beta),
                             acceleration);

  Vector<double> stiffness_times_predictor(dof_handler.n_dofs());
  stiffness_matrix.vmult(stiffness_times_predictor, displacement_predictor);

  system_rhs = force;
  system_rhs -= stiffness_times_predictor;

  Vector<double> new_acceleration(dof_handler.n_dofs());

  for (const auto &[dof_index, value] : boundary_values)
    {
      new_acceleration[dof_index] = value;
      system_rhs[dof_index]       = value;
    }

  solve_linear_system(newmark_matrix,
                      new_acceleration,
                      system_rhs,
                      "Newmark acceleration");

  // Newmark corrector:
  //
  //   u^{n+1} = u_predict + beta dt^2 a^{n+1}
  //   v^{n+1} = v^n + dt ((1-gamma) a^n + gamma a^{n+1})
  displacement = displacement_predictor;
  displacement.add(beta * time_step * time_step, new_acceleration);

  velocity.add(time_step * (1.0 - gamma), acceleration);
  velocity.add(time_step * gamma, new_acceleration);

  acceleration = new_acceleration;

  for (const auto &[dof_index, value] : boundary_values)
    {
      displacement[dof_index] = value;
      velocity[dof_index]     = 0.0;
      acceleration[dof_index] = 0.0;
    }
}

void
WaveEquation::solve_linear_system(const SparseMatrix<double> &matrix,
                                  Vector<double>             &unknown,
                                  const Vector<double>       &rhs,
                                  const std::string          &label) const
{
  ReductionControl solver_control(/* maxiter = */ 5000,
                                  /* tolerance = */ 1.0e-14,
                                  /* reduce = */ 1.0e-10);

  using Clock = std::chrono::steady_clock;
  double setup_seconds = 0.0;
  double solve_seconds = 0.0;
  bool converged = true;
  const auto solve_with = [&](auto &preconditioner) {
    const auto setup_start = Clock::now();
    preconditioner.initialize(matrix); // Both defaults use relaxation 1.
    setup_seconds = std::chrono::duration<double>(Clock::now() - setup_start).count();
    const auto solve_start = Clock::now();
    try
      {
        if (solver_name == "cg")
          {
            SolverCG<Vector<double>> solver(solver_control);
            solver.solve(matrix, unknown, rhs, preconditioner);
          }
        else
          {
            // Right preconditioning makes GMRES test the original residual,
            // as CG does. 32 temporary vectors give a restart basis of 30.
            SolverGMRES<Vector<double>> solver(
              solver_control,
              SolverGMRES<Vector<double>>::AdditionalData(32, true, true));
            solver.solve(matrix, unknown, rhs, preconditioner);
          }
      }
    catch (const SolverControl::NoConvergence &)
      {
        converged = false;
      }
    solve_seconds = std::chrono::duration<double>(Clock::now() - solve_start).count();
  };

  if (preconditioner_name == "ssor")
    {
      PreconditionSSOR<SparseMatrix<double>> preconditioner;
      solve_with(preconditioner);
    }
  else
    {
      PreconditionJacobi<SparseMatrix<double>> preconditioner;
      solve_with(preconditioner);
    }

  Vector<double> residual(rhs.size());
  matrix.vmult(residual, unknown);
  residual -= rhs;
  const double norm = residual.l2_norm();
  const double target = std::max(1.e-14, 1.e-10 * rhs.l2_norm());
  converged = converged && std::isfinite(norm) && norm <= 1.01 * target;
  std::ofstream statistics(output_directory + "/solver_stats.csv", std::ios::app);
  statistics << std::setprecision(16)
             << (label == "initial acceleration" ? 0 : timestep_number + 1)
             << ',' << label << ',' << solver_name << ',' << preconditioner_name
             << ',' << converged << ',' << solver_control.last_step()
             << ',' << norm << ',' << norm / std::max(rhs.l2_norm(), 1.e-300)
             << ',' << setup_seconds << ',' << solve_seconds << '\n';
  AssertThrow(converged, ExcMessage("Linear solve failed its true residual check"));
  std::cout << "  Solved " << label << " in " << solver_control.last_step()
            << ' ' << solver_name << " iterations" << std::endl;
}

void
WaveEquation::check_matrix(const SparseMatrix<double> &matrix,
                           const std::string &label) const
{
  double asymmetry = 0.0;
  double scale = 0.0;
  double min_diagonal = std::numeric_limits<double>::max();
  for (unsigned int i = 0; i < matrix.m(); ++i)
    {
      min_diagonal = std::min(min_diagonal, matrix.diag_element(i));
      for (auto entry = matrix.begin(i); entry != matrix.end(i); ++entry)
        {
          scale = std::max(scale, std::abs(entry->value()));
          asymmetry = std::max(asymmetry,
            std::abs(entry->value() - matrix.el(entry->column(), i)));
        }
    }
  std::mt19937 generator(42);
  std::uniform_real_distribution<double> distribution(-1.0, 1.0);
  Vector<double> x(matrix.m()), ax(matrix.m());
  double min_rayleigh = std::numeric_limits<double>::max();
  for (unsigned int trial = 0; trial < 8; ++trial)
    {
      for (auto &value : x)
        value = distribution(generator);
      matrix.vmult(ax, x);
      min_rayleigh = std::min(min_rayleigh, (x * ax) / (x * x));
    }
  std::ofstream checks(output_directory + "/matrix_checks.csv", std::ios::app);
  checks << std::setprecision(16) << label << ',' << matrix.m() << ','
         << asymmetry / scale << ',' << min_diagonal << ',' << min_rayleigh << '\n';
  AssertThrow(asymmetry <= 1.e-12 * scale && min_diagonal > 0 && min_rayleigh > 0,
              ExcMessage("Matrix sanity check failed (not a proof of SPD)"));
}

void
WaveEquation::output_results()
{
  DataOut<dim> data_out;

  data_out.add_data_vector(dof_handler, displacement, "u");
  data_out.add_data_vector(dof_handler, velocity, "velocity");
  data_out.add_data_vector(dof_handler, acceleration, "acceleration");

  data_out.build_patches(fe->degree);

  const std::string file_name = make_output_name(timestep_number);
  const std::string full_path = output_directory + "/" + file_name;

  std::ofstream output_file(full_path);
  data_out.write_vtu(output_file);

  pvd_entries.emplace_back(time, file_name);
  write_pvd_record();

  std::cout << "  Output written to " << full_path << std::endl;
}

void
WaveEquation::write_pvd_record() const
{
  std::ofstream pvd_file(output_directory + "/wave_solution.pvd");
  pvd_file << std::setprecision(16);

  pvd_file << "<?xml version=\"1.0\"?>\n";
  pvd_file << "<VTKFile type=\"Collection\" version=\"0.1\" "
              "byte_order=\"LittleEndian\">\n";
  pvd_file << "  <Collection>\n";

  for (const auto &[entry_time, file_name] : pvd_entries)
    {
      pvd_file << "    <DataSet timestep=\"" << entry_time
               << "\" group=\"\" part=\"0\" file=\"" << file_name
               << "\"/>\n";
    }

  pvd_file << "  </Collection>\n";
  pvd_file << "</VTKFile>\n";
}

std::pair<double, double>
WaveEquation::compute_errors() const
{
  ExactSolution exact_solution(exact_omega());
  exact_solution.set_time(time);

  const QGauss<dim> quadrature_error(fe_degree + 2);

  Vector<double> error_per_cell(mesh.n_active_cells());

  VectorTools::integrate_difference(dof_handler,
                                    displacement,
                                    exact_solution,
                                    error_per_cell,
                                    quadrature_error,
                                    VectorTools::L2_norm);
  const double l2_error = VectorTools::compute_global_error(mesh,
                                                           error_per_cell,
                                                           VectorTools::L2_norm);

  VectorTools::integrate_difference(dof_handler,
                                    displacement,
                                    exact_solution,
                                    error_per_cell,
                                    quadrature_error,
                                    VectorTools::H1_seminorm);
  const double h1_error =
    VectorTools::compute_global_error(mesh,
                                      error_per_cell,
                                      VectorTools::H1_seminorm);

  return {l2_error, h1_error};
}

double
WaveEquation::compute_energy() const
{
  Vector<double> matrix_times_vector(dof_handler.n_dofs());

  mass_matrix.vmult(matrix_times_vector, velocity);
  const double kinetic_energy = 0.5 * (velocity * matrix_times_vector);

  stiffness_matrix.vmult(matrix_times_vector, displacement);
  const double potential_energy = 0.5 * (displacement * matrix_times_vector);

  return kinetic_energy + potential_energy;
}

void
WaveEquation::append_diagnostics()
{
  const auto [l2_error, h1_error] = compute_errors();
  const double energy             = compute_energy();

  const bool write_header =
    !std::filesystem::exists(output_directory + "/diagnostics.csv");

  std::ofstream diagnostics(output_directory + "/diagnostics.csv",
                            std::ios::app);

  if (write_header)
    diagnostics << "step,time,l2_error,h1_seminorm,energy\n";

  diagnostics << timestep_number << "," << std::setprecision(16) << time
              << "," << l2_error << "," << h1_error << "," << energy
              << "\n";

  std::cout << "  L2 error = " << l2_error
            << ", H1 seminorm error = " << h1_error
            << ", energy = " << energy << std::endl;
}

void
WaveEquation::append_probe()
{
  const Point<dim> centre(0.5, 0.5);
  const double numerical_value =
    VectorTools::point_value(dof_handler, displacement, centre);
  ExactSolution exact_solution(exact_omega());
  exact_solution.set_time(time);
  const double exact_value = exact_solution.value(centre);

  const bool write_header =
    !std::filesystem::exists(output_directory + "/probe.csv");
  std::ofstream probe(output_directory + "/probe.csv", std::ios::app);
  if (write_header)
    probe << "step,time,numerical,exact,absolute_error\n";
  probe << timestep_number << ',' << std::setprecision(16) << time << ','
        << numerical_value << ',' << exact_value << ','
        << std::abs(numerical_value - exact_value) << '\n';
}

double
WaveEquation::exact_omega() const
{
  if (validation_case == "forced")
    return pi;

  return std::sqrt(2.0) * pi;
}

void
WaveEquation::assemble_force(const double &time, Vector<double> &force) const
{
  force = 0.0;

  if (validation_case == "free")
    return;

  // Manufactured forcing:
  //
  //   u(x,y,t) = sin(pi x) sin(pi y) cos(omega t)
  //
  // Substitution gives:
  //
  //   f = u_tt - Delta u
  //     = (2 pi^2 - omega^2) sin(pi x) sin(pi y) cos(omega t).
  //
  // For "free", omega = sqrt(2) pi, so the coefficient is zero.
  // For "forced", omega = pi, so the coefficient is pi^2.
  const double omega             = exact_omega();
  const double forcing_coefficient = 2.0 * pi * pi - omega * omega;

  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q           = quadrature->size();

  FEValues<dim> fe_values(*fe,
                          *quadrature,
                          update_values | update_quadrature_points |
                            update_JxW_values);

  Vector<double> cell_force(dofs_per_cell);
  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);
      cell_force = 0.0;

      for (unsigned int q = 0; q < n_q; ++q)
        {
          const Point<dim> &p = fe_values.quadrature_point(q);
          const double force_value =
            forcing_coefficient * std::sin(pi * p[0]) * std::sin(pi * p[1]) *
            std::cos(omega * time);

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            cell_force(i) += force_value * fe_values.shape_value(i, q) *
                             fe_values.JxW(q);
        }

      cell->get_dof_indices(dof_indices);
      force.add(dof_indices, cell_force);
    }
}

void
WaveEquation::run()
{
  std::cout << "Validation case = " << validation_case << std::endl;

  setup();
  std::ofstream(output_directory + "/solver_stats.csv")
    << "step,system,solver,preconditioner,converged,iterations,residual,relative_residual,setup_seconds,solve_seconds\n";
  std::ofstream(output_directory + "/matrix_checks.csv")
    << "matrix,size,relative_asymmetry,min_diagonal,min_sampled_rayleigh\n";
  std::ofstream parameters(output_directory + "/parameters.csv");
  parameters << "case,N,degree,cells,dofs,h,dt,T,solver,preconditioner\n"
             << std::setprecision(16) << validation_case << ',' << n_subdivisions
             << ',' << fe_degree << ',' << mesh.n_active_cells() << ','
             << dof_handler.n_dofs() << ',' << std::sqrt(2.0) / n_subdivisions
             << ',' << time_step << ',' << final_time << ',' << solver_name
             << ',' << preconditioner_name << '\n';
  assemble_matrices();
  initialize_solution();

  time            = 0.0;
  timestep_number = 0;

  std::filesystem::remove(output_directory + "/diagnostics.csv");
  std::filesystem::remove(output_directory + "/probe.csv");

  std::cout << "===============================================" << std::endl;
  std::cout << "Starting time loop" << std::endl;

  output_results();
  append_diagnostics();
  append_probe();

  const auto n_steps = static_cast<unsigned int>(std::llround(final_time / time_step));
  while (timestep_number < n_steps)
    {
      std::cout << "-----------------------------------------------" << std::endl;

      solve_time_step();

      ++timestep_number;
      time = timestep_number * time_step;

      std::cout << "Timestep " << timestep_number << ", time = " << time
                << std::endl;

      append_diagnostics();
      append_probe();

      if (timestep_number % output_frequency == 0 ||
          time >= final_time - 0.5 * time_step)
        output_results();
    }

  std::cout << "===============================================" << std::endl;
  std::cout << "Simulation completed" << std::endl;
}
