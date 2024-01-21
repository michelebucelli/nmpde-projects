#include "Problem.hpp"

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/mapping_fe.h>

double Cylinder2D::InletVelocity::value(const Point<dim> &p,
                                        const unsigned int component) const {
  if (component == 0) {
    return 4.0 * U_m * p[1] * (H - p[1]) / (H * H);
  } else {
    return 0.0;
  }
}

double Cylinder3D::InletVelocity::value(const Point<dim> &p,
                                        const unsigned int component) const {
  if (component == 0) {
    return 16.0 * U_m * p[1] * p[2] * (H - p[1]) * (H - p[2]) / (H * H * H * H);
  } else {
    return 0.0;
  }
}

void Cylinder2D::InletVelocity::vector_value(const Point<dim> &p,
                                             Vector<double> &values) const {
  values[0] = value(p, 0);

  for (unsigned int i = 1; i < dim + 1; ++i) values[i] = 0.0;
}

void Cylinder3D::InletVelocity::vector_value(const Point<dim> &p,
                                             Vector<double> &values) const {
  values[0] = value(p, 0);

  for (unsigned int i = 1; i < dim + 1; ++i) values[i] = 0.0;
}

double Cylinder2D::get_reynolds_number() const {
  return 2.0 * inlet_velocity.value(Point<dim>(0, H / 2.0), 0) / 3.0 * D / nu;
}

double Cylinder3D::get_reynolds_number() const {
  return 4.0 * inlet_velocity.value(Point<dim>(0, H / 2.0, H / 2.0), 0) / 9.0 *
         D / nu;
}

Cylinder2D::Cylinder2D(const std::string &mesh_file_name_,
                       const unsigned int &degree_velocity_,
                       const unsigned int &degree_pressure_, const double &T_,
                       const double &deltat_)
    : NavierStokes(mesh_file_name_, degree_velocity_, degree_pressure_, T_,
                   deltat_),
      inlet_velocity(U_m, H),
      zero_function(dim + 1) {
  ro = 1.0;
  nu = 1e-3;

  initial_conditions = std::make_shared<Functions::ZeroFunction<dim>>(dim + 1);

  dirichlet_boundary_functions[0] = &inlet_velocity;
  dirichlet_boundary_functions[1] = &zero_function;

  neumann_boundary_functions[2] = p_out;
}

Cylinder3D::Cylinder3D(const std::string &mesh_file_name_,
                       const unsigned int &degree_velocity_,
                       const unsigned int &degree_pressure_, const double &T_,
                       const double &deltat_)
    : NavierStokes(mesh_file_name_, degree_velocity_, degree_pressure_, T_,
                   deltat_),
      inlet_velocity(U_m, H),
      zero_function(dim + 1) {
  ro = 1.0;
  nu = 1e-3;

  initial_conditions = std::make_shared<Functions::ZeroFunction<dim>>(dim + 1);

  dirichlet_boundary_functions[4] = &inlet_velocity;

  neumann_boundary_functions[2] = p_out;

  dirichlet_boundary_functions[1] = &zero_function;
  dirichlet_boundary_functions[3] = &zero_function;
  dirichlet_boundary_functions[5] = &zero_function;
}

EthierSteinman::EthierSteinman(const std::string &mesh_file_name_,
                               const unsigned int &degree_velocity_,
                               const unsigned int &degree_pressure_,
                               const double &T_, const double &deltat_,
                               const double &nu_)
    : NavierStokes(mesh_file_name_, degree_velocity_, degree_pressure_, T_,
                   deltat_),
      exact_solution(nu_) {
  ro = 1.0;
  nu = nu_;

  initial_conditions = std::make_shared<ExactSolution>(nu);

  for (unsigned int i = 1U; i <= 6U; i++) {
    dirichlet_boundary_functions[i] = &exact_solution.exact_velocity;
  }
}

double EthierSteinman::ExactSolution::ExactVelocity::value(
    const Point<dim> &p, const unsigned int component) const {
  double result = -a * std::exp(-nu * b * b * get_time());
  if (component == 0) {
    return result * std::exp(a * p[0]) * std::sin(a * p[1] + b * p[2]) +
           std::exp(a * p[2]) * std::cos(a * p[0] + b * p[1]);
  } else if (component == 1) {
    return result * std::exp(a * p[1]) * std::sin(a * p[2] + b * p[0]) +
           std::exp(a * p[0]) * std::cos(a * p[1] + b * p[2]);
  } else if (component == 2) {
    return result * std::exp(a * p[2]) * std::sin(a * p[0] + b * p[1]) +
           std::exp(a * p[1]) * std::cos(a * p[2] + b * p[0]);
  } else {
    return 0.0;
  }
}

void EthierSteinman::ExactSolution::ExactVelocity::vector_value(
    const Point<dim> &p, Vector<double> &values) const {
  for (unsigned int i = 0; i < dim; i++) {
    values[i] = value(p, i);
  }
  values[dim] = 0.0;
}

double EthierSteinman::ExactSolution::ExactPressure::value(
    const Point<dim> &p, const unsigned int /*component*/) const {
  return -a * a / 2.0 * std::exp(-2 * nu * b * b * get_time()) *
         (2.0 * std::sin(a * p[0] + b * p[1]) * std::cos(a * p[2] + b * p[0]) *
              std::exp(a * (p[1] + p[2])) +
          2.0 * std::sin(a * p[1] + b * p[2]) * std::cos(a * p[0] + b * p[1]) *
              std::exp(a * (p[0] + p[2])) +
          2.0 * std::sin(a * p[2] + b * p[0]) * std::cos(a * p[1] + b * p[2]) *
              std::exp(a * (p[0] + p[1])) +
          std::exp(2.0 * a * p[0]) + std::exp(2.0 * a * p[1]) +
          std::exp(2.0 * a * p[2]));
}

void EthierSteinman::ExactSolution::ExactPressure::vector_value(
    const Point<dim> &p, Vector<double> &values) const {
  values[dim] = value(p);
  for (unsigned int i = 0; i < dim; i++) {
    values[i] = 0.0;
  }
}

double EthierSteinman::ExactSolution::value(
    const Point<dim> &p, const unsigned int component) const {
  if (component < dim) {
    return exact_velocity.value(p, component);
  } else if (component == dim) {
    return exact_pressure.value(p);
  } else {
    return 0.0;
  }
}

void EthierSteinman::ExactSolution::vector_value(const Point<dim> &p,
                                                 Vector<double> &values) const {
  for (unsigned int i = 0; i < dim + 1; i++) {
    values[i] = value(p, i);
  }
}

double EthierSteinman::compute_error(const VectorTools::NormType &norm_type) {
  FE_SimplexP<dim> fe_mapping(1);
  MappingFE mapping(fe_mapping);

  // The error is an integral, and we approximate that integral using a
  // quadrature formula. To make sure we are accurate enough, we use a
  // quadrature formula with one node more than what we used in assembly.
  const QGaussSimplex<dim> quadrature_error(degree_velocity + 2);

  // Set the time for the exact solution.
  exact_solution.set_time(time_step * deltat);

  // First we compute the norm on each element, and store it in a vector.
  Vector<double> error_per_cell(mesh.n_active_cells());
  VectorTools::integrate_difference(mapping, dof_handler, solution,
                                    exact_solution, error_per_cell,
                                    quadrature_error, norm_type);

  // Then, we add out all the cells.
  const double error =
      VectorTools::compute_global_error(mesh, error_per_cell, norm_type);

  return error;
}

void EthierSteinman::apply_initial_conditions() {
  NavierStokes<dim>::apply_initial_conditions();
  pcout << "L2 error: " << compute_error(VectorTools::L2_norm) << std::endl;
}

void EthierSteinman::solve_time_step() {
  NavierStokes<dim>::solve_time_step();
  pcout << "L2 error: " << compute_error(VectorTools::L2_norm) << std::endl;
}