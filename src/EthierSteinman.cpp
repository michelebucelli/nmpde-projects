#include "EthierSteinman.hpp"

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/mapping_fe.h>

EthierSteinman::EthierSteinman(const std::string &mesh_file_name_,
                               const unsigned int &degree_velocity_,
                               const unsigned int &degree_pressure_,
                               const double &T_, const double &deltat_,
                               const SolverOptions &solver_options_,
                               const double &nu_)
    : NavierStokes(mesh_file_name_, degree_velocity_, degree_pressure_, T_,
                   deltat_, solver_options_),
      exact_solution(nu_),
      neumann_function(nu_) {
  pcout << "Solving Ethier-Steinman problem" << std::endl;
  nu = nu_;
  ro = 1.0;
  initial_conditions = std::make_shared<ExactSolution>(nu);
  initial_conditions->set_time(0.0);

  dirichlet_boundary_functions[1U] = &exact_solution.exact_velocity;
  for (unsigned int i = 3U; i <= 6U; i++) {
    dirichlet_boundary_functions[i] = &exact_solution.exact_velocity;
  }

  neumann_boundary_functions[2U] = &neumann_function;
}

double EthierSteinman::ExactSolution::ExactVelocity::value(
    const Point<dim> &p, const unsigned int component) const {
  const double result = -a * std::exp(-nu * b * b * get_time());
  if (component == 0) {
    return result * (std::exp(a * p[0]) * std::sin(a * p[1] + b * p[2]) +
                     std::exp(a * p[2]) * std::cos(a * p[0] + b * p[1]));
  } else if (component == 1) {
    return result * (std::exp(a * p[1]) * std::sin(a * p[2] + b * p[0]) +
                     std::exp(a * p[0]) * std::cos(a * p[1] + b * p[2]));
  } else if (component == 2) {
    return result * (std::exp(a * p[2]) * std::sin(a * p[0] + b * p[1]) +
                     std::exp(a * p[1]) * std::cos(a * p[2] + b * p[0]));
  } else {
    return 0.0;
  }
}

void EthierSteinman::ExactSolution::ExactVelocity::vector_value(
    const Point<dim> &p, Vector<double> &values) const {
  for (unsigned int i = 0; i < dim + 1; i++) {
    values[i] = value(p, i);
  }
}

Tensor<1, EthierSteinman::dim>
EthierSteinman::ExactSolution::ExactVelocity::gradient(
    const Point<dim> &p, const unsigned int component) const {
  Tensor<1, dim> result;

  for (unsigned int i = 0; i < dim; i++) {
    result[i] = -a * std::exp(-nu * b * b * get_time());
  }

  if (component == 0) {
    result[0] *= (a * std::exp(a * p[0]) * std::sin(a * p[1] + b * p[2]) -
                  a * std::exp(a * p[2]) * std::sin(a * p[0] + b * p[1]));
    result[1] *= (a * std::exp(a * p[0]) * std::cos(a * p[1] + b * p[2]) -
                  b * std::exp(a * p[2]) * std::sin(a * p[0] + b * p[1]));
    result[2] *= (b * std::exp(a * p[0]) * std::cos(a * p[1] + b * p[2]) +
                  a * std::exp(a * p[2]) * std::cos(a * p[0] + b * p[1]));
  } else if (component == 1) {
    result[0] *= (b * std::exp(a * p[1]) * std::cos(a * p[2] + b * p[0]) +
                  a * std::exp(a * p[0]) * std::cos(a * p[1] + b * p[2]));
    result[1] *= (a * std::exp(a * p[1]) * std::sin(a * p[2] + b * p[0]) -
                  a * std::exp(a * p[0]) * std::sin(a * p[1] + b * p[2]));
    result[2] *= (a * std::exp(a * p[1]) * std::cos(a * p[2] + b * p[0]) -
                  b * std::exp(a * p[0]) * std::sin(a * p[1] + b * p[2]));
  } else if (component == 2) {
    result[0] *= (a * std::exp(a * p[2]) * std::cos(a * p[0] + b * p[1]) -
                  b * std::exp(a * p[1]) * std::sin(a * p[2] + b * p[0]));
    result[1] *= (b * std::exp(a * p[2]) * std::cos(a * p[0] + b * p[1]) +
                  a * std::exp(a * p[1]) * std::cos(a * p[2] + b * p[0]));
    result[2] *= (a * std::exp(a * p[2]) * std::sin(a * p[0] + b * p[1]) -
                  a * std::exp(a * p[1]) * std::sin(a * p[2] + b * p[0]));
  } else {
    for (unsigned int i = 0; i < dim; i++) {
      result[i] = 0.0;
    }
  }

  return result;
}

void EthierSteinman::ExactSolution::ExactVelocity::vector_gradient(
    const Point<dim> &p, std::vector<Tensor<1, dim>> &values) const {
  for (unsigned int i = 0; i < dim + 1; i++) {
    values[i] = gradient(p, i);
  }
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
  values[0] = value(p);
}

double EthierSteinman::ExactSolution::value(
    const Point<dim> &p, const unsigned int component) const {
  // Set the time for the exact velocity and pressure.
  exact_velocity.set_time(get_time());
  exact_pressure.set_time(get_time());

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
  // Set the time for the exact velocity and pressure.
  exact_velocity.set_time(get_time());
  exact_pressure.set_time(get_time());

  for (unsigned int i = 0; i < dim + 1; i++) {
    values[i] = value(p, i);
  }
}

double EthierSteinman::NeumannFunction::value(
    const Point<dim> &p, const unsigned int component) const {
  // Set the time for the exact solution.
  exact_solution.set_time(get_time());

  // Calculate the Neumann function.
  // This result was obtained by setting the normal vector to -j.
  if (component == 0 || component == 2) {
    Tensor<1, dim> velocity_gradient =
        exact_solution.exact_velocity.gradient(p, component);
    return -nu * velocity_gradient[1];
  } else if (component == 1) {
    Tensor<1, dim> velocity_gradient =
        exact_solution.exact_velocity.gradient(p, component);
    return -nu * velocity_gradient[1] + exact_solution.exact_pressure.value(p);
  } else {
    return 0.0;
  }
}

void EthierSteinman::NeumannFunction::vector_value(
    const Point<dim> &p, Vector<double> &values) const {
  for (unsigned int i = 0; i < dim; i++) {
    values[i] = value(p, i);
  }
  values[dim] = 0.0;
}

double EthierSteinman::compute_error(const VectorTools::NormType &norm_type,
                                     bool velocity) {
  FE_SimplexP<dim> fe_mapping(1);
  MappingFE mapping(fe_mapping);

  // Set the time for the exact solution.
  exact_solution.set_time(time_step * deltat);

  // First we compute the norm on each element, and store it in a vector.
  Vector<double> error_per_cell(mesh.n_active_cells());

  if (velocity) {
    // The error is an integral, and we approximate that integral using a
    // quadrature formula. To make sure we are accurate enough, we use a
    // quadrature formula with one node more than what we used in
    // assembly.
    QGaussSimplex<dim> quadrature_error(degree_velocity + 2);

    // Source: https://www.dealii.org/current/doxygen/deal.II/step_20.html.
    ComponentSelectFunction<dim> mask(std::make_pair(0, dim), dim + 1);

    VectorTools::integrate_difference(
        mapping, dof_handler, solution, exact_solution.exact_velocity,
        error_per_cell, quadrature_error, norm_type, &mask);
  } else {
    // Do the same for the pressure.
    QGaussSimplex<dim> quadrature_error(degree_pressure + 2);
    ComponentSelectFunction<dim> mask(dim, dim + 1);
    VectorTools::integrate_difference(mapping, dof_handler, solution,
                                      exact_solution, error_per_cell,
                                      quadrature_error, norm_type, &mask);
  }

  // Then, we add out all the cells.
  const double error =
      VectorTools::compute_global_error(mesh, error_per_cell, norm_type);

  return error;
}

void EthierSteinman::apply_initial_conditions() {
  NavierStokes<dim>::apply_initial_conditions();
  pcout << "  L2 error on the velocity: "
        << compute_error(VectorTools::L2_norm, true) << std::endl;
  pcout << "  H1 error on the velocity: "
        << compute_error(VectorTools::H1_norm, true) << std::endl;
}

void EthierSteinman::solve_time_step() {
  NavierStokes<dim>::solve_time_step();
  pcout << "  L2 error on the velocity: "
        << compute_error(VectorTools::L2_norm, true) << std::endl;
  pcout << "  H1 error on the velocity: "
        << compute_error(VectorTools::H1_norm, true) << std::endl;
  pcout << "  L2 error on the pressure: "
        << compute_error(VectorTools::L2_norm, false) << std::endl;
}