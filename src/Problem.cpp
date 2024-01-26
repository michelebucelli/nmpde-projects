#include "Problem.hpp"

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/mapping_fe.h>

template class Cylinder<2U>;
template class Cylinder<3U>;

template <unsigned int dim>
Cylinder<dim>::Cylinder(const std::string &mesh_file_name_,
                        const unsigned int &degree_velocity_,
                        const unsigned int &degree_pressure_, const double &T_,
                        const double &deltat_)
    : NavierStokes<dim>(mesh_file_name_, degree_velocity_, degree_pressure_, T_,
                        deltat_),
      zero_function(dim + 1) {
  NavierStokes<dim>::ro = 1.0;
  NavierStokes<dim>::nu = 1e-3,
  NavierStokes<dim>::initial_conditions =
      std::make_shared<Functions::ZeroFunction<dim>>(dim + 1);
}

template <unsigned int dim>
void Cylinder<dim>::update_lift_drag() {
  NavierStokes<dim>::pcout << "  Calculating lift and drag forces" << std::endl;

  const unsigned int dofs_per_cell = NavierStokes<dim>::fe->dofs_per_cell;
  const unsigned int n_q_face = NavierStokes<dim>::quadrature_face->size();

  FEFaceValues<dim> fe_face_values(*NavierStokes<dim>::fe,
                                   *NavierStokes<dim>::quadrature_face,
                                   update_values | update_quadrature_points |
                                       update_JxW_values | update_gradients);

  FEValuesExtractors::Vector velocity(0);
  FEValuesExtractors::Scalar pressure(dim);

  // Declare variables to store lift and drag forces.
  double local_lift_force = 0.0;
  double local_drag_force = 0.0;

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  // Declare a vector which will contain the values of the old solution at
  // quadrature points.
  std::vector<Tensor<1, dim>> velocity_values(n_q_face);
  std::vector<Tensor<2, dim>> velocity_gradients(n_q_face);
  std::vector<double> pressure_values(n_q_face);

  for (const auto &cell :
       NavierStokes<dim>::dof_handler.active_cell_iterators()) {
    if (!cell->is_locally_owned()) continue;

    for (unsigned int f = 0; f < cell->n_faces(); ++f) {
      if (cell->face(f)->at_boundary() &&
          cell->face(f)->boundary_id() == obstacle_tag) {
        fe_face_values.reinit(cell, f);

        // Calculate the value of the previous solution at quadrature points.
        // Source:
        // https://www.dealii.org/current/doxygen/deal.II/group__vector__valued.html
        fe_face_values[velocity].get_function_values(
            NavierStokes<dim>::solution_owned, velocity_values);

        fe_face_values[velocity].get_function_gradients(
            NavierStokes<dim>::solution_owned, velocity_gradients);

        fe_face_values[pressure].get_function_values(
            NavierStokes<dim>::solution_owned, pressure_values);

        for (unsigned int q = 0; q < n_q_face; ++q) {
          local_lift_force += -((NavierStokes<dim>::ro * NavierStokes<dim>::nu *
                                     velocity_gradients[q][1][0] +
                                 pressure_values[q] * velocity_values[q][1]) *
                                fe_face_values.JxW(q));

          local_drag_force += (NavierStokes<dim>::ro * NavierStokes<dim>::nu *
                                   velocity_gradients[q][1][0] -
                               pressure_values[q] * velocity_values[q][0]) *
                              fe_face_values.JxW(q);
        }
      }
    }
  }

  // Sum the drag and lift forces across all processes.
  lift_force = Utilities::MPI::sum(local_lift_force, MPI_COMM_WORLD);
  drag_force = Utilities::MPI::sum(local_drag_force, MPI_COMM_WORLD);
}

template <unsigned int dim>
double Cylinder<dim>::get_reynolds_number() const {
  return get_mean_velocity() * D / NavierStokes<dim>::nu;
}

template <unsigned int dim>
double Cylinder<dim>::get_drag() const {
  return 2.0 * drag_force /
         (NavierStokes<dim>::ro * get_mean_velocity() * get_mean_velocity() *
          D);
}

template <unsigned int dim>
double Cylinder<dim>::get_lift() const {
  return 2.0 * lift_force /
         (NavierStokes<dim>::ro * get_mean_velocity() * get_mean_velocity() *
          D);
}

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

double Cylinder2D::get_mean_velocity() const {
  return 2.0 * inlet_velocity.value(Point<dim>(0, H / 2.0), 0) / 3.0;
}

double Cylinder3D::get_mean_velocity() const {
  return 4.0 * inlet_velocity.value(Point<dim>(0, H / 2.0, H / 2.0), 0) / 9.0;
}

Cylinder2D::Cylinder2D(const std::string &mesh_file_name_,
                       const unsigned int &degree_velocity_,
                       const unsigned int &degree_pressure_, const double &T_,
                       const double &deltat_)
    : Cylinder<dim>(mesh_file_name_, degree_velocity_, degree_pressure_, T_,
                    deltat_),
      inlet_velocity(1.5, H) {
  U_m = 1.5;
  obstacle_tag = 5U;

  dirichlet_boundary_functions[4] = &inlet_velocity;

  neumann_boundary_functions[2] = &zero_function;

  dirichlet_boundary_functions[1] = &zero_function;
  dirichlet_boundary_functions[3] = &zero_function;
  dirichlet_boundary_functions[5] = &zero_function;
}

Cylinder3D::Cylinder3D(const std::string &mesh_file_name_,
                       const unsigned int &degree_velocity_,
                       const unsigned int &degree_pressure_, const double &T_,
                       const double &deltat_)
    : Cylinder<dim>(mesh_file_name_, degree_velocity_, degree_pressure_, T_,
                    deltat_),
      inlet_velocity(2.25, H) {
  U_m = 2.25;
  obstacle_tag = 7U;

  dirichlet_boundary_functions[5] = &inlet_velocity;

  neumann_boundary_functions[3] = &zero_function;

  dirichlet_boundary_functions[1] = &zero_function;
  dirichlet_boundary_functions[2] = &zero_function;
  dirichlet_boundary_functions[4] = &zero_function;
  dirichlet_boundary_functions[6] = &zero_function;
  dirichlet_boundary_functions[7] = &zero_function;
}

EthierSteinman::EthierSteinman(const std::string &mesh_file_name_,
                               const unsigned int &degree_velocity_,
                               const unsigned int &degree_pressure_,
                               const double &T_, const double &deltat_,
                               const double &nu_)
    : NavierStokes(mesh_file_name_, degree_velocity_, degree_pressure_, T_,
                   deltat_),
      exact_solution(nu_),
      neumann_function(nu_) {
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
  for (unsigned int i = 0; i < dim; i++) {
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
  for (unsigned int i = 0; i < dim; i++) {
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
  pcout << "L2 error on the velocity: "
        << compute_error(VectorTools::L2_norm, true) << std::endl;
  pcout << "H1 error on the velocity: "
        << compute_error(VectorTools::H1_norm, true) << std::endl;
}

void EthierSteinman::solve_time_step() {
  NavierStokes<dim>::solve_time_step();
  pcout << "L2 error on the velocity: "
        << compute_error(VectorTools::L2_norm, true) << std::endl;
  pcout << "H1 error on the velocity: "
        << compute_error(VectorTools::H1_norm, true) << std::endl;
  pcout << "L2 error on the pressure: "
        << compute_error(VectorTools::L2_norm, false) << std::endl;
}

Step::Step(const std::string &mesh_file_name_,
           const unsigned int &degree_velocity_,
           const unsigned int &degree_pressure_, const double &T_,
           const double &deltat_, const double &alpha_)
    : NavierStokes(mesh_file_name_, degree_velocity_, degree_pressure_, T_,
                   deltat_),
      inlet_velocity(alpha_),
      zero_function(dim + 1) {
  ro = 1.0;
  nu = 1.0;
  initial_conditions = std::make_shared<Functions::ZeroFunction<dim>>(dim + 1);

  dirichlet_boundary_functions[0] = &inlet_velocity;

  neumann_boundary_functions[2] = &neumann_function;

  dirichlet_boundary_functions[1] = &zero_function;
}

double Step::InletVelocity::value(const Point<dim> &p,
                                  const unsigned int component) const {
  if (component == 0)
    return -alpha * p[1] * (2.0 - p[1]) * (1.0 - p[2]) * (2.0 - p[2]);
  else
    return 0.0;
}

void Step::InletVelocity::vector_value(const Point<dim> &p,
                                       Vector<double> &values) const {
  values[0] = value(p, 0);

  for (unsigned int i = 1; i < dim + 1; ++i) values[i] = 0.0;
}

double Step::NeumannFunction::value(const Point<dim> & /*p*/,
                                    const unsigned int component) const {
  if (component == 0)
    return -p_out;
  else
    return 0.0;
}

void Step::NeumannFunction::vector_value(const Point<dim> & /*p*/,
                                         Vector<double> &values) const {
  values[0] = -p_out;

  for (unsigned int i = 1; i < dim + 1; ++i) values[i] = 0.0;
}