#include "Cylinder.hpp"

template class Cylinder<2U>;
template class Cylinder<3U>;

template <unsigned int dim>
Cylinder<dim>::Cylinder(const std::string &mesh_file_name_,
                        const unsigned int &degree_velocity_,
                        const unsigned int &degree_pressure_, const double &T_,
                        const double &deltat_,
                        const PreconditionerType &preconditioner_type_)
    : NavierStokes<dim>(mesh_file_name_, degree_velocity_, degree_pressure_, T_,
                        deltat_, preconditioner_type_),
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
                       const double &deltat_,
                       const PreconditionerType &preconditioner_type_)
    : Cylinder<dim>(mesh_file_name_, degree_velocity_, degree_pressure_, T_,
                    deltat_, preconditioner_type_),
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
                       const double &deltat_,
                       const PreconditionerType &preconditioner_type_)
    : Cylinder<dim>(mesh_file_name_, degree_velocity_, degree_pressure_, T_,
                    deltat_, preconditioner_type_),
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