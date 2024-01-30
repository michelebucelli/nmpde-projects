#include "Cylinder.hpp"

template class Cylinder<2U>;
template class Cylinder<3U>;

template <unsigned int dim>
void Cylinder<dim>::update_lift_drag() {
  NavierStokes<dim>::pcout << "  Calculating lift and drag forces" << std::endl;

  const unsigned int dofs_per_cell = NavierStokes<dim>::fe->dofs_per_cell;
  const unsigned int n_q_face = NavierStokes<dim>::quadrature_face->size();

  FEFaceValues<dim> fe_face_values(
      *NavierStokes<dim>::fe, *NavierStokes<dim>::quadrature_face,
      update_values | update_quadrature_points | update_JxW_values |
          update_gradients | update_normal_vectors);

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

        // Calculate the values and gradients of the solution at quadrature
        // points. Source:
        // https://www.dealii.org/current/doxygen/deal.II/group__vector__valued.html
        fe_face_values[velocity].get_function_values(
            NavierStokes<dim>::solution, velocity_values);

        fe_face_values[velocity].get_function_gradients(
            NavierStokes<dim>::solution, velocity_gradients);

        fe_face_values[pressure].get_function_values(
            NavierStokes<dim>::solution, pressure_values);

        for (unsigned int q = 0; q < n_q_face; ++q) {
          // Using as sources
          // https://www.sciencedirect.com/science/article/pii/S0997754619303838
          // and https://www.mate.polimi.it/biblioteca/add/qmox/mox84.pdf.

          // Get the normal vector.
          // In the formula in the paper, the normal vector is
          // defined in the opposite direction to the one in the mesh.
          Tensor<1, dim> normal = -fe_face_values.normal_vector(q);

          // Calculate the stress tensor.
          Tensor<2, dim> stress_tensor;
          stress_tensor = velocity_gradients[q];
          for (unsigned int i = 0; i < dim; i++) {
            for (unsigned int j = 0; j < dim; j++) {
              stress_tensor[i][j] += velocity_gradients[q][j][i];
            }
          }
          stress_tensor *= NavierStokes<dim>::nu;
          for (unsigned int i = 0; i < dim; i++) {
            stress_tensor[i][i] -= pressure_values[q];
          }

          // Calculate the force acting on the cylinder.
          Tensor<1, dim> force;
          force = NavierStokes<dim>::ro * stress_tensor * normal *
                  fe_face_values.JxW(q);

          // Update drag and lift forces.
          local_drag_force += force[0];
          local_lift_force += force[1];
        }
      }
    }
  }

  // Sum the drag and lift forces across all processes.
  lift_force = Utilities::MPI::sum(local_lift_force, MPI_COMM_WORLD);
  drag_force = Utilities::MPI::sum(local_drag_force, MPI_COMM_WORLD);

  // Print the results.
  NavierStokes<dim>::pcout << "  Lift coefficient: " << get_lift() << std::endl;
  NavierStokes<dim>::pcout << "  Drag coefficient: " << get_drag() << std::endl;
}

template <unsigned int dim>
double Cylinder<dim>::get_reynolds_number() const {
  return get_mean_velocity() * D / NavierStokes<dim>::nu;
}

double Cylinder2D::get_drag() const {
  const double mean_velocity = get_mean_velocity();
  return 2.0 * drag_force / (ro * mean_velocity * mean_velocity * D);
}

double Cylinder2D::get_lift() const {
  const double mean_velocity = get_mean_velocity();
  return 2.0 * lift_force / (ro * mean_velocity * mean_velocity * D * H);
}

double Cylinder3D::get_drag() const {
  const double mean_velocity = get_mean_velocity();
  return 2.0 * drag_force / (ro * mean_velocity * mean_velocity * D);
}

double Cylinder3D::get_lift() const {
  const double mean_velocity = get_mean_velocity();
  return 2.0 * lift_force / (ro * mean_velocity * mean_velocity * D * H);
}