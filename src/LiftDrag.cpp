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

        // Calculate the value of the previous solution at quadrature points.
        // Source:
        // https://www.dealii.org/current/doxygen/deal.II/group__vector__valued.html
        fe_face_values[velocity].get_function_values(
            NavierStokes<dim>::solution, velocity_values);

        fe_face_values[velocity].get_function_gradients(
            NavierStokes<dim>::solution, velocity_gradients);

        fe_face_values[pressure].get_function_values(
            NavierStokes<dim>::solution, pressure_values);

        for (unsigned int q = 0; q < n_q_face; ++q) {
          // Get normal and tangent vectors.
          // In formula in flow past a cylinder paper, the normal vector is
          // defined in the opposite direction to the one in the mesh.
          Tensor<1, dim> normal = -fe_face_values.normal_vector(q);
          Tensor<1, dim> tangent;
          tangent[0] = normal[1];
          tangent[1] = -normal[0];
          for (unsigned int i = 2; i < dim; ++i) {
            tangent[i] = 0.0;
          }

          // Get local pressure value and velocity gradient.
          // What we usually refer to as "pressure" is multiplied by ro, as it
          // actually represents the pressure divided by the density.
          double pressure = NavierStokes<dim>::ro * pressure_values[q];
          Tensor<2, dim> velocity_gradient = velocity_gradients[q];

          // Calculate the gradient of the tangential velocity.
          Tensor<1, dim> tangential_gradient;
          for (unsigned int j = 0; j < dim; ++j) {
            tangential_gradient[j] = 0.0;
            for (unsigned int i = 0; i < dim; ++i) {
              tangential_gradient[j] += tangent[i] * velocity_gradient[i][j];
            }
          }

          // Calculate the normal derivative of the tangential velocity.
          double der_velocity = scalar_product(tangential_gradient, normal);

          local_drag_force += (NavierStokes<dim>::ro * NavierStokes<dim>::nu *
                                   der_velocity * normal[1] -
                               pressure * normal[0]) *
                              fe_face_values.JxW(q);

          local_lift_force -= (NavierStokes<dim>::ro * NavierStokes<dim>::nu *
                                   der_velocity * normal[0] -
                               pressure * normal[1]) *
                              fe_face_values.JxW(q);
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