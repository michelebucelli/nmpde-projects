#include <deal.II/base/config.h>
#include <deal.II/fe/component_mask.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

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

template <unsigned int dim>
void Cylinder<dim>::update_lift_drag_weak() {
  // Find a phi_inf in V_h that satisfies the boundary conditions
  // - on the obstacle: phi_inf = -u_in
  // - on the other boundaries: phi_inf = 0
  // Since the constant part of the matrix is square, invertible and with the
  // right dimensions, phi_inf will be calculated solving a system with that
  // matrix and using a vector of 1s as right-hand side, then applying the
  // Dirichlet conditions.

  // Create matrices and vectors.
  TrilinosWrappers::SparseMatrix matrix;
  matrix.reinit(NavierStokes<dim>::constant_matrix.block(0, 0));
  matrix.copy_from(NavierStokes<dim>::constant_matrix.block(0, 0));
  TrilinosWrappers::MPI::Vector rhs;
  rhs.reinit(NavierStokes<dim>::system_rhs.block(0));
  rhs.add(1.0);
  TrilinosWrappers::MPI::Vector phi_inf_owned;
  TrilinosWrappers::MPI::Vector phi_inf;
  phi_inf_owned.reinit(NavierStokes<dim>::solution_owned.block(0));
  phi_inf.reinit(NavierStokes<dim>::solution.block(0));

  // Create a map with the correct Dirichlet boundary functions.
  std::map<types::boundary_id, const Function<dim> *>
      const_dirichlet_boundary_functions;
  for (auto iter = NavierStokes<dim>::dirichlet_boundary_functions.begin();
       iter != NavierStokes<dim>::dirichlet_boundary_functions.end(); ++iter) {
    if (iter->first != obstacle_tag) {
      Function<dim> *zero_function_ptr = &zero_function;
      const_dirichlet_boundary_functions[iter->first] =
          const_cast<const Function<dim> *>(zero_function_ptr);
    } else {
      const_dirichlet_boundary_functions[iter->first] =
          const_cast<const Function<dim> *>(&(*inlet_velocity));
    }
  }
  for (auto iter = NavierStokes<dim>::neumann_boundary_functions.begin();
       iter != NavierStokes<dim>::neumann_boundary_functions.end(); ++iter) {
    Function<dim> *zero_function_ptr = &zero_function;
    const_dirichlet_boundary_functions[iter->first] =
        const_cast<const Function<dim> *>(zero_function_ptr);
  }

  // Apply the boundary conditions.
  ComponentMask mask;
  if constexpr (dim == 2) {
    mask = ComponentMask({true, true, false});
  } else {
    mask = ComponentMask({true, true, true, false});
  }
  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(NavierStokes<dim>::dof_handler,
                                           const_dirichlet_boundary_functions,
                                           boundary_values, mask);
  MatrixTools::apply_boundary_values(boundary_values, matrix, phi_inf,
                                     rhs, false);

  // Solve the system.
  SolverControl solver_control(
      NavierStokes<dim>::solver_options.maxiter,
      NavierStokes<dim>::solver_options.tol * rhs.l2_norm());
  TrilinosWrappers::PreconditionAMG precondition;
  precondition.initialize(matrix);
  SolverGMRES<TrilinosWrappers::MPI::Vector> solver(solver_control);
  solver.solve(matrix, phi_inf_owned, rhs, precondition);
  phi_inf = phi_inf_owned;

  DataOut<dim> data_out;

  // Output (for debugging).
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(
          dim, DataComponentInterpretation::component_is_part_of_vector);
  std::vector<std::string> names;
  if constexpr (dim == 2) {
    names = {"velocity", "velocity"};
  } else {
    names = {"velocity", "velocity", "velocity"};
  }
  data_out.add_data_vector(NavierStokes<dim>::dof_handler, phi_inf, names,
                           data_component_interpretation);
  std::vector<unsigned int> partition_int(
      NavierStokes<dim>::mesh.n_active_cells());
  GridTools::get_subdomain_association(NavierStokes<dim>::mesh, partition_int);
  const Vector<double> partitioning(partition_int.begin(), partition_int.end());
  data_out.add_data_vector(partitioning, "partitioning");
  data_out.build_patches();
  const std::string output_file_name = "output-phi-inf";
  data_out.write_vtu_with_pvtu_record("../results/", output_file_name,
                                      NavierStokes<dim>::time_step,
                                      MPI_COMM_WORLD, 3);
}