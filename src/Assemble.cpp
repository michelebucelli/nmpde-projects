#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include "NavierStokes.hpp"

template class NavierStokes<2U>;
template class NavierStokes<3U>;

template <unsigned int dim>
void NavierStokes<dim>::assemble_constant() {
  pcout << "===============================================" << std::endl;
  pcout << "Assembling the constant terms of the system" << std::endl;

  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q = quadrature->size();

  FEValues<dim> fe_values(*fe, *quadrature,
                          update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> velocity_mass_cell_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> pressure_mass_cell_matrix(dofs_per_cell, dofs_per_cell);

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  constant_matrix = 0.0;
  velocity_mass = 0.0;
  pressure_mass = 0.0;

  FEValuesExtractors::Vector velocity(0);
  FEValuesExtractors::Scalar pressure(dim);

  for (const auto &cell : dof_handler.active_cell_iterators()) {
    if (!cell->is_locally_owned()) continue;

    fe_values.reinit(cell);

    cell_matrix = 0.0;
    velocity_mass_cell_matrix = 0.0;

    for (unsigned int q = 0; q < n_q; ++q) {
      for (unsigned int i = 0; i < dofs_per_cell; ++i) {
        for (unsigned int j = 0; j < dofs_per_cell; ++j) {
          // Viscosity term (A).
          cell_matrix(i, j) +=
              nu *
              scalar_product(fe_values[velocity].gradient(i, q),
                             fe_values[velocity].gradient(j, q)) *
              fe_values.JxW(q);

          // Mass term (M/deltat).
          const double mass_value =
              scalar_product(fe_values[velocity].value(j, q),
                             fe_values[velocity].value(i, q)) *
              fe_values.JxW(q) / deltat;
          cell_matrix(i, j) += mass_value;
          velocity_mass_cell_matrix(i, j) += mass_value;

          // Pressure term in the momentum equation (B^T).
          cell_matrix(i, j) -= fe_values[velocity].divergence(i, q) *
                               fe_values[pressure].value(j, q) *
                               fe_values.JxW(q);

          // Pressure term in the continuity equation (-B).
          cell_matrix(i, j) += fe_values[velocity].divergence(j, q) *
                               fe_values[pressure].value(i, q) *
                               fe_values.JxW(q);

          // Pressure mass term.
          pressure_mass_cell_matrix(i, j) += fe_values[pressure].value(j, q) *
                                             fe_values[pressure].value(i, q) *
                                             fe_values.JxW(q);
        }
      }
    }

    cell->get_dof_indices(dof_indices);

    constant_matrix.add(dof_indices, cell_matrix);
    velocity_mass.add(dof_indices, velocity_mass_cell_matrix);
    pressure_mass.add(dof_indices, pressure_mass_cell_matrix);
  }

  constant_matrix.compress(VectorOperation::add);
  velocity_mass.compress(VectorOperation::add);
  pressure_mass.compress(VectorOperation::add);
}

template <unsigned int dim>
void NavierStokes<dim>::assemble_time_dependent() {
  pcout << "  Assembling the nonlinear term" << std::endl;

  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q = quadrature->size();
  const unsigned int n_q_face = quadrature_face->size();

  FEValues<dim> fe_values(*fe, *quadrature,
                          update_values | update_quadrature_points |
                              update_JxW_values | update_gradients);
  FEFaceValues<dim> fe_face_values(
      *fe, *quadrature_face,
      update_values | update_normal_vectors | update_JxW_values);

  FEValuesExtractors::Vector velocity(0);
  FEValuesExtractors::Scalar pressure(dim);

  // Assemble the nonlinear term.
  {
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

    std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
    std::vector<Tensor<1, dim>> old_solution_values(n_q);
    // std::vector<Tensor<2, dim>> old_solution_gradients(n_q);

    // Copy the contant term.
    system_matrix.copy_from(constant_matrix);

    FEValuesExtractors::Vector velocity(0);
    FEValuesExtractors::Scalar pressure(dim);

    for (const auto &cell : dof_handler.active_cell_iterators()) {
      if (!cell->is_locally_owned()) continue;

      fe_values.reinit(cell);

      // Calculate the value and gradient of the previous solution.
      // Source:
      // https://www.dealii.org/current/doxygen/deal.II/group__vector__valued.html
      // https://www.dealii.org/current/doxygen/deal.II/classFEValuesViews_1_1Vector.html#ace19727c285e8035282a6f3f66ce7f18
      fe_values[velocity].get_function_values(solution, old_solution_values);
      // fe_values[velocity].get_function_gradients(solution,
      // old_solution_gradients);

      cell_matrix = 0.0;

      for (unsigned int q = 0; q < n_q; ++q) {
        for (unsigned int i = 0; i < dofs_per_cell; ++i) {
          for (unsigned int j = 0; j < dofs_per_cell; ++j) {
            // Calculate (u . nabla) u.
            Tensor<1, dim> local_old_solution_value = old_solution_values[q];
            // Tensor<2, dim> local_old_solution_gradient =
            // old_solution_gradients[q];
            Tensor<1, dim> nonlinear_term;

            for (unsigned int k = 0; k < dim; k++) {
              nonlinear_term[k] = 0.0;
              for (unsigned int l = 0; l < dim; l++) {
                // local_old_solution_gradient component order:
                // first: component of the vector valued function
                // second: derivative direction
                // local_old_solution_gradient[0][1] = du^n_x / dy
                nonlinear_term[k] += local_old_solution_value[l] *
                                     fe_values[velocity].gradient(j, q)[k][l];
              }
            }

            // Add the term (u . nabla) uv.
            cell_matrix(i, j) +=
                scalar_product(nonlinear_term,
                               fe_values[velocity].value(i, q)) *
                fe_values.JxW(q);
          }
        }
      }

      cell->get_dof_indices(dof_indices);

      system_matrix.add(dof_indices, cell_matrix);
    }

    system_matrix.compress(VectorOperation::add);
  }

  pcout << "  Assembling the right-hand side" << std::endl;

  // Assemble the right-hand side.
  {
    Vector<double> cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

    system_rhs = 0.0;

    for (const auto &cell : dof_handler.active_cell_iterators()) {
      if (!cell->is_locally_owned()) continue;

      // Boundary integral for Neumann BCs.
      if (cell->at_boundary()) {
        for (unsigned int f = 0; f < cell->n_faces(); ++f) {
          if (cell->face(f)->at_boundary() &&
              neumann_boundary_functions.count(cell->face(f)->boundary_id())) {
            fe_face_values.reinit(cell, f);

            for (unsigned int q = 0; q < n_q_face; ++q) {
              for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                cell_rhs(i) -=
                    neumann_boundary_functions[cell->face(f)->boundary_id()] *
                    scalar_product(fe_face_values.normal_vector(q),
                                   fe_face_values[velocity].value(i, q)) *
                    fe_face_values.JxW(q);
              }
            }
          }
        }
      }

      cell->get_dof_indices(dof_indices);
      system_rhs.add(dof_indices, cell_rhs);
    }

    // Add the term that comes from the old solution.
    velocity_mass.vmult_add(system_rhs, solution_owned);

    system_rhs.compress(VectorOperation::add);
  }

  pcout << "  Applying Dirichlet (BC)" << std::endl;

  // Dirichlet boundary conditions.
  {
    std::map<types::global_dof_index, double> boundary_values;

    ComponentMask mask;
    if constexpr (dim == 2) {
      mask = ComponentMask({true, true, false});
    } else {
      mask = ComponentMask({true, true, true, false});
    }

    VectorTools::interpolate_boundary_values(
        dof_handler, dirichlet_boundary_functions, boundary_values, mask);

    MatrixTools::apply_boundary_values(boundary_values, system_matrix, solution,
                                       system_rhs, false);
  }
}