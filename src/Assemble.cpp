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

  // Create the matrices as full matrices.
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> velocity_mass_cell_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> pressure_mass_cell_matrix(dofs_per_cell, dofs_per_cell);

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  // Reset the matrices.
  constant_matrix = 0.0;
  velocity_mass = 0.0;
  pressure_mass = 0.0;

  // Create extractors for velocity (u) and pressure (p).
  FEValuesExtractors::Vector velocity(0);
  FEValuesExtractors::Scalar pressure(dim);

  for (const auto &cell : dof_handler.active_cell_iterators()) {
    if (!cell->is_locally_owned()) continue;

    fe_values.reinit(cell);

    // Reset the cell matrices.
    cell_matrix = 0.0;
    velocity_mass_cell_matrix = 0.0;
    pressure_mass_cell_matrix = 0.0;

    for (unsigned int q = 0; q < n_q; ++q) {
      for (unsigned int i = 0; i < dofs_per_cell; ++i) {
        for (unsigned int j = 0; j < dofs_per_cell; ++j) {
          // Viscosity term (A).
          cell_matrix(i, j) +=
              nu *
              scalar_product(fe_values[velocity].gradient(i, q),
                             fe_values[velocity].gradient(j, q)) *
              fe_values.JxW(q);

          // Velocity mass term (M/deltat).
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

          // Pressure mass term (M_p/nu) (for preconditioning).
          pressure_mass_cell_matrix(i, j) += fe_values[pressure].value(j, q) *
                                             fe_values[pressure].value(i, q) /
                                             nu * fe_values.JxW(q);
        }
      }
    }

    cell->get_dof_indices(dof_indices);

    // Add the cell terms to the matrices.
    constant_matrix.add(dof_indices, cell_matrix);
    velocity_mass.add(dof_indices, velocity_mass_cell_matrix);
    pressure_mass.add(dof_indices, pressure_mass_cell_matrix);
  }

  // Each process might have written to some rows it does not own (for instance,
  // if it owns elements that are adjacent to elements owned by some other
  // process). Therefore, at the end of the assembly, processes need to exchange
  // information: the compress method allows to do this.
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
      update_values | update_quadrature_points | update_JxW_values);

  FEValuesExtractors::Vector velocity(0);
  FEValuesExtractors::Scalar pressure(dim);

  // Assemble the nonlinear term.
  {
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

    std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

    // Copy the constant term.
    system_matrix.copy_from(constant_matrix);

    // Declare a vector which will contain the values of the old solution at
    // quadrature points.
    std::vector<Tensor<1, dim>> old_solution_values(n_q);

    FEValuesExtractors::Vector velocity(0);
    FEValuesExtractors::Scalar pressure(dim);

    for (const auto &cell : dof_handler.active_cell_iterators()) {
      if (!cell->is_locally_owned()) continue;

      fe_values.reinit(cell);

      // Calculate the value of the previous solution at quadrature points.
      // Source:
      // https://www.dealii.org/current/doxygen/deal.II/group__vector__valued.html
      fe_values[velocity].get_function_values(solution_owned,
                                              old_solution_values);

      cell_matrix = 0.0;

      for (unsigned int q = 0; q < n_q; ++q) {
        // Declare tensors to store the old solution value at this quadrature
        // point and part of the nonlinear term.
        Tensor<1, dim> local_old_solution_value = old_solution_values[q];
        Tensor<1, dim> nonlinear_term;

        for (unsigned int i = 0; i < dofs_per_cell; ++i) {
          for (unsigned int j = 0; j < dofs_per_cell; ++j) {
            // Calculate (u . nabla) u.
            for (unsigned int k = 0; k < dim; k++) {
              nonlinear_term[k] = 0.0;
              for (unsigned int l = 0; l < dim; l++) {
                // Gradient component order:
                // first: component of the vector valued function
                // second: derivative direction
                // So gradient[0][1] = du_x / dy
                nonlinear_term[k] += local_old_solution_value[l] *
                                     fe_values[velocity].gradient(j, q)[k][l];
              }
            }

            // Add the term (u . nabla) uv to the matrix.
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

            // Set the correct time for the Neumann function.
            Function<dim> *neumann_function =
                neumann_boundary_functions[cell->face(f)->boundary_id()];
            neumann_function->set_time(time_step * deltat);

            // Declare vectors containing the local value of the Neumann
            // function.
            Vector<double> neumann_loc(dim + 1);
            Tensor<1, dim> neumann_loc_tensor;

            for (unsigned int q = 0; q < n_q_face; ++q) {
              // Get the local value of the Neumann function.
              neumann_function->vector_value(fe_face_values.quadrature_point(q),
                                             neumann_loc);
              for (unsigned int d = 0; d < dim; ++d)
                neumann_loc_tensor[d] = neumann_loc[d];

              for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                // Add the Neumann term.
                cell_rhs(i) +=
                    scalar_product(neumann_loc_tensor,
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

    system_rhs.compress(VectorOperation::add);

    // Add the term that comes from the old solution.
    velocity_mass.block(0, 0).vmult_add(system_rhs.block(0),
                                        solution_owned.block(0));
  }

  pcout << "  Applying Dirichlet (BC)" << std::endl;

  // Dirichlet boundary conditions.
  {
    std::map<types::global_dof_index, double> boundary_values;

    // Since only velocity has a Dirichlet (BC), create a mask.
    ComponentMask mask;
    if constexpr (dim == 2) {
      mask = ComponentMask({true, true, false});
    } else {
      mask = ComponentMask({true, true, true, false});
    }

    std::map<types::boundary_id, const Function<dim> *>
        const_dirichlet_boundary_functions;

    // Set the correct time for the Dirichlet functions and create a constant
    // copy of the map (needed because interpolate_boundary_values requires a
    // map with constant values).
    for (auto iter = dirichlet_boundary_functions.begin();
         iter != dirichlet_boundary_functions.end(); ++iter) {
      iter->second->set_time(time_step * deltat);
      const_dirichlet_boundary_functions[iter->first] =
          const_cast<const Function<dim> *>(iter->second);
    }

    VectorTools::interpolate_boundary_values(
        dof_handler, const_dirichlet_boundary_functions, boundary_values, mask);

    MatrixTools::apply_boundary_values(boundary_values, system_matrix, solution,
                                       system_rhs, false);
  }
}