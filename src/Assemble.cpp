#include "NavierStokes.hpp"

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/component_mask.h>

#include <deal.II/lac/full_matrix.h>

#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

void
NavierStokes::assemble_constant()
{
  pcout << "===============================================" << std::endl;
  pcout << "Assembling the system" << std::endl;

  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q           = quadrature->size();

  FEValues<dim>     fe_values(*fe,
                          *quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> mass_cell_matrix(dofs_per_cell, dofs_per_cell);

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  constant_matrix = 0.0;
  velocity_mass   = 0.0;

  FEValuesExtractors::Vector velocity(0);
  FEValuesExtractors::Scalar pressure(dim);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);

      cell_matrix               = 0.0;
      mass_cell_matrix          = 0.0;

      for (unsigned int q = 0; q < n_q; ++q)
        {
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                  // Viscosity term (A).
                  cell_matrix(i, j) += nu *
                                       scalar_product(fe_values[velocity].gradient(i, q),
                                       fe_values[velocity].gradient(j, q)) *
                                       fe_values.JxW(q);

                  // Mass term (M/deltat).
                  const double mass_value = scalar_product(
                                            fe_values[velocity].value(j,q),
                                            fe_values[velocity].value(i,q)) *
                                            fe_values.JxW(q) /
                                            deltat;
                  cell_matrix(i, j)      += mass_value;
                  mass_cell_matrix(i, j) += mass_value;

                  // Pressure term in the momentum equation (B^T).
                  cell_matrix(i, j) -= fe_values[velocity].divergence(i, q) *
                                       fe_values[pressure].value(j, q) *
                                       fe_values.JxW(q);

                  // Pressure term in the continuity equation (-B).
                  cell_matrix(i, j) += fe_values[velocity].divergence(j, q) *
                                       fe_values[pressure].value(i, q) *
                                       fe_values.JxW(q);
                }              
            }
        }

      cell->get_dof_indices(dof_indices);

      constant_matrix.add(dof_indices, cell_matrix);
      velocity_mass.add(dof_indices, mass_cell_matrix);
    }

  constant_matrix.compress(VectorOperation::add);
  velocity_mass.compress(VectorOperation::add);
}

void
NavierStokes::assemble_time_dependent()
{
  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q           = quadrature->size();
  const unsigned int n_q_face      = quadrature_face->size();

  FEValues<dim> fe_values(*fe,
                          *quadrature,
                          update_values | update_quadrature_points |
                            update_JxW_values);
  FEFaceValues<dim> fe_face_values(*fe,
                                   *quadrature_face,
                                   update_values | update_normal_vectors |
                                     update_JxW_values);

  FEValuesExtractors::Vector velocity(0);
  FEValuesExtractors::Scalar pressure(dim);

  // Assemble the nonlinear term.
  {
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

    std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
    std::vector<Tensor<2, dim>> old_solution_gradients(n_q);

    // Copy the contant term.
    system_matrix.copy_from(constant_matrix);

    FEValuesExtractors::Vector velocity(0);
    FEValuesExtractors::Scalar pressure(dim);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (!cell->is_locally_owned())
          continue;

        fe_values.reinit(cell);

        // Calculate the gradient of the previous solution.
        // Source: https://www.dealii.org/current/doxygen/deal.II/group__vector__valued.html
        fe_values[velocity].get_function_gradients(solution, old_solution_gradients);

        cell_matrix = 0.0;

        for (unsigned int q = 0; q < n_q; ++q)
          {
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  {
                    // Calculate (u . nabla) u.
                    Tensor<2, dim> local_old_solution_gradient = old_solution_gradients[q];
                    Tensor<1, dim> nonlinear_term;

                    for(unsigned int k=0; k<dim; k++) {
                      nonlinear_term[i] = 0.0;
                      for(unsigned int l=0; l<dim; l++) {
                        nonlinear_term[i] += fe_values[velocity].value(j, q)[l] *
                                             local_old_solution_gradient[k][l];
                      }
                    }

                    // Add the term (u . nabla) uv.
                    cell_matrix(i, j) += scalar_product(
                                         nonlinear_term,
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

  // Assemble the right-hand side.
  {
    Vector<double> cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

    system_rhs = 0.0;

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (!cell->is_locally_owned())
          continue;

        // Boundary integral for Neumann BCs.
        if (cell->at_boundary())
          {
            for (unsigned int f = 0; f < cell->n_faces(); ++f)
              {
                if (cell->face(f)->at_boundary() &&
                    cell->face(f)->boundary_id() == 2)
                  {
                    fe_face_values.reinit(cell, f);

                    for (unsigned int q = 0; q < n_q_face; ++q)
                      {
                        for (unsigned int i = 0; i < dofs_per_cell; ++i)
                          {
                            cell_rhs(i) +=
                              -p_out *
                              scalar_product(fe_face_values.normal_vector(q),
                                            fe_face_values[velocity].value(i,q)) *
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
    
    // Dirichlet boundary conditions.
    {
      std::map<types::global_dof_index, double>           boundary_values;
      std::map<types::boundary_id, const Function<dim> *> boundary_functions;

      ComponentMask mask;
      if constexpr(dim == 2) {
        mask = ComponentMask({true, true, false});
      } else {
        mask = ComponentMask({true, true, true, false});
      }

      // We interpolate first the dirichlet condition on the inlet alone, then on
      // the other boundaries alone, so that the latter "win" over the former where the two
      // boundaries touch.
      boundary_functions[0] = &inlet_velocity;
      VectorTools::interpolate_boundary_values(dof_handler,
                                              boundary_functions,
                                              boundary_values,
                                              mask);

      boundary_functions.clear();
      Functions::ZeroFunction<dim> zero_function(dim + 1);
      boundary_functions[1] = &zero_function;
      VectorTools::interpolate_boundary_values(dof_handler,
                                              boundary_functions,
                                              boundary_values,
                                              mask);

      MatrixTools::apply_boundary_values(
        boundary_values, constant_matrix, solution, system_rhs, false);
    }  
  }
}