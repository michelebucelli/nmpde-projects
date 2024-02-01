#include <deal.II/base/config.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/component_mask.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include "Cylinder.hpp"

template class Cylinder<2U>;
template class Cylinder<3U>;

// A function that returns the vector -i.
template <unsigned int dim>
class NegativeIFunction : public Function<dim> {
 public:
  NegativeIFunction() : Function<dim>(dim + 1) {}
  virtual double value(const Point<dim> & /*p*/,
                       const unsigned int component = 0) const override {
    if (component == 0) {
      return -1.0;
    } else {
      return 0.0;
    }
  }
  virtual void vector_value(const Point<dim> & /*p*/,
                            Vector<double> &values) const override {
    values[0] = -1.0;
    for (unsigned int i = 1; i < dim + 1; i++) {
      values[i] = 0.0;
    }
  }
};

// A function that returns the vector -j.
template <unsigned int dim>
class NegativeJFunction : public Function<dim> {
 public:
  NegativeJFunction() : Function<dim>(dim + 1) {}
  virtual double value(const Point<dim> & /*p*/,
                       const unsigned int component = 0) const override {
    if (component == 1) {
      return -1.0;
    } else {
      return 0.0;
    }
  }
  virtual void vector_value(const Point<dim> & /*p*/,
                            Vector<double> &values) const override {
    for (unsigned int i = 0; i < dim + 1; i++) {
      values[i] = 0.0;
    }
    values[1] = -1.0;
  }
};

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

  // Allocate space for the stress tensor and force.
  Tensor<2, dim> stress_tensor;
  Tensor<1, dim> force;

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
          // In the formula in the benchmark and report, the normal vector is
          // defined in the opposite direction to the one in the mesh.
          const Tensor<1, dim> &neg_normal = fe_face_values.normal_vector(q);

          // Calculate the stress tensor.
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
          force = -NavierStokes<dim>::ro * stress_tensor * neg_normal *
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
  NavierStokes<dim>::pcout << "  Strong lift coefficient: " << get_lift(false)
                           << std::endl;
  NavierStokes<dim>::pcout << "  Strong drag coefficient: " << get_drag(false)
                           << std::endl;
}

// This function calculates two functions phi_inf_lift and phi_inf_drag that
// satisfy the following conditions:
// For the drag force:
// - on the obstacle: phi_inf = -i
// - on the other boundaries: phi_inf = 0
// For the lift force:
// - on the obstacle: phi_inf = -j
// - on the other boundaries: phi_inf = 0
// A simple way to achieve this is setting the integral over Omega of phi_inf *
// v_h to 0 for all test functions v_h in V_h and the appling lifting. This way,
// the mass matrix for the velocity, which is already assembled and not badly
// conditioned, is used.
template <unsigned int dim>
void Cylinder<dim>::calculate_phi_inf() {
  NavierStokes<dim>::pcout << "==============================================="
                           << std::endl;
  NavierStokes<dim>::pcout << "Calculating phi inf" << std::endl;

  NavierStokes<dim>::pcout << "  Creating the sparsity pattern" << std::endl;

  // Create the sparsity pattern and initialize the phi_inf vectors.
  Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1);
  for (unsigned int c = 0; c < dim + 1; ++c) {
    for (unsigned int d = 0; d < dim + 1; ++d) {
      if ((c == dim && d != dim) ||
          (d == dim && c != dim))  // pressure-velocity terms
        coupling[c][d] = DoFTools::none;
      else  // other combinations
        coupling[c][d] = DoFTools::always;
    }
  }
  TrilinosWrappers::BlockSparsityPattern sparsity(
      NavierStokes<dim>::block_owned_dofs, MPI_COMM_WORLD);
  DoFTools::make_sparsity_pattern(NavierStokes<dim>::dof_handler, coupling,
                                  sparsity);
  sparsity.compress();
  phi_inf_lift.reinit(NavierStokes<dim>::solution);
  phi_inf_drag.reinit(NavierStokes<dim>::solution);

  // Run this for each of the lift and drag coefficients.
  for (unsigned int f = 0; f < 2; f++) {
    if (f == 0) {
      NavierStokes<dim>::pcout
          << "  Calculating phi inf for the lift coefficient" << std::endl;
    } else {
      NavierStokes<dim>::pcout
          << "  Calculating phi inf for the drag coefficient" << std::endl;
    }

    NavierStokes<dim>::pcout << "  Initializing matrices and vectors"
                             << std::endl;

    // Create needed matrices and vectors.
    // Matrix.
    TrilinosWrappers::BlockSparseMatrix matrix;
    matrix.reinit(sparsity);
    matrix.block(0, 0).copy_from(NavierStokes<dim>::velocity_mass.block(0, 0));
    matrix.block(0, 0) *= NavierStokes<dim>::deltat;
    matrix.block(1, 1) = 0.0;
    // Rhs.
    TrilinosWrappers::MPI::BlockVector rhs;
    rhs.reinit(NavierStokes<dim>::system_rhs);
    rhs = 0.0;
    // Solution.
    TrilinosWrappers::MPI::BlockVector phi_inf_owned;
    TrilinosWrappers::MPI::BlockVector *phi_inf =
        (f == 0 ? &phi_inf_lift : &phi_inf_drag);
    phi_inf_owned.reinit(NavierStokes<dim>::solution_owned);
    phi_inf->reinit(NavierStokes<dim>::solution);

    // Handle boundary conditions.
    NavierStokes<dim>::pcout << "  Setting boundary conditions" << std::endl;
    // Set the Dirichlet function on the obstacle boundary.
    std::shared_ptr<const Function<dim>> obstacle_function;
    if (f == 0) {
      obstacle_function = std::make_shared<const NegativeJFunction<dim>>();
    } else {
      obstacle_function = std::make_shared<const NegativeIFunction<dim>>();
    }
    // Create a map with the correct Dirichlet boundary functions.
    std::map<types::boundary_id, const Function<dim> *>
        const_dirichlet_boundary_functions;
    for (auto iter = NavierStokes<dim>::dirichlet_boundary_functions.begin();
         iter != NavierStokes<dim>::dirichlet_boundary_functions.end();
         ++iter) {
      if (iter->first != obstacle_tag) {
        Function<dim> *zero_function_ptr = &zero_function;
        const_dirichlet_boundary_functions[iter->first] =
            const_cast<const Function<dim> *>(zero_function_ptr);
      } else {
        const_dirichlet_boundary_functions[iter->first] = &(*obstacle_function);
      }
    }
    for (auto iter = NavierStokes<dim>::neumann_boundary_functions.begin();
         iter != NavierStokes<dim>::neumann_boundary_functions.end(); ++iter) {
      Function<dim> *zero_function_ptr = &zero_function;
      const_dirichlet_boundary_functions[iter->first] =
          const_cast<const Function<dim> *>(zero_function_ptr);
    }

    // Apply the boundary conditions.
    std::map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values(NavierStokes<dim>::dof_handler,
                                             const_dirichlet_boundary_functions,
                                             boundary_values);
    MatrixTools::apply_boundary_values(boundary_values, matrix, phi_inf_owned,
                                       rhs, false);

    // Solve the system.
    NavierStokes<dim>::pcout << "  Solving the system" << std::endl;
    constexpr double tolerance = 1e-12;
    SolverControl solver_control(NavierStokes<dim>::solver_options.maxiter,
                                 tolerance);
    TrilinosWrappers::PreconditionILU precondition;
    precondition.initialize(matrix.block(0, 0));
    SolverGMRES<TrilinosWrappers::MPI::Vector> solver(solver_control);
    solver.solve(matrix.block(0, 0), phi_inf_owned.block(0), rhs.block(0),
                 precondition);
    phi_inf_owned.block(1) = 0.0;
    *phi_inf = phi_inf_owned;

    // Output (for debugging).
    NavierStokes<dim>::pcout << "  Writing the result for debugging"
                             << std::endl;
    DataOut<dim> data_out;
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
        data_component_interpretation(
            dim, DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation.push_back(
        DataComponentInterpretation::component_is_scalar);
    std::vector<std::string> names;
    if constexpr (dim == 2) {
      names = {"velocity", "velocity", "placeholder"};
    } else {
      names = {"velocity", "velocity", "velocity", "placeholder"};
    }
    data_out.add_data_vector(NavierStokes<dim>::dof_handler, *phi_inf, names,
                             data_component_interpretation);
    std::vector<unsigned int> partition_int(
        NavierStokes<dim>::mesh.n_active_cells());
    GridTools::get_subdomain_association(NavierStokes<dim>::mesh,
                                         partition_int);
    const Vector<double> partitioning(partition_int.begin(),
                                      partition_int.end());
    data_out.add_data_vector(partitioning, "partitioning");
    data_out.build_patches();
    const std::string output_file_name = "output-phi-inf";
    data_out.write_vtu_with_pvtu_record("../results/", output_file_name, f,
                                        MPI_COMM_WORLD, 3);
  }
}

// Compute the lift and drag coefficients integrating over the domain.
template <unsigned int dim>
void Cylinder<dim>::update_lift_drag_weak() {
  const unsigned int dofs_per_cell = NavierStokes<dim>::fe->dofs_per_cell;
  const unsigned int n_q = NavierStokes<dim>::quadrature->size();

  FEValues<dim> fe_values(*NavierStokes<dim>::fe,
                          *NavierStokes<dim>::quadrature,
                          update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  // Declare vectors which will contain the values and gradients of the required
  // quantities at quadrature points.
  std::vector<Tensor<1, dim>> velocity_values(n_q);
  std::vector<Tensor<2, dim>> velocity_gradients(n_q);
  std::vector<double> pressure_values(n_q);
  std::vector<Tensor<1, dim>> phi_inf_values(n_q);
  std::vector<Tensor<2, dim>> phi_inf_gradients(n_q);

  // Create extractors for velocity (u) and pressure (p).
  FEValuesExtractors::Vector velocity(0);
  FEValuesExtractors::Scalar pressure(dim);

  // Allocate space for the transposed gradient of u and part of the nonlinear
  // term.
  Tensor<2, dim> transposed_gradient;
  Tensor<1, dim> nonlinear_term;

  for (unsigned int f = 0; f < 2; f++) {
    TrilinosWrappers::MPI::BlockVector *phi_inf =
        (f == 0 ? &phi_inf_lift : &phi_inf_drag);
    // Declare a variable to store the local force.
    double local_force = 0.0;

    for (const auto &cell :
         NavierStokes<dim>::dof_handler.active_cell_iterators()) {
      if (!cell->is_locally_owned()) continue;

      fe_values.reinit(cell);

      // Get the values of the needed quantities at the quadrature points for
      // this cell.
      fe_values[velocity].get_function_values(NavierStokes<dim>::solution,
                                              velocity_values);
      fe_values[velocity].get_function_gradients(NavierStokes<dim>::solution,
                                                 velocity_gradients);
      fe_values[pressure].get_function_values(NavierStokes<dim>::solution,
                                              pressure_values);
      fe_values[velocity].get_function_values(*phi_inf, phi_inf_values);
      fe_values[velocity].get_function_gradients(*phi_inf, phi_inf_gradients);

      for (unsigned int q = 0; q < n_q; ++q) {
        // Viscosity term a(u, phi_inf).
        for (unsigned int i = 0; i < dim; i++) {
          for (unsigned int j = 0; j < dim; j++) {
            transposed_gradient[i][j] = velocity_gradients[q][j][i];
          }
        }
        local_force +=
            NavierStokes<dim>::nu *
            scalar_product(velocity_gradients[q] + transposed_gradient,
                           phi_inf_gradients[q]) *
            fe_values.JxW(q);

        // Pressure term in the momentum equation b(p, phi_inf).
        double phi_inf_divergence = 0.0;
        for (unsigned int i = 0; i < dim; i++) {
          phi_inf_divergence += phi_inf_gradients[q][i][i];
        }
        local_force -=
            phi_inf_divergence * pressure_values[q] * fe_values.JxW(q);

        // Nonlinear term.
        // Calculate (u . nabla) u.
        for (unsigned int k = 0; k < dim; k++) {
          nonlinear_term[k] = 0.0;
          for (unsigned int l = 0; l < dim; l++) {
            nonlinear_term[k] +=
                velocity_values[q][l] * velocity_gradients[q][k][l];
          }
        }
        // Add the term (u . nabla) phi_inf.
        local_force += scalar_product(nonlinear_term, phi_inf_values[q]) *
                       fe_values.JxW(q);
      }
    }

    // Sum the lift and drag forces across all processes.
    if (f == 0) {
      lift_force_weak = Utilities::MPI::sum(local_force, MPI_COMM_WORLD);
    } else {
      drag_force_weak = Utilities::MPI::sum(local_force, MPI_COMM_WORLD);
    }
  }

  // Print the results.
  NavierStokes<dim>::pcout << "  Weak lift coefficient: " << get_lift(true)
                           << std::endl;
  NavierStokes<dim>::pcout << "  Weak drag coefficient: " << get_drag(true)
                           << std::endl;
}

template <unsigned int dim>
double Cylinder<dim>::get_reynolds_number() const {
  return get_mean_velocity() * D / NavierStokes<dim>::nu;
}

double Cylinder2D::get_drag(bool weak) const {
  const double mean_velocity = get_mean_velocity();
  const double force = weak ? drag_force_weak : drag_force;
  return 2.0 * force / (ro * mean_velocity * mean_velocity * D);
}

double Cylinder2D::get_lift(bool weak) const {
  const double mean_velocity = get_mean_velocity();
  const double force = weak ? lift_force_weak : lift_force;
  return 2.0 * force / (ro * mean_velocity * mean_velocity * D);
}

double Cylinder3D::get_drag(bool weak) const {
  const double mean_velocity = get_mean_velocity();
  const double force = weak ? drag_force_weak : drag_force;
  return 2.0 * force / (ro * mean_velocity * mean_velocity * D * H);
}

double Cylinder3D::get_lift(bool weak) const {
  const double mean_velocity = get_mean_velocity();
  const double force = weak ? lift_force_weak : lift_force;
  return 2.0 * force / (ro * mean_velocity * mean_velocity * D * H);
}