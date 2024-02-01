#include "Cylinder.hpp"

#include <fstream>

template class Cylinder<2U>;
template class Cylinder<3U>;

template <unsigned int dim>
Cylinder<dim>::Cylinder(const std::string &mesh_file_name_,
                        const unsigned int &degree_velocity_,
                        const unsigned int &degree_pressure_, const double &T_,
                        const double &deltat_, const double &U_m_,
                        const SolverOptions &solver_options_,
                        const bool &compute_lift_drag_)
    : NavierStokes<dim>(mesh_file_name_, degree_velocity_, degree_pressure_, T_,
                        deltat_, solver_options_),
      U_m(U_m_),
      compute_lift_drag(compute_lift_drag_),
      zero_function(dim + 1) {
  NavierStokes<dim>::ro = 1.0;   //[kg/m^3].
  NavierStokes<dim>::nu = 1e-3;  //[m^2/s].
  NavierStokes<dim>::initial_conditions =
      std::make_shared<Functions::ZeroFunction<dim>>(dim + 1);
}

template <unsigned int dim>
void Cylinder<dim>::apply_initial_conditions() {
  NavierStokes<dim>::apply_initial_conditions();

  // Set the header for the lift and drag file.
  if (NavierStokes<dim>::mpi_rank == 0 && compute_lift_drag) {
    std::ofstream file;
    file.open(lift_drag_output_file);
    file << "time_step(delta_t=" << NavierStokes<dim>::deltat
         << "s),lift_coefficient,drag_coefficient,lift_coefficient_weak,drag_"
            "coefficient_weak,reynolds_number\n";
    file.close();
  }
}

template <unsigned int dim>
void Cylinder<dim>::solve_time_step() {
  NavierStokes<dim>::solve_time_step();

  if (compute_lift_drag) {
    // Update lift and drag coefficients.
    if (NavierStokes<dim>::time_step == 1) {
      calculate_phi_inf();
    }
    update_lift_drag();
    update_lift_drag_weak();

    // Write the results to a file.
    if (NavierStokes<dim>::mpi_rank == 0) {
      std::ofstream file;
      file.open(lift_drag_output_file, std::ios::app);
      file << NavierStokes<dim>::time_step << "," << get_lift(false) << ","
           << get_drag(false) << "," << get_lift(true) << "," << get_drag(true)
           << "," << get_reynolds_number() << "\n";
      file.close();
    }
  }
}

double Cylinder2D::TimeIndependentInletVelocity::value(
    const Point<dim> &p, const unsigned int component) const {
  if (component == 0) {
    return 4.0 * U_m * p[1] * (H - p[1]) / (H * H);
  } else {
    return 0.0;
  }
}

double Cylinder3D::TimeDependentInletVelocity::value(
    const Point<dim> &p, const unsigned int component) const {
  if (component == 0) {
    return 16.0 * U_m * p[1] * p[2] * (H - p[1]) * (H - p[2]) /
           (H * H * H * H) * std::sin(M_PI * get_time() / 8.0);
  } else {
    return 0.0;
  }
}

double Cylinder3D::TimeIndependentInletVelocity::value(
    const Point<dim> &p, const unsigned int component) const {
  if (component == 0) {
    return 16.0 * U_m * p[1] * p[2] * (H - p[1]) * (H - p[2]) / (H * H * H * H);
  } else {
    return 0.0;
  }
}

double Cylinder2D::TimeDependentInletVelocity::value(
    const Point<dim> &p, const unsigned int component) const {
  if (component == 0) {
    return 4.0 * U_m * p[1] * (H - p[1]) / (H * H) *
           std::sin(M_PI * get_time() / 8.0);
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
  return 2.0 * inlet_velocity->value(Point<dim>(0, H / 2.0), 0) / 3.0;
}

double Cylinder3D::get_mean_velocity() const {
  return 4.0 * inlet_velocity->value(Point<dim>(0, H / 2.0, H / 2.0), 0) / 9.0;
}

Cylinder2D::Cylinder2D(const std::string &mesh_file_name_,
                       const unsigned int &degree_velocity_,
                       const unsigned int &degree_pressure_, const double &T_,
                       const double &deltat_, const double &U_m_,
                       const bool &time_dep_inlet,
                       const SolverOptions &solver_options_,
                       const bool &compute_lift_drag_)
    : Cylinder<dim>(mesh_file_name_, degree_velocity_, degree_pressure_, T_,
                    deltat_, U_m_, solver_options_, compute_lift_drag_) {
  pcout << "Solving 2D flow past a cylinder problem" << std::endl;

  if (time_dep_inlet) {
    inlet_velocity = std::make_shared<TimeDependentInletVelocity>(U_m, H);
  } else {
    inlet_velocity = std::make_shared<TimeIndependentInletVelocity>(U_m, H);
  }

  obstacle_tag = 5U;

  dirichlet_boundary_functions[4] = &(*inlet_velocity);

  neumann_boundary_functions[2] = &zero_function;

  dirichlet_boundary_functions[1] = &zero_function;
  dirichlet_boundary_functions[3] = &zero_function;
  dirichlet_boundary_functions[5] = &zero_function;
}

Cylinder3D::Cylinder3D(const std::string &mesh_file_name_,
                       const unsigned int &degree_velocity_,
                       const unsigned int &degree_pressure_, const double &T_,
                       const double &deltat_, const double &U_m_,
                       const bool &time_dep_inlet,
                       const SolverOptions &solver_options_,
                       const bool &compute_lift_drag_)
    : Cylinder<dim>(mesh_file_name_, degree_velocity_, degree_pressure_, T_,
                    deltat_, U_m_, solver_options_, compute_lift_drag_) {
  pcout << "Solving 3D flow past a cylinder problem" << std::endl;

  if (time_dep_inlet) {
    inlet_velocity = std::make_shared<TimeDependentInletVelocity>(U_m, H);
  } else {
    inlet_velocity = std::make_shared<TimeIndependentInletVelocity>(U_m, H);
  }

  obstacle_tag = 7U;

  dirichlet_boundary_functions[5] = &(*inlet_velocity);

  neumann_boundary_functions[3] = &zero_function;

  dirichlet_boundary_functions[1] = &zero_function;
  dirichlet_boundary_functions[2] = &zero_function;
  dirichlet_boundary_functions[4] = &zero_function;
  dirichlet_boundary_functions[6] = &zero_function;
  dirichlet_boundary_functions[7] = &zero_function;
}