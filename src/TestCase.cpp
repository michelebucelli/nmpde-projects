#include "TestCase.hpp"

void Cylinder2DNavierStokes::InletVelocity::vector_value(const Point<dim> &p, Vector<double> &values) const {
  values[0] = 4.0 * U_m * p[1] * (H - p[1]) / (H * H);

  for (unsigned int i = 1; i < dim + 1; ++i)
    values[i] = 0.0;
}

void Cylinder3DNavierStokes::InletVelocity::vector_value(const Point<dim> &p, Vector<double> &values) const {
  values[0] = 16.0 * U_m * p[1] * p[2] * (H - p[1]) * (H - p[2]) / (H * H * H * H);

  for (unsigned int i = 1; i < dim + 1; ++i)
    values[i] = 0.0;
}

double Cylinder2DNavierStokes::InletVelocity::value(const Point<dim> &p, const unsigned int component) const
{
  if (component == 0) {
    return 4.0 * U_m * p[1] * (H - p[1]) / (H * H);
  } else {
    return 0.0;
  }
}

double Cylinder3DNavierStokes::InletVelocity::value(const Point<dim> &p, const unsigned int component) const
{
  if (component == 0) {
    return 16.0 * U_m * p[1] * p[2] * (H - p[1]) * (H - p[2]) / (H * H * H * H);
  } else {
    return 0.0;
  }
}

double Cylinder2DNavierStokes::get_reynolds_number() const
{
  return 2.0 * inlet_velocity.value(Point<dim>(0, H/2.0), 0) / 3.0 * D / nu;
}

double Cylinder3DNavierStokes::get_reynolds_number() const
{
  return 4.0 * inlet_velocity.value(Point<dim>(0, H/2.0, H/2.0), 0) / 9.0 * D / nu;
}

Cylinder2DNavierStokes::Cylinder2DNavierStokes(
         const std::string  &mesh_file_name_,
         const unsigned int &degree_velocity_,
         const unsigned int &degree_pressure_,
         const double &T_,
         const double &deltat_)
    : NavierStokes(mesh_file_name_, degree_velocity_, degree_pressure_, T_, deltat_)
    , inlet_velocity(U_m, H)
    , zero_function(dim+1) {
  ro = 1.0;
  nu = 1e-3;

  initial_conditions = std::make_shared<Functions::ZeroFunction<dim>>(dim+1);

  dirichlet_boundary_functions[0] = &inlet_velocity;
  dirichlet_boundary_functions[1] = &zero_function;

  neumann_boundary_functions[2] = p_out;
}

Cylinder3DNavierStokes::Cylinder3DNavierStokes(
         const std::string  &mesh_file_name_,
         const unsigned int &degree_velocity_,
         const unsigned int &degree_pressure_,
         const double &T_,
         const double &deltat_)
    : NavierStokes(mesh_file_name_, degree_velocity_, degree_pressure_, T_, deltat_)
    , inlet_velocity(U_m, H)
    , zero_function(dim+1) {
  ro = 1.0;
  nu = 1e-3;

  initial_conditions = std::make_shared<Functions::ZeroFunction<dim>>(dim+1);
  
  dirichlet_boundary_functions[0] = &inlet_velocity;
  dirichlet_boundary_functions[1] = &zero_function;

  neumann_boundary_functions[2] = p_out;
}