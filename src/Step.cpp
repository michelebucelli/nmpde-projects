#include "Step.hpp"

Step::Step(const std::string &mesh_file_name_,
           const unsigned int &degree_velocity_,
           const unsigned int &degree_pressure_, const double &T_,
           const double &deltat_,
           const PreconditionerType &preconditioner_type_, const double &alpha_)
    : NavierStokes(mesh_file_name_, degree_velocity_, degree_pressure_, T_,
                   deltat_, preconditioner_type_),
      inlet_velocity(alpha_),
      zero_function(dim + 1) {
  ro = 1.0;
  nu = 1.0;
  initial_conditions = std::make_shared<Functions::ZeroFunction<dim>>(dim + 1);

  dirichlet_boundary_functions[0] = &inlet_velocity;

  neumann_boundary_functions[2] = &neumann_function;

  dirichlet_boundary_functions[1] = &zero_function;
}

double Step::InletVelocity::value(const Point<dim> &p,
                                  const unsigned int component) const {
  if (component == 0)
    return -alpha * p[1] * (2.0 - p[1]) * (1.0 - p[2]) * (2.0 - p[2]);
  else
    return 0.0;
}

void Step::InletVelocity::vector_value(const Point<dim> &p,
                                       Vector<double> &values) const {
  values[0] = value(p, 0);

  for (unsigned int i = 1; i < dim + 1; ++i) values[i] = 0.0;
}

double Step::NeumannFunction::value(const Point<dim> & /*p*/,
                                    const unsigned int component) const {
  if (component == 0)
    return -p_out;
  else
    return 0.0;
}

void Step::NeumannFunction::vector_value(const Point<dim> & /*p*/,
                                         Vector<double> &values) const {
  values[0] = -p_out;

  for (unsigned int i = 1; i < dim + 1; ++i) values[i] = 0.0;
}