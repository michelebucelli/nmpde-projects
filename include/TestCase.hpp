#ifndef TEST_CASE_HPP
#define TEST_CASE_HPP

#include "NavierStokes.hpp"

class Cylinder2DNavierStokes : public NavierStokes<2U> {
private:
  // Physical dimension.
  constexpr static unsigned int dim = 2;
public:
  // Function for inlet velocity in the 2D flow past a cyclinder test case. 
  // This actually returns an object with four
  // components (one for each velocity component, and one for the pressure), but
  // then only the first three are really used (see the component mask when
  // applying boundary conditions at the end of assembly). If we only return
  // three components, however, we may get an error message due to this function
  // being incompatible with the finite element space.
  class InletVelocity : public Function<dim>
  {
  public:
    InletVelocity(double U_m_, double H_)
      : Function<dim>(dim + 1), U_m(U_m_), H(H_)
    {}

    virtual void
    vector_value(const Point<dim> &p, Vector<double> &values) const override;

    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

  private:
    const double U_m;
    const double H;
  };

  // Reynolds number computation.
  double
  get_reynolds_number() const override;

  // Constructor.
  Cylinder2DNavierStokes(
        const std::string  &mesh_file_name_,
        const unsigned int &degree_velocity_,
        const unsigned int &degree_pressure_,
        const double &T_,
        const double &deltat_);

private:
  // Reference velocity.
  constexpr static double U_m = 1.5;
  // Pressure for Neumann BC.
  constexpr static double p_out = 10.0;
  // Cyclinder diameter [m].
  constexpr static double D = 0.1;
  // Inlet side length [m].
  constexpr static double H = 0.41;

  // Inlet velocity.
  const InletVelocity inlet_velocity;
  // Zero function.
  const Functions::ZeroFunction<dim> zero_function;
};



class Cylinder3DNavierStokes : public NavierStokes<3U> {
private:
  // Physical dimension.
  constexpr static unsigned int dim = 3;
public:
  // Function for inlet velocity in the 3D flow past a cyclinder test case. 
  class InletVelocity : public Function<dim>
  {
  public:
    InletVelocity(double U_m_, double H_)
      : Function<dim>(dim + 1), U_m(U_m_), H(H_)
    {}

    virtual void
    vector_value(const Point<dim> &p, Vector<double> &values) const override;

    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

  private:
    const double U_m;
    const double H;
  };

  // Reynolds number computation.
  double
  get_reynolds_number() const override;

  // Constructor.
  Cylinder3DNavierStokes(const std::string  &mesh_file_name_,
         const unsigned int &degree_velocity_,
         const unsigned int &degree_pressure_,
         const double &T_,
         const double &deltat_);

private:
  // Reference velocity.
  constexpr static double U_m = 2.25;
  // Pressure for Neumann BC.
  constexpr static double p_out = 10.0;
  // Cyclinder diameter [m].
  constexpr static double D = 0.1;
  // Inlet side length [m].
  constexpr static double H = 0.41;

  // Inlet velocity.
  const InletVelocity inlet_velocity;
  // Zero function.
  const Functions::ZeroFunction<dim> zero_function;
};

#endif