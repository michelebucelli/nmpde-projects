#ifndef CYLINDER_HPP
#define CYLINDER_HPP

#include "NavierStokes.hpp"

template <unsigned int dim>
class Cylinder : public NavierStokes<dim> {
 public:
  Cylinder(const std::string &mesh_file_name_,
           const unsigned int &degree_velocity_,
           const unsigned int &degree_pressure_, const double &T_,
           const double &deltat_,
           const PreconditionerType &preconditioner_type_);

  // Virtual destructor.
  virtual ~Cylinder() = default;

  // Update lift and drag forces.
  void update_lift_drag();

  // Compute and return the Reynolds number.
  double get_reynolds_number() const;

  // Return the lift coefficient.
  double get_lift() const;

  // Return the drag coefficient.
  double get_drag() const;

 protected:
  // Function to get the mean velocity.
  virtual double get_mean_velocity() const = 0;

  // Boundary tag for the obstacle.
  unsigned int obstacle_tag;
  // Reference velocity [m/s].
  double U_m;
  // Cyclinder diameter [m].
  const double D = 0.1;
  // Inlet side length [m].
  const double H = 0.41;
  // Lift force [N].
  double lift_force;
  // Drag force [N].
  double drag_force;
  // Zero function.
  Functions::ZeroFunction<dim> zero_function;
};

class Cylinder2D : public Cylinder<2U> {
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
  class InletVelocity : public Function<dim> {
   public:
    InletVelocity(double U_m_, double H_)
        : Function<dim>(dim + 1), U_m(U_m_), H(H_) {}

    virtual double value(const Point<dim> &p,
                         const unsigned int component = 0) const override;

    virtual void vector_value(const Point<dim> &p,
                              Vector<double> &values) const override;

   private:
    const double U_m;
    const double H;
  };

  // Constructor.
  Cylinder2D(const std::string &mesh_file_name_,
             const unsigned int &degree_velocity_,
             const unsigned int &degree_pressure_, const double &T_,
             const double &deltat_,
             const PreconditionerType &preconditioner_type_);

 private:
  // Function to get the mean velocity.
  double get_mean_velocity() const override;

  // Function for the inlet velocity.
  InletVelocity inlet_velocity;
};

class Cylinder3D : public Cylinder<3U> {
 private:
  // Physical dimension.
  constexpr static unsigned int dim = 3;

 public:
  // Function for inlet velocity in the 3D flow past a cyclinder test case.
  class InletVelocity : public Function<dim> {
   public:
    InletVelocity(double U_m_, double H_)
        : Function<dim>(dim + 1), U_m(U_m_), H(H_) {}

    virtual double value(const Point<dim> &p,
                         const unsigned int component = 0) const override;

    virtual void vector_value(const Point<dim> &p,
                              Vector<double> &values) const override;

   private:
    const double U_m;
    const double H;
  };

  // Constructor.
  Cylinder3D(const std::string &mesh_file_name_,
             const unsigned int &degree_velocity_,
             const unsigned int &degree_pressure_, const double &T_,
             const double &deltat_,
             const PreconditionerType &preconditioner_type_);

 private:
  // Function to get the mean velocity.
  double get_mean_velocity() const override;

  // Function for the inlet velocity.
  InletVelocity inlet_velocity;
};

#endif