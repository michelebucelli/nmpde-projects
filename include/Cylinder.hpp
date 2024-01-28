#ifndef CYLINDER_HPP
#define CYLINDER_HPP

#include "NavierStokes.hpp"

template <unsigned int dim>
class Cylinder : public NavierStokes<dim> {
 public:
  // Virtual destructor.
  virtual ~Cylinder() = default;

  // Update lift and drag forces.
  void update_lift_drag();

  // Apply initial conditions and set the header for the lift and drag file.
  void apply_initial_conditions() override;

  // Solve the problem for one time step and store the lift and drag
  // coefficients.
  void solve_time_step() override;

  // Compute and return the Reynolds number.
  double get_reynolds_number() const;

  // Return the lift coefficient.
  double get_lift() const;

  // Return the drag coefficient.
  double get_drag() const;

 protected:
  // Constructor.
  Cylinder(const std::string &mesh_file_name_,
           const unsigned int &degree_velocity_,
           const unsigned int &degree_pressure_, const double &T_,
           const double &deltat_, const double &U_m_,
           const SolverOptions &solver_options_);

  // Function to get the mean velocity.
  virtual double get_mean_velocity() const = 0;

  // Boundary tag for the obstacle.
  unsigned int obstacle_tag;
  // Reference velocity [m/s].
  const double U_m;
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
    virtual ~InletVelocity() = default;
    InletVelocity(double U_m_, double H_)
        : Function<dim>(dim + 1), U_m(U_m_), H(H_) {}
    virtual double value(const Point<dim> &p,
                         const unsigned int component = 0) const override = 0;
    virtual void vector_value(const Point<dim> &p,
                              Vector<double> &values) const override;

   protected:
    const double U_m;
    const double H;
  };

  // Class for the inlet velocity in tests 2D-1 and 2D-2.
  class TimeIndependentInletVelocity : public InletVelocity {
   public:
    TimeIndependentInletVelocity(double U_m_, double H_)
        : InletVelocity(U_m_, H_) {}
    double value(const Point<dim> &p,
                 const unsigned int component = 0) const override;
  };

  // Class for the inlet velocity in test 2D-3.
  class TimeDependentInletVelocity : public InletVelocity {
   public:
    TimeDependentInletVelocity(double U_m_, double H_)
        : InletVelocity(U_m_, H_) {}
    double value(const Point<dim> &p,
                 const unsigned int component = 0) const override;
  };

  // Constructor.
  Cylinder2D(const std::string &mesh_file_name_,
             const unsigned int &degree_velocity_,
             const unsigned int &degree_pressure_, const double &T_,
             const double &deltat_, const double &U_m_,
             const bool &time_dep_inlet, const SolverOptions &solver_options_);

 private:
  // Function to get the mean velocity.
  double get_mean_velocity() const override;

  // Function for the inlet velocity.
  std::shared_ptr<InletVelocity> inlet_velocity;
};

class Cylinder3D : public Cylinder<3U> {
 private:
  // Physical dimension.
  constexpr static unsigned int dim = 3;

 public:
  // Function for inlet velocity in the 3D flow past a cyclinder test case.
  class InletVelocity : public Function<dim> {
   public:
    virtual ~InletVelocity() = default;
    InletVelocity(double U_m_, double H_)
        : Function<dim>(dim + 1), U_m(U_m_), H(H_) {}
    virtual double value(const Point<dim> &p,
                         const unsigned int component = 0) const override = 0;
    virtual void vector_value(const Point<dim> &p,
                              Vector<double> &values) const override;

   protected:
    const double U_m;
    const double H;
  };

  // Class for the inlet velocity in tests 3D-1 and 3D-2.
  class TimeIndependentInletVelocity : public InletVelocity {
   public:
    TimeIndependentInletVelocity(double U_m_, double H_)
        : InletVelocity(U_m_, H_) {}
    double value(const Point<dim> &p,
                 const unsigned int component = 0) const override;
  };

  // Class for the inlet velocity in test 3D-3.
  class TimeDependentInletVelocity : public InletVelocity {
   public:
    TimeDependentInletVelocity(double U_m_, double H_)
        : InletVelocity(U_m_, H_) {}
    double value(const Point<dim> &p,
                 const unsigned int component = 0) const override;
  };

  // Constructor.
  Cylinder3D(const std::string &mesh_file_name_,
             const unsigned int &degree_velocity_,
             const unsigned int &degree_pressure_, const double &T_,
             const double &deltat_, const double &U_m_,
             const bool &time_dep_inlet, const SolverOptions &solver_options_);

 private:
  // Function to get the mean velocity.
  double get_mean_velocity() const override;

  // Function for the inlet velocity.
  std::shared_ptr<InletVelocity> inlet_velocity;
};

#endif