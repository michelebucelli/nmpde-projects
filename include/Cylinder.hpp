#ifndef CYLINDER_HPP
#define CYLINDER_HPP

#include "NavierStokes.hpp"

// In this file we define the classes for the flow past a cylinder test cases.
// To do this, we create a template class for the cylinder test case, which
// contains the common features of the test cases. This class is a template on
// the dimension, because it needs to be defined as constexpr for technical
// reasons. Then, we create two derived classes, one for the 2D case and one
// for the 3D case. These classes contain the specific features of the test
// cases, such as the inlet velocity and the initial conditions.

// This is a template class for the flow past a cylinder test case. It is
// templated on the dimension, which can be either 2 or 3. The reason for this
// is that we need the dimension to be known at compile time. The class is
// derived from the NavierStokes class, which is a template class itself, for
// the same reason. A similar approach is used in most of our other classes.
template <unsigned int dim>
class Cylinder : public NavierStokes<dim> {
 public:
  // Virtual destructor.
  virtual ~Cylinder() = default;

  // Update lift and drag forces. In particular, this function computes the
  // lift and drag coefficients and prints them to the console. These
  // coefficients will be added to a .csv file and plotted in a Python script.
  void update_lift_drag();

  // Like the previous function, but using the more precise method described in
  // https://www.mate.polimi.it/biblioteca/add/qmox/mox84.pdf.
  void update_lift_drag_weak();

  // This function applies the initial conditions for the problem, and also
  // sets the header for the lift and drag file if needed.
  void apply_initial_conditions() override;

  // This function performs the time stepping. It just performs one time step
  // and then stores the lift and drag coefficients if needed.
  void solve_time_step() override;

  // This function uses the definition of the Reynolds number to compute it.
  double get_reynolds_number() const;

  // Return the lift coefficient in the current configuration. If weak is true,
  // return the one calculated with the weak formulation.
  virtual double get_lift(bool weak) const = 0;

  // Return the drag coefficient in the current configuration. If weak is true,
  // return the one calculated with the weak formulation.
  virtual double get_drag(bool weak) const = 0;

 protected:
  // Constructor.
  Cylinder(const std::string &mesh_file_name_,
           const unsigned int &degree_velocity_,
           const unsigned int &degree_pressure_, const double &T_,
           const double &deltat_, const double &U_m_,
           const SolverOptions &solver_options_,
           const bool &compute_lift_drag_);

  // Function to get the mean velocity.
  virtual double get_mean_velocity() const = 0;

  // Function to calculate the two values of phi_inf.
  void calculate_phi_inf();

  // Path to the drag and lift coefficient logging file.
  const std::string lift_drag_output_file = "../results/lift_drag.csv";

  // When we apply the boundary conditions, we need to specify which boundary
  // is which. This is done by assigning a tag to each boundary. The tag is
  // just an integer, and we can assign it to each boundary in Gmsh. The
  // following constant is specifically used to tag the obstacle boundary.
  unsigned int obstacle_tag;

  // This is the reference velocity, expressed in meters per second.
  const double U_m;

  // This is the diameter of the cylinder, expressed in meters. Note that
  // by cylinder we mean the obstacle, which is a circle in 2D.
  const double D = 0.1;

  // This is the height of the domain, or the length of the inlet boundary,
  // expressed in meters.
  const double H = 0.41;

  // Whether to compute the lift and drag coefficients.
  const bool compute_lift_drag;

  // This is the lift force, expressed in Newtons. It is the component of the
  // force that is exerted on the cylinder in the direction perpendicular to
  // the flow. In the 3D case, "perpendicular" refers to the y direction.
  double lift_force;

  // This is the drag force, expressed in Newtons. It is the component of the
  // force that is exerted on the cylinder in the direction parallel to the
  // flow.
  double drag_force;

  // Same as the previous force, but calculated with the weak formulation.
  double lift_force_weak;

  // Same as the previous force, but calculated with the weak formulation.
  double drag_force_weak;

  // Zero function. This is a handful function that is used to set some boundary
  // conditions to zero.
  Functions::ZeroFunction<dim> zero_function;

  // Function for the inlet velocity.
  std::shared_ptr<Function<dim>> inlet_velocity;

  // Phi_inf (including ghost elements), used for weak calculation of the lift
  // coefficient.
  TrilinosWrappers::MPI::BlockVector phi_inf_lift;

  // Phi_inf (including ghost elements), used for weak calculation of the drag
  // coefficient.
  TrilinosWrappers::MPI::BlockVector phi_inf_drag;
};

class Cylinder2D : public Cylinder<2U> {
 private:
  // This is the physical dimension of the problem. We need to specify it here
  // because we need it to be a constant expression. If this were not the case,
  // we would have considered dropping the template approach.
  constexpr static unsigned int dim = 2;

 public:
  // This is a function that defines the inlet velocity in the 2D flow past a
  // cyclinder test case. The returned object has three components (one for each
  // dimensional component, and one for the pressure), but we only use the first
  // two (see the component mask when applying boundary conditions at the end
  // of assembly). If we only returned two components, however, we would get
  // an error message due to this function being incompatible with the finite
  // element space.
  class InletVelocity : public Function<dim> {
   public:
    // Virtual destructor.
    virtual ~InletVelocity() = default;
    // Constructor.
    InletVelocity(double U_m_, double H_)
        : Function<dim>(dim + 1), U_m(U_m_), H(H_) {}

    // When defining vector-valued functions, we need to define the value
    // function, which returns the value of the function at a given point and
    // component...
    virtual double value(const Point<dim> &p,
                         const unsigned int component = 0) const override = 0;

    // ... and the vector_value function, which returns the value of the
    // function at a given point for all components.
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
             const bool &time_dep_inlet, const SolverOptions &solver_options_,
             const bool &compute_lift_drag_);

 private:
  // Function to get the mean velocity.
  double get_mean_velocity() const override;

  // Function to get the lift coefficient.
  double get_lift(bool weak) const override;

  // Function to get the drag coefficient.
  double get_drag(bool weak) const override;
};

class Cylinder3D : public Cylinder<3U> {
 private:
  // This is the physical dimension of the problem. We need to specify it here
  // because we need it to be a constant expression. If this were not the case,
  // we would have considered dropping the template approach.
  constexpr static unsigned int dim = 3;

 public:
  // Function for inlet velocity in the 3D flow past a cyclinder test case.
  class InletVelocity : public Function<dim> {
   public:
    // Virtual destructor.
    virtual ~InletVelocity() = default;
    // Constructor.
    InletVelocity(double U_m_, double H_)
        : Function<dim>(dim + 1), U_m(U_m_), H(H_) {}

    // Value function (returns the value of the function at a given point and
    // component).
    virtual double value(const Point<dim> &p,
                         const unsigned int component = 0) const override = 0;

    // Vector value function (returns the value of the function at a given point
    // for all components).
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
             const bool &time_dep_inlet, const SolverOptions &solver_options_,
             const bool &compute_lift_drag_);

 private:
  // Function to get the mean velocity.
  double get_mean_velocity() const override;

  // Function to get the lift coefficient.
  double get_lift(bool weak) const override;

  // Function to get the drag coefficient.
  double get_drag(bool weak) const override;
};

#endif