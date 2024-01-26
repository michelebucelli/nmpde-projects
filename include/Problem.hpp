#ifndef PROBLEM_HPP
#define PROBLEM_HPP

#include <deal.II/numerics/vector_tools.h>

#include "NavierStokes.hpp"

template <unsigned int dim>
class Cylinder : public NavierStokes<dim> {
 public:
  Cylinder(const std::string &mesh_file_name_,
           const unsigned int &degree_velocity_,
           const unsigned int &degree_pressure_, const double &T_,
           const double &deltat_);

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
             const double &deltat_);

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
             const double &deltat_);

 private:
  // Function to get the mean velocity.
  double get_mean_velocity() const override;

  // Function for the inlet velocity.
  InletVelocity inlet_velocity;
};

class EthierSteinman : public NavierStokes<3U> {
 private:
  // Physical dimension.
  constexpr static unsigned int dim = 3;

 public:
  // Constructor.
  EthierSteinman(const std::string &mesh_file_name_,
                 const unsigned int &degree_velocity_,
                 const unsigned int &degree_pressure_, const double &T_,
                 const double &deltat_, const double &nu_);

  // Class for the exact solution, containing both velocity and pressure.
  class ExactSolution : public Function<dim> {
   public:
    // Constructor.
    ExactSolution(double nu_)
        : Function<dim>(dim + 1),
          nu(nu_),
          exact_velocity(nu_),
          exact_pressure(nu_) {}

    // Evaluation.
    virtual double value(const Point<dim> &p,
                         const unsigned int component) const override;
    virtual void vector_value(const Point<dim> &p,
                              Vector<double> &values) const override;

    // Class for the exact velocity.
    class ExactVelocity : public Function<dim> {
     public:
      // Constructor.
      ExactVelocity(double nu_) : Function<dim>(dim), nu(nu_) {}
      // Evaluation.
      virtual double value(const Point<dim> &p,
                           const unsigned int component) const override;
      virtual void vector_value(const Point<dim> &p,
                                Vector<double> &values) const override;
      virtual Tensor<1, dim> gradient(
          const Point<dim> &p, const unsigned int component) const override;
      virtual void vector_gradient(
          const Point<dim> &p,
          std::vector<Tensor<1, dim>> &values) const override;

     private:
      const double nu;
    };

    // Class for the exact pressure.
    // This returns a vector with 4 components, of which the first three are
    // empty.
    class ExactPressure : public Function<dim> {
     public:
      // Constructor.
      ExactPressure(double nu_) : Function<dim>(1), nu(nu_) {}
      // Evaluation.
      virtual double value(const Point<dim> &p,
                           const unsigned int /*component*/ = 0) const override;
      virtual void vector_value(const Point<dim> &p,
                                Vector<double> &values) const override;

     private:
      const double nu;
    };

   private:
    // Diffusion coefficient.
    const double nu;

   public:
    // Exact velocity, mutable to set the time.
    mutable ExactVelocity exact_velocity;
    // Exact pressure, mutable to set the time.
    mutable ExactPressure exact_pressure;
  };

  class NeumannFunction : public Function<dim> {
   public:
    // Constructor.
    NeumannFunction(double nu_)
        : Function<dim>(dim + 1), nu(nu_), exact_solution(nu_) {}
    // Evaluation.
    virtual double value(const Point<dim> &p,
                         const unsigned int component) const override;
    virtual void vector_value(const Point<dim> &p,
                              Vector<double> &values) const override;

   private:
    const double nu;

   public:
    // Exact solution, mutable to set the time.
    mutable ExactSolution exact_solution;
  };

  // Apply the initial conditions and print the error.
  void apply_initial_conditions() override;

  // Solve the problem for one time step and print the error.
  void solve_time_step() override;

  // Compute the error (if velocity is true, the error on the velocity is
  // computed, otherwise the error on the pressure is, the pressure does not
  // support the H1 norm).
  double compute_error(const VectorTools::NormType &norm_type, bool velocity);

 private:
  // Parameters of the problem.
  static constexpr double a = M_PI / 4.0;
  static constexpr double b = M_PI / 2.0;
  // Exact solution.
  ExactSolution exact_solution;
  // Neumann function on y=-1.
  NeumannFunction neumann_function;
};

class Step : public NavierStokes<3U> {
 private:
  // Physical dimension.
  constexpr static unsigned int dim = 3;

 public:
  // Function for the inlet velocity in the step test case.
  class InletVelocity : public Function<dim> {
   public:
    // Constructor.
    InletVelocity(double alpha_) : Function<dim>(dim + 1), alpha(alpha_) {}

    virtual double value(const Point<dim> &p,
                         const unsigned int component = 0) const override;

    virtual void vector_value(const Point<dim> &p,
                              Vector<double> &values) const override;

   private:
    double alpha;
  };

  // Function for the Neumann (BC) in the step test case.
  class NeumannFunction : public Function<dim> {
   public:
    // Constructor.
    NeumannFunction() : Function<dim>(dim + 1) {}

    virtual double value(const Point<dim> &p,
                         const unsigned int component = 0) const override;

    virtual void vector_value(const Point<dim> &p,
                              Vector<double> &values) const override;

   private:
    static constexpr double p_out = 10.0;
  };

  // Constructor.
  Step(const std::string &mesh_file_name_, const unsigned int &degree_velocity_,
       const unsigned int &degree_pressure_, const double &T_,
       const double &deltat_, const double &alpha_);

 private:
  // Inlet velocity.
  InletVelocity inlet_velocity;
  // Zero function.
  Functions::ZeroFunction<dim> zero_function;
  // Coefficient for inlet velocity.
  double alpha;
  // Function for Neumann (BC).
  NeumannFunction neumann_function;
};

#endif