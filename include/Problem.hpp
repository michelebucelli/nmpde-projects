#ifndef PROBLEM_HPP
#define PROBLEM_HPP

#include <deal.II/numerics/vector_tools.h>

#include "NavierStokes.hpp"

class Cylinder2D : public NavierStokes<2U> {
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

  // Reynolds number computation.
  double get_reynolds_number() const;
  void calc_lift_drag();

  // Constructor.
  Cylinder2D(const std::string &mesh_file_name_,
             const unsigned int &degree_velocity_,
             const unsigned int &degree_pressure_, const double &T_,
             const double &deltat_);

 private:
  constexpr static int OBSTACLE_ID = 5;
  // Reference velocity.
  constexpr static double U_m = 1.5;
  // Cyclinder diameter [m].
  constexpr static double D = 0.1;
  // Inlet side length [m].
  constexpr static double H = 0.41;

  double _lift = 0.0;
  double _drag = 0.0;

  // Inlet velocity.
  InletVelocity inlet_velocity;
  // Zero function.
  Functions::ZeroFunction<dim> zero_function;
};

class Cylinder3D : public NavierStokes<3U> {
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

  // Reynolds number computation.
  double get_reynolds_number() const;

  // Constructor.
  Cylinder3D(const std::string &mesh_file_name_,
             const unsigned int &degree_velocity_,
             const unsigned int &degree_pressure_, const double &T_,
             const double &deltat_);

 private:
  // Reference velocity.
  constexpr static double U_m = 2.25;
  // Cyclinder diameter [m].
  constexpr static double D = 0.1;
  // Inlet side length [m].
  constexpr static double H = 0.41;

  // Inlet velocity.
  InletVelocity inlet_velocity;
  // Zero function.
  Functions::ZeroFunction<dim> zero_function;
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

  // Compute the error (if block=0, the velocity error is computed, if block=1,
  // the pressure one is).
  double compute_error(const VectorTools::NormType &norm_type,
                       unsigned int block);

 private:
  // Parameters of the problem.
  static constexpr double a = M_PI / 4.0;
  static constexpr double b = M_PI / 2.0;
  // Exact solution.
  ExactSolution exact_solution;
  // Neumann function on y=-1.
  NeumannFunction neumann_function;
};

#endif