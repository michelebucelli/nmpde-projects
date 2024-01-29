#ifndef STEP_HPP
#define STEP_HPP

#include "NavierStokes.hpp"

class Step : public NavierStokes<3U> {
 private:
  // This is the physical dimension of the problem. We need to specify it here
  // because we need it to be a constant expression. If this were not the case,
  // we would have considered dropping the template approach.
  constexpr static unsigned int dim = 3;

 public:
  // Function for the inlet velocity in the step test case.
  class InletVelocity : public Function<dim> {
   public:
    // Constructor.
    InletVelocity(double alpha_) : Function<dim>(dim + 1), alpha(alpha_) {}

    // When defining vector-valued functions, we need to define the value
    // function, which returns the value of the function at a given point and
    // component...
    virtual double value(const Point<dim> &p,
                         const unsigned int component = 0) const override;

    // ... and the vector_value function, which returns the value of the
    // function at a given point for all components.
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

    // Evaluation.
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
       const double &deltat_, const SolverOptions &solver_options_,
       const double &alpha_);

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