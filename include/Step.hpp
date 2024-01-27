#ifndef STEP_HPP
#define STEP_HPP

#include "NavierStokes.hpp"

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