#ifndef ETHIER_STEINMAN_HPP
#define ETHIER_STEINMAN_HPP

#include "NavierStokes.hpp"
#include <deal.II/numerics/vector_tools.h>

class EthierSteinman : public NavierStokes<3U> {
 private:
  // Physical dimension.
  constexpr static unsigned int dim = 3;

 public:
  // Constructor.
  EthierSteinman(const std::string &mesh_file_name_,
                 const unsigned int &degree_velocity_,
                 const unsigned int &degree_pressure_, const double &T_,
                 const double &deltat_, const PreconditionerType & preconditioner_type_, const double &nu_);

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

#endif