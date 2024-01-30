#ifndef ETHIER_STEINMAN_HPP
#define ETHIER_STEINMAN_HPP

#include <deal.II/numerics/vector_tools.h>

#include "NavierStokes.hpp"

// We implement the Ethier-Steinman test case, as a means to test our code for
// correctness. This is because the exact solution of this specific problem is
// known. The problem is defined on the domain (-1,1)^3.

// This is a class for the Ethier-Steinman test case. The class is derived from
// the NavierStokes class, which is a template class on the physical dimension,
// which needs to be a constant expression. The Ethier-Steinman problem is in 3D
// and therefore inherits from NaverStokes<3U>.
class EthierSteinman : public NavierStokes<3U> {
 private:
  // This is the physical dimension of the problem. We need to specify it here
  // because we need it to be a constant expression. If this were not the case,
  // we would have considered dropping the template approach.
  constexpr static unsigned int dim = 3;

 public:
  // Constructor.
  EthierSteinman(const std::string &mesh_file_name_,
                 const unsigned int &degree_velocity_,
                 const unsigned int &degree_pressure_, const double &T_,
                 const double &deltat_, const SolverOptions &solver_options_,
                 const double &nu_);

  // Class for the exact solution, containing both velocity and pressure.
  class ExactSolution : public Function<dim> {
   public:
    // Constructor.
    ExactSolution(double nu_)
        : Function<dim>(dim + 1),
          nu(nu_),
          exact_velocity(nu_),
          exact_pressure(nu_) {}

    // When defining vector-valued functions, we need to define the value
    // function, which returns the value of the function at a given point and
    // component...
    virtual double value(const Point<dim> &p,
                         const unsigned int component) const override;

    // ... and the vector_value function, which returns the value of the
    // function at a given point for all components.
    virtual void vector_value(const Point<dim> &p,
                              Vector<double> &values) const override;

    // This is a function object, which defines the exact velocity. Since the
    // problem's exact solution is known, we can define it as a function object
    // and use it to compute the error of our numerical solution. To be able to
    // compute the H1 norm of the error, the exact gradient is computed as well.
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

    // Same as above, for the pressure. Note that the pressure is a scalar
    // function, but this function returns a vector with 4 components, of
    // which the first three are empty.
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
    // This is a function object, which defines the exact velocity. Since we
    // need to set the time, it's defined as mutable.
    mutable ExactVelocity exact_velocity;

    // This is a function object, which defines the exact pressure. Since we
    // need to set the time, it's defined as mutable.
    mutable ExactPressure exact_pressure;
  };

  // Class for the Neumann function, needed to impose the boundary conditions.
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
    // This is a function object, which defines the exact solution. Since we
    // need to set the time, it's defined as mutable.
    mutable ExactSolution exact_solution;
  };

  // Apply the initial conditions and print the error.
  void apply_initial_conditions() override;

  // Solve the problem for one time step and print the error.
  void solve_time_step() override;

  // Compute the error (if velocity is true, the error on the velocity is
  // computed, otherwise the error on the pressure is). Note that the pressure
  // does not support the H1 norm.
  double compute_error(const VectorTools::NormType &norm_type, bool velocity);

 private:
  // Parameters of the problem.
  static constexpr double a = M_PI / 4.0;
  static constexpr double b = M_PI / 2.0;

  // Exact solution.
  ExactSolution exact_solution;

  // Neumann function on y = -1.
  NeumannFunction neumann_function;
};

#endif