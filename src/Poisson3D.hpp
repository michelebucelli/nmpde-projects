#ifndef POISSON_3D_HPP
#define POISSON_3D_HPP

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/lac/eigen.h> //For power method
#include <deal.II/lac/vector_memory.h> // Power method

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <filesystem>
#include <fstream>
#include <iostream>

#include <vector>
#include <cmath> // For pow and sqrt functions
#include <random>

#include <unistd.h> // For command line arguments

using namespace dealii;

/**
 * Class managing the differential problem.
 */

class Poisson3D
{
public:
  // Physical dimension (1D, 2D, 3D)
  static constexpr unsigned int dim = 3;

  // Diffusion coefficient.
  // In deal.ii, functions are implemented by deriving the dealii::Function
  // class, which provides an interface for the computation of function values
  // and their derivatives.

  template <int dim>
  class DiffusionCoefficient : public Function<dim>
  {
  public:
      // Define a structure to hold sphere data
      struct Sphere
      {
          Point<dim> center;
          double radius;
      };
  
      // Constructor
      DiffusionCoefficient(const std::vector<Sphere>& spheres, double p)
      : spheres(spheres), p_value(p) {}

      DiffusionCoefficient(){}

      // Evaluation
      virtual double value(const Point<dim>& p, const unsigned int /*component*/ = 0) const override
      {
          for (const auto& sphere : spheres)
          {
              if (is_point_inside_sphere(p, sphere))
                  return std::pow(10.0, p_value);
          }
          return 1.0;
      }

  private:
      std::vector<Sphere> spheres;
      double p_value;

      bool is_point_inside_sphere(const Point<dim>& point, const Sphere& sphere) const
      {
          double distance_squared = 0.0;
          for (unsigned int i = 0; i < dim; ++i)
          {
              distance_squared += std::pow(point[i] - sphere.center[i], 2);
          }
          return distance_squared <= std::pow(sphere.radius, 2);
      }
  };

  // Forcing term.
  class ForcingTerm : public Function<dim>
  {
  public:
    // Constructor.
    ForcingTerm()
    {}

    // Evaluation.
    virtual double
    value(const Point<dim> & /*p*/, const unsigned int /*component*/ = 0) const
    {
      return 1.0;//We have to define the f we would want to use
    }
  };

  // Neumann boundary conditions.
  class FunctionH : public Function<dim>
  {
  public:
    // Constructor.
    FunctionH()
    {}

    // Evaluation:
    virtual double
    value(const Point<dim> /*&p*/, const unsigned int /*component*/ = 0) const
    {
      return 0.0;
    }
  };
// Dirichlet boundary conditions.
  class FunctionG : public Function<dim>
  {
  public:
    // Constructor.
    FunctionG()
    {}

    // Evaluation.
    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int /*component*/ = 0) const override
    {
      return 0.0;
    }
  };


  // Constructor.
  Poisson3D()
  {};



  bool isMatrixSingular(dealii::SparseMatrix<double> &matrix) {
      dealii::SolverControl solver_control(1000, 1e-12);
      dealii::SolverCG<> solver(solver_control);

      dealii::Vector<double> x(matrix.n()), b(matrix.m());
      b = 1; // Or any other non-zero vector

      try {
          solver.solve(matrix, x, b, dealii::PreconditionIdentity());
          return false;
      } catch (std::exception &e) {
          return true; // The matrix could be singular or poorly conditioned
      }
  }

  void initialize_diffusion_coefficient(double p_value) {
        // Create a vector of spheres
        std::vector<DiffusionCoefficient<dim>::Sphere> spheres;

        // Generate random spheres
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis_center(0.0, 1.0);
        std::uniform_real_distribution<> dis_radius(0.01, 0.1); // Small spokes to make sure they fit in the hub
        for (int i = 0; i < 10; ++i) {
            Point<dim> center = {dis_center(gen), dis_center(gen), dis_center(gen)};
            double radius = dis_radius(gen);
            spheres.push_back({center, radius});
        }



        // Initialize the diffusion coefficient
        diffusion_coefficient = DiffusionCoefficient<dim>(spheres, p_value);
  };

  void initialize_diffusion_coefficient_symmetric(double p_value) {
    // Create a vector of spheres
    std::vector<DiffusionCoefficient<dim>::Sphere> spheres;
    int n=2;
    // Define symmetrically distributed centers
    double step = 1.0 / (n+1);  // 'n' is the number of steps to divide each dimension
    double radius = 0.05;       // Fixed radius for each sphere

    // Generate symmetrically placed spheres
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= n; ++j) {
            for (int k = 1; k <= n; ++k) {
                if (spheres.size() < 10) {
                    Point<dim> center = {i * step, j * step, k * step};
                    spheres.push_back({center, radius});
                }
            }
        }
    }


    // Initialize the diffusion coefficient
    diffusion_coefficient = DiffusionCoefficient<dim>(spheres, p_value);
};

  void
  manage_flags(int argc, char ** argv);

  void
  write_csv(const long int elapsed_time, int iterations) const;

  static std::string
  extractFileName(const std::string& filePath);

  // Initialization.
  void
  setup();

  // System assembly.
  void
  assemble();

  // System solution.
  void
  solve();

  // Output.
  void
  output() const;

protected:
  // Path to the mesh file.
  std::string mesh_file_name = "../mesh/mesh-cube-20.msh";

  std::string preconditioner_name;

  bool symmetric;
  double p_value;
  // Polynomial degree.
  unsigned int r;

  // Forcing term.
  ForcingTerm forcing_term;

  DiffusionCoefficient<dim> diffusion_coefficient;


  //! Should we remove the members h, g?
  // h(x).
  FunctionH function_h;
  // g(x).
  FunctionG function_g;

  // Triangulation.
  Triangulation<dim> mesh;

  // Finite element space.
  // We use a unique_ptr here so that we can choose the type and degree of the
  // finite elements at runtime (the degree is a constructor parameter). The
  // class FiniteElement<dim> is an abstract class from which all types of
  // finite elements implemented by deal.ii inherit.
  std::unique_ptr<FiniteElement<dim>> fe;

  // Quadrature formula.
  // We use a unique_ptr here so that we can choose the type and order of the
  // quadrature formula at runtime (the order is a constructor parameter).
  std::unique_ptr<Quadrature<dim>> quadrature;

  // Quadrature formula used on boundary lines.
  std::unique_ptr<Quadrature<dim - 1>> quadrature_boundary;

  // DoF handler.
  DoFHandler<dim> dof_handler;

  // Sparsity pattern.
  SparsityPattern sparsity_pattern;

  // System matrix.
  SparseMatrix<double> system_matrix;

  // System right-hand side.
  Vector<double> system_rhs;

  // System solution.
  Vector<double> solution;
};

#endif