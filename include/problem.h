#ifndef PROBLEM_H
#define PROBLEM_H

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <filesystem>
#include <fstream>
#include <iostream>

using namespace dealii;

/**
 * Class managing the differential problem.
 */
template<int dim, int fe_degree>
class Problem
{
public:

    // Diffusion coefficient.
    class MuCoefficient : public Function<dim>
    {
    public:
        // Constructor.
        MuCoefficient()
        {}
        
        // Evaluation.
        virtual double
        value(const Point<dim> & /*p*/,
              const unsigned int /*component*/ = 0) const override
        {
            return 0.1; // -∇·(µ∇u)
        }
    };

    // Advection coefficient.
    class BetaCoefficient : public Function<dim>
    {
    public:
        // Constructor.
        BetaCoefficient()
        {}

        // Evaluation.
        virtual double
        value(const Point<dim> & /*p*/,
              const unsigned int component = 0) const override
        {
            if(component == 0) // ∇·(βu)
                return 0.0;
            else
                return 0.0;
        }
    };

    // Reaction coefficient.
    class GammaCoefficient : public Function<dim>
    {
    public:
        // Constructor.
        GammaCoefficient()
        {}

        // Evaluation.
        virtual double
        value(const Point<dim> & /*p*/,
              const unsigned int /*component*/ = 0) const override
        {
            return 0.1; // γu
        }
    };

    // Forcing term.
    class SourceFunction : public Function<dim>
    {
    public:
        // Constructor.
        SourceFunction()
        {}

        // Evaluation.
        virtual double
        value(const Point<dim> & p,
            const unsigned int /*component*/ = 0) const override
        {
            return p[0] * p[1]; // f
        }
    };

    // Dirichlet boundary conditions.
    class GFunction : public Function<dim>
    {
    public:
        // Constructor.
        GFunction()
        {}

        // Evaluation.
        virtual double
        value(const Point<dim> & /*p*/,
            const unsigned int /*component*/ = 0) const override
        {
            return 0.1; // g
        }
    };

    // Neumann boundary conditions.
    class HFunction : public Function<dim>
    {
    public:
        // Constructor.
        HFunction()
        {}

        // Evaluation.
        virtual double
        value(const Point<dim> & /*p*/,
            const unsigned int /*component*/ = 0) const override
        {
            return 0.1; // h
        }
    };

    // Constructor.
    Problem()
        : mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
        , mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
        , mesh(MPI_COMM_WORLD)
        , pcout(std::cout, mpi_rank == 0)
    {}

    // Initialization.
    void
    setup(const unsigned int ref_level);

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

    // Number of MPI processes.
    const unsigned int mpi_size;

    // This MPI process.
    const unsigned int mpi_rank;

    // Diffusion coefficient.
    MuCoefficient mu_coefficient;

    // Advection coefficient.
    BetaCoefficient beta_coefficient;

    // Reaction coefficient.
    GammaCoefficient gamma_coefficient;

    // Forcing term.
    SourceFunction source_function;

    // g(x).
    GFunction g_function;

    // h(x).
    HFunction h_function;

    // Triangulation.
    parallel::distributed::Triangulation<dim> mesh;

    // Finite element space.
    std::unique_ptr<FiniteElement<dim>> fe;

    // Quadrature formula.
    std::unique_ptr<Quadrature<dim>> quadrature;

    // Quadrature formula used on boundary lines.
    std::unique_ptr<Quadrature<dim - 1>> quadrature_boundary;

    // DoF handler.
    DoFHandler<dim> dof_handler;

    // System matrix.
    TrilinosWrappers::SparseMatrix system_matrix;

    // System right-hand side.
    TrilinosWrappers::MPI::Vector system_rhs;

    // System solution.
    TrilinosWrappers::MPI::Vector solution;

    // Parallel output stream.
    ConditionalOStream pcout;

    // DoFs owned by current process.
    IndexSet locally_owned_dofs;
};

#endif