#ifndef MATRIX_FREE_ADR
#define MATRIX_FREE_ADR

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/conditional_ostream.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>
#include <fstream>

#include "ADROperator.hpp"

using namespace dealii;

template<int dim, int fe_degree>
class MatrixFreeADR
{
public:
    MatrixFreeADR(const unsigned int N_,
                  const VectorizedFunction<dim, double> &diffusion_coeff_,
                  const VectorizedFunction<dim, double> &forcing_term_);

    void setup();
    void assemble_rhs();
    void solve();
    void output();

    double compute_error(const VectorTools::NormType &norm_type,
                         const Function<dim> &exact_solution) const;

private:
    const unsigned int mpi_rank;
    const unsigned int mpi_size;

    // Number of elements for each dimension
    const unsigned int N;

    // Diffusion coefficient
    const VectorizedFunction<dim, double> diffusion_coeff;
    // Forcing term
    const VectorizedFunction<dim, double> forcing_term;

    // TODO: Can I use fullydistributed?
    parallel::distributed::Triangulation<dim> mesh;

    const FE_Q<dim> fe;
    DoFHandler<dim> dof_handler;

    // TODO: What's the purpose of this?
    const MappingQ1<dim> mapping;

    ADROperator<dim, fe_degree, double> system_matrix;

    AffineConstraints<double> constraints;

    DistributedVector<double> system_rhs;
    DistributedVector<double> solution;

    ConditionalOStream pcout;
};

template<int dim, int fe_degree>
MatrixFreeADR<dim, fe_degree>::MatrixFreeADR(
    const unsigned int N_,
    const VectorizedFunction<dim, double> &diffusion_coeff_,
    const VectorizedFunction<dim, double> &forcing_term_
)
    : mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    , mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
    , N(N_)
    , diffusion_coeff(diffusion_coeff_)
    , forcing_term(forcing_term_)
    , mesh(MPI_COMM_WORLD)
    , fe(fe_degree)
    , system_matrix(diffusion_coeff_)
    , pcout(std::cout, mpi_rank == 0)
{}

template<int dim, int fe_degree>
void MatrixFreeADR<dim, fe_degree>::setup()
{
    {
        // TODO: Check matrix-free simplex mesh support
        GridGenerator::subdivided_hyper_cube(mesh, N, 0.0, 1.0, false);

        dof_handler.reinit(mesh);
        // TODO: Check DoFRenumbering::matrix_free_data_locality()
        dof_handler.distribute_dofs(fe);

        const IndexSet locally_relevant_dofs =
            DoFTools::extract_locally_relevant_dofs(dof_handler);

        constraints.reinit(locally_relevant_dofs);
        VectorTools::interpolate_boundary_values(mapping,
                                                 dof_handler,
                                                 0,
                                                 Functions::ZeroFunction<dim>(),
                                                 constraints);
        constraints.close();
    }

    {
        system_matrix.clear();

        typename MatrixFree<dim, double>::AdditionalData additional_data;
        // TODO: Allow threading within each MPI process
        additional_data.tasks_parallel_scheme =
            MatrixFree<dim, double>::AdditionalData::none;
        // Specifies what should be computed (and cached) on cells
        additional_data.mapping_update_flags =
            (update_gradients | update_JxW_values | update_quadrature_points);

        std::shared_ptr<MatrixFree<dim, double>> system_mf_storage(
            new MatrixFree<dim, double>());
        system_mf_storage->reinit(mapping,
                                  dof_handler,
                                  constraints,
                                  // dim=1 since basis functions are decomposed
                                  // as tensor products of 1D functions
                                  QGauss<1>(fe_degree + 1),
                                  additional_data);
        system_matrix.initialize(system_mf_storage);

        system_matrix.initialize_dof_vector(system_rhs);
        system_matrix.initialize_dof_vector(solution);
    }
}

template<int dim, int fe_degree>
void MatrixFreeADR<dim, fe_degree>::assemble_rhs()
{
    system_rhs = 0;

    FEEvaluation<dim, fe_degree> phi(*system_matrix.get_matrix_free());
    for (unsigned int cell = 0;
         cell < system_matrix.get_matrix_free()->n_cell_batches(); ++cell) {

        phi.reinit(cell);
        for (const unsigned int q : phi.quadrature_point_indices()) {
            phi.submit_value(forcing_term(phi.quadrature_point(q)), q);
        }
        phi.integrate(EvaluationFlags::values);
        phi.distribute_local_to_global(system_rhs);
    }
    system_rhs.compress(VectorOperation::add);
}

template<int dim, int fe_degree>
void MatrixFreeADR<dim, fe_degree>::solve()
{
    SolverControl solver_control(1000, 1e-12 * system_rhs.l2_norm());
    SolverCG<DistributedVector<double>> solver(solver_control);

    PreconditionIdentity preconditioner;

    constraints.set_zero(solution); // Is this needed?
    solver.solve(system_matrix, solution, system_rhs, preconditioner);
    // Solve for the constrained DoFs
    constraints.distribute(solution);

    pcout << solver_control.last_step() << " CG iterations" << std::endl;
}

template<int dim, int fe_degree>
void MatrixFreeADR<dim, fe_degree>::output()
{
    DataOut<dim> data_out;

    solution.update_ghost_values();
    data_out.add_data_vector(dof_handler, solution, "solution");
    data_out.build_patches(mapping);

    const std::string output_file_name =
        "solution-" + std::to_string(N) + ".vtk";
    std::ofstream output_file(output_file_name);
    data_out.write_vtk(output_file);

    pcout << "Output saved to " << output_file_name << std::endl;
}

template<int dim, int fe_degree>
double MatrixFreeADR<dim, fe_degree>::
    compute_error(const VectorTools::NormType &norm_type,
                  const Function<dim> &exact_solution) const
{
    Vector<double> error_per_cell(mesh.n_active_cells());

    VectorTools::integrate_difference(mapping,
                                      dof_handler,
                                      solution,
                                      exact_solution,
                                      error_per_cell,
                                      QGauss<dim>(fe_degree + 2), 
                                      norm_type);

    const double error =
        VectorTools::compute_global_error(mesh, error_per_cell, norm_type);

    return error;
}

#endif
