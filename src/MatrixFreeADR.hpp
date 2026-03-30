#ifndef MATRIX_FREE_ADR
#define MATRIX_FREE_ADR

#include <deal.II/base/timer.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/conditional_ostream.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/lac/solver_gmres.h>
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

    MatrixFreeADR(const unsigned int N_);

    void setup();
    void assemble_rhs();
    void solve();
    void output();

    double compute_error(const VectorTools::NormType &norm_type,
                         const Function<dim> &exact_solution) const;

private:
    const unsigned int mpi_rank;
    const unsigned int mpi_size;

    const unsigned int N;

    // TODO: Can I use fullydistributed?
    parallel::distributed::Triangulation<dim> mesh;

    const FE_Q<dim> fe;
    DoFHandler<dim> dof_handler;

    // TODO: What's the purpose of this?
    const MappingQ1<dim> mapping;

    std::unique_ptr<ADROperator<dim, fe_degree, double>> system_matrix_dev;

    AffineConstraints<double> constraints;

    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;

    // TODO: Make sure DEAL_II_MPI_WITH_DEVICE_SUPPORT is defined
    DistributedDeviceVector<double> system_rhs_dev;
    DistributedDeviceVector<double> solution_dev;
    DistributedVector<double> ghost_solution_host;

    ConditionalOStream pcout;
};

template<int dim, int fe_degree>
MatrixFreeADR<dim, fe_degree>::MatrixFreeADR(const unsigned int N_)
    : mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    , mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
    , N(N_)
    , mesh(MPI_COMM_WORLD)
    , fe(fe_degree)
    , pcout(std::cout, mpi_rank == 0)
{}

template<int dim, int fe_degree>
void MatrixFreeADR<dim, fe_degree>::setup()
{
    // TODO: Check matrix-free simplex mesh support
    GridGenerator::subdivided_hyper_cube(mesh, N, 0.0, 1.0, false);

    dof_handler.reinit(mesh);
    // TODO: Check DoFRenumbering::matrix_free_data_locality()
    dof_handler.distribute_dofs(fe);

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    locally_relevant_dofs =
        DoFTools::extract_locally_relevant_dofs(dof_handler);

    constraints.clear();
    constraints.reinit(locally_owned_dofs, locally_relevant_dofs);

    VectorTools::interpolate_boundary_values(mapping,
                                             dof_handler,
                                             0,
                                             Functions::ZeroFunction<dim>(),
                                             constraints);

    constraints.close();

    system_matrix_dev.reset(new ADROperator<dim, fe_degree, double>(
        dof_handler, mapping, constraints));

    ghost_solution_host.reinit(locally_owned_dofs,
                               locally_relevant_dofs,
                               MPI_COMM_WORLD);

    system_matrix_dev->initialize_dof_vector(solution_dev);
    system_rhs_dev.reinit(solution_dev);
}

template<int dim, int fe_degree>
void MatrixFreeADR<dim, fe_degree>::assemble_rhs()
{
#ifdef LOG_TIMERS
    Timer time;
#endif

    // TODO: Assemble rhs directly on device using cell_loop,
    // since Portable::FEEvaluation::read_dof_values_plain()
    // is missing, use two set of constraints to handle
    // inhomogeneous dirichlet boundary conditions

    DistributedVector<double> system_rhs_host(locally_owned_dofs,
                                              locally_relevant_dofs,
                                              MPI_COMM_WORLD);

    const QGauss<dim> quadrature(fe_degree + 1);

    FEValues<dim> fe_values(fe,
                            quadrature,
                            update_values | update_JxW_values |
                            update_quadrature_points);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q = quadrature.size();

    Vector<double> cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators()) {
        if (cell->is_locally_owned()) {

            cell_rhs = 0;
            fe_values.reinit(cell);

            for (unsigned int q = 0; q < n_q; ++q) {
                const double f_loc = 1.0;

                for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                    cell_rhs(i) += f_loc *
                                   fe_values.shape_value(i, q) *
                                   fe_values.JxW(q);
                }
            }

            cell->get_dof_indices(local_dof_indices);
            constraints.distribute_local_to_global(cell_rhs,
                                                   local_dof_indices,
                                                   system_rhs_host);
        }
    }

    system_rhs_host.compress(VectorOperation::add);

    // Copy rhs to device memory
    LinearAlgebra::ReadWriteVector<double> rw_vector(locally_owned_dofs);
    rw_vector.import_elements(system_rhs_host, VectorOperation::insert);
    system_rhs_dev.import_elements(rw_vector, VectorOperation::insert);

#ifdef LOG_TIMERS
    pcout << std::setw(20) << std::left << "assemble_rhs()"
          << std::fixed << time.wall_time() << "s" << std::endl;
#endif
}

template<int dim, int fe_degree>
void MatrixFreeADR<dim, fe_degree>::solve()
{
#ifdef LOG_TIMERS
    Timer time;
#endif

    // Not really needed currently
    system_matrix_dev->compute_diagonal();

    SolverControl solver_control(5000, 1e-12 * system_rhs_dev.l2_norm());
    SolverGMRES<DistributedDeviceVector<double>> solver(solver_control);

    // TODO: Figure out how to use other preconditioners
    PreconditionIdentity preconditioner;

    solver.solve(
        *system_matrix_dev, solution_dev, system_rhs_dev, preconditioner);

    pcout << solver_control.last_step() << " GMRES iterations" << std::endl;

    // Copy solution to host
    LinearAlgebra::ReadWriteVector<double> rw_vector(locally_owned_dofs);
    rw_vector.import_elements(solution_dev, VectorOperation::insert);
    ghost_solution_host.import_elements(rw_vector, VectorOperation::insert);

    constraints.distribute(ghost_solution_host);

    ghost_solution_host.update_ghost_values();

#ifdef LOG_TIMERS
    const double elapsed = time.wall_time();
    const double elapsed_per_iter = elapsed / solver_control.last_step();

    pcout << std::setw(20) << std::left << "solve() [tot/iter]"
          << std::fixed << elapsed << "s/" << elapsed_per_iter
          << "s" << std::endl;
#endif
}

template<int dim, int fe_degree>
void MatrixFreeADR<dim, fe_degree>::output()
{
    DataOut<dim> data_out;

    data_out.add_data_vector(dof_handler, ghost_solution_host, "solution");
    data_out.build_patches(mapping);

    data_out.write_vtu_with_pvtu_record("./",
                                       "solution-" + std::to_string(N),
                                        0,
                                        MPI_COMM_WORLD);
}

template<int dim, int fe_degree>
double MatrixFreeADR<dim, fe_degree>::
    compute_error(const VectorTools::NormType &norm_type,
                  const Function<dim> &exact_solution) const
{
    Vector<double> error_per_cell(mesh.n_active_cells());

    VectorTools::integrate_difference(mapping,
                                      dof_handler,
                                      ghost_solution_host,
                                      exact_solution,
                                      error_per_cell,
                                      QGauss<dim>(fe_degree + 2),
                                      norm_type);

    const double error =
        VectorTools::compute_global_error(mesh, error_per_cell, norm_type);

    return error;
}

#endif
