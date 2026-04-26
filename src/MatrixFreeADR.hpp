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
#include <functional>

#include "ADROperator.hpp"

using namespace dealii;

template<int dim, int fe_degree, typename ForcingTerm>
class MatrixFreeADR
{
public:

    MatrixFreeADR(const unsigned int N_,
                  std::unique_ptr<Function<dim>> dirichlet_bc_);

    void setup();
    void assemble_rhs();
    void solve();
    void output();

    const ConditionalOStream &get_pcout() { return pcout; }

    double compute_error(const VectorTools::NormType &norm_type,
                         const Function<dim> &exact_solution) const;

private:
    const unsigned int mpi_rank;
    const unsigned int mpi_size;

    const unsigned int N;

    std::unique_ptr<Function<dim>> dirichlet_bc;

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
    DistributedDeviceVector<double> solution_update_dev;
    DistributedVector<double> ghost_solution_host;

    ConditionalOStream pcout;
};

template<int dim, int fe_degree, typename ForcingTerm>
MatrixFreeADR<dim, fe_degree, ForcingTerm>::MatrixFreeADR(
    const unsigned int N_,
    std::unique_ptr<Function<dim>> dirichlet_bc_
)
    : mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    , mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
    , N(N_)
    , dirichlet_bc(std::move(dirichlet_bc_))
    , mesh(MPI_COMM_WORLD)
    , fe(fe_degree)
    , pcout(std::cout, mpi_rank == 0)
{}

template<int dim, int fe_degree, typename ForcingTerm>
void MatrixFreeADR<dim, fe_degree, ForcingTerm>::setup()
{
#ifdef LOG_TIMERS
    Timer time;
#endif

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
                                             *dirichlet_bc,
                                             constraints);

    constraints.close();

    system_matrix_dev.reset(new ADROperator<dim, fe_degree, double>(
        dof_handler, mapping, constraints));

    //system_matrix_dev->template precompute_coefs<ForcingTerm>(dof_handler);

    ghost_solution_host.reinit(locally_owned_dofs,
                               locally_relevant_dofs,
                               MPI_COMM_WORLD);

    system_matrix_dev->initialize_dof_vector(solution_dev);
    solution_update_dev.reinit(solution_dev);
    system_rhs_dev.reinit(solution_dev);

#ifdef LOG_TIMERS
    pcout << std::setw(20) << std::left << "setup()"
          << std::fixed << time.wall_time() << "s" << std::endl;
#endif
}

template<int dim, int fe_degree, typename ForcingTerm>
void MatrixFreeADR<dim, fe_degree, ForcingTerm>::assemble_rhs()
{
#ifdef LOG_TIMERS
    Timer time;
#endif

    ghost_solution_host = 0;
    constraints.distribute(ghost_solution_host);
    ghost_solution_host.update_ghost_values();

    // TODO: Can you avoid this transfer?
    LinearAlgebra::ReadWriteVector<double> rw_vector(locally_owned_dofs);
    rw_vector.import_elements(ghost_solution_host, VectorOperation::insert);
    solution_dev.import_elements(rw_vector, VectorOperation::insert);

    system_matrix_dev->template compute_rhs<ForcingTerm>(system_rhs_dev,
                                                         solution_dev);

#ifdef LOG_TIMERS
    pcout << std::setw(20) << std::left << "assemble_rhs()"
          << std::fixed << time.wall_time() << "s" << std::endl;
#endif
}

template<int dim, int fe_degree, typename ForcingTerm>
void MatrixFreeADR<dim, fe_degree, ForcingTerm>::solve()
{
#ifdef LOG_TIMERS
    Timer time;
#endif

    // Not really needed currently
    system_matrix_dev->compute_diagonal();

    SolverControl solver_control(10000, 1e-12 * system_rhs_dev.l2_norm());
    SolverCG<DistributedDeviceVector<double>> solver(solver_control);

    // TODO: Figure out how to use other preconditioners
    PreconditionIdentity preconditioner;

    solver.solve(*system_matrix_dev,
                 solution_update_dev,
                 system_rhs_dev,
                 preconditioner);

    solution_dev += solution_update_dev;

    pcout << solver_control.last_step() << " GMRES iterations" << std::endl;

    // Copy solution to host
    // TODO: Time d2h copy
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

template<int dim, int fe_degree, typename ForcingTerm>
void MatrixFreeADR<dim, fe_degree, ForcingTerm>::output()
{
    DataOut<dim> data_out;

    data_out.add_data_vector(dof_handler, ghost_solution_host, "solution");
    data_out.build_patches(mapping);

    data_out.write_vtu_with_pvtu_record("./",
                                       "solution-" + std::to_string(N),
                                        0,
                                        MPI_COMM_WORLD);
}

template<int dim, int fe_degree, typename ForcingTerm>
double MatrixFreeADR<dim, fe_degree, ForcingTerm>::
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
