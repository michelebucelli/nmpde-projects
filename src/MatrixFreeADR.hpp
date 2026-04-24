#ifndef MATRIX_FREE_ADR
#define MATRIX_FREE_ADR

#include <deal.II/base/timer.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/conditional_ostream.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/fe/fe_q.h>
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

    MatrixFreeADR(const unsigned int N_,
                  const VectorArrayFunction<dim, double> &advection_coeff,
                  const ScalarArrayFunction<dim, double> &diffusion_coeff,
                  const ScalarArrayFunction<dim, double> &reaction_coeff,
                  const ScalarArrayFunction<dim, double> &forcing_term_);

    void set_dirichlet_bc(const types::boundary_id boundary_id,
                          const std::shared_ptr<Function<dim>> &bc_function);

    void set_neumann_bc(const types::boundary_id boundary_id,
                        const ScalarArrayFunction<dim, double> &bc_function);

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

    const ScalarArrayFunction<dim, double> forcing_term;

    std::map<types::boundary_id, std::shared_ptr<Function<dim>>>
    dirichlet_bcs;

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
    DistributedVector<double> solution_update;

    ConditionalOStream pcout;
};

template<int dim, int fe_degree>
MatrixFreeADR<dim, fe_degree>::MatrixFreeADR(
    const unsigned int N_,
    const VectorArrayFunction<dim, double> &advection_coeff,
    const ScalarArrayFunction<dim, double> &diffusion_coeff,
    const ScalarArrayFunction<dim, double> &reaction_coeff,
    const ScalarArrayFunction<dim, double> &forcing_term_
)
    : mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    , mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
    , N(N_)
    , forcing_term(forcing_term_)
    , mesh(MPI_COMM_WORLD)
    , fe(fe_degree)
    , system_matrix(advection_coeff, diffusion_coeff, reaction_coeff)
    , pcout(std::cout, mpi_rank == 0)
{
    //for (unsigned int b_id = 0; b_id < dim * 2; ++b_id) {
    //    dirichlet_bcs[b_id] =
    //        std::make_shared<Functions::ZeroFunction<dim>>();
    //}
}

template<int dim, int fe_degree>
void MatrixFreeADR<dim, fe_degree>::set_dirichlet_bc(
    const types::boundary_id boundary_id,
    const std::shared_ptr<Function<dim>> &bc_function)
{
    dirichlet_bcs[boundary_id] = bc_function;
    system_matrix.neumann_bcs.erase(boundary_id);
}

template<int dim, int fe_degree>
void MatrixFreeADR<dim, fe_degree>::set_neumann_bc(
    const types::boundary_id boundary_id,
    const ScalarArrayFunction<dim, double> &bc_function)
{
    system_matrix.neumann_bcs[boundary_id] = bc_function;
    dirichlet_bcs.erase(boundary_id);
}

template<int dim, int fe_degree>
void MatrixFreeADR<dim, fe_degree>::setup()
{
    {
        // TODO: Check matrix-free simplex mesh support
        GridGenerator::subdivided_hyper_cube(mesh, N, 0.0, 1.0, true);

        dof_handler.reinit(mesh);
        // TODO: Check DoFRenumbering::matrix_free_data_locality()
        dof_handler.distribute_dofs(fe);

        const IndexSet locally_relevant_dofs =
            DoFTools::extract_locally_relevant_dofs(dof_handler);

        constraints.reinit(dof_handler.locally_owned_dofs(),
                           locally_relevant_dofs);

        std::map<types::boundary_id, const Function<dim> *> bc_functions;
        for (const auto &entry : dirichlet_bcs) {
            bc_functions[entry.first] = entry.second.get();
        }

        VectorTools::interpolate_boundary_values(mapping,
                                                 dof_handler,
                                                 bc_functions,
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
            (update_values | update_gradients | update_JxW_values |
             update_quadrature_points);
        // TODO: Investigate these flags, if for example, update_normal_vector
        // is omitted, the code still converges
        additional_data.mapping_update_flags_boundary_faces =
            (update_values | update_JxW_values | update_quadrature_points |
             update_normal_vectors);

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
        system_matrix.initialize_dof_vector(solution_update);
    }
}

template<int dim, int fe_degree>
void MatrixFreeADR<dim, fe_degree>::assemble_rhs()
{
    // NOTE: There's another method to handle inhomogeneous dirichlet
    // bcs which should be faster for high-order elements with hanging-nodes.
    // See https://dealii.org/current/doxygen/deal.II/group__matrixfree.html

#ifdef LOG_TIMERS
    Timer time;
#endif

    solution = 0;
    constraints.distribute(solution);
    solution.update_ghost_values();
    system_rhs = 0;

    FEEvaluation<dim, fe_degree> phi(*system_matrix.get_matrix_free());
    FEFaceEvaluation<dim, fe_degree> phi_face(
        *system_matrix.get_matrix_free());

    for (unsigned int cell = 0;
         cell < system_matrix.get_matrix_free()->n_cell_batches(); ++cell) {

        phi.reinit(cell);
        // Ignores constraints on solution, reading all DoFs
        phi.read_dof_values_plain(solution);
        phi.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
        for (const unsigned int q : phi.quadrature_point_indices()) {
            const ArrayPoint<dim, double> q_p = phi.quadrature_point(q);
            const VectorizedArray<double> q_val = phi.get_value(q);

            phi.submit_gradient(
                system_matrix.advection_coeff(q_p) * q_val +
                -system_matrix.diffusion_coeff(q_p) * phi.get_gradient(q), q);

            phi.submit_value(
                forcing_term(q_p) +
                -system_matrix.reaction_coeff(q_p) * q_val, q);

        }
        phi.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
        phi.distribute_local_to_global(system_rhs);
    }

    const unsigned int n_inner_face_batches =
        system_matrix.get_matrix_free()->n_inner_face_batches();
    const unsigned int n_boundary_face_batches =
        system_matrix.get_matrix_free()->n_boundary_face_batches();

    for (unsigned int face = n_inner_face_batches;
         face < n_inner_face_batches + n_boundary_face_batches;
         ++face) {

        phi_face.reinit(face);

        /* Faces in the same batch are on the same boundary. */
        const types::boundary_id b_id = phi_face.boundary_id();

        if (system_matrix.neumann_bcs.count(b_id)) {
            const ScalarArrayFunction<dim, double> &bc_function =
                system_matrix.neumann_bcs[b_id];

            phi_face.gather_evaluate(solution, EvaluationFlags::values);
            for (const unsigned int q : phi_face.quadrature_point_indices()) {
                const ArrayPoint<dim, double> q_p = phi_face.quadrature_point(q);

                phi_face.submit_value(
                    -system_matrix.advection_coeff(q_p) *
                    phi_face.normal_vector(q) * phi_face.get_value(q) +
                    bc_function(q_p), q);
            }
            phi_face.integrate_scatter(EvaluationFlags::values, system_rhs);
        }
    }

    system_rhs.compress(VectorOperation::add);

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

    SolverControl solver_control(5000, 1e-12 * system_rhs.l2_norm());
    SolverGMRES<DistributedVector<double>> solver(solver_control);

    // TODO: Figure out how to use other preconditioners
    PreconditionIdentity preconditioner;

    // No constraints stuff here?
    solver.solve(system_matrix, solution_update, system_rhs, preconditioner);
    solution += solution_update;

    //pcout << solver_control.last_step() << " GMRES iterations" << std::endl;

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

    solution.update_ghost_values();
    data_out.add_data_vector(dof_handler, solution, "solution");
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
