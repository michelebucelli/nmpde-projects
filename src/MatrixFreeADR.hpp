#ifndef MATRIX_FREE_ADR_HPP
#define MATRIX_FREE_ADR_HPP

#include <deal.II/base/timer.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>
#include <functional>

#include "ADROperator.hpp"

using namespace dealii;

// Neumann BCs evaluate purely to scalars, so std::function is fine for boundary (RHS)
template<int dim>
using ScalarFunc = std::function<double(const Point<dim> &)>;

template<int dim, int fe_degree, 
         typename AdvectionFunc, typename DiffusionFunc, typename ReactionFunc, typename ForcingFunc>
class MatrixFreeADR
{
public:
    MatrixFreeADR(const unsigned int N_,
                  const AdvectionFunc &advection_coeff,
                  const DiffusionFunc &diffusion_coeff,
                  const ReactionFunc &reaction_coeff,
                  const ForcingFunc &forcing_term_);

    void set_dirichlet_bc(const types::boundary_id boundary_id,
                          const std::shared_ptr<Function<dim>> &bc_function);

    void set_neumann_bc(const types::boundary_id boundary_id,
                        const ScalarFunc<dim> &bc_function);

    void setup();
    void assemble_rhs();
    void solve();
    void output();

    double compute_error(const VectorTools::NormType &norm_type,
                         const Function<dim> &exact_solution) const;

private:
    const unsigned int mpi_rank;
    const unsigned int N;

    const ForcingFunc forcing_term;

    std::map<types::boundary_id, std::shared_ptr<Function<dim>>> dirichlet_bcs;
    std::map<types::boundary_id, ScalarFunc<dim>> neumann_bcs;

    parallel::distributed::Triangulation<dim> mesh;
    const FE_Q<dim> fe;
    DoFHandler<dim> dof_handler;
    const MappingQ1<dim> mapping;

    ADROperator<dim, fe_degree, double, AdvectionFunc, DiffusionFunc, ReactionFunc> system_matrix;
    AffineConstraints<double> constraints;

    DistributedVector<double> system_rhs;
    DistributedVector<double> solution;
    DistributedVector<double> solution_update;

    ConditionalOStream pcout;
};

template<int dim, int fe_degree, typename A, typename D, typename R, typename F>
MatrixFreeADR<dim, fe_degree, A, D, R, F>::MatrixFreeADR(
    const unsigned int N_, const A &adv, const D &diff, const R &react, const F &forc)
    : mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    , N(N_)
    , forcing_term(forc)
    , mesh(MPI_COMM_WORLD)
    , fe(fe_degree)
    , system_matrix(adv, diff, react)
    , pcout(std::cout, mpi_rank == 0)
{}

template<int dim, int fe_degree, typename A, typename D, typename R, typename F>
void MatrixFreeADR<dim, fe_degree, A, D, R, F>::set_dirichlet_bc(
    const types::boundary_id boundary_id, const std::shared_ptr<Function<dim>> &bc_function)
{
    dirichlet_bcs[boundary_id] = bc_function;
    neumann_bcs.erase(boundary_id);
}

template<int dim, int fe_degree, typename A, typename D, typename R, typename F>
void MatrixFreeADR<dim, fe_degree, A, D, R, F>::set_neumann_bc(
    const types::boundary_id boundary_id, const ScalarFunc<dim> &bc_function)
{
    neumann_bcs[boundary_id] = bc_function;
    dirichlet_bcs.erase(boundary_id);
}

template<int dim, int fe_degree, typename A, typename D, typename R, typename F>
void MatrixFreeADR<dim, fe_degree, A, D, R, F>::setup()
{
    GridGenerator::subdivided_hyper_cube(mesh, N, 0.0, 1.0, true);
    dof_handler.reinit(mesh);
    dof_handler.distribute_dofs(fe);
    
    // CRITICAL: Cache locality for matrix-free vector operations
    DoFRenumbering::Cuthill_McKee(dof_handler);

    const IndexSet locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);
    constraints.reinit(locally_relevant_dofs);

    std::map<types::boundary_id, const Function<dim> *> bc_functions;
    for (const auto &entry : dirichlet_bcs) {
        bc_functions[entry.first] = entry.second.get();
    }

    VectorTools::interpolate_boundary_values(mapping, dof_handler, bc_functions, constraints);
    constraints.close();

    system_matrix.clear();

    typename MatrixFree<dim, double>::AdditionalData additional_data;
    additional_data.tasks_parallel_scheme = MatrixFree<dim, double>::AdditionalData::none;
    additional_data.mapping_update_flags = (update_values | update_gradients | update_JxW_values | update_quadrature_points);
    // Boundary faces only need JxW and quadrature points for Neumann RHS (no normal vector needed now)
    additional_data.mapping_update_flags_boundary_faces = (update_values | update_JxW_values | update_quadrature_points);

    std::shared_ptr<MatrixFree<dim, double>> system_mf_storage(new MatrixFree<dim, double>());
    system_mf_storage->reinit(mapping, dof_handler, constraints, QGauss<1>(fe_degree + 1), additional_data);
    system_matrix.initialize(system_mf_storage);

    system_matrix.initialize_dof_vector(system_rhs);
    system_matrix.initialize_dof_vector(solution);
    system_matrix.initialize_dof_vector(solution_update);
}

template<int dim, int fe_degree, typename A, typename D, typename R, typename F>
void MatrixFreeADR<dim, fe_degree, A, D, R, F>::assemble_rhs()
{
    solution = 0;
    constraints.distribute(solution);
    solution.update_ghost_values();
    system_rhs = 0;

    FEEvaluation<dim, fe_degree> phi(*system_matrix.get_matrix_free());
    FEFaceEvaluation<dim, fe_degree> phi_face(*system_matrix.get_matrix_free());

    for (unsigned int cell = 0; cell < system_matrix.get_matrix_free()->n_cell_batches(); ++cell) {
        phi.reinit(cell);
        phi.read_dof_values_plain(solution);
        phi.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
        
        for (const unsigned int q : phi.quadrature_point_indices()) {
            const ArrayPoint<dim, double> q_p = phi.quadrature_point(q);
            const VectorizedArray<double> q_val = phi.get_value(q);

            // Compute residual: f - A(u)
            phi.submit_gradient(-system_matrix.diffusion_coeff(q_p) * phi.get_gradient(q), q);
            phi.submit_value(forcing_term(q_p) 
                             - system_matrix.reaction_coeff(q_p) * q_val 
                             - system_matrix.advection_coeff(q_p) * phi.get_gradient(q), q);
        }
        phi.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
        phi.distribute_local_to_global(system_rhs);
    }

    const unsigned int n_inner = system_matrix.get_matrix_free()->n_inner_face_batches();
    const unsigned int n_bound = system_matrix.get_matrix_free()->n_boundary_face_batches();

    for (unsigned int face = n_inner; face < n_inner + n_bound; ++face) {
        phi_face.reinit(face);
        const types::boundary_id b_id = phi_face.boundary_id();

        if (neumann_bcs.count(b_id)) {
            const auto &bc_function = neumann_bcs[b_id];
            
            // Just apply the standard Neumann RHS formulation: \int h * v ds
            for (const unsigned int q : phi_face.quadrature_point_indices()) {
                // Notice we only need standard Point<dim> for std::function evaluation
                // We extract element 0 because all vectorized array entries on the boundary 
                // face batch have the same geometry for Cartesian grids. 
                // For fully deformed meshes, loop over VectorizedArray elements.
                Point<dim> p;
                for(unsigned int d=0; d<dim; ++d) p[d] = phi_face.quadrature_point(q)[d][0];
                
                phi_face.submit_value(make_vectorized_array(bc_function(p)), q);
            }
            phi_face.integrate_scatter(EvaluationFlags::values, system_rhs);
        }
    }
    system_rhs.compress(VectorOperation::add);
}

template<int dim, int fe_degree, typename A, typename D, typename R, typename F>
void MatrixFreeADR<dim, fe_degree, A, D, R, F>::solve()
{
    // Ensure diagonal is ready for the preconditioner
    system_matrix.compute_diagonal();

    SolverControl solver_control(5000, 1e-12 * system_rhs.l2_norm());
    SolverGMRES<DistributedVector<double>> solver(solver_control);

    // Apply Jacobi preconditioning to massively speed up GMRES
    PreconditionJacobi<ADROperator<dim, fe_degree, double, A, D, R>> preconditioner;
    preconditioner.initialize(system_matrix, 1.0); 

    solver.solve(system_matrix, solution_update, system_rhs, preconditioner);
    constraints.distribute(solution_update);
    solution += solution_update;
}

template<int dim, int fe_degree, typename A, typename D, typename R, typename F>
void MatrixFreeADR<dim, fe_degree, A, D, R, F>::output()
{
    DataOut<dim> data_out;
    solution.update_ghost_values();
    data_out.add_data_vector(dof_handler, solution, "solution");
    data_out.build_patches(mapping);
    data_out.write_vtu_with_pvtu_record("./", "solution-" + std::to_string(N), 0, MPI_COMM_WORLD);
}

template<int dim, int fe_degree, typename A, typename D, typename R, typename F>
double MatrixFreeADR<dim, fe_degree, A, D, R, F>::compute_error(
    const VectorTools::NormType &norm_type, const Function<dim> &exact_solution) const
{
    Vector<double> error_per_cell(mesh.n_active_cells());
    VectorTools::integrate_difference(mapping, dof_handler, solution, exact_solution,
                                      error_per_cell, QGauss<dim>(fe_degree + 2), norm_type);
    return VectorTools::compute_global_error(mesh, error_per_cell, norm_type);
}

#endif
