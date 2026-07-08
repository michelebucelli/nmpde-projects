#ifndef MATRIX_FREE_ADR_HPP
#define MATRIX_FREE_ADR_HPP

#include <deal.II/base/timer.h>
#include <deal.II/base/config.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/multithread_info.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>
#include <functional>

#include "ADR.hpp"
#include "ADRSpec.hpp"
#include "ADRFactory.hpp"
#include "ADROperator.hpp"

using namespace dealii;

template<unsigned int fe_degree, typename ADRSpec>
class MatrixFreeADR :
    public ADR<ADRSpec,
               typename ADRSpec::ScalarType::value_type,
               DistributedVector<typename ADRSpec::ScalarType::value_type>>
{
    using number = typename ADRSpec::ScalarType::value_type;
    // We need to explicitly bring into this scope the base template
    // class members due to the two-phase lookup (or use this->member)
    using ADR<ADRSpec, number, DistributedVector<number>>::dim;

    friend class ADRFactory;

    using ArrayScalar = typename ADRSpec::ScalarType;
    using ArrayTensor = typename ADRSpec::TensorType;
    using ArrayPoint = typename ADRSpec::PointType;
    using MatrixFreeType = MatrixFree<dim, number, ArrayScalar>;

    MatrixFreeADR(const parallel::distributed::Triangulation<dim> &mesh);

public:
    void setup() override;
    void setup(bool dof_reorder);

    void assemble() override;

    std::pair<unsigned int, double> solve(unsigned int max_iter,
                                          number rel_tol,
                                          bool apply_preconditioner) override;

    template <class Vector1_, class Vector2_>
    void vmult(Vector1_& dest, const Vector2_& src);

protected:
    ADROperator<fe_degree, ADRSpec> system_matrix;

    DistributedVector<number> solution_update;
};

template<unsigned int fe_degree, typename ADRSpec>
MatrixFreeADR<fe_degree, ADRSpec>::MatrixFreeADR(
    const parallel::distributed::Triangulation<dim> &mesh
)
    : ADR<ADRSpec, number, DistributedVector<number>>(mesh, fe_degree)
{}

template<unsigned int fe_degree, typename ADRSpec>
void MatrixFreeADR<fe_degree, ADRSpec>::setup()
{
    setup(true);
}

template<unsigned int fe_degree, typename ADRSpec>
void MatrixFreeADR<fe_degree, ADRSpec>::setup(bool dof_reorder)
{
    this->dof_handler.distribute_dofs(this->fe);

    IndexSet locally_relevant_dofs =
        DoFTools::extract_locally_relevant_dofs(this->dof_handler);

#if DEAL_II_VERSION_GTE(9, 6, 0)
    this->constraints.reinit(this->dof_handler.locally_owned_dofs(),
                             locally_relevant_dofs);
#else
    this->constraints.reinit(locally_relevant_dofs);
#endif

    DoFTools::make_hanging_node_constraints(this->dof_handler, this->constraints);

    BoundaryFuncMap<Function<dim>> tmp_dirichlet_map;
    std::map<types::boundary_id, const Function<dim> *> dirichlet_map;
    this->local_spec.fill_dirichlet_map(tmp_dirichlet_map);
    for (const auto& [boundary_id, dirichlet_function] : tmp_dirichlet_map) {
        dirichlet_map[boundary_id] = dirichlet_function.get();
    }

    VectorTools::interpolate_boundary_values(
        this->mapping, this->dof_handler, dirichlet_map, this->constraints);

    this->constraints.close();

    system_matrix.clear();

    typename MatrixFreeType::AdditionalData additional_data;
        // WARNING: This should be set to none if MPI is used
        // to parallelize also WITHIN each node
    additional_data.tasks_parallel_scheme =
        (MultithreadInfo::n_threads() > 1) ?
        MatrixFreeType::AdditionalData::partition_partition :
        MatrixFreeType::AdditionalData::none;
    // Specifies what should be computed (and cached) on cells
    additional_data.mapping_update_flags =
        (update_values | update_gradients | update_JxW_values |
         update_quadrature_points);
    additional_data.mapping_update_flags_boundary_faces =
        (update_values | update_JxW_values | update_quadrature_points);

    if (dof_reorder) {
        DoFRenumbering::matrix_free_data_locality(
            this->dof_handler, this->constraints, additional_data);
        // NOTE: compute_matrix_free_data_locality()
        // returns the reordered indices only, could be insightful

        // Constraints need to be setup again after reordering

        locally_relevant_dofs =
        DoFTools::extract_locally_relevant_dofs(this->dof_handler);

#if DEAL_II_VERSION_GTE(9, 6, 0)
        this->constraints.reinit(this->dof_handler.locally_owned_dofs(),
                                 locally_relevant_dofs);
#else
        this->constraints.reinit(locally_relevant_dofs);
#endif
        DoFTools::make_hanging_node_constraints(this->dof_handler, this->constraints);

        VectorTools::interpolate_boundary_values(
            this->mapping, this->dof_handler, dirichlet_map, this->constraints);

        this->constraints.close();
    }

    std::shared_ptr<MatrixFreeType> system_mf_storage(new MatrixFreeType());
    system_mf_storage->reinit(this->mapping,
                              this->dof_handler,
                              this->constraints,
                              // dim=1 since basis functions are decomposed
                              // as tensor products of 1D functions
                              QGauss<1>(fe_degree + 1),
                              additional_data);

    system_matrix.initialize(system_mf_storage);

    system_matrix.initialize_dof_vector(this->system_rhs);
    system_matrix.initialize_dof_vector(this->solution);
    system_matrix.initialize_dof_vector(solution_update);
}

template<unsigned int fe_degree, typename ADRSpec>
void MatrixFreeADR<fe_degree, ADRSpec>::assemble()
{
    this->solution = 0;
    this->constraints.distribute(this->solution);
    this->solution.update_ghost_values();
    this->system_rhs = 0;

    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number, ArrayScalar>
    phi(*system_matrix.get_matrix_free());

    FEFaceEvaluation<dim, fe_degree, fe_degree + 1, 1, number, ArrayScalar>
    phi_face(*system_matrix.get_matrix_free());

    for (unsigned int cell = 0;
         cell < system_matrix.get_matrix_free()->n_cell_batches(); ++cell) {

        phi.reinit(cell);
        // Ignores constraints on solution, reading all DoFs
        phi.read_dof_values_plain(this->solution);
        phi.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);

        for (const unsigned int q : phi.quadrature_point_indices()) {
            const ArrayPoint q_p = phi.quadrature_point(q);

            const ArrayScalar q_val = phi.get_value(q);
            const ArrayTensor q_grad = phi.get_gradient(q);

            phi.submit_value(this->local_spec.forcing(q_p) -
                             this->local_spec.reaction(q_p) * q_val -
                             this->local_spec.advection(q_p) * q_grad, q);

            phi.submit_gradient(-this->local_spec.diffusion(q_p) * q_grad, q);
        }

        phi.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
        phi.distribute_local_to_global(this->system_rhs);
    }

    BoundaryFuncMap<NeumannFunction<ArrayScalar, ArrayPoint>> neumann_map;
    this->local_spec.fill_neumann_map(neumann_map);

    const unsigned int n_inner =
        system_matrix.get_matrix_free()->n_inner_face_batches();
    const unsigned int n_bound =
        system_matrix.get_matrix_free()->n_boundary_face_batches();

    for (unsigned int face = n_inner; face < n_inner + n_bound; ++face) {
        phi_face.reinit(face);
        const types::boundary_id b_id = phi_face.boundary_id();

        if (neumann_map.count(b_id)) {
            const auto &bc_function = neumann_map[b_id];

            for (const unsigned int q : phi_face.quadrature_point_indices()) {
                const ArrayPoint q_p = phi_face.quadrature_point(q);

                const ArrayScalar mu_val = this->local_spec.diffusion(q_p);
                const ArrayScalar h_val  = bc_function->value(q_p);

                phi_face.submit_value(mu_val * h_val, q);
            }
            phi_face.integrate_scatter(EvaluationFlags::values, this->system_rhs);
        }
    }
    this->system_rhs.compress(VectorOperation::add);
}

template<unsigned int fe_degree, typename ADRSpec>
std::pair<unsigned int, double> MatrixFreeADR<fe_degree, ADRSpec>::solve(
    unsigned int max_iter, number rel_tol, bool apply_preconditioner)
{
    // we can reuse somewhere this piece of code
    // static constexpr number tolerance_factor =
    //     std::is_same_v<number, float> ? 1e-6f : 1e-12;
    SolverControl* solver_control;

    if (rel_tol == 0.0) {
        // same as in MatrixBasedADR::solve
        solver_control = new IterationNumberControl(max_iter, 0.0);
    } else {
        solver_control = new SolverControl(
            max_iter, rel_tol * this->system_rhs.l2_norm());
    }

    SolverGMRES<DistributedVector<number>> solver(*solver_control);

    auto run_solve = [&](const auto &preconditioner) {
        Timer solver_timer;
        solver.solve(system_matrix,
                     solution_update,
                     this->system_rhs,
                     preconditioner);
        return solver_timer.wall_time();
    };

    double solver_time = 0.0;
    // this should later be a switch based on an
    // ENUM value of supported preconditioners
    if (apply_preconditioner) {
        // Ensure diagonal is ready for the preconditioner
        system_matrix.compute_diagonal();

        PreconditionJacobi<ADROperator<fe_degree, ADRSpec>> preconditioner;
        preconditioner.initialize(system_matrix, 1.0);

        solver_time = run_solve(preconditioner);
    } else {
        solver_time = run_solve(PreconditionIdentity());
    }

    this->solution += solution_update;

    this->solution.update_ghost_values();

    const auto res = solver_control->last_step();
    delete solver_control;

    return {res, solver_time};
}

template<unsigned int fe_degree, typename ADRSpec>
template <class Vector1_, class Vector2_>
void MatrixFreeADR<fe_degree, ADRSpec>::vmult(
    Vector1_& dest, const Vector2_& src) {
    system_matrix.vmult(dest, src);
}

#endif
