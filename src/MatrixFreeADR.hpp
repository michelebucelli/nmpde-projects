#ifndef MATRIX_FREE_ADR_HPP
#define MATRIX_FREE_ADR_HPP

#include <deal.II/base/timer.h>
#include <deal.II/base/config.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

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
    // @DISCUSSION: in the setup function the mesh is initialized with
    // a specific grid generator...
    //
    // As users of the benchmark we would like to swap at runtime, without
    // recompilation, the type of mesh, which could be either generated or
    // user supplied. The point was testing structured and unstructured meshes.
    //
    // So I would change this interface and add a mesh generator function

    void setup() override;
    void setup(bool dof_reorder);

    void assemble() override;

    // @DISCUSSION: as users of the benchmark, we need some experiments:
    //   1. measure the time per iteration MB vs. MF without any preconditioner
    //   2. sweep the preconditioning strategies
    //
    // so this makes sense to have control parameters as function arguments
    //
    // As Tommaso pointed out, rtol should be specialized on number, otherwise
    // on floating point runs, it won't converge. Let's make it a responsibility
    // of the caller of this class
    //
    // Ideally this interface should be replicated for the matrix based solver
    std::pair<unsigned int, double> solve(unsigned int max_iter,
                                          number rel_tol,
                                          bool apply_preconditioner) override;

protected:
    AffineConstraints<number> constraints;

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
    constraints.reinit(this->dof_handler.locally_owned_dofs(),
                       locally_relevant_dofs);
#else
    constraints.reinit(locally_relevant_dofs);
#endif

    BoundaryFuncMap<Function<dim>> tmp_dirichlet_map;
    std::map<types::boundary_id, const Function<dim> *> dirichlet_map;
    this->local_spec.fill_dirichlet_map(tmp_dirichlet_map);
    for (const auto& [boundary_id, dirichlet_function] : tmp_dirichlet_map) {
        dirichlet_map[boundary_id] = dirichlet_function.get();
    }

    VectorTools::interpolate_boundary_values(
        this->mapping, this->dof_handler, dirichlet_map, constraints);

    constraints.close();

    system_matrix.clear();

    typename MatrixFreeType::AdditionalData additional_data;
    additional_data.tasks_parallel_scheme =
        // WARNING: This should be set to none if MPI is used
        // to parallelize also WITHIN each node
        MatrixFreeType::AdditionalData::partition_partition;
    // Specifies what should be computed (and cached) on cells
    additional_data.mapping_update_flags =
        (update_values | update_gradients | update_JxW_values |
         update_quadrature_points);
    additional_data.mapping_update_flags_boundary_faces =
        (update_values | update_JxW_values | update_quadrature_points);

    if (dof_reorder) {
        DoFRenumbering::matrix_free_data_locality(
            this->dof_handler, constraints, additional_data);
        // NOTE: compute_matrix_free_data_locality()
        // returns the reordered indices only, could be insightful

        // Constraints need to be setup again after reordering

        locally_relevant_dofs =
        DoFTools::extract_locally_relevant_dofs(this->dof_handler);

#if DEAL_II_VERSION_GTE(9, 6, 0)
        constraints.reinit(this->dof_handler.locally_owned_dofs(),
                           locally_relevant_dofs);
#else
        constraints.reinit(locally_relevant_dofs);
#endif
        VectorTools::interpolate_boundary_values(
            this->mapping, this->dof_handler, dirichlet_map, constraints);

        constraints.close();
    }

    std::shared_ptr<MatrixFreeType> system_mf_storage(new MatrixFreeType());
    system_mf_storage->reinit(this->mapping,
                              this->dof_handler,
                              constraints,
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
    // NOTE: There's another method to handle inhomogeneous dirichlet
    // bcs which should be faster for high-order elements with hanging-nodes.
    // See https://dealii.org/current/doxygen/deal.II/group__matrixfree.html

#ifdef LOG_TIMERS
    Timer time;
#endif

    this->solution = 0;
    constraints.distribute(this->solution);
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

#ifdef LOG_TIMERS
    this->pcout << std::setw(20) << std::left << "assemble()"
                << std::fixed << time.wall_time() << "s" << std::endl;
#endif
}

template<unsigned int fe_degree, typename ADRSpec>
std::pair<unsigned int, double> MatrixFreeADR<fe_degree, ADRSpec>::solve(
    unsigned int max_iter, number rel_tol, bool apply_preconditioner)
{
#ifdef LOG_TIMERS
    Timer time;
#endif

    // TODO: Find a way to implement solve() in the base ADR class

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
        // TODO: Consider performing compute_diagonal() in the assembly
        system_matrix.compute_diagonal();

        PreconditionJacobi<ADROperator<fe_degree, ADRSpec>> preconditioner;
        preconditioner.initialize(system_matrix, 1.0);

        solver_time = run_solve(preconditioner);
    } else {
        solver_time = run_solve(PreconditionIdentity());
    }

    this->solution += solution_update;

    this->solution.update_ghost_values();

    //pcout << solver_control.last_step() << " GMRES iterations" << std::endl;

#ifdef LOG_TIMERS
    const double elapsed = time.wall_time();
    const double elapsed_per_iter = elapsed / solver_control->last_step();

    this->pcout << std::setw(20) << std::left << "solve() [tot/iter]"
                << std::fixed << elapsed << "s/" << elapsed_per_iter
                << "s" << std::endl;

#endif

    const auto res = solver_control->last_step();
    delete solver_control;

    return {res, solver_time};
}

#endif
