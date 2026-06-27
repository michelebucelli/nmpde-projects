#ifndef ADR_OPERATOR_HPP
#define ADR_OPERATOR_HPP

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/portable_matrix_free.h>
#include <deal.II/matrix_free/portable_fe_evaluation.h>

#include "../ADR.hpp"

using namespace dealii;

template<typename number>
using DistributedVector =
    LinearAlgebra::distributed::Vector<number, MemorySpace::Host>;

template<typename number>
using DistributedDeviceVector =
    LinearAlgebra::distributed::Vector<number, MemorySpace::Default>;

template<unsigned int fe_degree, typename ADRSpec>
class QuadADROperator
{
    using number = typename ADRSpec::ScalarType;

public:
    DEAL_II_HOST_DEVICE QuadADROperator() {}

    DEAL_II_HOST_DEVICE void operator()(
        Portable::FEEvaluation<
            ADRSpec::dim, fe_degree, fe_degree + 1, 1, number> *fe_eval,
        const int q) const;
};

template<unsigned int fe_degree, typename ADRSpec> DEAL_II_HOST_DEVICE
void QuadADROperator<fe_degree, ADRSpec>::operator()(
    Portable::FEEvaluation<
        ADRSpec::dim, fe_degree, fe_degree + 1, 1, number> *fe_eval,
    const int q) const
{
    const unsigned int cell_idx = fe_eval->get_current_cell_index();
    const typename Portable::MatrixFree<ADRSpec::dim, number>::Data *data =
        fe_eval->get_matrix_free_data();

    const auto q_p = data->get_quadrature_point(cell_idx, q);
    const auto q_val = fe_eval->get_value(q);
    const auto q_grad = fe_eval->get_gradient(q);

    const ADRSpec spec;

    fe_eval->submit_value(spec.reaction(q_p) * q_val +
                          spec.advection(q_p) * q_grad, q);

    fe_eval->submit_gradient(spec.diffusion(q_p) * q_grad, q);
}

template<unsigned int fe_degree, typename ADRSpec>
class QuadRhsEvaluator
{
    using number = typename ADRSpec::ScalarType;

public:
    DEAL_II_HOST_DEVICE QuadRhsEvaluator() {}

    DEAL_II_HOST_DEVICE void operator()(
        Portable::FEEvaluation<
            ADRSpec::dim, fe_degree, fe_degree + 1, 1, number> *fe_eval,
        const int q) const;

    // TODO: Remove this if not used
    static constexpr unsigned int n_q_points =
        Utilities::pow(fe_degree + 1, ADRSpec::dim);
};

template<unsigned int fe_degree, typename ADRSpec>
DEAL_II_HOST_DEVICE void
QuadRhsEvaluator<fe_degree, ADRSpec>::operator()(
    Portable::FEEvaluation<
        ADRSpec::dim, fe_degree, fe_degree + 1, 1, number> *fe_eval,
    const int q) const
{
    // TODO: How do I get the boundary id here? Could be used to build the
    // neumann contribution on the device

    const unsigned int cell_idx = fe_eval->get_current_cell_index();
    const typename Portable::MatrixFree<ADRSpec::dim, number>::Data *data =
        fe_eval->get_matrix_free_data();

    const auto q_p = data->get_quadrature_point(cell_idx, q);
    const auto q_val = fe_eval->get_value(q);
    const auto q_grad = fe_eval->get_gradient(q);

    // TODO: Perhaps store this object as class member
    // to reduce the number of instantiations? Check of parallel_for
    // schedules the loop
    const ADRSpec spec;

    fe_eval->submit_value(spec.forcing(q_p) -
                          spec.reaction(q_p) * q_val -
                          spec.advection(q_p) * q_grad, q);

    fe_eval->submit_gradient(-spec.diffusion(q_p) * q_grad, q);
}

template<unsigned int fe_degree, typename ADRSpec>
class LocalADROperator
{
    using number = typename ADRSpec::ScalarType;

public:
    LocalADROperator() {}

    DEAL_II_HOST_DEVICE void operator()(
        const typename Portable::MatrixFree<ADRSpec::dim, number>::Data *data,
        const Portable::DeviceVector<number> &src,
        Portable::DeviceVector<number> &dst) const;

    static constexpr unsigned int n_q_points =
        Utilities::pow(fe_degree + 1, ADRSpec::dim);

    static constexpr unsigned int n_local_dofs = n_q_points;
};

template<unsigned int fe_degree, typename ADRSpec>
DEAL_II_HOST_DEVICE void LocalADROperator<fe_degree, ADRSpec>::operator()(
    const typename Portable::MatrixFree<ADRSpec::dim, number>::Data *data,
    const Portable::DeviceVector<number> &src,
    Portable::DeviceVector<number> &dst) const
{
    Portable::FEEvaluation<
        ADRSpec::dim, fe_degree, fe_degree + 1, 1, number> fe_eval(data);

    fe_eval.read_dof_values(src);
    fe_eval.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
    fe_eval.apply_for_each_quad_point(QuadADROperator<fe_degree, ADRSpec>());
    fe_eval.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
    fe_eval.distribute_local_to_global(dst);
}

template<unsigned int fe_degree, typename ADRSpec>
class LocalRhsEvaluator
{
    using number = typename ADRSpec::ScalarType;

public:
    LocalRhsEvaluator() {}

    DEAL_II_HOST_DEVICE void operator()(
        const typename Portable::MatrixFree<ADRSpec::dim, number>::Data *data,
        const Portable::DeviceVector<number> &src,
        Portable::DeviceVector<number> &dst) const;

    static constexpr unsigned int n_q_points =
        Utilities::pow(fe_degree + 1, ADRSpec::dim);

    static constexpr unsigned int n_local_dofs = n_q_points;
};

template<unsigned int fe_degree, typename ADRSpec>
DEAL_II_HOST_DEVICE void
LocalRhsEvaluator<fe_degree, ADRSpec>::operator()(
    const typename Portable::MatrixFree<ADRSpec::dim, number>::Data *data,
    const Portable::DeviceVector<number> &src,
    Portable::DeviceVector<number> &dst) const
{
    Portable::FEEvaluation<
        ADRSpec::dim, fe_degree, fe_degree + 1, 1, number> fe_eval(data);

    // NOTE: I guess there's no need for read_dof_values_plain()
    // since constraints are not applied with the Portable version of
    // read_dof_values(). That's why copy_constrained_values() is used, right?
    fe_eval.read_dof_values(src);
    fe_eval.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
    fe_eval.apply_for_each_quad_point(QuadRhsEvaluator<fe_degree, ADRSpec>());
    fe_eval.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
    fe_eval.distribute_local_to_global(dst);
}

template<unsigned int fe_degree, typename ADRSpec>
class ADROperator : public EnableObserverPointer
{
    static constexpr int dim = ADRSpec::dim;
    using number = typename ADRSpec::ScalarType;

public:
    ADROperator(const DoFHandler<dim> &dof_handler,
                const MappingQ<dim> &mapping,
                const AffineConstraints<number> &constraints);

    void precompute_coefs(const DoFHandler<dim> &dof_handler);

    void vmult(DistributedDeviceVector<number> &dst,
               const DistributedDeviceVector<number> &src) const;

    void initialize_dof_vector(DistributedDeviceVector<number> &vec) const;

    void compute_rhs(DistributedDeviceVector<number> &dst,
                     const DistributedDeviceVector<number> &src) const;

    void compute_diagonal();

    std::shared_ptr<DiagonalMatrix<DistributedDeviceVector<number>>>
    get_matrix_diagonal_inverse() const;

    types::global_dof_index m() const;

    types::global_dof_index n() const;

    number el(const types::global_dof_index row,
              const types::global_dof_index col) const;

private:
    Portable::MatrixFree<dim, number> mf_data;

    std::shared_ptr<DiagonalMatrix<DistributedDeviceVector<number>>>
    inverse_diagonal_entries;
};

template<unsigned int fe_degree, typename ADRSpec>
ADROperator<fe_degree, ADRSpec>::ADROperator(
    const DoFHandler<dim> &dof_handler,
    const MappingQ<dim> &mapping,
    const AffineConstraints<number> &constraints)
{
    typename Portable::MatrixFree<dim, double>::AdditionalData
    additional_data;
    // Specifies what should be computed (and cached) on cells
    additional_data.mapping_update_flags =
        (update_values | update_gradients |
         update_JxW_values | update_quadrature_points);

    mf_data.reinit(mapping,
                   dof_handler,
                   constraints,
                   QGauss<1>(fe_degree + 1),
                   additional_data);
}

template<unsigned int fe_degree, typename ADRSpec>
void ADROperator<fe_degree, ADRSpec>::precompute_coefs(
    const DoFHandler<dim> &dof_handler)
{
    /*
    const unsigned int n_owned_cells =
        dynamic_cast<const parallel::TriangulationBase<dim> *>(
            &dof_handler.get_triangulation())->n_locally_owned_active_cells();

    forcing_coefs.reinit(n_owned_cells *
        LocalADROperator<dim, fe_degree, number>::n_q_points);

    const ForcingEvaluator<dim, fe_degree, number, ForcingPointFunctor>
    forcing_eval(forcing_coefs.get_values());

    mf_data.evaluate_coefficients(forcing_eval);
    */
}

template<unsigned int fe_degree, typename ADRSpec>
void ADROperator<fe_degree, ADRSpec>::vmult(
    DistributedDeviceVector<number> &dst,
    const DistributedDeviceVector<number> &src) const
{
    dst = 0;
    LocalADROperator<fe_degree, ADRSpec> loc_adr_operator;
    mf_data.cell_loop(loc_adr_operator, src, dst);
    // cell_loop is evaluated for all dofs, but for constrained
    // ones, we want the action of the operator to be the identity,
    // we thus copy the constrained src values into the dst
    mf_data.copy_constrained_values(src, dst);
}

template<unsigned int fe_degree, typename ADRSpec>
void ADROperator<fe_degree, ADRSpec>::initialize_dof_vector(
    DistributedDeviceVector<number> &vec) const
{
    mf_data.initialize_dof_vector(vec);
}

template<unsigned int fe_degree, typename ADRSpec>
void ADROperator<fe_degree, ADRSpec>::compute_rhs(
    DistributedDeviceVector<number> &dst,
    const DistributedDeviceVector<number> &src) const
{
    dst = 0;
    LocalRhsEvaluator<fe_degree, ADRSpec> loc_rhs_operator;
    mf_data.cell_loop(loc_rhs_operator, src, dst);
    // Making sure that rhs is zeroed out for constrained dofs
    // (we want solution_update = 0 there, since the initial solution
    // has already the inhomogeneous dirichlet values)
    mf_data.set_constrained_values(0, dst);
}

template<unsigned int fe_degree, typename ADRSpec>
void ADROperator<fe_degree, ADRSpec>::compute_diagonal()
{
    inverse_diagonal_entries.reset(
        new DiagonalMatrix<DistributedDeviceVector<number>>());
    DistributedDeviceVector<double> &inverse_diagonal =
        inverse_diagonal_entries->get_vector();
    // Resizes inverse_diagonal to the number of DoFs
    initialize_dof_vector(inverse_diagonal);

    // Applies local_compute_diagonal to unit vectors
    // to compute the operator diagonal
    MatrixFreeTools::
    compute_diagonal<dim, fe_degree, fe_degree + 1, 1, number>(
        mf_data,
        inverse_diagonal,
        QuadADROperator<fe_degree, ADRSpec>(),
        EvaluationFlags::values | EvaluationFlags::gradients,
        EvaluationFlags::values | EvaluationFlags::gradients);

    // TODO: Is this not needed anymore?
    // Set Dirichlet constrained DoFs to 1 to allow invertibility
    // this->set_constrained_entries_to_one(inverse_diagonal);

    number *raw_diagonal = inverse_diagonal.get_values();

    Kokkos::parallel_for(
        inverse_diagonal.locally_owned_size(), KOKKOS_LAMBDA(int i) {
            raw_diagonal[i] = 1. / raw_diagonal[i];
        });

    // Kokkos:fence() here?
}

template<unsigned int fe_degree, typename ADRSpec>
std::shared_ptr<DiagonalMatrix<DistributedDeviceVector<
    typename ADROperator<fe_degree, ADRSpec>::number>>>
ADROperator<fe_degree, ADRSpec>::get_matrix_diagonal_inverse() const
{
    return inverse_diagonal_entries;
}

template<unsigned int fe_degree, typename ADRSpec>
types::global_dof_index ADROperator<fe_degree, ADRSpec>::m() const
{
    // TODO: This returns the local size, right? I have to understand
    // how vmult is distributed across MPI processes
    return mf_data.get_vector_partitioner()->size();
}

template<unsigned int fe_degree, typename ADRSpec>
types::global_dof_index ADROperator<fe_degree, ADRSpec>::n() const
{
    return m();
}

template<unsigned int fe_degree, typename ADRSpec>
typename ADROperator<fe_degree, ADRSpec>::number
ADROperator<fe_degree, ADRSpec>::el(
    const types::global_dof_index row,
    [[maybe_unused]] const types::global_dof_index col) const
{
    return 1.0 / (*inverse_diagonal_entries)(row, row);
}

#endif
