#ifndef ADR_OPERATOR_HPP
#define ADR_OPERATOR_HPP

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/portable_matrix_free.h>
#include <deal.II/matrix_free/portable_fe_evaluation.h>

using namespace dealii;

template<typename number>
using DistributedVector =
    LinearAlgebra::distributed::Vector<number, MemorySpace::Host>;

template<typename number>
using DistributedDeviceVector =
    LinearAlgebra::distributed::Vector<number, MemorySpace::Default>;

template<int dim, int fe_degree, typename number, typename PointFunctor>
class ForcingEvaluator
{
    // TODO: Try to find a way to compute the forcing
    // coefficient on the fly (while computing the rhs)
    // and compare the performances
public:
    ForcingEvaluator(number *coefs_) :
        coefs(coefs_)
    {}

    DEAL_II_HOST_DEVICE void operator()(
        const typename Portable::MatrixFree<dim, number>::Data *data,
        const unsigned int cell,
        const unsigned int q) const;

    static constexpr unsigned int n_q_points =
        Utilities::pow(fe_degree + 1, dim);

private:
    number *coefs;
};

template<int dim, int fe_degree, typename number, typename PointFunctor>
DEAL_II_HOST_DEVICE void
ForcingEvaluator<dim, fe_degree, number, PointFunctor>::operator()(
    const typename Portable::MatrixFree<dim, number>::Data *data,
    const unsigned int cell,
    const unsigned int q) const
{
    const unsigned int q_idx = data->local_q_point_id(cell, n_q_points, q);
    const Point<dim> q_point = data->get_quadrature_point(cell, q);

    const PointFunctor forcing_functor;
    // forcing_functor call can be inlined
    coefs[q_idx] = forcing_functor(q_point);
}

template<int dim, int fe_degree, typename number>
class ADROperatorQuad
{
public:
    DEAL_II_HOST_DEVICE ADROperatorQuad() {}

    DEAL_II_HOST_DEVICE void operator()(
        Portable::FEEvaluation<
            dim, fe_degree, fe_degree + 1, 1, number> *fe_eval,
        const int q_point) const;
};

template<int dim, int fe_degree, typename number> DEAL_II_HOST_DEVICE
void ADROperatorQuad<dim, fe_degree, number>::operator()(
    Portable::FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> *fe_eval,
    const int q_point) const
{
    auto q_grad = fe_eval->get_gradient(q_point);
    //fe_eval->submit_value(q_grad + fe_eval->get_value(q_point), q_point);
    fe_eval->submit_gradient(q_grad, q_point);
}

template<int dim, int fe_degree, typename number>
class QuadRhsEvaluator
{
public:
    DEAL_II_HOST_DEVICE QuadRhsEvaluator(number *forcing_) :
        forcing(forcing_)
    {}

    DEAL_II_HOST_DEVICE void operator()(
        Portable::FEEvaluation<
            dim, fe_degree, fe_degree + 1, 1, number> *fe_eval,
        const int q_point) const;

    static constexpr unsigned int n_q_points =
        Utilities::pow(fe_degree + 1, dim);

private:
    number *forcing;
};

template<int dim, int fe_degree, typename number> DEAL_II_HOST_DEVICE
void QuadRhsEvaluator<dim, fe_degree, number>::operator()(
    Portable::FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> *fe_eval,
    const int q_point) const
{
    const unsigned int cell_idx = fe_eval->get_current_cell_index();
    const typename Portable::MatrixFree<dim, number>::Data *data =
        fe_eval->get_matrix_free_data();

    const unsigned int q_idx =
        data->local_q_point_id(cell_idx, n_q_points, q_point);

    fe_eval->submit_value(forcing[q_idx], q_point);
}

template<int dim, int fe_degree, typename number>
class LocalADROperator
{
public:
    LocalADROperator() {}

    DEAL_II_HOST_DEVICE void operator()(
        const typename Portable::MatrixFree<dim, number>::Data *data,
        const Portable::DeviceVector<number> &src,
        Portable::DeviceVector<number> &dst) const;

    static constexpr unsigned int n_q_points =
        Utilities::pow(fe_degree + 1, dim);

    static constexpr unsigned int n_local_dofs = n_q_points;
};

template<int dim, int fe_degree, typename number>
DEAL_II_HOST_DEVICE void LocalADROperator<dim, fe_degree, number>::operator()(
    const typename Portable::MatrixFree<dim, number>::Data *data,
    const Portable::DeviceVector<number> &src,
    Portable::DeviceVector<number> &dst) const
{
    Portable::FEEvaluation<
        dim, fe_degree, fe_degree + 1, 1, number> fe_eval(data);

    fe_eval.read_dof_values(src);
    fe_eval.evaluate(EvaluationFlags::gradients);
    fe_eval.apply_for_each_quad_point(
        ADROperatorQuad<dim, fe_degree, number>());
    fe_eval.integrate(EvaluationFlags::gradients);
    fe_eval.distribute_local_to_global(dst);
}

template<int dim, int fe_degree, typename number>
class LocalRhsEvaluator
{
public:
    LocalRhsEvaluator(number *forcing_) :
        forcing(forcing_)
    {}

    DEAL_II_HOST_DEVICE void operator()(
        const typename Portable::MatrixFree<dim, number>::Data *data,
        const Portable::DeviceVector<number> &src,
        Portable::DeviceVector<number> &dst) const;

    static constexpr unsigned int n_q_points =
        Utilities::pow(fe_degree + 1, dim);

    static constexpr unsigned int n_local_dofs = n_q_points;

private:
    number *forcing;
};

template<int dim, int fe_degree, typename number>
DEAL_II_HOST_DEVICE void LocalRhsEvaluator<dim, fe_degree, number>::operator()(
    const typename Portable::MatrixFree<dim, number>::Data *data,
    const Portable::DeviceVector<number> &src,
    Portable::DeviceVector<number> &dst) const
{
    Portable::FEEvaluation<
        dim, fe_degree, fe_degree + 1, 1, number> fe_eval(data);

    fe_eval.apply_for_each_quad_point(
        QuadRhsEvaluator<dim, fe_degree, number>(forcing));
    fe_eval.integrate(EvaluationFlags::values);
    // FIX: I getting slightly worse convergence rates when
    // computing the rhs on the device (only for fe_degree=1)
    // Could distribute_local_to_global be the issue? This
    // call is supposed to also apply constraints, which is
    // not something we care about when building the rhs.
    // Perhaps building the rhs with cell_loop is not the
    // best option? I could try using evaluate_coefficients
    // or two sets of constraints (which is anyway needed for
    // inhomogeneous constraints)
    fe_eval.distribute_local_to_global(dst);
}

template<int dim, int fe_degree, typename number>
class ADROperator : public EnableObserverPointer
{
public:
    ADROperator(const DoFHandler<dim> &dof_handler,
                const MappingQ<dim> &mapping,
                const AffineConstraints<number> &constraints);

    template<typename ForcingPointFunctor>
    void precompute_coefs(const DoFHandler<dim> &dof_handler);

    void vmult(DistributedDeviceVector<number> &dst,
               const DistributedDeviceVector<number> &src) const;

    void initialize_dof_vector(DistributedDeviceVector<number> &vec) const;

    void compute_rhs(DistributedDeviceVector<number> &dst) const;

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

    DistributedDeviceVector<number> forcing_coefs;
};

template<int dim, int fe_degree, typename number>
ADROperator<dim, fe_degree, number>::ADROperator(
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

template<int dim, int fe_degree, typename number>
template<typename ForcingPointFunctor>
void ADROperator<dim, fe_degree, number>::precompute_coefs(
    const DoFHandler<dim> &dof_handler)
{
    const unsigned int n_owned_cells =
        dynamic_cast<const parallel::TriangulationBase<dim> *>(
            &dof_handler.get_triangulation())->n_locally_owned_active_cells();

    forcing_coefs.reinit(n_owned_cells *
        LocalADROperator<dim, fe_degree, number>::n_q_points);

    const ForcingEvaluator<dim, fe_degree, number, ForcingPointFunctor>
    forcing_eval(forcing_coefs.get_values());

    mf_data.evaluate_coefficients(forcing_eval);
}

template<int dim, int fe_degree, typename number>
void ADROperator<dim, fe_degree, number>::vmult(
    DistributedDeviceVector<number> &dst,
    const DistributedDeviceVector<number> &src) const
{
    dst = 0;
    LocalADROperator<dim, fe_degree, number> loc_adr_operator;
    mf_data.cell_loop(loc_adr_operator, src, dst);
    mf_data.copy_constrained_values(src, dst);
}

template<int dim, int fe_degree, typename number>
void ADROperator<dim, fe_degree, number>::initialize_dof_vector(
    DistributedDeviceVector<number> &vec) const
{
    mf_data.initialize_dof_vector(vec);
}

template<int dim, int fe_degree, typename number>
void ADROperator<dim, fe_degree, number>::compute_rhs(
    DistributedDeviceVector<number> &dst) const
{
    // TODO: Since Portable::FEEvaluation::read_dof_values_plain()
    // is missing, use two set of constraints to handle
    // inhomogeneous dirichlet boundary conditions

    dst = 0;
    LocalRhsEvaluator<dim, fe_degree, number> loc_rhs_operator(
        forcing_coefs.get_values());
    mf_data.cell_loop(loc_rhs_operator, dst, dst);
}

template<int dim, int fe_degree, typename number>
void ADROperator<dim, fe_degree, number>::compute_diagonal()
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
        ADROperatorQuad<dim, fe_degree, number>(),
        EvaluationFlags::gradients,
        EvaluationFlags::gradients);

    // TODO: Is this not needed anymore?
    // Set Dirichlet constrained DoFs to 1 to allow invertibility
    // this->set_constrained_entries_to_one(inverse_diagonal);

    number *raw_diagonal = inverse_diagonal.get_values();

    Kokkos::parallel_for(
        inverse_diagonal.locally_owned_size(), KOKKOS_LAMBDA(int i) {
            raw_diagonal[i] = 1. / raw_diagonal[i];
        });
}

template<int dim, int fe_degree, typename number>
std::shared_ptr<DiagonalMatrix<DistributedDeviceVector<number>>>
ADROperator<dim, fe_degree, number>::get_matrix_diagonal_inverse() const
{
    return inverse_diagonal_entries;
}

template<int dim, int fe_degree, typename number>
types::global_dof_index ADROperator<dim, fe_degree, number>::m() const
{
    // TODO: This returns the local size, right? I have to understand
    // how vmult is distributed across MPI processes
    return mf_data.get_vector_partitioner()->size();
}

template<int dim, int fe_degree, typename number>
types::global_dof_index ADROperator<dim, fe_degree, number>::n() const
{
    return m();
}

template<int dim, int fe_degree, typename number>
number ADROperator<dim, fe_degree, number>::el(
    const types::global_dof_index row,
    [[maybe_unused]] const types::global_dof_index col) const
{
    return 1.0 / (*inverse_diagonal_entries)(row, row);
}

#endif
