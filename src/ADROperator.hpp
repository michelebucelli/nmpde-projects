#ifndef ADR_OPERATOR_HPP
#define ADR_OPERATOR_HPP

#include <deal.II/base/point.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/fe_evaluation.h>

using namespace dealii;

template<typename number>
using DistributedVector = LinearAlgebra::distributed::Vector<number>;

template<int fe_degree, typename ADRSpec>
class ADROperator :
    public MatrixFreeOperators::Base<
        ADRSpec::dim,
        DistributedVector<typename ADRSpec::ScalarType::value_type>,
        typename ADRSpec::ScalarType>
{
public:
    ADROperator();

    virtual void compute_diagonal() override;

private:
    static constexpr int dim = ADRSpec::dim;

    using number = typename ADRSpec::ScalarType::value_type;
    // TODO: Consider using auto inside the quadrature loops
    // instead of using these aliases
    using ArrayScalar = typename ADRSpec::ScalarType;
    using ArrayPoint = typename ADRSpec::PointType;
    using ArrayTensor = typename ADRSpec::TensorType;
    using MatrixFreeType = MatrixFree<dim, number, ArrayScalar>;

    virtual void
    apply_add(DistributedVector<number> &dst,
              const DistributedVector<number> &src) const override;

    void
    local_apply(const MatrixFreeType &data,
                DistributedVector<number> &dst,
                const DistributedVector<number> &src,
                const std::pair<unsigned int,
                                unsigned int> &cell_range) const;

    void
    empty_face_apply(const MatrixFreeType &,
                     DistributedVector<number> &,
                     const DistributedVector<number> &,
                     const std::pair<unsigned int, unsigned int> &) const {}

    void local_compute_diagonal(
        FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number, ArrayScalar>
        &integrator) const;

    ADRSpec local_spec;
};

template<int fe_degree, typename ADRSpec>
ADROperator<fe_degree, ADRSpec>::ADROperator()
    : MatrixFreeOperators::Base<dim, DistributedVector<number>, ArrayScalar>()
{}

template<int fe_degree, typename ADRSpec>
void ADROperator<fe_degree, ADRSpec>::compute_diagonal()
{
    this->inverse_diagonal_entries.reset(
        new DiagonalMatrix<DistributedVector<number>>());
    DistributedVector<number> &inverse_diagonal =
        this->inverse_diagonal_entries->get_vector();
    // Resizes inverse_diagonal to the number of DoFs
    this->data->initialize_dof_vector(inverse_diagonal);

    // Applies local_compute_diagonal to unit vectors
    // to compute the operator diagonal
    MatrixFreeTools::compute_diagonal(*this->data,
                                      inverse_diagonal,
                                      &ADROperator::local_compute_diagonal,
                                      this);

    // Set Dirichlet constrained DoFs to 1 to allow invertibility
    this->set_constrained_entries_to_one(inverse_diagonal);

    for (unsigned int i = 0; i < inverse_diagonal.locally_owned_size(); ++i) {
        inverse_diagonal.local_element(i) =
            // cast to avoid warning
            static_cast<number>(1.) / inverse_diagonal.local_element(i);
    }
}

template<int fe_degree, typename ADRSpec>
void ADROperator<fe_degree, ADRSpec>::local_compute_diagonal(
    FEEvaluation<
        dim, fe_degree, fe_degree + 1, 1, number, ArrayScalar> &phi) const
{
    // No need to read any src value, we're multiplying by unit
    // vectors to compute the diagonal

    phi.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);

    // TODO: Move this loop inside a function to avoid duplicating
    // the weak form evaluation in local_apply
    for (const unsigned int q : phi.quadrature_point_indices()) {
        const ArrayPoint q_p = phi.quadrature_point(q);

        const ArrayScalar q_val = phi.get_value(q);
        const ArrayTensor q_grad = phi.get_gradient(q);

        phi.submit_value(local_spec.reaction(q_p) * q_val +
                         local_spec.advection(q_p) * q_grad, q);

        phi.submit_gradient(local_spec.diffusion(q_p) * q_grad, q);
    }

    phi.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
}

template<int fe_degree, typename ADRSpec>
void ADROperator<fe_degree, ADRSpec>::apply_add(
    DistributedVector<number> &dst,
    const DistributedVector<number> &src) const
{
    // Calls local_apply in batches to allow a degree of parallelization
    this->data->loop(&ADROperator::local_apply,
                     &ADROperator::empty_face_apply,
                     &ADROperator::empty_face_apply,
                     this,
                     dst,
                     src);
}

template<int fe_degree, typename ADRSpec>
void ADROperator<fe_degree, ADRSpec>::local_apply(
    const MatrixFreeType &data,
    DistributedVector<number> &dst,
    const DistributedVector<number> &src,
    const std::pair<unsigned int, unsigned int> &cell_range) const
{
    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number, ArrayScalar>
    phi(data);

    // TODO: Check if precomputing coefficients is faster
    // than re-evaluating them on the fly

    for (unsigned int cell = cell_range.first;
                      cell < cell_range.second; ++cell) {

        phi.reinit(cell);
        // Reads src values (u) on the DoFs of the current cell,
        // indirectly performing P_cell * u
        phi.read_dof_values(src);
        // Interpolates the src gradient over the quadrature points
        // by performing B_ref_cell * u. The evaluation exploits the tensor
        // product structure of the basis functions on the reference cell,
        // allowing for a product of the type B_grad,x * U * B_val,y where U
        // is a reshaped version of u and B_ref_cell = B_grad,x tdot B_val,y
        phi.evaluate(EvaluationFlags::gradients | EvaluationFlags::values);
        for (const unsigned int q : phi.quadrature_point_indices()) {
            const ArrayPoint q_p = phi.quadrature_point(q);

            const ArrayScalar q_val = phi.get_value(q);
            const ArrayTensor q_grad = phi.get_gradient(q);

            // Writes a contribution on the current quadrature point
            // (the gradient of u at q), which will be multiplied by the
            // gradient of the test functions. Note that get_gradient()
            // multiplies by J_cell^-T and the quadrature weight,
            // while submit_gradients() multiplies by J_cell^-1 and JxW.
            // TLDR: we're preparing to compute B_cell^T * D_cell * B_cell

            phi.submit_value(local_spec.reaction(q_p) * q_val +
                             local_spec.advection(q_p) * q_grad, q);

            phi.submit_gradient(local_spec.diffusion(q_p) * q_grad, q);
        }
        // Multiplies the values stored at the quadrature points by the
        // gradient of each basis function, and sums up the result for all
        // quadrature nodes
        phi.integrate(EvaluationFlags::gradients | EvaluationFlags::values);
        // Multiplies the result by P_cell^T and adds into dst
        phi.distribute_local_to_global(dst);
    }
}

#endif
