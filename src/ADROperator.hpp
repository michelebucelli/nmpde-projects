#ifndef ADR_OPERATOR_HPP
#define ADR_OPERATOR_HPP

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/fe_evaluation.h>

using namespace dealii;

template<typename number>
using DistributedVector = LinearAlgebra::distributed::Vector<number>;

// NOTE: I could probably remove the template on the number

template<int dim, int fe_degree, typename number>
class ADROperator :
    public MatrixFreeOperators::Base<dim, DistributedVector<number>>
{
public:
    ADROperator();

    virtual void compute_diagonal() override;

private:
    virtual void
    apply_add(DistributedVector<number> &dst,
              const DistributedVector<number> &src) const override;

    void
    local_apply(const MatrixFree<dim, number> &data,
                DistributedVector<number> &dst,
                const DistributedVector<number> &src,
                const std::pair<unsigned int,
                                unsigned int> &cell_range) const;

    void local_compute_diagonal(
        FEEvaluation<dim, fe_degree,
                          fe_degree + 1, 1, number> &integrator) const;
};

template<int dim, int fe_degree, typename number>
ADROperator<dim, fe_degree, number>::ADROperator()
    : MatrixFreeOperators::Base<dim, DistributedVector<number>>()
{}

template<int dim, int fe_degree, typename number>
void ADROperator<dim, fe_degree, number>::compute_diagonal()
{
    this->inverse_diagonal_entries.reset(
        new DiagonalMatrix<DistributedVector<number>>());
    DistributedVector<double> &inverse_diagonal =
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
            1. / inverse_diagonal.local_element(i);
    }
}

template<int dim, int fe_degree, typename number>
void ADROperator<dim, fe_degree, number>::local_compute_diagonal(
    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> &phi) const
{
    // No need to read any src value, we're multiplying by
    // unit vectors to compute the diagonal

    phi.evaluate(EvaluationFlags::gradients);

    for (const unsigned int q : phi.quadrature_point_indices()) {
        phi.submit_gradient(phi.get_gradient(q), q);
    }

    phi.integrate(EvaluationFlags::gradients);
}

template<int dim, int fe_degree, typename number>
void ADROperator<dim, fe_degree, number>::
    apply_add(DistributedVector<number> &dst,
              const DistributedVector<number> &src) const
{
    // Calls local_apply in batches to allow a degree of parallelization
    this->data->cell_loop(&ADROperator::local_apply, this, dst, src);
}

template<int dim, int fe_degree, typename number>
void ADROperator<dim, fe_degree, number>::local_apply(
    const MatrixFree<dim, number> &data,
    DistributedVector<number> &dst,
    const DistributedVector<number> &src,
    const std::pair<unsigned int, unsigned int> &cell_range) const
{
    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi(data);

    // WARNING: Make sure to compile with -march=native to vectorize
    // this loop. Additionally, any coefficient evaluation function
    // called within this loop should be templated on the data type,
    // so that it supports VectorizedArray values

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
        phi.evaluate(EvaluationFlags::gradients);
        for (const unsigned int q : phi.quadrature_point_indices()) {
            // Writes a contribution on the current quadrature point
            // (the gradient of u at q), which will be multiplied by the
            // gradient of the test functions. Note that get_gradient()
            // multiplies by J_cell^-T and the quadrature weight,
            // while submit_gradients() multiplies by J_cell^-1 and JxW.

            // TLDR: we're preparing to compute B_cell^T * D_cell * B_cell
            phi.submit_gradient(phi.get_gradient(q), q);
        }
        // Multiplies the values stored at the quadrature points by the
        // gradient of each basis function, and sums up the result for all
        // quadrature nodes, performing the 
        phi.integrate(EvaluationFlags::gradients);
        // Multiplies the result by P_cell^T and adds into dst
        phi.distribute_local_to_global(dst);
    }
}

#endif
