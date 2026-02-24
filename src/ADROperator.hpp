#ifndef ADR_OPERATOR_HPP
#define ADR_OPERATOR_HPP

#include <deal.II/base/point.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/fe_evaluation.h>

#include <functional>

using namespace dealii;

template<typename number>
using DistributedVector = LinearAlgebra::distributed::Vector<number>;

template<int dim, typename number>
using ArrayPoint = Point<dim, VectorizedArray<number>>;

template<int dim, typename number>
using ScalarArrayFunction = std::function<
    VectorizedArray<number>(const ArrayPoint<dim, number> &)>;

template<int dim, typename number>
using VectorArrayFunction = std::function<
    Tensor<1, dim, VectorizedArray<number>>(
        const ArrayPoint<dim, number> &)>;

// NOTE: I could probably remove the template on the number

template<int dim, int fe_degree, typename number>
class ADROperator :
    public MatrixFreeOperators::Base<dim, DistributedVector<number>>
{
public:
    ADROperator(const VectorArrayFunction<dim, number> &advection_coeff_,
                const ScalarArrayFunction<dim, number> &diffusion_coeff_,
                const ScalarArrayFunction<dim, number> &reaction_coeff_);

    virtual void compute_diagonal() override;

    const VectorArrayFunction<dim, number> advection_coeff;
    const ScalarArrayFunction<dim, number> diffusion_coeff;
    const ScalarArrayFunction<dim, number> reaction_coeff;

    // TODO: I don't really like this one here, since dirichlet
    // function map is defined inside MatrixFreeADR. Moreover
    // it would be better to define a getter
    std::map<types::boundary_id, ScalarArrayFunction<dim, number>>
    neumann_bcs;

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

    void local_apply_inner_faces(
        const MatrixFree<dim, number> &data,
        DistributedVector<number> &dst,
        const DistributedVector<number> &src,
        const std::pair<unsigned int, unsigned int> &face_range) const;

    void local_apply_boundary_faces(
        const MatrixFree<dim, number> &data,
        DistributedVector<number> &dst,
        const DistributedVector<number> &src,
        const std::pair<unsigned int, unsigned int> &face_range) const;

    void local_compute_diagonal(
        FEEvaluation<dim, fe_degree,
                          fe_degree + 1, 1, number> &integrator) const;
};

template<int dim, int fe_degree, typename number>
ADROperator<dim, fe_degree, number>::ADROperator(
    const VectorArrayFunction<dim, number> &advection_coeff_,
    const ScalarArrayFunction<dim, number> &diffusion_coeff_,
    const ScalarArrayFunction<dim, number> &reaction_coeff_
)
    : MatrixFreeOperators::Base<dim, DistributedVector<number>>()
    , advection_coeff(advection_coeff_)
    , diffusion_coeff(diffusion_coeff_)
    , reaction_coeff(reaction_coeff_)
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

    // TODO: If we actually need the diagonal of the operator (I guess
    // it's mandatory for multigrid) we need to specify an additional
    // function pointer to account for the contribution of the boundary
    // component of the advection term

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
    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> &) const
{
    // NOTE: The MS converges even if this function is not implemented
    // so I guess it's only needed with multigrid

    /*

    // No need to read any src value, we're multiplying by
    // unit vectors to compute the diagonal
    phi.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);

    // TODO: Avoid duplicating this loop, reuse the code in local_apply()

    for (const unsigned int q : phi.quadrature_point_indices()) {
        const ArrayPoint<dim, number> q_p = phi.quadrature_point(q);
        const VectorizedArray<double> q_val = phi.get_value(q);

        phi.submit_gradient(
            //-advection_coeff(q_p) * q_val +
            diffusion_coeff(q_p) * phi.get_gradient(q), q);

        phi.submit_value(
            reaction_coeff(q_p) * q_val, q);
    }

    phi.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
    */
}

template<int dim, int fe_degree, typename number>
void ADROperator<dim, fe_degree, number>::
    apply_add(DistributedVector<number> &dst,
              const DistributedVector<number> &src) const
{
    // Calls local_apply in batches to allow a degree of parallelization
    this->data->loop(&ADROperator::local_apply,
                     &ADROperator::local_apply_inner_faces,
                     &ADROperator::local_apply_boundary_faces,
                     this,
                     dst,
                     src);
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
        phi.evaluate(EvaluationFlags::gradients | EvaluationFlags::values);
        for (const unsigned int q : phi.quadrature_point_indices()) {
            // Writes a contribution on the current quadrature point
            // (the gradient of u at q), which will be multiplied by the
            // gradient of the test functions. Note that get_gradient()
            // multiplies by J_cell^-T and the quadrature weight,
            // while submit_gradients() multiplies by J_cell^-1 and JxW.
            // TLDR: we're preparing to compute B_cell^T * D_cell * B_cell

            const ArrayPoint<dim, number> q_p = phi.quadrature_point(q);
            const VectorizedArray<double> q_val = phi.get_value(q);

            // NOTE: The advection term in the weak formulation is split
            // using integration by parts into a domain integral and a
            // boundary integral to avoid computing the divergence of the
            // advection coefficient

            phi.submit_gradient(
                -advection_coeff(q_p) * q_val +
                diffusion_coeff(q_p) * phi.get_gradient(q), q);

            phi.submit_value(
                reaction_coeff(q_p) * q_val, q);
        }
        // Multiplies the values stored at the quadrature points by the
        // gradient of each basis function, and sums up the result for all
        // quadrature nodes
        phi.integrate(EvaluationFlags::gradients | EvaluationFlags::values);
        // Multiplies the result by P_cell^T and adds into dst
        phi.distribute_local_to_global(dst);
    }
}

template<int dim, int fe_degree, typename number>
void ADROperator<dim, fe_degree, number>::local_apply_inner_faces(
    const MatrixFree<dim, number> &,
    DistributedVector<number> &,
    const DistributedVector<number> &,
    const std::pair<unsigned int, unsigned int> &) const
{
    // TODO: Nothing to do here, check if we can avoid defining this
    // function and still be able to run MatrixFree::loop()
}

template<int dim, int fe_degree, typename number>
void ADROperator<dim, fe_degree, number>::local_apply_boundary_faces(
    const MatrixFree<dim, number> &data,
    DistributedVector<number> &dst,
    const DistributedVector<number> &src,
    const std::pair<unsigned int, unsigned int> &face_range) const
{
    FEFaceEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi_face(data);

    // Accounting for the contribution of the boundary component
    // of the advection term (split using integration by parts)

    for (unsigned int face = face_range.first;
                      face < face_range.second; ++face) {

        phi_face.reinit(face);
        const types::boundary_id b_id = phi_face.boundary_id();

        if (neumann_bcs.count(b_id)) {
            const ScalarArrayFunction<dim, double> &bc_function =
                neumann_bcs.at(b_id);

            phi_face.gather_evaluate(src, EvaluationFlags::values);

            for (const unsigned int q : phi_face.quadrature_point_indices()) {
                const ArrayPoint<dim, double> q_p = phi_face.quadrature_point(q);

                phi_face.submit_value(
                    advection_coeff(q_p) *
                    phi_face.normal_vector(q) * phi_face.get_value(q) +
                    bc_function(q_p), q);
            }
            phi_face.integrate_scatter(EvaluationFlags::values, dst);
        }
    }
}


#endif
