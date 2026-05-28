#ifndef ADR_OPERATOR_HPP
#define ADR_OPERATOR_HPP

#include <deal.II/base/point.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/fe_evaluation.h>

#include <functional>
// #include "bench_config.h"

using namespace dealii;


namespace BenchConfig {
    using RealNumber = double;
    constexpr unsigned int simd_width = 1;
};


template<typename number>
using DistributedVector = LinearAlgebra::distributed::Vector<number>;

using VectorizedArrayType = VectorizedArray<BenchConfig::RealNumber, BenchConfig::simd_width>;

template<int dim>
using ArrayPoint = Point<dim, VectorizedArrayType>;

template<int dim>
using MatrixFreeType = MatrixFree<dim, BenchConfig::RealNumber, VectorizedArrayType>;



template<int dim, int fe_degree, typename number, 
         typename AdvectionFunc, typename DiffusionFunc, typename ReactionFunc>
class ADROperator : public MatrixFreeOperators::Base<dim, DistributedVector<number>, VectorizedArrayType>
{
public:
    ADROperator(const AdvectionFunc &advection_coeff_,
                const DiffusionFunc &diffusion_coeff_,
                const ReactionFunc &reaction_coeff_);

    virtual void compute_diagonal() override;

    const AdvectionFunc advection_coeff;
    const DiffusionFunc diffusion_coeff;
    const ReactionFunc reaction_coeff;

private:
    virtual void apply_add(DistributedVector<number> &dst,
                           const DistributedVector<number> &src) const override;

    void
    local_apply(const MatrixFreeType<dim> &data,
                DistributedVector<number> &dst,
                const DistributedVector<number> &src,
                const std::pair<unsigned int,
                                unsigned int> &cell_range) const;

    void empty_face_apply(const MatrixFreeType<dim> &,
                          DistributedVector<number> &,
                          const DistributedVector<number> &,
                          const std::pair<unsigned int, unsigned int> &) const {}

    void local_compute_diagonal(
        FEEvaluation<dim, fe_degree,
                          fe_degree + 1, 1, number, VectorizedArrayType> &integrator) const;
};

template<int dim, int fe_degree, typename number, 
         typename AdvectionFunc, typename DiffusionFunc, typename ReactionFunc>
ADROperator<dim, fe_degree, number, AdvectionFunc, DiffusionFunc, ReactionFunc>::ADROperator(
    const AdvectionFunc &advection_coeff_,
    const DiffusionFunc &diffusion_coeff_,
    const ReactionFunc &reaction_coeff_
)
    : MatrixFreeOperators::Base<dim, DistributedVector<number>, VectorizedArrayType>()
    , advection_coeff(advection_coeff_)
    , diffusion_coeff(diffusion_coeff_)
    , reaction_coeff(reaction_coeff_)
{}

template<int dim, int fe_degree, typename number, 
         typename AdvectionFunc, typename DiffusionFunc, typename ReactionFunc>
void ADROperator<dim, fe_degree, number, AdvectionFunc, DiffusionFunc, ReactionFunc>::compute_diagonal()
{
    this->inverse_diagonal_entries.reset(
        new DiagonalMatrix<DistributedVector<number>>());
    DistributedVector<number> &inverse_diagonal =
        this->inverse_diagonal_entries->get_vector();
    // Resizes inverse_diagonal to the number of DoFs
    this->data->initialize_dof_vector(inverse_diagonal);

    MatrixFreeTools::compute_diagonal(*this->data, inverse_diagonal,
                                      &ADROperator::local_compute_diagonal, this);

    this->set_constrained_entries_to_one(inverse_diagonal);

    for (unsigned int i = 0; i < inverse_diagonal.locally_owned_size(); ++i) {
        inverse_diagonal.local_element(i) =
            //cast to avoid warning
            static_cast<number>(1.) / inverse_diagonal.local_element(i);
    }
}

template<int dim, int fe_degree, typename number, 
         typename AdvectionFunc, typename DiffusionFunc, typename ReactionFunc>
void ADROperator<dim, fe_degree, number, AdvectionFunc, DiffusionFunc, ReactionFunc>::local_compute_diagonal(
    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number, VectorizedArrayType> &phi) const
{
    // Evaluates against unit vectors implicitly handled by compute_diagonal
    phi.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);

    for (const unsigned int q : phi.quadrature_point_indices()) {
        const ArrayPoint<dim> q_p = phi.quadrature_point(q);
        const VectorizedArrayType q_val = phi.get_value(q);

        // Standard Galerkin evaluation (no Integration By Parts for advection)
        phi.submit_value(reaction_coeff(q_p) * q_val + 
                         advection_coeff(q_p) * phi.get_gradient(q), q);
                         
        phi.submit_gradient(diffusion_coeff(q_p) * phi.get_gradient(q), q);
    }

    phi.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
}

template<int dim, int fe_degree, typename number, 
         typename AdvectionFunc, typename DiffusionFunc, typename ReactionFunc>
void ADROperator<dim, fe_degree, number, AdvectionFunc, DiffusionFunc, ReactionFunc>::apply_add(
    DistributedVector<number> &dst, const DistributedVector<number> &src) const
{
    // Inner and boundary faces are empty because we are using Continuous Galerkin 
    // and standard weak form (no upwinding or IBP surface terms).
    this->data->loop(&ADROperator::local_apply,
                     &ADROperator::empty_face_apply,
                     &ADROperator::empty_face_apply,
                     this, dst, src);
}

template<int dim, int fe_degree, typename number, 
         typename AdvectionFunc, typename DiffusionFunc, typename ReactionFunc>
void ADROperator<dim, fe_degree, number, AdvectionFunc, DiffusionFunc, ReactionFunc>::local_apply(
    const MatrixFreeType<dim> &data,
    DistributedVector<number> &dst,
    const DistributedVector<number> &src,
    const std::pair<unsigned int, unsigned int> &cell_range) const
{
    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number, VectorizedArrayType> phi(data);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell) {
        phi.reinit(cell);
        phi.read_dof_values(src);
        phi.evaluate(EvaluationFlags::gradients | EvaluationFlags::values);
        
        for (const unsigned int q : phi.quadrature_point_indices()) {
            const ArrayPoint<dim> q_p = phi.quadrature_point(q);
            const VectorizedArrayType q_val = phi.get_value(q);

            phi.submit_value(reaction_coeff(q_p) * q_val + 
                             advection_coeff(q_p) * phi.get_gradient(q), q);
                             
            phi.submit_gradient(diffusion_coeff(q_p) * phi.get_gradient(q), q);
        }
        phi.integrate(EvaluationFlags::gradients | EvaluationFlags::values);
        phi.distribute_local_to_global(dst);
    }
}

#endif
