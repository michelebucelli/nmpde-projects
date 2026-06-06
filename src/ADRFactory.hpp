#ifndef ADR_FACTORY_HPP
#define ADR_FACTORY_HPP

#include "MatrixBasedADR.hpp"
#include "MatrixFreeADR.hpp"

using namespace dealii;

class ADRFactory
{
public:
    template<int dim,
             typename number,
             template <int dim_,
                   typename Scalar_,
                   typename Tensor_,
                   typename Point_> typename ADRSpec>
    static auto matrix_based(const Triangulation<dim> &mesh, int fe_degree)
    {
        using ScalarPoint = Point<dim, number>;
        using ScalarTensor = Tensor<1, dim, number>;

        return MatrixBasedADR<
            ADRSpec<dim, number, ScalarTensor, ScalarPoint>>(mesh, fe_degree);
    }

    template<int dim,
             int fe_degree,
             unsigned int simd_width,
             typename number,
             template <int dim_,
                       typename Scalar_,
                       typename Tensor_,
                       typename Point_> typename ADRSpec>
    static auto matrix_free(const unsigned int N)
    {
        using ArrayScalar = VectorizedArray<number, simd_width>;
        using ArrayPoint = Point<dim, ArrayScalar>;
        using ArrayTensor = Tensor<1, dim, ArrayScalar>;

        return MatrixFreeADR<
            fe_degree, ADRSpec<dim, ArrayScalar, ArrayTensor, ArrayPoint>>(N);
    }
};

#endif
