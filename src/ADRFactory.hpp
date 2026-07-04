#ifndef ADR_FACTORY_HPP
#define ADR_FACTORY_HPP

#include <deal.II/base/numbers.h>

#ifdef BUILD_WITH_KOKKOS
#include "portable/MatrixFreeADR.hpp"
#else
#include "MatrixBasedADR.hpp"
#include "MatrixFreeADR.hpp"
#endif

using namespace dealii;

class ADRFactory
{
public:
#ifndef BUILD_WITH_KOKKOS
    template<unsigned int dim,
             typename number,
             template <unsigned int dim_,
                   typename Scalar_,
                   typename Tensor_,
                   typename Point_> typename ADRSpec>
    static auto matrix_based(
        const parallel::distributed::Triangulation<dim> &mesh,
        unsigned int fe_degree)
    {
        using ScalarPoint = Point<dim, number>;
        using ScalarTensor = Tensor<1, dim, number>;

        using Spec = ADRSpec<dim, number, ScalarTensor, ScalarPoint>;

        static_assert(Spec::dim == dim);

        return MatrixBasedADR<
            ADRSpec<dim, number, ScalarTensor, ScalarPoint>>(mesh, fe_degree);
    }
#endif

    template<unsigned int dim,
             unsigned int fe_degree,
             typename number,
             template <unsigned int dim_,
                       typename Scalar_,
                       typename Tensor_,
                       typename Point_> typename ADRSpec,
             // Ignored when BUILD_WITH_KOKKOS is defined
             unsigned int simd_width =
                internal::VectorizedArrayWidthSpecifier<number>::max_width>
    static auto matrix_free(
        const parallel::distributed::Triangulation<dim> &mesh)
    {
#ifdef BUILD_WITH_KOKKOS
        using ArrayScalar = number;
#else
        using ArrayScalar = VectorizedArray<number, simd_width>;
#endif
        using ArrayPoint = Point<dim, ArrayScalar>;
        using ArrayTensor = Tensor<1, dim, ArrayScalar>;

        using Spec = ADRSpec<dim, ArrayScalar, ArrayTensor, ArrayPoint>;

        static_assert(Spec::dim == dim);

        return MatrixFreeADR<fe_degree, Spec>(mesh);
    }
};

#endif
