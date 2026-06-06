#ifndef ADR_SPEC_HPP
#define ADR_SPEC_HPP

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

using namespace dealii;

// This is a design proposal for a description of the ADR problem
// The design goals are (in decreasing order of importance):
//   1. centralize the ADR parameters (coefficients, boundary conditions...)
//      in a single entity
//
//   2. make this description stateless, so that it can be passed around
//      to other objects or devices without copy
//
//   3. use the same exact description for both matrix based and matrix free simulations
//
//   4. remove function call overheads by favouring inlining (we already solved that)
//
//
// Here are the bits that I still have to understand:
//   * vectorized implementations: the coefficient functions must be parametric on the
//     return type. In the matrix free context, we need them to yield VectorizedArray
//     I think we can template on the return type -> DONE
//
//   * how to pass functions to maps: we need to be careful about memory leaks. Since
//     the AdrSpec is stateless, we can't pass references to Function because they will
//     die after the `fill_.._map` call.
//
//     For the same reason we can't pass just pointers, because the memory will leak.
//     We must pass smart pointers, but then we will need to clone the map for using
//     `VectorTools::interpolate_boundary_values`. It's not a big deal

template<typename Function>
using BoundaryFuncMap = std::map<types::boundary_id,
                                 std::unique_ptr<Function>>;

template<typename Scalar, typename Point>
struct NeumannFunction
{
    virtual Scalar value(const Point &p) const = 0;
    virtual ~NeumannFunction() = default;
};

template <int dim_, typename Scalar, typename Tensor, typename Point>
struct ADRSpec
{
    static constexpr int dim = dim_;

    using ScalarType = Scalar;
    using TensorType = Tensor;
    using PointType = Point;

    virtual Tensor advection(const Point &p) const = 0; // beta
    virtual Scalar diffusion(const Point &p) const = 0; // mu
    virtual Scalar reaction(const Point &p) const = 0; // gamma
    virtual Scalar forcing(const Point &p) const = 0; // f

    virtual void fill_dirichlet_map(
        BoundaryFuncMap<Function<dim>> &function_map) const = 0;

    virtual void fill_neumann_map(
        BoundaryFuncMap<NeumannFunction<Scalar, Point>>
        &function_map) const = 0;

    // virtual constructor to avoid compiler warnings
    virtual ~ADRSpec() = default;
};

#endif // ADR_SPEC_HPP
