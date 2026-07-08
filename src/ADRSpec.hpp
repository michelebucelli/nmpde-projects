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

template<typename Function>
using BoundaryFuncMap = std::map<types::boundary_id,
                                 std::unique_ptr<Function>>;

template<typename Scalar, typename Point>
struct NeumannFunction
{
    virtual Scalar value(const Point &p) const = 0;
    virtual ~NeumannFunction() = default;
};

template <unsigned int dim_, typename Scalar, typename Tensor, typename Point>
struct ADRSpec
{
    static constexpr unsigned int dim = dim_;

    using ScalarType = Scalar;
    using TensorType = Tensor;
    using PointType = Point;

    DEAL_II_HOST_DEVICE virtual Tensor advection(const Point &p) const = 0; // beta
    DEAL_II_HOST_DEVICE virtual Scalar diffusion(const Point &p) const = 0; // mu
    DEAL_II_HOST_DEVICE virtual Scalar reaction(const Point &p) const = 0; // gamma
    DEAL_II_HOST_DEVICE virtual Scalar forcing(const Point &p) const = 0; // f

    virtual void fill_dirichlet_map(
        BoundaryFuncMap<Function<dim>> &function_map) const = 0;

    virtual void fill_neumann_map(
        BoundaryFuncMap<NeumannFunction<Scalar, Point>>
        &function_map) const = 0;

    // virtual constructor to avoid compiler warnings
    virtual ~ADRSpec() = default;
};

#endif // ADR_SPEC_HPP
