#ifndef CONCRETE_SPEC_HPP
#define CONCRETE_SPEC_HPP

#include <deal.II/base/point.h>
#include <deal.II/base/function.h>

#include <memory>

#include "ADRSpec.hpp"

using namespace dealii;

constexpr unsigned int dim = 2;
constexpr unsigned int fe_degree = 2;

class ExactSolution : public Function<dim>
{
public:
    virtual double value(const Point<dim> &p,
                         const unsigned int = 0) const override
    {
        return std::sin(2 * M_PI * p[0]) * std::cos(2 * M_PI * p[1]);
    }

    virtual Tensor<1, dim> gradient(const Point<dim> &p,
                                    const unsigned int = 0) const override
    {
        Tensor<1, dim> res;
        res[0] = 2 * M_PI * std::cos(2 * M_PI * p[0])
                          * std::cos(2 * M_PI * p[1]);
        res[1] = -2 * M_PI * std::sin(2 * M_PI * p[0])
                           * std::sin(2 * M_PI * p[1]);
        return res;
    }
};

template<typename Scalar, typename Point>
struct NeumannBC : public NeumannFunction<Scalar, Point>
{
    Scalar value(const Point &p) const override
    {
        return -2 * M_PI * std::cos(2 * M_PI * p[0])
                         * std::cos(2 * M_PI * p[1]);
    }
};

template<int dim, typename Scalar, typename Tensor, typename Point>
struct ConcreteSpec : public ADRSpec<dim, Scalar, Tensor, Point>
{
    ConcreteSpec() {}

    Tensor advection(const Point &p) const override
    {
        Tensor res;
        res[0] = p[1];
        res[1] = p[0];
        return res;
    }

    Scalar diffusion(const Point &p) const override
    {
        // WARNING: Avoid having mu = 0 on boundaries, the problem
        // could lose well-posedness if you want to apply neumann there
        return 1.0 + p[0] + p[1];
    }

    Scalar reaction(const Point &p) const override
    {
        return p[1];
    }

    Scalar forcing(const Point &p) const override
    {
        return 8 * M_PI * M_PI * (1.0 + p[0] + p[1])
                 * std::sin(2 * M_PI * p[0])
                 * std::cos(2 * M_PI * p[1]) -
               2 * M_PI * std::cos(2 * M_PI * p[0])
                        * std::cos(2 * M_PI * p[1]) +
               2 * M_PI * std::sin(2 * M_PI * p[0])
                        * std::sin(2 * M_PI * p[1]) +
                   p[1] * std::sin(2 * M_PI * p[0])
                        * std::cos(2 * M_PI * p[1]) +
               2 * M_PI * p[1] * std::cos(2 * M_PI * p[0])
                        * std::cos(2 * M_PI * p[1]) -
               2 * M_PI * p[0] * std::sin(2 * M_PI * p[0])
                        * std::sin(2 * M_PI * p[1]);
    }

    void fill_dirichlet_map(
        BoundaryFuncMap<Function<dim>> &function_map) const override
    {
        function_map[1] = std::make_unique<ExactSolution>();
        function_map[2] = std::make_unique<ExactSolution>();
        function_map[3] = std::make_unique<ExactSolution>();
    }

    void fill_neumann_map(BoundaryFuncMap<NeumannFunction<Scalar, Point>>
                          &function_map) const override
    {
        function_map[0] = std::make_unique<NeumannBC<Scalar, Point>>();
    }
};

#endif
