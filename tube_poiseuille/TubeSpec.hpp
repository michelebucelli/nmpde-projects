#ifndef TUBE_SPEC_HPP
#define TUBE_SPEC_HPP

#include <deal.II/base/point.h>
#include <deal.II/base/function.h>

#include <memory>

#include "ADRSpec.hpp"

using namespace dealii;
constexpr unsigned int dim = 3;

template<int dim, typename Scalar, typename Tensor, typename Point>
struct TubeSpec : public ADRSpec<dim, Scalar, Tensor, Point>
{

    static constexpr double tube_radius = 1.0;
    static constexpr double max_flow_velocity = 1.0;
    static constexpr double Peclet = 1.0;

    DEAL_II_HOST_DEVICE TubeSpec() {}

    DEAL_II_HOST_DEVICE Tensor advection(const Point &p) const override
    {
        Scalar r_squared = p[1]*p[1] + p[2]*p[2];
        Scalar poiseuille_velocity = max_flow_velocity * (1.0 - r_squared/(tube_radius*tube_radius));
        Tensor res;
        res[0] = poiseuille_velocity;
        res[1] = 0.0;
        res[2] = 0.0;
        return res;
    }

    DEAL_II_HOST_DEVICE Scalar diffusion(const Point &p) const override
    {
        return 2*tube_radius*max_flow_velocity/Peclet;
    }

    DEAL_II_HOST_DEVICE Scalar reaction(const Point &p) const override
    {
        return 0.0;
    }

    DEAL_II_HOST_DEVICE Scalar forcing(const Point &p) const override
    {
        return 0.0;
    }
    
    class ParabolicDistributionInlet : public Function<dim>
    {
    public:   //forced to specify dealii:: because the compiler was confused by the alias in the template
        virtual double value(const dealii::Point<dim, double> &p,
                            const unsigned int = 0) const override
        {
            double c_in = 1.0;
            double r_squared = p[1]*p[1] + p[2]*p[2];
            return c_in  * (1.0 - r_squared/(tube_radius*tube_radius));
        }
    };

    void fill_dirichlet_map(
        BoundaryFuncMap<Function<dim>> &function_map) const override
    {
        function_map[0] = std::make_unique<Functions::ZeroFunction<dim>>();
        function_map[1] = std::make_unique<ParabolicDistributionInlet>();
    }

    void fill_neumann_map(BoundaryFuncMap<NeumannFunction<Scalar, Point>>
                          &function_map) const override
    {}
};

#endif