#ifndef TUBE_SPEC_HPP
#define TUBE_SPEC_HPP

#include <deal.II/base/point.h>
#include <deal.II/base/function.h>

#include <memory>

#include "ADRSpec.hpp"

using namespace dealii;

constexpr unsigned int dim = 3;

// Provide default template arguments to allow querying
// TubeSpec<>::tube_radius
template<unsigned int dim = 3,
         typename Scalar = double,
         typename Tensor = dealii::Tensor<1, dim>,
         typename Point = Point<dim>>
struct TubeSpec : public ADRSpec<dim, Scalar, Tensor, Point>
{
    static constexpr double tube_radius = 1.0;
    static constexpr double max_flow_velocity = 1.0;
    static constexpr double peclet = 200.0;
    static constexpr double damkohler = 0.5;

    double photoreactor_starting_point = 0.0;
    double photoreactor_window = 2.0;
    double window_bound_sharpness = 20.0;
    double photoreaction_coefficient = (damkohler * max_flow_velocity) / photoreactor_window;

    class ParabolicDistributionInlet : public Function<dim>
    {
    public:
        virtual double value(const dealii::Point<dim> &p,
                             const unsigned int = 0) const override
        {
            double c_in = 1.0;
            double r_squared = p[1] * p[1] + p[2] * p[2];
            return c_in  * (1.0 - r_squared / (tube_radius * tube_radius));
        }
    };

    DEAL_II_HOST_DEVICE TubeSpec() {}

    DEAL_II_HOST_DEVICE Tensor advection(const Point &p) const override
    {
        Scalar r_squared = p[1] * p[1] + p[2] * p[2];
        Scalar poiseuille_velocity = max_flow_velocity *
            (1.0 - r_squared / (tube_radius * tube_radius));

        Tensor res;
        res[0] = poiseuille_velocity;
        res[1] = 0.0;
        res[2] = 0.0;
        return res;
    }

    DEAL_II_HOST_DEVICE Scalar diffusion(const Point &p) const override
    {
        return 2.0 * tube_radius * max_flow_velocity / peclet;
    }

    DEAL_II_HOST_DEVICE Scalar reaction(const Point &p) const override
    {
        Scalar window = 0.5 * (std::tanh(window_bound_sharpness * (p[0] - photoreactor_starting_point)) -
                               std::tanh(window_bound_sharpness * (p[0] - photoreactor_starting_point - photoreactor_window)));

        return 0.0;// photoreaction_coefficient * window;
    }

    DEAL_II_HOST_DEVICE Scalar forcing(const Point &p) const override
    {
        return 0.0;
    }

    void fill_dirichlet_map(
        BoundaryFuncMap<Function<dim>> &function_map) const override
    {
        function_map[0] = std::make_unique<Functions::ZeroFunction<dim>>();
        function_map[1] = std::make_unique<ParabolicDistributionInlet>();
    }

    void fill_neumann_map(BoundaryFuncMap<NeumannFunction<Scalar, Point>>
                          &function_map) const override
    {
        // Homogeneous Neumann on outlet
    }
};

#endif
