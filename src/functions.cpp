#include "functions.h"

using VectorType = VectorizedArray<double>;

template <int dim>
double MuFunction<dim>::value(const Point<dim> & p, const unsigned int component) const
{
    return value<double>(p, component);
}

template <int dim>
template <typename number>
number MuFunction<dim>::value(const Point<dim, number> & /*p*/, const unsigned int /*component*/) const
{
    return 0.1; // -∇·(µ∇u)
}

template <int dim>
Tensor<1, dim, double> BetaFunction<dim>::gradient(const Point<dim> & p, const unsigned int /*component*/) const
{
    return gradient(p);
}

template <int dim>
template <typename number>
Tensor<1, dim, number> BetaFunction<dim>::gradient(const Point<dim, number> & /*p*/, const unsigned int /*component*/) const
{
    Tensor<1, dim, number> beta; // ∇·(βu)

    beta[0] = 0.0;
    beta[1] = 0.0;
    if constexpr (dim == 3)
        beta[2] = 0.0;

    return beta;
}

template <int dim>
double GammaFunction<dim>::value(const Point<dim> & p, const unsigned int component) const
{
    return value<double>(p, component);
}

template <int dim>
template <typename number>
number GammaFunction<dim>::value(const Point<dim, number> & /*p*/, const unsigned int /*component*/) const
{
    return 0.1; // γu
}

template <int dim>
double SourceFunction<dim>::value(const Point<dim> & p, const unsigned int component) const
{
    return value<double>(p, component);
}

template <int dim>
template <typename number>
number SourceFunction<dim>::value(const Point<dim, number> & p, const unsigned int /*component*/) const
{
    return p[0] * p[1]; // f
}

template <int dim>
double GFunction<dim>::value(const Point<dim> & /*p*/, const unsigned int /*component*/) const
{
    return 0.1; // g
}

template <int dim>
double HFunction<dim>::value(const Point<dim> & p, const unsigned int component) const
{
    return value<double>(p, component);
}

template <int dim>
template <typename number>
number HFunction<dim>::value(const Point<dim, number> & /*p*/, const unsigned int /*component*/) const
{
    return 0.1; // h
}

template class MuFunction<2>;
template class MuFunction<3>;
template class BetaFunction<2>;
template class BetaFunction<3>;
template class GammaFunction<2>;
template class GammaFunction<3>;
template class SourceFunction<2>;
template class SourceFunction<3>;
template class GFunction<2>;
template class GFunction<3>;
template class HFunction<2>;
template class HFunction<3>;

template VectorType MuFunction<2>::value(const dealii::Point<2, VectorType> &, const unsigned int) const;
template VectorType MuFunction<3>::value(const dealii::Point<3, VectorType> &, const unsigned int) const;

template Tensor<1, 2, VectorType> BetaFunction<2>::gradient(const dealii::Point<2, VectorType> &, const unsigned int) const;
template Tensor<1, 3, VectorType> BetaFunction<3>::gradient(const dealii::Point<3, VectorType> &, const unsigned int) const;

template VectorType GammaFunction<2>::value(const dealii::Point<2, VectorType> &, const unsigned int) const;
template VectorType GammaFunction<3>::value(const dealii::Point<3, VectorType> &, const unsigned int) const;

template VectorType SourceFunction<2>::value(const dealii::Point<2, VectorType> &, const unsigned int) const;
template VectorType SourceFunction<3>::value(const dealii::Point<3, VectorType> &, const unsigned int) const;

template VectorType HFunction<2>::value(const dealii::Point<2, VectorType> &, const unsigned int) const;
template VectorType HFunction<3>::value(const dealii::Point<3, VectorType> &, const unsigned int) const;
