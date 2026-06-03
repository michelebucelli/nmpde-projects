#ifndef ADR_SPEC_HPP
#define ADR_SPEC_HPP


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
//     I think we can template on the return type
//
//   * how to pass functions to maps: we need to be careful about memory leaks. Since
//     the AdrSpec is stateless, we can't pass references to Function because they will
//     die after the `fill_.._map` call.
//
//     For the same reason we can't pass just pointers, because the memory will leak.
//     We must pass smart pointers, but then we will need to clone the map for using
//     `VectorTools::interpolate_boundary_values`. It's not a big deal
template <int dim>
struct AdrSpec {
	AdrSpec();

	double         mu    (Point<dim>& p) const;
	double         gamma (Point<dim>& p) const;
	Tensor<1, dim> beta  (Point<dim>& p) const;

	void fill_dirichlet_map(std::map<types::boundary_id, const Function<dim> *> &function_map) const;
	void fill_neumann_map  (std::map<types::boundary_id, const Function<dim> *> &function_map) const;
};


#endif // ADR_SPEC_HPP
