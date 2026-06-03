// The purpose of this file is to port the same computations
// of `main.cpp` with a matrix based backend

#include "AdrSpec.hpp"
#include "MatrixBasedAdr.hpp"


constexpr unsigned int dim = 2;


class ExactSolution : public Function<dim> {
	public:
		double value(const Point<dim> &p, const unsigned int = 0) const override
		{
			return sin(2 * M_PI * p[0]) * cos(2 * M_PI * p[1]);
		}

		Tensor<1, dim> gradient(const Point<dim> &p, const unsigned int = 0) const override
		{
			Tensor<1, dim> res;
			res[0] = 2 * M_PI * cos(2 * M_PI * p[0]) * cos(2 * M_PI * p[1]);
			res[1] = -2 * M_PI * sin(2 * M_PI * p[0]) * sin(2 * M_PI * p[1]);
			return res;
		}
};


struct ConcreteSpec {
	double mu(Point<dim>& p) const {
		return p[0];
	}

	double gamma(Point<dim>& p) const {
		return p[1];
	}

	Tensor<1,dim> beta(Point<dim>& p) const {
		Tensor<1,dim> res;
		res[0] = p[1];
		res[1] = p[0];

		return res;
	}

	double forcing(Point<dim>& p) const {
		return 8 * M_PI * M_PI * p[0] * sin(2 * M_PI * p[0]) * cos(2 * M_PI * p[1]) 
			- 2 * M_PI * cos(2 * M_PI * p[0]) * cos(2 * M_PI * p[1]) 
			+ p[1] * sin(2 * M_PI * p[0]) * cos(2 * M_PI * p[1]) 
			+ 2 * M_PI * p[1] * cos(2 * M_PI * p[0]) * cos(2 * M_PI * p[1]) 
			- 2 * M_PI * p[0] * sin(2 * M_PI * p[0]) * sin(2 * M_PI * p[1]);
	}

	void fill_dirichlet_map(std::map<types::boundary_id, std::unique_ptr<Function<dim>>> &function_map) const {
		function_map[0] = std::make_unique<ExactSolution>();
		function_map[1] = std::make_unique<ExactSolution>();
	}

	void fill_neumann_map(std::map<types::boundary_id, std::unique_ptr<Function<dim>>> &function_map) const {
		function_map[2] = std::make_unique<Functions::ZeroFunction<dim>>();
		function_map[3] = std::make_unique<Functions::ZeroFunction<dim>>();
	}
};


int main(void) {
	Triangulation<2> mesh;
	ConcreteSpec spec;
	return 0;
}
