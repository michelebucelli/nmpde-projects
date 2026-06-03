// The purpose of this file is to port the same computations
// of `main.cpp` with a matrix based backend
#include "MatrixBasedAdr.hpp"
#include <deal.II/grid/grid_generator.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/base/convergence_table.h>
using namespace dealii;


constexpr unsigned int dim = 2;
constexpr unsigned int fe_degree = 2;


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
	double mu(const Point<dim>& p) const {
		return p[0];
	}

	double gamma(const Point<dim>& p) const {
		return p[1];
	}

	Tensor<1,dim> beta(const Point<dim>& p) const {
		Tensor<1,dim> res;
		res[0] = p[1];
		res[1] = p[0];

		return res;
	}

	double forcing(const Point<dim>& p) const {
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
	auto exact_solution = std::make_shared<ExactSolution>();
	ConvergenceTable table;
	const std::vector<unsigned int> Ns = { 10, 20, 40, 80 };

	for (const unsigned int N : Ns) {
		Triangulation<dim> mesh;
		GridGenerator::subdivided_hyper_cube(mesh, N);

		MatrixBasedAdr<dim, ConcreteSpec> mb_adr(mesh, fe_degree);

		mb_adr.setup_system();
		mb_adr.assemble_system();
		mb_adr.solve(5000, 1e-12, true);

		const double error_L2 = mb_adr.compute_error(VectorTools::L2_norm, *exact_solution, mesh);
		const double error_H1 = mb_adr.compute_error(VectorTools::H1_norm, *exact_solution, mesh);

		table.add_value("h", 1.0 / N);
		table.add_value("eL2", error_L2);
		table.add_value("eH1", error_H1);
	}

	table.evaluate_all_convergence_rates(ConvergenceTable::reduction_rate_log2);
	table.set_scientific("eL2", true);
	table.set_scientific("eH1", true);
	table.write_text(std::cout);

	return 0;
}
