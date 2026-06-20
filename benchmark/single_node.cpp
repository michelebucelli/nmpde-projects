// The purpose of this file is to directly compare the matrix based and
// matrix free implementations down to the most primitive layer: the
// time per iteration of a GMRES solver
//
// The classes `MatrixFreeADR` and `MatrixBasedADR` are compared against
// an advection-diffusion-reaction problem on a mesh which could be either
// structured or unstructured
//
// By performing a fixed number of GMRES iterations without any preconditioning
// we get the raw time per iteration (per dof) metric. This metric is relevant
// for general runs as one could use a matrix based preconditioner while computing
// the action of A in a matrix free fashion
//
// For simplicity we restrict on using one node (one MPI rank) because it's
// easier to collect hardware performance counters. Other files will implement
// scalability tests on multiple nodes
//
// The acquisition campaign is challenging because of the many dimensions of
// the exploration space. We want to sweep:
//   * dofs: as the number of dofs on a single node is connected to the behaviour
//     of cache memories. We should see sudden drops in performance right
//     when the data can't fit in cache
//
//   * fe_degree: as it's connected to the arithmetic intensity of the matvec
//     operation. Matrix free methods should become more efficient for higher
//     order polynomials
//
//   * mesh type: structured meshes have a predictable access patterns, the matrix
//     free implementation could not be so dominant in the case of unstructured meshes
//
// We won't explore the cartesian product of those decision variables, but in principle
// we could. If we make this program parametric on the triplet (dofs, fe_degree, mesh type)
// we can use external scripts to sweep the parameter space.
// The problem is that matrix free implementations need to be templated on the fe_degree.
// To avoid recompilation, we could specialize the template for all degrees in [1..20].
// It will inflate the compilation time, but the cost is only paid once
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_generator.h>

#include <iostream>
#include <optional>

#include "ADRFactory.hpp"
#include "ConcreteSpec.hpp"


using namespace dealii;


enum class MeshType {
	Structured,
	Unstructured
};


enum class SolverType {
	MatrixBased,
	MatrixFree
};


struct RunParameters {
	int N;
	int degree;
	MeshType mesh_type;
	SolverType solver_type;
};


std::optional<RunParameters>
parse_arguments(int argc, char* argv[]) {
	MeshType mesh_type;
	SolverType solver_type;

	int N = -1;
	std::string mesh_name;

	if (argc != 4) return {};

	try {
		N = std::stoi(argv[1]);
		mesh_type = MeshType::Structured;
	} catch(std::exception const &) {
		mesh_type = MeshType::Unstructured;
	}

	const int degree = std::atoi(argv[2]);

	std::string solver_str {argv[3]};
	if (solver_str == "mf") {
		solver_type = SolverType::MatrixFree;
	} else if (solver_str == "mb") {
		solver_type = SolverType::MatrixBased;
	} else {
		std::cerr << "Error when defining solver type\n";
		return std::nullopt;
	}

	return RunParameters{N, degree, mesh_type, solver_type};
}


template <typename ADRSolver>
double run_benchmark(ADRSolver &adr_solver) {
	adr_solver.setup();
	adr_solver.assemble();

	const int maxit = 100;

	Timer timer;
	adr_solver.solve(maxit, 0.0, false);
	timer.stop();

	return timer.wall_time() / maxit;
}


constexpr int simd_width = 1;
constexpr int fe_degree = 2;


int main(int argc, char* argv[]) {
	Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

	const auto run_parameters = parse_arguments(argc, argv);

	if (not run_parameters.has_value()) {
		std::cerr
			<< "Need exactly 3 arguments:\n"
			<< "  * N (int or str)\n"
			<< "  * degree (int)\n"
			<< "  * implementation (str) which could be either `mb` or `mf`\n"
			<< std::endl;

		return 1;
	}

	RunParameters params = *run_parameters;
	parallel::distributed::Triangulation<dim> mesh(MPI_COMM_WORLD);

	double time_for_it;
	int ndofs;
	std::string mesh_str;

	if (params.mesh_type == MeshType::Structured) {
		GridGenerator::subdivided_hyper_cube(mesh, params.N, 0.0, 1.0, true);
		mesh_str = "structured";
	} else {
		const std::string mesh_filename = argv[1];

		GridIn<dim> gridin;
		gridin.attach_triangulation(mesh);
		std::ifstream f(mesh_filename);

		gridin.read_msh(f);

		mesh_str = "unstructured";
	}


	// I'm afraid that the benchmarking section requires some repetition
	if (params.solver_type == SolverType::MatrixBased) {
		auto adr_mb = ADRFactory::matrix_based<
			dim, double, ConcreteSpec>(mesh, params.degree);

		time_for_it = run_benchmark(adr_mb);
		ndofs = adr_mb.n_dofs();

	} else {
		// TODO: fe_degree is still harcoded, we need to come up with a solution
		auto adr_mf = ADRFactory::matrix_free<
			dim, fe_degree, simd_width, double, ConcreteSpec>(mesh);

		time_for_it = run_benchmark(adr_mf);
		ndofs = adr_mf.n_dofs();
	}

	std::cout << "ncells,ndofs,fe_degree,time_for_it,mesh_str" << std::endl;
	std::cout
		<< mesh.n_cells() << ","
		<< ndofs          << ","
		<< params.degree  << ","
		<< time_for_it    << ","
		<< mesh_str       << std::endl;

	return 0;
}
