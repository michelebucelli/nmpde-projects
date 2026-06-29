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
#include <deal.II/base/mpi.h>

#include <iostream>
#include <fstream>
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
	unsigned int N;
	unsigned int degree;
	MeshType mesh_type;
	SolverType solver_type;
	unsigned int num_iters;
};


std::optional<RunParameters>
parse_arguments(int argc, char* argv[]) {
	MeshType mesh_type;
	SolverType solver_type;

	unsigned int N = 0;
	std::string mesh_name;

	if (argc != 5) return {};

	try {
		int N_parsed = std::stoi(argv[1]);
		if (N_parsed <= 0) {
			std::cerr << "Structured mesh division N must be positive\n";
			return std::nullopt;
		}
		N = static_cast<unsigned int>(N_parsed);
		mesh_type = MeshType::Structured;
	} catch(std::exception const &) {
		mesh_type = MeshType::Unstructured;
	}

	int degree_parsed = std::atoi(argv[2]);
	if (degree_parsed <= 0) {
		std::cerr << "Degree must be positive\n";
		return std::nullopt;
	}
	const unsigned int degree = static_cast<unsigned int>(degree_parsed);

	std::string solver_str {argv[3]};
	if (solver_str == "mf") {
		solver_type = SolverType::MatrixFree;
	} else if (solver_str == "mb") {
		solver_type = SolverType::MatrixBased;
	} else {
		std::cerr << "Error when defining solver type\n";
		return std::nullopt;
	}

	int num_iters_parsed = -1;
	try {
		num_iters_parsed = std::stoi(argv[4]);
	} catch(std::exception const &) {
		std::cerr << "Error when parsing number of iterations\n";
		return std::nullopt;
	}
	if (num_iters_parsed <= 0) {
		std::cerr << "Number of iterations must be positive\n";
		return std::nullopt;
	}
	const unsigned int num_iters = static_cast<unsigned int>(num_iters_parsed);

	return RunParameters{N, degree, mesh_type, solver_type, num_iters};
}


struct RunTimers {
	double setup_time,
	       assemble_time,
	       solve_time,
	       solver_iteration_time;
};


template <typename ADRSolver>
RunTimers run_benchmark(ADRSolver &adr_solver,
                        unsigned int num_iters) {
	RunTimers res;
	Timer timer;

	MPI_Barrier(MPI_COMM_WORLD);
	timer.restart();
	adr_solver.setup();
	timer.stop();
	res.setup_time = Utilities::MPI::max(timer.wall_time(), MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
	timer.restart();
	adr_solver.assemble();
	timer.stop();
	res.assemble_time = Utilities::MPI::max(timer.wall_time(), MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
	timer.restart();
	const auto solve_result = adr_solver.solve(num_iters, 0.0, false);
	timer.stop();
	res.solve_time = Utilities::MPI::max(timer.wall_time(), MPI_COMM_WORLD);
	res.solver_iteration_time = Utilities::MPI::max(solve_result.second / num_iters, MPI_COMM_WORLD);

	return res;
}

constexpr int simd_width = VectorizedArray<double>::size();

int main(int argc, char* argv[]) {
	Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

	const auto run_parameters = parse_arguments(argc, argv);

	if (not run_parameters.has_value()) {
		std::cerr
			<< "Need exactly 4 arguments:\n"
			<< "  * N (int or str)\n"
			<< "  * degree (int)\n"
			<< "  * implementation (str) which could be either `mb` or `mf`\n"
			<< "  * iterations (int)\n"
			<< std::endl;

		return 1;
	}

	RunParameters params = *run_parameters;
	parallel::distributed::Triangulation<dim> mesh(MPI_COMM_WORLD);

	RunTimers run_timers;
	unsigned int ndofs = 0;
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
#ifndef BUILD_WITH_KOKKOS
		auto adr_mb = ADRFactory::matrix_based<
			dim, double, ConcreteSpec>(mesh, params.degree);

		run_timers = run_benchmark(adr_mb, params.num_iters);
		ndofs = adr_mb.n_dofs();
#else
        std::cout << "WARNING: matrix-based solver not currently "
                  << "available when BUILD_WITH_KOKKOS=ON" << std::endl;
        return 1;
#endif

	} else {
		#define BENCH_MF(degree) do {                                         \
			auto adr_mf = ADRFactory::matrix_free<                        \
				dim, degree, double, ConcreteSpec, simd_width>(mesh); \
                                                                                      \
		        run_timers = run_benchmark(adr_mf, params.num_iters);         \
			ndofs = adr_mf.n_dofs();                                      \
		} while (0)

		switch(params.degree) {
			case 1: BENCH_MF(1); break;
			case 2: BENCH_MF(2); break;
			case 3: BENCH_MF(3); break;
			case 4: BENCH_MF(4); break;
			case 5: BENCH_MF(5); break;
			case 6: BENCH_MF(6); break;
			case 7: BENCH_MF(7); break;
			case 8: BENCH_MF(8); break;
			case 9: BENCH_MF(9); break;
			case 10: BENCH_MF(10); break;
			case 11: BENCH_MF(11); break;
			case 12: BENCH_MF(12); break;
			case 13: BENCH_MF(13); break;
			case 14: BENCH_MF(14); break;
			default:
				 return 1;
		}
	}

	if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) {
		std::cout << "ncells,ndofs,fe_degree,iterations,setup_time,assemble_time,solve_time,solver_iteration_time,solver_iter_time_per_dof,mesh_str" << std::endl;
		std::cout
			<< mesh.n_cells()            << ","
			<< ndofs                     << ","
			<< params.degree             << ","
			<< params.num_iters          << ","
			<< run_timers.setup_time     << ","
			<< run_timers.assemble_time  << ","
			<< run_timers.solve_time     << ","
			<< run_timers.solver_iteration_time << ","
			<< run_timers.solver_iteration_time / ndofs << ","
			<< mesh_str
			<< std::endl;
	}

	return 0;
}
