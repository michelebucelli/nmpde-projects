#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/base/mpi.h>

#include <iostream>
#include <fstream>
#include <optional>

#include "ADRFactory.hpp"
#include "TubeSpec.hpp"

void cylindric_wall_adapted_mesh(parallel::distributed::Triangulation<dim> &mesh,
		double length,
		double radius,
		unsigned int length_subdivision,
		unsigned int global_refinement_steps,
		unsigned int wall_refinement_steps)
{
	GridGenerator::subdivided_cylinder(mesh, length_subdivision, radius, length / 2);
	//to have the inlet around the origin
	//Tensor<1, 3, double> shift_right{{length/2, 0.0, 0.0}};
	//GridTools::shift(shift_right, mesh);

	mesh.refine_global(global_refinement_steps);

	for (unsigned int step = 0; step < wall_refinement_steps; ++step)
	{
		for (const auto &cell : mesh.active_cell_iterators())
		{
			if (cell->at_boundary())
			{
				cell->set_refine_flag();
			}
		}
		mesh.execute_coarsening_and_refinement();
	}
}


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
	if (num_iters_parsed < 0) {
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
		unsigned int& num_iters) {
	RunTimers res;
	Timer timer(MPI_COMM_WORLD, true);

	timer.restart();
	adr_solver.setup();
	timer.stop();
	res.setup_time = timer.wall_time();

	timer.restart();
	adr_solver.assemble();
	timer.stop();
	res.assemble_time = timer.wall_time();

	if (num_iters != 0) {
		timer.restart();
		adr_solver.solve(num_iters, 0.0, false);
		timer.stop();
	} else {
		timer.restart();
		const auto [niters, _] = adr_solver.solve(1000, 1e-6, true);
		timer.stop();

		num_iters = niters;
	}

	res.solve_time = timer.wall_time();
	res.solver_iteration_time = res.solve_time / num_iters;

	return res;
}

constexpr int simd_width = VectorizedArray<double>::size();

int main(int argc, char* argv[])
{
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

	double radius = TubeSpec<>::tube_radius;
	double length = radius * 8.0;

	parallel::distributed::Triangulation<dim> mesh(MPI_COMM_WORLD);


	RunTimers run_timers;
	unsigned int ndofs = 0;
	std::string mesh_str;

	if (params.mesh_type == MeshType::Structured) {
		// GUARDARE QUI PER CONTROLLARE IL NUMERO DI DOFS
		cylindric_wall_adapted_mesh(mesh, length, radius, 10, params.N, 1);
		mesh_str = "structured";
	} else {
		return 1;
	}


	// I'm afraid that the benchmarking section requires some repetition
	if (params.solver_type == SolverType::MatrixBased) {
#ifndef BUILD_WITH_KOKKOS
		auto adr_mb = ADRFactory::matrix_based<
			dim, double, TubeSpec>(mesh, params.degree);

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
				dim, degree, double, TubeSpec, simd_width>(mesh); \
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
