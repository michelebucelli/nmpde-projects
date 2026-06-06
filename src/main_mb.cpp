// The purpose of this file is to port the same computations
// of `main.cpp` with a matrix based backend

#include <deal.II/base/convergence_table.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/numerics/vector_tools.h>

#include "ADRFactory.hpp"
#include "ConcreteSpec.hpp"

using namespace dealii;

int main(int argc, char *argv[])
{
    Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

    auto exact_solution = std::make_shared<ExactSolution>();
    ConvergenceTable table;
    const std::vector<unsigned int> Ns = { 10, 20, 40, 80 };

    for (const unsigned int N : Ns) {

        parallel::distributed::Triangulation<dim> mesh(MPI_COMM_WORLD);
        GridGenerator::subdivided_hyper_cube(mesh, N, 0.0, 1.0, true);

        auto mb_adr = ADRFactory::matrix_based<
            dim, double, ConcreteSpec>(mesh, fe_degree);

        mb_adr.setup_system();
        mb_adr.assemble_system();
        mb_adr.solve(5000, 1e-12, true);

        const double error_L2 = mb_adr.compute_error(
            VectorTools::L2_norm, *exact_solution, mesh);

        const double error_H1 = mb_adr.compute_error(
            VectorTools::H1_norm, *exact_solution, mesh);

        table.add_value("h", 1.0 / N);
        table.add_value("eL2", error_L2);
        table.add_value("eH1", error_H1);
    }

    table.evaluate_all_convergence_rates(
        ConvergenceTable::reduction_rate_log2);

    table.set_scientific("eL2", true);
    table.set_scientific("eH1", true);
    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) {
        table.write_text(std::cout);
    }

    return 0;
}
