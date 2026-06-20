#include <deal.II/base/convergence_table.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>

#include <memory>
#include <iostream>

#include "ADRFactory.hpp"
#include "ConcreteSpec.hpp"

using namespace dealii;

constexpr unsigned int simd_width = 1;
constexpr unsigned int fe_degree = 2;

template<typename ADRSolver>
void run_convergence_sample(
    ADRSolver &adr_solver,
    unsigned int N,
    const ExactSolution &exact_solution,
    ConvergenceTable &table)
{
    adr_solver.setup();
    adr_solver.assemble();
    adr_solver.solve(5000, 1e-12, true);

    const double error_L2 = adr_solver.compute_error(
        VectorTools::L2_norm, exact_solution);
    const double error_H1 = adr_solver.compute_error(
        VectorTools::H1_norm, exact_solution);

    table.add_value("h", 1.0 / N);
    table.add_value("eL2", error_L2);
    table.add_value("eH1", error_H1);
}

void log_convergence_table(ConvergenceTable &table,
                           const std::string &title)
{
    table.evaluate_all_convergence_rates(ConvergenceTable::
                                         reduction_rate_log2);
    table.set_scientific("eL2", true);
    table.set_scientific("eH1", true);

    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) {
        std::cout << title << std::endl;
        table.write_text(std::cout);
    }
}

int main(int argc, char *argv[])
{
    Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

    auto exact_solution = std::make_shared<ExactSolution>();

    ConvergenceTable table_mf;
    ConvergenceTable table_mb;

    const std::vector<unsigned int> Ns = { 10, 20, 40, 80 };

    for (const unsigned int N : Ns) {

        parallel::distributed::Triangulation<dim> mesh(MPI_COMM_WORLD);
        // TODO: Check matrix-free simplex mesh support
        GridGenerator::subdivided_hyper_cube(mesh, N, 0.0, 1.0, true);

        // ConcreteSpec could be passed as template template parameter
        // to MatrixFreeADR constructor, but to reduce the verbosity of the
        // MatrixFreeADR template signature, we explicitly instantiate
        // the ConcreteSpec template before passing it to MatrixFreeADR.
        // This step is hidden from the final user with the
        // MatrixFreeADRFactory class, that accept the ConcreteSpec as
        // template template argument and returns the MatrixFreeADR instance.
        // Internally, MatrixFreeADR can query ConcreteSpec to obtain
        // its template arguments (e.g. ConcreteSpec::ScalarType).
        // The same can be said for MatrixBasedADR.

        auto adr_mf = ADRFactory::matrix_free<
            dim, fe_degree, simd_width, double, ConcreteSpec>(mesh);

        auto adr_mb = ADRFactory::matrix_based<
            dim, double, ConcreteSpec>(mesh, fe_degree);

        run_convergence_sample(adr_mf, N, *exact_solution, table_mf);
        run_convergence_sample(adr_mb, N, *exact_solution, table_mb);
    }

    log_convergence_table(table_mf, "\nmatrix-free:");
    log_convergence_table(table_mb, "\nmatrix-based:");

    return 0;
}
