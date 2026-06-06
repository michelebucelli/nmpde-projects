#include <deal.II/base/convergence_table.h>

#include <iostream>

#include "ADRFactory.hpp"
#include "ConcreteSpec.hpp"

using namespace dealii;

constexpr unsigned int simd_width = 1;

int main(int argc, char *argv[])
{
    Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

    auto exact_solution = std::make_shared<ExactSolution>();
    ConvergenceTable table;
    const std::vector<unsigned int> Ns = { 10, 20, 40, 80 };

    for (const unsigned int N : Ns) {

        // ConcreteSpec could be passed as template template parameter,
        // but to reduce the verbosity of the MatrixFreeADR template class
        // signature, we explicitly instantiate the template before passing
        // it to MatrixFreeADR. This step is hidden from the final user with
        // the MatrixFreeADRFactory class, that accept the ConcreteSpec as
        // template template argument and returns the MatrixFreeADR instance.
        // Internally, MatrixFreeADR can query ConcreteSpec to obtain
        // its template arguments (e.g. ConcreteSpec::ScalarType).

        auto mf_adr = ADRFactory::matrix_free<
            dim, fe_degree, simd_width, double, ConcreteSpec>(N);

        mf_adr.setup();
        mf_adr.assemble_rhs();
        mf_adr.solve(5000, 1e-12, true);
        mf_adr.output();

        const double error_L2 = mf_adr.compute_error(VectorTools::L2_norm,
                                                     *exact_solution);
        const double error_H1 = mf_adr.compute_error(VectorTools::H1_norm,
                                                     *exact_solution);

        table.add_value("h", 1.0 / N);
        table.add_value("eL2", error_L2);
        table.add_value("eH1", error_H1);
    }

    table.evaluate_all_convergence_rates(ConvergenceTable::
                                         reduction_rate_log2);
    table.set_scientific("eL2", true);
    table.set_scientific("eH1", true);

    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) {
        table.write_text(std::cout);
    }

    return 0;
}
