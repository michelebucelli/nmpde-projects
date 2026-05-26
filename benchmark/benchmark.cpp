#include <MatrixFreeADR.hpp>
#include "bench_config.h"

#include <iostream>

#include <deal.II/base/convergence_table.h>


int main(int argc, char *argv[])
{
    Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

    const unsigned int n_vec_doubles = VectorizedArrayType::size();
    std::cout << "n_vec_doubles=" << n_vec_doubles << std::endl;

    constexpr unsigned int dim = BenchConfig::dim;

    std::cout << "Running with parameters d: " << dim << " p: " << BenchConfig::fe_degree << std::endl;

    VectorArrayFunction<dim> beta =
        [](const ArrayPoint<dim> &p) {
            Tensor<1, dim, VectorizedArrayType> res;
            // TODO: Test convergence using a coefficient with
            // non-const divergence
            res[0] = p[1];
            res[1] = p[0];
            return res;
        };

    ScalarArrayFunction<dim> mu =
        [](const ArrayPoint<dim> &p) { return p[0]; };

    ScalarArrayFunction<dim> gamma =
        [](const ArrayPoint<dim> &p) { return p[1]; };

    ScalarArrayFunction<dim> f =
        [](const ArrayPoint<dim> &p)
    {
               // Diffusion term
        return 8 * M_PI * M_PI * p[0] * sin(2 * M_PI * p[0])
                                      * cos(2 * M_PI * p[1]) -
               2 * M_PI * cos(2 * M_PI * p[0]) * cos(2 * M_PI * p[1]) +
               // Reaction term
               p[1] * sin(2 * M_PI * p[0]) * cos(2 * M_PI * p[1]) +
               // Advection term
               2 * M_PI * p[1] * cos(2 * M_PI * p[0]) * cos(2 * M_PI * p[1]) -
               2 * M_PI * p[0] * sin(2 * M_PI * p[0]) * sin(2 * M_PI * p[1]);
    };

    ScalarArrayFunction<dim> n_bc =
        [](const ArrayPoint<dim> &){ return 0; };

    auto exact_solution = std::make_shared<BenchConfig::ExactSolution>();

    ConvergenceTable table;
    const std::vector<unsigned int> Ns = { 10, 20, 40, 80 };

    for (const unsigned int N : Ns) {
        MatrixFreeADR<dim, BenchConfig::fe_degree, BenchConfig::RealNumber> mf_adr(N, beta, mu, gamma, f);

        mf_adr.set_dirichlet_bc(0, exact_solution);
        mf_adr.set_dirichlet_bc(1, exact_solution);
        //mf_adr.set_dirichlet_bc(2, exact_solution);
        //mf_adr.set_dirichlet_bc(3, exact_solution);

        mf_adr.set_neumann_bc(2, n_bc);
        mf_adr.set_neumann_bc(3, n_bc);

        mf_adr.setup();
        mf_adr.assemble_rhs();
        mf_adr.solve();

        const double error_L2 =
            mf_adr.compute_error(VectorTools::L2_norm, *exact_solution);
        const double error_H1 =
            mf_adr.compute_error(VectorTools::H1_norm, *exact_solution);

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