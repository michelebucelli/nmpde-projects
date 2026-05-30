#include "MatrixFreeADR.hpp"
#include <iostream>
#include <deal.II/base/convergence_table.h>


using namespace dealii;


constexpr unsigned int dim = 2;

// Highly recommended to set fe_degree to at least 2 or 3 to see the true 
// performance benefit of Matrix-Free over SparseMatrix assemblies.
constexpr int fe_degree = 2;

constexpr unsigned int simd_width = 1;
using number = double;
using VectorizedArrayType = VectorizedArray<number, simd_width>;
using ArrayPoint = Point<dim, VectorizedArrayType>;
using MatrixFreeType = MatrixFree<dim, number, VectorizedArrayType>;


class ExactSolution : public Function<dim>
{
public:
    virtual double value(const Point<dim> &p, const unsigned int = 0) const override
    {
        return sin(2 * M_PI * p[0]) * cos(2 * M_PI * p[1]);
    }

    virtual Tensor<1, dim> gradient(const Point<dim> &p, const unsigned int = 0) const override
    {
        Tensor<1, dim> res;
        res[0] = 2 * M_PI * cos(2 * M_PI * p[0]) * cos(2 * M_PI * p[1]);
        res[1] = -2 * M_PI * sin(2 * M_PI * p[0]) * sin(2 * M_PI * p[1]);
        return res;
    }
};

int main(int argc, char *argv[])
{
    Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

    // Define lambda functors. By passing these directly as template arguments,
    // the compiler will completely inline them, ensuring AVX2 vectorization.
    auto beta = [](const ArrayPoint &p) {
        Tensor<1, dim, VectorizedArrayType> res;
        res[0] = p[1];
        res[1] = p[0];
        return res;
    };

    auto mu = [](const ArrayPoint &p) { return p[0]; };
    auto gamma = [](const ArrayPoint &p) { return p[1]; };

    auto f = [](const ArrayPoint &p) {
        return 8 * M_PI * M_PI * p[0] * sin(2 * M_PI * p[0]) * cos(2 * M_PI * p[1]) 
               - 2 * M_PI * cos(2 * M_PI * p[0]) * cos(2 * M_PI * p[1]) 
               + p[1] * sin(2 * M_PI * p[0]) * cos(2 * M_PI * p[1]) 
               + 2 * M_PI * p[1] * cos(2 * M_PI * p[0]) * cos(2 * M_PI * p[1]) 
               - 2 * M_PI * p[0] * sin(2 * M_PI * p[0]) * sin(2 * M_PI * p[1]);
    };

    ScalarFunc<dim, double> n_bc = [](const Point<dim> &) { return 0.0; };

    auto exact_solution = std::make_shared<ExactSolution>();
    ConvergenceTable table;
    const std::vector<unsigned int> Ns = { 10, 20, 40, 80 };

    for (const unsigned int N : Ns) {
        // Instantiate using decltype to push the lambda types to the templates
        MatrixFreeADR<dim, fe_degree, number, simd_width, decltype(beta), decltype(mu), decltype(gamma), decltype(f)>
            mf_adr(N, beta, mu, gamma, f);

        mf_adr.set_dirichlet_bc(0, exact_solution);
        mf_adr.set_dirichlet_bc(1, exact_solution);
        mf_adr.set_neumann_bc(2, n_bc);
        mf_adr.set_neumann_bc(3, n_bc);

        mf_adr.setup();
        mf_adr.assemble_rhs();
        mf_adr.solve(5000, 1e-12, true);
        mf_adr.output();

        const double error_L2 = mf_adr.compute_error(VectorTools::L2_norm, *exact_solution);
        const double error_H1 = mf_adr.compute_error(VectorTools::H1_norm, *exact_solution);

        table.add_value("h", 1.0 / N);
        table.add_value("eL2", error_L2);
        table.add_value("eH1", error_H1);
    }

    table.evaluate_all_convergence_rates(ConvergenceTable::reduction_rate_log2);
    table.set_scientific("eL2", true);
    table.set_scientific("eH1", true);

    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) {
        table.write_text(std::cout);
    }

    return 0;
}
