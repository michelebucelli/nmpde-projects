#include "ADROperator.hpp"
#include "MatrixFreeADR.hpp"

#include <iostream>

#include <deal.II/base/convergence_table.h>

static const unsigned int dim = 2;

class ExactSolution : public Function<dim>
{
public:
    ExactSolution() {}

    virtual double value(const Point<dim> &p,
                         const unsigned int = 0) const override
    {
        return p[0] + p[1];//sin(2 * M_PI * p[0]) * cos(2 * M_PI * p[1]);
    }

    virtual Tensor<1, dim> gradient(const Point<dim> &p,
                                    const unsigned int = 0) const override
    {
        Tensor<1, dim> res;

        res[0] = 2 * M_PI * cos(2 * M_PI * p[0]) * cos(2 * M_PI * p[1]);
        res[1] = -2 * M_PI * sin(2 * M_PI * p[0]) * sin(2 * M_PI * p[1]);

        return res;
    }
};

/* TODO: Figure out how to enable AVX2 when compiling. */

int main(int argc, char *argv[])
{
    Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

    const unsigned int n_vec_doubles = VectorizedArray<double>::size();
    std::cout << "n_vec_doubles=" << n_vec_doubles << std::endl;

    const unsigned int fe_degree = 1;

    VectorizedFunction<dim, double> mu =
        [](const VectorizedPoint<dim, double> &p) { return 1.0; };

    VectorizedFunction<dim, double> f =
        [](const VectorizedPoint<dim, double> &p)
    {
        //return 8 * M_PI * M_PI * p[0] * sin(2 * M_PI * p[0])
        //                              * cos(2 * M_PI * p[1]) -
        //       2 * M_PI * cos(2 * M_PI * p[0]) * cos(2 * M_PI * p[1]);

        return -5;
    };

    VectorizedFunction<dim, double> n_bc =
        [](const VectorizedPoint<dim, double> &p){ return p[1]; };

    auto exact_solution = std::make_shared<ExactSolution>();

    ConvergenceTable table;
    const std::vector<unsigned int> Ns = { 10, 20, 40, 80 };

    for (const unsigned int N : Ns) {
        MatrixFreeADR<dim, fe_degree> mf_adr(N, mu, f);

        mf_adr.set_dirichlet_bc(0, exact_solution);
        mf_adr.set_dirichlet_bc(1, exact_solution);

        mf_adr.set_neumann_bc(2, n_bc);
        mf_adr.set_neumann_bc(3, n_bc);

        mf_adr.setup();
        mf_adr.assemble_rhs();
        mf_adr.solve();
        mf_adr.output();

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
