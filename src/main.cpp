#include "MatrixFreeADR.hpp"

#include <deal.II/base/convergence_table.h>

static const unsigned int dim = 2;

class ExactSolution : public Function<dim>
{
public:
    ExactSolution() {}

    virtual double value(const Point<dim> &p,
                         const unsigned int = 0) const override
    {
        return sin(2 * M_PI * p[0]) * cos(2 * M_PI * p[1] - M_PI / 2);
    }

    virtual Tensor<1, dim> gradient(const Point<dim> &p,
                                    const unsigned int = 0) const override
    {
        Tensor<1, dim> res;

        res[0] = 2 * M_PI * cos(2 * M_PI * p[0])
                          * cos(2 * M_PI * p[1] - M_PI / 2);

        res[1] = -2 * M_PI * sin(2 * M_PI * p[0])
                           * sin(2 * M_PI * p[1] - M_PI / 2);

        return res;
    }
};

class ForcingTerm
{
public:
    ForcingTerm() {}

    DEAL_II_HOST_DEVICE double operator()(const Point<dim> &p) const
    {
        return 8 * M_PI * M_PI * sin(2 * M_PI * p[0])
                               * cos(2 * M_PI * p[1] - M_PI / 2);
    }
};

int main(int argc, char *argv[]) {

    Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

    const unsigned int fe_degree = 2;
    const std::vector<unsigned int> Ns = { 10, 20, 40, 80 };
    ExactSolution exact_solution;

    ConvergenceTable table;

    for (const unsigned int N : Ns) {
        MatrixFreeADR<dim, fe_degree, ForcingTerm> mf_adr(N);

        const ConditionalOStream &pcout = mf_adr.get_pcout();
        pcout << "dim=" << dim << " fe_degree=" << fe_degree
              << " N=" << N << std::endl;

#ifdef LOG_TIMERS
        Timer time;
#endif
        mf_adr.setup();
        mf_adr.assemble_rhs();
        mf_adr.solve();
#ifdef LOG_TIMERS
        pcout << std::setw(20) << std::left << "total" << std::fixed
              << time.wall_time() << "s" << std::endl << std::endl;
#endif

        mf_adr.output();

        const double error_L2 =
            mf_adr.compute_error(VectorTools::L2_norm, exact_solution);
        const double error_H1 =
            mf_adr.compute_error(VectorTools::H1_norm, exact_solution);

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
