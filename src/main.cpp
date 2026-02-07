#include "MatrixFreeADR.hpp"

int main(int argc, char *argv[]) {

    const unsigned int dim = 2;
    const unsigned int fe_degree = 1;
    const unsigned int N = 20;

    Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

    MatrixFreeADR<dim, fe_degree> mf_adr(N);

    mf_adr.setup();
    mf_adr.assemble_rhs();
    mf_adr.solve();
    mf_adr.output();

    return 0;
}
