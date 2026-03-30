#include "MatrixFreeADR.hpp"

#include <iostream>

static const unsigned int dim = 2;

/* TODO: Figure out how to enable AVX2 when compiling. */

int main(int argc, char *argv[])
{
    Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

    const unsigned int fe_degree = 1;
    const unsigned int N = 40;

    MatrixFreeADR<dim, fe_degree> mf_adr(N);

    mf_adr.setup();
    mf_adr.assemble_rhs();
    mf_adr.solve();
    mf_adr.output();

    return 0;
}
