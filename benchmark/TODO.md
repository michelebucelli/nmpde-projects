## CPU Benchmark TODOs

Read this https://dealii.org/workshop-2020/slides-kronbichler-dealii-2020.pdf to get an idea.

- [ ] benchmark single node sweeps
- [ ] benchmark matrix_free_data_locality()
- [ ] benchmark different mf threading schemes
- [ ] measure memory usage
- [ ] MPI overhead with mpiP?

For these we should restrict the profiling to the SolverGMRES::solve() call
- [ ] measure cache behavior
- [ ] measure instruction throughput -> against peak FLOPS
- [ ] measure tma metrics
- [ ] compute roofline

