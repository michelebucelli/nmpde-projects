#ifndef ADR_HPP
#define ADR_HPP

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/conditional_ostream.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <string>
#include <utility>
#include <iostream>

using namespace dealii;

// TODO: Consider having an additional non-templated interface
// at the top of the hierarchy

template<typename ADRSpec, typename number, typename Vector_>
class ADR
{
protected:
    static constexpr int dim = ADRSpec::dim;

    ADR(const parallel::distributed::Triangulation<dim> &mesh_,
        unsigned int fe_degree);

    const ADRSpec local_spec;

    // WARNING: If named fe_degree, it would shadow
    // MatrixFreeADR fe_degree template parameter!
    // Alternatively, we could have a virtual get_fe_degree()
    const unsigned int fe_degree_;

    const MappingQ1<dim> mapping;
    const FE_Q<dim> fe;

    const parallel::distributed::Triangulation<dim> &mesh;

    DoFHandler<dim> dof_handler;

    Vector_ solution;

    ConditionalOStream pcout;

public:
    virtual void setup() = 0;

    virtual void assemble() = 0;

    virtual std::pair<unsigned int, double> solve(unsigned int max_iter,
                                                  number rel_tol,
                                                  bool apply_preconditioner) = 0;

    unsigned int n_dofs() const;

    void output(const std::string &output_name,
                const std::string &output_dir = "./");

    double compute_error(const VectorTools::NormType &norm_type,
                         const Function<dim, number> &exact_solution) const;

    virtual ~ADR() = default;
    Vector_ system_rhs;
};

template<typename ADRSpec, typename number, typename Vector_>
ADR<ADRSpec, number, Vector_>::ADR(
    const parallel::distributed::Triangulation<dim> &mesh_,
    unsigned int fe_degree
)
    : fe_degree_(fe_degree)
    , fe(fe_degree)
    , mesh(mesh_)
    , dof_handler(mesh)
    , pcout(std::cout,
            Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
{}

template<typename ADRSpec, typename number, typename Vector_>
void ADR<ADRSpec, number, Vector_>::output(const std::string &output_name,
                                           const std::string &output_dir)
{
    DataOut<dim> data_out;

    // Not really needed, it's called inside solve() by MatrixFreeADR
    solution.update_ghost_values();

    data_out.add_data_vector(dof_handler, solution, "solution");
    data_out.build_patches(mapping);

    data_out.write_vtu_with_pvtu_record(output_dir,
                                        output_name,
                                        0,
                                        MPI_COMM_WORLD);
}

template<typename ADRSpec, typename number, typename Vector_>
double ADR<ADRSpec, number, Vector_>::compute_error(
    const VectorTools::NormType &norm_type,
    const Function<dim, number> &exact_solution) const
{
    Vector<double> error_per_cell(mesh.n_active_cells());

    VectorTools::integrate_difference(mapping,
                                      dof_handler,
                                      solution,
                                      exact_solution,
                                      error_per_cell,
                                      QGauss<dim>(fe_degree_ + 2),
                                      norm_type);

    const double error =
        VectorTools::compute_global_error(mesh, error_per_cell, norm_type);

    return error;
}

template<typename ADRSpec, typename number, typename Vector_>
unsigned int ADR<ADRSpec, number, Vector_>::n_dofs() const
{
    // WARNING: n_dofs() returns types::global_dof_index
    return dof_handler.n_dofs();
}

#endif
