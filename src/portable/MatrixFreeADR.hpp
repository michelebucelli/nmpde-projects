#ifndef PORTABLE_MATRIX_FREE_ADR
#define PORTABLE_MATRIX_FREE_ADR

#include <deal.II/base/timer.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/conditional_ostream.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>
#include <fstream>
#include <functional>

#include "../ADR.hpp"
#include "../ADRSpec.hpp"
#include "ADROperator.hpp"

using namespace dealii;

template<unsigned int fe_degree, typename ADRSpec>
class MatrixFreeADR :
    public ADR<ADRSpec,
               typename ADRSpec::ScalarType,
               DistributedDeviceVector<typename ADRSpec::ScalarType>>
{
    using number = typename ADRSpec::ScalarType;
    using ADR<ADRSpec, number, DistributedDeviceVector<number>>::dim;

    friend class ADRFactory;

    using PointType = typename ADRSpec::PointType;

public:
    MatrixFreeADR(const parallel::distributed::Triangulation<dim> &mesh);

    void setup() override;

    void assemble() override;

    std::pair<unsigned int, double> solve(unsigned int max_iter,
                                          number rel_tol,
                                          bool apply_preconditioner) override;

    void output(const std::string &output_name,
                const std::string &output_dir = "./");

    double compute_error(const VectorTools::NormType &norm_type,
                         const Function<dim> &exact_solution) const;

    template <class Vector1_, class Vector2_>
    void vmult(Vector1_& dest, const Vector2_& src);

private:
    std::unique_ptr<ADROperator<fe_degree, ADRSpec>> system_matrix_dev;

    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;

    // TODO: Make sure DEAL_II_MPI_WITH_DEVICE_SUPPORT is defined
    DistributedDeviceVector<number> solution_update_dev;
    DistributedVector<number> ghost_solution_host;
};

template<unsigned int fe_degree, typename ADRSpec>
MatrixFreeADR<fe_degree, ADRSpec>::MatrixFreeADR(
    const parallel::distributed::Triangulation<dim> &mesh
)
    : ADR<ADRSpec, number, DistributedDeviceVector<number>>(mesh, fe_degree)
{}

template<unsigned int fe_degree, typename ADRSpec>
void MatrixFreeADR<fe_degree, ADRSpec>::setup()
{
#ifdef LOG_TIMERS
    Timer time;
#endif

    // TODO: Check DoFRenumbering::matrix_free_data_locality()
    this->dof_handler.distribute_dofs(this->fe);

    locally_owned_dofs = this->dof_handler.locally_owned_dofs();
    locally_relevant_dofs =
        DoFTools::extract_locally_relevant_dofs(this->dof_handler);

    this->constraints.clear();
    this->constraints.reinit(locally_owned_dofs, locally_relevant_dofs);

    DoFTools::make_hanging_node_constraints(this->dof_handler,
                                            this->constraints);

    BoundaryFuncMap<Function<dim>> tmp_dirichlet_map;
    std::map<types::boundary_id, const Function<dim> *> dirichlet_map;
    this->local_spec.fill_dirichlet_map(tmp_dirichlet_map);
    for (const auto& [boundary_id, dirichlet_function] : tmp_dirichlet_map) {
        dirichlet_map[boundary_id] = dirichlet_function.get();
    }

    VectorTools::interpolate_boundary_values(this->mapping,
                                             this->dof_handler,
                                             dirichlet_map,
                                             this->constraints);

    this->constraints.close();

    system_matrix_dev.reset(new ADROperator<fe_degree, ADRSpec>(
        this->dof_handler, this->mapping, this->constraints));

    //system_matrix_dev->template precompute_coefs<ForcingTerm>(dof_handler);

    ghost_solution_host.reinit(locally_owned_dofs,
                               locally_relevant_dofs,
                               MPI_COMM_WORLD);

    system_matrix_dev->initialize_dof_vector(this->solution);
    solution_update_dev.reinit(this->solution);
    this->system_rhs.reinit(this->solution);

#ifdef LOG_TIMERS
    this->pcout << std::setw(20) << std::left << "setup()"
                << std::fixed << time.wall_time() << "s" << std::endl;
#endif
}

template<unsigned int fe_degree, typename ADRSpec>
void MatrixFreeADR<fe_degree, ADRSpec>::assemble()
{
#ifdef LOG_TIMERS
    Timer time;
#endif

    ghost_solution_host = 0;
    this->constraints.distribute(ghost_solution_host);
    ghost_solution_host.update_ghost_values();

    // TODO: Can you avoid this transfer? I guess only if there's a way
    // to interpolate constraints directly on the device
    LinearAlgebra::ReadWriteVector<number> rw_vector(locally_owned_dofs);
    // TODO: Does this two-copy transfer have to do with the usage of
    // non-pinned memory?
    rw_vector.import_elements(ghost_solution_host, VectorOperation::insert);
    this->solution.import_elements(rw_vector, VectorOperation::insert);

    system_matrix_dev->compute_rhs(this->system_rhs, this->solution);

    // TODO: Neumann bcs rhs contribution is currently computed on the host.
    // Check wether overlapping this computation with compute_rhs is possible.

    DistributedVector<number> rhs_neumann_term(locally_owned_dofs,
                                               locally_relevant_dofs,
                                               MPI_COMM_WORLD);

    BoundaryFuncMap<NeumannFunction<number, PointType>> neumann_map;
    this->local_spec.fill_neumann_map(neumann_map);

    rhs_neumann_term = 0;

    const QGauss<dim - 1> quadrature_face(fe_degree + 1);

    FEFaceValues<dim> fe_face_values(this->mapping,
                                     this->fe,
                                     quadrature_face,
                                     update_values | update_JxW_values |
                                     update_quadrature_points);

    const unsigned int dofs_per_cell = this->fe.n_dofs_per_cell();
    const unsigned int n_q_face = quadrature_face.size();

    Vector<number> cell_neumann_term(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // TODO: Manually use kokkos to dispatch this loop on the host
    for (const auto &cell : this->dof_handler.active_cell_iterators()) {
        if (cell->is_locally_owned() && cell->at_boundary()) {

            cell_neumann_term = 0;

            for (unsigned int face = 0; face < cell->n_faces(); ++face) {
                const types::boundary_id b_id =
                    cell->face(face)->boundary_id();

                if (cell->face(face)->at_boundary() &&
                    neumann_map.count(b_id)) {

                    const auto &bc_function = neumann_map[b_id];

                    fe_face_values.reinit(cell, face);
                    for (unsigned int q = 0; q < n_q_face; ++q) {
                        const PointType q_p =
                            fe_face_values.quadrature_point(q); 

                        const number mu_val = this->local_spec.diffusion(q_p);
                        const number h_val = bc_function->value(q_p);

                        for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                            cell_neumann_term(i) +=
                                mu_val *
                                h_val *
                                fe_face_values.shape_value(i, q) *
                                fe_face_values.JxW(q);
                        }
                    }
                }
            }

            cell->get_dof_indices(local_dof_indices);
            // TODO: We could avoid using constraints object, right?
            this->constraints.distribute_local_to_global(cell_neumann_term,
                                                         local_dof_indices,
                                                         rhs_neumann_term);
        }
    }

    rhs_neumann_term.compress(VectorOperation::add);

    rw_vector.import_elements(rhs_neumann_term, VectorOperation::insert);
    this->system_rhs.import_elements(rw_vector, VectorOperation::add);

#ifdef LOG_TIMERS
    this->pcout << std::setw(20) << std::left << "assemble()"
                << std::fixed << time.wall_time() << "s" << std::endl;
#endif
}

template<unsigned int fe_degree, typename ADRSpec>
std::pair<unsigned int, double> MatrixFreeADR<fe_degree, ADRSpec>::solve(
    unsigned int max_iter, number rel_tol, bool apply_preconditioner)
{
#ifdef LOG_TIMERS
    Timer time;
#endif

    SolverControl *solver_control;

    if (rel_tol == 0.0) {
        // same as in MatrixBasedADR::solve
        solver_control = new IterationNumberControl(max_iter, 0.0);
    } else {
        solver_control = new SolverControl(
            max_iter, rel_tol * this->system_rhs.l2_norm());
    }

    SolverGMRES<DistributedDeviceVector<number>> solver(*solver_control);

    auto run_solve = [&](const auto &preconditioner) {
        Timer solver_timer;
        solver.solve(*system_matrix_dev,
                     solution_update_dev,
                     this->system_rhs,
                     preconditioner);
        return solver_timer.wall_time();
    };

    double solver_time = 0.0;

    if (apply_preconditioner) {
        // NOTE: This could probably be moved in assemble or setup
        system_matrix_dev->compute_diagonal();

        const auto &preconditioner =
            *system_matrix_dev->get_matrix_diagonal_inverse();

        solver_time = run_solve(preconditioner);
    } else {
        PreconditionIdentity preconditioner;

        solver_time = run_solve(preconditioner);
    }

    this->solution += solution_update_dev;

    // Copy solution to host
    // TODO: Time d2h copy
    LinearAlgebra::ReadWriteVector<number> rw_vector(locally_owned_dofs);
    rw_vector.import_elements(this->solution, VectorOperation::insert);
    ghost_solution_host.import_elements(rw_vector, VectorOperation::insert);

    this->constraints.distribute(ghost_solution_host);

    ghost_solution_host.update_ghost_values();

#ifdef LOG_TIMERS
    const double elapsed = time.wall_time();
    const double elapsed_per_iter = elapsed / solver_control->last_step();

    this->pcout << std::setw(20) << std::left << "solve() [tot/iter]"
                << std::fixed << elapsed << "s/" << elapsed_per_iter
                << "s" << std::endl;
#endif

    return {solver_control->last_step(), solver_time};
}

template<unsigned int fe_degree, typename ADRSpec>
void MatrixFreeADR<fe_degree, ADRSpec>::
    output(const std::string &output_name, const std::string &output_dir)
{
    DataOut<dim> data_out;

    data_out.add_data_vector(
        this->dof_handler, ghost_solution_host, "solution");

    data_out.build_patches(this->mapping);

    data_out.write_vtu_with_pvtu_record(output_dir,
                                        output_name,
                                        0,
                                        MPI_COMM_WORLD);
}

template<unsigned int fe_degree, typename ADRSpec>
double MatrixFreeADR<fe_degree, ADRSpec>::
    compute_error(const VectorTools::NormType &norm_type,
                  const Function<dim> &exact_solution) const
{
    Vector<double> error_per_cell(this->mesh.n_active_cells());

    // TODO: Check if you can use the base class compute_error()
    // which computes the error on this->solution instead of
    // ghost_solution_host, same for output()
    VectorTools::integrate_difference(this->mapping,
                                      this->dof_handler,
                                      ghost_solution_host,
                                      exact_solution,
                                      error_per_cell,
                                      QGauss<dim>(fe_degree + 2),
                                      norm_type);

    const double error =
        VectorTools::compute_global_error(this->mesh,
                                          error_per_cell,
                                          norm_type);

    return error;
}

#endif
