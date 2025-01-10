#include <fstream>

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/distributed/tria.h>

#include "operator.h"
#include "setup.h"
#include "functions.h"
#include "problem.h"

using namespace dealii;
using VectorType = LinearAlgebra::distributed::Vector<double>;

template <int dim, int fe_degree, int ref_level>
void run_simulation_matrix_free(ConditionalOStream &pcout){

    const unsigned int mpi_size = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

    parallel::distributed::Triangulation<dim> mesh(MPI_COMM_WORLD);
    FE_Q<dim> fe(fe_degree);
    DoFHandler<dim> dof_handler;
    MappingQ1<dim> mapping;
    QGauss<dim> quad(fe.degree + 1);
    IndexSet locally_dofs;
    AffineConstraints<double> constraints;
    VectorType rhs = 0.0;
    VectorType solution = 0.0;

    setup_problem(mesh, fe, dof_handler, ref_level);

    constraints.clear();
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_dofs);
    constraints.reinit(locally_dofs);
    const types::boundary_id dirichlet_boundary_id = 0;
    GFunction<dim> g_function;
    VectorTools::interpolate_boundary_values(dof_handler, dirichlet_boundary_id, g_function, constraints);
    constraints.close();
    
    CustomOperator<dim, fe_degree> custom_operator(constraints);

    typename MatrixFree<dim, double>::AdditionalData additional_data;
    additional_data.mapping_update_flags = (update_values | update_gradients | update_JxW_values | update_quadrature_points);
    additional_data.mapping_update_flags_boundary_faces = (update_values | update_JxW_values | update_quadrature_points);
    additional_data.mapping_update_flags_faces_by_cells  = (update_values | update_JxW_values | update_quadrature_points);

    std::shared_ptr<MatrixFree<dim, double>> matrix_free_ptr(new MatrixFree<dim, double>());

    matrix_free_ptr->reinit(mapping, dof_handler, constraints, quad, additional_data);
    custom_operator.initialize(matrix_free_ptr);

    custom_operator.evaluate_mu(MuFunction<dim>());
    custom_operator.evaluate_beta(BetaFunction<dim>());
    custom_operator.evaluate_gamma(GammaFunction<dim>());
    custom_operator.initialize_dof_vector(solution);
    custom_operator.initialize_dof_vector(rhs);

    SourceFunction<dim> source_function;
    HFunction<dim> h_function;
    const types::boundary_id neumann_boundary_id = 1;

    rhs = 0;
    FEEvaluation<dim, fe_degree> phi(*custom_operator.get_matrix_free());
    FEFaceEvaluation<dim, fe_degree> phi_face(*custom_operator.get_matrix_free(), true, 0, 0, 0, 0, 0, 0);

    for (unsigned int cell = 0; cell < custom_operator.get_matrix_free()->n_cell_batches(); ++cell)
    {
        phi.reinit(cell);
        for (unsigned int q = 0; q < phi.n_q_points; ++q)
            phi.submit_value(source_function.value(phi.quadrature_point(q)), q);
        phi.integrate(EvaluationFlags::values);
        phi.distribute_local_to_global(rhs);
    }

    for (unsigned int face = 0; face < custom_operator.get_matrix_free()->n_boundary_face_batches(); ++face)
    {
        phi_face.reinit(face);

        if (phi_face.boundary_id() == neumann_boundary_id)
        {
            for (unsigned int q = 0; q < phi_face.n_q_points; ++q)
                phi_face.submit_value(h_function.value(phi_face.quadrature_point(q)), q);

            phi_face.integrate(EvaluationFlags::values);
            phi_face.distribute_local_to_global(rhs);
        }
    }

    rhs.compress(VectorOperation::add);
    rhs.update_ghost_values();
    constraints.distribute(rhs);

    SolverControl solver_control(10000, 1e-12);
    SolverCG<VectorType> solver(solver_control);

    solver.solve(custom_operator, solution, rhs, PreconditionIdentity());

    constraints.distribute(solution);
    solution.update_ghost_values();
    rhs.update_ghost_values();

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "solution");
    data_out.add_data_vector(rhs, "rhs");
    data_out.build_patches(mapping);

    data_out.write_vtu_with_pvtu_record("./", "solution" + std::to_string(mpi_size), 0, MPI_COMM_WORLD, 5);

    pcout << "N-iterations: " << solver_control.last_step() << std::endl;

}

template <int dim, int fe_degree, int ref_level>
void run_simulation_classic(){

    Problem<dim, fe_degree> problem;
    
    problem.setup(ref_level);
    problem.assemble();
    problem.solve();
    problem.output();

    return;
}

int main(int argc, char *argv[]){

    const unsigned int dim       = 2;
    const unsigned int fe_degree = 1;
    const unsigned int ref_level = 9;

    Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
    const unsigned int mpi_size = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
    const unsigned int mpi_rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
    ConditionalOStream pcout(std::cout, mpi_rank == 0);

    pcout << "Running with vectorization on " << VectorizedArray<double>::size() << " cells(" << Utilities::System::get_current_vectorization_level() << ")" << std::endl;
    pcout << std::endl;

    pcout << "Running matrix-free simulation with " << mpi_size << " processors" << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    run_simulation_matrix_free<dim, fe_degree, ref_level>(pcout);
    auto end_time = std::chrono::high_resolution_clock::now();
    pcout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << " ms" << std::endl;
    pcout << std::endl;

    pcout << "Running classic simulation with " << mpi_size << " processors" << std::endl;
    start_time = std::chrono::high_resolution_clock::now();
    run_simulation_classic<dim, fe_degree, ref_level>();
    end_time = std::chrono::high_resolution_clock::now();
    pcout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << " ms" << std::endl;
    pcout << std::endl;

    return 0;
}