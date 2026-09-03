// Unit tests for the geometric multigrid preconditioner. They check that a
// V-cycle can be applied to double-precision vectors (the mixed-precision
// interface), that the cycle is linear, and that FGMRES iteration counts stay
// bounded when the mesh is refined.

#include "cdr_operator.hpp"
#include "mg_preconditioner.hpp"
#include "test_harness.hpp"

#include <cmath>

#include <deal.II/base/function.h>
#include <deal.II/base/mpi.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/numerics/vector_tools.h>

using namespace dealii;

namespace {

template <int dim, int fe_degree> class MGFixture {
  parallel::distributed::Triangulation<dim> triangulation;
  const MappingQ1<dim> mapping;
  const FE_Q<dim> fe;
  AffineConstraints<double> constraints;

public:
  using Number = double;
  using LevelNumber = float;
  using VectorType = LinearAlgebra::distributed::Vector<Number>;

  MGFixture(const unsigned int n_refinements)
      : triangulation(MPI_COMM_WORLD,
                      Triangulation<dim>::limit_level_difference_at_vertices,
                      parallel::distributed::Triangulation<
                          dim>::construct_multigrid_hierarchy),
        fe(fe_degree), dof_handler(triangulation) {
    GridGenerator::hyper_cube(triangulation, 0.0, 1.0, true);
    triangulation.refine_global(n_refinements);
    dof_handler.distribute_dofs(fe);
    dof_handler.distribute_mg_dofs();

    constraints.reinit(DoFTools::extract_locally_relevant_dofs(dof_handler));
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);

    for (const types::boundary_id id : triangulation.get_boundary_ids()) {
      if (id != neumann_boundary_id) {
        dirichlet_ids.insert(id);
        VectorTools::interpolate_boundary_values(mapping, dof_handler, id,
                                                 Functions::ZeroFunction<dim>(),
                                                 constraints);
      }
    }
    constraints.close();

    auto matrix_free = std::make_shared<MatrixFree<dim, Number>>();
    matrix_free->reinit(
        mapping, dof_handler, constraints, QGauss<1>(fe_degree + 1),
        matrix_free_data<dim, Number>(update_gradients | update_JxW_values));

    op.initialize(matrix_free);
    op.set_coefficients(Coefficients<dim>(1.0, 1.0));

    mg.build(dof_handler, mapping, Coefficients<dim>(1.0, 1.0), dirichlet_ids);
  }

  auto vector(const Number value = 0.0) const -> VectorType {
    VectorType v;
    op.initialize_dof_vector(v);
    v = value;
    constraints.set_zero(v);
    return v;
  }

  auto varying_vector(const Number seed) const -> VectorType {
    VectorType v = vector();
    for (unsigned int i = 0; i < v.locally_owned_size(); ++i) {
      const types::global_dof_index global =
          v.get_partitioner()->local_to_global(i);
      v.local_element(i) = std::sin(seed * (global + 1));
    }
    constraints.set_zero(v);
    return v;
  }

  auto apply_mg(const VectorType &source) const -> VectorType {
    VectorType destination = vector();
    mg.get().vmult(destination, source);
    return destination;
  }

  auto fgmres_iterations(const VectorType &rhs) const -> unsigned int {
    VectorType solution = vector();
    SolverControl solver_control(200, 1e-8 * rhs.l2_norm());
    typename SolverFGMRES<VectorType>::AdditionalData fgmres_data;
    fgmres_data.max_basis_size = 50;
    SolverFGMRES<VectorType> solver(solver_control, fgmres_data);
    solver.solve(op, solution, rhs, mg.get());
    return solver_control.last_step();
  }

  DoFHandler<dim> dof_handler;
  std::set<types::boundary_id> dirichlet_ids;
  Operator<dim, fe_degree, Number> op;
  MultigridPreconditioner<dim, fe_degree, LevelNumber> mg;
};

auto test_vmult_is_finite(TestSuite &suite) -> void {
  const MGFixture<2, 2> fixture(3);
  const auto result = fixture.apply_mg(fixture.varying_vector(0.41));
  suite.check(std::isfinite(result.l2_norm()) && result.l2_norm() > 0.0,
              "MG vmult on a double vector is finite and non-zero");
}

auto test_vmult_is_linear(TestSuite &suite) -> void {
  const MGFixture<2, 2> fixture(3);
  const auto x = fixture.varying_vector(0.29);
  const auto y = fixture.varying_vector(0.83);

  auto combination = x;
  combination.sadd(2.0, 3.0, y);

  auto expected = fixture.apply_mg(x);
  expected.sadd(2.0, 3.0, fixture.apply_mg(y));

  auto difference = fixture.apply_mg(combination);
  difference -= expected;

  // The V-cycle is linear in exact arithmetic, but the MG levels run in
  // float, so double-float round trips of 2x+3y versus 2 P(x)+3 P(y) only
  // agree to single-precision scale.
  suite.check(expected.l2_norm() > 1e-8,
              "linearity test is not comparing zeros");
  suite.close(difference.l2_norm() / expected.l2_norm(), 0.0, 1e-4,
              "mixed-precision V-cycle is linear");
}

template <int dim, int fe_degree>
auto test_fgmres_iteration_count(TestSuite &suite,
                                 const unsigned int n_refinements) -> void {
  const std::string tag = std::to_string(dim) + "D Q" +
                          std::to_string(fe_degree) + " refine " +
                          std::to_string(n_refinements) + " ";

  const MGFixture<dim, fe_degree> fixture(n_refinements);
  const auto rhs = fixture.varying_vector(0.17);
  const unsigned int iterations = fixture.fgmres_iterations(rhs);

  suite.check(iterations > 0, tag + "FGMRES took at least one iteration");
  suite.check(iterations < 25,
              tag + "FGMRES with MG stays under 25 iterations");
}

auto test_iterations_do_not_grow_with_refinement(TestSuite &suite) -> void {
  const MGFixture<2, 2> coarse(3);
  const MGFixture<2, 2> fine(4);
  const unsigned int coarse_its =
      coarse.fgmres_iterations(coarse.varying_vector(0.17));
  const unsigned int fine_its =
      fine.fgmres_iterations(fine.varying_vector(0.17));

  suite.check(fine_its <= coarse_its + 8,
              "MG iteration count stays bounded under one refinement");
}

} // namespace

int main(int argc, char *argv[]) {
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  TestSuite suite("mg_preconditioner");

  test_vmult_is_finite(suite);
  test_vmult_is_linear(suite);

  test_fgmres_iteration_count<2, 2>(suite, 3);
  test_fgmres_iteration_count<3, 2>(suite, 2);
  test_iterations_do_not_grow_with_refinement(suite);

  return suite.report();
}
