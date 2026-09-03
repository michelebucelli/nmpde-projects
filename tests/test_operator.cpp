// Unit tests for the matrix-free operator in cdr_operator.hpp. Each test
// isolates one term of a(u, v) = (mu grad u, grad v) + (beta . grad u, v) +
// (gamma u, v) by zeroing the coefficients of the other two, so a wrong sign or
// a missing term shows up on its own rather than as a convergence failure.

#include "cdr_operator.hpp"
#include "test_harness.hpp"

#include <deal.II/base/mpi.h>
#include <deal.II/grid/grid_generator.h>

using namespace dealii;

namespace {

template <int dim>
auto zero_convection(Coefficients<dim> c) -> Coefficients<dim> {
  c.beta = Tensor<1, dim>();
  return c;
}

// A cube meshed once, with an unconstrained operator on top of it. Leaving the
// constraints empty keeps vmult equal to the plain bilinear form, so identities
// such as "the constant is in the kernel" are not perturbed by boundary rows.
template <int dim, int fe_degree> class OperatorFixture {
  // Declared first so that the mesh outlives, and is constructed before, the
  // DoFHandler and the operator that refer to it.
  parallel::distributed::Triangulation<dim> triangulation;
  const MappingQ1<dim> mapping;
  const FE_Q<dim> fe;
  AffineConstraints<double> constraints;

public:
  using Number = double;
  using VectorType = LinearAlgebra::distributed::Vector<Number>;

  OperatorFixture(const Coefficients<dim> &coefficients,
                  const unsigned int n_refinements)
      : triangulation(MPI_COMM_WORLD), fe(fe_degree),
        dof_handler(triangulation) {
    GridGenerator::hyper_cube(triangulation, 0.0, 1.0, true);
    triangulation.refine_global(n_refinements);
    dof_handler.distribute_dofs(fe);

    constraints.reinit(DoFTools::extract_locally_relevant_dofs(dof_handler));
    constraints.close();

    auto matrix_free = std::make_shared<MatrixFree<dim, Number>>();
    matrix_free->reinit(
        mapping, dof_handler, constraints, QGauss<1>(fe_degree + 1),
        matrix_free_data<dim, Number>(update_gradients | update_JxW_values));

    op.initialize(matrix_free);
    op.set_coefficients(coefficients);
  }

  auto vector(const Number value = 0.0) const -> VectorType {
    VectorType v;
    op.initialize_dof_vector(v);
    v = value;
    return v;
  }

  // A reproducible non-constant vector, so tests do not accidentally probe the
  // operator only on its kernel.
  auto varying_vector(const Number seed) const -> VectorType {
    VectorType v = vector();
    for (unsigned int i = 0; i < v.locally_owned_size(); ++i) {
      const types::global_dof_index global =
          v.get_partitioner()->local_to_global(i);
      v.local_element(i) = std::sin(seed * (global + 1));
    }
    return v;
  }

  auto apply(const VectorType &source) const -> VectorType {
    VectorType destination = vector();
    op.vmult(destination, source);
    return destination;
  }

  DoFHandler<dim> dof_handler;
  Operator<dim, fe_degree, Number> op;
};

auto test_matrix_free_data(TestSuite &suite) -> void {
  const auto data =
      matrix_free_data<2, double>(update_gradients | update_JxW_values);
  suite.check(data.tasks_parallel_scheme ==
                  MatrixFree<2, double>::AdditionalData::partition_partition,
              "matrix_free_data enables the TBB partition_partition scheduler");
  suite.check(data.mapping_update_flags ==
                  (update_gradients | update_JxW_values),
              "matrix_free_data forwards the cell update flags");
}

// Diffusion and convection both differentiate their argument, so a constant
// vector has to land in the kernel of either term on its own.
template <int dim, int fe_degree>
auto test_constants_in_kernel(TestSuite &suite) -> void {
  const std::string tag =
      std::to_string(dim) + "D Q" + std::to_string(fe_degree) + " ";

  const OperatorFixture<dim, fe_degree> diffusion(
      zero_convection(Coefficients<dim>(1.0, 0.0)), 2);
  suite.close(diffusion.apply(diffusion.vector(1.0)).l2_norm(), 0.0, 1e-12,
              tag + "diffusion annihilates constants");

  Coefficients<dim> convection_only(0.0, 0.0);
  const OperatorFixture<dim, fe_degree> convection(convection_only, 2);
  suite.close(convection.apply(convection.vector(1.0)).l2_norm(), 0.0, 1e-12,
              tag + "convection annihilates constants");
}

// With only the reaction term left the operator is the mass matrix, whose
// entries sum to the measure of the unit cube.
template <int dim, int fe_degree>
auto test_reaction_is_the_mass_matrix(TestSuite &suite) -> void {
  const std::string tag =
      std::to_string(dim) + "D Q" + std::to_string(fe_degree) + " ";

  for (const double gamma : {1.0, 2.5}) {
    const OperatorFixture<dim, fe_degree> fixture(
        zero_convection(Coefficients<dim>(0.0, gamma)), 2);
    const auto ones = fixture.vector(1.0);
    suite.close(ones * fixture.apply(ones), gamma, 1e-12,
                tag + "mass matrix sums to gamma |Omega| at gamma = " +
                    std::to_string(gamma));
  }
}

// Diffusion and reaction are symmetric bilinear forms; the assembled operator
// has to reflect that.
template <int dim, int fe_degree>
auto test_symmetry_without_convection(TestSuite &suite) -> void {
  const std::string tag =
      std::to_string(dim) + "D Q" + std::to_string(fe_degree) + " ";

  const OperatorFixture<dim, fe_degree> fixture(
      zero_convection(Coefficients<dim>(1.5, 0.75)), 2);
  const auto x = fixture.varying_vector(0.37);
  const auto y = fixture.varying_vector(1.11);

  const double xAy = x * fixture.apply(y);
  const double yAx = y * fixture.apply(x);
  suite.close(xAy, yAx, 1e-10,
              tag + "diffusion-reaction operator is symmetric");
  suite.check(std::abs(xAy) > 1e-8,
              tag + "symmetry test uses vectors outside the kernel");

  // Convection is skew, so switching it on has to break that symmetry.
  const OperatorFixture<dim, fe_degree> with_convection(
      Coefficients<dim>(1.5, 0.75), 2);
  suite.check(std::abs(x * with_convection.apply(y) -
                       y * with_convection.apply(x)) > 1e-8,
              tag + "convection makes the operator non-symmetric");
}

// a(u, v) is affine in (mu, gamma, beta), so the full operator must equal the
// sum of its three single-term counterparts.
template <int dim, int fe_degree>
auto test_terms_superpose(TestSuite &suite) -> void {
  const std::string tag =
      std::to_string(dim) + "D Q" + std::to_string(fe_degree) + " ";

  const Coefficients<dim> full(1.7, 0.9);
  const OperatorFixture<dim, fe_degree> combined(full, 2);

  Coefficients<dim> diffusion_only(full.mu, 0.0);
  diffusion_only.beta = Tensor<1, dim>();
  Coefficients<dim> reaction_only(0.0, full.gamma);
  reaction_only.beta = Tensor<1, dim>();
  Coefficients<dim> convection_only(0.0, 0.0);
  convection_only.beta = full.beta;

  const OperatorFixture<dim, fe_degree> diffusion(diffusion_only, 2);
  const OperatorFixture<dim, fe_degree> reaction(reaction_only, 2);
  const OperatorFixture<dim, fe_degree> convection(convection_only, 2);

  const auto x = combined.varying_vector(0.53);
  auto expected = diffusion.apply(x);
  expected += reaction.apply(x);
  expected += convection.apply(x);

  auto difference = combined.apply(x);
  difference -= expected;
  suite.close(difference.l2_norm(), 0.0, 1e-12,
              tag + "the three terms superpose");
  suite.check(expected.l2_norm() > 1e-8,
              tag + "superposition test is not comparing zeros");
}

template <int dim, int fe_degree>
auto test_linearity(TestSuite &suite) -> void {
  const std::string tag =
      std::to_string(dim) + "D Q" + std::to_string(fe_degree) + " ";

  const OperatorFixture<dim, fe_degree> fixture(Coefficients<dim>(1.2, 0.4), 2);
  const auto x = fixture.varying_vector(0.29);
  const auto y = fixture.varying_vector(0.83);

  auto combination = x;
  combination.sadd(2.0, 3.0, y); // 2 x + 3 y

  auto expected = fixture.apply(x);
  expected.sadd(2.0, 3.0, fixture.apply(y));

  auto difference = fixture.apply(combination);
  difference -= expected;
  suite.close(difference.l2_norm(), 0.0, 1e-12, tag + "vmult is linear");
}

// compute_diagonal stores 1 / A_ii, which has to agree with the i-th entry of
// A e_i extracted one unit vector at a time.
template <int dim, int fe_degree>
auto test_compute_diagonal(TestSuite &suite) -> void {
  const std::string tag =
      std::to_string(dim) + "D Q" + std::to_string(fe_degree) + " ";

  OperatorFixture<dim, fe_degree> fixture(Coefficients<dim>(1.3, 0.6), 2);
  fixture.op.compute_diagonal();
  const auto &inverse_diagonal =
      fixture.op.get_matrix_diagonal_inverse()->get_vector();

  const types::global_dof_index n_dofs = fixture.dof_handler.n_dofs();
  unsigned int n_compared = 0;

  for (types::global_dof_index i = 0; i < n_dofs; ++i) {
    auto unit = fixture.vector();
    const bool owned = unit.get_partitioner()->in_local_range(i);
    if (owned) {
      unit(i) = 1.0;
    }

    const auto column = fixture.apply(unit);
    if (!owned) {
      continue;
    }

    const double entry = column(i);
    suite.check(std::abs(entry) > 1e-10,
                tag + "diagonal entry " + std::to_string(i) + " is non-zero");
    suite.close(inverse_diagonal(i), 1.0 / entry, 1e-10,
                tag + "inverse diagonal entry " + std::to_string(i));
    ++n_compared;
  }

  // A rank may legitimately own no degrees of freedom, so the count is checked
  // globally rather than per rank.
  suite.check(Utilities::MPI::sum(n_compared, MPI_COMM_WORLD) == n_dofs,
              tag + "every diagonal entry was compared");
}

} // namespace

int main(int argc, char *argv[]) {
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  TestSuite suite("operator");

  test_matrix_free_data(suite);

  test_constants_in_kernel<2, 1>(suite);
  test_constants_in_kernel<2, 2>(suite);
  test_constants_in_kernel<3, 2>(suite);

  test_reaction_is_the_mass_matrix<2, 1>(suite);
  test_reaction_is_the_mass_matrix<2, 2>(suite);
  test_reaction_is_the_mass_matrix<3, 2>(suite);

  test_symmetry_without_convection<2, 2>(suite);
  test_symmetry_without_convection<3, 2>(suite);

  test_terms_superpose<2, 2>(suite);
  test_terms_superpose<3, 2>(suite);

  test_linearity<2, 2>(suite);

  test_compute_diagonal<2, 1>(suite);
  test_compute_diagonal<2, 2>(suite);

  return suite.report();
}
