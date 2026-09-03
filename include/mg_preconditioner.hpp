#pragma once

#include "cdr_operator.hpp"

#include <deal.II/base/config.h>
#include <deal.II/base/mg_level_object.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/lac/precondition.h>

#include <deal.II/matrix_free/operators.h>

#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/multigrid.h>

#include <memory>

using namespace dealii;

// Geometric multigrid V-cycle for the matrix-free CDR operator. Level
// matrices, smoothers and transfers all run in LevelNumber (float in the
// solver), while PreconditionMG::vmult accepts the double-precision system
// vectors and converts at the fine-grid interface.
template <int dim, int fe_degree, typename LevelNumber>
class MultigridPreconditioner {
public:
  using LevelMatrixType = Operator<dim, fe_degree, LevelNumber>;
  using LevelVectorType = LinearAlgebra::distributed::Vector<LevelNumber>;
  using TransferType = MGTransferMatrixFree<dim, LevelNumber>;
  using SmootherType = PreconditionChebyshev<LevelMatrixType, LevelVectorType>;
  using Preconditioner = PreconditionMG<dim, LevelVectorType, TransferType>;

  auto build(const DoFHandler<dim> &dof_handler, const MappingQ1<dim> &mapping,
             const Coefficients<dim> &coefficients,
             const std::set<types::boundary_id> &dirichlet_boundary_ids)
      -> void {
    clear();

    const unsigned int n_levels =
        dof_handler.get_triangulation().n_global_levels();

    mg_constrained_dofs.initialize(dof_handler);
    mg_constrained_dofs.make_zero_boundary_constraints(dof_handler,
                                                       dirichlet_boundary_ids);

    mg_matrices.resize(0, n_levels - 1);
    const Coefficients<dim, LevelNumber> level_coefficients(coefficients);

    for (unsigned int level = 0; level < n_levels; ++level) {
      AffineConstraints<double> level_constraints;
      level_constraints.reinit(
          DoFTools::extract_locally_relevant_level_dofs(dof_handler, level));
      level_constraints.add_lines(
          mg_constrained_dofs.get_boundary_indices(level));
      level_constraints.close();

      auto additional_data = matrix_free_data<dim, LevelNumber>(
          update_gradients | update_JxW_values);
      additional_data.mg_level = level;

      auto matrix_free = std::make_shared<MatrixFree<dim, LevelNumber>>();
      matrix_free->reinit(mapping, dof_handler, level_constraints,
                          QGauss<1>(fe_degree + 1), additional_data);

      mg_matrices[level].initialize(matrix_free, mg_constrained_dofs, level);
      mg_matrices[level].set_coefficients(level_coefficients);
    }

    mg_transfer = std::make_unique<TransferType>(mg_constrained_dofs);
    mg_transfer->build(dof_handler);

    MGLevelObject<typename SmootherType::AdditionalData> smoother_data;
    smoother_data.resize(0, n_levels - 1);
    for (unsigned int level = 0; level < n_levels; ++level) {
      if (level > 0) {
        smoother_data[level].smoothing_range = 15.;
        smoother_data[level].degree = 5;
        smoother_data[level].eig_cg_n_iterations = 10;
#if DEAL_II_VERSION_GTE(9, 6, 0)
        // Convection makes A non-symmetric; the power iteration only needs
        // vmult, unlike the default Lanczos/CG estimate.
        smoother_data[level].eigenvalue_algorithm =
            SmootherType::AdditionalData::EigenvalueAlgorithm::power_iteration;
#endif
      } else {
        smoother_data[0].smoothing_range = 1e-3;
        smoother_data[0].degree = numbers::invalid_unsigned_int;
        smoother_data[0].eig_cg_n_iterations = mg_matrices[0].m();
      }

      mg_matrices[level].compute_diagonal();
      smoother_data[level].preconditioner =
          mg_matrices[level].get_matrix_diagonal_inverse();
    }
    mg_smoother.initialize(mg_matrices, smoother_data);

    mg_coarse.initialize(mg_smoother);

    mg_interface_matrices.resize(0, n_levels - 1);
    for (unsigned int level = 0; level < n_levels; ++level) {
      mg_interface_matrices[level].initialize(mg_matrices[level]);
    }

    mg_matrix.initialize(mg_matrices);
    mg_interface.initialize(mg_interface_matrices);

    mg = std::make_unique<Multigrid<LevelVectorType>>(
        mg_matrix, mg_coarse, *mg_transfer, mg_smoother, mg_smoother);
    mg->set_edge_matrices(mg_interface, mg_interface);

    preconditioner =
        std::make_unique<Preconditioner>(dof_handler, *mg, *mg_transfer);
  }

  auto clear() -> void {
    preconditioner.reset();
    mg.reset();
    mg_transfer.reset();
    mg_smoother.clear();
    mg_coarse.clear();
    mg_matrix.reset();
    mg_interface.reset();
    mg_interface_matrices.clear_elements();
    mg_matrices.clear_elements();
    mg_constrained_dofs.clear();
  }

  auto get() const -> const Preconditioner & { return *preconditioner; }

  // Bytes held by the level matrices and the transfers on this rank. Zero
  // when the V-cycle was never built, so the caller can add it blindly.
  auto memory_consumption() const -> std::size_t {
    if (preconditioner == nullptr) {
      return 0;
    }

    std::size_t bytes = mg_transfer->memory_consumption();
    for (unsigned int level = mg_matrices.min_level();
         level <= mg_matrices.max_level(); ++level) {
      bytes += mg_matrices[level].get_matrix_free()->memory_consumption();
    }
    return bytes;
  }

private:
  using Smoother = mg::SmootherRelaxation<SmootherType, LevelVectorType>;
  using CoarseGrid = MGCoarseGridApplySmoother<LevelVectorType>;
  using WrappedMatrix = mg::Matrix<LevelVectorType>;

  MGConstrainedDoFs mg_constrained_dofs;
  MGLevelObject<LevelMatrixType> mg_matrices;
  MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<LevelMatrixType>>
      mg_interface_matrices;
  std::unique_ptr<TransferType> mg_transfer;
  Smoother mg_smoother;
  CoarseGrid mg_coarse;
  WrappedMatrix mg_matrix;
  WrappedMatrix mg_interface;
  std::unique_ptr<Multigrid<LevelVectorType>> mg;
  std::unique_ptr<Preconditioner> preconditioner;
};
