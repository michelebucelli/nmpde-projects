#include "Precondition.hpp"

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>

void PreconditionBlockDiagonal::initialize(
    const TrilinosWrappers::SparseMatrix &velocity_stiffness_,
    const TrilinosWrappers::SparseMatrix &pressure_mass_) {
  velocity_stiffness = &velocity_stiffness_;
  pressure_mass = &pressure_mass_;

  preconditioner_velocity.initialize(velocity_stiffness_);
  preconditioner_pressure.initialize(pressure_mass_);
}

void PreconditionBlockDiagonal::vmult(
    TrilinosWrappers::MPI::BlockVector &dst,
    const TrilinosWrappers::MPI::BlockVector &src) const {
  SolverControl solver_control_velocity(1000, 1e-2 * src.block(0).l2_norm());
  SolverGMRES<TrilinosWrappers::MPI::Vector> solver_gmres_velocity(
      solver_control_velocity);
  solver_gmres_velocity.solve(*velocity_stiffness, dst.block(0), src.block(0),
                              preconditioner_velocity);

  SolverControl solver_control_pressure(1000, 1e-2 * src.block(1).l2_norm());
  SolverCG<TrilinosWrappers::MPI::Vector> solver_cg_pressure(
      solver_control_pressure);
  solver_cg_pressure.solve(*pressure_mass, dst.block(1), src.block(1),
                           preconditioner_pressure);
}

void PreconditionSIMPLE::initialize(
    const TrilinosWrappers::SparseMatrix &F_matrix_,
    const TrilinosWrappers::SparseMatrix &negB_matrix_,
    const TrilinosWrappers::SparseMatrix &Bt_matrix_,
    const TrilinosWrappers::MPI::BlockVector &vec, const double &alpha_) {
  F_matrix = &F_matrix_;
  negB_matrix = &negB_matrix_;
  Bt_matrix = &Bt_matrix_;
  alpha = alpha_;

  Dinv_vector.reinit(vec.block(0));
  tmp.reinit(vec);
  tmp.block(0).add(1.0);
  TrilinosWrappers::PreconditionJacobi precondition_jacobi;
  precondition_jacobi.initialize(*F_matrix);
  precondition_jacobi.vmult(Dinv_vector, tmp.block(0));

  negB_matrix->mmult(negS_matrix, *Bt_matrix, Dinv_vector);
}

void PreconditionSIMPLE::vmult(
    TrilinosWrappers::MPI::BlockVector &dst,
    const TrilinosWrappers::MPI::BlockVector &src) const {
  tmp.reinit(src);
  // Step 1: solve [F 0; B -S]sol1 = src.
  // Step 1.1: solve F*sol1_u = src_u.
  SolverControl solver_control_F(1000, 1e-2 * src.block(0).l2_norm());
  SolverGMRES<TrilinosWrappers::MPI::Vector> solver_gmres_F(solver_control_F);
  solver_gmres_F.solve(*F_matrix, tmp.block(0), src.block(0), preconditioner_F);
  // Step 1.2: solve -S*sol1_p = -B*sol1_u + src_p.
  tmp.block(1) = src.block(1);
  negB_matrix->vmult_add(tmp.block(1), tmp.block(0));
  SolverControl solver_control_S(1000, 1e-2 * src.block(0).l2_norm());
  SolverCG<TrilinosWrappers::MPI::Vector> solver_cg_S(solver_control_S);
  solver_cg_S.solve(negS_matrix, dst.block(1), tmp.block(1), preconditioner_S);

  // Step 2: solve [I D^-1*B^T; 0 alpha*I]dst = sol1.
  // Step 2.1: solve alpha*I*dst_p = sol1_p.
  dst.block(1) /= alpha;
  // Step 2.2: solve dst_u = sol1_u - D^-1*B^T*dst_p.
  dst.block(0) = tmp.block(0);
  tmp.block(0).scale(Dinv_vector);
  Bt_matrix->vmult(tmp.block(0), dst.block(1));
  dst.block(0) -= tmp.block(0);
}