#include "Precondition.hpp"

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>

void PreconditionBlockDiagonal::initialize(
    const TrilinosWrappers::SparseMatrix &velocity_stiffness_,
    const TrilinosWrappers::SparseMatrix &pressure_mass_,
    const unsigned int &maxit_, const double &tol_) {
  maxit = maxit_;
  tol = tol_;
  // Save a reference to the input matrices.
  velocity_stiffness = &velocity_stiffness_;
  pressure_mass = &pressure_mass_;

  // Initialize the preconditioners.
  preconditioner_velocity.initialize(velocity_stiffness_);
  preconditioner_pressure.initialize(pressure_mass_);
}

void PreconditionBlockDiagonal::vmult(
    TrilinosWrappers::MPI::BlockVector &dst,
    const TrilinosWrappers::MPI::BlockVector &src) const {
  // Solve the top-left block.
  SolverControl solver_control_velocity(maxit, tol * src.block(0).l2_norm());
  SolverGMRES<TrilinosWrappers::MPI::Vector> solver_gmres_velocity(
      solver_control_velocity);
  solver_gmres_velocity.solve(*velocity_stiffness, dst.block(0), src.block(0),
                              preconditioner_velocity);

  // Solve the bottom-right block.
  SolverControl solver_control_pressure(maxit, tol * src.block(1).l2_norm());
  SolverCG<TrilinosWrappers::MPI::Vector> solver_cg_pressure(
      solver_control_pressure);
  solver_cg_pressure.solve(*pressure_mass, dst.block(1), src.block(1),
                           preconditioner_pressure);
}

void PreconditionSIMPLE::initialize(
    const TrilinosWrappers::SparseMatrix &F_matrix_,
    const TrilinosWrappers::SparseMatrix &negB_matrix_,
    const TrilinosWrappers::SparseMatrix &Bt_matrix_,
    const TrilinosWrappers::MPI::BlockVector &vec, const double &alpha_,
    const unsigned int &maxit_, const double &tol_) {
  alpha = alpha_;
  maxit = maxit_;
  tol = tol_;
  // Save a reference to the input matrices.
  F_matrix = &F_matrix_;
  negB_matrix = &negB_matrix_;
  Bt_matrix = &Bt_matrix_;

  // Save the negated inverse diagonal of F.
  negDinv_vector.reinit(vec.block(0));
  for (unsigned int index : negDinv_vector.locally_owned_elements()) {
    negDinv_vector[index] = -1.0 / F_matrix->diag_element(index);
  }

  // Create the matrix S.
  negB_matrix->mmult(S_matrix, *Bt_matrix, negDinv_vector);

  // Initialize the preconditioners.
  preconditioner_F.initialize(*F_matrix);
  preconditioner_S.initialize(S_matrix);
}

void PreconditionSIMPLE::vmult(
    TrilinosWrappers::MPI::BlockVector &dst,
    const TrilinosWrappers::MPI::BlockVector &src) const {
  tmp.reinit(src);
  // Step 1: solve [F 0; B -S]sol1 = src.
  // Step 1.1: solve F*sol1_u = src_u.
  SolverControl solver_control_F(maxit, tol * src.block(0).l2_norm());
  SolverGMRES<TrilinosWrappers::MPI::Vector> solver_F(solver_control_F);
  solver_F.solve(*F_matrix, tmp.block(0), src.block(0), preconditioner_F);
  // Step 1.2: solve S*sol1_p = B*sol1_u - src_p.
  Bt_matrix->Tvmult(tmp.block(1), tmp.block(0));
  tmp.block(1) -= src.block(1);
  SolverControl solver_control_S(maxit, tol * tmp.block(1).l2_norm());
  SolverGMRES<TrilinosWrappers::MPI::Vector> solver_S(solver_control_S);
  solver_S.solve(S_matrix, dst.block(1), tmp.block(1), preconditioner_S);

  // Step 2: solve [I D^-1*B^T; 0 alpha*I]dst = sol1.
  // Step 2.1: solve alpha*I*dst_p = sol1_p.
  dst.block(1) /= alpha;
  // Step 2.2: solve dst_u = sol1_u - D^-1*B^T*dst_p.
  dst.block(0) = tmp.block(0);
  Bt_matrix->vmult(tmp.block(0), dst.block(1));
  tmp.block(0).scale(negDinv_vector);
  dst.block(0) += tmp.block(0);
}

void PreconditionaSIMPLE::initialize(
    const TrilinosWrappers::SparseMatrix &F_matrix_,
    const TrilinosWrappers::SparseMatrix &negB_matrix_,
    const TrilinosWrappers::SparseMatrix &Bt_matrix_,
    const TrilinosWrappers::MPI::BlockVector &vec, const double &alpha_,
    const bool &use_inner_solver_, const unsigned int &maxit_,
    const double &tol_) {
  alpha = alpha_;
  use_inner_solver = use_inner_solver_;
  maxit = maxit_;
  tol = tol_;
  // Save a reference to the input matrices.
  F_matrix = &F_matrix_;
  negB_matrix = &negB_matrix_;
  Bt_matrix = &Bt_matrix_;

  // Save the diagonal and inverse diagonal of F.
  D_vector.reinit(vec.block(0));
  Dinv_vector.reinit(vec.block(0));
  for (unsigned int index : D_vector.locally_owned_elements()) {
    const double value = F_matrix->diag_element(index);
    D_vector[index] = value;
    Dinv_vector[index] = 1.0 / value;
  }

  // Create the matrix -S.
  negB_matrix->mmult(negS_matrix, *Bt_matrix, Dinv_vector);

  // Initialize the preconditioners.
  preconditioner_F.initialize(*F_matrix);
  preconditioner_S.initialize(negS_matrix);
}

void PreconditionaSIMPLE::vmult(
    TrilinosWrappers::MPI::BlockVector &dst,
    const TrilinosWrappers::MPI::BlockVector &src) const {
  tmp.reinit(src);
  // Step 1: multiply src by [F^-1 0; 0 I].
  if (use_inner_solver) {
    SolverControl solver_control_F(maxit, tol * src.block(0).l2_norm());
    SolverGMRES<TrilinosWrappers::MPI::Vector> solver_F(solver_control_F);
    solver_F.solve(*F_matrix, dst.block(0), src.block(0), preconditioner_F);
  } else {
    preconditioner_F.vmult(dst.block(0), src.block(0));
  }
  tmp.block(1) = src.block(1);
  // Step 2: multiply the result by [I 0; -B I].
  negB_matrix->vmult_add(tmp.block(1), dst.block(0));
  // Step 3: multiply the result by [I 0; 0 -S^-1].
  if (use_inner_solver) {
    SolverControl solver_control_S(maxit, tol * tmp.block(1).l2_norm());
    SolverGMRES<TrilinosWrappers::MPI::Vector> solver_S(solver_control_S);
    solver_S.solve(negS_matrix, dst.block(1), tmp.block(1), preconditioner_S);
  } else {
    preconditioner_S.vmult(dst.block(1), tmp.block(1));
  }
  // Step 4: multiply the result by [D 0; 0 I/alpha].
  dst.block(0).scale(D_vector);
  dst.block(1) /= alpha;
  // Step 5: multiply the result by [I -B^T; 0 I].
  Bt_matrix->vmult(tmp.block(0), dst.block(1));
  dst.block(0) -= tmp.block(0);
  // Step 6: multiply the result by [D^-1 0; 0 I].
  dst.block(0).scale(Dinv_vector);
}

void PreconditionYoshida::initialize(
    const TrilinosWrappers::SparseMatrix &F_matrix_,
    const TrilinosWrappers::SparseMatrix &negB_matrix_,
    const TrilinosWrappers::SparseMatrix &Bt_matrix_,
    const TrilinosWrappers::SparseMatrix &M_dt_matrix_,
    const TrilinosWrappers::MPI::BlockVector &vec, const unsigned int &maxit_,
    const double &tol_) {
  maxit = maxit_;
  tol = tol_;
  // Save a reference to the input matrices.
  F_matrix = &F_matrix_;
  negB_matrix = &negB_matrix_;
  Bt_matrix = &Bt_matrix_;

  // Save the inverse diagonal of M_dt.
  Dinv_vector.reinit(vec.block(0));
  for (unsigned int index : Dinv_vector.locally_owned_elements()) {
    Dinv_vector[index] = 1.0 / M_dt_matrix_.diag_element(index);
  }

  // Create the matrix -S.
  negB_matrix->mmult(negS_matrix, *Bt_matrix, Dinv_vector);

  // Initialize the preconditioners.
  preconditioner_F.initialize(*F_matrix);
  preconditioner_S.initialize(negS_matrix);
}

void PreconditionYoshida::vmult(
    TrilinosWrappers::MPI::BlockVector &dst,
    const TrilinosWrappers::MPI::BlockVector &src) const {
  tmp.reinit(src);
  // Step 1: solve [F 0; B -S]sol1 = src.
  // Step 1.1: solve F*sol1_u = src_u.
  tmp.block(0) = dst.block(0);
  SolverControl solver_control_F(maxit, tol * src.block(0).l2_norm());
  SolverGMRES<TrilinosWrappers::MPI::Vector> solver_F(solver_control_F);
  solver_F.solve(*F_matrix, tmp.block(0), src.block(0), preconditioner_F);
  // Step 1.2: solve -S*sol1_p = -B*sol1_u + src_p.
  tmp.block(1) = src.block(1);
  negB_matrix->vmult_add(tmp.block(1), tmp.block(0));
  SolverControl solver_control_S(maxit, tol * tmp.block(1).l2_norm());
  SolverGMRES<TrilinosWrappers::MPI::Vector> solver_S(solver_control_S);
  solver_S.solve(negS_matrix, dst.block(1), tmp.block(1), preconditioner_S);

  // Step 2: solve [I F^-1*B^T; 0 I]dst = sol1.
  tmp_2 = src.block(0);
  dst.block(0) = tmp.block(0);
  Bt_matrix->vmult(tmp.block(0), dst.block(1));
  SolverControl solver_control_F2(maxit, tol * tmp.block(0).l2_norm());
  SolverGMRES<TrilinosWrappers::MPI::Vector> solver_gmres_F2(solver_control_F);
  solver_gmres_F2.solve(*F_matrix, tmp_2, tmp.block(0), preconditioner_F);
  dst.block(0) -= tmp_2;
}