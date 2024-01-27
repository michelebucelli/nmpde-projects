#ifndef PRECONDITION_HPP
#define PRECONDITION_HPP

#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

using namespace dealii;

// Abstract class implementing a preconditioner for the Navier-Stokes problem.
class BlockPrecondition {
 public:
  // Virtual destructor.
  virtual ~BlockPrecondition() = default;
  // Application of the preconditioner.
  virtual void vmult(TrilinosWrappers::MPI::BlockVector &dst,
                     const TrilinosWrappers::MPI::BlockVector &src) const = 0;
};

// Block-diagonal preconditioner [M/deltat + A + C  0; 0 1/nu*M_p].
// Adapted from the one proposed for the Stokes problem in laboratory 9.
class PreconditionBlockDiagonal : public BlockPrecondition {
 public:
  // Initialize the preconditioner.
  void initialize(const TrilinosWrappers::SparseMatrix &velocity_stiffness_,
                  const TrilinosWrappers::SparseMatrix &pressure_mass_);

  // Application of the preconditioner.
  void vmult(TrilinosWrappers::MPI::BlockVector &dst,
             const TrilinosWrappers::MPI::BlockVector &src) const override;

 private:
  // Velocity stiffness matrix.
  const TrilinosWrappers::SparseMatrix *velocity_stiffness;

  // Preconditioner used for the velocity block.
  TrilinosWrappers::PreconditionILU preconditioner_velocity;

  // Pressure mass matrix.
  const TrilinosWrappers::SparseMatrix *pressure_mass;

  // Preconditioner used for the pressure block.
  TrilinosWrappers::PreconditionILU preconditioner_pressure;
};

// SIMPLE preconditioner.
class PreconditionSIMPLE : public BlockPrecondition {
 public:
  // Initialize the preconditioner.
  void initialize(const TrilinosWrappers::SparseMatrix &F_matrix_,
                  const TrilinosWrappers::SparseMatrix &negB_matrix_,
                  const TrilinosWrappers::SparseMatrix &Bt_matrix_,
                  const TrilinosWrappers::MPI::BlockVector &vec,
                  const double &alpha_);

  // Application of the preconditioner.
  void vmult(TrilinosWrappers::MPI::BlockVector &dst,
             const TrilinosWrappers::MPI::BlockVector &src) const override;

 private:
  // Damping parameter (must be in (0,1]).
  double alpha;
  // Matrix F (top left block of the system matrix).
  const TrilinosWrappers::SparseMatrix *F_matrix;
  // Matrix -B (bottom left block of the system matrix).
  const TrilinosWrappers::SparseMatrix *negB_matrix;
  // Matrix B^T (top right block of the system matrix).
  const TrilinosWrappers::SparseMatrix *Bt_matrix;
  // Matrix -D^-1, negative inverse diagonal of F.
  TrilinosWrappers::MPI::Vector negDinv_vector;
  // Matrix S := B*D^-1*B^T.
  TrilinosWrappers::SparseMatrix S_matrix;
  // Preconditioner used for the block multiplied by F.
  TrilinosWrappers::PreconditionAMG preconditioner_F;
  // Preconditioner used for the block multiplied by S.
  TrilinosWrappers::PreconditionAMG preconditioner_S;
  // Temporary vector.
  mutable TrilinosWrappers::MPI::BlockVector tmp;
};

// aSIMPLE preconditioner.
class PreconditionaSIMPLE : public BlockPrecondition {
 public:
  // Initialize the preconditioner.
  void initialize(const TrilinosWrappers::SparseMatrix &F_matrix_,
                  const TrilinosWrappers::SparseMatrix &negB_matrix_,
                  const TrilinosWrappers::SparseMatrix &Bt_matrix_,
                  const TrilinosWrappers::MPI::BlockVector &vec,
                  const double &alpha_);

  // Application of the preconditioner.
  void vmult(TrilinosWrappers::MPI::BlockVector &dst,
             const TrilinosWrappers::MPI::BlockVector &src) const override;

 private:
  // Damping parameter (must be in (0,1]).
  double alpha;
  // Matrix F (top left block of the system matrix).
  const TrilinosWrappers::SparseMatrix *F_matrix;
  // Matrix -B (bottom left block of the system matrix).
  const TrilinosWrappers::SparseMatrix *negB_matrix;
  // Matrix B^T (top right block of the system matrix).
  const TrilinosWrappers::SparseMatrix *Bt_matrix;
  // Matrix D, diagonal of F.
  TrilinosWrappers::MPI::Vector D_vector;
  // Matrix D^-1.
  TrilinosWrappers::MPI::Vector Dinv_vector;
  // Matrix -S := -B*D^-1*B^T.
  TrilinosWrappers::SparseMatrix negS_matrix;
  // Preconditioner used for the block multiplied by F.
  TrilinosWrappers::PreconditionAMG preconditioner_F;
  // Preconditioner used for the block multiplied by S.
  TrilinosWrappers::PreconditionAMG preconditioner_S;
  // Temporary vector.
  mutable TrilinosWrappers::MPI::BlockVector tmp;
};

// Yoshida preconditioner.
class PreconditionYoshida : public BlockPrecondition {
 public:
  // Initialize the preconditioner.
  void initialize(const TrilinosWrappers::SparseMatrix &F_matrix_,
                  const TrilinosWrappers::SparseMatrix &negB_matrix_,
                  const TrilinosWrappers::SparseMatrix &Bt_matrix_,
                  const TrilinosWrappers::SparseMatrix &M_dt_matrix_,
                  const TrilinosWrappers::MPI::BlockVector &vec);

  // Application of the preconditioner.
  void vmult(TrilinosWrappers::MPI::BlockVector &dst,
             const TrilinosWrappers::MPI::BlockVector &src) const override;

 private:
  // Matrix F (top left block of the system matrix).
  const TrilinosWrappers::SparseMatrix *F_matrix;
  // Matrix -B (bottom left block of the system matrix).
  const TrilinosWrappers::SparseMatrix *negB_matrix;
  // Matrix B^T (top right block of the system matrix).
  const TrilinosWrappers::SparseMatrix *Bt_matrix;
  // Matrix D^-1, inverse diagonal of M/deltat.
  TrilinosWrappers::MPI::Vector Dinv_vector;
  // Matrix -S := -B*D^-1*B^T.
  TrilinosWrappers::SparseMatrix negS_matrix;
  // Preconditioner used for the block multiplied by F.
  TrilinosWrappers::PreconditionAMG preconditioner_F;
  // Preconditioner used for the block multiplied by S.
  TrilinosWrappers::PreconditionAMG preconditioner_S;
  // Temporary vectors.
  mutable TrilinosWrappers::MPI::BlockVector tmp;
  mutable TrilinosWrappers::MPI::Vector tmp_2;
};

// Yoshida preconditioner.
class PreconditionaYoshida : public BlockPrecondition {
 public:
  // Initialize the preconditioner.
  void initialize(const TrilinosWrappers::SparseMatrix &F_matrix_,
                  const TrilinosWrappers::SparseMatrix &negB_matrix_,
                  const TrilinosWrappers::SparseMatrix &Bt_matrix_,
                  const TrilinosWrappers::SparseMatrix &M_dt_matrix_,
                  const TrilinosWrappers::MPI::BlockVector &vec);

  // Application of the preconditioner.
  void vmult(TrilinosWrappers::MPI::BlockVector &dst,
             const TrilinosWrappers::MPI::BlockVector &src) const override;

 private:
  // Matrix F (top left block of the system matrix).
  const TrilinosWrappers::SparseMatrix *F_matrix;
  // Matrix -B (bottom left block of the system matrix).
  const TrilinosWrappers::SparseMatrix *negB_matrix;
  // Matrix B^T (top right block of the system matrix).
  const TrilinosWrappers::SparseMatrix *Bt_matrix;
  // Matrix deltat*Ml^-1, inverse of the the diagonal matrix whose elements are
  // the sum over all columns of M/deltat.
  TrilinosWrappers::MPI::Vector lumped_invM_dt_matrix;
  // Matrix -S := -B*deltat*Ml^-1*B^T.
  TrilinosWrappers::SparseMatrix negS_matrix;
  // Preconditioner used for the block multiplied by F.
  TrilinosWrappers::PreconditionAMG preconditioner_F;
  // Preconditioner used for the block multiplied by S.
  TrilinosWrappers::PreconditionAMG preconditioner_S;
  // Temporary vector.
  mutable TrilinosWrappers::MPI::BlockVector tmp;
};

#endif