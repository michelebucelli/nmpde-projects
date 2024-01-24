#ifndef PRECONDITION_HPP
#define PRECONDITION_HPP

#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

using namespace dealii;

// Block-diagonal preconditioner [M/deltat + A + C  0; 0 1/nu*M_p].
class PreconditionBlockDiagonal {
 public:
  // Initialize the preconditioner.
  void initialize(const TrilinosWrappers::SparseMatrix &velocity_stiffness_,
                  const TrilinosWrappers::SparseMatrix &pressure_mass_);

  // Application of the preconditioner.
  void vmult(TrilinosWrappers::MPI::BlockVector &dst,
             const TrilinosWrappers::MPI::BlockVector &src) const;

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

// SIMPLE preconditioner
class PreconditionSIMPLE {
 public:
  // Initialize the preconditioner.
  void initialize(const TrilinosWrappers::SparseMatrix &F_matrix_,
                  const TrilinosWrappers::SparseMatrix &negB_matrix_,
                  const TrilinosWrappers::SparseMatrix &Bt_matrix_,
                  const TrilinosWrappers::MPI::BlockVector &vec,
                  const double &alpha_);

  // Application of the preconditioner.
  void vmult(TrilinosWrappers::MPI::BlockVector &dst,
             const TrilinosWrappers::MPI::BlockVector &src) const;

 private:
  // Damping parameter (must be in (0,1]).
  double alpha;
  // Matrix F (top left block of the system matrix).
  const TrilinosWrappers::SparseMatrix *F_matrix;
  // Matrix -B (bottom left block of the system matrix).
  const TrilinosWrappers::SparseMatrix *negB_matrix;
  // Matrix B^T (top right block of the system matrix).
  const TrilinosWrappers::SparseMatrix *Bt_matrix;
  // Matrix D^-1, inverse diagonal of F.
  TrilinosWrappers::MPI::Vector Dinv_vector;
  // Matrix -S := -B*D^-1*B^T.
  TrilinosWrappers::SparseMatrix negS_matrix;
  // Preconditioner used for the block multiplied by F.
  TrilinosWrappers::PreconditionILU preconditioner_F;
  // Preconditioner used for the block multiplied by S.
  TrilinosWrappers::PreconditionILU preconditioner_S;
  // Temporary vector.
  mutable TrilinosWrappers::MPI::BlockVector tmp;
};

// aSIMPLE preconditioner
class PreconditionaSIMPLE {
 public:
  // Initialize the preconditioner.
  void initialize(const TrilinosWrappers::SparseMatrix &F_matrix_,
                  const TrilinosWrappers::SparseMatrix &negB_matrix_,
                  const TrilinosWrappers::SparseMatrix &Bt_matrix_,
                  const TrilinosWrappers::MPI::BlockVector &vec,
                  const double &alpha_);

  // Application of the preconditioner.
  void vmult(TrilinosWrappers::MPI::BlockVector &dst,
             const TrilinosWrappers::MPI::BlockVector &src) const;

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
  TrilinosWrappers::PreconditionILU preconditioner_F;
  // Preconditioner used for the block multiplied by S.
  TrilinosWrappers::PreconditionILU preconditioner_S;
  // Temporary vector.
  mutable TrilinosWrappers::MPI::BlockVector tmp;
};

#endif