#ifndef PRECONDITION_HPP
#define PRECONDITION_HPP

#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

using namespace dealii;

// Block-diagonal preconditioner.
class PreconditionBlockDiagonal  //[M/deltat + A B^T + C 0; 0 1/mu*M_p]
{
 public:
  // Initialize the preconditioner, given the velocity stiffness matrix, the
  // pressure mass matrix.
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

#endif