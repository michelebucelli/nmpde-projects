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