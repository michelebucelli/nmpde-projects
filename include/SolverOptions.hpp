#ifndef SOLVER_OPTIONS_HPP
#define SOLVER_OPTIONS_HPP

enum preconditioner_id { BLOCK_DIAGONAL, SIMPLE, ASIMPLE, YOSHIDA, AYOSHIDA };

struct SolverOptions {
  SolverOptions(const unsigned int &maxiter_, const double &tol_, const preconditioner_id &id_, const double &alpha_ = 1)
      : maxiter(maxiter_), tol(tol_), id(id_), alpha(alpha_) {
    switch (id) {
      case BLOCK_DIAGONAL:
        use_pressure_mass = true;
        use_lumped_mass = false;
        break;

      case SIMPLE:
      case ASIMPLE:
      case YOSHIDA:
        use_pressure_mass = false;
        use_lumped_mass = false;
        break;

      case AYOSHIDA:
        use_pressure_mass = false;
        use_lumped_mass = true;
        break;
    }
  }
  // Maximum number of iterations for the main linear solver.
  const unsigned int maxiter;
  // Relative tolerance for the main linear solver.
  const double tol;
  // Type of preconditioner.
  const preconditioner_id id;
  // Damping parameter (only relevant for SIMPLE and aSIMPLE preconditioners).
  const double alpha;
  // Whether the pressure mass matrix is used.
  bool use_pressure_mass;
  // Whether the lumped velocity mass matrix is used.
  bool use_lumped_mass;
};

#endif