#ifndef SOLVER_OPTIONS_HPP
#define SOLVER_OPTIONS_HPP

enum preconditioner_id { BLOCK_DIAGONAL, SIMPLE, ASIMPLE, YOSHIDA };

struct SolverOptions {
  SolverOptions(const preconditioner_id &id_, const unsigned int &maxiter_,
                const double &tol_, const bool &use_inner_solver_,
                const unsigned int &maxiter_inner_ = 1000,
                const double &tol_inner_ = 1e-2, const double &alpha_ = 1.0)
      : id(id_),
        maxiter(maxiter_),
        tol(tol_),
        use_inner_solver(use_inner_solver_),
        maxiter_inner(maxiter_inner_),
        tol_inner(tol_inner_),
        alpha(alpha_) {
    switch (id) {
      case BLOCK_DIAGONAL:
        use_pressure_mass = true;
        break;

      case SIMPLE:
      case ASIMPLE:
      case YOSHIDA:
        use_pressure_mass = false;
        break;
    }
  }
  // Type of preconditioner.
  const preconditioner_id id;
  // Maximum number of iterations for the main linear solver.
  const unsigned int maxiter;
  // Relative tolerance for the main linear solver.
  const double tol;
  // Whether to use an inner solver (only relevant for aSIMPLE and aYoshida
  // preconditioners).
  const bool use_inner_solver;
  // Maximum number of iterations for the inner linear solvers (if any).
  const unsigned int maxiter_inner;
  // Relative tolerance for the inner linear solvers (if any).
  const double tol_inner;
  // Damping parameter (only relevant for SIMPLE and aSIMPLE preconditioners).
  const double alpha;
  // Whether the pressure mass matrix is used.
  bool use_pressure_mass;
};

#endif