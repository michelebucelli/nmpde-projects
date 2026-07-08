#pragma once

#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

//#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_cg.h>

#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/solver_control.h>

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>

class SolverAdapter;
#include "args.hpp"

#include "Precondition.hpp"

using namespace dealii;
using std::literals::operator""s;

/**
 * Base class for solver adapters.
 */
class SolverAdapter {

protected:
  SolverAdapter() = default;

  virtual std::ostream& operator<<(std::ostream& os) const = 0;

public:
  virtual ~SolverAdapter() = default;

  static SolverAdapter *get_new(std::string solver_name);

  /**
   * Initialize the solver with the system matrix.
   */
  virtual void initialize(const TrilinosWrappers::SparseMatrix &matrix) = 0;

  /**
   * Solve the linear system.
   */
  virtual void solve(TrilinosWrappers::MPI::Vector       &solution,
                     const TrilinosWrappers::MPI::Vector &rhs) = 0;
  
  virtual unsigned int get_iterations() const = 0;
  virtual std::string get_name() const = 0;
  virtual std::string get_preconditioner_name() const { return "none"; }
  virtual std::unordered_map<std::string, std::function<char**(char**, char**)>> &add_extraoptions(std::unordered_map<std::string, std::function<char**(char**, char**)>> & args_options);

  friend std::ostream& operator<<(std::ostream& os, const SolverAdapter& solver);
};

/**
 * Implementation of SolverAdapter for Trilinos direct solver.
 */
class DirectSolver : public SolverAdapter {
protected:
  virtual std::ostream& operator<<(std::ostream& os) const override;

public:

  virtual void initialize(const TrilinosWrappers::SparseMatrix &matrix_) override;

  virtual void solve(TrilinosWrappers::MPI::Vector       &solution,
                     const TrilinosWrappers::MPI::Vector &rhs) override;

  virtual unsigned int get_iterations() const override;

  virtual std::string get_name() const override;

private:
  const TrilinosWrappers::SparseMatrix *matrix = nullptr;
  //TrilinosWrappers::SolverDirect solver;
  SparseDirectMUMPS solver;
};

/**
 * Template implementation of SolverAdapter for Trilinos iterative solvers.
 */
template <typename SolverType>
class IterativeSolver : public SolverAdapter {
protected:
  virtual std::ostream& operator<<(std::ostream& os) const override;

public:
  IterativeSolver()
    : solver_control(DEFAULT_MAX_ITER, DEFAULT_TOLERANCE, DEFAULT_REDUCE_FACTOR)
    , precondition(PreconditionAdapter::get_new(DEFAULT_PRECONDITIONER))
  {}

  virtual void initialize(const TrilinosWrappers::SparseMatrix &matrix_) override;

  virtual void solve(TrilinosWrappers::MPI::Vector       &solution,
                     const TrilinosWrappers::MPI::Vector &rhs) override;

  virtual unsigned int get_iterations() const override;

  virtual std::string get_name() const override;

  virtual std::string get_preconditioner_name() const override {
    return precondition ? precondition->get_name() : "none";
  }

  virtual std::unordered_map<std::string, std::function<char**(char**, char**)>> &add_extraoptions(std::unordered_map<std::string, std::function<char**(char**, char**)>> & args_options) override;

  char** set_max_iter(char** begin, char** end);
  char** set_tol(char** begin, char** end);
  char** set_recude_factor(char** begin, char** end);
  char** set_preconditioner_type(char** begin, char** end);

private:
  ReductionControl solver_control;
  const TrilinosWrappers::SparseMatrix *matrix = nullptr;
  std::unique_ptr<PreconditionAdapter> precondition;

  std::unordered_map<std::string, std::function<char **(char **, char **)>> *args_options = nullptr;
};

template <>
inline std::string IterativeSolver<SolverCG<TrilinosWrappers::MPI::Vector>>::get_name() const;

template <>
inline std::string IterativeSolver<SolverGMRES<TrilinosWrappers::MPI::Vector>>::get_name() const;

extern template class IterativeSolver<SolverGMRES<TrilinosWrappers::MPI::Vector>>;
extern template class IterativeSolver<SolverCG<TrilinosWrappers::MPI::Vector>>;

