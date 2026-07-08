#pragma once

#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>

using namespace dealii;

/**
 * Base class for preconditioner adapters.
 */
class PreconditionAdapter {
protected:
  PreconditionAdapter() = default;

  virtual std::ostream& operator<<(std::ostream& os) const = 0;

public:
  virtual ~PreconditionAdapter() = default;

  static PreconditionAdapter *get_new(const std::string &name);

  /**
   * Initialize the preconditioner with the system matrix.
   */
  virtual void initialize(const TrilinosWrappers::SparseMatrix &matrix) = 0;

  /**
   * Get the underlying Trilinos preconditioner.
   */
  virtual const TrilinosWrappers::PreconditionBase& get() const = 0;

  /**
   * Get the name of the preconditioner.
   */
  virtual std::string get_name() const = 0;

  /**
   * Add extra command-line options for the preconditioner.
   */
  virtual std::unordered_map<std::string, std::function<char**(char**, char**)>> &
  add_extraoptions(std::unordered_map<std::string, std::function<char**(char**, char**)>> & args_options);

  friend std::ostream& operator<<(std::ostream& os, const PreconditionAdapter& precondition);
};

/**
 * Implementation of PreconditionAdapter for Trilinos Jacobi preconditioner.
 */
class PreconditionJacobi : public PreconditionAdapter {
protected:
  virtual std::ostream& operator<<(std::ostream& os) const override;

public:
  virtual void initialize(const TrilinosWrappers::SparseMatrix &matrix) override;
  virtual const TrilinosWrappers::PreconditionBase& get() const override;
  virtual std::string get_name() const override;

private:
  TrilinosWrappers::PreconditionJacobi preconditioner;
};

/**
 * Implementation of PreconditionAdapter for Trilinos ILU preconditioner.
 */
class PreconditionILU : public PreconditionAdapter {
protected:
  virtual std::ostream& operator<<(std::ostream& os) const override;

public:
  PreconditionILU();
  virtual void initialize(const TrilinosWrappers::SparseMatrix &matrix) override;
  virtual const TrilinosWrappers::PreconditionBase& get() const override;
  virtual std::string get_name() const override;
  virtual std::unordered_map<std::string, std::function<char**(char**, char**)>> &
  add_extraoptions(std::unordered_map<std::string, std::function<char**(char**, char**)>> & args_options) override;

  char** set_fill_in(char** begin, char** end);

private:
  TrilinosWrappers::PreconditionILU preconditioner;
  unsigned int fill_in;
};

/**
 * Implementation of PreconditionAdapter for Trilinos SSOR preconditioner.
 */
class PreconditionSSOR : public PreconditionAdapter {
protected:
  virtual std::ostream& operator<<(std::ostream& os) const override;

public:
  PreconditionSSOR();
  virtual void initialize(const TrilinosWrappers::SparseMatrix &matrix) override;
  virtual const TrilinosWrappers::PreconditionBase& get() const override;
  virtual std::string get_name() const override;
  virtual std::unordered_map<std::string, std::function<char**(char**, char**)>> &
  add_extraoptions(std::unordered_map<std::string, std::function<char**(char**, char**)>> & args_options) override;

  char** set_omega(char** begin, char** end);

  friend std::ostream& operator<<(std::ostream& os, const PreconditionAdapter& precondition);

private:
  TrilinosWrappers::PreconditionSSOR preconditioner;
  double omega;
};
