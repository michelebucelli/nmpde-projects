#include "Solver.hpp"
#include <deal.II/lac/solver_gmres.h>

template <typename SolverType>
SolverAdapter* get_new(){
  return static_cast<SolverAdapter*>(new SolverType());
}

const std::unordered_map<std::string, SolverAdapter*(*)()> solvers{
  {"direct", get_new<DirectSolver>},
  {"cg",     get_new<IterativeSolver<SolverCG<TrilinosWrappers::MPI::Vector>>>},
  {"gmres",  get_new<IterativeSolver<SolverGMRES<TrilinosWrappers::MPI::Vector>>>}
};

std::unordered_map<std::string, std::function<char **(char **, char **)>> &SolverAdapter::add_extraoptions(std::unordered_map<std::string, std::function<char **(char **, char **)>> &args_options) {
  return args_options;
};



template <typename SolverType>
std::unordered_map<std::string, std::function<char **(char **, char **)>> &IterativeSolver<SolverType>::add_extraoptions(std::unordered_map<std::string, std::function<char **(char **, char **)>> &args_options)
{
  SolverAdapter::add_extraoptions(args_options);
  std::unordered_map<std::string, std::function<char **(char **, char **)>> extra_options{
    {"--max_iter",       [this](char** begin, char** end)->char** {return this->set_max_iter(begin, end);}},
    {"--tol",            [this](char** begin, char** end)->char** {return this->set_tol(begin, end);}},
    {"--reduce",         [this](char** begin, char** end)->char** {return this->set_recude_factor(begin, end);}},
    {"-p",               [this](char** begin, char** end)->char** {return this->set_preconditioner_type(begin, end);}},
    {"--preconditioner", [this](char** begin, char** end)->char** {return this->set_preconditioner_type(begin, end);}}
  };
  if (precondition) precondition->add_extraoptions(extra_options);
  this->args_options=&args_options;
  for (const auto &pair : extra_options) {
    args_options[pair.first] = pair.second;
  }
  return args_options;
}

template <typename SolverType>
char** IterativeSolver<SolverType>::set_max_iter(char** begin, char** end) {
  if (begin != end) {
    solver_control.set_max_steps(std::atoi(*begin));
    return ++begin;
  }
  return begin;
}

template <typename SolverType>
char** IterativeSolver<SolverType>::set_tol(char** begin, char** end) {
  if (begin != end) {
    solver_control.set_tolerance(std::atof(*begin));
    return ++begin;
  }
  return begin;
}

template <typename SolverType>
char** IterativeSolver<SolverType>::set_recude_factor(char** begin, char** end) {
  if (begin != end) {
    solver_control.set_reduction(std::atof(*begin));
    return ++begin;
  }
  return begin;
}

template <typename SolverType>
char** IterativeSolver<SolverType>::set_preconditioner_type(char** begin, char** end) {
  if (begin != end) {
    std::string type = std::string(*begin);
    PreconditionAdapter *new_preconditioner = PreconditionAdapter::get_new(type);
    new_preconditioner->add_extraoptions(*args_options);
    precondition.reset(new_preconditioner);
    return ++begin;
  }
  return begin;
}

SolverAdapter *SolverAdapter::get_new(std::string solver_name)
{
  auto s = solvers.find(solver_name);
  if (s == solvers.end())
    return nullptr;
  return s->second();
}

void DirectSolver::initialize(const TrilinosWrappers::SparseMatrix &matrix_){
  matrix = &matrix_;
  solver.initialize(matrix_);
  // solver.initialize(matrix_, TrilinosWrappers::SolverDirect::AdditionalData{false, "Amesos_Mumps"});
}

template <typename SolverType>
void IterativeSolver<SolverType>::initialize(const TrilinosWrappers::SparseMatrix &matrix_){
  matrix = &matrix_;
  if (precondition) precondition->initialize(matrix_);
}

void DirectSolver::solve(TrilinosWrappers::MPI::Vector &solution, const TrilinosWrappers::MPI::Vector &rhs){
  solver.vmult(solution, rhs);
}

template <typename SolverType>
void IterativeSolver<SolverType>::solve(TrilinosWrappers::MPI::Vector &solution, const TrilinosWrappers::MPI::Vector &rhs){

  // TrilinosWrappers::PreconditionSSOR preconditioner;
  // preconditioner.initialize(
  //   *matrix, TrilinosWrappers::PreconditionSSOR::AdditionalData(1.0));

  // solver_control = ReductionControl(/* maxiter = */ 10000,
  //                                   /* tolerance = */ 1.0e-16,
  //                                   /* reduce = */ 1.0e-6);
 
  SolverType solver(solver_control);
  solver.solve(*matrix, solution, rhs, precondition->get());
}

unsigned int DirectSolver::get_iterations() const{
  return 1;
}

template <typename SolverType>
unsigned int IterativeSolver<SolverType>::get_iterations() const{
  return solver_control.last_step();
}

std::string DirectSolver::get_name() const{
  return "Direct";
}

template <>
inline std::string IterativeSolver<SolverCG<TrilinosWrappers::MPI::Vector>>::get_name() const {
  return "CG";
}

template <>
inline std::string IterativeSolver<SolverGMRES<TrilinosWrappers::MPI::Vector>>::get_name() const {
  return "GMRES";
}

std::ostream& operator<<(std::ostream& os, const SolverAdapter& solver) {
  return solver.operator<<(os);
}

std::ostream& DirectSolver::operator<<(std::ostream& os) const {
  os << "Direct";
  return os;
}

template <typename SolverType>
std::ostream& IterativeSolver<SolverType>::operator<<(std::ostream& os) const {
  os << get_name() << std::endl
     << "\tMax iter:      " << solver_control.max_steps() << std::endl
     << "\tTolerance:     " << solver_control.tolerance() << std::endl
     << "\tReduction:     " << solver_control.reduction();
  if (precondition) {
    os << std::endl
       << "\tPrecond:       " << *precondition;
  } else {
    os << std::endl
       << "\tPrecond:       none";
  }
  return os;
}

template class IterativeSolver<SolverGMRES<TrilinosWrappers::MPI::Vector>>;
template class IterativeSolver<SolverCG<TrilinosWrappers::MPI::Vector>>;
