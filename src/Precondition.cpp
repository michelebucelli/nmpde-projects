#include "Precondition.hpp"

template <typename PrecondType>
PreconditionAdapter* get_new_precond(){
  return static_cast<PreconditionAdapter*>(new PrecondType());
}

const std::unordered_map<std::string, PreconditionAdapter*(*)()> preconditioners{
  {"jacobi", get_new_precond<PreconditionJacobi>},
  {"ilu",    get_new_precond<PreconditionILU>},
  {"ssor",   get_new_precond<PreconditionSSOR>}
};

std::unordered_map<std::string, std::function<char **(char **, char **)>> &PreconditionAdapter::add_extraoptions(std::unordered_map<std::string, std::function<char **(char **, char **)>> &args_options) {
  return args_options;
}

PreconditionAdapter *PreconditionAdapter::get_new(const std::string &name) {
  auto it = preconditioners.find(name);
  if (it != preconditioners.end())
    return it->second();
  return nullptr;
}

// --- PreconditionJacobi ---

void PreconditionJacobi::initialize(const TrilinosWrappers::SparseMatrix &matrix) {
  preconditioner.initialize(matrix);
}

const TrilinosWrappers::PreconditionBase& PreconditionJacobi::get() const {
  return preconditioner;
}

std::string PreconditionJacobi::get_name() const {
  return "Jacobi";
}

// --- PreconditionILU ---

PreconditionILU::PreconditionILU() : fill_in(DEFAULT_ILU_FILL_IN) {}

void PreconditionILU::initialize(const TrilinosWrappers::SparseMatrix &matrix) {
  preconditioner.initialize(matrix, TrilinosWrappers::PreconditionILU::AdditionalData(fill_in));
}

const TrilinosWrappers::PreconditionBase& PreconditionILU::get() const {
  return preconditioner;
}

std::string PreconditionILU::get_name() const {
  return "ILU";
}

std::unordered_map<std::string, std::function<char **(char **, char **)>> &
PreconditionILU::add_extraoptions(std::unordered_map<std::string, std::function<char **(char **, char **)>> & args_options) {
  args_options["--ilu_fill"] = [this](char** begin, char** end) { return this->set_fill_in(begin, end); };
  return args_options;
}

char** PreconditionILU::set_fill_in(char** begin, char** end) {
  if (begin != end) {
    this->fill_in = std::atoi(*begin);
    return ++begin;
  }
  return begin;
}

// --- PreconditionSSOR ---

PreconditionSSOR::PreconditionSSOR() : omega(DEFAULT_SSOR_OMEGA) {}

void PreconditionSSOR::initialize(const TrilinosWrappers::SparseMatrix &matrix) {
  preconditioner.initialize(matrix, TrilinosWrappers::PreconditionSSOR::AdditionalData(omega));
}

const TrilinosWrappers::PreconditionBase& PreconditionSSOR::get() const {
  return preconditioner;
}

std::string PreconditionSSOR::get_name() const {
  return "SSOR";
}

std::unordered_map<std::string, std::function<char **(char **, char **)>> &
PreconditionSSOR::add_extraoptions(std::unordered_map<std::string, std::function<char **(char **, char **)>> & args_options) {
  args_options["--ssor_omega"] = [this](char** begin, char** end) { return this->set_omega(begin, end); };
  return args_options;
}

char** PreconditionSSOR::set_omega(char** begin, char** end) {
  if (begin != end) {
    this->omega = std::atof(*begin);
    return ++begin;
  }
  return begin;
}

std::ostream& operator<<(std::ostream& os, const PreconditionAdapter& precondition) {
  return precondition.operator<<(os);
}

std::ostream& PreconditionJacobi::operator<<(std::ostream& os) const {
  os << "Jacobi";
  return os;
}

std::ostream& PreconditionILU::operator<<(std::ostream& os) const {
  os << "ILU" << std::endl
     << "\t\tFill-in:       " << fill_in;
  return os;
}

std::ostream& PreconditionSSOR::operator<<(std::ostream& os) const {
  os << "SSOR" << std::endl
     << "\t\tOmega:         " << omega;
  return os;
}
