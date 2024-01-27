#ifndef PRECONDITIONER_TYPE_HPP
#define PRECONDITIONER_TYPE_HPP

enum preconditioner_id { BLOCK_DIAGONAL, SIMPLE, ASIMPLE, YOSHIDA, AYOSHIDA };

struct PreconditionerType {
  PreconditionerType(const preconditioner_id &id_, const double &alpha_ = 1)
      : id(id_), alpha(alpha_) {
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
  const preconditioner_id id;
  const double alpha;
  bool use_pressure_mass;
  bool use_lumped_mass;
};

#endif