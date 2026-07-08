#ifndef CURRENT_HPP
#define CURRENT_HPP

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <cmath>
#include <memory>
#include "Solver.hpp"

using namespace dealii;

/**
 * Class managing the differential problem.
 */
class Current
{
public:
  // Physical dimension (1D, 2D, 3D)
  static constexpr unsigned int dim = 3;

  // Constructor.
  Current(const std::string                            &mesh_file_name_,
       const std::string                               &output_file_name_,
       const unsigned int                              &r_,
       SolverAdapter                                   *solver_,
       const double                                    &T_,
       const double                                    &theta_,
       const double                                    &delta_t_)
    : mesh_file_name(mesh_file_name_)
    , output_file_name(output_file_name_)
    , r(r_)
    , T(T_)
    , theta(theta_)
    , delta_t(delta_t_)
    , mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
    , mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    , mesh(MPI_COMM_WORLD)
    , solver(solver_)
    , pcout(std::cout, mpi_rank == 0)
  {}

  void set_model(const std::string& model);

  // Run the time-dependent simulation.
  void
  run();

protected:
  // Initialization.
  void
  setup();

  // System assembly.
  void
  assemble();

  // System solution.
  void
  solve_linear_system();

  void integrate_auxiliar_variables();

  void compute_ionic_currents();

  void check_activation_time();

  // Output.
  void
  output() const;

  void
  output_activation_time() const;

  // Name of the mesh.
  const std::string mesh_file_name;
  const std::string output_file_name;

  // Polynomial degree.
  const unsigned int r;

  // Final time.
  const double T;

  // Theta parameter for the theta method.
  const double theta;

  // Time step.
  const double delta_t;

  // Current time in ms.
  double time = 0.0;

  bool has_non_activated_cells = true;

  // Current timestep number.
  unsigned int timestep_number = 0;

  // // Forcing term.
  // std::function<double(const Point<dim> &, const double &)> f;

  // Number of MPI processes.
  const unsigned int mpi_size;

  // Rank of the current MPI process.
  const unsigned int mpi_rank;

  // Triangulation.
  parallel::fullydistributed::Triangulation<dim> mesh;

  // Finite element space.
  std::unique_ptr<FiniteElement<dim>> fe;

  // Quadrature formula.
  std::unique_ptr<Quadrature<dim>> quadrature;

  // DoF handler.
  DoFHandler<dim> dof_handler;

  // System matrix.
  TrilinosWrappers::SparseMatrix system_matrix;

  //Sistem solver.
  std::unique_ptr<SolverAdapter> solver;

  // System right-hand side.
  TrilinosWrappers::MPI::Vector system_rhs;

  // System solution, without ghost elements.
  TrilinosWrappers::MPI::Vector v;

  // System solution, with ghost elements.
  TrilinosWrappers::MPI::Vector v_ghost;

  //
  TrilinosWrappers::MPI::Vector activation_time;

  // Currents, without ghost elements.
  TrilinosWrappers::MPI::Vector J_fi;
  TrilinosWrappers::MPI::Vector J_so;
  TrilinosWrappers::MPI::Vector J_si;

  // Currents, with ghost elements.
  TrilinosWrappers::MPI::Vector J_fi_ghost;
  TrilinosWrappers::MPI::Vector J_so_ghost;
  TrilinosWrappers::MPI::Vector J_si_ghost;

  // Auxiliar variables, without ghost elements.
  TrilinosWrappers::MPI::Vector w1;
  TrilinosWrappers::MPI::Vector w2;
  TrilinosWrappers::MPI::Vector w3;

  // Output stream for process 0.
  ConditionalOStream pcout;

  
  // --- Valori limite del voltaggio ---
  double v_o;           // voltaggio a riposo (adimensionale)
  double v_u;           // voltaggio massimo upstroke

  // --- Soglie per le funzioni di Heaviside ---
  double theta_w1;       // soglia per J_fi e gate w1
  double theta_w2;       // soglia per J_so, J_si, gate w2
  double theta_w1_m;     // soglia per tau_w1- (quale ramo)
  double theta_o;       // soglia per tau_o e w2_inf

  // --- Parametri per tau_w1- (costante di tempo w1 in chiusura) ---
  double tau_w1_1_m;       // tau_w1 - quando u < theta_w1_m
  double tau_w1_2_m;       // tau_w1 - quando u > theta_w1_m

  // --- Parametri per tau_w1+ (costante di tempo w1 in apertura) ---
  double tau_w1_p;

  // --- Parametri per tau_w2- (costante di tempo w2 in chiusura) ---
  double tau_w2_1_m;      // valore minimo di tau_w2-
  double tau_w2_2_m;      // valore massimo di tau_w2-
  double k_w2_m;         // slope della sigmoide per tau_w2-
  double v_w2_m;         // punto di mezzo della sigmoide

  // --- Parametri per tau_w2+ (costante di tempo w2 in apertura) ---
  double tau_w2_p;

  // --- Parametri per J_fi (fast inward - sodio) ---
  double tau_fi;

  // --- Parametri per J_so (slow outward - potassio) ---
  double tau_o1;        // tau_o quando u < theta_o
  double tau_o2;        // tau_o quando u > theta_o
  double tau_so1;       // tau_so minimo
  double tau_so2;       // tau_so massimo
  double k_so;          // slope della sigmoide per tau_so
  double v_so;          // punto di mezzo della sigmoide

  // --- Parametri per w3 (quarta variabile, morfologia AP) ---
  double tau_w3_1;        // tau_w3 quando v < theta_w2
  double tau_w3_2;        // tau_w3 quando v > theta_w2
  double k_w3;           // slope della tanh per w3_inf
  double v_w3;           // punto di mezzo della tanh

  // --- Parametri per J_si (slow inward - calcio) ---
  double tau_si;

  // --- Parametri per w2_inf ---
  double tau_w2_inf;     // usato nel calcolo di w2_inf
  double w2_inf_star;    // valore di w2_inf quando v > theta_o

  // --- Diffusione ---
  // // D = 1.171 cm^2/s dal paper (Appendice A)
  // // Qui in unità adimensionali del modello
  // double D = 0.1171;        // adattato alle unità del problema

  
  static constexpr double beta_mm = 140.0; // surface-to-volume ratio in mm^-1
  static constexpr double Cm = 0.01e-6;    // membrane capacitance in F/mm^2 (0.01 μF/mm^2)

  static constexpr double sigma_i_long = 0.17;  // intracellular longitudinal conductivity S/m
  static constexpr double sigma_e_long = 0.62;  // extracellular longitudinal conductivity S/m
  static constexpr double sigma_i_trans = 0.019; // intracellular transverse conductivity S/m
  static constexpr double sigma_e_trans = 0.24;  // extracellular transverse conductivity S/m
  
  // per modello anisotropico
  Tensor<2, 3> diffusion_tensor;
  
  // ============================================================
};

#endif
