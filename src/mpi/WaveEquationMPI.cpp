#include "WaveEquationMPI.hpp"

#include <deal.II/base/function_lib.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/numbers.h>
#include <deal.II/grid/grid_generator.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace { constexpr double pi = numbers::PI; }

class MPIExact : public Function<2>
{
public:
  explicit MPIExact(double w) : omega(w) {}
  double value(const Point<2> &p, unsigned int = 0) const override
  { return std::sin(pi*p[0])*std::sin(pi*p[1])*std::cos(omega*get_time()); }
  Tensor<1,2> gradient(const Point<2> &p, unsigned int = 0) const override
  {
    Tensor<1,2> g;
    g[0]=pi*std::cos(pi*p[0])*std::sin(pi*p[1])*std::cos(omega*get_time());
    g[1]=pi*std::sin(pi*p[0])*std::cos(pi*p[1])*std::cos(omega*get_time());
    return g;
  }
private: double omega;
};

WaveEquationMPI::WaveEquationMPI(unsigned int n_, unsigned int degree_,
                                 double T, double dt_, unsigned int frequency,
                                 std::string directory, std::string case_name_)
  : n(n_), degree(degree_), output_frequency(frequency), final_time(T), dt(dt_),
    output_directory(std::move(directory)), case_name(std::move(case_name_)),
    rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)),
    n_ranks(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)), mesh(MPI_COMM_WORLD),
    fe(degree), quadrature(degree+1), dof_handler(mesh),
    pcout(std::cout, rank == 0)
{}

double WaveEquationMPI::omega() const
{ return case_name == "forced" ? pi : std::sqrt(2.)*pi; }

void WaveEquationMPI::setup()
{
  // deal.II partitions this structured grid directly across MPI ranks.
  GridGenerator::subdivided_hyper_cube(mesh, n, 0., 1.);
  dof_handler.distribute_dofs(fe);
  owned = dof_handler.locally_owned_dofs();
  relevant = DoFTools::extract_locally_relevant_dofs(dof_handler);
  constraints.clear();
  VectorTools::interpolate_boundary_values(dof_handler, 0,
                                           Functions::ZeroFunction<dim>(),
                                           constraints);
  constraints.close();

  sparsity.reinit(owned, MPI_COMM_WORLD);
  DoFTools::make_sparsity_pattern(dof_handler, sparsity, constraints, false);
  sparsity.compress();
  mass.reinit(sparsity); stiffness.reinit(sparsity); effective.reinit(sparsity);
  u_owned.reinit(owned, MPI_COMM_WORLD); v_owned.reinit(owned, MPI_COMM_WORLD);
  a_owned.reinit(owned, MPI_COMM_WORLD); rhs.reinit(owned, MPI_COMM_WORLD);
  u.reinit(owned, relevant, MPI_COMM_WORLD);
  v.reinit(owned, relevant, MPI_COMM_WORLD);
  a.reinit(owned, relevant, MPI_COMM_WORLD);
}

void WaveEquationMPI::assemble_matrices()
{
  FEValues<dim> values(fe, quadrature, update_values|update_gradients|update_JxW_values);
  FullMatrix<double> cm(fe.dofs_per_cell, fe.dofs_per_cell);
  FullMatrix<double> ck(fe.dofs_per_cell, fe.dofs_per_cell);
  std::vector<types::global_dof_index> indices(fe.dofs_per_cell);
  mass=0.; stiffness=0.;
  for (const auto &cell : dof_handler.active_cell_iterators()) if (cell->is_locally_owned())
    {
      values.reinit(cell); cm=0.; ck=0.;
      for (unsigned int q=0; q<quadrature.size(); ++q)
        for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
          for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
            {
              cm(i,j)+=values.shape_value(i,q)*values.shape_value(j,q)*values.JxW(q);
              ck(i,j)+=scalar_product(values.shape_grad(i,q),values.shape_grad(j,q))*values.JxW(q);
            }
      cell->get_dof_indices(indices);
      constraints.distribute_local_to_global(cm, indices, mass);
      constraints.distribute_local_to_global(ck, indices, stiffness);
    }
  mass.compress(VectorOperation::add); stiffness.compress(VectorOperation::add);
  effective.copy_from(mass); effective.add(.25*dt*dt, stiffness);
  for (const auto index : owned)
    if (constraints.is_constrained(index))
      { mass.set(index, index, 1.); effective.set(index, index, 1.); }
  mass.compress(VectorOperation::insert); effective.compress(VectorOperation::insert);
}

void WaveEquationMPI::assemble_force(double t, TrilinosWrappers::MPI::Vector &force) const
{
  force=0.; if (case_name=="free") return;
  FEValues<dim> values(fe, quadrature, update_values|update_quadrature_points|update_JxW_values);
  Vector<double> local(fe.dofs_per_cell);
  std::vector<types::global_dof_index> indices(fe.dofs_per_cell);
  const double c=2.*pi*pi-omega()*omega();
  for (const auto &cell : dof_handler.active_cell_iterators()) if (cell->is_locally_owned())
    {
      values.reinit(cell); local=0.;
      for (unsigned int q=0;q<quadrature.size();++q)
        {
          const auto &p=values.quadrature_point(q);
          const double f=c*std::sin(pi*p[0])*std::sin(pi*p[1])*std::cos(omega()*t);
          for (unsigned int i=0;i<fe.dofs_per_cell;++i)
            local(i)+=f*values.shape_value(i,q)*values.JxW(q);
        }
      cell->get_dof_indices(indices); constraints.distribute_local_to_global(local,indices,force);
    }
  force.compress(VectorOperation::add);
}

void WaveEquationMPI::initialize()
{
  MPIExact exact(omega()); exact.set_time(0.);
  VectorTools::interpolate(dof_handler, exact, u_owned); constraints.distribute(u_owned);
  v_owned=0.; solve_initial_acceleration(); u=u_owned; v=v_owned; a=a_owned;
}

void WaveEquationMPI::solve_initial_acceleration()
{
  TrilinosWrappers::MPI::Vector force(owned,MPI_COMM_WORLD), ku(owned,MPI_COMM_WORLD);
  assemble_force(0.,force); stiffness.vmult(ku,u_owned); rhs=force; rhs-=ku;
  SolverControl control(5000,1.e-12); SolverCG<TrilinosWrappers::MPI::Vector> solver(control);
  TrilinosWrappers::PreconditionSSOR preconditioner; preconditioner.initialize(mass);
  solver.solve(mass,a_owned,rhs,preconditioner); constraints.distribute(a_owned);
}

void WaveEquationMPI::solve_step(double t)
{
  TrilinosWrappers::MPI::Vector force(owned,MPI_COMM_WORLD), predictor(u_owned), ka(owned,MPI_COMM_WORLD), next(owned,MPI_COMM_WORLD);
  assemble_force(t, force);
  predictor.add(dt,v_owned); predictor.add(dt*dt*(.5-.25),a_owned);
  stiffness.vmult(ka,predictor); rhs=force; rhs-=ka;
  SolverControl control(5000,1.e-12); SolverCG<TrilinosWrappers::MPI::Vector> solver(control);
  TrilinosWrappers::PreconditionSSOR preconditioner; preconditioner.initialize(effective);
  solver.solve(effective,next,rhs,preconditioner);
  u_owned=predictor; u_owned.add(.25*dt*dt,next);
  v_owned.add(.5*dt,a_owned); v_owned.add(.5*dt,next); a_owned=next;
  constraints.distribute(u_owned); constraints.distribute(v_owned); constraints.distribute(a_owned);
  u=u_owned; v=v_owned; a=a_owned;
}

double WaveEquationMPI::energy() const
{
  TrilinosWrappers::MPI::Vector mv(owned,MPI_COMM_WORLD), kv(owned,MPI_COMM_WORLD);
  mass.vmult(mv,v_owned); stiffness.vmult(kv,u_owned);
  return .5*(v_owned*mv+u_owned*kv);
}

std::pair<double,double> WaveEquationMPI::errors() const
{
  MPIExact exact(omega()); exact.set_time(current_time);
  QGauss<dim> q(degree+2); Vector<float> cell_errors(mesh.n_active_cells());
  VectorTools::integrate_difference(dof_handler,u,exact,cell_errors,q,VectorTools::L2_norm);
  const double l2=VectorTools::compute_global_error(mesh,cell_errors,VectorTools::L2_norm);
  VectorTools::integrate_difference(dof_handler,u,exact,cell_errors,q,VectorTools::H1_seminorm);
  return {l2,VectorTools::compute_global_error(mesh,cell_errors,VectorTools::H1_seminorm)};
}

void WaveEquationMPI::output()
{
  DataOut<dim> data; data.add_data_vector(dof_handler,u,"u"); data.add_data_vector(dof_handler,v,"velocity"); data.add_data_vector(dof_handler,a,"acceleration"); data.build_patches(fe.degree);
  data.write_vtu_with_pvtu_record(output_directory + "/","wave_solution",current_step,MPI_COMM_WORLD,fe.degree);
}

void WaveEquationMPI::diagnostics()
{
  const auto e=errors();
  const double total_energy = energy();
  if(rank==0) { std::ofstream file(output_directory+"/diagnostics.csv",std::ios::app); file<<std::setprecision(16)<<current_step<<','<<current_time<<','<<e.first<<','<<e.second<<','<<total_energy<<'\n'; }
}

void WaveEquationMPI::run()
{
  const double start_time = MPI_Wtime();
  std::filesystem::create_directories(output_directory);
  setup(); assemble_matrices(); initialize();
  const auto steps=static_cast<unsigned int>(std::llround(final_time/dt));
  current_step=0; current_time=0.;
  if(rank==0) std::ofstream(output_directory+"/diagnostics.csv")
    << "step,time,l2_error,h1_seminorm,energy\n";
  diagnostics();
  output();
  const double setup_elapsed = MPI_Wtime() - start_time;
  const double loop_start = MPI_Wtime();
  for(;current_step<steps;++current_step)
    {
      current_time=(current_step+1)*dt;
      solve_step(current_time);
      diagnostics();
    }
  const double loop_elapsed = MPI_Wtime() - loop_start;
  const double postprocess_start = MPI_Wtime();
  output();
  const double postprocess_elapsed = MPI_Wtime() - postprocess_start;
  const double elapsed = Utilities::MPI::max(MPI_Wtime() - start_time, MPI_COMM_WORLD);
  const double setup_max = Utilities::MPI::max(setup_elapsed, MPI_COMM_WORLD);
  const double loop_max = Utilities::MPI::max(loop_elapsed, MPI_COMM_WORLD);
  const double postprocess_max = Utilities::MPI::max(postprocess_elapsed, MPI_COMM_WORLD);
  if (rank == 0)
    {
      std::ofstream timing(output_directory + "/timing.csv");
      timing << "n_ranks,n_subdivisions,dofs,cells,setup_seconds,time_loop_seconds,"
                "postprocess_seconds,elapsed_seconds\n";
      timing << n_ranks << ',' << n << ',' << dof_handler.n_dofs() << ','
             << mesh.n_global_active_cells() << ',' << std::setprecision(16)
             << setup_max << ',' << loop_max << ',' << postprocess_max << ','
             << elapsed << '\n';
      pcout << "MPI run completed with " << n_ranks << " processes in "
            << elapsed << " seconds\n";
    }
}
