#ifndef CARDIAC_PROBLEM_HPP
#define CARDIAC_PROBLEM_HPP

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/timer.h>

#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/affine_constraints.h>

#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/grid/grid_in.h> // Per leggere le mesh!

#include "IonicModel.hpp"

#include <fstream>
#include <iostream>

namespace CardiacProject {
    using namespace dealii;
    // Le due modalità di gestione della mesh:
    // - InternalGeneration: la mesh viene generata internamente con GridGenerator.
    // - ExternalLoad: la mesh viene caricata da un file esterno.
    enum class MeshMode {
        InternalGeneration,
        ExternalLoad
    };

    // Driver numerico completo del problema elettrofisiologico.
    // La classe incapsula mesh, DoF, assembly, solve PDE/ODE e output parallelo.
    template <int dim>
    class CardiacProblem {
    public:
        // Costruisce il problema e inizializza i dati temporali di default.
        CardiacProblem();

        // Costruisce il problema specificando la modalità di mesh e il nome del file se necessario.
        CardiacProblem(MeshMode mode = MeshMode::InternalGeneration, const std::string &mesh_filename = "");

        // Esegue l'intero ciclo di simulazione: mesh, setup, time stepping e output.
        void run();

    private:
        // Prepara la mesh parallela e marca la regione di stimolo.
        void setup_grid();
        // Distribuisce i DoF e inizializza matrici, vettori e vincoli.
        void setup_system();
        // Assembla il sistema lineare della parte diffusive + termine di reazione.
        void assemble_system();
        // Integra localmente la parte ODE del modello ionico.
        void solve_ode();
        // Risolve il sistema lineare per il potenziale transmembrana.
        void solve_pde();
        // Scrive file VTU/PVTU per visualizzazione in ParaView.
        void output_results(const unsigned int time_step);

        // Variabili per salvare la configurazione scelta
        MeshMode current_mode;
        std::string external_mesh_filename;

        // Precondizionatore AMG per il solver lineare, gestito con smart pointer per sicurezza.
        std::shared_ptr<TrilinosWrappers::PreconditionAMG> amg_preconditioner;

        // Comunicazione MPI e partizionamento della mesh.
        MPI_Comm mpi_communicator;
        parallel::distributed::Triangulation<dim, dim> triangulation;
        FE_Q<dim> fe;
        DoFHandler<dim> dof_handler;

        // DoF posseduti localmente e DoF rilevanti con ghost.
        IndexSet locally_owned_dofs;
        IndexSet locally_relevant_dofs;

        // Vincoli di continuità e gestione dei nodi hanging.
        AffineConstraints<double> constraints;

        // Matrici sparse globali del problema lineare.
        TrilinosWrappers::SparseMatrix system_matrix;
        TrilinosWrappers::SparseMatrix mass_matrix; 
        
        // Vettori MPI: solution contiene solo gli owned DoF, mentre
        // locally_relevant_solution include anche i ghost per la lettura.
        TrilinosWrappers::MPI::Vector solution;  // Solo locally_owned (per solver)
        TrilinosWrappers::MPI::Vector locally_relevant_solution;  // Con ghost elements (per assembly)
        TrilinosWrappers::MPI::Vector system_rhs;
        TrilinosWrappers::MPI::Vector activation_time;

        // Strumenti di logging e profilazione per output parallelo.
        ConditionalOStream pcout;
        TimerOutput computing_timer;

        // Parametri temporali della simulazione.
        double time;
        double time_step;
        double final_time;
        
        // Gate ionici memorizzati localmente per efficienza di accesso.
        std::vector<double> gate_v;
        std::vector<double> gate_w;
        std::vector<double> gate_s;
        IonicModel ionic_model;
    };
}

#endif
