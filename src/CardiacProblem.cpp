#include "../include/CardiacProblem.hpp"
#include <memory> // Inserito per supportare la gestione degli smart pointers (std::shared_ptr)

namespace CardiacProject {
    using namespace dealii;

    // =========================================================================
    // COSTRUTTORE DELLA CLASSE
    // Inizializza l'infrastruttura di calcolo distribuito MPI e i parametri temporali
    // =========================================================================
    template <int dim>
    CardiacProblem<dim>::CardiacProblem(MeshMode mode, const std::string &mesh_filename)
        : current_mode(mode),                       
          external_mesh_filename(mesh_filename),    
          mpi_communicator(MPI_COMM_WORLD), // Inizializzazione del contesto di comunicazione globale MPI
          triangulation(mpi_communicator),  // La mesh parallela distribuita è legata al comunicatore
          fe(1),                            // Elementi finiti lineari di Lagrange (Q1)
          dof_handler(triangulation),       // Gestore dei Gradi di Libertà legato alla mesh parallela
          pcout(std::cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)), // Stream di log attivo SOLO sul nodo Master (Processo 0)
          computing_timer(mpi_communicator, pcout, TimerOutput::summary,
                          TimerOutput::wall_times) // Profiler HPC per il tracciamento dei colli di bottiglia temporali
    {
        // Parametri di integrazione temporale del problema accoppiato reazione-diffusione
        time = 0.0;
        time_step = 0.001; // ms - Passo temporale ridotto per garantire la stabilità numerica delle ODE rigide
        final_time = 80.0; // ms - Orizzonte temporale sufficiente per osservare l'intero fronte di propagazione
    }

    // =========================================================================
    // SETUP DELLA GRIGLIA COMPUTAZIONALE
    // Gestisce il binario di produzione (mesh esterna .msh) o il benchmark (mesh interna)
    // =========================================================================
    template <int dim>
    void CardiacProblem<dim>::setup_grid() {
        TimerOutput::Scope t(computing_timer, "Setup Grid");

        if (current_mode == MeshMode::ExternalLoad) {
            // -----------------------------------------------------------------
            // MODALITÀ PRODUZIONE: Caricamento mesh anatomica da file esterno
            // -----------------------------------------------------------------
            pcout << "=> Modalità Produzione: Caricamento mesh ottimizzata dal file " 
                  << external_mesh_filename << " ..." << std::endl;
            
            GridIn<dim> grid_in;
            grid_in.attach_triangulation(triangulation);
            std::ifstream input_file(external_mesh_filename);
            
            AssertThrow(input_file.is_open(), 
                        ExcMessage("Errore: impossibile aprire il file della mesh!"));
            
            grid_in.read_msh(input_file); // Lettura del formato Gmsh in modalità parallela distribuita
            pcout << "Mesh esterna caricata con successo!" << std::endl;

        } else {
            // -----------------------------------------------------------------
            // MODALITÀ BENCHMARK: Generazione analitica del dominio rettangolare
            // -----------------------------------------------------------------
            pcout << "=> Modalità Benchmark: Generazione interna della mesh in corso..." << std::endl;
            
            // h definisce il passo di discretizzazione spaziale iniziale (0.2 mm = risoluzione standard)
            double h = 0.2; // h=0.05 fa crashare il programma

            std::vector<unsigned int> repetitions(dim);
            Point<dim> p1, p2;

            // Definizione dei confini fisici del tessuto cardiaco in millimetri (20x7x3 mm)
            p1[0] = 0.0; p1[1] = 0.0; 
            p2[0] = 20.0; p2[1] = 7.0;
            if (dim == 3) {
                p1[2] = 0.0; p2[2] = 3.0;
            }

            // Calcolo automatico delle suddivisioni dei blocchi FEM in base ad h
            for(unsigned int d=0; d<dim; ++d)
                repetitions[d] = static_cast<unsigned int>((p2[d]-p1[d])/h);

            GridGenerator::subdivided_hyper_rectangle(triangulation, repetitions, p1, p2);

            // Assegnazione dei tag del materiale per identificare geometricamente l'area dello stimolo elettrico
            for (const auto &cell : triangulation.active_cell_iterators()) {
                if (cell->is_locally_owned()) {
                    Point<dim> center = cell->center();
                    bool in_stimulus_zone = true;
                    // Lo stimolo iniziale viene applicato in un volume confinato nell'origine (cubo di lato 1.5 mm)
                    for(unsigned int d=0; d<dim; ++d)
                        if(center[d] > 1.5) in_stimulus_zone = false;

                    if (in_stimulus_zone) cell->set_material_id(1); // ID 1: Tessuto stimolato/attivo
                    else cell->set_material_id(0);                  // ID 0: Tessuto passivo a riposo
                }
            }
        }
    }

    // =========================================================================
    // SETUP DELLE STRUTTURE ALGEBRICHE DISTRIBUITE
    // Alloca lo spazio in memoria RAM per i vettori e le matrici MPI parallele
    // =========================================================================
    template <int dim>
    void CardiacProblem<dim>::setup_system() {
        TimerOutput::Scope t(computing_timer, "Setup System");
            
            // Distribuzione globale dei gradi di libertà sui nodi della mesh parallela
            dof_handler.distribute_dofs(fe);

            // Estrazione degli insiemi di DoF per la segmentazione della memoria MPI
            locally_owned_dofs = dof_handler.locally_owned_dofs(); // DoF fisicamente gestiti da questo core
            locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler); // DoF locali + Nodi Ghost di confine
            
            // Log diagnostico sul partizionamento del bilanciamento del carico tra i processi MPI
            pcout << "Process " << Utilities::MPI::this_mpi_process(mpi_communicator) 
                  << ": owned=" << locally_owned_dofs.n_elements() 
                  << ", relevant=" << locally_relevant_dofs.n_elements()
                  << ", ghost=" << (locally_relevant_dofs.n_elements() - locally_owned_dofs.n_elements())
                  << std::endl;

            pcout << "Number of active cells: " << triangulation.n_global_active_cells() << std::endl;
            pcout << "Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl;

            // Inizializzazione delle matrici di vincolo per garantire la continuità (es. nodi hanging)
            constraints.clear();
            constraints.reinit(locally_relevant_dofs);
            DoFTools::make_hanging_node_constraints(dof_handler, constraints);
            constraints.close();

            // ==============================================================
            // SOLUZIONE DEL BUG DI ISOLAMENTO PARALLELO: ALLOCAZIONE NATIVA TRILINOS
            // Inizializziamo lo SparsityPattern legandolo direttamente ai DoF owned
            // e al comunicatore MPI. Evitiamo classi seriali intermedie che 
            // troncherebbero i canali di interscambio e accoppiamento inter-processore.
            // ==============================================================
            TrilinosWrappers::SparsityPattern sp(locally_owned_dofs, mpi_communicator);
            DoFTools::make_sparsity_pattern(dof_handler, sp, constraints, false);
            sp.compress(); // Sincronizzazione di rete per congelare i blocchi off-diagonal dei nodi ghost

            // Allocazione delle matrici sparse parallele basate sulla topologia MPI validata
            system_matrix.reinit(sp);
            mass_matrix.reinit(sp);
            
            // Allocazione dei vettori distribuiti. solution e system_rhs allocano solo i DoF locali (owned)
            solution.reinit(locally_owned_dofs, mpi_communicator);
            system_rhs.reinit(locally_owned_dofs, mpi_communicator);
            
            // locally_relevant_solution alloca anche i nodi ghost necessari per scambiare i dati di bordo
            locally_relevant_solution.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

            // Allocazione delle strutture di memoria per le variabili ioniche contigue della ODE
            unsigned int n_local_dofs = locally_owned_dofs.n_elements();
            gate_v.resize(n_local_dofs, 1.0); // v inizia a 1.0 (stato di riposo/pronto all'attivazione)
            gate_w.resize(n_local_dofs, 1.0); // w inizia a 1.0
            gate_s.resize(n_local_dofs, 0.0); // s inizia a 0.0 (canali del calcio chiusi)

            // Vettore diagnostico di post-processing distribuito per mappare i tempi di attivazione
            activation_time.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
            activation_time = -1.0; // Inizializzazione a -1.0 (nodo non ancora attivato)
    }

    // =========================================================================
    // ASSEMBLAGGIO DEI SISTEMI LINEARI GLOBALO
    // Sfrutta il Mass Lumping per ottimizzare l'algebra lineare parallela
    // =========================================================================
    template <int dim>
    void CardiacProblem<dim>::assemble_system()
    {
        TimerOutput::Scope t(computing_timer, "Assembly");

        // Lucchetti statici per evitare di ricalcolare le matrici ad ogni istante temporale
        static bool matrices_are_assembled = false;
        static TrilinosWrappers::MPI::Vector stimulus_vector;

        // ---------------------------------------------------------------------
        // FASE 1: SETUP UNA-TANTUM (Eseguito esclusivamente a t = 0)
        // ---------------------------------------------------------------------
        if (!matrices_are_assembled) {
            system_matrix = 0.0;
            mass_matrix = 0.0;
            
            stimulus_vector.reinit(locally_owned_dofs, mpi_communicator);
            stimulus_vector = 0.0;

            QGauss<dim> quadrature_formula(fe.degree + 1); // Punti di quadratura gaussiani
            FEValues<dim> fe_values(fe, quadrature_formula,
                                    update_values | update_gradients | update_quadrature_points | update_JxW_values);

            const unsigned int dofs_per_cell = fe.dofs_per_cell;
            const unsigned int n_q_points = quadrature_formula.size();

            FullMatrix<double> cell_system_matrix(dofs_per_cell, dofs_per_cell);
            FullMatrix<double> cell_mass_matrix(dofs_per_cell, dofs_per_cell);
            Vector<double> cell_stimulus(dofs_per_cell);
            std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

            // Costanti fisiche conducibilità del tessuto anisotropo (mm^2/ms)
            const double diff_l = 0.12;  // Conducibilità lungo la fibra cardidaca
            const double diff_t = 0.013; // Conducibilità trasversale alla fibra
            const double Stim_val = 0.5; // Intensità normalizzata della corrente di stimolo

            Tensor<2, dim> diffusion;
            diffusion[0][0] = diff_l;
            diffusion[1][1] = diff_t;
            if (dim == 3) diffusion[2][2] = diff_t;

            unsigned int n_cells_assembled = 0;

            // Loop globale sugli elementi finiti locali assegnati al core corrente
            for (const auto &cell : dof_handler.active_cell_iterators()) {
                if (cell->is_locally_owned()) {
                    n_cells_assembled++;
                    cell_system_matrix = 0;
                    cell_mass_matrix = 0;
                    cell_stimulus = 0;
                    
                    fe_values.reinit(cell);

                    for (unsigned int q = 0; q < n_q_points; ++q) {
                        const double JxW = fe_values.JxW(q); // Determinante Jacobiano * Peso di quadratura
                        Point<dim> q_point = fe_values.quadrature_point(q);
                        
                        bool in_stim = true;
                        for(unsigned int d=0; d<dim; ++d) 
                            if(q_point[d] > 1.5) in_stim = false;

                        for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                            const auto phi_i = fe_values.shape_value(i, q); // Funzione di forma i
                            
                            // -------------------------------------------------------------
                            // IMPLEMENTAZIONE MASS LUMPING DIAGONALE
                            // Ancoriamo la massa esclusivamente sull'indice principale (i, i)
                            // Questo elimina numericamente gli undershoots e stabilizza l'onda
                            // -------------------------------------------------------------
                            double mass_lumped = (1.0 / time_step) * phi_i * JxW;
                            cell_system_matrix(i, i) += mass_lumped;
                            cell_mass_matrix(i, i) += mass_lumped;

                            const auto grad_phi_i = fe_values.shape_grad(i, q);
                            for (unsigned int j = 0; j < dofs_per_cell; ++j) {
                                const auto grad_phi_j = fe_values.shape_grad(j, q);
                                double stiff_term = grad_phi_i * (diffusion * grad_phi_j);
                                
                                // La matrice di rigidezza spaziale (diffusione) rimane consistente
                                cell_system_matrix(i, j) += stiff_term * JxW;
                            }

                            if (in_stim) {
                                cell_stimulus(i) += Stim_val * phi_i * JxW;
                            }
                        }
                    }
                    cell->get_dof_indices(local_dof_indices);
                    
                    // Trasferimento dei contributi locali nei contenitori globali distribuiti MPI
                    constraints.distribute_local_to_global(cell_system_matrix, local_dof_indices, system_matrix);
                    constraints.distribute_local_to_global(cell_mass_matrix, local_dof_indices, mass_matrix);
                    constraints.distribute_local_to_global(cell_stimulus, local_dof_indices, stimulus_vector);
                }
            }
            // Compressione algebrica finale e somme parziali inter-processore per le matrici statiche
            system_matrix.compress(VectorOperation::add);
            mass_matrix.compress(VectorOperation::add);
            stimulus_vector.compress(VectorOperation::add);
            
            matrices_are_assembled = true;
            if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0) {
                pcout << "=> Matrici Pre-Assemblate. MASS LUMPING e HPC OPTIMIZATION!" << std::endl;
            }
        }

        // ---------------------------------------------------------------------
        // FASE 2: EFFETTIVA PRODUZIONE TEMPORALE (Eseguita 40.000 volte)
        // Ottimizzazione vettoriale: zero cicli interni e accessi RAM sequenziali
        // ---------------------------------------------------------------------
        system_rhs = 0.0;
        
        // Operazione parallela bloccante ottimizzata: system_rhs = mass_matrix * solution
        // Essendo mass_matrix lumpata (diagonale pura), questa vmult è fulminea
        mass_matrix.vmult(system_rhs, solution); 
        
        // Applicazione temporale del vettore di stimolo (durata transitoria = 2.0 ms)
        if (time <= 2.0) {
            system_rhs.add(1.0, stimulus_vector);
        }
    }
    
    // =========================================================================
    // RISOLUTORE DELLA EQUAZIONE ALLE DERIVATE PARZIALI (PDE)
    // =========================================================================
    template <int dim>
    void CardiacProblem<dim>::solve_pde()
    {
        TimerOutput::Scope t(computing_timer, "Solve PDE");

        // Solutore iterativo a Gradienti Coniugati (CG), ottimale per matrici simmetriche definite positive
        SolverControl solver_control(10000, 1e-6 * system_rhs.l2_norm());
        TrilinosWrappers::SolverCG solver(solver_control);
        
        // ==============================================================
        // SCELTA STRATEGICA PRECONDIZIONATORE: JACOBI DIAGONALE PURE
        // Grazie al Mass Lumping e al dt ridotto, la matrice è fortemente
        // a diagonale dominante. Jacobi non richiede alcuna comunicazione 
        // di rete inter-processore, azzerando la latenza MPI latente dell'AMG.
        // ==============================================================
        static std::shared_ptr<TrilinosWrappers::PreconditionJacobi> preconditioner;
        
        if (!preconditioner) {
            if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
                pcout << "=> Inizializzazione Precondizionatore Jacobi..." << std::endl;

            preconditioner = std::make_shared<TrilinosWrappers::PreconditionJacobi>();
            preconditioner->initialize(system_matrix);
        }
        
        // Risoluzione parallela distribuita del sistema lineare PDE
        solver.solve(system_matrix, solution, system_rhs, *preconditioner);

        // Distribuzione e sincronizzazione dei vincoli tra le partizioni di memoria
        constraints.distribute(solution);
        solution.compress(VectorOperation::insert);
    }

    // =========================================================================
    // RISOLUTORE DELLE EQUAZIONI DIFFERENZIALI ORDINARIE (ODE - CANALI IONICI)
    // Sfrutta il metodo di Forward Euler VETTORIALIZZATO (AVX2/SIMD)
    // =========================================================================
    template <int dim>
    void CardiacProblem<dim>::solve_ode()
    {
        TimerOutput::Scope t(computing_timer, "Solve ODE");
                
        // LARGHEZZA DEL REGISTRO SIMD (es. 4 per processori AVX2)
        constexpr unsigned int vector_width = VectorizedArray<double>::size();
        
        unsigned int local_idx = 0;
        IndexSet::ElementIterator idx = locally_owned_dofs.begin();
        IndexSet::ElementIterator endIdx = locally_owned_dofs.end();
        unsigned int n_local = locally_owned_dofs.n_elements();

        // ---------------------------------------------------------------------
        // 1. CICLO VETTORIALE PRINCIPALE (Avanza a blocchi di 'vector_width')
        // ---------------------------------------------------------------------
        for(; local_idx + vector_width <= n_local; local_idx += vector_width) 
        {
            VectorizedArray<double> u_vec, v_vec, w_vec, s_vec;
            std::vector<types::global_dof_index> global_indices(vector_width);

            // FASE DI LOAD: Riempiamo i registri SIMD
            for (unsigned int v = 0; v < vector_width; ++v) {
                global_indices[v] = *idx;
                u_vec[v] = solution(*idx);
                v_vec[v] = gate_v[local_idx + v];
                w_vec[v] = gate_w[local_idx + v];
                s_vec[v] = gate_s[local_idx + v];
                ++idx;
            }

            IonicState<VectorizedArray<double>> state = {u_vec, v_vec, w_vec, s_vec};

            // FASE DI COMPUTE (Tutta la chimica di Bueno-Orovio per 4 nodi in un colpo solo)
            VectorizedArray<double> J_ion = ionic_model.compute_total_current(state);

            // Salvavita numerico SIMD (clipping tra -100 e 100 usando max/min vettoriali)
            J_ion = std::max(VectorizedArray<double>(-100.0), std::min(VectorizedArray<double>(100.0), J_ion));

            state.u -= VectorizedArray<double>(time_step) * J_ion; 
            state.u = std::max(VectorizedArray<double>(0.0), std::min(VectorizedArray<double>(1.6), state.u)); 

            ionic_model.evolve_gates(state, time_step);
            
            state.v = std::max(VectorizedArray<double>(0.0), std::min(VectorizedArray<double>(1.0), state.v));
            state.w = std::max(VectorizedArray<double>(0.0), std::min(VectorizedArray<double>(1.0), state.w));
            state.s = std::max(VectorizedArray<double>(0.0), std::min(VectorizedArray<double>(1.0), state.s));

            // FASE DI STORE: Scarichiamo i risultati in RAM
            for (unsigned int v = 0; v < vector_width; ++v) {
                solution(global_indices[v]) = state.u[v];
                gate_v[local_idx + v] = state.v[v];
                gate_w[local_idx + v] = state.w[v];
                gate_s[local_idx + v] = state.s[v];
            }
        }

        // ---------------------------------------------------------------------
        // 2. TAIL-LOOP SEQUENZIALE
        // Calcola in modo scalare gli ultimissimi nodi (resto della divisione intera)
        // ---------------------------------------------------------------------
        for(; local_idx < n_local; ++local_idx, ++idx) 
        {
            types::global_dof_index global_i = *idx;
            
            IonicState<double> state;
            state.u = solution(global_i);
            state.v = gate_v[local_idx];
            state.w = gate_w[local_idx];
            state.s = gate_s[local_idx];

            double J_ion = ionic_model.compute_total_current(state);
            if (std::abs(J_ion) > 100.0) J_ion = (J_ion > 0 ? 100.0 : -100.0);

            state.u -= time_step * J_ion; 
            state.u = std::max(0.0, std::min(1.6, state.u)); 

            ionic_model.evolve_gates(state, time_step);
            
            state.v = std::max(0.0, std::min(1.0, state.v));
            state.w = std::max(0.0, std::min(1.0, state.w));
            state.s = std::max(0.0, std::min(1.0, state.s));

            solution(global_i) = state.u;
            gate_v[local_idx] = state.v;
            gate_w[local_idx] = state.w;
            gate_s[local_idx] = state.s;
        }
        
        solution.compress(VectorOperation::insert);
    }

    // =========================================================================
    // CICLO DI EVOLUZIONE TEMPORALE (IL DRIVER PRINCIPALE)
    // Coordinatore globale della pipeline numerica parallela
    // =========================================================================
    template <int dim>
    void CardiacProblem<dim>::run()
    {
        pcout << "Running Cardiac Problem on " 
            << Utilities::MPI::n_mpi_processes(mpi_communicator) 
            << " processes." << std::endl;

        setup_grid();
        setup_system();

        // Condizioni iniziali a potenziale nullo ovunque (interpolazione analitica)
        VectorTools::interpolate(dof_handler, Functions::ZeroFunction<dim>(), solution);
        locally_relevant_solution = solution;  

        output_results(0); // Scrittura del frame di baseline a t=0

        unsigned int time_step_number = 0;
        const unsigned int total_steps = static_cast<unsigned int>(final_time / time_step);
        unsigned int last_progress = 0;
        
        if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
            std::cout << "Progress: 0%" << std::flush;
        
        // ---------------------------------------------------------------------
        // MAIN LOOPS DI SIMULAZIONE
        // ---------------------------------------------------------------------
        while (time < final_time)
        {
            time_step_number++;
            time += time_step;
            
            // Catena dell'Operator Splitting spaziotemporale
            assemble_system();
            solve_pde(); // Fase diffusiva macroscopica
            solve_ode(); // Fase reattiva microscopica cellulare
            
            // NOTA OPTIMIZATION: Sincronizzazione barriera locally_relevant_solution rimossa da qui 
            // per evitare di intasare la rete 40.000 volte inutilmente.
            
            // Logging testuale periodico dello stato dell'upstroke
            if (time_step_number % 100 == 0 && Utilities::MPI::this_mpi_process(mpi_communicator) == 0) {
                std::cout << "Step " << time_step_number << ", t=" << time << "ms" << std::endl;
            }
            
            // Aggiornamento asincrono della barra di progressione visiva su terminale master
            if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0) {
                unsigned int progress = (time_step_number * 100) / total_steps;
                if (progress >= last_progress + 5) {
                    std::cout << "Progress: " << progress << "%" << std::endl;
                    last_progress = progress;
                }
            }
            
            // -----------------------------------------------------------------
            // HACK CACHE MEMORIA PER ACTIVATION TIME
            // Scansione diretta basata sull'array continuo dei DoF di proprietà,
            // evitando l'oneroso ciclo sugli iteratori attivi delle celle geometriche
            // -----------------------------------------------------------------
            const double threshold = 0.5; // Soglia convenzionale di superamento della depolarizzazione
            for (IndexSet::ElementIterator it = locally_owned_dofs.begin(); it != locally_owned_dofs.end(); ++it) {
                types::global_dof_index idx = *it;
                if (activation_time(idx) == -1.0 && solution(idx) > threshold) {
                    activation_time(idx) = time; // Registrazione cronometrica invariabile
                }
            }
            activation_time.compress(VectorOperation::insert);

            // -----------------------------------------------------------------
            // EFFETTIVA ESPORTAZIONE DATI SU DISCO (GUEST SYNC ONLY FOR PARAVIEW)
            // Limitiamo l'I/O ad intervalli di 100 passi per non strozzare la CPU
            // -----------------------------------------------------------------
            if (time_step_number % 100 == 0) {
                locally_relevant_solution = solution; // Unica sincronizzazione barriera ammessa
                output_results(time_step_number);
            }
        } // <-- Fine del ciclo temporale while

        if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
            std::cout << "\rProgress: 100% - Simulation completed!" << std::endl;
    }
    
    // =========================================================================
    // EXPORT RISULTATI PARALLELI IN FORMATO VTK/PARAVIEW (.vtu / .pvtu)
    // =========================================================================
    template <int dim>
    void CardiacProblem<dim>::output_results(const unsigned int step) {
        DataOut<dim> data_out;
        data_out.attach_dof_handler(dof_handler);

        // Vettori esportati completi di accoppiamento dei nodi ghost per prevenire artefatti grafici
        data_out.add_data_vector(locally_relevant_solution, "V");
        data_out.add_data_vector(activation_time, "ActivationTime");
        data_out.build_patches();

        // Generazione del file record master .pvtu che cuce insieme i frammenti .vtu locali dei processi
        data_out.write_vtu_with_pvtu_record("./", "solution", step, mpi_communicator);
    }

    // Istanziazioni esplicite per le dimensioni geometriche supportate dal binario
    template class CardiacProblem<2>;
    template class CardiacProblem<3>;
}
