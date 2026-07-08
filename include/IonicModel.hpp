#ifndef IONIC_MODEL_HPP
#define IONIC_MODEL_HPP

#include <algorithm>
#include <cmath>
#include <deal.II/base/vectorization.h> // Indispensabile per le operazioni SIMD

/**
 * @file IonicModel.hpp
 * @brief Implementazione del Modello Minimo Fenomenologico per il potenziale d'azione ventricolare umano.
 * * Basato ESATTAMENTE sulle equazioni e sui parametri (Set Epicardico) della Tabella 1 del paper:
 * "Minimal model for human ventricular action potentials in tissue"
 * Alfonso Bueno-Orovio, Elizabeth M. Cherry, Flavio H. Fenton.
 * Journal of Theoretical Biology 253 (2008) 544-560.
 * * Il modello riduce la complessità dei canali ionici a 3 correnti principali (J_fi, J_so, J_si)
 * governate dal potenziale u e da 3 variabili di gating (v, w, s).
 */

using namespace dealii; // Aggiunto per permettere l'uso diretto di SIMDComparison

// Stato locale templatizzato per accettare sia 'double' singolo sia 'VectorizedArray<double>'
template <typename Number>
struct IonicState {
    Number u; // Potenziale transmembrana (adimensionale)
    Number v; // Gate v (governa la corrente rapida del sodio)
    Number w; // Gate w (governa le correnti lente, recupero)
    Number s; // Gate s (governa la corrente lenta entrante, plateau)
};

class IonicModel {
public:
    IonicModel() {
        // Parametri per cellula EPICARDICA (Rif: Tabella 1, Bueno-Orovio et al. 2008)
        // Soglie di voltaggio
        u_o = 0.0; 
        u_u = 1.55;             // Picco massimo del potenziale
        theta_v = 0.3;          // Soglia di ATTIVAZIONE per la corrente rapida (Sodio)
        theta_w = 0.13;         // Soglia di ATTIVAZIONE per le correnti lente
        theta_v_minus = 0.006;  // Soglia di pre-allarme/inattivazione per il gate v
        theta_o = 0.006;        // Soglia per l'attivazione della corrente uscente di base

        // Costanti di tempo per il gate v (Sodio)
        tau_v1_minus = 60.0; 
        tau_v2_minus = 1150.0;
        tau_v_plus = 1.4506; 

        // Costanti di tempo per il gate w (Potassio/Recupero)
        tau_w1_minus = 60.0; 
        tau_w2_minus = 15.0;
        k_w_minus = 65.0; 
        u_w_minus = 0.03; 
        tau_w_plus = 200.0;
        tau_w_inf = 0.07;       // Parametro w_inf corretto da Tabella 1
        w_inf_star = 0.94;

        // Costanti di tempo per le correnti (J_fi, J_so)
        tau_fi = 0.11;          // Velocità upstroke (più è piccolo, più è ripido)
        tau_o1 = 400.0; 
        tau_o2 = 6.0;
        tau_so1 = 30.0181;      // Valore esatto da Tabella 1
        tau_so2 = 0.9957; 
        k_so = 2.0458; 
        u_so = 0.65;

        // Costanti di tempo per il gate s (Calcio/Plateau)
        tau_s1 = 2.7342; 
        tau_s2 = 16.0; 
        k_s = 2.0994; 
        u_s = 0.9087;
        tau_si = 1.8875; 
    }

    /**
     * @brief Calcola la somma delle correnti ioniche (RHS dell'equazione Monodominio).
     * Funziona ora in modo agnostico per scalari e vettori SIMD.
     */
    template <typename Number>
    Number compute_total_current(const IonicState<Number> &state) const {
        Number u = state.u;
        Number v = state.v;
        Number w = state.w;
        Number s = state.s;

        // 1. J_fi: Maschera SIMD al posto del branching ternario 
        Number H_u_thv = compare_and_apply_mask<SIMDComparison::greater_than>(u, Number(theta_v), Number(1.0), Number(0.0));
        Number J_fi = -(v * H_u_thv * (u - Number(theta_v)) * (Number(u_u) - u)) / Number(tau_fi);

        // 2. J_so
        Number H_u_thw = compare_and_apply_mask<SIMDComparison::greater_than>(u, Number(theta_w), Number(1.0), Number(0.0));
        Number J_so = ((u - Number(u_o)) * (Number(1.0) - H_u_thw) / tau_o(u)) + (H_u_thw / tau_so(u));

        // 3. J_si
        Number J_si = - (H_u_thw * w * s) / Number(tau_si);

        return J_fi + J_so + J_si;
    }

    /**
     * @brief Evolve le variabili di gating nel tempo usando Forward Euler esplicito.
     */
    template <typename Number>
    void evolve_gates(IonicState<Number> &state, double dt) const {
        Number u = state.u;
        Number time_step = Number(dt); // Promozione di dt a Number
        
        // --- Evoluzione GATE V (Sodio) ---
        Number H_u_thv_m = compare_and_apply_mask<SIMDComparison::greater_than>(u, Number(theta_v_minus), Number(1.0), Number(0.0));
        Number v_inf = compare_and_apply_mask<SIMDComparison::less_than>(u, Number(theta_v_minus), Number(1.0), Number(0.0));
        Number tau_v_m = (Number(1.0) - H_u_thv_m) * Number(tau_v1_minus) + H_u_thv_m * Number(tau_v2_minus);
        
        Number H_u_thv = compare_and_apply_mask<SIMDComparison::greater_than>(u, Number(theta_v), Number(1.0), Number(0.0));
        Number dv = (Number(1.0) - H_u_thv) * (v_inf - state.v) / tau_v_m - H_u_thv * state.v / Number(tau_v_plus);
        state.v += time_step * dv;

        // --- Evoluzione GATE W (Potassio) ---
        Number H_u_tho = compare_and_apply_mask<SIMDComparison::greater_than>(u, Number(theta_o), Number(1.0), Number(0.0));
        Number w_inf = (Number(1.0) - H_u_tho) * (Number(1.0) - (u / Number(tau_w_inf))) + H_u_tho * Number(w_inf_star);
        
        // Sostituzione di std::tanh con la nostra helper SIMD-compatibile
        Number tau_w_m = Number(tau_w1_minus) + (Number(tau_w2_minus) - Number(tau_w1_minus)) * (Number(1.0) + custom_tanh(Number(k_w_minus) * (u - Number(u_w_minus)))) / Number(2.0);
        
        Number H_u_thw = compare_and_apply_mask<SIMDComparison::greater_than>(u, Number(theta_w), Number(1.0), Number(0.0));
        Number dw = (Number(1.0) - H_u_thw) * (w_inf - state.w) / tau_w_m - H_u_thw * state.w / Number(tau_w_plus);
        state.w += time_step * dw;

        // --- Evoluzione GATE S (Calcio) ---
        Number tau_s = (Number(1.0) - H_u_thw) * Number(tau_s1) + H_u_thw * Number(tau_s2); 
        Number ds = (((Number(1.0) + custom_tanh(Number(k_s) * (u - Number(u_s)))) * Number(0.5)) - state.s) / tau_s;
        state.s += time_step * ds;
    }

private:
    double u_o, u_u, theta_v, theta_w, theta_v_minus, theta_o;
    double tau_v1_minus, tau_v2_minus, tau_v_plus;
    double tau_w1_minus, tau_w2_minus, k_w_minus, u_w_minus, tau_w_plus;
    double tau_fi, tau_o1, tau_o2, tau_so1, tau_so2, k_so, u_so;
    double tau_s1, tau_s2, k_s, u_s, tau_si;
    double w_inf_star;
    double tau_w_inf; 

    // Funzioni helper aggiornate con i template
    template <typename Number>
    Number tau_o(Number u) const {
        Number H_u_tho = compare_and_apply_mask<SIMDComparison::greater_than>(u, Number(theta_o), Number(1.0), Number(0.0));
        return (Number(1.0) - H_u_tho) * Number(tau_o1) + H_u_tho * Number(tau_o2);
    }

    template <typename Number>
    Number tau_so(Number u) const {
        return Number(tau_so1) + (Number(tau_so2) - Number(tau_so1)) * (Number(1.0) + custom_tanh(Number(k_so)*(u - Number(u_so)))) * Number(0.5);
    }

    // HELPER MATEMATICO PER IL SIMD (Sostituisce std::tanh)
    template <typename Number>
    Number custom_tanh(Number x) const {
        Number exp_pos = std::exp(x);
        Number exp_neg = std::exp(-x);
        return (exp_pos - exp_neg) / (exp_pos + exp_neg);
    }
};

#endif
