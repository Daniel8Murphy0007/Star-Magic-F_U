/**
 * ================================================================================================
 * Header: MagnetarSGR0501_4516.h
 *
 * Description: C++ Module for SGR 0501+4516 Magnetar Class
 *              This is the first module in a series of 500+ code files for the Universal Quantum
 *              Field Framework (UQFF) simulations, focusing on magnetar evolution and gravity
 *              equations derived from Hubble datasets, high-energy lab simulations, and UQFF
 *              refinements (dated May 08, 2025, updated for full term inclusion on October 08, 2025).
 *
 * Purpose: Encapsulates the Master Universal Gravity Equation (MUGE) for SGR 0501+4516 magnetar
 *          evolution. Now includes ALL terms (no omissions): base gravity, cosmic expansion,
 *          magnetic decay, UQFF Ug components with f_TRZ, Lambda, scaled EM, GW, quantum uncertainty,
 *          fluid dynamics (effective acceleration; note: may overlap with base self-gravity),
 *          oscillatory waves (real part; assumes consistent units), and DM/density perturbations
 *          (converted to acceleration via /M). Supports dynamic variable updates for all parameters.
 *
 * Integration: Designed for inclusion in base program 'ziqn233h.cpp' (not present here).
 *              Instantiate class in main: MagnetarSGR0501_4516 mag;
 *              Compute: double g = mag.compute_g_Magnetar(t);
 *
 * Key Features:
 *   - Default values from UQFF document, with approximations for previously "negligible" terms.
 *   - Setter methods for updates: setVar(double new_val) or addToVar(double delta)/subtractFromVar(double delta).
 *   - Computes g_Magnetar(r, t) with every term explicitly included.
 *
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: October 08, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

#ifndef MAGNETAR_SGR0501_4516_H
#define MAGNETAR_SGR0501_4516_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <map>
#include <string>
#include <functional>
#include <random>
#include <algorithm>
#include <sstream>
#include <vector>

class MagnetarSGR0501_4516 {
private:
    // Core parameters (mutable for updates)
    double G;               // Gravitational constant
    double M;               // Magnetar mass
    double r;               // Radius
    double H0;              // Hubble constant (s^-1)
    double B0;              // Initial magnetic field
    double tau_B;           // B decay timescale (s)
    double B_crit;          // Critical B field
    double Lambda;          // Cosmological constant
    double c_light;         // Speed of light
    double q_charge;        // Charge (proton)
    double v_surf;          // Surface velocity
    double f_TRZ;           // Time-reversal factor
    double rho_vac_UA;      // UA vacuum density
    double rho_vac_SCm;     // SCm vacuum density
    double P_init;          // Initial rotation period (s)
    double tau_Omega;       // Omega decay timescale (s)
    double scale_EM;        // EM scaling factor
    double proton_mass;     // Proton mass for EM acceleration

    // Additional parameters for full inclusion of terms
    double hbar;            // Reduced Planck's constant
    double t_Hubble;        // Hubble time (s)
    double delta_x;         // Position uncertainty (m)
    double delta_p;         // Momentum uncertainty (kg m/s)
    double integral_psi;    // Wavefunction integral approximation
    double rho_fluid;       // Fluid density (kg/m^3)
    double A_osc;           // Oscillatory amplitude (assumed m/s^2)
    double k_osc;           // Wave number (1/m)
    double omega_osc;       // Angular frequency (rad/s)
    double x_pos;           // Position for oscillation (m)
    double t_Hubble_gyr;    // Hubble time in Gyr (for oscillatory prefactor)
    double M_DM_factor;     // Dark matter mass fraction
    double delta_rho_over_rho; // Density perturbation fraction

    // Computed caches (updated on demand)
    double ug1_base;        // Cached Ug1 = G*M/r^2

public:
    // Constructor with default UQFF values
    MagnetarSGR0501_4516() {
        initializeDefaults();
    }

    // Destructor (empty)
    ~MagnetarSGR0501_4516() {}

    // Initialization method (called in constructor)
    void initializeDefaults() {
        G = 6.6743e-11;
        M = 1.4 * 1.989e30;
        r = 20e3;
        H0 = 2.184e-18;
        B0 = 1e10;
        tau_B = 4000 * 3.156e7;
        B_crit = 1e11;
        Lambda = 1.1e-52;
        c_light = 3e8;
        q_charge = 1.602e-19;
        v_surf = 1e6;
        f_TRZ = 0.1;
        rho_vac_UA = 7.09e-36;
        rho_vac_SCm = 7.09e-37;
        P_init = 5.0;
        tau_Omega = 10000 * 3.156e7;
        scale_EM = 1e-12;
        proton_mass = 1.673e-27;

        // Full terms defaults
        hbar = 1.0546e-34;
        t_Hubble = 13.8e9 * 3.156e7;
        delta_x = 1e-10;  // Arbitrary for uncertainty principle
        delta_p = hbar / delta_x;
        integral_psi = 1.0;
        rho_fluid = 1e17;
        A_osc = 1e10;     // Arbitrary amplitude to scale ~1e10 m/s^2
        k_osc = 1.0 / r;
        omega_osc = 2 * M_PI / P_init;
        x_pos = r;
        t_Hubble_gyr = 13.8;
        M_DM_factor = 0.1;
        delta_rho_over_rho = 1e-5;

        updateCache();
    }

    // Cache update for efficiency (call after parameter changes)
    void updateCache() {
        ug1_base = (G * M) / (r * r);
        // Update dependent params if needed, e.g., delta_p = hbar / delta_x; but kept independent
    }

    // Universal setter for any variable (by name, for flexibility)
    bool setVariable(const std::string& varName, double newValue) {
        if (varName == "G") { G = newValue; }
        else if (varName == "M") { M = newValue; }
        else if (varName == "r") { r = newValue; }
        else if (varName == "H0") { H0 = newValue; }
        else if (varName == "B0") { B0 = newValue; }
        else if (varName == "tau_B") { tau_B = newValue; }
        else if (varName == "B_crit") { B_crit = newValue; }
        else if (varName == "Lambda") { Lambda = newValue; }
        else if (varName == "c_light") { c_light = newValue; }
        else if (varName == "q_charge") { q_charge = newValue; }
        else if (varName == "v_surf") { v_surf = newValue; }
        else if (varName == "f_TRZ") { f_TRZ = newValue; }
        else if (varName == "rho_vac_UA") { rho_vac_UA = newValue; }
        else if (varName == "rho_vac_SCm") { rho_vac_SCm = newValue; }
        else if (varName == "P_init") { P_init = newValue; }
        else if (varName == "tau_Omega") { tau_Omega = newValue; }
        else if (varName == "scale_EM") { scale_EM = newValue; }
        else if (varName == "proton_mass") { proton_mass = newValue; }
        // Full terms
        else if (varName == "hbar") { hbar = newValue; }
        else if (varName == "t_Hubble") { t_Hubble = newValue; }
        else if (varName == "delta_x") { delta_x = newValue; }
        else if (varName == "delta_p") { delta_p = newValue; }
        else if (varName == "integral_psi") { integral_psi = newValue; }
        else if (varName == "rho_fluid") { rho_fluid = newValue; }
        else if (varName == "A_osc") { A_osc = newValue; }
        else if (varName == "k_osc") { k_osc = newValue; }
        else if (varName == "omega_osc") { omega_osc = newValue; }
        else if (varName == "x_pos") { x_pos = newValue; }
        else if (varName == "t_Hubble_gyr") { t_Hubble_gyr = newValue; }
        else if (varName == "M_DM_factor") { M_DM_factor = newValue; }
        else if (varName == "delta_rho_over_rho") { delta_rho_over_rho = newValue; }
        else {
            std::cerr << "Error: Unknown variable '" << varName << "'." << std::endl;
            return false;
        }
        updateCache();
        return true;
    }

    // Addition method for variables
    bool addToVariable(const std::string& varName, double delta) {
        if (!setVariable(varName, getVariable(varName) + delta)) {
            return false;
        }
        updateCache();
        return true;
    }

    // Subtraction method for variables
    bool subtractFromVariable(const std::string& varName, double delta) {
        return addToVariable(varName, -delta);
    }

    // Getter for any variable (helper for add/subtract)
    double getVariable(const std::string& varName) const {
        if (varName == "G") return G;
        else if (varName == "M") return M;
        else if (varName == "r") return r;
        else if (varName == "H0") return H0;
        else if (varName == "B0") return B0;
        else if (varName == "tau_B") return tau_B;
        else if (varName == "B_crit") return B_crit;
        else if (varName == "Lambda") return Lambda;
        else if (varName == "c_light") return c_light;
        else if (varName == "q_charge") return q_charge;
        else if (varName == "v_surf") return v_surf;
        else if (varName == "f_TRZ") return f_TRZ;
        else if (varName == "rho_vac_UA") return rho_vac_UA;
        else if (varName == "rho_vac_SCm") return rho_vac_SCm;
        else if (varName == "P_init") return P_init;
        else if (varName == "tau_Omega") return tau_Omega;
        else if (varName == "scale_EM") return scale_EM;
        else if (varName == "proton_mass") return proton_mass;
        // Full terms
        else if (varName == "hbar") return hbar;
        else if (varName == "t_Hubble") return t_Hubble;
        else if (varName == "delta_x") return delta_x;
        else if (varName == "delta_p") return delta_p;
        else if (varName == "integral_psi") return integral_psi;
        else if (varName == "rho_fluid") return rho_fluid;
        else if (varName == "A_osc") return A_osc;
        else if (varName == "k_osc") return k_osc;
        else if (varName == "omega_osc") return omega_osc;
        else if (varName == "x_pos") return x_pos;
        else if (varName == "t_Hubble_gyr") return t_Hubble_gyr;
        else if (varName == "M_DM_factor") return M_DM_factor;
        else if (varName == "delta_rho_over_rho") return delta_rho_over_rho;
        else {
            std::cerr << "Error: Unknown variable '" << varName << "'." << std::endl;
            return 0.0;
        }
    }

    // B(t) computation
    double B_t(double t) const {
        return B0 * exp(-t / tau_B);
    }

    // Omega(t) computation
    double Omega_t(double t) const {
        return (2 * M_PI / P_init) * exp(-t / tau_Omega);
    }

    // dOmega/dt computation
    double dOmega_dt(double t) const {
        double omega0 = 2 * M_PI / P_init;
        return omega0 * (-1.0 / tau_Omega) * exp(-t / tau_Omega);
    }

    // Ug terms computation
    double compute_Ug(double Bt) const {
        double Ug1 = ug1_base;
        double Ug2 = 0.0;
        double Ug3 = 0.0;
        double Ug4 = Ug1 * (1 - Bt / B_crit);
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1 + f_TRZ);
    }

    // Volume computation for fluid
    double compute_V() const {
        return (4.0 / 3.0) * M_PI * r * r * r;
    }

    // Main MUGE computation (now includes ALL terms)
    double compute_g_Magnetar(double t) const {
        if (t < 0) {
            std::cerr << "Error: Time t must be non-negative." << std::endl;
            return 0.0;
        }

        double Bt = B_t(t);
        double dOdt = dOmega_dt(t);

        // Term 1: Base + H0 + B corrections
        double corr_H = 1 + H0 * t;
        double corr_B = 1 - Bt / B_crit;
        double term1 = ug1_base * corr_H * corr_B;

        // Term 2: UQFF Ug with f_TRZ
        double term2 = compute_Ug(Bt);

        // Term 3: Lambda
        double term3 = (Lambda * c_light * c_light) / 3.0;

        // Term 4: Scaled EM
        double cross_vB = v_surf * Bt;  // Magnitude, assuming perpendicular
        double em_base = (q_charge * cross_vB) / proton_mass;
        double corr_UA = 1 + (rho_vac_UA / rho_vac_SCm);
        double term4 = (em_base * corr_UA) * scale_EM;

        // Term 5: GW
        double gw_prefactor = (G * M * M) / (pow(c_light, 4) * r);
        double term5 = gw_prefactor * (dOdt * dOdt);

        // Quantum uncertainty term
        double sqrt_unc = sqrt(delta_x * delta_p);
        double term_q = (hbar / sqrt_unc) * integral_psi * (2 * M_PI / t_Hubble);

        // Fluid term (effective acceleration; note: may overlap with base self-gravity)
        double V = compute_V();
        double term_fluid = (rho_fluid * V * ug1_base) / M;

        // Oscillatory terms (real parts)
        double term_osc1 = 2 * A_osc * cos(k_osc * x_pos) * cos(omega_osc * t);
        double term_osc2 = (2 * M_PI / t_Hubble) * A_osc * cos(k_osc * x_pos - omega_osc * t);  // Adjusted for units consistency
        double term_osc = term_osc1 + term_osc2;

        // DM and density perturbation term (converted to acceleration)
        double M_dm = M * M_DM_factor;
        double pert1 = delta_rho_over_rho;
        double pert2 = 3 * G * M / (r * r * r);
        double term_dm_force_like = (M + M_dm) * (pert1 + pert2);
        double term_DM = term_dm_force_like / M;

        // Total g_Magnetar (all terms summed)
        return term1 + term2 + term3 + term4 + term5 + term_q + term_fluid + term_osc + term_DM;
    }

    // Debug/Output method (for transparency in base program)
    void printParameters(std::ostream& os = std::cout) const {
        os << std::fixed << std::setprecision(3);
        os << "SGR 0501+4516 Parameters:" << std::endl;
        os << "G: " << G << ", M: " << M << ", r: " << r << std::endl;
        os << "H0: " << H0 << ", B0: " << B0 << ", tau_B: " << tau_B << std::endl;
        os << "f_TRZ: " << f_TRZ << ", rho_fluid: " << rho_fluid << ", M_DM_factor: " << M_DM_factor << std::endl;
        os << "A_osc: " << A_osc << ", delta_rho_over_rho: " << delta_rho_over_rho << std::endl;
        os << "ug1_base: " << ug1_base << std::endl;
    }

    // Example computation at t=5000 years (for testing; now includes all terms)
    double exampleAt5000Years() const {
        double t_example = 5000 * 3.156e7;
        return compute_g_Magnetar(t_example);
    }

    // ========== ENHANCED DYNAMIC CAPABILITIES (25 methods) ==========

    // --- Variable Management (5 methods) ---
    bool createVariable(const std::string& name, double value);
    bool removeVariable(const std::string& name);
    bool cloneVariable(const std::string& src, const std::string& dest);
    std::vector<std::string> listVariables() const;
    std::string getSystemName() const;

    // --- Batch Operations (2 methods) ---
    bool transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func);
    bool scaleVariableGroup(const std::vector<std::string>& names, double factor);

    // --- Self-Expansion (4 methods) ---
    void expandParameterSpace(double factor);
    void expandMagneticScale(double B_factor, double tau_factor);
    void expandDecayScale(double omega_factor, double tau_omega_factor);
    void expandFluidDMScale(double rho_factor, double DM_factor);

    // --- Self-Refinement (3 methods) ---
    void autoRefineParameters(const std::vector<std::pair<double, double>>& observations);
    void calibrateToObservations(const std::vector<double>& times, const std::vector<double>& g_obs);
    double optimizeForMetric(std::function<double(double)> metric, double t_start, double t_end, int steps);

    // --- Parameter Exploration (1 method) ---
    std::vector<std::map<std::string, double>> generateVariations(int count, double variation_pct);

    // --- Adaptive Evolution (2 methods) ---
    void mutateParameters(double mutation_rate, std::mt19937& rng);
    void evolveSystem(int generations, std::function<double(const MagnetarSGR0501_4516&)> fitness);

    // --- State Management (4 methods) ---
    bool saveState(const std::string& stateName);
    bool restoreState(const std::string& stateName);
    std::vector<std::string> listSavedStates() const;
    std::string exportState() const;

    // --- System Analysis (4 methods) ---
    std::map<std::string, double> sensitivityAnalysis(double t, double delta_pct);
    std::string generateReport(double t) const;
    bool validateConsistency() const;
    bool autoCorrectAnomalies();
};

#endif // MAGNETAR_SGR0501_4516_Halue

// ========== IMPLEMENTATION OF ENHANCED METHODS (Outside class) ==========

// Anonymous namespace for state storage (external to class since class uses member variables)
namespace {
    std::map<std::string, std::map<std::string, double>> magnetar_sgr0501_saved_states;
}

// --- Variable Management (5 methods) ---
bool MagnetarSGR0501_4516::createVariable(const std::string& name, double value) {
    // For magnetar class with member variables, we sync to member if it exists
    return setVariable(name, value);
}

bool MagnetarSGR0501_4516::removeVariable(const std::string& name) {
    // Cannot remove built-in member variables, return false
    std::cerr << "Warning: Cannot remove built-in variable '" << name << "' in magnetar class." << std::endl;
    return false;
}

bool MagnetarSGR0501_4516::cloneVariable(const std::string& src, const std::string& dest) {
    double val = getVariable(src);
    return setVariable(dest, val);
}

std::vector<std::string> MagnetarSGR0501_4516::listVariables() const {
    return {"G", "M", "r", "H0", "B0", "tau_B", "B_crit", "Lambda", "c_light", "q_charge",
            "v_surf", "f_TRZ", "rho_vac_UA", "rho_vac_SCm", "P_init", "tau_Omega", "scale_EM", "proton_mass",
            "hbar", "t_Hubble", "delta_x", "delta_p", "integral_psi", "rho_fluid", "A_osc", "k_osc",
            "omega_osc", "x_pos", "t_Hubble_gyr", "M_DM_factor", "delta_rho_over_rho"};
}

std::string MagnetarSGR0501_4516::getSystemName() const {
    return "MagnetarSGR0501_4516";
}

// --- Batch Operations (2 methods) ---
bool MagnetarSGR0501_4516::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        double val = getVariable(name);
        if (!setVariable(name, func(val))) return false;
    }
    return true;
}

bool MagnetarSGR0501_4516::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    return transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// --- Self-Expansion (4 methods) ---
void MagnetarSGR0501_4516::expandParameterSpace(double factor) {
    std::vector<std::string> expandable = {"M", "r", "B0", "tau_B", "rho_fluid", "A_osc", "M_DM_factor"};
    scaleVariableGroup(expandable, factor);
}

void MagnetarSGR0501_4516::expandMagneticScale(double B_factor, double tau_factor) {
    setVariable("B0", getVariable("B0") * B_factor);
    setVariable("B_crit", getVariable("B_crit") * B_factor);
    setVariable("tau_B", getVariable("tau_B") * tau_factor);
}

void MagnetarSGR0501_4516::expandDecayScale(double omega_factor, double tau_omega_factor) {
    setVariable("P_init", getVariable("P_init") / omega_factor);  // Period inversely related
    setVariable("tau_Omega", getVariable("tau_Omega") * tau_omega_factor);
}

void MagnetarSGR0501_4516::expandFluidDMScale(double rho_factor, double DM_factor) {
    setVariable("rho_fluid", getVariable("rho_fluid") * rho_factor);
    setVariable("M_DM_factor", getVariable("M_DM_factor") * DM_factor);
    setVariable("delta_rho_over_rho", getVariable("delta_rho_over_rho") * DM_factor);
}

// --- Self-Refinement (3 methods) ---
void MagnetarSGR0501_4516::autoRefineParameters(const std::vector<std::pair<double, double>>& observations) {
    if (observations.empty()) return;
    
    double sum_error = 0.0;
    for (const auto& obs : observations) {
        double t = obs.first;
        double g_obs = obs.second;
        double g_calc = compute_g_Magnetar(t);
        sum_error += std::abs(g_calc - g_obs);
    }
    double avg_error = sum_error / observations.size();
    
    // Simple refinement: adjust B0 and tau_B based on error
    if (avg_error > 1e-6) {
        double adj_factor = 1.0 - std::min(0.1, avg_error / 1e6);
        setVariable("B0", getVariable("B0") * adj_factor);
        setVariable("tau_B", getVariable("tau_B") * (2.0 - adj_factor));
    }
}

void MagnetarSGR0501_4516::calibrateToObservations(const std::vector<double>& times, const std::vector<double>& g_obs) {
    if (times.size() != g_obs.size() || times.empty()) return;
    
    std::vector<std::pair<double, double>> obs;
    for (size_t i = 0; i < times.size(); ++i) {
        obs.push_back({times[i], g_obs[i]});
    }
    
    // Iterative refinement (5 passes)
    for (int iter = 0; iter < 5; ++iter) {
        autoRefineParameters(obs);
    }
}

double MagnetarSGR0501_4516::optimizeForMetric(std::function<double(double)> metric, double t_start, double t_end, int steps) {
    double best_score = -1e100;
    double dt = (t_end - t_start) / steps;
    
    for (int i = 0; i <= steps; ++i) {
        double t = t_start + i * dt;
        double g = compute_g_Magnetar(t);
        double score = metric(g);
        if (score > best_score) best_score = score;
    }
    return best_score;
}

// --- Parameter Exploration (1 method) ---
std::vector<std::map<std::string, double>> MagnetarSGR0501_4516::generateVariations(int count, double variation_pct) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-variation_pct/100.0, variation_pct/100.0);
    
    auto vars = listVariables();
    for (int i = 0; i < count; ++i) {
        std::map<std::string, double> variant;
        for (const auto& v : vars) {
            double val = getVariable(v);
            double variation = val * (1.0 + dis(gen));
            variant[v] = variation;
        }
        variations.push_back(variant);
    }
    return variations;
}

// --- Adaptive Evolution (2 methods) ---
void MagnetarSGR0501_4516::mutateParameters(double mutation_rate, std::mt19937& rng) {
    std::uniform_real_distribution<> dis(-mutation_rate, mutation_rate);
    auto vars = listVariables();
    
    for (const auto& v : vars) {
        if (v == "c_light" || v == "G" || v == "hbar") continue; // Skip constants
        double val = getVariable(v);
        double delta = val * dis(rng);
        setVariable(v, val + delta);
    }
}

void MagnetarSGR0501_4516::evolveSystem(int generations, std::function<double(const MagnetarSGR0501_4516&)> fitness) {
    std::random_device rd;
    std::mt19937 rng(rd());
    
    double best_fitness = fitness(*this);
    saveState("evolution_best");
    
    for (int gen = 0; gen < generations; ++gen) {
        saveState("evolution_temp");
        mutateParameters(0.05, rng);
        
        double new_fitness = fitness(*this);
        if (new_fitness > best_fitness) {
            best_fitness = new_fitness;
            saveState("evolution_best");
        } else {
            restoreState("evolution_temp");
        }
    }
    
    restoreState("evolution_best");
}

// --- State Management (4 methods) ---
bool MagnetarSGR0501_4516::saveState(const std::string& stateName) {
    std::map<std::string, double> state;
    auto vars = listVariables();
    for (const auto& v : vars) {
        state[v] = getVariable(v);
    }
    magnetar_sgr0501_saved_states[stateName] = state;
    return true;
}

bool MagnetarSGR0501_4516::restoreState(const std::string& stateName) {
    auto it = magnetar_sgr0501_saved_states.find(stateName);
    if (it == magnetar_sgr0501_saved_states.end()) return false;
    
    for (const auto& pair : it->second) {
        setVariable(pair.first, pair.second);
    }
    return true;
}

std::vector<std::string> MagnetarSGR0501_4516::listSavedStates() const {
    std::vector<std::string> names;
    for (const auto& pair : magnetar_sgr0501_saved_states) {
        names.push_back(pair.first);
    }
    return names;
}

std::string MagnetarSGR0501_4516::exportState() const {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(6);
    oss << "MagnetarSGR0501_4516 State Export:\n";
    auto vars = listVariables();
    for (const auto& v : vars) {
        oss << v << " = " << getVariable(v) << "\n";
    }
    return oss.str();
}

// --- System Analysis (4 methods) ---
std::map<std::string, double> MagnetarSGR0501_4516::sensitivityAnalysis(double t, double delta_pct) {
    std::map<std::string, double> sensitivities;
    double g_base = compute_g_Magnetar(t);
    
    auto vars = listVariables();
    for (const auto& v : vars) {
        if (v == "c_light" || v == "G" || v == "hbar") continue; // Skip constants
        
        double original = getVariable(v);
        double delta = original * delta_pct / 100.0;
        
        setVariable(v, original + delta);
        double g_plus = compute_g_Magnetar(t);
        setVariable(v, original);
        
        double sensitivity = (g_base != 0.0) ? std::abs((g_plus - g_base) / g_base) : 0.0;
        sensitivities[v] = sensitivity;
    }
    
    return sensitivities;
}

std::string MagnetarSGR0501_4516::generateReport(double t) const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6);
    oss << "============================================\n";
    oss << "SGR 0501+4516 MAGNETAR SYSTEM REPORT\n";
    oss << "============================================\n";
    oss << "Time: t = " << t << " s (" << (t/3.156e7) << " years)\n\n";
    
    oss << "Physical Parameters:\n";
    oss << "  Mass M = " << M << " kg (" << (M/1.989e30) << " M_sun)\n";
    oss << "  Radius r = " << r << " m (" << (r/1e3) << " km)\n";
    oss << "  B-field B0 = " << B0 << " T (B_crit = " << B_crit << " T)\n";
    oss << "  B(t) = " << B_t(t) << " T\n";
    oss << "  Period P_init = " << P_init << " s\n";
    oss << "  Omega(t) = " << Omega_t(t) << " rad/s\n";
    oss << "  Fluid density rho = " << rho_fluid << " kg/m^3\n";
    oss << "  DM factor = " << M_DM_factor << "\n\n";
    
    oss << "Computed Acceleration:\n";
    oss << "  g_Magnetar(t) = " << compute_g_Magnetar(t) << " m/s^2\n\n";
    
    oss << "UQFF Terms:\n";
    double Bt = B_t(t);
    double corr_H = 1 + H0 * t;
    double corr_B = 1 - Bt / B_crit;
    oss << "  Base (with H0, B): " << (ug1_base * corr_H * corr_B) << " m/s^2\n";
    oss << "  Ug total: " << compute_Ug(Bt) << " m/s^2\n";
    oss << "  Lambda: " << ((Lambda * c_light * c_light) / 3.0) << " m/s^2\n";
    
    double cross_vB = v_surf * Bt;
    double em_base = (q_charge * cross_vB) / proton_mass;
    double corr_UA = 1 + (rho_vac_UA / rho_vac_SCm);
    oss << "  EM (scaled): " << (em_base * corr_UA * scale_EM) << " m/s^2\n";
    
    double dOdt = dOmega_dt(t);
    double gw_prefactor = (G * M * M) / (pow(c_light, 4) * r);
    oss << "  GW: " << (gw_prefactor * dOdt * dOdt) << " m/s^2\n";
    
    double sqrt_unc = sqrt(delta_x * delta_p);
    oss << "  Quantum: " << ((hbar / sqrt_unc) * integral_psi * (2 * M_PI / t_Hubble)) << " m/s^2\n";
    
    double V = compute_V();
    oss << "  Fluid: " << ((rho_fluid * V * ug1_base) / M) << " m/s^2\n";
    
    oss << "  Oscillatory: (combined real parts)\n";
    
    double M_dm = M * M_DM_factor;
    double pert1 = delta_rho_over_rho;
    double pert2 = 3 * G * M / (r * r * r);
    double term_dm_force_like = (M + M_dm) * (pert1 + pert2);
    oss << "  DM: " << (term_dm_force_like / M) << " m/s^2\n";
    
    oss << "============================================\n";
    return oss.str();
}

bool MagnetarSGR0501_4516::validateConsistency() const {
    bool valid = true;
    
    if (M <= 0 || r <= 0) { std::cerr << "Error: M and r must be positive.\n"; valid = false; }
    if (B0 < 0 || B_crit <= 0) { std::cerr << "Error: B0, B_crit must be non-negative/positive.\n"; valid = false; }
    if (tau_B <= 0 || tau_Omega <= 0) { std::cerr << "Error: Decay timescales must be positive.\n"; valid = false; }
    if (rho_fluid < 0) { std::cerr << "Error: Fluid density must be non-negative.\n"; valid = false; }
    if (M_DM_factor < 0 || M_DM_factor > 1.0) { std::cerr << "Warning: DM factor outside [0,1].\n"; }
    
    return valid;
}

bool MagnetarSGR0501_4516::autoCorrectAnomalies() {
    bool corrected = false;
    
    if (M <= 0) { M = 1.4 * 1.989e30; corrected = true; }
    if (r <= 0) { r = 20e3; corrected = true; }
    if (B0 < 0) { B0 = 1e10; corrected = true; }
    if (B_crit <= 0) { B_crit = 1e11; corrected = true; }
    if (tau_B <= 0) { tau_B = 4000 * 3.156e7; corrected = true; }
    if (tau_Omega <= 0) { tau_Omega = 10000 * 3.156e7; corrected = true; }
    if (rho_fluid < 0) { rho_fluid = 1e17; corrected = true; }
    if (M_DM_factor < 0) { M_DM_factor = 0.1; corrected = true; }
    if (M_DM_factor > 1.0) { M_DM_factor = 1.0; corrected = true; }
    
    if (corrected) updateCache();
    return corrected;
}

// ========== ENHANCED EXAMPLE FUNCTION ==========
void enhancedMagnetarSGR0501Example() {
    std::cout << "\n========== ENHANCED SGR 0501+4516 MAGNETAR UQFF EXAMPLE ==========\n\n";
    
    MagnetarSGR0501_4516 mag;
    
    // Step 1: Initial state
    std::cout << "Step 1: Initial Configuration\n";
    mag.printParameters();
    double t0 = 0.0;
    std::cout << "g_Magnetar(t=0) = " << mag.compute_g_Magnetar(t0) << " m/s^2\n\n";
    
    // Step 2: Time evolution
    std::cout << "Step 2: Time Evolution (0, 1000, 3000, 5000 years)\n";
    for (double t_yr : {0.0, 1000.0, 3000.0, 5000.0}) {
        double t = t_yr * 3.156e7;
        std::cout << "  t = " << t_yr << " yr: g = " << mag.compute_g_Magnetar(t) 
                  << " m/s^2, B(t) = " << mag.B_t(t) << " T, Omega(t) = " << mag.Omega_t(t) << " rad/s\n";
    }
    std::cout << "\n";
    
    // Step 3: Magnetic field scaling
    std::cout << "Step 3: Magnetic Field Scaling (B0 x1.5, tau_B x0.8)\n";
    mag.expandMagneticScale(1.5, 0.8);
    double t_test = 2000 * 3.156e7;
    std::cout << "After expansion: B0 = " << mag.getVariable("B0") << " T, tau_B = " << mag.getVariable("tau_B") << " s\n";
    std::cout << "g_Magnetar(t=2000 yr) = " << mag.compute_g_Magnetar(t_test) << " m/s^2\n\n";
    
    // Step 4: Rotation decay scaling
    std::cout << "Step 4: Rotation Decay Scaling (Omega x1.2, tau_Omega x1.5)\n";
    mag.expandDecayScale(1.2, 1.5);
    std::cout << "After expansion: P_init = " << mag.getVariable("P_init") << " s, tau_Omega = " << mag.getVariable("tau_Omega") << " s\n";
    std::cout << "g_Magnetar(t=2000 yr) = " << mag.compute_g_Magnetar(t_test) << " m/s^2\n\n";
    
    // Step 5: Fluid & DM scaling
    std::cout << "Step 5: Fluid & DM Scaling (rho x2.0, DM factor x1.3)\n";
    mag.expandFluidDMScale(2.0, 1.3);
    std::cout << "After expansion: rho_fluid = " << mag.getVariable("rho_fluid") << " kg/m^3, M_DM_factor = " << mag.getVariable("M_DM_factor") << "\n";
    std::cout << "g_Magnetar(t=2000 yr) = " << mag.compute_g_Magnetar(t_test) << " m/s^2\n\n";
    
    // Step 6: State save/restore
    std::cout << "Step 6: State Management\n";
    mag.saveState("expanded_state");
    mag.setVariable("B0", 5e9);
    std::cout << "Modified B0 to 5e9 T: g = " << mag.compute_g_Magnetar(t_test) << " m/s^2\n";
    mag.restoreState("expanded_state");
    std::cout << "Restored state: B0 = " << mag.getVariable("B0") << " T, g = " << mag.compute_g_Magnetar(t_test) << " m/s^2\n";
    std::cout << "Saved states: ";
    for (const auto& s : mag.listSavedStates()) std::cout << s << " ";
    std::cout << "\n\n";
    
    // Step 7: Sensitivity analysis
    std::cout << "Step 7: Sensitivity Analysis at t=2000 yr (top 5 parameters)\n";
    auto sens = mag.sensitivityAnalysis(t_test, 1.0);
    std::vector<std::pair<std::string, double>> sens_vec(sens.begin(), sens.end());
    std::sort(sens_vec.begin(), sens_vec.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
    for (int i = 0; i < std::min(5, (int)sens_vec.size()); ++i) {
        std::cout << "  " << sens_vec[i].first << ": " << sens_vec[i].second << "\n";
    }
    std::cout << "\n";
    
    // Step 8: Generate variations
    std::cout << "Step 8: Generate Parameter Variations (3 variants, 10% variation)\n";
    auto variations = mag.generateVariations(3, 10.0);
    for (size_t i = 0; i < variations.size(); ++i) {
        std::cout << "  Variant " << (i+1) << ": B0 = " << variations[i]["B0"] 
                  << " T, M = " << variations[i]["M"] << " kg\n";
    }
    std::cout << "\n";
    
    // Step 9: Batch transformation
    std::cout << "Step 9: Batch Transform (scale mass parameters by 1.1)\n";
    mag.transformVariableGroup({"M", "rho_fluid"}, [](double v) { return v * 1.1; });
    std::cout << "After transform: M = " << mag.getVariable("M") << " kg, rho_fluid = " << mag.getVariable("rho_fluid") << " kg/m^3\n";
    std::cout << "g_Magnetar(t=2000 yr) = " << mag.compute_g_Magnetar(t_test) << " m/s^2\n\n";
    
    // Step 10: Consistency validation
    std::cout << "Step 10: Consistency Validation\n";
    bool valid = mag.validateConsistency();
    std::cout << "System is " << (valid ? "VALID" : "INVALID") << "\n\n";
    
    // Step 11: Metric optimization
    std::cout << "Step 11: Optimize for Maximum g (t=0 to 10000 yr, 100 steps)\n";
    double max_g = mag.optimizeForMetric([](double g) { return g; }, 0.0, 10000 * 3.156e7, 100);
    std::cout << "Maximum g found: " << max_g << " m/s^2\n\n";
    
    // Step 12: Full system report
    std::cout << "Step 12: Full System Report at t=3000 yr\n";
    double t_report = 3000 * 3.156e7;
    std::cout << mag.generateReport(t_report) << "\n";
    
    // Step 13: B-field sweep
    std::cout << "Step 13: B-Field Sweep (B0 = 0.5e10, 1.0e10, 1.5e10, 2.0e10 T)\n";
    mag.saveState("before_sweep");
    for (double B0_val : {0.5e10, 1.0e10, 1.5e10, 2.0e10}) {
        mag.setVariable("B0", B0_val);
        std::cout << "  B0 = " << B0_val << " T: g(t=2000yr) = " << mag.compute_g_Magnetar(t_test) << " m/s^2\n";
    }
    mag.restoreState("before_sweep");
    std::cout << "\n";
    
    // Step 14: Decay timescale sweep
    std::cout << "Step 14: Magnetic Decay Timescale Sweep (tau_B = 2000, 4000, 6000, 8000 yr)\n";
    for (double tau_yr : {2000.0, 4000.0, 6000.0, 8000.0}) {
        mag.setVariable("tau_B", tau_yr * 3.156e7);
        std::cout << "  tau_B = " << tau_yr << " yr: g(t=2000yr) = " << mag.compute_g_Magnetar(t_test) << " m/s^2\n";
    }
    mag.restoreState("before_sweep");
    std::cout << "\n";
    
    // Step 15: Rotation period sweep
    std::cout << "Step 15: Initial Rotation Period Sweep (P_init = 3, 5, 7, 10 s)\n";
    for (double P_val : {3.0, 5.0, 7.0, 10.0}) {
        mag.setVariable("P_init", P_val);
        std::cout << "  P_init = " << P_val << " s: g(t=2000yr) = " << mag.compute_g_Magnetar(t_test) 
                  << " m/s^2, Omega(t=2000yr) = " << mag.Omega_t(t_test) << " rad/s\n";
    }
    mag.restoreState("before_sweep");
    std::cout << "\n";
    
    // Step 16: DM factor sweep
    std::cout << "Step 16: Dark Matter Factor Sweep (M_DM_factor = 0.0, 0.1, 0.2, 0.3)\n";
    for (double dm : {0.0, 0.1, 0.2, 0.3}) {
        mag.setVariable("M_DM_factor", dm);
        std::cout << "  M_DM_factor = " << dm << ": g(t=2000yr) = " << mag.compute_g_Magnetar(t_test) << " m/s^2\n";
    }
    mag.restoreState("before_sweep");
    std::cout << "\n";
    
    // Step 17: Export final state
    std::cout << "Step 17: Export Final State\n";
    std::cout << mag.exportState() << "\n";
    
    std::cout << "========== ENHANCED EXAMPLE COMPLETE ==========\n\n";
}