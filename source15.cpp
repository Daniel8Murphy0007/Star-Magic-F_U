/**
 * ================================================================================================
 * Header: SMBHSgrAStar.h
 *
 * Description: C++ Module for Sagittarius A* (Sgr A*) Supermassive Black Hole Class
 *              This is the third module in a series of 500+ code files for the Universal Quantum
 *              Field Framework (UQFF) simulations, focusing on SMBH evolution and gravity
 *              equations derived from Hubble datasets, high-energy lab simulations, and UQFF
 *              refinements (dated May 09, 2025, updated for full term inclusion on October 08, 2025).
 *
 * Purpose: Encapsulates the Master Universal Gravity Equation (MUGE) for Sgr A* evolution.
 *          Includes ALL terms: base gravity with mass growth M(t), cosmic expansion (H_0), magnetic decay,
 *          UQFF Ug components with f_TRZ, Lambda, quantum uncertainty, EM (with B(t)), fluid dynamics,
 *          oscillatory waves, DM/density perturbations with precession sin(30�), and GW term.
 *          Supports dynamic variable updates for all parameters.
 *
 * Integration: Designed for inclusion in base program 'ziqn233h.cpp' (not present here).
 *              Instantiate class in main: SMBHSgrAStar sgrA;
 *              Compute: double g = sgrA.compute_g_SgrA(t);
 *
 * Key Features:
 *   - Default values from UQFF document, with approximations for all terms.
 *   - Units handled: B(t) converted to T (1 G = 10^-4 T); energy terms converted to effective acceleration.
 *   - Setter methods for updates: setVar(double new_val) or addToVar(double delta)/subtractFromVar(double delta).
 *   - Computes g_SgrA(r, t) with every term explicitly included.
 *
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: October 08, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

#ifndef SMBH_SGR_A_STAR_H
#define SMBH_SGR_A_STAR_H

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

class SMBHSgrAStar {
private:
    // Core parameters (mutable for updates)
    double G;               // Gravitational constant
    double M_initial;       // Initial SMBH mass
    double r;               // Schwarzschild radius
    double H0;              // Hubble constant (s^-1)
    double B0_G;            // Initial magnetic field (G)
    double tau_B;           // B decay timescale (s)
    double B_crit;          // Critical B field (T)
    double Lambda;          // Cosmological constant
    double c_light;         // Speed of light
    double q_charge;        // Charge (proton)
    double v_surf;          // Surface velocity (arbitrary for BH)
    double f_TRZ;           // Time-reversal factor
    double M_dot_0;         // Initial mass accretion rate factor
    double tau_acc;         // Accretion timescale (s)
    double spin_factor;     // Spin factor (0.3)
    double tau_Omega;       // Omega decay timescale (s)

    // Additional parameters for full inclusion of terms
    double hbar;            // Reduced Planck's constant
    double t_Hubble;        // Hubble time (s)
    double delta_x;         // Position uncertainty (m)
    double delta_p;         // Momentum uncertainty (kg m/s)
    double integral_psi;    // Wavefunction integral approximation
    double rho_fluid;       // Fluid density (kg/m^3)
    double A_osc;           // Oscillatory amplitude (m/s^2)
    double k_osc;           // Wave number (1/m)
    double omega_osc;       // Angular frequency (rad/s)
    double x_pos;           // Position for oscillation (m)
    double t_Hubble_gyr;    // Hubble time in Gyr
    double M_DM_factor;     // Dark matter mass fraction
    double delta_rho_over_rho; // Density perturbation fraction
    double precession_angle_deg; // Precession angle (degrees)

    // Computed caches (updated on demand)
    double ug1_base;        // Cached Ug1 for initial M (will recompute with M(t))

public:
    // Constructor with default UQFF values
    SMBHSgrAStar() {
        initializeDefaults();
    }

    // Destructor (empty)
    ~SMBHSgrAStar() {}

    // Initialization method (called in constructor)
    void initializeDefaults() {
        G = 6.6743e-11;
        M_initial = 4.3e6 * 1.989e30;
        r = 1.27e10;
        H0 = 2.184e-18;
        B0_G = 1e4;  // G
        tau_B = 1e6 * 3.156e7;
        B_crit = 1e11;  // T
        Lambda = 1.1e-52;
        c_light = 3e8;
        q_charge = 1.602e-19;
        v_surf = 1e6;  // Arbitrary
        f_TRZ = 0.1;
        M_dot_0 = 0.01;
        tau_acc = 9e9 * 3.156e7;
        spin_factor = 0.3;
        tau_Omega = 9e9 * 3.156e7;

        // Full terms defaults
        hbar = 1.0546e-34;
        t_Hubble = 13.8e9 * 3.156e7;
        t_Hubble_gyr = 13.8;
        delta_x = 1e-10;
        delta_p = hbar / delta_x;
        integral_psi = 1.0;
        rho_fluid = 1e17;  // Arbitrary for accretion disk
        A_osc = 1e6;       // Scaled down for BH
        k_osc = 1.0 / r;
        omega_osc = 2 * M_PI / (r / c_light);  // Orbital-like
        x_pos = r;
        M_DM_factor = 0.1;
        delta_rho_over_rho = 1e-5;
        precession_angle_deg = 30.0;

        updateCache();
    }

    // Cache update for efficiency (call after parameter changes)
    void updateCache() {
        ug1_base = (G * M_initial) / (r * r);
    }

    // Universal setter for any variable (by name, for flexibility)
    bool setVariable(const std::string& varName, double newValue) {
        if (varName == "G") { G = newValue; }
        else if (varName == "M_initial") { M_initial = newValue; }
        else if (varName == "r") { r = newValue; }
        else if (varName == "H0") { H0 = newValue; }
        else if (varName == "B0_G") { B0_G = newValue; }
        else if (varName == "tau_B") { tau_B = newValue; }
        else if (varName == "B_crit") { B_crit = newValue; }
        else if (varName == "Lambda") { Lambda = newValue; }
        else if (varName == "c_light") { c_light = newValue; }
        else if (varName == "q_charge") { q_charge = newValue; }
        else if (varName == "v_surf") { v_surf = newValue; }
        else if (varName == "f_TRZ") { f_TRZ = newValue; }
        else if (varName == "M_dot_0") { M_dot_0 = newValue; }
        else if (varName == "tau_acc") { tau_acc = newValue; }
        else if (varName == "spin_factor") { spin_factor = newValue; }
        else if (varName == "tau_Omega") { tau_Omega = newValue; }
        // Full terms
        else if (varName == "hbar") { hbar = newValue; }
        else if (varName == "t_Hubble") { t_Hubble = newValue; }
        else if (varName == "t_Hubble_gyr") { t_Hubble_gyr = newValue; }
        else if (varName == "delta_x") { delta_x = newValue; }
        else if (varName == "delta_p") { delta_p = newValue; }
        else if (varName == "integral_psi") { integral_psi = newValue; }
        else if (varName == "rho_fluid") { rho_fluid = newValue; }
        else if (varName == "A_osc") { A_osc = newValue; }
        else if (varName == "k_osc") { k_osc = newValue; }
        else if (varName == "omega_osc") { omega_osc = newValue; }
        else if (varName == "x_pos") { x_pos = newValue; }
        else if (varName == "M_DM_factor") { M_DM_factor = newValue; }
        else if (varName == "delta_rho_over_rho") { delta_rho_over_rho = newValue; }
        else if (varName == "precession_angle_deg") { precession_angle_deg = newValue; }
        else {
            std::cerr << "Error: Unknown variable '" << varName << "'." << std::endl;
            return false;
        }
        updateCache();
        return true;
    }

    // Addition method for variables
    bool addToVariable(const std::string& varName, double delta) {
        return setVariable(varName, getVariable(varName) + delta);
    }

    // Subtraction method for variables
    bool subtractFromVariable(const std::string& varName, double delta) {
        return addToVariable(varName, -delta);
    }

    // Getter for any variable (helper for add/subtract)
    double getVariable(const std::string& varName) const {
        if (varName == "G") return G;
        else if (varName == "M_initial") return M_initial;
        else if (varName == "r") return r;
        else if (varName == "H0") return H0;
        else if (varName == "B0_G") return B0_G;
        else if (varName == "tau_B") return tau_B;
        else if (varName == "B_crit") return B_crit;
        else if (varName == "Lambda") return Lambda;
        else if (varName == "c_light") return c_light;
        else if (varName == "q_charge") return q_charge;
        else if (varName == "v_surf") return v_surf;
        else if (varName == "f_TRZ") return f_TRZ;
        else if (varName == "M_dot_0") return M_dot_0;
        else if (varName == "tau_acc") return tau_acc;
        else if (varName == "spin_factor") return spin_factor;
        else if (varName == "tau_Omega") return tau_Omega;
        // Full terms
        else if (varName == "hbar") return hbar;
        else if (varName == "t_Hubble") return t_Hubble;
        else if (varName == "t_Hubble_gyr") return t_Hubble_gyr;
        else if (varName == "delta_x") return delta_x;
        else if (varName == "delta_p") return delta_p;
        else if (varName == "integral_psi") return integral_psi;
        else if (varName == "rho_fluid") return rho_fluid;
        else if (varName == "A_osc") return A_osc;
        else if (varName == "k_osc") return k_osc;
        else if (varName == "omega_osc") return omega_osc;
        else if (varName == "x_pos") return x_pos;
        else if (varName == "M_DM_factor") return M_DM_factor;
        else if (varName == "delta_rho_over_rho") return delta_rho_over_rho;
        else if (varName == "precession_angle_deg") return precession_angle_deg;
        else {
            std::cerr << "Error: Unknown variable '" << varName << "'." << std::endl;
            return 0.0;
        }
    }

    // M(t) computation
    double M_t(double t) const {
        double M_dot = M_dot_0 * exp(-t / tau_acc);
        return M_initial * (1 + M_dot);
    }

    // B(t) in T
    double B_t(double t) const {
        double B_G = B0_G * exp(-t / tau_B);
        return B_G * 1e-4;  // G to T
    }

    // Omega(t) computation
    double Omega_t(double t) const {
        double omega0 = spin_factor * c_light / r;
        return omega0 * exp(-t / tau_Omega);
    }

    // dOmega/dt computation
    double dOmega_dt(double t) const {
        double omega0 = spin_factor * c_light / r;
        return omega0 * (-1.0 / tau_Omega) * exp(-t / tau_Omega);
    }

    // Ug terms computation
    double compute_Ug(double Mt, double Bt) const {
        double Ug1 = (G * Mt) / (r * r);
        double Ug2 = 0.0;
        double Ug3 = 0.0;
        double corr_B = 1 - Bt / B_crit;
        double Ug4 = Ug1 * corr_B;
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1 + f_TRZ);
    }

    // Volume computation for fluid
    double compute_V() const {
        return (4.0 / 3.0) * M_PI * r * r * r;
    }

    // Main MUGE computation (includes ALL terms)
    double compute_g_SgrA(double t) const {
        if (t < 0) {
            std::cerr << "Error: Time t must be non-negative." << std::endl;
            return 0.0;
        }

        double Mt = M_t(t);
        double Bt = B_t(t);
        double dOdt = dOmega_dt(t);
        double ug1_t = (G * Mt) / (r * r);

        // Term 1: Base + H0 + B corrections
        double corr_H = 1 + H0 * t;
        double corr_B = 1 - Bt / B_crit;
        double term1 = ug1_t * corr_H * corr_B;

        // Term 2: UQFF Ug with f_TRZ
        double term2 = compute_Ug(Mt, Bt);

        // Term 3: Lambda
        double term3 = (Lambda * c_light * c_light) / 3.0;

        // Term 4: EM (v x B, no scaling or UA here)
        double cross_vB = v_surf * Bt;  // Magnitude
        double em_base = q_charge * cross_vB / 1.673e-27;  // Acceleration
        double term4 = em_base;

        // Term 5: GW
        double gw_prefactor = (G * Mt * Mt) / (pow(c_light, 4) * r);
        double term5 = gw_prefactor * (dOdt * dOdt);

        // Quantum uncertainty term
        double sqrt_unc = sqrt(delta_x * delta_p);
        double term_q = (hbar / sqrt_unc) * integral_psi * (2 * M_PI / t_Hubble);

        // Fluid term (effective acceleration)
        double V = compute_V();
        double term_fluid = (rho_fluid * V * ug1_t) / Mt;

        // Oscillatory terms (real parts)
        double term_osc1 = 2 * A_osc * cos(k_osc * x_pos) * cos(omega_osc * t);
        double arg = k_osc * x_pos - omega_osc * t;
        double term_osc2 = (2 * M_PI / t_Hubble_gyr) * A_osc * cos(arg);
        double term_osc = term_osc1 + term_osc2;

        // DM and density perturbation term with precession (converted to acceleration)
        double M_dm = Mt * M_DM_factor;
        double sin_prec = sin(precession_angle_deg * M_PI / 180.0);
        double pert1 = delta_rho_over_rho;
        double pert2 = 3 * G * Mt / (r * r * r);
        double term_dm_force_like = (Mt + M_dm) * (pert1 + pert2 * sin_prec);
        double term_DM = term_dm_force_like / Mt;

        // Total g_SgrA (all terms summed)
        return term1 + term2 + term3 + term4 + term5 + term_q + term_fluid + term_osc + term_DM;
    }

    // Debug/Output method (for transparency in base program)
    void printParameters(std::ostream& os = std::cout) const {
        os << std::fixed << std::setprecision(3);
        os << "Sgr A* Parameters:" << std::endl;
        os << "G: " << G << ", M_initial: " << M_initial << ", r: " << r << std::endl;
        os << "H0: " << H0 << ", B0_G: " << B0_G << ", tau_B: " << tau_B << std::endl;
        os << "f_TRZ: " << f_TRZ << ", M_dot_0: " << M_dot_0 << ", tau_acc: " << tau_acc << std::endl;
        os << "rho_fluid: " << rho_fluid << ", M_DM_factor: " << M_DM_factor << std::endl;
        os << "A_osc: " << A_osc << ", precession_angle_deg: " << precession_angle_deg << std::endl;
        os << "ug1_base: " << ug1_base << std::endl;
    }

    // Example computation at t=4.5 Gyr (for testing)
    double exampleAt4_5Gyr() const {
        double t_example = 4.5e9 * 3.156e7;
        return compute_g_SgrA(t_example);
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
    void expandAccretionScale(double M_dot_factor, double tau_acc_factor);
    void expandMagneticScale(double B_factor, double tau_B_factor);
    void expandDMPrecessionScale(double DM_factor, double angle_factor);

    // --- Self-Refinement (3 methods) ---
    void autoRefineParameters(const std::vector<std::pair<double, double>>& observations);
    void calibrateToObservations(const std::vector<double>& times, const std::vector<double>& g_obs);
    double optimizeForMetric(std::function<double(double)> metric, double t_start, double t_end, int steps);

    // --- Parameter Exploration (1 method) ---
    std::vector<std::map<std::string, double>> generateVariations(int count, double variation_pct);

    // --- Adaptive Evolution (2 methods) ---
    void mutateParameters(double mutation_rate, std::mt19937& rng);
    void evolveSystem(int generations, std::function<double(const SMBHSgrAStar&)> fitness);

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

#endif // SMBH_SGR_A_STAR_H

// ========== IMPLEMENTATION OF ENHANCED METHODS (Outside class) ==========

// Anonymous namespace for state storage (external to class since class uses member variables)
namespace {
    std::map<std::string, std::map<std::string, double>> smbh_sgra_saved_states;
}

// --- Variable Management (5 methods) ---
bool SMBHSgrAStar::createVariable(const std::string& name, double value) {
    // For SMBH class with member variables, we sync to member if it exists
    return setVariable(name, value);
}

bool SMBHSgrAStar::removeVariable(const std::string& name) {
    // Cannot remove built-in member variables, return false
    std::cerr << "Warning: Cannot remove built-in variable '" << name << "' in SMBH class." << std::endl;
    return false;
}

bool SMBHSgrAStar::cloneVariable(const std::string& src, const std::string& dest) {
    double val = getVariable(src);
    return setVariable(dest, val);
}

std::vector<std::string> SMBHSgrAStar::listVariables() const {
    return {"G", "M_initial", "r", "H0", "B0_G", "tau_B", "B_crit", "Lambda", "c_light", "q_charge",
            "v_surf", "f_TRZ", "M_dot_0", "tau_acc", "spin_factor", "tau_Omega",
            "hbar", "t_Hubble", "t_Hubble_gyr", "delta_x", "delta_p", "integral_psi", "rho_fluid",
            "A_osc", "k_osc", "omega_osc", "x_pos", "M_DM_factor", "delta_rho_over_rho", "precession_angle_deg"};
}

std::string SMBHSgrAStar::getSystemName() const {
    return "SMBHSgrAStar";
}

// --- Batch Operations (2 methods) ---
bool SMBHSgrAStar::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        double val = getVariable(name);
        if (!setVariable(name, func(val))) return false;
    }
    return true;
}

bool SMBHSgrAStar::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    return transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// --- Self-Expansion (4 methods) ---
void SMBHSgrAStar::expandParameterSpace(double factor) {
    std::vector<std::string> expandable = {"M_initial", "r", "B0_G", "tau_B", "rho_fluid", "A_osc", "M_DM_factor"};
    scaleVariableGroup(expandable, factor);
}

void SMBHSgrAStar::expandAccretionScale(double M_dot_factor, double tau_acc_factor) {
    setVariable("M_dot_0", getVariable("M_dot_0") * M_dot_factor);
    setVariable("tau_acc", getVariable("tau_acc") * tau_acc_factor);
}

void SMBHSgrAStar::expandMagneticScale(double B_factor, double tau_B_factor) {
    setVariable("B0_G", getVariable("B0_G") * B_factor);
    setVariable("B_crit", getVariable("B_crit") * B_factor);
    setVariable("tau_B", getVariable("tau_B") * tau_B_factor);
}

void SMBHSgrAStar::expandDMPrecessionScale(double DM_factor, double angle_factor) {
    setVariable("M_DM_factor", getVariable("M_DM_factor") * DM_factor);
    setVariable("delta_rho_over_rho", getVariable("delta_rho_over_rho") * DM_factor);
    setVariable("precession_angle_deg", getVariable("precession_angle_deg") * angle_factor);
}

// --- Self-Refinement (3 methods) ---
void SMBHSgrAStar::autoRefineParameters(const std::vector<std::pair<double, double>>& observations) {
    if (observations.empty()) return;
    
    double sum_error = 0.0;
    for (const auto& obs : observations) {
        double t = obs.first;
        double g_obs = obs.second;
        double g_calc = compute_g_SgrA(t);
        sum_error += std::abs(g_calc - g_obs);
    }
    double avg_error = sum_error / observations.size();
    
    // Simple refinement: adjust M_dot_0 and tau_acc based on error
    if (avg_error > 1e-6) {
        double adj_factor = 1.0 - std::min(0.1, avg_error / 1e6);
        setVariable("M_dot_0", getVariable("M_dot_0") * adj_factor);
        setVariable("tau_acc", getVariable("tau_acc") * (2.0 - adj_factor));
    }
}

void SMBHSgrAStar::calibrateToObservations(const std::vector<double>& times, const std::vector<double>& g_obs) {
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

double SMBHSgrAStar::optimizeForMetric(std::function<double(double)> metric, double t_start, double t_end, int steps) {
    double best_score = -1e100;
    double dt = (t_end - t_start) / steps;
    
    for (int i = 0; i <= steps; ++i) {
        double t = t_start + i * dt;
        double g = compute_g_SgrA(t);
        double score = metric(g);
        if (score > best_score) best_score = score;
    }
    return best_score;
}

// --- Parameter Exploration (1 method) ---
std::vector<std::map<std::string, double>> SMBHSgrAStar::generateVariations(int count, double variation_pct) {
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
void SMBHSgrAStar::mutateParameters(double mutation_rate, std::mt19937& rng) {
    std::uniform_real_distribution<> dis(-mutation_rate, mutation_rate);
    auto vars = listVariables();
    
    for (const auto& v : vars) {
        if (v == "c_light" || v == "G" || v == "hbar") continue; // Skip constants
        double val = getVariable(v);
        double delta = val * dis(rng);
        setVariable(v, val + delta);
    }
}

void SMBHSgrAStar::evolveSystem(int generations, std::function<double(const SMBHSgrAStar&)> fitness) {
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
bool SMBHSgrAStar::saveState(const std::string& stateName) {
    std::map<std::string, double> state;
    auto vars = listVariables();
    for (const auto& v : vars) {
        state[v] = getVariable(v);
    }
    smbh_sgra_saved_states[stateName] = state;
    return true;
}

bool SMBHSgrAStar::restoreState(const std::string& stateName) {
    auto it = smbh_sgra_saved_states.find(stateName);
    if (it == smbh_sgra_saved_states.end()) return false;
    
    for (const auto& pair : it->second) {
        setVariable(pair.first, pair.second);
    }
    return true;
}

std::vector<std::string> SMBHSgrAStar::listSavedStates() const {
    std::vector<std::string> names;
    for (const auto& pair : smbh_sgra_saved_states) {
        names.push_back(pair.first);
    }
    return names;
}

std::string SMBHSgrAStar::exportState() const {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(6);
    oss << "SMBHSgrAStar State Export:\n";
    auto vars = listVariables();
    for (const auto& v : vars) {
        oss << v << " = " << getVariable(v) << "\n";
    }
    return oss.str();
}

// --- System Analysis (4 methods) ---
std::map<std::string, double> SMBHSgrAStar::sensitivityAnalysis(double t, double delta_pct) {
    std::map<std::string, double> sensitivities;
    double g_base = compute_g_SgrA(t);
    
    auto vars = listVariables();
    for (const auto& v : vars) {
        if (v == "c_light" || v == "G" || v == "hbar") continue; // Skip constants
        
        double original = getVariable(v);
        double delta = original * delta_pct / 100.0;
        
        setVariable(v, original + delta);
        double g_plus = compute_g_SgrA(t);
        setVariable(v, original);
        
        double sensitivity = (g_base != 0.0) ? std::abs((g_plus - g_base) / g_base) : 0.0;
        sensitivities[v] = sensitivity;
    }
    
    return sensitivities;
}

std::string SMBHSgrAStar::generateReport(double t) const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6);
    oss << "============================================\n";
    oss << "SGR A* SUPERMASSIVE BLACK HOLE REPORT\n";
    oss << "============================================\n";
    oss << "Time: t = " << t << " s (" << (t/3.156e7/1e9) << " Gyr)\n\n";
    
    oss << "Physical Parameters:\n";
    oss << "  Initial Mass M_initial = " << M_initial << " kg (" << (M_initial/1.989e30) << " M_sun)\n";
    oss << "  M(t) = " << M_t(t) << " kg (" << (M_t(t)/1.989e30) << " M_sun)\n";
    oss << "  Schwarzschild radius r = " << r << " m (" << (r/1e3) << " km)\n";
    oss << "  B-field B0 = " << B0_G << " G (B_crit = " << B_crit << " T)\n";
    oss << "  B(t) = " << B_t(t) << " T\n";
    oss << "  Accretion M_dot_0 = " << M_dot_0 << ", tau_acc = " << tau_acc << " s\n";
    oss << "  Spin factor = " << spin_factor << ", Omega(t) = " << Omega_t(t) << " rad/s\n";
    oss << "  Fluid density rho = " << rho_fluid << " kg/m^3\n";
    oss << "  DM factor = " << M_DM_factor << ", precession angle = " << precession_angle_deg << " deg\n\n";
    
    oss << "Computed Acceleration:\n";
    oss << "  g_SgrA(t) = " << compute_g_SgrA(t) << " m/s^2\n\n";
    
    oss << "UQFF Terms:\n";
    double Mt = M_t(t);
    double Bt = B_t(t);
    double ug1_t = (G * Mt) / (r * r);
    double corr_H = 1 + H0 * t;
    double corr_B = 1 - Bt / B_crit;
    oss << "  Base (with H0, B, M(t)): " << (ug1_t * corr_H * corr_B) << " m/s^2\n";
    oss << "  Ug total: " << compute_Ug(Mt, Bt) << " m/s^2\n";
    oss << "  Lambda: " << ((Lambda * c_light * c_light) / 3.0) << " m/s^2\n";
    
    double cross_vB = v_surf * Bt;
    double em_base = q_charge * cross_vB / 1.673e-27;
    oss << "  EM: " << em_base << " m/s^2\n";
    
    double dOdt = dOmega_dt(t);
    double gw_prefactor = (G * Mt * Mt) / (pow(c_light, 4) * r);
    oss << "  GW: " << (gw_prefactor * dOdt * dOdt) << " m/s^2\n";
    
    double sqrt_unc = sqrt(delta_x * delta_p);
    oss << "  Quantum: " << ((hbar / sqrt_unc) * integral_psi * (2 * M_PI / t_Hubble)) << " m/s^2\n";
    
    double V = compute_V();
    oss << "  Fluid: " << ((rho_fluid * V * ug1_t) / Mt) << " m/s^2\n";
    
    oss << "  Oscillatory: (combined real parts)\n";
    
    double M_dm = Mt * M_DM_factor;
    double sin_prec = sin(precession_angle_deg * M_PI / 180.0);
    double pert1 = delta_rho_over_rho;
    double pert2 = 3 * G * Mt / (r * r * r);
    double term_dm_force_like = (Mt + M_dm) * (pert1 + pert2 * sin_prec);
    oss << "  DM (with precession): " << (term_dm_force_like / Mt) << " m/s^2\n";
    
    oss << "============================================\n";
    return oss.str();
}

bool SMBHSgrAStar::validateConsistency() const {
    bool valid = true;
    
    if (M_initial <= 0 || r <= 0) { std::cerr << "Error: M_initial and r must be positive.\n"; valid = false; }
    if (B0_G < 0 || B_crit <= 0) { std::cerr << "Error: B0_G, B_crit must be non-negative/positive.\n"; valid = false; }
    if (tau_B <= 0 || tau_acc <= 0 || tau_Omega <= 0) { std::cerr << "Error: Decay/accretion timescales must be positive.\n"; valid = false; }
    if (rho_fluid < 0) { std::cerr << "Error: Fluid density must be non-negative.\n"; valid = false; }
    if (M_DM_factor < 0 || M_DM_factor > 1.0) { std::cerr << "Warning: DM factor outside [0,1].\n"; }
    if (spin_factor < 0 || spin_factor > 1.0) { std::cerr << "Warning: Spin factor outside [0,1].\n"; }
    
    return valid;
}

bool SMBHSgrAStar::autoCorrectAnomalies() {
    bool corrected = false;
    
    if (M_initial <= 0) { M_initial = 4.3e6 * 1.989e30; corrected = true; }
    if (r <= 0) { r = 1.27e10; corrected = true; }
    if (B0_G < 0) { B0_G = 1e4; corrected = true; }
    if (B_crit <= 0) { B_crit = 1e11; corrected = true; }
    if (tau_B <= 0) { tau_B = 1e6 * 3.156e7; corrected = true; }
    if (tau_acc <= 0) { tau_acc = 9e9 * 3.156e7; corrected = true; }
    if (tau_Omega <= 0) { tau_Omega = 9e9 * 3.156e7; corrected = true; }
    if (rho_fluid < 0) { rho_fluid = 1e17; corrected = true; }
    if (M_DM_factor < 0) { M_DM_factor = 0.1; corrected = true; }
    if (M_DM_factor > 1.0) { M_DM_factor = 1.0; corrected = true; }
    if (spin_factor < 0) { spin_factor = 0.0; corrected = true; }
    if (spin_factor > 1.0) { spin_factor = 1.0; corrected = true; }
    
    if (corrected) updateCache();
    return corrected;
}

// ========== ENHANCED EXAMPLE FUNCTION ==========
void enhancedSgrAExample() {
    std::cout << "\n========== ENHANCED SGR A* SUPERMASSIVE BLACK HOLE UQFF EXAMPLE ==========\n\n";
    
    SMBHSgrAStar sgrA;
    
    // Step 1: Initial state
    std::cout << "Step 1: Initial Configuration\n";
    sgrA.printParameters();
    double t0 = 0.0;
    std::cout << "g_SgrA(t=0) = " << sgrA.compute_g_SgrA(t0) << " m/s^2\n\n";
    
    // Step 2: Time evolution (0, 1, 3, 5, 9 Gyr)
    std::cout << "Step 2: Time Evolution (0, 1, 3, 5, 9 Gyr)\n";
    for (double t_gyr : {0.0, 1.0, 3.0, 5.0, 9.0}) {
        double t = t_gyr * 1e9 * 3.156e7;
        std::cout << "  t = " << t_gyr << " Gyr: g = " << sgrA.compute_g_SgrA(t) 
                  << " m/s^2, M(t) = " << (sgrA.M_t(t)/1.989e30) << " M_sun, B(t) = " << sgrA.B_t(t) << " T\n";
    }
    std::cout << "\n";
    
    // Step 3: Accretion scaling
    std::cout << "Step 3: Accretion Scaling (M_dot x1.5, tau_acc x0.8)\n";
    sgrA.expandAccretionScale(1.5, 0.8);
    double t_test = 4.5e9 * 3.156e7;
    std::cout << "After expansion: M_dot_0 = " << sgrA.getVariable("M_dot_0") << ", tau_acc = " << sgrA.getVariable("tau_acc") << " s\n";
    std::cout << "g_SgrA(t=4.5 Gyr) = " << sgrA.compute_g_SgrA(t_test) << " m/s^2\n\n";
    
    // Step 4: Magnetic field scaling
    std::cout << "Step 4: Magnetic Field Scaling (B0 x1.2, tau_B x1.3)\n";
    sgrA.expandMagneticScale(1.2, 1.3);
    std::cout << "After expansion: B0_G = " << sgrA.getVariable("B0_G") << " G, tau_B = " << sgrA.getVariable("tau_B") << " s\n";
    std::cout << "g_SgrA(t=4.5 Gyr) = " << sgrA.compute_g_SgrA(t_test) << " m/s^2\n\n";
    
    // Step 5: DM & precession scaling
    std::cout << "Step 5: DM & Precession Scaling (DM factor x1.4, angle x1.1)\n";
    sgrA.expandDMPrecessionScale(1.4, 1.1);
    std::cout << "After expansion: M_DM_factor = " << sgrA.getVariable("M_DM_factor") 
              << ", precession_angle_deg = " << sgrA.getVariable("precession_angle_deg") << " deg\n";
    std::cout << "g_SgrA(t=4.5 Gyr) = " << sgrA.compute_g_SgrA(t_test) << " m/s^2\n\n";
    
    // Step 6: State save/restore
    std::cout << "Step 6: State Management\n";
    sgrA.saveState("expanded_state");
    sgrA.setVariable("M_initial", 5e6 * 1.989e30);
    std::cout << "Modified M_initial to 5e6 M_sun: g = " << sgrA.compute_g_SgrA(t_test) << " m/s^2\n";
    sgrA.restoreState("expanded_state");
    std::cout << "Restored state: M_initial = " << (sgrA.getVariable("M_initial")/1.989e30) 
              << " M_sun, g = " << sgrA.compute_g_SgrA(t_test) << " m/s^2\n";
    std::cout << "Saved states: ";
    for (const auto& s : sgrA.listSavedStates()) std::cout << s << " ";
    std::cout << "\n\n";
    
    // Step 7: Sensitivity analysis
    std::cout << "Step 7: Sensitivity Analysis at t=4.5 Gyr (top 5 parameters)\n";
    auto sens = sgrA.sensitivityAnalysis(t_test, 1.0);
    std::vector<std::pair<std::string, double>> sens_vec(sens.begin(), sens.end());
    std::sort(sens_vec.begin(), sens_vec.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
    for (int i = 0; i < std::min(5, (int)sens_vec.size()); ++i) {
        std::cout << "  " << sens_vec[i].first << ": " << sens_vec[i].second << "\n";
    }
    std::cout << "\n";
    
    // Step 8: Generate variations
    std::cout << "Step 8: Generate Parameter Variations (3 variants, 10% variation)\n";
    auto variations = sgrA.generateVariations(3, 10.0);
    for (size_t i = 0; i < variations.size(); ++i) {
        std::cout << "  Variant " << (i+1) << ": M_initial = " << (variations[i]["M_initial"]/1.989e30) 
                  << " M_sun, B0_G = " << variations[i]["B0_G"] << " G\n";
    }
    std::cout << "\n";
    
    // Step 9: Batch transformation
    std::cout << "Step 9: Batch Transform (scale accretion parameters by 1.15)\n";
    sgrA.transformVariableGroup({"M_dot_0", "rho_fluid"}, [](double v) { return v * 1.15; });
    std::cout << "After transform: M_dot_0 = " << sgrA.getVariable("M_dot_0") 
              << ", rho_fluid = " << sgrA.getVariable("rho_fluid") << " kg/m^3\n";
    std::cout << "g_SgrA(t=4.5 Gyr) = " << sgrA.compute_g_SgrA(t_test) << " m/s^2\n\n";
    
    // Step 10: Consistency validation
    std::cout << "Step 10: Consistency Validation\n";
    bool valid = sgrA.validateConsistency();
    std::cout << "System is " << (valid ? "VALID" : "INVALID") << "\n\n";
    
    // Step 11: Metric optimization
    std::cout << "Step 11: Optimize for Maximum g (t=0 to 10 Gyr, 100 steps)\n";
    double max_g = sgrA.optimizeForMetric([](double g) { return g; }, 0.0, 10e9 * 3.156e7, 100);
    std::cout << "Maximum g found: " << max_g << " m/s^2\n\n";
    
    // Step 12: Full system report
    std::cout << "Step 12: Full System Report at t=5 Gyr\n";
    double t_report = 5e9 * 3.156e7;
    std::cout << sgrA.generateReport(t_report) << "\n";
    
    // Step 13: Mass sweep
    std::cout << "Step 13: Initial Mass Sweep (M_initial = 3, 4.3, 5, 6 x 10^6 M_sun)\n";
    sgrA.saveState("before_sweep");
    for (double M_factor : {3.0, 4.3, 5.0, 6.0}) {
        sgrA.setVariable("M_initial", M_factor * 1e6 * 1.989e30);
        std::cout << "  M_initial = " << M_factor << " x 10^6 M_sun: g(t=4.5Gyr) = " << sgrA.compute_g_SgrA(t_test) << " m/s^2\n";
    }
    sgrA.restoreState("before_sweep");
    std::cout << "\n";
    
    // Step 14: Accretion rate sweep
    std::cout << "Step 14: Accretion Rate Sweep (M_dot_0 = 0.005, 0.01, 0.02, 0.03)\n";
    for (double M_dot : {0.005, 0.01, 0.02, 0.03}) {
        sgrA.setVariable("M_dot_0", M_dot);
        std::cout << "  M_dot_0 = " << M_dot << ": M(t=4.5Gyr) = " << (sgrA.M_t(t_test)/1.989e30) 
                  << " M_sun, g = " << sgrA.compute_g_SgrA(t_test) << " m/s^2\n";
    }
    sgrA.restoreState("before_sweep");
    std::cout << "\n";
    
    // Step 15: Spin factor sweep
    std::cout << "Step 15: Spin Factor Sweep (spin_factor = 0.1, 0.3, 0.5, 0.7)\n";
    for (double spin : {0.1, 0.3, 0.5, 0.7}) {
        sgrA.setVariable("spin_factor", spin);
        std::cout << "  spin_factor = " << spin << ": Omega(t=4.5Gyr) = " << sgrA.Omega_t(t_test) 
                  << " rad/s, g = " << sgrA.compute_g_SgrA(t_test) << " m/s^2\n";
    }
    sgrA.restoreState("before_sweep");
    std::cout << "\n";
    
    // Step 16: Precession angle sweep
    std::cout << "Step 16: Precession Angle Sweep (angle = 0, 15, 30, 45, 60 deg)\n";
    for (double angle : {0.0, 15.0, 30.0, 45.0, 60.0}) {
        sgrA.setVariable("precession_angle_deg", angle);
        std::cout << "  angle = " << angle << " deg: g(t=4.5Gyr) = " << sgrA.compute_g_SgrA(t_test) << " m/s^2\n";
    }
    sgrA.restoreState("before_sweep");
    std::cout << "\n";
    
    // Step 17: DM factor sweep
    std::cout << "Step 17: Dark Matter Factor Sweep (M_DM_factor = 0.0, 0.1, 0.2, 0.3)\n";
    for (double dm : {0.0, 0.1, 0.2, 0.3}) {
        sgrA.setVariable("M_DM_factor", dm);
        std::cout << "  M_DM_factor = " << dm << ": g(t=4.5Gyr) = " << sgrA.compute_g_SgrA(t_test) << " m/s^2\n";
    }
    sgrA.restoreState("before_sweep");
    std::cout << "\n";
    
    // Step 18: Export final state
    std::cout << "Step 18: Export Final State\n";
    std::cout << sgrA.exportState() << "\n";
    
    std::cout << "========== ENHANCED EXAMPLE COMPLETE ==========\n\n";
}