/**
 * ================================================================================================
 * Header: MagnetarSGR1745_2900.h
 *
 * Description: C++ Module for SGR 1745-2900 Magnetar Class
 *              This is the second module in a series of 500+ code files for the Universal Quantum
 *              Field Framework (UQFF) simulations, focusing on magnetar evolution and gravity
 *              equations derived from Chandra X-ray Observatory datasets, high-energy lab simulations,
 *              and UQFF refinements (dated May 11, 2025, updated for full term inclusion on October 08, 2025).
 *
 * Purpose: Encapsulates the Master Universal Gravity Equation (MUGE) for SGR 1745-2900 magnetar
 *          evolution, including black hole proximity (Sgr A*), magnetic energy, and outburst decay.
 *          Includes ALL terms: base gravity, cosmic expansion (H(z)), magnetic decay, BH influence,
 *          UQFF Ug components with f_sc, Lambda, quantum uncertainty, EM, fluid, oscillatory waves,
 *          DM/density perturbations, magnetic energy (effective g), and decay power (cumulative energy effective g).
 *          Supports dynamic variable updates for all parameters.
 *
 * Integration: Designed for inclusion in base program 'ziqn233h.cpp' (not present here).
 *              Instantiate class in main: MagnetarSGR1745_2900 mag;
 *              Compute: double g = mag.compute_g_Magnetar(t);
 *
 * Key Features:
 *   - Default values from UQFF document, with numerical computations for H(z), M_mag, etc.
 *   - Units handled: Energy terms converted to effective acceleration via / (M * r); decay as cumulative energy / (M * r).
 *   - Setter methods for updates: setVar(double new_val) or addToVar(double delta)/subtractFromVar(double delta).
 *   - Computes g_Magnetar(r, t) with every term explicitly included.
 *
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: October 08, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

#ifndef MAGNETAR_SGR1745_2900_H
#define MAGNETAR_SGR1745_2900_H

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

class MagnetarSGR1745_2900 {
private:
    // Core parameters (mutable for updates)
    double G;               // Gravitational constant
    double M;               // Magnetar mass
    double r;               // Radius
    double Hz;              // Hubble parameter at z (s^-1)
    double B0;              // Initial magnetic field
    double tau_B;           // B decay timescale (s) - not used in this eq, but for consistency
    double B_crit;          // Critical B field
    double Lambda;          // Cosmological constant
    double c_light;         // Speed of light
    double q_charge;        // Charge (proton)
    double v_surf;          // Surface velocity
    double f_sc;            // Superconductive factor (computed as 1 - B/B_crit)
    double rho_vac_UA;      // UA vacuum density - not used here
    double rho_vac_SCm;     // SCm vacuum density - not used here
    double P_init;          // Initial rotation period (s)
    double tau_Omega;       // Omega decay timescale (s)
    double scale_EM;        // EM scaling factor
    double proton_mass;     // Proton mass for EM acceleration
    double M_BH;            // Black hole mass
    double r_BH;            // Distance to black hole
    double mu0;             // Vacuum permeability
    double L0_W;            // Initial luminosity (W)
    double tau_decay;       // Decay timescale (s)

    // Additional parameters for full inclusion of terms
    double hbar;            // Reduced Planck's constant
    double t_Hubble;        // Hubble time (s)
    double t_Hubble_gyr;    // Hubble time in Gyr
    double delta_x;         // Position uncertainty (m)
    double delta_p;         // Momentum uncertainty (kg m/s)
    double integral_psi;    // Wavefunction integral approximation
    double rho_fluid;       // Fluid density (kg/m^3)
    double A_osc;           // Oscillatory amplitude (m/s^2)
    double k_osc;           // Wave number (1/m)
    double omega_osc;       // Angular frequency (rad/s)
    double x_pos;           // Position for oscillation (m)
    double M_DM_factor;     // Dark matter mass fraction
    double delta_rho_over_rho; // Density perturbation fraction

    // Computed caches (updated on demand)
    double ug1_base;        // Cached Ug1 = G*M/r^2
    double B;               // Current B (static in this model)

public:
    // Constructor with default UQFF values
    MagnetarSGR1745_2900() {
        initializeDefaults();
    }

    // Destructor (empty)
    ~MagnetarSGR1745_2900() {}

    // Initialization method (called in constructor)
    void initializeDefaults() {
        G = 6.6743e-11;
        M = 1.4 * 1.989e30;
        r = 1e4;
        Hz = 2.269e-18;  // Computed H(z)
        B0 = 2e10;
        B = B0;  // Static for this model
        tau_B = 4000 * 3.15576e7;  // Default, not used
        B_crit = 1e11;
        Lambda = 1.1e-52;
        c_light = 3e8;
        q_charge = 1.602e-19;
        v_surf = 1e6;
        f_sc = 1 - (B / B_crit);  // Initial
        rho_vac_UA = 7.09e-36;
        rho_vac_SCm = 7.09e-37;
        P_init = 3.76;  // Pulse period
        tau_Omega = 10000 * 3.15576e7;
        scale_EM = 1e-12;
        proton_mass = 1.673e-27;
        M_BH = 4e6 * 1.989e30;
        r_BH = 2.83e16;
        mu0 = 4 * M_PI * 1e-7;
        L0_W = 5e28;  // 5e35 erg/s = 5e28 W
        tau_decay = 3.5 * 365.25 * 24 * 3600;  // 3.5 years in s

        // Full terms defaults
        hbar = 1.0546e-34;
        t_Hubble = 13.8e9 * 3.15576e7;
        t_Hubble_gyr = 13.8;
        delta_x = 1e-10;
        delta_p = hbar / delta_x;
        integral_psi = 1.0;
        rho_fluid = 1e17;
        A_osc = 1e10;
        k_osc = 1.0 / r;
        omega_osc = 2 * M_PI / P_init;
        x_pos = r;
        M_DM_factor = 0.1;
        delta_rho_over_rho = 1e-5;

        updateCache();
    }

    // Cache update for efficiency (call after parameter changes)
    void updateCache() {
        ug1_base = (G * M) / (r * r);
        f_sc = 1 - (B / B_crit);
    }

    // Universal setter for any variable (by name, for flexibility)
    bool setVariable(const std::string& varName, double newValue) {
        if (varName == "G") { G = newValue; }
        else if (varName == "M") { M = newValue; }
        else if (varName == "r") { r = newValue; }
        else if (varName == "Hz") { Hz = newValue; }
        else if (varName == "B0") { B0 = newValue; B = newValue; }
        else if (varName == "tau_B") { tau_B = newValue; }
        else if (varName == "B_crit") { B_crit = newValue; }
        else if (varName == "Lambda") { Lambda = newValue; }
        else if (varName == "c_light") { c_light = newValue; }
        else if (varName == "q_charge") { q_charge = newValue; }
        else if (varName == "v_surf") { v_surf = newValue; }
        else if (varName == "f_sc") { f_sc = newValue; }
        else if (varName == "rho_vac_UA") { rho_vac_UA = newValue; }
        else if (varName == "rho_vac_SCm") { rho_vac_SCm = newValue; }
        else if (varName == "P_init") { P_init = newValue; }
        else if (varName == "tau_Omega") { tau_Omega = newValue; }
        else if (varName == "scale_EM") { scale_EM = newValue; }
        else if (varName == "proton_mass") { proton_mass = newValue; }
        else if (varName == "M_BH") { M_BH = newValue; }
        else if (varName == "r_BH") { r_BH = newValue; }
        else if (varName == "mu0") { mu0 = newValue; }
        else if (varName == "L0_W") { L0_W = newValue; }
        else if (varName == "tau_decay") { tau_decay = newValue; }
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
        else if (varName == "M") return M;
        else if (varName == "r") return r;
        else if (varName == "Hz") return Hz;
        else if (varName == "B0") return B0;
        else if (varName == "tau_B") return tau_B;
        else if (varName == "B_crit") return B_crit;
        else if (varName == "Lambda") return Lambda;
        else if (varName == "c_light") return c_light;
        else if (varName == "q_charge") return q_charge;
        else if (varName == "v_surf") return v_surf;
        else if (varName == "f_sc") return f_sc;
        else if (varName == "rho_vac_UA") return rho_vac_UA;
        else if (varName == "rho_vac_SCm") return rho_vac_SCm;
        else if (varName == "P_init") return P_init;
        else if (varName == "tau_Omega") return tau_Omega;
        else if (varName == "scale_EM") return scale_EM;
        else if (varName == "proton_mass") return proton_mass;
        else if (varName == "M_BH") return M_BH;
        else if (varName == "r_BH") return r_BH;
        else if (varName == "mu0") return mu0;
        else if (varName == "L0_W") return L0_W;
        else if (varName == "tau_decay") return tau_decay;
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
        else {
            std::cerr << "Error: Unknown variable '" << varName << "'." << std::endl;
            return 0.0;
        }
    }

    // B(t) - static B for this model
    double B_t(double /*t*/) const {
        return B;
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
    double compute_Ug() const {
        double Ug1 = ug1_base;
        double Ug2 = 0.0;
        double Ug3 = 0.0;
        double Ug4 = Ug1 * f_sc;
        return Ug1 + Ug2 + Ug3 + Ug4;
    }

    // Volume computation
    double compute_V() const {
        return (4.0 / 3.0) * M_PI * r * r * r;
    }

    // Magnetic energy M_mag (J)
    double compute_M_mag() const {
        double V = compute_V();
        return (B_t(0) * B_t(0) / (2 * mu0)) * V;
    }

    // Cumulative decay energy up to t (J)
    double compute_cumulative_D(double t) const {
        double exp_term = exp(-t / tau_decay);
        return L0_W * tau_decay * (1 - exp_term);
    }

    // Main MUGE computation (includes ALL terms)
    double compute_g_Magnetar(double t) const {
        if (t < 0) {
            std::cerr << "Error: Time t must be non-negative." << std::endl;
            return 0.0;
        }

        double Bt = B_t(t);
        double dOdt = dOmega_dt(t);

        // f_sc update
        double current_f_sc = 1 - (Bt / B_crit);

        // Term 1: Base + H(z) + B corrections
        double corr_H = 1 + Hz * t;
        double corr_B = current_f_sc;
        double term1 = ug1_base * corr_H * corr_B;

        // BH term
        double term_BH = (G * M_BH) / (r_BH * r_BH);

        // Term 2: UQFF Ug 
        double term2 = compute_Ug();

        // Term 3: Lambda
        double term3 = (Lambda * c_light * c_light) / 3.0;

        // Term 4: Scaled EM (v x B magnitude)
        double cross_vB = v_surf * Bt;
        double em_base = (q_charge * cross_vB) / proton_mass;
        double term4 = em_base * scale_EM;  // UA not used here

        // Term 5: GW (assumed same as previous)
        double gw_prefactor = (G * M * M) / (pow(c_light, 4) * r);
        double term5 = gw_prefactor * (dOdt * dOdt);

        // Quantum uncertainty term
        double sqrt_unc = sqrt(delta_x * delta_p);
        double term_q = (hbar / sqrt_unc) * integral_psi * (2 * M_PI / t_Hubble);

        // Fluid term (effective acceleration)
        double V = compute_V();
        double term_fluid = (rho_fluid * V * ug1_base) / M;

        // Oscillatory terms (real parts)
        double term_osc1 = 2 * A_osc * cos(k_osc * x_pos) * cos(omega_osc * t);
        double arg = k_osc * x_pos - omega_osc * t;
        double term_osc2 = (2 * M_PI / t_Hubble_gyr) * A_osc * cos(arg);
        double term_osc = term_osc1 + term_osc2;

        // DM and density perturbation term (converted to acceleration)
        double M_dm = M * M_DM_factor;
        double pert1 = delta_rho_over_rho;
        double pert2 = 3 * G * M / (r * r * r);
        double term_dm_force_like = (M + M_dm) * (pert1 + pert2);
        double term_DM = term_dm_force_like / M;

        // Magnetic energy term (effective g)
        double M_mag = compute_M_mag();
        double term_mag = M_mag / (M * r);

        // Decay term (cumulative energy effective g)
        double cum_D = compute_cumulative_D(t);
        double term_decay = cum_D / (M * r);

        // Total g_Magnetar (all terms summed)
        return term1 + term_BH + term2 + term3 + term4 + term5 + term_q + term_fluid + term_osc + term_DM + term_mag + term_decay;
    }

    // Debug/Output method (for transparency in base program)
    void printParameters(std::ostream& os = std::cout) const {
        os << std::fixed << std::setprecision(3);
        os << "SGR 1745-2900 Parameters:" << std::endl;
        os << "G: " << G << ", M: " << M << ", r: " << r << std::endl;
        os << "Hz: " << Hz << ", B: " << B << ", M_BH: " << M_BH << ", r_BH: " << r_BH << std::endl;
        os << "L0_W: " << L0_W << ", tau_decay: " << tau_decay << std::endl;
        os << "f_sc: " << f_sc << ", rho_fluid: " << rho_fluid << ", M_DM_factor: " << M_DM_factor << std::endl;
        os << "A_osc: " << A_osc << ", delta_rho_over_rho: " << delta_rho_over_rho << std::endl;
        double M_mag = compute_M_mag();
        os << "M_mag (J): " << M_mag << ", ug1_base: " << ug1_base << std::endl;
    }

    // Example computation at t=1 year (for testing)
    double exampleAtOneYear() const {
        double t_example = 1.0 * 365.25 * 24 * 3600;
        return compute_g_Magnetar(t_example);
    }

    // ===== ENHANCED DYNAMIC CAPABILITIES =====
    // Variable Management
    void createVariable(const std::string& name, double value);
    void removeVariable(const std::string& name);
    void cloneVariable(const std::string& source, const std::string& dest);
    std::vector<std::string> listVariables();
    std::string getSystemName();

    // Batch Operations
    void transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func);
    void scaleVariableGroup(const std::vector<std::string>& names, double factor);

    // Self-Expansion (domain-specific)
    void expandParameterSpace(double scale_factor);
    void expandMagneticScale(double scale_factor);
    void expandDecayScale(double scale_factor);
    void expandBlackHoleScale(double scale_factor);

    // Self-Refinement
    void autoRefineParameters(double tolerance);
    void calibrateToObservations(const std::map<std::string, double>& obs_data);
    void optimizeForMetric(const std::string& var_name, double target_value, int iterations);

    // Parameter Exploration
    std::vector<std::map<std::string, double>> generateVariations(int n_variations);

    // Adaptive Evolution
    void mutateParameters(double mutation_rate);
    void evolveSystem(int generations, std::function<double()> fitness_function);

    // State Management
    void saveState(const std::string& label);
    void restoreState(const std::string& label);
    std::vector<std::string> listSavedStates();
    std::string exportState(double t);

    // System Analysis
    std::map<std::string, double> sensitivityAnalysis(const std::string& param, double t, double delta);
    std::string generateReport(double t);
    bool validateConsistency();
    void autoCorrectAnomalies();

private:
    std::map<std::string, double> extra_variables;
};

#endif // MAGNETAR_SGR1745_2900_H

// ===== IMPLEMENTATION OF ENHANCED DYNAMIC CAPABILITIES =====

namespace {
    std::map<std::string, std::map<std::string, double>> magnetar_saved_states;
}

// Variable Management
void MagnetarSGR1745_2900::createVariable(const std::string& name, double value) {
    extra_variables[name] = value;
}

void MagnetarSGR1745_2900::removeVariable(const std::string& name) {
    extra_variables.erase(name);
}

void MagnetarSGR1745_2900::cloneVariable(const std::string& source, const std::string& dest) {
    if (extra_variables.find(source) != extra_variables.end()) {
        extra_variables[dest] = extra_variables[source];
    } else {
        // Try built-in variable
        double val = getVariable(source);
        extra_variables[dest] = val;
    }
}

std::vector<std::string> MagnetarSGR1745_2900::listVariables() {
    std::vector<std::string> names = {
        "G", "M", "r", "Hz", "B0", "tau_B", "B_crit", "Lambda", "c_light",
        "q_charge", "v_surf", "f_sc", "rho_vac_UA", "rho_vac_SCm", "P_init",
        "tau_Omega", "scale_EM", "proton_mass", "M_BH", "r_BH", "mu0", "L0_W",
        "tau_decay", "hbar", "t_Hubble", "t_Hubble_gyr", "delta_x", "delta_p",
        "integral_psi", "rho_fluid", "A_osc", "k_osc", "omega_osc", "x_pos",
        "M_DM_factor", "delta_rho_over_rho"
    };
    for (const auto& pair : extra_variables) {
        names.push_back(pair.first);
    }
    return names;
}

std::string MagnetarSGR1745_2900::getSystemName() {
    return "SGR 1745-2900 Magnetar near Sgr A* - Full UQFF & SM Integration";
}

// Batch Operations
void MagnetarSGR1745_2900::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        double current = getVariable(name);
        setVariable(name, func(current));
    }
}

void MagnetarSGR1745_2900::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// Self-Expansion (domain-specific for Magnetar)
void MagnetarSGR1745_2900::expandParameterSpace(double scale_factor) {
    // Scale core parameters
    r *= scale_factor;
    M *= scale_factor;
    L0_W *= scale_factor;
    updateCache();
}

void MagnetarSGR1745_2900::expandMagneticScale(double scale_factor) {
    // Scale magnetic parameters (B0, B, tau_B)
    B0 *= scale_factor;
    B *= scale_factor;
    tau_B *= scale_factor;
    updateCache();
}

void MagnetarSGR1745_2900::expandDecayScale(double scale_factor) {
    // Scale decay parameters (L0_W, tau_decay)
    L0_W *= scale_factor;
    tau_decay *= scale_factor;
}

void MagnetarSGR1745_2900::expandBlackHoleScale(double scale_factor) {
    // Scale black hole interaction parameters (M_BH, r_BH)
    M_BH *= scale_factor;
    r_BH *= scale_factor;
}

// Self-Refinement
void MagnetarSGR1745_2900::autoRefineParameters(double tolerance) {
    // Enforce physical constraints
    if (M <= 0) M = 1.4 * 1.989e30;
    if (r <= 0) r = 1e4;
    if (B0 <= 0) B0 = 2e10;
    if (B <= 0) B = B0;
    if (L0_W <= 0) L0_W = 5e28;
    if (tau_decay <= 0) tau_decay = 3.5 * 365.25 * 24 * 3600;
    if (M_BH <= 0) M_BH = 4e6 * 1.989e30;
    if (r_BH <= 0) r_BH = 2.83e16;
    if (rho_fluid <= 0) rho_fluid = 1e17;
    updateCache();
}

void MagnetarSGR1745_2900::calibrateToObservations(const std::map<std::string, double>& obs_data) {
    for (const auto& obs : obs_data) {
        setVariable(obs.first, obs.second);
    }
    updateCache();
}

void MagnetarSGR1745_2900::optimizeForMetric(const std::string& var_name, double target_value, int iterations) {
    double current = getVariable(var_name);
    if (current == 0.0 && var_name != "unknown") {
        return;  // Variable exists but is zero or unknown
    }
    
    double best_value = current;
    double best_error = std::abs(current - target_value);
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.9, 1.1);
    
    for (int i = 0; i < iterations; ++i) {
        double test_value = current * dis(gen);
        double error = std::abs(test_value - target_value);
        if (error < best_error) {
            best_error = error;
            best_value = test_value;
        }
    }
    setVariable(var_name, best_value);
}

// Parameter Exploration
std::vector<std::map<std::string, double>> MagnetarSGR1745_2900::generateVariations(int n_variations) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.8, 1.2);
    
    for (int i = 0; i < n_variations; ++i) {
        std::map<std::string, double> variation;
        variation["B0"] = B0 * dis(gen);
        variation["L0_W"] = L0_W * dis(gen);
        variation["tau_decay"] = tau_decay * dis(gen);
        variation["M_BH"] = M_BH * dis(gen);
        variation["r_BH"] = r_BH * dis(gen);
        variations.push_back(variation);
    }
    return variations;
}

// Adaptive Evolution
void MagnetarSGR1745_2900::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-mutation_rate, mutation_rate);
    
    B0 *= (1.0 + dis(gen));
    B *= (1.0 + dis(gen));
    L0_W *= (1.0 + dis(gen));
    tau_decay *= (1.0 + dis(gen));
    updateCache();
}

void MagnetarSGR1745_2900::evolveSystem(int generations, std::function<double()> fitness_function) {
    double best_fitness = fitness_function();
    
    // Save current state
    double saved_B0 = B0, saved_B = B, saved_L0W = L0_W, saved_tau = tau_decay;
    
    for (int gen = 0; gen < generations; ++gen) {
        mutateParameters(0.1);
        double fitness = fitness_function();
        if (fitness > best_fitness) {
            best_fitness = fitness;
            saved_B0 = B0;
            saved_B = B;
            saved_L0W = L0_W;
            saved_tau = tau_decay;
        } else {
            // Revert
            B0 = saved_B0;
            B = saved_B;
            L0_W = saved_L0W;
            tau_decay = saved_tau;
            updateCache();
        }
    }
}

// State Management
void MagnetarSGR1745_2900::saveState(const std::string& label) {
    std::map<std::string, double> state;
    auto var_names = listVariables();
    for (const auto& name : var_names) {
        state[name] = getVariable(name);
    }
    magnetar_saved_states[label] = state;
}

void MagnetarSGR1745_2900::restoreState(const std::string& label) {
    if (magnetar_saved_states.find(label) != magnetar_saved_states.end()) {
        const auto& state = magnetar_saved_states[label];
        for (const auto& pair : state) {
            setVariable(pair.first, pair.second);
        }
    }
}

std::vector<std::string> MagnetarSGR1745_2900::listSavedStates() {
    std::vector<std::string> labels;
    for (const auto& pair : magnetar_saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string MagnetarSGR1745_2900::exportState(double t) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(6);
    oss << "SGR 1745-2900 State Export at t=" << t << " s (" << t/(365.25*24*3600) << " yr):\n";
    oss << "M=" << M << " kg (" << M/1.989e30 << " M_sun)\n";
    oss << "r=" << r << " m (" << r/1e3 << " km)\n";
    oss << "B0=" << B0 << " T, B(t)=" << B_t(t) << " T\n";
    oss << "L0_W=" << L0_W << " W, tau_decay=" << tau_decay << " s (" << tau_decay/(365.25*24*3600) << " yr)\n";
    oss << "M_BH=" << M_BH << " kg (" << M_BH/1.989e30 << " M_sun), r_BH=" << r_BH << " m\n";
    oss << "rho_fluid=" << rho_fluid << " kg/m³, A_osc=" << A_osc << " m/s²\n";
    oss << "M_mag=" << compute_M_mag() << " J, cumulative_D=" << compute_cumulative_D(t) << " J\n";
    oss << "g_total=" << compute_g_Magnetar(t) << " m/s²\n";
    return oss.str();
}

// System Analysis
std::map<std::string, double> MagnetarSGR1745_2900::sensitivityAnalysis(const std::string& param, double t, double delta) {
    std::map<std::string, double> result;
    double original = getVariable(param);
    if (original == 0.0 && param != "unknown") {
        result["error"] = -1;
        return result;
    }
    
    double g_original = compute_g_Magnetar(t);
    
    setVariable(param, original * (1.0 + delta));
    double g_plus = compute_g_Magnetar(t);
    
    setVariable(param, original * (1.0 - delta));
    double g_minus = compute_g_Magnetar(t);
    
    setVariable(param, original);
    
    result["dg/d" + param] = (g_plus - g_minus) / (2.0 * delta * original);
    result["g_original"] = g_original;
    result["g_plus"] = g_plus;
    result["g_minus"] = g_minus;
    
    return result;
}

std::string MagnetarSGR1745_2900::generateReport(double t) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(6);
    oss << "===== SGR 1745-2900 MAGNETAR UQFF Module Report (t=" << t << " s) =====\n\n";
    oss << "System: " << getSystemName() << "\n\n";
    
    oss << "Core Parameters:\n";
    oss << "  M = " << M << " kg (" << M/1.989e30 << " M_sun)\n";
    oss << "  r = " << r << " m (" << r/1e3 << " km)\n";
    oss << "  B0 = " << B0 << " T, B(t) = " << B_t(t) << " T\n";
    oss << "  B_crit = " << B_crit << " T\n";
    oss << "  P_init = " << P_init << " s (pulse period)\n";
    oss << "  L0_W = " << L0_W << " W (initial luminosity)\n";
    oss << "  tau_decay = " << tau_decay << " s (" << tau_decay/(365.25*24*3600) << " yr)\n";
    oss << "  M_BH = " << M_BH << " kg (" << M_BH/1.989e30 << " M_sun, Sgr A*)\n";
    oss << "  r_BH = " << r_BH << " m (distance to Sgr A*)\n\n";
    
    double Bt = B_t(t);
    double dOdt = dOmega_dt(t);
    double current_f_sc = 1 - (Bt / B_crit);
    double corr_H = 1 + Hz * t;
    double corr_B = current_f_sc;
    double term1 = ug1_base * corr_H * corr_B;
    double term_BH = (G * M_BH) / (r_BH * r_BH);
    double term2 = compute_Ug();
    double term3 = (Lambda * c_light * c_light) / 3.0;
    double cross_vB = v_surf * Bt;
    double em_base = (q_charge * cross_vB) / proton_mass;
    double term4 = em_base * scale_EM;
    double gw_prefactor = (G * M * M) / (pow(c_light, 4) * r);
    double term5 = gw_prefactor * (dOdt * dOdt);
    double sqrt_unc = sqrt(delta_x * delta_p);
    double term_q = (hbar / sqrt_unc) * integral_psi * (2 * M_PI / t_Hubble);
    double V = compute_V();
    double term_fluid = (rho_fluid * V * ug1_base) / M;
    double term_osc1 = 2 * A_osc * cos(k_osc * x_pos) * cos(omega_osc * t);
    double arg = k_osc * x_pos - omega_osc * t;
    double term_osc2 = (2 * M_PI / t_Hubble_gyr) * A_osc * cos(arg);
    double term_osc = term_osc1 + term_osc2;
    double M_dm = M * M_DM_factor;
    double pert1 = delta_rho_over_rho;
    double pert2 = 3 * G * M / (r * r * r);
    double term_dm_force_like = (M + M_dm) * (pert1 + pert2);
    double term_DM = term_dm_force_like / M;
    double M_mag = compute_M_mag();
    double term_mag = M_mag / (M * r);
    double cum_D = compute_cumulative_D(t);
    double term_decay = cum_D / (M * r);
    double g_total = compute_g_Magnetar(t);
    
    oss << "Term Breakdown:\n";
    oss << "  Base gravity (H(z), B corrections) = " << term1 << " m/s²\n";
    oss << "  BH term (Sgr A*) = " << term_BH << " m/s²\n";
    oss << "  Ug_sum (Ug1+Ug2+Ug3+Ug4) = " << term2 << " m/s²\n";
    oss << "  Lambda_term = " << term3 << " m/s²\n";
    oss << "  EM_term (v×B) = " << term4 << " m/s²\n";
    oss << "  GW_term (dΩ/dt)² = " << term5 << " m/s²\n";
    oss << "  Quantum_term = " << term_q << " m/s²\n";
    oss << "  Fluid_term = " << term_fluid << " m/s²\n";
    oss << "  Oscillatory_term = " << term_osc << " m/s²\n";
    oss << "  DM_term = " << term_DM << " m/s²\n";
    oss << "  Magnetic_energy_term = " << term_mag << " m/s² [M_mag=" << M_mag << " J]\n";
    oss << "  Decay_energy_term = " << term_decay << " m/s² [cumulative_D=" << cum_D << " J]\n\n";
    oss << "TOTAL g = " << g_total << " m/s²\n\n";
    
    oss << "Physics Notes:\n";
    oss << "- SGR 1745-2900 is a magnetar near Sgr A* supermassive black hole\n";
    oss << "- Ultra-strong magnetic field B~" << B0 << " T with superconductivity correction f_sc=" << current_f_sc << "\n";
    oss << "- X-ray outburst decay L(t)=L0*exp(-t/tau), cumulative energy affects gravity\n";
    oss << "- Black hole proximity (r_BH=" << r_BH << " m) adds tidal term\n";
    oss << "- Full UQFF+SM: gravity, BH, Ug1-4, Lambda, EM, GW, quantum, fluid, oscillatory, DM, magnetic energy, decay\n";
    oss << "- Rotation period P=" << P_init << " s, slowing with tau_Omega=" << tau_Omega << " s\n";
    oss << "- Dense NS interior: rho_fluid=" << rho_fluid << " kg/m³\n\n";
    
    return oss.str();
}

bool MagnetarSGR1745_2900::validateConsistency() {
    bool valid = true;
    if (M <= 0) valid = false;
    if (r <= 0) valid = false;
    if (B0 <= 0) valid = false;
    if (B <= 0) valid = false;
    if (L0_W <= 0) valid = false;
    if (tau_decay <= 0) valid = false;
    if (M_BH <= 0) valid = false;
    if (r_BH <= 0) valid = false;
    if (rho_fluid <= 0) valid = false;
    return valid;
}

void MagnetarSGR1745_2900::autoCorrectAnomalies() {
    if (M <= 0) M = 1.4 * 1.989e30;
    if (r <= 0) r = 1e4;
    if (B0 <= 0) B0 = 2e10;
    if (B <= 0) B = B0;
    if (L0_W <= 0) L0_W = 5e28;
    if (tau_decay <= 0) tau_decay = 3.5 * 365.25 * 24 * 3600;
    if (M_BH <= 0) M_BH = 4e6 * 1.989e30;
    if (r_BH <= 0) r_BH = 2.83e16;
    if (rho_fluid <= 0) rho_fluid = 1e17;
    updateCache();
}

// Enhanced example usage demonstration
void enhanced_magnetar_example() {
    MagnetarSGR1745_2900 mag;
    double t_1yr = 365.25 * 24 * 3600;  // 1 year in seconds
    
    std::cout << "===== ENHANCED SGR 1745-2900 MAGNETAR MODULE DEMONSTRATION =====\n\n";
    
    // Step 1: Variable management
    std::cout << "Step 1: Variable Management\n";
    mag.createVariable("custom_outburst_factor", 1.5);
    mag.cloneVariable("B0", "B0_backup");
    std::vector<std::string> vars = mag.listVariables();
    std::cout << "Total variables: " << vars.size() << "\n";
    std::cout << "System: " << mag.getSystemName() << "\n\n";
    
    // Step 2: Batch scaling
    std::cout << "Step 2: Batch Scaling (Magnetic parameters)\n";
    mag.scaleVariableGroup({"B0", "L0_W", "tau_decay"}, 1.1);
    std::cout << "Scaled B0, L0_W, tau_decay by 1.1\n\n";
    
    // Step 3: Self-expansion (different physics domains)
    std::cout << "Step 3: Self-Expansion\n";
    mag.expandMagneticScale(1.05);  // Magnetic +5%
    std::cout << "Expanded magnetic scale +5%\n";
    mag.expandDecayScale(1.03);  // Decay +3%
    std::cout << "Expanded decay scale +3%\n";
    mag.expandBlackHoleScale(1.02);  // BH interaction +2%
    std::cout << "Expanded black hole scale +2%\n\n";
    
    // Step 4: Self-refinement
    std::cout << "Step 4: Self-Refinement\n";
    mag.autoRefineParameters(1e-10);
    std::cout << "Auto-refined parameters\n";
    std::map<std::string, double> obs_data = {
        {"B0", 2.2e10},
        {"L0_W", 5.5e28},
        {"tau_decay", 3.6 * 365.25 * 24 * 3600}
    };
    mag.calibrateToObservations(obs_data);
    std::cout << "Calibrated to observations\n\n";
    
    // Step 5: Optimize for specific metric
    std::cout << "Step 5: Optimize for B0~2e10 T\n";
    mag.optimizeForMetric("B0", 2e10, 50);
    std::cout << "Optimization complete\n\n";
    
    // Step 6: Generate variations
    std::cout << "Step 6: Generate 15 Parameter Variations\n";
    auto variations = mag.generateVariations(15);
    std::cout << "Generated " << variations.size() << " variations\n\n";
    
    // Step 7: State management
    std::cout << "Step 7: State Management\n";
    mag.saveState("initial");
    mag.scaleVariableGroup({"B0", "L0_W"}, 1.2);
    mag.saveState("enhanced_magnetic");
    mag.expandDecayScale(0.8);
    mag.saveState("reduced_decay");
    std::cout << "Saved 3 states\n\n";
    
    // Step 8: Sensitivity analysis
    std::cout << "Step 8: Sensitivity Analysis (B0 at t=1yr)\n";
    mag.restoreState("initial");
    auto sensitivity = mag.sensitivityAnalysis("B0", t_1yr, 0.1);
    std::cout << "dg/dB0 = " << std::scientific << sensitivity["dg/dB0"] << " (m/s²)/T\n\n";
    
    // Step 9: System validation
    std::cout << "Step 9: System Validation\n";
    bool valid = mag.validateConsistency();
    std::cout << "System consistency: " << (valid ? "VALID" : "INVALID") << "\n";
    if (!valid) {
        mag.autoCorrectAnomalies();
        std::cout << "Auto-corrected anomalies\n";
    }
    std::cout << "\n";
    
    // Step 10: Comprehensive report
    std::cout << "Step 10: Comprehensive Report (t=1yr)\n";
    std::string report = mag.generateReport(t_1yr);
    std::cout << report << "\n";
    
    // Step 11: Adaptive evolution
    std::cout << "Step 11: Adaptive Evolution (25 generations)\n";
    auto fitness_fn = [&mag, t_1yr]() -> double {
        double g = mag.compute_g_Magnetar(t_1yr);
        return -std::abs(std::log10(std::abs(g)) + 5.0);  // Target g~1e-5 m/s²
    };
    mag.evolveSystem(25, fitness_fn);
    std::cout << "Evolution complete\n\n";
    
    // Step 12: Time evolution comparison
    std::cout << "Step 12: Time Evolution (0 to 10 years)\n";
    std::vector<double> times = {0.0, 1*t_1yr, 2*t_1yr, 3.5*t_1yr, 5*t_1yr, 10*t_1yr};
    for (double t : times) {
        double g = mag.compute_g_Magnetar(t);
        std::cout << "t=" << std::scientific << t << " s (" << t/(365.25*24*3600) << " yr): g=" << g << " m/s²\n";
    }
    std::cout << "\n";
    
    // Step 13: Magnetic field evolution
    std::cout << "Step 13: Magnetic Field B(t) Evolution\n";
    for (double t : times) {
        double Bt = mag.B_t(t);
        std::cout << "t=" << t/(365.25*24*3600) << " yr: B(t)=" << std::scientific << Bt << " T\n";
    }
    std::cout << "\n";
    
    // Step 14: Decay energy accumulation
    std::cout << "Step 14: Cumulative Decay Energy\n";
    for (double t : times) {
        double cum_D = mag.compute_cumulative_D(t);
        std::cout << "t=" << t/(365.25*24*3600) << " yr: cumulative_D=" << std::scientific << cum_D << " J\n";
    }
    std::cout << "\n";
    
    // Step 15: Multi-parameter sensitivity
    std::cout << "Step 15: Multi-Parameter Sensitivity (t=1yr)\n";
    std::vector<std::string> params = {"B0", "L0_W", "tau_decay", "M_BH", "r_BH"};
    for (const auto& param : params) {
        auto sens = mag.sensitivityAnalysis(param, t_1yr, 0.05);
        if (sens.find("error") == sens.end()) {
            std::cout << "dg/d" << param << " = " << std::scientific << sens["dg/d" + param] << "\n";
        }
    }
    std::cout << "\n";
    
    // Step 16: State restoration
    std::cout << "Step 16: State Restoration\n";
    mag.restoreState("initial");
    std::cout << "Restored initial state\n";
    double g_initial = mag.compute_g_Magnetar(t_1yr);
    std::cout << "Initial g = " << std::scientific << g_initial << " m/s²\n\n";
    
    // Step 17: Final state export
    std::cout << "Step 17: Final State Export\n";
    std::cout << mag.exportState(t_1yr) << "\n";
    
    std::cout << "===== DEMONSTRATION COMPLETE =====\n";
}

#endif // MAGNETAR_SGR1745_2900_H