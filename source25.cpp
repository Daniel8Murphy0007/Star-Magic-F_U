/**
 * ================================================================================================
 * Header: NGC1275.h
 *
 * Description: C++ Module for NGC 1275 (Magnetic Monster Perseus A) Class
 *              This is the sixteenth module in a series of 500+ code files for the Universal Quantum
 *              Field Framework (UQFF) simulations, focusing on active galactic nucleus evolution
 *              and gravity equations derived from Hubble datasets, high-energy lab simulations, and
 *              UQFF refinements (dated May 09, 2025, updated for full term inclusion on October 08, 2025).
 *
 * Purpose: Encapsulates the Master Universal Gravity Equation (MUGE) for NGC 1275 evolution.
 *          Includes ALL terms: base gravity (static M), cosmic expansion (H(z)), magnetic field B(t),
 *          filament support F(t), black hole influence, UQFF Ug components with f_TRZ, Lambda,
 *          quantum uncertainty, scaled EM with [UA], fluid dynamics, oscillatory waves,
 *          DM/density perturbations, and cooling flow term. Supports dynamic variable updates.
 *
 * Integration: Designed for inclusion in base program 'ziqn233h.cpp' (not present here).
 *              Instantiate class in main: NGC1275 ngc1275;
 *              Compute: double g = ngc1275.compute_g_NGC1275(t);
 *
 * Key Features:
 *   - Default values from UQFF document: M = 1e11 Msun, r = 1.893e21 m (200k ly), M_BH = 8e8 Msun,
 *     z = 0.0176, Hz ? 2.20e-18 s^-1, B0 = 5e-9 T, tau_B = 100 Myr, F0 = 0.1, tau_fil = 100 Myr,
 *     rho_cool = 1e-20 kg/m^3, v_cool = 3e3 m/s.
 *   - Units handled: Msun to kg, ly to m; cooling term as (rho * v_cool^2) / rho_fluid for acceleration.
 *   - Setter methods for updates: setVar(double new_val) or addToVar(double delta)/subtractFromVar(double delta).
 *   - Computes g_NGC1275(r, t) with every term explicitly included.
 *
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: October 08, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

#ifndef NGC_1275_H
#define NGC_1275_H

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

class NGC1275 {
private:
    // Core parameters (mutable for updates)
    double G;               // Gravitational constant
    double M;               // Total galaxy mass (kg)
    double r;               // Radius (m)
    double Hz;              // Hubble parameter at z (s^-1)
    double B0;              // Initial magnetic field (T)
    double tau_B;           // B decay timescale (s)
    double B_crit;          // Critical B field (T)
    double Lambda;          // Cosmological constant
    double c_light;         // Speed of light
    double q_charge;        // Charge (proton)
    double gas_v;           // Gas velocity for EM (m/s)
    double f_TRZ;           // Time-reversal factor
    double M_BH;            // Black hole mass (kg)
    double r_BH;            // Black hole radius (m)
    double F0;              // Initial filament factor
    double tau_fil;         // Filament timescale (s)
    double rho_cool;        // Cooling flow density (kg/m^3)
    double v_cool;          // Cooling flow velocity (m/s)
    double rho_fluid;       // Fluid density (kg/m^3)
    double rho_vac_UA;      // UA vacuum density (J/m^3)
    double rho_vac_SCm;     // SCm vacuum density (J/m^3)
    double scale_EM;        // EM scaling factor
    double proton_mass;     // Proton mass for EM acceleration
    double z_gal;           // Galaxy redshift

    // Additional parameters for full inclusion of terms
    double hbar;            // Reduced Planck's constant
    double t_Hubble;        // Hubble time (s)
    double delta_x;         // Position uncertainty (m)
    double delta_p;         // Momentum uncertainty (kg m/s)
    double integral_psi;    // Wavefunction integral approximation
    double A_osc;           // Oscillatory amplitude (m/s^2)
    double k_osc;           // Wave number (1/m)
    double omega_osc;       // Angular frequency (rad/s)
    double x_pos;           // Position for oscillation (m)
    double t_Hubble_gyr;    // Hubble time in Gyr
    double M_DM_factor;     // Dark matter mass fraction
    double delta_rho_over_rho; // Density perturbation fraction

    // Computed caches (updated on demand)
    double ug1_base;        // Cached Ug1 = G*M/r^2
    double g_BH;            // Cached BH acceleration

public:
    // Constructor with default UQFF values
    NGC1275() {
        initializeDefaults();
    }

    // Destructor (empty)
    ~NGC1275() {}

    // Initialization method (called in constructor)
    void initializeDefaults() {
        G = 6.6743e-11;
        double M_sun = 1.989e30;
        M = 1e11 * M_sun;
        double ly_to_m = 9.461e15;
        r = 200000.0 * ly_to_m;
        z_gal = 0.0176;
        double Hz_kms = 70 * sqrt(0.3 * pow(1 + z_gal, 3) + 0.7);  // km/s/Mpc
        Hz = (Hz_kms * 1000 / 3.086e19);  // s^-1
        B0 = 5e-9;
        tau_B = 100e6 * 3.156e7;
        B_crit = 1e11;
        Lambda = 1.1e-52;
        c_light = 3e8;
        q_charge = 1.602e-19;
        gas_v = 1e5;
        f_TRZ = 0.1;
        M_BH = 8e8 * M_sun;
        r_BH = 1e18;  // Approximate influence radius
        F0 = 0.1;
        tau_fil = 100e6 * 3.156e7;
        rho_cool = 1e-20;
        v_cool = 3e3;
        rho_fluid = 1e-20;
        rho_vac_UA = 7.09e-36;
        rho_vac_SCm = 7.09e-37;
        scale_EM = 1e-12;
        proton_mass = 1.673e-27;

        // Full terms defaults
        hbar = 1.0546e-34;
        t_Hubble = 13.8e9 * 3.156e7;
        t_Hubble_gyr = 13.8;
        delta_x = 1e-10;
        delta_p = hbar / delta_x;
        integral_psi = 1.0;
        A_osc = 1e-10;
        k_osc = 1.0 / r;
        omega_osc = 2 * M_PI / (r / c_light);
        x_pos = r;
        M_DM_factor = 0.1;
        delta_rho_over_rho = 1e-5;

        updateCache();
    }

    // Cache update for efficiency (call after parameter changes)
    void updateCache() {
        ug1_base = (G * M) / (r * r);
        g_BH = (G * M_BH) / (r_BH * r_BH);
    }

    // Universal setter for any variable (by name, for flexibility)
    bool setVariable(const std::string& varName, double newValue) {
        if (varName == "G") { G = newValue; }
        else if (varName == "M") { M = newValue; }
        else if (varName == "r") { r = newValue; }
        else if (varName == "Hz") { Hz = newValue; }
        else if (varName == "B0") { B0 = newValue; }
        else if (varName == "tau_B") { tau_B = newValue; }
        else if (varName == "B_crit") { B_crit = newValue; }
        else if (varName == "Lambda") { Lambda = newValue; }
        else if (varName == "c_light") { c_light = newValue; }
        else if (varName == "q_charge") { q_charge = newValue; }
        else if (varName == "gas_v") { gas_v = newValue; }
        else if (varName == "f_TRZ") { f_TRZ = newValue; }
        else if (varName == "M_BH") { M_BH = newValue; }
        else if (varName == "r_BH") { r_BH = newValue; }
        else if (varName == "F0") { F0 = newValue; }
        else if (varName == "tau_fil") { tau_fil = newValue; }
        else if (varName == "rho_cool") { rho_cool = newValue; }
        else if (varName == "v_cool") { v_cool = newValue; }
        else if (varName == "rho_fluid") { rho_fluid = newValue; }
        else if (varName == "rho_vac_UA") { rho_vac_UA = newValue; }
        else if (varName == "rho_vac_SCm") { rho_vac_SCm = newValue; }
        else if (varName == "scale_EM") { scale_EM = newValue; }
        else if (varName == "proton_mass") { proton_mass = newValue; }
        else if (varName == "z_gal") { z_gal = newValue; }
        // Full terms
        else if (varName == "hbar") { hbar = newValue; }
        else if (varName == "t_Hubble") { t_Hubble = newValue; }
        else if (varName == "t_Hubble_gyr") { t_Hubble_gyr = newValue; }
        else if (varName == "delta_x") { delta_x = newValue; }
        else if (varName == "delta_p") { delta_p = newValue; }
        else if (varName == "integral_psi") { integral_psi = newValue; }
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
        else if (varName == "gas_v") return gas_v;
        else if (varName == "f_TRZ") return f_TRZ;
        else if (varName == "M_BH") return M_BH;
        else if (varName == "r_BH") return r_BH;
        else if (varName == "F0") return F0;
        else if (varName == "tau_fil") return tau_fil;
        else if (varName == "rho_cool") return rho_cool;
        else if (varName == "v_cool") return v_cool;
        else if (varName == "rho_fluid") return rho_fluid;
        else if (varName == "rho_vac_UA") return rho_vac_UA;
        else if (varName == "rho_vac_SCm") return rho_vac_SCm;
        else if (varName == "scale_EM") return scale_EM;
        else if (varName == "proton_mass") return proton_mass;
        else if (varName == "z_gal") return z_gal;
        // Full terms
        else if (varName == "hbar") return hbar;
        else if (varName == "t_Hubble") return t_Hubble;
        else if (varName == "t_Hubble_gyr") return t_Hubble_gyr;
        else if (varName == "delta_x") return delta_x;
        else if (varName == "delta_p") return delta_p;
        else if (varName == "integral_psi") return integral_psi;
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

    // B(t) computation
    double B_t(double t) const {
        return B0 * exp(-t / tau_B);
    }

    // F(t) computation
    double F_t(double t) const {
        return F0 * exp(-t / tau_fil);
    }

    // Ug terms computation
    double compute_Ug(double Bt, double Ft) const {
        double Ug1 = ug1_base;
        double Ug2 = 0.0;
        double Ug3 = 0.0;
        double corr_B = 1 - Bt / B_crit;
        double Ug4 = Ug1 * corr_B;
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1 + f_TRZ) * (1 + Ft);
    }

    // Volume computation for fluid
    double compute_V() const {
        return (4.0 / 3.0) * M_PI * r * r * r;
    }

    // Main MUGE computation (includes ALL terms)
    double compute_g_NGC1275(double t) const {
        if (t < 0) {
            std::cerr << "Error: Time t must be non-negative." << std::endl;
            return 0.0;
        }

        double Bt = B_t(t);
        double Ft = F_t(t);

        // Term 1: Base + Hz + B + F corrections
        double corr_H = 1 + Hz * t;
        double corr_B = 1 - Bt / B_crit;
        double corr_F = 1 + Ft;
        double term1 = ug1_base * corr_H * corr_B * corr_F;

        // BH term
        double term_BH = g_BH;

        // Term 2: UQFF Ug with f_TRZ, B, F
        double term2 = compute_Ug(Bt, Ft);

        // Term 3: Lambda
        double term3 = (Lambda * c_light * c_light) / 3.0;

        // Term 4: Scaled EM with UA
        double cross_vB = gas_v * Bt;  // Magnitude, assuming perpendicular
        double em_base = (q_charge * cross_vB) / proton_mass;
        double corr_UA = 1 + (rho_vac_UA / rho_vac_SCm);
        double term4 = (em_base * corr_UA) * scale_EM;

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

        // Cooling flow term (pressure / density for acceleration)
        double cool_pressure = rho_cool * v_cool * v_cool;
        double term_cool = cool_pressure / rho_fluid;

        // Total g_NGC1275 (all terms summed)
        return term1 + term_BH + term2 + term3 + term4 + term_q + term_fluid + term_osc + term_DM + term_cool;
    }

    // Debug/Output method (for transparency in base program)
    void printParameters(std::ostream& os = std::cout) const {
        os << std::fixed << std::setprecision(3);
        os << "NGC 1275 Parameters:" << std::endl;
        os << "G: " << G << ", M: " << M << ", r: " << r << std::endl;
        os << "Hz: " << Hz << ", B0: " << B0 << ", tau_B: " << tau_B << std::endl;
        os << "f_TRZ: " << f_TRZ << ", M_BH: " << M_BH << ", r_BH: " << r_BH << std::endl;
        os << "F0: " << F0 << ", tau_fil: " << tau_fil << std::endl;
        os << "rho_fluid: " << rho_fluid << ", rho_cool: " << rho_cool << ", v_cool: " << v_cool << std::endl;
        os << "gas_v: " << gas_v << ", M_DM_factor: " << M_DM_factor << std::endl;
        os << "A_osc: " << A_osc << ", delta_rho_over_rho: " << delta_rho_over_rho << std::endl;
        os << "ug1_base: " << ug1_base << ", g_BH: " << g_BH << std::endl;
    }

    // Example computation at t=50 Myr (for testing)
    double exampleAt50Myr() const {
        double t_example = 50e6 * 3.156e7;
        return compute_g_NGC1275(t_example);
    }

    // ========== ENHANCED DYNAMIC CAPABILITIES (25 methods) ==========
    
    // Variable Management (5 methods)
    void createVariable(const std::string& name, double value);
    bool removeVariable(const std::string& name);
    void cloneVariable(const std::string& source, const std::string& dest);
    std::vector<std::string> listVariables() const;
    std::string getSystemName() const { return "NGC1275"; }
    
    // Batch Operations (2 methods)
    void transformVariableGroup(const std::vector<std::string>& vars, std::function<double(double)> func);
    void scaleVariableGroup(const std::vector<std::string>& vars, double factor);
    
    // Self-Expansion (4 methods - domain-specific for AGN)
    void expandParameterSpace(double scale_factor);
    void expandGalaxyScale(double M_scale, double r_scale);
    void expandBlackHoleScale(double M_BH_scale, double r_BH_scale);
    void expandMagneticFilamentScale(double B0_scale, double F0_scale);
    
    // Self-Refinement (3 methods)
    void autoRefineParameters(const std::vector<std::pair<double, double>>& observations);
    void calibrateToObservations(const std::vector<std::pair<double, double>>& obs_data);
    double optimizeForMetric(std::function<double(double)> metric, double t_start, double t_end, int steps);
    
    // Parameter Exploration (1 method)
    std::vector<std::map<std::string, double>> generateVariations(int count, double variation_percent);
    
    // Adaptive Evolution (2 methods)
    void mutateParameters(double mutation_rate);
    void evolveSystem(int generations, std::function<double(const NGC1275&)> fitness);
    
    // State Management (4 methods)
    void saveState(const std::string& label);
    bool restoreState(const std::string& label);
    std::vector<std::string> listSavedStates() const;
    std::string exportState() const;
    
    // System Analysis (4 methods)
    std::map<std::string, double> sensitivityAnalysis(double t, double perturbation);
    std::string generateReport(double t) const;
    bool validateConsistency() const;
    bool autoCorrectAnomalies();
};

#endif // NGC_1275_H

// ========== IMPLEMENTATION OF ENHANCED METHODS ==========

namespace {
    // Storage for saved states (anonymous namespace for encapsulation)
    std::map<std::string, std::map<std::string, double>> ngc1275_saved_states;
}

// Variable Management
void NGC1275::createVariable(const std::string& name, double value) {
    setVariable(name, value);
}

bool NGC1275::removeVariable(const std::string& name) {
    // Cannot truly remove core variables, but can reset to defaults
    initializeDefaults();
    return true;
}

void NGC1275::cloneVariable(const std::string& source, const std::string& dest) {
    double value = getVariable(source);
    setVariable(dest, value);
}

std::vector<std::string> NGC1275::listVariables() const {
    return {"G", "M", "r", "Hz", "B0", "tau_B", "B_crit", "Lambda", "c_light", "q_charge", 
            "gas_v", "f_TRZ", "M_BH", "r_BH", "F0", "tau_fil", "rho_cool", "v_cool", 
            "rho_fluid", "rho_vac_UA", "rho_vac_SCm", "scale_EM", "proton_mass", "z_gal",
            "hbar", "t_Hubble", "t_Hubble_gyr", "delta_x", "delta_p", "integral_psi",
            "A_osc", "k_osc", "omega_osc", "x_pos", "M_DM_factor", "delta_rho_over_rho"};
}

// Batch Operations
void NGC1275::transformVariableGroup(const std::vector<std::string>& vars, std::function<double(double)> func) {
    for (const auto& var : vars) {
        double current = getVariable(var);
        setVariable(var, func(current));
    }
}

void NGC1275::scaleVariableGroup(const std::vector<std::string>& vars, double factor) {
    transformVariableGroup(vars, [factor](double v) { return v * factor; });
}

// Self-Expansion (domain-specific for AGN)
void NGC1275::expandParameterSpace(double scale_factor) {
    std::vector<std::string> scalable = {"M", "r", "M_BH", "r_BH", "B0", "F0", "tau_B", "tau_fil", 
                                          "rho_cool", "v_cool"};
    scaleVariableGroup(scalable, scale_factor);
}

void NGC1275::expandGalaxyScale(double M_scale, double r_scale) {
    setVariable("M", getVariable("M") * M_scale);
    setVariable("r", getVariable("r") * r_scale);
}

void NGC1275::expandBlackHoleScale(double M_BH_scale, double r_BH_scale) {
    setVariable("M_BH", getVariable("M_BH") * M_BH_scale);
    setVariable("r_BH", getVariable("r_BH") * r_BH_scale);
}

void NGC1275::expandMagneticFilamentScale(double B0_scale, double F0_scale) {
    setVariable("B0", getVariable("B0") * B0_scale);
    setVariable("F0", getVariable("F0") * F0_scale);
}

// Self-Refinement
void NGC1275::autoRefineParameters(const std::vector<std::pair<double, double>>& observations) {
    calibrateToObservations(observations);
}

void NGC1275::calibrateToObservations(const std::vector<std::pair<double, double>>& obs_data) {
    if (obs_data.empty()) return;
    
    double total_error = 0.0;
    for (const auto& obs : obs_data) {
        double t = obs.first;
        double g_obs = obs.second;
        double g_model = compute_g_NGC1275(t);
        total_error += std::abs(g_model - g_obs);
    }
    
    double avg_error = total_error / obs_data.size();
    if (avg_error > 1e-9) {
        double correction = 0.95;
        setVariable("B0", getVariable("B0") * correction);
    }
}

double NGC1275::optimizeForMetric(std::function<double(double)> metric, double t_start, double t_end, int steps) {
    double best_metric = -1e100;
    double dt = (t_end - t_start) / steps;
    
    for (int i = 0; i <= steps; ++i) {
        double t = t_start + i * dt;
        double g = compute_g_NGC1275(t);
        double m = metric(g);
        if (m > best_metric) best_metric = m;
    }
    return best_metric;
}

// Parameter Exploration
std::vector<std::map<std::string, double>> NGC1275::generateVariations(int count, double variation_percent) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-variation_percent / 100.0, variation_percent / 100.0);
    
    auto vars = listVariables();
    for (int i = 0; i < count; ++i) {
        std::map<std::string, double> variant;
        for (const auto& var : vars) {
            double base = getVariable(var);
            double variation = base * (1.0 + dis(gen));
            variant[var] = variation;
        }
        variations.push_back(variant);
    }
    return variations;
}

// Adaptive Evolution
void NGC1275::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-mutation_rate, mutation_rate);
    
    auto vars = listVariables();
    for (const auto& var : vars) {
        double current = getVariable(var);
        double mutated = current * (1.0 + dis(gen));
        setVariable(var, mutated);
    }
}

void NGC1275::evolveSystem(int generations, std::function<double(const NGC1275&)> fitness) {
    double best_fitness = fitness(*this);
    saveState("evolution_best");
    
    for (int gen = 0; gen < generations; ++gen) {
        saveState("evolution_temp");
        mutateParameters(0.1);
        
        double current_fitness = fitness(*this);
        if (current_fitness > best_fitness) {
            best_fitness = current_fitness;
            saveState("evolution_best");
        } else {
            restoreState("evolution_temp");
        }
    }
    restoreState("evolution_best");
}

// State Management
void NGC1275::saveState(const std::string& label) {
    std::map<std::string, double> state;
    auto vars = listVariables();
    for (const auto& var : vars) {
        state[var] = getVariable(var);
    }
    ngc1275_saved_states[label] = state;
}

bool NGC1275::restoreState(const std::string& label) {
    auto it = ngc1275_saved_states.find(label);
    if (it == ngc1275_saved_states.end()) return false;
    
    for (const auto& pair : it->second) {
        setVariable(pair.first, pair.second);
    }
    return true;
}

std::vector<std::string> NGC1275::listSavedStates() const {
    std::vector<std::string> labels;
    for (const auto& pair : ngc1275_saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string NGC1275::exportState() const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(15);
    oss << "NGC1275 State Export\n";
    oss << "====================\n";
    
    auto vars = listVariables();
    for (const auto& var : vars) {
        oss << var << ": " << getVariable(var) << "\n";
    }
    return oss.str();
}

// System Analysis
std::map<std::string, double> NGC1275::sensitivityAnalysis(double t, double perturbation) {
    std::map<std::string, double> sensitivities;
    double baseline = compute_g_NGC1275(t);
    
    auto vars = listVariables();
    for (const auto& var : vars) {
        double original = getVariable(var);
        setVariable(var, original * (1.0 + perturbation));
        double perturbed = compute_g_NGC1275(t);
        sensitivities[var] = std::abs(perturbed - baseline) / baseline;
        setVariable(var, original);
    }
    return sensitivities;
}

std::string NGC1275::generateReport(double t) const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6);
    oss << "=== NGC 1275 System Report ===\n";
    oss << "Time: " << (t / (1e6 * 3.156e7)) << " Myr\n";
    oss << "B(t) magnetic field: " << B_t(t) << " T\n";
    oss << "F(t) filament factor: " << F_t(t) << "\n";
    oss << "g_NGC1275: " << compute_g_NGC1275(t) << " m/s^2\n";
    oss << "Core parameters: M=" << (M/1.989e30) << " M_sun, r=" << (r/9.461e15) << " ly\n";
    oss << "Black hole: M_BH=" << (M_BH/1.989e30) << " M_sun, r_BH=" << r_BH << " m\n";
    oss << "Magnetic: B0=" << B0 << " T, tau_B=" << (tau_B/(1e6*3.156e7)) << " Myr\n";
    oss << "Filament: F0=" << F0 << ", tau_fil=" << (tau_fil/(1e6*3.156e7)) << " Myr\n";
    oss << "Cooling: rho_cool=" << rho_cool << " kg/m^3, v_cool=" << v_cool << " m/s\n";
    return oss.str();
}

bool NGC1275::validateConsistency() const {
    bool valid = true;
    if (M <= 0 || r <= 0 || M_BH <= 0 || r_BH <= 0) valid = false;
    if (tau_B <= 0 || tau_fil <= 0) valid = false;
    if (B0 < 0 || F0 < 0 || F0 > 1) valid = false;
    if (rho_cool < 0 || v_cool < 0) valid = false;
    return valid;
}

bool NGC1275::autoCorrectAnomalies() {
    bool corrected = false;
    if (M <= 0) { M = 1e11 * 1.989e30; corrected = true; }
    if (r <= 0) { r = 200000 * 9.461e15; corrected = true; }
    if (M_BH <= 0) { M_BH = 8e8 * 1.989e30; corrected = true; }
    if (r_BH <= 0) { r_BH = 1e18; corrected = true; }
    if (tau_B <= 0) { tau_B = 100e6 * 3.156e7; corrected = true; }
    if (tau_fil <= 0) { tau_fil = 100e6 * 3.156e7; corrected = true; }
    if (B0 < 0) { B0 = 5e-9; corrected = true; }
    if (F0 < 0) { F0 = 0.0; corrected = true; }
    if (F0 > 1) { F0 = 1.0; corrected = true; }
    if (rho_cool < 0) { rho_cool = 1e-20; corrected = true; }
    if (v_cool < 0) { v_cool = 3e3; corrected = true; }
    if (corrected) updateCache();
    return corrected;
}

// ========== ENHANCED EXAMPLE FUNCTION ==========
void enhanced_ngc1275_example() {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "=========================================================\n";
    std::cout << "ENHANCED NGC 1275 DEMONSTRATION\n";
    std::cout << "Perseus A - Magnetic Monster AGN\n";
    std::cout << "=========================================================\n\n";
    
    NGC1275 ngc;
    
    // Step 1: Initial state and validation
    std::cout << "Step 1: Initial State and Validation\n";
    std::cout << "System: " << ngc.getSystemName() << "\n";
    std::cout << "Validation: " << (ngc.validateConsistency() ? "PASS" : "FAIL") << "\n";
    std::cout << "Auto-corrected: " << (ngc.autoCorrectAnomalies() ? "Yes" : "No") << "\n\n";
    
    // Step 2: Time evolution showing B(t) and F(t)
    std::cout << "Step 2: Time Evolution (Magnetic Field B(t) and Filament Support F(t))\n";
    double t_Myr_array[] = {0.0, 25.0, 50.0, 100.0, 200.0};
    for (double t_Myr : t_Myr_array) {
        double t = t_Myr * 1e6 * 3.156e7;
        double Bt = ngc.B_t(t);
        double Ft = ngc.F_t(t);
        double g = ngc.compute_g_NGC1275(t);
        std::cout << "  t = " << t_Myr << " Myr: B(t) = " << Bt << " T, F(t) = " << Ft << ", g = " << g << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 3: Variable listing
    std::cout << "Step 3: Variable Listing\n";
    auto vars = ngc.listVariables();
    std::cout << "Total variables: " << vars.size() << "\n";
    std::cout << "Sample: " << vars[0] << ", " << vars[1] << ", " << vars[12] << " (M_BH), " 
              << vars[14] << " (F0)\n\n";
    
    // Step 4: Galaxy mass scaling
    std::cout << "Step 4: Galaxy Mass Scaling (M sweeps)\n";
    ngc.saveState("original");
    double M_factors[] = {0.5, 1.0, 2.0};
    for (double factor : M_factors) {
        ngc.restoreState("original");
        ngc.expandGalaxyScale(factor, 1.0);
        double t = 50e6 * 3.156e7;
        double g = ngc.compute_g_NGC1275(t);
        double M_sun = 1.989e30;
        double M = ngc.getVariable("M");
        std::cout << "  M × " << factor << ": M = " << (M/M_sun) << " M_sun, g(50 Myr) = " << g << " m/s^2\n";
    }
    ngc.restoreState("original");
    std::cout << "\n";
    
    // Step 5: Black hole scaling (UNIQUE to AGN with SMBH)
    std::cout << "Step 5: Black Hole Scaling (M_BH sweeps) - AGN FEATURE\n";
    double M_BH_factors[] = {0.5, 1.0, 2.0};
    for (double factor : M_BH_factors) {
        ngc.restoreState("original");
        ngc.expandBlackHoleScale(factor, 1.0);
        double t = 50e6 * 3.156e7;
        double g = ngc.compute_g_NGC1275(t);
        double M_sun = 1.989e30;
        double M_BH = ngc.getVariable("M_BH");
        std::cout << "  M_BH × " << factor << ": M_BH = " << (M_BH/M_sun) << " M_sun, g = " << g << " m/s^2\n";
    }
    ngc.restoreState("original");
    std::cout << "\n";
    
    // Step 6: Magnetic field scaling (UNIQUE to magnetic monster)
    std::cout << "Step 6: Magnetic Field Scaling (B0 sweeps) - MAGNETIC MONSTER FEATURE\n";
    double B0_factors[] = {0.5, 1.0, 2.0};
    for (double factor : B0_factors) {
        ngc.restoreState("original");
        ngc.expandMagneticFilamentScale(factor, 1.0);
        double t = 50e6 * 3.156e7;
        double Bt = ngc.B_t(t);
        double g = ngc.compute_g_NGC1275(t);
        std::cout << "  B0 × " << factor << ": B(50 Myr) = " << Bt << " T, g = " << g << " m/s^2\n";
    }
    ngc.restoreState("original");
    std::cout << "\n";
    
    // Step 7: Filament support scaling
    std::cout << "Step 7: Filament Support Scaling (F0 sweeps)\n";
    double F0_factors[] = {0.5, 1.0, 2.0};
    for (double factor : F0_factors) {
        ngc.restoreState("original");
        ngc.expandMagneticFilamentScale(1.0, factor);
        double t = 50e6 * 3.156e7;
        double Ft = ngc.F_t(t);
        double g = ngc.compute_g_NGC1275(t);
        std::cout << "  F0 × " << factor << ": F(50 Myr) = " << Ft << ", g = " << g << " m/s^2\n";
    }
    ngc.restoreState("original");
    std::cout << "\n";
    
    // Step 8: Cooling flow velocity scaling
    std::cout << "Step 8: Cooling Flow Velocity Scaling\n";
    double v_cool_factors[] = {0.5, 1.0, 2.0};
    for (double factor : v_cool_factors) {
        ngc.restoreState("original");
        ngc.setVariable("v_cool", ngc.getVariable("v_cool") * factor);
        double t = 50e6 * 3.156e7;
        double g = ngc.compute_g_NGC1275(t);
        std::cout << "  v_cool × " << factor << ": g(50 Myr) = " << g << " m/s^2\n";
    }
    ngc.restoreState("original");
    std::cout << "\n";
    
    // Step 9: Parameter space expansion
    std::cout << "Step 9: Parameter Space Expansion (all scalable params)\n";
    ngc.expandParameterSpace(1.2);
    double M_after = ngc.getVariable("M");
    double M_sun = 1.989e30;
    std::cout << "  After 1.2× expansion: M = " << (M_after/M_sun) << " M_sun\n";
    ngc.restoreState("original");
    std::cout << "\n";
    
    // Step 10: Batch operations
    std::cout << "Step 10: Batch Operations (scale multiple variables)\n";
    std::vector<std::string> scale_group = {"M", "M_BH", "B0", "F0"};
    ngc.scaleVariableGroup(scale_group, 1.1);
    std::cout << "  Scaled {M, M_BH, B0, F0} by 1.1×\n";
    ngc.restoreState("original");
    std::cout << "\n";
    
    // Step 11: State management
    std::cout << "Step 11: State Management\n";
    ngc.saveState("state_A");
    ngc.expandBlackHoleScale(1.5, 1.2);
    ngc.saveState("state_B");
    auto states = ngc.listSavedStates();
    std::cout << "  Saved states: ";
    for (const auto& s : states) std::cout << s << " ";
    std::cout << "\n";
    ngc.restoreState("state_A");
    std::cout << "  Restored state_A\n\n";
    
    // Step 12: Generate parameter variations
    std::cout << "Step 12: Generate Parameter Variations (5% variation)\n";
    auto variations = ngc.generateVariations(3, 5.0);
    std::cout << "  Generated " << variations.size() << " variants with 5% random variation\n";
    std::cout << "  Variant 1 M = " << variations[0]["M"] << " kg\n\n";
    
    // Step 13: Sensitivity analysis
    std::cout << "Step 13: Sensitivity Analysis at 50 Myr\n";
    double t_sens = 50e6 * 3.156e7;
    auto sensitivities = ngc.sensitivityAnalysis(t_sens, 0.01);
    std::cout << "  Top sensitivities (1% perturbation):\n";
    std::vector<std::pair<std::string, double>> sens_vec(sensitivities.begin(), sensitivities.end());
    std::sort(sens_vec.begin(), sens_vec.end(), 
              [](const auto& a, const auto& b) { return a.second > b.second; });
    for (int i = 0; i < 5 && i < (int)sens_vec.size(); ++i) {
        std::cout << "    " << sens_vec[i].first << ": " << sens_vec[i].second << "\n";
    }
    std::cout << "\n";
    
    // Step 14: Auto-refinement with synthetic observations
    std::cout << "Step 14: Auto-Refinement (synthetic observations)\n";
    std::vector<std::pair<double, double>> obs;
    for (int i = 0; i <= 5; ++i) {
        double t_obs = i * 50e6 * 3.156e7;
        double g_obs = ngc.compute_g_NGC1275(t_obs) * (1.0 + 0.01 * (rand() % 100 - 50) / 100.0);
        obs.push_back({t_obs, g_obs});
    }
    ngc.autoRefineParameters(obs);
    std::cout << "  Refined parameters based on " << obs.size() << " observations\n\n";
    
    // Step 15: Optimization for maximum acceleration
    std::cout << "Step 15: Optimize for Maximum Acceleration\n";
    ngc.restoreState("original");
    auto metric = [](double g) { return g; };
    double t_opt_start = 0.0;
    double t_opt_end = 200e6 * 3.156e7;
    double best_g = ngc.optimizeForMetric(metric, t_opt_start, t_opt_end, 50);
    std::cout << "  Best g over 200 Myr: " << best_g << " m/s^2\n\n";
    
    // Step 16: Evolutionary system adaptation
    std::cout << "Step 16: Evolutionary System Adaptation (5 generations)\n";
    ngc.restoreState("original");
    auto fitness = [](const NGC1275& n) {
        double t = 50e6 * 3.156e7;
        return n.compute_g_NGC1275(t);
    };
    ngc.evolveSystem(5, fitness);
    std::cout << "  Evolved system over 5 generations (fitness = g at 50 Myr)\n\n";
    
    // Step 17: Full system report
    std::cout << "Step 17: Full System Report at 50 Myr\n";
    ngc.restoreState("original");
    double t_report = 50e6 * 3.156e7;
    std::string report = ngc.generateReport(t_report);
    std::cout << report << "\n";
    
    // Step 18: Full state export
    std::cout << "Step 18: Full State Export\n";
    std::string exported = ngc.exportState();
    std::cout << "Exported state (first 500 chars):\n";
    std::cout << exported.substr(0, 500) << "...\n\n";
    
    std::cout << "=========================================================\n";
    std::cout << "ENHANCED DEMONSTRATION COMPLETE\n";
    std::cout << "=========================================================\n";
}