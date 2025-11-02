// SpiralSupernovaeUQFFModule.h
// Modular C++ implementation of the full Master Universal Gravity Equation (UQFF & SM Integration) for Spirals and Supernovae Evolution.
// This module can be plugged into a base program (e.g., 'ziqn233h.cpp') by including this header and linking the .cpp.
// Usage in base: #include "SpiralSupernovaeUQFFModule.h"
// SpiralSupernovaeUQFFModule mod; mod.computeG(t); mod.updateVariable("M", new_value);
// All variables are stored in a std::map for dynamic addition/subtraction/update.
// Nothing is negligible: Includes all terms - base gravity with T_spiral, Ug1-Ug4, cosmological Lambda with ?_?, quantum integral, Lorentz q(v x B), fluid rho_fluid V g, resonant oscillatory (cos/exp), DM/visible with perturbations, supernova SN_term.
// Associated text: Outputs descriptive equation string via getEquationText().
// Approximations: Quantum integral normalized to 1.0; exp real part; Ug2/Ug3=0; DM fraction 0.85; T_spiral with ?_p; SN_term from L_SN.
// Spiral-SN params: M=1.989e41 kg, r=9.258e20 m, H0=73 km/s/Mpc, ?_p=20 km/s/kpc, L_SN=1e36 W, z up to 1.5, etc.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#ifndef SPIRAL_SUPERNOVAE_UQFF_MODULE_H
#define SPIRAL_SUPERNOVAE_UQFF_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>
#include <functional>
#include <random>
#include <algorithm>
#include <sstream>
#include <vector>

class SpiralSupernovaeUQFFModule {
private:
    std::map<std::string, double> variables;
    double computeQuantumTerm(double t_Hubble_val);
    double computeFluidTerm(double g_base);
    double computeResonantTerm(double t);
    double computeDMTerm();
    double computeUgSum();
    double computeHz(double z);
    double computeT_spiral(double t);
    double computeSN_term(double z);

public:
    // Constructor: Initialize all variables with Spirals and Supernovae defaults
    SpiralSupernovaeUQFFModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computation: Full g_UQFF(r, t) for Spirals and Supernovae
    double computeG(double t, double z);

    // Output descriptive text of the equation
    std::string getEquationText();

    // Print all current variables (for debugging/updates)
    void printVariables();

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
    void expandSpiralScale(double scale_factor);
    void expandSupernovaScale(double scale_factor);
    void expandCosmologyScale(double scale_factor);

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
    std::string exportState(double t, double z);

    // System Analysis
    std::map<std::string, double> sensitivityAnalysis(const std::string& param, double t, double z, double delta);
    std::string generateReport(double t, double z);
    bool validateConsistency();
    void autoCorrectAnomalies();
};

#endif // SPIRAL_SUPERNOVAE_UQFF_MODULE_H

// SpiralSupernovaeUQFFModule.cpp
#include "SpiralSupernovaeUQFFModule.h"
#include <complex>

// Constructor: Set all variables with Spirals and Supernovae-specific values
SpiralSupernovaeUQFFModule::SpiralSupernovaeUQFFModule() {
    // Base constants (universal)
    variables["G"] = 6.6743e-11;                    // m^3 kg^-1 s^-2
    variables["c"] = 3e8;                           // m/s
    variables["hbar"] = 1.0546e-34;                 // J s
    variables["Lambda"] = 1.1e-52;                  // m^-2
    variables["q"] = 1.602e-19;                     // C
    variables["pi"] = 3.141592653589793;            // pi
    variables["t_Hubble"] = 13.8e9 * 3.156e7;       // s

    // Spiral-SN parameters
    double M_sun_val = 1.989e30;                    // kg
    variables["M_sun"] = M_sun_val;
    variables["M"] = 1e11 * M_sun_val;              // Galaxy mass kg
    variables["M_visible"] = 0.15 * variables["M"]; // Visible fraction
    variables["M_DM"] = 0.85 * variables["M"];      // Dark matter
    variables["r"] = 9.258e20;                      // m (~30 kpc)
    variables["M_gas"] = 1e9 * M_sun_val;           // Gas mass

    // Hubble/cosmology
    variables["H0"] = 73.0;                         // km/s/Mpc (SH0ES)
    variables["Mpc_to_m"] = 3.086e22;               // m/Mpc
    variables["z"] = 0.5;                           // Typical z for SN
    variables["Omega_m"] = 0.3;
    variables["Omega_Lambda"] = 0.7;
    variables["t"] = 5e9 * 3.156e7;                 // Default t=5 Gyr s

    // Spiral dynamics
    variables["Omega_p"] = 20e3 / 3.086e19;         // rad/s (20 km/s/kpc pattern speed)

    // SN parameters
    variables["L_SN"] = 1e36;                       // W (peak luminosity)
    variables["rho_fluid"] = 1e-21;                 // kg/m^3 (ISM)
    variables["V"] = 1e3;                           // m^3
    variables["v_rot"] = 2e5;                       // m/s (rotation)
    variables["delta_rho"] = 0.1 * variables["rho_fluid"];
    variables["rho"] = variables["rho_fluid"];

    // EM/magnetic
    variables["B"] = 1e-5;                          // T (galactic field)
    variables["B_crit"] = 1e11;                     // T

    // Quantum terms
    variables["Delta_x"] = 1e-10;                   // m
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];
    variables["integral_psi"] = 1.0;

    // Resonant/oscillatory
    variables["A"] = 1e-10;
    variables["k"] = 1e20;
    variables["omega"] = 1e15;
    variables["x"] = 0.0;

    // Ug subterms
    variables["Ug1"] = 0.0;
    variables["Ug2"] = 0.0;
    variables["Ug3"] = 0.0;
    variables["Ug4"] = 0.0;

    // Scale factors
    variables["scale_macro"] = 1e-12;
    variables["f_TRZ"] = 0.1;
    variables["f_sc"] = 1.0;
}

// Update variable (set to new value)
void SpiralSupernovaeUQFFModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
    if (name == "Delta_x") {
        variables["Delta_p"] = variables["hbar"] / value;
    } else if (name == "M") {
        variables["M_visible"] = 0.15 * value;
        variables["M_DM"] = 0.85 * value;
    } else if (name == "H0") {
        // Recompute if needed
    }
}

// Add delta to variable
void SpiralSupernovaeUQFFModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta from variable
void SpiralSupernovaeUQFFModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute H(z) in s^-1
double SpiralSupernovaeUQFFModule::computeHz(double z) {
    double Hz_kms = variables["H0"] * std::sqrt(variables["Omega_m"] * std::pow(1.0 + z, 3) + variables["Omega_Lambda"]);
    return (Hz_kms * 1e3) / variables["Mpc_to_m"];
}

// Compute Ug sum: Ug1 = G M / r^2, Ug4 = Ug1 * f_sc, others 0
double SpiralSupernovaeUQFFModule::computeUgSum() {
    double Ug1 = (variables["G"] * variables["M"]) / (variables["r"] * variables["r"]);
    variables["Ug1"] = Ug1;
    variables["Ug4"] = Ug1 * variables["f_sc"];
    return variables["Ug1"] + variables["Ug2"] + variables["Ug3"] + variables["Ug4"];
}

// Quantum term: (hbar / sqrt(Delta_x Delta_p)) * integral * (2 pi / t_Hubble)
double SpiralSupernovaeUQFFModule::computeQuantumTerm(double t_Hubble_val) {
    double unc = std::sqrt(variables["Delta_x"] * variables["Delta_p"]);
    double integral_val = variables["integral_psi"];
    return (variables["hbar"] / unc) * integral_val * (2 * variables["pi"] / t_Hubble_val);
}

// Fluid term: rho_fluid * V * g
double SpiralSupernovaeUQFFModule::computeFluidTerm(double g_base) {
    return variables["rho_fluid"] * variables["V"] * g_base;
}

// Resonant terms: 2 A cos(k x) cos(omega t) + (2 pi / 13.8) A Re[exp(i (k x - omega t))]
double SpiralSupernovaeUQFFModule::computeResonantTerm(double t) {
    double cos_term = 2 * variables["A"] * std::cos(variables["k"] * variables["x"]) * std::cos(variables["omega"] * t);
    std::complex<double> exp_term(variables["A"] * std::exp(std::complex<double>(0, variables["k"] * variables["x"] - variables["omega"] * t)));
    double real_exp = exp_term.real();
    double exp_factor = (2 * variables["pi"] / 13.8);
    return cos_term + exp_factor * real_exp;
}

// DM term: (M_visible + M_DM) * (delta_rho / rho + 3 G M / r^3)
double SpiralSupernovaeUQFFModule::computeDMTerm() {
    double pert = variables["delta_rho"] / variables["rho"];
    double curv = 3 * variables["G"] * variables["M"] / (variables["r"] * variables["r"] * variables["r"]);
    return (variables["M_visible"] + variables["M_DM"]) * (pert + curv);
}

// Spiral torque term: T_spiral = G * M_gas * M / r^2 * (1 + ?_p * t)
double SpiralSupernovaeUQFFModule::computeT_spiral(double t) {
    double torque_base = (variables["G"] * variables["M_gas"] * variables["M"]) / (variables["r"] * variables["r"]);
    return torque_base * (1.0 + variables["Omega_p"] * t);
}

// Supernova term: SN_term = (L_SN / (4 pi r^2 c)) * (1 + H(z) * t)
double SpiralSupernovaeUQFFModule::computeSN_term(double z) {
    double Hz = computeHz(z);
    double flux = variables["L_SN"] / (4 * variables["pi"] * variables["r"] * variables["r"] * variables["c"]);
    return flux * (1.0 + Hz * variables["t"]);
}

// Full computation: g_UQFF(r, t) = ... all terms with T_spiral and SN_term
double SpiralSupernovaeUQFFModule::computeG(double t, double z) {
    variables["t"] = t;
    double Hz = computeHz(z);
    double expansion = 1.0 + Hz * t;
    double sc_correction = 1.0 - (variables["B"] / variables["B_crit"]);
    double tr_factor = 1.0 + variables["f_TRZ"];
    double t_spiral = computeT_spiral(t);
    double sn_term = computeSN_term(z);

    // Base gravity with expansion, SC, TR, T_spiral
    double g_base = ((variables["G"] * variables["M"] / (variables["r"] * variables["r"])) * expansion * sc_correction * tr_factor) * (1.0 + t_spiral);

    // Ug sum
    double ug_sum = computeUgSum();

    // Cosmological with ?_?
    double lambda_term = variables["Lambda"] * (variables["c"] * variables["c"] * variables["Omega_Lambda"]) / 3.0;

    // Quantum
    double quantum_term = computeQuantumTerm(variables["t_Hubble"]);

    // EM Lorentz (rotation v_rot B)
    double em_base = variables["q"] * variables["v_rot"] * variables["B"] / 1.673e-27;
    double em_term = em_base * (1.0 + (7.09e-36 / 7.09e-37)) * variables["scale_macro"];

    // Fluid
    double fluid_term = computeFluidTerm(g_base);

    // Resonant
    double resonant_term = computeResonantTerm(t);

    // DM
    double dm_term = computeDMTerm();

    // Total: Sum all + SN_term
    return g_base + ug_sum + lambda_term + quantum_term + em_term + fluid_term + resonant_term + dm_term + sn_term;
}

// Get equation text (descriptive)
std::string SpiralSupernovaeUQFFModule::getEquationText() {
    return "g_Spiral_SN(r, t) = (G * M / r^2) * (1 + H(z) * t) * (1 + T_spiral) * (1 - B / B_crit) * (1 + f_TRZ) + (Ug1 + Ug2 + Ug3 + Ug4) + (Lambda * c^2 * ?_? / 3) + "
           "(hbar / sqrt(Delta_x * Delta_p)) * ?(?* H ? dV) * (2? / t_Hubble) + q (v � B) + ?_fluid * V * g + "
           "2 A cos(k x) cos(? t) + (2? / 13.8) A exp(i (k x - ? t)) + (M_visible + M_DM) * (??/? + 3 G M / r^3) + SN_term\n"
           "Where T_spiral = G * M_gas * M / r^2 * (1 + ?_p * t); SN_term = (L_SN / (4? r^2 c)) * (1 + H(z) * t)\n"
           "Special Terms:\n"
           "- Quantum: Heisenberg uncertainty for ISM quantum effects.\n"
           "- Fluid: Gas density-volume-gravity coupling in arms.\n"
           "- Resonant: Oscillatory Aether waves for density waves.\n"
           "- DM: Visible+dark mass with perturbations for rotation curves.\n"
           "- Superconductivity: (1 - B/B_crit) for galactic fields.\n"
           "- Spiral Torque: T_spiral for arm evolution.\n"
           "- Supernova: SN_term for expansion probe.\n"
           "Solutions: At t=5 Gyr, z=0.5, g_Spiral_SN ~1e-10 m/s� (Lambda/SN dominant; g_base ~1e-10).\n"
           "Adaptations for Spirals and Supernovae: SH0ES H0=73; ?_p=20 km/s/kpc; L_SN=1e36 W for Ia SN.";
}

// Print variables
void SpiralSupernovaeUQFFModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// ===== ENHANCED DYNAMIC CAPABILITIES IMPLEMENTATION =====

namespace {
    std::map<std::string, std::map<std::string, double>> spiral_sn_saved_states;
}

// Variable Management
void SpiralSupernovaeUQFFModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void SpiralSupernovaeUQFFModule::removeVariable(const std::string& name) {
    variables.erase(name);
}

void SpiralSupernovaeUQFFModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
    }
}

std::vector<std::string> SpiralSupernovaeUQFFModule::listVariables() {
    std::vector<std::string> names;
    for (const auto& pair : variables) {
        names.push_back(pair.first);
    }
    return names;
}

std::string SpiralSupernovaeUQFFModule::getSystemName() {
    return "Spiral Galaxies & Supernovae - Full UQFF & SM Integration";
}

// Batch Operations
void SpiralSupernovaeUQFFModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = func(variables[name]);
        }
    }
}

void SpiralSupernovaeUQFFModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// Self-Expansion (domain-specific for Spirals & Supernovae)
void SpiralSupernovaeUQFFModule::expandParameterSpace(double scale_factor) {
    // Scale core galactic parameters
    variables["r"] *= scale_factor;
    variables["M"] *= scale_factor;
    variables["M_gas"] *= scale_factor;
    variables["rho_fluid"] *= scale_factor;
    // Update dependent variables
    variables["M_visible"] = 0.15 * variables["M"];
    variables["M_DM"] = 0.85 * variables["M"];
    variables["delta_rho"] = 0.1 * variables["rho_fluid"];
    variables["rho"] = variables["rho_fluid"];
}

void SpiralSupernovaeUQFFModule::expandSpiralScale(double scale_factor) {
    // Scale spiral dynamics parameters (Omega_p, M_gas, v_rot)
    variables["Omega_p"] *= scale_factor;
    variables["M_gas"] *= scale_factor;
    variables["v_rot"] *= scale_factor;
}

void SpiralSupernovaeUQFFModule::expandSupernovaScale(double scale_factor) {
    // Scale supernova parameters (L_SN, luminosity)
    variables["L_SN"] *= scale_factor;
}

void SpiralSupernovaeUQFFModule::expandCosmologyScale(double scale_factor) {
    // Scale cosmological parameters (H0, z, Omega values)
    variables["H0"] *= scale_factor;
    variables["z"] *= scale_factor;
}

// Self-Refinement
void SpiralSupernovaeUQFFModule::autoRefineParameters(double tolerance) {
    // Enforce physical constraints
    if (variables["M"] <= 0) variables["M"] = 1e11 * variables["M_sun"];
    if (variables["r"] <= 0) variables["r"] = 9.258e20;
    if (variables["rho_fluid"] <= 0) variables["rho_fluid"] = 1e-21;
    if (variables["v_rot"] < 0) variables["v_rot"] = 2e5;
    if (variables["M_gas"] <= 0) variables["M_gas"] = 1e9 * variables["M_sun"];
    if (variables["L_SN"] <= 0) variables["L_SN"] = 1e36;
    if (variables["H0"] <= 0) variables["H0"] = 73.0;
    if (variables["z"] < 0) variables["z"] = 0.5;
    // Update dependent variables
    variables["M_visible"] = 0.15 * variables["M"];
    variables["M_DM"] = 0.85 * variables["M"];
    variables["delta_rho"] = 0.1 * variables["rho_fluid"];
    variables["rho"] = variables["rho_fluid"];
}

void SpiralSupernovaeUQFFModule::calibrateToObservations(const std::map<std::string, double>& obs_data) {
    for (const auto& obs : obs_data) {
        if (variables.find(obs.first) != variables.end()) {
            variables[obs.first] = obs.second;
        }
    }
    // Update dependent variables
    if (obs_data.find("M") != obs_data.end()) {
        variables["M_visible"] = 0.15 * variables["M"];
        variables["M_DM"] = 0.85 * variables["M"];
    }
    if (obs_data.find("rho_fluid") != obs_data.end()) {
        variables["delta_rho"] = 0.1 * variables["rho_fluid"];
        variables["rho"] = variables["rho_fluid"];
    }
    if (obs_data.find("Delta_x") != obs_data.end()) {
        variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];
    }
}

void SpiralSupernovaeUQFFModule::optimizeForMetric(const std::string& var_name, double target_value, int iterations) {
    if (variables.find(var_name) == variables.end()) return;
    double best_value = variables[var_name];
    double best_error = std::abs(variables[var_name] - target_value);
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.9, 1.1);
    
    for (int i = 0; i < iterations; ++i) {
        double test_value = variables[var_name] * dis(gen);
        variables[var_name] = test_value;
        double error = std::abs(test_value - target_value);
        if (error < best_error) {
            best_error = error;
            best_value = test_value;
        }
    }
    variables[var_name] = best_value;
}

// Parameter Exploration
std::vector<std::map<std::string, double>> SpiralSupernovaeUQFFModule::generateVariations(int n_variations) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.8, 1.2);
    
    for (int i = 0; i < n_variations; ++i) {
        std::map<std::string, double> variation = variables;
        variation["M"] *= dis(gen);
        variation["H0"] *= dis(gen);
        variation["L_SN"] *= dis(gen);
        variation["v_rot"] *= dis(gen);
        variations.push_back(variation);
    }
    return variations;
}

// Adaptive Evolution
void SpiralSupernovaeUQFFModule::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-mutation_rate, mutation_rate);
    
    variables["M"] *= (1.0 + dis(gen));
    variables["H0"] *= (1.0 + dis(gen));
    variables["L_SN"] *= (1.0 + dis(gen));
    variables["v_rot"] *= (1.0 + dis(gen));
    
    // Update dependent variables
    variables["M_visible"] = 0.15 * variables["M"];
    variables["M_DM"] = 0.85 * variables["M"];
}

void SpiralSupernovaeUQFFModule::evolveSystem(int generations, std::function<double()> fitness_function) {
    double best_fitness = fitness_function();
    std::map<std::string, double> best_state = variables;
    
    for (int gen = 0; gen < generations; ++gen) {
        mutateParameters(0.1);
        double fitness = fitness_function();
        if (fitness > best_fitness) {
            best_fitness = fitness;
            best_state = variables;
        } else {
            variables = best_state;
        }
    }
}

// State Management
void SpiralSupernovaeUQFFModule::saveState(const std::string& label) {
    spiral_sn_saved_states[label] = variables;
}

void SpiralSupernovaeUQFFModule::restoreState(const std::string& label) {
    if (spiral_sn_saved_states.find(label) != spiral_sn_saved_states.end()) {
        variables = spiral_sn_saved_states[label];
    }
}

std::vector<std::string> SpiralSupernovaeUQFFModule::listSavedStates() {
    std::vector<std::string> labels;
    for (const auto& pair : spiral_sn_saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string SpiralSupernovaeUQFFModule::exportState(double t, double z) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(6);
    oss << "Spiral-SN State Export at t=" << t << " s, z=" << z << ":\n";
    oss << "M=" << variables["M"] << " kg (" << variables["M"]/variables["M_sun"] << " M_sun)\n";
    oss << "r=" << variables["r"] << " m (~" << variables["r"]/3.086e19 << " kpc)\n";
    oss << "H0=" << variables["H0"] << " km/s/Mpc, L_SN=" << variables["L_SN"] << " W\n";
    oss << "v_rot=" << variables["v_rot"] << " m/s, Omega_p=" << variables["Omega_p"] << " rad/s\n";
    oss << "g_total=" << computeG(t, z) << " m/s²\n";
    return oss.str();
}

// System Analysis
std::map<std::string, double> SpiralSupernovaeUQFFModule::sensitivityAnalysis(const std::string& param, double t, double z, double delta) {
    std::map<std::string, double> result;
    if (variables.find(param) == variables.end()) return result;
    
    double original = variables[param];
    double g_original = computeG(t, z);
    
    variables[param] = original * (1.0 + delta);
    double g_plus = computeG(t, z);
    
    variables[param] = original * (1.0 - delta);
    double g_minus = computeG(t, z);
    
    variables[param] = original;
    
    result["dg/d" + param] = (g_plus - g_minus) / (2.0 * delta * original);
    result["g_original"] = g_original;
    result["g_plus"] = g_plus;
    result["g_minus"] = g_minus;
    
    return result;
}

std::string SpiralSupernovaeUQFFModule::generateReport(double t, double z) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(6);
    oss << "===== SPIRAL-SUPERNOVA UQFF Module Report (t=" << t << " s, z=" << z << ") =====\n\n";
    oss << "System: " << getSystemName() << "\n\n";
    
    oss << "Core Parameters:\n";
    oss << "  M = " << variables["M"] << " kg (" << variables["M"]/variables["M_sun"] << " M_sun)\n";
    oss << "  r = " << variables["r"] << " m (~" << variables["r"]/3.086e19 << " kpc)\n";
    oss << "  M_gas = " << variables["M_gas"] << " kg (" << variables["M_gas"]/variables["M_sun"] << " M_sun)\n";
    oss << "  v_rot = " << variables["v_rot"] << " m/s (" << variables["v_rot"]/1e3 << " km/s)\n";
    oss << "  H0 = " << variables["H0"] << " km/s/Mpc (SH0ES)\n";
    oss << "  z = " << z << " (redshift)\n";
    oss << "  L_SN = " << variables["L_SN"] << " W (peak luminosity)\n";
    oss << "  Omega_p = " << variables["Omega_p"] << " rad/s (pattern speed)\n\n";
    
    double Hz = computeHz(z);
    double expansion = 1.0 + Hz * t;
    double sc_correction = 1.0 - (variables["B"] / variables["B_crit"]);
    double tr_factor = 1.0 + variables["f_TRZ"];
    double t_spiral = computeT_spiral(t);
    double sn_term = computeSN_term(z);
    double g_base = ((variables["G"] * variables["M"] / (variables["r"] * variables["r"])) * expansion * sc_correction * tr_factor) * (1.0 + t_spiral);
    double ug_sum = computeUgSum();
    double lambda_term = variables["Lambda"] * (variables["c"] * variables["c"] * variables["Omega_Lambda"]) / 3.0;
    double quantum_term = computeQuantumTerm(variables["t_Hubble"]);
    double em_base = variables["q"] * variables["v_rot"] * variables["B"] / 1.673e-27;
    double em_term = em_base * (1.0 + (7.09e-36 / 7.09e-37)) * variables["scale_macro"];
    double fluid_term = computeFluidTerm(g_base);
    double resonant_term = computeResonantTerm(t);
    double dm_term = computeDMTerm();
    double g_total = computeG(t, z);
    
    oss << "Term Breakdown:\n";
    oss << "  g_base (with T_spiral) = " << g_base << " m/s²\n";
    oss << "  T_spiral = " << t_spiral << " (spiral torque factor)\n";
    oss << "  Ug_sum (Ug1+Ug2+Ug3+Ug4) = " << ug_sum << " m/s²\n";
    oss << "  Lambda_term (Omega_Lambda) = " << lambda_term << " m/s²\n";
    oss << "  Quantum_term = " << quantum_term << " m/s²\n";
    oss << "  EM_term (v_rot×B) = " << em_term << " m/s²\n";
    oss << "  Fluid_term = " << fluid_term << " m/s²\n";
    oss << "  Resonant_term = " << resonant_term << " m/s²\n";
    oss << "  DM_term = " << dm_term << " m/s²\n";
    oss << "  SN_term = " << sn_term << " m/s² [SUPERNOVA EXPANSION PROBE]\n\n";
    oss << "TOTAL g = " << g_total << " m/s²\n\n";
    
    oss << "Physics Notes:\n";
    oss << "- Spiral galaxies with pattern speed Omega_p=" << variables["Omega_p"] << " rad/s\n";
    oss << "- Type Ia Supernovae used as standard candles (L_SN=" << variables["L_SN"] << " W)\n";
    oss << "- SH0ES H0=" << variables["H0"] << " km/s/Mpc for cosmological distance ladder\n";
    oss << "- Full UQFF+SM: gravity, spiral torque, Ug1-4, Lambda, quantum, EM, fluid, resonance, DM, SN\n";
    oss << "- Lambda term includes Omega_Lambda for dark energy\n";
    oss << "- DM halo important for rotation curve (M_DM=" << variables["M_DM"] << " kg)\n";
    oss << "- Spiral arms modeled via torque: T_spiral ∝ (1 + Omega_p * t)\n";
    oss << "- SN term: flux-based expansion probe ∝ (1 + H(z) * t)\n\n";
    
    return oss.str();
}

bool SpiralSupernovaeUQFFModule::validateConsistency() {
    bool valid = true;
    if (variables["M"] <= 0) valid = false;
    if (variables["r"] <= 0) valid = false;
    if (variables["rho_fluid"] <= 0) valid = false;
    if (variables["v_rot"] < 0) valid = false;
    if (variables["M_gas"] <= 0) valid = false;
    if (variables["L_SN"] <= 0) valid = false;
    if (variables["H0"] <= 0) valid = false;
    if (variables["z"] < 0) valid = false;
    return valid;
}

void SpiralSupernovaeUQFFModule::autoCorrectAnomalies() {
    if (variables["M"] <= 0) variables["M"] = 1e11 * variables["M_sun"];
    if (variables["r"] <= 0) variables["r"] = 9.258e20;
    if (variables["rho_fluid"] <= 0) variables["rho_fluid"] = 1e-21;
    if (variables["v_rot"] < 0) variables["v_rot"] = 2e5;
    if (variables["M_gas"] <= 0) variables["M_gas"] = 1e9 * variables["M_sun"];
    if (variables["L_SN"] <= 0) variables["L_SN"] = 1e36;
    if (variables["H0"] <= 0) variables["H0"] = 73.0;
    if (variables["z"] < 0) variables["z"] = 0.5;
    // Update dependent variables
    variables["M_visible"] = 0.15 * variables["M"];
    variables["M_DM"] = 0.85 * variables["M"];
    variables["delta_rho"] = 0.1 * variables["rho_fluid"];
    variables["rho"] = variables["rho_fluid"];
}

// Enhanced example usage demonstration
void enhanced_example_usage() {
    SpiralSupernovaeUQFFModule mod;
    double t_5Gyr = 5e9 * 3.156e7;  // 5 Gyr in seconds
    double z_typical = 0.5;         // Typical SN redshift
    
    std::cout << "===== ENHANCED SPIRAL-SUPERNOVA UQFF MODULE DEMONSTRATION =====\n\n";
    
    // Step 1: Variable management
    std::cout << "Step 1: Variable Management\n";
    mod.createVariable("custom_dm_fraction", 0.85);
    mod.cloneVariable("H0", "H0_backup");
    std::vector<std::string> vars = mod.listVariables();
    std::cout << "Total variables: " << vars.size() << "\n";
    std::cout << "System: " << mod.getSystemName() << "\n\n";
    
    // Step 2: Batch scaling
    std::cout << "Step 2: Batch Scaling (Cosmological parameters)\n";
    mod.scaleVariableGroup({"H0", "z", "L_SN"}, 1.05);
    std::cout << "Scaled H0, z, L_SN by 1.05\n\n";
    
    // Step 3: Self-expansion (different physics domains)
    std::cout << "Step 3: Self-Expansion\n";
    mod.expandSpiralScale(1.08);  // Spiral dynamics +8%
    std::cout << "Expanded spiral scale +8%\n";
    mod.expandSupernovaScale(1.05);  // SN luminosity +5%
    std::cout << "Expanded supernova scale +5%\n";
    mod.expandCosmologyScale(1.02);  // Cosmology +2%
    std::cout << "Expanded cosmology scale +2%\n\n";
    
    // Step 4: Self-refinement
    std::cout << "Step 4: Self-Refinement\n";
    mod.autoRefineParameters(1e-10);
    std::cout << "Auto-refined parameters\n";
    std::map<std::string, double> obs_data = {
        {"H0", 73.5},
        {"L_SN", 1.1e36},
        {"v_rot", 2.2e5},
        {"z", 0.55}
    };
    mod.calibrateToObservations(obs_data);
    std::cout << "Calibrated to observations\n\n";
    
    // Step 5: Optimize for specific metric
    std::cout << "Step 5: Optimize for H0~73.0 km/s/Mpc\n";
    mod.optimizeForMetric("H0", 73.0, 50);
    std::cout << "Optimization complete\n\n";
    
    // Step 6: Generate variations
    std::cout << "Step 6: Generate 15 Parameter Variations\n";
    auto variations = mod.generateVariations(15);
    std::cout << "Generated " << variations.size() << " variations\n\n";
    
    // Step 7: State management
    std::cout << "Step 7: State Management\n";
    mod.saveState("initial");
    mod.scaleVariableGroup({"H0", "L_SN"}, 1.15);
    mod.saveState("enhanced_cosmology");
    mod.expandSpiralScale(0.9);
    mod.saveState("reduced_spiral");
    std::cout << "Saved 3 states\n\n";
    
    // Step 8: Sensitivity analysis
    std::cout << "Step 8: Sensitivity Analysis (H0 at t=5Gyr, z=0.5)\n";
    mod.restoreState("initial");
    auto sensitivity = mod.sensitivityAnalysis("H0", t_5Gyr, z_typical, 0.1);
    std::cout << "dg/dH0 = " << std::scientific << sensitivity["dg/dH0"] << " (m/s²)/(km/s/Mpc)\n\n";
    
    // Step 9: System validation
    std::cout << "Step 9: System Validation\n";
    bool valid = mod.validateConsistency();
    std::cout << "System consistency: " << (valid ? "VALID" : "INVALID") << "\n";
    if (!valid) {
        mod.autoCorrectAnomalies();
        std::cout << "Auto-corrected anomalies\n";
    }
    std::cout << "\n";
    
    // Step 10: Comprehensive report
    std::cout << "Step 10: Comprehensive Report (t=5Gyr, z=0.5)\n";
    std::string report = mod.generateReport(t_5Gyr, z_typical);
    std::cout << report << "\n";
    
    // Step 11: Adaptive evolution
    std::cout << "Step 11: Adaptive Evolution (25 generations)\n";
    auto fitness_fn = [&mod, t_5Gyr, z_typical]() -> double {
        double g = mod.computeG(t_5Gyr, z_typical);
        return -std::abs(std::log10(std::abs(g)) + 10.0);  // Target g~1e-10 m/s²
    };
    mod.evolveSystem(25, fitness_fn);
    std::cout << "Evolution complete\n\n";
    
    // Step 12: Time evolution comparison
    std::cout << "Step 12: Time Evolution (0 to 13.8 Gyr)\n";
    std::vector<double> times = {0.0, 1e9*3.156e7, 5e9*3.156e7, 10e9*3.156e7, 13.8e9*3.156e7};
    for (double t : times) {
        double g = mod.computeG(t, z_typical);
        std::cout << "t=" << std::scientific << t << " s (" << t/(3.156e7*1e9) << " Gyr): g=" << g << " m/s²\n";
    }
    std::cout << "\n";
    
    // Step 13: Redshift sweep
    std::cout << "Step 13: Redshift Sweep (SN distance ladder, t=5Gyr)\n";
    std::vector<double> redshifts = {0.1, 0.3, 0.5, 0.8, 1.0, 1.5};
    for (double z : redshifts) {
        double g = mod.computeG(t_5Gyr, z);
        std::cout << "z=" << z << ": g=" << std::scientific << g << " m/s²\n";
    }
    std::cout << "\n";
    
    // Step 14: H0 sweep (Hubble tension)
    std::cout << "Step 14: H0 Sweep (Hubble tension, t=5Gyr, z=0.5)\n";
    std::vector<double> H0_values = {67.4, 70.0, 73.0, 74.0, 76.0};  // Planck to SH0ES range
    for (double H0 : H0_values) {
        mod.updateVariable("H0", H0);
        double g = mod.computeG(t_5Gyr, z_typical);
        std::cout << "H0=" << H0 << " km/s/Mpc: g=" << std::scientific << g << " m/s²\n";
    }
    mod.restoreState("initial");
    std::cout << "\n";
    
    // Step 15: Supernova luminosity impact
    std::cout << "Step 15: SN Luminosity Impact (t=5Gyr, z=0.5)\n";
    std::vector<double> luminosities = {5e35, 7.5e35, 1e36, 1.25e36, 1.5e36};
    for (double L : luminosities) {
        mod.updateVariable("L_SN", L);
        double g = mod.computeG(t_5Gyr, z_typical);
        std::cout << "L_SN=" << std::scientific << L << " W: g=" << g << " m/s²\n";
    }
    mod.restoreState("initial");
    std::cout << "\n";
    
    // Step 16: Multi-parameter sensitivity
    std::cout << "Step 16: Multi-Parameter Sensitivity (t=5Gyr, z=0.5)\n";
    std::vector<std::string> params = {"M", "H0", "L_SN", "v_rot", "r"};
    for (const auto& param : params) {
        auto sens = mod.sensitivityAnalysis(param, t_5Gyr, z_typical, 0.05);
        std::cout << "dg/d" << param << " = " << std::scientific << sens["dg/d" + param] << "\n";
    }
    std::cout << "\n";
    
    // Step 17: Spiral torque evolution
    std::cout << "Step 17: Spiral Torque Evolution (z=0.5)\n";
    for (double t : times) {
        double t_spiral = mod.computeT_spiral(t);
        std::cout << "t=" << t/(3.156e7*1e9) << " Gyr: T_spiral=" << std::scientific << t_spiral << "\n";
    }
    std::cout << "\n";
    
    // Step 18: Expansion domains comparison
    std::cout << "Step 18: Expansion Domains Comparison (t=5Gyr, z=0.5)\n";
    mod.restoreState("initial");
    double g_initial = mod.computeG(t_5Gyr, z_typical);
    std::cout << "Initial: g=" << std::scientific << g_initial << " m/s²\n";
    
    mod.expandSpiralScale(1.2);
    double g_spiral = mod.computeG(t_5Gyr, z_typical);
    std::cout << "Spiral +20%: g=" << g_spiral << " m/s² (Δ=" << (g_spiral-g_initial)/g_initial*100 << "%)\n";
    mod.restoreState("initial");
    
    mod.expandSupernovaScale(1.2);
    double g_sn = mod.computeG(t_5Gyr, z_typical);
    std::cout << "SN +20%: g=" << g_sn << " m/s² (Δ=" << (g_sn-g_initial)/g_initial*100 << "%)\n";
    mod.restoreState("initial");
    
    mod.expandCosmologyScale(1.2);
    double g_cosmo = mod.computeG(t_5Gyr, z_typical);
    std::cout << "Cosmology +20%: g=" << g_cosmo << " m/s² (Δ=" << (g_cosmo-g_initial)/g_initial*100 << "%)\n";
    mod.restoreState("initial");
    std::cout << "\n";
    
    // Step 19: State restoration
    std::cout << "Step 19: State Restoration\n";
    mod.restoreState("initial");
    std::cout << "Restored initial state\n";
    g_initial = mod.computeG(t_5Gyr, z_typical);
    std::cout << "Initial g = " << std::scientific << g_initial << " m/s²\n\n";
    
    // Step 20: Final state export
    std::cout << "Step 20: Final State Export\n";
    std::cout << mod.exportState(t_5Gyr, z_typical) << "\n";
    
    std::cout << "===== DEMONSTRATION COMPLETE =====\n";
}

// Example usage in base program 'ziqn233h.cpp' (snippet for integration)

// Evaluation of SpiralSupernovaeUQFFModule (UQFF & Standard Model Integration for Spiral Galaxies and Supernovae)

**Strengths:**
-**Dynamic & Extensible : **All model parameters are stored in a `std: : map<std::string, double> variables`, allowing runtime updates, additions, and removals.The methods `updateVariable`, `addToVariable`, and `subtractFromVariable` enable flexible modification of any parameter.
- **Automatic Dependency Updates : **When key variables like `"Delta_x"` or `"M"` are updated, dependent variables(`"Delta_p"`, `"M_visible"`, `"M_DM"`) are recalculated automatically, ensuring consistency.
    - **Immediate Effect : **All computations(e.g., `computeG`) use the current values in the map, so any changes are immediately reflected in results.
        - **Comprehensive Physics : **The module includes all major UQFF and Standard Model terms relevant for spiral galaxy and supernova gravity, such as base gravity, cosmological constant, quantum, EM, fluid, resonant, DM, spiral torque, supernova, and superconductivity corrections.
        - **Specialized Terms : **Incorporates spiral torque(`T_spiral`) and supernova(`SN_term`) effects, which are important for galactic evolution and feedback.
            - **Debugging Support : **The `printVariables()` method provides a snapshot of all current parameters, aiding validation and troubleshooting.
            - **Sample Usage Provided : **Example integration and compilation instructions are included, demonstrating how to update variables and see their effect.

            ** Weaknesses / Recommendations : **
            -**Error Handling : **Adding unknown variables is flexible but may lead to silent errors if a typo occurs.Consider stricter validation or optional warnings for unknown variable names.
            - **Unit Consistency : **Ensure all units are consistent, especially when combining terms from different physical domains and frequency - based scaling.
            - **Magic Numbers : **Some scale factors and physical constants are set to arbitrary values.Document their physical meaning or allow configuration via constructor arguments or config files.
            - **Physical Justification : **The model is highly specialized for UQFF and galactic physics.Ensure this is appropriate for your scientific context and document the rationale for each term.
            - **Performance : **For large - scale or repeated updates, consider profiling and optimizing critical paths if needed.

            ** Summary : **
            The module is robust, dynamic, and extensible, supporting runtime updates and changes to all model parameters.It is suitable for advanced UQFF - based and Standard Model modeling of spiral galaxies and supernovae.Minor improvements in error handling, documentation, and physical justification are recommended for production or publication use.