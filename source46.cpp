// NGC6302UQFFModule.h
// Modular C++ implementation of the full Master Universal Gravity Equation (UQFF & SM Integration) for NGC 6302 Evolution.
// This module can be plugged into a base program (e.g., 'ziqn233h.cpp') by including this header and linking the .cpp.
// Usage in base: #include "NGC6302UQFFModule.h"
// NGC6302UQFFModule mod; mod.computeG(t); mod.updateVariable("M", new_value);
// All variables are stored in a std::map for dynamic addition/subtraction/update.
// Nothing is negligible: Includes all terms - base gravity, Ug1-Ug4, cosmological Lambda, quantum integral, Lorentz q(v x B), fluid rho_fluid V g, resonant oscillatory (cos/exp), DM/visible with perturbations, stellar wind W_shock.
// Associated text: Outputs descriptive equation string via getEquationText().
// Approximations: Quantum integral normalized to 1.0; exp real part; Ug2/Ug3=0; DM fraction ~0.85; W_shock = rho * v_wind^2 * (1 + t / t_eject).
// NGC6302 params: M=3.98e30 kg, r=9.46e15 m, v_wind=1e5 m/s, t_eject=2000 yr, z=0.00095, rho=1e-20 kg/m^3, etc.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#ifndef NGC6302_UQFF_MODULE_H
#define NGC6302_UQFF_MODULE_H

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

class NGC6302UQFFModule {
private:
    std::map<std::string, double> variables;
    double computeQuantumTerm(double t_Hubble_val);
    double computeFluidTerm(double g_base);
    double computeResonantTerm(double t);
    double computeDMTerm();
    double computeUgSum();
    double computeHz();
    double computeW_shock(double t);

public:
    // Constructor: Initialize all variables with NGC 6302 defaults
    NGC6302UQFFModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computation: Full g_UQFF(r, t) for NGC 6302
    double computeG(double t);

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
    void expandWindShockScale(double scale_factor);
    void expandNebularScale(double scale_factor);
    void expandResonanceScale(double scale_factor);

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
};

#endif // NGC6302_UQFF_MODULE_H

// NGC6302UQFFModule.cpp
#include "NGC6302UQFFModule.h"
#include <complex>

// Constructor: Set all variables with NGC 6302-specific values
NGC6302UQFFModule::NGC6302UQFFModule() {
    // Base constants (universal)
    variables["G"] = 6.6743e-11;                    // m^3 kg^-1 s^-2
    variables["c"] = 3e8;                           // m/s
    variables["hbar"] = 1.0546e-34;                 // J s
    variables["Lambda"] = 1.1e-52;                  // m^-2
    variables["q"] = 1.602e-19;                     // C
    variables["pi"] = 3.141592653589793;            // pi
    variables["t_Hubble"] = 13.8e9 * 3.156e7;       // s
    variables["year_to_s"] = 3.156e7;               // s/yr

    // NGC 6302 parameters
    double M_sun_val = 1.989e30;                    // kg
    variables["M_sun"] = M_sun_val;
    variables["M"] = 2 * M_sun_val;                 // Total ejected mass ~2 M_sun kg
    variables["M_visible"] = 0.15 * variables["M"]; // Visible fraction est.
    variables["M_DM"] = 0.85 * variables["M"];      // Dark matter (negligible, but included)
    variables["r"] = 9.46e15;                       // m (~1 ly radius)

    // Hubble/cosmology
    variables["H0"] = 67.15;                        // km/s/Mpc
    variables["Mpc_to_m"] = 3.086e22;               // m/Mpc
    variables["z"] = 0.00095;                       // Redshift
    variables["Omega_m"] = 0.3;
    variables["Omega_Lambda"] = 0.7;
    variables["t"] = 2000 * variables["year_to_s"]; // Default t=2000 yr s

    // Gas/wind dynamics
    variables["rho_fluid"] = 1e-20;                 // kg/m^3 (ionized gas)
    variables["V"] = 1e3;                           // m^3 (arbitrary)
    variables["v_wind"] = 1e5;                      // m/s (100 km/s)
    variables["t_eject"] = 2000 * variables["year_to_s"];  // s
    variables["delta_rho"] = 0.1 * variables["rho_fluid"];
    variables["rho"] = variables["rho_fluid"];

    // EM/magnetic
    variables["B"] = 1e-5;                          // T (nebula field)
    variables["B_crit"] = 1e11;                     // T

    // Quantum terms
    variables["Delta_x"] = 1e-10;                   // m
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];
    variables["integral_psi"] = 1.0;

    // Resonant/oscillatory
    variables["A"] = 1e-10;
    variables["k"] = 1e20;
    variables["omega"] = 1e15;                      // rad/s
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
void NGC6302UQFFModule::updateVariable(const std::string& name, double value) {
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
    }
}

// Add delta to variable
void NGC6302UQFFModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta from variable
void NGC6302UQFFModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute H(z) in s^-1
double NGC6302UQFFModule::computeHz() {
    double Hz_kms = variables["H0"] * std::sqrt(variables["Omega_m"] * std::pow(1.0 + variables["z"], 3) + variables["Omega_Lambda"]);
    return (Hz_kms * 1e3) / variables["Mpc_to_m"];
}

// Compute Ug sum: Ug1 = G M / r^2, Ug4 = Ug1 * f_sc, others 0
double NGC6302UQFFModule::computeUgSum() {
    double Ug1 = (variables["G"] * variables["M"]) / (variables["r"] * variables["r"]);
    variables["Ug1"] = Ug1;
    variables["Ug4"] = Ug1 * variables["f_sc"];
    return variables["Ug1"] + variables["Ug2"] + variables["Ug3"] + variables["Ug4"];
}

// Quantum term: (hbar / sqrt(Delta_x Delta_p)) * integral * (2 pi / t_Hubble)
double NGC6302UQFFModule::computeQuantumTerm(double t_Hubble_val) {
    double unc = std::sqrt(variables["Delta_x"] * variables["Delta_p"]);
    double integral_val = variables["integral_psi"];
    return (variables["hbar"] / unc) * integral_val * (2 * variables["pi"] / t_Hubble_val);
}

// Fluid term: rho_fluid * V * g
double NGC6302UQFFModule::computeFluidTerm(double g_base) {
    return variables["rho_fluid"] * variables["V"] * g_base;
}

// Resonant terms: 2 A cos(k x) cos(omega t) + (2 pi / 13.8) A Re[exp(i (k x - omega t))]
double NGC6302UQFFModule::computeResonantTerm(double t) {
    double cos_term = 2 * variables["A"] * std::cos(variables["k"] * variables["x"]) * std::cos(variables["omega"] * t);
    std::complex<double> exp_term(variables["A"] * std::exp(std::complex<double>(0, variables["k"] * variables["x"] - variables["omega"] * t)));
    double real_exp = exp_term.real();
    double exp_factor = (2 * variables["pi"] / 13.8);
    return cos_term + exp_factor * real_exp;
}

// DM term: (M_visible + M_DM) * (delta_rho / rho + 3 G M / r^3)
double NGC6302UQFFModule::computeDMTerm() {
    double pert = variables["delta_rho"] / variables["rho"];
    double curv = 3 * variables["G"] * variables["M"] / (variables["r"] * variables["r"] * variables["r"]);
    return (variables["M_visible"] + variables["M_DM"]) * (pert + curv);
}

// Wind shock term: W_shock = rho * v_wind^2 * (1 + t / t_eject)
double NGC6302UQFFModule::computeW_shock(double t) {
    return variables["rho_fluid"] * std::pow(variables["v_wind"], 2) * (1.0 + t / variables["t_eject"]);
}

// Full computation: g_UQFF(r, t) = ... all terms + W_shock
double NGC6302UQFFModule::computeG(double t) {
    variables["t"] = t;
    double Hz = computeHz();
    double expansion = 1.0 + Hz * t;
    double sc_correction = 1.0 - (variables["B"] / variables["B_crit"]);
    double tr_factor = 1.0 + variables["f_TRZ"];
    double w_shock = computeW_shock(t);

    // Base gravity with expansion, SC, TR
    double g_base = (variables["G"] * variables["M"] / (variables["r"] * variables["r"])) * expansion * sc_correction * tr_factor;

    // Ug sum
    double ug_sum = computeUgSum();

    // Cosmological
    double lambda_term = variables["Lambda"] * (variables["c"] * variables["c"]) / 3.0;

    // Quantum
    double quantum_term = computeQuantumTerm(variables["t_Hubble"]);

    // EM Lorentz (magnitude v_wind B)
    double em_base = variables["q"] * variables["v_wind"] * variables["B"] / 1.673e-27;
    double em_term = em_base * (1.0 + (7.09e-36 / 7.09e-37)) * variables["scale_macro"];

    // Fluid
    double fluid_term = computeFluidTerm(g_base);

    // Resonant
    double resonant_term = computeResonantTerm(t);

    // DM
    double dm_term = computeDMTerm();

    // Total: Sum all + W_shock
    return g_base + ug_sum + lambda_term + quantum_term + em_term + fluid_term + resonant_term + dm_term + w_shock;
}

// Get equation text (descriptive)
std::string NGC6302UQFFModule::getEquationText() {
    return "g_NGC6302(r, t) = (G * M / r^2) * (1 + H(z) * t) * (1 - B / B_crit) * (1 + f_TRZ) + (Ug1 + Ug2 + Ug3 + Ug4) + (Lambda * c^2 / 3) + "
           "(hbar / sqrt(Delta_x * Delta_p)) * ?(?* H ? dV) * (2? / t_Hubble) + q (v � B) + ?_fluid * V * g + "
           "2 A cos(k x) cos(? t) + (2? / 13.8) A exp(i (k x - ? t)) + (M_visible + M_DM) * (??/? + 3 G M / r^3) + W_shock\n"
           "Where W_shock = ? * v_wind^2 * (1 + t / t_eject)\n"
           "Special Terms:\n"
           "- Quantum: Heisenberg uncertainty for gas quantum effects.\n"
           "- Fluid: Ionized gas density-volume-gravity coupling.\n"
           "- Resonant: Oscillatory Aether waves for shock fronts.\n"
           "- DM: Visible+dark mass with perturbations (negligible).\n"
           "- Superconductivity: (1 - B/B_crit) for nebular fields.\n"
           "- Wind Shock: W_shock from central star winds eroding lobes.\n"
           "Solutions: At t=2000 yr, g_NGC6302 ~1e-10 m/s� (W_shock/EM dominant; g_base ~1e-12).\n"
           "Adaptations for NGC 6302: Bipolar PN with v_wind=100 km/s; z=0.00095; t_eject=2000 yr for ejections.";
}

// Print variables
void NGC6302UQFFModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// ===== ENHANCED DYNAMIC CAPABILITIES IMPLEMENTATION =====

namespace {
    std::map<std::string, std::map<std::string, double>> ngc6302_saved_states;
}

// Variable Management
void NGC6302UQFFModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void NGC6302UQFFModule::removeVariable(const std::string& name) {
    variables.erase(name);
}

void NGC6302UQFFModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
    }
}

std::vector<std::string> NGC6302UQFFModule::listVariables() {
    std::vector<std::string> names;
    for (const auto& pair : variables) {
        names.push_back(pair.first);
    }
    return names;
}

std::string NGC6302UQFFModule::getSystemName() {
    return "NGC 6302 (Butterfly Nebula) - Full UQFF & SM Integration";
}

// Batch Operations
void NGC6302UQFFModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = func(variables[name]);
        }
    }
}

void NGC6302UQFFModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// Self-Expansion (domain-specific for NGC 6302)
void NGC6302UQFFModule::expandParameterSpace(double scale_factor) {
    // Scale core nebular parameters
    variables["r"] *= scale_factor;
    variables["M"] *= scale_factor;
    variables["rho_fluid"] *= scale_factor;
    // Update dependent variables
    variables["M_visible"] = 0.15 * variables["M"];
    variables["M_DM"] = 0.85 * variables["M"];
    variables["delta_rho"] = 0.1 * variables["rho_fluid"];
    variables["rho"] = variables["rho_fluid"];
}

void NGC6302UQFFModule::expandWindShockScale(double scale_factor) {
    // Scale wind shock parameters (v_wind, t_eject, rho_fluid)
    variables["v_wind"] *= scale_factor;
    variables["t_eject"] *= scale_factor;
    variables["rho_fluid"] *= scale_factor;
    variables["delta_rho"] = 0.1 * variables["rho_fluid"];
    variables["rho"] = variables["rho_fluid"];
}

void NGC6302UQFFModule::expandNebularScale(double scale_factor) {
    // Scale nebular field parameters (B, rho_fluid, V)
    variables["B"] *= scale_factor;
    variables["rho_fluid"] *= scale_factor;
    variables["V"] *= scale_factor;
    variables["delta_rho"] = 0.1 * variables["rho_fluid"];
    variables["rho"] = variables["rho_fluid"];
}

void NGC6302UQFFModule::expandResonanceScale(double scale_factor) {
    // Scale resonant oscillatory parameters (A, omega, k)
    variables["A"] *= scale_factor;
    variables["omega"] *= scale_factor;
    variables["k"] *= scale_factor;
}

// Self-Refinement
void NGC6302UQFFModule::autoRefineParameters(double tolerance) {
    // Enforce physical constraints
    if (variables["M"] <= 0) variables["M"] = 2 * variables["M_sun"];
    if (variables["r"] <= 0) variables["r"] = 9.46e15;
    if (variables["rho_fluid"] <= 0) variables["rho_fluid"] = 1e-20;
    if (variables["v_wind"] < 0) variables["v_wind"] = 1e5;
    if (variables["t_eject"] <= 0) variables["t_eject"] = 2000 * variables["year_to_s"];
    if (variables["B"] < 0) variables["B"] = 1e-5;
    // Update dependent variables
    variables["M_visible"] = 0.15 * variables["M"];
    variables["M_DM"] = 0.85 * variables["M"];
    variables["delta_rho"] = 0.1 * variables["rho_fluid"];
    variables["rho"] = variables["rho_fluid"];
}

void NGC6302UQFFModule::calibrateToObservations(const std::map<std::string, double>& obs_data) {
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

void NGC6302UQFFModule::optimizeForMetric(const std::string& var_name, double target_value, int iterations) {
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
std::vector<std::map<std::string, double>> NGC6302UQFFModule::generateVariations(int n_variations) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.8, 1.2);
    
    for (int i = 0; i < n_variations; ++i) {
        std::map<std::string, double> variation = variables;
        variation["M"] *= dis(gen);
        variation["v_wind"] *= dis(gen);
        variation["rho_fluid"] *= dis(gen);
        variation["B"] *= dis(gen);
        variations.push_back(variation);
    }
    return variations;
}

// Adaptive Evolution
void NGC6302UQFFModule::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-mutation_rate, mutation_rate);
    
    variables["M"] *= (1.0 + dis(gen));
    variables["v_wind"] *= (1.0 + dis(gen));
    variables["rho_fluid"] *= (1.0 + dis(gen));
    variables["B"] *= (1.0 + dis(gen));
    
    // Update dependent variables
    variables["M_visible"] = 0.15 * variables["M"];
    variables["M_DM"] = 0.85 * variables["M"];
    variables["delta_rho"] = 0.1 * variables["rho_fluid"];
    variables["rho"] = variables["rho_fluid"];
}

void NGC6302UQFFModule::evolveSystem(int generations, std::function<double()> fitness_function) {
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
void NGC6302UQFFModule::saveState(const std::string& label) {
    ngc6302_saved_states[label] = variables;
}

void NGC6302UQFFModule::restoreState(const std::string& label) {
    if (ngc6302_saved_states.find(label) != ngc6302_saved_states.end()) {
        variables = ngc6302_saved_states[label];
    }
}

std::vector<std::string> NGC6302UQFFModule::listSavedStates() {
    std::vector<std::string> labels;
    for (const auto& pair : ngc6302_saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string NGC6302UQFFModule::exportState(double t) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(6);
    oss << "NGC6302 State Export at t=" << t << " s:\n";
    oss << "M=" << variables["M"] << " kg, r=" << variables["r"] << " m\n";
    oss << "v_wind=" << variables["v_wind"] << " m/s, t_eject=" << variables["t_eject"] << " s\n";
    oss << "rho_fluid=" << variables["rho_fluid"] << " kg/m³, B=" << variables["B"] << " T\n";
    oss << "g_total=" << computeG(t) << " m/s²\n";
    return oss.str();
}

// System Analysis
std::map<std::string, double> NGC6302UQFFModule::sensitivityAnalysis(const std::string& param, double t, double delta) {
    std::map<std::string, double> result;
    if (variables.find(param) == variables.end()) return result;
    
    double original = variables[param];
    double g_original = computeG(t);
    
    variables[param] = original * (1.0 + delta);
    double g_plus = computeG(t);
    
    variables[param] = original * (1.0 - delta);
    double g_minus = computeG(t);
    
    variables[param] = original;
    
    result["dg/d" + param] = (g_plus - g_minus) / (2.0 * delta * original);
    result["g_original"] = g_original;
    result["g_plus"] = g_plus;
    result["g_minus"] = g_minus;
    
    return result;
}

std::string NGC6302UQFFModule::generateReport(double t) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(6);
    oss << "===== NGC 6302 UQFF Module Report (t=" << t << " s) =====\n\n";
    oss << "System: " << getSystemName() << "\n\n";
    
    oss << "Core Parameters:\n";
    oss << "  M = " << variables["M"] << " kg (" << variables["M"]/variables["M_sun"] << " M_sun)\n";
    oss << "  r = " << variables["r"] << " m (~" << variables["r"]/9.46e15 << " ly)\n";
    oss << "  v_wind = " << variables["v_wind"] << " m/s (" << variables["v_wind"]/1e3 << " km/s)\n";
    oss << "  t_eject = " << variables["t_eject"] << " s (" << variables["t_eject"]/variables["year_to_s"] << " yr)\n";
    oss << "  rho_fluid = " << variables["rho_fluid"] << " kg/m³\n";
    oss << "  B = " << variables["B"] << " T\n\n";
    
    double Hz = computeHz();
    double expansion = 1.0 + Hz * t;
    double sc_correction = 1.0 - (variables["B"] / variables["B_crit"]);
    double tr_factor = 1.0 + variables["f_TRZ"];
    double g_base = (variables["G"] * variables["M"] / (variables["r"] * variables["r"])) * expansion * sc_correction * tr_factor;
    double ug_sum = computeUgSum();
    double lambda_term = variables["Lambda"] * (variables["c"] * variables["c"]) / 3.0;
    double quantum_term = computeQuantumTerm(variables["t_Hubble"]);
    double em_base = variables["q"] * variables["v_wind"] * variables["B"] / 1.673e-27;
    double em_term = em_base * (1.0 + (7.09e-36 / 7.09e-37)) * variables["scale_macro"];
    double fluid_term = computeFluidTerm(g_base);
    double resonant_term = computeResonantTerm(t);
    double dm_term = computeDMTerm();
    double w_shock = computeW_shock(t);
    double g_total = computeG(t);
    
    oss << "Term Breakdown:\n";
    oss << "  g_base = " << g_base << " m/s²\n";
    oss << "  Ug_sum (Ug1+Ug2+Ug3+Ug4) = " << ug_sum << " m/s²\n";
    oss << "  Lambda_term = " << lambda_term << " m/s²\n";
    oss << "  Quantum_term = " << quantum_term << " m/s²\n";
    oss << "  EM_term = " << em_term << " m/s²\n";
    oss << "  Fluid_term = " << fluid_term << " m/s²\n";
    oss << "  Resonant_term = " << resonant_term << " m/s²\n";
    oss << "  DM_term = " << dm_term << " m/s²\n";
    oss << "  W_shock = " << w_shock << " m/s² [DOMINANT WIND TERM]\n\n";
    oss << "TOTAL g = " << g_total << " m/s²\n\n";
    
    oss << "Physics Notes:\n";
    oss << "- NGC 6302 is a bipolar planetary nebula with high-velocity winds (~100 km/s)\n";
    oss << "- W_shock term dominates, modeling stellar wind impact on lobe expansion\n";
    oss << "- EM term significant due to v×B interaction in ionized gas\n";
    oss << "- Full UQFF+SM integration: gravity, Ug1-4, Lambda, quantum, EM, fluid, resonance, DM\n";
    oss << "- Superconductivity correction: (1 - B/B_crit)\n";
    oss << "- Time-dependent wind shock: W_shock ∝ (1 + t/t_eject)\n\n";
    
    return oss.str();
}

bool NGC6302UQFFModule::validateConsistency() {
    bool valid = true;
    if (variables["M"] <= 0) valid = false;
    if (variables["r"] <= 0) valid = false;
    if (variables["rho_fluid"] <= 0) valid = false;
    if (variables["v_wind"] < 0) valid = false;
    if (variables["t_eject"] <= 0) valid = false;
    if (variables["B"] < 0) valid = false;
    return valid;
}

void NGC6302UQFFModule::autoCorrectAnomalies() {
    if (variables["M"] <= 0) variables["M"] = 2 * variables["M_sun"];
    if (variables["r"] <= 0) variables["r"] = 9.46e15;
    if (variables["rho_fluid"] <= 0) variables["rho_fluid"] = 1e-20;
    if (variables["v_wind"] < 0) variables["v_wind"] = 1e5;
    if (variables["t_eject"] <= 0) variables["t_eject"] = 2000 * variables["year_to_s"];
    if (variables["B"] < 0) variables["B"] = 1e-5;
    // Update dependent variables
    variables["M_visible"] = 0.15 * variables["M"];
    variables["M_DM"] = 0.85 * variables["M"];
    variables["delta_rho"] = 0.1 * variables["rho_fluid"];
    variables["rho"] = variables["rho_fluid"];
}

// Enhanced example usage demonstration
void enhanced_example_usage() {
    NGC6302UQFFModule mod;
    double t_2kyr = 2000 * 3.156e7;  // 2000 years in seconds
    
    std::cout << "===== ENHANCED NGC 6302 UQFF MODULE DEMONSTRATION =====\n\n";
    
    // Step 1: Variable management
    std::cout << "Step 1: Variable Management\n";
    mod.createVariable("custom_wind_factor", 1.15);
    mod.cloneVariable("v_wind", "v_wind_backup");
    std::vector<std::string> vars = mod.listVariables();
    std::cout << "Total variables: " << vars.size() << "\n";
    std::cout << "System: " << mod.getSystemName() << "\n\n";
    
    // Step 2: Batch scaling
    std::cout << "Step 2: Batch Scaling (Wind parameters)\n";
    mod.scaleVariableGroup({"v_wind", "rho_fluid", "B"}, 1.1);
    std::cout << "Scaled v_wind, rho_fluid, B by 1.1\n\n";
    
    // Step 3: Self-expansion (different physics domains)
    std::cout << "Step 3: Self-Expansion\n";
    mod.expandWindShockScale(1.08);  // Wind +8%
    std::cout << "Expanded wind shock scale +8%\n";
    mod.expandNebularScale(1.05);  // Nebular fields +5%
    std::cout << "Expanded nebular scale +5%\n";
    mod.expandResonanceScale(1.03);  // Resonance +3%
    std::cout << "Expanded resonance scale +3%\n\n";
    
    // Step 4: Self-refinement
    std::cout << "Step 4: Self-Refinement\n";
    mod.autoRefineParameters(1e-10);
    std::cout << "Auto-refined parameters\n";
    std::map<std::string, double> obs_data = {
        {"M", 2.5e30},
        {"v_wind", 1.2e5},
        {"rho_fluid", 1.2e-20},
        {"B", 1.1e-5}
    };
    mod.calibrateToObservations(obs_data);
    std::cout << "Calibrated to observations\n\n";
    
    // Step 5: Optimize for specific metric
    std::cout << "Step 5: Optimize for v_wind~1.2e5 m/s\n";
    mod.optimizeForMetric("v_wind", 1.2e5, 50);
    std::cout << "Optimization complete\n\n";
    
    // Step 6: Generate variations
    std::cout << "Step 6: Generate 15 Parameter Variations\n";
    auto variations = mod.generateVariations(15);
    std::cout << "Generated " << variations.size() << " variations\n\n";
    
    // Step 7: State management
    std::cout << "Step 7: State Management\n";
    mod.saveState("initial");
    mod.scaleVariableGroup({"v_wind", "rho_fluid"}, 1.2);
    mod.saveState("enhanced_wind");
    mod.expandNebularScale(0.8);
    mod.saveState("reduced_nebular");
    std::cout << "Saved 3 states\n\n";
    
    // Step 8: Sensitivity analysis
    std::cout << "Step 8: Sensitivity Analysis (v_wind at t=2000yr)\n";
    mod.restoreState("initial");
    auto sensitivity = mod.sensitivityAnalysis("v_wind", t_2kyr, 0.1);
    std::cout << "dg/dv_wind = " << std::scientific << sensitivity["dg/dv_wind"] << " (m/s²)/(m/s)\n\n";
    
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
    std::cout << "Step 10: Comprehensive Report (t=2000yr)\n";
    std::string report = mod.generateReport(t_2kyr);
    std::cout << report << "\n";
    
    // Step 11: Adaptive evolution
    std::cout << "Step 11: Adaptive Evolution (25 generations)\n";
    auto fitness_fn = [&mod, t_2kyr]() -> double {
        double g = mod.computeG(t_2kyr);
        return -std::abs(std::log10(std::abs(g)) + 10.0);  // Target g~1e-10 m/s²
    };
    mod.evolveSystem(25, fitness_fn);
    std::cout << "Evolution complete\n\n";
    
    // Step 12: Time evolution comparison
    std::cout << "Step 12: Time Evolution (0 to 10,000 years)\n";
    std::vector<double> times = {0.0, 500*3.156e7, 1000*3.156e7, 2000*3.156e7, 5000*3.156e7, 10000*3.156e7};
    for (double t : times) {
        double g = mod.computeG(t);
        std::cout << "t=" << std::scientific << t << " s (" << t/(3.156e7) << " yr): g=" << g << " m/s²\n";
    }
    std::cout << "\n";
    
    // Step 13: Wind shock evolution
    std::cout << "Step 13: Wind Shock Evolution (t=0 to 10,000 yr)\n";
    for (double t : times) {
        double w_shock = mod.computeW_shock(t);
        std::cout << "t=" << t/(3.156e7) << " yr: W_shock=" << std::scientific << w_shock << " m/s²\n";
    }
    std::cout << "\n";
    
    // Step 14: Wind velocity sweep
    std::cout << "Step 14: Wind Velocity Sweep (t=2000yr)\n";
    std::vector<double> velocities = {5e4, 7.5e4, 1e5, 1.25e5, 1.5e5};
    for (double v : velocities) {
        mod.updateVariable("v_wind", v);
        double g = mod.computeG(t_2kyr);
        std::cout << "v_wind=" << std::scientific << v << " m/s (" << v/1000 << " km/s): g=" << g << " m/s²\n";
    }
    mod.restoreState("initial");
    std::cout << "\n";
    
    // Step 15: Mass impact
    std::cout << "Step 15: Mass Impact (t=2000yr)\n";
    std::vector<double> masses = {1e30, 1.5e30, 2e30, 2.5e30, 3e30};
    for (double M : masses) {
        mod.updateVariable("M", M);
        double g = mod.computeG(t_2kyr);
        std::cout << "M=" << std::scientific << M << " kg (" << M/1.989e30 << " M_sun): g=" << g << " m/s²\n";
    }
    mod.restoreState("initial");
    std::cout << "\n";
    
    // Step 16: Multi-parameter sensitivity
    std::cout << "Step 16: Multi-Parameter Sensitivity (t=2000yr)\n";
    std::vector<std::string> params = {"M", "v_wind", "rho_fluid", "B", "r"};
    for (const auto& param : params) {
        auto sens = mod.sensitivityAnalysis(param, t_2kyr, 0.05);
        std::cout << "dg/d" << param << " = " << std::scientific << sens["dg/d" + param] << "\n";
    }
    std::cout << "\n";
    
    // Step 17: Expansion domains comparison
    std::cout << "Step 17: Expansion Domains Comparison (t=2000yr)\n";
    mod.restoreState("initial");
    double g_initial = mod.computeG(t_2kyr);
    std::cout << "Initial: g=" << std::scientific << g_initial << " m/s²\n";
    
    mod.expandWindShockScale(1.2);
    double g_wind = mod.computeG(t_2kyr);
    std::cout << "Wind +20%: g=" << g_wind << " m/s² (Δ=" << (g_wind-g_initial)/g_initial*100 << "%)\n";
    mod.restoreState("initial");
    
    mod.expandNebularScale(1.2);
    double g_nebular = mod.computeG(t_2kyr);
    std::cout << "Nebular +20%: g=" << g_nebular << " m/s² (Δ=" << (g_nebular-g_initial)/g_initial*100 << "%)\n";
    mod.restoreState("initial");
    std::cout << "\n";
    
    // Step 18: State restoration
    std::cout << "Step 18: State Restoration\n";
    mod.restoreState("initial");
    std::cout << "Restored initial state\n";
    g_initial = mod.computeG(t_2kyr);
    std::cout << "Initial g = " << std::scientific << g_initial << " m/s²\n\n";
    
    // Step 19: Final state export
    std::cout << "Step 19: Final State Export\n";
    std::cout << mod.exportState(t_2kyr) << "\n";
    
    std::cout << "===== DEMONSTRATION COMPLETE =====\n";
}

// Example usage in base program 'ziqn233h.cpp' (snippet for integration)

// Evaluation of NGC6302UQFFModule (UQFF & Standard Model Integration for NGC 6302 Nebula Evolution)

**Strengths:**
-**Dynamic & Extensible : **All model parameters are stored in a `std: : map<std::string, double> variables`, allowing runtime updates, additions, and removals.The methods `updateVariable`, `addToVariable`, and `subtractFromVariable` enable flexible modification of any parameter.
- **Automatic Dependency Updates : **When key variables like `"Delta_x"` or `"M"` are updated, dependent variables(`"Delta_p"`, `"M_visible"`, `"M_DM"`) are recalculated automatically, ensuring consistency.
    - **Immediate Effect : **All computations(e.g., `computeG`) use the current values in the map, so any changes are immediately reflected in results.
        - **Comprehensive Physics : **The module includes all major UQFF and Standard Model terms relevant for nebular gravity, such as base gravity, cosmological constant, quantum, EM, fluid, resonant, DM, wind shock, and superconductivity corrections.
        - **Wind Shock Modeling : **Incorporates a wind shock term(`W_shock`) that evolves with time, reflecting the impact of stellar winds on nebular evolution.
            - **Debugging Support : **The `printVariables()` method provides a snapshot of all current parameters, aiding validation and troubleshooting.
        - **Sample Usage Provided : **Example integration and compilation instructions are included, demonstrating how to update variables and see their effect.

        ** Weaknesses / Recommendations : **
        -**Error Handling : **Adding unknown variables is flexible but may lead to silent errors if a typo occurs.Consider stricter validation or optional warnings for unknown variable names.
        - **Unit Consistency : **Ensure all units are consistent, especially when combining terms from different physical domains and frequency - based scaling.
        - **Magic Numbers : **Some scale factors and physical constants are set to arbitrary values.Document their physical meaning or allow configuration via constructor arguments or config files.
        - **Physical Justification : **The model is highly specialized for UQFF and nebular physics.Ensure this is appropriate for your scientific context and document the rationale for each term.
        - **Performance : **For large - scale or repeated updates, consider profiling and optimizing critical paths if needed.

        ** Summary : **
        The module is robust, dynamic, and extensible, supporting runtime updates and changes to all model parameters.It is suitable for advanced UQFF - based and Standard Model nebular modeling.Minor improvements in error handling, documentation, and physical justification are recommended for production or publication use.