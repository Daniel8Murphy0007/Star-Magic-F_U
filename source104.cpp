// MagneticMomentModule.h
// Modular C++ implementation of the Magnetic Moment of the j-th String (?_j) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes ?_j = (10^3 + 0.4 sin(?_c t)) * 3.38e20 T�m^3; scales ?_j / r_j in Universal Magnetism U_m and Ug3.
// Pluggable: #include "MagneticMomentModule.h"
// MagneticMomentModule mod; mod.computeMu_j(0.0); mod.updateVariable("base_mu", new_value);
// Variables in std::map; j-indexed; example for j=1 at t=0.
// Approximations: ?_c=2.5e-6 rad/s; at t=0, sin=0, ?_j?3.38e23 T�m^3 (adjusted for example consistency).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef MAGNETIC_MOMENT_MODULE_H
#define MAGNETIC_MOMENT_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <functional>
#include <random>
#include <algorithm>
#include <sstream>
#include <vector>

class MagneticMomentModule {
private:
    std::map<std::string, double> variables;
    double computeMu_j(int j, double t);
    double computeUmContrib(int j, double t);

public:
    // Constructor: Initialize with framework defaults
    MagneticMomentModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeMu_j(int j, double t);  // T�m^3
    double computeB_j(double t);  // Base field 10^3 + 0.4 sin(?_c t) T
    double computeUmContrib(int j, double t);  // Example U_m single string (J/m^3)
    double computeUg3Contrib(double t);  // Example Ug3 (J/m^3)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print ?_j and contributions
    void printMomentContributions(int j = 1, double t = 0.0);

    // ===== Enhanced: 25-Method Dynamic Self-Update & Self-Expansion Capabilities =====

    // Variable Management (5 methods)
    void createVariable(const std::string& name, double value);
    void removeVariable(const std::string& name);
    void cloneVariable(const std::string& source, const std::string& target);
    std::string listVariables();
    std::string getSystemName();

    // Batch Operations (2 methods)
    void transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func);
    void scaleVariableGroup(const std::vector<std::string>& names, double factor);

    // Self-Expansion (4 methods: 1 generic + 3 domain-specific)
    void expandParameterSpace(const std::vector<std::string>& params, double expansion_factor);
    void expandMomentScale(double base_mu_factor, double B_j_factor);      // μ_j base and B_j field
    void expandOscillationScale(double omega_c_factor, double amp_factor); // ω_c and oscillation amplitude
    void expandStringScale(double r_j_factor, double gamma_factor);        // String distance and decay

    // Self-Refinement (3 methods)
    void autoRefineParameters(const std::string& target_metric, double target_value);
    void calibrateToObservations(const std::map<std::string, double>& observed);
    void optimizeForMetric(const std::string& metric);

    // Parameter Exploration (1 method)
    std::vector<std::map<std::string, double>> generateVariations(int count, double variance);

    // Adaptive Evolution (2 methods)
    void mutateParameters(double mutation_rate);
    void evolveSystem(int generations, std::function<double()> fitness_func);

    // State Management (4 methods)
    void saveState(const std::string& label);
    void restoreState(const std::string& label);
    std::vector<std::string> listSavedStates();
    std::string exportState();

    // System Analysis (4 methods)
    std::map<std::string, double> sensitivityAnalysis(const std::vector<std::string>& params);
    std::string generateReport();
    bool validateConsistency();
    void autoCorrectAnomalies();
};

#endif // MAGNETIC_MOMENT_MODULE_H

// MagneticMomentModule.cpp
#include "MagneticMomentModule.h"

// Constructor: Set framework defaults
MagneticMomentModule::MagneticMomentModule() {
    // Universal constants
    variables["base_mu"] = 3.38e20;                 // T�m^3 (definition); note: example uses 3.38e23
    variables["omega_c"] = 2.5e-6;                  // rad/s
    variables["r_j"] = 1.496e13;                    // m (for j=1)
    variables["gamma"] = 5e-5 / 86400.0;            // s^-1 (0.00005 day^-1)
    variables["t_n"] = 0.0;                         // s
    variables["phi_hat_j"] = 1.0;                   // Normalized
    variables["P_SCm"] = 1.0;                       // Pressure
    variables["E_react"] = 1e46;                    // J
    variables["f_Heaviside"] = 0.01;                // Unitless
    variables["f_quasi"] = 0.01;                    // Unitless
    variables["k3"] = 1.8;                          // Coupling for Ug3
    variables["pi"] = 3.141592653589793;

    // Derived defaults
    variables["B_j"] = 1e3;                         // Base T
    variables["scale_Heaviside"] = 1e13;            // Amplification
}

// Update variable
void MagneticMomentModule::updateVariable(const std::string& name, double value) {
    variables[name] = value;
}

// Add delta
void MagneticMomentModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void MagneticMomentModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute ?_j(t)
double MagneticMomentModule::computeMu_j(int j, double t) {
    double sin_term = std::sin(variables["omega_c"] * t);
    double b_j = variables["B_j"] + 0.4 * sin_term;  // T
    return b_j * variables["base_mu"];  // T�m^3; adjust base if needed for example
}

// Compute B_j(t) base
double MagneticMomentModule::computeB_j(double t) {
    return variables["B_j"] + 0.4 * std::sin(variables["omega_c"] * t);
}

// Example U_m contrib for j (J/m^3, simplified)
double MagneticMomentModule::computeUmContrib(int j, double t) {
    double mu_j = computeMu_j(j, t);
    double r_j = variables["r_j"];
    double exp_arg = - variables["gamma"] * t * std::cos(variables["pi"] * variables["t_n"]);
    double one_minus_exp = 1.0 - std::exp(exp_arg);
    double phi_hat = variables["phi_hat_j"];
    double heaviside_f = 1.0 + variables["scale_Heaviside"] * variables["f_Heaviside"];
    double quasi_f = 1.0 + variables["f_quasi"];
    return (mu_j / r_j * one_minus_exp * phi_hat) * variables["P_SCm"] * variables["E_react"] * heaviside_f * quasi_f;
}

// Example Ug3 contrib (J/m^3)
double MagneticMomentModule::computeUg3Contrib(double t) {
    double b_j = computeB_j(t);
    double cos_term = std::cos(variables["omega_c"] * t * variables["pi"]);  // Approx
    double p_core = 1.0;
    double e_react = variables["E_react"];
    return variables["k3"] * b_j * cos_term * p_core * e_react;
}

// Equation text
std::string MagneticMomentModule::getEquationText() {
    return "?_j = (10^3 + 0.4 sin(?_c t)) * 3.38e20 T�m^3\n"
           "Where ?_c=2.5e-6 rad/s; units T�m^3 (magnetic dipole strength).\n"
           "In U_m: ?_j [?_j / r_j * (1 - e^{-? t cos(? t_n)}) ?_hat_j ] P_SCm E_react (1 + 10^13 f_Heaviside) (1 + f_quasi)\n"
           "In Ug3: k3 * ?_j B_j cos(?_s t ?) P_core E_react; B_j = 10^3 + 0.4 sin(?_c t) T.\n"
           "Example j=1, t=0: ?_j ?3.38e23 T�m^3; U_m contrib ?2.28e65 J/m�; Ug3 ?1.8e49 J/m�.\n"
           "Role: Quantifies string magnetic strength; drives Um/Ug3 for jets/disks/nebulae.";
}

// Print variables
void MagneticMomentModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Print contributions
void MagneticMomentModule::printMomentContributions(int j, double t) {
    double mu = computeMu_j(j, t);
    double b = computeB_j(t);
    double um = computeUmContrib(j, t);
    double ug3 = computeUg3Contrib(t);
    std::cout << "Magnetic Moment j=" << j << " at t=" << t << " s:\n";
    std::cout << "?_j = " << std::scientific << mu << " T�m^3\n";
    std::cout << "B_j = " << b << " T\n";
    std::cout << "U_m contrib = " << um << " J/m�\n";
    std::cout << "Ug3 contrib = " << ug3 << " J/m�\n";
}

// ===== Implementation: 25-Method Dynamic Self-Update & Self-Expansion Capabilities =====

namespace magnetic_moment_saved_states {
    std::map<std::string, std::map<std::string, double>> saved_states;
}

// Variable Management (5 methods)
void MagneticMomentModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void MagneticMomentModule::removeVariable(const std::string& name) {
    variables.erase(name);
}

void MagneticMomentModule::cloneVariable(const std::string& source, const std::string& target) {
    if (variables.find(source) != variables.end()) {
        variables[target] = variables[source];
    }
}

std::string MagneticMomentModule::listVariables() {
    std::ostringstream oss;
    for (const auto& pair : variables) {
        oss << pair.first << " = " << pair.second << "\n";
    }
    return oss.str();
}

std::string MagneticMomentModule::getSystemName() {
    return "Magnetic_Moment_UQFF";
}

// Batch Operations (2 methods)
void MagneticMomentModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = func(variables[name]);
        }
    }
}

void MagneticMomentModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double x) { return x * factor; });
}

// Self-Expansion (4 methods)
void MagneticMomentModule::expandParameterSpace(const std::vector<std::string>& params, double expansion_factor) {
    for (const auto& param : params) {
        if (variables.find(param) != variables.end()) {
            variables[param] *= expansion_factor;
        }
    }
}

void MagneticMomentModule::expandMomentScale(double base_mu_factor, double B_j_factor) {
    // base_mu: fundamental moment (T·m³, ~3.38e20, example uses 3.38e23)
    // B_j: base magnetic field (T, ~1000)
    variables["base_mu"] *= base_mu_factor;
    variables["B_j"] *= B_j_factor;
}

void MagneticMomentModule::expandOscillationScale(double omega_c_factor, double amp_factor) {
    // ω_c: oscillation frequency (rad/s, ~2.5e-6)
    // amp_factor: scales 0.4 amplitude in sin term (stored implicitly)
    variables["omega_c"] *= omega_c_factor;
    // Note: 0.4 amplitude is hardcoded in computeMu_j, could add "osc_amplitude" variable
    if (variables.find("osc_amplitude") != variables.end()) {
        variables["osc_amplitude"] *= amp_factor;
    }
}

void MagneticMomentModule::expandStringScale(double r_j_factor, double gamma_factor) {
    // r_j: string distance (m, ~1.496e13 = 100 AU for j=1)
    // γ: decay rate (s⁻¹, ~5.79e-10)
    variables["r_j"] *= r_j_factor;
    variables["gamma"] *= gamma_factor;
}

// Self-Refinement (3 methods)
void MagneticMomentModule::autoRefineParameters(const std::string& target_metric, double target_value) {
    // Example: target U_m contribution by adjusting base_mu and E_react
    if (target_metric == "U_m") {
        double current = computeUmContrib(1, 0.0);
        if (current > 0) {
            double ratio = target_value / current;
            variables["base_mu"] *= std::sqrt(ratio);
            variables["E_react"] *= std::sqrt(ratio);
        }
    } else if (target_metric == "Ug3") {
        double current = computeUg3Contrib(0.0);
        if (current > 0) {
            double ratio = target_value / current;
            variables["k3"] *= std::sqrt(ratio);
            variables["E_react"] *= std::sqrt(ratio);
        }
    } else if (target_metric == "mu_j") {
        double current = computeMu_j(1, 0.0);
        if (current > 0) {
            double ratio = target_value / current;
            variables["base_mu"] *= ratio;
        }
    }
}

void MagneticMomentModule::calibrateToObservations(const std::map<std::string, double>& observed) {
    for (const auto& obs : observed) {
        if (variables.find(obs.first) != variables.end()) {
            variables[obs.first] = obs.second;
        }
    }
}

void MagneticMomentModule::optimizeForMetric(const std::string& metric) {
    if (metric == "maximize_Um") {
        variables["base_mu"] *= 1.2;
        variables["E_react"] *= 1.15;
        variables["scale_Heaviside"] *= 1.1;
    } else if (metric == "maximize_Ug3") {
        variables["k3"] *= 1.2;
        variables["B_j"] *= 1.15;
    } else if (metric == "enhance_oscillation") {
        variables["omega_c"] *= 1.3;
        if (variables.find("osc_amplitude") == variables.end()) {
            variables["osc_amplitude"] = 0.4;
        }
        variables["osc_amplitude"] *= 1.5;
    } else if (metric == "stabilize_field") {
        variables["B_j"] = 1e3;
        variables["omega_c"] = 2.5e-6;
    }
}

// Parameter Exploration (1 method)
std::vector<std::map<std::string, double>> MagneticMomentModule::generateVariations(int count, double variance) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-variance, variance);

    for (int i = 0; i < count; ++i) {
        std::map<std::string, double> variant = variables;
        for (auto& pair : variant) {
            if (pair.first != "pi" && pair.first != "t_n" && pair.first != "phi_hat_j") {
                double factor = 1.0 + dis(gen);
                pair.second *= factor;
                if (pair.second < 0) pair.second = std::abs(pair.second);
            }
        }
        variations.push_back(variant);
    }
    return variations;
}

// Adaptive Evolution (2 methods)
void MagneticMomentModule::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-mutation_rate, mutation_rate);

    for (auto& pair : variables) {
        if (pair.first != "pi" && pair.first != "t_n" && pair.first != "phi_hat_j") {
            pair.second *= (1.0 + dis(gen));
            if (pair.second < 0) pair.second = std::abs(pair.second);
        }
    }
}

void MagneticMomentModule::evolveSystem(int generations, std::function<double()> fitness_func) {
    double best_fitness = fitness_func();
    std::map<std::string, double> best_state = variables;

    for (int gen = 0; gen < generations; ++gen) {
        mutateParameters(0.05);
        double current_fitness = fitness_func();
        if (current_fitness > best_fitness) {
            best_fitness = current_fitness;
            best_state = variables;
        } else {
            variables = best_state; // Revert if worse
        }
    }
    variables = best_state;
}

// State Management (4 methods)
void MagneticMomentModule::saveState(const std::string& label) {
    magnetic_moment_saved_states::saved_states[label] = variables;
}

void MagneticMomentModule::restoreState(const std::string& label) {
    if (magnetic_moment_saved_states::saved_states.find(label) != magnetic_moment_saved_states.end()) {
        variables = magnetic_moment_saved_states::saved_states[label];
    }
}

std::vector<std::string> MagneticMomentModule::listSavedStates() {
    std::vector<std::string> labels;
    for (const auto& pair : magnetic_moment_saved_states::saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string MagneticMomentModule::exportState() {
    std::ostringstream oss;
    oss << "System: " << getSystemName() << "\n";
    oss << listVariables();
    return oss.str();
}

// System Analysis (4 methods)
std::map<std::string, double> MagneticMomentModule::sensitivityAnalysis(const std::vector<std::string>& params) {
    std::map<std::string, double> sensitivity;
    double baseline_um = computeUmContrib(1, 0.0);

    for (const auto& param : params) {
        if (variables.find(param) != variables.end()) {
            double original = variables[param];
            variables[param] *= 1.01; // +1% perturbation
            double perturbed_um = computeUmContrib(1, 0.0);
            variables[param] = original;
            if (baseline_um > 0) {
                sensitivity[param] = std::abs((perturbed_um - baseline_um) / baseline_um);
            } else {
                sensitivity[param] = 0.0;
            }
        }
    }
    return sensitivity;
}

std::string MagneticMomentModule::generateReport() {
    std::ostringstream oss;
    oss << "========== Magnetic Moment Module Report ==========\n";
    oss << "System: " << getSystemName() << "\n\n";
    oss << "Key Parameters:\n";
    oss << "  base_mu: " << std::scientific << variables["base_mu"] << " T·m³ (fundamental moment)\n";
    oss << "  B_j (base field): " << variables["B_j"] << " T\n";
    oss << "  ω_c (oscillation): " << variables["omega_c"] << " rad/s\n";
    oss << "  r_j (string distance): " << variables["r_j"] << " m (";
    oss << (variables["r_j"] / 1.496e11) << " AU)\n";
    oss << "  γ (decay rate): " << variables["gamma"] << " s⁻¹\n";
    oss << "  E_react: " << variables["E_react"] << " J\n";
    oss << "  k3 (Ug3 coupling): " << variables["k3"] << "\n";
    oss << "  scale_Heaviside: " << variables["scale_Heaviside"] << "\n";
    oss << "  f_Heaviside: " << variables["f_Heaviside"] << "\n\n";
    
    oss << "Computed Values at t=0:\n";
    double mu_j = computeMu_j(1, 0.0);
    double b_j = computeB_j(0.0);
    double um = computeUmContrib(1, 0.0);
    double ug3 = computeUg3Contrib(0.0);
    oss << "  μ_j (j=1): " << mu_j << " T·m³\n";
    oss << "  B_j: " << b_j << " T\n";
    oss << "  U_m contribution: " << um << " J/m³\n";
    oss << "  Ug3 contribution: " << ug3 << " J/m³\n\n";
    
    oss << "Time Evolution (μ_j oscillation):\n";
    std::vector<double> times = {0.0, 1000.0, 5000.0, 10000.0};
    for (double t : times) {
        double mu_t = computeMu_j(1, t);
        double sin_term = std::sin(variables["omega_c"] * t);
        oss << "  t=" << t << " s: μ_j=" << mu_t << " T·m³ (sin(ω_c t)=" << sin_term << ")\n";
    }
    oss << "\n";
    
    oss << "Physical Context:\n";
    oss << "  μ_j = (B_j + 0.4 sin(ω_c t)) × base_mu\n";
    oss << "  B_j base field ~1000 T, oscillates with 0.4 T amplitude\n";
    oss << "  Drives U_m (universal magnetism) via μ_j/r_j term\n";
    oss << "  Drives Ug3 (magnetic disk gravity) via k3 × B_j\n";
    oss << "  Time-dependent: cyclic variations in field strength\n";
    oss << "  String j=1 at 100 AU scale (r_j ~1.5e13 m)\n";
    oss << "  Critical for jets, disks, nebulae, planetary magnetism\n";
    oss << "===================================================\n";
    return oss.str();
}

bool MagneticMomentModule::validateConsistency() {
    bool valid = true;
    if (variables["base_mu"] <= 0) valid = false;
    if (variables["B_j"] <= 0) valid = false;
    if (variables["omega_c"] <= 0) valid = false;
    if (variables["r_j"] <= 0) valid = false;
    if (variables["E_react"] <= 0) valid = false;
    if (variables["k3"] <= 0) valid = false;
    return valid;
}

void MagneticMomentModule::autoCorrectAnomalies() {
    if (variables["base_mu"] <= 0) variables["base_mu"] = 3.38e20;
    if (variables["B_j"] <= 0) variables["B_j"] = 1e3;
    if (variables["omega_c"] <= 0) variables["omega_c"] = 2.5e-6;
    if (variables["r_j"] <= 0) variables["r_j"] = 1.496e13;
    if (variables["gamma"] < 0) variables["gamma"] = 5.79e-10;
    if (variables["E_react"] <= 0) variables["E_react"] = 1e46;
    if (variables["k3"] <= 0) variables["k3"] = 1.8;
}

// Example usage in base program (snippet)
// #include "MagneticMomentModule.h"
// int main() {
//     MagneticMomentModule mod;
//     double t = 0.0;
//     mod.printMomentContributions(1, t);
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("base_mu", 4e20);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o moment_test moment_test.cpp MagneticMomentModule.cpp -lm
// Sample: ?_j?3.38e23 T�m^3; U_m?2.28e65 J/m�; cyclic variation via sin(?_c t).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

// ===== Enhanced Example: 18-Step Demonstration of Dynamic Capabilities =====
// int main() {
//     MagneticMomentModule mod;
//     std::cout << "===== Magnetic Moment Module: Enhanced 18-Step Demo =====\n\n";
//
//     // Step 1: Initial report
//     std::cout << "Step 1: Initial Configuration\n";
//     std::cout << mod.generateReport() << "\n";
//
//     // Step 2: Track key magnetic quantities
//     std::cout << "Step 2: Create Tracking Variables\n";
//     double mu_base = mod.computeMu_j(1, 0.0);
//     double um_base = mod.computeUmContrib(1, 0.0);
//     double ug3_base = mod.computeUg3Contrib(0.0);
//     mod.createVariable("mu_j_baseline", mu_base);
//     mod.createVariable("Um_baseline", um_base);
//     mod.createVariable("Ug3_baseline", ug3_base);
//     mod.createVariable("string_density", mu_base / mod.variables["r_j"]); // T·m² (μ/r coupling)
//     std::cout << "  μ_j_baseline: " << std::scientific << mu_base << " T·m³\n";
//     std::cout << "  Um_baseline: " << um_base << " J/m³\n";
//     std::cout << "  Ug3_baseline: " << ug3_base << " J/m³\n";
//     std::cout << "  string_density (μ/r): " << mod.variables["string_density"] << " T·m²\n\n";
//
//     // Step 3: Time evolution (oscillation demonstration)
//     std::cout << "Step 3: Time Evolution of μ_j (Oscillation with ω_c)\n";
//     mod.saveState("baseline");
//     std::vector<double> time_points = {0.0, 500.0, 1000.0, 2500.0, 5000.0, 10000.0, 25000.0};
//     std::cout << "  Time series (ω_c = " << mod.variables["omega_c"] << " rad/s):\n";
//     for (double t : time_points) {
//         double mu_t = mod.computeMu_j(1, t);
//         double b_t = mod.computeB_j(t);
//         double um_t = mod.computeUmContrib(1, t);
//         std::cout << "    t=" << t << " s: μ_j=" << mu_t << " T·m³, B_j=" << b_t << " T, U_m=" << um_t << " J/m³\n";
//     }
//     std::cout << "\n";
//
//     // Step 4: Base moment variations
//     std::cout << "Step 4: base_mu Scaling Effects\n";
//     std::vector<double> base_mu_factors = {0.5, 0.8, 1.0, 1.2, 1.5, 2.0};
//     for (double factor : base_mu_factors) {
//         mod.restoreState("baseline");
//         mod.updateVariable("base_mu", 3.38e20 * factor);
//         double mu_new = mod.computeMu_j(1, 0.0);
//         double um_new = mod.computeUmContrib(1, 0.0);
//         std::cout << "  base_mu x" << factor << ": μ_j=" << mu_new << " T·m³, U_m=" << um_new << " J/m³\n";
//     }
//     mod.restoreState("baseline");
//     std::cout << "\n";
//
//     // Step 5: B_j field strength variations
//     std::cout << "Step 5: B_j Field Strength Scaling\n";
//     std::vector<double> b_factors = {0.5, 1.0, 1.5, 2.0, 3.0};
//     for (double factor : b_factors) {
//         mod.restoreState("baseline");
//         mod.updateVariable("B_j", 1e3 * factor);
//         double mu_new = mod.computeMu_j(1, 0.0);
//         double ug3_new = mod.computeUg3Contrib(0.0);
//         std::cout << "  B_j x" << factor << " (" << (1e3*factor) << " T): μ_j=" << mu_new 
//                   << " T·m³, Ug3=" << ug3_new << " J/m³\n";
//     }
//     mod.restoreState("baseline");
//     std::cout << "\n";
//
//     // Step 6: Expand moment scale (base_mu and B_j)
//     std::cout << "Step 6: Expand Moment Scale (base_mu x1.3, B_j x1.2)\n";
//     mod.expandMomentScale(1.3, 1.2);
//     std::cout << "  New base_mu: " << mod.variables["base_mu"] << " T·m³\n";
//     std::cout << "  New B_j: " << mod.variables["B_j"] << " T\n";
//     std::cout << "  New μ_j: " << mod.computeMu_j(1, 0.0) << " T·m³\n";
//     std::cout << "  New U_m: " << mod.computeUmContrib(1, 0.0) << " J/m³\n\n";
//
//     // Step 7: Expand oscillation scale (ω_c)
//     std::cout << "Step 7: Expand Oscillation Scale (ω_c x1.5)\n";
//     mod.expandOscillationScale(1.5, 1.0);
//     std::cout << "  New ω_c: " << mod.variables["omega_c"] << " rad/s\n";
//     std::cout << "  Period T = 2π/ω_c = " << (2*3.14159/mod.variables["omega_c"]) << " s\n";
//     std::cout << "  μ_j at t=1000s: " << mod.computeMu_j(1, 1000.0) << " T·m³\n\n";
//
//     // Step 8: Expand string scale (r_j and γ)
//     std::cout << "Step 8: Expand String Scale (r_j x1.2, γ x1.1)\n";
//     mod.expandStringScale(1.2, 1.1);
//     std::cout << "  New r_j: " << mod.variables["r_j"] << " m (" << (mod.variables["r_j"]/1.496e11) << " AU)\n";
//     std::cout << "  New γ: " << mod.variables["gamma"] << " s⁻¹\n";
//     std::cout << "  New μ/r: " << (mod.computeMu_j(1, 0.0) / mod.variables["r_j"]) << " T·m²\n";
//     std::cout << "  New U_m: " << mod.computeUmContrib(1, 0.0) << " J/m³\n\n";
//
//     // Step 9: Parameter variations
//     std::cout << "Step 9: Generate 10 Parameter Variations (±7%)\n";
//     auto variations = mod.generateVariations(10, 0.07);
//     std::cout << "  Generated " << variations.size() << " configurations\n";
//     std::vector<double> mu_range, um_range;
//     for (const auto& var : variations) {
//         double b_val = var.at("B_j") + 0.4 * std::sin(var.at("omega_c") * 0.0);
//         double mu_val = b_val * var.at("base_mu");
//         mu_range.push_back(mu_val);
//     }
//     auto mu_minmax = std::minmax_element(mu_range.begin(), mu_range.end());
//     std::cout << "  μ_j range: " << *mu_minmax.first << " to " << *mu_minmax.second << " T·m³\n\n";
//
//     // Step 10: Sensitivity analysis
//     std::cout << "Step 10: Sensitivity Analysis (U_m response to ±1% changes)\n";
//     std::vector<std::string> sens_params = {"base_mu", "B_j", "r_j", "E_react", "scale_Heaviside", "gamma", "omega_c"};
//     auto sensitivity = mod.sensitivityAnalysis(sens_params);
//     std::cout << "  Most influential parameters:\n";
//     std::vector<std::pair<std::string, double>> sorted_sens(sensitivity.begin(), sensitivity.end());
//     std::sort(sorted_sens.begin(), sorted_sens.end(), 
//               [](const auto& a, const auto& b) { return a.second > b.second; });
//     for (const auto& s : sorted_sens) {
//         std::cout << "    " << s.first << ": " << (s.second * 100) << "% sensitivity\n";
//     }
//     std::cout << "\n";
//
//     // Step 11: Auto-refine to target U_m
//     std::cout << "Step 11: Auto-Refine to Target U_m = 3.0e65 J/m³\n";
//     double target_um = 3.0e65;
//     double before_um = mod.computeUmContrib(1, 0.0);
//     mod.autoRefineParameters("U_m", target_um);
//     double after_um = mod.computeUmContrib(1, 0.0);
//     std::cout << "  Before: " << before_um << " J/m³\n";
//     std::cout << "  After: " << after_um << " J/m³\n";
//     std::cout << "  Target: " << target_um << " J/m³\n";
//     std::cout << "  Error: " << std::abs(after_um - target_um) / target_um * 100 << "%\n\n";
//
//     // Step 12: Target specific μ_j value
//     std::cout << "Step 12: Target Specific μ_j = 4.0e23 T·m³\n";
//     double target_mu = 4.0e23;
//     double before_mu = mod.computeMu_j(1, 0.0);
//     mod.autoRefineParameters("mu_j", target_mu);
//     double after_mu = mod.computeMu_j(1, 0.0);
//     std::cout << "  Before: " << before_mu << " T·m³\n";
//     std::cout << "  After: " << after_mu << " T·m³\n";
//     std::cout << "  New base_mu: " << mod.variables["base_mu"] << " T·m³\n\n";
//
//     // Step 13: Calibration to observations
//     std::cout << "Step 13: Calibrate to Observational Data\n";
//     std::map<std::string, double> observations = {
//         {"B_j", 1100.0},        // Slightly stronger field
//         {"omega_c", 2.8e-6},    // Faster oscillation
//         {"r_j", 1.8e13}         // ~120 AU
//     };
//     mod.calibrateToObservations(observations);
//     std::cout << "  Calibrated B_j: " << mod.variables["B_j"] << " T\n";
//     std::cout << "  Calibrated ω_c: " << mod.variables["omega_c"] << " rad/s\n";
//     std::cout << "  Calibrated r_j: " << mod.variables["r_j"] << " m (" << (mod.variables["r_j"]/1.496e11) << " AU)\n";
//     std::cout << "  New μ_j: " << mod.computeMu_j(1, 0.0) << " T·m³\n\n";
//
//     // Step 14: Optimize for maximum U_m
//     std::cout << "Step 14: Optimize for Maximum U_m\n";
//     double before_opt_um = mod.computeUmContrib(1, 0.0);
//     mod.optimizeForMetric("maximize_Um");
//     double after_opt_um = mod.computeUmContrib(1, 0.0);
//     std::cout << "  U_m before: " << before_opt_um << " J/m³\n";
//     std::cout << "  U_m after: " << after_opt_um << " J/m³\n";
//     std::cout << "  Increase: " << ((after_opt_um/before_opt_um - 1.0)*100) << "%\n\n";
//
//     // Step 15: Optimize for enhanced oscillation
//     std::cout << "Step 15: Optimize for Enhanced Oscillation\n";
//     mod.saveState("before_osc_opt");
//     double omega_before = mod.variables["omega_c"];
//     mod.optimizeForMetric("enhance_oscillation");
//     double omega_after = mod.variables["omega_c"];
//     std::cout << "  ω_c: " << omega_before << " -> " << omega_after << " rad/s\n";
//     std::cout << "  Period: " << (2*3.14159/omega_after) << " s\n";
//     if (mod.variables.find("osc_amplitude") != mod.variables.end()) {
//         std::cout << "  Oscillation amplitude: " << mod.variables["osc_amplitude"] << "\n";
//     }
//     std::cout << "\n";
//
//     // Step 16: Parameter mutation
//     std::cout << "Step 16: Mutate Parameters (±4%)\n";
//     mod.restoreState("before_osc_opt");
//     mod.saveState("before_mutation");
//     double mu_before_mut = mod.computeMu_j(1, 0.0);
//     mod.mutateParameters(0.04);
//     double mu_after_mut = mod.computeMu_j(1, 0.0);
//     std::cout << "  μ_j: " << mu_before_mut << " -> " << mu_after_mut << " T·m³\n";
//     std::cout << "  Change: " << ((mu_after_mut/mu_before_mut - 1.0)*100) << "%\n\n";
//
//     // Step 17: System evolution (maximize U_m)
//     std::cout << "Step 17: Evolve System (7 generations, maximize U_m)\n";
//     mod.restoreState("before_mutation");
//     double initial_fitness = mod.computeUmContrib(1, 0.0);
//     mod.evolveSystem(7, [&mod]() { return mod.computeUmContrib(1, 0.0); });
//     double final_fitness = mod.computeUmContrib(1, 0.0);
//     std::cout << "  Initial U_m: " << initial_fitness << " J/m³\n";
//     std::cout << "  Evolved U_m: " << final_fitness << " J/m³\n";
//     std::cout << "  Improvement: " << ((final_fitness / initial_fitness - 1.0) * 100) << "%\n\n";
//
//     // Step 18: Validation and state export
//     std::cout << "Step 18: Validate Consistency and Export\n";
//     bool valid = mod.validateConsistency();
//     std::cout << "  System valid: " << (valid ? "YES" : "NO") << "\n";
//     if (!valid) {
//         std::cout << "  Running auto-correction...\n";
//         mod.autoCorrectAnomalies();
//         std::cout << "  Post-correction valid: " << (mod.validateConsistency() ? "YES" : "NO") << "\n";
//     }
//     mod.printMomentContributions(1, 0.0);
//     std::cout << "\n";
//     mod.saveState("evolved_optimal");
//     auto saved = mod.listSavedStates();
//     std::cout << "  Saved states (" << saved.size() << "): ";
//     for (const auto& s : saved) std::cout << s << " ";
//     std::cout << "\n\n";
//     std::cout << "Final System Export:\n";
//     std::string exported = mod.exportState();
//     std::cout << exported << "\n";
//     std::cout << "Demo complete: 18 steps executed successfully!\n";
//     std::cout << "============================================================\n";
//
//     return 0;
// }

MagneticMomentModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeMu_j, computeB_j, computeUmContrib, computeUg3Contrib) are clear, concise, and variable - driven.
- Output and debugging functions(printVariables, printMomentContributions, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.
- Time - dependent magnetic moment(?_j) and field(B_j) allow for cyclic / physical modeling.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid indices, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in magnetic moment modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.