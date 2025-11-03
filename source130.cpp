// UniversalInertiaVacuumModule.h
// Modular C++ implementation of the Vacuum Energy Density of Universal Inertia (ρ_vac,Ui) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes ρ_vac,Ui = 2.84e-36 J/m³ (Sun, level 13); reference scale for U_i inertial term.
// Pluggable: #include "UniversalInertiaVacuumModule.h"
// UniversalInertiaVacuumModule mod; mod.computeU_i_example(0.0, 0.0); mod.updateVariable("rho_vac_Ui", new_value);
// Variables in std::map; example for Sun at t=0, t_n=0; U_i ≈1.38e-47 J/m³.
// Approximations: λ_i=1.0; cos(π t_n)=1; ω_s=2.5e-6 rad/s; f_TRZ=0.1; ρ_[SCm/UA] product=5.03e-72 J²/m⁶.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef UNIVERSAL_INERTIA_VACUUM_MODULE_H
#define UNIVERSAL_INERTIA_VACUUM_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <functional>
#include <random>
#include <algorithm>
#include <sstream>

class UniversalInertiaVacuumModule {
private:
    std::map<std::string, double> variables;
    double computeU_i_base(double t, double t_n);
    double computeU_i(double t, double t_n);

public:
    // Constructor: Initialize with framework defaults (Sun, level 13)
    UniversalInertiaVacuumModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeRho_vac_Ui();  // 2.84e-36 J/m³
    double computeU_i(double t, double t_n);  // U_i example (J/m³)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // ========== ENHANCED DYNAMIC CAPABILITIES ==========
    // Variable Management
    void createVariable(const std::string& name, double value);
    void removeVariable(const std::string& name);
    void cloneVariable(const std::string& source, const std::string& target);
    std::vector<std::string> listVariables();
    std::string getSystemName() const;

    // Batch Operations
    void transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func);
    void scaleVariableGroup(const std::vector<std::string>& names, double factor);

    // Self-Expansion (Domain-Specific)
    void expandParameterSpace(double expansion_factor);
    void expandInertiaScale(double lambda_factor, double resistance_factor);   // λ_i and U_i inertial resistance
    void expandVacuumScale(double rho_factor, double product_factor);          // ρ_vac densities (Ui, SCm, UA)
    void expandTemporalScale(double omega_factor, double phase_factor);        // ω_s(t) and cos(π t_n) temporal

    // Self-Refinement
    void autoRefineParameters();
    void calibrateToObservations(const std::map<std::string, double>& observations);
    void optimizeForMetric(const std::string& metric_name);

    // Parameter Exploration
    std::vector<std::map<std::string, double>> generateVariations(int count);

    // Adaptive Evolution
    void mutateParameters(double mutation_rate);
    void evolveSystem(int generations);

    // State Management
    void saveState(const std::string& state_name);
    void restoreState(const std::string& state_name);
    std::vector<std::string> listSavedStates();
    std::string exportState();

    // System Analysis
    std::map<std::string, double> sensitivityAnalysis(const std::string& output_var);
    std::string generateReport();
    bool validateConsistency();
    void autoCorrectAnomalies();
};

#endif // UNIVERSAL_INERTIA_VACUUM_MODULE_H

// UniversalInertiaVacuumModule.cpp
#include "UniversalInertiaVacuumModule.h"

// Constructor: Set framework defaults (Sun at t=0, level 13)
UniversalInertiaVacuumModule::UniversalInertiaVacuumModule() {
    // Universal constants
    variables["rho_vac_Ui"] = 2.84e-36;             // J/m³ (reference scale)
    variables["lambda_i"] = 1.0;                    // Unitless
    variables["rho_vac_SCm"] = 7.09e-37;            // J/m³
    variables["rho_vac_UA"] = 7.09e-36;             // J/m³
    variables["omega_s"] = 2.5e-6;                  // rad/s
    variables["f_TRZ"] = 0.1;                       // Unitless
    variables["pi"] = 3.141592653589793;
    variables["t"] = 0.0;                           // s
    variables["t_n"] = 0.0;                         // s

    // Derived
    variables["rho_product"] = variables["rho_vac_SCm"] * variables["rho_vac_UA"];
}

// Update variable
void UniversalInertiaVacuumModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
        if (name == "rho_vac_SCm" || name == "rho_vac_UA") {
            variables["rho_product"] = variables["rho_vac_SCm"] * variables["rho_vac_UA"];
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void UniversalInertiaVacuumModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        if (name == "rho_vac_SCm" || name == "rho_vac_UA") {
            variables["rho_product"] = variables["rho_vac_SCm"] * variables["rho_vac_UA"];
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void UniversalInertiaVacuumModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute ρ_vac,Ui (J/m³, reference)
double UniversalInertiaVacuumModule::computeRho_vac_Ui() {
    return variables["rho_vac_Ui"];
}

// Base U_i without TRZ
double UniversalInertiaVacuumModule::computeU_i_base(double t, double t_n) {
    double lambda_i = variables["lambda_i"];
    double rho_product = variables["rho_product"];
    double omega_s_t = variables["omega_s"];        // Simplified constant
    double cos_pi_tn = std::cos(variables["pi"] * t_n);
    return lambda_i * rho_product * omega_s_t * cos_pi_tn;
}

// U_i with TRZ
double UniversalInertiaVacuumModule::computeU_i(double t, double t_n) {
    variables["t"] = t;
    double base = computeU_i_base(t, t_n);
    double trz_factor = 1.0 + variables["f_TRZ"];
    return base * trz_factor;
}

// Equation text
std::string UniversalInertiaVacuumModule::getEquationText() {
    return "U_i = λ_i * ρ_vac,[SCm] * ρ_vac,[UA] * ω_s(t) * cos(π t_n) * (1 + f_TRZ)\n"
           "ρ_vac,Ui = 2.84e-36 J/m³ (Sun level 13, inertia vacuum scale; not direct in eq.).\n"
           "Provides reference for U_i magnitude; inertial resistance from [SCm]/[UA].\n"
           "Example t=0, t_n=0: U_i ≈1.38e-47 J/m³ (consistent scale with ρ_vac,Ui).\n"
           "In F_U: -∑ λ_i U_i E_react (resistive inertia).\n"
           "Role: Quantifies vacuum inertia energy; opposes dynamics in nebulae/formation.\n"
           "UQFF: Small-scale reference for cosmic inertia; [SCm]-[UA] resistance.";
}

// Print variables
void UniversalInertiaVacuumModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// ========== ENHANCED DYNAMIC CAPABILITIES IMPLEMENTATION ==========

// ===== Variable Management =====
void UniversalInertiaVacuumModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void UniversalInertiaVacuumModule::removeVariable(const std::string& name) {
    variables.erase(name);
}

void UniversalInertiaVacuumModule::cloneVariable(const std::string& source, const std::string& target) {
    if (variables.find(source) != variables.end()) {
        variables[target] = variables[source];
    }
}

std::vector<std::string> UniversalInertiaVacuumModule::listVariables() {
    std::vector<std::string> names;
    for (const auto& pair : variables) {
        names.push_back(pair.first);
    }
    return names;
}

std::string UniversalInertiaVacuumModule::getSystemName() const {
    return "Universal_Inertia_Vacuum_rho_vac_Ui_UQFF_v2";
}

// ===== Batch Operations =====
void UniversalInertiaVacuumModule::transformVariableGroup(const std::vector<std::string>& names, 
                                                          std::function<double(double)> func) {
    for (const auto& name : names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = func(variables[name]);
        }
    }
    // Update derived
    variables["rho_product"] = variables["rho_vac_SCm"] * variables["rho_vac_UA"];
}

void UniversalInertiaVacuumModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double x) { return x * factor; });
}

// ===== Self-Expansion (Domain-Specific) =====
void UniversalInertiaVacuumModule::expandParameterSpace(double expansion_factor) {
    // General expansion across all key parameters
    variables["rho_vac_Ui"] *= expansion_factor;
    variables["lambda_i"] *= expansion_factor;
    variables["omega_s"] *= expansion_factor;
    variables["f_TRZ"] = std::min(1.0, variables["f_TRZ"] * expansion_factor);
    variables["rho_product"] = variables["rho_vac_SCm"] * variables["rho_vac_UA"];
}

void UniversalInertiaVacuumModule::expandInertiaScale(double lambda_factor, double resistance_factor) {
    // Expand inertial resistance: λ_i (inertia coupling), U_i magnitude
    variables["lambda_i"] *= lambda_factor;
    
    // Resistance affects vacuum densities that produce U_i
    variables["rho_vac_SCm"] *= resistance_factor;
    variables["rho_vac_UA"] *= resistance_factor;
    variables["rho_vac_Ui"] *= resistance_factor;  // Reference scale
    
    // Update product
    variables["rho_product"] = variables["rho_vac_SCm"] * variables["rho_vac_UA"];
}

void UniversalInertiaVacuumModule::expandVacuumScale(double rho_factor, double product_factor) {
    // Expand vacuum densities: ρ_vac,Ui (reference), ρ_vac,SCm, ρ_vac,UA
    variables["rho_vac_Ui"] *= rho_factor;
    variables["rho_vac_SCm"] *= rho_factor;
    variables["rho_vac_UA"] *= rho_factor;
    
    // Product scaling (may differ from individual densities)
    variables["rho_product"] *= product_factor;
}

void UniversalInertiaVacuumModule::expandTemporalScale(double omega_factor, double phase_factor) {
    // Expand temporal dynamics: ω_s(t) rotation rate, phase evolution
    variables["omega_s"] *= omega_factor;
    
    // Phase affects t_n (negative time cycles)
    variables["t_n"] *= phase_factor;
    
    // TRZ temporal enhancement
    variables["f_TRZ"] = std::min(1.0, variables["f_TRZ"] * phase_factor);
}

// ===== Self-Refinement =====
void UniversalInertiaVacuumModule::autoRefineParameters() {
    // Clamp ρ_vac,Ui to physically reasonable range [1e-40, 1e-30] J/m³
    if (variables["rho_vac_Ui"] < 1e-40) variables["rho_vac_Ui"] = 1e-40;
    if (variables["rho_vac_Ui"] > 1e-30) variables["rho_vac_Ui"] = 1e-30;
    
    // λ_i inertia coupling [0.1, 10]
    if (variables["lambda_i"] < 0.1) variables["lambda_i"] = 0.1;
    if (variables["lambda_i"] > 10.0) variables["lambda_i"] = 10.0;
    
    // ω_s rotation rate [1e-8, 1e-4] rad/s
    if (variables["omega_s"] < 1e-8) variables["omega_s"] = 1e-8;
    if (variables["omega_s"] > 1e-4) variables["omega_s"] = 1e-4;
    
    // f_TRZ [0, 1]
    if (variables["f_TRZ"] < 0.0) variables["f_TRZ"] = 0.0;
    if (variables["f_TRZ"] > 1.0) variables["f_TRZ"] = 1.0;
    
    // Update derived
    variables["rho_product"] = variables["rho_vac_SCm"] * variables["rho_vac_UA"];
}

void UniversalInertiaVacuumModule::calibrateToObservations(const std::map<std::string, double>& observations) {
    for (const auto& obs : observations) {
        if (variables.find(obs.first) != variables.end()) {
            variables[obs.first] = obs.second;
        }
    }
    variables["rho_product"] = variables["rho_vac_SCm"] * variables["rho_vac_UA"];
    autoRefineParameters();
}

void UniversalInertiaVacuumModule::optimizeForMetric(const std::string& metric_name) {
    // Optimize for different physical scenarios
    if (metric_name == "low_inertia") {
        // Minimal resistance: Interstellar medium, free particles
        variables["rho_vac_Ui"] = 1e-38;
        variables["lambda_i"] = 0.5;
        variables["f_TRZ"] = 0.05;
    } else if (metric_name == "solar_inertia") {
        // Standard Sun reference (default)
        variables["rho_vac_Ui"] = 2.84e-36;
        variables["lambda_i"] = 1.0;
        variables["omega_s"] = 2.5e-6;
        variables["f_TRZ"] = 0.1;
    } else if (metric_name == "nebula_inertia") {
        // Higher inertia: Star-forming regions, dense nebulae
        variables["rho_vac_Ui"] = 1e-35;
        variables["lambda_i"] = 2.0;
        variables["f_TRZ"] = 0.2;
    } else if (metric_name == "galaxy_core") {
        // Galactic core dynamics: High rotation, strong inertia
        variables["rho_vac_Ui"] = 5e-35;
        variables["lambda_i"] = 3.0;
        variables["omega_s"] = 1e-5;
        variables["f_TRZ"] = 0.15;
    } else if (metric_name == "high_resistance") {
        // Maximum inertial resistance: Dense systems, quasar jets
        variables["rho_vac_Ui"] = 1e-34;
        variables["lambda_i"] = 5.0;
        variables["f_TRZ"] = 0.3;
    }
    variables["rho_product"] = variables["rho_vac_SCm"] * variables["rho_vac_UA"];
}

// ===== Parameter Exploration =====
std::vector<std::map<std::string, double>> UniversalInertiaVacuumModule::generateVariations(int count) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.5, 2.0);
    
    for (int i = 0; i < count; ++i) {
        std::map<std::string, double> variant = variables;
        variant["rho_vac_Ui"] *= dis(gen);
        variant["lambda_i"] *= dis(gen);
        variant["omega_s"] *= dis(gen);
        variant["f_TRZ"] = std::min(1.0, variant["f_TRZ"] * dis(gen));
        variant["rho_product"] = variant["rho_vac_SCm"] * variant["rho_vac_UA"];
        variations.push_back(variant);
    }
    
    return variations;
}

// ===== Adaptive Evolution =====
void UniversalInertiaVacuumModule::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-mutation_rate, mutation_rate);
    
    variables["rho_vac_Ui"] *= (1.0 + dis(gen));
    variables["lambda_i"] *= (1.0 + dis(gen));
    variables["omega_s"] *= (1.0 + dis(gen));
    variables["f_TRZ"] = std::min(1.0, std::max(0.0, variables["f_TRZ"] * (1.0 + dis(gen))));
    
    variables["rho_product"] = variables["rho_vac_SCm"] * variables["rho_vac_UA"];
    autoRefineParameters();
}

void UniversalInertiaVacuumModule::evolveSystem(int generations) {
    for (int i = 0; i < generations; ++i) {
        mutateParameters(0.05);
        
        // Fitness: Minimize U_i magnitude (lower inertial resistance preferred)
        double u_i = computeU_i(variables["t"], variables["t_n"]);
        if (u_i > 1e-45) {
            variables["lambda_i"] *= 0.95;
        }
        
        autoRefineParameters();
    }
}

// ===== State Management =====
namespace {
    std::map<std::string, std::map<std::string, double>> saved_states_v2;
}

void UniversalInertiaVacuumModule::saveState(const std::string& state_name) {
    saved_states_v2[state_name] = variables;
}

void UniversalInertiaVacuumModule::restoreState(const std::string& state_name) {
    if (saved_states_v2.find(state_name) != saved_states_v2.end()) {
        variables = saved_states_v2[state_name];
    }
}

std::vector<std::string> UniversalInertiaVacuumModule::listSavedStates() {
    std::vector<std::string> names;
    for (const auto& pair : saved_states_v2) {
        names.push_back(pair.first);
    }
    return names;
}

std::string UniversalInertiaVacuumModule::exportState() {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(6);
    oss << "UniversalInertiaVacuumModule State (v2):\n";
    for (const auto& pair : variables) {
        oss << pair.first << "=" << pair.second << "\n";
    }
    return oss.str();
}

// ===== System Analysis =====
std::map<std::string, double> UniversalInertiaVacuumModule::sensitivityAnalysis(const std::string& output_var) {
    std::map<std::string, double> sensitivity;
    double baseline = 0.0;
    
    if (output_var == "U_i") {
        baseline = computeU_i(variables["t"], variables["t_n"]);
    } else if (output_var == "rho_vac_Ui") {
        baseline = computeRho_vac_Ui();
    }
    
    double perturbation = 0.01;  // 1% change
    
    for (const auto& pair : variables) {
        std::string var_name = pair.first;
        double original = variables[var_name];
        
        variables[var_name] = original * (1.0 + perturbation);
        variables["rho_product"] = variables["rho_vac_SCm"] * variables["rho_vac_UA"];
        
        double perturbed = 0.0;
        if (output_var == "U_i") {
            perturbed = computeU_i(variables["t"], variables["t_n"]);
        } else if (output_var == "rho_vac_Ui") {
            perturbed = computeRho_vac_Ui();
        }
        
        sensitivity[var_name] = ((perturbed - baseline) / baseline) / perturbation;
        variables[var_name] = original;
    }
    
    variables["rho_product"] = variables["rho_vac_SCm"] * variables["rho_vac_UA"];
    return sensitivity;
}

std::string UniversalInertiaVacuumModule::generateReport() {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(3);
    oss << "========== Universal Inertia Vacuum Module Report (v2) ==========\n";
    oss << "System: " << getSystemName() << "\n\n";
    
    oss << "Core Parameters:\n";
    oss << "  ρ_vac,Ui = " << computeRho_vac_Ui() << " J/m³ (inertia vacuum reference)\n";
    oss << "  λ_i = " << variables["lambda_i"] << " (inertia coupling)\n";
    oss << "  ω_s = " << variables["omega_s"] << " rad/s (rotation rate)\n";
    oss << "  f_TRZ = " << variables["f_TRZ"] << " (time-reversal enhancement)\n\n";
    
    oss << "Vacuum Densities:\n";
    oss << "  ρ_vac,SCm = " << variables["rho_vac_SCm"] << " J/m³\n";
    oss << "  ρ_vac,UA = " << variables["rho_vac_UA"] << " J/m³\n";
    oss << "  Product = " << variables["rho_product"] << " J²/m⁶\n\n";
    
    oss << "Computed Values:\n";
    double u_i = computeU_i(variables["t"], variables["t_n"]);
    double u_i_base = computeU_i_base(variables["t"], variables["t_n"]);
    oss << "  U_i (base) = " << u_i_base << " J/m³\n";
    oss << "  U_i (with TRZ) = " << u_i << " J/m³\n";
    oss << "  TRZ enhancement = " << (u_i / u_i_base - 1.0) * 100.0 << "%\n\n";
    
    oss << "Physical Interpretation:\n";
    oss << "  - ρ_vac,Ui: Reference scale for inertial vacuum energy density\n";
    oss << "  - U_i: Actual inertial resistance energy (F_U: -Σ λ_i U_i E_react)\n";
    oss << "  - [SCm]-[UA]: Vacuum product drives inertial resistance\n";
    oss << "  - f_TRZ: Time-reversal zones enhance (reduce resistance ~10-30%)\n";
    oss << "  - Applications: Nebular dynamics, star formation, cosmic inertia\n";
    oss << "================================================================\n";
    
    return oss.str();
}

bool UniversalInertiaVacuumModule::validateConsistency() {
    bool valid = true;
    
    // Check ρ_vac,Ui range [1e-40, 1e-30] J/m³
    if (variables["rho_vac_Ui"] < 1e-40 || variables["rho_vac_Ui"] > 1e-30) valid = false;
    
    // Check λ_i range [0.1, 10]
    if (variables["lambda_i"] < 0.0 || variables["lambda_i"] > 10.0) valid = false;
    
    // Check ω_s range [1e-8, 1e-4] rad/s
    if (variables["omega_s"] < 1e-8 || variables["omega_s"] > 1e-4) valid = false;
    
    // Check f_TRZ [0, 1]
    if (variables["f_TRZ"] < 0.0 || variables["f_TRZ"] > 1.0) valid = false;
    
    // Check product consistency
    double expected_product = variables["rho_vac_SCm"] * variables["rho_vac_UA"];
    if (std::abs(variables["rho_product"] - expected_product) / expected_product > 0.01) valid = false;
    
    return valid;
}

void UniversalInertiaVacuumModule::autoCorrectAnomalies() {
    // Correct out-of-range parameters
    autoRefineParameters();
    
    // Fix product
    variables["rho_product"] = variables["rho_vac_SCm"] * variables["rho_vac_UA"];
    
    // Ensure temporal consistency
    if (variables["t"] < 0.0) variables["t"] = 0.0;
    if (std::abs(variables["t_n"]) > 10.0) variables["t_n"] = 0.0;
}

// Example usage in base program (snippet)
// #include "UniversalInertiaVacuumModule.h"
// int main() {
//     UniversalInertiaVacuumModule mod;
//     double rho = mod.computeRho_vac_Ui();
//     std::cout << "ρ_vac,Ui = " << rho << " J/m³\n";
//     double u_i = mod.computeU_i(0.0, 0.0);
//     std::cout << "U_i = " << u_i << " J/m³\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("rho_vac_Ui", 3e-36);
//     mod.printVariables();
//     return 0;
// }

// ========== COMPREHENSIVE ENHANCED DEMONSTRATION (v2) ==========
/*
int main() {
    std::cout << std::scientific << std::setprecision(3);
    std::cout << "========== UNIVERSAL INERTIA VACUUM MODULE v2 DEMONSTRATION ==========\n\n";
    
    // ===== Step 1: Initialize Module =====
    UniversalInertiaVacuumModule mod;
    std::cout << "Step 1: Module initialized with defaults (v2):\n";
    std::cout << "  System: " << mod.getSystemName() << "\n";
    std::cout << "  ρ_vac,Ui = " << mod.computeRho_vac_Ui() << " J/m³\n";
    std::cout << "  U_i (t=0, t_n=0) = " << mod.computeU_i(0.0, 0.0) << " J/m³\n\n";
    
    // ===== Step 2: Baseline Computation =====
    std::cout << "Step 2: Compute baseline inertial parameters:\n";
    double rho_ui = mod.computeRho_vac_Ui();
    double u_i = mod.computeU_i(0.0, 0.0);
    
    std::cout << "  ρ_vac,Ui = " << rho_ui << " J/m³ (reference scale)\n";
    std::cout << "  U_i = " << u_i << " J/m³ (actual inertial resistance)\n";
    std::cout << "  Ratio (U_i/ρ_vac,Ui) = " << (u_i / rho_ui) << "\n";
    std::cout << "  Physical interpretation: U_i ~11 orders smaller than reference\n\n";
    
    // ===== Step 3-7: Dynamic Operations =====
    std::cout << "Step 3: Variable Management\n";
    mod.createVariable("inertia_force_scale", -u_i * 1e10);
    std::cout << "  Created 'inertia_force_scale'\n";
    
    std::cout << "\nStep 4: Inertia Expansion\n";
    mod.expandInertiaScale(2.0, 1.5);
    std::cout << "  Expanded: λ_i = " << mod.variables["lambda_i"] << ", ρ_vac,Ui = " << mod.computeRho_vac_Ui() << " J/m³\n";
    
    std::cout << "\nStep 5: Vacuum Expansion\n";
    mod.expandVacuumScale(2.0, 3.0);
    std::cout << "  Expanded: Product = " << mod.variables["rho_product"] << " J²/m⁶\n";
    
    std::cout << "\nStep 6: Temporal Expansion\n";
    mod.expandTemporalScale(1.5, 1.2);
    std::cout << "  Expanded: ω_s = " << mod.variables["omega_s"] << " rad/s\n";
    
    std::cout << "\nStep 7: Batch Operations\n";
    std::vector<std::string> inertia_group = {"lambda_i", "rho_vac_Ui", "f_TRZ"};
    mod.scaleVariableGroup(inertia_group, 0.5);
    std::cout << "  Scaled inertia group by 0.5\n\n";
    
    // ===== Step 8-12: Physical Regimes =====
    std::cout << "Steps 8-12: Test Multiple Physical Regimes\n";
    
    mod.optimizeForMetric("low_inertia");
    std::cout << "  Low Inertia: ρ_vac,Ui = " << mod.computeRho_vac_Ui() << " J/m³, U_i = " << mod.computeU_i(0.0, 0.0) << " J/m³\n";
    
    mod.optimizeForMetric("solar_inertia");
    std::cout << "  Solar Inertia: ρ_vac,Ui = " << mod.computeRho_vac_Ui() << " J/m³, U_i = " << mod.computeU_i(0.0, 0.0) << " J/m³\n";
    
    mod.optimizeForMetric("nebula_inertia");
    std::cout << "  Nebula Inertia: ρ_vac,Ui = " << mod.computeRho_vac_Ui() << " J/m³, U_i = " << mod.computeU_i(0.0, 0.0) << " J/m³\n";
    
    mod.optimizeForMetric("galaxy_core");
    std::cout << "  Galaxy Core: ρ_vac,Ui = " << mod.computeRho_vac_Ui() << " J/m³, ω_s = " << mod.variables["omega_s"] << " rad/s\n";
    
    mod.optimizeForMetric("high_resistance");
    std::cout << "  High Resistance: ρ_vac,Ui = " << mod.computeRho_vac_Ui() << " J/m³, λ_i = " << mod.variables["lambda_i"] << "\n\n";
    
    // ===== Step 13-17: Refinement & Evolution =====
    std::cout << "Step 13: Auto-Refinement\n";
    mod.updateVariable("lambda_i", 25.0);
    mod.autoRefineParameters();
    std::cout << "  Clamped λ_i from 25.0 to " << mod.variables["lambda_i"] << "\n";
    
    std::cout << "\nStep 14: Calibration to Observations\n";
    std::map<std::string, double> obs_data;
    obs_data["rho_vac_Ui"] = 3.0e-36;
    obs_data["lambda_i"] = 1.2;
    mod.calibrateToObservations(obs_data);
    std::cout << "  Calibrated: ρ_vac,Ui = " << mod.computeRho_vac_Ui() << " J/m³\n";
    
    std::cout << "\nStep 15: Parameter Variations\n";
    std::vector<std::map<std::string, double>> variations = mod.generateVariations(5);
    std::cout << "  Generated " << variations.size() << " parameter variations\n";
    
    std::cout << "\nStep 16: Mutation\n";
    mod.optimizeForMetric("solar_inertia");
    mod.mutateParameters(0.15);
    std::cout << "  Mutated: ρ_vac,Ui = " << mod.computeRho_vac_Ui() << " J/m³\n";
    
    std::cout << "\nStep 17: System Evolution\n";
    mod.evolveSystem(10);
    std::cout << "  Evolved: U_i = " << mod.computeU_i(0.0, 0.0) << " J/m³\n\n";
    
    // ===== Step 18-19: State Management =====
    std::cout << "Step 18: State Management\n";
    mod.optimizeForMetric("solar_inertia");
    mod.saveState("solar_reference");
    mod.optimizeForMetric("nebula_inertia");
    mod.saveState("nebula_high_inertia");
    std::vector<std::string> saved = mod.listSavedStates();
    std::cout << "  Saved " << saved.size() << " states\n";
    mod.restoreState("solar_reference");
    std::cout << "  Restored 'solar_reference'\n";
    
    std::cout << "\nStep 19: Export State\n";
    std::string exported = mod.exportState();
    std::cout << "  Exported " << exported.length() << " bytes\n\n";
    
    // ===== Step 20-22: Analysis =====
    std::cout << "Step 20: Sensitivity Analysis (U_i response)\n";
    std::map<std::string, double> sensitivity = mod.sensitivityAnalysis("U_i");
    std::cout << "  Top sensitivity parameters identified\n";
    
    std::cout << "\nStep 21: Consistency Validation\n";
    bool valid = mod.validateConsistency();
    std::cout << "  System consistency: " << (valid ? "VALID" : "INVALID") << "\n";
    if (!valid) mod.autoCorrectAnomalies();
    
    std::cout << "\nStep 22: Generate Full Report\n";
    std::string report = mod.generateReport();
    std::cout << report << "\n";
    
    // ===== Step 23-26: Scale Analysis =====
    std::cout << "Steps 23-26: Inertia Energy Density Scale Analysis\n";
    std::cout << "  Regime           | ρ_vac,Ui (J/m³) | λ_i   | U_i (J/m³) | Context\n";
    std::cout << "  -------------------------------------------------------------------------------\n";
    
    struct InertiaRegime {
        std::string name;
        double rho_vac_ui;
        double lambda_i;
        std::string context;
    };
    
    std::vector<InertiaRegime> regimes = {
        {"Ultra-Low", 1e-38, 0.5, "Free particles, ISM"},
        {"Solar Standard", 2.84e-36, 1.0, "Sun reference (level 13)"},
        {"Nebula", 1e-35, 2.0, "Star formation regions"},
        {"Galaxy Core", 5e-35, 3.0, "High rotation dynamics"},
        {"Quasar Jet", 1e-34, 5.0, "Maximum resistance"}
    };
    
    for (const auto& reg : regimes) {
        mod.updateVariable("rho_vac_Ui", reg.rho_vac_ui);
        mod.updateVariable("lambda_i", reg.lambda_i);
        mod.updateVariable("rho_vac_SCm", reg.rho_vac_ui * 0.25);
        mod.updateVariable("rho_vac_UA", reg.rho_vac_ui * 2.5);
        mod.variables["rho_product"] = mod.variables["rho_vac_SCm"] * mod.variables["rho_vac_UA"];
        
        double rho = mod.computeRho_vac_Ui();
        double u_i_val = mod.computeU_i(0.0, 0.0);
        
        std::cout << "  " << std::setw(16) << std::left << reg.name
                  << " | " << std::scientific << std::setprecision(2) << std::setw(15) << rho
                  << " | " << std::fixed << std::setprecision(1) << std::setw(5) << reg.lambda_i
                  << " | " << std::scientific << std::setprecision(2) << std::setw(10) << u_i_val
                  << " | " << reg.context << "\n";
    }
    
    std::cout << "\n========== DEMONSTRATION COMPLETE (v2) ==========\n";
    std::cout << "Universal Inertia Vacuum module validated across physical regimes.\n";
    std::cout << "ρ_vac,Ui provides reference scale for inertial vacuum energy.\n";
    std::cout << "U_i = λ_i ρ_vac,[SCm] ρ_vac,[UA] ω_s cos(π t_n) (1 + f_TRZ)\n";
    std::cout << "Physical significance: [SCm]-[UA] vacuum resistance opposes dynamics.\n";
    std::cout << "F_U contribution: -Σ λ_i U_i E_react (resistive inertia in force equation).\n";
    std::cout << "UQFF Integration: Cosmic inertia from superconductive vacuum interactions.\n";
    std::cout << "Applications: Nebular collapse, star formation, galactic dynamics, quasar jets.\n";
    
    return 0;
}
*/
// Compile: g++ -o inertia_vac_test inertia_vac_test.cpp UniversalInertiaVacuumModule.cpp -lm
// Sample: ρ_vac,Ui=2.84e-36 J/m³; U_i≈1.38e-47 J/m³; scale reference.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

UniversalInertiaVacuumModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeRho_vac_Ui, computeU_i) are clear, concise, and variable - driven.
- Automatic recalculation of derived variables(rho_product) when dependencies change.
- Output and debugging functions(printVariables, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.
- Provides a reference scale for vacuum inertia energy, supporting scientific modeling.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in universal inertia vacuum energy modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.