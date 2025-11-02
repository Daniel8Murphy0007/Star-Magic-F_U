// UgCouplingModule.h
// Modular C++ implementation of the Coupling Constants for Ug Ranges (k_i) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes scaled Universal Gravity terms k_i * U_gi for i=1-4 (Ug1-Ug4), with k1=1.5, k2=1.2, k3=1.8, k4=1.0 (unitless).
// Pluggable: #include "UgCouplingModule.h"
// UgCouplingModule mod; mod.computeSumK_Ugi(); mod.updateVariable("U_g1", new_value);
// Variables in std::map; sum contributes to F_U; placeholders for full U_gi equations.
// Approximations: t_n=0, cos(? t_n)=1; ?_def=0, etc.; example values from Sun at t=0.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef UG_COUPLING_MODULE_H
#define UG_COUPLING_MODULE_H

#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <functional>
#include <random>
#include <algorithm>
#include <sstream>

class UgCouplingModule {
private:
    std::map<std::string, double> variables;
    std::vector<double> k_values;  // [k1, k2, k3, k4]
    std::vector<double> computeAllK_Ugi();

public:
    // Constructor: Initialize with framework defaults
    UgCouplingModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeK_i(int i);  // k_i for specific i (1-4)
    double computeU_gi(int i);  // Placeholder U_gi (J/m^3)
    double computeK_Ugi(int i);  // k_i * U_gi
    std::vector<double> computeAllK_Ugi();  // All four k_i * U_gi
    double computeSumK_Ugi();  // Sum for F_U contribution

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print all k_i * U_gi
    void printK_Ugi();

    // ===== Dynamic Self-Update & Self-Expansion Capabilities =====
    
    // 1. Variable Management (4 methods)
    void createVariable(const std::string& name, double value);
    void removeVariable(const std::string& name);
    void cloneVariable(const std::string& source, const std::string& dest);
    std::vector<std::string> listVariables();

    // 2. Batch Operations (2 methods)
    void transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func);
    void scaleVariableGroup(const std::vector<std::string>& names, double factor);

    // 3. Self-Expansion (4 methods: parameter space + 3 domain-specific scales)
    void expandParameterSpace(const std::vector<std::string>& new_params);
    void expandCouplingScale(double factor);
    void expandGravityScale(double factor);
    void expandEnergyScale(double factor);

    // 4. Self-Refinement (3 methods)
    void autoRefineParameters(double tolerance);
    void calibrateToObservations(const std::map<std::string, double>& obs_data);
    void optimizeForMetric(std::function<double(UgCouplingModule&)> metric);

    // 5. Parameter Exploration (2 methods)
    std::vector<std::map<std::string, double>> generateVariations(int n_variations);
    std::map<std::string, double> findOptimalParameters(std::function<double(UgCouplingModule&)> objective, int iterations);

    // 6. Adaptive Evolution (2 methods)
    void mutateParameters(double mutation_rate);
    void evolveSystem(int generations, std::function<double(UgCouplingModule&)> fitness);

    // 7. State Management (4 methods)
    void saveState(const std::string& label);
    void restoreState(const std::string& label);
    std::vector<std::string> listSavedStates();
    std::map<std::string, double> exportState();

    // 8. System Analysis (4 methods)
    std::map<std::string, double> sensitivityAnalysis(const std::string& var_name, double delta);
    std::string generateReport();
    bool validateConsistency();
    void autoCorrectAnomalies();
};

#endif // UG_COUPLING_MODULE_H

// UgCouplingModule.cpp
#include "UgCouplingModule.h"

// Constructor: Set framework defaults
UgCouplingModule::UgCouplingModule() {
    // Coupling constants (unitless)
    k_values = {1.5, 1.2, 1.8, 1.0};               // k1=1.5, k2=1.2, k3=1.8, k4=1.0

    // U_gi defaults (example from Sun at t=0, J/m^3)
    variables["U_g1"] = 1.39e26;                    // Internal Dipole
    variables["U_g2"] = 1.18e53;                    // Outer Field Bubble
    variables["U_g3"] = 1.8e49;                     // Magnetic Strings Disk
    variables["U_g4"] = 2.50e-20;                   // Star-Black Hole Interactions

    // Shared params (placeholders)
    variables["mu_s"] = 1.0;                        // Magnetic moment
    variables["M_s"] = 1.989e30;                    // Stellar mass kg
    variables["r"] = 1e11;                          // m
    variables["alpha"] = 1e-10;                     // Decay rate s^-1
    variables["t_n"] = 0.0;                         // s
    variables["pi"] = 3.141592653589793;
    variables["delta_def"] = 0.0;                   // Deformation
    variables["rho_vac_UA"] = 7.09e-36;             // J/m^3
    variables["rho_vac_SCm"] = 7.09e-37;            // J/m^3
    variables["S_r_Rb"] = 1.0;                      // Step function
    variables["delta_sw"] = 0.0;                    // Swirl deformation
    variables["v_sw"] = 0.0;                        // Solar wind velocity
    variables["H_SCm"] = 1.0;                       // Heaviside SCm
    variables["E_react"] = 1.0;                     // Reactive energy
    variables["M_bh"] = 8.15e36;                    // kg
    variables["d_g"] = 2.55e20;                     // m
    variables["f_feedback"] = 0.0;                  // Feedback factor
}

// Update variable
void UgCouplingModule::updateVariable(const std::string& name, double value) {
    variables[name] = value;
}

// Add delta
void UgCouplingModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void UgCouplingModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute k_i (1-based index)
double UgCouplingModule::computeK_i(int i) {
    if (i < 1 || i > 4) {
        std::cerr << "Invalid i: " << i << ". Using k1." << std::endl;
        return k_values[0];
    }
    return k_values[i-1];
}

// Placeholder compute U_gi (simplified; full eqs require more params)
double UgCouplingModule::computeU_gi(int i) {
    std::string key = "U_g" + std::to_string(i);
    if (variables.find(key) != variables.end()) {
        return variables[key];
    }
    std::cerr << "U_g" << i << " not defined. Returning 0." << std::endl;
    return 0.0;
}

// Compute k_i * U_gi
double UgCouplingModule::computeK_Ugi(int i) {
    return computeK_i(i) * computeU_gi(i);
}

// Compute all k_i * U_gi
std::vector<double> UgCouplingModule::computeAllK_Ugi() {
    std::vector<double> k_ugi(4);
    for (int i = 1; i <= 4; ++i) {
        k_ugi[i-1] = computeK_Ugi(i);
    }
    return k_ugi;
}

// Sum k_i * U_gi for F_U
double UgCouplingModule::computeSumK_Ugi() {
    auto all = computeAllK_Ugi();
    double sum = 0.0;
    for (double val : all) {
        sum += val;
    }
    return sum;
}

// Equation text
std::string UgCouplingModule::getEquationText() {
    return "F_U = ? [k_i * U_gi(r,t,M_s,?_s,T_s,B_s,?_vac,[SCm],?_vac,[UA],t_n) - ?_i * ... ] + other terms\n"
           "k_i (unitless): k1=1.5 (Ug1 Internal Dipole), k2=1.2 (Ug2 Outer Bubble),\n"
           "k3=1.8 (Ug3 Magnetic Disk), k4=1.0 (Ug4 Star-BH).\n"
           "U_g1 = k1 * ?_s ?(M_s/r) e^{-? t} cos(? t_n) (1+?_def);\n"
           "U_g2 = k2 * (?_UA + ?_SCm) M_s / r^2 * S(r-R_b) (1+?_sw v_sw) H_SCm E_react;\n"
           "U_g3 = k3 * (?_SCm + ?_UA) ?_g M_s / d_g * cos(? t_n) (1 + f_mag_str);\n"
           "U_g4 = k4 * ?_SCm M_bh / d_g * e^{-? t} cos(? t_n) (1 + f_feedback).\n"
           "Example Sun t=0: ? k_i U_gi ?1.42e53 J/m� (Ug2 dominant).\n"
           "Role: Scales Ug strengths; normalizes contributions in F_U.";
}

// Print variables
void UgCouplingModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
    std::cout << "k_i values: ";
    for (size_t j = 0; j < k_values.size(); ++j) {
        std::cout << "k" << (j+1) << "=" << k_values[j] << " ";
    }
    std::cout << std::endl;
}

// Print k_i * U_gi
void UgCouplingModule::printK_Ugi() {
    auto all = computeAllK_Ugi();
    std::cout << "Scaled Ug Terms k_i * U_gi (J/m³):\n";
    for (int i = 1; i <= 4; ++i) {
        std::cout << "k" << i << " * U_g" << i << " = " << std::scientific << all[i-1] << std::endl;
    }
    std::cout << "Sum Σ k_i U_gi = " << std::scientific << computeSumK_Ugi() << std::endl;
}

// ===== Implementation of Dynamic Self-Update & Self-Expansion Capabilities =====

namespace {
    // Static storage for saved states
    std::map<std::string, std::map<std::string, double>> ug_coupling_saved_states;
    std::map<std::string, std::vector<double>> ug_coupling_saved_k_values;
}

// 1. Variable Management

void UgCouplingModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void UgCouplingModule::removeVariable(const std::string& name) {
    auto it = variables.find(name);
    if (it != variables.end()) {
        variables.erase(it);
    }
}

void UgCouplingModule::cloneVariable(const std::string& source, const std::string& dest) {
    auto it = variables.find(source);
    if (it != variables.end()) {
        variables[dest] = it->second;
    }
}

std::vector<std::string> UgCouplingModule::listVariables() {
    std::vector<std::string> var_names;
    for (const auto& pair : variables) {
        var_names.push_back(pair.first);
    }
    return var_names;
}

// 2. Batch Operations

void UgCouplingModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        auto it = variables.find(name);
        if (it != variables.end()) {
            it->second = func(it->second);
        }
    }
}

void UgCouplingModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// 3. Self-Expansion

void UgCouplingModule::expandParameterSpace(const std::vector<std::string>& new_params) {
    for (const auto& param : new_params) {
        if (variables.find(param) == variables.end()) {
            variables[param] = 0.0;
        }
    }
}

void UgCouplingModule::expandCouplingScale(double factor) {
    // Scale coupling constants k_i
    for (size_t i = 0; i < k_values.size(); ++i) {
        k_values[i] *= factor;
    }
}

void UgCouplingModule::expandGravityScale(double factor) {
    // Scale gravity-related terms: U_g1 through U_g4, M_s, M_bh
    std::vector<std::string> gravity_vars = {"U_g1", "U_g2", "U_g3", "U_g4", "M_s", "M_bh"};
    scaleVariableGroup(gravity_vars, factor);
}

void UgCouplingModule::expandEnergyScale(double factor) {
    // Scale energy-related terms: U_gi, E_react, vacuum densities
    std::vector<std::string> energy_vars = {"U_g1", "U_g2", "U_g3", "U_g4", "E_react", 
                                             "rho_vac_UA", "rho_vac_SCm"};
    scaleVariableGroup(energy_vars, factor);
}

// 4. Self-Refinement

void UgCouplingModule::autoRefineParameters(double tolerance) {
    // Ensure coupling constants are positive
    for (size_t i = 0; i < k_values.size(); ++i) {
        if (k_values[i] < 0) {
            k_values[i] = std::abs(k_values[i]);
        }
    }
    
    // Ensure E_react is normalized
    if (variables["E_react"] < 0) {
        variables["E_react"] = 1.0;
    }
    
    // Ensure physical positivity of masses and densities
    if (variables["M_s"] <= 0) {
        variables["M_s"] = 1.989e30;
    }
    if (variables["M_bh"] <= 0) {
        variables["M_bh"] = 8.15e36;
    }
    if (variables["rho_vac_UA"] < 0) {
        variables["rho_vac_UA"] = 7.09e-36;
    }
    if (variables["rho_vac_SCm"] < 0) {
        variables["rho_vac_SCm"] = 7.09e-37;
    }
    
    // Ensure U_gi positivity
    for (int i = 1; i <= 4; ++i) {
        std::string key = "U_g" + std::to_string(i);
        if (variables.find(key) != variables.end() && variables[key] < 0) {
            variables[key] = std::abs(variables[key]);
        }
    }
}

void UgCouplingModule::calibrateToObservations(const std::map<std::string, double>& obs_data) {
    for (const auto& obs : obs_data) {
        if (variables.find(obs.first) != variables.end()) {
            variables[obs.first] = obs.second;
        }
    }
    // Auto-sync dependencies
    autoRefineParameters(1e-10);
}

void UgCouplingModule::optimizeForMetric(std::function<double(UgCouplingModule&)> metric) {
    double best_score = metric(*this);
    std::map<std::string, double> best_state = variables;
    std::vector<double> best_k = k_values;
    
    std::default_random_engine gen;
    std::uniform_real_distribution<double> dist(0.9, 1.1);
    
    for (int iter = 0; iter < 100; iter++) {
        // Mutate key coupling parameters
        std::vector<std::string> key_params = {"U_g1", "U_g2", "U_g3", "U_g4", "M_s", "M_bh", "E_react"};
        for (const auto& param : key_params) {
            if (variables.find(param) != variables.end()) {
                variables[param] *= dist(gen);
            }
        }
        
        // Mutate coupling constants
        for (size_t i = 0; i < k_values.size(); ++i) {
            k_values[i] *= dist(gen);
        }
        
        autoRefineParameters(1e-10);
        
        double score = metric(*this);
        if (score > best_score) {
            best_score = score;
            best_state = variables;
            best_k = k_values;
        } else {
            variables = best_state;
            k_values = best_k;
        }
    }
}

// 5. Parameter Exploration

std::vector<std::map<std::string, double>> UgCouplingModule::generateVariations(int n_variations) {
    std::vector<std::map<std::string, double>> variations;
    std::default_random_engine gen;
    std::uniform_real_distribution<double> dist(0.8, 1.2);
    
    std::map<std::string, double> original = variables;
    std::vector<double> original_k = k_values;
    std::vector<std::string> vary_params = {"U_g1", "U_g2", "U_g3", "U_g4", "M_s", "M_bh", "E_react", "mu_s"};
    
    for (int i = 0; i < n_variations; i++) {
        for (const auto& param : vary_params) {
            if (variables.find(param) != variables.end()) {
                variables[param] = original[param] * dist(gen);
            }
        }
        
        // Vary coupling constants
        for (size_t j = 0; j < k_values.size(); ++j) {
            k_values[j] = original_k[j] * dist(gen);
        }
        
        autoRefineParameters(1e-10);
        
        // Store with k_values encoded in special variables
        std::map<std::string, double> variation = variables;
        for (size_t j = 0; j < k_values.size(); ++j) {
            variation["k_" + std::to_string(j+1)] = k_values[j];
        }
        variations.push_back(variation);
    }
    
    variables = original;
    k_values = original_k;
    return variations;
}

std::map<std::string, double> UgCouplingModule::findOptimalParameters(std::function<double(UgCouplingModule&)> objective, int iterations) {
    double best_score = objective(*this);
    std::map<std::string, double> best_params = variables;
    std::vector<double> best_k = k_values;
    
    std::default_random_engine gen;
    std::uniform_real_distribution<double> dist(0.5, 1.5);
    
    for (int iter = 0; iter < iterations; iter++) {
        std::vector<std::string> opt_params = {"U_g1", "U_g2", "U_g3", "U_g4", "M_s", "M_bh", "E_react"};
        for (const auto& param : opt_params) {
            if (variables.find(param) != variables.end()) {
                variables[param] *= dist(gen);
            }
        }
        
        // Vary k values
        for (size_t j = 0; j < k_values.size(); ++j) {
            k_values[j] *= dist(gen);
        }
        
        autoRefineParameters(1e-10);
        
        double score = objective(*this);
        if (score > best_score) {
            best_score = score;
            best_params = variables;
            best_k = k_values;
        }
    }
    
    variables = best_params;
    k_values = best_k;
    
    // Encode k_values into return map
    for (size_t j = 0; j < k_values.size(); ++j) {
        best_params["k_" + std::to_string(j+1)] = k_values[j];
    }
    return best_params;
}

// 6. Adaptive Evolution

void UgCouplingModule::mutateParameters(double mutation_rate) {
    std::default_random_engine gen;
    std::uniform_real_distribution<double> dist(-mutation_rate, mutation_rate);
    
    std::vector<std::string> mutable_params = {"U_g1", "U_g2", "U_g3", "U_g4", 
                                                 "M_s", "M_bh", "mu_s", "E_react", 
                                                 "alpha", "d_g", "rho_vac_UA", "rho_vac_SCm"};
    for (const auto& param : mutable_params) {
        if (variables.find(param) != variables.end()) {
            double mutation = 1.0 + dist(gen);
            variables[param] *= mutation;
        }
    }
    
    // Mutate k values
    for (size_t i = 0; i < k_values.size(); ++i) {
        double mutation = 1.0 + dist(gen);
        k_values[i] *= mutation;
    }
    
    autoRefineParameters(1e-10);
}

void UgCouplingModule::evolveSystem(int generations, std::function<double(UgCouplingModule&)> fitness) {
    for (int gen = 0; gen < generations; gen++) {
        double current_fitness = fitness(*this);
        std::map<std::string, double> current_state = variables;
        std::vector<double> current_k = k_values;
        
        mutateParameters(0.1);
        
        double new_fitness = fitness(*this);
        if (new_fitness < current_fitness) {
            variables = current_state;  // Revert if fitness decreased
            k_values = current_k;
        }
    }
}

// 7. State Management

void UgCouplingModule::saveState(const std::string& label) {
    ug_coupling_saved_states[label] = variables;
    ug_coupling_saved_k_values[label] = k_values;
}

void UgCouplingModule::restoreState(const std::string& label) {
    auto it = ug_coupling_saved_states.find(label);
    if (it != ug_coupling_saved_states.end()) {
        variables = it->second;
    }
    auto it_k = ug_coupling_saved_k_values.find(label);
    if (it_k != ug_coupling_saved_k_values.end()) {
        k_values = it_k->second;
    }
}

std::vector<std::string> UgCouplingModule::listSavedStates() {
    std::vector<std::string> state_labels;
    for (const auto& pair : ug_coupling_saved_states) {
        state_labels.push_back(pair.first);
    }
    return state_labels;
}

std::map<std::string, double> UgCouplingModule::exportState() {
    std::map<std::string, double> state = variables;
    // Include k_values in export
    for (size_t i = 0; i < k_values.size(); ++i) {
        state["k_" + std::to_string(i+1)] = k_values[i];
    }
    return state;
}

// 8. System Analysis

std::map<std::string, double> UgCouplingModule::sensitivityAnalysis(const std::string& var_name, double delta) {
    std::map<std::string, double> sensitivity;
    
    auto it = variables.find(var_name);
    if (it == variables.end()) {
        return sensitivity;
    }
    
    double original_val = it->second;
    
    // Test sensitivity for each scaled term and sum
    for (int i = 1; i <= 4; ++i) {
        variables[var_name] = original_val * (1.0 + delta);
        autoRefineParameters(1e-10);
        double k_ugi_plus = computeK_Ugi(i);
        
        variables[var_name] = original_val * (1.0 - delta);
        autoRefineParameters(1e-10);
        double k_ugi_minus = computeK_Ugi(i);
        
        double sens = (k_ugi_plus - k_ugi_minus) / (2.0 * delta * original_val);
        sensitivity["k" + std::to_string(i) + "_Ug" + std::to_string(i)] = sens;
    }
    
    // Test sensitivity for sum
    variables[var_name] = original_val * (1.0 + delta);
    autoRefineParameters(1e-10);
    double sum_plus = computeSumK_Ugi();
    
    variables[var_name] = original_val * (1.0 - delta);
    autoRefineParameters(1e-10);
    double sum_minus = computeSumK_Ugi();
    
    double sum_sens = (sum_plus - sum_minus) / (2.0 * delta * original_val);
    sensitivity["sum_k_Ugi"] = sum_sens;
    
    variables[var_name] = original_val;
    autoRefineParameters(1e-10);
    return sensitivity;
}

std::string UgCouplingModule::generateReport() {
    std::ostringstream report;
    report << "===== UQFF Ug Coupling Constants Module Report =====\n";
    report << std::scientific;
    report << "Coupling Constants k_i (unitless):\n";
    for (size_t i = 0; i < k_values.size(); ++i) {
        report << "  k" << (i+1) << " = " << k_values[i] << " (";
        if (i == 0) report << "Internal Dipole)\n";
        else if (i == 1) report << "Outer Bubble)\n";
        else if (i == 2) report << "Magnetic Disk)\n";
        else if (i == 3) report << "Star-BH)\n";
    }
    report << "\n";
    
    report << "Universal Gravity Terms U_gi (J/m³):\n";
    for (int i = 1; i <= 4; ++i) {
        std::string key = "U_g" + std::to_string(i);
        report << "  U_g" << i << " = " << variables[key] << "\n";
    }
    report << "\n";
    
    auto all_k_ugi = computeAllK_Ugi();
    report << "Scaled Terms k_i × U_gi (J/m³):\n";
    for (int i = 1; i <= 4; ++i) {
        report << "  k" << i << " × U_g" << i << " = " << all_k_ugi[i-1] << "\n";
    }
    report << "\n";
    
    double sum = computeSumK_Ugi();
    report << "Sum Σ(k_i × U_gi) = " << sum << " J/m³\n";
    report << "Dominant term: U_g" << (std::max_element(all_k_ugi.begin(), all_k_ugi.end()) - all_k_ugi.begin() + 1) << "\n\n";
    
    report << "Physical Parameters:\n";
    report << "  M_s = " << variables["M_s"] << " kg (" << variables["M_s"]/1.989e30 << " Msun)\n";
    report << "  M_bh = " << variables["M_bh"] << " kg (" << variables["M_bh"]/1.989e30 << " Msun)\n";
    report << "  E_react = " << variables["E_react"] << "\n";
    report << "  ρ_vac,UA = " << variables["rho_vac_UA"] << " J/m³\n";
    report << "  ρ_vac,SCm = " << variables["rho_vac_SCm"] << " J/m³\n";
    report << "Saved states: " << ug_coupling_saved_states.size() << "\n";
    report << "===============================================\n";
    return report.str();
}

bool UgCouplingModule::validateConsistency() {
    bool valid = true;
    
    // Check k_i positivity
    for (size_t i = 0; i < k_values.size(); ++i) {
        if (k_values[i] < 0) {
            valid = false;
        }
    }
    
    // Check physical positivity
    if (variables["M_s"] <= 0 || variables["M_bh"] <= 0) {
        valid = false;
    }
    
    // Check U_gi positivity (gravity terms should be positive)
    for (int i = 1; i <= 4; ++i) {
        std::string key = "U_g" + std::to_string(i);
        if (variables.find(key) != variables.end() && variables[key] < 0) {
            valid = false;
        }
    }
    
    // Check E_react non-negativity
    if (variables["E_react"] < 0) {
        valid = false;
    }
    
    // Check vacuum densities
    if (variables["rho_vac_UA"] < 0 || variables["rho_vac_SCm"] < 0) {
        valid = false;
    }
    
    return valid;
}

void UgCouplingModule::autoCorrectAnomalies() {
    // Correct k_i to positive
    for (size_t i = 0; i < k_values.size(); ++i) {
        if (k_values[i] < 0) {
            k_values[i] = std::abs(k_values[i]);
        }
        if (k_values[i] == 0.0) {
            // Reset to defaults
            std::vector<double> defaults = {1.5, 1.2, 1.8, 1.0};
            k_values[i] = defaults[i];
        }
    }
    
    // Enforce physical defaults for critical parameters
    if (variables["M_s"] <= 0) {
        variables["M_s"] = 1.989e30;
    }
    if (variables["M_bh"] <= 0) {
        variables["M_bh"] = 8.15e36;
    }
    if (variables["E_react"] < 0) {
        variables["E_react"] = 1.0;
    }
    
    // Ensure U_gi are positive
    std::vector<double> ug_defaults = {1.39e26, 1.18e53, 1.8e49, 2.50e-20};
    for (int i = 1; i <= 4; ++i) {
        std::string key = "U_g" + std::to_string(i);
        if (variables.find(key) != variables.end() && variables[key] < 0) {
            variables[key] = std::abs(variables[key]);
        }
        if (variables[key] == 0.0) {
            variables[key] = ug_defaults[i-1];
        }
    }
    
    // Ensure vacuum densities are physical
    if (variables["rho_vac_UA"] < 0) {
        variables["rho_vac_UA"] = 7.09e-36;
    }
    if (variables["rho_vac_SCm"] < 0) {
        variables["rho_vac_SCm"] = 7.09e-37;
    }
}

// Example usage in base program (snippet)
// #include "UgCouplingModule.h"
// int main() {
//     UgCouplingModule mod;
//     double sum = mod.computeSumK_Ugi();
//     std::cout << "Σ k_i U_gi = " << sum << " J/m³\n";
//     mod.printK_Ugi();
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("U_g3", 2e49);
//     mod.printVariables();
//
//     // ===== Demonstrate Dynamic Self-Update & Self-Expansion =====
//     
//     // 1. Variable management
//     mod.createVariable("k_custom", 1.4);
//     mod.cloneVariable("U_g1", "U_g1_backup");
//     std::cout << "Variables: " << mod.listVariables().size() << " total\n";
//     
//     // 2. Batch operations on gravity terms
//     std::vector<std::string> gravity_group = {"U_g1", "U_g2", "U_g3", "U_g4"};
//     mod.scaleVariableGroup(gravity_group, 1.1);  // 10% gravity increase
//     
//     // 3. Self-expansion
//     mod.expandCouplingScale(1.05);  // 5% stronger coupling across all k_i
//     mod.expandGravityScale(1.08);  // 8% gravity enhancement
//     mod.expandEnergyScale(1.12);  // 12% energy boost
//     std::cout << "After expansion: k1 = " << mod.computeK_i(1) << "\n";
//     
//     // 4. Self-refinement
//     mod.autoRefineParameters(1e-10);
//     std::map<std::string, double> obs = {{"U_g1", 1.5e26}, {"M_s", 2e30}};
//     mod.calibrateToObservations(obs);
//     
//     // 5. Parameter exploration (optimize coupling balance)
//     auto balance_objective = [](UgCouplingModule& m) {
//         double sum = m.computeSumK_Ugi();
//         return -std::abs(sum - 1.5e53);  // Target specific sum
//     };
//     mod.optimizeForMetric(balance_objective);
//     
//     // 6. Generate coupling scenario variations
//     auto variations = mod.generateVariations(5);
//     std::cout << "Generated " << variations.size() << " coupling scenarios\n";
//     
//     // 7. State management
//     mod.saveState("optimal_coupling");
//     mod.expandCouplingScale(0.9);  // Test weaker coupling
//     mod.restoreState("optimal_coupling");  // Revert
//     std::cout << "Saved states: " << mod.listSavedStates().size() << "\n";
//     
//     // 8. Sensitivity analysis for U_g2 (dominant term)
//     auto u_g2_sensitivity = mod.sensitivityAnalysis("U_g2", 0.1);
//     std::cout << "U_g2 sensitivity:\n";
//     for (const auto& s : u_g2_sensitivity) {
//         std::cout << "  " << s.first << ": " << s.second << "\n";
//     }
//     
//     // 9. System validation
//     bool valid = mod.validateConsistency();
//     std::cout << "System consistency: " << (valid ? "VALID" : "INVALID") << "\n";
//     if (!valid) mod.autoCorrectAnomalies();
//     
//     // 10. Comprehensive report
//     std::cout << mod.generateReport();
//     
//     // 11. Adaptive evolution (optimize coupling distribution over 12 generations)
//     auto coupling_fitness = [](UgCouplingModule& m) {
//         auto all_k_ugi = m.computeAllK_Ugi();
//         // Minimize variance for balanced contributions
//         double mean = 0.0;
//         for (double val : all_k_ugi) mean += val;
//         mean /= all_k_ugi.size();
//         double variance = 0.0;
//         for (double val : all_k_ugi) variance += (val - mean) * (val - mean);
//         return -variance;
//     };
//     mod.evolveSystem(12, coupling_fitness);
//     std::cout << "Evolved system over 12 generations\n";
//     
//     // 12. Multi-term sensitivity comparison
//     std::cout << "Sensitivity analysis for all U_gi:\n";
//     for (int i = 1; i <= 4; ++i) {
//         std::string param = "U_g" + std::to_string(i);
//         auto sens = mod.sensitivityAnalysis(param, 0.05);
//         std::cout << "  " << param << ": sum sensitivity = " << sens["sum_k_Ugi"] << "\n";
//     }
//     
//     // 13. Coupling constant exploration
//     std::cout << "Sum Σ(k_i × U_gi) at various coupling scales:\n";
//     for (double scale : {0.8, 1.0, 1.2, 1.5}) {
//         mod.expandCouplingScale(scale / mod.computeK_i(1) * 1.5);  // Normalize to k1
//         double sum_scaled = mod.computeSumK_Ugi();
//         std::cout << "  Scale " << scale << ": Σ = " << sum_scaled << " J/m³\n";
//     }
//     
//     // 14. Final state export
//     auto final_state = mod.exportState();
//     std::cout << "Final k1 = " << final_state["k_1"] << "\n";
//     std::cout << "Final U_g1 = " << final_state["U_g1"] << " J/m³\n";
//     std::cout << "Final M_s = " << final_state["M_s"]/1.989e30 << " Msun\n";
//     std::cout << "Final Sum = " << mod.computeSumK_Ugi() << " J/m³\n";
//
//     return 0;
// }
// Compile: g++ -o ug_coupling_test ug_coupling_test.cpp UgCouplingModule.cpp -lm
// Sample Output: Sum ?1.42e53 J/m�; k3 amplifies Ug3 by 1.8.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

UgCouplingModule Evaluation

Strengths :
-Modular, extensible design for computing scaled universal gravity terms(k_i * U_gi) in the UQFF framework.
- Clear encapsulation of variables and coupling constants using std::map and std::vector, supporting dynamic updates and easy extension.
- Implements core physical concepts : scaling of Ug terms via k_i, contribution to the unified field(F_U), and separation of physical roles for Ug1 - Ug4.
- Approximations and physical meaning are well - documented in comments and equation text.
- Output functions for variable state and scaled Ug terms support debugging and transparency.
- Handles dynamic updates to variables and recalculates dependent terms as needed.
- Example values and equation text provide context for scientific use and validation.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for flexibility and scalability.
- Some calculations use magic numbers or lack explanatory comments; define named constants and clarify logic.
- Minimal error handling(e.g., invalid index for k_i / U_gi, division by zero); add validation for robustness.
- Unit consistency should be checked and documented for all physical quantities.
- For large - scale or performance - critical simulations, optimize data structures and reduce redundant calculations.
- std::map and std::vector are flexible but may be less efficient than structured types for very large models.
- Expand documentation for function purposes and physical meaning.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in gravity coupling modeling.It implements the UQFF coupling concept faithfully and adapts to various scenarios.For production or high - performance applications, address the recommendations for improved robustness, maintainability, and scalability.