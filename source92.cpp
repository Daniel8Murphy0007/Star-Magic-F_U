// BuoyancyCouplingModule.h
// Modular C++ implementation of the Buoyancy Coupling Constants (?_i) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes the Universal Buoyancy terms U_bi = -?_i * U_gi * ?_g * (M_bh / d_g) * E_react for i=1 to 4 (Ug1-Ug4).
// Pluggable: #include "BuoyancyCouplingModule.h"
// BuoyancyCouplingModule mod; mod.computeU_bi(1); mod.updateVariable("beta", new_value);
// Variables in std::map; ?_i=0.6 uniform (unitless); opposes gravity with 60% scaling.
// Approximations: cos(? t_n)=1 at t_n=0; ?_sw * ?_vac,sw ?0; U_UA=1; computes per i or sum.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef BUOYANCY_COUPLING_MODULE_H
#define BUOYANCY_COUPLING_MODULE_H

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

class BuoyancyCouplingModule {
private:
    std::map<std::string, double> variables;
    std::vector<double> computeAllU_bi();

public:
    // Constructor: Initialize with framework defaults
    BuoyancyCouplingModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeBeta(int i);  // ?_i = 0.6 for all i
    double computeU_bi(int i);  // U_bi for specific i (Ug1-4)
    std::vector<double> computeAllU_bi();  // All four U_bi
    double computeF_U_contribution();  // Sum ?_i terms in F_U

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print all U_bi
    void printU_bi();

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
    void expandBuoyancyScale(double factor);
    void expandGravityScale(double factor);
    void expandEnergyScale(double factor);

    // 4. Self-Refinement (3 methods)
    void autoRefineParameters(double tolerance);
    void calibrateToObservations(const std::map<std::string, double>& obs_data);
    void optimizeForMetric(std::function<double(BuoyancyCouplingModule&)> metric);

    // 5. Parameter Exploration (2 methods)
    std::vector<std::map<std::string, double>> generateVariations(int n_variations);
    std::map<std::string, double> findOptimalParameters(std::function<double(BuoyancyCouplingModule&)> objective, int iterations);

    // 6. Adaptive Evolution (2 methods)
    void mutateParameters(double mutation_rate);
    void evolveSystem(int generations, std::function<double(BuoyancyCouplingModule&)> fitness);

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

#endif // BUOYANCY_COUPLING_MODULE_H

// BuoyancyCouplingModule.cpp
#include "BuoyancyCouplingModule.h"

// Constructor: Set framework defaults
BuoyancyCouplingModule::BuoyancyCouplingModule() {
    // Universal constants
    variables["beta"] = 0.6;                        // ?_i uniform (unitless)
    variables["Omega_g"] = 7.3e-16;                 // rad/s (galactic spin)
    variables["M_bh"] = 8.15e36;                    // kg (black hole mass)
    variables["d_g"] = 2.55e20;                     // m (galactic distance)
    variables["E_react"] = 1.0;                     // Reactive energy (normalized)
    variables["epsilon_sw"] = 0.001;                // Swirl factor
    variables["rho_vac_sw"] = 8e-21;                // J/m^3
    variables["U_UA"] = 1.0;                        // Universal Aether factor
    variables["t_n"] = 0.0;                         // Time node (s)
    variables["pi"] = 3.141592653589793;

    // U_gi defaults (example from doc for Ug1; others placeholder)
    variables["U_g1"] = 1.39e26;                    // J/m^3
    variables["U_g2"] = 1e25;                       // Placeholder J/m^3
    variables["U_g3"] = 1e24;                       // Placeholder J/m^3
    variables["U_g4"] = 1e23;                       // Placeholder J/m^3
}

// Update variable
void BuoyancyCouplingModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
    // Dependency: If U_gi changes, U_bi updates in compute
}

// Add delta
void BuoyancyCouplingModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void BuoyancyCouplingModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute ?_i (uniform 0.6)
double BuoyancyCouplingModule::computeBeta(int i) {
    return variables["beta"];  // Uniform for all i
}

// Compute U_bi for specific i
double BuoyancyCouplingModule::computeU_bi(int i) {
    std::string ug_key = "U_g" + std::to_string(i);
    if (variables.find(ug_key) == variables.end()) {
        std::cerr << "U_g" << i << " not found. Using 1e26 default." << std::endl;
        return -0.6 * 1e26 * variables["Omega_g"] * (variables["M_bh"] / variables["d_g"]) * variables["E_react"];
    }
    double U_gi = variables[ug_key];
    double beta_i = computeBeta(i);
    double M_bh_over_d_g = variables["M_bh"] / variables["d_g"];
    double swirl_factor = 1.0 + variables["epsilon_sw"] * variables["rho_vac_sw"];
    double cos_term = std::cos(variables["pi"] * variables["t_n"]);
    return -beta_i * U_gi * variables["Omega_g"] * M_bh_over_d_g * variables["E_react"] * swirl_factor * variables["U_UA"] * cos_term;
}

// Compute all U_bi
std::vector<double> BuoyancyCouplingModule::computeAllU_bi() {
    std::vector<double> u_bi(4);
    for (int i = 1; i <= 4; ++i) {
        u_bi[i-1] = computeU_bi(i);
    }
    return u_bi;
}

// Contribution to F_U (sum of -?_i terms)
double BuoyancyCouplingModule::computeF_U_contribution() {
    auto all_u_bi = computeAllU_bi();
    double sum = 0.0;
    for (double val : all_u_bi) {
        sum += val;
    }
    return sum;
}

// Equation text
std::string BuoyancyCouplingModule::getEquationText() {
    return "U_bi = -?_i * U_gi * ?_g * (M_bh / d_g) * E_react * (1 + ?_sw * ?_vac,sw) * U_UA * cos(? t_n)\n"
           "Where ?_i = 0.6 (unitless, uniform for i=1-4: Ug1-Ug4);\n"
           "Opposes gravity: 60% scaling of gravitational term.\n"
           "In F_U: ? [k_i U_gi - ?_i U_gi ?_g (M_bh/d_g) E_react] + other terms.\n"
           "Role: Stabilizes systems (e.g., molecular clouds, nebulae); counteracts Ug collapse.\n"
           "Example Ug1: U_b1 ? -1.94e27 J/m� (at t_n=0, Sun params)."
           "UQFF: Uniform buoyancy across scales; tunable for refinements.";
}

// Print variables
void BuoyancyCouplingModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Print U_bi
void BuoyancyCouplingModule::printU_bi() {
    auto all_u_bi = computeAllU_bi();
    std::cout << "Universal Buoyancy Terms U_bi (J/m³):\n";
    for (int i = 1; i <= 4; ++i) {
        std::cout << "U_b" << i << " = " << std::scientific << all_u_bi[i-1] << std::endl;
    }
    std::cout << "F_U Buoyancy Contribution (sum): " << std::scientific << computeF_U_contribution() << std::endl;
}

// ===== Implementation of Dynamic Self-Update & Self-Expansion Capabilities =====

namespace {
    // Static storage for saved states
    std::map<std::string, std::map<std::string, double>> buoyancy_saved_states;
}

// 1. Variable Management

void BuoyancyCouplingModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void BuoyancyCouplingModule::removeVariable(const std::string& name) {
    auto it = variables.find(name);
    if (it != variables.end()) {
        variables.erase(it);
    }
}

void BuoyancyCouplingModule::cloneVariable(const std::string& source, const std::string& dest) {
    auto it = variables.find(source);
    if (it != variables.end()) {
        variables[dest] = it->second;
    }
}

std::vector<std::string> BuoyancyCouplingModule::listVariables() {
    std::vector<std::string> var_names;
    for (const auto& pair : variables) {
        var_names.push_back(pair.first);
    }
    return var_names;
}

// 2. Batch Operations

void BuoyancyCouplingModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        auto it = variables.find(name);
        if (it != variables.end()) {
            it->second = func(it->second);
        }
    }
}

void BuoyancyCouplingModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// 3. Self-Expansion

void BuoyancyCouplingModule::expandParameterSpace(const std::vector<std::string>& new_params) {
    for (const auto& param : new_params) {
        if (variables.find(param) == variables.end()) {
            variables[param] = 0.0;
        }
    }
}

void BuoyancyCouplingModule::expandBuoyancyScale(double factor) {
    // Scale buoyancy coefficient and related terms
    variables["beta"] *= factor;
    
    // Scale swirl factor (affects buoyancy)
    if (variables.find("epsilon_sw") != variables.end()) {
        variables["epsilon_sw"] *= factor;
    }
    
    // Scale vacuum swirl density
    if (variables.find("rho_vac_sw") != variables.end()) {
        variables["rho_vac_sw"] *= factor;
    }
    
    // Scale U_UA factor
    if (variables.find("U_UA") != variables.end()) {
        variables["U_UA"] *= factor;
    }
}

void BuoyancyCouplingModule::expandGravityScale(double factor) {
    // Scale gravity-related terms: U_g1 through U_g4
    std::vector<std::string> gravity_vars = {"U_g1", "U_g2", "U_g3", "U_g4"};
    scaleVariableGroup(gravity_vars, factor);
    
    // Scale galactic parameters
    if (variables.find("Omega_g") != variables.end()) {
        variables["Omega_g"] *= factor;
    }
    if (variables.find("M_bh") != variables.end()) {
        variables["M_bh"] *= factor;
    }
}

void BuoyancyCouplingModule::expandEnergyScale(double factor) {
    // Scale energy-related terms
    std::vector<std::string> energy_vars = {"E_react", "U_g1", "U_g2", "U_g3", "U_g4", "rho_vac_sw"};
    scaleVariableGroup(energy_vars, factor);
}

// 4. Self-Refinement

void BuoyancyCouplingModule::autoRefineParameters(double tolerance) {
    // Ensure beta is within physical bounds [0, 1]
    if (variables["beta"] < 0.0) {
        variables["beta"] = 0.0;
    } else if (variables["beta"] > 1.0) {
        variables["beta"] = 1.0;
    }
    
    // Ensure E_react is normalized
    if (std::abs(variables["E_react"] - 1.0) > tolerance && variables["E_react"] > 0) {
        // Keep as-is if non-standard scenario, but validate positivity
        if (variables["E_react"] < 0) {
            variables["E_react"] = 1.0;
        }
    }
    
    // Ensure U_UA is positive
    if (variables["U_UA"] <= 0) {
        variables["U_UA"] = 1.0;
    }
    
    // Ensure physical positivity of masses and densities
    if (variables["M_bh"] <= 0) {
        variables["M_bh"] = 8.15e36;
    }
    if (variables["rho_vac_sw"] < 0) {
        variables["rho_vac_sw"] = 8e-21;
    }
}

void BuoyancyCouplingModule::calibrateToObservations(const std::map<std::string, double>& obs_data) {
    for (const auto& obs : obs_data) {
        if (variables.find(obs.first) != variables.end()) {
            variables[obs.first] = obs.second;
        }
    }
    // Auto-refine to maintain consistency
    autoRefineParameters(1e-10);
}

void BuoyancyCouplingModule::optimizeForMetric(std::function<double(BuoyancyCouplingModule&)> metric) {
    double best_score = metric(*this);
    std::map<std::string, double> best_state = variables;
    
    std::default_random_engine gen;
    std::uniform_real_distribution<double> dist(0.9, 1.1);
    
    for (int iter = 0; iter < 100; iter++) {
        // Mutate key buoyancy parameters
        std::vector<std::string> key_params = {"beta", "U_g1", "U_g2", "U_g3", "U_g4", "Omega_g", "E_react"};
        for (const auto& param : key_params) {
            if (variables.find(param) != variables.end()) {
                variables[param] *= dist(gen);
            }
        }
        autoRefineParameters(1e-10);
        
        double score = metric(*this);
        if (score > best_score) {
            best_score = score;
            best_state = variables;
        } else {
            variables = best_state;
        }
    }
}

// 5. Parameter Exploration

std::vector<std::map<std::string, double>> BuoyancyCouplingModule::generateVariations(int n_variations) {
    std::vector<std::map<std::string, double>> variations;
    std::default_random_engine gen;
    std::uniform_real_distribution<double> dist(0.8, 1.2);
    
    std::map<std::string, double> original = variables;
    std::vector<std::string> vary_params = {"beta", "U_g1", "U_g2", "U_g3", "U_g4", "Omega_g", "M_bh", "E_react", "epsilon_sw"};
    
    for (int i = 0; i < n_variations; i++) {
        for (const auto& param : vary_params) {
            if (variables.find(param) != variables.end()) {
                variables[param] = original[param] * dist(gen);
            }
        }
        autoRefineParameters(1e-10);
        variations.push_back(variables);
    }
    
    variables = original;
    return variations;
}

std::map<std::string, double> BuoyancyCouplingModule::findOptimalParameters(std::function<double(BuoyancyCouplingModule&)> objective, int iterations) {
    double best_score = objective(*this);
    std::map<std::string, double> best_params = variables;
    
    std::default_random_engine gen;
    std::uniform_real_distribution<double> dist(0.5, 1.5);
    
    for (int iter = 0; iter < iterations; iter++) {
        std::vector<std::string> opt_params = {"beta", "U_g1", "U_g2", "U_g3", "U_g4", "Omega_g", "E_react"};
        for (const auto& param : opt_params) {
            if (variables.find(param) != variables.end()) {
                variables[param] *= dist(gen);
            }
        }
        autoRefineParameters(1e-10);
        
        double score = objective(*this);
        if (score > best_score) {
            best_score = score;
            best_params = variables;
        }
    }
    
    variables = best_params;
    return best_params;
}

// 6. Adaptive Evolution

void BuoyancyCouplingModule::mutateParameters(double mutation_rate) {
    std::default_random_engine gen;
    std::uniform_real_distribution<double> dist(-mutation_rate, mutation_rate);
    
    std::vector<std::string> mutable_params = {"beta", "U_g1", "U_g2", "U_g3", "U_g4", 
                                                 "Omega_g", "M_bh", "d_g", "E_react", 
                                                 "epsilon_sw", "rho_vac_sw", "U_UA"};
    for (const auto& param : mutable_params) {
        if (variables.find(param) != variables.end()) {
            double mutation = 1.0 + dist(gen);
            variables[param] *= mutation;
        }
    }
    autoRefineParameters(1e-10);
}

void BuoyancyCouplingModule::evolveSystem(int generations, std::function<double(BuoyancyCouplingModule&)> fitness) {
    for (int gen = 0; gen < generations; gen++) {
        double current_fitness = fitness(*this);
        std::map<std::string, double> current_state = variables;
        
        mutateParameters(0.1);
        
        double new_fitness = fitness(*this);
        if (new_fitness < current_fitness) {
            variables = current_state;  // Revert if fitness decreased
        }
    }
}

// 7. State Management

void BuoyancyCouplingModule::saveState(const std::string& label) {
    buoyancy_saved_states[label] = variables;
}

void BuoyancyCouplingModule::restoreState(const std::string& label) {
    auto it = buoyancy_saved_states.find(label);
    if (it != buoyancy_saved_states.end()) {
        variables = it->second;
    }
}

std::vector<std::string> BuoyancyCouplingModule::listSavedStates() {
    std::vector<std::string> state_labels;
    for (const auto& pair : buoyancy_saved_states) {
        state_labels.push_back(pair.first);
    }
    return state_labels;
}

std::map<std::string, double> BuoyancyCouplingModule::exportState() {
    return variables;
}

// 8. System Analysis

std::map<std::string, double> BuoyancyCouplingModule::sensitivityAnalysis(const std::string& var_name, double delta) {
    std::map<std::string, double> sensitivity;
    
    auto it = variables.find(var_name);
    if (it == variables.end()) {
        return sensitivity;
    }
    
    double original_val = it->second;
    
    // Test sensitivity for all four buoyancy terms
    for (int i = 1; i <= 4; i++) {
        variables[var_name] = original_val * (1.0 + delta);
        autoRefineParameters(1e-10);
        double u_bi_plus = computeU_bi(i);
        
        variables[var_name] = original_val * (1.0 - delta);
        autoRefineParameters(1e-10);
        double u_bi_minus = computeU_bi(i);
        
        double sens = (u_bi_plus - u_bi_minus) / (2.0 * delta * original_val);
        sensitivity["U_b" + std::to_string(i)] = sens;
    }
    
    // Test sensitivity for total F_U contribution
    variables[var_name] = original_val * (1.0 + delta);
    autoRefineParameters(1e-10);
    double f_u_plus = computeF_U_contribution();
    
    variables[var_name] = original_val * (1.0 - delta);
    autoRefineParameters(1e-10);
    double f_u_minus = computeF_U_contribution();
    
    double f_u_sens = (f_u_plus - f_u_minus) / (2.0 * delta * original_val);
    sensitivity["F_U_total"] = f_u_sens;
    
    variables[var_name] = original_val;
    autoRefineParameters(1e-10);
    return sensitivity;
}

std::string BuoyancyCouplingModule::generateReport() {
    std::ostringstream report;
    report << "===== UQFF Buoyancy Coupling Module Report =====\n";
    report << std::scientific;
    report << "Buoyancy Coefficient: β = " << variables["beta"] << " (unitless)\n";
    report << "Galactic Parameters:\n";
    report << "  Ω_g = " << variables["Omega_g"] << " rad/s\n";
    report << "  M_bh = " << variables["M_bh"] << " kg (" << variables["M_bh"]/1.989e30 << " Msun)\n";
    report << "  d_g = " << variables["d_g"] << " m\n";
    report << "Energy: E_react = " << variables["E_react"] << "\n";
    report << "Swirl: ε_sw = " << variables["epsilon_sw"] << ", ρ_vac,sw = " << variables["rho_vac_sw"] << " J/m³\n";
    report << "Universal Aether: U_UA = " << variables["U_UA"] << "\n";
    report << "Time node: t_n = " << variables["t_n"] << " s\n\n";
    
    report << "Gravity Terms (U_gi):\n";
    for (int i = 1; i <= 4; i++) {
        std::string key = "U_g" + std::to_string(i);
        report << "  U_g" << i << " = " << variables[key] << " J/m³\n";
    }
    report << "\n";
    
    auto all_u_bi = computeAllU_bi();
    report << "Buoyancy Terms (U_bi):\n";
    for (int i = 1; i <= 4; i++) {
        report << "  U_b" << i << " = " << all_u_bi[i-1] << " J/m³\n";
    }
    report << "\n";
    
    double f_u_contrib = computeF_U_contribution();
    report << "F_U Buoyancy Contribution (sum): " << f_u_contrib << " J/m³\n";
    report << "Opposition to gravity: " << (variables["beta"] * 100.0) << "%\n";
    report << "Saved states: " << buoyancy_saved_states.size() << "\n";
    report << "===============================================\n";
    return report.str();
}

bool BuoyancyCouplingModule::validateConsistency() {
    bool valid = true;
    
    // Check beta in [0, 1]
    if (variables["beta"] < 0.0 || variables["beta"] > 1.0) {
        valid = false;
    }
    
    // Check physical positivity
    if (variables["M_bh"] <= 0 || variables["d_g"] <= 0) {
        valid = false;
    }
    
    // Check U_gi positivity (gravity terms should be positive)
    for (int i = 1; i <= 4; i++) {
        std::string key = "U_g" + std::to_string(i);
        if (variables.find(key) != variables.end() && variables[key] <= 0) {
            valid = false;
        }
    }
    
    // Check U_UA positivity
    if (variables["U_UA"] <= 0) {
        valid = false;
    }
    
    // Check E_react non-negativity
    if (variables["E_react"] < 0) {
        valid = false;
    }
    
    // Check Omega_g physical range
    if (variables["Omega_g"] < 0 || variables["Omega_g"] > 1e-10) {  // Galactic rotation ~10^-15 to 10^-16
        valid = false;
    }
    
    return valid;
}

void BuoyancyCouplingModule::autoCorrectAnomalies() {
    // Correct beta to [0, 1]
    if (variables["beta"] < 0.0) {
        variables["beta"] = 0.0;
    } else if (variables["beta"] > 1.0) {
        variables["beta"] = 1.0;
    }
    
    // Enforce physical defaults for critical parameters
    if (variables["M_bh"] <= 0) {
        variables["M_bh"] = 8.15e36;
    }
    if (variables["d_g"] <= 0) {
        variables["d_g"] = 2.55e20;
    }
    if (variables["U_UA"] <= 0) {
        variables["U_UA"] = 1.0;
    }
    if (variables["E_react"] < 0) {
        variables["E_react"] = 1.0;
    }
    if (variables["Omega_g"] < 0) {
        variables["Omega_g"] = 7.3e-16;
    }
    
    // Ensure U_gi are positive
    for (int i = 1; i <= 4; i++) {
        std::string key = "U_g" + std::to_string(i);
        if (variables.find(key) != variables.end() && variables[key] <= 0) {
            variables[key] = 1e26 / std::pow(10, i-1);  // Default scaling
        }
    }
    
    // Ensure vacuum density is physical
    if (variables["rho_vac_sw"] < 0) {
        variables["rho_vac_sw"] = 8e-21;
    }
}

// Example usage in base program (snippet)
// #include "BuoyancyCouplingModule.h"
// int main() {
//     BuoyancyCouplingModule mod;
//     double u_b1 = mod.computeU_bi(1);
//     std::cout << "U_b1 = " << u_b1 << " J/m³\n";
//     mod.printU_bi();
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("beta", 0.7);
//     mod.printVariables();
//
//     // ===== Demonstrate Dynamic Self-Update & Self-Expansion =====
//     
//     // 1. Variable management
//     mod.createVariable("beta_refined", 0.65);
//     mod.cloneVariable("U_g1", "U_g1_backup");
//     std::cout << "Variables: " << mod.listVariables().size() << " total\n";
//     
//     // 2. Batch operations on gravity terms
//     std::vector<std::string> gravity_group = {"U_g1", "U_g2", "U_g3", "U_g4"};
//     mod.scaleVariableGroup(gravity_group, 1.15);  // 15% gravity increase
//     
//     // 3. Self-expansion
//     mod.expandBuoyancyScale(1.1);  // 10% stronger buoyancy opposition
//     mod.expandGravityScale(1.05);  // 5% gravity enhancement
//     mod.expandEnergyScale(1.08);  // 8% energy boost
//     std::cout << "After expansion: β = " << mod.exportState()["beta"] << "\n";
//     
//     // 4. Self-refinement
//     mod.autoRefineParameters(1e-10);
//     std::map<std::string, double> obs = {{"beta", 0.58}, {"Omega_g", 7.5e-16}};
//     mod.calibrateToObservations(obs);
//     
//     // 5. Parameter exploration (optimize buoyancy balance)
//     auto balance_objective = [](BuoyancyCouplingModule& m) {
//         double f_u = m.computeF_U_contribution();
//         return -std::abs(f_u + 1e27);  // Target specific opposition strength
//     };
//     mod.optimizeForMetric(balance_objective);
//     
//     // 6. Generate buoyancy scenario variations
//     auto variations = mod.generateVariations(5);
//     std::cout << "Generated " << variations.size() << " buoyancy scenarios\n";
//     
//     // 7. State management
//     mod.saveState("optimal_buoyancy");
//     mod.expandBuoyancyScale(0.8);  // Test weaker buoyancy
//     mod.restoreState("optimal_buoyancy");  // Revert
//     std::cout << "Saved states: " << mod.listSavedStates().size() << "\n";
//     
//     // 8. Sensitivity analysis for beta
//     auto beta_sensitivity = mod.sensitivityAnalysis("beta", 0.1);
//     std::cout << "Beta sensitivity across U_bi terms: " << beta_sensitivity.size() << " values\n";
//     
//     // 9. System validation
//     bool valid = mod.validateConsistency();
//     std::cout << "System consistency: " << (valid ? "VALID" : "INVALID") << "\n";
//     if (!valid) mod.autoCorrectAnomalies();
//     
//     // 10. Comprehensive report
//     std::cout << mod.generateReport();
//     
//     // 11. Adaptive evolution (optimize buoyancy over 15 generations)
//     auto buoyancy_fitness = [](BuoyancyCouplingModule& m) {
//         auto all_u_bi = m.computeAllU_bi();
//         double variance = 0.0;
//         double mean = 0.0;
//         for (double val : all_u_bi) mean += val;
//         mean /= all_u_bi.size();
//         for (double val : all_u_bi) variance += (val - mean) * (val - mean);
//         return -variance;  // Minimize variance (uniform opposition)
//     };
//     mod.evolveSystem(15, buoyancy_fitness);
//     std::cout << "Evolved system over 15 generations\n";
//     
//     // 12. Multi-term sensitivity comparison
//     std::cout << "Sensitivity analysis:\n";
//     for (const std::string& param : {"beta", "Omega_g", "E_react"}) {
//         auto sens = mod.sensitivityAnalysis(param, 0.05);
//         std::cout << "  " << param << ": F_U sensitivity = " << sens["F_U_total"] << "\n";
//     }
//     
//     // 13. Final state export
//     auto final_state = mod.exportState();
//     std::cout << "Final β = " << final_state["beta"] << "\n";
//     std::cout << "Final U_g1 = " << final_state["U_g1"] << " J/m³\n";
//     std::cout << "Final F_U contribution = " << mod.computeF_U_contribution() << " J/m³\n";
//
//     return 0;
// }
// Compile: g++ -o buoyancy_test buoyancy_test.cpp BuoyancyCouplingModule.cpp -lm
// Sample Output: U_b1 ? -1.94e27 J/m�; sum opposes gravity by ~60% scaled.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

BuoyancyCouplingModule Evaluation

Strengths :
-Modular, extensible design for computing buoyancy coupling constants(?_i) and universal buoyancy terms(U_bi) in the UQFF framework.
- Clear encapsulation of variables and buoyancy terms using std::map and std::vector, supporting dynamic updates and easy extension.
- Implements core physical concepts : ?_i scaling, opposition to gravity, and contributions to the unified field(F_U).
- Approximations and physical meaning are well - documented in comments and equation text.
- Output functions for variable state and U_bi terms support debugging and transparency.
- Uniform ?_i(0.6) simplifies analysis and tuning; supports per - term and summed contributions.
- Handles dynamic updates to variables and recalculates dependent terms as needed.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for flexibility and scalability.
- Some calculations use magic numbers or lack explanatory comments; define named constants and clarify logic.
- Minimal error handling(e.g., division by zero, invalid variable names); add validation for robustness.
- Unit consistency should be checked and documented for all physical quantities.
- For large - scale or performance - critical simulations, optimize data structures and reduce redundant calculations.
- std::map and std::vector are flexible but may be less efficient than structured types for very large models.
- Expand documentation for function purposes and physical meaning.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in buoyancy coupling modeling.It implements the UQFF buoyancy concept faithfully and adapts to various scenarios.For production or high - performance applications, address the recommendations for improved robustness, maintainability, and scalability.