// SolarWindBuoyancyModule.h
// Modular C++ implementation of the Buoyancy Modulation by Solar Wind Density (?_sw) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes the modulation factor (1 + ?_sw * ?_vac,sw) in the Universal Buoyancy term U_bi, with ?_sw=0.001 (unitless).
// Pluggable: #include "SolarWindBuoyancyModule.h"
// SolarWindBuoyancyModule mod; mod.computeModulationFactor(); mod.updateVariable("epsilon_sw", new_value);
// Variables in std::map; negligible correction ~8e-24; integrates into U_b1 example computation.
// Approximations: cos(? t_n)=1; U_UA=1; ?_vac,sw as energy density (8e-21 J/m^3).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef SOLAR_WIND_BUOYANCY_MODULE_H
#define SOLAR_WIND_BUOYANCY_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <functional>
#include <vector>
#include <random>
#include <algorithm>
#include <sstream>

class SolarWindBuoyancyModule {
private:
    std::map<std::string, double> variables;
    double computeModulationFactor();
    double computeU_b1();  // Example U_b1 integration

public:
    // Constructor: Initialize with framework defaults
    SolarWindBuoyancyModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeEpsilon_sw();  // ?_sw = 0.001 (unitless)
    double computeModulationFactor();  // 1 + ?_sw * ?_vac,sw
    double computeU_b1();  // Full U_b1 with modulation

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

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
    void expandModulationScale(double factor);
    void expandEnergyScale(double factor);
    void expandGravityScale(double factor);

    // 4. Self-Refinement (3 methods)
    void autoRefineParameters(double tolerance);
    void calibrateToObservations(const std::map<std::string, double>& obs_data);
    void optimizeForMetric(std::function<double(SolarWindBuoyancyModule&)> metric);

    // 5. Parameter Exploration (2 methods)
    std::vector<std::map<std::string, double>> generateVariations(int n_variations);
    std::map<std::string, double> findOptimalParameters(std::function<double(SolarWindBuoyancyModule&)> objective, int iterations);

    // 6. Adaptive Evolution (2 methods)
    void mutateParameters(double mutation_rate);
    void evolveSystem(int generations, std::function<double(SolarWindBuoyancyModule&)> fitness);

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

#endif // SOLAR_WIND_BUOYANCY_MODULE_H

// SolarWindBuoyancyModule.cpp
#include "SolarWindBuoyancyModule.h"

// Constructor: Set framework defaults
SolarWindBuoyancyModule::SolarWindBuoyancyModule() {
    // Universal constants
    variables["epsilon_sw"] = 0.001;                // Buoyancy modulation (unitless)
    variables["rho_vac_sw"] = 8e-21;                // J/m^3 (solar wind energy density)
    variables["beta_1"] = 0.6;                      // From buoyancy coupling
    variables["U_g1"] = 1.39e26;                    // J/m^3 (Ug1 example)
    variables["Omega_g"] = 7.3e-16;                 // rad/s
    variables["M_bh"] = 8.15e36;                    // kg
    variables["d_g"] = 2.55e20;                     // m
    variables["E_react"] = 1.0;                     // Normalized
    variables["U_UA"] = 1.0;                        // Universal Aether factor
    variables["t_n"] = 0.0;                         // s
    variables["pi"] = 3.141592653589793;

    // Derived defaults
    variables["modulation_factor"] = computeModulationFactor();
}

// Update variable
void SolarWindBuoyancyModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
        if (name == "epsilon_sw" || name == "rho_vac_sw") {
            variables["modulation_factor"] = computeModulationFactor();
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void SolarWindBuoyancyModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        if (name == "epsilon_sw" || name == "rho_vac_sw") {
            variables["modulation_factor"] = computeModulationFactor();
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void SolarWindBuoyancyModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute ?_sw (fixed 0.001)
double SolarWindBuoyancyModule::computeEpsilon_sw() {
    return variables["epsilon_sw"];
}

// Compute modulation factor 1 + ?_sw * ?_vac,sw
double SolarWindBuoyancyModule::computeModulationFactor() {
    return 1.0 + variables["epsilon_sw"] * variables["rho_vac_sw"];
}

// Compute example U_b1 with modulation
double SolarWindBuoyancyModule::computeU_b1() {
    double beta_1 = variables["beta_1"];
    double U_g1 = variables["U_g1"];
    double Omega_g = variables["Omega_g"];
    double M_bh_over_d_g = variables["M_bh"] / variables["d_g"];
    double E_react = variables["E_react"];
    double mod_factor = computeModulationFactor();
    double U_UA = variables["U_UA"];
    double cos_term = std::cos(variables["pi"] * variables["t_n"]);
    return -beta_1 * U_g1 * Omega_g * M_bh_over_d_g * E_react * mod_factor * U_UA * cos_term;
}

// Equation text
std::string SolarWindBuoyancyModule::getEquationText() {
    return "Modulation Factor = 1 + ?_sw * ?_vac,sw\n"
           "Where ?_sw = 0.001 (unitless); ?_vac,sw = 8e-21 J/m�.\n"
           "In U_bi: ... * (1 + ?_sw * ?_vac,sw) * ... ?1 (negligible correction ~8e-24).\n"
           "U_b1 = -?_1 U_g1 ?_g (M_bh / d_g) * modulation * U_UA * cos(? t_n)\n"
           "? -1.94e27 J/m� (at t_n=0, Sun params; modulation ?1).\n"
           "Role: Minor solar wind density effect on buoyancy; stabilizes heliosphere/nebulae.\n"
           "UQFF: Scales counterforce to Ug; negligible but flexible for variations.";
}

// Print variables
void SolarWindBuoyancyModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// ===== Implementation of Dynamic Self-Update & Self-Expansion Capabilities =====

namespace {
    // Static storage for saved states
    std::map<std::string, std::map<std::string, double>> solar_wind_saved_states;
}

// 1. Variable Management

void SolarWindBuoyancyModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void SolarWindBuoyancyModule::removeVariable(const std::string& name) {
    auto it = variables.find(name);
    if (it != variables.end()) {
        variables.erase(it);
    }
}

void SolarWindBuoyancyModule::cloneVariable(const std::string& source, const std::string& dest) {
    auto it = variables.find(source);
    if (it != variables.end()) {
        variables[dest] = it->second;
    }
}

std::vector<std::string> SolarWindBuoyancyModule::listVariables() {
    std::vector<std::string> var_names;
    for (const auto& pair : variables) {
        var_names.push_back(pair.first);
    }
    return var_names;
}

// 2. Batch Operations

void SolarWindBuoyancyModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        auto it = variables.find(name);
        if (it != variables.end()) {
            it->second = func(it->second);
        }
    }
    // Update modulation factor if affected
    if (std::find(names.begin(), names.end(), "epsilon_sw") != names.end() ||
        std::find(names.begin(), names.end(), "rho_vac_sw") != names.end()) {
        variables["modulation_factor"] = computeModulationFactor();
    }
}

void SolarWindBuoyancyModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// 3. Self-Expansion

void SolarWindBuoyancyModule::expandParameterSpace(const std::vector<std::string>& new_params) {
    for (const auto& param : new_params) {
        if (variables.find(param) == variables.end()) {
            variables[param] = 0.0;
        }
    }
}

void SolarWindBuoyancyModule::expandModulationScale(double factor) {
    // Scale modulation parameters: epsilon_sw, rho_vac_sw
    variables["epsilon_sw"] *= factor;
    variables["rho_vac_sw"] *= factor;
    
    // Auto-sync modulation factor
    variables["modulation_factor"] = computeModulationFactor();
}

void SolarWindBuoyancyModule::expandEnergyScale(double factor) {
    // Scale energy-related terms: U_g1, rho_vac_sw, E_react
    std::vector<std::string> energy_vars = {"U_g1", "rho_vac_sw", "E_react"};
    scaleVariableGroup(energy_vars, factor);
    
    // Auto-sync modulation factor
    variables["modulation_factor"] = computeModulationFactor();
}

void SolarWindBuoyancyModule::expandGravityScale(double factor) {
    // Scale gravity-related terms: U_g1, Omega_g, M_bh
    std::vector<std::string> gravity_vars = {"U_g1", "Omega_g", "M_bh"};
    scaleVariableGroup(gravity_vars, factor);
}

// 4. Self-Refinement

void SolarWindBuoyancyModule::autoRefineParameters(double tolerance) {
    // Ensure beta_1 is within physical bounds [0, 1]
    if (variables["beta_1"] < 0.0) {
        variables["beta_1"] = 0.0;
    } else if (variables["beta_1"] > 1.0) {
        variables["beta_1"] = 1.0;
    }
    
    // Ensure epsilon_sw is small (typically << 1)
    if (variables["epsilon_sw"] < 0.0) {
        variables["epsilon_sw"] = 0.0;
    } else if (variables["epsilon_sw"] > 0.1) {
        variables["epsilon_sw"] = 0.1;  // Cap at 10%
    }
    
    // Ensure E_react is normalized
    if (variables["E_react"] < 0) {
        variables["E_react"] = 1.0;
    }
    
    // Ensure U_UA is positive
    if (variables["U_UA"] <= 0) {
        variables["U_UA"] = 1.0;
    }
    
    // Ensure physical positivity
    if (variables["M_bh"] <= 0) {
        variables["M_bh"] = 8.15e36;
    }
    if (variables["rho_vac_sw"] < 0) {
        variables["rho_vac_sw"] = 8e-21;
    }
    if (variables["U_g1"] <= 0) {
        variables["U_g1"] = 1.39e26;
    }
    
    // Recompute modulation factor
    double expected_mod = computeModulationFactor();
    if (std::abs(variables["modulation_factor"] - expected_mod) > tolerance) {
        variables["modulation_factor"] = expected_mod;
    }
}

void SolarWindBuoyancyModule::calibrateToObservations(const std::map<std::string, double>& obs_data) {
    for (const auto& obs : obs_data) {
        if (variables.find(obs.first) != variables.end()) {
            variables[obs.first] = obs.second;
        }
    }
    // Auto-sync dependencies
    variables["modulation_factor"] = computeModulationFactor();
    autoRefineParameters(1e-10);
}

void SolarWindBuoyancyModule::optimizeForMetric(std::function<double(SolarWindBuoyancyModule&)> metric) {
    double best_score = metric(*this);
    std::map<std::string, double> best_state = variables;
    
    std::default_random_engine gen;
    std::uniform_real_distribution<double> dist(0.9, 1.1);
    
    for (int iter = 0; iter < 100; iter++) {
        // Mutate key solar wind parameters
        std::vector<std::string> key_params = {"epsilon_sw", "rho_vac_sw", "beta_1", "U_g1", "Omega_g", "E_react"};
        for (const auto& param : key_params) {
            if (variables.find(param) != variables.end()) {
                variables[param] *= dist(gen);
            }
        }
        variables["modulation_factor"] = computeModulationFactor();
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

std::vector<std::map<std::string, double>> SolarWindBuoyancyModule::generateVariations(int n_variations) {
    std::vector<std::map<std::string, double>> variations;
    std::default_random_engine gen;
    std::uniform_real_distribution<double> dist(0.8, 1.2);
    
    std::map<std::string, double> original = variables;
    std::vector<std::string> vary_params = {"epsilon_sw", "rho_vac_sw", "beta_1", "U_g1", "Omega_g", "M_bh", "E_react"};
    
    for (int i = 0; i < n_variations; i++) {
        for (const auto& param : vary_params) {
            if (variables.find(param) != variables.end()) {
                variables[param] = original[param] * dist(gen);
            }
        }
        variables["modulation_factor"] = computeModulationFactor();
        autoRefineParameters(1e-10);
        variations.push_back(variables);
    }
    
    variables = original;
    return variations;
}

std::map<std::string, double> SolarWindBuoyancyModule::findOptimalParameters(std::function<double(SolarWindBuoyancyModule&)> objective, int iterations) {
    double best_score = objective(*this);
    std::map<std::string, double> best_params = variables;
    
    std::default_random_engine gen;
    std::uniform_real_distribution<double> dist(0.5, 1.5);
    
    for (int iter = 0; iter < iterations; iter++) {
        std::vector<std::string> opt_params = {"epsilon_sw", "rho_vac_sw", "beta_1", "U_g1", "Omega_g", "E_react"};
        for (const auto& param : opt_params) {
            if (variables.find(param) != variables.end()) {
                variables[param] *= dist(gen);
            }
        }
        variables["modulation_factor"] = computeModulationFactor();
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

void SolarWindBuoyancyModule::mutateParameters(double mutation_rate) {
    std::default_random_engine gen;
    std::uniform_real_distribution<double> dist(-mutation_rate, mutation_rate);
    
    std::vector<std::string> mutable_params = {"epsilon_sw", "rho_vac_sw", "beta_1", "U_g1", 
                                                 "Omega_g", "M_bh", "d_g", "E_react", "U_UA"};
    for (const auto& param : mutable_params) {
        if (variables.find(param) != variables.end()) {
            double mutation = 1.0 + dist(gen);
            variables[param] *= mutation;
        }
    }
    variables["modulation_factor"] = computeModulationFactor();
    autoRefineParameters(1e-10);
}

void SolarWindBuoyancyModule::evolveSystem(int generations, std::function<double(SolarWindBuoyancyModule&)> fitness) {
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

void SolarWindBuoyancyModule::saveState(const std::string& label) {
    solar_wind_saved_states[label] = variables;
}

void SolarWindBuoyancyModule::restoreState(const std::string& label) {
    auto it = solar_wind_saved_states.find(label);
    if (it != solar_wind_saved_states.end()) {
        variables = it->second;
    }
}

std::vector<std::string> SolarWindBuoyancyModule::listSavedStates() {
    std::vector<std::string> state_labels;
    for (const auto& pair : solar_wind_saved_states) {
        state_labels.push_back(pair.first);
    }
    return state_labels;
}

std::map<std::string, double> SolarWindBuoyancyModule::exportState() {
    return variables;
}

// 8. System Analysis

std::map<std::string, double> SolarWindBuoyancyModule::sensitivityAnalysis(const std::string& var_name, double delta) {
    std::map<std::string, double> sensitivity;
    
    auto it = variables.find(var_name);
    if (it == variables.end()) {
        return sensitivity;
    }
    
    double original_val = it->second;
    
    // Test sensitivity for modulation factor
    variables[var_name] = original_val * (1.0 + delta);
    variables["modulation_factor"] = computeModulationFactor();
    autoRefineParameters(1e-10);
    double mod_plus = computeModulationFactor();
    
    variables[var_name] = original_val * (1.0 - delta);
    variables["modulation_factor"] = computeModulationFactor();
    autoRefineParameters(1e-10);
    double mod_minus = computeModulationFactor();
    
    double mod_sens = (mod_plus - mod_minus) / (2.0 * delta * original_val);
    sensitivity["modulation_factor"] = mod_sens;
    
    // Test sensitivity for U_b1
    variables[var_name] = original_val * (1.0 + delta);
    variables["modulation_factor"] = computeModulationFactor();
    autoRefineParameters(1e-10);
    double u_b1_plus = computeU_b1();
    
    variables[var_name] = original_val * (1.0 - delta);
    variables["modulation_factor"] = computeModulationFactor();
    autoRefineParameters(1e-10);
    double u_b1_minus = computeU_b1();
    
    double u_b1_sens = (u_b1_plus - u_b1_minus) / (2.0 * delta * original_val);
    sensitivity["U_b1"] = u_b1_sens;
    
    variables[var_name] = original_val;
    variables["modulation_factor"] = computeModulationFactor();
    autoRefineParameters(1e-10);
    return sensitivity;
}

std::string SolarWindBuoyancyModule::generateReport() {
    std::ostringstream report;
    report << "===== UQFF Solar Wind Buoyancy Modulation Report =====\n";
    report << std::scientific;
    report << "Solar Wind Modulation: ε_sw = " << variables["epsilon_sw"] << " (unitless)\n";
    report << "Vacuum Energy Density: ρ_vac,sw = " << variables["rho_vac_sw"] << " J/m³\n";
    report << "Modulation Factor: (1 + ε_sw × ρ_vac,sw) = " << variables["modulation_factor"] << "\n";
    report << "Effect: " << ((variables["modulation_factor"] - 1.0) * 100.0) << "% correction\n\n";
    
    report << "Buoyancy Coefficient: β₁ = " << variables["beta_1"] << " (unitless)\n";
    report << "Galactic Parameters:\n";
    report << "  Ω_g = " << variables["Omega_g"] << " rad/s\n";
    report << "  M_bh = " << variables["M_bh"] << " kg (" << variables["M_bh"]/1.989e30 << " Msun)\n";
    report << "  d_g = " << variables["d_g"] << " m\n";
    report << "Gravity Term: U_g1 = " << variables["U_g1"] << " J/m³\n";
    report << "Energy: E_react = " << variables["E_react"] << "\n";
    report << "Universal Aether: U_UA = " << variables["U_UA"] << "\n";
    report << "Time node: t_n = " << variables["t_n"] << " s\n\n";
    
    double u_b1 = computeU_b1();
    report << "Computed U_b1 (with modulation): " << u_b1 << " J/m³\n";
    report << "Saved states: " << solar_wind_saved_states.size() << "\n";
    report << "===============================================\n";
    return report.str();
}

bool SolarWindBuoyancyModule::validateConsistency() {
    bool valid = true;
    
    // Check beta_1 in [0, 1]
    if (variables["beta_1"] < 0.0 || variables["beta_1"] > 1.0) {
        valid = false;
    }
    
    // Check epsilon_sw is reasonable (small correction)
    if (variables["epsilon_sw"] < 0.0 || variables["epsilon_sw"] > 0.1) {
        valid = false;
    }
    
    // Check physical positivity
    if (variables["M_bh"] <= 0 || variables["d_g"] <= 0 || variables["U_g1"] <= 0) {
        valid = false;
    }
    
    // Check U_UA positivity
    if (variables["U_UA"] <= 0) {
        valid = false;
    }
    
    // Check E_react non-negativity
    if (variables["E_react"] < 0) {
        valid = false;
    }
    
    // Check modulation factor consistency
    double expected_mod = computeModulationFactor();
    if (std::abs(variables["modulation_factor"] - expected_mod) > 1e-12) {
        valid = false;
    }
    
    // Check rho_vac_sw positivity
    if (variables["rho_vac_sw"] < 0) {
        valid = false;
    }
    
    return valid;
}

void SolarWindBuoyancyModule::autoCorrectAnomalies() {
    // Correct beta_1 to [0, 1]
    if (variables["beta_1"] < 0.0) {
        variables["beta_1"] = 0.0;
    } else if (variables["beta_1"] > 1.0) {
        variables["beta_1"] = 1.0;
    }
    
    // Correct epsilon_sw to [0, 0.1]
    if (variables["epsilon_sw"] < 0.0) {
        variables["epsilon_sw"] = 0.0;
    } else if (variables["epsilon_sw"] > 0.1) {
        variables["epsilon_sw"] = 0.001;  // Reset to default
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
    if (variables["U_g1"] <= 0) {
        variables["U_g1"] = 1.39e26;
    }
    if (variables["rho_vac_sw"] < 0) {
        variables["rho_vac_sw"] = 8e-21;
    }
    if (variables["Omega_g"] < 0) {
        variables["Omega_g"] = 7.3e-16;
    }
    
    // Recompute modulation factor
    variables["modulation_factor"] = computeModulationFactor();
}

// Example usage in base program (snippet)
// #include "SolarWindBuoyancyModule.h"
// int main() {
//     SolarWindBuoyancyModule mod;
//     double mod_factor = mod.computeModulationFactor();
//     std::cout << "Modulation Factor = " << mod_factor << std::endl;
//     double u_b1 = mod.computeU_b1();
//     std::cout << "U_b1 = " << u_b1 << " J/m³\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("epsilon_sw", 0.002);
//     mod.printVariables();
//
//     // ===== Demonstrate Dynamic Self-Update & Self-Expansion =====
//     
//     // 1. Variable management
//     mod.createVariable("epsilon_sw_alt", 0.0015);
//     mod.cloneVariable("rho_vac_sw", "rho_vac_sw_backup");
//     std::cout << "Variables: " << mod.listVariables().size() << " total\n";
//     
//     // 2. Batch operations on modulation parameters
//     std::vector<std::string> mod_group = {"epsilon_sw", "rho_vac_sw"};
//     mod.scaleVariableGroup(mod_group, 1.2);  // 20% stronger modulation
//     
//     // 3. Self-expansion
//     mod.expandModulationScale(1.15);  // 15% modulation enhancement
//     mod.expandEnergyScale(1.05);  // 5% energy boost
//     mod.expandGravityScale(1.08);  // 8% gravity increase
//     std::cout << "After expansion: ε_sw = " << mod.exportState()["epsilon_sw"] << "\n";
//     
//     // 4. Self-refinement
//     mod.autoRefineParameters(1e-10);
//     std::map<std::string, double> obs = {{"epsilon_sw", 0.0012}, {"rho_vac_sw", 9e-21}};
//     mod.calibrateToObservations(obs);
//     
//     // 5. Parameter exploration (optimize modulation effect)
//     auto mod_objective = [](SolarWindBuoyancyModule& m) {
//         double mod_factor = m.computeModulationFactor();
//         return -std::abs(mod_factor - 1.00001);  // Target 0.001% effect
//     };
//     mod.optimizeForMetric(mod_objective);
//     
//     // 6. Generate solar wind scenario variations
//     auto variations = mod.generateVariations(5);
//     std::cout << "Generated " << variations.size() << " solar wind scenarios\n";
//     
//     // 7. State management
//     mod.saveState("optimal_modulation");
//     mod.expandModulationScale(0.8);  // Test weaker modulation
//     mod.restoreState("optimal_modulation");  // Revert
//     std::cout << "Saved states: " << mod.listSavedStates().size() << "\n";
//     
//     // 8. Sensitivity analysis for epsilon_sw
//     auto epsilon_sensitivity = mod.sensitivityAnalysis("epsilon_sw", 0.1);
//     std::cout << "ε_sw sensitivity:\n";
//     std::cout << "  Modulation factor: " << epsilon_sensitivity["modulation_factor"] << "\n";
//     std::cout << "  U_b1: " << epsilon_sensitivity["U_b1"] << "\n";
//     
//     // 9. System validation
//     bool valid = mod.validateConsistency();
//     std::cout << "System consistency: " << (valid ? "VALID" : "INVALID") << "\n";
//     if (!valid) mod.autoCorrectAnomalies();
//     
//     // 10. Comprehensive report
//     std::cout << mod.generateReport();
//     
//     // 11. Adaptive evolution (optimize modulation over 10 generations)
//     auto modulation_fitness = [](SolarWindBuoyancyModule& m) {
//         double mod_factor = m.computeModulationFactor();
//         double u_b1 = m.computeU_b1();
//         return -(std::abs(mod_factor - 1.00001) + std::abs(u_b1 + 1.94e27) / 1e25);
//     };
//     mod.evolveSystem(10, modulation_fitness);
//     std::cout << "Evolved system over 10 generations\n";
//     
//     // 12. Multi-parameter sensitivity comparison
//     std::cout << "Sensitivity analysis:\n";
//     for (const std::string& param : {"epsilon_sw", "rho_vac_sw", "beta_1"}) {
//         auto sens = mod.sensitivityAnalysis(param, 0.05);
//         std::cout << "  " << param << ": mod_factor sens = " << sens["modulation_factor"] 
//                   << ", U_b1 sens = " << sens["U_b1"] << "\n";
//     }
//     
//     // 13. Modulation effect exploration
//     std::cout << "Modulation effect at various ε_sw:\n";
//     for (double eps : {0.0001, 0.0005, 0.001, 0.002, 0.005}) {
//         mod.updateVariable("epsilon_sw", eps);
//         double mod_f = mod.computeModulationFactor();
//         double effect = (mod_f - 1.0) * 100.0;
//         std::cout << "  ε_sw=" << eps << ": " << effect << "% correction\n";
//     }
//     
//     // 14. Final state export
//     auto final_state = mod.exportState();
//     std::cout << "Final ε_sw = " << final_state["epsilon_sw"] << "\n";
//     std::cout << "Final ρ_vac,sw = " << final_state["rho_vac_sw"] << " J/m³\n";
//     std::cout << "Final modulation = " << final_state["modulation_factor"] << "\n";
//     std::cout << "Final U_b1 = " << mod.computeU_b1() << " J/m³\n";
//
//     return 0;
// }
// Compile: g++ -o sw_mod_test sw_mod_test.cpp SolarWindBuoyancyModule.cpp -lm
// Sample Output: Modulation ?1.000000000000008; U_b1 ? -1.94e27 J/m� (unchanged).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

SolarWindBuoyancyModule Evaluation

Strengths :
-Modular, extensible design for modeling solar wind density modulation in the UQFF buoyancy framework.
- Clear encapsulation of variables using std::map, supporting dynamic updates and easy extension.
- Implements core physical concepts : modulation factor(1 + ?_sw * ?_vac, sw), integration into U_b1 computation, and solar wind effects.
- Approximations and physical meaning are well - documented in comments and equation text.
- Output functions for variable state and equation text support debugging and transparency.
- Handles dynamic updates to ?_sw and ?_vac, sw, recalculating the modulation factor as needed.
- Negligible correction is correctly handled, but the structure allows for future flexibility.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for flexibility and scalability.
- Some calculations use magic numbers or lack explanatory comments; define named constants and clarify logic.
- Minimal error handling(e.g., division by zero, invalid variable names); add validation for robustness.
- Unit consistency should be checked and documented for all physical quantities.
- For large - scale or performance - critical simulations, optimize data structures and reduce redundant calculations.
- std::map is flexible but may be less efficient than structured types for very large models.
- Expand documentation for function purposes and physical meaning.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in solar wind buoyancy modeling.It implements the UQFF modulation concept faithfully and adapts to various scenarios.For production or high - performance applications, address the recommendations for improved robustness, maintainability, and scalability.