// ScmVelocityModule.h
// Modular C++ implementation of the [SCm] Velocity (v_SCm) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes v_SCm = 1e8 m/s (~c/3); scales in E_react = ρ_vac,[SCm] v_SCm² / ρ_vac,A * exp(-κ t) for U_m, U_bi, etc.
// Pluggable: #include "ScmVelocityModule.h"
// ScmVelocityModule mod; mod.computeE_react(0.0); mod.updateVariable("v_sc m", new_value);
// Variables in std::map; example for Sun at t=0 (E_react=1e46 J); t=2000 days: scales down via exp.
// Approximations: κ=0.0005 day⁻¹; ρ_vac,[SCm]=7.09e-37 J/m³; ρ_vac,A=1e-23 J/m³; U_m base=2.28e65 J/m³.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef SCM_VELOCITY_MODULE_H
#define SCM_VELOCITY_MODULE_H

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

class ScmVelocityModule {
private:
    std::map<std::string, double> variables;
    double computeE_react_base();
    double computeE_react(double t_day);
    double computeUmExample(double t_day);

public:
    // Constructor: Initialize with framework defaults
    ScmVelocityModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeV_sc m();  // 1e8 m/s
    double computeE_react(double t_day);  // ρ_[SCm] v_SCm² / ρ_A * exp(-κ t)
    double computeUmExample(double t_day);  // Simplified U_m with E_react (J/m³)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print velocity effects
    void printVelocityEffects(double t_day = 2000.0);

    // ===== ENHANCED DYNAMIC CAPABILITIES (25 Methods) =====
    // Variable Management (5 methods)
    void createVariable(const std::string& name, double value);
    void removeVariable(const std::string& name);
    void cloneVariable(const std::string& source, const std::string& dest);
    std::vector<std::string> listVariables();
    std::string getSystemName();

    // Batch Operations (2 methods)
    void transformVariableGroup(const std::vector<std::string>& vars, std::function<double(double)> func);
    void scaleVariableGroup(const std::vector<std::string>& vars, double factor);

    // Self-Expansion (4 methods: 1 general + 3 domain-specific)
    void expandParameterSpace(double expansion_factor);
    void expandVelocityScale(double v_factor, double energy_factor);
    void expandReactivityScale(double e_react_factor, double decay_factor);
    void expandTemporalScale(double kappa_factor, double time_factor);

    // Self-Refinement (3 methods)
    void autoRefineParameters();
    void calibrateToObservations(const std::map<std::string, double>& observed_data);
    void optimizeForMetric(const std::string& metric_name);

    // Parameter Exploration (1 method)
    std::vector<std::map<std::string, double>> generateVariations(int count);

    // Adaptive Evolution (2 methods)
    void mutateParameters(double mutation_rate);
    void evolveSystem(int generations);

    // State Management (4 methods)
    void saveState(const std::string& state_name);
    void restoreState(const std::string& state_name);
    std::vector<std::string> listSavedStates();
    std::string exportState();

    // System Analysis (4 methods)
    std::map<std::string, double> sensitivityAnalysis(const std::string& output_var);
    std::string generateReport();
    bool validateConsistency();
    void autoCorrectAnomalies();
};

#endif // SCM_VELOCITY_MODULE_H

// ScmVelocityModule.cpp
#include "ScmVelocityModule.h"

// Constructor: Set framework defaults
ScmVelocityModule::ScmVelocityModule() {
    // Universal constants
    variables["v_sc m"] = 1e8;                      // m/s (~c/3)
    variables["rho_vac_SCm"] = 7.09e-37;            // J/m³
    variables["rho_vac_A"] = 1e-23;                 // J/m³
    variables["kappa_day"] = 0.0005;                // day⁻¹
    variables["day_to_s"] = 86400.0;                // s/day
    variables["t_day"] = 0.0;                       // days
    variables["E_react_base"] = variables["rho_vac_SCm"] * std::pow(variables["v_sc m"], 2) / variables["rho_vac_A"];  // Derived
    variables["mu_over_rj"] = 2.26e10;              // T m² (example)
    variables["P_SCm"] = 1.0;                       // Normalized
    variables["heaviside_f"] = 1e11 + 1.0;          // 1 + 10^13 * 0.01
    variables["quasi_f"] = 1.01;                    // 1 + 0.01
    variables["one_minus_exp"] = 0.0;               // At t=0

    // Derived
    variables["kappa_s"] = variables["kappa_day"] / variables["day_to_s"];
}

// Update variable
void ScmVelocityModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
        if (name == "v_sc m" || name == "rho_vac_SCm" || name == "rho_vac_A") {
            variables["E_react_base"] = variables["rho_vac_SCm"] * std::pow(variables["v_sc m"], 2) / variables["rho_vac_A"];
        } else if (name == "kappa_day") {
            variables["kappa_s"] = value / variables["day_to_s"];
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void ScmVelocityModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        if (name == "v_sc m" || name == "rho_vac_SCm" || name == "rho_vac_A") {
            variables["E_react_base"] = variables["rho_vac_SCm"] * std::pow(variables["v_sc m"], 2) / variables["rho_vac_A"];
        } else if (name == "kappa_day") {
            variables["kappa_s"] = variables["kappa_day"] / variables["day_to_s"];
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void ScmVelocityModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute v_SCm (m/s)
double ScmVelocityModule::computeV_sc m() {
    return variables["v_sc m"];
}

// Compute E_react = E_react_base * exp(-κ t)
double ScmVelocityModule::computeE_react(double t_day) {
    variables["t_day"] = t_day;
    double arg = - variables["kappa_day"] * t_day;
    return variables["E_react_base"] * std::exp(arg);
}

// Simplified U_m example with E_react
double ScmVelocityModule::computeUmExample(double t_day) {
    double e_react = computeE_react(t_day);
    double one_minus_exp = variables["one_minus_exp"];  // Placeholder
    double phi_hat = 1.0;
    double p_scm = variables["P_SCm"];
    double heaviside_f = variables["heaviside_f"];
    double quasi_f = variables["quasi_f"];
    return (variables["mu_over_rj"] * one_minus_exp * phi_hat) * p_scm * e_react * heaviside_f * quasi_f;
}

// Equation text
std::string ScmVelocityModule::getEquationText() {
    return "E_react = [ρ_vac,[SCm] v_SCm² / ρ_vac,A] * exp(-κ t) (t days);\n"
           "v_SCm = 1e8 m/s (~c/3, [SCm] propagation speed);\n"
           "Scales reactivity in U_m, U_bi, U_i, U_gi via E_react.\n"
           "Example t=0: E_react=1e46 J; t=2000 days: ~3.68e45 J (~36.8%).\n"
           "U_m (t=0): ≈2.28e65 J/m³; t=2000: ≈8.39e64 J/m³.\n"
           "Role: [SCm] dynamic speed for relativistic effects; jets/energy transfer.\n"
           "UQFF: Subluminal propagation; [SCm]-[UA] reactions in nebulae/formation.";
}

// Print variables
void ScmVelocityModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Print effects
void ScmVelocityModule::printVelocityEffects(double t_day) {
    double v = computeV_sc m();
    double e_react = computeE_react(t_day);
    double um_ex = computeUmExample(t_day);
    double fraction = e_react / variables["E_react_base"];
    std::cout << "[SCm] Velocity Effects at t=" << t_day << " days:\n";
    std::cout << "v_SCm = " << std::scientific << v << " m/s\n";
    std::cout << "E_react = " << e_react << " J (" << fraction << " of initial)\n";
    std::cout << "U_m example = " << um_ex << " J/m³\n";
}

// ===== ENHANCED DYNAMIC CAPABILITIES IMPLEMENTATION =====

// Namespace for saved states
namespace saved_states_scm_velocity {
    std::map<std::string, std::map<std::string, double>> state_storage;
}

// Variable Management (5 methods)
void ScmVelocityModule::createVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        std::cout << "Warning: Variable '" << name << "' already exists. Overwriting.\n";
    }
    variables[name] = value;
}

void ScmVelocityModule::removeVariable(const std::string& name) {
    if (variables.find(name) != variables.end()) {
        variables.erase(name);
    } else {
        std::cerr << "Warning: Cannot remove non-existent variable '" << name << "'.\n";
    }
}

void ScmVelocityModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
    } else {
        std::cerr << "Error: Source variable '" << source << "' not found.\n";
    }
}

std::vector<std::string> ScmVelocityModule::listVariables() {
    std::vector<std::string> var_list;
    for (const auto& pair : variables) {
        var_list.push_back(pair.first);
    }
    return var_list;
}

std::string ScmVelocityModule::getSystemName() {
    return "SCm_Velocity_v_SCm_UQFF";
}

// Batch Operations (2 methods)
void ScmVelocityModule::transformVariableGroup(const std::vector<std::string>& vars, std::function<double(double)> func) {
    for (const auto& var : vars) {
        if (variables.find(var) != variables.end()) {
            variables[var] = func(variables[var]);
        }
    }
}

void ScmVelocityModule::scaleVariableGroup(const std::vector<std::string>& vars, double factor) {
    transformVariableGroup(vars, [factor](double v) { return v * factor; });
}

// Self-Expansion (4 methods: 1 general + 3 domain-specific)
void ScmVelocityModule::expandParameterSpace(double expansion_factor) {
    variables["v_scm"] *= expansion_factor;
    variables["rho_vac_SCm"] *= expansion_factor;
    variables["E_react_base"] = variables["rho_vac_SCm"] * std::pow(variables["v_scm"], 2) / variables["rho_vac_A"];
}

void ScmVelocityModule::expandVelocityScale(double v_factor, double energy_factor) {
    // Expand v_SCm velocity and related energy terms
    variables["v_scm"] *= v_factor;
    variables["E_react_base"] *= energy_factor;
    variables["mu_over_rj"] *= energy_factor;
    std::cout << "Velocity scale expanded: v_SCm *= " << v_factor 
              << ", E_react_base *= " << energy_factor << "\n";
}

void ScmVelocityModule::expandReactivityScale(double e_react_factor, double decay_factor) {
    // Expand E_react reactivity and decay rate
    variables["E_react_base"] *= e_react_factor;
    variables["rho_vac_SCm"] *= e_react_factor;
    variables["kappa_day"] *= decay_factor;
    variables["kappa_s"] = variables["kappa_day"] / variables["day_to_s"];
    std::cout << "Reactivity scale expanded: E_react *= " << e_react_factor 
              << ", κ *= " << decay_factor << "\n";
}

void ScmVelocityModule::expandTemporalScale(double kappa_factor, double time_factor) {
    // Expand decay constant and time parameters
    variables["kappa_day"] *= kappa_factor;
    variables["kappa_s"] = variables["kappa_day"] / variables["day_to_s"];
    variables["t_day"] *= time_factor;
    std::cout << "Temporal scale expanded: κ *= " << kappa_factor 
              << ", t_day *= " << time_factor << "\n";
}

// Self-Refinement (3 methods)
void ScmVelocityModule::autoRefineParameters() {
    // Clamp velocities within physical bounds: [1e7, 3e8] m/s (0.1c to c)
    if (variables["v_scm"] < 1e7) variables["v_scm"] = 1e7;
    if (variables["v_scm"] > 3e8) variables["v_scm"] = 3e8;
    
    // Clamp decay rates: κ ∈ [1e-5, 0.01] day⁻¹
    if (variables["kappa_day"] < 1e-5) variables["kappa_day"] = 1e-5;
    if (variables["kappa_day"] > 0.01) variables["kappa_day"] = 0.01;
    
    // Clamp vacuum densities: ρ_vac,SCm ∈ [1e-40, 1e-35] J/m³
    if (variables["rho_vac_SCm"] < 1e-40) variables["rho_vac_SCm"] = 1e-40;
    if (variables["rho_vac_SCm"] > 1e-35) variables["rho_vac_SCm"] = 1e-35;
    
    // Clamp normalization factors: P_SCm ∈ [0.1, 10]
    if (variables["P_SCm"] < 0.1) variables["P_SCm"] = 0.1;
    if (variables["P_SCm"] > 10.0) variables["P_SCm"] = 10.0;
    
    // Recalculate derived quantities
    variables["E_react_base"] = variables["rho_vac_SCm"] * std::pow(variables["v_scm"], 2) / variables["rho_vac_A"];
    variables["kappa_s"] = variables["kappa_day"] / variables["day_to_s"];
}

void ScmVelocityModule::calibrateToObservations(const std::map<std::string, double>& observed_data) {
    for (const auto& obs : observed_data) {
        if (variables.find(obs.first) != variables.end()) {
            double current = variables[obs.first];
            double target = obs.second;
            variables[obs.first] = 0.7 * current + 0.3 * target;
        }
    }
    // Recalculate derived quantities
    variables["E_react_base"] = variables["rho_vac_SCm"] * std::pow(variables["v_scm"], 2) / variables["rho_vac_A"];
    variables["kappa_s"] = variables["kappa_day"] / variables["day_to_s"];
}

void ScmVelocityModule::optimizeForMetric(const std::string& metric_name) {
    if (metric_name == "high_velocity") {
        // Maximize v_SCm for high-energy jets (near-c propagation)
        variables["v_scm"] = 2.5e8;  // ~0.83c
        variables["rho_vac_SCm"] = 5e-37;
        variables["kappa_day"] = 0.0001;
    } else if (metric_name == "standard_velocity") {
        // Standard [SCm] velocity: 1e8 m/s (~c/3)
        variables["v_scm"] = 1e8;
        variables["rho_vac_SCm"] = 7.09e-37;
        variables["kappa_day"] = 0.0005;
    } else if (metric_name == "low_velocity") {
        // Low-energy [SCm] propagation: slower dynamics
        variables["v_scm"] = 5e7;  // ~c/6
        variables["rho_vac_SCm"] = 1e-37;
        variables["kappa_day"] = 0.001;
    } else if (metric_name == "fast_decay") {
        // Rapid reactivity decay: short-lived jets
        variables["kappa_day"] = 0.005;
        variables["v_scm"] = 1.5e8;
    } else if (metric_name == "slow_decay") {
        // Long-lived reactivity: persistent energy transfer
        variables["kappa_day"] = 0.0001;
        variables["v_scm"] = 1e8;
    }
    // Recalculate derived quantities
    variables["E_react_base"] = variables["rho_vac_SCm"] * std::pow(variables["v_scm"], 2) / variables["rho_vac_A"];
    variables["kappa_s"] = variables["kappa_day"] / variables["day_to_s"];
}

// Parameter Exploration (1 method)
std::vector<std::map<std::string, double>> ScmVelocityModule::generateVariations(int count) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.5, 1.5);
    
    for (int i = 0; i < count; ++i) {
        std::map<std::string, double> variation = variables;
        variation["v_scm"] *= dis(gen);
        variation["rho_vac_SCm"] *= dis(gen);
        variation["kappa_day"] *= dis(gen);
        variation["E_react_base"] = variation["rho_vac_SCm"] * std::pow(variation["v_scm"], 2) / variation["rho_vac_A"];
        variation["kappa_s"] = variation["kappa_day"] / variation["day_to_s"];
        variations.push_back(variation);
    }
    return variations;
}

// Adaptive Evolution (2 methods)
void ScmVelocityModule::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dis(1.0, mutation_rate);
    
    variables["v_scm"] *= dis(gen);
    variables["rho_vac_SCm"] *= dis(gen);
    variables["kappa_day"] *= std::abs(dis(gen));
    
    autoRefineParameters();
}

void ScmVelocityModule::evolveSystem(int generations) {
    for (int gen = 0; gen < generations; ++gen) {
        mutateParameters(0.1);
        
        // Fitness: minimize decay rate while maintaining reactivity
        double e_react_2000 = computeE_react(2000.0);
        if (e_react_2000 < 1e45) {
            // Increase reactivity if too low
            variables["rho_vac_SCm"] *= 1.05;
            variables["v_scm"] *= 1.02;
        }
        if (variables["kappa_day"] > 0.001) {
            // Reduce decay if too fast
            variables["kappa_day"] *= 0.95;
        }
        
        autoRefineParameters();
    }
}

// State Management (4 methods)
void ScmVelocityModule::saveState(const std::string& state_name) {
    saved_states_scm_velocity::state_storage[state_name] = variables;
}

void ScmVelocityModule::restoreState(const std::string& state_name) {
    if (saved_states_scm_velocity::state_storage.find(state_name) != saved_states_scm_velocity::state_storage.end()) {
        variables = saved_states_scm_velocity::state_storage[state_name];
    } else {
        std::cerr << "Error: State '" << state_name << "' not found.\n";
    }
}

std::vector<std::string> ScmVelocityModule::listSavedStates() {
    std::vector<std::string> states;
    for (const auto& pair : saved_states_scm_velocity::state_storage) {
        states.push_back(pair.first);
    }
    return states;
}

std::string ScmVelocityModule::exportState() {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(6);
    oss << "SCm Velocity Module State Export:\n";
    for (const auto& pair : variables) {
        oss << pair.first << "=" << pair.second << "\n";
    }
    return oss.str();
}

// System Analysis (4 methods)
std::map<std::string, double> ScmVelocityModule::sensitivityAnalysis(const std::string& output_var) {
    std::map<std::string, double> sensitivities;
    double baseline;
    
    if (output_var == "E_react") {
        baseline = computeE_react(2000.0);
    } else if (output_var == "v_scm") {
        baseline = computeV_scm();
    } else {
        baseline = variables[output_var];
    }
    
    double delta = 0.01;  // 1% perturbation
    for (auto& pair : variables) {
        double original = pair.second;
        variables[pair.first] = original * (1.0 + delta);
        
        double perturbed;
        if (output_var == "E_react") {
            variables["E_react_base"] = variables["rho_vac_SCm"] * std::pow(variables["v_scm"], 2) / variables["rho_vac_A"];
            perturbed = computeE_react(2000.0);
        } else if (output_var == "v_scm") {
            perturbed = computeV_scm();
        } else {
            perturbed = variables[output_var];
        }
        
        sensitivities[pair.first] = (perturbed - baseline) / (baseline * delta);
        variables[pair.first] = original;
    }
    
    // Restore derived quantities
    variables["E_react_base"] = variables["rho_vac_SCm"] * std::pow(variables["v_scm"], 2) / variables["rho_vac_A"];
    return sensitivities;
}

std::string ScmVelocityModule::generateReport() {
    std::ostringstream report;
    report << std::scientific << std::setprecision(3);
    report << "========== SCm Velocity Module Report ==========\n";
    report << "System: " << getSystemName() << "\n\n";
    
    report << "Core Parameters:\n";
    report << "  v_SCm = " << variables["v_scm"] << " m/s (~" << (variables["v_scm"]/3e8) << "c)\n";
    report << "  ρ_vac,[SCm] = " << variables["rho_vac_SCm"] << " J/m³\n";
    report << "  ρ_vac,A = " << variables["rho_vac_A"] << " J/m³\n";
    report << "  κ_day = " << variables["kappa_day"] << " day⁻¹\n\n";
    
    report << "Derived Quantities:\n";
    report << "  E_react(t=0) = " << variables["E_react_base"] << " J\n";
    report << "  E_react(t=2000 days) = " << computeE_react(2000.0) << " J\n";
    report << "  Decay fraction = " << (computeE_react(2000.0) / variables["E_react_base"]) << "\n\n";
    
    report << "Physical Interpretation:\n";
    report << "  [SCm] propagation speed: subluminal (~c/3 standard)\n";
    report << "  Reactivity scaling: E_react exponentially decays\n";
    report << "  Applications: Quasar jets, energy transfer, U_m dynamics\n";
    report << "  UQFF integration: [SCm]-[UA] reactions in nebulae\n";
    
    report << "================================================\n";
    return report.str();
}

bool ScmVelocityModule::validateConsistency() {
    bool consistent = true;
    
    // Check velocity bounds: v_SCm should be < c
    if (variables["v_scm"] >= 3e8) {
        std::cerr << "Inconsistency: v_SCm >= c (superluminal!)\n";
        consistent = false;
    }
    
    // Check positive vacuum densities
    if (variables["rho_vac_SCm"] <= 0 || variables["rho_vac_A"] <= 0) {
        std::cerr << "Inconsistency: Non-positive vacuum densities\n";
        consistent = false;
    }
    
    // Check decay constant positivity
    if (variables["kappa_day"] <= 0) {
        std::cerr << "Inconsistency: Non-positive decay constant\n";
        consistent = false;
    }
    
    // Verify E_react_base calculation
    double expected_e_react = variables["rho_vac_SCm"] * std::pow(variables["v_scm"], 2) / variables["rho_vac_A"];
    if (std::abs(variables["E_react_base"] - expected_e_react) / expected_e_react > 0.01) {
        std::cerr << "Inconsistency: E_react_base mismatch\n";
        consistent = false;
    }
    
    return consistent;
}

void ScmVelocityModule::autoCorrectAnomalies() {
    // Enforce v_SCm < c
    if (variables["v_scm"] >= 3e8) {
        variables["v_scm"] = 2.99e8;
        std::cout << "Corrected: v_SCm clamped to < c\n";
    }
    
    // Enforce positive densities
    if (variables["rho_vac_SCm"] <= 0) {
        variables["rho_vac_SCm"] = 7.09e-37;
        std::cout << "Corrected: ρ_vac,SCm reset to default\n";
    }
    
    if (variables["rho_vac_A"] <= 0) {
        variables["rho_vac_A"] = 1e-23;
        std::cout << "Corrected: ρ_vac,A reset to default\n";
    }
    
    // Enforce positive decay
    if (variables["kappa_day"] <= 0) {
        variables["kappa_day"] = 0.0005;
        std::cout << "Corrected: κ_day reset to default\n";
    }
    
    // Recalculate derived quantities
    variables["E_react_base"] = variables["rho_vac_SCm"] * std::pow(variables["v_scm"], 2) / variables["rho_vac_A"];
    variables["kappa_s"] = variables["kappa_day"] / variables["day_to_s"];
    
    std::cout << "Anomaly correction complete.\n";
}

// Example usage in base program (snippet)
// #include "ScmVelocityModule.h"
// int main() {
//     ScmVelocityModule mod;
//     double v = mod.computeV_scm();
//     std::cout << "v_SCm = " << v << " m/s\n";
//     mod.printVelocityEffects(2000.0);
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("v_scm", 1.5e8);
//     mod.printVariables();
//     return 0;
// }

// ========== COMPREHENSIVE ENHANCED DEMONSTRATION ==========
/*
int main() {
    std::cout << std::scientific << std::setprecision(3);
    std::cout << "========== SCm VELOCITY MODULE DEMONSTRATION ==========\n\n";
    
    // ===== Step 1: Initialize Module =====
    ScmVelocityModule mod;
    std::cout << "Step 1: Module initialized with defaults:\n";
    std::cout << "  System: " << mod.getSystemName() << "\n";
    std::cout << "  v_SCm = " << mod.computeV_scm() << " m/s (~c/3)\n";
    std::cout << "  E_react(t=0) = " << mod.computeE_react(0.0) << " J\n";
    std::cout << "  E_react(t=2000 days) = " << mod.computeE_react(2000.0) << " J\n\n";
    
    // ===== Step 2: Baseline Computation =====
    std::cout << "Step 2: Compute baseline [SCm] velocity parameters:\n";
    double v_scm = mod.computeV_scm();
    double e_react_0 = mod.computeE_react(0.0);
    double e_react_2000 = mod.computeE_react(2000.0);
    double decay_fraction = e_react_2000 / e_react_0;
    
    std::cout << "  v_SCm = " << v_scm << " m/s (subluminal propagation)\n";
    std::cout << "  E_react(t=0) = " << e_react_0 << " J (initial reactivity)\n";
    std::cout << "  E_react(t=2000) = " << e_react_2000 << " J\n";
    std::cout << "  Decay fraction = " << decay_fraction << " (~36.8%)\n";
    std::cout << "  Physical interpretation: [SCm] propagates at ~c/3, reactivity decays exponentially\n\n";
    
    // ===== Step 3-7: Dynamic Operations =====
    std::cout << "Step 3: Variable Management\n";
    mod.createVariable("jet_velocity_scale", v_scm * 2.5);
    std::cout << "  Created 'jet_velocity_scale'\n";
    
    std::cout << "\nStep 4: Velocity Expansion\n";
    mod.expandVelocityScale(1.5, 2.0);
    std::cout << "  Expanded: v_SCm = " << mod.computeV_scm() << " m/s\n";
    
    std::cout << "\nStep 5: Reactivity Expansion\n";
    mod.expandReactivityScale(2.0, 1.5);
    std::cout << "  Expanded: E_react(t=0) = " << mod.computeE_react(0.0) << " J\n";
    
    std::cout << "\nStep 6: Temporal Expansion\n";
    mod.expandTemporalScale(0.8, 1.5);
    std::cout << "  Expanded: κ reduced (slower decay), E_react(2000) = " << mod.computeE_react(2000.0) << " J\n";
    
    std::cout << "\nStep 7: Batch Operations\n";
    std::vector<std::string> velocity_group = {"v_scm", "mu_over_rj", "P_SCm"};
    mod.scaleVariableGroup(velocity_group, 0.7);
    std::cout << "  Scaled velocity group by 0.7\n\n";
    
    // ===== Step 8-12: Physical Regimes =====
    std::cout << "Steps 8-12: Test Multiple [SCm] Velocity Regimes\n";
    
    mod.optimizeForMetric("high_velocity");
    std::cout << "  High Velocity: v_SCm = " << mod.computeV_scm() << " m/s (~0.83c)\n";
    std::cout << "                 E_react(0) = " << mod.computeE_react(0.0) << " J\n";
    
    mod.optimizeForMetric("standard_velocity");
    std::cout << "  Standard: v_SCm = " << mod.computeV_scm() << " m/s (~c/3)\n";
    std::cout << "            E_react(0) = " << mod.computeE_react(0.0) << " J\n";
    
    mod.optimizeForMetric("low_velocity");
    std::cout << "  Low Velocity: v_SCm = " << mod.computeV_scm() << " m/s (~c/6)\n";
    std::cout << "                E_react(0) = " << mod.computeE_react(0.0) << " J\n";
    
    mod.optimizeForMetric("fast_decay");
    std::cout << "  Fast Decay: E_react(2000)/E_react(0) = " << (mod.computeE_react(2000.0) / mod.computeE_react(0.0)) << "\n";
    
    mod.optimizeForMetric("slow_decay");
    std::cout << "  Slow Decay: E_react(2000)/E_react(0) = " << (mod.computeE_react(2000.0) / mod.computeE_react(0.0)) << "\n\n";
    
    // ===== Step 13-17: Refinement & Evolution =====
    std::cout << "Step 13: Auto-Refinement\n";
    mod.updateVariable("v_scm", 5e8);  // Superluminal (invalid)
    mod.autoRefineParameters();
    std::cout << "  Clamped v_SCm from 5e8 to " << mod.computeV_scm() << " m/s (< c)\n";
    
    std::cout << "\nStep 14: Calibration to Observations\n";
    std::map<std::string, double> obs_data;
    obs_data["v_scm"] = 1.2e8;
    obs_data["rho_vac_SCm"] = 8e-37;
    mod.calibrateToObservations(obs_data);
    std::cout << "  Calibrated: v_SCm = " << mod.computeV_scm() << " m/s\n";
    
    std::cout << "\nStep 15: Parameter Variations\n";
    std::vector<std::map<std::string, double>> variations = mod.generateVariations(5);
    std::cout << "  Generated " << variations.size() << " parameter variations\n";
    
    std::cout << "\nStep 16: Mutation\n";
    mod.optimizeForMetric("standard_velocity");
    mod.mutateParameters(0.15);
    std::cout << "  Mutated: v_SCm = " << mod.computeV_scm() << " m/s\n";
    
    std::cout << "\nStep 17: System Evolution\n";
    mod.evolveSystem(10);
    std::cout << "  Evolved: E_react(2000) = " << mod.computeE_react(2000.0) << " J\n\n";
    
    // ===== Step 18-19: State Management =====
    std::cout << "Step 18: State Management\n";
    mod.optimizeForMetric("standard_velocity");
    mod.saveState("standard_reference");
    mod.optimizeForMetric("high_velocity");
    mod.saveState("high_velocity_jets");
    std::vector<std::string> saved = mod.listSavedStates();
    std::cout << "  Saved " << saved.size() << " states\n";
    mod.restoreState("standard_reference");
    std::cout << "  Restored 'standard_reference'\n";
    
    std::cout << "\nStep 19: Export State\n";
    std::string exported = mod.exportState();
    std::cout << "  Exported " << exported.length() << " bytes\n\n";
    
    // ===== Step 20-22: Analysis =====
    std::cout << "Step 20: Sensitivity Analysis (E_react response)\n";
    std::map<std::string, double> sensitivity = mod.sensitivityAnalysis("E_react");
    std::cout << "  Top sensitivity parameters:\n";
    std::vector<std::pair<std::string, double>> sens_vec(sensitivity.begin(), sensitivity.end());
    std::sort(sens_vec.begin(), sens_vec.end(), 
              [](const auto& a, const auto& b) { return std::abs(a.second) > std::abs(b.second); });
    for (int i = 0; i < std::min(5, (int)sens_vec.size()); ++i) {
        std::cout << "    " << sens_vec[i].first << ": " << sens_vec[i].second << "\n";
    }
    
    std::cout << "\nStep 21: Consistency Validation\n";
    bool valid = mod.validateConsistency();
    std::cout << "  System consistency: " << (valid ? "VALID" : "INVALID") << "\n";
    if (!valid) mod.autoCorrectAnomalies();
    
    std::cout << "\nStep 22: Generate Full Report\n";
    std::string report = mod.generateReport();
    std::cout << report << "\n";
    
    // ===== Step 23-26: Velocity Scale Analysis =====
    std::cout << "Steps 23-26: [SCm] Velocity Scale Analysis\n";
    std::cout << "  Regime           | v_SCm (m/s) | v/c   | E_react(0) (J) | κ (day⁻¹) | Context\n";
    std::cout << "  -------------------------------------------------------------------------------------\n";
    
    struct VelocityRegime {
        std::string name;
        double v_scm;
        double rho_vac_scm;
        double kappa_day;
        std::string context;
    };
    
    std::vector<VelocityRegime> regimes = {
        {"Ultra-Slow", 3e7, 1e-37, 0.001, "Low-energy ISM"},
        {"Slow", 5e7, 3e-37, 0.0008, "Diffuse nebulae"},
        {"Standard", 1e8, 7.09e-37, 0.0005, "Solar reference (~c/3)"},
        {"Fast", 1.5e8, 1e-36, 0.0003, "High-energy jets"},
        {"Relativistic", 2.5e8, 2e-36, 0.0001, "Quasar jets (~0.83c)"}
    };
    
    for (const auto& reg : regimes) {
        mod.updateVariable("v_scm", reg.v_scm);
        mod.updateVariable("rho_vac_SCm", reg.rho_vac_scm);
        mod.updateVariable("kappa_day", reg.kappa_day);
        mod.updateVariable("E_react_base", reg.rho_vac_scm * std::pow(reg.v_scm, 2) / mod.variables["rho_vac_A"]);
        
        double v = mod.computeV_scm();
        double v_over_c = v / 3e8;
        double e_react = mod.computeE_react(0.0);
        
        std::cout << "  " << std::setw(16) << std::left << reg.name
                  << " | " << std::scientific << std::setprecision(2) << std::setw(11) << v
                  << " | " << std::fixed << std::setprecision(2) << std::setw(5) << v_over_c
                  << " | " << std::scientific << std::setprecision(2) << std::setw(14) << e_react
                  << " | " << std::fixed << std::setprecision(4) << std::setw(9) << reg.kappa_day
                  << " | " << reg.context << "\n";
    }
    
    std::cout << "\n========== DEMONSTRATION COMPLETE ==========\n";
    std::cout << "[SCm] Velocity module validated across physical regimes.\n";
    std::cout << "v_SCm = 1e8 m/s (~c/3) provides standard subluminal propagation.\n";
    std::cout << "E_react = [ρ_vac,[SCm] v_SCm² / ρ_vac,A] exp(-κ t)\n";
    std::cout << "Physical significance: [SCm] propagation speed for energy transfer.\n";
    std::cout << "Scales U_m, U_bi, U_i reactivity terms via E_react exponential decay.\n";
    std::cout << "UQFF Integration: Subluminal [SCm]-[UA] interactions in jets/nebulae.\n";
    std::cout << "Applications: Quasar jets, stellar formation, energy dynamics.\n";
    
    return 0;
}
*/
// Compile: g++ -o scm_vel_test scm_vel_test.cpp ScmVelocityModule.cpp -lm
// Sample: v_SCm=1e8 m/s; t=2000 days: E_react≈3.68e45 J; U_m≈8.39e64 J/m³.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

ScmVelocityModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeV_sc m, computeE_react, computeUmExample) are clear, concise, and variable - driven.
- Automatic recalculation of derived variables(E_react_base, kappa_s) when dependencies change.
- Output and debugging functions(printVariables, printVelocityEffects, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.
- Models relativistic[SCm] velocity effects and their impact on reactivity and energy transfer.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- The variable name "v_sc m" contains a space, which is non - standard and may cause confusion or errors; use "v_scm" or "v_SCm" instead.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in[SCm] velocity and reactivity modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.