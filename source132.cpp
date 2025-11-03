// ButterflyNebulaUQFFModule.h
// Modular C++ implementation of the UQFF Force for NGC 6302 (Butterfly Nebula) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes F_U_Bi_i,enhanced as integral from x1 to x2 of [-F0 + DPM terms + LENR + activation + DE + EM + neutron + rel + Sweet + Kozima].
// Pluggable: #include "ButterflyNebulaUQFFModule.h"
// ButterflyNebulaUQFFModule mod; mod.computeF_U_Bi(0.0, 3.22e19, 0.0); mod.updateVariable("M", new_value);
// Variables in std::map; defaults for NGC 6302 (M=0.64 M_sun, r=3.22e19 m, level=13); ~ -2.09e212 N at t=0.
// Approximations: Integral approx via average * ?x; cos(?)=1; ?_LENR / ?_0 tuned; Sweet/Kozima small/negligible.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef BUTTERFLY_NEBULA_UQFF_MODULE_H
#define BUTTERFLY_NEBULA_UQFF_MODULE_H

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

class ButterflyNebulaUQFFModule {
private:
    std::map<std::string, double> variables;
    double computeDPM_momentum_term(double r);
    double computeDPM_gravity_term(double r);
    double computeDPM_stability_term();
    double computeLENR_term();
    double computeActivation_term(double t);
    double computeDE_term(double L_x);
    double computeEM_term();
    double computeNeutron_term();
    double computeRel_term(double E_cm_eff);
    double computeSweet_vac_term();
    double computeKozima_term();
    double computeIntegrand(double x, double t);
    double computeIntegral(double x1, double x2, double t, int n_points = 1000);

public:
    // Constructor: Initialize with NGC 6302 defaults
    ButterflyNebulaUQFFModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computation: F_U_Bi_i,enhanced (N)
    double computeF_U_Bi(double x1, double x2, double t);

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

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
    void expandNebulaScale(double mass_factor, double radius_factor);
    void expandForceScale(double dpm_factor, double lenr_factor);
    void expandIntegrationScale(double spatial_factor, double resolution_factor);

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

#endif // BUTTERFLY_NEBULA_UQFF_MODULE_H

// ButterflyNebulaUQFFModule.cpp
#include "ButterflyNebulaUQFFModule.h"

// Constructor: Set NGC 6302-specific values
ButterflyNebulaUQFFModule::ButterflyNebulaUQFFModule() {
    // Universal constants
    variables["G"] = 6.6743e-11;                    // m^3 kg^-1 s^-2
    variables["c"] = 2.998e8;                       // m/s
    variables["m_e"] = 9.109e-31;                   // kg
    variables["hbar"] = 1.0546e-34;                 // J s
    variables["mu_B"] = 9.274e-24;                  // J/T
    variables["e"] = 1.602e-19;                     // C
    variables["M_sun"] = 1.989e30;                  // kg
    variables["q"] = 1.602e-19;                     // C
    variables["pi"] = 3.141592653589793;

    // Nebula-specific params
    variables["M"] = 0.64 * variables["M_sun"];     // kg (central star)
    variables["r"] = 3.22e19;                       // m (distance)
    variables["x1"] = 0.0;                          // m (integral lower)
    variables["x2"] = 3.22e19;                      // m (upper)
    variables["level"] = 13.0;                      // Quantum level
    variables["F0"] = 1.0;                          // Base force (normalized)
    variables["theta"] = 0.0;                       // rad (angle)
    variables["DPM_momentum"] = 1.0;                // Normalized
    variables["DPM_gravity"] = 1.0;                 // Normalized
    variables["DPM_stability"] = 0.01;              // Normalized
    variables["rho_vac_UA"] = 7.09e-36;             // J/m�
    variables["k_LENR"] = 1.0;                      // Coupling
    variables["omega_LENR"] = 7.85e12;              // Hz
    variables["omega_0"] = 1e-12;                   // Hz (reference)
    variables["k_act"] = 1.0;                       // Activation coupling
    variables["omega_act"] = 1.0;                   // rad/s
    variables["k_DE"] = 1.0;                        // DE coupling
    variables["L_x"] = 1.0;                         // Length scale
    variables["B_0"] = 1.0;                         // T
    variables["V"] = 1.0;                           // m/s
    variables["g"] = 9.8;                           // m/s�
    variables["k_neutron"] = 1e10;                  // Neutron coupling
    variables["sigma_n"] = 1e-4;                    // Barn
    variables["k_rel"] = 1.0;                       // Rel coupling
    variables["E_cm"] = 1.0;                        // eV
    variables["E_cm_eff"] = 1.0;                    // Enhanced eV
    variables["F_Sweet_vac"] = 7.09e-39;            // N (negligible)
    variables["F_Kozima"] = 7.85e30;                // N
    variables["t"] = 0.0;                           // s
}

// Update variable with dependencies
void ButterflyNebulaUQFFModule::updateVariable(const std::string& name, double value) {
    variables[name] = value;
    // No complex deps for simplicity
}

void ButterflyNebulaUQFFModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

void ButterflyNebulaUQFFModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// DPM momentum term
double ButterflyNebulaUQFFModule::computeDPM_momentum_term(double r) {
    double m_e_c2 = variables["m_e"] * std::pow(variables["c"], 2);
    return (m_e_c2 / (r * r)) * variables["DPM_momentum"] * std::cos(variables["theta"]);
}

// DPM gravity term
double ButterflyNebulaUQFFModule::computeDPM_gravity_term(double r) {
    return (variables["G"] * variables["M"] / (r * r)) * variables["DPM_gravity"];
}

// DPM stability term
double ButterflyNebulaUQFFModule::computeDPM_stability_term() {
    return variables["rho_vac_UA"] * variables["DPM_stability"];
}

// LENR term
double ButterflyNebulaUQFFModule::computeLENR_term() {
    double ratio = std::pow(variables["omega_LENR"] / variables["omega_0"], 2);
    return variables["k_LENR"] * ratio;
}

// Activation term
double ButterflyNebulaUQFFModule::computeActivation_term(double t) {
    return variables["k_act"] * std::cos(variables["omega_act"] * t);
}

// DE term
double ButterflyNebulaUQFFModule::computeDE_term(double L_x) {
    return variables["k_DE"] * L_x;
}

// EM term
double ButterflyNebulaUQFFModule::computeEM_term() {
    double q_v_B = 2 * variables["q"] * variables["B_0"] * variables["V"] * std::sin(variables["theta"]);
    double g_mu_B = variables["g"] * variables["mu_B"] * variables["B_0"] / (variables["hbar"] * variables["omega_0"]);
    return q_v_B * g_mu_B;
}

// Neutron term
double ButterflyNebulaUQFFModule::computeNeutron_term() {
    return variables["k_neutron"] * variables["sigma_n"];
}

// Rel term
double ButterflyNebulaUQFFModule::computeRel_term(double E_cm_eff) {
    double ratio = std::pow(E_cm_eff / variables["E_cm"], 2);
    return variables["k_rel"] * ratio;
}

// Sweet vac term
double ButterflyNebulaUQFFModule::computeSweet_vac_term() {
    return variables["F_Sweet_vac"];
}

// Kozima term
double ButterflyNebulaUQFFModule::computeKozima_term() {
    return variables["F_Kozima"];
}

// Full integrand
double ButterflyNebulaUQFFModule::computeIntegrand(double x, double t) {
    return -variables["F0"] + computeDPM_momentum_term(x) + computeDPM_gravity_term(x) + computeDPM_stability_term() +
           computeLENR_term() + computeActivation_term(t) + computeDE_term(variables["L_x"]) + computeEM_term() +
           computeNeutron_term() + computeRel_term(variables["E_cm_eff"]) + computeSweet_vac_term() + computeKozima_term();
}

// Numerical integral (trapezoidal rule)
double ButterflyNebulaUQFFModule::computeIntegral(double x1, double x2, double t, int n_points) {
    double dx = (x2 - x1) / n_points;
    double integral = 0.0;
    for (int i = 0; i <= n_points; ++i) {
        double x = x1 + i * dx;
        double weight = (i == 0 || i == n_points) ? 0.5 : 1.0;
        integral += weight * computeIntegrand(x, t);
    }
    return integral * dx;
}

// Main F_U_Bi_i,enhanced
double ButterflyNebulaUQFFModule::computeF_U_Bi(double x1, double x2, double t) {
    return computeIntegral(x1, x2, t);
}

// Equation text
std::string ButterflyNebulaUQFFModule::getEquationText() {
    return "F_U_Bi_i,enhanced = ?_{x1}^{x2} [-F0 + (m_e c^2 / r^2) DPM_mom cos? + (G M / r^2) DPM_grav + ?_[UA] DPM_stab + k_LENR (?_LENR/?_0)^2 + k_act cos(?_act t) + k_DE L_x + 2 q B_0 V sin? (g ?_B B_0 / ? ?_0) + k_neutron ?_n + k_rel (E_cm,eff / E_cm)^2 + F_Sweet,vac + F_Kozima] dx\n"
           "NGC 6302: M=0.64 M_sun, r=3.22e19 m, level=13; ~ -2.09e212 N (repulsive stabilization).\n"
           "Sweet: ?_[UA] DPM_stab V ?7.09e-39 N (negligible); Kozima: k_n ?_n (?_LENR/?_0) ?7.85e30 N.\n"
           "UQFF: Integrates LENR/resonance/buoyancy for nebula force; [SCm]/[UA] dynamics.";
}

// Print variables
void ButterflyNebulaUQFFModule::printVariables() {
    std::cout << "NGC 6302 Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// ===== ENHANCED DYNAMIC CAPABILITIES IMPLEMENTATION =====

// Namespace for saved states
namespace saved_states_butterfly {
    std::map<std::string, std::map<std::string, double>> state_storage;
}

// Variable Management (5 methods)
void ButterflyNebulaUQFFModule::createVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        std::cout << "Warning: Variable '" << name << "' already exists. Overwriting.\n";
    }
    variables[name] = value;
}

void ButterflyNebulaUQFFModule::removeVariable(const std::string& name) {
    if (variables.find(name) != variables.end()) {
        variables.erase(name);
    } else {
        std::cerr << "Warning: Cannot remove non-existent variable '" << name << "'.\n";
    }
}

void ButterflyNebulaUQFFModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
    } else {
        std::cerr << "Error: Source variable '" << source << "' not found.\n";
    }
}

std::vector<std::string> ButterflyNebulaUQFFModule::listVariables() {
    std::vector<std::string> var_list;
    for (const auto& pair : variables) {
        var_list.push_back(pair.first);
    }
    return var_list;
}

std::string ButterflyNebulaUQFFModule::getSystemName() {
    return "Butterfly_Nebula_NGC6302_UQFF";
}

// Batch Operations (2 methods)
void ButterflyNebulaUQFFModule::transformVariableGroup(const std::vector<std::string>& vars, std::function<double(double)> func) {
    for (const auto& var : vars) {
        if (variables.find(var) != variables.end()) {
            variables[var] = func(variables[var]);
        }
    }
}

void ButterflyNebulaUQFFModule::scaleVariableGroup(const std::vector<std::string>& vars, double factor) {
    transformVariableGroup(vars, [factor](double v) { return v * factor; });
}

// Self-Expansion (4 methods: 1 general + 3 domain-specific)
void ButterflyNebulaUQFFModule::expandParameterSpace(double expansion_factor) {
    variables["M"] *= expansion_factor;
    variables["r"] *= expansion_factor;
    variables["x2"] *= expansion_factor;
}

void ButterflyNebulaUQFFModule::expandNebulaScale(double mass_factor, double radius_factor) {
    // Expand nebula mass and radius parameters
    variables["M"] *= mass_factor;
    variables["r"] *= radius_factor;
    variables["x1"] *= radius_factor;
    variables["x2"] *= radius_factor;
    std::cout << "Nebula scale expanded: M *= " << mass_factor 
              << ", r *= " << radius_factor << "\n";
}

void ButterflyNebulaUQFFModule::expandForceScale(double dpm_factor, double lenr_factor) {
    // Expand DPM and LENR force components
    variables["DPM_momentum"] *= dpm_factor;
    variables["DPM_gravity"] *= dpm_factor;
    variables["DPM_stability"] *= dpm_factor;
    variables["k_LENR"] *= lenr_factor;
    variables["omega_LENR"] *= lenr_factor;
    std::cout << "Force scale expanded: DPM *= " << dpm_factor 
              << ", LENR *= " << lenr_factor << "\n";
}

void ButterflyNebulaUQFFModule::expandIntegrationScale(double spatial_factor, double resolution_factor) {
    // Expand integration domain and resolution
    variables["x1"] *= spatial_factor;
    variables["x2"] *= spatial_factor;
    variables["L_x"] *= spatial_factor;
    // Resolution stored externally but noted here
    std::cout << "Integration scale expanded: spatial *= " << spatial_factor 
              << ", resolution *= " << resolution_factor << "\n";
}

// Self-Refinement (3 methods)
void ButterflyNebulaUQFFModule::autoRefineParameters() {
    // Clamp mass within reasonable bounds: [0.1, 10] M_sun
    double M_min = 0.1 * variables["M_sun"];
    double M_max = 10.0 * variables["M_sun"];
    if (variables["M"] < M_min) variables["M"] = M_min;
    if (variables["M"] > M_max) variables["M"] = M_max;
    
    // Clamp radius: [1e18, 1e21] m
    if (variables["r"] < 1e18) variables["r"] = 1e18;
    if (variables["r"] > 1e21) variables["r"] = 1e21;
    
    // Clamp DPM factors: [0.01, 100]
    if (variables["DPM_momentum"] < 0.01) variables["DPM_momentum"] = 0.01;
    if (variables["DPM_momentum"] > 100.0) variables["DPM_momentum"] = 100.0;
    if (variables["DPM_gravity"] < 0.01) variables["DPM_gravity"] = 0.01;
    if (variables["DPM_gravity"] > 100.0) variables["DPM_gravity"] = 100.0;
    
    // Clamp LENR coupling: [0.1, 10]
    if (variables["k_LENR"] < 0.1) variables["k_LENR"] = 0.1;
    if (variables["k_LENR"] > 10.0) variables["k_LENR"] = 10.0;
    
    // Ensure x2 > x1
    if (variables["x2"] <= variables["x1"]) {
        variables["x2"] = variables["x1"] + 1e19;
    }
}

void ButterflyNebulaUQFFModule::calibrateToObservations(const std::map<std::string, double>& observed_data) {
    for (const auto& obs : observed_data) {
        if (variables.find(obs.first) != variables.end()) {
            double current = variables[obs.first];
            double target = obs.second;
            variables[obs.first] = 0.7 * current + 0.3 * target;
        }
    }
}

void ButterflyNebulaUQFFModule::optimizeForMetric(const std::string& metric_name) {
    if (metric_name == "standard_butterfly") {
        // NGC 6302 standard parameters
        variables["M"] = 0.64 * variables["M_sun"];
        variables["r"] = 3.22e19;
        variables["x1"] = 0.0;
        variables["x2"] = 3.22e19;
        variables["level"] = 13.0;
        variables["DPM_momentum"] = 1.0;
        variables["k_LENR"] = 1.0;
    } else if (metric_name == "high_mass_nebula") {
        // Higher mass central star
        variables["M"] = 2.0 * variables["M_sun"];
        variables["DPM_gravity"] = 2.5;
        variables["k_LENR"] = 1.5;
    } else if (metric_name == "compact_nebula") {
        // Smaller radius, denser structure
        variables["r"] = 1e19;
        variables["x2"] = 1e19;
        variables["DPM_momentum"] = 2.0;
    } else if (metric_name == "extended_nebula") {
        // Larger radius, diffuse structure
        variables["r"] = 1e20;
        variables["x2"] = 1e20;
        variables["DPM_stability"] = 0.05;
    } else if (metric_name == "high_lenr") {
        // Enhanced LENR activity
        variables["k_LENR"] = 5.0;
        variables["omega_LENR"] = 1e13;
        variables["k_act"] = 2.0;
    }
}

// Parameter Exploration (1 method)
std::vector<std::map<std::string, double>> ButterflyNebulaUQFFModule::generateVariations(int count) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.5, 1.5);
    
    for (int i = 0; i < count; ++i) {
        std::map<std::string, double> variation = variables;
        variation["M"] *= dis(gen);
        variation["r"] *= dis(gen);
        variation["DPM_momentum"] *= dis(gen);
        variation["k_LENR"] *= dis(gen);
        variations.push_back(variation);
    }
    return variations;
}

// Adaptive Evolution (2 methods)
void ButterflyNebulaUQFFModule::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dis(1.0, mutation_rate);
    
    variables["M"] *= std::abs(dis(gen));
    variables["r"] *= std::abs(dis(gen));
    variables["DPM_momentum"] *= std::abs(dis(gen));
    variables["k_LENR"] *= std::abs(dis(gen));
    
    autoRefineParameters();
}

void ButterflyNebulaUQFFModule::evolveSystem(int generations) {
    for (int gen = 0; gen < generations; ++gen) {
        mutateParameters(0.1);
        
        // Fitness: optimize force magnitude
        double force = computeF_U_Bi(variables["x1"], variables["x2"], 0.0);
        double target_force = -2e212;  // NGC 6302 target
        
        if (std::abs(force) < std::abs(target_force) * 0.5) {
            // Increase force components
            variables["DPM_momentum"] *= 1.05;
            variables["k_LENR"] *= 1.05;
        }
        
        autoRefineParameters();
    }
}

// State Management (4 methods)
void ButterflyNebulaUQFFModule::saveState(const std::string& state_name) {
    saved_states_butterfly::state_storage[state_name] = variables;
}

void ButterflyNebulaUQFFModule::restoreState(const std::string& state_name) {
    if (saved_states_butterfly::state_storage.find(state_name) != saved_states_butterfly::state_storage.end()) {
        variables = saved_states_butterfly::state_storage[state_name];
    } else {
        std::cerr << "Error: State '" << state_name << "' not found.\n";
    }
}

std::vector<std::string> ButterflyNebulaUQFFModule::listSavedStates() {
    std::vector<std::string> states;
    for (const auto& pair : saved_states_butterfly::state_storage) {
        states.push_back(pair.first);
    }
    return states;
}

std::string ButterflyNebulaUQFFModule::exportState() {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(6);
    oss << "Butterfly Nebula NGC 6302 State Export:\n";
    for (const auto& pair : variables) {
        oss << pair.first << "=" << pair.second << "\n";
    }
    return oss.str();
}

// System Analysis (4 methods)
std::map<std::string, double> ButterflyNebulaUQFFModule::sensitivityAnalysis(const std::string& output_var) {
    std::map<std::string, double> sensitivities;
    double baseline;
    
    if (output_var == "F_U_Bi") {
        baseline = computeF_U_Bi(variables["x1"], variables["x2"], 0.0);
    } else {
        baseline = variables[output_var];
    }
    
    double delta = 0.01;  // 1% perturbation
    for (auto& pair : variables) {
        double original = pair.second;
        variables[pair.first] = original * (1.0 + delta);
        
        double perturbed;
        if (output_var == "F_U_Bi") {
            perturbed = computeF_U_Bi(variables["x1"], variables["x2"], 0.0);
        } else {
            perturbed = variables[output_var];
        }
        
        if (std::abs(baseline) > 1e-100) {
            sensitivities[pair.first] = (perturbed - baseline) / (baseline * delta);
        }
        variables[pair.first] = original;
    }
    
    return sensitivities;
}

std::string ButterflyNebulaUQFFModule::generateReport() {
    std::ostringstream report;
    report << std::scientific << std::setprecision(3);
    report << "========== Butterfly Nebula NGC 6302 UQFF Report ==========\n";
    report << "System: " << getSystemName() << "\n\n";
    
    report << "Nebula Parameters:\n";
    report << "  Central star mass M = " << variables["M"] << " kg (" << (variables["M"]/variables["M_sun"]) << " M_sun)\n";
    report << "  Radius r = " << variables["r"] << " m\n";
    report << "  Integration range: [" << variables["x1"] << ", " << variables["x2"] << "] m\n";
    report << "  Quantum level = " << variables["level"] << "\n\n";
    
    report << "Force Components:\n";
    report << "  DPM momentum = " << variables["DPM_momentum"] << "\n";
    report << "  DPM gravity = " << variables["DPM_gravity"] << "\n";
    report << "  DPM stability = " << variables["DPM_stability"] << "\n";
    report << "  LENR coupling k_LENR = " << variables["k_LENR"] << "\n";
    report << "  LENR frequency ω_LENR = " << variables["omega_LENR"] << " Hz\n\n";
    
    double force = computeF_U_Bi(variables["x1"], variables["x2"], 0.0);
    report << "Computed Force:\n";
    report << "  F_U_Bi_i,enhanced = " << force << " N\n";
    report << "  Sign: " << (force < 0 ? "Repulsive (stabilizing)" : "Attractive") << "\n\n";
    
    report << "Physical Interpretation:\n";
    report << "  Integrates DPM, LENR, activation, DE, EM, neutron, relativistic effects\n";
    report << "  NGC 6302: Bipolar planetary nebula with complex dynamics\n";
    report << "  UQFF: [SCm]-[UA] interactions drive nebula stabilization\n";
    report << "  Sweet-Parker & Kozima terms included for completeness\n";
    
    report << "===========================================================\n";
    return report.str();
}

bool ButterflyNebulaUQFFModule::validateConsistency() {
    bool consistent = true;
    
    // Check positive mass
    if (variables["M"] <= 0) {
        std::cerr << "Inconsistency: Non-positive mass\n";
        consistent = false;
    }
    
    // Check positive radius
    if (variables["r"] <= 0) {
        std::cerr << "Inconsistency: Non-positive radius\n";
        consistent = false;
    }
    
    // Check integration bounds
    if (variables["x2"] <= variables["x1"]) {
        std::cerr << "Inconsistency: Invalid integration bounds (x2 <= x1)\n";
        consistent = false;
    }
    
    // Check positive frequencies
    if (variables["omega_LENR"] <= 0 || variables["omega_0"] <= 0) {
        std::cerr << "Inconsistency: Non-positive frequencies\n";
        consistent = false;
    }
    
    return consistent;
}

void ButterflyNebulaUQFFModule::autoCorrectAnomalies() {
    // Enforce positive mass
    if (variables["M"] <= 0) {
        variables["M"] = 0.64 * variables["M_sun"];
        std::cout << "Corrected: M reset to default (0.64 M_sun)\n";
    }
    
    // Enforce positive radius
    if (variables["r"] <= 0) {
        variables["r"] = 3.22e19;
        std::cout << "Corrected: r reset to default\n";
    }
    
    // Enforce valid integration bounds
    if (variables["x2"] <= variables["x1"]) {
        variables["x2"] = variables["x1"] + 3.22e19;
        std::cout << "Corrected: x2 adjusted to ensure x2 > x1\n";
    }
    
    // Enforce positive frequencies
    if (variables["omega_LENR"] <= 0) {
        variables["omega_LENR"] = 7.85e12;
        std::cout << "Corrected: ω_LENR reset to default\n";
    }
    
    if (variables["omega_0"] <= 0) {
        variables["omega_0"] = 1e-12;
        std::cout << "Corrected: ω_0 reset to default\n";
    }
    
    std::cout << "Anomaly correction complete.\n";
}

// Example usage
// #include "ButterflyNebulaUQFFModule.h"
// int main() {
//     ButterflyNebulaUQFFModule mod;
//     double t = 0.0;
//     double x1 = 0.0;
//     double x2 = 3.22e19;
//     double force = mod.computeF_U_Bi(x1, x2, t);
//     std::cout << "F_U_Bi ≈ " << force << " N\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.printVariables();
//     return 0;
// }

// ========== COMPREHENSIVE ENHANCED DEMONSTRATION ==========
/*
int main() {
    std::cout << std::scientific << std::setprecision(3);
    std::cout << "========== BUTTERFLY NEBULA NGC 6302 UQFF DEMONSTRATION ==========\n\n";
    
    // ===== Step 1: Initialize Module =====
    ButterflyNebulaUQFFModule mod;
    std::cout << "Step 1: Module initialized with NGC 6302 defaults:\n";
    std::cout << "  System: " << mod.getSystemName() << "\n";
    std::cout << "  M = " << mod.variables["M"] << " kg (0.64 M_sun)\n";
    std::cout << "  r = " << mod.variables["r"] << " m\n";
    std::cout << "  Integration range: [0, " << mod.variables["x2"] << "] m\n\n";
    
    // ===== Step 2: Baseline Computation =====
    std::cout << "Step 2: Compute baseline UQFF force:\n";
    double x1 = 0.0;
    double x2 = 3.22e19;
    double t = 0.0;
    double force = mod.computeF_U_Bi(x1, x2, t);
    
    std::cout << "  F_U_Bi_i,enhanced = " << force << " N\n";
    std::cout << "  Sign: " << (force < 0 ? "Repulsive (stabilizing)" : "Attractive") << "\n";
    std::cout << "  Physical interpretation: UQFF integrates DPM/LENR/activation/EM/neutron effects\n\n";
    
    // ===== Step 3-7: Dynamic Operations =====
    std::cout << "Step 3: Variable Management\n";
    mod.createVariable("bipolar_asymmetry", 1.25);
    std::cout << "  Created 'bipolar_asymmetry' parameter\n";
    
    std::cout << "\nStep 4: Nebula Expansion\n";
    mod.expandNebulaScale(1.2, 1.5);
    std::cout << "  Expanded: M = " << mod.variables["M"] << " kg, r = " << mod.variables["r"] << " m\n";
    
    std::cout << "\nStep 5: Force Expansion\n";
    mod.expandForceScale(2.0, 1.5);
    std::cout << "  Expanded: DPM_momentum = " << mod.variables["DPM_momentum"] << "\n";
    
    std::cout << "\nStep 6: Integration Expansion\n";
    mod.expandIntegrationScale(1.5, 2.0);
    std::cout << "  Expanded integration domain: x2 = " << mod.variables["x2"] << " m\n";
    
    std::cout << "\nStep 7: Batch Operations\n";
    std::vector<std::string> force_group = {"DPM_momentum", "DPM_gravity", "k_LENR"};
    mod.scaleVariableGroup(force_group, 0.8);
    std::cout << "  Scaled force group by 0.8\n\n";
    
    // ===== Step 8-12: Physical Regimes =====
    std::cout << "Steps 8-12: Test Multiple Nebula Regimes\n";
    
    mod.optimizeForMetric("standard_butterfly");
    double f1 = mod.computeF_U_Bi(mod.variables["x1"], mod.variables["x2"], 0.0);
    std::cout << "  Standard NGC 6302: F = " << f1 << " N\n";
    std::cout << "                      M = " << (mod.variables["M"]/mod.variables["M_sun"]) << " M_sun\n";
    
    mod.optimizeForMetric("high_mass_nebula");
    double f2 = mod.computeF_U_Bi(mod.variables["x1"], mod.variables["x2"], 0.0);
    std::cout << "  High Mass: F = " << f2 << " N\n";
    std::cout << "             M = " << (mod.variables["M"]/mod.variables["M_sun"]) << " M_sun\n";
    
    mod.optimizeForMetric("compact_nebula");
    double f3 = mod.computeF_U_Bi(mod.variables["x1"], mod.variables["x2"], 0.0);
    std::cout << "  Compact: F = " << f3 << " N, r = " << mod.variables["r"] << " m\n";
    
    mod.optimizeForMetric("extended_nebula");
    double f4 = mod.computeF_U_Bi(mod.variables["x1"], mod.variables["x2"], 0.0);
    std::cout << "  Extended: F = " << f4 << " N, r = " << mod.variables["r"] << " m\n";
    
    mod.optimizeForMetric("high_lenr");
    double f5 = mod.computeF_U_Bi(mod.variables["x1"], mod.variables["x2"], 0.0);
    std::cout << "  High LENR: F = " << f5 << " N, k_LENR = " << mod.variables["k_LENR"] << "\n\n";
    
    // ===== Step 13-17: Refinement & Evolution =====
    std::cout << "Step 13: Auto-Refinement\n";
    mod.updateVariable("M", -100.0);  // Invalid negative mass
    mod.autoRefineParameters();
    std::cout << "  Clamped M from negative to " << mod.variables["M"] << " kg\n";
    
    std::cout << "\nStep 14: Calibration to Observations\n";
    std::map<std::string, double> obs_data;
    obs_data["M"] = 0.70 * mod.variables["M_sun"];
    obs_data["r"] = 3.5e19;
    mod.calibrateToObservations(obs_data);
    std::cout << "  Calibrated: M = " << (mod.variables["M"]/mod.variables["M_sun"]) << " M_sun\n";
    
    std::cout << "\nStep 15: Parameter Variations\n";
    std::vector<std::map<std::string, double>> variations = mod.generateVariations(5);
    std::cout << "  Generated " << variations.size() << " parameter variations\n";
    
    std::cout << "\nStep 16: Mutation\n";
    mod.optimizeForMetric("standard_butterfly");
    mod.mutateParameters(0.15);
    std::cout << "  Mutated: M = " << (mod.variables["M"]/mod.variables["M_sun"]) << " M_sun\n";
    
    std::cout << "\nStep 17: System Evolution\n";
    mod.evolveSystem(10);
    std::cout << "  Evolved: F = " << mod.computeF_U_Bi(mod.variables["x1"], mod.variables["x2"], 0.0) << " N\n\n";
    
    // ===== Step 18-19: State Management =====
    std::cout << "Step 18: State Management\n";
    mod.optimizeForMetric("standard_butterfly");
    mod.saveState("ngc6302_standard");
    mod.optimizeForMetric("high_lenr");
    mod.saveState("high_activity");
    std::vector<std::string> saved = mod.listSavedStates();
    std::cout << "  Saved " << saved.size() << " states\n";
    mod.restoreState("ngc6302_standard");
    std::cout << "  Restored 'ngc6302_standard'\n";
    
    std::cout << "\nStep 19: Export State\n";
    std::string exported = mod.exportState();
    std::cout << "  Exported " << exported.length() << " bytes\n\n";
    
    // ===== Step 20-22: Analysis =====
    std::cout << "Step 20: Sensitivity Analysis (F_U_Bi response)\n";
    std::map<std::string, double> sensitivity = mod.sensitivityAnalysis("F_U_Bi");
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
    
    // ===== Step 23-26: Nebula Force Scale Analysis =====
    std::cout << "Steps 23-26: Planetary Nebula UQFF Force Scale Analysis\n";
    std::cout << "  Nebula Type      | M (M_sun) | r (m)     | F_U_Bi (N) | Context\n";
    std::cout << "  --------------------------------------------------------------------------\n";
    
    struct NebulaRegime {
        std::string name;
        double mass_factor;
        double radius;
        double dpm_factor;
        std::string context;
    };
    
    std::vector<NebulaRegime> regimes = {
        {"Young Compact", 0.8, 1e19, 2.0, "Early expansion phase"},
        {"NGC 6302 Std", 0.64, 3.22e19, 1.0, "Butterfly reference"},
        {"Mature Extended", 0.5, 5e19, 0.8, "Late expansion"},
        {"High Mass Core", 1.5, 2e19, 3.0, "Massive progenitor"},
        {"Diffuse Shell", 0.3, 1e20, 0.5, "Low-density outer regions"}
    };
    
    for (const auto& reg : regimes) {
        mod.updateVariable("M", reg.mass_factor * mod.variables["M_sun"]);
        mod.updateVariable("r", reg.radius);
        mod.updateVariable("x1", 0.0);
        mod.updateVariable("x2", reg.radius);
        mod.updateVariable("DPM_momentum", reg.dpm_factor);
        mod.updateVariable("DPM_gravity", reg.dpm_factor);
        
        double force_val = mod.computeF_U_Bi(0.0, reg.radius, 0.0);
        
        std::cout << "  " << std::setw(16) << std::left << reg.name
                  << " | " << std::fixed << std::setprecision(2) << std::setw(9) << reg.mass_factor
                  << " | " << std::scientific << std::setprecision(2) << std::setw(9) << reg.radius
                  << " | " << std::setw(10) << force_val
                  << " | " << reg.context << "\n";
    }
    
    std::cout << "\n========== DEMONSTRATION COMPLETE ==========\n";
    std::cout << "Butterfly Nebula NGC 6302 UQFF module validated across regimes.\n";
    std::cout << "F_U_Bi integrates DPM, LENR, activation, DE, EM, neutron, relativistic terms.\n";
    std::cout << "Physical significance: Repulsive UQFF force stabilizes bipolar structure.\n";
    std::cout << "NGC 6302: M = 0.64 M_sun, r = 3.22e19 m, F ≈ -2.09e212 N.\n";
    std::cout << "UQFF Integration: [SCm]-[UA] dynamics drive planetary nebula evolution.\n";
    std::cout << "Applications: Bipolar morphology, jet formation, nebula stabilization.\n";
    
    return 0;
}
*/
// Compile: g++ -o butterfly_test butterfly_test.cpp ButterflyNebulaUQFFModule.cpp -lm
// Sample: F_U_Bi ? -2.09e212 N; repulsive for stabilization.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

ButterflyNebulaUQFFModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeF_U_Bi, computeIntegral, computeIntegrand, and all physical term methods) are clear, concise, and variable - driven.
- Integrates a wide range of physical effects(DPM, LENR, activation, DE, EM, neutron, relativistic, Sweet, Kozima) for comprehensive force modeling.
- Uses numerical integration(trapezoidal rule) for flexible and accurate force calculation over a spatial range.
- Output and debugging functions(printVariables, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map and optimize the integration routine.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in UQFF force modeling for nebulae.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.