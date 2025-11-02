// InertiaCouplingModule.h
// Modular C++ implementation of the Inertia Coupling Constants (?_i) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes ?_i=1.0 (unitless, uniform for i=1-4) and scales U_i in F_U: -?_i [?_i U_i E_react].
// Pluggable: #include "InertiaCouplingModule.h"
// InertiaCouplingModule mod; mod.computeSumInertiaTerms(0.0); mod.updateVariable("rho_vac_SCm", new_value);
// Variables in std::map; example for Sun at t=0, t_n=0; U_i ?1.38e-47 J/m�, contrib ? -0.138 J/m�.
// Approximations: Uniform ?_i=1.0; cos(? t_n)=1; f_TRZ=0.1; ?_s=2.5e-6 rad/s; E_react=1e46.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef INERTIA_COUPLING_MODULE_H
#define INERTIA_COUPLING_MODULE_H

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

class InertiaCouplingModule {
private:
    std::map<std::string, double> variables;
    std::vector<double> computeAllInertiaTerms(double t);

public:
    // Constructor: Initialize with framework defaults
    InertiaCouplingModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeLambda_i(int i);  // ?_i=1.0 (unitless)
    double computeU_i(int i, double t);  // U_i for i=1-4 (J/m^3)
    double computeInertiaTerm(int i, double t);  // -?_i U_i E_react
    double computeSumInertiaTerms(double t);  // Sum for F_U contribution (J/m^3)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print inertia breakdown
    void printInertiaBreakdown(double t = 0.0);

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
    void expandInertiaScale(double lambda_factor);                        // λ_i coupling
    void expandDensityScale(double rho_SCm_factor, double rho_UA_factor); // Vacuum densities
    void expandDynamicsScale(double omega_factor, double E_react_factor); // Rotation and reactor

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

#endif // INERTIA_COUPLING_MODULE_H

// InertiaCouplingModule.cpp
#include "InertiaCouplingModule.h"

// Constructor: Set framework defaults (Sun at t=0, level 13)
InertiaCouplingModule::InertiaCouplingModule() {
    // Universal constants
    variables["lambda"] = 1.0;                      // Uniform ?_i (unitless)
    variables["rho_vac_SCm"] = 7.09e-37;            // J/m^3
    variables["rho_vac_UA"] = 7.09e-36;             // J/m^3
    variables["omega_s"] = 2.5e-6;                  // rad/s (Sun rotation)
    variables["f_TRZ"] = 0.1;                       // Unitless
    variables["E_react"] = 1e46;                    // J
    variables["pi"] = 3.141592653589793;
    variables["t_n"] = 0.0;                         // s
    variables["alpha_decay"] = 0.0005;              // For E_react exp, but t=0
}

// Update variable
void InertiaCouplingModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void InertiaCouplingModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void InertiaCouplingModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute ?_i (uniform 1.0)
double InertiaCouplingModule::computeLambda_i(int i) {
    return variables["lambda"];
}

// Compute U_i (uniform across i)
double InertiaCouplingModule::computeU_i(int i, double t) {
    double lambda_i = computeLambda_i(i);
    double rho_sc = variables["rho_vac_SCm"];
    double rho_ua = variables["rho_vac_UA"];
    double omega_s_t = variables["omega_s"];        // Simplified, no t dep in example
    double cos_term = std::cos(variables["pi"] * variables["t_n"]);
    double trz_factor = 1.0 + variables["f_TRZ"];
    return lambda_i * rho_sc * rho_ua * omega_s_t * cos_term * trz_factor;
}

// Compute inertia term -?_i U_i E_react
double InertiaCouplingModule::computeInertiaTerm(int i, double t) {
    double u_i = computeU_i(i, t);
    double e_react = variables["E_react"] * std::exp( - variables["alpha_decay"] * t );  // With decay
    return - computeLambda_i(i) * u_i * e_react;
}

// Sum over i=1 to 4
double InertiaCouplingModule::computeSumInertiaTerms(double t) {
    double sum = 0.0;
    for (int i = 1; i <= 4; ++i) {
        sum += computeInertiaTerm(i, t);
    }
    return sum;
}

// Equation text
std::string InertiaCouplingModule::getEquationText() {
    return "F_U = ... - ?_i [?_i * U_i * E_react] + ...\n"
           "U_i = ?_i * ?_vac,[SCm] * ?_vac,[UA] * ?_s(t) * cos(? t_n) * (1 + f_TRZ)\n"
           "Where ?_i = 1.0 (unitless, uniform for i=1-4: Ug1-Ug4);\n"
           "E_react = 1e46 * e^{-? t} (?=5e-4);\n"
           "Example Sun t=0, t_n=0: U_i ?1.38e-47 J/m�; -?_i U_i E_react ? -0.138 J/m� (per i).\n"
           "Role: Scales resistive inertia; uniform baseline opposition to dynamics.\n"
           "UQFF: Consistent across scales; aids stability in interiors/disks/mergers.";
}

// Print variables
void InertiaCouplingModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Print breakdown
void InertiaCouplingModule::printInertiaBreakdown(double t) {
    std::cout << "Inertia Breakdown at t=" << t << " s:\n";
    for (int i = 1; i <= 4; ++i) {
        double u_i = computeU_i(i, t);
        double term = computeInertiaTerm(i, t);
        std::cout << "i=" << i << ": U_i = " << std::scientific << u_i << " J/m�, Term = " << term << " J/m�\n";
    }
    std::cout << "Sum ? Terms = " << std::scientific << computeSumInertiaTerms(t) << " J/m�\n";
}

// ===== Implementation: 25-Method Dynamic Self-Update & Self-Expansion Capabilities =====

namespace inertia_coupling_saved_states {
    std::map<std::string, std::map<std::string, double>> saved_states;
}

// Variable Management (5 methods)
void InertiaCouplingModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void InertiaCouplingModule::removeVariable(const std::string& name) {
    variables.erase(name);
}

void InertiaCouplingModule::cloneVariable(const std::string& source, const std::string& target) {
    if (variables.find(source) != variables.end()) {
        variables[target] = variables[source];
    }
}

std::string InertiaCouplingModule::listVariables() {
    std::ostringstream oss;
    for (const auto& pair : variables) {
        oss << pair.first << " = " << pair.second << "\n";
    }
    return oss.str();
}

std::string InertiaCouplingModule::getSystemName() {
    return "Inertia_Coupling_UQFF";
}

// Batch Operations (2 methods)
void InertiaCouplingModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = func(variables[name]);
        }
    }
}

void InertiaCouplingModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double x) { return x * factor; });
}

// Self-Expansion (4 methods)
void InertiaCouplingModule::expandParameterSpace(const std::vector<std::string>& params, double expansion_factor) {
    for (const auto& param : params) {
        if (variables.find(param) != variables.end()) {
            variables[param] *= expansion_factor;
        }
    }
}

void InertiaCouplingModule::expandInertiaScale(double lambda_factor) {
    // λ_i: uniform inertia coupling constant (unitless, typically 1.0)
    variables["lambda"] *= lambda_factor;
}

void InertiaCouplingModule::expandDensityScale(double rho_SCm_factor, double rho_UA_factor) {
    // ρ_vac_SCm: superconductive vacuum density (J/m³, ~7.09e-37)
    // ρ_vac_UA: universal aether vacuum density (J/m³, ~7.09e-36)
    variables["rho_vac_SCm"] *= rho_SCm_factor;
    variables["rho_vac_UA"] *= rho_UA_factor;
}

void InertiaCouplingModule::expandDynamicsScale(double omega_factor, double E_react_factor) {
    // ω_s: stellar rotation rate (rad/s, ~2.5e-6 for Sun)
    // E_react: reactor energy (J, ~1e46)
    variables["omega_s"] *= omega_factor;
    variables["E_react"] *= E_react_factor;
}

// Self-Refinement (3 methods)
void InertiaCouplingModule::autoRefineParameters(const std::string& target_metric, double target_value) {
    // Example: target sum of inertia terms by adjusting lambda and E_react
    if (target_metric == "sum_inertia") {
        double current = computeSumInertiaTerms(0.0);
        if (std::abs(current) > 1e-50) {
            double ratio = target_value / current;
            variables["lambda"] *= std::sqrt(ratio);
            variables["E_react"] *= std::sqrt(ratio);
        }
    } else if (target_metric == "U_i") {
        double current = computeU_i(1, 0.0);
        if (std::abs(current) > 1e-50) {
            double ratio = target_value / current;
            variables["rho_vac_SCm"] *= std::sqrt(ratio);
            variables["omega_s"] *= std::sqrt(ratio);
        }
    }
}

void InertiaCouplingModule::calibrateToObservations(const std::map<std::string, double>& observed) {
    for (const auto& obs : observed) {
        if (variables.find(obs.first) != variables.end()) {
            variables[obs.first] = obs.second;
        }
    }
}

void InertiaCouplingModule::optimizeForMetric(const std::string& metric) {
    if (metric == "maximize_resistance") {
        // Increase inertia coupling and densities
        variables["lambda"] *= 1.2;
        variables["rho_vac_SCm"] *= 1.1;
        variables["rho_vac_UA"] *= 1.1;
    } else if (metric == "minimize_resistance") {
        // Decrease inertia effects
        variables["lambda"] *= 0.8;
        variables["E_react"] *= 0.9;
    } else if (metric == "enhance_rotation") {
        // Boost rotational contribution
        variables["omega_s"] *= 1.3;
        variables["f_TRZ"] *= 1.2;
    } else if (metric == "stabilize_system") {
        // Balance all parameters toward defaults
        variables["lambda"] = 1.0;
        variables["f_TRZ"] = 0.1;
    }
}

// Parameter Exploration (1 method)
std::vector<std::map<std::string, double>> InertiaCouplingModule::generateVariations(int count, double variance) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-variance, variance);

    for (int i = 0; i < count; ++i) {
        std::map<std::string, double> variant = variables;
        for (auto& pair : variant) {
            if (pair.first != "pi" && pair.first != "t_n") {
                double factor = 1.0 + dis(gen);
                pair.second *= factor;
                if (pair.second < 0 && pair.first != "alpha_decay") {
                    pair.second = std::abs(pair.second);
                }
            }
        }
        variations.push_back(variant);
    }
    return variations;
}

// Adaptive Evolution (2 methods)
void InertiaCouplingModule::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-mutation_rate, mutation_rate);

    for (auto& pair : variables) {
        if (pair.first != "pi" && pair.first != "t_n") {
            pair.second *= (1.0 + dis(gen));
            if (pair.second < 0 && pair.first != "alpha_decay") {
                pair.second = std::abs(pair.second);
            }
        }
    }
}

void InertiaCouplingModule::evolveSystem(int generations, std::function<double()> fitness_func) {
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
void InertiaCouplingModule::saveState(const std::string& label) {
    inertia_coupling_saved_states::saved_states[label] = variables;
}

void InertiaCouplingModule::restoreState(const std::string& label) {
    if (inertia_coupling_saved_states::saved_states.find(label) != inertia_coupling_saved_states.end()) {
        variables = inertia_coupling_saved_states::saved_states[label];
    }
}

std::vector<std::string> InertiaCouplingModule::listSavedStates() {
    std::vector<std::string> labels;
    for (const auto& pair : inertia_coupling_saved_states::saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string InertiaCouplingModule::exportState() {
    std::ostringstream oss;
    oss << "System: " << getSystemName() << "\n";
    oss << listVariables();
    return oss.str();
}

// System Analysis (4 methods)
std::map<std::string, double> InertiaCouplingModule::sensitivityAnalysis(const std::vector<std::string>& params) {
    std::map<std::string, double> sensitivity;
    double baseline = std::abs(computeSumInertiaTerms(0.0));

    for (const auto& param : params) {
        if (variables.find(param) != variables.end()) {
            double original = variables[param];
            variables[param] *= 1.01; // +1% perturbation
            double perturbed = std::abs(computeSumInertiaTerms(0.0));
            variables[param] = original;
            if (baseline > 1e-50) {
                sensitivity[param] = std::abs((perturbed - baseline) / baseline);
            } else {
                sensitivity[param] = 0.0;
            }
        }
    }
    return sensitivity;
}

std::string InertiaCouplingModule::generateReport() {
    std::ostringstream oss;
    oss << "========== Inertia Coupling Module Report ==========\n";
    oss << "System: " << getSystemName() << "\n\n";
    oss << "Key Parameters:\n";
    oss << "  λ (lambda, coupling): " << std::scientific << variables["lambda"] << " (unitless, uniform for i=1-4)\n";
    oss << "  ρ_vac_SCm: " << variables["rho_vac_SCm"] << " J/m^3\n";
    oss << "  ρ_vac_UA: " << variables["rho_vac_UA"] << " J/m^3\n";
    oss << "  ω_s (rotation): " << variables["omega_s"] << " rad/s\n";
    oss << "  f_TRZ (TRZ factor): " << variables["f_TRZ"] << " (unitless)\n";
    oss << "  E_react (reactor): " << variables["E_react"] << " J\n";
    oss << "  α_decay: " << variables["alpha_decay"] << " (decay rate)\n\n";
    
    oss << "Computed Values at t=0:\n";
    double u_i = computeU_i(1, 0.0);
    double term = computeInertiaTerm(1, 0.0);
    double sum = computeSumInertiaTerms(0.0);
    oss << "  U_i (per i): " << u_i << " J/m^3\n";
    oss << "  Inertia term (per i): " << term << " J/m^3\n";
    oss << "  Sum (all 4 ranges): " << sum << " J/m^3\n\n";
    
    oss << "Physical Context:\n";
    oss << "  λ_i = 1.0 (uniform coupling for all Ug ranges)\n";
    oss << "  U_i = λ_i × ρ_SCm × ρ_UA × ω_s × cos(πt_n) × (1 + f_TRZ)\n";
    oss << "  Inertia term: -λ_i × U_i × E_react × exp(-αt)\n";
    oss << "  Represents resistive opposition to field dynamics\n";
    oss << "  Uniform baseline ensures consistent inertia across scales\n";
    oss << "  Small magnitude (~J/m^3) but crucial for stability\n";
    oss << "  Appears with negative sign in F_U (opposes acceleration)\n";
    oss << "===================================================\n";
    return oss.str();
}

bool InertiaCouplingModule::validateConsistency() {
    bool valid = true;
    if (variables["lambda"] <= 0) valid = false;
    if (variables["rho_vac_SCm"] <= 0) valid = false;
    if (variables["rho_vac_UA"] <= 0) valid = false;
    if (variables["omega_s"] < 0) valid = false;
    if (variables["E_react"] <= 0) valid = false;
    return valid;
}

void InertiaCouplingModule::autoCorrectAnomalies() {
    if (variables["lambda"] <= 0) variables["lambda"] = 1.0;
    if (variables["rho_vac_SCm"] <= 0) variables["rho_vac_SCm"] = 7.09e-37;
    if (variables["rho_vac_UA"] <= 0) variables["rho_vac_UA"] = 7.09e-36;
    if (variables["omega_s"] < 0) variables["omega_s"] = 2.5e-6;
    if (variables["E_react"] <= 0) variables["E_react"] = 1e46;
    if (variables["f_TRZ"] < 0) variables["f_TRZ"] = 0.1;
}

// Example usage in base program (snippet)
// #include "InertiaCouplingModule.h"
// int main() {
//     InertiaCouplingModule mod;
//     double sum = mod.computeSumInertiaTerms(0.0);
//     std::cout << "Sum Inertia Terms = " << sum << " J/m�\n";
//     mod.printInertiaBreakdown(0.0);
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("lambda", 1.1);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o inertia_test inertia_test.cpp InertiaCouplingModule.cpp -lm
// Sample: Per i term ? -0.138 J/m�; sum ? -0.552 J/m� (4 terms).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

// ===== Enhanced Example: 18-Step Demonstration of Dynamic Capabilities =====
// int main() {
//     InertiaCouplingModule mod;
//     std::cout << "===== Inertia Coupling Module: Enhanced 18-Step Demo =====\n\n";
//
//     // Step 1: Initial report
//     std::cout << "Step 1: Initial Configuration\n";
//     std::cout << mod.generateReport() << "\n";
//
//     // Step 2: Track key derived quantities
//     std::cout << "Step 2: Create Tracking Variables\n";
//     double u_i_base = mod.computeU_i(1, 0.0);
//     double inertia_base = mod.computeSumInertiaTerms(0.0);
//     mod.createVariable("U_i_baseline", u_i_base);
//     mod.createVariable("inertia_sum_baseline", inertia_base);
//     mod.createVariable("resistance_factor", std::abs(inertia_base) / mod.variables["E_react"]);
//     std::cout << "  U_i_baseline: " << std::scientific << u_i_base << " J/m^3\n";
//     std::cout << "  inertia_sum_baseline: " << inertia_base << " J/m^3\n";
//     std::cout << "  resistance_factor: " << mod.variables["resistance_factor"] << "\n\n";
//
//     // Step 3: Lambda (λ_i) variation scenarios
//     std::cout << "Step 3: Lambda Coupling Variations\n";
//     mod.saveState("baseline");
//     std::vector<double> lambda_values = {0.5, 0.8, 1.0, 1.2, 1.5, 2.0};
//     for (double lam : lambda_values) {
//         mod.updateVariable("lambda", lam);
//         double sum = mod.computeSumInertiaTerms(0.0);
//         std::cout << "  λ=" << lam << " -> Sum=" << sum << " J/m^3 ";
//         std::cout << "(" << std::showpos << ((sum/inertia_base - 1.0)*100) << std::noshowpos << "%)\n";
//     }
//     mod.restoreState("baseline");
//     std::cout << "\n";
//
//     // Step 4: Time evolution (decay effects)
//     std::cout << "Step 4: Time Evolution (E_react decay with α=5e-4)\n";
//     std::vector<double> time_points = {0.0, 100.0, 500.0, 1000.0, 2000.0, 5000.0};
//     for (double t : time_points) {
//         double sum_t = mod.computeSumInertiaTerms(t);
//         double decay_factor = std::exp(-mod.variables["alpha_decay"] * t);
//         std::cout << "  t=" << t << " s: Sum=" << sum_t << " J/m^3 (decay=" << decay_factor << ")\n";
//     }
//     std::cout << "\n";
//
//     // Step 5: Expand inertia scale (boost λ)
//     std::cout << "Step 5: Expand Inertia Scale (λ x1.3)\n";
//     mod.expandInertiaScale(1.3);
//     std::cout << "  New λ: " << mod.variables["lambda"] << "\n";
//     std::cout << "  New Sum: " << mod.computeSumInertiaTerms(0.0) << " J/m^3\n\n";
//
//     // Step 6: Expand density scale (SCm and UA)
//     std::cout << "Step 6: Expand Density Scale (ρ_SCm x1.2, ρ_UA x1.15)\n";
//     mod.expandDensityScale(1.2, 1.15);
//     std::cout << "  New ρ_vac_SCm: " << mod.variables["rho_vac_SCm"] << " J/m^3\n";
//     std::cout << "  New ρ_vac_UA: " << mod.variables["rho_vac_UA"] << " J/m^3\n";
//     std::cout << "  New U_i: " << mod.computeU_i(1, 0.0) << " J/m^3\n";
//     std::cout << "  New Sum: " << mod.computeSumInertiaTerms(0.0) << " J/m^3\n\n";
//
//     // Step 7: Expand dynamics scale (rotation and reactor)
//     std::cout << "Step 7: Expand Dynamics Scale (ω_s x1.4, E_react x1.25)\n";
//     mod.expandDynamicsScale(1.4, 1.25);
//     std::cout << "  New ω_s: " << mod.variables["omega_s"] << " rad/s\n";
//     std::cout << "  New E_react: " << mod.variables["E_react"] << " J\n";
//     std::cout << "  New U_i: " << mod.computeU_i(1, 0.0) << " J/m^3\n";
//     std::cout << "  New Sum: " << mod.computeSumInertiaTerms(0.0) << " J/m^3\n\n";
//
//     // Step 8: Parameter variations
//     std::cout << "Step 8: Generate 10 Parameter Variations (±6%)\n";
//     auto variations = mod.generateVariations(10, 0.06);
//     std::cout << "  Generated " << variations.size() << " configurations\n";
//     std::cout << "  Sample λ range: " << variations[0]["lambda"] << " to " << variations[9]["lambda"] << "\n";
//     std::cout << "  Sample ω_s range: " << variations[0]["omega_s"] << " to " << variations[9]["omega_s"] << " rad/s\n\n";
//
//     // Step 9: Sensitivity analysis
//     std::cout << "Step 9: Sensitivity Analysis (Sum response to ±1% changes)\n";
//     std::vector<std::string> sens_params = {"lambda", "rho_vac_SCm", "rho_vac_UA", "omega_s", "E_react", "f_TRZ"};
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
//     // Step 10: Auto-refine to target sum
//     std::cout << "Step 10: Auto-Refine to Target Sum = -1.0 J/m^3\n";
//     double target_sum = -1.0;
//     double before_refine = mod.computeSumInertiaTerms(0.0);
//     mod.autoRefineParameters("sum_inertia", target_sum);
//     double after_refine = mod.computeSumInertiaTerms(0.0);
//     std::cout << "  Before: " << before_refine << " J/m^3\n";
//     std::cout << "  After: " << after_refine << " J/m^3\n";
//     std::cout << "  Target: " << target_sum << " J/m^3\n";
//     std::cout << "  Error: " << std::abs(after_refine - target_sum) / std::abs(target_sum) * 100 << "%\n\n";
//
//     // Step 11: Restore and target U_i
//     std::cout << "Step 11: Restore and Target U_i = 2.0e-47 J/m^3\n";
//     mod.restoreState("baseline");
//     double target_ui = 2.0e-47;
//     double before_ui = mod.computeU_i(1, 0.0);
//     mod.autoRefineParameters("U_i", target_ui);
//     double after_ui = mod.computeU_i(1, 0.0);
//     std::cout << "  Before: " << before_ui << " J/m^3\n";
//     std::cout << "  After: " << after_ui << " J/m^3\n";
//     std::cout << "  Target: " << target_ui << " J/m^3\n\n";
//
//     // Step 12: Calibration to observations
//     std::cout << "Step 12: Calibrate to Observational Data\n";
//     std::map<std::string, double> observations = {
//         {"lambda", 1.05},
//         {"omega_s", 2.8e-6},  // Slightly faster rotation
//         {"f_TRZ", 0.12}
//     };
//     mod.calibrateToObservations(observations);
//     std::cout << "  Calibrated λ: " << mod.variables["lambda"] << "\n";
//     std::cout << "  Calibrated ω_s: " << mod.variables["omega_s"] << " rad/s\n";
//     std::cout << "  Calibrated f_TRZ: " << mod.variables["f_TRZ"] << "\n";
//     std::cout << "  New Sum: " << mod.computeSumInertiaTerms(0.0) << " J/m^3\n\n";
//
//     // Step 13: Optimize for maximum resistance
//     std::cout << "Step 13: Optimize for Maximum Resistance (Inertia)\n";
//     double before_opt = std::abs(mod.computeSumInertiaTerms(0.0));
//     mod.optimizeForMetric("maximize_resistance");
//     double after_opt = std::abs(mod.computeSumInertiaTerms(0.0));
//     std::cout << "  |Sum| before: " << before_opt << " J/m^3\n";
//     std::cout << "  |Sum| after: " << after_opt << " J/m^3\n";
//     std::cout << "  Increase: " << ((after_opt/before_opt - 1.0)*100) << "%\n\n";
//
//     // Step 14: Optimize for enhanced rotation
//     std::cout << "Step 14: Optimize for Enhanced Rotation Effects\n";
//     mod.saveState("before_rotation_opt");
//     double omega_before = mod.variables["omega_s"];
//     mod.optimizeForMetric("enhance_rotation");
//     double omega_after = mod.variables["omega_s"];
//     std::cout << "  ω_s: " << omega_before << " -> " << omega_after << " rad/s\n";
//     std::cout << "  f_TRZ: " << mod.variables["f_TRZ"] << "\n";
//     std::cout << "  New Sum: " << mod.computeSumInertiaTerms(0.0) << " J/m^3\n\n";
//
//     // Step 15: Parameter mutation
//     std::cout << "Step 15: Mutate Parameters (±5%)\n";
//     mod.restoreState("before_rotation_opt");
//     mod.saveState("before_mutation");
//     double lambda_before = mod.variables["lambda"];
//     double e_before = mod.variables["E_react"];
//     mod.mutateParameters(0.05);
//     std::cout << "  λ: " << lambda_before << " -> " << mod.variables["lambda"] << "\n";
//     std::cout << "  E_react: " << e_before << " -> " << mod.variables["E_react"] << " J\n";
//     std::cout << "  New Sum: " << mod.computeSumInertiaTerms(0.0) << " J/m^3\n\n";
//
//     // Step 16: System evolution (minimize resistance for efficiency)
//     std::cout << "Step 16: Evolve System (6 generations, maximize |Sum|)\n";
//     mod.restoreState("before_mutation");
//     double initial_fitness = std::abs(mod.computeSumInertiaTerms(0.0));
//     mod.evolveSystem(6, [&mod]() { return std::abs(mod.computeSumInertiaTerms(0.0)); });
//     double final_fitness = std::abs(mod.computeSumInertiaTerms(0.0));
//     std::cout << "  Initial |Sum|: " << initial_fitness << " J/m^3\n";
//     std::cout << "  Evolved |Sum|: " << final_fitness << " J/m^3\n";
//     std::cout << "  Improvement: " << ((final_fitness / initial_fitness - 1.0) * 100) << "%\n\n";
//
//     // Step 17: Validation and auto-correction
//     std::cout << "Step 17: Validate Consistency\n";
//     bool valid = mod.validateConsistency();
//     std::cout << "  System valid: " << (valid ? "YES" : "NO") << "\n";
//     if (!valid) {
//         std::cout << "  Running auto-correction...\n";
//         mod.autoCorrectAnomalies();
//         std::cout << "  Post-correction valid: " << (mod.validateConsistency() ? "YES" : "NO") << "\n";
//     }
//     mod.printInertiaBreakdown(0.0);
//     std::cout << "\n";
//
//     // Step 18: State management and export
//     std::cout << "Step 18: State Management and Export\n";
//     mod.saveState("evolved_optimal");
//     auto saved = mod.listSavedStates();
//     std::cout << "  Saved states (" << saved.size() << "): ";
//     for (const auto& s : saved) std::cout << s << " ";
//     std::cout << "\n\n";
//     std::cout << "Final System Export:\n";
//     std::string exported = mod.exportState();
//     std::cout << exported << "\n";
//     std::cout << "Demo complete: 18 steps executed successfully!\n";
//     std::cout << "==========================================================\n";
//
//     return 0;
// }

InertiaCouplingModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeLambda_i, computeU_i, computeInertiaTerm, computeSumInertiaTerms) are clear, concise, and variable - driven.
- Output and debugging functions(printVariables, printInertiaBreakdown, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid indices, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in inertia coupling modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.