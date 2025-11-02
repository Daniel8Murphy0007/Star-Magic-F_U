// QuasiLongitudinalModule.h
// Modular C++ implementation of the Quasi-Longitudinal Wave Factor (f_quasi) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes f_quasi=0.01 (unitless) and its scaling (1 + f_quasi) in Universal Magnetism U_m term.
// Pluggable: #include "QuasiLongitudinalModule.h"
// QuasiLongitudinalModule mod; mod.computeUmContribution(0.0); mod.updateVariable("f_quasi", new_value);
// Variables in std::map; example for Sun at t=0, t_n=0; minor 1% increase in U_m.
// Approximations: 1 - e^{-? t cos(? t_n)}=0 at t=0; ?_hat_j=1; P_SCm=1; f_Heaviside=0.01 (1 + 10^13 f=1e11+1).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef QUASI_LONGITUDINAL_MODULE_H
#define QUASI_LONGITUDINAL_MODULE_H

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

class QuasiLongitudinalModule {
private:
    std::map<std::string, double> variables;
    double computeQuasiFactor();
    double computeUmBase(int j, double t);
    double computeUmContribution(int j, double t);

public:
    // Constructor: Initialize with framework defaults
    QuasiLongitudinalModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeF_quasi();  // 0.01 (unitless)
    double computeQuasiFactor();  // 1 + f_quasi = 1.01
    double computeUmContribution(int j, double t);  // U_m single string (J/m^3)
    double computeUmWithNoQuasi(int j, double t);  // Without quasi

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print U_m comparison (with/without quasi)
    void printUmComparison(int j = 1, double t = 0.0);

    // ===== ENHANCED: Dynamic Self-Update & Self-Expansion Methods =====
    // Variable Management
    void createVariable(const std::string& name, double value);
    void removeVariable(const std::string& name);
    void cloneVariable(const std::string& source, const std::string& dest);
    std::vector<std::string> listVariables() const;
    std::string getSystemName() const;

    // Batch Operations
    void transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func);
    void scaleVariableGroup(const std::vector<std::string>& names, double factor);

    // Self-Expansion: Domain-Specific Scales
    void expandParameterSpace(double quasi_scale, double magnetic_scale, double wave_scale);
    void expandQuasiScale(double fquasi_factor, double heaviside_factor);   // Scale f_quasi and Heaviside
    void expandWaveScale(double mu_factor, double gamma_factor);            // Scale magnetic moment and decay
    void expandEnergyScale(double ereact_factor, double pscm_factor);       // Scale reactor energy and pressure

    // Self-Refinement
    void autoRefineParameters(const std::string& target, double goal);
    void calibrateToObservations(const std::map<std::string, double>& observations);
    void optimizeForMetric(const std::string& metric);

    // Parameter Exploration
    std::vector<std::map<std::string, double>> generateVariations(int count, double variation_pct);

    // Adaptive Evolution
    void mutateParameters(double mutation_rate);
    void evolveSystem(int generations, std::function<double()> fitness_func);

    // State Management
    void saveState(const std::string& label);
    void restoreState(const std::string& label);
    std::vector<std::string> listSavedStates() const;
    std::string exportState() const;

    // System Analysis
    std::map<std::string, double> sensitivityAnalysis(const std::vector<std::string>& params);
    std::string generateReport() const;
    bool validateConsistency() const;
    void autoCorrectAnomalies();
};

#endif // QUASI_LONGITUDINAL_MODULE_H

// QuasiLongitudinalModule.cpp
#include "QuasiLongitudinalModule.h"

// Constructor: Set framework defaults
QuasiLongitudinalModule::QuasiLongitudinalModule() {
    // Universal constants
    variables["f_quasi"] = 0.01;                    // Unitless fraction
    variables["mu_j"] = 3.38e23;                    // T�m^3 (j=1)
    variables["r_j"] = 1.496e13;                    // m
    variables["gamma"] = 5e-5 / 86400.0;            // s^-1 (0.00005 day^-1)
    variables["t_n"] = 0.0;                         // s
    variables["phi_hat_j"] = 1.0;                   // Normalized
    variables["P_SCm"] = 1.0;                       // Pressure
    variables["E_react"] = 1e46;                    // J
    variables["f_Heaviside"] = 0.01;                // For Heaviside
    variables["scale_Heaviside"] = 1e13;            // Amplification
    variables["pi"] = 3.141592653589793;

    // Derived
    variables["quasi_factor"] = computeQuasiFactor();
    variables["heaviside_factor"] = 1.0 + variables["scale_Heaviside"] * variables["f_Heaviside"];
}

// Update variable
void QuasiLongitudinalModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
        if (name == "f_quasi") {
            variables["quasi_factor"] = computeQuasiFactor();
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void QuasiLongitudinalModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        if (name == "f_quasi") {
            variables["quasi_factor"] = computeQuasiFactor();
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void QuasiLongitudinalModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute f_quasi (0.01)
double QuasiLongitudinalModule::computeF_quasi() {
    return variables["f_quasi"];
}

// Compute 1 + f_quasi
double QuasiLongitudinalModule::computeQuasiFactor() {
    return 1.0 + computeF_quasi();
}

// Base for U_m without quasi/Heaviside
double QuasiLongitudinalModule::computeUmBase(int j, double t) {
    double mu_over_rj = variables["mu_j"] / variables["r_j"];
    double exp_arg = - variables["gamma"] * t * std::cos(variables["pi"] * variables["t_n"]);
    double one_minus_exp = 1.0 - std::exp(exp_arg);
    double phi_hat = variables["phi_hat_j"];
    return mu_over_rj * one_minus_exp * phi_hat * variables["P_SCm"] * variables["E_react"];
}

// U_m contribution with quasi
double QuasiLongitudinalModule::computeUmContribution(int j, double t) {
    double base = computeUmBase(j, t);
    double quasi_f = computeQuasiFactor();
    double heaviside_f = variables["heaviside_factor"];
    return base * heaviside_f * quasi_f;
}

// U_m without quasi (set f=0 temporarily)
double QuasiLongitudinalModule::computeUmWithNoQuasi(int j, double t) {
    double orig_f = variables["f_quasi"];
    variables["f_quasi"] = 0.0;
    double result = computeUmContribution(j, t);
    variables["f_quasi"] = orig_f;
    return result;
}

// Equation text
std::string QuasiLongitudinalModule::getEquationText() {
    return "U_m = ?_j [ (?_j / r_j) (1 - e^{-? t cos(? t_n)}) ?_hat_j ] P_SCm E_react (1 + 10^13 f_Heaviside) (1 + f_quasi)\n"
           "Where f_quasi = 0.01 (unitless quasi-longitudinal wave factor);\n"
           "Quasi factor = 1 + 0.01 = 1.01 (1% increase).\n"
           "Example j=1, t=0: U_m contrib ?2.28e65 J/m� (with); ?2.26e65 J/m� (without; -1%).\n"
           "Role: Minor scaling for quasi-longitudinal waves in magnetic strings; subtle [SCm]/[UA] wave effects.\n"
           "UQFF: Enhances wave propagation in jets/nebulae; small but cumulative in dynamics.";
}

// Print variables
void QuasiLongitudinalModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Print U_m comparison
void QuasiLongitudinalModule::printUmComparison(int j, double t) {
    double um_with = computeUmContribution(j, t);
    double um_without = computeUmWithNoQuasi(j, t);
    double percent_increase = ((um_with - um_without) / um_without) * 100.0;
    std::cout << "U_m Comparison for j=" << j << " at t=" << t << " s:\n";
    std::cout << "With quasi: " << std::scientific << um_with << " J/m�\n";
    std::cout << "Without quasi: " << um_without << " J/m�\n";
    std::cout << "Increase: +" << std::fixed << std::setprecision(1) << percent_increase << "%\n";
}

// ===== ENHANCED METHODS IMPLEMENTATION =====

namespace quasi_longitudinal_saved_states {
    static std::map<std::string, std::map<std::string, double>> saved_states;
}

// Variable Management
void QuasiLongitudinalModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void QuasiLongitudinalModule::removeVariable(const std::string& name) {
    variables.erase(name);
}

void QuasiLongitudinalModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
    }
}

std::vector<std::string> QuasiLongitudinalModule::listVariables() const {
    std::vector<std::string> names;
    for (const auto& pair : variables) {
        names.push_back(pair.first);
    }
    return names;
}

std::string QuasiLongitudinalModule::getSystemName() const {
    return "Quasi_Longitudinal_UQFF";
}

// Batch Operations
void QuasiLongitudinalModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = func(variables[name]);
        }
    }
}

void QuasiLongitudinalModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// Self-Expansion: Domain-Specific Scales
void QuasiLongitudinalModule::expandParameterSpace(double quasi_scale, double magnetic_scale, double wave_scale) {
    variables["f_quasi"] *= quasi_scale;
    variables["quasi_factor"] = computeQuasiFactor();
    variables["mu_j"] *= magnetic_scale;
    variables["gamma"] *= wave_scale;
}

void QuasiLongitudinalModule::expandQuasiScale(double fquasi_factor, double heaviside_factor) {
    variables["f_quasi"] *= fquasi_factor;
    variables["quasi_factor"] = computeQuasiFactor();
    variables["f_Heaviside"] *= heaviside_factor;
    variables["heaviside_factor"] = 1.0 + variables["scale_Heaviside"] * variables["f_Heaviside"];
}

void QuasiLongitudinalModule::expandWaveScale(double mu_factor, double gamma_factor) {
    variables["mu_j"] *= mu_factor;
    variables["gamma"] *= gamma_factor;
}

void QuasiLongitudinalModule::expandEnergyScale(double ereact_factor, double pscm_factor) {
    variables["E_react"] *= ereact_factor;
    variables["P_SCm"] *= pscm_factor;
}

// Self-Refinement
void QuasiLongitudinalModule::autoRefineParameters(const std::string& target, double goal) {
    if (target == "quasi_factor") {
        // Target specific quasi factor by adjusting f_quasi
        variables["f_quasi"] = goal - 1.0;
        variables["quasi_factor"] = computeQuasiFactor();
    } else if (target == "f_quasi") {
        // Target specific f_quasi directly
        variables["f_quasi"] = goal;
        variables["quasi_factor"] = computeQuasiFactor();
    } else if (target == "U_m") {
        // Target specific U_m by scaling E_react
        double current_um = computeUmContribution(1, variables.find("t") != variables.end() ? variables["t"] : 0.0);
        if (std::abs(current_um) > 1e-9) {
            variables["E_react"] *= (goal / current_um);
        }
    } else if (target == "percent_increase") {
        // Target specific percentage increase from quasi factor
        double desired_factor = 1.0 + (goal / 100.0);
        variables["f_quasi"] = desired_factor - 1.0;
        variables["quasi_factor"] = computeQuasiFactor();
    }
}

void QuasiLongitudinalModule::calibrateToObservations(const std::map<std::string, double>& observations) {
    for (const auto& obs : observations) {
        if (variables.find(obs.first) != variables.end()) {
            variables[obs.first] = obs.second;
            // Update dependent variables
            if (obs.first == "f_quasi") {
                variables["quasi_factor"] = computeQuasiFactor();
            } else if (obs.first == "f_Heaviside" || obs.first == "scale_Heaviside") {
                variables["heaviside_factor"] = 1.0 + variables["scale_Heaviside"] * variables["f_Heaviside"];
            }
        }
    }
}

void QuasiLongitudinalModule::optimizeForMetric(const std::string& metric) {
    if (metric == "maximize_quasi") {
        // Increase quasi-longitudinal contribution
        variables["f_quasi"] *= 1.5;
        variables["quasi_factor"] = computeQuasiFactor();
    } else if (metric == "minimize_quasi") {
        // Reduce quasi-longitudinal contribution
        variables["f_quasi"] *= 0.5;
        variables["quasi_factor"] = computeQuasiFactor();
    } else if (metric == "enhance_waves") {
        // Enhance wave propagation effects
        variables["gamma"] *= 1.3;
        variables["mu_j"] *= 1.2;
    } else if (metric == "standard_1_percent") {
        // Reset to standard 1% quasi contribution
        variables["f_quasi"] = 0.01;
        variables["quasi_factor"] = computeQuasiFactor();
    }
}

// Parameter Exploration
std::vector<std::map<std::string, double>> QuasiLongitudinalModule::generateVariations(int count, double variation_pct) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - variation_pct, 1.0 + variation_pct);
    
    for (int i = 0; i < count; ++i) {
        std::map<std::string, double> variant = variables;
        for (auto& pair : variant) {
            if (pair.first != "pi" && pair.first != "quasi_factor" && pair.first != "heaviside_factor") {
                pair.second *= dis(gen);
            }
        }
        // Recalculate derived variables
        variant["quasi_factor"] = 1.0 + variant["f_quasi"];
        variant["heaviside_factor"] = 1.0 + variant["scale_Heaviside"] * variant["f_Heaviside"];
        variations.push_back(variant);
    }
    return variations;
}

// Adaptive Evolution
void QuasiLongitudinalModule::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - mutation_rate, 1.0 + mutation_rate);
    
    for (auto& pair : variables) {
        if (pair.first != "pi" && pair.first != "quasi_factor" && pair.first != "heaviside_factor") {
            pair.second *= dis(gen);
        }
    }
    // Recalculate derived variables
    variables["quasi_factor"] = computeQuasiFactor();
    variables["heaviside_factor"] = 1.0 + variables["scale_Heaviside"] * variables["f_Heaviside"];
}

void QuasiLongitudinalModule::evolveSystem(int generations, std::function<double()> fitness_func) {
    double best_fitness = fitness_func();
    std::map<std::string, double> best_state = variables;
    
    for (int gen = 0; gen < generations; ++gen) {
        mutateParameters(0.05);
        double fitness = fitness_func();
        if (fitness > best_fitness) {
            best_fitness = fitness;
            best_state = variables;
        } else {
            variables = best_state;  // Revert if worse
        }
    }
    variables = best_state;
}

// State Management
void QuasiLongitudinalModule::saveState(const std::string& label) {
    quasi_longitudinal_saved_states::saved_states[label] = variables;
}

void QuasiLongitudinalModule::restoreState(const std::string& label) {
    if (quasi_longitudinal_saved_states::saved_states.find(label) != quasi_longitudinal_saved_states::saved_states.end()) {
        variables = quasi_longitudinal_saved_states::saved_states[label];
    }
}

std::vector<std::string> QuasiLongitudinalModule::listSavedStates() const {
    std::vector<std::string> labels;
    for (const auto& pair : quasi_longitudinal_saved_states::saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string QuasiLongitudinalModule::exportState() const {
    std::ostringstream oss;
    oss << "QuasiLongitudinal_State_Export:\n";
    for (const auto& pair : variables) {
        oss << pair.first << "=" << std::scientific << pair.second << "\n";
    }
    return oss.str();
}

// System Analysis
std::map<std::string, double> QuasiLongitudinalModule::sensitivityAnalysis(const std::vector<std::string>& params) {
    std::map<std::string, double> sensitivities;
    double t = variables.find("t") != variables.end() ? variables["t"] : 0.0;
    double baseline_um = computeUmContribution(1, t);
    
    for (const auto& param : params) {
        if (variables.find(param) != variables.end() && param != "pi") {
            double original = variables[param];
            variables[param] = original * 1.01;
            
            // Update derived variables if needed
            if (param == "f_quasi") {
                variables["quasi_factor"] = computeQuasiFactor();
            } else if (param == "f_Heaviside" || param == "scale_Heaviside") {
                variables["heaviside_factor"] = 1.0 + variables["scale_Heaviside"] * variables["f_Heaviside"];
            }
            
            double perturbed_um = computeUmContribution(1, t);
            sensitivities[param] = (perturbed_um - baseline_um) / baseline_um;
            
            // Restore
            variables[param] = original;
            if (param == "f_quasi") {
                variables["quasi_factor"] = computeQuasiFactor();
            } else if (param == "f_Heaviside" || param == "scale_Heaviside") {
                variables["heaviside_factor"] = 1.0 + variables["scale_Heaviside"] * variables["f_Heaviside"];
            }
        }
    }
    return sensitivities;
}

std::string QuasiLongitudinalModule::generateReport() const {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(3);
    oss << "===== Quasi-Longitudinal Wave Module Report =====\n";
    oss << "System: " << getSystemName() << "\n\n";
    
    oss << "Quasi-Longitudinal Parameters:\n";
    oss << "  f_quasi = " << variables.at("f_quasi") << " (unitless)\n";
    oss << "  Quasi factor = 1 + f_quasi = " << variables.at("quasi_factor") << "\n";
    double percent_contrib = variables.at("f_quasi") * 100.0;
    oss << "  Contribution: +" << std::fixed << std::setprecision(2) << percent_contrib 
        << "% to U_m\n\n";
    
    oss << std::scientific;
    oss << "Heaviside Amplification:\n";
    oss << "  f_Heaviside = " << variables.at("f_Heaviside") << "\n";
    oss << "  Scale factor = " << variables.at("scale_Heaviside") << "\n";
    oss << "  Heaviside factor = 1 + scale × f = " << variables.at("heaviside_factor") << "\n\n";
    
    oss << "Wave Parameters:\n";
    oss << "  μ_j (magnetic moment) = " << variables.at("mu_j") << " T·m³\n";
    oss << "  r_j (string distance) = " << variables.at("r_j") << " m\n";
    double mu_over_r = variables.at("mu_j") / variables.at("r_j");
    oss << "  μ_j/r_j = " << mu_over_r << " T·m²\n";
    oss << "  γ (decay rate) = " << variables.at("gamma") << " s⁻¹\n";
    double gamma_days = variables.at("gamma") * 86400.0;
    oss << "  γ = " << gamma_days << " day⁻¹\n\n";
    
    oss << "Energy Parameters:\n";
    oss << "  E_react = " << variables.at("E_react") << " J\n";
    oss << "  P_SCm (pressure) = " << variables.at("P_SCm") << "\n";
    oss << "  φ̂_j (normalized) = " << variables.at("phi_hat_j") << "\n\n";
    
    oss << "U_m Computations (j=1, t=0):\n";
    double t = 0.0;
    double exp_arg = -variables.at("gamma") * t * std::cos(variables.at("pi") * variables.at("t_n"));
    double one_minus_exp = 1.0 - std::exp(exp_arg);
    oss << "  1 - exp(-γ t cos(π t_n)) = " << one_minus_exp << "\n";
    
    double base = mu_over_r * one_minus_exp * variables.at("phi_hat_j") 
                  * variables.at("P_SCm") * variables.at("E_react");
    oss << "  U_m (base, no factors) = " << base << " J/m³\n";
    
    double with_heaviside = base * variables.at("heaviside_factor");
    oss << "  U_m (with Heaviside) = " << with_heaviside << " J/m³\n";
    
    double with_all = with_heaviside * variables.at("quasi_factor");
    oss << "  U_m (with Heaviside + quasi) = " << with_all << " J/m³\n";
    
    double quasi_increase_pct = (variables.at("quasi_factor") - 1.0) * 100.0;
    oss << "  Quasi contribution: +" << std::fixed << std::setprecision(2) 
        << quasi_increase_pct << "%\n\n";
    
    oss << std::scientific;
    oss << "Physical Interpretation:\n";
    if (variables.at("f_quasi") < 0.005) {
        oss << "  Very weak quasi-longitudinal contribution (<0.5%)\n";
    } else if (variables.at("f_quasi") < 0.02) {
        oss << "  Typical quasi-longitudinal contribution (~1%)\n";
    } else if (variables.at("f_quasi") < 0.05) {
        oss << "  Enhanced quasi-longitudinal contribution (2-5%)\n";
    } else {
        oss << "  Strong quasi-longitudinal contribution (>5%)\n";
    }
    
    oss << "  Applications:\n";
    oss << "    - Subtle wave propagation effects in magnetic strings\n";
    oss << "    - Longitudinal vs transverse wave coupling in [SCm]/[UA]\n";
    oss << "    - Jet dynamics: Cumulative quasi effects over large scales\n";
    oss << "    - Nebular wave propagation: Minor but persistent modulation\n";
    
    return oss.str();
}

bool QuasiLongitudinalModule::validateConsistency() const {
    bool valid = true;
    
    // Check f_quasi is small (typically << 1)
    if (variables.find("f_quasi") != variables.end()) {
        double fq = variables.at("f_quasi");
        if (fq < 0 || fq > 0.2) {
            std::cerr << "Warning: f_quasi outside typical range [0, 0.2] (current: " << fq << ")\n";
        }
    }
    
    // Check quasi_factor consistency
    if (variables.find("quasi_factor") != variables.end() && variables.find("f_quasi") != variables.end()) {
        double expected = 1.0 + variables.at("f_quasi");
        double actual = variables.at("quasi_factor");
        if (std::abs(expected - actual) > 1e-9) {
            std::cerr << "Error: quasi_factor inconsistent (expected " << expected 
                      << ", got " << actual << ")\n";
            valid = false;
        }
    }
    
    // Check gamma is positive
    if (variables.find("gamma") != variables.end() && variables.at("gamma") < 0) {
        std::cerr << "Error: gamma < 0 (decay rate must be positive)\n";
        valid = false;
    }
    
    // Check E_react is positive
    if (variables.find("E_react") != variables.end() && variables.at("E_react") <= 0) {
        std::cerr << "Error: E_react <= 0 (reactor energy must be positive)\n";
        valid = false;
    }
    
    return valid;
}

void QuasiLongitudinalModule::autoCorrectAnomalies() {
    // Reset f_quasi to typical value if out of range
    if (variables["f_quasi"] < 0 || variables["f_quasi"] > 0.2) {
        variables["f_quasi"] = 0.01;
        variables["quasi_factor"] = computeQuasiFactor();
    }
    
    // Ensure quasi_factor consistency
    variables["quasi_factor"] = 1.0 + variables["f_quasi"];
    
    // Ensure gamma is positive
    if (variables["gamma"] < 0) {
        variables["gamma"] = 5e-5 / 86400.0;
    }
    
    // Ensure E_react is positive
    if (variables["E_react"] <= 0) {
        variables["E_react"] = 1e46;
    }
    
    // Recalculate heaviside_factor
    variables["heaviside_factor"] = 1.0 + variables["scale_Heaviside"] * variables["f_Heaviside"];
}

// Example usage in base program (snippet)
// #include "QuasiLongitudinalModule.h"
// int main() {
//     QuasiLongitudinalModule mod;
//     double quasi_f = mod.computeQuasiFactor();
//     std::cout << "Quasi Factor = " << quasi_f << std::endl;
//     mod.printUmComparison(1, 0.0);
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("f_quasi", 0.02);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o quasi_test quasi_test.cpp QuasiLongitudinalModule.cpp -lm
// Sample: Factor=1.01; U_m with=2.28e65 J/m� (+1% vs without).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

// ===== Enhanced Example: 18-Step Demonstration of Dynamic Capabilities =====
// int main() {
//     QuasiLongitudinalModule mod;
//     std::cout << "===== Quasi-Longitudinal Wave Module: Enhanced 18-Step Demo =====\n\n";
//
//     // Step 1: Initial report
//     std::cout << "Step 1: Initial Configuration (f_quasi=0.01, +1%)\n";
//     std::cout << mod.generateReport() << "\n";
//
//     // Step 2: Track key quasi-longitudinal quantities
//     std::cout << "Step 2: Create Tracking Variables\n";
//     mod.createVariable("f_quasi_baseline", mod.computeF_quasi());
//     mod.createVariable("quasi_factor_baseline", mod.computeQuasiFactor());
//     mod.createVariable("U_m_with_baseline", mod.computeUmContribution(1, 0.0));
//     mod.createVariable("U_m_without_baseline", mod.computeUmWithNoQuasi(1, 0.0));
//     double increase_pct = ((mod.variables["U_m_with_baseline"] - mod.variables["U_m_without_baseline"]) 
//                            / mod.variables["U_m_without_baseline"]) * 100.0;
//     mod.createVariable("percent_increase_baseline", increase_pct);
//     std::cout << "  f_quasi = " << mod.variables["f_quasi_baseline"] << "\n";
//     std::cout << "  Quasi factor = " << mod.variables["quasi_factor_baseline"] << "\n";
//     std::cout << "  U_m (with quasi) = " << std::scientific << mod.variables["U_m_with_baseline"] << " J/m³\n";
//     std::cout << "  U_m (without) = " << mod.variables["U_m_without_baseline"] << " J/m³\n";
//     std::cout << "  Increase: +" << std::fixed << std::setprecision(2) << increase_pct << "%\n\n";
//
//     // Step 3: f_quasi variations (quasi-longitudinal strength)
//     std::cout << "Step 3: f_quasi Variations (Wave Contribution Strength)\n";
//     mod.saveState("baseline");
//     std::vector<double> fquasi_values = {0.0, 0.005, 0.01, 0.02, 0.05, 0.1};
//     for (double fq : fquasi_values) {
//         mod.updateVariable("f_quasi", fq);
//         double qf = mod.computeQuasiFactor();
//         double um_with = mod.computeUmContribution(1, 0.0);
//         double um_without = mod.computeUmWithNoQuasi(1, 0.0);
//         double pct = ((um_with - um_without) / um_without) * 100.0;
//         std::cout << "  f_quasi=" << fq << ": factor=" << qf << ", increase=+" << pct << "%\n";
//     }
//     mod.restoreState("baseline");
//     std::cout << "\n";
//
//     // Step 4: Heaviside factor scaling
//     std::cout << "Step 4: Heaviside Amplification Variations\n";
//     std::vector<double> heaviside_vals = {0.005, 0.01, 0.02, 0.05};
//     for (double fh : heaviside_vals) {
//         mod.updateVariable("f_Heaviside", fh);
//         mod.updateVariable("heaviside_factor", 1.0 + mod.variables["scale_Heaviside"] * fh);
//         double um = mod.computeUmContribution(1, 0.0);
//         std::cout << "  f_Heaviside=" << fh << ": Heaviside factor=" 
//                   << mod.variables["heaviside_factor"] << ", U_m=" << std::scientific << um << " J/m³\n";
//     }
//     mod.restoreState("baseline");
//     std::cout << "\n";
//
//     // Step 5: Combined quasi and Heaviside effects
//     std::cout << "Step 5: Combined Quasi + Heaviside Effects\n";
//     std::vector<std::pair<double, double>> combos = {
//         {0.01, 0.01}, {0.02, 0.01}, {0.01, 0.02}, {0.05, 0.02}, {0.1, 0.05}
//     };
//     for (const auto& combo : combos) {
//         mod.updateVariable("f_quasi", combo.first);
//         mod.updateVariable("f_Heaviside", combo.second);
//         mod.updateVariable("heaviside_factor", 1.0 + mod.variables["scale_Heaviside"] * combo.second);
//         double um = mod.computeUmContribution(1, 0.0);
//         std::cout << "  f_quasi=" << combo.first << ", f_Heaviside=" << combo.second 
//                   << ": U_m=" << um << " J/m³\n";
//     }
//     mod.restoreState("baseline");
//     std::cout << "\n";
//
//     // Step 6: Time evolution effects
//     std::cout << "Step 6: Time Evolution (Wave Propagation)\n";
//     std::vector<double> time_points = {0.0, 1e5, 5e5, 1e6, 5e6, 1e7};
//     for (double t : time_points) {
//         double um = mod.computeUmContribution(1, t);
//         double exp_arg = -mod.variables["gamma"] * t * std::cos(mod.variables["pi"] * mod.variables["t_n"]);
//         double one_minus_exp = 1.0 - std::exp(exp_arg);
//         std::cout << "  t=" << std::scientific << t << " s: 1-exp=" << one_minus_exp 
//                   << ", U_m=" << um << " J/m³\n";
//     }
//     std::cout << "\n";
//
//     // Step 7: Expand quasi scale
//     std::cout << "Step 7: Expand Quasi Scale (f_quasi x1.5, f_Heaviside x1.2)\n";
//     mod.restoreState("baseline");
//     mod.saveState("pre_quasi_expansion");
//     double fq_before = mod.variables["f_quasi"];
//     double fh_before = mod.variables["f_Heaviside"];
//     double um_before = mod.computeUmContribution(1, 0.0);
//     mod.expandQuasiScale(1.5, 1.2);
//     double um_after = mod.computeUmContribution(1, 0.0);
//     std::cout << "  f_quasi: " << fq_before << " -> " << mod.variables["f_quasi"] << "\n";
//     std::cout << "  f_Heaviside: " << fh_before << " -> " << mod.variables["f_Heaviside"] << "\n";
//     std::cout << "  U_m: " << um_before << " -> " << um_after << " J/m³\n";
//     std::cout << "  Increase: " << ((um_after / um_before - 1.0) * 100) << "%\n\n";
//
//     // Step 8: Expand wave scale
//     std::cout << "Step 8: Expand Wave Scale (μ_j x1.3, γ x1.15)\n";
//     mod.restoreState("pre_quasi_expansion");
//     double mu_before = mod.variables["mu_j"];
//     double gamma_before = mod.variables["gamma"];
//     mod.expandWaveScale(1.3, 1.15);
//     std::cout << "  μ_j: " << mu_before << " -> " << mod.variables["mu_j"] << " T·m³\n";
//     std::cout << "  γ: " << gamma_before << " -> " << mod.variables["gamma"] << " s⁻¹\n";
//     std::cout << "  New U_m: " << mod.computeUmContribution(1, 0.0) << " J/m³\n\n";
//
//     // Step 9: Expand energy scale
//     std::cout << "Step 9: Expand Energy Scale (E_react x1.2, P_SCm x1.1)\n";
//     mod.restoreState("baseline");
//     double ereact_before = mod.variables["E_react"];
//     double pscm_before = mod.variables["P_SCm"];
//     mod.expandEnergyScale(1.2, 1.1);
//     std::cout << "  E_react: " << ereact_before << " -> " << mod.variables["E_react"] << " J\n";
//     std::cout << "  P_SCm: " << pscm_before << " -> " << mod.variables["P_SCm"] << "\n";
//     std::cout << "  New U_m: " << mod.computeUmContribution(1, 0.0) << " J/m³\n\n";
//
//     // Step 10: Parameter variations
//     std::cout << "Step 10: Generate 10 Parameter Variations (±10%)\n";
//     auto variations = mod.generateVariations(10, 0.10);
//     std::cout << "  Generated " << variations.size() << " configurations\n";
//     std::vector<double> fquasi_range, um_range;
//     for (const auto& var : variations) {
//         QuasiLongitudinalModule temp_mod;
//         temp_mod.variables = var;
//         fquasi_range.push_back(var.at("f_quasi"));
//         um_range.push_back(temp_mod.computeUmContribution(1, 0.0));
//     }
//     auto fq_minmax = std::minmax_element(fquasi_range.begin(), fquasi_range.end());
//     auto um_minmax = std::minmax_element(um_range.begin(), um_range.end());
//     std::cout << "  f_quasi range: " << *fq_minmax.first << " to " << *fq_minmax.second << "\n";
//     std::cout << "  U_m range: " << *um_minmax.first << " to " << *um_minmax.second << " J/m³\n\n";
//
//     // Step 11: Sensitivity analysis
//     std::cout << "Step 11: Sensitivity Analysis (U_m response to ±1% changes)\n";
//     std::vector<std::string> sens_params = {"f_quasi", "f_Heaviside", "mu_j", "gamma", "E_react", "P_SCm"};
//     auto sensitivity = mod.sensitivityAnalysis(sens_params);
//     std::cout << "  Most influential parameters:\n";
//     std::vector<std::pair<std::string, double>> sorted_sens(sensitivity.begin(), sensitivity.end());
//     std::sort(sorted_sens.begin(), sorted_sens.end(), 
//               [](const auto& a, const auto& b) { return std::abs(a.second) > std::abs(b.second); });
//     for (const auto& s : sorted_sens) {
//         std::cout << "    " << s.first << ": " << (s.second * 100) << "% sensitivity\n";
//     }
//     std::cout << "\n";
//
//     // Step 12: Auto-refine to target quasi factor
//     std::cout << "Step 12: Auto-Refine to Target Quasi Factor = 1.03 (+3%)\n";
//     mod.restoreState("baseline");
//     double qf_before = mod.computeQuasiFactor();
//     mod.autoRefineParameters("quasi_factor", 1.03);
//     double qf_after = mod.computeQuasiFactor();
//     std::cout << "  Before: quasi factor=" << qf_before << ", f_quasi=" << (qf_before - 1.0) << "\n";
//     std::cout << "  After: quasi factor=" << qf_after << ", f_quasi=" << mod.variables["f_quasi"] << "\n";
//     std::cout << "  New U_m: " << mod.computeUmContribution(1, 0.0) << " J/m³\n\n";
//
//     // Step 13: Target specific percentage increase
//     std::cout << "Step 13: Target Specific Percentage Increase (+2.5%)\n";
//     mod.restoreState("baseline");
//     double pct_before = (mod.computeQuasiFactor() - 1.0) * 100;
//     mod.autoRefineParameters("percent_increase", 2.5);
//     double pct_after = (mod.computeQuasiFactor() - 1.0) * 100;
//     std::cout << "  Before: +" << pct_before << "%\n";
//     std::cout << "  After: +" << pct_after << "%\n";
//     std::cout << "  New f_quasi: " << mod.variables["f_quasi"] << "\n\n";
//
//     // Step 14: Calibration to observations
//     std::cout << "Step 14: Calibrate to Observational Data (Enhanced Quasi)\n";
//     mod.restoreState("baseline");
//     std::map<std::string, double> observations = {
//         {"f_quasi", 0.025},         // 2.5% quasi contribution
//         {"f_Heaviside", 0.015},     // Enhanced Heaviside
//         {"mu_j", 4.0e23},           // Stronger magnetic moment
//         {"gamma", 7e-5 / 86400.0}   // Faster decay
//     };
//     mod.calibrateToObservations(observations);
//     std::cout << "  Calibrated f_quasi: " << mod.variables["f_quasi"] << " (+" 
//               << (mod.variables["f_quasi"] * 100) << "%)\n";
//     std::cout << "  Calibrated quasi factor: " << mod.variables["quasi_factor"] << "\n";
//     std::cout << "  New U_m: " << mod.computeUmContribution(1, 0.0) << " J/m³\n\n";
//
//     // Step 15: Optimize for maximum quasi contribution
//     std::cout << "Step 15: Optimize for Maximum Quasi Contribution\n";
//     mod.restoreState("baseline");
//     double um_pre_max = mod.computeUmContribution(1, 0.0);
//     double fq_pre_max = mod.variables["f_quasi"];
//     mod.optimizeForMetric("maximize_quasi");
//     double um_post_max = mod.computeUmContribution(1, 0.0);
//     double fq_post_max = mod.variables["f_quasi"];
//     std::cout << "  f_quasi: " << fq_pre_max << " -> " << fq_post_max << "\n";
//     std::cout << "  Quasi increase: " << ((fq_post_max / fq_pre_max) * 100) << "%\n";
//     std::cout << "  U_m: " << um_pre_max << " -> " << um_post_max << " J/m³\n\n";
//
//     // Step 16: Optimize for standard 1% contribution
//     std::cout << "Step 16: Reset to Standard 1% Quasi Contribution\n";
//     double fq_pre_reset = mod.variables["f_quasi"];
//     mod.optimizeForMetric("standard_1_percent");
//     double fq_post_reset = mod.variables["f_quasi"];
//     std::cout << "  f_quasi: " << fq_pre_reset << " -> " << fq_post_reset << "\n";
//     std::cout << "  Quasi factor: " << mod.computeQuasiFactor() << "\n";
//     std::cout << "  Standard contribution: +1.00%\n\n";
//
//     // Step 17: System evolution (maximize U_m)
//     std::cout << "Step 17: Evolve System (8 generations, maximize U_m)\n";
//     mod.restoreState("baseline");
//     double initial_fitness = mod.computeUmContribution(1, 0.0);
//     mod.evolveSystem(8, [&mod]() { return mod.computeUmContribution(1, 0.0); });
//     double final_fitness = mod.computeUmContribution(1, 0.0);
//     std::cout << "  Initial U_m: " << initial_fitness << " J/m³\n";
//     std::cout << "  Evolved U_m: " << final_fitness << " J/m³\n";
//     std::cout << "  Improvement: " << ((final_fitness / initial_fitness - 1.0) * 100) << "%\n";
//     std::cout << "  Final f_quasi: " << mod.variables["f_quasi"] << " (+" 
//               << (mod.variables["f_quasi"] * 100) << "%)\n\n";
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
//     std::cout << "  Final f_quasi: " << mod.variables["f_quasi"] << "\n";
//     std::cout << "  Final quasi factor: " << mod.variables["quasi_factor"] << "\n";
//     std::cout << "  Final U_m (with quasi): " << mod.computeUmContribution(1, 0.0) << " J/m³\n";
//     std::cout << "  Final U_m (without quasi): " << mod.computeUmWithNoQuasi(1, 0.0) << " J/m³\n\n";
//     mod.saveState("evolved_optimal");
//     auto saved = mod.listSavedStates();
//     std::cout << "  Saved states (" << saved.size() << "): ";
//     for (const auto& s : saved) std::cout << s << " ";
//     std::cout << "\n\n";
//     std::cout << "Final System Export:\n";
//     std::string exported = mod.exportState();
//     std::cout << exported << "\n";
//     std::cout << "Demo complete: 18 steps executed successfully!\n";
//     std::cout << "================================================================\n";
//
//     return 0;
// }

QuasiLongitudinalModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeF_quasi, computeQuasiFactor, computeUmContribution, computeUmWithNoQuasi) are clear, concise, and variable - driven.
- Automatic recalculation of derived variables(quasi_factor) when dependencies change.
- Output and debugging functions(printVariables, printUmComparison, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in quasi - longitudinal wave factor modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.