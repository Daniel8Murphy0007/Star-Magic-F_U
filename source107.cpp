// PiConstantModule.h
// Modular C++ implementation of the Mathematical Constant Pi (?) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes ? ?3.14159 (unitless) and its role in oscillatory terms like cos(? t_n), sin(?_c t), with ?_c=2? / period.
// Pluggable: #include "PiConstantModule.h"
// PiConstantModule mod; mod.computeCosPiTn(0.0); mod.updateVariable("t_n", new_value);
// Variables in std::map; examples for U_m ?_j and U_g1 cos(? t_n) at t=0, t_n=0.
// Approximations: ?=3.141592653589793; sin(?_c * 0)=0; cos(? * 0)=1.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef PI_CONSTANT_MODULE_H
#define PI_CONSTANT_MODULE_H

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

class PiConstantModule {
private:
    std::map<std::string, double> variables;
    double computeCosPiTn(double t_n);
    double computeSinOmegaCT(double t);
    double computeMuJExample(double t);
    double computeUg1CosTerm(double t_n);

public:
    // Constructor: Initialize with framework defaults
    PiConstantModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computePi();  // ?3.141592653589793 (unitless)
    double computeCosPiTn(double t_n);  // cos(? t_n)
    double computeSinOmegaCT(double t);  // sin(?_c t), ?_c=2? / period
    double computeMuJExample(double t);  // Example ?_j with sin(?_c t)
    double computeUg1CosTerm(double t_n);  // cos(? t_n) in U_g1

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

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
    void expandParameterSpace(double pi_scale, double period_scale, double amplitude_scale);
    void expandOscillationScale(double omega_factor, double period_factor);   // Scale ω_c and period
    void expandMagneticScale(double mu_factor, double B_factor);              // Scale μ_j components
    void expandPhaseScale(double t_factor, double tn_factor);                 // Scale time variables

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

#endif // PI_CONSTANT_MODULE_H

// PiConstantModule.cpp
#include "PiConstantModule.h"

// Constructor: Set framework defaults
PiConstantModule::PiConstantModule() {
    // Mathematical constants
    variables["pi"] = 3.141592653589793;            // Unitless
    variables["t_n"] = 0.0;                         // Negative time factor
    variables["t"] = 0.0;                           // Time
    variables["period"] = 3.96e8;                   // s (example solar cycle)
    variables["omega_c"] = 2.0 * variables["pi"] / variables["period"];  // rad/s
    variables["base_mu"] = 3.38e20;                 // T�m^3
    variables["B_j"] = 1e3;                         // Base T
}

// Update variable
void PiConstantModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
        if (name == "period") {
            variables["omega_c"] = 2.0 * variables["pi"] / value;
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void PiConstantModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        if (name == "period") {
            variables["omega_c"] = 2.0 * variables["pi"] / variables[name];
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void PiConstantModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute ?
double PiConstantModule::computePi() {
    return variables["pi"];
}

// Compute cos(? t_n)
double PiConstantModule::computeCosPiTn(double t_n) {
    variables["t_n"] = t_n;
    return std::cos(computePi() * t_n);
}

// Compute sin(?_c t)
double PiConstantModule::computeSinOmegaCT(double t) {
    variables["t"] = t;
    return std::sin(variables["omega_c"] * t);
}

// Example ?_j = (10^3 + 0.4 sin(?_c t)) * 3.38e20 T�m^3
double PiConstantModule::computeMuJExample(double t) {
    double sin_omega = computeSinOmegaCT(t);
    double b_j = variables["B_j"] + 0.4 * sin_omega;
    return b_j * variables["base_mu"];
}

// Example cos(? t_n) in U_g1
double PiConstantModule::computeUg1CosTerm(double t_n) {
    return computeCosPiTn(t_n);
}

// Equation text
std::string PiConstantModule::getEquationText() {
    return "? ? 3.141592653589793 (unitless mathematical constant).\n"
           "Role: Defines periodicity in oscillations; C=2? r; trig args (sin/cos with 2? cycle).\n"
           "In U_m: ?_j = (10^3 + 0.4 sin(?_c t)) * 3.38e20; ?_c = 2? / period.\n"
           "In U_g1: ... cos(? t_n) ... (time-reversal oscillations).\n"
           "Example t=0, t_n=0: sin(?_c t)=0 ? ?_j=3.38e23 T�m^3; cos(? t_n)=1.\n"
           "UQFF: Ensures cyclic/TRZ dynamics; solar cycles, rotations in nebulae/quasars.";
}

// Print variables
void PiConstantModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// ===== ENHANCED METHODS IMPLEMENTATION =====

namespace pi_constant_saved_states {
    static std::map<std::string, std::map<std::string, double>> saved_states;
}

// Variable Management
void PiConstantModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void PiConstantModule::removeVariable(const std::string& name) {
    variables.erase(name);
}

void PiConstantModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
    }
}

std::vector<std::string> PiConstantModule::listVariables() const {
    std::vector<std::string> names;
    for (const auto& pair : variables) {
        names.push_back(pair.first);
    }
    return names;
}

std::string PiConstantModule::getSystemName() const {
    return "Pi_Constant_UQFF";
}

// Batch Operations
void PiConstantModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = func(variables[name]);
        }
    }
}

void PiConstantModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// Self-Expansion: Domain-Specific Scales
void PiConstantModule::expandParameterSpace(double pi_scale, double period_scale, double amplitude_scale) {
    // Note: pi typically stays constant (mathematical), but allow for framework flexibility
    variables["pi"] *= pi_scale;
    variables["period"] *= period_scale;
    variables["omega_c"] = 2.0 * variables["pi"] / variables["period"];
    variables["base_mu"] *= amplitude_scale;
}

void PiConstantModule::expandOscillationScale(double omega_factor, double period_factor) {
    variables["period"] *= period_factor;
    variables["omega_c"] = 2.0 * variables["pi"] / variables["period"];
    variables["omega_c"] *= omega_factor;  // Additional adjustment beyond period
}

void PiConstantModule::expandMagneticScale(double mu_factor, double B_factor) {
    variables["base_mu"] *= mu_factor;
    variables["B_j"] *= B_factor;
}

void PiConstantModule::expandPhaseScale(double t_factor, double tn_factor) {
    variables["t"] *= t_factor;
    variables["t_n"] *= tn_factor;
}

// Self-Refinement
void PiConstantModule::autoRefineParameters(const std::string& target, double goal) {
    if (target == "cos_pi_tn") {
        // Target specific cos(π t_n) by adjusting t_n
        if (std::abs(goal) <= 1.0) {
            variables["t_n"] = std::acos(goal) / variables["pi"];
        }
    } else if (target == "sin_omega_t") {
        // Target specific sin(ω_c t) by adjusting t
        if (std::abs(goal) <= 1.0) {
            variables["t"] = std::asin(goal) / variables["omega_c"];
        }
    } else if (target == "mu_j") {
        // Target specific μ_j by scaling base_mu
        double current_mu = computeMuJExample(variables["t"]);
        if (std::abs(current_mu) > 1e-9) {
            variables["base_mu"] *= (goal / current_mu);
        }
    } else if (target == "omega_c") {
        // Target specific ω_c by adjusting period
        variables["omega_c"] = goal;
        variables["period"] = 2.0 * variables["pi"] / goal;
    } else if (target == "period") {
        // Target specific period
        variables["period"] = goal;
        variables["omega_c"] = 2.0 * variables["pi"] / goal;
    }
}

void PiConstantModule::calibrateToObservations(const std::map<std::string, double>& observations) {
    for (const auto& obs : observations) {
        if (variables.find(obs.first) != variables.end()) {
            variables[obs.first] = obs.second;
            // Update dependent variables
            if (obs.first == "period") {
                variables["omega_c"] = 2.0 * variables["pi"] / obs.second;
            }
        }
    }
}

void PiConstantModule::optimizeForMetric(const std::string& metric) {
    if (metric == "maximize_mu_j") {
        // Maximize μ_j: increase base_mu and B_j
        variables["base_mu"] *= 1.5;
        variables["B_j"] *= 1.3;
    } else if (metric == "faster_oscillation") {
        // Increase oscillation frequency (decrease period)
        variables["period"] *= 0.7;
        variables["omega_c"] = 2.0 * variables["pi"] / variables["period"];
    } else if (metric == "slower_oscillation") {
        // Decrease oscillation frequency (increase period)
        variables["period"] *= 1.5;
        variables["omega_c"] = 2.0 * variables["pi"] / variables["period"];
    } else if (metric == "solar_cycle_alignment") {
        // Align to 11-year solar cycle
        variables["period"] = 11.0 * 365.25 * 86400.0;  // 11 years in seconds
        variables["omega_c"] = 2.0 * variables["pi"] / variables["period"];
    }
}

// Parameter Exploration
std::vector<std::map<std::string, double>> PiConstantModule::generateVariations(int count, double variation_pct) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - variation_pct, 1.0 + variation_pct);
    
    for (int i = 0; i < count; ++i) {
        std::map<std::string, double> variant = variables;
        for (auto& pair : variant) {
            // Don't vary pi (mathematical constant) significantly
            if (pair.first == "pi") {
                pair.second *= (1.0 + dis(gen) * 0.001);  // Tiny variation for numerical exploration
            } else {
                pair.second *= dis(gen);
            }
        }
        // Recalculate omega_c for consistency
        variant["omega_c"] = 2.0 * variant["pi"] / variant["period"];
        variations.push_back(variant);
    }
    return variations;
}

// Adaptive Evolution
void PiConstantModule::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - mutation_rate, 1.0 + mutation_rate);
    
    for (auto& pair : variables) {
        if (pair.first == "pi") {
            pair.second *= (1.0 + dis(gen) * 0.001 - 0.0005);  // Minimal mutation
        } else if (pair.first != "omega_c") {  // Don't mutate omega_c directly
            pair.second *= dis(gen);
        }
    }
    // Recalculate omega_c
    variables["omega_c"] = 2.0 * variables["pi"] / variables["period"];
}

void PiConstantModule::evolveSystem(int generations, std::function<double()> fitness_func) {
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
void PiConstantModule::saveState(const std::string& label) {
    pi_constant_saved_states::saved_states[label] = variables;
}

void PiConstantModule::restoreState(const std::string& label) {
    if (pi_constant_saved_states::saved_states.find(label) != pi_constant_saved_states::saved_states.end()) {
        variables = pi_constant_saved_states::saved_states[label];
    }
}

std::vector<std::string> PiConstantModule::listSavedStates() const {
    std::vector<std::string> labels;
    for (const auto& pair : pi_constant_saved_states::saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string PiConstantModule::exportState() const {
    std::ostringstream oss;
    oss << "PiConstant_State_Export:\n";
    for (const auto& pair : variables) {
        oss << pair.first << "=" << std::scientific << pair.second << "\n";
    }
    return oss.str();
}

// System Analysis
std::map<std::string, double> PiConstantModule::sensitivityAnalysis(const std::vector<std::string>& params) {
    std::map<std::string, double> sensitivities;
    double baseline_mu = computeMuJExample(variables["t"]);
    
    for (const auto& param : params) {
        if (variables.find(param) != variables.end()) {
            double original = variables[param];
            variables[param] = original * 1.01;
            
            // Recalculate omega_c if period changed
            if (param == "period") {
                variables["omega_c"] = 2.0 * variables["pi"] / variables[param];
            }
            
            double perturbed_mu = computeMuJExample(variables["t"]);
            sensitivities[param] = (perturbed_mu - baseline_mu) / baseline_mu;
            
            // Restore
            variables[param] = original;
            if (param == "period") {
                variables["omega_c"] = 2.0 * variables["pi"] / original;
            }
        }
    }
    return sensitivities;
}

std::string PiConstantModule::generateReport() const {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(6);
    oss << "===== Pi Constant Module Report =====\n";
    oss << "System: " << getSystemName() << "\n\n";
    
    oss << "Mathematical Constant:\n";
    oss << "  π = " << variables.at("pi") << " (unitless)\n";
    oss << "  Exact value: 3.141592653589793\n";
    double pi_error = std::abs(variables.at("pi") - 3.141592653589793);
    oss << "  Deviation: " << pi_error << " (" << (pi_error / 3.141592653589793 * 100) << "%)\n\n";
    
    oss << "Oscillation Parameters:\n";
    oss << "  Period = " << variables.at("period") << " s\n";
    double period_years = variables.at("period") / (365.25 * 86400.0);
    oss << "  Period = " << period_years << " years\n";
    oss << "  ω_c = 2π/T = " << variables.at("omega_c") << " rad/s\n";
    double freq_hz = variables.at("omega_c") / (2.0 * variables.at("pi"));
    oss << "  Frequency = " << freq_hz << " Hz\n\n";
    
    oss << "Time Variables:\n";
    oss << "  t = " << variables.at("t") << " s\n";
    oss << "  t_n = " << variables.at("t_n") << " (negative time factor)\n\n";
    
    oss << "Trigonometric Functions:\n";
    double cos_pi_tn = std::cos(variables.at("pi") * variables.at("t_n"));
    double sin_omega_t = std::sin(variables.at("omega_c") * variables.at("t"));
    oss << "  cos(π t_n) = " << cos_pi_tn << "\n";
    oss << "  sin(ω_c t) = " << sin_omega_t << "\n";
    double phase_pi = variables.at("pi") * variables.at("t_n");
    double phase_omega = variables.at("omega_c") * variables.at("t");
    oss << "  Phase (π t_n): " << phase_pi << " rad = " << (phase_pi * 180.0 / variables.at("pi")) << "°\n";
    oss << "  Phase (ω_c t): " << phase_omega << " rad = " << (phase_omega * 180.0 / variables.at("pi")) << "°\n\n";
    
    oss << "Magnetic Moment Parameters:\n";
    oss << "  base_μ = " << variables.at("base_mu") << " T·m³\n";
    oss << "  B_j (base) = " << variables.at("B_j") << " T\n";
    double b_j_current = variables.at("B_j") + 0.4 * sin_omega_t;
    double mu_j_current = b_j_current * variables.at("base_mu");
    oss << "  B_j (current) = B_j + 0.4·sin(ω_c t) = " << b_j_current << " T\n";
    oss << "  μ_j (current) = " << mu_j_current << " T·m³\n\n";
    
    oss << "Physical Interpretation:\n";
    oss << "  π role: Defines 2π periodicity in oscillations\n";
    oss << "  cos(π t_n): Time-reversal oscillations in U_g1, U_g2, etc.\n";
    oss << "  sin(ω_c t): Magnetic moment modulation in μ_j\n";
    oss << "  Period ≈ " << period_years << " years: ";
    if (period_years > 10.0 && period_years < 12.0) {
        oss << "Near solar cycle (11 years)\n";
    } else if (period_years > 200.0 && period_years < 250.0) {
        oss << "Galactic orbital scale\n";
    } else {
        oss << "Custom oscillation timescale\n";
    }
    oss << "  Applications: Solar cycles, stellar rotation, nebular dynamics, quasar jets\n";
    
    return oss.str();
}

bool PiConstantModule::validateConsistency() const {
    bool valid = true;
    
    // Check pi is close to mathematical constant
    double pi_val = variables.at("pi");
    if (std::abs(pi_val - 3.141592653589793) > 0.01) {
        std::cerr << "Warning: π deviates significantly from 3.14159... (current: " << pi_val << ")\n";
    }
    
    // Check period is positive
    if (variables.find("period") != variables.end() && variables.at("period") <= 0) {
        std::cerr << "Error: period <= 0 (period must be positive)\n";
        valid = false;
    }
    
    // Check omega_c consistency with period
    if (variables.find("omega_c") != variables.end() && variables.find("period") != variables.end()) {
        double expected_omega = 2.0 * variables.at("pi") / variables.at("period");
        double actual_omega = variables.at("omega_c");
        if (std::abs(expected_omega - actual_omega) / expected_omega > 0.01) {
            std::cerr << "Warning: ω_c inconsistent with period (expected " << expected_omega 
                      << ", got " << actual_omega << ")\n";
        }
    }
    
    // Check base_mu is positive
    if (variables.find("base_mu") != variables.end() && variables.at("base_mu") <= 0) {
        std::cerr << "Error: base_mu <= 0 (magnetic moment base must be positive)\n";
        valid = false;
    }
    
    return valid;
}

void PiConstantModule::autoCorrectAnomalies() {
    // Reset pi to mathematical constant
    if (std::abs(variables["pi"] - 3.141592653589793) > 0.1) {
        variables["pi"] = 3.141592653589793;
    }
    
    // Ensure period is positive
    if (variables["period"] <= 0) {
        variables["period"] = 3.96e8;  // Default solar cycle
    }
    
    // Recalculate omega_c for consistency
    variables["omega_c"] = 2.0 * variables["pi"] / variables["period"];
    
    // Ensure base_mu is positive
    if (variables["base_mu"] <= 0) {
        variables["base_mu"] = 3.38e20;
    }
}

// Example usage in base program (snippet)
// #include "PiConstantModule.h"
// int main() {
//     PiConstantModule mod;
//     double pi_val = mod.computePi();
//     std::cout << "? ? " << pi_val << std::endl;
//     double mu = mod.computeMuJExample(0.0);
//     std::cout << "?_j (t=0) = " << mu << " T�m^3\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("t_n", 1.0);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o pi_test pi_test.cpp PiConstantModule.cpp -lm
// Sample: ?=3.14159; ?_j?3.38e23 T�m^3; cos(?*0)=1.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

// ===== Enhanced Example: 18-Step Demonstration of Dynamic Capabilities =====
// int main() {
//     PiConstantModule mod;
//     std::cout << "===== Pi Constant Module: Enhanced 18-Step Demo =====\n\n";
//
//     // Step 1: Initial report
//     std::cout << "Step 1: Initial Configuration (π = 3.14159...)\n";
//     std::cout << mod.generateReport() << "\n";
//
//     // Step 2: Track key π-dependent quantities
//     std::cout << "Step 2: Create Tracking Variables\n";
//     mod.createVariable("pi_baseline", mod.computePi());
//     mod.createVariable("omega_baseline", mod.variables["omega_c"]);
//     mod.createVariable("period_baseline", mod.variables["period"]);
//     mod.createVariable("mu_j_baseline", mod.computeMuJExample(0.0));
//     std::cout << "  π = " << std::setprecision(15) << mod.variables["pi_baseline"] << "\n";
//     std::cout << "  ω_c = " << std::scientific << mod.variables["omega_baseline"] << " rad/s\n";
//     std::cout << "  Period = " << mod.variables["period_baseline"] << " s\n";
//     std::cout << "  μ_j (t=0) = " << mod.variables["mu_j_baseline"] << " T·m³\n\n";
//
//     // Step 3: cos(π t_n) oscillation analysis
//     std::cout << "Step 3: cos(π t_n) Oscillation Analysis\n";
//     mod.saveState("baseline");
//     std::vector<double> tn_values = {0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0};
//     for (double tn : tn_values) {
//         double cos_val = mod.computeCosPiTn(tn);
//         double phase_deg = tn * 180.0;
//         std::cout << "  t_n=" << tn << ": cos(π t_n)=" << cos_val 
//                   << ", phase=" << phase_deg << "°\n";
//     }
//     mod.restoreState("baseline");
//     std::cout << "\n";
//
//     // Step 4: sin(ω_c t) magnetic modulation
//     std::cout << "Step 4: sin(ω_c t) Magnetic Moment Modulation\n";
//     double period = mod.variables["period"];
//     std::vector<double> phase_fractions = {0.0, 0.125, 0.25, 0.375, 0.5, 0.75, 1.0};
//     for (double frac : phase_fractions) {
//         double t = frac * period;
//         double sin_val = mod.computeSinOmegaCT(t);
//         double mu_j = mod.computeMuJExample(t);
//         std::cout << "  t=" << (frac * 100) << "% period: sin(ω_c t)=" << sin_val 
//                   << ", μ_j=" << mu_j << " T·m³\n";
//     }
//     std::cout << "\n";
//
//     // Step 5: Period variations (frequency sweep)
//     std::cout << "Step 5: Period Scaling (Frequency Variations)\n";
//     std::vector<double> period_factors = {0.5, 0.8, 1.0, 1.5, 2.0, 5.0};
//     double base_period = 3.96e8;
//     for (double factor : period_factors) {
//         mod.updateVariable("period", base_period * factor);
//         double omega = mod.variables["omega_c"];
//         double freq = omega / (2.0 * mod.variables["pi"]);
//         double years = mod.variables["period"] / (365.25 * 86400.0);
//         std::cout << "  Period x" << factor << " (" << years << " years): "
//                   << "ω_c=" << omega << " rad/s, f=" << freq << " Hz\n";
//     }
//     mod.restoreState("baseline");
//     std::cout << "\n";
//
//     // Step 6: Magnetic field amplitude variations
//     std::cout << "Step 6: Magnetic Field Scaling (B_j and base_μ)\n";
//     std::vector<double> mag_factors = {0.5, 0.8, 1.0, 1.2, 1.5, 2.0};
//     for (double factor : mag_factors) {
//         mod.updateVariable("base_mu", 3.38e20 * factor);
//         mod.updateVariable("B_j", 1e3 * factor);
//         double mu_j = mod.computeMuJExample(0.0);
//         std::cout << "  Magnetic x" << factor << ": base_μ=" << mod.variables["base_mu"]
//                   << " T·m³, B_j=" << mod.variables["B_j"] << " T, μ_j=" << mu_j << " T·m³\n";
//     }
//     mod.restoreState("baseline");
//     std::cout << "\n";
//
//     // Step 7: Expand oscillation scale
//     std::cout << "Step 7: Expand Oscillation Scale (ω_c x1.2, period x0.9)\n";
//     mod.saveState("pre_oscillation");
//     double omega_before = mod.variables["omega_c"];
//     double period_before = mod.variables["period"];
//     mod.expandOscillationScale(1.2, 0.9);
//     double omega_after = mod.variables["omega_c"];
//     double period_after = mod.variables["period"];
//     std::cout << "  ω_c: " << omega_before << " -> " << omega_after << " rad/s\n";
//     std::cout << "  Period: " << period_before << " -> " << period_after << " s\n";
//     std::cout << "  Period: " << (period_before / (365.25 * 86400)) << " -> " 
//               << (period_after / (365.25 * 86400)) << " years\n\n";
//
//     // Step 8: Expand magnetic scale
//     std::cout << "Step 8: Expand Magnetic Scale (base_μ x1.3, B_j x1.15)\n";
//     mod.restoreState("pre_oscillation");
//     double mu_before = mod.variables["base_mu"];
//     double bj_before = mod.variables["B_j"];
//     mod.expandMagneticScale(1.3, 1.15);
//     std::cout << "  base_μ: " << mu_before << " -> " << mod.variables["base_mu"] << " T·m³\n";
//     std::cout << "  B_j: " << bj_before << " -> " << mod.variables["B_j"] << " T\n";
//     std::cout << "  μ_j (t=0): " << mod.computeMuJExample(0.0) << " T·m³\n\n";
//
//     // Step 9: Expand phase scale
//     std::cout << "Step 9: Expand Phase Scale (t x1.5, t_n x1.2)\n";
//     mod.updateVariable("t", 1000.0);
//     mod.updateVariable("t_n", 0.5);
//     double t_before = mod.variables["t"];
//     double tn_before = mod.variables["t_n"];
//     mod.expandPhaseScale(1.5, 1.2);
//     std::cout << "  t: " << t_before << " -> " << mod.variables["t"] << " s\n";
//     std::cout << "  t_n: " << tn_before << " -> " << mod.variables["t_n"] << "\n";
//     std::cout << "  cos(π t_n): " << mod.computeCosPiTn(mod.variables["t_n"]) << "\n\n";
//
//     // Step 10: Parameter variations
//     std::cout << "Step 10: Generate 10 Parameter Variations (±8%)\n";
//     mod.restoreState("baseline");
//     auto variations = mod.generateVariations(10, 0.08);
//     std::cout << "  Generated " << variations.size() << " configurations\n";
//     std::vector<double> omega_range, mu_range;
//     for (const auto& var : variations) {
//         omega_range.push_back(var.at("omega_c"));
//         double sin_val = std::sin(var.at("omega_c") * var.at("t"));
//         double mu = (var.at("B_j") + 0.4 * sin_val) * var.at("base_mu");
//         mu_range.push_back(mu);
//     }
//     auto omega_minmax = std::minmax_element(omega_range.begin(), omega_range.end());
//     auto mu_minmax = std::minmax_element(mu_range.begin(), mu_range.end());
//     std::cout << "  ω_c range: " << *omega_minmax.first << " to " << *omega_minmax.second << " rad/s\n";
//     std::cout << "  μ_j range: " << *mu_minmax.first << " to " << *mu_minmax.second << " T·m³\n\n";
//
//     // Step 11: Sensitivity analysis
//     std::cout << "Step 11: Sensitivity Analysis (μ_j response to ±1% changes)\n";
//     std::vector<std::string> sens_params = {"base_mu", "B_j", "period", "omega_c", "t", "pi"};
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
//     // Step 12: Auto-refine to target cos(π t_n)
//     std::cout << "Step 12: Auto-Refine to Target cos(π t_n) = 0.707 (45°)\n";
//     double cos_before = mod.computeCosPiTn(mod.variables["t_n"]);
//     mod.autoRefineParameters("cos_pi_tn", 0.707);
//     double cos_after = mod.computeCosPiTn(mod.variables["t_n"]);
//     std::cout << "  Before: cos(π t_n)=" << cos_before << "\n";
//     std::cout << "  After: cos(π t_n)=" << cos_after << "\n";
//     std::cout << "  New t_n: " << mod.variables["t_n"] << "\n";
//     std::cout << "  Phase: " << (mod.variables["t_n"] * 180.0) << "°\n\n";
//
//     // Step 13: Target specific period (solar cycle)
//     std::cout << "Step 13: Target Period = 11 Years (Solar Cycle)\n";
//     mod.restoreState("baseline");
//     double period_before = mod.variables["period"];
//     double years_before = period_before / (365.25 * 86400.0);
//     mod.autoRefineParameters("period", 11.0 * 365.25 * 86400.0);
//     double period_after = mod.variables["period"];
//     double years_after = period_after / (365.25 * 86400.0);
//     std::cout << "  Before: " << years_before << " years (" << period_before << " s)\n";
//     std::cout << "  After: " << years_after << " years (" << period_after << " s)\n";
//     std::cout << "  New ω_c: " << mod.variables["omega_c"] << " rad/s\n\n";
//
//     // Step 14: Calibration to observations
//     std::cout << "Step 14: Calibrate to Observational Data\n";
//     mod.restoreState("baseline");
//     std::map<std::string, double> observations = {
//         {"period", 12.5 * 365.25 * 86400.0},    // 12.5 year cycle
//         {"base_mu", 4.0e20},                     // Higher magnetic base
//         {"B_j", 1.2e3},                          // Stronger field
//         {"t", 5000.0}                            // Time point
//     };
//     mod.calibrateToObservations(observations);
//     std::cout << "  Calibrated period: " << (mod.variables["period"] / (365.25 * 86400.0)) << " years\n";
//     std::cout << "  Calibrated ω_c: " << mod.variables["omega_c"] << " rad/s\n";
//     std::cout << "  Calibrated base_μ: " << mod.variables["base_mu"] << " T·m³\n";
//     std::cout << "  New μ_j: " << mod.computeMuJExample(mod.variables["t"]) << " T·m³\n\n";
//
//     // Step 15: Optimize for faster oscillation
//     std::cout << "Step 15: Optimize for Faster Oscillation\n";
//     mod.restoreState("baseline");
//     double freq_before = mod.variables["omega_c"] / (2.0 * mod.variables["pi"]);
//     mod.optimizeForMetric("faster_oscillation");
//     double freq_after = mod.variables["omega_c"] / (2.0 * mod.variables["pi"]);
//     std::cout << "  Frequency: " << freq_before << " -> " << freq_after << " Hz\n";
//     std::cout << "  Period: " << (mod.variables["period"] / (365.25 * 86400.0)) << " years\n";
//     std::cout << "  Speedup: " << (freq_after / freq_before) << "x\n\n";
//
//     // Step 16: Optimize for solar cycle alignment
//     std::cout << "Step 16: Optimize for Solar Cycle Alignment (11 years)\n";
//     mod.restoreState("baseline");
//     double period_pre_solar = mod.variables["period"];
//     mod.optimizeForMetric("solar_cycle_alignment");
//     double period_post_solar = mod.variables["period"];
//     std::cout << "  Period: " << (period_pre_solar / (365.25 * 86400.0)) << " -> " 
//               << (period_post_solar / (365.25 * 86400.0)) << " years\n";
//     std::cout << "  Solar cycle period: 11.0 years\n";
//     std::cout << "  New ω_c: " << mod.variables["omega_c"] << " rad/s\n\n";
//
//     // Step 17: System evolution (maximize μ_j)
//     std::cout << "Step 17: Evolve System (8 generations, maximize μ_j)\n";
//     mod.restoreState("baseline");
//     mod.updateVariable("t", 1000.0);
//     double initial_fitness = mod.computeMuJExample(mod.variables["t"]);
//     mod.evolveSystem(8, [&mod]() { return mod.computeMuJExample(mod.variables["t"]); });
//     double final_fitness = mod.computeMuJExample(mod.variables["t"]);
//     std::cout << "  Initial μ_j: " << initial_fitness << " T·m³\n";
//     std::cout << "  Evolved μ_j: " << final_fitness << " T·m³\n";
//     std::cout << "  Improvement: " << ((final_fitness / initial_fitness - 1.0) * 100) << "%\n";
//     std::cout << "  Final period: " << (mod.variables["period"] / (365.25 * 86400.0)) << " years\n\n";
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
//     std::cout << "  Final π: " << std::setprecision(15) << mod.variables["pi"] << "\n";
//     std::cout << "  Final period: " << (mod.variables["period"] / (365.25 * 86400.0)) << " years\n";
//     std::cout << "  Final ω_c: " << std::scientific << mod.variables["omega_c"] << " rad/s\n";
//     std::cout << "  Final μ_j: " << mod.computeMuJExample(mod.variables["t"]) << " T·m³\n\n";
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

PiConstantModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computePi, computeCosPiTn, computeSinOmegaCT, computeMuJExample, computeUg1CosTerm) are clear, concise, and variable - driven.
- Handles updates to dependent variables(e.g., period, omega_c) automatically for consistency.
- Output and debugging functions(printVariables, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in mathematical constant and oscillatory modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.