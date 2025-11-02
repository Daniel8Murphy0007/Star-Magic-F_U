// NegativeTimeModule.h
// Modular C++ implementation of the Negative Time Factor (t_n) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes t_n = t - t_0 (s or days, allows t_n < 0); used in cos(? t_n) for oscillations and exp(-? t cos(? t_n)) for growth/decay.
// Pluggable: #include "NegativeTimeModule.h"
// NegativeTimeModule mod; mod.computeCosPiTn(1000.0); mod.updateVariable("t_0", new_value);
// Variables in std::map; defaults t_0=0, t=0 (t_n=0); example for U_m term with t_n negative.
// Approximations: cos even function; ?=5e-5 day^-1; at t_n=-1, exp term negative (growth).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef NEGATIVE_TIME_MODULE_H
#define NEGATIVE_TIME_MODULE_H

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

class NegativeTimeModule {
private:
    std::map<std::string, double> variables;
    double computeCosPiTn(double t_n);
    double computeExpTerm(double gamma, double t, double t_n);

public:
    // Constructor: Initialize with framework defaults
    NegativeTimeModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeT_n(double t);  // t_n = t - t_0 (s/days)
    double computeCosPiTn(double t);  // cos(? t_n)
    double computeExpTerm(double gamma, double t);  // exp(-? t cos(? t_n))
    double computeOneMinusExp(double gamma, double t);  // 1 - exp(-? t cos(? t_n))
    double computeUmExample(double t, double mu_over_rj = 2.26e10);  // Simplified U_m contrib

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print t_n effects (positive/negative)
    void printTnEffects(double t, double gamma = 5e-5);

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
    void expandParameterSpace(double time_scale, double gamma_scale, double energy_scale);
    void expandTimeScale(double t_factor, double t0_factor);          // Scale t and t_0
    void expandDecayScale(double gamma_factor, double exp_factor);     // Scale gamma (decay rate)
    void expandOscillationScale(double pi_factor, double cos_factor);  // Scale π and oscillations

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

#endif // NEGATIVE_TIME_MODULE_H

// NegativeTimeModule.cpp
#include "NegativeTimeModule.h"

// Constructor: Set framework defaults
NegativeTimeModule::NegativeTimeModule() {
    // Universal constants
    variables["t_0"] = 0.0;                         // Reference time (s/days)
    variables["t"] = 0.0;                           // Current time
    variables["gamma"] = 5e-5;                      // day^-1 (example)
    variables["pi"] = 3.141592653589793;
    variables["mu_over_rj"] = 2.26e10;              // T m^2 (example)
    variables["P_SCm"] = 1.0;                       // Normalized
    variables["E_react"] = 1e46;                    // J
    variables["heaviside_f"] = 1e11 + 1.0;          // 1 + 10^13 * 0.01
    variables["quasi_f"] = 1.01;                    // 1 + 0.01
}

// Update variable
void NegativeTimeModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void NegativeTimeModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void NegativeTimeModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute t_n = t - t_0
double NegativeTimeModule::computeT_n(double t) {
    variables["t"] = t;
    return t - variables["t_0"];
}

// Compute cos(? t_n)
double NegativeTimeModule::computeCosPiTn(double t) {
    double t_n = computeT_n(t);
    return std::cos(variables["pi"] * t_n);
}

// Compute exp(-? t cos(? t_n))
double NegativeTimeModule::computeExpTerm(double gamma, double t) {
    double cos_pi_tn = computeCosPiTn(t);
    double arg = - gamma * t * cos_pi_tn;
    return std::exp(arg);
}

// Compute 1 - exp(-? t cos(? t_n))
double NegativeTimeModule::computeOneMinusExp(double gamma, double t) {
    return 1.0 - computeExpTerm(gamma, t);
}

// Simplified U_m example contrib
double NegativeTimeModule::computeUmExample(double t, double mu_over_rj) {
    double gamma = variables["gamma"];
    double one_minus_exp = computeOneMinusExp(gamma, t);
    double phi_hat = 1.0;
    double p_scm = variables["P_SCm"];
    double e_react = variables["E_react"];
    double heaviside_f = variables["heaviside_f"];
    double quasi_f = variables["quasi_f"];
    return (mu_over_rj * one_minus_exp * phi_hat) * p_scm * e_react * heaviside_f * quasi_f;
}

// Equation text
std::string NegativeTimeModule::getEquationText() {
    return "t_n = t - t_0 (s/days, allows t_n < 0 for time-reversal);\n"
           "Used in: cos(? t_n) for oscillations; exp(-? t cos(? t_n)) for decay/growth.\n"
           "In U_m: ... (1 - exp(-? t cos(? t_n))) ...;\n"
           "Negative t_n: e.g., t_n=-1 ? cos(-?)=-1 ? exp(? t) >1 (growth, negentropic).\n"
           "Example t=1000 days, ?=5e-5 day^-1, t_0=0: 1-exp ?0.049, U_m ?1.12e66 J/m�.\n"
           "t_n=-1000: same (cos even); t_n=-1: 1-exp ? -0.051 (growth phase).\n"
           "Role: Models cyclic/TRZ dynamics; forward/reverse time in nebulae/mergers/jets.";
}

// Print variables
void NegativeTimeModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Print effects for positive/negative t_n
void NegativeTimeModule::printTnEffects(double t, double gamma) {
    double t_n_pos = computeT_n(t);  // Positive example
    double cos_pos = computeCosPiTn(t);
    double exp_pos = computeExpTerm(gamma, t);
    double one_minus_pos = computeOneMinusExp(gamma, t);
    double um_pos = computeUmExample(t);

    // Negative t_n: adjust t_0 to make t_n negative
    double orig_t0 = variables["t_0"];
    variables["t_0"] = t + 1.0;  // t_n = t - (t+1) = -1
    double t_n_neg = computeT_n(t);
    double cos_neg = computeCosPiTn(t);
    double exp_neg = computeExpTerm(gamma, t);
    double one_minus_neg = computeOneMinusExp(gamma, t);
    double um_neg = computeUmExample(t);

    variables["t_0"] = orig_t0;  // Restore

    std::cout << "t_n Effects at t=" << t << " (?=" << gamma << "):\n";
    std::cout << "Positive t_n (" << t_n_pos << "): cos(? t_n)=" << cos_pos << ", 1-exp=" << one_minus_pos << ", U_m?" << um_pos << " J/m�\n";
    std::cout << "Negative t_n (" << t_n_neg << "): cos(? t_n)=" << cos_neg << ", 1-exp=" << one_minus_neg << ", U_m?" << um_neg << " J/m�\n";
}

// ===== ENHANCED METHODS IMPLEMENTATION =====

namespace negative_time_saved_states {
    static std::map<std::string, std::map<std::string, double>> saved_states;
}

// Variable Management
void NegativeTimeModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void NegativeTimeModule::removeVariable(const std::string& name) {
    variables.erase(name);
}

void NegativeTimeModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
    }
}

std::vector<std::string> NegativeTimeModule::listVariables() const {
    std::vector<std::string> names;
    for (const auto& pair : variables) {
        names.push_back(pair.first);
    }
    return names;
}

std::string NegativeTimeModule::getSystemName() const {
    return "Negative_Time_UQFF";
}

// Batch Operations
void NegativeTimeModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = func(variables[name]);
        }
    }
}

void NegativeTimeModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// Self-Expansion: Domain-Specific Scales
void NegativeTimeModule::expandParameterSpace(double time_scale, double gamma_scale, double energy_scale) {
    variables["t"] *= time_scale;
    variables["t_0"] *= time_scale;
    variables["gamma"] *= gamma_scale;
    variables["E_react"] *= energy_scale;
}

void NegativeTimeModule::expandTimeScale(double t_factor, double t0_factor) {
    variables["t"] *= t_factor;
    variables["t_0"] *= t0_factor;
}

void NegativeTimeModule::expandDecayScale(double gamma_factor, double exp_factor) {
    variables["gamma"] *= gamma_factor;
    if (variables.find("heaviside_f") != variables.end()) {
        variables["heaviside_f"] *= exp_factor;
    }
}

void NegativeTimeModule::expandOscillationScale(double pi_factor, double cos_factor) {
    // pi_factor typically 1.0 (mathematical constant), cos_factor adjusts amplitudes
    variables["pi"] *= pi_factor;
    if (variables.find("quasi_f") != variables.end()) {
        variables["quasi_f"] *= cos_factor;
    }
}

// Self-Refinement
void NegativeTimeModule::autoRefineParameters(const std::string& target, double goal) {
    if (target == "t_n") {
        // Target specific t_n value by adjusting t_0
        double current_tn = computeT_n(variables["t"]);
        variables["t_0"] += (current_tn - goal);
    } else if (target == "cos_pi_tn") {
        // Target specific cos(π t_n) by adjusting t_0
        double desired_tn = std::acos(goal) / variables["pi"];
        variables["t_0"] = variables["t"] - desired_tn;
    } else if (target == "exp_term") {
        // Target exp(-γ t cos(π t_n)) by adjusting gamma
        double t = variables["t"];
        double cos_val = computeCosPiTn(t);
        if (std::abs(cos_val * t) > 1e-9) {
            variables["gamma"] = -std::log(goal) / (t * cos_val);
        }
    } else if (target == "U_m") {
        // Target U_m by scaling E_react
        double current_um = computeUmExample(variables["t"]);
        if (std::abs(current_um) > 1e-9) {
            variables["E_react"] *= (goal / current_um);
        }
    }
}

void NegativeTimeModule::calibrateToObservations(const std::map<std::string, double>& observations) {
    for (const auto& obs : observations) {
        if (variables.find(obs.first) != variables.end()) {
            variables[obs.first] = obs.second;
        }
    }
}

void NegativeTimeModule::optimizeForMetric(const std::string& metric) {
    if (metric == "maximize_growth") {
        // Maximize growth: minimize t_0 to make t_n positive and large
        variables["t_0"] = 0.0;
        variables["gamma"] *= 1.2;
    } else if (metric == "maximize_decay") {
        // Maximize decay: increase gamma
        variables["gamma"] *= 1.5;
    } else if (metric == "enhance_oscillation") {
        // Enhance oscillation amplitude via quasi_f
        if (variables.find("quasi_f") != variables.end()) {
            variables["quasi_f"] *= 1.3;
        }
    } else if (metric == "time_reversal") {
        // Create negative t_n scenario
        variables["t_0"] = variables["t"] + 1.0;
    }
}

// Parameter Exploration
std::vector<std::map<std::string, double>> NegativeTimeModule::generateVariations(int count, double variation_pct) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - variation_pct, 1.0 + variation_pct);
    
    for (int i = 0; i < count; ++i) {
        std::map<std::string, double> variant = variables;
        for (auto& pair : variant) {
            if (pair.first != "pi") {  // Don't vary mathematical constant
                pair.second *= dis(gen);
            }
        }
        variations.push_back(variant);
    }
    return variations;
}

// Adaptive Evolution
void NegativeTimeModule::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - mutation_rate, 1.0 + mutation_rate);
    
    for (auto& pair : variables) {
        if (pair.first != "pi") {
            pair.second *= dis(gen);
        }
    }
}

void NegativeTimeModule::evolveSystem(int generations, std::function<double()> fitness_func) {
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
void NegativeTimeModule::saveState(const std::string& label) {
    negative_time_saved_states::saved_states[label] = variables;
}

void NegativeTimeModule::restoreState(const std::string& label) {
    if (negative_time_saved_states::saved_states.find(label) != negative_time_saved_states::saved_states.end()) {
        variables = negative_time_saved_states::saved_states[label];
    }
}

std::vector<std::string> NegativeTimeModule::listSavedStates() const {
    std::vector<std::string> labels;
    for (const auto& pair : negative_time_saved_states::saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string NegativeTimeModule::exportState() const {
    std::ostringstream oss;
    oss << "NegativeTime_State_Export:\n";
    for (const auto& pair : variables) {
        oss << pair.first << "=" << std::scientific << pair.second << "\n";
    }
    return oss.str();
}

// System Analysis
std::map<std::string, double> NegativeTimeModule::sensitivityAnalysis(const std::vector<std::string>& params) {
    std::map<std::string, double> sensitivities;
    double baseline_um = computeUmExample(variables["t"]);
    
    for (const auto& param : params) {
        if (variables.find(param) != variables.end() && param != "pi") {
            double original = variables[param];
            variables[param] = original * 1.01;
            double perturbed_um = computeUmExample(variables["t"]);
            sensitivities[param] = (perturbed_um - baseline_um) / baseline_um;
            variables[param] = original;
        }
    }
    return sensitivities;
}

std::string NegativeTimeModule::generateReport() const {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(3);
    oss << "===== Negative Time Module Report =====\n";
    oss << "System: " << getSystemName() << "\n\n";
    
    oss << "Time Parameters:\n";
    oss << "  t = " << variables.at("t") << " days\n";
    oss << "  t_0 = " << variables.at("t_0") << " days\n";
    double t_n = variables.at("t") - variables.at("t_0");
    oss << "  t_n = t - t_0 = " << t_n << " days\n";
    oss << "  t_n sign: " << (t_n >= 0 ? "Positive (forward time)" : "Negative (time reversal)") << "\n\n";
    
    oss << "Decay/Oscillation Parameters:\n";
    oss << "  gamma = " << variables.at("gamma") << " day^-1\n";
    oss << "  pi = " << variables.at("pi") << "\n";
    double cos_val = std::cos(variables.at("pi") * t_n);
    oss << "  cos(π t_n) = " << cos_val << "\n\n";
    
    oss << "Exponential Terms:\n";
    double t = variables.at("t");
    double gamma = variables.at("gamma");
    double exp_arg = -gamma * t * cos_val;
    double exp_term = std::exp(exp_arg);
    double one_minus_exp = 1.0 - exp_term;
    oss << "  exp(-γ t cos(π t_n)) = " << exp_term << "\n";
    oss << "  1 - exp(...) = " << one_minus_exp << "\n";
    oss << "  Regime: " << (one_minus_exp > 0 ? "Decay/Growth" : "Enhanced Growth (negentropic)") << "\n\n";
    
    oss << "Energy Scaling:\n";
    oss << "  E_react = " << variables.at("E_react") << " J\n";
    oss << "  Heaviside factor = " << variables.at("heaviside_f") << "\n";
    oss << "  Quasi factor = " << variables.at("quasi_f") << "\n\n";
    
    oss << "U_m Contribution (example):\n";
    double mu_rj = variables.at("mu_over_rj");
    double um = mu_rj * one_minus_exp * variables.at("P_SCm") * variables.at("E_react") 
                * variables.at("heaviside_f") * variables.at("quasi_f");
    oss << "  μ_j/r_j = " << mu_rj << " T·m²\n";
    oss << "  U_m ≈ " << um << " J/m³\n\n";
    
    oss << "Physical Interpretation:\n";
    if (t_n < 0) {
        oss << "  Negative t_n: Time-reversal regime\n";
        oss << "  Applications: Nebular collapse, merger rewind, TRZ dynamics\n";
    } else {
        oss << "  Positive t_n: Forward time regime\n";
        oss << "  Applications: Stellar evolution, jet propagation, standard dynamics\n";
    }
    oss << "  cos(π t_n) oscillation: " << (cos_val > 0 ? "Positive phase" : "Negative phase") << "\n";
    oss << "  Growth/decay: " << (exp_arg < 0 ? "Growth (negentropic)" : "Decay") << "\n";
    
    return oss.str();
}

bool NegativeTimeModule::validateConsistency() const {
    bool valid = true;
    
    // Check critical parameters exist and are reasonable
    if (variables.find("gamma") != variables.end() && variables.at("gamma") < 0) {
        std::cerr << "Warning: gamma < 0 (decay rate should be positive)\n";
        valid = false;
    }
    if (variables.find("E_react") != variables.end() && variables.at("E_react") <= 0) {
        std::cerr << "Warning: E_react <= 0 (energy should be positive)\n";
        valid = false;
    }
    if (variables.find("pi") != variables.end()) {
        double pi_val = variables.at("pi");
        if (std::abs(pi_val - 3.141592653589793) > 0.01) {
            std::cerr << "Warning: pi deviates significantly from mathematical constant\n";
        }
    }
    
    return valid;
}

void NegativeTimeModule::autoCorrectAnomalies() {
    // Reset to safe defaults if anomalies detected
    if (variables["gamma"] < 0) {
        variables["gamma"] = 5e-5;
    }
    if (variables["E_react"] <= 0) {
        variables["E_react"] = 1e46;
    }
    if (std::abs(variables["pi"] - 3.141592653589793) > 0.1) {
        variables["pi"] = 3.141592653589793;
    }
}

// Example usage in base program (snippet)
// #include "NegativeTimeModule.h"
// int main() {
//     NegativeTimeModule mod;
//     double t = 1000.0;  // days
//     mod.printTnEffects(t);
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("t_0", 500.0);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o tn_test tn_test.cpp NegativeTimeModule.cpp -lm
// Sample: Positive: 1-exp?0.049; Negative t_n=-1: 1-exp?-0.051 (growth); U_m scales accordingly.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

// ===== Enhanced Example: 18-Step Demonstration of Dynamic Capabilities =====
// int main() {
//     NegativeTimeModule mod;
//     std::cout << "===== Negative Time Module: Enhanced 18-Step Demo =====\n\n";
//
//     // Step 1: Initial report
//     std::cout << "Step 1: Initial Configuration (t=0, t_0=0)\n";
//     std::cout << mod.generateReport() << "\n";
//
//     // Step 2: Track key time quantities
//     std::cout << "Step 2: Create Tracking Variables\n";
//     mod.createVariable("t_n_baseline", mod.computeT_n(mod.variables["t"]));
//     mod.createVariable("cos_baseline", mod.computeCosPiTn(mod.variables["t"]));
//     mod.createVariable("exp_baseline", mod.computeExpTerm(mod.variables["gamma"], mod.variables["t"]));
//     mod.createVariable("U_m_baseline", mod.computeUmExample(mod.variables["t"]));
//     std::cout << "  t_n = " << mod.variables["t_n_baseline"] << " days\n";
//     std::cout << "  cos(π t_n) = " << mod.variables["cos_baseline"] << "\n";
//     std::cout << "  exp term = " << mod.variables["exp_baseline"] << "\n";
//     std::cout << "  U_m = " << std::scientific << mod.variables["U_m_baseline"] << " J/m³\n\n";
//
//     // Step 3: Time evolution (positive t_n)
//     std::cout << "Step 3: Time Evolution (Forward Time, t_n > 0)\n";
//     mod.saveState("baseline");
//     std::vector<double> time_points = {0, 100, 500, 1000, 2000, 5000};
//     for (double t : time_points) {
//         double tn = mod.computeT_n(t);
//         double cos_val = mod.computeCosPiTn(t);
//         double exp_val = mod.computeExpTerm(mod.variables["gamma"], t);
//         double um = mod.computeUmExample(t);
//         std::cout << "  t=" << t << " days: t_n=" << tn << ", cos(π t_n)=" << cos_val 
//                   << ", exp=" << exp_val << ", U_m=" << um << " J/m³\n";
//     }
//     mod.restoreState("baseline");
//     std::cout << "\n";
//
//     // Step 4: Negative time scenarios (t_n < 0)
//     std::cout << "Step 4: Negative Time Scenarios (Time Reversal, t_n < 0)\n";
//     mod.updateVariable("t", 1000.0);
//     std::vector<double> t0_vals = {1001, 1100, 1500, 2000};
//     for (double t0 : t0_vals) {
//         mod.updateVariable("t_0", t0);
//         double tn = mod.computeT_n(mod.variables["t"]);
//         double cos_val = mod.computeCosPiTn(mod.variables["t"]);
//         double one_minus = mod.computeOneMinusExp(mod.variables["gamma"], mod.variables["t"]);
//         std::cout << "  t_0=" << t0 << ": t_n=" << tn << " (negative), cos(π t_n)=" << cos_val 
//                   << ", 1-exp=" << one_minus << "\n";
//     }
//     mod.restoreState("baseline");
//     std::cout << "\n";
//
//     // Step 5: Gamma variations (decay rate)
//     std::cout << "Step 5: Gamma Scaling (Decay Rate Variations)\n";
//     mod.updateVariable("t", 1000.0);
//     std::vector<double> gamma_factors = {0.5, 0.8, 1.0, 1.5, 2.0, 5.0};
//     for (double factor : gamma_factors) {
//         mod.updateVariable("gamma", 5e-5 * factor);
//         double exp_val = mod.computeExpTerm(mod.variables["gamma"], mod.variables["t"]);
//         double one_minus = mod.computeOneMinusExp(mod.variables["gamma"], mod.variables["t"]);
//         std::cout << "  γ x" << factor << " (γ=" << mod.variables["gamma"] << " day^-1): "
//                   << "exp=" << exp_val << ", 1-exp=" << one_minus << "\n";
//     }
//     mod.restoreState("baseline");
//     std::cout << "\n";
//
//     // Step 6: Oscillation phase analysis
//     std::cout << "Step 6: cos(π t_n) Oscillation Phase Analysis\n";
//     mod.updateVariable("t_0", 0.0);
//     std::vector<double> tn_phases = {-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0};
//     for (double tn : tn_phases) {
//         mod.updateVariable("t", tn);
//         double cos_val = mod.computeCosPiTn(mod.variables["t"]);
//         double phase_deg = std::acos(cos_val) * 180.0 / mod.variables["pi"];
//         std::cout << "  t_n=" << tn << ": cos(π t_n)=" << cos_val 
//                   << ", phase≈" << phase_deg << "°\n";
//     }
//     mod.restoreState("baseline");
//     std::cout << "\n";
//
//     // Step 7: Expand time scale
//     std::cout << "Step 7: Expand Time Scale (t x1.5, t_0 x1.2)\n";
//     mod.updateVariable("t", 1000.0);
//     mod.updateVariable("t_0", 100.0);
//     mod.saveState("pre_time_scale");
//     double tn_before = mod.computeT_n(mod.variables["t"]);
//     mod.expandTimeScale(1.5, 1.2);
//     double tn_after = mod.computeT_n(mod.variables["t"]);
//     std::cout << "  Before: t=" << (1000.0) << ", t_0=" << (100.0) << ", t_n=" << tn_before << "\n";
//     std::cout << "  After: t=" << mod.variables["t"] << ", t_0=" << mod.variables["t_0"] 
//               << ", t_n=" << tn_after << "\n\n";
//
//     // Step 8: Expand decay scale
//     std::cout << "Step 8: Expand Decay Scale (γ x1.3, heaviside x1.1)\n";
//     mod.restoreState("pre_time_scale");
//     double gamma_before = mod.variables["gamma"];
//     double heaviside_before = mod.variables["heaviside_f"];
//     mod.expandDecayScale(1.3, 1.1);
//     std::cout << "  γ: " << gamma_before << " -> " << mod.variables["gamma"] << " day^-1\n";
//     std::cout << "  Heaviside factor: " << heaviside_before << " -> " << mod.variables["heaviside_f"] << "\n\n";
//
//     // Step 9: Expand oscillation scale
//     std::cout << "Step 9: Expand Oscillation Scale (quasi_f x1.2)\n";
//     double quasi_before = mod.variables["quasi_f"];
//     mod.expandOscillationScale(1.0, 1.2);
//     std::cout << "  quasi_f: " << quasi_before << " -> " << mod.variables["quasi_f"] << "\n";
//     std::cout << "  New U_m: " << mod.computeUmExample(mod.variables["t"]) << " J/m³\n\n";
//
//     // Step 10: Parameter variations
//     std::cout << "Step 10: Generate 10 Parameter Variations (±10%)\n";
//     mod.restoreState("baseline");
//     mod.updateVariable("t", 1000.0);
//     auto variations = mod.generateVariations(10, 0.10);
//     std::cout << "  Generated " << variations.size() << " configurations\n";
//     std::vector<double> um_range;
//     for (const auto& var : variations) {
//         NegativeTimeModule temp_mod;
//         temp_mod.variables = var;
//         um_range.push_back(temp_mod.computeUmExample(temp_mod.variables["t"]));
//     }
//     auto um_minmax = std::minmax_element(um_range.begin(), um_range.end());
//     std::cout << "  U_m range: " << *um_minmax.first << " to " << *um_minmax.second << " J/m³\n\n";
//
//     // Step 11: Sensitivity analysis
//     std::cout << "Step 11: Sensitivity Analysis (U_m response to ±1% changes)\n";
//     mod.updateVariable("t", 1000.0);
//     std::vector<std::string> sens_params = {"t", "t_0", "gamma", "E_react", "mu_over_rj", "heaviside_f"};
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
//     // Step 12: Auto-refine to target t_n
//     std::cout << "Step 12: Auto-Refine to Target t_n = 500 days\n";
//     mod.updateVariable("t", 1000.0);
//     mod.updateVariable("t_0", 100.0);
//     double tn_before_refine = mod.computeT_n(mod.variables["t"]);
//     mod.autoRefineParameters("t_n", 500.0);
//     double tn_after_refine = mod.computeT_n(mod.variables["t"]);
//     std::cout << "  Before: t_n=" << tn_before_refine << " days\n";
//     std::cout << "  After: t_n=" << tn_after_refine << " days\n";
//     std::cout << "  New t_0: " << mod.variables["t_0"] << " days\n\n";
//
//     // Step 13: Target specific cos(π t_n)
//     std::cout << "Step 13: Target cos(π t_n) = 0.5\n";
//     mod.restoreState("baseline");
//     mod.updateVariable("t", 1000.0);
//     double cos_before = mod.computeCosPiTn(mod.variables["t"]);
//     mod.autoRefineParameters("cos_pi_tn", 0.5);
//     double cos_after = mod.computeCosPiTn(mod.variables["t"]);
//     std::cout << "  Before: cos(π t_n)=" << cos_before << "\n";
//     std::cout << "  After: cos(π t_n)=" << cos_after << "\n";
//     std::cout << "  New t_0: " << mod.variables["t_0"] << " days\n\n";
//
//     // Step 14: Calibration to observations
//     std::cout << "Step 14: Calibrate to Observational Data\n";
//     mod.restoreState("baseline");
//     std::map<std::string, double> observations = {
//         {"t", 2000.0},          // 2000 days elapsed
//         {"t_0", 500.0},         // Reference at 500 days
//         {"gamma", 8e-5},        // Faster decay
//         {"E_react", 1.2e46}     // Higher reactor energy
//     };
//     mod.calibrateToObservations(observations);
//     std::cout << "  Calibrated t: " << mod.variables["t"] << " days\n";
//     std::cout << "  Calibrated t_n: " << mod.computeT_n(mod.variables["t"]) << " days\n";
//     std::cout << "  Calibrated γ: " << mod.variables["gamma"] << " day^-1\n";
//     std::cout << "  New U_m: " << mod.computeUmExample(mod.variables["t"]) << " J/m³\n\n";
//
//     // Step 15: Optimize for maximum growth
//     std::cout << "Step 15: Optimize for Maximum Growth (maximize 1-exp)\n";
//     mod.restoreState("baseline");
//     mod.updateVariable("t", 1000.0);
//     double growth_before = mod.computeOneMinusExp(mod.variables["gamma"], mod.variables["t"]);
//     mod.optimizeForMetric("maximize_growth");
//     double growth_after = mod.computeOneMinusExp(mod.variables["gamma"], mod.variables["t"]);
//     std::cout << "  1-exp before: " << growth_before << "\n";
//     std::cout << "  1-exp after: " << growth_after << "\n";
//     std::cout << "  New t_0: " << mod.variables["t_0"] << " days\n";
//     std::cout << "  New γ: " << mod.variables["gamma"] << " day^-1\n\n";
//
//     // Step 16: Time reversal scenario
//     std::cout << "Step 16: Optimize for Time Reversal (t_n < 0)\n";
//     mod.restoreState("baseline");
//     mod.updateVariable("t", 1000.0);
//     double tn_pre_reversal = mod.computeT_n(mod.variables["t"]);
//     mod.optimizeForMetric("time_reversal");
//     double tn_post_reversal = mod.computeT_n(mod.variables["t"]);
//     std::cout << "  t_n before: " << tn_pre_reversal << " days (forward)\n";
//     std::cout << "  t_n after: " << tn_post_reversal << " days (reversal)\n";
//     std::cout << "  New t_0: " << mod.variables["t_0"] << " days\n";
//     std::cout << "  Interpretation: Negentropic/TRZ regime\n\n";
//
//     // Step 17: System evolution (maximize U_m)
//     std::cout << "Step 17: Evolve System (6 generations, maximize U_m)\n";
//     mod.restoreState("baseline");
//     mod.updateVariable("t", 1000.0);
//     double initial_fitness = mod.computeUmExample(mod.variables["t"]);
//     mod.evolveSystem(6, [&mod]() { return mod.computeUmExample(mod.variables["t"]); });
//     double final_fitness = mod.computeUmExample(mod.variables["t"]);
//     std::cout << "  Initial U_m: " << initial_fitness << " J/m³\n";
//     std::cout << "  Evolved U_m: " << final_fitness << " J/m³\n";
//     std::cout << "  Improvement: " << ((final_fitness / initial_fitness - 1.0) * 100) << "%\n";
//     std::cout << "  Final t_n: " << mod.computeT_n(mod.variables["t"]) << " days\n\n";
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
//     std::cout << "  Final t: " << mod.variables["t"] << " days\n";
//     std::cout << "  Final t_n: " << mod.computeT_n(mod.variables["t"]) << " days\n";
//     std::cout << "  Final cos(π t_n): " << mod.computeCosPiTn(mod.variables["t"]) << "\n";
//     std::cout << "  Final U_m: " << mod.computeUmExample(mod.variables["t"]) << " J/m³\n\n";
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

NegativeTimeModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeT_n, computeCosPiTn, computeExpTerm, computeOneMinusExp, computeUmExample) are clear, concise, and variable - driven.
- Handles negative time values(t_n < 0) robustly, enabling modeling of time - reversal and cyclic effects.
    - Output and debugging functions(printVariables, printTnEffects, getEquationText) provide transparency and aid validation.
    - Well - documented physical meaning and example calculations in comments and equation text.

    Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in negative time factor modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.