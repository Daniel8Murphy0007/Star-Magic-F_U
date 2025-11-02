// ReciprocationDecayModule.h
// Modular C++ implementation of the Reciprocation Decay Rate (?) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes ?=0.00005 day?� (~5.8e-10 s?�); used in exp(-? t cos(? t_n)) for U_m decay.
// Pluggable: #include "ReciprocationDecayModule.h"
// ReciprocationDecayModule mod; mod.computeOneMinusExp(1000.0, 0.0); mod.updateVariable("gamma_day", new_value);
// Variables in std::map; example for t=1000 days, t_n=0; 1-exp ?0.049.
// Approximations: cos(? t_n)=1; timescale ~55 years; ?_j / r_j=2.26e10 T m�.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef RECIPROCATION_DECAY_MODULE_H
#define RECIPROCATION_DECAY_MODULE_H

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

class ReciprocationDecayModule {
private:
    std::map<std::string, double> variables;
    double computeGamma_s();  // ? in s?�
    double computeCosPiTn(double t_n);
    double computeExpTerm(double t_day, double t_n);
    double computeOneMinusExp(double t_day, double t_n);

public:
    // Constructor: Initialize with framework defaults
    ReciprocationDecayModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeGamma_day();  // 0.00005 day?�
    double computeGamma_s();    // ~5.8e-10 s?�
    double computeCosPiTn(double t_n);  // cos(? t_n)
    double computeExpTerm(double t_day, double t_n);  // exp(-? t cos(? t_n))
    double computeOneMinusExp(double t_day, double t_n);  // 1 - exp(...)
    double computeUmExample(double t_day, double t_n, double mu_over_rj = 2.26e10);  // Simplified U_m contrib

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print decay effects
    void printDecayEffects(double t_day = 1000.0, double t_n = 0.0);

    // ===== ENHANCED METHODS =====
    // Variable Management
    void createVariable(const std::string& name, double value);
    void removeVariable(const std::string& name);
    void cloneVariable(const std::string& source, const std::string& dest);
    std::vector<std::string> listVariables() const;
    std::string getSystemName() const;
    double getVariable(const std::string& name) const { return variables.at(name); }

    // Batch Operations
    void transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func);
    void scaleVariableGroup(const std::vector<std::string>& names, double factor);

    // Self-Expansion: Domain-Specific Scales
    void expandParameterSpace(double decay_scale, double time_scale, double energy_scale);
    void expandDecayScale(double gamma_factor, double timescale_factor);  // γ and timescales
    void expandOscillationScale(double pi_factor, double cos_factor);    // π and cos(π t_n)
    void expandMagneticScale(double mu_factor, double heaviside_factor); // μ_j/r_j and Heaviside

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

#endif // RECIPROCATION_DECAY_MODULE_H

// ReciprocationDecayModule.cpp
#include "ReciprocationDecayModule.h"

// Constructor: Set framework defaults
ReciprocationDecayModule::ReciprocationDecayModule() {
    // Universal constants
    variables["gamma_day"] = 0.00005;               // day?�
    variables["day_to_s"] = 86400.0;                // s/day
    variables["t_n"] = 0.0;                         // days
    variables["t_day"] = 0.0;                       // days
    variables["pi"] = 3.141592653589793;
    variables["mu_over_rj"] = 2.26e10;              // T m�
    variables["P_SCm"] = 1.0;                       // Normalized
    variables["E_react"] = 1e46;                    // J
    variables["heaviside_f"] = 1e11 + 1.0;          // 1 + 10^13 * 0.01
    variables["quasi_f"] = 1.01;                    // 1 + 0.01

    // Derived
    variables["gamma_s"] = computeGamma_s();
}

// Update variable
void ReciprocationDecayModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
        if (name == "gamma_day") {
            variables["gamma_s"] = computeGamma_s();
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void ReciprocationDecayModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        if (name == "gamma_day") {
            variables["gamma_s"] = computeGamma_s();
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void ReciprocationDecayModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute ? (day?�)
double ReciprocationDecayModule::computeGamma_day() {
    return variables["gamma_day"];
}

// Compute ? in s?�
double ReciprocationDecayModule::computeGamma_s() {
    return computeGamma_day() / variables["day_to_s"];
}

// Compute cos(? t_n)
double ReciprocationDecayModule::computeCosPiTn(double t_n) {
    variables["t_n"] = t_n;
    return std::cos(variables["pi"] * t_n);
}

// Compute exp(-? t cos(? t_n)) (t in days)
double ReciprocationDecayModule::computeExpTerm(double t_day, double t_n) {
    variables["t_day"] = t_day;
    double cos_pi_tn = computeCosPiTn(t_n);
    double arg = - computeGamma_day() * t_day * cos_pi_tn;
    return std::exp(arg);
}

// Compute 1 - exp(-? t cos(? t_n))
double ReciprocationDecayModule::computeOneMinusExp(double t_day, double t_n) {
    return 1.0 - computeExpTerm(t_day, t_n);
}

// Simplified U_m example (J/m�)
double ReciprocationDecayModule::computeUmExample(double t_day, double t_n, double mu_over_rj) {
    double one_minus_exp = computeOneMinusExp(t_day, t_n);
    double phi_hat = 1.0;
    double p_scm = variables["P_SCm"];
    double e_react = variables["E_react"];
    double heaviside_f = variables["heaviside_f"];
    double quasi_f = variables["quasi_f"];
    return (mu_over_rj * one_minus_exp * phi_hat) * p_scm * e_react * heaviside_f * quasi_f;
}

// Equation text
std::string ReciprocationDecayModule::getEquationText() {
    return "? = 0.00005 day?� (~5.8e-10 s?�; timescale ~55 years);\n"
           "In U_m: ... (1 - exp(-? t cos(? t_n))) ... (t days, reciprocating decay/growth).\n"
           "Negative cos(? t_n): exp(+? t) >1 (growth, negentropic TRZ).\n"
           "Example t=1000 days, t_n=0: 1-exp ?0.049, U_m ?1.12e66 J/m�.\n"
           "UQFF: Slow decay for magnetic strings; cyclic via cos(? t_n) in jets/nebulae/mergers.";
}

// Print variables
void ReciprocationDecayModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Print effects
void ReciprocationDecayModule::printDecayEffects(double t_day, double t_n) {
    double cos_pi = computeCosPiTn(t_n);
    double exp_val = computeExpTerm(t_day, t_n);
    double one_minus = computeOneMinusExp(t_day, t_n);
    double um_ex = computeUmExample(t_day, t_n);
    std::cout << "Decay Effects at t=" << t_day << " days, t_n=" << t_n << ":\n";
    std::cout << "cos(? t_n) = " << cos_pi << "\n";
    std::cout << "exp(-? t cos(? t_n)) = " << exp_val << "\n";
    std::cout << "1 - exp(...) = " << one_minus << "\n";
    std::cout << "U_m example contrib = " << um_ex << " J/m�\n";
}

// ===== ENHANCED METHODS IMPLEMENTATION =====

namespace reciprocation_decay_saved_states {
    static std::map<std::string, std::map<std::string, double>> saved_states;
}

// Variable Management
void ReciprocationDecayModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void ReciprocationDecayModule::removeVariable(const std::string& name) {
    variables.erase(name);
}

void ReciprocationDecayModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
    }
}

std::vector<std::string> ReciprocationDecayModule::listVariables() const {
    std::vector<std::string> names;
    for (const auto& pair : variables) {
        names.push_back(pair.first);
    }
    return names;
}

std::string ReciprocationDecayModule::getSystemName() const {
    return "Reciprocation_Decay_UQFF";
}

// Batch Operations
void ReciprocationDecayModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = func(variables[name]);
        }
    }
    // Update derived
    variables["gamma_s"] = computeGamma_s();
}

void ReciprocationDecayModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// Self-Expansion: Domain-Specific Scales
void ReciprocationDecayModule::expandParameterSpace(double decay_scale, double time_scale, double energy_scale) {
    variables["gamma_day"] *= decay_scale;
    variables["gamma_s"] = computeGamma_s();
    variables["t_day"] *= time_scale;
    variables["t_n"] *= time_scale;
    variables["E_react"] *= energy_scale;
}

void ReciprocationDecayModule::expandDecayScale(double gamma_factor, double timescale_factor) {
    variables["gamma_day"] *= gamma_factor;
    variables["gamma_s"] = computeGamma_s();
    // Effective timescale = 1/γ, so inverse scaling
    double effective_timescale = 1.0 / variables["gamma_day"];
    effective_timescale *= timescale_factor;
    variables["gamma_day"] = 1.0 / effective_timescale;
    variables["gamma_s"] = computeGamma_s();
}

void ReciprocationDecayModule::expandOscillationScale(double pi_factor, double cos_factor) {
    // Note: π is a mathematical constant, but we can scale t_n to affect cos(π t_n)
    variables["t_n"] *= cos_factor;
    // If we want to scale π itself (non-physical but for testing):
    // variables["pi"] *= pi_factor;
}

void ReciprocationDecayModule::expandMagneticScale(double mu_factor, double heaviside_factor) {
    variables["mu_over_rj"] *= mu_factor;
    variables["heaviside_f"] *= heaviside_factor;
}

// Self-Refinement
void ReciprocationDecayModule::autoRefineParameters(const std::string& target, double goal) {
    if (target == "gamma_day") {
        variables["gamma_day"] = goal;
        variables["gamma_s"] = computeGamma_s();
    } else if (target == "gamma_s") {
        variables["gamma_s"] = goal;
        variables["gamma_day"] = goal * variables["day_to_s"];
    } else if (target == "timescale_years") {
        // Timescale ≈ 1/γ (in days), convert to years: timescale_days / 365.25
        double timescale_days = goal * 365.25;
        variables["gamma_day"] = 1.0 / timescale_days;
        variables["gamma_s"] = computeGamma_s();
    } else if (target == "decay_fraction") {
        // Target specific decay fraction at t_day=1000, t_n=0
        // 1 - exp(-γ t) = goal
        // exp(-γ t) = 1 - goal
        // -γ t = ln(1 - goal)
        // γ = -ln(1 - goal) / t
        double t = 1000.0;  // reference time
        if (goal < 1.0 && goal > 0.0) {
            variables["gamma_day"] = -std::log(1.0 - goal) / t;
            variables["gamma_s"] = computeGamma_s();
        }
    }
}

void ReciprocationDecayModule::calibrateToObservations(const std::map<std::string, double>& observations) {
    for (const auto& obs : observations) {
        if (variables.find(obs.first) != variables.end()) {
            variables[obs.first] = obs.second;
            if (obs.first == "gamma_day") {
                variables["gamma_s"] = computeGamma_s();
            }
        }
    }
}

void ReciprocationDecayModule::optimizeForMetric(const std::string& metric) {
    if (metric == "fast_decay") {
        // Increase decay rate (shorter timescale)
        variables["gamma_day"] *= 2.0;
        variables["gamma_s"] = computeGamma_s();
    } else if (metric == "slow_decay") {
        // Decrease decay rate (longer timescale)
        variables["gamma_day"] *= 0.5;
        variables["gamma_s"] = computeGamma_s();
    } else if (metric == "standard_55years") {
        // Reset to standard ~55 year timescale (γ = 0.00005 day⁻¹)
        variables["gamma_day"] = 0.00005;
        variables["gamma_s"] = computeGamma_s();
    } else if (metric == "century_scale") {
        // 100 year timescale: γ = 1/(100*365.25) ≈ 2.74e-5 day⁻¹
        variables["gamma_day"] = 1.0 / (100.0 * 365.25);
        variables["gamma_s"] = computeGamma_s();
    } else if (metric == "decade_scale") {
        // 10 year timescale: γ = 1/(10*365.25) ≈ 2.74e-4 day⁻¹
        variables["gamma_day"] = 1.0 / (10.0 * 365.25);
        variables["gamma_s"] = computeGamma_s();
    }
}

// Parameter Exploration
std::vector<std::map<std::string, double>> ReciprocationDecayModule::generateVariations(int count, double variation_pct) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - variation_pct, 1.0 + variation_pct);
    
    for (int i = 0; i < count; ++i) {
        std::map<std::string, double> variant = variables;
        for (auto& pair : variant) {
            if (pair.first != "day_to_s" && pair.first != "pi" && pair.first != "gamma_s") {
                pair.second *= dis(gen);
            }
        }
        // Recalculate derived
        variant["gamma_s"] = variant["gamma_day"] / variant["day_to_s"];
        variations.push_back(variant);
    }
    return variations;
}

// Adaptive Evolution
void ReciprocationDecayModule::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - mutation_rate, 1.0 + mutation_rate);
    
    for (auto& pair : variables) {
        if (pair.first != "day_to_s" && pair.first != "pi" && pair.first != "gamma_s") {
            pair.second *= dis(gen);
        }
    }
    // Recalculate derived
    variables["gamma_s"] = computeGamma_s();
}

void ReciprocationDecayModule::evolveSystem(int generations, std::function<double()> fitness_func) {
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
void ReciprocationDecayModule::saveState(const std::string& label) {
    reciprocation_decay_saved_states::saved_states[label] = variables;
}

void ReciprocationDecayModule::restoreState(const std::string& label) {
    if (reciprocation_decay_saved_states::saved_states.find(label) != reciprocation_decay_saved_states::saved_states.end()) {
        variables = reciprocation_decay_saved_states::saved_states[label];
    }
}

std::vector<std::string> ReciprocationDecayModule::listSavedStates() const {
    std::vector<std::string> labels;
    for (const auto& pair : reciprocation_decay_saved_states::saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string ReciprocationDecayModule::exportState() const {
    std::ostringstream oss;
    oss << "ReciprocationDecay_State_Export:\n";
    for (const auto& pair : variables) {
        oss << pair.first << "=" << std::scientific << pair.second << "\n";
    }
    return oss.str();
}

// System Analysis
std::map<std::string, double> ReciprocationDecayModule::sensitivityAnalysis(const std::vector<std::string>& params) {
    std::map<std::string, double> sensitivities;
    double t_test = 1000.0;  // Test at 1000 days
    double t_n_test = 0.0;
    double baseline = computeOneMinusExp(t_test, t_n_test);
    
    for (const auto& param : params) {
        if (variables.find(param) != variables.end() && param != "day_to_s" && param != "pi") {
            double original = variables[param];
            variables[param] = original * 1.01;
            
            if (param == "gamma_day") {
                variables["gamma_s"] = computeGamma_s();
            }
            
            double perturbed = computeOneMinusExp(t_test, t_n_test);
            sensitivities[param] = (perturbed - baseline) / baseline;
            
            // Restore
            variables[param] = original;
            if (param == "gamma_day") {
                variables["gamma_s"] = computeGamma_s();
            }
        }
    }
    return sensitivities;
}

std::string ReciprocationDecayModule::generateReport() const {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(3);
    oss << "===== Reciprocation Decay Module Report =====\n";
    oss << "System: " << getSystemName() << "\n\n";
    
    oss << "Decay Rate Parameters:\n";
    oss << "  γ (day⁻¹) = " << variables.at("gamma_day") << " day⁻¹\n";
    oss << "  γ (s⁻¹) = " << variables.at("gamma_s") << " s⁻¹\n";
    double timescale_days = 1.0 / variables.at("gamma_day");
    double timescale_years = timescale_days / 365.25;
    oss << "  Characteristic timescale = " << std::fixed << std::setprecision(1) 
        << timescale_days << " days ≈ " << timescale_years << " years\n\n";
    
    oss << std::scientific;
    oss << "Time Parameters:\n";
    oss << "  t (current time) = " << variables.at("t_day") << " days\n";
    oss << "  t_n (negative time) = " << variables.at("t_n") << " days\n";
    oss << "  π = " << std::fixed << std::setprecision(15) << variables.at("pi") << "\n\n";
    
    oss << std::scientific;
    oss << "Magnetic and Energy Parameters:\n";
    oss << "  μ_j/r_j = " << variables.at("mu_over_rj") << " T·m²\n";
    oss << "  P_SCm = " << variables.at("P_SCm") << "\n";
    oss << "  E_react = " << variables.at("E_react") << " J\n";
    oss << "  Heaviside factor = " << variables.at("heaviside_f") << "\n";
    oss << "  Quasi factor = " << variables.at("quasi_f") << "\n\n";
    
    oss << "Decay Function Examples:\n";
    std::vector<double> test_times = {100.0, 500.0, 1000.0, 5000.0, 10000.0};
    for (double t : test_times) {
        double cos_val = std::cos(variables.at("pi") * 0.0);  // t_n = 0
        double exp_val = std::exp(-variables.at("gamma_day") * t * cos_val);
        double one_minus = 1.0 - exp_val;
        oss << "  t=" << std::fixed << std::setprecision(0) << t << " days: ";
        oss << "exp(-γt) = " << std::scientific << exp_val << ", ";
        oss << "1-exp = " << one_minus << "\n";
    }
    oss << "\n";
    
    oss << "Oscillation via cos(π t_n):\n";
    std::vector<double> test_tn = {-1.0, -0.5, 0.0, 0.5, 1.0};
    for (double tn : test_tn) {
        double cos_val = std::cos(variables.at("pi") * tn);
        oss << "  t_n=" << std::fixed << std::setprecision(1) << tn << ": ";
        oss << "cos(π t_n) = " << std::fixed << std::setprecision(3) << cos_val;
        if (cos_val < 0) {
            oss << " (negative → exp(+γt) growth, negentropic)\n";
        } else {
            oss << " (positive → exp(-γt) decay)\n";
        }
    }
    oss << "\n" << std::scientific;
    
    oss << "U_m Contribution Example (t=1000 days, t_n=0):\n";
    double t_test = 1000.0;
    double tn_test = 0.0;
    double cos_pi = std::cos(variables.at("pi") * tn_test);
    double exp_term = std::exp(-variables.at("gamma_day") * t_test * cos_pi);
    double one_minus = 1.0 - exp_term;
    double um_contrib = variables.at("mu_over_rj") * one_minus * 1.0 * 
                        variables.at("P_SCm") * variables.at("E_react") * 
                        variables.at("heaviside_f") * variables.at("quasi_f");
    oss << "  cos(π t_n) = " << cos_pi << "\n";
    oss << "  exp(-γ t cos(π t_n)) = " << exp_term << "\n";
    oss << "  1 - exp(...) = " << one_minus << "\n";
    oss << "  U_m ≈ " << um_contrib << " J/m³\n\n";
    
    oss << "Physical Interpretation:\n";
    if (timescale_years < 10.0) {
        oss << "  Fast decay timescale (<10 years)\n";
    } else if (timescale_years < 100.0) {
        oss << "  Moderate decay timescale (10-100 years)\n";
    } else {
        oss << "  Slow decay timescale (>100 years)\n";
    }
    
    oss << "  Applications:\n";
    oss << "    - Magnetic string decay: Long-term energy dissipation\n";
    oss << "    - Reciprocation: cos(π t_n) creates cyclic decay/growth\n";
    oss << "    - Negentropic regimes: Negative cos → growth (TRZ dynamics)\n";
    oss << "    - Jet/nebula evolution: 55-year timescale aligns with observed structures\n";
    oss << "    - Merger dynamics: Time-reversal scenarios with negative time\n";
    
    return oss.str();
}

bool ReciprocationDecayModule::validateConsistency() const {
    bool valid = true;
    
    // Check γ is positive
    if (variables.find("gamma_day") != variables.end() && variables.at("gamma_day") <= 0) {
        std::cerr << "Error: gamma_day <= 0 (decay rate must be positive)\n";
        valid = false;
    }
    
    // Check γ_s consistency
    if (variables.find("gamma_s") != variables.end()) {
        double expected = variables.at("gamma_day") / variables.at("day_to_s");
        double actual = variables.at("gamma_s");
        if (std::abs(expected - actual) / expected > 1e-9) {
            std::cerr << "Error: gamma_s inconsistent (expected " << expected << ", got " << actual << ")\n";
            valid = false;
        }
    }
    
    // Check timescale is reasonable (1-1000 years)
    double timescale_years = (1.0 / variables.at("gamma_day")) / 365.25;
    if (timescale_years < 1.0 || timescale_years > 1000.0) {
        std::cerr << "Warning: Timescale outside typical range [1, 1000] years (current: " 
                  << timescale_years << " years)\n";
    }
    
    return valid;
}

void ReciprocationDecayModule::autoCorrectAnomalies() {
    // Reset γ to typical value if out of range
    double timescale_years = (1.0 / variables["gamma_day"]) / 365.25;
    if (variables["gamma_day"] <= 0 || timescale_years > 1000.0) {
        variables["gamma_day"] = 0.00005;  // Standard 55 years
        variables["gamma_s"] = computeGamma_s();
    }
    
    // Ensure day_to_s is correct
    if (variables["day_to_s"] != 86400.0) {
        variables["day_to_s"] = 86400.0;
        variables["gamma_s"] = computeGamma_s();
    }
    
    // Ensure π is correct
    if (std::abs(variables["pi"] - 3.141592653589793) > 1e-12) {
        variables["pi"] = 3.141592653589793;
    }
}

// Example usage in base program (snippet)
int main() {
    ReciprocationDecayModule module;
    std::cout << "===== Reciprocation Decay Module Enhanced Demonstration =====\n\n";
    
    // Step 1: Report initial state
    std::cout << "STEP 1: Initial Configuration (γ = 0.00005 day⁻¹, ~55 years)\n";
    std::cout << module.generateReport() << "\n";
    
    // Step 2: Compute decay at various times
    std::cout << "STEP 2: Decay Function at Different Times (t_n = 0)\n";
    std::vector<double> test_times = {100.0, 500.0, 1000.0, 5000.0, 10000.0, 20000.0};
    for (double t : test_times) {
        double one_minus_exp = module.computeOneMinusExp(t, 0.0);
        std::cout << "  t=" << std::fixed << std::setprecision(0) << t << " days: ";
        std::cout << "1-exp(-γt) = " << std::scientific << one_minus_exp << "\n";
    }
    std::cout << "\n";
    
    // Step 3: Save initial state
    std::cout << "STEP 3: Save Initial State\n";
    module.saveState("standard_55years");
    std::cout << "State saved as 'standard_55years'\n\n";
    
    // Step 4: Test oscillation via cos(π t_n)
    std::cout << "STEP 4: Test Oscillation via cos(π t_n) at t=1000 days\n";
    std::vector<double> test_tn = {-1.0, -0.5, 0.0, 0.5, 1.0};
    for (double tn : test_tn) {
        double cos_val = module.computeCosPiTn(tn);
        double exp_val = module.computeExpTerm(1000.0, tn);
        double one_minus = module.computeOneMinusExp(1000.0, tn);
        std::cout << "  t_n=" << std::fixed << std::setprecision(1) << tn << ": ";
        std::cout << "cos(π t_n)=" << std::fixed << std::setprecision(3) << cos_val << ", ";
        std::cout << "1-exp=" << std::scientific << one_minus;
        if (cos_val < 0) {
            std::cout << " (GROWTH, negentropic)\n";
        } else {
            std::cout << " (decay)\n";
        }
    }
    std::cout << "\n";
    
    // Step 5: Test fast decay (10 year timescale)
    std::cout << "STEP 5: Test Fast Decay (10 year timescale)\n";
    module.optimizeForMetric("decade_scale");
    double gamma_fast = module.getVariable("gamma_day");
    double timescale_fast = (1.0 / gamma_fast) / 365.25;
    double decay_fast = module.computeOneMinusExp(1000.0, 0.0);
    std::cout << "Fast γ = " << std::scientific << gamma_fast << " day⁻¹\n";
    std::cout << "Timescale = " << std::fixed << std::setprecision(1) << timescale_fast << " years\n";
    std::cout << "1-exp at t=1000 days = " << std::scientific << decay_fast << "\n\n";
    
    // Step 6: Test slow decay (100 year timescale)
    std::cout << "STEP 6: Test Slow Decay (100 year timescale)\n";
    module.restoreState("standard_55years");
    module.optimizeForMetric("century_scale");
    double gamma_slow = module.getVariable("gamma_day");
    double timescale_slow = (1.0 / gamma_slow) / 365.25;
    double decay_slow = module.computeOneMinusExp(1000.0, 0.0);
    std::cout << "Slow γ = " << std::scientific << gamma_slow << " day⁻¹\n";
    std::cout << "Timescale = " << std::fixed << std::setprecision(1) << timescale_slow << " years\n";
    std::cout << "1-exp at t=1000 days = " << std::scientific << decay_slow << "\n\n";
    
    // Step 7: Restore and expand decay scale
    std::cout << "STEP 7: Expand Decay Scale (γ x2, faster decay)\n";
    module.restoreState("standard_55years");
    module.expandDecayScale(2.0, 1.0);
    double gamma_expanded = module.getVariable("gamma_day");
    double decay_expanded = module.computeOneMinusExp(1000.0, 0.0);
    std::cout << "Expanded γ = " << std::scientific << gamma_expanded << " day⁻¹\n";
    std::cout << "New timescale = " << std::fixed << std::setprecision(1) 
              << ((1.0/gamma_expanded)/365.25) << " years\n";
    std::cout << "1-exp at t=1000 days = " << std::scientific << decay_expanded << "\n\n";
    
    // Step 8: Restore and expand magnetic scale
    std::cout << "STEP 8: Expand Magnetic Scale (μ/r x2, Heaviside x1.5)\n";
    module.restoreState("standard_55years");
    module.expandMagneticScale(2.0, 1.5);
    double mu_new = module.getVariable("mu_over_rj");
    double hf_new = module.getVariable("heaviside_f");
    double um_new = module.computeUmExample(1000.0, 0.0);
    std::cout << "New μ_j/r_j = " << std::scientific << mu_new << " T·m²\n";
    std::cout << "New Heaviside = " << hf_new << "\n";
    std::cout << "U_m at t=1000 days = " << um_new << " J/m³\n\n";
    
    // Step 9: Sensitivity analysis
    std::cout << "STEP 9: Sensitivity Analysis (at t=1000 days, t_n=0)\n";
    module.restoreState("standard_55years");
    std::vector<std::string> params = {"gamma_day", "t_day", "t_n", "mu_over_rj", "E_react"};
    auto sensitivities = module.sensitivityAnalysis(params);
    for (const auto& pair : sensitivities) {
        std::cout << "  ∂(1-exp)/∂" << pair.first << " ≈ " << std::scientific << pair.second << " (normalized)\n";
    }
    std::cout << "\n";
    
    // Step 10: Generate variations
    std::cout << "STEP 10: Generate Parameter Variations (5 variants, ±10%)\n";
    auto variations = module.generateVariations(5, 0.1);
    for (int i = 0; i < variations.size(); ++i) {
        double var_gamma = variations[i]["gamma_day"];
        double var_timescale = (1.0 / var_gamma) / 365.25;
        std::cout << "  Variant " << (i+1) << ": γ=" << std::scientific << var_gamma 
                  << " day⁻¹ (timescale=" << std::fixed << std::setprecision(1) 
                  << var_timescale << " years)\n";
    }
    std::cout << "\n";
    
    // Step 11: Auto-refine to target decay fraction
    std::cout << "STEP 11: Auto-Refine to Target Decay Fraction (0.1 at t=1000 days)\n";
    module.restoreState("standard_55years");
    module.autoRefineParameters("decay_fraction", 0.1);
    double refined_gamma = module.getVariable("gamma_day");
    double refined_decay = module.computeOneMinusExp(1000.0, 0.0);
    std::cout << "Refined γ = " << std::scientific << refined_gamma << " day⁻¹\n";
    std::cout << "Achieved decay fraction = " << refined_decay << " (target: 0.1)\n";
    std::cout << "Timescale = " << std::fixed << std::setprecision(1) 
              << ((1.0/refined_gamma)/365.25) << " years\n\n";
    
    // Step 12: Auto-refine to target timescale
    std::cout << "STEP 12: Auto-Refine to Target Timescale (30 years)\n";
    module.restoreState("standard_55years");
    module.autoRefineParameters("timescale_years", 30.0);
    double gamma_30y = module.getVariable("gamma_day");
    double timescale_30y = (1.0 / gamma_30y) / 365.25;
    std::cout << "Refined γ = " << std::scientific << gamma_30y << " day⁻¹\n";
    std::cout << "Achieved timescale = " << std::fixed << std::setprecision(1) 
              << timescale_30y << " years (target: 30.0)\n\n";
    
    // Step 13: Calibrate to observations
    std::cout << "STEP 13: Calibrate to Observational Data\n";
    module.restoreState("standard_55years");
    std::map<std::string, double> observations;
    observations["gamma_day"] = 0.0001;  // Observed faster decay
    observations["mu_over_rj"] = 3.0e10;  // Observed stronger field
    module.calibrateToObservations(observations);
    std::cout << "Calibrated γ = " << std::scientific << module.getVariable("gamma_day") << " day⁻¹\n";
    std::cout << "Calibrated μ/r = " << module.getVariable("mu_over_rj") << " T·m²\n";
    std::cout << "Timescale = " << std::fixed << std::setprecision(1) 
              << ((1.0/module.getVariable("gamma_day"))/365.25) << " years\n\n";
    
    // Step 14: Optimize for metric
    std::cout << "STEP 14: Optimize for 'standard_55years' Metric\n";
    module.optimizeForMetric("standard_55years");
    std::cout << "Optimized γ = " << std::scientific << module.getVariable("gamma_day") << " day⁻¹\n";
    std::cout << "Timescale = " << std::fixed << std::setprecision(1) 
              << ((1.0/module.getVariable("gamma_day"))/365.25) << " years\n\n";
    
    // Step 15: Mutate parameters
    std::cout << "STEP 15: Mutate Parameters (5% random variation)\n";
    module.mutateParameters(0.05);
    std::cout << "Mutated γ = " << std::scientific << module.getVariable("gamma_day") << " day⁻¹\n";
    std::cout << "Mutated timescale = " << std::fixed << std::setprecision(1) 
              << ((1.0/module.getVariable("gamma_day"))/365.25) << " years\n\n";
    
    // Step 16: Validate consistency
    std::cout << "STEP 16: Validate Consistency\n";
    bool valid = module.validateConsistency();
    std::cout << "Consistency check: " << (valid ? "PASSED" : "FAILED") << "\n\n";
    
    // Step 17: Introduce anomaly and auto-correct
    std::cout << "STEP 17: Introduce Anomaly and Auto-Correct\n";
    module.createVariable("gamma_anomaly", -0.001);  // Invalid negative γ
    module.removeVariable("gamma_day");
    module.createVariable("gamma_day", -0.001);
    std::cout << "Introduced invalid γ = " << module.getVariable("gamma_day") << " day⁻¹\n";
    module.autoCorrectAnomalies();
    std::cout << "Auto-corrected γ = " << std::scientific << module.getVariable("gamma_day") << " day⁻¹\n";
    std::cout << "Timescale = " << std::fixed << std::setprecision(1) 
              << ((1.0/module.getVariable("gamma_day"))/365.25) << " years\n";
    bool valid_after = module.validateConsistency();
    std::cout << "Consistency after correction: " << (valid_after ? "PASSED" : "FAILED") << "\n\n";
    
    // Step 18: List saved states and export
    std::cout << "STEP 18: List Saved States and Export Final State\n";
    auto states = module.listSavedStates();
    std::cout << "Saved states:\n";
    for (const auto& state : states) {
        std::cout << "  - " << state << "\n";
    }
    std::cout << "\n";
    std::string exported = module.exportState();
    std::cout << exported << "\n";
    
    std::cout << "===== Demonstration Complete =====\n";
    return 0;
}

// Original commented example preserved below for reference:
// #include "ReciprocationDecayModule.h"
// int main() {
//     ReciprocationDecayModule mod;
//     double gamma = mod.computeGamma_day();
//     std::cout << "γ = " << gamma << " day⁻¹\n";
//     mod.printDecayEffects(1000.0, 0.0);
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("gamma_day", 0.0001);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o decay_test decay_test.cpp ReciprocationDecayModule.cpp -lm
// Sample: γ=5e-5 day⁻¹; t=1000 days: 1-exp≈0.049; U_m≈1.12e66 J/m³.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

ReciprocationDecayModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeGamma_day, computeGamma_s, computeCosPiTn, computeExpTerm, computeOneMinusExp, computeUmExample) are clear, concise, and variable - driven.
- Automatic recalculation of derived variables(gamma_s) when dependencies change.
- Output and debugging functions(printVariables, printDecayEffects, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in reciprocation decay modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.