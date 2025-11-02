// ScmReactivityDecayModule.h
// Modular C++ implementation of the [SCm] Reactivity Decay Rate (?) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes ?=0.0005 day?� (~5.8e-6 s?�); used in E_react = 10^46 * exp(-? t) for decay in U_m, U_bi, etc.
// Pluggable: #include "ScmReactivityDecayModule.h"
// ScmReactivityDecayModule mod; mod.computeE_react(0.0); mod.updateVariable("kappa_day", new_value);
// Variables in std::map; example for Sun at t=0 (E_react=1e46); t=2000 days: ~3.68e45.
// Approximations: t in days; timescale ~5.5 years; integrates into U_m example.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef SCM_REACTIVITY_DECAY_MODULE_H
#define SCM_REACTIVITY_DECAY_MODULE_H

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

class ScmReactivityDecayModule {
private:
    std::map<std::string, double> variables;
    double computeKappa_s();  // ? in s?�
    double computeE_react(double t_day);
    double computeUmExample(double t_day);

public:
    // Constructor: Initialize with framework defaults
    ScmReactivityDecayModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeKappa_day();  // 0.0005 day?�
    double computeKappa_s();    // ~5.8e-6 s?�
    double computeE_react(double t_day);  // 1e46 * exp(-? t)
    double computeUmExample(double t_day);  // Simplified U_m with E_react (J/m�)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print decay effects
    void printDecayEffects(double t_day = 2000.0);

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
    void expandParameterSpace(double decay_scale, double energy_scale, double time_scale);
    void expandDecayScale(double kappa_factor, double timescale_factor);  // κ and timescales
    void expandEnergyScale(double ereact_factor, double base_factor);     // E_react and base energy
    void expandMagneticScale(double mu_factor, double amplification_factor); // μ and amplifications

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

#endif // SCM_REACTIVITY_DECAY_MODULE_H

// ScmReactivityDecayModule.cpp
#include "ScmReactivityDecayModule.h"

// Constructor: Set framework defaults
ScmReactivityDecayModule::ScmReactivityDecayModule() {
    // Universal constants
    variables["kappa_day"] = 0.0005;                // day?�
    variables["day_to_s"] = 86400.0;                // s/day
    variables["E_react_base"] = 1e46;               // J
    variables["t_day"] = 0.0;                       // days
    variables["mu_over_rj"] = 2.26e10;              // T m� (example)
    variables["P_SCm"] = 1.0;                       // Normalized
    variables["heaviside_f"] = 1e11 + 1.0;          // 1 + 10^13 * 0.01
    variables["quasi_f"] = 1.01;                    // 1 + 0.01
    variables["one_minus_exp"] = 1.0;               // At t=0

    // Derived
    variables["kappa_s"] = computeKappa_s();
}

// Update variable
void ScmReactivityDecayModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
        if (name == "kappa_day") {
            variables["kappa_s"] = computeKappa_s();
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void ScmReactivityDecayModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        if (name == "kappa_day") {
            variables["kappa_s"] = computeKappa_s();
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void ScmReactivityDecayModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute ? (day?�)
double ScmReactivityDecayModule::computeKappa_day() {
    return variables["kappa_day"];
}

// Compute ? in s?�
double ScmReactivityDecayModule::computeKappa_s() {
    return computeKappa_day() / variables["day_to_s"];
}

// Compute E_react = 1e46 * exp(-? t)
double ScmReactivityDecayModule::computeE_react(double t_day) {
    variables["t_day"] = t_day;
    double arg = - computeKappa_day() * t_day;
    return variables["E_react_base"] * std::exp(arg);
}

// Simplified U_m example with E_react
double ScmReactivityDecayModule::computeUmExample(double t_day) {
    double e_react = computeE_react(t_day);
    double one_minus_exp = variables["one_minus_exp"];  // Placeholder; full would compute
    double phi_hat = 1.0;
    double p_scm = variables["P_SCm"];
    double heaviside_f = variables["heaviside_f"];
    double quasi_f = variables["quasi_f"];
    return (variables["mu_over_rj"] * one_minus_exp * phi_hat) * p_scm * e_react * heaviside_f * quasi_f;
}

// Equation text
std::string ScmReactivityDecayModule::getEquationText() {
    return "E_react = 10^46 * exp(-? t) (t days); ?=0.0005 day?� (~5.8e-6 s?�, timescale ~5.5 years).\n"
           "In U_m, U_bi, U_i, U_gi: ... * E_react * ... (decays [SCm] reactivity).\n"
           "Example t=0: E_react=1e46 J; t=2000 days: ~3.68e45 J (~36.8%).\n"
           "U_m (t=0): ?2.28e65 J/m�; t=2000: ?8.39e64 J/m�.\n"
           "Role: Gradual [SCm]-[UA] interaction loss; temporal evolution in jets/nebulae/mergers.\n"
           "UQFF: Models reactivity decay; energy dissipation over cosmic time.";
}

// Print variables
void ScmReactivityDecayModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Print effects
void ScmReactivityDecayModule::printDecayEffects(double t_day) {
    double e_react = computeE_react(t_day);
    double um_ex = computeUmExample(t_day);
    double fraction = e_react / variables["E_react_base"];
    std::cout << "[SCm] Decay Effects at t=" << t_day << " days:\n";
    std::cout << "E_react = " << std::scientific << e_react << " J (" << fraction << " of initial)\n";
    std::cout << "U_m example = " << um_ex << " J/m�\n";
}

// ===== ENHANCED METHODS IMPLEMENTATION =====

namespace scm_reactivity_decay_saved_states {
    static std::map<std::string, std::map<std::string, double>> saved_states;
}

// Variable Management
void ScmReactivityDecayModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void ScmReactivityDecayModule::removeVariable(const std::string& name) {
    variables.erase(name);
}

void ScmReactivityDecayModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
    }
}

std::vector<std::string> ScmReactivityDecayModule::listVariables() const {
    std::vector<std::string> names;
    for (const auto& pair : variables) {
        names.push_back(pair.first);
    }
    return names;
}

std::string ScmReactivityDecayModule::getSystemName() const {
    return "SCm_Reactivity_Decay_UQFF";
}

// Batch Operations
void ScmReactivityDecayModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = func(variables[name]);
        }
    }
    // Update derived
    variables["kappa_s"] = computeKappa_s();
}

void ScmReactivityDecayModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// Self-Expansion: Domain-Specific Scales
void ScmReactivityDecayModule::expandParameterSpace(double decay_scale, double energy_scale, double time_scale) {
    variables["kappa_day"] *= decay_scale;
    variables["kappa_s"] = computeKappa_s();
    variables["E_react_base"] *= energy_scale;
    variables["t_day"] *= time_scale;
}

void ScmReactivityDecayModule::expandDecayScale(double kappa_factor, double timescale_factor) {
    variables["kappa_day"] *= kappa_factor;
    variables["kappa_s"] = computeKappa_s();
    // Effective timescale = 1/κ, so inverse scaling
    double effective_timescale = 1.0 / variables["kappa_day"];
    effective_timescale *= timescale_factor;
    variables["kappa_day"] = 1.0 / effective_timescale;
    variables["kappa_s"] = computeKappa_s();
}

void ScmReactivityDecayModule::expandEnergyScale(double ereact_factor, double base_factor) {
    variables["E_react_base"] *= base_factor;
    // Note: Current E_react at current time also scales
}

void ScmReactivityDecayModule::expandMagneticScale(double mu_factor, double amplification_factor) {
    variables["mu_over_rj"] *= mu_factor;
    variables["heaviside_f"] *= amplification_factor;
    variables["quasi_f"] *= amplification_factor;
}

// Self-Refinement
void ScmReactivityDecayModule::autoRefineParameters(const std::string& target, double goal) {
    if (target == "kappa_day") {
        variables["kappa_day"] = goal;
        variables["kappa_s"] = computeKappa_s();
    } else if (target == "kappa_s") {
        variables["kappa_s"] = goal;
        variables["kappa_day"] = goal * variables["day_to_s"];
    } else if (target == "timescale_years") {
        // Timescale ≈ 1/κ (in days), convert to years: timescale_days / 365.25
        double timescale_days = goal * 365.25;
        variables["kappa_day"] = 1.0 / timescale_days;
        variables["kappa_s"] = computeKappa_s();
    } else if (target == "E_react_at_time") {
        // Target specific E_react at t_day=2000
        // E_react = E_base * exp(-κ t) = goal
        // E_base = goal / exp(-κ t)
        double t = 2000.0;
        double exp_factor = std::exp(-variables["kappa_day"] * t);
        if (std::abs(exp_factor) > 1e-9) {
            variables["E_react_base"] = goal / exp_factor;
        }
    } else if (target == "decay_fraction") {
        // Target specific decay fraction at t_day=2000
        // E_react / E_base = exp(-κ t) = goal
        // -κ t = ln(goal)
        // κ = -ln(goal) / t
        double t = 2000.0;
        if (goal > 0.0 && goal < 1.0) {
            variables["kappa_day"] = -std::log(goal) / t;
            variables["kappa_s"] = computeKappa_s();
        }
    }
}

void ScmReactivityDecayModule::calibrateToObservations(const std::map<std::string, double>& observations) {
    for (const auto& obs : observations) {
        if (variables.find(obs.first) != variables.end()) {
            variables[obs.first] = obs.second;
            if (obs.first == "kappa_day") {
                variables["kappa_s"] = computeKappa_s();
            }
        }
    }
}

void ScmReactivityDecayModule::optimizeForMetric(const std::string& metric) {
    if (metric == "fast_decay") {
        // Increase decay rate (shorter timescale)
        variables["kappa_day"] *= 2.0;
        variables["kappa_s"] = computeKappa_s();
    } else if (metric == "slow_decay") {
        // Decrease decay rate (longer timescale)
        variables["kappa_day"] *= 0.5;
        variables["kappa_s"] = computeKappa_s();
    } else if (metric == "standard_5_5years") {
        // Reset to standard ~5.5 year timescale (κ = 0.0005 day⁻¹)
        variables["kappa_day"] = 0.0005;
        variables["kappa_s"] = computeKappa_s();
    } else if (metric == "decade_scale") {
        // 10 year timescale: κ = 1/(10*365.25) ≈ 2.74e-4 day⁻¹
        variables["kappa_day"] = 1.0 / (10.0 * 365.25);
        variables["kappa_s"] = computeKappa_s();
    } else if (metric == "year_scale") {
        // 1 year timescale: κ = 1/365.25 ≈ 2.74e-3 day⁻¹
        variables["kappa_day"] = 1.0 / 365.25;
        variables["kappa_s"] = computeKappa_s();
    } else if (metric == "century_scale") {
        // 100 year timescale: κ = 1/(100*365.25) ≈ 2.74e-5 day⁻¹
        variables["kappa_day"] = 1.0 / (100.0 * 365.25);
        variables["kappa_s"] = computeKappa_s();
    }
}

// Parameter Exploration
std::vector<std::map<std::string, double>> ScmReactivityDecayModule::generateVariations(int count, double variation_pct) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - variation_pct, 1.0 + variation_pct);
    
    for (int i = 0; i < count; ++i) {
        std::map<std::string, double> variant = variables;
        for (auto& pair : variant) {
            if (pair.first != "day_to_s" && pair.first != "kappa_s") {
                pair.second *= dis(gen);
            }
        }
        // Recalculate derived
        variant["kappa_s"] = variant["kappa_day"] / variant["day_to_s"];
        variations.push_back(variant);
    }
    return variations;
}

// Adaptive Evolution
void ScmReactivityDecayModule::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - mutation_rate, 1.0 + mutation_rate);
    
    for (auto& pair : variables) {
        if (pair.first != "day_to_s" && pair.first != "kappa_s") {
            pair.second *= dis(gen);
        }
    }
    // Recalculate derived
    variables["kappa_s"] = computeKappa_s();
}

void ScmReactivityDecayModule::evolveSystem(int generations, std::function<double()> fitness_func) {
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
void ScmReactivityDecayModule::saveState(const std::string& label) {
    scm_reactivity_decay_saved_states::saved_states[label] = variables;
}

void ScmReactivityDecayModule::restoreState(const std::string& label) {
    if (scm_reactivity_decay_saved_states::saved_states.find(label) != scm_reactivity_decay_saved_states.end()) {
        variables = scm_reactivity_decay_saved_states::saved_states[label];
    }
}

std::vector<std::string> ScmReactivityDecayModule::listSavedStates() const {
    std::vector<std::string> labels;
    for (const auto& pair : scm_reactivity_decay_saved_states::saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string ScmReactivityDecayModule::exportState() const {
    std::ostringstream oss;
    oss << "ScmReactivityDecay_State_Export:\n";
    for (const auto& pair : variables) {
        oss << pair.first << "=" << std::scientific << pair.second << "\n";
    }
    return oss.str();
}

// System Analysis
std::map<std::string, double> ScmReactivityDecayModule::sensitivityAnalysis(const std::vector<std::string>& params) {
    std::map<std::string, double> sensitivities;
    double t_test = 2000.0;  // Test at 2000 days
    double baseline = computeE_react(t_test);
    
    for (const auto& param : params) {
        if (variables.find(param) != variables.end() && param != "day_to_s") {
            double original = variables[param];
            variables[param] = original * 1.01;
            
            if (param == "kappa_day") {
                variables["kappa_s"] = computeKappa_s();
            }
            
            double perturbed = computeE_react(t_test);
            sensitivities[param] = (perturbed - baseline) / baseline;
            
            // Restore
            variables[param] = original;
            if (param == "kappa_day") {
                variables["kappa_s"] = computeKappa_s();
            }
        }
    }
    return sensitivities;
}

std::string ScmReactivityDecayModule::generateReport() const {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(3);
    oss << "===== [SCm] Reactivity Decay Module Report =====\n";
    oss << "System: " << getSystemName() << "\n\n";
    
    oss << "Decay Rate Parameters:\n";
    oss << "  κ (day⁻¹) = " << variables.at("kappa_day") << " day⁻¹\n";
    oss << "  κ (s⁻¹) = " << variables.at("kappa_s") << " s⁻¹\n";
    double timescale_days = 1.0 / variables.at("kappa_day");
    double timescale_years = timescale_days / 365.25;
    oss << "  Characteristic timescale = " << std::fixed << std::setprecision(1) 
        << timescale_days << " days ≈ " << timescale_years << " years\n\n";
    
    oss << std::scientific;
    oss << "Energy Parameters:\n";
    oss << "  E_react_base = " << variables.at("E_react_base") << " J\n";
    oss << "  Current t = " << variables.at("t_day") << " days\n\n";
    
    oss << "E_react Decay Examples:\n";
    std::vector<double> test_times = {0.0, 500.0, 1000.0, 2000.0, 5000.0, 10000.0};
    for (double t : test_times) {
        double e_react = variables.at("E_react_base") * std::exp(-variables.at("kappa_day") * t);
        double fraction = e_react / variables.at("E_react_base");
        oss << "  t=" << std::fixed << std::setprecision(0) << t << " days: ";
        oss << "E_react = " << std::scientific << e_react << " J (";
        oss << std::fixed << std::setprecision(1) << (fraction * 100.0) << "% of initial)\n";
    }
    oss << "\n";
    
    oss << std::scientific;
    oss << "Magnetic and Amplification Parameters:\n";
    oss << "  μ_j/r_j = " << variables.at("mu_over_rj") << " T·m²\n";
    oss << "  P_SCm = " << variables.at("P_SCm") << "\n";
    oss << "  Heaviside factor = " << variables.at("heaviside_f") << "\n";
    oss << "  Quasi factor = " << variables.at("quasi_f") << "\n";
    oss << "  1 - exp (placeholder) = " << variables.at("one_minus_exp") << "\n\n";
    
    oss << "U_m Example Computation (with E_react decay):\n";
    double t_example = 2000.0;
    double e_react_2000 = variables.at("E_react_base") * std::exp(-variables.at("kappa_day") * t_example);
    double um_2000 = variables.at("mu_over_rj") * variables.at("one_minus_exp") * 1.0 * 
                     variables.at("P_SCm") * e_react_2000 * variables.at("heaviside_f") * variables.at("quasi_f");
    oss << "  t=" << std::fixed << std::setprecision(0) << t_example << " days:\n";
    oss << "    E_react = " << std::scientific << e_react_2000 << " J\n";
    oss << "    U_m ≈ " << um_2000 << " J/m³\n\n";
    
    // Compare to t=0
    double e_react_0 = variables.at("E_react_base");
    double um_0 = variables.at("mu_over_rj") * variables.at("one_minus_exp") * 1.0 * 
                  variables.at("P_SCm") * e_react_0 * variables.at("heaviside_f") * variables.at("quasi_f");
    oss << "  Comparison (t=2000 vs t=0):\n";
    oss << "    E_react ratio = " << std::fixed << std::setprecision(3) << (e_react_2000 / e_react_0) << "\n";
    oss << "    U_m ratio = " << (um_2000 / um_0) << "\n\n";
    
    oss << std::scientific;
    oss << "Physical Interpretation:\n";
    if (timescale_years < 1.0) {
        oss << "  Very fast reactivity decay (<1 year timescale)\n";
    } else if (timescale_years < 10.0) {
        oss << "  Fast reactivity decay (1-10 year timescale)\n";
    } else if (timescale_years < 100.0) {
        oss << "  Moderate reactivity decay (10-100 year timescale)\n";
    } else {
        oss << "  Slow reactivity decay (>100 year timescale)\n";
    }
    
    oss << "  Applications:\n";
    oss << "    - [SCm]-[UA] interaction decay: Gradual energy loss over cosmic time\n";
    oss << "    - Jet evolution: E_react decay affects magnetic string energy\n";
    oss << "    - Nebular aging: ~5.5 year timescale matches observed structures\n";
    oss << "    - Star formation: Reactivity decay modulates interior dynamics\n";
    oss << "    - Merger dynamics: Long-term energy dissipation\n";
    oss << "    - Universal terms: E_react appears in U_m, U_bi, U_i, U_gi\n";
    
    return oss.str();
}

bool ScmReactivityDecayModule::validateConsistency() const {
    bool valid = true;
    
    // Check κ is positive
    if (variables.find("kappa_day") != variables.end() && variables.at("kappa_day") <= 0) {
        std::cerr << "Error: kappa_day <= 0 (decay rate must be positive)\n";
        valid = false;
    }
    
    // Check κ_s consistency
    if (variables.find("kappa_s") != variables.end()) {
        double expected = variables.at("kappa_day") / variables.at("day_to_s");
        double actual = variables.at("kappa_s");
        if (std::abs(expected - actual) / expected > 1e-9) {
            std::cerr << "Error: kappa_s inconsistent (expected " << expected << ", got " << actual << ")\n";
            valid = false;
        }
    }
    
    // Check timescale is reasonable (0.1-1000 years)
    double timescale_years = (1.0 / variables.at("kappa_day")) / 365.25;
    if (timescale_years < 0.1 || timescale_years > 1000.0) {
        std::cerr << "Warning: Timescale outside typical range [0.1, 1000] years (current: " 
                  << timescale_years << " years)\n";
    }
    
    // Check E_react_base is positive
    if (variables.find("E_react_base") != variables.end() && variables.at("E_react_base") <= 0) {
        std::cerr << "Error: E_react_base <= 0 (base energy must be positive)\n";
        valid = false;
    }
    
    return valid;
}

void ScmReactivityDecayModule::autoCorrectAnomalies() {
    // Reset κ to typical value if out of range
    double timescale_years = (1.0 / variables["kappa_day"]) / 365.25;
    if (variables["kappa_day"] <= 0 || timescale_years > 1000.0 || timescale_years < 0.1) {
        variables["kappa_day"] = 0.0005;  // Standard 5.5 years
        variables["kappa_s"] = computeKappa_s();
    }
    
    // Ensure E_react_base is positive and reasonable
    if (variables["E_react_base"] <= 0 || variables["E_react_base"] > 1e50) {
        variables["E_react_base"] = 1e46;  // Standard value
    }
    
    // Ensure day_to_s is correct
    if (variables["day_to_s"] != 86400.0) {
        variables["day_to_s"] = 86400.0;
        variables["kappa_s"] = computeKappa_s();
    }
}

// Example usage in base program (snippet)
int main() {
    ScmReactivityDecayModule module;
    std::cout << "===== [SCm] Reactivity Decay Module Enhanced Demonstration =====\n\n";
    
    // Step 1: Report initial state
    std::cout << "STEP 1: Initial Configuration (κ = 0.0005 day⁻¹, ~5.5 years)\n";
    std::cout << module.generateReport() << "\n";
    
    // Step 2: Compute E_react at various times
    std::cout << "STEP 2: E_react Decay at Different Times\n";
    std::vector<double> test_times = {0.0, 500.0, 1000.0, 2000.0, 5000.0, 10000.0};
    for (double t : test_times) {
        double e_react = module.computeE_react(t);
        double fraction = e_react / module.getVariable("E_react_base");
        std::cout << "  t=" << std::fixed << std::setprecision(0) << t << " days: ";
        std::cout << "E_react = " << std::scientific << e_react << " J (";
        std::cout << std::fixed << std::setprecision(1) << (fraction * 100.0) << "%)\n";
    }
    std::cout << "\n";
    
    // Step 3: Save initial state
    std::cout << "STEP 3: Save Initial State\n";
    module.saveState("standard_5_5years");
    std::cout << "State saved as 'standard_5_5years'\n\n";
    
    // Step 4: Test fast decay (1 year timescale)
    std::cout << "STEP 4: Test Fast Decay (1 year timescale)\n";
    module.optimizeForMetric("year_scale");
    double kappa_fast = module.getVariable("kappa_day");
    double timescale_fast = (1.0 / kappa_fast) / 365.25;
    double ereact_fast = module.computeE_react(2000.0);
    std::cout << "Fast κ = " << std::scientific << kappa_fast << " day⁻¹\n";
    std::cout << "Timescale = " << std::fixed << std::setprecision(1) << timescale_fast << " years\n";
    std::cout << "E_react at t=2000 days = " << std::scientific << ereact_fast << " J\n";
    std::cout << "Fraction remaining: " << std::fixed << std::setprecision(3) 
              << (ereact_fast / module.getVariable("E_react_base")) << "\n\n";
    
    // Step 5: Test slow decay (100 year timescale)
    std::cout << "STEP 5: Test Slow Decay (100 year timescale)\n";
    module.restoreState("standard_5_5years");
    module.optimizeForMetric("century_scale");
    double kappa_slow = module.getVariable("kappa_day");
    double timescale_slow = (1.0 / kappa_slow) / 365.25;
    double ereact_slow = module.computeE_react(2000.0);
    std::cout << "Slow κ = " << std::scientific << kappa_slow << " day⁻¹\n";
    std::cout << "Timescale = " << std::fixed << std::setprecision(1) << timescale_slow << " years\n";
    std::cout << "E_react at t=2000 days = " << std::scientific << ereact_slow << " J\n";
    std::cout << "Fraction remaining: " << std::fixed << std::setprecision(3) 
              << (ereact_slow / module.getVariable("E_react_base")) << "\n\n";
    
    // Step 6: Test decade scale
    std::cout << "STEP 6: Test Decade Scale (10 year timescale)\n";
    module.restoreState("standard_5_5years");
    module.optimizeForMetric("decade_scale");
    double kappa_decade = module.getVariable("kappa_day");
    double timescale_decade = (1.0 / kappa_decade) / 365.25;
    double ereact_decade = module.computeE_react(2000.0);
    std::cout << "Decade κ = " << std::scientific << kappa_decade << " day⁻¹\n";
    std::cout << "Timescale = " << std::fixed << std::setprecision(1) << timescale_decade << " years\n";
    std::cout << "E_react at t=2000 days = " << std::scientific << ereact_decade << " J\n\n";
    
    // Step 7: Restore and expand decay scale
    std::cout << "STEP 7: Expand Decay Scale (κ x2, faster decay)\n";
    module.restoreState("standard_5_5years");
    module.expandDecayScale(2.0, 1.0);
    double kappa_expanded = module.getVariable("kappa_day");
    double ereact_expanded = module.computeE_react(2000.0);
    std::cout << "Expanded κ = " << std::scientific << kappa_expanded << " day⁻¹\n";
    std::cout << "New timescale = " << std::fixed << std::setprecision(1) 
              << ((1.0/kappa_expanded)/365.25) << " years\n";
    std::cout << "E_react at t=2000 days = " << std::scientific << ereact_expanded << " J\n\n";
    
    // Step 8: Restore and expand energy scale
    std::cout << "STEP 8: Expand Energy Scale (E_react_base x2)\n";
    module.restoreState("standard_5_5years");
    module.expandEnergyScale(1.0, 2.0);
    double ebase_new = module.getVariable("E_react_base");
    double ereact_scaled = module.computeE_react(2000.0);
    std::cout << "New E_react_base = " << std::scientific << ebase_new << " J\n";
    std::cout << "E_react at t=2000 days = " << ereact_scaled << " J\n";
    std::cout << "Scaling factor: " << std::fixed << std::setprecision(1) << (ebase_new / 1e46) << "x\n\n";
    
    // Step 9: Restore and expand magnetic scale
    std::cout << "STEP 9: Expand Magnetic Scale (μ/r x2, amplifications x1.5)\n";
    module.restoreState("standard_5_5years");
    module.expandMagneticScale(2.0, 1.5);
    double mu_new = module.getVariable("mu_over_rj");
    double hf_new = module.getVariable("heaviside_f");
    double um_new = module.computeUmExample(2000.0);
    std::cout << "New μ_j/r_j = " << std::scientific << mu_new << " T·m²\n";
    std::cout << "New Heaviside = " << hf_new << "\n";
    std::cout << "U_m at t=2000 days = " << um_new << " J/m³\n\n";
    
    // Step 10: Sensitivity analysis
    std::cout << "STEP 10: Sensitivity Analysis (at t=2000 days)\n";
    module.restoreState("standard_5_5years");
    std::vector<std::string> params = {"kappa_day", "E_react_base", "t_day", "mu_over_rj"};
    auto sensitivities = module.sensitivityAnalysis(params);
    for (const auto& pair : sensitivities) {
        std::cout << "  ∂E_react/∂" << pair.first << " ≈ " << std::scientific << pair.second << " (normalized)\n";
    }
    std::cout << "\n";
    
    // Step 11: Generate variations
    std::cout << "STEP 11: Generate Parameter Variations (5 variants, ±10%)\n";
    auto variations = module.generateVariations(5, 0.1);
    for (int i = 0; i < variations.size(); ++i) {
        double var_kappa = variations[i]["kappa_day"];
        double var_timescale = (1.0 / var_kappa) / 365.25;
        std::cout << "  Variant " << (i+1) << ": κ=" << std::scientific << var_kappa 
                  << " day⁻¹ (timescale=" << std::fixed << std::setprecision(1) 
                  << var_timescale << " years)\n";
    }
    std::cout << "\n";
    
    // Step 12: Auto-refine to target decay fraction
    std::cout << "STEP 12: Auto-Refine to Target Decay Fraction (50% at t=2000 days)\n";
    module.restoreState("standard_5_5years");
    module.autoRefineParameters("decay_fraction", 0.5);
    double refined_kappa = module.getVariable("kappa_day");
    double refined_ereact = module.computeE_react(2000.0);
    double refined_fraction = refined_ereact / module.getVariable("E_react_base");
    std::cout << "Refined κ = " << std::scientific << refined_kappa << " day⁻¹\n";
    std::cout << "E_react at t=2000 days = " << refined_ereact << " J\n";
    std::cout << "Achieved fraction = " << std::fixed << std::setprecision(3) 
              << refined_fraction << " (target: 0.5)\n";
    std::cout << "Timescale = " << std::fixed << std::setprecision(1) 
              << ((1.0/refined_kappa)/365.25) << " years\n\n";
    
    // Step 13: Auto-refine to target timescale
    std::cout << "STEP 13: Auto-Refine to Target Timescale (3 years)\n";
    module.restoreState("standard_5_5years");
    module.autoRefineParameters("timescale_years", 3.0);
    double kappa_3y = module.getVariable("kappa_day");
    double timescale_3y = (1.0 / kappa_3y) / 365.25;
    std::cout << "Refined κ = " << std::scientific << kappa_3y << " day⁻¹\n";
    std::cout << "Achieved timescale = " << std::fixed << std::setprecision(1) 
              << timescale_3y << " years (target: 3.0)\n\n";
    
    // Step 14: Auto-refine to target E_react at specific time
    std::cout << "STEP 14: Auto-Refine to Target E_react (5e45 J at t=2000 days)\n";
    module.restoreState("standard_5_5years");
    module.autoRefineParameters("E_react_at_time", 5e45);
    double refined_ebase = module.getVariable("E_react_base");
    double achieved_ereact = module.computeE_react(2000.0);
    std::cout << "Refined E_react_base = " << std::scientific << refined_ebase << " J\n";
    std::cout << "Achieved E_react at t=2000 = " << achieved_ereact << " J (target: 5e45)\n\n";
    
    // Step 15: Calibrate to observations
    std::cout << "STEP 15: Calibrate to Observational Data\n";
    module.restoreState("standard_5_5years");
    std::map<std::string, double> observations;
    observations["kappa_day"] = 0.001;  // Observed faster decay
    observations["E_react_base"] = 1.5e46;  // Observed higher base energy
    module.calibrateToObservations(observations);
    std::cout << "Calibrated κ = " << std::scientific << module.getVariable("kappa_day") << " day⁻¹\n";
    std::cout << "Calibrated E_base = " << module.getVariable("E_react_base") << " J\n";
    std::cout << "Timescale = " << std::fixed << std::setprecision(1) 
              << ((1.0/module.getVariable("kappa_day"))/365.25) << " years\n\n";
    
    // Step 16: Optimize for metric
    std::cout << "STEP 16: Optimize for 'standard_5_5years' Metric\n";
    module.optimizeForMetric("standard_5_5years");
    std::cout << "Optimized κ = " << std::scientific << module.getVariable("kappa_day") << " day⁻¹\n";
    std::cout << "Timescale = " << std::fixed << std::setprecision(1) 
              << ((1.0/module.getVariable("kappa_day"))/365.25) << " years\n\n";
    
    // Step 17: Mutate parameters
    std::cout << "STEP 17: Mutate Parameters (5% random variation)\n";
    module.mutateParameters(0.05);
    std::cout << "Mutated κ = " << std::scientific << module.getVariable("kappa_day") << " day⁻¹\n";
    std::cout << "Mutated timescale = " << std::fixed << std::setprecision(1) 
              << ((1.0/module.getVariable("kappa_day"))/365.25) << " years\n\n";
    
    // Step 18: Validate consistency
    std::cout << "STEP 18: Validate Consistency\n";
    bool valid = module.validateConsistency();
    std::cout << "Consistency check: " << (valid ? "PASSED" : "FAILED") << "\n\n";
    
    // Step 19: Introduce anomaly and auto-correct
    std::cout << "STEP 19: Introduce Anomaly and Auto-Correct\n";
    module.createVariable("kappa_anomaly", -0.001);  // Invalid negative κ
    module.removeVariable("kappa_day");
    module.createVariable("kappa_day", -0.001);
    std::cout << "Introduced invalid κ = " << module.getVariable("kappa_day") << " day⁻¹\n";
    module.autoCorrectAnomalies();
    std::cout << "Auto-corrected κ = " << std::scientific << module.getVariable("kappa_day") << " day⁻¹\n";
    std::cout << "Timescale = " << std::fixed << std::setprecision(1) 
              << ((1.0/module.getVariable("kappa_day"))/365.25) << " years\n";
    bool valid_after = module.validateConsistency();
    std::cout << "Consistency after correction: " << (valid_after ? "PASSED" : "FAILED") << "\n\n";
    
    // Step 20: List saved states and export
    std::cout << "STEP 20: List Saved States and Export Final State\n";
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
// #include "ScmReactivityDecayModule.h"
// int main() {
//     ScmReactivityDecayModule mod;
//     double kappa = mod.computeKappa_day();
//     std::cout << "κ = " << kappa << " day⁻¹\n";
//     mod.printDecayEffects(2000.0);
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("kappa_day", 0.001);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o scm_decay_test scm_decay_test.cpp ScmReactivityDecayModule.cpp -lm
// Sample: κ=5e-4 day⁻¹; t=2000 days: E_react≈3.68e45 J; U_m≈8.39e64 J/m³.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

ScmReactivityDecayModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeKappa_day, computeKappa_s, computeE_react, computeUmExample) are clear, concise, and variable - driven.
- Automatic recalculation of derived variables(kappa_s) when dependencies change.
- Output and debugging functions(printVariables, printDecayEffects, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in[SCm] reactivity decay modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.