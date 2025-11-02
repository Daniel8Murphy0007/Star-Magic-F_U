// CorePenetrationModule.h
// Modular C++ implementation of the Planetary Core Penetration Factor (P_core) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes P_core ?1 (unitless for Sun, ~1e-3 for planets); scales P_core in Universal Gravity U_g3 term.
// Pluggable: #include "CorePenetrationModule.h"
// CorePenetrationModule mod; mod.computeU_g3(0.0); mod.updateVariable("P_core", new_value);
// Variables in std::map; example for Sun at t=0; planet mode with P_core=1e-3.
// Approximations: cos(?_s t ?)=1 at t=0; E_react=1e46; B_j=1e3 T.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef CORE_PENETRATION_MODULE_H
#define CORE_PENETRATION_MODULE_H

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

class CorePenetrationModule {
private:
    std::map<std::string, double> variables;
    double computeU_g3(double t);

public:
    // Constructor: Initialize with framework defaults (Sun)
    CorePenetrationModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeP_core();  // ?1 for Sun (unitless)
    double computeU_g3(double t);  // U_g3 with P_core (J/m^3)
    double computeU_g3_planet(double t);  // For planet P_core=1e-3

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
    void expandParameterSpace(double pcore_scale, double magnetic_scale, double energy_scale);
    void expandCoreScale(double pcore_factor, double planet_factor);        // Scale P_core (stellar vs planetary)
    void expandMagneticDiskScale(double bj_factor, double k3_factor);       // Scale B_j and k_3
    void expandOscillationScale(double omega_factor, double cos_factor);    // Scale ω_s and oscillations

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

#endif // CORE_PENETRATION_MODULE_H

// CorePenetrationModule.cpp
#include "CorePenetrationModule.h"

// Constructor: Set framework defaults (Sun at t=0)
CorePenetrationModule::CorePenetrationModule() {
    // Universal constants
    variables["P_core"] = 1.0;                      // Unitless ?1 for Sun
    variables["k_3"] = 1.8;                         // Coupling
    variables["B_j"] = 1e3;                         // T (base)
    variables["omega_s"] = 2.5e-6;                  // rad/s
    variables["P_core_planet"] = 1e-3;              // For planets
    variables["E_react"] = 1e46;                    // J
    variables["pi"] = 3.141592653589793;
    variables["t"] = 0.0;                           // s
}

// Update variable
void CorePenetrationModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void CorePenetrationModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void CorePenetrationModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute P_core ?1
double CorePenetrationModule::computeP_core() {
    return variables["P_core"];
}

// Compute U_g3 with P_core
double CorePenetrationModule::computeU_g3(double t) {
    variables["t"] = t;
    double k_3 = variables["k_3"];
    double b_j = variables["B_j"];  // Simplified, no sin term in example
    double cos_term = std::cos(variables["omega_s"] * t * variables["pi"]);
    double p_core = computeP_core();
    double e_react = variables["E_react"];
    return k_3 * b_j * cos_term * p_core * e_react;
}

// U_g3 for planet (P_core=1e-3)
double CorePenetrationModule::computeU_g3_planet(double t) {
    double orig_p = variables["P_core"];
    variables["P_core"] = variables["P_core_planet"];
    double result = computeU_g3(t);
    variables["P_core"] = orig_p;
    return result;
}

// Equation text
std::string CorePenetrationModule::getEquationText() {
    return "U_g3 = k_3 * ?_j B_j(r,?,t,?_vac,[SCm]) * cos(?_s(t) t ?) * P_core * E_react\n"
           "Where P_core ?1 (unitless for Sun, ~1e-3 for planets; core penetration).\n"
           "Scales magnetic disk gravity for core [SCm] influence.\n"
           "Example Sun t=0: U_g3 ?1.8e49 J/m� (P_core=1);\n"
           "Planet: ?1.8e46 J/m� (P_core=1e-3, -3 orders).\n"
           "Role: Adjusts core interactions; full for stellar plasma, reduced for solid cores.\n"
           "UQFF: Models penetration in nebulae/star formation/disks.";
}

// Print variables
void CorePenetrationModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// ===== ENHANCED METHODS IMPLEMENTATION =====

namespace core_penetration_saved_states {
    static std::map<std::string, std::map<std::string, double>> saved_states;
}

// Variable Management
void CorePenetrationModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void CorePenetrationModule::removeVariable(const std::string& name) {
    variables.erase(name);
}

void CorePenetrationModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
    }
}

std::vector<std::string> CorePenetrationModule::listVariables() const {
    std::vector<std::string> names;
    for (const auto& pair : variables) {
        names.push_back(pair.first);
    }
    return names;
}

std::string CorePenetrationModule::getSystemName() const {
    return "Core_Penetration_UQFF";
}

// Batch Operations
void CorePenetrationModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = func(variables[name]);
        }
    }
}

void CorePenetrationModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// Self-Expansion: Domain-Specific Scales
void CorePenetrationModule::expandParameterSpace(double pcore_scale, double magnetic_scale, double energy_scale) {
    variables["P_core"] *= pcore_scale;
    variables["P_core_planet"] *= pcore_scale;
    variables["B_j"] *= magnetic_scale;
    variables["k_3"] *= magnetic_scale;
    variables["E_react"] *= energy_scale;
}

void CorePenetrationModule::expandCoreScale(double pcore_factor, double planet_factor) {
    variables["P_core"] *= pcore_factor;
    variables["P_core_planet"] *= planet_factor;
}

void CorePenetrationModule::expandMagneticDiskScale(double bj_factor, double k3_factor) {
    variables["B_j"] *= bj_factor;
    variables["k_3"] *= k3_factor;
}

void CorePenetrationModule::expandOscillationScale(double omega_factor, double cos_factor) {
    variables["omega_s"] *= omega_factor;
    // cos_factor could adjust phase or amplitude in future extensions
}

// Self-Refinement
void CorePenetrationModule::autoRefineParameters(const std::string& target, double goal) {
    if (target == "U_g3") {
        // Target specific U_g3 by scaling E_react
        double current_ug3 = computeU_g3(variables["t"]);
        if (std::abs(current_ug3) > 1e-9) {
            variables["E_react"] *= (goal / current_ug3);
        }
    } else if (target == "P_core") {
        // Target specific P_core directly
        variables["P_core"] = goal;
    } else if (target == "stellar_to_planetary_ratio") {
        // Target specific ratio between stellar and planetary P_core
        variables["P_core_planet"] = variables["P_core"] / goal;
    } else if (target == "cos_term") {
        // Target specific cos(ω_s t π) by adjusting t
        double desired_cos = goal;
        if (std::abs(desired_cos) <= 1.0) {
            double arg = std::acos(desired_cos);
            variables["t"] = arg / (variables["omega_s"] * variables["pi"]);
        }
    }
}

void CorePenetrationModule::calibrateToObservations(const std::map<std::string, double>& observations) {
    for (const auto& obs : observations) {
        if (variables.find(obs.first) != variables.end()) {
            variables[obs.first] = obs.second;
        }
    }
}

void CorePenetrationModule::optimizeForMetric(const std::string& metric) {
    if (metric == "maximize_stellar") {
        // Maximize stellar core penetration
        variables["P_core"] = 1.0;
    } else if (metric == "maximize_planetary") {
        // Maximize planetary core penetration (still << stellar)
        variables["P_core_planet"] = 5e-3;
    } else if (metric == "enhance_magnetic_disk") {
        // Enhance magnetic disk contribution
        variables["k_3"] *= 1.3;
        variables["B_j"] *= 1.2;
    } else if (metric == "reduce_planetary_gap") {
        // Reduce the gap between stellar and planetary
        double stellar = variables["P_core"];
        double planetary = variables["P_core_planet"];
        variables["P_core_planet"] = std::sqrt(stellar * planetary);  // Geometric mean
    }
}

// Parameter Exploration
std::vector<std::map<std::string, double>> CorePenetrationModule::generateVariations(int count, double variation_pct) {
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
void CorePenetrationModule::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - mutation_rate, 1.0 + mutation_rate);
    
    for (auto& pair : variables) {
        if (pair.first != "pi") {
            pair.second *= dis(gen);
        }
    }
}

void CorePenetrationModule::evolveSystem(int generations, std::function<double()> fitness_func) {
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
void CorePenetrationModule::saveState(const std::string& label) {
    core_penetration_saved_states::saved_states[label] = variables;
}

void CorePenetrationModule::restoreState(const std::string& label) {
    if (core_penetration_saved_states::saved_states.find(label) != core_penetration_saved_states::saved_states.end()) {
        variables = core_penetration_saved_states::saved_states[label];
    }
}

std::vector<std::string> CorePenetrationModule::listSavedStates() const {
    std::vector<std::string> labels;
    for (const auto& pair : core_penetration_saved_states::saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string CorePenetrationModule::exportState() const {
    std::ostringstream oss;
    oss << "CorePenetration_State_Export:\n";
    for (const auto& pair : variables) {
        oss << pair.first << "=" << std::scientific << pair.second << "\n";
    }
    return oss.str();
}

// System Analysis
std::map<std::string, double> CorePenetrationModule::sensitivityAnalysis(const std::vector<std::string>& params) {
    std::map<std::string, double> sensitivities;
    double baseline_ug3 = computeU_g3(variables["t"]);
    
    for (const auto& param : params) {
        if (variables.find(param) != variables.end() && param != "pi") {
            double original = variables[param];
            variables[param] = original * 1.01;
            double perturbed_ug3 = computeU_g3(variables["t"]);
            sensitivities[param] = (perturbed_ug3 - baseline_ug3) / baseline_ug3;
            variables[param] = original;
        }
    }
    return sensitivities;
}

std::string CorePenetrationModule::generateReport() const {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(3);
    oss << "===== Core Penetration Module Report =====\n";
    oss << "System: " << getSystemName() << "\n\n";
    
    oss << "Core Penetration Parameters:\n";
    oss << "  P_core (stellar) = " << variables.at("P_core") << " (unitless)\n";
    oss << "  P_core (planetary) = " << variables.at("P_core_planet") << " (unitless)\n";
    double ratio = variables.at("P_core") / variables.at("P_core_planet");
    oss << "  Stellar/Planetary ratio = " << ratio << ":1\n";
    oss << "  Orders of magnitude difference: " << std::log10(ratio) << "\n\n";
    
    oss << "Magnetic Disk Parameters:\n";
    oss << "  k_3 (coupling) = " << variables.at("k_3") << "\n";
    oss << "  B_j (magnetic field) = " << variables.at("B_j") << " T\n";
    oss << "  ω_s (solar angular frequency) = " << variables.at("omega_s") << " rad/s\n";
    double period_s = 2.0 * variables.at("pi") / variables.at("omega_s");
    double period_days = period_s / 86400.0;
    oss << "  Period = 2π/ω_s = " << period_days << " days\n\n";
    
    oss << "Energy Parameters:\n";
    oss << "  E_react = " << variables.at("E_react") << " J\n";
    oss << "  Time t = " << variables.at("t") << " s\n\n";
    
    oss << "U_g3 Computations (Magnetic Disk Gravity):\n";
    double t = variables.at("t");
    double cos_term = std::cos(variables.at("omega_s") * t * variables.at("pi"));
    oss << "  cos(ω_s t π) = " << cos_term << "\n";
    
    // Stellar U_g3
    double k_3 = variables.at("k_3");
    double b_j = variables.at("B_j");
    double p_core_stellar = variables.at("P_core");
    double e_react = variables.at("E_react");
    double ug3_stellar = k_3 * b_j * cos_term * p_core_stellar * e_react;
    oss << "  U_g3 (stellar, P_core=" << p_core_stellar << ") = " << ug3_stellar << " J/m³\n";
    
    // Planetary U_g3
    double p_core_planet = variables.at("P_core_planet");
    double ug3_planet = k_3 * b_j * cos_term * p_core_planet * e_react;
    oss << "  U_g3 (planetary, P_core=" << p_core_planet << ") = " << ug3_planet << " J/m³\n";
    oss << "  U_g3 stellar/planetary ratio = " << (ug3_stellar / ug3_planet) << ":1\n\n";
    
    oss << "Physical Interpretation:\n";
    if (p_core_stellar >= 0.9) {
        oss << "  Stellar mode: Full core penetration (plasma core)\n";
    } else if (p_core_stellar >= 0.5) {
        oss << "  Hybrid mode: Partial core penetration\n";
    } else {
        oss << "  Reduced mode: Limited core penetration\n";
    }
    
    if (p_core_planet >= 1e-2) {
        oss << "  Planetary mode: Enhanced solid core coupling\n";
    } else if (p_core_planet >= 1e-4) {
        oss << "  Planetary mode: Typical solid core (~1e-3)\n";
    } else {
        oss << "  Planetary mode: Very weak core coupling\n";
    }
    
    oss << "  Applications:\n";
    oss << "    - Stellar: Core [SCm] in Sun, fully coupled to magnetic disk\n";
    oss << "    - Planetary: Solid cores (Earth, Jupiter) with reduced [SCm] influence\n";
    oss << "    - Star formation: Nebular core penetration modeling\n";
    oss << "    - AGN/Quasars: Central engine core coupling\n";
    
    return oss.str();
}

bool CorePenetrationModule::validateConsistency() const {
    bool valid = true;
    
    // Check P_core is in reasonable range [0, 1]
    if (variables.find("P_core") != variables.end()) {
        double pcore = variables.at("P_core");
        if (pcore < 0 || pcore > 1.5) {
            std::cerr << "Warning: P_core outside typical range [0, 1.5] (current: " << pcore << ")\n";
        }
    }
    
    // Check P_core_planet << P_core
    if (variables.find("P_core") != variables.end() && variables.find("P_core_planet") != variables.end()) {
        double stellar = variables.at("P_core");
        double planetary = variables.at("P_core_planet");
        if (planetary >= stellar) {
            std::cerr << "Error: P_core_planet >= P_core (planetary should be much less than stellar)\n";
            valid = false;
        }
    }
    
    // Check k_3 is positive
    if (variables.find("k_3") != variables.end() && variables.at("k_3") <= 0) {
        std::cerr << "Error: k_3 <= 0 (coupling constant must be positive)\n";
        valid = false;
    }
    
    // Check E_react is positive
    if (variables.find("E_react") != variables.end() && variables.at("E_react") <= 0) {
        std::cerr << "Error: E_react <= 0 (reactor energy must be positive)\n";
        valid = false;
    }
    
    return valid;
}

void CorePenetrationModule::autoCorrectAnomalies() {
    // Reset P_core to safe defaults
    if (variables["P_core"] < 0 || variables["P_core"] > 1.5) {
        variables["P_core"] = 1.0;
    }
    
    // Ensure planetary << stellar
    if (variables["P_core_planet"] >= variables["P_core"]) {
        variables["P_core_planet"] = variables["P_core"] * 1e-3;
    }
    
    // Ensure positive k_3
    if (variables["k_3"] <= 0) {
        variables["k_3"] = 1.8;
    }
    
    // Ensure positive E_react
    if (variables["E_react"] <= 0) {
        variables["E_react"] = 1e46;
    }
}

// Example usage in base program (snippet)
// #include "CorePenetrationModule.h"
// int main() {
//     CorePenetrationModule mod;
//     double p = mod.computeP_core();
//     std::cout << "P_core ? " << p << std::endl;
//     double u_g3_sun = mod.computeU_g3(0.0);
//     std::cout << "U_g3 (Sun) = " << u_g3_sun << " J/m�\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("P_core", 1e-3);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o core_test core_test.cpp CorePenetrationModule.cpp -lm
// Sample: P_core=1; U_g3?1.8e49 J/m� (Sun); scales for planetary cores.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

// ===== Enhanced Example: 18-Step Demonstration of Dynamic Capabilities =====
// int main() {
//     CorePenetrationModule mod;
//     std::cout << "===== Core Penetration Module: Enhanced 18-Step Demo =====\n\n";
//
//     // Step 1: Initial report
//     std::cout << "Step 1: Initial Configuration (Stellar Mode, P_core=1)\n";
//     std::cout << mod.generateReport() << "\n";
//
//     // Step 2: Track key core penetration quantities
//     std::cout << "Step 2: Create Tracking Variables\n";
//     mod.createVariable("P_core_baseline", mod.computeP_core());
//     mod.createVariable("U_g3_stellar_baseline", mod.computeU_g3(0.0));
//     mod.createVariable("U_g3_planet_baseline", mod.computeU_g3_planet(0.0));
//     double ratio = mod.variables["U_g3_stellar_baseline"] / mod.variables["U_g3_planet_baseline"];
//     mod.createVariable("stellar_planet_ratio", ratio);
//     std::cout << "  P_core (stellar) = " << mod.variables["P_core_baseline"] << "\n";
//     std::cout << "  U_g3 (stellar) = " << std::scientific << mod.variables["U_g3_stellar_baseline"] << " J/m³\n";
//     std::cout << "  U_g3 (planetary) = " << mod.variables["U_g3_planet_baseline"] << " J/m³\n";
//     std::cout << "  Ratio = " << ratio << ":1 (stellar " << std::log10(ratio) << " orders higher)\n\n";
//
//     // Step 3: P_core variations (stellar scenarios)
//     std::cout << "Step 3: P_core Variations (Stellar Scenarios)\n";
//     mod.saveState("baseline");
//     std::vector<double> pcore_stellar = {0.5, 0.7, 0.85, 1.0, 1.1, 1.2};
//     for (double pc : pcore_stellar) {
//         mod.updateVariable("P_core", pc);
//         double ug3 = mod.computeU_g3(0.0);
//         std::cout << "  P_core=" << pc << ": U_g3=" << ug3 << " J/m³\n";
//     }
//     mod.restoreState("baseline");
//     std::cout << "\n";
//
//     // Step 4: P_core variations (planetary scenarios)
//     std::cout << "Step 4: P_core Variations (Planetary Scenarios)\n";
//     std::vector<double> pcore_planet = {1e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2};
//     for (double pc : pcore_planet) {
//         mod.updateVariable("P_core_planet", pc);
//         double ug3 = mod.computeU_g3_planet(0.0);
//         std::cout << "  P_core_planet=" << pc << ": U_g3=" << ug3 << " J/m³\n";
//     }
//     mod.restoreState("baseline");
//     std::cout << "\n";
//
//     // Step 5: Stellar vs planetary comparison
//     std::cout << "Step 5: Stellar vs Planetary Comparison (Various Objects)\n";
//     std::vector<std::pair<std::string, std::pair<double, double>>> objects = {
//         {"Sun", {1.0, 0.0}},
//         {"Jupiter", {0.1, 5e-3}},
//         {"Earth", {0.0, 1e-3}},
//         {"Neutron Star", {1.2, 0.0}},
//         {"Proto-planet", {0.5, 8e-3}}
//     };
//     for (const auto& obj : objects) {
//         mod.updateVariable("P_core", obj.second.first);
//         mod.updateVariable("P_core_planet", obj.second.second);
//         double ug3 = (obj.second.first > 0.01) ? mod.computeU_g3(0.0) : mod.computeU_g3_planet(0.0);
//         std::cout << "  " << obj.first << " (P_core=" << obj.second.first 
//                   << ", planet=" << obj.second.second << "): U_g3=" << ug3 << " J/m³\n";
//     }
//     mod.restoreState("baseline");
//     std::cout << "\n";
//
//     // Step 6: Time evolution (oscillation analysis)
//     std::cout << "Step 6: Time Evolution (cos(ω_s t π) Oscillation)\n";
//     std::vector<double> time_points = {0.0, 1e5, 2e5, 5e5, 1e6, 2e6};
//     for (double t : time_points) {
//         double cos_val = std::cos(mod.variables["omega_s"] * t * mod.variables["pi"]);
//         double ug3 = mod.computeU_g3(t);
//         std::cout << "  t=" << t << " s: cos(ω_s t π)=" << cos_val << ", U_g3=" << ug3 << " J/m³\n";
//     }
//     std::cout << "\n";
//
//     // Step 7: Expand core scale
//     std::cout << "Step 7: Expand Core Scale (P_core x1.1, P_core_planet x1.5)\n";
//     mod.restoreState("baseline");
//     mod.saveState("pre_core_expansion");
//     double pcore_before = mod.variables["P_core"];
//     double pplanet_before = mod.variables["P_core_planet"];
//     mod.expandCoreScale(1.1, 1.5);
//     std::cout << "  P_core: " << pcore_before << " -> " << mod.variables["P_core"] << "\n";
//     std::cout << "  P_core_planet: " << pplanet_before << " -> " << mod.variables["P_core_planet"] << "\n";
//     std::cout << "  New U_g3 (stellar): " << mod.computeU_g3(0.0) << " J/m³\n";
//     std::cout << "  New U_g3 (planetary): " << mod.computeU_g3_planet(0.0) << " J/m³\n\n";
//
//     // Step 8: Expand magnetic disk scale
//     std::cout << "Step 8: Expand Magnetic Disk Scale (B_j x1.3, k_3 x1.2)\n";
//     mod.restoreState("pre_core_expansion");
//     double bj_before = mod.variables["B_j"];
//     double k3_before = mod.variables["k_3"];
//     mod.expandMagneticDiskScale(1.3, 1.2);
//     std::cout << "  B_j: " << bj_before << " -> " << mod.variables["B_j"] << " T\n";
//     std::cout << "  k_3: " << k3_before << " -> " << mod.variables["k_3"] << "\n";
//     std::cout << "  New U_g3: " << mod.computeU_g3(0.0) << " J/m³\n\n";
//
//     // Step 9: Expand oscillation scale
//     std::cout << "Step 9: Expand Oscillation Scale (ω_s x1.5)\n";
//     mod.restoreState("baseline");
//     double omega_before = mod.variables["omega_s"];
//     double period_before = 2.0 * mod.variables["pi"] / omega_before / 86400.0;
//     mod.expandOscillationScale(1.5, 1.0);
//     double omega_after = mod.variables["omega_s"];
//     double period_after = 2.0 * mod.variables["pi"] / omega_after / 86400.0;
//     std::cout << "  ω_s: " << omega_before << " -> " << omega_after << " rad/s\n";
//     std::cout << "  Period: " << period_before << " -> " << period_after << " days\n\n";
//
//     // Step 10: Parameter variations
//     std::cout << "Step 10: Generate 10 Parameter Variations (±10%)\n";
//     auto variations = mod.generateVariations(10, 0.10);
//     std::cout << "  Generated " << variations.size() << " configurations\n";
//     std::vector<double> ug3_range, pcore_range;
//     for (const auto& var : variations) {
//         CorePenetrationModule temp_mod;
//         temp_mod.variables = var;
//         ug3_range.push_back(temp_mod.computeU_g3(0.0));
//         pcore_range.push_back(var.at("P_core"));
//     }
//     auto ug3_minmax = std::minmax_element(ug3_range.begin(), ug3_range.end());
//     auto pcore_minmax = std::minmax_element(pcore_range.begin(), pcore_range.end());
//     std::cout << "  U_g3 range: " << *ug3_minmax.first << " to " << *ug3_minmax.second << " J/m³\n";
//     std::cout << "  P_core range: " << *pcore_minmax.first << " to " << *pcore_minmax.second << "\n\n";
//
//     // Step 11: Sensitivity analysis
//     std::cout << "Step 11: Sensitivity Analysis (U_g3 response to ±1% changes)\n";
//     std::vector<std::string> sens_params = {"P_core", "B_j", "k_3", "E_react", "omega_s"};
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
//     // Step 12: Auto-refine to target U_g3
//     std::cout << "Step 12: Auto-Refine to Target U_g3 = 2.5e49 J/m³\n";
//     mod.restoreState("baseline");
//     double ug3_before = mod.computeU_g3(0.0);
//     mod.autoRefineParameters("U_g3", 2.5e49);
//     double ug3_after = mod.computeU_g3(0.0);
//     std::cout << "  Before: " << ug3_before << " J/m³\n";
//     std::cout << "  After: " << ug3_after << " J/m³\n";
//     std::cout << "  New E_react: " << mod.variables["E_react"] << " J\n\n";
//
//     // Step 13: Target specific P_core
//     std::cout << "Step 13: Target P_core = 0.8 (Partially Penetrating Core)\n";
//     mod.restoreState("baseline");
//     double pc_before = mod.variables["P_core"];
//     mod.autoRefineParameters("P_core", 0.8);
//     double pc_after = mod.variables["P_core"];
//     std::cout << "  P_core: " << pc_before << " -> " << pc_after << "\n";
//     std::cout << "  New U_g3: " << mod.computeU_g3(0.0) << " J/m³\n\n";
//
//     // Step 14: Calibration to observations
//     std::cout << "Step 14: Calibrate to Observational Data (Giant Planet Core)\n";
//     mod.restoreState("baseline");
//     std::map<std::string, double> observations = {
//         {"P_core", 0.2},            // Weak stellar component
//         {"P_core_planet", 3e-3},    // Enhanced planetary
//         {"B_j", 1.5e3},             // Stronger field
//         {"k_3", 2.0}                // Higher coupling
//     };
//     mod.calibrateToObservations(observations);
//     std::cout << "  Calibrated P_core: " << mod.variables["P_core"] << "\n";
//     std::cout << "  Calibrated P_core_planet: " << mod.variables["P_core_planet"] << "\n";
//     std::cout << "  Calibrated B_j: " << mod.variables["B_j"] << " T\n";
//     std::cout << "  New U_g3 (mixed): " << mod.computeU_g3(0.0) << " J/m³\n\n";
//
//     // Step 15: Optimize for stellar dominance
//     std::cout << "Step 15: Optimize for Maximum Stellar Core Penetration\n";
//     mod.restoreState("baseline");
//     double pc_pre_stellar = mod.variables["P_core"];
//     mod.optimizeForMetric("maximize_stellar");
//     double pc_post_stellar = mod.variables["P_core"];
//     std::cout << "  P_core: " << pc_pre_stellar << " -> " << pc_post_stellar << "\n";
//     std::cout << "  New U_g3: " << mod.computeU_g3(0.0) << " J/m³\n\n";
//
//     // Step 16: Optimize for enhanced magnetic disk
//     std::cout << "Step 16: Optimize for Enhanced Magnetic Disk\n";
//     mod.restoreState("baseline");
//     double k3_pre = mod.variables["k_3"];
//     double bj_pre = mod.variables["B_j"];
//     double ug3_pre = mod.computeU_g3(0.0);
//     mod.optimizeForMetric("enhance_magnetic_disk");
//     double ug3_post = mod.computeU_g3(0.0);
//     std::cout << "  k_3: " << k3_pre << " -> " << mod.variables["k_3"] << "\n";
//     std::cout << "  B_j: " << bj_pre << " -> " << mod.variables["B_j"] << " T\n";
//     std::cout << "  U_g3: " << ug3_pre << " -> " << ug3_post << " J/m³\n";
//     std::cout << "  Enhancement: " << ((ug3_post / ug3_pre - 1.0) * 100) << "%\n\n";
//
//     // Step 17: System evolution (maximize U_g3)
//     std::cout << "Step 17: Evolve System (8 generations, maximize U_g3)\n";
//     mod.restoreState("baseline");
//     double initial_fitness = mod.computeU_g3(0.0);
//     mod.evolveSystem(8, [&mod]() { return mod.computeU_g3(0.0); });
//     double final_fitness = mod.computeU_g3(0.0);
//     std::cout << "  Initial U_g3: " << initial_fitness << " J/m³\n";
//     std::cout << "  Evolved U_g3: " << final_fitness << " J/m³\n";
//     std::cout << "  Improvement: " << ((final_fitness / initial_fitness - 1.0) * 100) << "%\n";
//     std::cout << "  Final P_core: " << mod.variables["P_core"] << "\n\n";
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
//     std::cout << "  Final P_core: " << mod.variables["P_core"] << "\n";
//     std::cout << "  Final P_core_planet: " << mod.variables["P_core_planet"] << "\n";
//     std::cout << "  Final U_g3 (stellar): " << mod.computeU_g3(0.0) << " J/m³\n";
//     std::cout << "  Final U_g3 (planetary): " << mod.computeU_g3_planet(0.0) << " J/m³\n\n";
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

CorePenetrationModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeP_core, computeU_g3, computeU_g3_planet) are clear, concise, and variable - driven.
- Output and debugging functions(printVariables, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.
- Supports both stellar and planetary scenarios by switching P_core values.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in core penetration modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.