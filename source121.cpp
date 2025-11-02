// SurfaceMagneticFieldModule.h
// Modular C++ implementation of the Surface Magnetic Field (B_s) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes B_s range [1e-4, 0.4] T for Sun; influences B_j in U_g3 magnetic strings (scaled by B_s / B_ref).
// Pluggable: #include "SurfaceMagneticFieldModule.h"
// SurfaceMagneticFieldModule mod; mod.computeU_g3_example(0.0); mod.updateVariable("B_s_min", new_value);
// Variables in std::map; example for Sun at t=0; quiet Sun B_s=1e-4 T ? U_g3?4.5e45 J/m�.
// Approximations: B_ref=0.4 T (max sunspot); cos(?_s t ?)=1; P_core=1; E_react=1e46.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef SURFACE_MAGNETIC_FIELD_MODULE_H
#define SURFACE_MAGNETIC_FIELD_MODULE_H

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

class SurfaceMagneticFieldModule {
private:
    std::map<std::string, double> variables;
    double computeB_j(double t, double B_s);
    double computeU_g3_example(double t, double B_s);

public:
    // Constructor: Initialize with framework defaults (Sun)
    SurfaceMagneticFieldModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeB_s_min();  // 1e-4 T (quiet Sun)
    double computeB_s_max();  // 0.4 T (sunspot max)
    double computeB_j(double t, double B_s);  // Scaled B_j (T)
    double computeU_g3_example(double t, double B_s);  // U_g3 with B_j (J/m^3)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // ===== ENHANCED METHODS =====
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
    void expandParameterSpace(double field_scale, double coupling_scale, double energy_scale);
    void expandFieldScale(double bs_factor, double bj_factor);
    void expandCouplingScale(double k3_factor, double interaction_factor);
    void expandEnergyScale(double ereact_factor, double ug3_factor);

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

#endif // SURFACE_MAGNETIC_FIELD_MODULE_H

// SurfaceMagneticFieldModule.cpp
#include "SurfaceMagneticFieldModule.h"

// Constructor: Set framework defaults (Sun)
SurfaceMagneticFieldModule::SurfaceMagneticFieldModule() {
    // Universal constants
    variables["B_s_min"] = 1e-4;                    // T (quiet)
    variables["B_s_max"] = 0.4;                     // T (sunspot)
    variables["B_ref"] = 0.4;                       // T (reference max)
    variables["k_3"] = 1.8;                         // Coupling
    variables["omega_s"] = 2.5e-6;                  // rad/s
    variables["P_core"] = 1.0;                      // Unitless
    variables["E_react"] = 1e46;                    // J
    variables["pi"] = 3.141592653589793;
    variables["t"] = 0.0;                           // s
}

// Update variable
void SurfaceMagneticFieldModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void SurfaceMagneticFieldModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void SurfaceMagneticFieldModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute B_s min (T)
double SurfaceMagneticFieldModule::computeB_s_min() {
    return variables["B_s_min"];
}

// Compute B_s max (T)
double SurfaceMagneticFieldModule::computeB_s_max() {
    return variables["B_s_max"];
}

// Compute B_j scaled by B_s / B_ref
double SurfaceMagneticFieldModule::computeB_j(double t, double B_s) {
    variables["t"] = t;
    double base_b = variables["B_ref"] + 0.4 * std::sin(variables["omega_s"] * t);  // Hypothetical cycle
    return base_b * (B_s / variables["B_ref"]);
}

// U_g3 example with B_j
double SurfaceMagneticFieldModule::computeU_g3_example(double t, double B_s) {
    double k_3 = variables["k_3"];
    double b_j = computeB_j(t, B_s);
    double cos_term = std::cos(variables["omega_s"] * t * variables["pi"]);
    double p_core = variables["P_core"];
    double e_react = variables["E_react"];
    return k_3 * b_j * cos_term * p_core * e_react;
}

// Equation text
std::string SurfaceMagneticFieldModule::getEquationText() {
    return "B_j ? (10^3 + 0.4 sin(?_s t)) * (B_s / 0.4) T (hypothetical scaling);\n"
           "U_g3 = k_3 * ? B_j * cos(?_s t ?) * P_core * E_react\n"
           "Where B_s = [1e-4, 0.4] T (Sun surface; quiet to sunspot).\n"
           "B_ref=0.4 T (max); scales string fields by surface B_s.\n"
           "Example t=0, B_s=0.4 T: B_j?1e3 T, U_g3?1.8e49 J/m�;\n"
           "B_s=1e-4 T: B_j?0.25 T, U_g3?4.5e45 J/m� (-4 orders).\n"
           "Role: Baseline magnetic strength for strings; variability in U_g3/disks.\n"
           "UQFF: Surface fields drive cosmic magnetism; extensible for planets.";
}

// Print variables
void SurfaceMagneticFieldModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// ===== ENHANCED METHODS IMPLEMENTATION =====

namespace surface_magnetic_field_saved_states {
    static std::map<std::string, std::map<std::string, double>> saved_states;
}

// Variable Management
void SurfaceMagneticFieldModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void SurfaceMagneticFieldModule::removeVariable(const std::string& name) {
    variables.erase(name);
}

void SurfaceMagneticFieldModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
    }
}

std::vector<std::string> SurfaceMagneticFieldModule::listVariables() const {
    std::vector<std::string> names;
    for (const auto& pair : variables) {
        names.push_back(pair.first);
    }
    return names;
}

std::string SurfaceMagneticFieldModule::getSystemName() const {
    return "Surface_Magnetic_Field_Bs_UQFF";
}

// Batch Operations
void SurfaceMagneticFieldModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = func(variables[name]);
        }
    }
}

void SurfaceMagneticFieldModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// Self-Expansion: Domain-Specific Scales
void SurfaceMagneticFieldModule::expandParameterSpace(double field_scale, double coupling_scale, double energy_scale) {
    // Scale magnetic fields
    variables["B_s_min"] *= field_scale;
    variables["B_s_max"] *= field_scale;
    variables["B_ref"] *= field_scale;
    
    // Scale coupling
    variables["k_3"] *= coupling_scale;
    
    // Scale energy
    variables["E_react"] *= energy_scale;
}

void SurfaceMagneticFieldModule::expandFieldScale(double bs_factor, double bj_factor) {
    variables["B_s_min"] *= bs_factor;
    variables["B_s_max"] *= bs_factor;
    variables["B_ref"] *= bs_factor;
    
    // Advanced field characteristics
    if (variables.find("B_s_range_T") == variables.end()) {
        // Range: max - min
        variables["B_s_range_T"] = variables["B_s_max"] - variables["B_s_min"];
    }
    variables["B_s_range_T"] *= bs_factor;
    
    // Field variability
    if (variables.find("field_variability") == variables.end()) {
        // Ratio: max/min
        variables["field_variability"] = variables["B_s_max"] / variables["B_s_min"];
    }
    // Variability stays constant with uniform scaling
    
    // Magnetic energy density
    if (variables.find("magnetic_energy_density_Jm3") == variables.end()) {
        // B²/(2μ₀), μ₀ = 4π×10⁻⁷
        double mu_0 = 4.0 * variables["pi"] * 1e-7;
        double B_avg = (variables["B_s_min"] + variables["B_s_max"]) / 2.0;
        variables["magnetic_energy_density_Jm3"] = (B_avg * B_avg) / (2.0 * mu_0);
    }
    variables["magnetic_energy_density_Jm3"] *= bs_factor * bs_factor;
    
    // Field strength scale
    if (variables.find("field_strength_scale_T") == variables.end()) {
        // Typical field strength
        variables["field_strength_scale_T"] = std::sqrt(variables["B_s_min"] * variables["B_s_max"]);
    }
    variables["field_strength_scale_T"] *= bs_factor;
    
    // String field B_j scaling
    if (variables.find("B_j_typical_T") == variables.end()) {
        // Typical B_j at reference
        variables["B_j_typical_T"] = variables["B_ref"];
    }
    variables["B_j_typical_T"] *= bj_factor;
}

void SurfaceMagneticFieldModule::expandCouplingScale(double k3_factor, double interaction_factor) {
    variables["k_3"] *= k3_factor;
    
    // Coupling characteristics
    if (variables.find("coupling_strength") == variables.end()) {
        variables["coupling_strength"] = variables["k_3"];
    }
    variables["coupling_strength"] *= k3_factor;
    
    // Interaction scale (how B_j affects U_g3)
    if (variables.find("interaction_scale") == variables.end()) {
        // k_3 * B_ref * E_react
        variables["interaction_scale"] = variables["k_3"] * variables["B_ref"] * variables["E_react"];
    }
    variables["interaction_scale"] *= k3_factor * interaction_factor;
    
    // Effective coupling
    if (variables.find("effective_coupling") == variables.end()) {
        variables["effective_coupling"] = variables["k_3"];
    }
    variables["effective_coupling"] *= k3_factor * interaction_factor;
}

void SurfaceMagneticFieldModule::expandEnergyScale(double ereact_factor, double ug3_factor) {
    variables["E_react"] *= ereact_factor;
    
    // Energy characteristics
    if (variables.find("reactor_energy_J") == variables.end()) {
        variables["reactor_energy_J"] = variables["E_react"];
    }
    variables["reactor_energy_J"] *= ereact_factor;
    
    // U_g3 scale
    if (variables.find("U_g3_scale_Jm3") == variables.end()) {
        // Typical U_g3 magnitude
        double t_test = 0.0;
        double B_s_test = variables["B_s_max"];
        variables["U_g3_scale_Jm3"] = computeU_g3_example(t_test, B_s_test);
    }
    variables["U_g3_scale_Jm3"] *= ug3_factor;
    
    // Power output (hypothetical)
    if (variables.find("power_output_W") == variables.end()) {
        // E_react / cycle_time
        double cycle_time = 2.0 * variables["pi"] / variables["omega_s"];
        variables["power_output_W"] = variables["E_react"] / cycle_time;
    }
    variables["power_output_W"] *= ereact_factor;
}

// Self-Refinement
void SurfaceMagneticFieldModule::autoRefineParameters(const std::string& target, double goal) {
    if (target == "B_s_min") {
        variables["B_s_min"] = goal;
        // Update range
        if (variables.find("B_s_range_T") != variables.end()) {
            variables["B_s_range_T"] = variables["B_s_max"] - variables["B_s_min"];
        }
        // Update variability
        if (variables.find("field_variability") != variables.end()) {
            variables["field_variability"] = variables["B_s_max"] / variables["B_s_min"];
        }
    } else if (target == "B_s_max") {
        variables["B_s_max"] = goal;
        if (variables.find("B_s_range_T") != variables.end()) {
            variables["B_s_range_T"] = variables["B_s_max"] - variables["B_s_min"];
        }
        if (variables.find("field_variability") != variables.end()) {
            variables["field_variability"] = variables["B_s_max"] / variables["B_s_min"];
        }
    } else if (target == "B_ref") {
        variables["B_ref"] = goal;
    } else if (target == "k_3") {
        variables["k_3"] = goal;
    } else if (target == "E_react") {
        variables["E_react"] = goal;
    } else if (target == "omega_s") {
        variables["omega_s"] = goal;
    } else if (target == "U_g3_at_max") {
        // Target specific U_g3 at B_s_max, t=0
        double current_U_g3 = computeU_g3_example(0.0, variables["B_s_max"]);
        if (current_U_g3 > 0) {
            double scale = goal / current_U_g3;
            variables["k_3"] *= scale;
        }
    } else if (target == "magnetic_energy_density_Jm3") {
        if (variables.find("magnetic_energy_density_Jm3") == variables.end()) {
            double mu_0 = 4.0 * variables["pi"] * 1e-7;
            double B_avg = (variables["B_s_min"] + variables["B_s_max"]) / 2.0;
            variables["magnetic_energy_density_Jm3"] = (B_avg * B_avg) / (2.0 * mu_0);
        }
        variables["magnetic_energy_density_Jm3"] = goal;
        // Back-calculate B_avg
        double mu_0 = 4.0 * variables["pi"] * 1e-7;
        double B_avg = std::sqrt(goal * 2.0 * mu_0);
        // Scale B_s to match
        double current_avg = (variables["B_s_min"] + variables["B_s_max"]) / 2.0;
        if (current_avg > 0) {
            double scale = B_avg / current_avg;
            variables["B_s_min"] *= scale;
            variables["B_s_max"] *= scale;
        }
    }
}

void SurfaceMagneticFieldModule::calibrateToObservations(const std::map<std::string, double>& observations) {
    for (const auto& obs : observations) {
        if (variables.find(obs.first) != variables.end()) {
            variables[obs.first] = obs.second;
        }
    }
    // Update derived quantities
    if (variables.find("B_s_range_T") != variables.end()) {
        variables["B_s_range_T"] = variables["B_s_max"] - variables["B_s_min"];
    }
    if (variables.find("field_variability") != variables.end() && variables["B_s_min"] > 0) {
        variables["field_variability"] = variables["B_s_max"] / variables["B_s_min"];
    }
}

void SurfaceMagneticFieldModule::optimizeForMetric(const std::string& metric) {
    if (metric == "sun_quiet") {
        // Quiet Sun (solar minimum)
        variables["B_s_min"] = 1e-4;
        variables["B_s_max"] = 0.01;  // Low activity
        variables["B_ref"] = 0.4;
    } else if (metric == "sun_active") {
        // Active Sun (solar maximum)
        variables["B_s_min"] = 1e-3;
        variables["B_s_max"] = 0.4;  // Sunspot max
        variables["B_ref"] = 0.4;
    } else if (metric == "sun_standard") {
        // Standard Sun
        variables["B_s_min"] = 1e-4;
        variables["B_s_max"] = 0.4;
        variables["B_ref"] = 0.4;
    } else if (metric == "weak_field") {
        // Weak magnetic field (old star, low activity)
        variables["B_s_min"] = 1e-5;
        variables["B_s_max"] = 0.01;
        variables["B_ref"] = 0.01;
    } else if (metric == "strong_field") {
        // Strong magnetic field (young star, high activity)
        variables["B_s_min"] = 1e-3;
        variables["B_s_max"] = 1.0;  // Stronger than Sun
        variables["B_ref"] = 1.0;
    } else if (metric == "neutron_star") {
        // Neutron star (extreme)
        variables["B_s_min"] = 1e8;  // 100 MT
        variables["B_s_max"] = 1e11; // 100 GT
        variables["B_ref"] = 1e11;
    } else if (metric == "white_dwarf") {
        // White dwarf magnetic
        variables["B_s_min"] = 1e2;   // 100 T
        variables["B_s_max"] = 1e6;   // 1 MT
        variables["B_ref"] = 1e6;
    } else if (metric == "jupiter") {
        // Jupiter-like planet
        variables["B_s_min"] = 4e-4;  // ~4 Gauss
        variables["B_s_max"] = 1e-3;  // ~10 Gauss
        variables["B_ref"] = 1e-3;
    } else if (metric == "earth") {
        // Earth-like planet
        variables["B_s_min"] = 3e-5;  // ~0.3 Gauss
        variables["B_s_max"] = 6e-5;  // ~0.6 Gauss
        variables["B_ref"] = 6e-5;
    } else if (metric == "magnetar") {
        // Magnetar (most extreme known)
        variables["B_s_min"] = 1e10;  // 10 GT
        variables["B_s_max"] = 1e11;  // 100 GT
        variables["B_ref"] = 1e11;
    }
}

// Parameter Exploration
std::vector<std::map<std::string, double>> SurfaceMagneticFieldModule::generateVariations(int count, double variation_pct) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - variation_pct, 1.0 + variation_pct);
    
    for (int i = 0; i < count; ++i) {
        std::map<std::string, double> variant = variables;
        for (auto& pair : variant) {
            if (pair.first != "pi") {  // Don't vary π
                pair.second *= dis(gen);
            }
        }
        variations.push_back(variant);
    }
    return variations;
}

// Adaptive Evolution
void SurfaceMagneticFieldModule::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - mutation_rate, 1.0 + mutation_rate);
    
    for (auto& pair : variables) {
        if (pair.first != "pi") {  // Don't mutate π
            pair.second *= dis(gen);
        }
    }
}

void SurfaceMagneticFieldModule::evolveSystem(int generations, std::function<double()> fitness_func) {
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
void SurfaceMagneticFieldModule::saveState(const std::string& label) {
    surface_magnetic_field_saved_states::saved_states[label] = variables;
}

void SurfaceMagneticFieldModule::restoreState(const std::string& label) {
    if (surface_magnetic_field_saved_states::saved_states.find(label) != surface_magnetic_field_saved_states::saved_states.end()) {
        variables = surface_magnetic_field_saved_states::saved_states[label];
    }
}

std::vector<std::string> SurfaceMagneticFieldModule::listSavedStates() const {
    std::vector<std::string> labels;
    for (const auto& pair : surface_magnetic_field_saved_states::saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string SurfaceMagneticFieldModule::exportState() const {
    std::ostringstream oss;
    oss << "SurfaceMagneticField_State_Export:\n";
    for (const auto& pair : variables) {
        oss << pair.first << "=" << std::scientific << pair.second << "\n";
    }
    return oss.str();
}

// System Analysis
std::map<std::string, double> SurfaceMagneticFieldModule::sensitivityAnalysis(const std::vector<std::string>& params) {
    std::map<std::string, double> sensitivities;
    double t_test = 0.0;
    double B_s_test = variables["B_s_max"];
    double baseline = computeU_g3_example(t_test, B_s_test);
    
    for (const auto& param : params) {
        if (variables.find(param) != variables.end() && param != "pi") {
            double original = variables[param];
            variables[param] = original * 1.01;
            
            double perturbed = computeU_g3_example(t_test, B_s_test);
            if (baseline > 0) {
                sensitivities[param] = (perturbed - baseline) / baseline;
            } else {
                sensitivities[param] = 0.0;
            }
            
            // Restore
            variables[param] = original;
        }
    }
    return sensitivities;
}

std::string SurfaceMagneticFieldModule::generateReport() const {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(3);
    oss << "===== Surface Magnetic Field (B_s) Module Report =====\n";
    oss << "System: " << getSystemName() << "\n\n";
    
    oss << "Surface Magnetic Field Range:\n";
    oss << "  B_s_min = " << variables.at("B_s_min") << " T (quiet/minimum)\n";
    oss << "  B_s_max = " << variables.at("B_s_max") << " T (active/maximum)\n";
    oss << "  B_ref = " << variables.at("B_ref") << " T (reference)\n";
    
    double range = variables.at("B_s_max") - variables.at("B_s_min");
    oss << "  Range = " << range << " T\n";
    
    if (variables.at("B_s_min") > 0) {
        double variability = variables.at("B_s_max") / variables.at("B_s_min");
        oss << "  Variability (max/min) = " << std::fixed << std::setprecision(1) 
            << variability << "x\n";
    }
    oss << "\n";
    
    oss << std::scientific;
    if (variables.find("B_s_range_T") != variables.end()) {
        oss << "Field Characteristics:\n";
        oss << "  B_s range = " << variables.at("B_s_range_T") << " T\n";
        if (variables.find("field_variability") != variables.end()) {
            oss << "  Field variability = " << std::fixed << std::setprecision(1) 
                << variables.at("field_variability") << "x\n";
        }
        oss << std::scientific;
        if (variables.find("magnetic_energy_density_Jm3") != variables.end()) {
            oss << "  Magnetic energy density = " << variables.at("magnetic_energy_density_Jm3") << " J/m³\n";
        }
        if (variables.find("field_strength_scale_T") != variables.end()) {
            oss << "  Field strength scale = " << variables.at("field_strength_scale_T") << " T\n";
        }
        oss << "\n";
    }
    
    oss << "Coupling & Energy:\n";
    oss << "  k_3 = " << variables.at("k_3") << " (coupling constant)\n";
    oss << "  E_react = " << variables.at("E_react") << " J (reactor energy)\n";
    oss << "  ω_s = " << variables.at("omega_s") << " rad/s (rotation)\n";
    oss << "  P_core = " << variables.at("P_core") << " (core power)\n";
    
    double cycle_period = 2.0 * variables.at("pi") / variables.at("omega_s");
    oss << "  Rotation period = " << cycle_period << " s (";
    oss << std::fixed << std::setprecision(1) << (cycle_period / 86400.0) << " days)\n";
    oss << "\n";
    
    oss << std::scientific;
    oss << "Magnetic String Field B_j:\n";
    double t_test = 0.0;
    std::vector<double> B_s_tests = {variables.at("B_s_min"), 
                                       (variables.at("B_s_min") + variables.at("B_s_max")) / 2.0,
                                       variables.at("B_s_max")};
    std::vector<std::string> labels = {"min", "avg", "max"};
    
    for (size_t i = 0; i < B_s_tests.size(); ++i) {
        double B_j = const_cast<SurfaceMagneticFieldModule*>(this)->computeB_j(t_test, B_s_tests[i]);
        oss << "  B_j at B_s=" << labels[i] << " (" << B_s_tests[i] << " T): " << B_j << " T\n";
    }
    oss << "\n";
    
    oss << "U_g3 Energy Density (t=0):\n";
    for (size_t i = 0; i < B_s_tests.size(); ++i) {
        double U_g3 = const_cast<SurfaceMagneticFieldModule*>(this)->computeU_g3_example(t_test, B_s_tests[i]);
        oss << "  U_g3 at B_s=" << labels[i] << ": " << U_g3 << " J/m³\n";
    }
    oss << "\n";
    
    oss << "Physical Interpretation:\n";
    double B_max = variables.at("B_s_max");
    if (B_max > 1e9) {
        oss << "  Magnetar regime (B_s > 1 GT, extreme magnetic field)\n";
    } else if (B_max > 1e6) {
        oss << "  White dwarf magnetic regime (B_s > 1 MT)\n";
    } else if (B_max > 1e3) {
        oss << "  Neutron star regime (B_s > 1 kT)\n";
    } else if (B_max > 0.1) {
        oss << "  Strong stellar field (B_s > 0.1 T, active star)\n";
    } else if (B_max > 1e-3) {
        oss << "  Moderate stellar field (B_s 1-1000 mT, solar-type)\n";
    } else if (B_max > 1e-5) {
        oss << "  Planetary field (B_s ~10-1000 μT, Earth/Jupiter-like)\n";
    } else {
        oss << "  Weak field (B_s < 10 μT, low activity)\n";
    }
    
    double B_min = variables.at("B_s_min");
    if (B_min < 1e-6) {
        oss << "  Very quiet baseline (B_s_min < 1 μT)\n";
    } else if (B_min < 1e-4) {
        oss << "  Quiet baseline (B_s_min < 100 μT)\n";
    } else if (B_min < 1e-2) {
        oss << "  Moderate baseline (B_s_min 100 μT - 10 mT)\n";
    } else {
        oss << "  High baseline (B_s_min > 10 mT, always active)\n";
    }
    
    oss << "\n  Applications:\n";
    oss << "    - Solar cycle modeling: B_s variability drives 11-year cycle\n";
    oss << "    - Magnetic string dynamics: B_j scales with surface B_s\n";
    oss << "    - U_g3 field strength: Orbital stability influenced by B_s\n";
    oss << "    - Stellar activity: Sunspot/starspot coverage ∝ B_s\n";
    oss << "    - Planetary magnetism: Dipole fields in UQFF framework\n";
    oss << "    - Compact objects: Neutron stars, white dwarfs, magnetars\n";
    oss << "    - Magnetic reconnection: Flare/CME energy from B_s gradients\n";
    
    return oss.str();
}

bool SurfaceMagneticFieldModule::validateConsistency() const {
    bool valid = true;
    
    // Check B_s_min < B_s_max
    if (variables.at("B_s_min") >= variables.at("B_s_max")) {
        std::cerr << "Error: B_s_min must be less than B_s_max\n";
        valid = false;
    }
    
    // Check positive fields
    if (variables.at("B_s_min") <= 0 || variables.at("B_s_max") <= 0 || variables.at("B_ref") <= 0) {
        std::cerr << "Error: Magnetic field strengths must be positive\n";
        valid = false;
    }
    
    // Check k_3 is positive
    if (variables.at("k_3") <= 0) {
        std::cerr << "Error: Coupling constant k_3 must be positive\n";
        valid = false;
    }
    
    // Check E_react is positive
    if (variables.at("E_react") <= 0) {
        std::cerr << "Error: Reactor energy E_react must be positive\n";
        valid = false;
    }
    
    // Check omega_s is positive
    if (variables.at("omega_s") <= 0) {
        std::cerr << "Error: Rotation rate omega_s must be positive\n";
        valid = false;
    }
    
    // Check reasonable B_s range
    if (variables.at("B_s_max") / variables.at("B_s_min") > 1e10) {
        std::cerr << "Warning: Very large B_s variability (max/min > 1e10)\n";
    }
    
    return valid;
}

void SurfaceMagneticFieldModule::autoCorrectAnomalies() {
    // Ensure B_s_min < B_s_max
    if (variables["B_s_min"] >= variables["B_s_max"]) {
        variables["B_s_min"] = 1e-4;
        variables["B_s_max"] = 0.4;
    }
    
    // Ensure positive fields
    if (variables["B_s_min"] <= 0) {
        variables["B_s_min"] = 1e-4;
    }
    if (variables["B_s_max"] <= 0) {
        variables["B_s_max"] = 0.4;
    }
    if (variables["B_ref"] <= 0) {
        variables["B_ref"] = 0.4;
    }
    
    // Ensure k_3 is positive and reasonable
    if (variables["k_3"] <= 0 || variables["k_3"] > 100.0) {
        variables["k_3"] = 1.8;
    }
    
    // Ensure E_react is positive
    if (variables["E_react"] <= 0) {
        variables["E_react"] = 1e46;
    }
    
    // Ensure omega_s is positive
    if (variables["omega_s"] <= 0) {
        variables["omega_s"] = 2.5e-6;
    }
    
    // Ensure P_core is reasonable
    if (variables["P_core"] <= 0 || variables["P_core"] > 10.0) {
        variables["P_core"] = 1.0;
    }
}

// Example usage in base program (snippet)
// #include "SurfaceMagneticFieldModule.h"
// int main() {
//     SurfaceMagneticFieldModule mod;
//     double b_min = mod.computeB_s_min();
//     std::cout << "B_s range: " << b_min << " to " << mod.computeB_s_max() << " T\n";
//     double u_g3 = mod.computeU_g3_example(0.0, 1e-4);
//     std::cout << "U_g3 (quiet Sun) = " << u_g3 << " J/m�\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("B_s_min", 5e-5);
//     mod.printVariables();
//     return 0;
// }

// ========== COMPREHENSIVE ENHANCED DEMONSTRATION ==========
/*
int main() {
    std::cout << std::scientific << std::setprecision(3);
    std::cout << "========== SURFACE MAGNETIC FIELD MODULE DEMONSTRATION ==========\n\n";
    
    // ===== Step 1: Initialize Module =====
    SurfaceMagneticFieldModule mod;
    std::cout << "Step 1: Module initialized with defaults (Sun):\n";
    std::cout << "  System: " << mod.getSystemName() << "\n";
    std::cout << "  B_s range: [" << mod.computeB_s_min() << ", " << mod.computeB_s_max() << "] T\n";
    std::cout << "  B_ref = " << mod.variables["B_ref"] << " T\n\n";
    
    // ===== Step 2: Compute Baseline Fields =====
    std::cout << "Step 2: Compute baseline B_j and U_g3:\n";
    double t_test = 0.0;
    
    double B_j_quiet = mod.computeB_j(t_test, mod.variables["B_s_min"]);
    double U_g3_quiet = mod.computeU_g3_example(t_test, mod.variables["B_s_min"]);
    std::cout << "  Quiet Sun (B_s=" << mod.variables["B_s_min"] << " T):\n";
    std::cout << "    B_j = " << B_j_quiet << " T\n";
    std::cout << "    U_g3 = " << U_g3_quiet << " J/m³\n";
    
    double B_j_active = mod.computeB_j(t_test, mod.variables["B_s_max"]);
    double U_g3_active = mod.computeU_g3_example(t_test, mod.variables["B_s_max"]);
    std::cout << "  Active Sun (B_s=" << mod.variables["B_s_max"] << " T):\n";
    std::cout << "    B_j = " << B_j_active << " T\n";
    std::cout << "    U_g3 = " << U_g3_active << " J/m³\n";
    std::cout << "  Amplification: " << std::fixed << std::setprecision(1) 
              << (U_g3_active / U_g3_quiet) << "x\n\n";
    
    // ===== Step 3: Variable Management =====
    std::cout << std::scientific;
    std::cout << "Step 3: Variable Management\n";
    mod.createVariable("B_s_range_T", mod.variables["B_s_max"] - mod.variables["B_s_min"]);
    std::cout << "  Created 'B_s_range_T' = " << mod.variables["B_s_range_T"] << " T\n";
    
    std::vector<std::string> all_vars = mod.listVariables();
    std::cout << "  Total variables: " << all_vars.size() << "\n";
    
    mod.cloneVariable("B_s_max", "B_s_max_backup");
    std::cout << "  Cloned B_s_max to backup\n\n";
    
    // ===== Step 4: Batch Operations =====
    std::cout << "Step 4: Batch Operations (scale field strengths)\n";
    std::vector<std::string> field_params = {"B_s_min", "B_s_max", "B_ref"};
    mod.scaleVariableGroup(field_params, 2.0);  // Double fields
    std::cout << "  Scaled field parameters by 2.0x:\n";
    std::cout << "    B_s_min = " << mod.variables["B_s_min"] << " T\n";
    std::cout << "    B_s_max = " << mod.variables["B_s_max"] << " T\n";
    std::cout << "    B_ref = " << mod.variables["B_ref"] << " T\n";
    
    double U_g3_scaled = mod.computeU_g3_example(t_test, mod.variables["B_s_max"]);
    std::cout << "    U_g3 (scaled) = " << U_g3_scaled << " J/m³\n\n";
    
    // Restore
    mod.scaleVariableGroup(field_params, 0.5);
    
    // ===== Step 5: Self-Expansion - Field Scale =====
    std::cout << "Step 5: Self-Expansion - Field Scale\n";
    mod.saveState("before_field_expansion");
    std::cout << "  Initial B_s range: [" << mod.variables["B_s_min"] << ", " 
              << mod.variables["B_s_max"] << "] T\n";
    
    mod.expandFieldScale(3.0, 2.0);  // 3x B_s, 2x B_j
    std::cout << "  After expandFieldScale(3.0, 2.0):\n";
    std::cout << "    B_s_min = " << mod.variables["B_s_min"] << " T\n";
    std::cout << "    B_s_max = " << mod.variables["B_s_max"] << " T\n";
    if (mod.variables.find("magnetic_energy_density_Jm3") != mod.variables.end()) {
        std::cout << "    Magnetic energy density = " << mod.variables["magnetic_energy_density_Jm3"] << " J/m³\n";
    }
    std::cout << "\n";
    
    mod.restoreState("before_field_expansion");
    
    // ===== Step 6-25 continued... =====
    std::cout << "========== DEMONSTRATION COMPLETE ==========\n";
    
    return 0;
}
*/

// Compile: g++ -o surface_b_test surface_b_test.cpp SurfaceMagneticFieldModule.cpp -lm
// Sample: B_s [1e-4, 0.4] T; U_g3 quiet?4.5e45 J/m�; scales magnetic influence.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

SurfaceMagneticFieldModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeB_s_min, computeB_s_max, computeB_j, computeU_g3_example) are clear, concise, and variable - driven.
- Output and debugging functions(printVariables, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.
- Supports a wide range of surface magnetic field strengths, enabling both quiet and active Sun scenarios.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in surface magnetic field modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.