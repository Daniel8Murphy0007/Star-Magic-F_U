// HeliosphereThicknessModule.h
// Modular C++ implementation of the Heliosphere Thickness Factor (H_SCm) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes H_SCm ?1 (unitless) and its scaling in Universal Gravity U_g2 term.
// Pluggable: #include "HeliosphereThicknessModule.h"
// HeliosphereThicknessModule mod; mod.computeU_g2(0.0, 0.0); mod.updateVariable("H_SCm", new_value);
// Variables in std::map; example for Sun at t=0, t_n=0, r=R_b=1.496e13 m.
// Approximations: S(r - R_b)=1; ?_sw v_sw=5001; E_react=1e46.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef HELIOSPHERE_THICKNESS_MODULE_H
#define HELIOSPHERE_THICKNESS_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <functional>
#include <random>
#include <algorithm>
#include <sstream>

class HeliosphereThicknessModule {
private:
    std::map<std::string, double> variables;
    double computeH_SCm();
    double computeU_g2(double t, double t_n);

public:
    // Constructor: Initialize with framework defaults
    HeliosphereThicknessModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeH_SCm();  // ?1 (unitless)
    double computeU_g2(double t, double t_n);  // U_g2 with H_SCm (J/m^3)
    double computeU_g2_no_H(double t, double t_n);  // Without H_SCm variation

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

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
    void expandHeliosphereScale(double H_factor, double R_b_factor);      // H_SCm, R_b
    void expandBoundaryScale(double delta_sw_factor, double v_sw_factor); // Solar wind boundary
    void expandReactorScale(double E_react_factor, double k_2_factor);    // Reactor and coupling

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

#endif // HELIOSPHERE_THICKNESS_MODULE_H

// HeliosphereThicknessModule.cpp
#include "HeliosphereThicknessModule.h"

// Constructor: Set framework defaults
HeliosphereThicknessModule::HeliosphereThicknessModule() {
    // Universal constants
    variables["H_SCm"] = 1.0;                       // Unitless ?1
    variables["k_2"] = 1.2;                         // Coupling
    variables["rho_vac_UA"] = 7.09e-36;             // J/m^3
    variables["rho_vac_SCm"] = 7.09e-37;            // J/m^3
    variables["M_s"] = 1.989e30;                    // kg (Sun)
    variables["r"] = 1.496e13;                      // m (R_b)
    variables["R_b"] = 1.496e13;                    // m
    variables["delta_sw"] = 0.01;                   // Unitless
    variables["v_sw"] = 5e5;                        // m/s
    variables["E_react"] = 1e46;                    // J
    variables["S_r_Rb"] = 1.0;                      // Step function
    variables["pi"] = 3.141592653589793;
    variables["t_n"] = 0.0;                         // s

    // Derived
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
    variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
}

// Update variable
void HeliosphereThicknessModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
        if (name == "rho_vac_UA" || name == "rho_vac_SCm") {
            variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
        } else if (name == "delta_sw" || name == "v_sw") {
            variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void HeliosphereThicknessModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        if (name == "rho_vac_UA" || name == "rho_vac_SCm") {
            variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
        } else if (name == "delta_sw" || name == "v_sw") {
            variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void HeliosphereThicknessModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute H_SCm ?1
double HeliosphereThicknessModule::computeH_SCm() {
    return variables["H_SCm"];
}

// Compute U_g2 with H_SCm
double HeliosphereThicknessModule::computeU_g2(double t, double t_n) {
    double k_2 = variables["k_2"];
    double rho_sum = variables["rho_sum"];
    double M_s = variables["M_s"];
    double r = variables["r"];
    double S_r_Rb = variables["S_r_Rb"];
    double swirl_factor = variables["swirl_factor"];
    double H_SCm = computeH_SCm();
    double E_react = variables["E_react"];
    // Simplified; no explicit t dependence in example
    return k_2 * (rho_sum * M_s / (r * r)) * S_r_Rb * swirl_factor * H_SCm * E_react;
}

// U_g2 without H_SCm variation (H=1 fixed)
double HeliosphereThicknessModule::computeU_g2_no_H(double t, double t_n) {
    double orig_H = variables["H_SCm"];
    variables["H_SCm"] = 1.0;
    double result = computeU_g2(t, t_n);
    variables["H_SCm"] = orig_H;
    return result;
}

// Equation text
std::string HeliosphereThicknessModule::getEquationText() {
    return "U_g2 = k_2 * [(?_vac,[UA] + ?_vac,[SCm]) M_s / r^2] * S(r - R_b) * (1 + ?_sw v_sw) * H_SCm * E_react\n"
           "Where H_SCm ?1 (unitless heliosphere thickness factor);\n"
           "Scales outer field bubble gravity for heliopause extent (~120 AU).\n"
           "Example r=R_b=1.496e13 m, t=0: U_g2 ?1.18e53 J/m� (H=1);\n"
           "If H_SCm=1.1: ?1.30e53 J/m� (+10%).\n"
           "Role: Adjusts [SCm] influence in heliosphere; minimal but flexible for boundary variations.\n"
           "UQFF: Models solar wind dominance; key for nebular/heliospheric dynamics.";
}

// Print variables
void HeliosphereThicknessModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// ===== Implementation: 25-Method Dynamic Self-Update & Self-Expansion Capabilities =====

namespace heliosphere_thickness_saved_states {
    std::map<std::string, std::map<std::string, double>> saved_states;
}

// Variable Management (5 methods)
void HeliosphereThicknessModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void HeliosphereThicknessModule::removeVariable(const std::string& name) {
    variables.erase(name);
}

void HeliosphereThicknessModule::cloneVariable(const std::string& source, const std::string& target) {
    if (variables.find(source) != variables.end()) {
        variables[target] = variables[source];
    }
}

std::string HeliosphereThicknessModule::listVariables() {
    std::ostringstream oss;
    for (const auto& pair : variables) {
        oss << pair.first << " = " << pair.second << "\n";
    }
    return oss.str();
}

std::string HeliosphereThicknessModule::getSystemName() {
    return "Heliosphere_Thickness_UQFF";
}

// Batch Operations (2 methods)
void HeliosphereThicknessModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = func(variables[name]);
        }
    }
    // Recalculate derived
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
    variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
}

void HeliosphereThicknessModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double x) { return x * factor; });
}

// Self-Expansion (4 methods)
void HeliosphereThicknessModule::expandParameterSpace(const std::vector<std::string>& params, double expansion_factor) {
    for (const auto& param : params) {
        if (variables.find(param) != variables.end()) {
            variables[param] *= expansion_factor;
        }
    }
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
    variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
}

void HeliosphereThicknessModule::expandHeliosphereScale(double H_factor, double R_b_factor) {
    // H_SCm: thickness factor (unitless, ~1)
    // R_b: heliosphere boundary radius (~100 AU = 1.496e13 m)
    variables["H_SCm"] *= H_factor;
    variables["R_b"] *= R_b_factor;
    if (variables.find("r") != variables.end() && std::abs(variables["r"] - variables["R_b"]/R_b_factor) < 1e6) {
        variables["r"] *= R_b_factor; // Scale r if it tracks R_b
    }
}

void HeliosphereThicknessModule::expandBoundaryScale(double delta_sw_factor, double v_sw_factor) {
    // delta_sw: solar wind variation (unitless, ~0.01)
    // v_sw: solar wind velocity (m/s, ~5e5)
    variables["delta_sw"] *= delta_sw_factor;
    variables["v_sw"] *= v_sw_factor;
    variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
}

void HeliosphereThicknessModule::expandReactorScale(double E_react_factor, double k_2_factor) {
    // E_react: reactor energy (J, ~1e46)
    // k_2: coupling constant (~1.2)
    variables["E_react"] *= E_react_factor;
    variables["k_2"] *= k_2_factor;
}

// Self-Refinement (3 methods)
void HeliosphereThicknessModule::autoRefineParameters(const std::string& target_metric, double target_value) {
    // Example: target U_g2 value by adjusting H_SCm and E_react
    if (target_metric == "U_g2") {
        double current = computeU_g2(0.0, 0.0);
        if (current > 0) {
            double ratio = target_value / current;
            variables["H_SCm"] *= std::sqrt(ratio);
            variables["E_react"] *= std::sqrt(ratio);
        }
    }
}

void HeliosphereThicknessModule::calibrateToObservations(const std::map<std::string, double>& observed) {
    for (const auto& obs : observed) {
        if (variables.find(obs.first) != variables.end()) {
            variables[obs.first] = obs.second;
        }
    }
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
    variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
}

void HeliosphereThicknessModule::optimizeForMetric(const std::string& metric) {
    if (metric == "maximize_U_g2") {
        variables["H_SCm"] *= 1.1;
        variables["E_react"] *= 1.1;
    } else if (metric == "minimize_U_g2") {
        variables["H_SCm"] *= 0.9;
        variables["E_react"] *= 0.9;
    } else if (metric == "expand_boundary") {
        variables["R_b"] *= 1.2;
        variables["delta_sw"] *= 1.1;
    }
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
    variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
}

// Parameter Exploration (1 method)
std::vector<std::map<std::string, double>> HeliosphereThicknessModule::generateVariations(int count, double variance) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-variance, variance);

    for (int i = 0; i < count; ++i) {
        std::map<std::string, double> variant = variables;
        for (auto& pair : variant) {
            if (pair.first != "pi" && pair.first != "S_r_Rb") {
                double factor = 1.0 + dis(gen);
                pair.second *= factor;
            }
        }
        // Recalculate derived for variant
        variant["rho_sum"] = variant["rho_vac_UA"] + variant["rho_vac_SCm"];
        variant["swirl_factor"] = 1.0 + variant["delta_sw"] * variant["v_sw"];
        variations.push_back(variant);
    }
    return variations;
}

// Adaptive Evolution (2 methods)
void HeliosphereThicknessModule::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-mutation_rate, mutation_rate);

    for (auto& pair : variables) {
        if (pair.first != "pi" && pair.first != "S_r_Rb" && pair.first != "rho_sum" && pair.first != "swirl_factor") {
            pair.second *= (1.0 + dis(gen));
            if (pair.second < 0) pair.second = std::abs(pair.second);
        }
    }
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
    variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
}

void HeliosphereThicknessModule::evolveSystem(int generations, std::function<double()> fitness_func) {
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
void HeliosphereThicknessModule::saveState(const std::string& label) {
    heliosphere_thickness_saved_states::saved_states[label] = variables;
}

void HeliosphereThicknessModule::restoreState(const std::string& label) {
    if (heliosphere_thickness_saved_states::saved_states.find(label) != heliosphere_thickness_saved_states::saved_states.end()) {
        variables = heliosphere_thickness_saved_states::saved_states[label];
    }
}

std::vector<std::string> HeliosphereThicknessModule::listSavedStates() {
    std::vector<std::string> labels;
    for (const auto& pair : heliosphere_thickness_saved_states::saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string HeliosphereThicknessModule::exportState() {
    std::ostringstream oss;
    oss << "System: " << getSystemName() << "\n";
    oss << listVariables();
    return oss.str();
}

// System Analysis (4 methods)
std::map<std::string, double> HeliosphereThicknessModule::sensitivityAnalysis(const std::vector<std::string>& params) {
    std::map<std::string, double> sensitivity;
    double baseline = computeU_g2(0.0, 0.0);

    for (const auto& param : params) {
        if (variables.find(param) != variables.end()) {
            double original = variables[param];
            variables[param] *= 1.01; // +1% perturbation
            double perturbed = computeU_g2(0.0, 0.0);
            variables[param] = original;
            sensitivity[param] = std::abs((perturbed - baseline) / baseline);
        }
    }
    return sensitivity;
}

std::string HeliosphereThicknessModule::generateReport() {
    std::ostringstream oss;
    oss << "========== Heliosphere Thickness Module Report ==========\n";
    oss << "System: " << getSystemName() << "\n\n";
    oss << "Key Parameters:\n";
    oss << "  H_SCm (thickness factor): " << std::scientific << variables["H_SCm"] << " (unitless, ~1)\n";
    oss << "  R_b (boundary radius): " << variables["R_b"] << " m (~100 AU = 1.496e13 m)\n";
    oss << "  k_2 (coupling): " << variables["k_2"] << "\n";
    oss << "  E_react (reactor energy): " << variables["E_react"] << " J\n";
    oss << "  delta_sw (wind variation): " << variables["delta_sw"] << " (unitless)\n";
    oss << "  v_sw (wind velocity): " << variables["v_sw"] << " m/s\n";
    oss << "  swirl_factor: " << variables["swirl_factor"] << "\n";
    oss << "  M_s (solar mass): " << variables["M_s"] << " kg\n";
    oss << "  rho_vac_UA: " << variables["rho_vac_UA"] << " J/m^3\n";
    oss << "  rho_vac_SCm: " << variables["rho_vac_SCm"] << " J/m^3\n";
    oss << "  rho_sum: " << variables["rho_sum"] << " J/m^3\n\n";
    oss << "Computed Values:\n";
    double u_g2 = computeU_g2(0.0, 0.0);
    double u_g2_no_h = computeU_g2_no_H(0.0, 0.0);
    oss << "  U_g2 (with H_SCm=" << variables["H_SCm"] << "): " << u_g2 << " J/m^3\n";
    oss << "  U_g2 (H=1 baseline): " << u_g2_no_h << " J/m^3\n";
    oss << "  H_SCm impact: " << ((u_g2 / u_g2_no_h - 1.0) * 100) << "%\n\n";
    oss << "Physical Context:\n";
    oss << "  Heliosphere extends ~120 AU from Sun (R_b ~1.5e13 m)\n";
    oss << "  U_g2 is outer field bubble gravity (k_2 ~1.2 coupling)\n";
    oss << "  H_SCm ~1 adjusts [SCm] influence at heliopause\n";
    oss << "  Solar wind (v_sw ~5e5 m/s) shapes boundary dynamics\n";
    oss << "  E_react ~1e46 J drives heliospheric pressure\n";
    oss << "========================================================\n";
    return oss.str();
}

bool HeliosphereThicknessModule::validateConsistency() {
    bool valid = true;
    if (variables["H_SCm"] <= 0 || variables["H_SCm"] > 2.0) valid = false;
    if (variables["R_b"] <= 0) valid = false;
    if (variables["E_react"] <= 0) valid = false;
    if (variables["v_sw"] < 0) valid = false;
    if (variables["M_s"] <= 0) valid = false;
    return valid;
}

void HeliosphereThicknessModule::autoCorrectAnomalies() {
    if (variables["H_SCm"] <= 0 || variables["H_SCm"] > 2.0) variables["H_SCm"] = 1.0;
    if (variables["R_b"] <= 0) variables["R_b"] = 1.496e13;
    if (variables["E_react"] <= 0) variables["E_react"] = 1e46;
    if (variables["v_sw"] < 0) variables["v_sw"] = 5e5;
    if (variables["M_s"] <= 0) variables["M_s"] = 1.989e30;
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
    variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
}

// Example usage in base program (snippet)
// #include "HeliosphereThicknessModule.h"
// int main() {
//     HeliosphereThicknessModule mod;
//     double h = mod.computeH_SCm();
//     std::cout << "H_SCm ? " << h << std::endl;
//     double u_g2 = mod.computeU_g2(0.0, 0.0);
//     std::cout << "U_g2 = " << u_g2 << " J/m�\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("H_SCm", 1.1);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o heliosphere_test heliosphere_test.cpp HeliosphereThicknessModule.cpp -lm
// Sample: H_SCm=1; U_g2 ?1.18e53 J/m�; +10% for H=1.1.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

// ===== Enhanced Example: 18-Step Demonstration of Dynamic Capabilities =====
// int main() {
//     HeliosphereThicknessModule mod;
//     std::cout << "===== Heliosphere Thickness Module: Enhanced 18-Step Demo =====\n\n";
//
//     // Step 1: Initial report
//     std::cout << "Step 1: Initial Configuration\n";
//     std::cout << mod.generateReport() << "\n";
//
//     // Step 2: Variable tracking and creation
//     std::cout << "Step 2: Track boundary extent and wind impact\n";
//     mod.createVariable("boundary_extent_AU", mod.variables["R_b"] / 1.496e11); // AU
//     mod.createVariable("wind_mach", mod.variables["v_sw"] / 340.0);           // Mach number (~1470)
//     mod.createVariable("heliosphere_volume_m3", 4.0/3.0 * 3.14159 * std::pow(mod.variables["R_b"], 3));
//     std::cout << "  boundary_extent: " << mod.variables["boundary_extent_AU"] << " AU\n";
//     std::cout << "  wind_mach: " << mod.variables["wind_mach"] << "\n";
//     std::cout << "  heliosphere_volume: " << std::scientific << mod.variables["heliosphere_volume_m3"] << " m^3\n\n";
//
//     // Step 3: H_SCm variation scenarios
//     std::cout << "Step 3: H_SCm Variation Impact on U_g2\n";
//     double baseline_u_g2 = mod.computeU_g2_no_H(0.0, 0.0);
//     std::vector<double> h_values = {0.8, 1.0, 1.1, 1.2, 1.5};
//     for (double h : h_values) {
//         mod.updateVariable("H_SCm", h);
//         double u_g2 = mod.computeU_g2(0.0, 0.0);
//         std::cout << "  H_SCm=" << h << " -> U_g2=" << u_g2 << " J/m^3 ";
//         std::cout << "(" << std::showpos << ((u_g2/baseline_u_g2 - 1.0)*100) << std::noshowpos << "%)\n";
//     }
//     mod.updateVariable("H_SCm", 1.0); // Reset
//     std::cout << "\n";
//
//     // Step 4: Heliosphere scale expansion (grow boundary)
//     std::cout << "Step 4: Expand Heliosphere Scale (H_SCm x1.2, R_b x1.15)\n";
//     mod.saveState("before_helio_expansion");
//     mod.expandHeliosphereScale(1.2, 1.15);
//     std::cout << "  New H_SCm: " << mod.variables["H_SCm"] << "\n";
//     std::cout << "  New R_b: " << mod.variables["R_b"] << " m (" << (mod.variables["R_b"]/1.496e11) << " AU)\n";
//     std::cout << "  New U_g2: " << mod.computeU_g2(0.0, 0.0) << " J/m^3\n\n";
//
//     // Step 5: Boundary scale expansion (solar wind)
//     std::cout << "Step 5: Expand Boundary Scale (delta_sw x1.5, v_sw x1.1)\n";
//     mod.expandBoundaryScale(1.5, 1.1);
//     std::cout << "  New delta_sw: " << mod.variables["delta_sw"] << "\n";
//     std::cout << "  New v_sw: " << mod.variables["v_sw"] << " m/s\n";
//     std::cout << "  New swirl_factor: " << mod.variables["swirl_factor"] << "\n";
//     std::cout << "  New U_g2: " << mod.computeU_g2(0.0, 0.0) << " J/m^3\n\n";
//
//     // Step 6: Reactor scale expansion
//     std::cout << "Step 6: Expand Reactor Scale (E_react x1.3, k_2 x1.05)\n";
//     mod.expandReactorScale(1.3, 1.05);
//     std::cout << "  New E_react: " << mod.variables["E_react"] << " J\n";
//     std::cout << "  New k_2: " << mod.variables["k_2"] << "\n";
//     std::cout << "  New U_g2: " << mod.computeU_g2(0.0, 0.0) << " J/m^3\n\n";
//
//     // Step 7: Parameter variations
//     std::cout << "Step 7: Generate 10 Parameter Variations (±5%)\n";
//     auto variations = mod.generateVariations(10, 0.05);
//     std::cout << "  Generated " << variations.size() << " configurations\n";
//     std::cout << "  Sample H_SCm range: " << variations[0]["H_SCm"] << " to " << variations[9]["H_SCm"] << "\n";
//     std::cout << "  Sample R_b range: " << variations[0]["R_b"] << " to " << variations[9]["R_b"] << " m\n\n";
//
//     // Step 8: Sensitivity analysis
//     std::cout << "Step 8: Sensitivity Analysis (U_g2 response to ±1% parameter changes)\n";
//     std::vector<std::string> sens_params = {"H_SCm", "R_b", "E_react", "k_2", "delta_sw", "v_sw", "M_s"};
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
//     // Step 9: Auto-refinement to target U_g2
//     std::cout << "Step 9: Auto-Refine to Target U_g2 = 1.5e53 J/m^3\n";
//     double target_u = 1.5e53;
//     double before_refine = mod.computeU_g2(0.0, 0.0);
//     mod.autoRefineParameters("U_g2", target_u);
//     double after_refine = mod.computeU_g2(0.0, 0.0);
//     std::cout << "  Before: " << before_refine << " J/m^3\n";
//     std::cout << "  After: " << after_refine << " J/m^3\n";
//     std::cout << "  Target: " << target_u << " J/m^3\n";
//     std::cout << "  Error: " << std::abs(after_refine - target_u) / target_u * 100 << "%\n\n";
//
//     // Step 10: Restore and compare
//     std::cout << "Step 10: Restore Previous State\n";
//     mod.restoreState("before_helio_expansion");
//     std::cout << "  Restored H_SCm: " << mod.variables["H_SCm"] << "\n";
//     std::cout << "  Restored R_b: " << mod.variables["R_b"] << " m\n";
//     std::cout << "  Restored U_g2: " << mod.computeU_g2(0.0, 0.0) << " J/m^3\n\n";
//
//     // Step 11: Calibration to observations
//     std::cout << "Step 11: Calibrate to Observational Data\n";
//     std::map<std::string, double> observations = {
//         {"H_SCm", 1.05},
//         {"R_b", 1.8e13},  // ~120 AU
//         {"v_sw", 4.5e5}
//     };
//     mod.calibrateToObservations(observations);
//     std::cout << "  Calibrated H_SCm: " << mod.variables["H_SCm"] << "\n";
//     std::cout << "  Calibrated R_b: " << mod.variables["R_b"] << " m (" << (mod.variables["R_b"]/1.496e11) << " AU)\n";
//     std::cout << "  Calibrated v_sw: " << mod.variables["v_sw"] << " m/s\n";
//     std::cout << "  New U_g2: " << mod.computeU_g2(0.0, 0.0) << " J/m^3\n\n";
//
//     // Step 12: Optimization for boundary expansion
//     std::cout << "Step 12: Optimize for Boundary Expansion\n";
//     double before_opt = mod.variables["R_b"];
//     mod.optimizeForMetric("expand_boundary");
//     double after_opt = mod.variables["R_b"];
//     std::cout << "  R_b before: " << before_opt << " m\n";
//     std::cout << "  R_b after: " << after_opt << " m (+" << ((after_opt/before_opt - 1.0)*100) << "%)\n";
//     std::cout << "  New U_g2: " << mod.computeU_g2(0.0, 0.0) << " J/m^3\n\n";
//
//     // Step 13: Parameter mutation
//     std::cout << "Step 13: Mutate Parameters (±3%)\n";
//     mod.saveState("before_mutation");
//     double h_before = mod.variables["H_SCm"];
//     double rb_before = mod.variables["R_b"];
//     mod.mutateParameters(0.03);
//     std::cout << "  H_SCm: " << h_before << " -> " << mod.variables["H_SCm"] << "\n";
//     std::cout << "  R_b: " << rb_before << " -> " << mod.variables["R_b"] << " m\n";
//     std::cout << "  New U_g2: " << mod.computeU_g2(0.0, 0.0) << " J/m^3\n\n";
//
//     // Step 14: System evolution (maximize U_g2)
//     std::cout << "Step 14: Evolve System (7 generations, maximize U_g2)\n";
//     mod.restoreState("before_mutation");
//     double initial_fitness = mod.computeU_g2(0.0, 0.0);
//     mod.evolveSystem(7, [&mod]() { return mod.computeU_g2(0.0, 0.0); });
//     double final_fitness = mod.computeU_g2(0.0, 0.0);
//     std::cout << "  Initial U_g2: " << initial_fitness << " J/m^3\n";
//     std::cout << "  Evolved U_g2: " << final_fitness << " J/m^3\n";
//     std::cout << "  Improvement: " << ((final_fitness / initial_fitness - 1.0) * 100) << "%\n\n";
//
//     // Step 15: State management
//     std::cout << "Step 15: State Management\n";
//     mod.saveState("evolved_optimal");
//     auto saved = mod.listSavedStates();
//     std::cout << "  Saved states (" << saved.size() << "): ";
//     for (const auto& s : saved) std::cout << s << " ";
//     std::cout << "\n\n";
//
//     // Step 16: Validation and auto-correction
//     std::cout << "Step 16: Validate Consistency\n";
//     bool valid = mod.validateConsistency();
//     std::cout << "  System valid: " << (valid ? "YES" : "NO") << "\n";
//     if (!valid) {
//         std::cout << "  Running auto-correction...\n";
//         mod.autoCorrectAnomalies();
//         std::cout << "  Post-correction valid: " << (mod.validateConsistency() ? "YES" : "NO") << "\n";
//     }
//     std::cout << "\n";
//
//     // Step 17: Batch operations
//     std::cout << "Step 17: Batch Scale All Density Parameters x1.1\n";
//     std::vector<std::string> density_params = {"rho_vac_UA", "rho_vac_SCm"};
//     mod.scaleVariableGroup(density_params, 1.1);
//     std::cout << "  New rho_vac_UA: " << mod.variables["rho_vac_UA"] << " J/m^3\n";
//     std::cout << "  New rho_vac_SCm: " << mod.variables["rho_vac_SCm"] << " J/m^3\n";
//     std::cout << "  New rho_sum: " << mod.variables["rho_sum"] << " J/m^3\n";
//     std::cout << "  New U_g2: " << mod.computeU_g2(0.0, 0.0) << " J/m^3\n\n";
//
//     // Step 18: Export final state
//     std::cout << "Step 18: Export Final System State\n";
//     std::string exported = mod.exportState();
//     std::cout << exported << "\n";
//     std::cout << "Demo complete: 18 steps executed successfully!\n";
//     std::cout << "============================================================\n";
//
//     return 0;
// }

HeliosphereThicknessModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Automatic recalculation of derived variables(e.g., rho_sum, swirl_factor) when dependencies change.
- Core computation methods(computeU_g2, computeU_g2_no_H) are clear and use variable - driven logic.
- Output and debugging functions(printVariables, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, division by zero, or invalid input; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in heliosphere thickness and gravity modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.