// OuterFieldBubbleModule.h
// Modular C++ implementation of the Radius of the Outer Field Bubble (R_b) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes R_b=1.496e13 m (100 AU); defines S(r - R_b) step function in Universal Gravity U_g2 term.
// Pluggable: #include "OuterFieldBubbleModule.h"
// OuterFieldBubbleModule mod; mod.computeU_g2(1.5e13); mod.updateVariable("R_b", new_value);
// Variables in std::map; example for Sun at t=0; S=1 for r >= R_b, 0 otherwise.
// Approximations: S step=1 at r=R_b; ?_sw v_sw=5001; E_react=1e46.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef OUTER_FIELD_BUBBLE_MODULE_H
#define OUTER_FIELD_BUBBLE_MODULE_H

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

class OuterFieldBubbleModule {
private:
    std::map<std::string, double> variables;
    double computeS_r_Rb(double r);
    double computeU_g2(double r);

public:
    // Constructor: Initialize with framework defaults (Sun)
    OuterFieldBubbleModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeR_b();  // 1.496e13 m (100 AU)
    double computeR_bInAU();  // 100 AU
    double computeS_r_Rb(double r);  // Step function
    double computeU_g2(double r);  // U_g2 with S(r - R_b) (J/m^3)

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
    void expandParameterSpace(double bubble_scale, double density_scale, double energy_scale);
    void expandBubbleScale(double rb_factor, double k2_factor);               // Scale R_b and k_2
    void expandDensityScale(double rho_ua_factor, double rho_scm_factor);    // Scale vacuum densities
    void expandSolarWindScale(double vsw_factor, double delta_factor);        // Scale solar wind parameters

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

#endif // OUTER_FIELD_BUBBLE_MODULE_H

// OuterFieldBubbleModule.cpp
#include "OuterFieldBubbleModule.h"

// Constructor: Set framework defaults (Sun at t=0)
OuterFieldBubbleModule::OuterFieldBubbleModule() {
    // Universal constants
    variables["R_b"] = 1.496e13;                    // m (100 AU)
    variables["AU_to_m"] = 1.496e11;                // m/AU
    variables["k_2"] = 1.2;                         // Coupling
    variables["rho_vac_UA"] = 7.09e-36;             // J/m^3
    variables["rho_vac_SCm"] = 7.09e-37;            // J/m^3
    variables["M_s"] = 1.989e30;                    // kg
    variables["r"] = 1.496e13;                      // m (default = R_b)
    variables["delta_sw"] = 0.01;                   // Unitless
    variables["v_sw"] = 5e5;                        // m/s
    variables["H_SCm"] = 1.0;                       // Unitless
    variables["E_react"] = 1e46;                    // J

    // Derived
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
    variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
}

// Update variable
void OuterFieldBubbleModule::updateVariable(const std::string& name, double value) {
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
void OuterFieldBubbleModule::addToVariable(const std::string& name, double delta) {
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
void OuterFieldBubbleModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute R_b (m)
double OuterFieldBubbleModule::computeR_b() {
    return variables["R_b"];
}

// R_b in AU
double OuterFieldBubbleModule::computeR_bInAU() {
    return computeR_b() / variables["AU_to_m"];
}

// Step function S(r - R_b): 1 if r >= R_b, 0 otherwise
double OuterFieldBubbleModule::computeS_r_Rb(double r) {
    return (r >= computeR_b()) ? 1.0 : 0.0;
}

// Compute U_g2 with S(r - R_b)
double OuterFieldBubbleModule::computeU_g2(double r) {
    variables["r"] = r;
    double k_2 = variables["k_2"];
    double rho_sum = variables["rho_sum"];
    double M_s = variables["M_s"];
    double s_step = computeS_r_Rb(r);
    double swirl_factor = variables["swirl_factor"];
    double h_scm = variables["H_SCm"];
    double e_react = variables["E_react"];
    return k_2 * (rho_sum * M_s / (r * r)) * s_step * swirl_factor * h_scm * e_react;
}

// Equation text
std::string OuterFieldBubbleModule::getEquationText() {
    return "U_g2 = k_2 * [(?_vac,[UA] + ?_vac,[SCm]) M_s / r^2] * S(r - R_b) * (1 + ?_sw v_sw) * H_SCm * E_react\n"
           "Where R_b = 1.496e13 m (100 AU, outer bubble radius);\n"
           "S(r - R_b) = 1 (r >= R_b), 0 otherwise (step function).\n"
           "Example r=R_b: U_g2 ?1.18e53 J/m�; r < R_b (e.g., 1 AU): U_g2=0.\n"
           "Role: Defines external gravity boundary (~heliopause); activates U_g2 beyond R_b.\n"
           "UQFF: Separates internal/external fields; models heliosphere/nebular extent.";
}

// Print variables
void OuterFieldBubbleModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// ===== ENHANCED METHODS IMPLEMENTATION =====

namespace outer_field_bubble_saved_states {
    static std::map<std::string, std::map<std::string, double>> saved_states;
}

// Variable Management
void OuterFieldBubbleModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void OuterFieldBubbleModule::removeVariable(const std::string& name) {
    variables.erase(name);
}

void OuterFieldBubbleModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
    }
}

std::vector<std::string> OuterFieldBubbleModule::listVariables() const {
    std::vector<std::string> names;
    for (const auto& pair : variables) {
        names.push_back(pair.first);
    }
    return names;
}

std::string OuterFieldBubbleModule::getSystemName() const {
    return "Outer_Field_Bubble_UQFF";
}

// Batch Operations
void OuterFieldBubbleModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = func(variables[name]);
        }
    }
}

void OuterFieldBubbleModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// Self-Expansion: Domain-Specific Scales
void OuterFieldBubbleModule::expandParameterSpace(double bubble_scale, double density_scale, double energy_scale) {
    variables["R_b"] *= bubble_scale;
    variables["rho_vac_UA"] *= density_scale;
    variables["rho_vac_SCm"] *= density_scale;
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
    variables["E_react"] *= energy_scale;
}

void OuterFieldBubbleModule::expandBubbleScale(double rb_factor, double k2_factor) {
    variables["R_b"] *= rb_factor;
    variables["k_2"] *= k2_factor;
}

void OuterFieldBubbleModule::expandDensityScale(double rho_ua_factor, double rho_scm_factor) {
    variables["rho_vac_UA"] *= rho_ua_factor;
    variables["rho_vac_SCm"] *= rho_scm_factor;
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
}

void OuterFieldBubbleModule::expandSolarWindScale(double vsw_factor, double delta_factor) {
    variables["v_sw"] *= vsw_factor;
    variables["delta_sw"] *= delta_factor;
    variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
}

// Self-Refinement
void OuterFieldBubbleModule::autoRefineParameters(const std::string& target, double goal) {
    if (target == "R_b") {
        // Target specific R_b directly
        variables["R_b"] = goal;
    } else if (target == "R_b_AU") {
        // Target R_b in AU
        variables["R_b"] = goal * variables["AU_to_m"];
    } else if (target == "U_g2") {
        // Target specific U_g2 at R_b by scaling E_react
        double current_ug2 = computeU_g2(variables["R_b"]);
        if (std::abs(current_ug2) > 1e-9) {
            variables["E_react"] *= (goal / current_ug2);
        }
    } else if (target == "step_boundary") {
        // Adjust R_b so that step function activates at specific radius
        variables["R_b"] = goal;
    }
}

void OuterFieldBubbleModule::calibrateToObservations(const std::map<std::string, double>& observations) {
    for (const auto& obs : observations) {
        if (variables.find(obs.first) != variables.end()) {
            variables[obs.first] = obs.second;
            // Update dependent variables
            if (obs.first == "rho_vac_UA" || obs.first == "rho_vac_SCm") {
                variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
            } else if (obs.first == "delta_sw" || obs.first == "v_sw") {
                variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
            }
        }
    }
}

void OuterFieldBubbleModule::optimizeForMetric(const std::string& metric) {
    if (metric == "maximize_bubble") {
        // Increase bubble radius (larger heliosphere)
        variables["R_b"] *= 1.3;
    } else if (metric == "minimize_bubble") {
        // Decrease bubble radius (compressed heliosphere)
        variables["R_b"] *= 0.7;
    } else if (metric == "enhance_external_field") {
        // Enhance U_g2 strength
        variables["k_2"] *= 1.3;
        variables["E_react"] *= 1.2;
    } else if (metric == "standard_heliopause") {
        // Reset to standard 100 AU heliopause
        variables["R_b"] = 1.496e13;
    } else if (metric == "voyager_crossing") {
        // Set to Voyager 1 heliopause crossing (~122 AU)
        variables["R_b"] = 122.0 * variables["AU_to_m"];
    }
}

// Parameter Exploration
std::vector<std::map<std::string, double>> OuterFieldBubbleModule::generateVariations(int count, double variation_pct) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - variation_pct, 1.0 + variation_pct);
    
    for (int i = 0; i < count; ++i) {
        std::map<std::string, double> variant = variables;
        for (auto& pair : variant) {
            if (pair.first != "AU_to_m" && pair.first != "rho_sum" && pair.first != "swirl_factor") {
                pair.second *= dis(gen);
            }
        }
        // Recalculate derived variables
        variant["rho_sum"] = variant["rho_vac_UA"] + variant["rho_vac_SCm"];
        variant["swirl_factor"] = 1.0 + variant["delta_sw"] * variant["v_sw"];
        variations.push_back(variant);
    }
    return variations;
}

// Adaptive Evolution
void OuterFieldBubbleModule::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - mutation_rate, 1.0 + mutation_rate);
    
    for (auto& pair : variables) {
        if (pair.first != "AU_to_m" && pair.first != "rho_sum" && pair.first != "swirl_factor") {
            pair.second *= dis(gen);
        }
    }
    // Recalculate derived variables
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
    variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
}

void OuterFieldBubbleModule::evolveSystem(int generations, std::function<double()> fitness_func) {
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
void OuterFieldBubbleModule::saveState(const std::string& label) {
    outer_field_bubble_saved_states::saved_states[label] = variables;
}

void OuterFieldBubbleModule::restoreState(const std::string& label) {
    if (outer_field_bubble_saved_states::saved_states.find(label) != outer_field_bubble_saved_states::saved_states.end()) {
        variables = outer_field_bubble_saved_states::saved_states[label];
    }
}

std::vector<std::string> OuterFieldBubbleModule::listSavedStates() const {
    std::vector<std::string> labels;
    for (const auto& pair : outer_field_bubble_saved_states::saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string OuterFieldBubbleModule::exportState() const {
    std::ostringstream oss;
    oss << "OuterFieldBubble_State_Export:\n";
    for (const auto& pair : variables) {
        oss << pair.first << "=" << std::scientific << pair.second << "\n";
    }
    return oss.str();
}

// System Analysis
std::map<std::string, double> OuterFieldBubbleModule::sensitivityAnalysis(const std::vector<std::string>& params) {
    std::map<std::string, double> sensitivities;
    double r = variables["R_b"];  // Test at bubble boundary
    double baseline_ug2 = computeU_g2(r);
    
    for (const auto& param : params) {
        if (variables.find(param) != variables.end() && param != "AU_to_m") {
            double original = variables[param];
            variables[param] = original * 1.01;
            
            // Update derived variables if needed
            if (param == "rho_vac_UA" || param == "rho_vac_SCm") {
                variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
            } else if (param == "delta_sw" || param == "v_sw") {
                variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
            }
            
            double perturbed_ug2 = computeU_g2(r);
            sensitivities[param] = (perturbed_ug2 - baseline_ug2) / baseline_ug2;
            
            // Restore
            variables[param] = original;
            if (param == "rho_vac_UA" || param == "rho_vac_SCm") {
                variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
            } else if (param == "delta_sw" || param == "v_sw") {
                variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
            }
        }
    }
    return sensitivities;
}

std::string OuterFieldBubbleModule::generateReport() const {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(3);
    oss << "===== Outer Field Bubble Module Report =====\n";
    oss << "System: " << getSystemName() << "\n\n";
    
    oss << "Bubble Parameters:\n";
    oss << "  R_b (bubble radius) = " << variables.at("R_b") << " m\n";
    double rb_au = variables.at("R_b") / variables.at("AU_to_m");
    oss << "  R_b = " << std::fixed << std::setprecision(1) << rb_au << " AU\n";
    oss << "  k_2 (coupling) = " << std::scientific << variables.at("k_2") << "\n\n";
    
    oss << "Vacuum Densities:\n";
    oss << "  ρ_vac,[UA] = " << variables.at("rho_vac_UA") << " J/m³\n";
    oss << "  ρ_vac,[SCm] = " << variables.at("rho_vac_SCm") << " J/m³\n";
    oss << "  ρ_sum = ρ_UA + ρ_SCm = " << variables.at("rho_sum") << " J/m³\n\n";
    
    oss << "Solar Wind Parameters:\n";
    oss << "  v_sw (solar wind speed) = " << variables.at("v_sw") << " m/s\n";
    oss << "  δ_sw (swirl parameter) = " << variables.at("delta_sw") << "\n";
    oss << "  Swirl factor = 1 + δ_sw × v_sw = " << variables.at("swirl_factor") << "\n\n";
    
    oss << "Additional Parameters:\n";
    oss << "  M_s (solar mass) = " << variables.at("M_s") << " kg\n";
    oss << "  H_SCm (heliosphere thickness) = " << variables.at("H_SCm") << "\n";
    oss << "  E_react = " << variables.at("E_react") << " J\n\n";
    
    oss << "Step Function S(r - R_b):\n";
    std::vector<double> test_radii = {0.5 * variables.at("R_b"), 0.9 * variables.at("R_b"), 
                                       variables.at("R_b"), 1.1 * variables.at("R_b"), 
                                       1.5 * variables.at("R_b")};
    for (double r : test_radii) {
        double r_au = r / variables.at("AU_to_m");
        double s_val = (r >= variables.at("R_b")) ? 1.0 : 0.0;
        oss << "  r=" << std::fixed << std::setprecision(1) << r_au << " AU: S=" << s_val;
        if (s_val == 0.0) {
            oss << " (inside bubble, U_g2=0)\n";
        } else {
            oss << " (outside bubble, U_g2 active)\n";
        }
    }
    oss << "\n" << std::scientific;
    
    oss << "U_g2 Computations:\n";
    // Inside bubble (r < R_b)
    double r_inside = 0.5 * variables.at("R_b");
    double ug2_inside = 0.0;  // Step function = 0
    oss << "  r=" << std::fixed << std::setprecision(1) << (r_inside / variables.at("AU_to_m")) 
        << " AU (inside): U_g2=" << std::scientific << ug2_inside << " J/m³\n";
    
    // At boundary (r = R_b)
    double r_boundary = variables.at("R_b");
    double k_2 = variables.at("k_2");
    double rho_sum = variables.at("rho_sum");
    double M_s = variables.at("M_s");
    double swirl_factor = variables.at("swirl_factor");
    double h_scm = variables.at("H_SCm");
    double e_react = variables.at("E_react");
    double ug2_boundary = k_2 * (rho_sum * M_s / (r_boundary * r_boundary)) * 1.0 * swirl_factor * h_scm * e_react;
    oss << "  r=" << std::fixed << std::setprecision(1) << (r_boundary / variables.at("AU_to_m"))
        << " AU (boundary): U_g2=" << std::scientific << ug2_boundary << " J/m³\n";
    
    // Outside bubble (r > R_b)
    double r_outside = 1.5 * variables.at("R_b");
    double ug2_outside = k_2 * (rho_sum * M_s / (r_outside * r_outside)) * 1.0 * swirl_factor * h_scm * e_react;
    oss << "  r=" << std::fixed << std::setprecision(1) << (r_outside / variables.at("AU_to_m"))
        << " AU (outside): U_g2=" << std::scientific << ug2_outside << " J/m³\n\n";
    
    oss << "Physical Interpretation:\n";
    if (rb_au < 80.0) {
        oss << "  Compressed heliosphere (<80 AU)\n";
    } else if (rb_au < 110.0) {
        oss << "  Typical heliosphere (~100 AU)\n";
    } else if (rb_au < 150.0) {
        oss << "  Extended heliosphere (110-150 AU, like Voyager crossing)\n";
    } else {
        oss << "  Very extended heliosphere (>150 AU)\n";
    }
    
    oss << "  Applications:\n";
    oss << "    - Heliopause boundary: Separates solar wind from interstellar medium\n";
    oss << "    - U_g2 activation: External gravity field beyond R_b\n";
    oss << "    - Stellar winds: Defines stellar influence sphere\n";
    oss << "    - Nebular boundaries: Gas cloud extent modeling\n";
    
    return oss.str();
}

bool OuterFieldBubbleModule::validateConsistency() const {
    bool valid = true;
    
    // Check R_b is positive and reasonable
    if (variables.find("R_b") != variables.end()) {
        double rb = variables.at("R_b");
        double rb_au = rb / variables.at("AU_to_m");
        if (rb <= 0) {
            std::cerr << "Error: R_b <= 0 (bubble radius must be positive)\n";
            valid = false;
        }
        if (rb_au < 10.0 || rb_au > 1000.0) {
            std::cerr << "Warning: R_b outside typical range [10, 1000] AU (current: " << rb_au << " AU)\n";
        }
    }
    
    // Check k_2 is positive
    if (variables.find("k_2") != variables.end() && variables.at("k_2") <= 0) {
        std::cerr << "Error: k_2 <= 0 (coupling constant must be positive)\n";
        valid = false;
    }
    
    // Check vacuum densities are positive
    if (variables.find("rho_vac_UA") != variables.end() && variables.at("rho_vac_UA") <= 0) {
        std::cerr << "Error: rho_vac_UA <= 0 (vacuum density must be positive)\n";
        valid = false;
    }
    
    // Check rho_sum consistency
    if (variables.find("rho_sum") != variables.end()) {
        double expected = variables.at("rho_vac_UA") + variables.at("rho_vac_SCm");
        double actual = variables.at("rho_sum");
        if (std::abs(expected - actual) / expected > 1e-9) {
            std::cerr << "Error: rho_sum inconsistent (expected " << expected << ", got " << actual << ")\n";
            valid = false;
        }
    }
    
    return valid;
}

void OuterFieldBubbleModule::autoCorrectAnomalies() {
    // Reset R_b to typical value if out of range
    if (variables["R_b"] <= 0 || variables["R_b"] / variables["AU_to_m"] > 1000.0) {
        variables["R_b"] = 1.496e13;  // 100 AU
    }
    
    // Ensure k_2 is positive
    if (variables["k_2"] <= 0) {
        variables["k_2"] = 1.2;
    }
    
    // Ensure vacuum densities are positive
    if (variables["rho_vac_UA"] <= 0) {
        variables["rho_vac_UA"] = 7.09e-36;
    }
    if (variables["rho_vac_SCm"] <= 0) {
        variables["rho_vac_SCm"] = 7.09e-37;
    }
    
    // Recalculate derived variables
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
    variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
}

// Example usage in base program (snippet)
int main() {
    OuterFieldBubbleModule module;
    std::cout << "===== Outer Field Bubble Module Enhanced Demonstration =====\n\n";
    
    // Step 1: Report initial state
    std::cout << "STEP 1: Initial Configuration (Standard 100 AU Heliosphere)\n";
    std::cout << module.generateReport() << "\n";
    
    // Step 2: Compute U_g2 at various radii
    std::cout << "STEP 2: U_g2 Computation at Different Radii\n";
    double AU_to_m = 1.496e11;
    std::vector<double> test_radii_au = {50.0, 90.0, 100.0, 120.0, 150.0};
    for (double r_au : test_radii_au) {
        double r = r_au * AU_to_m;
        double ug2 = module.computeU_g2(r);
        std::cout << "  r=" << r_au << " AU: U_g2=" << std::scientific << ug2 << " J/m³";
        if (r < module.getVariable("R_b")) {
            std::cout << " (inside bubble, S=0)\n";
        } else {
            std::cout << " (outside bubble, S=1)\n";
        }
    }
    std::cout << "\n";
    
    // Step 3: Save initial state
    std::cout << "STEP 3: Save Initial State\n";
    module.saveState("standard_100AU");
    std::cout << "State saved as 'standard_100AU'\n\n";
    
    // Step 4: Test compressed heliosphere (50 AU)
    std::cout << "STEP 4: Test Compressed Heliosphere (50 AU)\n";
    module.createVariable("R_b_new", 50.0 * AU_to_m);
    module.expandBubbleScale(0.5, 1.0);  // Compress R_b by 50%
    double ug2_compressed = module.computeU_g2(50.0 * AU_to_m);
    std::cout << "Compressed R_b = 50 AU\n";
    std::cout << "U_g2 at new boundary = " << std::scientific << ug2_compressed << " J/m³\n\n";
    
    // Step 5: Test extended heliosphere (Voyager scale, 122 AU)
    std::cout << "STEP 5: Test Extended Heliosphere (Voyager 1 Crossing, 122 AU)\n";
    module.restoreState("standard_100AU");
    module.optimizeForMetric("voyager_crossing");
    double ug2_voyager = module.computeU_g2(module.getVariable("R_b"));
    std::cout << "Extended R_b = " << (module.getVariable("R_b") / AU_to_m) << " AU\n";
    std::cout << "U_g2 at Voyager boundary = " << std::scientific << ug2_voyager << " J/m³\n\n";
    
    // Step 6: Test very extended heliosphere (150 AU)
    std::cout << "STEP 6: Test Very Extended Heliosphere (150 AU)\n";
    module.restoreState("standard_100AU");
    module.expandBubbleScale(1.5, 1.0);  // Expand R_b by 50%
    double ug2_extended = module.computeU_g2(module.getVariable("R_b"));
    std::cout << "Very extended R_b = " << (module.getVariable("R_b") / AU_to_m) << " AU\n";
    std::cout << "U_g2 at extended boundary = " << std::scientific << ug2_extended << " J/m³\n\n";
    
    // Step 7: Restore and expand density scale
    std::cout << "STEP 7: Expand Vacuum Density Scale (UA and SCm x2)\n";
    module.restoreState("standard_100AU");
    module.expandDensityScale(2.0, 2.0);
    double ug2_dense = module.computeU_g2(100.0 * AU_to_m);
    std::cout << "New rho_sum = " << std::scientific << module.getVariable("rho_sum") << " J/m³\n";
    std::cout << "U_g2 at 100 AU (doubled density) = " << ug2_dense << " J/m³\n\n";
    
    // Step 8: Restore and expand solar wind scale
    std::cout << "STEP 8: Expand Solar Wind Scale (v_sw x1.5, δ_sw x2)\n";
    module.restoreState("standard_100AU");
    module.expandSolarWindScale(1.5, 2.0);
    double new_swirl = module.getVariable("swirl_factor");
    double ug2_wind = module.computeU_g2(100.0 * AU_to_m);
    std::cout << "New swirl_factor = " << std::scientific << new_swirl << "\n";
    std::cout << "U_g2 at 100 AU (enhanced wind) = " << ug2_wind << " J/m³\n\n";
    
    // Step 9: Sensitivity analysis
    std::cout << "STEP 9: Sensitivity Analysis\n";
    module.restoreState("standard_100AU");
    std::vector<std::string> params = {"R_b", "k_2", "rho_vac_UA", "v_sw", "E_react"};
    auto sensitivities = module.sensitivityAnalysis(params);
    for (const auto& pair : sensitivities) {
        std::cout << "  ∂U_g2/∂" << pair.first << " ≈ " << std::scientific << pair.second << " (normalized)\n";
    }
    std::cout << "\n";
    
    // Step 10: Generate variations
    std::cout << "STEP 10: Generate Parameter Variations (5 variants, ±10%)\n";
    auto variations = module.generateVariations(5, 0.1);
    for (int i = 0; i < variations.size(); ++i) {
        double var_rb = variations[i]["R_b"] / AU_to_m;
        std::cout << "  Variant " << (i+1) << ": R_b=" << std::fixed << std::setprecision(1) 
                  << var_rb << " AU\n";
    }
    std::cout << "\n";
    
    // Step 11: Auto-refine to target R_b
    std::cout << "STEP 11: Auto-Refine to Target R_b = 110 AU\n";
    module.restoreState("standard_100AU");
    module.autoRefineParameters("R_b_AU", 110.0);
    double refined_rb = module.getVariable("R_b") / AU_to_m;
    double ug2_refined = module.computeU_g2(module.getVariable("R_b"));
    std::cout << "Refined R_b = " << std::fixed << std::setprecision(1) << refined_rb << " AU\n";
    std::cout << "U_g2 at refined boundary = " << std::scientific << ug2_refined << " J/m³\n\n";
    
    // Step 12: Calibrate to observations
    std::cout << "STEP 12: Calibrate to Observational Data\n";
    module.restoreState("standard_100AU");
    std::map<std::string, double> observations;
    observations["R_b"] = 95.0 * AU_to_m;  // Observed slightly smaller than 100 AU
    observations["v_sw"] = 4.5e5;  // Slightly slower solar wind
    module.calibrateToObservations(observations);
    std::cout << "Calibrated R_b = " << (module.getVariable("R_b") / AU_to_m) << " AU\n";
    std::cout << "Calibrated v_sw = " << std::scientific << module.getVariable("v_sw") << " m/s\n\n";
    
    // Step 13: Optimize for metric
    std::cout << "STEP 13: Optimize for 'standard_heliopause' Metric\n";
    module.optimizeForMetric("standard_heliopause");
    std::cout << "Optimized R_b = " << (module.getVariable("R_b") / AU_to_m) << " AU\n\n";
    
    // Step 14: Mutate parameters
    std::cout << "STEP 14: Mutate Parameters (5% random variation)\n";
    module.mutateParameters(0.05);
    std::cout << "Mutated R_b = " << (module.getVariable("R_b") / AU_to_m) << " AU\n\n";
    
    // Step 15: Validate consistency
    std::cout << "STEP 15: Validate Consistency\n";
    bool valid = module.validateConsistency();
    std::cout << "Consistency check: " << (valid ? "PASSED" : "FAILED") << "\n\n";
    
    // Step 16: Introduce anomaly and auto-correct
    std::cout << "STEP 16: Introduce Anomaly and Auto-Correct\n";
    module.createVariable("R_b_anomaly", -100.0);  // Invalid negative R_b
    module.removeVariable("R_b");
    module.createVariable("R_b", -100.0);
    std::cout << "Introduced invalid R_b = " << module.getVariable("R_b") << " m\n";
    module.autoCorrectAnomalies();
    std::cout << "Auto-corrected R_b = " << (module.getVariable("R_b") / AU_to_m) << " AU\n";
    bool valid_after = module.validateConsistency();
    std::cout << "Consistency after correction: " << (valid_after ? "PASSED" : "FAILED") << "\n\n";
    
    // Step 17: List saved states
    std::cout << "STEP 17: List Saved States\n";
    auto states = module.listSavedStates();
    for (const auto& state : states) {
        std::cout << "  - " << state << "\n";
    }
    std::cout << "\n";
    
    // Step 18: Export final state
    std::cout << "STEP 18: Export Final State\n";
    std::string exported = module.exportState();
    std::cout << exported << "\n";
    
    std::cout << "===== Demonstration Complete =====\n";
    return 0;
}

// Original commented example preserved below for reference:
// #include "OuterFieldBubbleModule.h"
// int main() {
//     OuterFieldBubbleModule mod;
//     double rb = mod.computeR_b();
//     std::cout << "R_b = " << rb << " m (" << mod.computeR_bInAU() << " AU)\n";
//     double u_g2 = mod.computeU_g2(1.5e13);  // r > R_b
//     std::cout << "U_g2 (r=1.5e13 m) = " << u_g2 << " J/m³\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("R_b", 2e13);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o bubble_test bubble_test.cpp OuterFieldBubbleModule.cpp -lm
// Sample: R_b=1.496e13 m (100 AU); U_g2≈1.18e53 J/m³ (r>=R_b); 0 inside.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

OuterFieldBubbleModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeR_b, computeR_bInAU, computeS_r_Rb, computeU_g2) are clear, concise, and variable - driven.
- Automatic recalculation of derived variables(rho_sum, swirl_factor) when dependencies change.
- Output and debugging functions(printVariables, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.
- Step function S(r - R_b) cleanly separates internal and external field regions.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in outer field bubble modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.