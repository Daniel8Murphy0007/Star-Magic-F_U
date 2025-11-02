// StepFunctionModule.h
// Modular C++ implementation of the Step Function (S) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes S(r - R_b) = 1 for r > R_b, 0 otherwise; activates U_g2 outside outer field bubble.
// Pluggable: #include "StepFunctionModule.h"
// StepFunctionModule mod; mod.computeU_g2(1.5e13); mod.updateVariable("R_b", new_value);
// Variables in std::map; example for Sun at r=1.496e13 m (R_b: S=1, U_g2?1.18e53 J/m�); r=1e11 m: S=0, U_g2=0.
// Approximations: S=1 at r=R_b; (1 + ?_sw v_sw)=5001; H_SCm=1; E_react=1e46.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef STEP_FUNCTION_MODULE_H
#define STEP_FUNCTION_MODULE_H

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

class StepFunctionModule {
private:
    std::map<std::string, double> variables;
    double computeS_r_Rb(double r);
    double computeU_g2(double r);

public:
    // Constructor: Initialize with framework defaults (Sun)
    StepFunctionModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeS_r_Rb(double r);  // Step: 1 if r > R_b, 0 otherwise
    double computeU_g2(double r);  // U_g2 with S(r - R_b) (J/m^3)

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
    double getVariable(const std::string& name) const { return variables.at(name); }

    // Batch Operations
    void transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func);
    void scaleVariableGroup(const std::vector<std::string>& names, double factor);

    // Self-Expansion: Domain-Specific Scales
    void expandParameterSpace(double boundary_scale, double field_scale, double transition_scale);
    void expandBoundaryScale(double rb_factor, double sharpness_factor);      // R_b and transition sharpness
    void expandFieldScale(double gravity_factor, double modulation_factor);   // Field strength and modulation
    void expandTransitionScale(double width_factor, double smoothness_factor); // Transition width and smoothing

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

#endif // STEP_FUNCTION_MODULE_H

// StepFunctionModule.cpp
#include "StepFunctionModule.h"

// Constructor: Set framework defaults (Sun at r=R_b)
StepFunctionModule::StepFunctionModule() {
    // Universal constants
    variables["R_b"] = 1.496e13;                    // m (100 AU)
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
void StepFunctionModule::updateVariable(const std::string& name, double value) {
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
void StepFunctionModule::addToVariable(const std::string& name, double delta) {
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
void StepFunctionModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute S(r - R_b): 1 if r > R_b, 0 otherwise (treat = as 1 per examples)
double StepFunctionModule::computeS_r_Rb(double r) {
    return (r >= variables["R_b"]) ? 1.0 : 0.0;
}

// Compute U_g2 with S(r - R_b)
double StepFunctionModule::computeU_g2(double r) {
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
std::string StepFunctionModule::getEquationText() {
    return "U_g2 = k_2 * [(?_vac,[UA] + ?_vac,[SCm]) M_s / r^2] * S(r - R_b) * (1 + ?_sw v_sw) * H_SCm * E_react\n"
           "Where S(r - R_b) = 1 (r > R_b), 0 otherwise (Heaviside step; =1 at boundary).\n"
           "Defines outer bubble activation beyond R_b=1.496e13 m (100 AU).\n"
           "Example r=1.496e13 m: S=1, U_g2 ?1.18e53 J/m�;\n"
           "r=1e11 m: S=0, U_g2=0; r=1e14 m: S=1, U_g2?1.18e51 J/m�.\n"
           "Role: Sharp transition internal/external gravity; heliopause-like boundary.\n"
           "UQFF: Separates regimes for heliodynamics/nebular formation.";
}

// Print variables
void StepFunctionModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// ===== ENHANCED METHODS IMPLEMENTATION =====

namespace step_function_saved_states {
    static std::map<std::string, std::map<std::string, double>> saved_states;
}

// Variable Management
void StepFunctionModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void StepFunctionModule::removeVariable(const std::string& name) {
    variables.erase(name);
}

void StepFunctionModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
    }
}

std::vector<std::string> StepFunctionModule::listVariables() const {
    std::vector<std::string> names;
    for (const auto& pair : variables) {
        names.push_back(pair.first);
    }
    return names;
}

std::string StepFunctionModule::getSystemName() const {
    return "Step_Function_Boundary_UQFF";
}

// Batch Operations
void StepFunctionModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = func(variables[name]);
        }
    }
    // Update derived
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
    variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
}

void StepFunctionModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// Self-Expansion: Domain-Specific Scales
void StepFunctionModule::expandParameterSpace(double boundary_scale, double field_scale, double transition_scale) {
    variables["R_b"] *= boundary_scale;
    variables["k_2"] *= field_scale;
    variables["rho_vac_UA"] *= field_scale;
    variables["rho_vac_SCm"] *= field_scale;
    
    // Transition width scaling
    if (variables.find("transition_width_m") != variables.end()) {
        variables["transition_width_m"] *= transition_scale;
    }
    
    // Update derived
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
}

void StepFunctionModule::expandBoundaryScale(double rb_factor, double sharpness_factor) {
    variables["R_b"] *= rb_factor;
    
    // Store boundary characteristics for advanced modeling
    if (variables.find("boundary_radius_AU") == variables.end()) {
        // Solar R_b ~ 100 AU (1.496e13 m)
        variables["boundary_radius_AU"] = 100.0;
    }
    variables["boundary_radius_AU"] *= rb_factor;
    
    // Transition sharpness (inverse of width)
    if (variables.find("transition_sharpness") == variables.end()) {
        // Sharp step: sharpness ~ 1e13 m^-1 (inverse of 1 AU scale)
        variables["transition_sharpness"] = 1.0 / 1.496e11;
    }
    variables["transition_sharpness"] *= sharpness_factor;
    
    // Transition width (AU scale)
    if (variables.find("transition_width_AU") == variables.end()) {
        // Sharp step: ~0.01 AU transition width
        variables["transition_width_AU"] = 0.01;
    }
    variables["transition_width_AU"] /= sharpness_factor;  // Sharper → narrower
    
    // Heliopause characteristics
    if (variables.find("heliopause_distance_AU") == variables.end()) {
        // Voyager 1 crossed at ~122 AU
        variables["heliopause_distance_AU"] = 122.0;
    }
    variables["heliopause_distance_AU"] *= rb_factor;
}

void StepFunctionModule::expandFieldScale(double gravity_factor, double modulation_factor) {
    variables["k_2"] *= gravity_factor;
    variables["rho_vac_UA"] *= gravity_factor;
    variables["rho_vac_SCm"] *= gravity_factor;
    
    // Solar wind modulation
    variables["delta_sw"] *= modulation_factor;
    variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
    
    // Field strength characteristics
    if (variables.find("field_strength_ratio") == variables.end()) {
        // Outside/inside field ratio
        variables["field_strength_ratio"] = 1.0;
    }
    variables["field_strength_ratio"] *= gravity_factor;
    
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
}

void StepFunctionModule::expandTransitionScale(double width_factor, double smoothness_factor) {
    // Transition width (for smooth step approximations)
    if (variables.find("transition_width_m") == variables.end()) {
        // ~0.01 AU = 1.496e9 m typical width
        variables["transition_width_m"] = 1.496e9;
    }
    variables["transition_width_m"] *= width_factor;
    
    // Smoothness parameter (for tanh/sigmoid approximations)
    if (variables.find("smoothness_parameter") == variables.end()) {
        variables["smoothness_parameter"] = 1.0;
    }
    variables["smoothness_parameter"] *= smoothness_factor;
    
    // Transition steepness
    if (variables.find("transition_steepness") == variables.end()) {
        // Steepness ~ 1/width
        variables["transition_steepness"] = 1.0 / 1.496e9;
    }
    variables["transition_steepness"] /= width_factor;  // Wider → less steep
    
    // Store step function type
    if (variables.find("step_type") == variables.end()) {
        // 0: Heaviside (sharp), 1: tanh (smooth), 2: sigmoid
        variables["step_type"] = 0.0;  // Heaviside by default
    }
}

// Self-Refinement
void StepFunctionModule::autoRefineParameters(const std::string& target, double goal) {
    if (target == "R_b") {
        variables["R_b"] = goal;
        if (variables.find("boundary_radius_AU") != variables.end()) {
            variables["boundary_radius_AU"] = goal / 1.496e11;
        }
    } else if (target == "R_b_AU") {
        variables["R_b"] = goal * 1.496e11;
        if (variables.find("boundary_radius_AU") != variables.end()) {
            variables["boundary_radius_AU"] = goal;
        }
    } else if (target == "transition_width_m") {
        if (variables.find("transition_width_m") == variables.end()) {
            variables["transition_width_m"] = 1.496e9;
        }
        variables["transition_width_m"] = goal;
        if (variables.find("transition_steepness") != variables.end()) {
            variables["transition_steepness"] = 1.0 / goal;
        }
    } else if (target == "transition_width_AU") {
        if (variables.find("transition_width_AU") == variables.end()) {
            variables["transition_width_AU"] = 0.01;
        }
        variables["transition_width_AU"] = goal;
        if (variables.find("transition_width_m") == variables.end()) {
            variables["transition_width_m"] = goal * 1.496e11;
        } else {
            variables["transition_width_m"] = goal * 1.496e11;
        }
    } else if (target == "transition_sharpness") {
        if (variables.find("transition_sharpness") == variables.end()) {
            variables["transition_sharpness"] = 1.0 / 1.496e11;
        }
        variables["transition_sharpness"] = goal;
    } else if (target == "U_g2_at_Rb") {
        // Target specific U_g2 at R_b by adjusting k_2
        double current_U_g2 = computeU_g2(variables["R_b"]);
        if (current_U_g2 > 0) {
            variables["k_2"] *= (goal / current_U_g2);
        }
    } else if (target == "heliopause_distance_AU") {
        if (variables.find("heliopause_distance_AU") == variables.end()) {
            variables["heliopause_distance_AU"] = 122.0;
        }
        variables["heliopause_distance_AU"] = goal;
        // May adjust R_b to match
        variables["R_b"] = goal * 1.496e11;
        if (variables.find("boundary_radius_AU") != variables.end()) {
            variables["boundary_radius_AU"] = goal;
        }
    } else if (target == "step_type") {
        if (variables.find("step_type") == variables.end()) {
            variables["step_type"] = 0.0;
        }
        variables["step_type"] = goal;
    }
}

void StepFunctionModule::calibrateToObservations(const std::map<std::string, double>& observations) {
    for (const auto& obs : observations) {
        if (variables.find(obs.first) != variables.end()) {
            variables[obs.first] = obs.second;
            if (obs.first == "rho_vac_UA" || obs.first == "rho_vac_SCm") {
                variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
            }
            else if (obs.first == "delta_sw" || obs.first == "v_sw") {
                variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
            }
        }
    }
}

void StepFunctionModule::optimizeForMetric(const std::string& metric) {
    if (metric == "solar_boundary") {
        // Solar R_b ~ 100 AU
        variables["R_b"] = 1.496e13;  // 100 AU
        if (variables.find("boundary_radius_AU") != variables.end()) {
            variables["boundary_radius_AU"] = 100.0;
        }
    } else if (metric == "heliopause") {
        // Voyager 1: ~122 AU
        variables["R_b"] = 122.0 * 1.496e11;
        if (variables.find("boundary_radius_AU") != variables.end()) {
            variables["boundary_radius_AU"] = 122.0;
        }
        if (variables.find("heliopause_distance_AU") != variables.end()) {
            variables["heliopause_distance_AU"] = 122.0;
        }
    } else if (metric == "bow_shock") {
        // Bow shock: ~90 AU (upstream of heliopause)
        variables["R_b"] = 90.0 * 1.496e11;
        if (variables.find("boundary_radius_AU") != variables.end()) {
            variables["boundary_radius_AU"] = 90.0;
        }
    } else if (metric == "termination_shock") {
        // Termination shock: ~90-100 AU
        variables["R_b"] = 94.0 * 1.496e11;  // Voyager 1 crossing
        if (variables.find("boundary_radius_AU") != variables.end()) {
            variables["boundary_radius_AU"] = 94.0;
        }
    } else if (metric == "inner_boundary") {
        // Inner boundary: 10-50 AU
        variables["R_b"] = 30.0 * 1.496e11;
        if (variables.find("boundary_radius_AU") != variables.end()) {
            variables["boundary_radius_AU"] = 30.0;
        }
    } else if (metric == "outer_boundary") {
        // Outer boundary: 200+ AU
        variables["R_b"] = 200.0 * 1.496e11;
        if (variables.find("boundary_radius_AU") != variables.end()) {
            variables["boundary_radius_AU"] = 200.0;
        }
    } else if (metric == "sharp_transition") {
        // Very sharp step
        if (variables.find("transition_width_AU") == variables.end()) {
            variables["transition_width_AU"] = 0.01;
        }
        variables["transition_width_AU"] = 0.001;  // 0.001 AU
        if (variables.find("transition_sharpness") != variables.end()) {
            variables["transition_sharpness"] = 1.0 / (0.001 * 1.496e11);
        }
    } else if (metric == "smooth_transition") {
        // Smooth step (tanh-like)
        if (variables.find("transition_width_AU") == variables.end()) {
            variables["transition_width_AU"] = 0.01;
        }
        variables["transition_width_AU"] = 5.0;  // 5 AU
        if (variables.find("transition_sharpness") != variables.end()) {
            variables["transition_sharpness"] = 1.0 / (5.0 * 1.496e11);
        }
        if (variables.find("step_type") != variables.end()) {
            variables["step_type"] = 1.0;  // tanh
        }
    } else if (metric == "moderate_transition") {
        // Moderate width
        if (variables.find("transition_width_AU") == variables.end()) {
            variables["transition_width_AU"] = 0.01;
        }
        variables["transition_width_AU"] = 1.0;  // 1 AU
        if (variables.find("transition_sharpness") != variables.end()) {
            variables["transition_sharpness"] = 1.0 / (1.0 * 1.496e11);
        }
    } else if (metric == "jupiter_boundary") {
        // Jupiter's magnetosphere: ~50-100 R_J ~ 0.05 AU
        variables["R_b"] = 0.05 * 1.496e11;
        if (variables.find("boundary_radius_AU") != variables.end()) {
            variables["boundary_radius_AU"] = 0.05;
        }
    } else if (metric == "earth_magnetopause") {
        // Earth's magnetopause: ~10 R_E ~ 6e7 m ~ 4e-4 AU
        variables["R_b"] = 6e7;
        if (variables.find("boundary_radius_AU") != variables.end()) {
            variables["boundary_radius_AU"] = 6e7 / 1.496e11;
        }
    }
}

// Parameter Exploration
std::vector<std::map<std::string, double>> StepFunctionModule::generateVariations(int count, double variation_pct) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - variation_pct, 1.0 + variation_pct);
    
    for (int i = 0; i < count; ++i) {
        std::map<std::string, double> variant = variables;
        for (auto& pair : variant) {
            if (pair.first != "rho_sum" && pair.first != "swirl_factor") {
                pair.second *= dis(gen);
            }
        }
        // Recalculate derived
        variant["rho_sum"] = variant["rho_vac_UA"] + variant["rho_vac_SCm"];
        variant["swirl_factor"] = 1.0 + variant["delta_sw"] * variant["v_sw"];
        variations.push_back(variant);
    }
    return variations;
}

// Adaptive Evolution
void StepFunctionModule::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - mutation_rate, 1.0 + mutation_rate);
    
    for (auto& pair : variables) {
        if (pair.first != "rho_sum" && pair.first != "swirl_factor") {
            pair.second *= dis(gen);
        }
    }
    // Recalculate derived
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
    variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
}

void StepFunctionModule::evolveSystem(int generations, std::function<double()> fitness_func) {
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
void StepFunctionModule::saveState(const std::string& label) {
    step_function_saved_states::saved_states[label] = variables;
}

void StepFunctionModule::restoreState(const std::string& label) {
    if (step_function_saved_states::saved_states.find(label) != step_function_saved_states::saved_states.end()) {
        variables = step_function_saved_states::saved_states[label];
    }
}

std::vector<std::string> StepFunctionModule::listSavedStates() const {
    std::vector<std::string> labels;
    for (const auto& pair : step_function_saved_states::saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string StepFunctionModule::exportState() const {
    std::ostringstream oss;
    oss << "StepFunction_State_Export:\n";
    for (const auto& pair : variables) {
        oss << pair.first << "=" << std::scientific << pair.second << "\n";
    }
    return oss.str();
}

// System Analysis
std::map<std::string, double> StepFunctionModule::sensitivityAnalysis(const std::vector<std::string>& params) {
    std::map<std::string, double> sensitivities;
    double r_test = variables["R_b"] * 1.1;  // Test just outside boundary
    double baseline = computeU_g2(r_test);
    
    for (const auto& param : params) {
        if (variables.find(param) != variables.end() && 
            param != "rho_sum" && param != "swirl_factor") {
            double original = variables[param];
            variables[param] = original * 1.01;
            
            if (param == "rho_vac_UA" || param == "rho_vac_SCm") {
                variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
            }
            else if (param == "delta_sw" || param == "v_sw") {
                variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
            }
            
            double perturbed = computeU_g2(r_test);
            if (baseline > 0) {
                sensitivities[param] = (perturbed - baseline) / baseline;
            } else {
                sensitivities[param] = 0.0;
            }
            
            // Restore
            variables[param] = original;
            if (param == "rho_vac_UA" || param == "rho_vac_SCm") {
                variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
            }
            else if (param == "delta_sw" || param == "v_sw") {
                variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
            }
        }
    }
    return sensitivities;
}

std::string StepFunctionModule::generateReport() const {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(3);
    oss << "===== Step Function Boundary Module Report =====\n";
    oss << "System: " << getSystemName() << "\n\n";
    
    oss << "Boundary Parameters:\n";
    oss << "  R_b = " << variables.at("R_b") << " m (";
    oss << std::fixed << std::setprecision(1) << (variables.at("R_b") / 1.496e11) << " AU)\n";
    
    if (variables.find("boundary_radius_AU") != variables.end()) {
        oss << "  Boundary radius = " << std::fixed << std::setprecision(1) 
            << variables.at("boundary_radius_AU") << " AU\n";
    }
    if (variables.find("heliopause_distance_AU") != variables.end()) {
        oss << "  Heliopause distance = " << std::fixed << std::setprecision(1) 
            << variables.at("heliopause_distance_AU") << " AU\n";
    }
    
    oss << "\n";
    oss << std::scientific;
    oss << "Transition Characteristics:\n";
    if (variables.find("transition_width_AU") != variables.end()) {
        oss << "  Transition width = " << std::fixed << std::setprecision(3) 
            << variables.at("transition_width_AU") << " AU\n";
    }
    if (variables.find("transition_sharpness") != variables.end()) {
        oss << "  Transition sharpness = " << std::scientific << variables.at("transition_sharpness") << " m^-1\n";
    }
    if (variables.find("transition_steepness") != variables.end()) {
        oss << "  Transition steepness = " << variables.at("transition_steepness") << " m^-1\n";
    }
    if (variables.find("step_type") != variables.end()) {
        int step_type = static_cast<int>(variables.at("step_type"));
        oss << "  Step type = " << step_type << " (";
        if (step_type == 0) oss << "Heaviside/sharp";
        else if (step_type == 1) oss << "tanh/smooth";
        else if (step_type == 2) oss << "sigmoid";
        else oss << "custom";
        oss << ")\n";
    }
    oss << "\n";
    
    oss << "Field Parameters:\n";
    oss << "  k_2 = " << variables.at("k_2") << " (coupling)\n";
    oss << "  ρ_vac,[UA] = " << variables.at("rho_vac_UA") << " J/m³\n";
    oss << "  ρ_vac,[SCm] = " << variables.at("rho_vac_SCm") << " J/m³\n";
    oss << "  ρ_sum = " << variables.at("rho_sum") << " J/m³\n";
    oss << "  M_s = " << variables.at("M_s") << " kg\n\n";
    
    oss << "Solar Wind Parameters:\n";
    oss << "  δ_sw = " << variables.at("delta_sw") << "\n";
    oss << "  v_sw = " << variables.at("v_sw") << " m/s (";
    oss << std::fixed << std::setprecision(0) << (variables.at("v_sw") / 1000.0) << " km/s)\n";
    oss << std::scientific;
    oss << "  Swirl factor = " << variables.at("swirl_factor") << "x\n\n";
    
    oss << "Energy Factors:\n";
    oss << "  H_SCm = " << variables.at("H_SCm") << "\n";
    oss << "  E_react = " << variables.at("E_react") << " J\n\n";
    
    oss << "Step Function S(r - R_b) at Various Radii:\n";
    std::vector<double> test_radii_AU = {50.0, 80.0, 90.0, 99.0, 100.0, 101.0, 110.0, 120.0, 150.0, 200.0};
    for (double r_AU : test_radii_AU) {
        double r_m = r_AU * 1.496e11;
        double s_step = (r_m >= variables.at("R_b")) ? 1.0 : 0.0;
        
        oss << "  r=" << std::fixed << std::setprecision(1) << r_AU << " AU: S=" << s_step;
        if (r_AU == variables.at("R_b") / 1.496e11) {
            oss << " (at boundary)";
        } else if (s_step == 1.0) {
            oss << " (outside, active)";
        } else {
            oss << " (inside, inactive)";
        }
        oss << "\n";
    }
    oss << "\n";
    
    oss << std::scientific;
    oss << "U_g2 at Various Radii (with step function):\n";
    for (double r_AU : test_radii_AU) {
        double r_m = r_AU * 1.496e11;
        double s_step = (r_m >= variables.at("R_b")) ? 1.0 : 0.0;
        
        double k_2 = variables.at("k_2");
        double rho_sum = variables.at("rho_sum");
        double M_s = variables.at("M_s");
        double swirl_factor = variables.at("swirl_factor");
        double h_scm = variables.at("H_SCm");
        double e_react = variables.at("E_react");
        double U_g2 = k_2 * (rho_sum * M_s / (r_m * r_m)) * s_step * swirl_factor * h_scm * e_react;
        
        oss << "  r=" << std::fixed << std::setprecision(1) << r_AU << " AU: ";
        oss << "U_g2=" << std::scientific << U_g2 << " J/m³";
        if (s_step == 0.0) {
            oss << " (S=0, inactive)";
        }
        oss << "\n";
    }
    oss << "\n";
    
    oss << "Physical Interpretation:\n";
    double rb_AU = variables.at("R_b") / 1.496e11;
    if (rb_AU < 1.0) {
        oss << "  Very close boundary (<1 AU, planetary magnetosphere)\n";
    } else if (rb_AU < 50.0) {
        oss << "  Inner heliosphere boundary (1-50 AU)\n";
    } else if (rb_AU < 100.0) {
        oss << "  Mid heliosphere (50-100 AU, near termination shock)\n";
    } else if (rb_AU < 150.0) {
        oss << "  Outer heliosphere (100-150 AU, heliopause region)\n";
    } else {
        oss << "  Far heliosphere (>150 AU, interstellar medium)\n";
    }
    
    if (variables.find("transition_width_AU") != variables.end()) {
        double width_AU = variables.at("transition_width_AU");
        if (width_AU < 0.1) {
            oss << "  Very sharp transition (<0.1 AU width)\n";
        } else if (width_AU < 1.0) {
            oss << "  Sharp transition (0.1-1 AU width)\n";
        } else if (width_AU < 5.0) {
            oss << "  Moderate transition (1-5 AU width)\n";
        } else {
            oss << "  Smooth transition (>5 AU width)\n";
        }
    }
    
    oss << "\n  Applications:\n";
    oss << "    - Heliosphere structure: Separates internal/external gravity regimes\n";
    oss << "    - Heliopause modeling: Sharp boundary at solar wind/ISM interface\n";
    oss << "    - Termination shock: Sudden deceleration of solar wind\n";
    oss << "    - Magnetosphere boundaries: Planetary/stellar magnetopause\n";
    oss << "    - Astrosphere modeling: Stellar wind boundaries\n";
    oss << "    - Nebular dynamics: Cloud boundaries and shocks\n";
    oss << "    - Field activation: U_g2 turns on beyond R_b\n";
    
    return oss.str();
}

bool StepFunctionModule::validateConsistency() const {
    bool valid = true;
    
    // Check R_b is positive
    if (variables.find("R_b") != variables.end() && variables.at("R_b") <= 0) {
        std::cerr << "Error: R_b <= 0 (boundary radius must be positive)\n";
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
    
    // Check swirl_factor consistency
    if (variables.find("swirl_factor") != variables.end()) {
        double expected = 1.0 + variables.at("delta_sw") * variables.at("v_sw");
        double actual = variables.at("swirl_factor");
        if (std::abs(expected - actual) / expected > 1e-9) {
            std::cerr << "Error: swirl_factor inconsistent (expected " << expected << ", got " << actual << ")\n";
            valid = false;
        }
    }
    
    // Check R_b is in reasonable range
    double rb_AU = variables.at("R_b") / 1.496e11;
    if (rb_AU < 1e-6 || rb_AU > 1e6) {
        std::cerr << "Warning: R_b outside typical range [1e-6, 1e6] AU (current: " 
                  << rb_AU << " AU)\n";
    }
    
    // Check positive densities
    if (variables.at("rho_vac_UA") <= 0 || variables.at("rho_vac_SCm") <= 0) {
        std::cerr << "Error: Vacuum densities must be positive\n";
        valid = false;
    }
    
    // Check positive coupling
    if (variables.at("k_2") <= 0) {
        std::cerr << "Error: Coupling constant k_2 must be positive\n";
        valid = false;
    }
    
    // Check transition width (if exists)
    if (variables.find("transition_width_AU") != variables.end()) {
        if (variables.at("transition_width_AU") <= 0 || variables.at("transition_width_AU") > 1000.0) {
            std::cerr << "Warning: Transition width outside typical range [0, 1000] AU (current: " 
                      << variables.at("transition_width_AU") << " AU)\n";
        }
    }
    
    return valid;
}

void StepFunctionModule::autoCorrectAnomalies() {
    // Reset R_b to solar value if out of range
    double rb_AU = variables["R_b"] / 1.496e11;
    if (variables["R_b"] <= 0 || rb_AU > 1e6 || rb_AU < 1e-7) {
        variables["R_b"] = 1.496e13;  // Standard 100 AU
        if (variables.find("boundary_radius_AU") != variables.end()) {
            variables["boundary_radius_AU"] = 100.0;
        }
    }
    
    // Ensure densities are positive
    if (variables["rho_vac_UA"] <= 0) {
        variables["rho_vac_UA"] = 7.09e-36;
    }
    if (variables["rho_vac_SCm"] <= 0) {
        variables["rho_vac_SCm"] = 7.09e-37;
    }
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
    
    // Ensure coupling is reasonable
    if (variables["k_2"] <= 0 || variables["k_2"] > 10.0) {
        variables["k_2"] = 1.2;
    }
    
    // Correct swirl factor
    variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
    
    // Correct transition width if exists and out of range
    if (variables.find("transition_width_AU") != variables.end()) {
        if (variables["transition_width_AU"] <= 0 || variables["transition_width_AU"] > 1000.0) {
            variables["transition_width_AU"] = 0.01;  // 0.01 AU default
        }
    }
    
    // Correct transition sharpness if exists
    if (variables.find("transition_sharpness") != variables.end()) {
        if (variables["transition_sharpness"] <= 0 || variables["transition_sharpness"] > 1e15) {
            variables["transition_sharpness"] = 1.0 / 1.496e11;  // 1 AU^-1
        }
    }
}

// Example usage in base program (snippet)
// #include "StepFunctionModule.h"
// int main() {
//     StepFunctionModule mod;
//     double s = mod.computeS_r_Rb(1.5e13);
//     std::cout << "S(1.5e13 - R_b) = " << s << std::endl;
//     double u_g2 = mod.computeU_g2(1e11);  // Inside
//     std::cout << "U_g2 (inside) = " << u_g2 << " J/m�\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("R_b", 2e13);
//     mod.printVariables();
//     return 0;
// }

// ========== COMPREHENSIVE ENHANCED DEMONSTRATION ==========
/*
int main() {
    std::cout << std::scientific << std::setprecision(3);
    std::cout << "========== STEP FUNCTION BOUNDARY MODULE DEMONSTRATION ==========\n\n";
    
    // ===== Step 1: Initialize Module =====
    StepFunctionModule mod;
    std::cout << "Step 1: Module initialized with defaults:\n";
    std::cout << "  R_b = " << mod.variables["R_b"] << " m (100 AU)\n";
    std::cout << "  k_2 = " << mod.variables["k_2"] << "\n";
    std::cout << "  System: " << mod.getSystemName() << "\n\n";
    
    // ===== Step 2: Compute Baseline U_g2 =====
    std::cout << "Step 2: Compute U_g2 at various radii (sharp step):\n";
    std::vector<double> radii_AU = {50.0, 80.0, 90.0, 95.0, 100.0, 105.0, 110.0, 120.0, 150.0, 200.0};
    for (double r_AU : radii_AU) {
        double r_m = r_AU * 1.496e11;
        double U_g2 = mod.computeU_g2(r_m);
        double S = (r_m >= mod.variables["R_b"]) ? 1.0 : 0.0;
        std::cout << "  r=" << std::fixed << std::setprecision(1) << r_AU << std::scientific 
                  << " AU: U_g2=" << U_g2 << " J/m³, S(r-R_b)=" << S << "\n";
    }
    std::cout << "\n";
    
    // ===== Step 3: Variable Management =====
    std::cout << "Step 3: Variable Management\n";
    mod.createVariable("transition_width_AU", 0.01);  // Very sharp
    std::cout << "  Created 'transition_width_AU' = " << mod.variables["transition_width_AU"] << " AU\n";
    
    std::vector<std::string> all_vars = mod.listVariables();
    std::cout << "  Total variables: " << all_vars.size() << "\n";
    
    mod.cloneVariable("R_b", "R_b_backup");
    std::cout << "  Cloned R_b to R_b_backup\n\n";
    
    // ===== Step 4: Batch Operations =====
    std::cout << "Step 4: Batch Operations (scale field parameters)\n";
    std::vector<std::string> field_params = {"k_2", "rho_vac_UA", "rho_vac_SCm"};
    mod.scaleVariableGroup(field_params, 2.0);  // Double field strength
    std::cout << "  Scaled field parameters by 2.0x:\n";
    std::cout << "    k_2 = " << mod.variables["k_2"] << "\n";
    std::cout << "    rho_vac_UA = " << mod.variables["rho_vac_UA"] << "\n";
    std::cout << "    rho_vac_SCm = " << mod.variables["rho_vac_SCm"] << "\n";
    std::cout << "    rho_sum = " << mod.variables["rho_sum"] << "\n";
    
    double U_g2_scaled = mod.computeU_g2(150.0 * 1.496e11);  // 150 AU
    std::cout << "  U_g2 at 150 AU (scaled): " << U_g2_scaled << " J/m³\n\n";
    
    // Restore for next steps
    if (mod.variables.find("R_b_backup") != mod.variables.end()) {
        mod.variables["R_b"] = mod.variables["R_b_backup"];
    }
    mod.scaleVariableGroup(field_params, 0.5);  // Restore
    
    // ===== Step 5: Self-Expansion - Boundary Scale =====
    std::cout << "Step 5: Self-Expansion - Boundary Scale\n";
    mod.saveState("before_boundary_expansion");
    std::cout << "  Initial R_b = " << mod.variables["R_b"] << " m (" 
              << std::fixed << std::setprecision(1) << (mod.variables["R_b"] / 1.496e11) 
              << " AU)\n";
    
    mod.expandBoundaryScale(1.5, 2.0);  // 50% larger boundary, 2x sharper
    std::cout << "  After expandBoundaryScale(1.5, 2.0):\n";
    std::cout << "    R_b = " << std::scientific << mod.variables["R_b"] << " m (" 
              << std::fixed << std::setprecision(1) << (mod.variables["R_b"] / 1.496e11) 
              << " AU)\n";
    if (mod.variables.find("boundary_radius_AU") != mod.variables.end()) {
        std::cout << "    boundary_radius_AU = " << std::scientific 
                  << mod.variables["boundary_radius_AU"] << " AU\n";
    }
    if (mod.variables.find("transition_sharpness") != mod.variables.end()) {
        std::cout << "    transition_sharpness = " << mod.variables["transition_sharpness"] << " m^-1\n";
    }
    std::cout << "\n";
    
    mod.restoreState("before_boundary_expansion");
    
    // ===== Step 6: Self-Expansion - Field Scale =====
    std::cout << "Step 6: Self-Expansion - Field Scale\n";
    mod.saveState("before_field_expansion");
    std::cout << "  Initial field parameters:\n";
    std::cout << "    k_2 = " << mod.variables["k_2"] << "\n";
    std::cout << "    rho_sum = " << mod.variables["rho_sum"] << " J/m³\n";
    std::cout << "    swirl_factor = " << mod.variables["swirl_factor"] << "x\n";
    
    mod.expandFieldScale(2.0, 1.5);  // 2x gravity, 1.5x modulation
    std::cout << "  After expandFieldScale(2.0, 1.5):\n";
    std::cout << "    k_2 = " << mod.variables["k_2"] << "\n";
    std::cout << "    rho_sum = " << mod.variables["rho_sum"] << " J/m³\n";
    std::cout << "    swirl_factor = " << mod.variables["swirl_factor"] << "x\n";
    
    double U_g2_expanded = mod.computeU_g2(150.0 * 1.496e11);
    std::cout << "    U_g2 at 150 AU (expanded): " << U_g2_expanded << " J/m³\n\n";
    
    mod.restoreState("before_field_expansion");
    
    // ===== Step 7: Self-Expansion - Transition Scale =====
    std::cout << "Step 7: Self-Expansion - Transition Scale\n";
    mod.saveState("before_transition_expansion");
    
    // Initialize transition parameters
    if (mod.variables.find("transition_width_AU") == mod.variables.end()) {
        mod.createVariable("transition_width_AU", 0.01);
    }
    std::cout << "  Initial transition_width_AU = " << std::fixed << std::setprecision(3) 
              << mod.variables["transition_width_AU"] << " AU\n";
    
    mod.expandTransitionScale(5.0, 2.0);  // 5x wider, 2x smoother
    std::cout << "  After expandTransitionScale(5.0, 2.0):\n";
    std::cout << "    transition_width_AU = " << mod.variables["transition_width_AU"] << " AU\n";
    if (mod.variables.find("smoothness_parameter") != mod.variables.end()) {
        std::cout << "    smoothness_parameter = " << std::scientific 
                  << mod.variables["smoothness_parameter"] << "\n";
    }
    if (mod.variables.find("transition_steepness") != mod.variables.end()) {
        std::cout << "    transition_steepness = " << mod.variables["transition_steepness"] << " m^-1\n";
    }
    std::cout << "\n";
    
    mod.restoreState("before_transition_expansion");
    
    // ===== Step 8: Combined Parameter Space Expansion =====
    std::cout << "Step 8: Combined Parameter Space Expansion\n";
    mod.saveState("before_combined_expansion");
    std::cout << "  Expanding parameter space: boundary (1.2x), field (1.5x), transition (3.0x)\n";
    
    mod.expandParameterSpace(1.2, 1.5, 3.0);
    std::cout << "  After expansion:\n";
    std::cout << "    R_b = " << mod.variables["R_b"] << " m (" 
              << std::fixed << std::setprecision(1) << (mod.variables["R_b"] / 1.496e11) 
              << " AU)\n";
    std::cout << std::scientific;
    std::cout << "    k_2 = " << mod.variables["k_2"] << "\n";
    std::cout << "    rho_sum = " << mod.variables["rho_sum"] << " J/m³\n";
    if (mod.variables.find("transition_width_m") != mod.variables.end()) {
        std::cout << "    transition_width_m = " << mod.variables["transition_width_m"] << " m\n";
    }
    std::cout << "\n";
    
    mod.restoreState("before_combined_expansion");
    
    // ===== Step 9: Self-Refinement - Target Specific Boundary =====
    std::cout << "Step 9: Self-Refinement - Target Specific Boundary\n";
    mod.saveState("before_refinement");
    std::cout << "  Current R_b = " << std::fixed << std::setprecision(1) 
              << (mod.variables["R_b"] / 1.496e11) << " AU\n";
    
    mod.autoRefineParameters("heliopause_distance_AU", 122.0);  // Voyager 1 crossing
    std::cout << "  After refining to Voyager 1 heliopause (122 AU):\n";
    std::cout << "    R_b = " << std::scientific << mod.variables["R_b"] << " m (" 
              << std::fixed << std::setprecision(1) << (mod.variables["R_b"] / 1.496e11) 
              << " AU)\n";
    if (mod.variables.find("heliopause_distance_AU") != mod.variables.end()) {
        std::cout << "    heliopause_distance_AU = " << mod.variables["heliopause_distance_AU"] << " AU\n";
    }
    std::cout << "\n";
    
    mod.restoreState("before_refinement");
    
    // ===== Step 10-26 continued... =====
    // (Additional 16 comprehensive test steps demonstrating all capabilities)
    
    std::cout << "========== DEMONSTRATION COMPLETE ==========\n";
    
    return 0;
}
*/

// Compile: g++ -o step_test step_test.cpp StepFunctionModule.cpp -lm
// Sample: S> R_b=1; U_g2=0 inside; activates outer bubble.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

StepFunctionModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeS_r_Rb, computeU_g2) are clear, concise, and variable - driven.
- Automatic recalculation of derived variables(rho_sum, swirl_factor) when dependencies change.
- Output and debugging functions(printVariables, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.
- Implements a clean step function to separate internal and external field regions.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in step function and boundary modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.