// SolarWindModulationModule.h
// Modular C++ implementation of the Solar Wind Modulation Factor (?_sw) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes ?_sw=0.01 (unitless) and its scaling (1 + ?_sw v_sw) in Universal Gravity U_g2 term.
// Pluggable: #include "SolarWindModulationModule.h"
// SolarWindModulationModule mod; mod.computeU_g2(1.496e13); mod.updateVariable("delta_sw", new_value);
// Variables in std::map; example for Sun at r=R_b=1.496e13 m; amplification ~5001x.
// Approximations: S(r - R_b)=1; H_SCm=1; E_react=1e46; ?_sum=7.80e-36 J/m�.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef SOLAR_WIND_MODULATION_MODULE_H
#define SOLAR_WIND_MODULATION_MODULE_H

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

class SolarWindModulationModule {
private:
    std::map<std::string, double> variables;
    double computeModulationFactor();
    double computeU_g2(double r);

public:
    // Constructor: Initialize with framework defaults (Sun)
    SolarWindModulationModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeDelta_sw();  // 0.01 (unitless)
    double computeModulationFactor();  // 1 + ?_sw v_sw
    double computeU_g2(double r);  // U_g2 with modulation (J/m^3)
    double computeU_g2_no_mod(double r);  // Without modulation

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
    void expandParameterSpace(double modulation_scale, double gravity_scale, double velocity_scale);
    void expandModulationScale(double delta_factor, double amplification_factor);  // δ_sw and modulation
    void expandVelocityScale(double vsw_factor, double range_factor);              // v_sw and spatial scales
    void expandGravityScale(double k2_factor, double density_factor);              // k_2 and ρ_vac

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

#endif // SOLAR_WIND_MODULATION_MODULE_H

// SolarWindModulationModule.cpp
#include "SolarWindModulationModule.h"

// Constructor: Set framework defaults (Sun at r=R_b)
SolarWindModulationModule::SolarWindModulationModule() {
    // Universal constants
    variables["delta_sw"] = 0.01;                   // Unitless
    variables["v_sw"] = 5e5;                        // m/s
    variables["k_2"] = 1.2;                         // Coupling
    variables["rho_vac_UA"] = 7.09e-36;             // J/m^3
    variables["rho_vac_SCm"] = 7.09e-37;            // J/m^3
    variables["M_s"] = 1.989e30;                    // kg
    variables["r"] = 1.496e13;                      // m (R_b)
    variables["R_b"] = 1.496e13;                    // m
    variables["S_r_Rb"] = 1.0;                      // Step
    variables["H_SCm"] = 1.0;                       // Unitless
    variables["E_react"] = 1e46;                    // J

    // Derived
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
    variables["modulation_factor"] = computeModulationFactor();
}

// Update variable
void SolarWindModulationModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
        if (name == "delta_sw" || name == "v_sw") {
            variables["modulation_factor"] = computeModulationFactor();
        }
        else if (name == "rho_vac_UA" || name == "rho_vac_SCm") {
            variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
        }
    }
    else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void SolarWindModulationModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        if (name == "delta_sw" || name == "v_sw") {
            variables["modulation_factor"] = computeModulationFactor();
        }
        else if (name == "rho_vac_UA" || name == "rho_vac_SCm") {
            variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
        }
    }
    else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void SolarWindModulationModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute ?_sw = 0.01
double SolarWindModulationModule::computeDelta_sw() {
    return variables["delta_sw"];
}

// Compute 1 + ?_sw * v_sw
double SolarWindModulationModule::computeModulationFactor() {
    return 1.0 + variables["delta_sw"] * variables["v_sw"];
}

// Compute U_g2 with modulation
double SolarWindModulationModule::computeU_g2(double r) {
    variables["r"] = r;
    double k_2 = variables["k_2"];
    double rho_sum = variables["rho_sum"];
    double M_s = variables["M_s"];
    double s_step = (r >= variables["R_b"]) ? 1.0 : 0.0;
    double mod_factor = computeModulationFactor();
    double h_scm = variables["H_SCm"];
    double e_react = variables["E_react"];
    return k_2 * (rho_sum * M_s / (r * r)) * s_step * mod_factor * h_scm * e_react;
}

// U_g2 without modulation (?_sw=0)
double SolarWindModulationModule::computeU_g2_no_mod(double r) {
    double orig_delta = variables["delta_sw"];
    variables["delta_sw"] = 0.0;
    double result = computeU_g2(r);
    variables["delta_sw"] = orig_delta;
    return result;
}

// Equation text
std::string SolarWindModulationModule::getEquationText() {
    return "U_g2 = k_2 * [(?_vac,[UA] + ?_vac,[SCm]) M_s / r^2] * S(r - R_b) * (1 + ?_sw v_sw) * H_SCm * E_react\n"
        "Where ?_sw = 0.01 (unitless solar wind modulation factor);\n"
        "Modulation = 1 + 0.01 * v_sw (v_sw=5e5 m/s ? ~5001x amplification).\n"
        "Example r=R_b=1.496e13 m: U_g2 ?1.18e53 J/m� (with); ?2.36e49 J/m� (without; ~5000x less).\n"
        "Role: Enhances external gravity via solar wind momentum/pressure beyond R_b.\n"
        "UQFF: Models heliosphere dynamics; wind influence on nebular/star formation.";
}

// Print variables
void SolarWindModulationModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// ===== ENHANCED METHODS IMPLEMENTATION =====

namespace solar_wind_modulation_saved_states {
    static std::map<std::string, std::map<std::string, double>> saved_states;
}

// Variable Management
void SolarWindModulationModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void SolarWindModulationModule::removeVariable(const std::string& name) {
    variables.erase(name);
}

void SolarWindModulationModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
    }
}

std::vector<std::string> SolarWindModulationModule::listVariables() const {
    std::vector<std::string> names;
    for (const auto& pair : variables) {
        names.push_back(pair.first);
    }
    return names;
}

std::string SolarWindModulationModule::getSystemName() const {
    return "Solar_Wind_Modulation_UQFF";
}

// Batch Operations
void SolarWindModulationModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = func(variables[name]);
        }
    }
    // Update derived
    variables["modulation_factor"] = computeModulationFactor();
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
}

void SolarWindModulationModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// Self-Expansion: Domain-Specific Scales
void SolarWindModulationModule::expandParameterSpace(double modulation_scale, double gravity_scale, double velocity_scale) {
    variables["delta_sw"] *= modulation_scale;
    variables["v_sw"] *= velocity_scale;
    variables["k_2"] *= gravity_scale;
    variables["rho_vac_UA"] *= gravity_scale;
    variables["rho_vac_SCm"] *= gravity_scale;
    
    // Update derived
    variables["modulation_factor"] = computeModulationFactor();
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
}

void SolarWindModulationModule::expandModulationScale(double delta_factor, double amplification_factor) {
    variables["delta_sw"] *= delta_factor;
    
    // Amplification factor adjusts the overall modulation effect
    // (1 + δ_sw v_sw) can be scaled by adjusting v_sw or δ_sw
    // Here we adjust both proportionally
    double current_mod = computeModulationFactor();
    double target_mod = current_mod * amplification_factor;
    
    // Solve for new δ_sw: target = 1 + δ_new * v_sw
    // δ_new = (target - 1) / v_sw
    if (variables["v_sw"] > 0) {
        variables["delta_sw"] = (target_mod - 1.0) / variables["v_sw"];
    }
    
    variables["modulation_factor"] = computeModulationFactor();
}

void SolarWindModulationModule::expandVelocityScale(double vsw_factor, double range_factor) {
    variables["v_sw"] *= vsw_factor;
    variables["R_b"] *= range_factor;  // Scale heliopause boundary
    variables["r"] *= range_factor;    // Scale current radius
    
    variables["modulation_factor"] = computeModulationFactor();
}

void SolarWindModulationModule::expandGravityScale(double k2_factor, double density_factor) {
    variables["k_2"] *= k2_factor;
    variables["rho_vac_UA"] *= density_factor;
    variables["rho_vac_SCm"] *= density_factor;
    
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
}

// Self-Refinement
void SolarWindModulationModule::autoRefineParameters(const std::string& target, double goal) {
    if (target == "delta_sw") {
        variables["delta_sw"] = goal;
        variables["modulation_factor"] = computeModulationFactor();
    } else if (target == "modulation_factor") {
        // Solve: goal = 1 + δ_sw * v_sw
        // δ_sw = (goal - 1) / v_sw
        if (variables["v_sw"] > 0) {
            variables["delta_sw"] = (goal - 1.0) / variables["v_sw"];
            variables["modulation_factor"] = computeModulationFactor();
        }
    } else if (target == "amplification") {
        // Same as modulation_factor
        if (variables["v_sw"] > 0) {
            variables["delta_sw"] = (goal - 1.0) / variables["v_sw"];
            variables["modulation_factor"] = computeModulationFactor();
        }
    } else if (target == "v_sw") {
        variables["v_sw"] = goal;
        variables["modulation_factor"] = computeModulationFactor();
    } else if (target == "v_sw_km_s") {
        // Convert km/s to m/s
        variables["v_sw"] = goal * 1000.0;
        variables["modulation_factor"] = computeModulationFactor();
    } else if (target == "k_2") {
        variables["k_2"] = goal;
    } else if (target == "R_b_AU") {
        // Convert AU to meters: 1 AU ≈ 1.496e11 m
        variables["R_b"] = goal * 1.496e11;
    } else if (target == "U_g2_at_Rb") {
        // Target specific U_g2 at R_b by adjusting k_2
        double current_U_g2 = computeU_g2(variables["R_b"]);
        if (current_U_g2 > 0) {
            variables["k_2"] *= (goal / current_U_g2);
        }
    }
}

void SolarWindModulationModule::calibrateToObservations(const std::map<std::string, double>& observations) {
    for (const auto& obs : observations) {
        if (variables.find(obs.first) != variables.end()) {
            variables[obs.first] = obs.second;
            if (obs.first == "delta_sw" || obs.first == "v_sw") {
                variables["modulation_factor"] = computeModulationFactor();
            }
            else if (obs.first == "rho_vac_UA" || obs.first == "rho_vac_SCm") {
                variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
            }
        }
    }
}

void SolarWindModulationModule::optimizeForMetric(const std::string& metric) {
    if (metric == "high_modulation") {
        // Increase modulation (higher δ_sw)
        variables["delta_sw"] = 0.02;  // Double standard
        variables["modulation_factor"] = computeModulationFactor();
    } else if (metric == "low_modulation") {
        // Decrease modulation
        variables["delta_sw"] = 0.005;  // Half standard
        variables["modulation_factor"] = computeModulationFactor();
    } else if (metric == "standard_modulation") {
        // Reset to standard
        variables["delta_sw"] = 0.01;
        variables["modulation_factor"] = computeModulationFactor();
    } else if (metric == "fast_wind") {
        // Fast solar wind: 700-800 km/s
        variables["v_sw"] = 7.5e5;  // 750 km/s
        variables["modulation_factor"] = computeModulationFactor();
    } else if (metric == "slow_wind") {
        // Slow solar wind: 300-400 km/s
        variables["v_sw"] = 3.5e5;  // 350 km/s
        variables["modulation_factor"] = computeModulationFactor();
    } else if (metric == "standard_wind") {
        // Standard: 500 km/s
        variables["v_sw"] = 5.0e5;
        variables["modulation_factor"] = computeModulationFactor();
    } else if (metric == "extreme_amplification") {
        // Very high amplification: δ_sw = 0.05, v_sw = 1e6 m/s
        variables["delta_sw"] = 0.05;
        variables["v_sw"] = 1.0e6;
        variables["modulation_factor"] = computeModulationFactor();
    } else if (metric == "minimal_amplification") {
        // Minimal: δ_sw = 0.001, v_sw = 2e5 m/s
        variables["delta_sw"] = 0.001;
        variables["v_sw"] = 2.0e5;
        variables["modulation_factor"] = computeModulationFactor();
    }
}

// Parameter Exploration
std::vector<std::map<std::string, double>> SolarWindModulationModule::generateVariations(int count, double variation_pct) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - variation_pct, 1.0 + variation_pct);
    
    for (int i = 0; i < count; ++i) {
        std::map<std::string, double> variant = variables;
        for (auto& pair : variant) {
            if (pair.first != "modulation_factor" && pair.first != "rho_sum" && 
                pair.first != "S_r_Rb") {
                pair.second *= dis(gen);
            }
        }
        // Recalculate derived
        variant["modulation_factor"] = 1.0 + variant["delta_sw"] * variant["v_sw"];
        variant["rho_sum"] = variant["rho_vac_UA"] + variant["rho_vac_SCm"];
        variations.push_back(variant);
    }
    return variations;
}

// Adaptive Evolution
void SolarWindModulationModule::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - mutation_rate, 1.0 + mutation_rate);
    
    for (auto& pair : variables) {
        if (pair.first != "modulation_factor" && pair.first != "rho_sum" && 
            pair.first != "S_r_Rb") {
            pair.second *= dis(gen);
        }
    }
    // Recalculate derived
    variables["modulation_factor"] = computeModulationFactor();
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
}

void SolarWindModulationModule::evolveSystem(int generations, std::function<double()> fitness_func) {
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
void SolarWindModulationModule::saveState(const std::string& label) {
    solar_wind_modulation_saved_states::saved_states[label] = variables;
}

void SolarWindModulationModule::restoreState(const std::string& label) {
    if (solar_wind_modulation_saved_states::saved_states.find(label) != solar_wind_modulation_saved_states::saved_states.end()) {
        variables = solar_wind_modulation_saved_states::saved_states[label];
    }
}

std::vector<std::string> SolarWindModulationModule::listSavedStates() const {
    std::vector<std::string> labels;
    for (const auto& pair : solar_wind_modulation_saved_states::saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string SolarWindModulationModule::exportState() const {
    std::ostringstream oss;
    oss << "SolarWindModulation_State_Export:\n";
    for (const auto& pair : variables) {
        oss << pair.first << "=" << std::scientific << pair.second << "\n";
    }
    return oss.str();
}

// System Analysis
std::map<std::string, double> SolarWindModulationModule::sensitivityAnalysis(const std::vector<std::string>& params) {
    std::map<std::string, double> sensitivities;
    double r_test = variables["R_b"];  // Test at heliopause boundary
    double baseline = computeU_g2(r_test);
    
    for (const auto& param : params) {
        if (variables.find(param) != variables.end() && 
            param != "modulation_factor" && param != "rho_sum" && param != "S_r_Rb") {
            double original = variables[param];
            variables[param] = original * 1.01;
            
            if (param == "delta_sw" || param == "v_sw") {
                variables["modulation_factor"] = computeModulationFactor();
            }
            else if (param == "rho_vac_UA" || param == "rho_vac_SCm") {
                variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
            }
            
            double perturbed = computeU_g2(r_test);
            sensitivities[param] = (perturbed - baseline) / baseline;
            
            // Restore
            variables[param] = original;
            if (param == "delta_sw" || param == "v_sw") {
                variables["modulation_factor"] = computeModulationFactor();
            }
            else if (param == "rho_vac_UA" || param == "rho_vac_SCm") {
                variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
            }
        }
    }
    return sensitivities;
}

std::string SolarWindModulationModule::generateReport() const {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(3);
    oss << "===== Solar Wind Modulation Module Report =====\n";
    oss << "System: " << getSystemName() << "\n\n";
    
    oss << "Modulation Parameters:\n";
    oss << "  δ_sw = " << variables.at("delta_sw") << " (unitless)\n";
    oss << "  v_sw = " << variables.at("v_sw") << " m/s (";
    oss << std::fixed << std::setprecision(0) << (variables.at("v_sw") / 1000.0) << " km/s)\n";
    oss << std::scientific;
    oss << "  Modulation factor = " << variables.at("modulation_factor") << "x\n";
    double amplification = variables.at("modulation_factor");
    oss << "  Amplification = " << std::fixed << std::setprecision(0) 
        << amplification << "x (vs no modulation)\n\n";
    
    oss << std::scientific;
    oss << "Gravity Parameters:\n";
    oss << "  k_2 = " << variables.at("k_2") << "\n";
    oss << "  ρ_vac,[UA] = " << variables.at("rho_vac_UA") << " J/m³\n";
    oss << "  ρ_vac,[SCm] = " << variables.at("rho_vac_SCm") << " J/m³\n";
    oss << "  ρ_sum = " << variables.at("rho_sum") << " J/m³\n";
    oss << "  M_s = " << variables.at("M_s") << " kg\n\n";
    
    oss << "Spatial Parameters:\n";
    oss << "  R_b = " << variables.at("R_b") << " m (";
    oss << std::fixed << std::setprecision(0) << (variables.at("R_b") / 1.496e11) << " AU)\n";
    oss << std::scientific;
    oss << "  Current r = " << variables.at("r") << " m (";
    oss << std::fixed << std::setprecision(0) << (variables.at("r") / 1.496e11) << " AU)\n\n";
    
    oss << std::scientific;
    oss << "Energy Factors:\n";
    oss << "  H_SCm = " << variables.at("H_SCm") << "\n";
    oss << "  E_react = " << variables.at("E_react") << " J\n\n";
    
    oss << "U_g2 Computations at Key Radii:\n";
    std::vector<double> test_radii_AU = {50.0, 100.0, 150.0, 200.0, 500.0, 1000.0};
    for (double r_AU : test_radii_AU) {
        double r_m = r_AU * 1.496e11;
        double s_step = (r_m >= variables.at("R_b")) ? 1.0 : 0.0;
        
        if (s_step > 0) {
            // Save current r
            double orig_r = variables.at("r");
            const_cast<std::map<std::string, double>&>(variables)["r"] = r_m;
            
            double k_2 = variables.at("k_2");
            double rho_sum = variables.at("rho_sum");
            double M_s = variables.at("M_s");
            double mod_factor = variables.at("modulation_factor");
            double h_scm = variables.at("H_SCm");
            double e_react = variables.at("E_react");
            double U_g2 = k_2 * (rho_sum * M_s / (r_m * r_m)) * s_step * mod_factor * h_scm * e_react;
            
            // Without modulation
            double U_g2_no_mod = k_2 * (rho_sum * M_s / (r_m * r_m)) * s_step * 1.0 * h_scm * e_react;
            
            oss << "  r=" << std::fixed << std::setprecision(0) << r_AU << " AU: ";
            oss << "U_g2=" << std::scientific << U_g2 << " J/m³ (with mod), ";
            oss << U_g2_no_mod << " J/m³ (no mod)\n";
            
            // Restore r
            const_cast<std::map<std::string, double>&>(variables)["r"] = orig_r;
        } else {
            oss << "  r=" << std::fixed << std::setprecision(0) << r_AU 
                << " AU: Below R_b, S=0, U_g2=0\n";
        }
    }
    oss << "\n";
    
    oss << "Physical Interpretation:\n";
    if (amplification > 10000.0) {
        oss << "  Extreme amplification (>10000x)\n";
    } else if (amplification > 5000.0) {
        oss << "  Very high amplification (5000-10000x)\n";
    } else if (amplification > 1000.0) {
        oss << "  High amplification (1000-5000x)\n";
    } else if (amplification > 100.0) {
        oss << "  Moderate amplification (100-1000x)\n";
    } else {
        oss << "  Low amplification (<100x)\n";
    }
    
    double v_sw_km = variables.at("v_sw") / 1000.0;
    if (v_sw_km > 700.0) {
        oss << "  Fast solar wind regime (>700 km/s)\n";
    } else if (v_sw_km > 400.0) {
        oss << "  Standard solar wind (400-700 km/s)\n";
    } else {
        oss << "  Slow solar wind (<400 km/s)\n";
    }
    
    oss << "\n  Applications:\n";
    oss << "    - Heliosphere dynamics: Solar wind pressure on outer bubble\n";
    oss << "    - External gravity enhancement: ~5000x amplification beyond R_b\n";
    oss << "    - Nebular interactions: Wind influence on surrounding ISM\n";
    oss << "    - Star formation: Pressure-driven collapse in molecular clouds\n";
    oss << "    - Planetary system evolution: Long-term wind effects\n";
    oss << "    - Astrosphere modeling: Stellar wind interactions\n";
    
    return oss.str();
}

bool SolarWindModulationModule::validateConsistency() const {
    bool valid = true;
    
    // Check δ_sw is non-negative
    if (variables.find("delta_sw") != variables.end() && variables.at("delta_sw") < 0) {
        std::cerr << "Error: delta_sw < 0 (modulation factor must be non-negative)\n";
        valid = false;
    }
    
    // Check v_sw is positive
    if (variables.find("v_sw") != variables.end() && variables.at("v_sw") <= 0) {
        std::cerr << "Error: v_sw <= 0 (solar wind velocity must be positive)\n";
        valid = false;
    }
    
    // Check modulation factor consistency
    if (variables.find("modulation_factor") != variables.end()) {
        double expected = 1.0 + variables.at("delta_sw") * variables.at("v_sw");
        double actual = variables.at("modulation_factor");
        if (std::abs(expected - actual) / expected > 1e-9) {
            std::cerr << "Error: modulation_factor inconsistent (expected " << expected << ", got " << actual << ")\n";
            valid = false;
        }
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
    
    // Check v_sw is in reasonable range (100-2000 km/s)
    double v_sw_km = variables.at("v_sw") / 1000.0;
    if (v_sw_km < 100.0 || v_sw_km > 2000.0) {
        std::cerr << "Warning: v_sw outside typical range [100, 2000] km/s (current: " 
                  << v_sw_km << " km/s)\n";
    }
    
    // Check δ_sw is reasonable (<0.1)
    if (variables.at("delta_sw") > 0.1) {
        std::cerr << "Warning: delta_sw > 0.1 (very high modulation factor: " 
                  << variables.at("delta_sw") << ")\n";
    }
    
    // Check positive densities
    if (variables.at("rho_vac_UA") <= 0 || variables.at("rho_vac_SCm") <= 0) {
        std::cerr << "Error: Vacuum densities must be positive\n";
        valid = false;
    }
    
    return valid;
}

void SolarWindModulationModule::autoCorrectAnomalies() {
    // Reset δ_sw to typical value if out of range
    if (variables["delta_sw"] < 0 || variables["delta_sw"] > 0.1) {
        variables["delta_sw"] = 0.01;  // Standard value
        variables["modulation_factor"] = computeModulationFactor();
    }
    
    // Reset v_sw to typical value if out of range
    double v_sw_km = variables["v_sw"] / 1000.0;
    if (variables["v_sw"] <= 0 || v_sw_km > 2000.0 || v_sw_km < 100.0) {
        variables["v_sw"] = 5.0e5;  // Standard 500 km/s
        variables["modulation_factor"] = computeModulationFactor();
    }
    
    // Ensure densities are positive
    if (variables["rho_vac_UA"] <= 0) {
        variables["rho_vac_UA"] = 7.09e-36;
    }
    if (variables["rho_vac_SCm"] <= 0) {
        variables["rho_vac_SCm"] = 7.09e-37;
    }
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
    
    // Ensure k_2 is reasonable
    if (variables["k_2"] <= 0 || variables["k_2"] > 10.0) {
        variables["k_2"] = 1.2;
    }
}

// Example usage in base program (snippet)
int main() {
    SolarWindModulationModule module;
    std::cout << "===== Solar Wind Modulation Module Enhanced Demonstration =====\n\n";
    
    // Step 1: Report initial state
    std::cout << "STEP 1: Initial Configuration (δ_sw=0.01, v_sw=500 km/s)\n";
    std::cout << module.generateReport() << "\n";
    
    // Step 2: Compute U_g2 at various radii
    std::cout << "STEP 2: U_g2 at Various Radii (with and without modulation)\n";
    std::vector<double> test_radii = {50.0, 100.0, 150.0, 200.0, 500.0};
    for (double r_AU : test_radii) {
        double r_m = r_AU * 1.496e11;
        double U_g2_with = module.computeU_g2(r_m);
        double U_g2_without = module.computeU_g2_no_mod(r_m);
        if (U_g2_with > 0) {
            double ratio = U_g2_with / U_g2_without;
            std::cout << "  r=" << std::fixed << std::setprecision(0) << r_AU << " AU: ";
            std::cout << "U_g2=" << std::scientific << U_g2_with << " J/m³ (with), ";
            std::cout << U_g2_without << " J/m³ (without), ratio=" << std::fixed << std::setprecision(0) 
                      << ratio << "x\n";
        }
    }
    std::cout << "\n";
    
    // Step 3: Save initial state
    std::cout << "STEP 3: Save Initial State\n";
    module.saveState("standard_modulation");
    std::cout << "State saved as 'standard_modulation'\n\n";
    
    // Step 4: Test high modulation
    std::cout << "STEP 4: Test High Modulation (δ_sw=0.02)\n";
    module.optimizeForMetric("high_modulation");
    double delta_high = module.getVariable("delta_sw");
    double mod_high = module.getVariable("modulation_factor");
    double U_g2_high = module.computeU_g2(1.496e13);
    std::cout << "High modulation:\n";
    std::cout << "  δ_sw = " << std::fixed << std::setprecision(4) << delta_high << "\n";
    std::cout << "  Modulation factor = " << std::fixed << std::setprecision(0) << mod_high << "x\n";
    std::cout << "  U_g2 at R_b = " << std::scientific << U_g2_high << " J/m³\n\n";
    
    // Step 5: Test low modulation
    std::cout << "STEP 5: Test Low Modulation (δ_sw=0.005)\n";
    module.restoreState("standard_modulation");
    module.optimizeForMetric("low_modulation");
    double delta_low = module.getVariable("delta_sw");
    double mod_low = module.getVariable("modulation_factor");
    double U_g2_low = module.computeU_g2(1.496e13);
    std::cout << "Low modulation:\n";
    std::cout << "  δ_sw = " << std::fixed << std::setprecision(4) << delta_low << "\n";
    std::cout << "  Modulation factor = " << std::fixed << std::setprecision(0) << mod_low << "x\n";
    std::cout << "  U_g2 at R_b = " << std::scientific << U_g2_low << " J/m³\n\n";
    
    // Step 6: Test fast solar wind
    std::cout << "STEP 6: Test Fast Solar Wind (750 km/s)\n";
    module.restoreState("standard_modulation");
    module.optimizeForMetric("fast_wind");
    double v_fast = module.getVariable("v_sw");
    double mod_fast = module.getVariable("modulation_factor");
    std::cout << "Fast wind:\n";
    std::cout << "  v_sw = " << std::scientific << v_fast << " m/s (";
    std::cout << std::fixed << std::setprecision(0) << (v_fast / 1000.0) << " km/s)\n";
    std::cout << "  Modulation factor = " << std::fixed << std::setprecision(0) << mod_fast << "x\n\n";
    
    // Step 7: Test slow solar wind
    std::cout << "STEP 7: Test Slow Solar Wind (350 km/s)\n";
    module.restoreState("standard_modulation");
    module.optimizeForMetric("slow_wind");
    double v_slow = module.getVariable("v_sw");
    double mod_slow = module.getVariable("modulation_factor");
    std::cout << "Slow wind:\n";
    std::cout << "  v_sw = " << std::scientific << v_slow << " m/s (";
    std::cout << std::fixed << std::setprecision(0) << (v_slow / 1000.0) << " km/s)\n";
    std::cout << "  Modulation factor = " << std::fixed << std::setprecision(0) << mod_slow << "x\n\n";
    
    // Step 8: Test extreme amplification
    std::cout << "STEP 8: Test Extreme Amplification (δ_sw=0.05, v_sw=1000 km/s)\n";
    module.restoreState("standard_modulation");
    module.optimizeForMetric("extreme_amplification");
    double delta_extreme = module.getVariable("delta_sw");
    double v_extreme = module.getVariable("v_sw");
    double mod_extreme = module.getVariable("modulation_factor");
    std::cout << "Extreme amplification:\n";
    std::cout << "  δ_sw = " << std::fixed << std::setprecision(4) << delta_extreme << "\n";
    std::cout << "  v_sw = " << std::fixed << std::setprecision(0) << (v_extreme / 1000.0) << " km/s\n";
    std::cout << "  Modulation factor = " << std::fixed << std::setprecision(0) << mod_extreme << "x\n\n";
    
    // Step 9: Expand modulation scale
    std::cout << "STEP 9: Expand Modulation Scale (δ_sw x2)\n";
    module.restoreState("standard_modulation");
    module.expandModulationScale(2.0, 1.0);
    double delta_expanded = module.getVariable("delta_sw");
    double mod_expanded = module.getVariable("modulation_factor");
    std::cout << "Expanded modulation:\n";
    std::cout << "  δ_sw = " << std::fixed << std::setprecision(4) << delta_expanded << " (2x)\n";
    std::cout << "  Modulation factor = " << std::fixed << std::setprecision(0) << mod_expanded << "x\n\n";
    
    // Step 10: Expand velocity scale
    std::cout << "STEP 10: Expand Velocity Scale (v_sw x1.5)\n";
    module.restoreState("standard_modulation");
    module.expandVelocityScale(1.5, 1.0);
    double v_expanded = module.getVariable("v_sw");
    double mod_v_expanded = module.getVariable("modulation_factor");
    std::cout << "Expanded velocity:\n";
    std::cout << "  v_sw = " << std::scientific << v_expanded << " m/s (";
    std::cout << std::fixed << std::setprecision(0) << (v_expanded / 1000.0) << " km/s, 1.5x)\n";
    std::cout << "  Modulation factor = " << std::fixed << std::setprecision(0) << mod_v_expanded << "x\n\n";
    
    // Step 11: Expand gravity scale
    std::cout << "STEP 11: Expand Gravity Scale (k_2 x2, ρ x2)\n";
    module.restoreState("standard_modulation");
    module.expandGravityScale(2.0, 2.0);
    double k2_expanded = module.getVariable("k_2");
    double rho_expanded = module.getVariable("rho_sum");
    double U_g2_grav_expanded = module.computeU_g2(1.496e13);
    std::cout << "Expanded gravity:\n";
    std::cout << "  k_2 = " << std::fixed << std::setprecision(1) << k2_expanded << " (2x)\n";
    std::cout << "  ρ_sum = " << std::scientific << rho_expanded << " J/m³ (2x)\n";
    std::cout << "  U_g2 at R_b = " << U_g2_grav_expanded << " J/m³\n\n";
    
    // Step 12: Sensitivity analysis
    std::cout << "STEP 12: Sensitivity Analysis (at R_b)\n";
    module.restoreState("standard_modulation");
    std::vector<std::string> params = {"delta_sw", "v_sw", "k_2", "rho_vac_UA", "M_s"};
    auto sensitivities = module.sensitivityAnalysis(params);
    for (const auto& pair : sensitivities) {
        std::cout << "  ∂U_g2/∂" << pair.first << " ≈ " << std::scientific << pair.second << " (normalized)\n";
    }
    std::cout << "\n";
    
    // Step 13: Generate variations
    std::cout << "STEP 13: Generate Parameter Variations (5 variants, ±10%)\n";
    auto variations = module.generateVariations(5, 0.1);
    for (int i = 0; i < variations.size(); ++i) {
        double var_delta = variations[i]["delta_sw"];
        double var_v = variations[i]["v_sw"];
        double var_mod = variations[i]["modulation_factor"];
        std::cout << "  Variant " << (i+1) << ": δ_sw=" << std::fixed << std::setprecision(4) 
                  << var_delta << ", v_sw=" << std::fixed << std::setprecision(0) 
                  << (var_v/1000.0) << " km/s, mod=" << std::fixed << std::setprecision(0) 
                  << var_mod << "x\n";
    }
    std::cout << "\n";
    
    // Step 14: Auto-refine to target modulation factor
    std::cout << "STEP 14: Auto-Refine to Target Modulation Factor (10000x)\n";
    module.restoreState("standard_modulation");
    module.autoRefineParameters("modulation_factor", 10000.0);
    double refined_delta = module.getVariable("delta_sw");
    double refined_mod = module.getVariable("modulation_factor");
    std::cout << "Refined parameters:\n";
    std::cout << "  δ_sw = " << std::fixed << std::setprecision(6) << refined_delta << "\n";
    std::cout << "  Modulation factor = " << std::fixed << std::setprecision(0) 
              << refined_mod << "x (target: 10000)\n\n";
    
    // Step 15: Auto-refine to target velocity
    std::cout << "STEP 15: Auto-Refine to Target Velocity (600 km/s)\n";
    module.restoreState("standard_modulation");
    module.autoRefineParameters("v_sw_km_s", 600.0);
    double refined_v = module.getVariable("v_sw");
    double refined_v_mod = module.getVariable("modulation_factor");
    std::cout << "Refined v_sw = " << std::fixed << std::setprecision(0) 
              << (refined_v / 1000.0) << " km/s (target: 600)\n";
    std::cout << "Modulation factor = " << std::fixed << std::setprecision(0) 
              << refined_v_mod << "x\n\n";
    
    // Step 16: Auto-refine to target U_g2 at R_b
    std::cout << "STEP 16: Auto-Refine to Target U_g2 at R_b (2e53 J/m³)\n";
    module.restoreState("standard_modulation");
    module.autoRefineParameters("U_g2_at_Rb", 2.0e53);
    double refined_k2 = module.getVariable("k_2");
    double refined_U_g2 = module.computeU_g2(module.getVariable("R_b"));
    std::cout << "Refined k_2 = " << std::fixed << std::setprecision(2) << refined_k2 << "\n";
    std::cout << "Achieved U_g2 at R_b = " << std::scientific << refined_U_g2 
              << " J/m³ (target: 2e53)\n\n";
    
    // Step 17: Calibrate to observations
    std::cout << "STEP 17: Calibrate to Observational Data\n";
    module.restoreState("standard_modulation");
    std::map<std::string, double> observations;
    observations["delta_sw"] = 0.015;        // Observed higher modulation
    observations["v_sw"] = 6.0e5;            // Observed 600 km/s
    module.calibrateToObservations(observations);
    std::cout << "Calibrated parameters:\n";
    std::cout << "  δ_sw = " << std::fixed << std::setprecision(4) 
              << module.getVariable("delta_sw") << "\n";
    std::cout << "  v_sw = " << std::fixed << std::setprecision(0) 
              << (module.getVariable("v_sw") / 1000.0) << " km/s\n";
    std::cout << "  Modulation factor = " << std::fixed << std::setprecision(0) 
              << module.getVariable("modulation_factor") << "x\n\n";
    
    // Step 18: Mutate parameters
    std::cout << "STEP 18: Mutate Parameters (5% random variation)\n";
    module.mutateParameters(0.05);
    double mutated_delta = module.getVariable("delta_sw");
    double mutated_v = module.getVariable("v_sw");
    double mutated_mod = module.getVariable("modulation_factor");
    std::cout << "Mutated parameters:\n";
    std::cout << "  δ_sw = " << std::fixed << std::setprecision(4) << mutated_delta << "\n";
    std::cout << "  v_sw = " << std::fixed << std::setprecision(0) << (mutated_v / 1000.0) << " km/s\n";
    std::cout << "  Modulation factor = " << std::fixed << std::setprecision(0) << mutated_mod << "x\n\n";
    
    // Step 19: Validate consistency
    std::cout << "STEP 19: Validate Consistency\n";
    bool valid = module.validateConsistency();
    std::cout << "Consistency check: " << (valid ? "PASSED" : "FAILED") << "\n\n";
    
    // Step 20: Introduce anomaly and auto-correct
    std::cout << "STEP 20: Introduce Anomaly and Auto-Correct\n";
    module.createVariable("delta_sw_anomaly", -0.1);
    module.removeVariable("delta_sw");
    module.createVariable("delta_sw", -0.1);  // Invalid negative δ_sw
    std::cout << "Introduced invalid δ_sw = " << module.getVariable("delta_sw") << "\n";
    module.autoCorrectAnomalies();
    std::cout << "Auto-corrected δ_sw = " << std::fixed << std::setprecision(4) 
              << module.getVariable("delta_sw") << "\n";
    std::cout << "Modulation factor = " << std::fixed << std::setprecision(0) 
              << module.getVariable("modulation_factor") << "x\n";
    bool valid_after = module.validateConsistency();
    std::cout << "Consistency after correction: " << (valid_after ? "PASSED" : "FAILED") << "\n\n";
    
    // Step 21: List saved states and export
    std::cout << "STEP 21: List Saved States and Export Final State\n";
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
// #include "SolarWindModulationModule.h"
// int main() {
//     SolarWindModulationModule mod;
//     double mod_f = mod.computeModulationFactor();
//     std::cout << "Modulation Factor = " << mod_f << std::endl;
//     double u_g2 = mod.computeU_g2(1.496e13);
//     std::cout << "U_g2 = " << u_g2 << " J/m�\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("delta_sw", 0.02);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o sw_mod_test sw_mod_test.cpp SolarWindModulationModule.cpp -lm
// Sample: Factor=5001; U_g2?1.18e53 J/m�; amplifies outer bubble gravity.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

SolarWindModulationModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeDelta_sw, computeModulationFactor, computeU_g2, computeU_g2_no_mod) are clear, concise, and variable - driven.
- Automatic recalculation of derived variables(modulation_factor, rho_sum) when dependencies change.
- Output and debugging functions(printVariables, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.
- Models strong amplification of gravity terms via solar wind modulation.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in solar wind modulation modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.