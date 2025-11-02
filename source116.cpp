// SolarWindVelocityModule.h
// Modular C++ implementation of the Solar Wind Velocity (v_sw) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes v_sw=5e5 m/s (500 km/s); scales (1 + δ_sw v_sw) in Universal Gravity U_g2 term.
// Pluggable: #include "SolarWindVelocityModule.h"
// SolarWindVelocityModule mod; mod.computeU_g2(1.496e13); mod.updateVariable("v_sw", new_value);
// Variables in std::map; example for Sun at r=R_b=1.496e13 m; amplification ~5001x.
// Approximations: S(r - R_b)=1; H_SCm=1; E_react=1e46; ρ_sum=7.80e-36 J/m³.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef SOLAR_WIND_VELOCITY_MODULE_H
#define SOLAR_WIND_VELOCITY_MODULE_H

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

class SolarWindVelocityModule {
private:
    std::map<std::string, double> variables;
    double computeModulationFactor();
    double computeU_g2(double r);

public:
    // Constructor: Initialize with framework defaults (Sun)
    SolarWindVelocityModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeV_sw();  // 5e5 m/s
    double computeV_swKmS();  // 500 km/s
    double computeModulationFactor();  // 1 + δ_sw v_sw
    double computeU_g2(double r);  // U_g2 with modulation (J/m^3)
    double computeU_g2_no_sw(double r);  // Without v_sw (set=0)

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
    void expandParameterSpace(double velocity_scale, double modulation_scale, double gravity_scale);
    void expandVelocityScale(double vsw_factor, double mach_factor);       // v_sw and Mach characteristics
    void expandModulationScale(double delta_factor, double coupling_factor); // δ_sw and modulation coupling
    void expandWindScale(double range_factor, double pressure_factor);      // Spatial and pressure effects

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

#endif // SOLAR_WIND_VELOCITY_MODULE_H

// SolarWindVelocityModule.cpp
#include "SolarWindVelocityModule.h"

// Constructor: Set framework defaults (Sun at r=R_b)
SolarWindVelocityModule::SolarWindVelocityModule() {
    // Universal constants
    variables["v_sw"] = 5e5;                        // m/s
    variables["delta_sw"] = 0.01;                   // Unitless
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
void SolarWindVelocityModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
        if (name == "v_sw" || name == "delta_sw") {
            variables["modulation_factor"] = computeModulationFactor();
        } else if (name == "rho_vac_UA" || name == "rho_vac_SCm") {
            variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void SolarWindVelocityModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        if (name == "v_sw" || name == "delta_sw") {
            variables["modulation_factor"] = computeModulationFactor();
        } else if (name == "rho_vac_UA" || name == "rho_vac_SCm") {
            variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void SolarWindVelocityModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute v_sw (m/s)
double SolarWindVelocityModule::computeV_sw() {
    return variables["v_sw"];
}

// v_sw in km/s
double SolarWindVelocityModule::computeV_swKmS() {
    return computeV_sw() / 1e3;
}

// Compute 1 + δ_sw * v_sw
double SolarWindVelocityModule::computeModulationFactor() {
    return 1.0 + variables["delta_sw"] * computeV_sw();
}

// Compute U_g2 with modulation
double SolarWindVelocityModule::computeU_g2(double r) {
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

// U_g2 without solar wind (v_sw=0)
double SolarWindVelocityModule::computeU_g2_no_sw(double r) {
    double orig_v = variables["v_sw"];
    variables["v_sw"] = 0.0;
    double result = computeU_g2(r);
    variables["v_sw"] = orig_v;
    return result;
}

// Equation text
std::string SolarWindVelocityModule::getEquationText() {
    return "U_g2 = k_2 * [(ρ_vac,[UA] + ρ_vac,[SCm]) M_s / r^2] * S(r - R_b) * (1 + δ_sw v_sw) * H_SCm * E_react\n"
           "Where v_sw = 5e5 m/s (500 km/s, typical solar wind speed at 1 AU+);\n"
           "Modulation = 1 + 0.01 * v_sw ≈5001 (amplifies ~5000x).\n"
           "Example r=R_b=1.496e13 m: U_g2 ≈1.18e53 J/m³ (with); ≈2.36e49 J/m³ (without v_sw; ~5000x less).\n"
           "Role: Solar wind momentum/pressure enhances external gravity beyond R_b (heliosphere).\n"
           "UQFF: Models wind shaping of fields; key for heliodynamics/nebular formation.";
}

// Print variables
void SolarWindVelocityModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// ===== ENHANCED METHODS IMPLEMENTATION =====

namespace solar_wind_velocity_saved_states {
    static std::map<std::string, std::map<std::string, double>> saved_states;
}

// Variable Management
void SolarWindVelocityModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void SolarWindVelocityModule::removeVariable(const std::string& name) {
    variables.erase(name);
}

void SolarWindVelocityModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
    }
}

std::vector<std::string> SolarWindVelocityModule::listVariables() const {
    std::vector<std::string> names;
    for (const auto& pair : variables) {
        names.push_back(pair.first);
    }
    return names;
}

std::string SolarWindVelocityModule::getSystemName() const {
    return "Solar_Wind_Velocity_UQFF";
}

// Batch Operations
void SolarWindVelocityModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = func(variables[name]);
        }
    }
    // Update derived
    variables["modulation_factor"] = computeModulationFactor();
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
}

void SolarWindVelocityModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// Self-Expansion: Domain-Specific Scales
void SolarWindVelocityModule::expandParameterSpace(double velocity_scale, double modulation_scale, double gravity_scale) {
    variables["v_sw"] *= velocity_scale;
    variables["delta_sw"] *= modulation_scale;
    variables["k_2"] *= gravity_scale;
    variables["rho_vac_UA"] *= gravity_scale;
    variables["rho_vac_SCm"] *= gravity_scale;
    
    // Update derived
    variables["modulation_factor"] = computeModulationFactor();
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
}

void SolarWindVelocityModule::expandVelocityScale(double vsw_factor, double mach_factor) {
    variables["v_sw"] *= vsw_factor;
    
    // Mach number considerations (for supersonic/subsonic regimes)
    // Store Mach-related parameters for advanced modeling
    if (variables.find("mach_number") == variables.end()) {
        // Typical solar wind is supersonic: M ~ 8-10
        variables["mach_number"] = 9.0;
    }
    variables["mach_number"] *= mach_factor;
    
    variables["modulation_factor"] = computeModulationFactor();
}

void SolarWindVelocityModule::expandModulationScale(double delta_factor, double coupling_factor) {
    variables["delta_sw"] *= delta_factor;
    
    // Coupling factor affects how velocity translates to modulation
    // Store for advanced coupling models
    if (variables.find("coupling_efficiency") == variables.end()) {
        variables["coupling_efficiency"] = 1.0;
    }
    variables["coupling_efficiency"] *= coupling_factor;
    
    variables["modulation_factor"] = computeModulationFactor();
}

void SolarWindVelocityModule::expandWindScale(double range_factor, double pressure_factor) {
    variables["R_b"] *= range_factor;  // Scale heliopause boundary
    variables["r"] *= range_factor;    // Scale current radius
    
    // Pressure factor affects wind pressure (stored for advanced models)
    if (variables.find("wind_pressure") == variables.end()) {
        // Dynamic pressure: ρ v²
        // For v_sw = 5e5 m/s, n ~ 5e6 m⁻³, ρ ~ 8.4e-21 kg/m³
        // P_dyn ~ 2.1e-9 Pa
        variables["wind_pressure"] = 2.1e-9;
    }
    variables["wind_pressure"] *= pressure_factor;
}

// Self-Refinement
void SolarWindVelocityModule::autoRefineParameters(const std::string& target, double goal) {
    if (target == "v_sw") {
        variables["v_sw"] = goal;
        variables["modulation_factor"] = computeModulationFactor();
    } else if (target == "v_sw_km_s") {
        // Convert km/s to m/s
        variables["v_sw"] = goal * 1000.0;
        variables["modulation_factor"] = computeModulationFactor();
    } else if (target == "modulation_factor") {
        // Solve: goal = 1 + δ_sw * v_sw
        // v_sw = (goal - 1) / δ_sw
        if (variables["delta_sw"] > 0) {
            variables["v_sw"] = (goal - 1.0) / variables["delta_sw"];
            variables["modulation_factor"] = computeModulationFactor();
        }
    } else if (target == "amplification") {
        // Same as modulation_factor
        if (variables["delta_sw"] > 0) {
            variables["v_sw"] = (goal - 1.0) / variables["delta_sw"];
            variables["modulation_factor"] = computeModulationFactor();
        }
    } else if (target == "mach_number") {
        if (variables.find("mach_number") == variables.end()) {
            variables["mach_number"] = 9.0;
        }
        variables["mach_number"] = goal;
    } else if (target == "wind_pressure") {
        if (variables.find("wind_pressure") == variables.end()) {
            variables["wind_pressure"] = 2.1e-9;
        }
        variables["wind_pressure"] = goal;
    } else if (target == "dynamic_pressure_Pa") {
        // P_dyn = ρ v² / 2 (for solar wind)
        // Approximate ρ from n ~ 5e6 m⁻³, m_p = 1.67e-27 kg
        // ρ ~ n * m_p ~ 8.4e-21 kg/m³
        // v = sqrt(2 * P / ρ)
        double rho_sw = 8.4e-21;  // kg/m³
        variables["v_sw"] = std::sqrt(2.0 * goal / rho_sw);
        variables["modulation_factor"] = computeModulationFactor();
    } else if (target == "U_g2_at_Rb") {
        // Target specific U_g2 at R_b by adjusting k_2
        double current_U_g2 = computeU_g2(variables["R_b"]);
        if (current_U_g2 > 0) {
            variables["k_2"] *= (goal / current_U_g2);
        }
    }
}

void SolarWindVelocityModule::calibrateToObservations(const std::map<std::string, double>& observations) {
    for (const auto& obs : observations) {
        if (variables.find(obs.first) != variables.end()) {
            variables[obs.first] = obs.second;
            if (obs.first == "v_sw" || obs.first == "delta_sw") {
                variables["modulation_factor"] = computeModulationFactor();
            }
            else if (obs.first == "rho_vac_UA" || obs.first == "rho_vac_SCm") {
                variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
            }
        }
    }
}

void SolarWindVelocityModule::optimizeForMetric(const std::string& metric) {
    if (metric == "fast_wind") {
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
    } else if (metric == "very_fast_wind") {
        // Extreme fast: 900-1000 km/s (coronal holes)
        variables["v_sw"] = 9.5e5;  // 950 km/s
        variables["modulation_factor"] = computeModulationFactor();
    } else if (metric == "very_slow_wind") {
        // Very slow: 200-250 km/s
        variables["v_sw"] = 2.25e5;  // 225 km/s
        variables["modulation_factor"] = computeModulationFactor();
    } else if (metric == "high_mach") {
        // High Mach number (supersonic)
        if (variables.find("mach_number") == variables.end()) {
            variables["mach_number"] = 9.0;
        }
        variables["mach_number"] = 15.0;
    } else if (metric == "low_mach") {
        // Lower Mach number
        if (variables.find("mach_number") == variables.end()) {
            variables["mach_number"] = 9.0;
        }
        variables["mach_number"] = 5.0;
    } else if (metric == "high_pressure") {
        // High dynamic pressure
        if (variables.find("wind_pressure") == variables.end()) {
            variables["wind_pressure"] = 2.1e-9;
        }
        variables["wind_pressure"] = 5.0e-9;  // Pa
    } else if (metric == "low_pressure") {
        // Low dynamic pressure
        if (variables.find("wind_pressure") == variables.end()) {
            variables["wind_pressure"] = 2.1e-9;
        }
        variables["wind_pressure"] = 1.0e-9;  // Pa
    }
}

// Parameter Exploration
std::vector<std::map<std::string, double>> SolarWindVelocityModule::generateVariations(int count, double variation_pct) {
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
void SolarWindVelocityModule::mutateParameters(double mutation_rate) {
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

void SolarWindVelocityModule::evolveSystem(int generations, std::function<double()> fitness_func) {
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
void SolarWindVelocityModule::saveState(const std::string& label) {
    solar_wind_velocity_saved_states::saved_states[label] = variables;
}

void SolarWindVelocityModule::restoreState(const std::string& label) {
    if (solar_wind_velocity_saved_states::saved_states.find(label) != solar_wind_velocity_saved_states::saved_states.end()) {
        variables = solar_wind_velocity_saved_states::saved_states[label];
    }
}

std::vector<std::string> SolarWindVelocityModule::listSavedStates() const {
    std::vector<std::string> labels;
    for (const auto& pair : solar_wind_velocity_saved_states::saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string SolarWindVelocityModule::exportState() const {
    std::ostringstream oss;
    oss << "SolarWindVelocity_State_Export:\n";
    for (const auto& pair : variables) {
        oss << pair.first << "=" << std::scientific << pair.second << "\n";
    }
    return oss.str();
}

// System Analysis
std::map<std::string, double> SolarWindVelocityModule::sensitivityAnalysis(const std::vector<std::string>& params) {
    std::map<std::string, double> sensitivities;
    double r_test = variables["R_b"];  // Test at heliopause boundary
    double baseline = computeU_g2(r_test);
    
    for (const auto& param : params) {
        if (variables.find(param) != variables.end() && 
            param != "modulation_factor" && param != "rho_sum" && param != "S_r_Rb") {
            double original = variables[param];
            variables[param] = original * 1.01;
            
            if (param == "v_sw" || param == "delta_sw") {
                variables["modulation_factor"] = computeModulationFactor();
            }
            else if (param == "rho_vac_UA" || param == "rho_vac_SCm") {
                variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
            }
            
            double perturbed = computeU_g2(r_test);
            sensitivities[param] = (perturbed - baseline) / baseline;
            
            // Restore
            variables[param] = original;
            if (param == "v_sw" || param == "delta_sw") {
                variables["modulation_factor"] = computeModulationFactor();
            }
            else if (param == "rho_vac_UA" || param == "rho_vac_SCm") {
                variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
            }
        }
    }
    return sensitivities;
}

std::string SolarWindVelocityModule::generateReport() const {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(3);
    oss << "===== Solar Wind Velocity Module Report =====\n";
    oss << "System: " << getSystemName() << "\n\n";
    
    oss << "Velocity Parameters:\n";
    oss << "  v_sw = " << variables.at("v_sw") << " m/s (";
    oss << std::fixed << std::setprecision(0) << (variables.at("v_sw") / 1000.0) << " km/s)\n";
    oss << std::scientific;
    oss << "  δ_sw = " << variables.at("delta_sw") << " (unitless)\n";
    oss << "  Modulation factor = " << variables.at("modulation_factor") << "x\n";
    double amplification = variables.at("modulation_factor");
    oss << "  Amplification = " << std::fixed << std::setprecision(0) 
        << amplification << "x (vs no wind)\n\n";
    
    // Mach number and pressure (if available)
    oss << std::scientific;
    if (variables.find("mach_number") != variables.end()) {
        oss << "Wind Characteristics:\n";
        oss << "  Mach number = " << std::fixed << std::setprecision(1) 
            << variables.at("mach_number") << " (supersonic)\n";
    }
    if (variables.find("wind_pressure") != variables.end()) {
        oss << "  Dynamic pressure = " << std::scientific << variables.at("wind_pressure") << " Pa\n";
    }
    oss << "\n";
    
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
            double k_2 = variables.at("k_2");
            double rho_sum = variables.at("rho_sum");
            double M_s = variables.at("M_s");
            double mod_factor = variables.at("modulation_factor");
            double h_scm = variables.at("H_SCm");
            double e_react = variables.at("E_react");
            double U_g2 = k_2 * (rho_sum * M_s / (r_m * r_m)) * s_step * mod_factor * h_scm * e_react;
            
            // Without v_sw
            double U_g2_no_sw = k_2 * (rho_sum * M_s / (r_m * r_m)) * s_step * 1.0 * h_scm * e_react;
            
            oss << "  r=" << std::fixed << std::setprecision(0) << r_AU << " AU: ";
            oss << "U_g2=" << std::scientific << U_g2 << " J/m³ (with wind), ";
            oss << U_g2_no_sw << " J/m³ (no wind)\n";
        } else {
            oss << "  r=" << std::fixed << std::setprecision(0) << r_AU 
                << " AU: Below R_b, S=0, U_g2=0\n";
        }
    }
    oss << "\n";
    
    oss << "Physical Interpretation:\n";
    double v_sw_km = variables.at("v_sw") / 1000.0;
    if (v_sw_km > 800.0) {
        oss << "  Very fast solar wind regime (>800 km/s, coronal holes)\n";
    } else if (v_sw_km > 600.0) {
        oss << "  Fast solar wind (600-800 km/s)\n";
    } else if (v_sw_km > 400.0) {
        oss << "  Standard solar wind (400-600 km/s)\n";
    } else if (v_sw_km > 250.0) {
        oss << "  Slow solar wind (250-400 km/s)\n";
    } else {
        oss << "  Very slow solar wind (<250 km/s)\n";
    }
    
    if (amplification > 10000.0) {
        oss << "  Extreme amplification (>10000x)\n";
    } else if (amplification > 5000.0) {
        oss << "  Very high amplification (5000-10000x)\n";
    } else if (amplification > 1000.0) {
        oss << "  High amplification (1000-5000x)\n";
    } else {
        oss << "  Moderate amplification (<1000x)\n";
    }
    
    oss << "\n  Applications:\n";
    oss << "    - Heliosphere dynamics: Wind velocity shapes outer bubble\n";
    oss << "    - Momentum transfer: Kinetic energy → gravity enhancement\n";
    oss << "    - Pressure balance: Wind ram pressure vs ISM\n";
    oss << "    - Heliopause structure: Termination shock formation\n";
    oss << "    - Stellar astrospheres: General wind-driven boundaries\n";
    oss << "    - Space weather: CME/flare-driven velocity variations\n";
    
    return oss.str();
}

bool SolarWindVelocityModule::validateConsistency() const {
    bool valid = true;
    
    // Check v_sw is positive
    if (variables.find("v_sw") != variables.end() && variables.at("v_sw") <= 0) {
        std::cerr << "Error: v_sw <= 0 (solar wind velocity must be positive)\n";
        valid = false;
    }
    
    // Check δ_sw is non-negative
    if (variables.find("delta_sw") != variables.end() && variables.at("delta_sw") < 0) {
        std::cerr << "Error: delta_sw < 0 (modulation factor must be non-negative)\n";
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
    
    // Check positive densities
    if (variables.at("rho_vac_UA") <= 0 || variables.at("rho_vac_SCm") <= 0) {
        std::cerr << "Error: Vacuum densities must be positive\n";
        valid = false;
    }
    
    // Check Mach number is reasonable (if exists)
    if (variables.find("mach_number") != variables.end()) {
        if (variables.at("mach_number") < 1.0 || variables.at("mach_number") > 20.0) {
            std::cerr << "Warning: Mach number outside typical range [1, 20] (current: " 
                      << variables.at("mach_number") << ")\n";
        }
    }
    
    return valid;
}

void SolarWindVelocityModule::autoCorrectAnomalies() {
    // Reset v_sw to typical value if out of range
    double v_sw_km = variables["v_sw"] / 1000.0;
    if (variables["v_sw"] <= 0 || v_sw_km > 2000.0 || v_sw_km < 100.0) {
        variables["v_sw"] = 5.0e5;  // Standard 500 km/s
        variables["modulation_factor"] = computeModulationFactor();
    }
    
    // Reset δ_sw if out of range
    if (variables["delta_sw"] < 0 || variables["delta_sw"] > 0.1) {
        variables["delta_sw"] = 0.01;  // Standard value
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
    
    // Correct Mach number if exists and out of range
    if (variables.find("mach_number") != variables.end()) {
        if (variables["mach_number"] < 1.0 || variables["mach_number"] > 20.0) {
            variables["mach_number"] = 9.0;  // Typical supersonic value
        }
    }
}

// Example usage in base program (snippet)
int main() {
    SolarWindVelocityModule module;
    std::cout << "===== Solar Wind Velocity Module Enhanced Demonstration =====\n\n";
    
    // Step 1: Report initial state
    std::cout << "STEP 1: Initial Configuration (v_sw=500 km/s)\n";
    std::cout << module.generateReport() << "\n";
    
    // Step 2: Compute U_g2 at various radii
    std::cout << "STEP 2: U_g2 at Various Radii (with and without wind)\n";
    std::vector<double> test_radii = {50.0, 100.0, 150.0, 200.0, 500.0};
    for (double r_AU : test_radii) {
        double r_m = r_AU * 1.496e11;
        double U_g2_with = module.computeU_g2(r_m);
        double U_g2_without = module.computeU_g2_no_sw(r_m);
        if (U_g2_with > 0) {
            double ratio = U_g2_with / U_g2_without;
            std::cout << "  r=" << std::fixed << std::setprecision(0) << r_AU << " AU: ";
            std::cout << "U_g2=" << std::scientific << U_g2_with << " J/m³ (with wind), ";
            std::cout << U_g2_without << " J/m³ (no wind), ratio=" << std::fixed << std::setprecision(0) 
                      << ratio << "x\n";
        }
    }
    std::cout << "\n";
    
    // Step 3: Save initial state
    std::cout << "STEP 3: Save Initial State\n";
    module.saveState("standard_500km_s");
    std::cout << "State saved as 'standard_500km_s'\n\n";
    
    // Step 4: Test fast wind
    std::cout << "STEP 4: Test Fast Solar Wind (750 km/s)\n";
    module.optimizeForMetric("fast_wind");
    double v_fast = module.getVariable("v_sw");
    double mod_fast = module.getVariable("modulation_factor");
    double U_g2_fast = module.computeU_g2(1.496e13);
    std::cout << "Fast wind:\n";
    std::cout << "  v_sw = " << std::fixed << std::setprecision(0) << (v_fast / 1000.0) << " km/s\n";
    std::cout << "  Modulation factor = " << std::fixed << std::setprecision(0) << mod_fast << "x\n";
    std::cout << "  U_g2 at R_b = " << std::scientific << U_g2_fast << " J/m³\n\n";
    
    // Step 5: Test slow wind
    std::cout << "STEP 5: Test Slow Solar Wind (350 km/s)\n";
    module.restoreState("standard_500km_s");
    module.optimizeForMetric("slow_wind");
    double v_slow = module.getVariable("v_sw");
    double mod_slow = module.getVariable("modulation_factor");
    double U_g2_slow = module.computeU_g2(1.496e13);
    std::cout << "Slow wind:\n";
    std::cout << "  v_sw = " << std::fixed << std::setprecision(0) << (v_slow / 1000.0) << " km/s\n";
    std::cout << "  Modulation factor = " << std::fixed << std::setprecision(0) << mod_slow << "x\n";
    std::cout << "  U_g2 at R_b = " << std::scientific << U_g2_slow << " J/m³\n\n";
    
    // Step 6: Test very fast wind
    std::cout << "STEP 6: Test Very Fast Wind (950 km/s, coronal holes)\n";
    module.restoreState("standard_500km_s");
    module.optimizeForMetric("very_fast_wind");
    double v_vfast = module.getVariable("v_sw");
    double mod_vfast = module.getVariable("modulation_factor");
    std::cout << "Very fast wind:\n";
    std::cout << "  v_sw = " << std::fixed << std::setprecision(0) << (v_vfast / 1000.0) << " km/s\n";
    std::cout << "  Modulation factor = " << std::fixed << std::setprecision(0) << mod_vfast << "x\n\n";
    
    // Step 7: Test very slow wind
    std::cout << "STEP 7: Test Very Slow Wind (225 km/s)\n";
    module.restoreState("standard_500km_s");
    module.optimizeForMetric("very_slow_wind");
    double v_vslow = module.getVariable("v_sw");
    double mod_vslow = module.getVariable("modulation_factor");
    std::cout << "Very slow wind:\n";
    std::cout << "  v_sw = " << std::fixed << std::setprecision(0) << (v_vslow / 1000.0) << " km/s\n";
    std::cout << "  Modulation factor = " << std::fixed << std::setprecision(0) << mod_vslow << "x\n\n";
    
    // Step 8: Expand velocity scale
    std::cout << "STEP 8: Expand Velocity Scale (v_sw x2)\n";
    module.restoreState("standard_500km_s");
    module.expandVelocityScale(2.0, 1.0);
    double v_expanded = module.getVariable("v_sw");
    double mod_expanded = module.getVariable("modulation_factor");
    std::cout << "Expanded velocity:\n";
    std::cout << "  v_sw = " << std::fixed << std::setprecision(0) << (v_expanded / 1000.0) << " km/s (2x)\n";
    std::cout << "  Modulation factor = " << std::fixed << std::setprecision(0) << mod_expanded << "x\n\n";
    
    // Step 9: Expand modulation scale
    std::cout << "STEP 9: Expand Modulation Scale (δ_sw x1.5)\n";
    module.restoreState("standard_500km_s");
    module.expandModulationScale(1.5, 1.0);
    double delta_expanded = module.getVariable("delta_sw");
    double mod_delta_expanded = module.getVariable("modulation_factor");
    std::cout << "Expanded modulation:\n";
    std::cout << "  δ_sw = " << std::fixed << std::setprecision(4) << delta_expanded << " (1.5x)\n";
    std::cout << "  Modulation factor = " << std::fixed << std::setprecision(0) << mod_delta_expanded << "x\n\n";
    
    // Step 10: Expand wind scale (range and pressure)
    std::cout << "STEP 10: Expand Wind Scale (R_b x1.2, pressure x2)\n";
    module.restoreState("standard_500km_s");
    module.expandWindScale(1.2, 2.0);
    double rb_expanded = module.getVariable("R_b");
    double press_expanded = module.getVariable("wind_pressure");
    std::cout << "Expanded wind scale:\n";
    std::cout << "  R_b = " << std::scientific << rb_expanded << " m (";
    std::cout << std::fixed << std::setprecision(0) << (rb_expanded / 1.496e11) << " AU)\n";
    std::cout << "  Wind pressure = " << std::scientific << press_expanded << " Pa\n\n";
    
    // Step 11: Test high Mach number
    std::cout << "STEP 11: Optimize for High Mach Number (M=15)\n";
    module.restoreState("standard_500km_s");
    module.optimizeForMetric("high_mach");
    double mach_high = module.getVariable("mach_number");
    std::cout << "High Mach regime:\n";
    std::cout << "  Mach number = " << std::fixed << std::setprecision(1) << mach_high << "\n\n";
    
    // Step 12: Sensitivity analysis
    std::cout << "STEP 12: Sensitivity Analysis (at R_b)\n";
    module.restoreState("standard_500km_s");
    std::vector<std::string> params = {"v_sw", "delta_sw", "k_2", "rho_vac_UA", "M_s"};
    auto sensitivities = module.sensitivityAnalysis(params);
    for (const auto& pair : sensitivities) {
        std::cout << "  ∂U_g2/∂" << pair.first << " ≈ " << std::scientific << pair.second << " (normalized)\n";
    }
    std::cout << "\n";
    
    // Step 13: Generate variations
    std::cout << "STEP 13: Generate Parameter Variations (5 variants, ±10%)\n";
    auto variations = module.generateVariations(5, 0.1);
    for (int i = 0; i < variations.size(); ++i) {
        double var_v = variations[i]["v_sw"];
        double var_mod = variations[i]["modulation_factor"];
        std::cout << "  Variant " << (i+1) << ": v_sw=" << std::fixed << std::setprecision(0) 
                  << (var_v/1000.0) << " km/s, mod=" << std::fixed << std::setprecision(0) 
                  << var_mod << "x\n";
    }
    std::cout << "\n";
    
    // Step 14: Auto-refine to target velocity
    std::cout << "STEP 14: Auto-Refine to Target Velocity (600 km/s)\n";
    module.restoreState("standard_500km_s");
    module.autoRefineParameters("v_sw_km_s", 600.0);
    double refined_v = module.getVariable("v_sw");
    double refined_mod = module.getVariable("modulation_factor");
    std::cout << "Refined v_sw = " << std::fixed << std::setprecision(0) 
              << (refined_v / 1000.0) << " km/s (target: 600)\n";
    std::cout << "Modulation factor = " << std::fixed << std::setprecision(0) 
              << refined_mod << "x\n\n";
    
    // Step 15: Auto-refine to target modulation factor
    std::cout << "STEP 15: Auto-Refine to Target Modulation Factor (8000x)\n";
    module.restoreState("standard_500km_s");
    module.autoRefineParameters("modulation_factor", 8000.0);
    double refined_v_mod = module.getVariable("v_sw");
    double refined_mod_target = module.getVariable("modulation_factor");
    std::cout << "Refined parameters:\n";
    std::cout << "  v_sw = " << std::fixed << std::setprecision(0) 
              << (refined_v_mod / 1000.0) << " km/s\n";
    std::cout << "  Modulation factor = " << std::fixed << std::setprecision(0) 
              << refined_mod_target << "x (target: 8000)\n\n";
    
    // Step 16: Auto-refine to target Mach number
    std::cout << "STEP 16: Auto-Refine to Target Mach Number (12.0)\n";
    module.restoreState("standard_500km_s");
    module.autoRefineParameters("mach_number", 12.0);
    double refined_mach = module.getVariable("mach_number");
    std::cout << "Refined Mach number = " << std::fixed << std::setprecision(1) 
              << refined_mach << " (target: 12.0)\n\n";
    
    // Step 17: Auto-refine to target dynamic pressure
    std::cout << "STEP 17: Auto-Refine to Target Dynamic Pressure (3e-9 Pa)\n";
    module.restoreState("standard_500km_s");
    module.autoRefineParameters("dynamic_pressure_Pa", 3.0e-9);
    double refined_v_press = module.getVariable("v_sw");
    std::cout << "Refined v_sw = " << std::fixed << std::setprecision(0) 
              << (refined_v_press / 1000.0) << " km/s\n";
    std::cout << "For dynamic pressure = " << std::scientific << 3.0e-9 << " Pa\n\n";
    
    // Step 18: Calibrate to observations
    std::cout << "STEP 18: Calibrate to Observational Data\n";
    module.restoreState("standard_500km_s");
    std::map<std::string, double> observations;
    observations["v_sw"] = 6.5e5;            // Observed 650 km/s
    observations["delta_sw"] = 0.012;        // Observed slightly higher
    module.calibrateToObservations(observations);
    std::cout << "Calibrated parameters:\n";
    std::cout << "  v_sw = " << std::fixed << std::setprecision(0) 
              << (module.getVariable("v_sw") / 1000.0) << " km/s\n";
    std::cout << "  δ_sw = " << std::fixed << std::setprecision(4) 
              << module.getVariable("delta_sw") << "\n";
    std::cout << "  Modulation factor = " << std::fixed << std::setprecision(0) 
              << module.getVariable("modulation_factor") << "x\n\n";
    
    // Step 19: Mutate parameters
    std::cout << "STEP 19: Mutate Parameters (5% random variation)\n";
    module.mutateParameters(0.05);
    double mutated_v = module.getVariable("v_sw");
    double mutated_mod = module.getVariable("modulation_factor");
    std::cout << "Mutated parameters:\n";
    std::cout << "  v_sw = " << std::fixed << std::setprecision(0) << (mutated_v / 1000.0) << " km/s\n";
    std::cout << "  Modulation factor = " << std::fixed << std::setprecision(0) << mutated_mod << "x\n\n";
    
    // Step 20: Validate consistency
    std::cout << "STEP 20: Validate Consistency\n";
    bool valid = module.validateConsistency();
    std::cout << "Consistency check: " << (valid ? "PASSED" : "FAILED") << "\n\n";
    
    // Step 21: Introduce anomaly and auto-correct
    std::cout << "STEP 21: Introduce Anomaly and Auto-Correct\n";
    module.createVariable("v_sw_anomaly", -1000.0);
    module.removeVariable("v_sw");
    module.createVariable("v_sw", -1000.0);  // Invalid negative v_sw
    std::cout << "Introduced invalid v_sw = " << module.getVariable("v_sw") << " m/s\n";
    module.autoCorrectAnomalies();
    std::cout << "Auto-corrected v_sw = " << std::fixed << std::setprecision(0) 
              << (module.getVariable("v_sw") / 1000.0) << " km/s\n";
    std::cout << "Modulation factor = " << std::fixed << std::setprecision(0) 
              << module.getVariable("modulation_factor") << "x\n";
    bool valid_after = module.validateConsistency();
    std::cout << "Consistency after correction: " << (valid_after ? "PASSED" : "FAILED") << "\n\n";
    
    // Step 22: List saved states and export
    std::cout << "STEP 22: List Saved States and Export Final State\n";
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
// #include "SolarWindVelocityModule.h"
// int main() {
//     SolarWindVelocityModule mod;
//     double v = mod.computeV_sw();
//     std::cout << "v_sw = " << v << " m/s (" << mod.computeV_swKmS() << " km/s)\n";
//     double u_g2 = mod.computeU_g2(1.496e13);
//     std::cout << "U_g2 = " << u_g2 << " J/m³\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("v_sw", 4e5);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o sw_vel_test sw_vel_test.cpp SolarWindVelocityModule.cpp -lm
// Sample: v_sw=5e5 m/s (500 km/s); U_g2≈1.18e53 J/m³; amplifies outer bubble.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

SolarWindVelocityModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeV_sw, computeV_swKmS, computeModulationFactor, computeU_g2, computeU_g2_no_sw) are clear, concise, and variable - driven.
- Automatic recalculation of derived variables(modulation_factor, rho_sum) when dependencies change.
- Output and debugging functions(printVariables, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.
- Models strong amplification of gravity terms via solar wind velocity.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in solar wind velocity modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.