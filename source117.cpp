// StellarMassModule.h
// Modular C++ implementation of the Stellar/Planetary Mass (M_s) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes M_s=1.989e30 kg (1 M_sun for Sun); scales M_s / r^2 in Universal Gravity U_g1 and U_g2 terms.
// Pluggable: #include "StellarMassModule.h"
// StellarMassModule mod; mod.computeU_g2(1.496e13); mod.updateVariable("M_s", new_value);
// Variables in std::map; example for Sun at r=1.496e13 m; U_g2 ?1.18e53 J/m�.
// Approximations: S(r - R_b)=1; (1 + ?_sw v_sw)=5001; H_SCm=1; E_react=1e46.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef STELLAR_MASS_MODULE_H
#define STELLAR_MASS_MODULE_H

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

class StellarMassModule {
private:
    std::map<std::string, double> variables;
    double computeM_sOverR2(double r);
    double computeU_g1(double r);
    double computeU_g2(double r);

public:
    // Constructor: Initialize with framework defaults (Sun)
    StellarMassModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeM_s();  // 1.989e30 kg
    double computeM_sInMsun();  // 1 M_sun
    double computeM_sOverR2(double r);  // M_s / r^2 (kg/m�)
    double computeU_g1(double r);  // U_g1 example (J/m^3)
    double computeU_g2(double r);  // U_g2 example (J/m^3)

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
    void expandParameterSpace(double mass_scale, double gravity_scale, double spatial_scale);
    void expandMassScale(double ms_factor, double density_factor);        // M_s and mass characteristics
    void expandGravityScale(double coupling_factor, double field_factor);  // k_1, k_2 and field strength
    void expandSpatialScale(double radius_factor, double boundary_factor); // Spatial extent and boundaries

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

#endif // STELLAR_MASS_MODULE_H

// StellarMassModule.cpp
#include "StellarMassModule.h"

// Constructor: Set framework defaults (Sun at r=R_b)
StellarMassModule::StellarMassModule() {
    // Universal constants
    variables["M_s"] = 1.989e30;                    // kg (Sun)
    variables["M_sun"] = 1.989e30;                  // kg
    variables["k_1"] = 1.5;                         // Coupling for U_g1
    variables["k_2"] = 1.2;                         // Coupling for U_g2
    variables["rho_vac_UA"] = 7.09e-36;             // J/m^3
    variables["rho_vac_SCm"] = 7.09e-37;            // J/m^3
    variables["r"] = 1.496e13;                      // m (example R_b)
    variables["R_b"] = 1.496e13;                    // m
    variables["S_r_Rb"] = 1.0;                      // Step
    variables["delta_sw"] = 0.01;                   // Unitless
    variables["v_sw"] = 5e5;                        // m/s
    variables["H_SCm"] = 1.0;                       // Unitless
    variables["E_react"] = 1e46;                    // J

    // Derived
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
    variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
}

// Update variable
void StellarMassModule::updateVariable(const std::string& name, double value) {
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
void StellarMassModule::addToVariable(const std::string& name, double delta) {
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
void StellarMassModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute M_s (kg)
double StellarMassModule::computeM_s() {
    return variables["M_s"];
}

// M_s in M_sun
double StellarMassModule::computeM_sInMsun() {
    return computeM_s() / variables["M_sun"];
}

// M_s / r^2 (kg/m�)
double StellarMassModule::computeM_sOverR2(double r) {
    variables["r"] = r;
    if (r == 0.0) return 0.0;
    return computeM_s() / (r * r);
}

// U_g1 example (internal, simplified)
double StellarMassModule::computeU_g1(double r) {
    double k_1 = variables["k_1"];
    double rho_sum = variables["rho_sum"];
    double m_over_r2 = computeM_sOverR2(r);
    double e_react = variables["E_react"];
    return k_1 * rho_sum * m_over_r2 * e_react;  // Simplified
}

// U_g2 example (outer bubble)
double StellarMassModule::computeU_g2(double r) {
    variables["r"] = r;
    double k_2 = variables["k_2"];
    double rho_sum = variables["rho_sum"];
    double s_step = (r >= variables["R_b"]) ? 1.0 : 0.0;
    double swirl_factor = variables["swirl_factor"];
    double h_scm = variables["H_SCm"];
    double e_react = variables["E_react"];
    return k_2 * rho_sum * computeM_sOverR2(r) * s_step * swirl_factor * h_scm * e_react;
}

// Equation text
std::string StellarMassModule::getEquationText() {
    return "U_g1 = k_1 * ?_vac,[UA/SCm] * (M_s / r^2) * ... E_react (internal dipole);\n"
           "U_g2 = k_2 * ?_vac,[UA/SCm] * (M_s / r^2) * S(r - R_b) * (1 + ?_sw v_sw) * H_SCm * E_react (outer bubble).\n"
           "Where M_s = 1.989e30 kg (1 M_sun for Sun).\n"
           "Scales gravity by mass; M_s / r^2 ?8.89e3 kg/m� at r=1.496e13 m.\n"
           "Example U_g2 (r=R_b): ?1.18e53 J/m�.\n"
           "Role: Central mass drives internal/external gravity; stellar/planetary dynamics.\n"
           "UQFF: Mass-dependent fields for nebulae/formation/mergers.";
}

// Print variables
void StellarMassModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// ===== ENHANCED METHODS IMPLEMENTATION =====

namespace stellar_mass_saved_states {
    static std::map<std::string, std::map<std::string, double>> saved_states;
}

// Variable Management
void StellarMassModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void StellarMassModule::removeVariable(const std::string& name) {
    variables.erase(name);
}

void StellarMassModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
    }
}

std::vector<std::string> StellarMassModule::listVariables() const {
    std::vector<std::string> names;
    for (const auto& pair : variables) {
        names.push_back(pair.first);
    }
    return names;
}

std::string StellarMassModule::getSystemName() const {
    return "Stellar_Planetary_Mass_UQFF";
}

// Batch Operations
void StellarMassModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = func(variables[name]);
        }
    }
    // Update derived
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
    variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
}

void StellarMassModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// Self-Expansion: Domain-Specific Scales
void StellarMassModule::expandParameterSpace(double mass_scale, double gravity_scale, double spatial_scale) {
    variables["M_s"] *= mass_scale;
    variables["k_1"] *= gravity_scale;
    variables["k_2"] *= gravity_scale;
    variables["rho_vac_UA"] *= gravity_scale;
    variables["rho_vac_SCm"] *= gravity_scale;
    variables["R_b"] *= spatial_scale;
    variables["r"] *= spatial_scale;
    
    // Update derived
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
}

void StellarMassModule::expandMassScale(double ms_factor, double density_factor) {
    variables["M_s"] *= ms_factor;
    
    // Store mass characteristics for advanced modeling
    if (variables.find("mass_density_kg_m3") == variables.end()) {
        // Solar average density: ρ_sun ~ 1408 kg/m³
        // R_sun = 6.96e8 m, V = 4/3 π R³ ~ 1.41e27 m³
        // ρ = M_s / V ~ 1408 kg/m³
        variables["mass_density_kg_m3"] = 1408.0;
    }
    variables["mass_density_kg_m3"] *= density_factor;
    
    // Surface gravity scaling
    if (variables.find("surface_gravity_m_s2") == variables.end()) {
        // g = G M_s / R_sun² ~ 274 m/s²
        variables["surface_gravity_m_s2"] = 274.0;
    }
    // g scales with M and 1/R²; if density changes, adjust accordingly
    variables["surface_gravity_m_s2"] *= ms_factor;
}

void StellarMassModule::expandGravityScale(double coupling_factor, double field_factor) {
    variables["k_1"] *= coupling_factor;
    variables["k_2"] *= coupling_factor;
    
    // Field strength parameters
    variables["rho_vac_UA"] *= field_factor;
    variables["rho_vac_SCm"] *= field_factor;
    
    // Store gravity field characteristics
    if (variables.find("gravity_field_strength") == variables.end()) {
        variables["gravity_field_strength"] = 1.0;
    }
    variables["gravity_field_strength"] *= field_factor;
    
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
}

void StellarMassModule::expandSpatialScale(double radius_factor, double boundary_factor) {
    variables["r"] *= radius_factor;
    variables["R_b"] *= boundary_factor;  // Scale heliopause/astrosphere boundary
    
    // Store spatial characteristics
    if (variables.find("stellar_radius_m") == variables.end()) {
        // R_sun = 6.96e8 m
        variables["stellar_radius_m"] = 6.96e8;
    }
    variables["stellar_radius_m"] *= radius_factor;
    
    // Hill sphere / Roche lobe considerations
    if (variables.find("hill_sphere_AU") == variables.end()) {
        // For Sun-Earth: r_H ~ 0.01 AU (1.5e9 m) at 1 AU
        variables["hill_sphere_AU"] = 0.01;
    }
    variables["hill_sphere_AU"] *= boundary_factor;
}

// Self-Refinement
void StellarMassModule::autoRefineParameters(const std::string& target, double goal) {
    if (target == "M_s") {
        variables["M_s"] = goal;
    } else if (target == "M_sun_multiple") {
        // Set M_s as multiple of solar mass
        variables["M_s"] = goal * variables["M_sun"];
    } else if (target == "mass_density_kg_m3") {
        if (variables.find("mass_density_kg_m3") == variables.end()) {
            variables["mass_density_kg_m3"] = 1408.0;
        }
        variables["mass_density_kg_m3"] = goal;
    } else if (target == "surface_gravity_m_s2") {
        if (variables.find("surface_gravity_m_s2") == variables.end()) {
            variables["surface_gravity_m_s2"] = 274.0;
        }
        variables["surface_gravity_m_s2"] = goal;
    } else if (target == "M_s_over_r2") {
        // Target specific M_s/r² value by adjusting M_s
        double r_current = variables["r"];
        if (r_current > 0) {
            variables["M_s"] = goal * r_current * r_current;
        }
    } else if (target == "U_g2_at_Rb") {
        // Target specific U_g2 at R_b by adjusting k_2
        double current_U_g2 = computeU_g2(variables["R_b"]);
        if (current_U_g2 > 0) {
            variables["k_2"] *= (goal / current_U_g2);
        }
    } else if (target == "U_g1_at_r") {
        // Target specific U_g1 at current r by adjusting k_1
        double current_U_g1 = computeU_g1(variables["r"]);
        if (current_U_g1 > 0) {
            variables["k_1"] *= (goal / current_U_g1);
        }
    } else if (target == "stellar_radius_m") {
        if (variables.find("stellar_radius_m") == variables.end()) {
            variables["stellar_radius_m"] = 6.96e8;
        }
        variables["stellar_radius_m"] = goal;
    }
}

void StellarMassModule::calibrateToObservations(const std::map<std::string, double>& observations) {
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

void StellarMassModule::optimizeForMetric(const std::string& metric) {
    if (metric == "solar_mass") {
        // Standard Sun: M_s = 1.989e30 kg
        variables["M_s"] = 1.989e30;
    } else if (metric == "jupiter_mass") {
        // Jupiter: M_J ~ 1.898e27 kg (0.000955 M_sun)
        variables["M_s"] = 1.898e27;
    } else if (metric == "earth_mass") {
        // Earth: M_E ~ 5.972e24 kg (3e-6 M_sun)
        variables["M_s"] = 5.972e24;
    } else if (metric == "low_mass_star") {
        // 0.1 M_sun (red dwarf threshold)
        variables["M_s"] = 0.1 * variables["M_sun"];
    } else if (metric == "high_mass_star") {
        // 10 M_sun (massive star)
        variables["M_s"] = 10.0 * variables["M_sun"];
    } else if (metric == "supermassive_star") {
        // 100 M_sun (very massive)
        variables["M_s"] = 100.0 * variables["M_sun"];
    } else if (metric == "neutron_star") {
        // 1.4 M_sun (typical neutron star)
        variables["M_s"] = 1.4 * variables["M_sun"];
    } else if (metric == "white_dwarf") {
        // 0.6 M_sun (typical white dwarf)
        variables["M_s"] = 0.6 * variables["M_sun"];
    } else if (metric == "brown_dwarf") {
        // 0.05 M_sun (50 M_J, brown dwarf)
        variables["M_s"] = 0.05 * variables["M_sun"];
    } else if (metric == "sgr_a_star") {
        // Sagittarius A*: M_BH ~ 4.15e6 M_sun
        variables["M_s"] = 4.15e6 * variables["M_sun"];
    } else if (metric == "high_density") {
        // High density stellar object
        if (variables.find("mass_density_kg_m3") == variables.end()) {
            variables["mass_density_kg_m3"] = 1408.0;
        }
        variables["mass_density_kg_m3"] = 1e5;  // 100,000 kg/m³ (white dwarf regime)
    } else if (metric == "low_density") {
        // Low density object
        if (variables.find("mass_density_kg_m3") == variables.end()) {
            variables["mass_density_kg_m3"] = 1408.0;
        }
        variables["mass_density_kg_m3"] = 100.0;  // Gas giant regime
    }
}

// Parameter Exploration
std::vector<std::map<std::string, double>> StellarMassModule::generateVariations(int count, double variation_pct) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - variation_pct, 1.0 + variation_pct);
    
    for (int i = 0; i < count; ++i) {
        std::map<std::string, double> variant = variables;
        for (auto& pair : variant) {
            if (pair.first != "rho_sum" && pair.first != "swirl_factor" && 
                pair.first != "S_r_Rb" && pair.first != "M_sun") {
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
void StellarMassModule::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - mutation_rate, 1.0 + mutation_rate);
    
    for (auto& pair : variables) {
        if (pair.first != "rho_sum" && pair.first != "swirl_factor" && 
            pair.first != "S_r_Rb" && pair.first != "M_sun") {
            pair.second *= dis(gen);
        }
    }
    // Recalculate derived
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
    variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
}

void StellarMassModule::evolveSystem(int generations, std::function<double()> fitness_func) {
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
void StellarMassModule::saveState(const std::string& label) {
    stellar_mass_saved_states::saved_states[label] = variables;
}

void StellarMassModule::restoreState(const std::string& label) {
    if (stellar_mass_saved_states::saved_states.find(label) != stellar_mass_saved_states::saved_states.end()) {
        variables = stellar_mass_saved_states::saved_states[label];
    }
}

std::vector<std::string> StellarMassModule::listSavedStates() const {
    std::vector<std::string> labels;
    for (const auto& pair : stellar_mass_saved_states::saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string StellarMassModule::exportState() const {
    std::ostringstream oss;
    oss << "StellarMass_State_Export:\n";
    for (const auto& pair : variables) {
        oss << pair.first << "=" << std::scientific << pair.second << "\n";
    }
    return oss.str();
}

// System Analysis
std::map<std::string, double> StellarMassModule::sensitivityAnalysis(const std::vector<std::string>& params) {
    std::map<std::string, double> sensitivities;
    double r_test = variables["R_b"];  // Test at boundary
    double baseline_g1 = computeU_g1(r_test);
    double baseline_g2 = computeU_g2(r_test);
    double baseline = baseline_g2;  // Use U_g2 as primary metric
    
    for (const auto& param : params) {
        if (variables.find(param) != variables.end() && 
            param != "rho_sum" && param != "swirl_factor" && 
            param != "S_r_Rb" && param != "M_sun") {
            double original = variables[param];
            variables[param] = original * 1.01;
            
            if (param == "rho_vac_UA" || param == "rho_vac_SCm") {
                variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
            }
            else if (param == "delta_sw" || param == "v_sw") {
                variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
            }
            
            double perturbed = computeU_g2(r_test);
            sensitivities[param] = (perturbed - baseline) / baseline;
            
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

std::string StellarMassModule::generateReport() const {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(3);
    oss << "===== Stellar/Planetary Mass Module Report =====\n";
    oss << "System: " << getSystemName() << "\n\n";
    
    oss << "Mass Parameters:\n";
    oss << "  M_s = " << variables.at("M_s") << " kg (";
    oss << std::fixed << std::setprecision(3) << (variables.at("M_s") / variables.at("M_sun")) << " M_sun)\n";
    oss << std::scientific;
    oss << "  M_sun = " << variables.at("M_sun") << " kg (reference)\n";
    
    if (variables.find("mass_density_kg_m3") != variables.end()) {
        oss << "  Mass density = " << std::fixed << std::setprecision(0) 
            << variables.at("mass_density_kg_m3") << " kg/m³\n";
    }
    if (variables.find("surface_gravity_m_s2") != variables.end()) {
        oss << "  Surface gravity = " << std::fixed << std::setprecision(1) 
            << variables.at("surface_gravity_m_s2") << " m/s²\n";
    }
    oss << "\n";
    
    oss << std::scientific;
    oss << "Gravity Parameters:\n";
    oss << "  k_1 = " << variables.at("k_1") << " (internal coupling)\n";
    oss << "  k_2 = " << variables.at("k_2") << " (outer coupling)\n";
    oss << "  ρ_vac,[UA] = " << variables.at("rho_vac_UA") << " J/m³\n";
    oss << "  ρ_vac,[SCm] = " << variables.at("rho_vac_SCm") << " J/m³\n";
    oss << "  ρ_sum = " << variables.at("rho_sum") << " J/m³\n\n";
    
    oss << "Spatial Parameters:\n";
    oss << "  R_b = " << variables.at("R_b") << " m (";
    oss << std::fixed << std::setprecision(0) << (variables.at("R_b") / 1.496e11) << " AU)\n";
    oss << std::scientific;
    oss << "  Current r = " << variables.at("r") << " m (";
    oss << std::fixed << std::setprecision(0) << (variables.at("r") / 1.496e11) << " AU)\n";
    
    if (variables.find("stellar_radius_m") != variables.end()) {
        oss << std::scientific;
        oss << "  Stellar radius = " << variables.at("stellar_radius_m") << " m (";
        oss << std::fixed << std::setprecision(2) << (variables.at("stellar_radius_m") / 6.96e8) << " R_sun)\n";
    }
    if (variables.find("hill_sphere_AU") != variables.end()) {
        oss << "  Hill sphere = " << std::fixed << std::setprecision(4) 
            << variables.at("hill_sphere_AU") << " AU\n";
    }
    oss << "\n";
    
    oss << std::scientific;
    oss << "Solar Wind Parameters:\n";
    oss << "  δ_sw = " << variables.at("delta_sw") << "\n";
    oss << "  v_sw = " << variables.at("v_sw") << " m/s (";
    oss << std::fixed << std::setprecision(0) << (variables.at("v_sw") / 1000.0) << " km/s)\n";
    oss << std::scientific;
    oss << "  Swirl factor = " << variables.at("swirl_factor") << "x\n\n";
    
    oss << "Energy Factors:\n";
    oss << "  H_SCm = " << variables.at("H_SCm") << "\n";
    oss << "  E_react = " << variables.at("E_react") << " J\n\n";
    
    // M_s/r² at current position
    double m_over_r2 = 0.0;
    if (variables.at("r") > 0) {
        m_over_r2 = variables.at("M_s") / (variables.at("r") * variables.at("r"));
    }
    oss << "M_s / r² at current r:\n";
    oss << "  M_s/r² = " << m_over_r2 << " kg/m²\n\n";
    
    oss << "U_g1 and U_g2 Computations at Key Radii:\n";
    std::vector<double> test_radii_AU = {0.1, 0.5, 1.0, 5.0, 10.0, 50.0, 100.0, 500.0};
    for (double r_AU : test_radii_AU) {
        double r_m = r_AU * 1.496e11;
        double s_step = (r_m >= variables.at("R_b")) ? 1.0 : 0.0;
        
        // U_g1 (always active)
        double k_1 = variables.at("k_1");
        double rho_sum = variables.at("rho_sum");
        double e_react = variables.at("E_react");
        double m_over_r2_test = variables.at("M_s") / (r_m * r_m);
        double U_g1 = k_1 * rho_sum * m_over_r2_test * e_react;
        
        // U_g2 (only beyond R_b)
        double k_2 = variables.at("k_2");
        double swirl_factor = variables.at("swirl_factor");
        double h_scm = variables.at("H_SCm");
        double U_g2 = k_2 * rho_sum * m_over_r2_test * s_step * swirl_factor * h_scm * e_react;
        
        oss << "  r=" << std::fixed << std::setprecision(1) << r_AU << " AU: ";
        oss << "U_g1=" << std::scientific << U_g1 << " J/m³";
        if (s_step > 0) {
            oss << ", U_g2=" << U_g2 << " J/m³";
        } else {
            oss << ", U_g2=0 (r < R_b)";
        }
        oss << "\n";
    }
    oss << "\n";
    
    oss << "Physical Interpretation:\n";
    double mass_ratio = variables.at("M_s") / variables.at("M_sun");
    if (mass_ratio > 100.0) {
        oss << "  Supermassive object (>100 M_sun, black hole regime)\n";
    } else if (mass_ratio > 10.0) {
        oss << "  Very massive star (10-100 M_sun)\n";
    } else if (mass_ratio > 2.0) {
        oss << "  Massive star (2-10 M_sun)\n";
    } else if (mass_ratio > 0.5) {
        oss << "  Solar-mass star (0.5-2 M_sun)\n";
    } else if (mass_ratio > 0.08) {
        oss << "  Low-mass star / red dwarf (0.08-0.5 M_sun)\n";
    } else if (mass_ratio > 0.01) {
        oss << "  Brown dwarf (0.01-0.08 M_sun, 10-80 M_J)\n";
    } else if (mass_ratio > 1e-4) {
        oss << "  Gas giant planet (0.0001-0.01 M_sun, ~0.1-10 M_J)\n";
    } else {
        oss << "  Rocky planet / asteroid (<0.0001 M_sun)\n";
    }
    
    oss << "\n  Applications:\n";
    oss << "    - Stellar dynamics: Mass determines gravity field strength\n";
    oss << "    - Planetary systems: M_s/r² governs orbital mechanics\n";
    oss << "    - Internal structure: U_g1 drives dipole/core effects\n";
    oss << "    - External fields: U_g2 shapes heliosphere/astrosphere\n";
    oss << "    - Binary systems: Mass ratio determines orbital evolution\n";
    oss << "    - Mergers/formation: Mass accumulation and field coupling\n";
    oss << "    - Nebular dynamics: Mass concentration and collapse\n";
    
    return oss.str();
}

bool StellarMassModule::validateConsistency() const {
    bool valid = true;
    
    // Check M_s is positive
    if (variables.find("M_s") != variables.end() && variables.at("M_s") <= 0) {
        std::cerr << "Error: M_s <= 0 (stellar/planetary mass must be positive)\n";
        valid = false;
    }
    
    // Check M_sun is positive
    if (variables.find("M_sun") != variables.end() && variables.at("M_sun") <= 0) {
        std::cerr << "Error: M_sun <= 0 (solar mass reference must be positive)\n";
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
    
    // Check M_s is in reasonable range (0.01 M_sun to 1e8 M_sun, covers planets to supermassive BH)
    double mass_ratio = variables.at("M_s") / variables.at("M_sun");
    if (mass_ratio < 1e-6 || mass_ratio > 1e9) {
        std::cerr << "Warning: M_s outside typical range [1e-6, 1e9] M_sun (current: " 
                  << mass_ratio << " M_sun)\n";
    }
    
    // Check positive densities
    if (variables.at("rho_vac_UA") <= 0 || variables.at("rho_vac_SCm") <= 0) {
        std::cerr << "Error: Vacuum densities must be positive\n";
        valid = false;
    }
    
    // Check positive couplings
    if (variables.at("k_1") <= 0 || variables.at("k_2") <= 0) {
        std::cerr << "Error: Coupling constants k_1, k_2 must be positive\n";
        valid = false;
    }
    
    // Check mass density (if exists)
    if (variables.find("mass_density_kg_m3") != variables.end()) {
        if (variables.at("mass_density_kg_m3") <= 0 || variables.at("mass_density_kg_m3") > 1e18) {
            std::cerr << "Warning: Mass density outside typical range [0, 1e18] kg/m³ (current: " 
                      << variables.at("mass_density_kg_m3") << " kg/m³)\n";
        }
    }
    
    return valid;
}

void StellarMassModule::autoCorrectAnomalies() {
    // Reset M_s to solar mass if out of range
    double mass_ratio = variables["M_s"] / variables["M_sun"];
    if (variables["M_s"] <= 0 || mass_ratio > 1e9 || mass_ratio < 1e-7) {
        variables["M_s"] = 1.989e30;  // Standard solar mass
    }
    
    // Reset M_sun if corrupted
    if (variables["M_sun"] <= 0 || variables["M_sun"] > 1e32) {
        variables["M_sun"] = 1.989e30;
    }
    
    // Ensure densities are positive
    if (variables["rho_vac_UA"] <= 0) {
        variables["rho_vac_UA"] = 7.09e-36;
    }
    if (variables["rho_vac_SCm"] <= 0) {
        variables["rho_vac_SCm"] = 7.09e-37;
    }
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
    
    // Ensure couplings are reasonable
    if (variables["k_1"] <= 0 || variables["k_1"] > 10.0) {
        variables["k_1"] = 1.5;
    }
    if (variables["k_2"] <= 0 || variables["k_2"] > 10.0) {
        variables["k_2"] = 1.2;
    }
    
    // Correct swirl factor
    variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
    
    // Correct mass density if exists and out of range
    if (variables.find("mass_density_kg_m3") != variables.end()) {
        if (variables["mass_density_kg_m3"] <= 0 || variables["mass_density_kg_m3"] > 1e18) {
            variables["mass_density_kg_m3"] = 1408.0;  // Solar density
        }
    }
    
    // Correct surface gravity if exists
    if (variables.find("surface_gravity_m_s2") != variables.end()) {
        if (variables["surface_gravity_m_s2"] <= 0 || variables["surface_gravity_m_s2"] > 1e15) {
            variables["surface_gravity_m_s2"] = 274.0;  // Solar surface gravity
        }
    }
}

// Example usage in base program (snippet)
int main() {
    StellarMassModule module;
    std::cout << "===== Stellar/Planetary Mass Module Enhanced Demonstration =====\n\n";
    
    // Step 1: Report initial state (solar mass)
    std::cout << "STEP 1: Initial Configuration (M_s = 1 M_sun)\n";
    std::cout << module.generateReport() << "\n";
    
    // Step 2: Compute M_s/r² at various radii
    std::cout << "STEP 2: M_s/r² at Various Radii\n";
    std::vector<double> test_radii = {0.1, 0.5, 1.0, 5.0, 10.0, 50.0, 100.0};
    for (double r_AU : test_radii) {
        double r_m = r_AU * 1.496e11;
        double m_over_r2 = module.computeM_sOverR2(r_m);
        std::cout << "  r=" << std::fixed << std::setprecision(1) << r_AU << " AU: ";
        std::cout << "M_s/r² = " << std::scientific << m_over_r2 << " kg/m²\n";
    }
    std::cout << "\n";
    
    // Step 3: Compute U_g1 and U_g2 at key radii
    std::cout << "STEP 3: U_g1 and U_g2 at Key Radii\n";
    for (double r_AU : test_radii) {
        double r_m = r_AU * 1.496e11;
        double U_g1 = module.computeU_g1(r_m);
        double U_g2 = module.computeU_g2(r_m);
        std::cout << "  r=" << std::fixed << std::setprecision(1) << r_AU << " AU: ";
        std::cout << "U_g1=" << std::scientific << U_g1 << " J/m³";
        if (U_g2 > 0) {
            std::cout << ", U_g2=" << U_g2 << " J/m³";
        } else {
            std::cout << ", U_g2=0 (r < R_b)";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
    
    // Step 4: Save solar mass state
    std::cout << "STEP 4: Save Solar Mass State\n";
    module.saveState("solar_mass_1Msun");
    std::cout << "State saved as 'solar_mass_1Msun'\n\n";
    
    // Step 5: Test Jupiter mass
    std::cout << "STEP 5: Test Jupiter Mass (M_J ~ 0.000955 M_sun)\n";
    module.optimizeForMetric("jupiter_mass");
    double m_jupiter = module.getVariable("M_s");
    double m_jupiter_ratio = module.computeM_sInMsun();
    double U_g2_jupiter = module.computeU_g2(1.496e13);
    std::cout << "Jupiter mass:\n";
    std::cout << "  M_s = " << std::scientific << m_jupiter << " kg (";
    std::cout << std::fixed << std::setprecision(6) << m_jupiter_ratio << " M_sun)\n";
    std::cout << "  U_g2 at R_b = " << std::scientific << U_g2_jupiter << " J/m³\n\n";
    
    // Step 6: Test Earth mass
    std::cout << "STEP 6: Test Earth Mass (M_E ~ 3e-6 M_sun)\n";
    module.restoreState("solar_mass_1Msun");
    module.optimizeForMetric("earth_mass");
    double m_earth = module.getVariable("M_s");
    double m_earth_ratio = module.computeM_sInMsun();
    double U_g2_earth = module.computeU_g2(1.496e13);
    std::cout << "Earth mass:\n";
    std::cout << "  M_s = " << std::scientific << m_earth << " kg (";
    std::cout << std::scientific << std::setprecision(2) << m_earth_ratio << " M_sun)\n";
    std::cout << "  U_g2 at R_b = " << std::scientific << U_g2_earth << " J/m³\n\n";
    
    // Step 7: Test low-mass star (0.1 M_sun, red dwarf)
    std::cout << "STEP 7: Test Low-Mass Star (0.1 M_sun, red dwarf)\n";
    module.restoreState("solar_mass_1Msun");
    module.optimizeForMetric("low_mass_star");
    double m_low = module.getVariable("M_s");
    double m_low_ratio = module.computeM_sInMsun();
    std::cout << "Low-mass star:\n";
    std::cout << "  M_s = " << std::scientific << m_low << " kg (";
    std::cout << std::fixed << std::setprecision(2) << m_low_ratio << " M_sun)\n\n";
    
    // Step 8: Test high-mass star (10 M_sun)
    std::cout << "STEP 8: Test High-Mass Star (10 M_sun)\n";
    module.restoreState("solar_mass_1Msun");
    module.optimizeForMetric("high_mass_star");
    double m_high = module.getVariable("M_s");
    double m_high_ratio = module.computeM_sInMsun();
    std::cout << "High-mass star:\n";
    std::cout << "  M_s = " << std::scientific << m_high << " kg (";
    std::cout << std::fixed << std::setprecision(1) << m_high_ratio << " M_sun)\n\n";
    
    // Step 9: Test supermassive star (100 M_sun)
    std::cout << "STEP 9: Test Supermassive Star (100 M_sun)\n";
    module.restoreState("solar_mass_1Msun");
    module.optimizeForMetric("supermassive_star");
    double m_super = module.getVariable("M_s");
    double m_super_ratio = module.computeM_sInMsun();
    std::cout << "Supermassive star:\n";
    std::cout << "  M_s = " << std::scientific << m_super << " kg (";
    std::cout << std::fixed << std::setprecision(0) << m_super_ratio << " M_sun)\n\n";
    
    // Step 10: Test neutron star (1.4 M_sun)
    std::cout << "STEP 10: Test Neutron Star (1.4 M_sun)\n";
    module.restoreState("solar_mass_1Msun");
    module.optimizeForMetric("neutron_star");
    double m_ns = module.getVariable("M_s");
    double m_ns_ratio = module.computeM_sInMsun();
    std::cout << "Neutron star:\n";
    std::cout << "  M_s = " << std::scientific << m_ns << " kg (";
    std::cout << std::fixed << std::setprecision(1) << m_ns_ratio << " M_sun)\n\n";
    
    // Step 11: Test Sagittarius A* (4.15e6 M_sun)
    std::cout << "STEP 11: Test Sagittarius A* (4.15e6 M_sun, supermassive black hole)\n";
    module.restoreState("solar_mass_1Msun");
    module.optimizeForMetric("sgr_a_star");
    double m_sgra = module.getVariable("M_s");
    double m_sgra_ratio = module.computeM_sInMsun();
    std::cout << "Sagittarius A*:\n";
    std::cout << "  M_s = " << std::scientific << m_sgra << " kg (";
    std::cout << std::scientific << std::setprecision(2) << m_sgra_ratio << " M_sun)\n\n";
    
    // Step 12: Expand mass scale
    std::cout << "STEP 12: Expand Mass Scale (M_s x2, density x1.5)\n";
    module.restoreState("solar_mass_1Msun");
    module.expandMassScale(2.0, 1.5);
    double m_expanded = module.getVariable("M_s");
    double density_expanded = module.getVariable("mass_density_kg_m3");
    std::cout << "Expanded mass:\n";
    std::cout << "  M_s = " << std::scientific << m_expanded << " kg (";
    std::cout << std::fixed << std::setprecision(1) << (m_expanded / 1.989e30) << " M_sun)\n";
    std::cout << "  Mass density = " << std::fixed << std::setprecision(0) << density_expanded << " kg/m³\n\n";
    
    // Step 13: Expand gravity scale
    std::cout << "STEP 13: Expand Gravity Scale (coupling x1.5, field x2)\n";
    module.restoreState("solar_mass_1Msun");
    module.expandGravityScale(1.5, 2.0);
    double k1_expanded = module.getVariable("k_1");
    double k2_expanded = module.getVariable("k_2");
    double field_strength = module.getVariable("gravity_field_strength");
    std::cout << "Expanded gravity:\n";
    std::cout << "  k_1 = " << std::fixed << std::setprecision(2) << k1_expanded << " (1.5x)\n";
    std::cout << "  k_2 = " << std::fixed << std::setprecision(2) << k2_expanded << " (1.5x)\n";
    std::cout << "  Field strength = " << std::fixed << std::setprecision(1) << field_strength << "x\n\n";
    
    // Step 14: Expand spatial scale
    std::cout << "STEP 14: Expand Spatial Scale (radius x2, boundary x1.5)\n";
    module.restoreState("solar_mass_1Msun");
    module.expandSpatialScale(2.0, 1.5);
    double r_expanded = module.getVariable("r");
    double rb_expanded = module.getVariable("R_b");
    double stellar_radius = module.getVariable("stellar_radius_m");
    std::cout << "Expanded spatial:\n";
    std::cout << "  Current r = " << std::scientific << r_expanded << " m (";
    std::cout << std::fixed << std::setprecision(0) << (r_expanded / 1.496e11) << " AU)\n";
    std::cout << std::scientific;
    std::cout << "  R_b = " << rb_expanded << " m (";
    std::cout << std::fixed << std::setprecision(0) << (rb_expanded / 1.496e11) << " AU)\n";
    std::cout << std::scientific;
    std::cout << "  Stellar radius = " << stellar_radius << " m (";
    std::cout << std::fixed << std::setprecision(1) << (stellar_radius / 6.96e8) << " R_sun)\n\n";
    
    // Step 15: Sensitivity analysis
    std::cout << "STEP 15: Sensitivity Analysis (at R_b)\n";
    module.restoreState("solar_mass_1Msun");
    std::vector<std::string> params = {"M_s", "k_1", "k_2", "rho_vac_UA", "delta_sw", "v_sw"};
    auto sensitivities = module.sensitivityAnalysis(params);
    for (const auto& pair : sensitivities) {
        std::cout << "  ∂U_g2/∂" << pair.first << " ≈ " << std::scientific << pair.second << " (normalized)\n";
    }
    std::cout << "\n";
    
    // Step 16: Generate variations
    std::cout << "STEP 16: Generate Parameter Variations (5 variants, ±10%)\n";
    auto variations = module.generateVariations(5, 0.1);
    for (int i = 0; i < variations.size(); ++i) {
        double var_m = variations[i]["M_s"];
        double var_ratio = var_m / variations[i]["M_sun"];
        std::cout << "  Variant " << (i+1) << ": M_s=" << std::scientific << var_m 
                  << " kg (" << std::fixed << std::setprecision(3) << var_ratio << " M_sun)\n";
    }
    std::cout << "\n";
    
    // Step 17: Auto-refine to target mass
    std::cout << "STEP 17: Auto-Refine to Target Mass (2.5 M_sun)\n";
    module.restoreState("solar_mass_1Msun");
    module.autoRefineParameters("M_sun_multiple", 2.5);
    double refined_m = module.getVariable("M_s");
    double refined_ratio = module.computeM_sInMsun();
    std::cout << "Refined M_s = " << std::scientific << refined_m << " kg (";
    std::cout << std::fixed << std::setprecision(1) << refined_ratio << " M_sun, target: 2.5)\n\n";
    
    // Step 18: Auto-refine to target M_s/r²
    std::cout << "STEP 18: Auto-Refine to Target M_s/r² (1e4 kg/m²)\n";
    module.restoreState("solar_mass_1Msun");
    module.autoRefineParameters("M_s_over_r2", 1e4);
    double refined_m_r2 = module.getVariable("M_s");
    double refined_m_over_r2 = module.computeM_sOverR2(module.getVariable("r"));
    std::cout << "Refined parameters:\n";
    std::cout << "  M_s = " << std::scientific << refined_m_r2 << " kg\n";
    std::cout << "  M_s/r² = " << refined_m_over_r2 << " kg/m² (target: 1e4)\n\n";
    
    // Step 19: Auto-refine to target U_g2
    std::cout << "STEP 19: Auto-Refine to Target U_g2 (1e54 J/m³ at R_b)\n";
    module.restoreState("solar_mass_1Msun");
    module.autoRefineParameters("U_g2_at_Rb", 1e54);
    double refined_k2 = module.getVariable("k_2");
    double refined_U_g2 = module.computeU_g2(module.getVariable("R_b"));
    std::cout << "Refined parameters:\n";
    std::cout << "  k_2 = " << std::fixed << std::setprecision(2) << refined_k2 << "\n";
    std::cout << "  U_g2 at R_b = " << std::scientific << refined_U_g2 << " J/m³ (target: 1e54)\n\n";
    
    // Step 20: Calibrate to observations
    std::cout << "STEP 20: Calibrate to Observational Data\n";
    module.restoreState("solar_mass_1Msun");
    std::map<std::string, double> observations;
    observations["M_s"] = 2.2e30;         // Observed massive star
    observations["k_1"] = 1.6;            // Observed coupling
    observations["rho_vac_UA"] = 8.0e-36; // Observed vacuum density
    module.calibrateToObservations(observations);
    std::cout << "Calibrated parameters:\n";
    std::cout << "  M_s = " << std::scientific << module.getVariable("M_s") << " kg (";
    std::cout << std::fixed << std::setprecision(2) << module.computeM_sInMsun() << " M_sun)\n";
    std::cout << "  k_1 = " << std::fixed << std::setprecision(2) << module.getVariable("k_1") << "\n";
    std::cout << "  ρ_vac,[UA] = " << std::scientific << module.getVariable("rho_vac_UA") << " J/m³\n\n";
    
    // Step 21: Mutate parameters
    std::cout << "STEP 21: Mutate Parameters (5% random variation)\n";
    module.mutateParameters(0.05);
    double mutated_m = module.getVariable("M_s");
    double mutated_ratio = module.computeM_sInMsun();
    std::cout << "Mutated parameters:\n";
    std::cout << "  M_s = " << std::scientific << mutated_m << " kg (";
    std::cout << std::fixed << std::setprecision(3) << mutated_ratio << " M_sun)\n\n";
    
    // Step 22: Validate consistency
    std::cout << "STEP 22: Validate Consistency\n";
    bool valid = module.validateConsistency();
    std::cout << "Consistency check: " << (valid ? "PASSED" : "FAILED") << "\n\n";
    
    // Step 23: Introduce anomaly and auto-correct
    std::cout << "STEP 23: Introduce Anomaly and Auto-Correct\n";
    module.createVariable("M_s_anomaly", -1e30);
    module.removeVariable("M_s");
    module.createVariable("M_s", -1e30);  // Invalid negative M_s
    std::cout << "Introduced invalid M_s = " << std::scientific << module.getVariable("M_s") << " kg\n";
    module.autoCorrectAnomalies();
    std::cout << "Auto-corrected M_s = " << module.getVariable("M_s") << " kg (";
    std::cout << std::fixed << std::setprecision(1) << module.computeM_sInMsun() << " M_sun)\n";
    bool valid_after = module.validateConsistency();
    std::cout << "Consistency after correction: " << (valid_after ? "PASSED" : "FAILED") << "\n\n";
    
    // Step 24: Test high and low density regimes
    std::cout << "STEP 24: Test High and Low Density Regimes\n";
    module.restoreState("solar_mass_1Msun");
    module.optimizeForMetric("high_density");
    double high_density = module.getVariable("mass_density_kg_m3");
    std::cout << "High density regime:\n";
    std::cout << "  Mass density = " << std::fixed << std::setprecision(0) << high_density 
              << " kg/m³ (white dwarf regime)\n";
    
    module.restoreState("solar_mass_1Msun");
    module.optimizeForMetric("low_density");
    double low_density = module.getVariable("mass_density_kg_m3");
    std::cout << "Low density regime:\n";
    std::cout << "  Mass density = " << std::fixed << std::setprecision(0) << low_density 
              << " kg/m³ (gas giant regime)\n\n";
    
    // Step 25: List saved states and export
    std::cout << "STEP 25: List Saved States and Export Final State\n";
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

StellarMassModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeM_s, computeM_sInMsun, computeM_sOverR2, computeU_g1, computeU_g2) are clear, concise, and variable - driven.
- Automatic recalculation of derived variables(rho_sum, swirl_factor) when dependencies change.
- Output and debugging functions(printVariables, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.
- Supports both internal and external gravity calculations(U_g1, U_g2) with mass scaling.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in stellar / planetary mass modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.