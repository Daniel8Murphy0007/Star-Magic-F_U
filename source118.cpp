// StellarRotationModule.h
// Modular C++ implementation of the Stellar/Planetary Rotation Rate (?_s) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes ?_s=2.5e-6 rad/s (~29-day Sun period); scales ?_s(t) in U_g3 cos(?_s t ?) and U_i ?_s cos(? t_n).
// Pluggable: #include "StellarRotationModule.h"
// StellarRotationModule mod; mod.computeU_g3(0.0); mod.updateVariable("omega_s", new_value);
// Variables in std::map; example for Sun at t=0, t_n=0; U_g3 ?1.8e49 J/m�, U_i ?1.38e-47 J/m�.
// Approximations: cos(? t_n)=1; f_TRZ=0.1; ?_i=1.0; ?_vac sum=7.80e-36 J/m�.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef STELLAR_ROTATION_MODULE_H
#define STELLAR_ROTATION_MODULE_H

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

class StellarRotationModule {
private:
    std::map<std::string, double> variables;
    double computeOmega_s_t(double t);  // ?_s(t), simplified constant
    double computeU_g3(double t);
    double computeU_i(double t, double t_n);

public:
    // Constructor: Initialize with framework defaults (Sun)
    StellarRotationModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeOmega_s();  // 2.5e-6 rad/s
    double computeOmega_s_t(double t);  // ?_s(t) (rad/s)
    double computePeriod_days();  // ~29 days
    double computeU_g3(double t);  // U_g3 example (J/m^3)
    double computeU_i(double t, double t_n);  // U_i example (J/m^3)

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
    void expandParameterSpace(double rotation_scale, double field_scale, double time_scale);
    void expandRotationScale(double omega_factor, double period_factor);     // ω_s and rotation period
    void expandFieldScale(double magnetic_factor, double coupling_factor);   // B_j and field coupling
    void expandTimeScale(double temporal_factor, double oscillation_factor); // Time evolution and oscillations

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

#endif // STELLAR_ROTATION_MODULE_H

// StellarRotationModule.cpp
#include "StellarRotationModule.h"

// Constructor: Set framework defaults (Sun at t=0)
StellarRotationModule::StellarRotationModule() {
    // Universal constants
    variables["omega_s"] = 2.5e-6;                  // rad/s
    variables["k_3"] = 1.8;                         // Coupling U_g3
    variables["B_j"] = 1e3;                         // T
    variables["P_core"] = 1.0;                      // Unitless
    variables["E_react"] = 1e46;                    // J
    variables["lambda_i"] = 1.0;                    // Unitless U_i
    variables["rho_vac_SCm"] = 7.09e-37;            // J/m^3
    variables["rho_vac_UA"] = 7.09e-36;             // J/m^3
    variables["f_TRZ"] = 0.1;                       // Unitless
    variables["pi"] = 3.141592653589793;
    variables["t"] = 0.0;                           // s
    variables["t_n"] = 0.0;                         // s

    // Derived
    variables["rho_sum"] = variables["rho_vac_SCm"] + variables["rho_vac_UA"];
    variables["day_to_s"] = 86400.0;                // s/day
}

// Update variable
void StellarRotationModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
        if (name == "rho_vac_SCm" || name == "rho_vac_UA") {
            variables["rho_sum"] = variables["rho_vac_SCm"] + variables["rho_vac_UA"];
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void StellarRotationModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        if (name == "rho_vac_SCm" || name == "rho_vac_UA") {
            variables["rho_sum"] = variables["rho_vac_SCm"] + variables["rho_vac_UA"];
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void StellarRotationModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute ?_s (rad/s)
double StellarRotationModule::computeOmega_s() {
    return variables["omega_s"];
}

// ?_s(t), simplified as constant (no t dep in example)
double StellarRotationModule::computeOmega_s_t(double t) {
    variables["t"] = t;
    return computeOmega_s();  // Constant for Sun
}

// Rotation period in days
double StellarRotationModule::computePeriod_days() {
    double period_s = 2.0 * M_PI / computeOmega_s();
    return period_s / variables["day_to_s"];
}

// U_g3 example
double StellarRotationModule::computeU_g3(double t) {
    double k_3 = variables["k_3"];
    double b_j = variables["B_j"];
    double cos_term = std::cos(variables["omega_s"] * t * variables["pi"]);
    double p_core = variables["P_core"];
    double e_react = variables["E_react"];
    return k_3 * b_j * cos_term * p_core * e_react;
}

// U_i example
double StellarRotationModule::computeU_i(double t, double t_n) {
    double lambda_i = variables["lambda_i"];
    double rho_sc = variables["rho_vac_SCm"];
    double rho_ua = variables["rho_vac_UA"];
    double omega_s_t = computeOmega_s_t(t);
    double cos_pi_tn = std::cos(variables["pi"] * t_n);
    double trz_factor = 1.0 + variables["f_TRZ"];
    return lambda_i * rho_sc * rho_ua * omega_s_t * cos_pi_tn * trz_factor;
}

// Equation text
std::string StellarRotationModule::getEquationText() {
    return "U_g3 = k_3 * ? B_j * cos(?_s(t) t ?) * P_core * E_react\n"
           "U_i = ?_i * ?_vac,[SCm] * ?_vac,[UA] * ?_s(t) * cos(? t_n) * (1 + f_TRZ)\n"
           "Where ?_s = 2.5e-6 rad/s (~29-day Sun equatorial rotation);\n"
           "Scales rotational oscillations/inertia.\n"
           "Example t=0, t_n=0: U_g3 ?1.8e49 J/m�; U_i ?1.38e-47 J/m�.\n"
           "Role: Introduces spin in disk gravity/inertia; stellar/planetary dynamics.\n"
           "UQFF: Rotational effects in nebulae/disks/formation/mergers.";
}

// Print variables
void StellarRotationModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// ===== ENHANCED METHODS IMPLEMENTATION =====

namespace stellar_rotation_saved_states {
    static std::map<std::string, std::map<std::string, double>> saved_states;
}

// Variable Management
void StellarRotationModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void StellarRotationModule::removeVariable(const std::string& name) {
    variables.erase(name);
}

void StellarRotationModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
    }
}

std::vector<std::string> StellarRotationModule::listVariables() const {
    std::vector<std::string> names;
    for (const auto& pair : variables) {
        names.push_back(pair.first);
    }
    return names;
}

std::string StellarRotationModule::getSystemName() const {
    return "Stellar_Planetary_Rotation_UQFF";
}

// Batch Operations
void StellarRotationModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = func(variables[name]);
        }
    }
    // Update derived
    variables["rho_sum"] = variables["rho_vac_SCm"] + variables["rho_vac_UA"];
}

void StellarRotationModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// Self-Expansion: Domain-Specific Scales
void StellarRotationModule::expandParameterSpace(double rotation_scale, double field_scale, double time_scale) {
    variables["omega_s"] *= rotation_scale;
    variables["B_j"] *= field_scale;
    variables["k_3"] *= field_scale;
    variables["t"] *= time_scale;
    variables["t_n"] *= time_scale;
}

void StellarRotationModule::expandRotationScale(double omega_factor, double period_factor) {
    variables["omega_s"] *= omega_factor;
    
    // Store rotation characteristics for advanced modeling
    if (variables.find("rotation_period_days") == variables.end()) {
        // Solar rotation: ~27 days (equator), ~29 days (typical), ~31 days (poles)
        variables["rotation_period_days"] = 29.0;
    }
    // Period scales inversely with omega
    variables["rotation_period_days"] /= omega_factor;
    variables["rotation_period_days"] *= period_factor;
    
    // Differential rotation modeling
    if (variables.find("differential_rotation_factor") == variables.end()) {
        // Sun: poles rotate ~25% slower than equator
        variables["differential_rotation_factor"] = 1.25;
    }
    variables["differential_rotation_factor"] *= period_factor;
    
    // Rossby number (rotation vs convection timescale)
    if (variables.find("rossby_number") == variables.end()) {
        // Ro ~ τ_conv / τ_rot; Sun: Ro ~ 1-2
        variables["rossby_number"] = 1.5;
    }
    variables["rossby_number"] /= omega_factor;  // Faster rotation → smaller Ro
}

void StellarRotationModule::expandFieldScale(double magnetic_factor, double coupling_factor) {
    variables["B_j"] *= magnetic_factor;
    variables["k_3"] *= coupling_factor;
    
    // Magnetic field characteristics
    if (variables.find("magnetic_field_strength_T") == variables.end()) {
        // Solar magnetic field: ~1000 T (internal), ~1e-4 T (surface)
        variables["magnetic_field_strength_T"] = 1000.0;
    }
    variables["magnetic_field_strength_T"] *= magnetic_factor;
    
    // Dynamo efficiency
    if (variables.find("dynamo_efficiency") == variables.end()) {
        variables["dynamo_efficiency"] = 1.0;
    }
    variables["dynamo_efficiency"] *= coupling_factor;
    
    // Magnetic energy density
    if (variables.find("magnetic_energy_density_J_m3") == variables.end()) {
        // B²/(2μ₀) ~ (1000 T)² / (2 × 4π×10⁻⁷) ~ 4e11 J/m³
        variables["magnetic_energy_density_J_m3"] = 4e11;
    }
    variables["magnetic_energy_density_J_m3"] *= (magnetic_factor * magnetic_factor);
}

void StellarRotationModule::expandTimeScale(double temporal_factor, double oscillation_factor) {
    variables["t"] *= temporal_factor;
    variables["t_n"] *= temporal_factor;
    
    // Oscillation characteristics
    if (variables.find("oscillation_amplitude") == variables.end()) {
        variables["oscillation_amplitude"] = 1.0;
    }
    variables["oscillation_amplitude"] *= oscillation_factor;
    
    // Phase evolution
    if (variables.find("phase_offset_rad") == variables.end()) {
        variables["phase_offset_rad"] = 0.0;
    }
    // Phase accumulates with time scaling
    variables["phase_offset_rad"] += (temporal_factor - 1.0) * variables["omega_s"] * variables["t"];
    
    // Rotation cycles completed
    if (variables.find("rotation_cycles") == variables.end()) {
        variables["rotation_cycles"] = 0.0;
    }
    variables["rotation_cycles"] += (temporal_factor - 1.0) * variables["omega_s"] * variables["t"] / (2.0 * M_PI);
}

// Self-Refinement
void StellarRotationModule::autoRefineParameters(const std::string& target, double goal) {
    if (target == "omega_s") {
        variables["omega_s"] = goal;
    } else if (target == "period_days") {
        // ω = 2π / (period × 86400)
        variables["omega_s"] = 2.0 * M_PI / (goal * variables["day_to_s"]);
        if (variables.find("rotation_period_days") != variables.end()) {
            variables["rotation_period_days"] = goal;
        }
    } else if (target == "rotation_period_days") {
        variables["omega_s"] = 2.0 * M_PI / (goal * variables["day_to_s"]);
        if (variables.find("rotation_period_days") == variables.end()) {
            variables["rotation_period_days"] = goal;
        } else {
            variables["rotation_period_days"] = goal;
        }
    } else if (target == "frequency_Hz") {
        // ω = 2π f
        variables["omega_s"] = 2.0 * M_PI * goal;
    } else if (target == "angular_velocity_deg_s") {
        // Convert deg/s to rad/s
        variables["omega_s"] = goal * M_PI / 180.0;
    } else if (target == "B_j") {
        variables["B_j"] = goal;
        if (variables.find("magnetic_field_strength_T") != variables.end()) {
            variables["magnetic_field_strength_T"] = goal;
        }
    } else if (target == "magnetic_field_strength_T") {
        variables["B_j"] = goal;
        if (variables.find("magnetic_field_strength_T") == variables.end()) {
            variables["magnetic_field_strength_T"] = goal;
        } else {
            variables["magnetic_field_strength_T"] = goal;
        }
    } else if (target == "U_g3_at_t0") {
        // Target specific U_g3 at t=0 by adjusting k_3
        double current_U_g3 = computeU_g3(0.0);
        if (current_U_g3 > 0) {
            variables["k_3"] *= (goal / current_U_g3);
        }
    } else if (target == "U_i_at_t0") {
        // Target specific U_i at t=0, t_n=0 by adjusting lambda_i
        double current_U_i = computeU_i(0.0, 0.0);
        if (current_U_i > 0) {
            variables["lambda_i"] *= (goal / current_U_i);
        }
    } else if (target == "rossby_number") {
        if (variables.find("rossby_number") == variables.end()) {
            variables["rossby_number"] = 1.5;
        }
        variables["rossby_number"] = goal;
    } else if (target == "differential_rotation_factor") {
        if (variables.find("differential_rotation_factor") == variables.end()) {
            variables["differential_rotation_factor"] = 1.25;
        }
        variables["differential_rotation_factor"] = goal;
    }
}

void StellarRotationModule::calibrateToObservations(const std::map<std::string, double>& observations) {
    for (const auto& obs : observations) {
        if (variables.find(obs.first) != variables.end()) {
            variables[obs.first] = obs.second;
            if (obs.first == "rho_vac_SCm" || obs.first == "rho_vac_UA") {
                variables["rho_sum"] = variables["rho_vac_SCm"] + variables["rho_vac_UA"];
            }
        }
    }
}

void StellarRotationModule::optimizeForMetric(const std::string& metric) {
    if (metric == "solar_rotation") {
        // Sun equatorial: ~27 days, typical: ~29 days
        variables["omega_s"] = 2.5e-6;  // rad/s, ~29 days
        if (variables.find("rotation_period_days") != variables.end()) {
            variables["rotation_period_days"] = 29.0;
        }
    } else if (metric == "fast_rotation") {
        // Fast rotators: 1-5 days
        double period_s = 3.0 * variables["day_to_s"];  // 3 days
        variables["omega_s"] = 2.0 * M_PI / period_s;
        if (variables.find("rotation_period_days") != variables.end()) {
            variables["rotation_period_days"] = 3.0;
        }
    } else if (metric == "very_fast_rotation") {
        // Very fast: hours to 1 day (young stars, pulsars)
        double period_s = 0.5 * variables["day_to_s"];  // 12 hours
        variables["omega_s"] = 2.0 * M_PI / period_s;
        if (variables.find("rotation_period_days") != variables.end()) {
            variables["rotation_period_days"] = 0.5;
        }
    } else if (metric == "slow_rotation") {
        // Slow rotators: 50-100 days
        double period_s = 70.0 * variables["day_to_s"];  // 70 days
        variables["omega_s"] = 2.0 * M_PI / period_s;
        if (variables.find("rotation_period_days") != variables.end()) {
            variables["rotation_period_days"] = 70.0;
        }
    } else if (metric == "very_slow_rotation") {
        // Very slow: months to years
        double period_s = 200.0 * variables["day_to_s"];  // ~200 days
        variables["omega_s"] = 2.0 * M_PI / period_s;
        if (variables.find("rotation_period_days") != variables.end()) {
            variables["rotation_period_days"] = 200.0;
        }
    } else if (metric == "jupiter_rotation") {
        // Jupiter: ~10 hours (very fast)
        double period_s = (10.0 / 24.0) * variables["day_to_s"];
        variables["omega_s"] = 2.0 * M_PI / period_s;
        if (variables.find("rotation_period_days") != variables.end()) {
            variables["rotation_period_days"] = 10.0 / 24.0;
        }
    } else if (metric == "earth_rotation") {
        // Earth: 24 hours (1 day)
        double period_s = variables["day_to_s"];
        variables["omega_s"] = 2.0 * M_PI / period_s;
        if (variables.find("rotation_period_days") != variables.end()) {
            variables["rotation_period_days"] = 1.0;
        }
    } else if (metric == "venus_rotation") {
        // Venus: ~243 days (retrograde, very slow)
        double period_s = 243.0 * variables["day_to_s"];
        variables["omega_s"] = 2.0 * M_PI / period_s;
        if (variables.find("rotation_period_days") != variables.end()) {
            variables["rotation_period_days"] = 243.0;
        }
    } else if (metric == "pulsar_rotation") {
        // Millisecond pulsar: ~0.001 s period
        double period_s = 0.001;
        variables["omega_s"] = 2.0 * M_PI / period_s;
        if (variables.find("rotation_period_days") != variables.end()) {
            variables["rotation_period_days"] = 0.001 / variables["day_to_s"];
        }
    } else if (metric == "high_magnetic_field") {
        // Strong internal field
        variables["B_j"] = 1e4;  // 10,000 T
        if (variables.find("magnetic_field_strength_T") != variables.end()) {
            variables["magnetic_field_strength_T"] = 1e4;
        }
    } else if (metric == "low_magnetic_field") {
        // Weak field
        variables["B_j"] = 100.0;  // 100 T
        if (variables.find("magnetic_field_strength_T") != variables.end()) {
            variables["magnetic_field_strength_T"] = 100.0;
        }
    } else if (metric == "high_rossby") {
        // Slow rotation regime (weak dynamo)
        if (variables.find("rossby_number") == variables.end()) {
            variables["rossby_number"] = 1.5;
        }
        variables["rossby_number"] = 5.0;
    } else if (metric == "low_rossby") {
        // Fast rotation regime (strong dynamo)
        if (variables.find("rossby_number") == variables.end()) {
            variables["rossby_number"] = 1.5;
        }
        variables["rossby_number"] = 0.1;
    }
}

// Parameter Exploration
std::vector<std::map<std::string, double>> StellarRotationModule::generateVariations(int count, double variation_pct) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - variation_pct, 1.0 + variation_pct);
    
    for (int i = 0; i < count; ++i) {
        std::map<std::string, double> variant = variables;
        for (auto& pair : variant) {
            if (pair.first != "rho_sum" && pair.first != "day_to_s" && pair.first != "pi") {
                pair.second *= dis(gen);
            }
        }
        // Recalculate derived
        variant["rho_sum"] = variant["rho_vac_SCm"] + variant["rho_vac_UA"];
        variations.push_back(variant);
    }
    return variations;
}

// Adaptive Evolution
void StellarRotationModule::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - mutation_rate, 1.0 + mutation_rate);
    
    for (auto& pair : variables) {
        if (pair.first != "rho_sum" && pair.first != "day_to_s" && pair.first != "pi") {
            pair.second *= dis(gen);
        }
    }
    // Recalculate derived
    variables["rho_sum"] = variables["rho_vac_SCm"] + variables["rho_vac_UA"];
}

void StellarRotationModule::evolveSystem(int generations, std::function<double()> fitness_func) {
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
void StellarRotationModule::saveState(const std::string& label) {
    stellar_rotation_saved_states::saved_states[label] = variables;
}

void StellarRotationModule::restoreState(const std::string& label) {
    if (stellar_rotation_saved_states::saved_states.find(label) != stellar_rotation_saved_states::saved_states.end()) {
        variables = stellar_rotation_saved_states::saved_states[label];
    }
}

std::vector<std::string> StellarRotationModule::listSavedStates() const {
    std::vector<std::string> labels;
    for (const auto& pair : stellar_rotation_saved_states::saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string StellarRotationModule::exportState() const {
    std::ostringstream oss;
    oss << "StellarRotation_State_Export:\n";
    for (const auto& pair : variables) {
        oss << pair.first << "=" << std::scientific << pair.second << "\n";
    }
    return oss.str();
}

// System Analysis
std::map<std::string, double> StellarRotationModule::sensitivityAnalysis(const std::vector<std::string>& params) {
    std::map<std::string, double> sensitivities;
    double t_test = 0.0;  // Test at t=0
    double baseline_g3 = computeU_g3(t_test);
    double baseline_i = computeU_i(t_test, 0.0);
    double baseline = baseline_g3;  // Use U_g3 as primary metric
    
    for (const auto& param : params) {
        if (variables.find(param) != variables.end() && 
            param != "rho_sum" && param != "day_to_s" && param != "pi") {
            double original = variables[param];
            variables[param] = original * 1.01;
            
            if (param == "rho_vac_SCm" || param == "rho_vac_UA") {
                variables["rho_sum"] = variables["rho_vac_SCm"] + variables["rho_vac_UA"];
            }
            
            double perturbed = computeU_g3(t_test);
            if (baseline > 0) {
                sensitivities[param] = (perturbed - baseline) / baseline;
            } else {
                sensitivities[param] = 0.0;
            }
            
            // Restore
            variables[param] = original;
            if (param == "rho_vac_SCm" || param == "rho_vac_UA") {
                variables["rho_sum"] = variables["rho_vac_SCm"] + variables["rho_vac_UA"];
            }
        }
    }
    return sensitivities;
}

std::string StellarRotationModule::generateReport() const {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(3);
    oss << "===== Stellar/Planetary Rotation Module Report =====\n";
    oss << "System: " << getSystemName() << "\n\n";
    
    oss << "Rotation Parameters:\n";
    oss << "  ω_s = " << variables.at("omega_s") << " rad/s\n";
    double period_days = 2.0 * M_PI / (variables.at("omega_s") * variables.at("day_to_s"));
    oss << "  Rotation period = " << std::fixed << std::setprecision(2) << period_days << " days\n";
    double frequency_Hz = variables.at("omega_s") / (2.0 * M_PI);
    oss << "  Frequency = " << std::scientific << frequency_Hz << " Hz\n";
    double angular_deg_s = variables.at("omega_s") * 180.0 / M_PI;
    oss << "  Angular velocity = " << std::scientific << angular_deg_s << " deg/s\n";
    
    if (variables.find("rotation_period_days") != variables.end()) {
        oss << "  Stored period = " << std::fixed << std::setprecision(2) 
            << variables.at("rotation_period_days") << " days\n";
    }
    if (variables.find("differential_rotation_factor") != variables.end()) {
        oss << "  Differential rotation = " << std::fixed << std::setprecision(2) 
            << variables.at("differential_rotation_factor") << "x (pole/equator)\n";
    }
    if (variables.find("rossby_number") != variables.end()) {
        oss << "  Rossby number = " << std::fixed << std::setprecision(2) 
            << variables.at("rossby_number") << " (rotation regime)\n";
    }
    oss << "\n";
    
    oss << std::scientific;
    oss << "Magnetic Field Parameters:\n";
    oss << "  B_j = " << variables.at("B_j") << " T (internal field)\n";
    if (variables.find("magnetic_field_strength_T") != variables.end()) {
        oss << "  Magnetic field strength = " << variables.at("magnetic_field_strength_T") << " T\n";
    }
    if (variables.find("magnetic_energy_density_J_m3") != variables.end()) {
        oss << "  Magnetic energy density = " << variables.at("magnetic_energy_density_J_m3") << " J/m³\n";
    }
    if (variables.find("dynamo_efficiency") != variables.end()) {
        oss << "  Dynamo efficiency = " << std::fixed << std::setprecision(2) 
            << variables.at("dynamo_efficiency") << "\n";
    }
    oss << "\n";
    
    oss << std::scientific;
    oss << "Coupling Parameters:\n";
    oss << "  k_3 = " << variables.at("k_3") << " (U_g3 coupling)\n";
    oss << "  λ_i = " << variables.at("lambda_i") << " (U_i coupling)\n";
    oss << "  P_core = " << variables.at("P_core") << "\n";
    oss << "  f_TRZ = " << variables.at("f_TRZ") << " (tachocline/radiative zone)\n\n";
    
    oss << "Vacuum Densities:\n";
    oss << "  ρ_vac,[SCm] = " << variables.at("rho_vac_SCm") << " J/m³\n";
    oss << "  ρ_vac,[UA] = " << variables.at("rho_vac_UA") << " J/m³\n";
    oss << "  ρ_sum = " << variables.at("rho_sum") << " J/m³\n\n";
    
    oss << "Energy Factor:\n";
    oss << "  E_react = " << variables.at("E_react") << " J\n\n";
    
    oss << "Time Parameters:\n";
    oss << "  Current t = " << variables.at("t") << " s\n";
    oss << "  Current t_n = " << variables.at("t_n") << " s\n";
    if (variables.find("rotation_cycles") != variables.end()) {
        oss << "  Rotation cycles completed = " << std::fixed << std::setprecision(1) 
            << variables.at("rotation_cycles") << "\n";
    }
    if (variables.find("phase_offset_rad") != variables.end()) {
        oss << "  Phase offset = " << std::scientific << variables.at("phase_offset_rad") << " rad\n";
    }
    oss << "\n";
    
    oss << "U_g3 and U_i at Various Times:\n";
    std::vector<double> test_times_days = {0.0, 7.25, 14.5, 21.75, 29.0, 58.0, 87.0};
    for (double t_days : test_times_days) {
        double t_s = t_days * variables.at("day_to_s");
        
        // U_g3
        double k_3 = variables.at("k_3");
        double b_j = variables.at("B_j");
        double cos_term = std::cos(variables.at("omega_s") * t_s * variables.at("pi"));
        double p_core = variables.at("P_core");
        double e_react = variables.at("E_react");
        double U_g3 = k_3 * b_j * cos_term * p_core * e_react;
        
        // U_i at t_n=0
        double lambda_i = variables.at("lambda_i");
        double rho_sc = variables.at("rho_vac_SCm");
        double rho_ua = variables.at("rho_vac_UA");
        double omega_s_t = variables.at("omega_s");
        double cos_pi_tn = 1.0;  // t_n=0
        double trz_factor = 1.0 + variables.at("f_TRZ");
        double U_i = lambda_i * rho_sc * rho_ua * omega_s_t * cos_pi_tn * trz_factor;
        
        oss << "  t=" << std::fixed << std::setprecision(2) << t_days << " days: ";
        oss << "U_g3=" << std::scientific << U_g3 << " J/m³, ";
        oss << "U_i=" << U_i << " J/m³\n";
    }
    oss << "\n";
    
    oss << "Physical Interpretation:\n";
    if (period_days < 1.0) {
        oss << "  Very fast rotator (<1 day, young stars/gas giants/pulsars)\n";
    } else if (period_days < 10.0) {
        oss << "  Fast rotator (1-10 days)\n";
    } else if (period_days < 40.0) {
        oss << "  Moderate rotator (10-40 days, solar-type)\n";
    } else if (period_days < 100.0) {
        oss << "  Slow rotator (40-100 days)\n";
    } else {
        oss << "  Very slow rotator (>100 days, evolved stars)\n";
    }
    
    if (variables.find("rossby_number") != variables.end()) {
        double ro = variables.at("rossby_number");
        if (ro < 0.5) {
            oss << "  Strong dynamo regime (Ro < 0.5, saturated magnetic activity)\n";
        } else if (ro < 2.0) {
            oss << "  Active dynamo regime (0.5 < Ro < 2, solar-type activity)\n";
        } else {
            oss << "  Weak dynamo regime (Ro > 2, reduced magnetic activity)\n";
        }
    }
    
    oss << "\n  Applications:\n";
    oss << "    - Stellar dynamics: Rotation drives magnetic dynamo\n";
    oss << "    - Disk gravity: U_g3 oscillates with ω_s, shapes magnetic strings\n";
    oss << "    - Inertial effects: U_i couples rotation to vacuum fields\n";
    oss << "    - Differential rotation: Latitude-dependent ω creates shear\n";
    oss << "    - Angular momentum: Rotation evolution in stellar/planetary systems\n";
    oss << "    - Activity cycles: Magnetic field modulation (e.g., 11-year solar)\n";
    oss << "    - Tachocline dynamics: Interface between radiative and convective zones\n";
    
    return oss.str();
}

bool StellarRotationModule::validateConsistency() const {
    bool valid = true;
    
    // Check ω_s is positive
    if (variables.find("omega_s") != variables.end() && variables.at("omega_s") <= 0) {
        std::cerr << "Error: omega_s <= 0 (rotation rate must be positive)\n";
        valid = false;
    }
    
    // Check rho_sum consistency
    if (variables.find("rho_sum") != variables.end()) {
        double expected = variables.at("rho_vac_SCm") + variables.at("rho_vac_UA");
        double actual = variables.at("rho_sum");
        if (std::abs(expected - actual) / expected > 1e-9) {
            std::cerr << "Error: rho_sum inconsistent (expected " << expected << ", got " << actual << ")\n";
            valid = false;
        }
    }
    
    // Check ω_s is in reasonable range
    double period_days = 2.0 * M_PI / (variables.at("omega_s") * variables.at("day_to_s"));
    if (period_days < 1e-6 || period_days > 1e4) {
        std::cerr << "Warning: Rotation period outside typical range [1e-6, 1e4] days (current: " 
                  << period_days << " days)\n";
    }
    
    // Check positive magnetic field
    if (variables.at("B_j") <= 0) {
        std::cerr << "Error: B_j <= 0 (magnetic field must be positive)\n";
        valid = false;
    }
    
    // Check positive densities
    if (variables.at("rho_vac_SCm") <= 0 || variables.at("rho_vac_UA") <= 0) {
        std::cerr << "Error: Vacuum densities must be positive\n";
        valid = false;
    }
    
    // Check positive couplings
    if (variables.at("k_3") <= 0 || variables.at("lambda_i") <= 0) {
        std::cerr << "Error: Coupling constants k_3, lambda_i must be positive\n";
        valid = false;
    }
    
    // Check Rossby number (if exists)
    if (variables.find("rossby_number") != variables.end()) {
        if (variables.at("rossby_number") <= 0 || variables.at("rossby_number") > 100.0) {
            std::cerr << "Warning: Rossby number outside typical range [0, 100] (current: " 
                      << variables.at("rossby_number") << ")\n";
        }
    }
    
    return valid;
}

void StellarRotationModule::autoCorrectAnomalies() {
    // Reset ω_s to solar value if out of range
    double period_days = 2.0 * M_PI / (variables["omega_s"] * variables["day_to_s"]);
    if (variables["omega_s"] <= 0 || period_days > 1e4 || period_days < 1e-6) {
        variables["omega_s"] = 2.5e-6;  // Standard solar rotation
        if (variables.find("rotation_period_days") != variables.end()) {
            variables["rotation_period_days"] = 29.0;
        }
    }
    
    // Ensure magnetic field is positive
    if (variables["B_j"] <= 0) {
        variables["B_j"] = 1e3;  // Standard 1000 T
        if (variables.find("magnetic_field_strength_T") != variables.end()) {
            variables["magnetic_field_strength_T"] = 1e3;
        }
    }
    
    // Ensure densities are positive
    if (variables["rho_vac_SCm"] <= 0) {
        variables["rho_vac_SCm"] = 7.09e-37;
    }
    if (variables["rho_vac_UA"] <= 0) {
        variables["rho_vac_UA"] = 7.09e-36;
    }
    variables["rho_sum"] = variables["rho_vac_SCm"] + variables["rho_vac_UA"];
    
    // Ensure couplings are reasonable
    if (variables["k_3"] <= 0 || variables["k_3"] > 10.0) {
        variables["k_3"] = 1.8;
    }
    if (variables["lambda_i"] <= 0 || variables["lambda_i"] > 10.0) {
        variables["lambda_i"] = 1.0;
    }
    
    // Correct Rossby number if exists and out of range
    if (variables.find("rossby_number") != variables.end()) {
        if (variables["rossby_number"] <= 0 || variables["rossby_number"] > 100.0) {
            variables["rossby_number"] = 1.5;  // Solar-type value
        }
    }
    
    // Correct differential rotation if exists
    if (variables.find("differential_rotation_factor") != variables.end()) {
        if (variables["differential_rotation_factor"] < 1.0 || variables["differential_rotation_factor"] > 2.0) {
            variables["differential_rotation_factor"] = 1.25;  // Solar value
        }
    }
}

// Example usage in base program (snippet)
int main() {
    StellarRotationModule module;
    std::cout << "===== Stellar/Planetary Rotation Module Enhanced Demonstration =====\n\n";
    
    // Step 1: Report initial state (solar rotation)
    std::cout << "STEP 1: Initial Configuration (ω_s = 2.5e-6 rad/s, ~29 days)\n";
    std::cout << module.generateReport() << "\n";
    
    // Step 2: Compute U_g3 at various times
    std::cout << "STEP 2: U_g3 at Various Times (one rotation period)\n";
    std::vector<double> test_times = {0.0, 7.25, 14.5, 21.75, 29.0};
    for (double t_days : test_times) {
        double t_s = t_days * 86400.0;
        double U_g3 = module.computeU_g3(t_s);
        std::cout << "  t=" << std::fixed << std::setprecision(2) << t_days << " days: ";
        std::cout << "U_g3=" << std::scientific << U_g3 << " J/m³\n";
    }
    std::cout << "\n";
    
    // Step 3: Compute U_i at various times
    std::cout << "STEP 3: U_i at Various Times (t_n=0)\n";
    for (double t_days : test_times) {
        double t_s = t_days * 86400.0;
        double U_i = module.computeU_i(t_s, 0.0);
        std::cout << "  t=" << std::fixed << std::setprecision(2) << t_days << " days: ";
        std::cout << "U_i=" << std::scientific << U_i << " J/m³\n";
    }
    std::cout << "\n";
    
    // Step 4: Save solar rotation state
    std::cout << "STEP 4: Save Solar Rotation State\n";
    module.saveState("solar_rotation_29days");
    std::cout << "State saved as 'solar_rotation_29days'\n\n";
    
    // Step 5: Test fast rotation (Jupiter, ~10 hours)
    std::cout << "STEP 5: Test Fast Rotation (Jupiter, ~10 hours)\n";
    module.optimizeForMetric("jupiter_rotation");
    double omega_jupiter = module.getVariable("omega_s");
    double period_jupiter = module.computePeriod_days();
    double U_g3_jupiter = module.computeU_g3(0.0);
    std::cout << "Jupiter rotation:\n";
    std::cout << "  ω_s = " << std::scientific << omega_jupiter << " rad/s\n";
    std::cout << "  Period = " << std::fixed << std::setprecision(2) << period_jupiter << " days (";
    std::cout << std::fixed << std::setprecision(1) << (period_jupiter * 24.0) << " hours)\n";
    std::cout << "  U_g3 at t=0 = " << std::scientific << U_g3_jupiter << " J/m³\n\n";
    
    // Step 6: Test Earth rotation (24 hours)
    std::cout << "STEP 6: Test Earth Rotation (24 hours)\n";
    module.restoreState("solar_rotation_29days");
    module.optimizeForMetric("earth_rotation");
    double omega_earth = module.getVariable("omega_s");
    double period_earth = module.computePeriod_days();
    std::cout << "Earth rotation:\n";
    std::cout << "  ω_s = " << std::scientific << omega_earth << " rad/s\n";
    std::cout << "  Period = " << std::fixed << std::setprecision(2) << period_earth << " day\n\n";
    
    // Step 7: Test very fast rotation (3 days)
    std::cout << "STEP 7: Test Very Fast Rotation (3 days, young star)\n";
    module.restoreState("solar_rotation_29days");
    module.optimizeForMetric("fast_rotation");
    double omega_fast = module.getVariable("omega_s");
    double period_fast = module.computePeriod_days();
    std::cout << "Fast rotation:\n";
    std::cout << "  ω_s = " << std::scientific << omega_fast << " rad/s\n";
    std::cout << "  Period = " << std::fixed << std::setprecision(1) << period_fast << " days\n\n";
    
    // Step 8: Test slow rotation (70 days)
    std::cout << "STEP 8: Test Slow Rotation (70 days)\n";
    module.restoreState("solar_rotation_29days");
    module.optimizeForMetric("slow_rotation");
    double omega_slow = module.getVariable("omega_s");
    double period_slow = module.computePeriod_days();
    std::cout << "Slow rotation:\n";
    std::cout << "  ω_s = " << std::scientific << omega_slow << " rad/s\n";
    std::cout << "  Period = " << std::fixed << std::setprecision(1) << period_slow << " days\n\n";
    
    // Step 9: Test Venus rotation (243 days, retrograde)
    std::cout << "STEP 9: Test Venus Rotation (243 days, very slow)\n";
    module.restoreState("solar_rotation_29days");
    module.optimizeForMetric("venus_rotation");
    double omega_venus = module.getVariable("omega_s");
    double period_venus = module.computePeriod_days();
    std::cout << "Venus rotation:\n";
    std::cout << "  ω_s = " << std::scientific << omega_venus << " rad/s\n";
    std::cout << "  Period = " << std::fixed << std::setprecision(0) << period_venus << " days\n\n";
    
    // Step 10: Test pulsar rotation (millisecond)
    std::cout << "STEP 10: Test Pulsar Rotation (1 ms period)\n";
    module.restoreState("solar_rotation_29days");
    module.optimizeForMetric("pulsar_rotation");
    double omega_pulsar = module.getVariable("omega_s");
    double period_pulsar = module.computePeriod_days();
    std::cout << "Pulsar rotation:\n";
    std::cout << "  ω_s = " << std::scientific << omega_pulsar << " rad/s\n";
    std::cout << "  Period = " << std::scientific << period_pulsar << " days (";
    std::cout << std::fixed << std::setprecision(3) << (period_pulsar * 86400.0) << " s)\n\n";
    
    // Step 11: Expand rotation scale
    std::cout << "STEP 11: Expand Rotation Scale (ω x2, period /2)\n";
    module.restoreState("solar_rotation_29days");
    module.expandRotationScale(2.0, 1.0);
    double omega_expanded = module.getVariable("omega_s");
    double period_expanded = module.computePeriod_days();
    std::cout << "Expanded rotation:\n";
    std::cout << "  ω_s = " << std::scientific << omega_expanded << " rad/s (2x)\n";
    std::cout << "  Period = " << std::fixed << std::setprecision(1) << period_expanded << " days (~half)\n\n";
    
    // Step 12: Expand field scale
    std::cout << "STEP 12: Expand Field Scale (B x3, coupling x1.5)\n";
    module.restoreState("solar_rotation_29days");
    module.expandFieldScale(3.0, 1.5);
    double b_expanded = module.getVariable("B_j");
    double k3_expanded = module.getVariable("k_3");
    double mag_energy = module.getVariable("magnetic_energy_density_J_m3");
    std::cout << "Expanded field:\n";
    std::cout << "  B_j = " << std::scientific << b_expanded << " T (3x)\n";
    std::cout << "  k_3 = " << std::fixed << std::setprecision(2) << k3_expanded << " (1.5x)\n";
    std::cout << "  Magnetic energy density = " << std::scientific << mag_energy << " J/m³ (9x)\n\n";
    
    // Step 13: Expand time scale
    std::cout << "STEP 13: Expand Time Scale (time x5, oscillation x2)\n";
    module.restoreState("solar_rotation_29days");
    module.expandTimeScale(5.0, 2.0);
    double t_expanded = module.getVariable("t");
    double cycles = module.getVariable("rotation_cycles");
    std::cout << "Expanded time:\n";
    std::cout << "  Current t = " << std::scientific << t_expanded << " s (5x)\n";
    std::cout << "  Rotation cycles = " << std::fixed << std::setprecision(1) << cycles << "\n\n";
    
    // Step 14: Test high magnetic field
    std::cout << "STEP 14: Test High Magnetic Field (10,000 T)\n";
    module.restoreState("solar_rotation_29days");
    module.optimizeForMetric("high_magnetic_field");
    double b_high = module.getVariable("B_j");
    double U_g3_high_B = module.computeU_g3(0.0);
    std::cout << "High magnetic field:\n";
    std::cout << "  B_j = " << std::scientific << b_high << " T\n";
    std::cout << "  U_g3 at t=0 = " << U_g3_high_B << " J/m³\n\n";
    
    // Step 15: Sensitivity analysis
    std::cout << "STEP 15: Sensitivity Analysis (at t=0)\n";
    module.restoreState("solar_rotation_29days");
    std::vector<std::string> params = {"omega_s", "B_j", "k_3", "lambda_i", "P_core"};
    auto sensitivities = module.sensitivityAnalysis(params);
    for (const auto& pair : sensitivities) {
        std::cout << "  ∂U_g3/∂" << pair.first << " ≈ " << std::scientific << pair.second << " (normalized)\n";
    }
    std::cout << "\n";
    
    // Step 16: Generate variations
    std::cout << "STEP 16: Generate Parameter Variations (5 variants, ±10%)\n";
    auto variations = module.generateVariations(5, 0.1);
    for (int i = 0; i < variations.size(); ++i) {
        double var_omega = variations[i]["omega_s"];
        double var_period = 2.0 * M_PI / (var_omega * 86400.0);
        std::cout << "  Variant " << (i+1) << ": ω_s=" << std::scientific << var_omega 
                  << " rad/s (period=" << std::fixed << std::setprecision(1) << var_period << " days)\n";
    }
    std::cout << "\n";
    
    // Step 17: Auto-refine to target period
    std::cout << "STEP 17: Auto-Refine to Target Period (20 days)\n";
    module.restoreState("solar_rotation_29days");
    module.autoRefineParameters("period_days", 20.0);
    double refined_omega = module.getVariable("omega_s");
    double refined_period = module.computePeriod_days();
    std::cout << "Refined parameters:\n";
    std::cout << "  ω_s = " << std::scientific << refined_omega << " rad/s\n";
    std::cout << "  Period = " << std::fixed << std::setprecision(1) << refined_period << " days (target: 20)\n\n";
    
    // Step 18: Auto-refine to target frequency
    std::cout << "STEP 18: Auto-Refine to Target Frequency (1e-5 Hz)\n";
    module.restoreState("solar_rotation_29days");
    module.autoRefineParameters("frequency_Hz", 1e-5);
    double refined_freq_omega = module.getVariable("omega_s");
    double refined_freq = refined_freq_omega / (2.0 * M_PI);
    std::cout << "Refined parameters:\n";
    std::cout << "  ω_s = " << std::scientific << refined_freq_omega << " rad/s\n";
    std::cout << "  Frequency = " << refined_freq << " Hz (target: 1e-5)\n\n";
    
    // Step 19: Auto-refine to target U_g3
    std::cout << "STEP 19: Auto-Refine to Target U_g3 (1e50 J/m³ at t=0)\n";
    module.restoreState("solar_rotation_29days");
    module.autoRefineParameters("U_g3_at_t0", 1e50);
    double refined_k3 = module.getVariable("k_3");
    double refined_U_g3 = module.computeU_g3(0.0);
    std::cout << "Refined parameters:\n";
    std::cout << "  k_3 = " << std::fixed << std::setprecision(2) << refined_k3 << "\n";
    std::cout << "  U_g3 at t=0 = " << std::scientific << refined_U_g3 << " J/m³ (target: 1e50)\n\n";
    
    // Step 20: Auto-refine Rossby number
    std::cout << "STEP 20: Auto-Refine Rossby Number (0.5, strong dynamo)\n";
    module.restoreState("solar_rotation_29days");
    module.autoRefineParameters("rossby_number", 0.5);
    double refined_ro = module.getVariable("rossby_number");
    std::cout << "Refined Rossby number = " << std::fixed << std::setprecision(2) 
              << refined_ro << " (target: 0.5, strong dynamo)\n\n";
    
    // Step 21: Calibrate to observations
    std::cout << "STEP 21: Calibrate to Observational Data\n";
    module.restoreState("solar_rotation_29days");
    std::map<std::string, double> observations;
    observations["omega_s"] = 3.0e-6;     // Observed faster rotation
    observations["B_j"] = 1.2e3;          // Observed field
    observations["k_3"] = 2.0;            // Observed coupling
    module.calibrateToObservations(observations);
    std::cout << "Calibrated parameters:\n";
    std::cout << "  ω_s = " << std::scientific << module.getVariable("omega_s") << " rad/s (";
    std::cout << std::fixed << std::setprecision(1) << module.computePeriod_days() << " days)\n";
    std::cout << "  B_j = " << std::scientific << module.getVariable("B_j") << " T\n";
    std::cout << "  k_3 = " << std::fixed << std::setprecision(1) << module.getVariable("k_3") << "\n\n";
    
    // Step 22: Mutate parameters
    std::cout << "STEP 22: Mutate Parameters (5% random variation)\n";
    module.mutateParameters(0.05);
    double mutated_omega = module.getVariable("omega_s");
    double mutated_period = module.computePeriod_days();
    std::cout << "Mutated parameters:\n";
    std::cout << "  ω_s = " << std::scientific << mutated_omega << " rad/s (";
    std::cout << std::fixed << std::setprecision(1) << mutated_period << " days)\n\n";
    
    // Step 23: Validate consistency
    std::cout << "STEP 23: Validate Consistency\n";
    bool valid = module.validateConsistency();
    std::cout << "Consistency check: " << (valid ? "PASSED" : "FAILED") << "\n\n";
    
    // Step 24: Introduce anomaly and auto-correct
    std::cout << "STEP 24: Introduce Anomaly and Auto-Correct\n";
    module.createVariable("omega_s_anomaly", -1e-6);
    module.removeVariable("omega_s");
    module.createVariable("omega_s", -1e-6);  // Invalid negative ω_s
    std::cout << "Introduced invalid ω_s = " << std::scientific << module.getVariable("omega_s") << " rad/s\n";
    module.autoCorrectAnomalies();
    std::cout << "Auto-corrected ω_s = " << module.getVariable("omega_s") << " rad/s (";
    std::cout << std::fixed << std::setprecision(1) << module.computePeriod_days() << " days)\n";
    bool valid_after = module.validateConsistency();
    std::cout << "Consistency after correction: " << (valid_after ? "PASSED" : "FAILED") << "\n\n";
    
    // Step 25: Test Rossby regimes
    std::cout << "STEP 25: Test High and Low Rossby Number Regimes\n";
    module.restoreState("solar_rotation_29days");
    module.optimizeForMetric("low_rossby");
    double ro_low = module.getVariable("rossby_number");
    std::cout << "Low Rossby regime (strong dynamo):\n";
    std::cout << "  Rossby number = " << std::fixed << std::setprecision(2) << ro_low << "\n";
    
    module.restoreState("solar_rotation_29days");
    module.optimizeForMetric("high_rossby");
    double ro_high = module.getVariable("rossby_number");
    std::cout << "High Rossby regime (weak dynamo):\n";
    std::cout << "  Rossby number = " << std::fixed << std::setprecision(1) << ro_high << "\n\n";
    
    // Step 26: List saved states and export
    std::cout << "STEP 26: List Saved States and Export Final State\n";
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

StellarRotationModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeOmega_s, computeOmega_s_t, computePeriod_days, computeU_g3, computeU_i) are clear, concise, and variable - driven.
- Automatic recalculation of derived variables(rho_sum) when dependencies change.
- Output and debugging functions(printVariables, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.
- Models both gravity and inertia effects of stellar / planetary rotation.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in stellar / planetary rotation modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.