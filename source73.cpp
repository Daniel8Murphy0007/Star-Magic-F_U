// NGC1300UQFFModule.h
// Modular C++ implementation of the Master Universal Gravity Equation (MUGE & UQFF Integration) for Barred Spiral Galaxy NGC 1300 Evolution.
// This module models NGC 1300's gravitational dynamics, incorporating bar-driven gas funneling, spiral arm density waves, star formation, dust lanes, and dark matter.
// Usage: #include "NGC1300UQFFModule.h" in base program; NGC1300UQFFModule mod; mod.computeG(t); mod.updateVariable("SFR", new_value);
// Variables in std::map for dynamic updates; supports F_env(t) with bar and wave terms.
// Approximations: psi_integral normalized to 1.0; H(t,z) with Omega_m=0.3, Omega_Lambda=0.7; no AGN terms.
// NGC 1300 params: M=1e11 Msun, r=11.79 kpc, SFR=1 Msun/yr, v_arm=200 km/s, B=1e-5 T, z=0.005, etc.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef NGC1300_UQFF_MODULE_H
#define NGC1300_UQFF_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>
#include <functional>
#include <vector>

class NGC1300UQFFModule {
private:
    std::map<std::string, double> variables;
    double computeHtz(double z_val);
    double computeFenv(double t);
    double computeUg1(double t);
    double computeUg2(double t);
    double computeUg3prime(double t);
    double computeUg4(double t);
    double computeUi(double t);
    double computePsiIntegral(double r, double t);
    double computeQuantumTerm(double t_Hubble_val, double r);
    double computeFluidTerm(double g_base);
    double computeDMTerm(double r);
    double computeUgSum(double r);
    double computeMsfFactor(double t);
    double computeRt(double t);

public:
    // Constructor: Initialize with NGC 1300 defaults
    NGC1300UQFFModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computation: g_NGC1300(r, t)
    double computeG(double t, double r);

    // Output descriptive text of the equation
    std::string getEquationText();

    // Print all current variables (for debugging)
    void printVariables();

    // ===== DYNAMIC SELF-UPDATE & SELF-EXPANSION CAPABILITIES =====
    
    // Dynamic variable management
    void createDynamicVariable(const std::string& name, double value);
    void removeDynamicVariable(const std::string& name);
    void cloneVariable(const std::string& source, const std::string& dest);
    void listAllVariables();
    
    // Batch operations
    void applyTransformToGroup(const std::vector<std::string>& varNames, std::function<double(double)> transform);
    void scaleVariableGroup(const std::vector<std::string>& varNames, double scale_factor);
    
    // Self-expansion capabilities
    void autoExpandParameterSpace(double scale_factor);
    void expandMassScale(double mass_multiplier);
    void expandSpatialScale(double spatial_multiplier);
    void expandTimeScale(double time_multiplier);
    
    // Self-refinement
    void autoRefineParameters(double tolerance);
    void calibrateToObservations(const std::map<std::string, double>& observed_values);
    void optimizeForMetric(const std::string& metric_name, double target_value);
    
    // Parameter exploration
    void generateVariations(int num_variations, double variation_range);
    void findOptimalParameters(const std::string& objective, int iterations);
    
    // Adaptive evolution
    void mutateParameters(double mutation_rate, double mutation_strength);
    void evolveSystem(int generations);
    
    // State management
    void saveState(const std::string& label);
    void restoreState(const std::string& label);
    void listSavedStates();
    void exportState(const std::string& filename);
    
    // System analysis
    void analyzeParameterSensitivity(const std::string& param_name);
    void generateSystemReport();
    void validatePhysicalConsistency();
    void autoCorrectAnomalies();
};

#endif // NGC1300_UQFF_MODULE_H

// NGC1300UQFFModule.cpp
#include "NGC1300UQFFModule.h"
#include <complex>

// Constructor: NGC 1300-specific values
NGC1300UQFFModule::NGC1300UQFFModule() {
    // Universal constants
    variables["G"] = 6.6743e-11;                    // m^3 kg^-1 s^-2
    variables["c"] = 3e8;                           // m/s
    variables["hbar"] = 1.0546e-34;                 // J s
    variables["Lambda"] = 1.1e-52;                  // m^-2
    variables["q"] = 1.602e-19;                     // C
    variables["pi"] = 3.141592653589793;            // pi
    variables["t_Hubble"] = 13.8e9 * 3.156e7;       // s
    variables["year_to_s"] = 3.156e7;               // s/yr
    variables["H0"] = 70.0;                         // km/s/Mpc
    variables["Mpc_to_m"] = 3.086e22;               // m/Mpc
    variables["Omega_m"] = 0.3;
    variables["Omega_Lambda"] = 0.7;
    double M_sun_val = 1.989e30;                    // kg
    double kpc_val = 3.086e19;                      // m

    // NGC 1300 parameters
    variables["M_visible"] = 7e10 * M_sun_val;      // kg
    variables["M_DM"] = 3e10 * M_sun_val;           // kg
    variables["M"] = variables["M_visible"] + variables["M_DM"];  // Total initial
    variables["M0"] = variables["M"];
    variables["SFR"] = 1 * M_sun_val / variables["year_to_s"];    // kg/s
    variables["r"] = 11.79e3 * kpc_val;             // m
    variables["z"] = 0.005;                         // Redshift
    variables["v_arm"] = 200e3;                     // m/s (gas velocity)
    variables["t"] = 1e9 * variables["year_to_s"];  // Default t=1 Gyr s

    // Dynamics
    variables["rho_fluid"] = 1e-21;                 // kg/m^3
    variables["V"] = 1e50;                          // m^3
    variables["B"] = 1e-5;                          // T
    variables["B_crit"] = 1e11;                     // T
    variables["Delta_x"] = 1e-10;                   // m
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];
    variables["integral_psi"] = 1.0;                // Normalized

    // Wave/oscillatory for spiral arms
    variables["A"] = 1e-10;
    variables["k"] = 1e20;
    variables["omega"] = 1e-15;                     // rad/s for density waves
    variables["x"] = 0.0;
    variables["v"] = variables["v_arm"];            // m/s
    variables["sigma"] = 1e3 * kpc_val;             // m for Gaussian

    // Ug subterms & Ui
    variables["Ug1"] = 0.0;                         // Dipole
    variables["Ug2"] = 0.0;                         // Superconductor
    variables["Ug3"] = 0.0;                         // External
    variables["Ug4"] = 0.0;                         // Reaction
    variables["Ui"] = 0.0;
    variables["mu_0"] = 4 * variables["pi"] * 1e-7; // H/m
    variables["rho_vac_SCm"] = 7.09e-37;            // J/m^3
    variables["rho_vac_UA"] = 7.09e-36;             // J/m^3
    variables["lambda_I"] = 1.0;
    variables["omega_i"] = 1e-8;                    // rad/s
    variables["t_n"] = 0.0;
    variables["F_RZ"] = 0.01;
    variables["k_4"] = 1.0;
    variables["k_SF"] = 1e-10;                      // N/Msun, adjusted to m/s^2
    variables["omega_spin"] = 1e-4;                 // rad/s (bar rotation)
    variables["I_dipole"] = 1e20;                   // A
    variables["A_dipole"] = 1e15;                   // m^2
    variables["H_aether"] = 1e-6;                   // A/m
    variables["delta_rho_over_rho"] = 1e-5;

    // Scales
    variables["scale_macro"] = 1e-12;
    variables["f_TRZ"] = 0.1;
    variables["f_sc"] = 1.0;
    variables["v_r"] = 1e3;                         // m/s radial velocity
    variables["rho"] = variables["rho_fluid"];
}

// Update variable (with dependents)
void NGC1300UQFFModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding." << std::endl;
        variables[name] = value;
    }
    if (name == "Delta_x") {
        variables["Delta_p"] = variables["hbar"] / value;
    } else if (name == "M") {
        variables["M_visible"] = 0.7 * value;
        variables["M_DM"] = 0.3 * value;
        variables["M0"] = value;
    } else if (name == "SFR") {
        // Adjust units if needed
    }
}

// Add/subtract
void NGC1300UQFFModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        variables[name] = delta;
    }
}
void NGC1300UQFFModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute H(t, z)
double NGC1300UQFFModule::computeHtz(double z_val) {
    double Hz_kms = variables["H0"] * std::sqrt(variables["Omega_m"] * std::pow(1.0 + z_val, 3) + variables["Omega_Lambda"]);
    return (Hz_kms * 1e3) / variables["Mpc_to_m"];
}

// M(t)
double NGC1300UQFFModule::computeMsfFactor(double t) {
    return variables["SFR"] * t / variables["M0"];
}

// r(t)
double NGC1300UQFFModule::computeRt(double t) {
    return variables["r"] + variables["v_r"] * t;
}

// F_env(t)
double NGC1300UQFFModule::computeFenv(double t) {
    double F_bar = 0.1 * (variables["G"] * variables["M"]) / (variables["r"] * variables["r"]);  // Bar funneling
    double F_SF = variables["k_SF"] * variables["SFR"] / 1.989e30;  // Normalize to m/s^2
    double F_wave = variables["rho_fluid"] * std::pow(variables["v_arm"], 2);  // Density wave
    return F_bar + F_SF + F_wave;
}

// Ug1: dipole
double NGC1300UQFFModule::computeUg1(double t) {
    double mu_dipole = variables["I_dipole"] * variables["A_dipole"] * variables["omega_spin"];
    return mu_dipole * variables["B"];
}

// Ug2: superconductor
double NGC1300UQFFModule::computeUg2(double t) {
    double B_super = variables["mu_0"] * variables["H_aether"];
    return (B_super * B_super) / (2 * variables["mu_0"]);
}

// Ug3': external (bar as external)
double NGC1300UQFFModule::computeUg3prime(double t) {
    double M_bar = 0.2 * variables["M"];  // Bar mass fraction
    double r_bar = 0.3 * variables["r"];  // Bar radius
    return (variables["G"] * M_bar) / (r_bar * r_bar);
}

// Ug4: reaction
double NGC1300UQFFModule::computeUg4(double t) {
    double E_react = 1e46 * std::exp(-0.0005 * t);
    return variables["k_4"] * E_react;
}

// Ui
double NGC1300UQFFModule::computeUi(double t) {
    return variables["lambda_I"] * (variables["rho_vac_SCm"] / variables["rho_vac_UA"]) * variables["omega_i"] * std::cos(variables["pi"] * variables["t_n"]) * (1 + variables["F_RZ"]);
}

// Psi integral (simplified)
double NGC1300UQFFModule::computePsiIntegral(double r, double t) {
    double A = variables["A"];
    double m = 2.0;  // m-mode for spiral
    double omega = variables["omega"];
    double sigma = variables["sigma"];
    std::complex<double> psi_spiral(A * std::exp(-r*r / (2 * sigma * sigma)) * std::exp(std::complex<double>(0, m * 0 - omega * t)));
    return std::norm(psi_spiral);  // |psi|^2
}

// Quantum term
double NGC1300UQFFModule::computeQuantumTerm(double t_Hubble_val, double r) {
    double unc = std::sqrt(variables["Delta_x"] * variables["Delta_p"]);
    double psi_int = computePsiIntegral(r, variables["t"]);
    return (variables["hbar"] / unc) * variables["integral_psi"] * (2 * variables["pi"] / t_Hubble_val) * psi_int;
}

// Fluid
double NGC1300UQFFModule::computeFluidTerm(double g_base) {
    return variables["rho_fluid"] * variables["V"] * g_base;
}

// DM
double NGC1300UQFFModule::computeDMTerm(double r) {
    double pert = variables["delta_rho_over_rho"];
    double curv = 3 * variables["G"] * variables["M"] / (r * r * r);
    return (variables["M_visible"] + variables["M_DM"]) * (pert + curv);
}

// Ug sum
double NGC1300UQFFModule::computeUgSum(double r) {
    double Ug_base = (variables["G"] * variables["M"]) / (r * r);
    variables["Ug1"] = computeUg1(variables["t"]);
    variables["Ug2"] = computeUg2(variables["t"]);
    variables["Ug3"] = computeUg3prime(variables["t"]);
    variables["Ug4"] = computeUg4(variables["t"]);
    return Ug_base + variables["Ug1"] + variables["Ug2"] + variables["Ug3"] + variables["Ug4"];
}

// Full g_NGC1300
double NGC1300UQFFModule::computeG(double t, double r) {
    variables["t"] = t;
    double msf_factor = computeMsfFactor(t);
    double m_factor = 1.0 + msf_factor;
    double Hz = computeHtz(variables["z"]);
    double expansion = 1.0 + Hz * t;
    double sc_correction = 1.0 - (variables["B"] / variables["B_crit"]);
    double f_env = computeFenv(t);
    double tr_factor = 1.0 + variables["f_TRZ"];
    double rt = computeRt(t);  // But use input r for profile

    // Base gravity
    double g_base = (variables["G"] * variables["M"] * m_factor / (r * r)) * expansion * sc_correction * (1.0 + f_env) * tr_factor;

    // Ug sum (includes base? Adjust: Ug sum without base)
    double ug_sum = computeUgSum(r) - g_base;  // Subtract to avoid double-count

    // Cosmological
    double lambda_term = variables["Lambda"] * (variables["c"] * variables["c"]) / 3.0;

    // Ui
    double ui_term = computeUi(t);

    // Quantum
    double quantum_term = computeQuantumTerm(variables["t_Hubble"], r);

    // Fluid
    double fluid_term = computeFluidTerm(g_base);

    // DM
    double dm_term = computeDMTerm(r);

    // Total
    return g_base + ug_sum + lambda_term + ui_term + quantum_term + fluid_term + dm_term;
}

// Equation text
std::string NGC1300UQFFModule::getEquationText() {
    return "g_NGC1300(r, t) = (G * M(t) / r(t)^2) * (1 + H(t, z)) * (1 - B(t) / B_crit) * (1 + F_env(t)) + "
           "(U_g1 + U_g2 + U_g3' + U_g4) + U_i + (Lambda * c^2 / 3) + "
           "(hbar / sqrt(?x * ?p)) * ? (?_total * H * ?_total dV) * (2? / t_Hubble) + "
           "?_fluid * V * g + (M_visible + M_DM) * (??/? + 3 G M / r^3)\n"
           "Where: M(t) = M * (1 + M_SF(t)); M_SF(t) = SFR * t; r(t) = r0 + v_r t;\n"
           "H(t, z) = H0 * sqrt(?m (1+z)^3 + ??); F_env(t) = F_bar + F_SF + F_wave;\n"
           "F_bar = 0.1 G M / r^2; F_wave = ? v_arm^2; U_g1 = ?_dipole * B; U_g2 = B_super^2 / (2 ?0);\n"
           "U_g3' = G M_bar / r_bar^2; U_g4 = k4 * E_react(t); U_i = ?_I * (?_SCm/?_UA) * ?_i * cos(? t_n) * (1 + F_RZ);\n"
           "?_total = A exp(-r^2/(2?^2)) exp(i(m? - ? t)) + bar terms; Insights: Attractive (g_base, Ug1, Ug3') vs. Repulsive (U_g2, ?) advance UQFF.\n"
           "Adaptations: Hubble ACS 2004 data; SFR=1 Msun/yr; M=1e11 Msun. Solutions: g ~2e36 m/s� at t=1 Gyr (DM/fluid dominant).";
}

// Print
void NGC1300UQFFModule::printVariables() {
    std::cout << "NGC 1300 Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// ===== DYNAMIC SELF-UPDATE & SELF-EXPANSION IMPLEMENTATIONS =====

// Static storage for saved states
static std::map<std::string, std::map<std::string, double>> ngc1300_saved_states;

// 1. Dynamic variable management
void NGC1300UQFFModule::createDynamicVariable(const std::string& name, double value) {
    variables[name] = value;
    std::cout << "Created dynamic variable: " << name << " = " << value << std::endl;
}

void NGC1300UQFFModule::removeDynamicVariable(const std::string& name) {
    if (variables.find(name) != variables.end()) {
        variables.erase(name);
        std::cout << "Removed dynamic variable: " << name << std::endl;
    } else {
        std::cerr << "Variable '" << name << "' not found for removal." << std::endl;
    }
}

void NGC1300UQFFModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
        std::cout << "Cloned " << source << " to " << dest << std::endl;
    } else {
        std::cerr << "Source variable '" << source << "' not found." << std::endl;
    }
}

void NGC1300UQFFModule::listAllVariables() {
    std::cout << "=== All NGC 1300 Variables (Total: " << variables.size() << ") ===" << std::endl;
    for (const auto& pair : variables) {
        std::cout << "  " << pair.first << " = " << pair.second << std::endl;
    }
}

// 2. Batch operations
void NGC1300UQFFModule::applyTransformToGroup(const std::vector<std::string>& varNames,
                                                std::function<double(double)> transform) {
    for (const auto& name : varNames) {
        if (variables.find(name) != variables.end()) {
            variables[name] = transform(variables[name]);
            std::cout << "Transformed " << name << " to " << variables[name] << std::endl;
        }
    }
}

void NGC1300UQFFModule::scaleVariableGroup(const std::vector<std::string>& varNames, double scale_factor) {
    applyTransformToGroup(varNames, [scale_factor](double val) { return val * scale_factor; });
}

// 3. Self-expansion capabilities
void NGC1300UQFFModule::autoExpandParameterSpace(double scale_factor) {
    std::cout << "Auto-expanding NGC 1300 parameter space by factor " << scale_factor << std::endl;
    std::vector<std::string> expandable = {"M", "SFR", "r", "v_arm", "rho_fluid"};
    scaleVariableGroup(expandable, scale_factor);
    
    // Auto-sync dependent variables
    variables["M_visible"] = 0.7 * variables["M"];
    variables["M_DM"] = 0.3 * variables["M"];
    variables["M0"] = variables["M"];
    
    std::cout << "  Parameter space expanded" << std::endl;
}

void NGC1300UQFFModule::expandMassScale(double mass_multiplier) {
    std::cout << "Expanding mass scale by " << mass_multiplier << std::endl;
    variables["M"] *= mass_multiplier;
    variables["M_visible"] = 0.7 * variables["M"];
    variables["M_DM"] = 0.3 * variables["M"];
    variables["M0"] = variables["M"];
    variables["SFR"] *= mass_multiplier;
    std::cout << "  M: " << variables["M"] << " kg" << std::endl;
    std::cout << "  M_visible: " << variables["M_visible"] << " kg" << std::endl;
    std::cout << "  M_DM: " << variables["M_DM"] << " kg" << std::endl;
    std::cout << "  SFR: " << variables["SFR"] << " kg/s" << std::endl;
}

void NGC1300UQFFModule::expandSpatialScale(double spatial_multiplier) {
    std::cout << "Expanding spatial scale by " << spatial_multiplier << std::endl;
    variables["r"] *= spatial_multiplier;
    variables["Delta_x"] *= spatial_multiplier;
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];
    variables["sigma"] *= spatial_multiplier;
    std::cout << "  r: " << variables["r"] << " m" << std::endl;
    std::cout << "  Delta_x: " << variables["Delta_x"] << " m" << std::endl;
    std::cout << "  sigma: " << variables["sigma"] << " m" << std::endl;
}

void NGC1300UQFFModule::expandTimeScale(double time_multiplier) {
    std::cout << "Expanding time scale by " << time_multiplier << std::endl;
    variables["t"] *= time_multiplier;
    variables["t_Hubble"] *= time_multiplier;
    variables["omega"] /= time_multiplier;
    variables["omega_i"] /= time_multiplier;
    variables["omega_spin"] /= time_multiplier;
    std::cout << "  t: " << variables["t"] << " s" << std::endl;
    std::cout << "  t_Hubble: " << variables["t_Hubble"] << " s" << std::endl;
    std::cout << "  omega: " << variables["omega"] << " rad/s" << std::endl;
}

// 4. Self-refinement
void NGC1300UQFFModule::autoRefineParameters(double tolerance) {
    std::cout << "Auto-refining NGC 1300 parameters with tolerance " << tolerance << std::endl;
    
    // Ensure M = M_visible + M_DM
    double M_total_computed = variables["M_visible"] + variables["M_DM"];
    if (std::abs(M_total_computed - variables["M"]) / variables["M"] > tolerance) {
        std::cout << "  Correcting total mass consistency" << std::endl;
        variables["M"] = M_total_computed;
        variables["M0"] = variables["M"];
    }
    
    // Ensure Delta_p = hbar / Delta_x
    double delta_p_computed = variables["hbar"] / variables["Delta_x"];
    if (std::abs(delta_p_computed - variables["Delta_p"]) / variables["Delta_p"] > tolerance) {
        std::cout << "  Correcting uncertainty relation" << std::endl;
        variables["Delta_p"] = delta_p_computed;
    }
    
    // Ensure v = v_arm consistency
    if (std::abs(variables["v"] - variables["v_arm"]) / variables["v_arm"] > tolerance) {
        std::cout << "  Synchronizing velocity parameters" << std::endl;
        variables["v"] = variables["v_arm"];
    }
    
    std::cout << "Refinement complete." << std::endl;
}

void NGC1300UQFFModule::calibrateToObservations(const std::map<std::string, double>& observed_values) {
    std::cout << "Calibrating to " << observed_values.size() << " NGC 1300 observations..." << std::endl;
    for (const auto& obs : observed_values) {
        if (variables.find(obs.first) != variables.end()) {
            double old_val = variables[obs.first];
            updateVariable(obs.first, obs.second);
            std::cout << "  " << obs.first << ": " << old_val << " -> " << obs.second << std::endl;
        }
    }
    std::cout << "Calibration complete." << std::endl;
}

void NGC1300UQFFModule::optimizeForMetric(const std::string& metric_name, double target_value) {
    std::cout << "Optimizing for metric: " << metric_name << " = " << target_value << std::endl;
    
    if (metric_name == "g_total") {
        double current_g = computeG(variables["t"], variables["r"]);
        double ratio = target_value / std::max(current_g, 1e-100);
        
        // Adjust M to reach target gravity
        variables["M"] *= ratio;
        variables["M_visible"] = 0.7 * variables["M"];
        variables["M_DM"] = 0.3 * variables["M"];
        variables["M0"] = variables["M"];
        std::cout << "  Adjusted M by " << ratio << std::endl;
    } else if (metric_name == "SFR") {
        updateVariable("SFR", target_value);
        std::cout << "  Set SFR to " << target_value << std::endl;
    }
    
    std::cout << "Optimization complete." << std::endl;
}

// 5. Parameter exploration
void NGC1300UQFFModule::generateVariations(int num_variations, double variation_range) {
    std::cout << "Generating " << num_variations << " NGC 1300 variations with range ±" 
              << (variation_range * 100) << "%" << std::endl;
    
    std::vector<std::string> key_params = {"M", "SFR", "r", "v_arm", "B", "rho_fluid"};
    
    for (int i = 0; i < num_variations; ++i) {
        std::cout << "  Variation " << (i+1) << ":" << std::endl;
        for (const auto& param : key_params) {
            if (variables.find(param) != variables.end()) {
                double base = variables[param];
                double variation = base * (1.0 + variation_range * (2.0 * (rand() / (double)RAND_MAX) - 1.0));
                std::cout << "    " << param << ": " << base << " -> " << variation << std::endl;
            }
        }
    }
}

void NGC1300UQFFModule::findOptimalParameters(const std::string& objective, int iterations) {
    std::cout << "Finding optimal NGC 1300 parameters for: " << objective 
              << " (" << iterations << " iterations)" << std::endl;
    
    double best_score = -1e100;
    std::map<std::string, double> best_params;
    
    for (int i = 0; i < iterations; ++i) {
        mutateParameters(0.7, 0.1);
        
        double score = computeG(variables["t"], variables["r"]);
        
        if (objective == "maximize_gravity") {
            if (score > best_score) {
                best_score = score;
                best_params = variables;
            }
        } else if (objective == "target_2e36") {
            if (std::abs(score - 2e36) < std::abs(best_score - 2e36)) {
                best_score = score;
                best_params = variables;
            }
        }
    }
    
    variables = best_params;
    std::cout << "Optimal g_NGC1300: " << best_score << " m/s^2" << std::endl;
}

// 6. Adaptive evolution
void NGC1300UQFFModule::mutateParameters(double mutation_rate, double mutation_strength) {
    std::vector<std::string> mutable_params = {"M", "SFR", "r", "v_arm", "B", "rho_fluid", "omega_spin"};
    
    for (const auto& param : mutable_params) {
        if (variables.find(param) != variables.end()) {
            if ((rand() / (double)RAND_MAX) < mutation_rate) {
                double mutation = 1.0 + mutation_strength * (2.0 * (rand() / (double)RAND_MAX) - 1.0);
                variables[param] *= mutation;
            }
        }
    }
    
    // Auto-sync dependencies
    variables["M_visible"] = 0.7 * variables["M"];
    variables["M_DM"] = 0.3 * variables["M"];
    variables["M0"] = variables["M"];
    variables["v"] = variables["v_arm"];
}

void NGC1300UQFFModule::evolveSystem(int generations) {
    std::cout << "Evolving NGC 1300 system over " << generations << " generations..." << std::endl;
    
    for (int gen = 0; gen < generations; ++gen) {
        mutateParameters(0.3, 0.08);
        
        double fitness = computeG(variables["t"], variables["r"]);
        
        if (gen % 10 == 0) {
            std::cout << "  Gen " << gen << ": g_NGC1300 = " << fitness << " m/s^2" << std::endl;
        }
    }
    
    std::cout << "Evolution complete." << std::endl;
}

// 7. State management
void NGC1300UQFFModule::saveState(const std::string& label) {
    ngc1300_saved_states[label] = variables;
    std::cout << "Saved NGC 1300 state: " << label << " (" << variables.size() << " variables)" << std::endl;
}

void NGC1300UQFFModule::restoreState(const std::string& label) {
    if (ngc1300_saved_states.find(label) != ngc1300_saved_states.end()) {
        variables = ngc1300_saved_states[label];
        std::cout << "Restored NGC 1300 state: " << label << std::endl;
    } else {
        std::cerr << "State '" << label << "' not found." << std::endl;
    }
}

void NGC1300UQFFModule::listSavedStates() {
    std::cout << "=== Saved NGC 1300 States (Total: " << ngc1300_saved_states.size() << ") ===" << std::endl;
    for (const auto& state : ngc1300_saved_states) {
        std::cout << "  " << state.first << " (" << state.second.size() << " variables)" << std::endl;
    }
}

void NGC1300UQFFModule::exportState(const std::string& filename) {
    std::cout << "Exporting NGC 1300 state to " << filename << " (not implemented - placeholder)" << std::endl;
    // In real implementation: write variables to file
}

// 8. System analysis
void NGC1300UQFFModule::analyzeParameterSensitivity(const std::string& param_name) {
    if (variables.find(param_name) == variables.end()) {
        std::cerr << "Parameter '" << param_name << "' not found." << std::endl;
        return;
    }
    
    std::cout << "=== NGC 1300 Sensitivity Analysis: " << param_name << " ===" << std::endl;
    
    double base_value = variables[param_name];
    double base_output = computeG(variables["t"], variables["r"]);
    
    std::vector<double> perturbations = {0.7, 0.85, 1.0, 1.15, 1.3};
    
    for (double factor : perturbations) {
        updateVariable(param_name, base_value * factor);
        
        double new_output = computeG(variables["t"], variables["r"]);
        double sensitivity = (new_output - base_output) / std::max(std::abs(base_output), 1e-100);
        
        std::cout << "  " << param_name << " * " << factor << " -> g_NGC1300 change: " 
                  << (sensitivity * 100) << "%" << std::endl;
    }
    
    updateVariable(param_name, base_value);  // Restore
}

void NGC1300UQFFModule::generateSystemReport() {
    std::cout << "\n========== NGC 1300 Barred Spiral Galaxy System Report ==========" << std::endl;
    std::cout << "Total Variables: " << variables.size() << std::endl;
    
    // Galaxy parameters
    std::cout << "\nGalaxy Parameters:" << std::endl;
    std::cout << "M (total): " << variables["M"] << " kg (" << (variables["M"]/1.989e30) << " M_sun)" << std::endl;
    std::cout << "M_visible: " << variables["M_visible"] << " kg" << std::endl;
    std::cout << "M_DM: " << variables["M_DM"] << " kg" << std::endl;
    std::cout << "r (radius): " << variables["r"] << " m (" << (variables["r"]/3.086e19) << " kpc)" << std::endl;
    std::cout << "z (redshift): " << variables["z"] << std::endl;
    
    std::cout << "\nDynamics:" << std::endl;
    std::cout << "SFR: " << variables["SFR"] << " kg/s (" << (variables["SFR"]*3.156e7/1.989e30) << " M_sun/yr)" << std::endl;
    std::cout << "v_arm: " << variables["v_arm"] << " m/s (" << (variables["v_arm"]/1e3) << " km/s)" << std::endl;
    std::cout << "omega_spin (bar): " << variables["omega_spin"] << " rad/s" << std::endl;
    std::cout << "B (magnetic): " << variables["B"] << " T" << std::endl;
    
    std::cout << "\nEnvironmental Forces:" << std::endl;
    double f_env = computeFenv(variables["t"]);
    std::cout << "F_env (total): " << f_env << " N" << std::endl;
    
    std::cout << "\nGravity Components:" << std::endl;
    double Ug1 = computeUg1(variables["t"]);
    double Ug2 = computeUg2(variables["t"]);
    double Ug3 = computeUg3prime(variables["t"]);
    double Ug4 = computeUg4(variables["t"]);
    double Ui = computeUi(variables["t"]);
    std::cout << "Ug1 (dipole): " << Ug1 << " J" << std::endl;
    std::cout << "Ug2 (superconductor): " << Ug2 << " J" << std::endl;
    std::cout << "Ug3' (bar external): " << Ug3 << " m/s^2" << std::endl;
    std::cout << "Ug4 (reaction): " << Ug4 << " J" << std::endl;
    std::cout << "Ui (internal): " << Ui << std::endl;
    
    // Current computation
    double g_total = computeG(variables["t"], variables["r"]);
    
    std::cout << "\nCurrent Computation:" << std::endl;
    std::cout << "t: " << variables["t"] << " s (" << (variables["t"]/(3.156e7*1e9)) << " Gyr)" << std::endl;
    std::cout << "g_NGC1300: " << g_total << " m/s^2" << std::endl;
    
    std::cout << "\nPhysics Regime:" << std::endl;
    if (variables["B"] / variables["B_crit"] < 0.01) {
        std::cout << "Weak magnetic regime (B << B_crit)" << std::endl;
    } else if (variables["B"] / variables["B_crit"] < 0.5) {
        std::cout << "Moderate magnetic regime" << std::endl;
    } else {
        std::cout << "Strong magnetic regime (approaching B_crit)" << std::endl;
    }
    
    std::cout << "==================================================\n" << std::endl;
}

void NGC1300UQFFModule::validatePhysicalConsistency() {
    std::cout << "Validating NGC 1300 physical consistency..." << std::endl;
    bool consistent = true;
    
    // Check for NaN/Inf
    for (const auto& pair : variables) {
        if (std::isnan(pair.second) || std::isinf(pair.second)) {
            std::cerr << "  ERROR: " << pair.first << " is NaN/Inf" << std::endl;
            consistent = false;
        }
    }
    
    // Check mass consistency
    double M_total = variables["M_visible"] + variables["M_DM"];
    if (std::abs(M_total - variables["M"]) / variables["M"] > 0.01) {
        std::cerr << "  WARNING: M != M_visible + M_DM (discrepancy: " 
                  << ((M_total - variables["M"]) / variables["M"] * 100) << "%)" << std::endl;
    }
    
    // Check uncertainty relation
    double delta_p_expected = variables["hbar"] / variables["Delta_x"];
    if (std::abs(delta_p_expected - variables["Delta_p"]) / variables["Delta_p"] > 0.01) {
        std::cerr << "  WARNING: Delta_p != hbar / Delta_x" << std::endl;
    }
    
    // Positive values
    if (variables["M"] <= 0) {
        std::cerr << "  ERROR: M must be positive" << std::endl;
        consistent = false;
    }
    
    if (variables["r"] <= 0) {
        std::cerr << "  ERROR: r must be positive" << std::endl;
        consistent = false;
    }
    
    if (variables["SFR"] < 0) {
        std::cerr << "  ERROR: SFR cannot be negative" << std::endl;
        consistent = false;
    }
    
    // Physical ranges
    if (variables["z"] < 0) {
        std::cerr << "  WARNING: Negative redshift (blueshift)" << std::endl;
    }
    
    if (variables["B"] > variables["B_crit"]) {
        std::cerr << "  WARNING: B exceeds B_crit (superconductor breakdown)" << std::endl;
    }
    
    if (consistent) {
        std::cout << "  All checks passed. NGC 1300 system is physically consistent." << std::endl;
    }
}

void NGC1300UQFFModule::autoCorrectAnomalies() {
    std::cout << "Auto-correcting NGC 1300 anomalies..." << std::endl;
    
    // Fix NaN/Inf
    for (auto& pair : variables) {
        if (std::isnan(pair.second) || std::isinf(pair.second)) {
            std::cout << "  Correcting " << pair.first << " (was NaN/Inf)" << std::endl;
            pair.second = 1.0;
        }
    }
    
    // Fix mass consistency
    double M_total = variables["M_visible"] + variables["M_DM"];
    if (std::abs(M_total - variables["M"]) / variables["M"] > 0.01) {
        std::cout << "  Correcting total mass to match components" << std::endl;
        variables["M"] = M_total;
        variables["M0"] = M_total;
    }
    
    // Fix uncertainty relation
    double delta_p_expected = variables["hbar"] / variables["Delta_x"];
    if (std::abs(delta_p_expected - variables["Delta_p"]) / variables["Delta_p"] > 0.01) {
        std::cout << "  Correcting Delta_p = hbar / Delta_x" << std::endl;
        variables["Delta_p"] = delta_p_expected;
    }
    
    // Ensure positive values
    if (variables["M"] <= 0) {
        std::cout << "  Correcting M to 1e11 M_sun" << std::endl;
        variables["M"] = 1e11 * 1.989e30;
        variables["M_visible"] = 0.7 * variables["M"];
        variables["M_DM"] = 0.3 * variables["M"];
        variables["M0"] = variables["M"];
    }
    
    if (variables["r"] <= 0) {
        std::cout << "  Correcting r to 11.79 kpc" << std::endl;
        variables["r"] = 11.79e3 * 3.086e19;
    }
    
    if (variables["SFR"] < 0) {
        std::cout << "  Correcting SFR to 1 M_sun/yr" << std::endl;
        variables["SFR"] = 1.989e30 / 3.156e7;
    }
    
    // Sync velocity
    if (variables["v"] != variables["v_arm"]) {
        std::cout << "  Synchronizing v = v_arm" << std::endl;
        variables["v"] = variables["v_arm"];
    }
    
    std::cout << "Auto-correction complete." << std::endl;
}

// Example usage
// Uncomment the following code to test the enhanced NGC 1300 module with dynamic capabilities
/*
#include "NGC1300UQFFModule.h"
int main() {
    NGC1300UQFFModule mod;
    
    std::cout << "===== NGC 1300 Barred Spiral Galaxy with Dynamic Capabilities =====" << std::endl;
    std::cout << "  UQFF Integration with bar-driven gas funneling and spiral arm dynamics\n" << std::endl;
    
    // Initial state
    std::cout << "1. Initial gravity calculation:" << std::endl;
    double t = 1e9 * 3.156e7;  // 1 Gyr
    double r = 5e3 * 3.086e19;  // 5 kpc
    double g = mod.computeG(t, r);
    std::cout << "g_NGC1300 = " << g << " m/s^2\n" << std::endl;
    
    // List all variables
    std::cout << "2. Variable inventory:" << std::endl;
    mod.listAllVariables();
    std::cout << std::endl;
    
    // Generate system report
    std::cout << "3. System report:" << std::endl;
    mod.generateSystemReport();
    std::cout << std::endl;
    
    // Save initial state
    std::cout << "4. Saving initial state:" << std::endl;
    mod.saveState("initial");
    std::cout << std::endl;
    
    // Double star formation rate
    std::cout << "5. Doubling SFR:" << std::endl;
    mod.updateVariable("SFR", 2 * mod.variables["SFR"]);
    g = mod.computeG(t, r);
    std::cout << "New g_NGC1300 = " << g << " m/s^2\n" << std::endl;
    
    // Expand mass scale
    std::cout << "6. Expanding mass scale by 1.5x:" << std::endl;
    mod.expandMassScale(1.5);
    g = mod.computeG(t, r);
    std::cout << "New g_NGC1300 = " << g << " m/s^2\n" << std::endl;
    
    // Save evolved state
    std::cout << "7. Saving evolved state:" << std::endl;
    mod.saveState("evolved_1.5x");
    std::cout << std::endl;
    
    // Analyze sensitivity to M
    std::cout << "8. Sensitivity analysis for M:" << std::endl;
    mod.analyzeParameterSensitivity("M");
    std::cout << std::endl;
    
    // Create dynamic variable
    std::cout << "9. Creating dynamic bar strength variable:" << std::endl;
    mod.createDynamicVariable("bar_strength", 0.25);
    std::cout << std::endl;
    
    // Optimize for target gravity
    std::cout << "10. Optimizing for target g = 2e36 m/s^2:" << std::endl;
    mod.optimizeForMetric("g_total", 2e36);
    g = mod.computeG(t, r);
    std::cout << "Optimized g_NGC1300 = " << g << " m/s^2\n" << std::endl;
    
    // Validate consistency
    std::cout << "11. Validating physical consistency:" << std::endl;
    mod.validatePhysicalConsistency();
    std::cout << std::endl;
    
    // Generate variations
    std::cout << "12. Generating 3 parameter variations (±15%):" << std::endl;
    mod.generateVariations(3, 0.15);
    std::cout << std::endl;
    
    // Restore initial state
    std::cout << "13. Restoring initial state:" << std::endl;
    mod.restoreState("initial");
    g = mod.computeG(t, r);
    std::cout << "Restored g_NGC1300 = " << g << " m/s^2\n" << std::endl;
    
    // List saved states
    std::cout << "14. List of saved states:" << std::endl;
    mod.listSavedStates();
    std::cout << std::endl;
    
    std::cout << "End of NGC 1300 demonstration with dynamic capabilities.\n" << std::endl;
    std::cout << mod.getEquationText() << std::endl;
    return 0;
}
*/
// Compile: g++ -o ngc1300_sim base.cpp NGC1300UQFFModule.cpp -lm
// Sample Output: g_NGC1300 ~ 2e36 m/s^2 (env/fluid dominant; repulsive terms advance framework).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

NGC1300UQFFModule Evaluation

Strengths :
-Modular, extensible design for modeling NGC 1300 galaxy gravity, including bar - driven gas funneling, spiral arm density waves, star formation, dust lanes, and dark matter.
- Comprehensive physics : gravity, cosmological expansion, magnetic fields, environmental / bar / wave effects, quantum, fluid, and DM terms.
- Dynamic variable management via std::map enables runtime updates and system adaptation.
- Clear separation of computation functions(e.g., Ug1�Ug4, F_env, quantum, fluid, DM), aiding maintainability.
- NGC 1300 - specific parameters are initialized for realistic simulation; supports easy modification.
- Output functions for equation text and variable state support debugging and documentation.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for flexibility and scalability.
- Some calculations use magic numbers or lack explanatory comments; define named constants and clarify logic.
- Minimal error handling(e.g., division by zero, invalid variable names); add validation for robustness.
- Unit consistency should be checked and documented for all physical quantities.
- For large - scale or performance - critical simulations, optimize data structures and reduce redundant calculations.
- std::map is flexible but may be less efficient than structured types for very large models.
- Expand documentation for function purposes and physical meaning.

Summary:
The code is well - structured, flexible, and suitable for scientific prototyping and educational use in galactic dynamics modeling.It implements a broad set of physical effects and adapts to various scenarios.For production or high - performance applications, address the recommendations for improved robustness, maintainability, and scalability.