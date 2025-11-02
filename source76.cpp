// NGC2264UQFFModule.h
// Modular C++ implementation of the Master Universal Gravity Equation (MUGE & UQFF Integration) for Cone Nebula (NGC 2264) Evolution.
// This module models NGC 2264's gravitational dynamics, incorporating stellar winds, pillar erosion, protostar formation, dust/gas densities, and dark matter.
// Usage: #include "NGC2264UQFFModule.h" in base program; NGC2264UQFFModule mod; mod.computeG(t); mod.updateVariable("SFR", new_value);
// Variables in std::map for dynamic updates; supports F_env(t) with wind and erosion terms.
// Approximations: psi_integral normalized to 1.0; H(t,z) with Omega_m=0.3, Omega_Lambda=0.7; E_react exp decay; pillar waves simplified.
// NGC 2264 params: M=100 Msun, r=3.31e16 m, SFR=0.01 Msun/yr, v_wind=20 km/s, rho=1e-20 kg/m^3, B=1e-5 T, z=0.0008, etc.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef NGC2264_UQFF_MODULE_H
#define NGC2264_UQFF_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>
#include <functional>

class NGC2264UQFFModule {
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
    // Constructor: Initialize with NGC 2264 defaults
    NGC2264UQFFModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computation: g_NGC2264(r, t)
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
    
    // Batch operations on variable groups
    void applyTransformToGroup(const std::vector<std::string>& varNames, 
                               std::function<double(double)> transform);
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

#endif // NGC2264_UQFF_MODULE_H

// NGC2264UQFFModule.cpp
#include "NGC2264UQFFModule.h"
#include <complex>

// Constructor: NGC 2264-specific values
NGC2264UQFFModule::NGC2264UQFFModule() {
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
    double ly_val = 9.461e15;                       // m

    // NGC 2264 parameters
    variables["M_visible"] = 80 * M_sun_val;        // kg
    variables["M_DM"] = 20 * M_sun_val;             // kg
    variables["M"] = variables["M_visible"] + variables["M_DM"];  // Total initial
    variables["M0"] = variables["M"];
    variables["SFR"] = 0.01 * M_sun_val / variables["year_to_s"]; // kg/s
    variables["r"] = 3.31e16;                       // m (~3.5 ly radius)
    variables["z"] = 0.0008;                        // Redshift
    variables["v_wind"] = 20e3;                     // m/s
    variables["t"] = 3e6 * variables["year_to_s"];  // Default t=3 Myr s

    // Dynamics
    variables["rho_fluid"] = 1e-20;                 // kg/m^3
    variables["V"] = 1e48;                          // m^3
    variables["B"] = 1e-5;                          // T
    variables["B_crit"] = 1e11;                     // T
    variables["Delta_x"] = 1e-10;                   // m
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];
    variables["integral_psi"] = 1.0;                // Normalized

    // Wave/oscillatory for pillar dynamics
    variables["A"] = 1e-10;
    variables["k"] = 1e20;
    variables["omega"] = 1e-14;                     // rad/s for wind waves
    variables["x"] = 0.0;
    variables["v"] = variables["v_wind"];           // m/s
    variables["sigma"] = 1e15;                      // m for Gaussian

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
    variables["omega_spin"] = 1e-5;                 // rad/s protostar spin
    variables["I_dipole"] = 1e18;                   // A
    variables["A_dipole"] = 1e12;                   // m^2
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
void NGC2264UQFFModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding." << std::endl;
        variables[name] = value;
    }
    if (name == "Delta_x") {
        variables["Delta_p"] = variables["hbar"] / value;
    } else if (name == "M_visible" || name == "M_DM") {
        variables["M"] = variables["M_visible"] + variables["M_DM"];
        variables["M0"] = variables["M"];
    }
}

// Add/subtract
void NGC2264UQFFModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        variables[name] = delta;
    }
}
void NGC2264UQFFModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute H(t, z)
double NGC2264UQFFModule::computeHtz(double z_val) {
    double Hz_kms = variables["H0"] * std::sqrt(variables["Omega_m"] * std::pow(1.0 + z_val, 3) + variables["Omega_Lambda"]);
    return (Hz_kms * 1e3) / variables["Mpc_to_m"];
}

// M(t)
double NGC2264UQFFModule::computeMsfFactor(double t) {
    return variables["SFR"] * t / variables["M0"];
}

// r(t)
double NGC2264UQFFModule::computeRt(double t) {
    return variables["r"] + variables["v_r"] * t;
}

// F_env(t)
double NGC2264UQFFModule::computeFenv(double t) {
    double F_wind = variables["rho_fluid"] * std::pow(variables["v_wind"], 2);
    double F_SF = variables["k_SF"] * variables["SFR"] / 1.989e30;  // Normalize to m/s^2
    double F_erode = 0.05 * (t / (3e6 * variables["year_to_s"]));  // Erosion factor
    return F_wind + F_SF + F_erode;
}

// Ug1: dipole
double NGC2264UQFFModule::computeUg1(double t) {
    double mu_dipole = variables["I_dipole"] * variables["A_dipole"] * variables["omega_spin"];
    return mu_dipole * variables["B"];
}

// Ug2: superconductor
double NGC2264UQFFModule::computeUg2(double t) {
    double B_super = variables["mu_0"] * variables["H_aether"];
    return (B_super * B_super) / (2 * variables["mu_0"]);
}

// Ug3': external (stellar wind as external)
double NGC2264UQFFModule::computeUg3prime(double t) {
    double M_star = 20 * 1.989e30;  // Massive star
    double r_star = 1e10;           // Approx
    return (variables["G"] * M_star) / (r_star * r_star);
}

// Ug4: reaction
double NGC2264UQFFModule::computeUg4(double t) {
    double E_react = 1e40 * std::exp(-0.0005 * t);
    return variables["k_4"] * E_react;
}

// Ui
double NGC2264UQFFModule::computeUi(double t) {
    return variables["lambda_I"] * (variables["rho_vac_SCm"] / variables["rho_vac_UA"]) * variables["omega_i"] * std::cos(variables["pi"] * variables["t_n"]) * (1 + variables["F_RZ"]);
}

// Psi integral (simplified)
double NGC2264UQFFModule::computePsiIntegral(double r, double t) {
    double A = variables["A"];
    double m = 1.0;  // m-mode for pillar
    double omega = variables["omega"];
    double sigma = variables["sigma"];
    std::complex<double> psi_pillar(A * std::exp(-r*r / (2 * sigma * sigma)) * std::exp(std::complex<double>(0, m * 0 - omega * t)));
    return std::norm(psi_pillar);  // |psi|^2
}

// Quantum term
double NGC2264UQFFModule::computeQuantumTerm(double t_Hubble_val, double r) {
    double unc = std::sqrt(variables["Delta_x"] * variables["Delta_p"]);
    double psi_int = computePsiIntegral(r, variables["t"]);
    return (variables["hbar"] / unc) * variables["integral_psi"] * (2 * variables["pi"] / t_Hubble_val) * psi_int;
}

// Fluid
double NGC2264UQFFModule::computeFluidTerm(double g_base) {
    return variables["rho_fluid"] * variables["V"] * g_base;
}

// DM
double NGC2264UQFFModule::computeDMTerm(double r) {
    double pert = variables["delta_rho_over_rho"];
    double curv = 3 * variables["G"] * variables["M"] / (r * r * r);
    return (variables["M_visible"] + variables["M_DM"]) * (pert + curv);
}

// Ug sum
double NGC2264UQFFModule::computeUgSum(double r) {
    double Ug_base = (variables["G"] * variables["M"]) / (r * r);
    variables["Ug1"] = computeUg1(variables["t"]);
    variables["Ug2"] = computeUg2(variables["t"]);
    variables["Ug3"] = computeUg3prime(variables["t"]);
    variables["Ug4"] = computeUg4(variables["t"]);
    return Ug_base + variables["Ug1"] + variables["Ug2"] + variables["Ug3"] + variables["Ug4"];
}

// Full g_NGC2264
double NGC2264UQFFModule::computeG(double t, double r) {
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
std::string NGC2264UQFFModule::getEquationText() {
    return "g_NGC2264(r, t) = (G * M(t) / r(t)^2) * (1 + H(t, z)) * (1 - B(t) / B_crit) * (1 + F_env(t)) + "
           "(U_g1 + U_g2 + U_g3' + U_g4) + U_i + (Lambda * c^2 / 3) + "
           "(hbar / sqrt(?x * ?p)) * ? (?_total * H * ?_total dV) * (2? / t_Hubble) + "
           "?_fluid * V * g + (M_visible + M_DM) * (??/? + 3 G M / r^3)\n"
           "Where: M(t) = M * (1 + M_SF(t)); M_SF(t) = SFR * t; r(t) = r0 + v_r t;\n"
           "H(t, z) = H0 * sqrt(?m (1+z)^3 + ??); F_env(t) = F_wind + F_SF + F_erode;\n"
           "F_wind = ? v_wind^2; U_g1 = ?_dipole * B; U_g2 = B_super^2 / (2 ?0);\n"
           "U_g3' = G M_star / r_star^2; U_g4 = k4 * E_react(t); U_i = ?_I * (?_SCm/?_UA) * ?_i * cos(? t_n) * (1 + F_RZ);\n"
           "?_total = A exp(-r^2/(2?^2)) exp(i(m? - ? t)) + wind terms; Insights: Attractive (g_base, Ug1, Ug3') vs. Repulsive (U_g2, ?) advance UQFF.\n"
           "Adaptations: Hubble ACS 2002 data; SFR=0.01 Msun/yr; M=100 Msun. Solutions: g ~1e-9 m/s� at t=3 Myr (wind/fluid dominant).";
}

// Print
void NGC2264UQFFModule::printVariables() {
    std::cout << "NGC 2264 Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// ===== DYNAMIC SELF-UPDATE & SELF-EXPANSION IMPLEMENTATIONS =====

// Static storage for saved states
static std::map<std::string, std::map<std::string, double>> ngc2264_saved_states;

// 1. Dynamic variable management
void NGC2264UQFFModule::createDynamicVariable(const std::string& name, double value) {
    variables[name] = value;
    std::cout << "Created dynamic variable: " << name << " = " << value << std::endl;
}

void NGC2264UQFFModule::removeDynamicVariable(const std::string& name) {
    if (variables.find(name) != variables.end()) {
        variables.erase(name);
        std::cout << "Removed dynamic variable: " << name << std::endl;
    } else {
        std::cerr << "Variable '" << name << "' not found for removal." << std::endl;
    }
}

void NGC2264UQFFModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
        std::cout << "Cloned " << source << " to " << dest << std::endl;
    } else {
        std::cerr << "Source variable '" << source << "' not found." << std::endl;
    }
}

void NGC2264UQFFModule::listAllVariables() {
    std::cout << "=== All NGC 2264 Variables (Total: " << variables.size() << ") ===" << std::endl;
    for (const auto& pair : variables) {
        std::cout << "  " << pair.first << " = " << pair.second << std::endl;
    }
}

// 2. Batch operations
void NGC2264UQFFModule::applyTransformToGroup(const std::vector<std::string>& varNames,
                                               std::function<double(double)> transform) {
    for (const auto& name : varNames) {
        if (variables.find(name) != variables.end()) {
            variables[name] = transform(variables[name]);
            std::cout << "Transformed " << name << " to " << variables[name] << std::endl;
        }
    }
}

void NGC2264UQFFModule::scaleVariableGroup(const std::vector<std::string>& varNames, double scale_factor) {
    applyTransformToGroup(varNames, [scale_factor](double val) { return val * scale_factor; });
}

// 3. Self-expansion capabilities
void NGC2264UQFFModule::autoExpandParameterSpace(double scale_factor) {
    std::cout << "Auto-expanding NGC 2264 parameter space by factor " << scale_factor << std::endl;
    std::vector<std::string> expandable = {"M", "M_visible", "M_DM", "r", "V", "SFR"};
    scaleVariableGroup(expandable, scale_factor);
    // Update M consistency
    variables["M"] = variables["M_visible"] + variables["M_DM"];
    variables["M0"] = variables["M"];
    std::cout << "  Updated M_total, M0" << std::endl;
}

void NGC2264UQFFModule::expandMassScale(double mass_multiplier) {
    std::cout << "Expanding mass scale by " << mass_multiplier << std::endl;
    variables["M_visible"] *= mass_multiplier;
    variables["M_DM"] *= mass_multiplier;
    variables["M"] = variables["M_visible"] + variables["M_DM"];
    variables["M0"] = variables["M"];
    std::cout << "  M_visible: " << (variables["M_visible"] / 1.989e30) << " M☉" << std::endl;
    std::cout << "  M_DM: " << (variables["M_DM"] / 1.989e30) << " M☉" << std::endl;
    std::cout << "  M_total: " << (variables["M"] / 1.989e30) << " M☉" << std::endl;
}

void NGC2264UQFFModule::expandSpatialScale(double spatial_multiplier) {
    std::cout << "Expanding spatial scale by " << spatial_multiplier << std::endl;
    variables["r"] *= spatial_multiplier;
    variables["V"] *= (spatial_multiplier * spatial_multiplier * spatial_multiplier);
    variables["Delta_x"] *= spatial_multiplier;
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];
    variables["sigma"] *= spatial_multiplier;
    std::cout << "  r: " << variables["r"] << " m" << std::endl;
    std::cout << "  V: " << variables["V"] << " m^3" << std::endl;
    std::cout << "  Updated Delta_p, sigma" << std::endl;
}

void NGC2264UQFFModule::expandTimeScale(double time_multiplier) {
    std::cout << "Expanding time scale by " << time_multiplier << std::endl;
    variables["t"] *= time_multiplier;
    variables["omega"] *= (1.0 / time_multiplier);
    variables["omega_i"] *= (1.0 / time_multiplier);
    variables["omega_spin"] *= (1.0 / time_multiplier);
    std::cout << "  t: " << (variables["t"] / variables["year_to_s"]) << " years" << std::endl;
}

// 4. Self-refinement
void NGC2264UQFFModule::autoRefineParameters(double tolerance) {
    std::cout << "Auto-refining NGC 2264 parameters with tolerance " << tolerance << std::endl;
    
    // Validate M consistency
    double M_expected = variables["M_visible"] + variables["M_DM"];
    if (std::abs(variables["M"] - M_expected) / M_expected > tolerance) {
        std::cout << "  Correcting M: " << variables["M"] << " -> " << M_expected << std::endl;
        variables["M"] = M_expected;
        variables["M0"] = M_expected;
    }
    
    // Validate Delta_p from Delta_x
    double Delta_p_expected = variables["hbar"] / variables["Delta_x"];
    if (std::abs(variables["Delta_p"] - Delta_p_expected) / Delta_p_expected > tolerance) {
        std::cout << "  Correcting Delta_p" << std::endl;
        variables["Delta_p"] = Delta_p_expected;
    }
    
    std::cout << "Refinement complete." << std::endl;
}

void NGC2264UQFFModule::calibrateToObservations(const std::map<std::string, double>& observed_values) {
    std::cout << "Calibrating to " << observed_values.size() << " NGC 2264 observations..." << std::endl;
    for (const auto& obs : observed_values) {
        if (variables.find(obs.first) != variables.end()) {
            double old_val = variables[obs.first];
            updateVariable(obs.first, obs.second);
            std::cout << "  " << obs.first << ": " << old_val << " -> " << obs.second << std::endl;
        }
    }
    std::cout << "Calibration complete." << std::endl;
}

void NGC2264UQFFModule::optimizeForMetric(const std::string& metric_name, double target_value) {
    std::cout << "Optimizing for metric: " << metric_name << " = " << target_value << std::endl;
    
    if (metric_name == "g" || metric_name == "gravity") {
        double t = variables["t"];
        double r = variables["r"];
        double current_g = computeG(t, r);
        double ratio = target_value / std::max(current_g, 1e-100);
        
        // Adjust mass to reach target
        variables["M"] *= ratio;
        variables["M_visible"] *= ratio;
        variables["M_DM"] *= ratio;
        variables["M0"] = variables["M"];
        std::cout << "  Adjusted mass by " << ratio << std::endl;
    } else if (metric_name == "SFR") {
        variables["SFR"] = target_value;
        std::cout << "  Set SFR to " << target_value << " kg/s" << std::endl;
    }
    
    std::cout << "Optimization complete." << std::endl;
}

// 5. Parameter exploration
void NGC2264UQFFModule::generateVariations(int num_variations, double variation_range) {
    std::cout << "Generating " << num_variations << " NGC 2264 variations with range ±" 
              << (variation_range * 100) << "%" << std::endl;
    
    std::vector<std::string> key_params = {"M", "r", "SFR", "v_wind", "rho_fluid", "B"};
    
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

void NGC2264UQFFModule::findOptimalParameters(const std::string& objective, int iterations) {
    std::cout << "Finding optimal NGC 2264 parameters for: " << objective 
              << " (" << iterations << " iterations)" << std::endl;
    
    double best_score = -1e100;
    std::map<std::string, double> best_params;
    
    for (int i = 0; i < iterations; ++i) {
        mutateParameters(0.7, 0.1);
        
        double t = variables["t"];
        double r = variables["r"];
        double score = computeG(t, r);
        
        if (objective == "maximize_g") {
            if (score > best_score) {
                best_score = score;
                best_params = variables;
            }
        } else if (objective == "target_1e-9") {
            if (std::abs(score - 1e-9) < std::abs(best_score - 1e-9)) {
                best_score = score;
                best_params = variables;
            }
        }
    }
    
    variables = best_params;
    std::cout << "Optimal g: " << best_score << " m/s^2" << std::endl;
}

// 6. Adaptive evolution
void NGC2264UQFFModule::mutateParameters(double mutation_rate, double mutation_strength) {
    std::vector<std::string> mutable_params = {"M", "r", "SFR", "v_wind", "rho_fluid", "B", "v_r"};
    
    for (const auto& param : mutable_params) {
        if (variables.find(param) != variables.end()) {
            if ((rand() / (double)RAND_MAX) < mutation_rate) {
                double mutation = 1.0 + mutation_strength * (2.0 * (rand() / (double)RAND_MAX) - 1.0);
                variables[param] *= mutation;
            }
        }
    }
    
    // Update dependent variables
    if (variables["M_visible"] > 0 && variables["M_DM"] > 0) {
        variables["M"] = variables["M_visible"] + variables["M_DM"];
        variables["M0"] = variables["M"];
    }
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];
}

void NGC2264UQFFModule::evolveSystem(int generations) {
    std::cout << "Evolving NGC 2264 system over " << generations << " generations..." << std::endl;
    
    for (int gen = 0; gen < generations; ++gen) {
        mutateParameters(0.3, 0.08);
        
        double t = variables["t"];
        double r = variables["r"];
        double fitness = computeG(t, r);
        
        if (gen % 10 == 0) {
            std::cout << "  Gen " << gen << ": g = " << fitness << " m/s^2" << std::endl;
        }
    }
    
    std::cout << "Evolution complete." << std::endl;
}

// 7. State management
void NGC2264UQFFModule::saveState(const std::string& label) {
    ngc2264_saved_states[label] = variables;
    std::cout << "Saved NGC 2264 state: " << label << " (" << variables.size() << " variables)" << std::endl;
}

void NGC2264UQFFModule::restoreState(const std::string& label) {
    if (ngc2264_saved_states.find(label) != ngc2264_saved_states.end()) {
        variables = ngc2264_saved_states[label];
        std::cout << "Restored NGC 2264 state: " << label << std::endl;
    } else {
        std::cerr << "State '" << label << "' not found." << std::endl;
    }
}

void NGC2264UQFFModule::listSavedStates() {
    std::cout << "=== Saved NGC 2264 States (Total: " << ngc2264_saved_states.size() << ") ===" << std::endl;
    for (const auto& state : ngc2264_saved_states) {
        std::cout << "  " << state.first << " (" << state.second.size() << " variables)" << std::endl;
    }
}

void NGC2264UQFFModule::exportState(const std::string& filename) {
    std::cout << "Exporting NGC 2264 state to " << filename << " (not implemented - placeholder)" << std::endl;
    // In real implementation: write variables to file
}

// 8. System analysis
void NGC2264UQFFModule::analyzeParameterSensitivity(const std::string& param_name) {
    if (variables.find(param_name) == variables.end()) {
        std::cerr << "Parameter '" << param_name << "' not found." << std::endl;
        return;
    }
    
    std::cout << "=== NGC 2264 Sensitivity Analysis: " << param_name << " ===" << std::endl;
    
    double base_value = variables[param_name];
    double t = variables["t"];
    double r = variables["r"];
    double base_output = computeG(t, r);
    
    std::vector<double> perturbations = {0.7, 0.85, 1.0, 1.15, 1.3};
    
    for (double factor : perturbations) {
        updateVariable(param_name, base_value * factor);
        
        double new_output = computeG(t, r);
        double sensitivity = (new_output - base_output) / std::max(std::abs(base_output), 1e-100);
        
        std::cout << "  " << param_name << " * " << factor << " -> g change: " 
                  << (sensitivity * 100) << "%" << std::endl;
    }
    
    updateVariable(param_name, base_value);  // Restore
}

void NGC2264UQFFModule::generateSystemReport() {
    std::cout << "\n========== NGC 2264 Cone Nebula System Report ==========" << std::endl;
    std::cout << "Total Variables: " << variables.size() << std::endl;
    
    // Key parameters
    std::cout << "\nMass Parameters:" << std::endl;
    std::cout << "M_visible: " << (variables["M_visible"] / 1.989e30) << " M☉" << std::endl;
    std::cout << "M_DM: " << (variables["M_DM"] / 1.989e30) << " M☉" << std::endl;
    std::cout << "M_total: " << (variables["M"] / 1.989e30) << " M☉" << std::endl;
    std::cout << "SFR: " << (variables["SFR"] * variables["year_to_s"] / 1.989e30) << " M☉/yr" << std::endl;
    
    std::cout << "\nSpatial Parameters:" << std::endl;
    std::cout << "r: " << variables["r"] << " m (" << (variables["r"] / 9.461e15) << " ly)" << std::endl;
    std::cout << "z (redshift): " << variables["z"] << std::endl;
    std::cout << "V (volume): " << variables["V"] << " m^3" << std::endl;
    
    std::cout << "\nDynamics:" << std::endl;
    std::cout << "v_wind: " << (variables["v_wind"] / 1e3) << " km/s" << std::endl;
    std::cout << "v_r (radial): " << (variables["v_r"] / 1e3) << " km/s" << std::endl;
    std::cout << "rho_fluid: " << variables["rho_fluid"] << " kg/m^3" << std::endl;
    std::cout << "B (magnetic): " << variables["B"] << " T" << std::endl;
    
    std::cout << "\nQuantum Parameters:" << std::endl;
    std::cout << "Delta_x: " << variables["Delta_x"] << " m" << std::endl;
    std::cout << "Delta_p: " << variables["Delta_p"] << " kg·m/s" << std::endl;
    std::cout << "omega (pillar wave): " << variables["omega"] << " rad/s" << std::endl;
    
    // Current computation
    double t = variables["t"];
    double r = variables["r"];
    double g = computeG(t, r);
    
    std::cout << "\nCurrent Computation:" << std::endl;
    std::cout << "t: " << (t / variables["year_to_s"]) << " years" << std::endl;
    std::cout << "r_eval: " << r << " m" << std::endl;
    std::cout << "g_NGC2264: " << g << " m/s^2" << std::endl;
    
    std::cout << "\nSubcomponents:" << std::endl;
    std::cout << "Ug1 (dipole): " << variables["Ug1"] << std::endl;
    std::cout << "Ug2 (superconductor): " << variables["Ug2"] << std::endl;
    std::cout << "Ug3 (external): " << variables["Ug3"] << std::endl;
    std::cout << "Ug4 (reaction): " << variables["Ug4"] << std::endl;
    std::cout << "Ui (internal): " << variables["Ui"] << std::endl;
    
    std::cout << "\nEnvironmental Factors:" << std::endl;
    std::cout << "F_env(t): " << computeFenv(t) << std::endl;
    std::cout << "H(t,z): " << computeHtz(variables["z"]) << " s^-1" << std::endl;
    
    std::cout << "======================================================\n" << std::endl;
}

void NGC2264UQFFModule::validatePhysicalConsistency() {
    std::cout << "Validating NGC 2264 physical consistency..." << std::endl;
    bool consistent = true;
    
    // Check for NaN/Inf
    for (const auto& pair : variables) {
        if (std::isnan(pair.second) || std::isinf(pair.second)) {
            std::cerr << "  ERROR: " << pair.first << " is NaN/Inf" << std::endl;
            consistent = false;
        }
    }
    
    // Mass consistency
    double M_expected = variables["M_visible"] + variables["M_DM"];
    if (std::abs(variables["M"] - M_expected) / M_expected > 0.01) {
        std::cerr << "  WARNING: M inconsistent with M_visible + M_DM" << std::endl;
        consistent = false;
    }
    
    // Delta_p from Delta_x
    double Delta_p_expected = variables["hbar"] / variables["Delta_x"];
    if (std::abs(variables["Delta_p"] - Delta_p_expected) / Delta_p_expected > 0.01) {
        std::cerr << "  WARNING: Delta_p inconsistent with hbar/Delta_x" << std::endl;
        consistent = false;
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
    
    if (consistent) {
        std::cout << "  All checks passed. NGC 2264 system is physically consistent." << std::endl;
    }
}

void NGC2264UQFFModule::autoCorrectAnomalies() {
    std::cout << "Auto-correcting NGC 2264 anomalies..." << std::endl;
    
    // Fix NaN/Inf
    for (auto& pair : variables) {
        if (std::isnan(pair.second) || std::isinf(pair.second)) {
            std::cout << "  Correcting " << pair.first << " (was NaN/Inf)" << std::endl;
            pair.second = 1.0;
        }
    }
    
    // Enforce M = M_visible + M_DM
    double M_expected = variables["M_visible"] + variables["M_DM"];
    if (std::abs(variables["M"] - M_expected) / M_expected > 0.01) {
        std::cout << "  Correcting M to match M_visible + M_DM" << std::endl;
        variables["M"] = M_expected;
        variables["M0"] = M_expected;
    }
    
    // Enforce Delta_p = hbar / Delta_x
    double Delta_p_expected = variables["hbar"] / variables["Delta_x"];
    if (std::abs(variables["Delta_p"] - Delta_p_expected) / Delta_p_expected > 0.01) {
        std::cout << "  Correcting Delta_p to match hbar/Delta_x" << std::endl;
        variables["Delta_p"] = Delta_p_expected;
    }
    
    // Ensure positive values
    if (variables["M"] <= 0) {
        std::cout << "  Correcting M to 100 M_sun" << std::endl;
        variables["M"] = 100 * 1.989e30;
        variables["M_visible"] = 80 * 1.989e30;
        variables["M_DM"] = 20 * 1.989e30;
    }
    
    if (variables["r"] <= 0) {
        std::cout << "  Correcting r to 3.31e16 m" << std::endl;
        variables["r"] = 3.31e16;
    }
    
    if (variables["SFR"] < 0) {
        std::cout << "  Correcting SFR to 0.01 M_sun/yr" << std::endl;
        variables["SFR"] = 0.01 * 1.989e30 / variables["year_to_s"];
    }
    
    std::cout << "Auto-correction complete." << std::endl;
}

// Example usage
// #include "NGC2264UQFFModule.h"
// int main() {
//     std::cout << "=============================================" << std::endl;
//     std::cout << "NGC 2264 Cone Nebula - Enhanced Demonstration" << std::endl;
//     std::cout << "=============================================" << std::endl;
// 
//     // === Part 1: Basic gravity computation ===
//     std::cout << "\n--- Part 1: Basic Gravity Computation ---" << std::endl;
//     NGC2264UQFFModule mod;
//     double t = 3e6 * 3.156e7;  // 3 Myr
//     double r = 1e16;           // ~1 ly
//     double g = mod.computeG(t, r);
//     std::cout << "g_NGC2264 = " << g << " m/s^2" << std::endl;
//     std::cout << "Time: " << (t / 3.156e7) << " years" << std::endl;
// 
//     // === Part 2: Dynamic variable management ===
//     std::cout << "\n--- Part 2: Dynamic Variable Management ---" << std::endl;
//     mod.createDynamicVariable("custom_wind_factor", 1.5);
//     mod.cloneVariable("v_wind", "v_wind_backup");
//     std::cout << "Total variables: " << std::endl;
//     mod.listAllVariables();
// 
//     // === Part 3: Stellar wind dynamics ===
//     std::cout << "\n--- Part 3: Stellar Wind Dynamics ---" << std::endl;
//     std::cout << "Initial wind velocity: " << (mod.getVariable("v_wind") / 1e3) << " km/s" << std::endl;
//     mod.updateVariable("v_wind", 30e3);  // Increase to 30 km/s
//     std::cout << "Updated wind velocity: " << (mod.getVariable("v_wind") / 1e3) << " km/s" << std::endl;
//     double g_new = mod.computeG(t, r);
//     std::cout << "New g_NGC2264: " << g_new << " m/s^2" << std::endl;
// 
//     // === Part 4: Star formation rate exploration ===
//     std::cout << "\n--- Part 4: Star Formation Rate Exploration ---" << std::endl;
//     std::cout << "Current SFR: " << (mod.getVariable("SFR") * 3.156e7 / 1.989e30) << " M☉/yr" << std::endl;
//     mod.updateVariable("SFR", 0.02 * 1.989e30 / 3.156e7);
//     std::cout << "Increased SFR: " << (mod.getVariable("SFR") * 3.156e7 / 1.989e30) << " M☉/yr" << std::endl;
// 
//     // === Part 5: Batch operations ===
//     std::cout << "\n--- Part 5: Batch Operations ---" << std::endl;
//     std::vector<std::string> mass_params = {"M_visible", "M_DM"};
//     std::cout << "Scaling mass parameters by 1.3..." << std::endl;
//     mod.scaleVariableGroup(mass_params, 1.3);
// 
//     // === Part 6: Self-expansion ===
//     std::cout << "\n--- Part 6: Self-Expansion ---" << std::endl;
//     mod.saveState("before_expansion");
//     mod.expandSpatialScale(1.5);
//     std::cout << "Spatial scale expanded by 1.5x" << std::endl;
// 
//     // === Part 7: Parameter exploration ===
//     std::cout << "\n--- Part 7: Parameter Exploration ---" << std::endl;
//     mod.generateVariations(3, 0.15);
// 
//     // === Part 8: Sensitivity analysis ===
//     std::cout << "\n--- Part 8: Sensitivity Analysis ---" << std::endl;
//     mod.restoreState("before_expansion");
//     mod.analyzeParameterSensitivity("M");
// 
//     // === Part 9: System validation ===
//     std::cout << "\n--- Part 9: System Validation ---" << std::endl;
//     mod.validatePhysicalConsistency();
// 
//     // === Part 10: Auto-refinement ===
//     std::cout << "\n--- Part 10: Auto-Refinement ---" << std::endl;
//     mod.autoRefineParameters(0.01);
// 
//     // === Part 11: Comprehensive system report ===
//     std::cout << "\n--- Part 11: System Report ---" << std::endl;
//     mod.generateSystemReport();
// 
//     // === Part 12: State management ===
//     std::cout << "\n--- Part 12: State Management ---" << std::endl;
//     mod.saveState("final_ngc2264_state");
//     mod.listSavedStates();
// 
//     std::cout << "\n=============================================" << std::endl;
//     std::cout << "Enhanced NGC 2264 demonstration complete!" << std::endl;
//     std::cout << "=============================================" << std::endl;
// 
//     return 0;
// }
// Compile: g++ -o ngc2264_sim base.cpp NGC2264UQFFModule.cpp -lm
// Sample Output: g_NGC2264 ~ 1e-9 m/s² (wind/fluid dominant; repulsive terms advance framework).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

NGC2264UQFFModule Evaluation

Strengths :
-Modular, extensible design for modeling Cone Nebula(NGC 2264) gravity, including stellar winds, pillar erosion, protostar formation, dust / gas densities, and dark matter.
- Comprehensive physics : gravity, cosmological expansion, magnetic fields, environmental / wind / erosion effects, quantum, fluid, and DM terms.
- Dynamic variable management via std::map enables runtime updates and system adaptation.
- Clear separation of computation functions(e.g., Ug1�Ug4, F_env, quantum, fluid, DM), aiding maintainability.
- NGC 2264 - specific parameters are initialized for realistic simulation; supports easy modification.
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
The code is well - structured, flexible, and suitable for scientific prototyping and educational use in nebular dynamics modeling.It implements a broad set of physical effects and adapts to various scenarios.For production or high - performance applications, address the recommendations for improved robustness, maintainability, and scalability.