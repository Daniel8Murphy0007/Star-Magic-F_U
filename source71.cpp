// NGC1316UQFFModule.h
// Modular C++ implementation of the Master Universal Gravity Equation (MUGE & UQFF Integration) for NGC 1316 (Hubble Spies Cosmic Dust Bunnies) Evolution.
// This module models NGC 1316's gravitational dynamics, incorporating merger history, tidal forces, star cluster disruption, dust lanes, AGN jets/radio lobes, and dark matter.
// Usage: #include "NGC1316UQFFModule.h" in base program; NGC1316UQFFModule mod; mod.computeG(t); mod.updateVariable("M_spiral", new_value);
// Variables in std::map for dynamic updates; supports F_env(t) with tidal and cluster terms; uses rho_dust for fluid.
// Approximations: psi_integral normalized to 1.0; H(t,z) with Omega_m=0.3, Omega_Lambda=0.7; E_react exp decay; dust waves simplified.
// NGC 1316 params: M=5e11 Msun, r=46 kpc, M_spiral=1e10 Msun, d=50 kpc, M_BH=1e8 Msun, M_cluster=1e6 Msun, rho_dust=1e-21 kg/m^3, B=1e-4 T, z=0.005, etc.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef NGC1316_UQFF_MODULE_H
#define NGC1316_UQFF_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>
#include <functional>
#include <vector>
#include <random>
#include <algorithm>

class NGC1316UQFFModule {
private:
    std::map<std::string, double> variables;
    double computeHtz(double z_val);
    double computeFenv(double t);
    double computeMmerge(double t);
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
    double computeRt(double t);

public:
    // Constructor: Initialize with NGC 1316 defaults
    NGC1316UQFFModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computation: g_NGC1316(r, t)
    double computeG(double t, double r);

    // Output descriptive text of the equation
    std::string getEquationText();

    // Print all current variables (for debugging)
    void printVariables();

    // ===== Dynamic Self-Update & Self-Expansion Capabilities =====
    
    // 1. Variable Management (4 methods)
    void createVariable(const std::string& name, double value);
    void removeVariable(const std::string& name);
    void cloneVariable(const std::string& source, const std::string& dest);
    std::vector<std::string> listVariables();

    // 2. Batch Operations (2 methods)
    void transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func);
    void scaleVariableGroup(const std::vector<std::string>& names, double factor);

    // 3. Self-Expansion (4 methods: parameter space + 3 domain-specific scales)
    void expandParameterSpace(const std::vector<std::string>& new_params);
    void expandMassScale(double factor);
    void expandSpatialScale(double factor);
    void expandTimeScale(double factor);

    // 4. Self-Refinement (3 methods)
    void autoRefineParameters(double tolerance);
    void calibrateToObservations(const std::map<std::string, double>& obs_data);
    void optimizeForMetric(std::function<double(NGC1316UQFFModule&)> metric);

    // 5. Parameter Exploration (2 methods)
    std::vector<std::map<std::string, double>> generateVariations(int n_variations);
    std::map<std::string, double> findOptimalParameters(std::function<double(NGC1316UQFFModule&)> objective, int iterations);

    // 6. Adaptive Evolution (2 methods)
    void mutateParameters(double mutation_rate);
    void evolveSystem(int generations, std::function<double(NGC1316UQFFModule&)> fitness);

    // 7. State Management (4 methods)
    void saveState(const std::string& label);
    void restoreState(const std::string& label);
    std::vector<std::string> listSavedStates();
    std::map<std::string, double> exportState();

    // 8. System Analysis (4 methods)
    std::map<std::string, double> sensitivityAnalysis(const std::string& var_name, double delta);
    std::string generateReport();
    bool validateConsistency();
    void autoCorrectAnomalies();
};

#endif // NGC1316_UQFF_MODULE_H

// NGC1316UQFFModule.cpp
#include "NGC1316UQFFModule.h"
#include <complex>

// Constructor: NGC 1316-specific values
NGC1316UQFFModule::NGC1316UQFFModule() {
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

    // NGC 1316 parameters
    variables["M_visible"] = 3.5e11 * M_sun_val;    // kg
    variables["M_DM"] = 1.5e11 * M_sun_val;         // kg
    variables["M"] = variables["M_visible"] + variables["M_DM"];  // Total initial
    variables["M0"] = variables["M"];
    variables["M_spiral"] = 1e10 * M_sun_val;       // kg (merger progenitor)
    variables["d_spiral"] = 50e3 * kpc_val;         // m
    variables["M_BH"] = 1e8 * M_sun_val;            // kg (AGN BH)
    variables["M_cluster"] = 1e6 * M_sun_val;       // kg
    variables["r"] = 46e3 * kpc_val;                // m
    variables["z"] = 0.005;                         // Redshift
    variables["tau_merge"] = 1e9 * variables["year_to_s"];  // s
    variables["t"] = 2e9 * variables["year_to_s"];  // Default t=2 Gyr s

    // Dynamics
    variables["rho_dust"] = 1e-21;                  // kg/m^3
    variables["V"] = 1e51;                          // m^3
    variables["B"] = 1e-4;                          // T (AGN)
    variables["B_crit"] = 1e11;                     // T
    variables["Delta_x"] = 1e-10;                   // m
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];
    variables["integral_psi"] = 1.0;                // Normalized

    // Wave/oscillatory for dust lanes
    variables["A"] = 1e-10;
    variables["k"] = 1e20;
    variables["omega"] = 1e-16;                     // rad/s for dust dynamics
    variables["x"] = 0.0;
    variables["v"] = 1e3;                           // m/s
    variables["sigma"] = 2e3 * kpc_val;             // m for Gaussian

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
    variables["k_cluster"] = 1e-12;                 // N/Msun, adjusted
    variables["omega_spin"] = 1e-3;                 // rad/s BH spin
    variables["I_dipole"] = 1e20;                   // A
    variables["A_dipole"] = 1e15;                   // m^2
    variables["H_aether"] = 1e-5;                   // A/m
    variables["delta_rho_over_rho"] = 1e-5;

    // Scales
    variables["scale_macro"] = 1e-12;
    variables["f_TRZ"] = 0.1;
    variables["f_sc"] = 1.0;
    variables["v_r"] = 1e3;                         // m/s radial velocity
    variables["rho"] = variables["rho_dust"];
}

// Update variable (with dependents)
void NGC1316UQFFModule::updateVariable(const std::string& name, double value) {
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
void NGC1316UQFFModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        variables[name] = delta;
    }
}
void NGC1316UQFFModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute H(t, z)
double NGC1316UQFFModule::computeHtz(double z_val) {
    double Hz_kms = variables["H0"] * std::sqrt(variables["Omega_m"] * std::pow(1.0 + z_val, 3) + variables["Omega_Lambda"]);
    return (Hz_kms * 1e3) / variables["Mpc_to_m"];
}

// M_merge(t)
double NGC1316UQFFModule::computeMmerge(double t) {
    return 1e10 * 1.989e30 * std::exp(-t / variables["tau_merge"]);
}

// r(t)
double NGC1316UQFFModule::computeRt(double t) {
    return variables["r"] + variables["v_r"] * t;
}

// F_env(t)
double NGC1316UQFFModule::computeFenv(double t) {
    double F_tidal = (variables["G"] * variables["M_spiral"]) / (variables["d_spiral"] * variables["d_spiral"]);
    double F_cluster = variables["k_cluster"] * (variables["M_cluster"] / 1.989e30);  // Normalize to m/s^2
    return F_tidal + F_cluster;
}

// Ug1: dipole
double NGC1316UQFFModule::computeUg1(double t) {
    double mu_dipole = variables["I_dipole"] * variables["A_dipole"] * variables["omega_spin"];
    return mu_dipole * variables["B"];
}

// Ug2: superconductor
double NGC1316UQFFModule::computeUg2(double t) {
    double B_super = variables["mu_0"] * variables["H_aether"];
    return (B_super * B_super) / (2 * variables["mu_0"]);
}

// Ug3': external
double NGC1316UQFFModule::computeUg3prime(double t) {
    return (variables["G"] * variables["M_spiral"]) / (variables["d_spiral"] * variables["d_spiral"]);
}

// Ug4: reaction
double NGC1316UQFFModule::computeUg4(double t) {
    double E_react = 1e46 * std::exp(-0.0005 * t);
    return variables["k_4"] * E_react;
}

// Ui
double NGC1316UQFFModule::computeUi(double t) {
    return variables["lambda_I"] * (variables["rho_vac_SCm"] / variables["rho_vac_UA"]) * variables["omega_i"] * std::cos(variables["pi"] * variables["t_n"]) * (1 + variables["F_RZ"]);
}

// Psi integral (simplified)
double NGC1316UQFFModule::computePsiIntegral(double r, double t) {
    double A = variables["A"];
    double m = 2.0;
    double omega = variables["omega"];
    double sigma = variables["sigma"];
    std::complex<double> psi_dust(A * std::exp(-r*r / (2 * sigma * sigma)) * std::exp(std::complex<double>(0, m * 0 - omega * t)));
    return std::norm(psi_dust);  // |psi|^2
}

// Quantum term
double NGC1316UQFFModule::computeQuantumTerm(double t_Hubble_val, double r) {
    double unc = std::sqrt(variables["Delta_x"] * variables["Delta_p"]);
    double psi_int = computePsiIntegral(r, variables["t"]);
    return (variables["hbar"] / unc) * variables["integral_psi"] * (2 * variables["pi"] / t_Hubble_val) * psi_int;
}

// Fluid with rho_dust
double NGC1316UQFFModule::computeFluidTerm(double g_base) {
    return variables["rho_dust"] * variables["V"] * g_base;
}

// DM
double NGC1316UQFFModule::computeDMTerm(double r) {
    double pert = variables["delta_rho_over_rho"];
    double curv = 3 * variables["G"] * variables["M"] / (r * r * r);
    return (variables["M_visible"] + variables["M_DM"]) * (pert + curv);
}

// Ug sum
double NGC1316UQFFModule::computeUgSum(double r) {
    double Ug_base = (variables["G"] * variables["M"]) / (r * r);
    variables["Ug1"] = computeUg1(variables["t"]);
    variables["Ug2"] = computeUg2(variables["t"]);
    variables["Ug3"] = computeUg3prime(variables["t"]);
    variables["Ug4"] = computeUg4(variables["t"]);
    return Ug_base + variables["Ug1"] + variables["Ug2"] + variables["Ug3"] + variables["Ug4"];
}

// Full g_NGC1316
double NGC1316UQFFModule::computeG(double t, double r) {
    variables["t"] = t;
    double m_merge = computeMmerge(t);
    double m_factor = 1.0 + m_merge / variables["M0"];
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
std::string NGC1316UQFFModule::getEquationText() {
    return "g_NGC1316(r, t) = (G * M(t) / r(t)^2) * (1 + H(t, z)) * (1 - B(t) / B_crit) * (1 + F_env(t)) + "
           "(U_g1 + U_g2 + U_g3' + U_g4) + U_i + (Lambda * c^2 / 3) + "
           "(hbar / sqrt(?x * ?p)) * ? (?_total * H * ?_total dV) * (2? / t_Hubble) + "
           "?_dust * V * g + (M_visible + M_DM) * (??/? + 3 G M / r^3)\n"
           "Where: M(t) = M * (1 + M_merge(t)); M_merge(t) = 1e10 Msun * exp(-t/?); r(t) = r0 + v_r t;\n"
           "H(t, z) = H0 * sqrt(?m (1+z)^3 + ??); F_env(t) = F_tidal + F_cluster;\n"
           "F_tidal = G M_spiral / d^2; F_cluster = k_cluster * M_cluster; U_g1 = ?_dipole * B; U_g2 = B_super^2 / (2 ?0);\n"
           "U_g3' = G M_spiral / d^2; U_g4 = k4 * E_react(t); U_i = ?_I * (?_SCm/?_UA) * ?_i * cos(? t_n) * (1 + F_RZ);\n"
           "?_total = A exp(-r^2/(2?^2)) exp(i(m? - ? t)) + BH terms; Insights: Attractive (g_base, Ug1, Ug3') vs. Repulsive (U_g1, U_g2, ?) advance UQFF.\n"
           "Adaptations: Hubble ACS 2003 data; M=5e11 Msun; rho_dust=1e-21 kg/m^3. Solutions: g ~2e37 m/s� at t=2 Gyr (DM/dust dominant).";
}

// Print
void NGC1316UQFFModule::printVariables() {
    std::cout << "NGC 1316 Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// ===== Implementation of Dynamic Self-Update & Self-Expansion Capabilities =====

namespace {
    // Static storage for saved states
    std::map<std::string, std::map<std::string, double>> ngc1316_saved_states;
}

// 1. Variable Management

void NGC1316UQFFModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void NGC1316UQFFModule::removeVariable(const std::string& name) {
    auto it = variables.find(name);
    if (it != variables.end()) {
        variables.erase(it);
    }
}

void NGC1316UQFFModule::cloneVariable(const std::string& source, const std::string& dest) {
    auto it = variables.find(source);
    if (it != variables.end()) {
        variables[dest] = it->second;
    }
}

std::vector<std::string> NGC1316UQFFModule::listVariables() {
    std::vector<std::string> var_names;
    for (const auto& pair : variables) {
        var_names.push_back(pair.first);
    }
    return var_names;
}

// 2. Batch Operations

void NGC1316UQFFModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        auto it = variables.find(name);
        if (it != variables.end()) {
            it->second = func(it->second);
        }
    }
}

void NGC1316UQFFModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// 3. Self-Expansion

void NGC1316UQFFModule::expandParameterSpace(const std::vector<std::string>& new_params) {
    for (const auto& param : new_params) {
        if (variables.find(param) == variables.end()) {
            variables[param] = 0.0;
        }
    }
}

void NGC1316UQFFModule::expandMassScale(double factor) {
    // Scale masses: M, M_visible, M_DM, M_spiral, M_BH, M_cluster
    std::vector<std::string> mass_vars = {"M", "M_visible", "M_DM", "M_spiral", "M_BH", "M_cluster", "M0"};
    scaleVariableGroup(mass_vars, factor);
    // Auto-sync M = M_visible + M_DM
    variables["M"] = variables["M_visible"] + variables["M_DM"];
    variables["M0"] = variables["M"];
}

void NGC1316UQFFModule::expandSpatialScale(double factor) {
    // Scale spatial: r, d_spiral, Delta_x, sigma, V
    std::vector<std::string> spatial_vars = {"r", "d_spiral", "Delta_x", "sigma"};
    scaleVariableGroup(spatial_vars, factor);
    // Auto-sync Delta_p and Volume scales cubically
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];
    variables["V"] *= (factor * factor * factor);
}

void NGC1316UQFFModule::expandTimeScale(double factor) {
    // Scale temporal: t, tau_merge, t_n, t_Hubble
    std::vector<std::string> time_vars = {"t", "tau_merge", "t_n", "t_Hubble"};
    scaleVariableGroup(time_vars, factor);
    // Auto-sync omega values (inverse time)
    if (variables.find("omega") != variables.end()) {
        variables["omega"] /= factor;
    }
    if (variables.find("omega_i") != variables.end()) {
        variables["omega_i"] /= factor;
    }
    if (variables.find("omega_spin") != variables.end()) {
        variables["omega_spin"] /= factor;
    }
}

// 4. Self-Refinement

void NGC1316UQFFModule::autoRefineParameters(double tolerance) {
    // Refine M to match M_visible + M_DM
    double expected_M = variables["M_visible"] + variables["M_DM"];
    if (std::abs(variables["M"] - expected_M) > tolerance) {
        variables["M"] = expected_M;
        variables["M0"] = expected_M;
    }
    
    // Refine Delta_p from Delta_x
    double expected_Delta_p = variables["hbar"] / variables["Delta_x"];
    if (std::abs(variables["Delta_p"] - expected_Delta_p) > tolerance) {
        variables["Delta_p"] = expected_Delta_p;
    }
    
    // Refine rho from rho_dust
    if (std::abs(variables["rho"] - variables["rho_dust"]) > tolerance) {
        variables["rho"] = variables["rho_dust"];
    }
}

void NGC1316UQFFModule::calibrateToObservations(const std::map<std::string, double>& obs_data) {
    for (const auto& obs : obs_data) {
        if (variables.find(obs.first) != variables.end()) {
            variables[obs.first] = obs.second;
        }
    }
    // Auto-sync dependencies
    autoRefineParameters(1e-10);
}

void NGC1316UQFFModule::optimizeForMetric(std::function<double(NGC1316UQFFModule&)> metric) {
    double best_score = metric(*this);
    std::map<std::string, double> best_state = variables;
    
    std::default_random_engine gen;
    std::uniform_real_distribution<double> dist(0.9, 1.1);
    
    for (int iter = 0; iter < 100; iter++) {
        // Mutate key NGC 1316 parameters
        std::vector<std::string> key_params = {"M_visible", "M_DM", "M_spiral", "M_BH", "rho_dust", "B", "tau_merge"};
        for (const auto& param : key_params) {
            if (variables.find(param) != variables.end()) {
                variables[param] *= dist(gen);
            }
        }
        autoRefineParameters(1e-10);
        
        double score = metric(*this);
        if (score > best_score) {
            best_score = score;
            best_state = variables;
        } else {
            variables = best_state;
        }
    }
}

// 5. Parameter Exploration

std::vector<std::map<std::string, double>> NGC1316UQFFModule::generateVariations(int n_variations) {
    std::vector<std::map<std::string, double>> variations;
    std::default_random_engine gen;
    std::uniform_real_distribution<double> dist(0.8, 1.2);
    
    std::map<std::string, double> original = variables;
    std::vector<std::string> vary_params = {"M_visible", "M_DM", "M_spiral", "M_BH", "M_cluster", "rho_dust", "B", "tau_merge", "z"};
    
    for (int i = 0; i < n_variations; i++) {
        for (const auto& param : vary_params) {
            if (variables.find(param) != variables.end()) {
                variables[param] = original[param] * dist(gen);
            }
        }
        autoRefineParameters(1e-10);
        variations.push_back(variables);
    }
    
    variables = original;
    return variations;
}

std::map<std::string, double> NGC1316UQFFModule::findOptimalParameters(std::function<double(NGC1316UQFFModule&)> objective, int iterations) {
    double best_score = objective(*this);
    std::map<std::string, double> best_params = variables;
    
    std::default_random_engine gen;
    std::uniform_real_distribution<double> dist(0.5, 1.5);
    
    for (int iter = 0; iter < iterations; iter++) {
        std::vector<std::string> opt_params = {"M_visible", "M_DM", "M_spiral", "M_BH", "rho_dust", "B", "omega_spin"};
        for (const auto& param : opt_params) {
            if (variables.find(param) != variables.end()) {
                variables[param] *= dist(gen);
            }
        }
        autoRefineParameters(1e-10);
        
        double score = objective(*this);
        if (score > best_score) {
            best_score = score;
            best_params = variables;
        }
    }
    
    variables = best_params;
    return best_params;
}

// 6. Adaptive Evolution

void NGC1316UQFFModule::mutateParameters(double mutation_rate) {
    std::default_random_engine gen;
    std::uniform_real_distribution<double> dist(-mutation_rate, mutation_rate);
    
    std::vector<std::string> mutable_params = {"M_visible", "M_DM", "M_spiral", "M_BH", "M_cluster", 
                                                 "rho_dust", "B", "tau_merge", "z", "omega_spin", 
                                                 "k_cluster", "I_dipole", "v_r"};
    for (const auto& param : mutable_params) {
        if (variables.find(param) != variables.end()) {
            double mutation = 1.0 + dist(gen);
            variables[param] *= mutation;
        }
    }
    autoRefineParameters(1e-10);
}

void NGC1316UQFFModule::evolveSystem(int generations, std::function<double(NGC1316UQFFModule&)> fitness) {
    for (int gen = 0; gen < generations; gen++) {
        double current_fitness = fitness(*this);
        std::map<std::string, double> current_state = variables;
        
        mutateParameters(0.1);
        
        double new_fitness = fitness(*this);
        if (new_fitness < current_fitness) {
            variables = current_state;  // Revert if fitness decreased
        }
    }
}

// 7. State Management

void NGC1316UQFFModule::saveState(const std::string& label) {
    ngc1316_saved_states[label] = variables;
}

void NGC1316UQFFModule::restoreState(const std::string& label) {
    auto it = ngc1316_saved_states.find(label);
    if (it != ngc1316_saved_states.end()) {
        variables = it->second;
    }
}

std::vector<std::string> NGC1316UQFFModule::listSavedStates() {
    std::vector<std::string> state_labels;
    for (const auto& pair : ngc1316_saved_states) {
        state_labels.push_back(pair.first);
    }
    return state_labels;
}

std::map<std::string, double> NGC1316UQFFModule::exportState() {
    return variables;
}

// 8. System Analysis

std::map<std::string, double> NGC1316UQFFModule::sensitivityAnalysis(const std::string& var_name, double delta) {
    std::map<std::string, double> sensitivity;
    
    auto it = variables.find(var_name);
    if (it == variables.end()) {
        return sensitivity;
    }
    
    double original_val = it->second;
    double base_g = computeG(variables["t"], variables["r"]);
    
    // Test sensitivity at multiple radii
    std::vector<double> radii = {10e3 * 3.086e19, 20e3 * 3.086e19, 30e3 * 3.086e19, 50e3 * 3.086e19};
    for (size_t i = 0; i < radii.size(); i++) {
        variables[var_name] = original_val * (1.0 + delta);
        autoRefineParameters(1e-10);
        double g_plus = computeG(variables["t"], radii[i]);
        
        variables[var_name] = original_val * (1.0 - delta);
        autoRefineParameters(1e-10);
        double g_minus = computeG(variables["t"], radii[i]);
        
        double sens = (g_plus - g_minus) / (2.0 * delta * original_val);
        sensitivity["r_" + std::to_string(i)] = sens;
    }
    
    variables[var_name] = original_val;
    autoRefineParameters(1e-10);
    return sensitivity;
}

std::string NGC1316UQFFModule::generateReport() {
    std::ostringstream report;
    report << "===== NGC 1316 (Fornax A) UQFF Module Report =====\n";
    report << std::scientific;
    report << "Total Mass: M = " << variables["M"] << " kg (" << variables["M"]/1.989e30 << " Msun)\n";
    report << "  M_visible = " << variables["M_visible"] << " kg\n";
    report << "  M_DM = " << variables["M_DM"] << " kg\n";
    report << "Merger Progenitor: M_spiral = " << variables["M_spiral"] << " kg (" << variables["M_spiral"]/1.989e30 << " Msun)\n";
    report << "  Distance: d_spiral = " << variables["d_spiral"] << " m (" << variables["d_spiral"]/3.086e19 << " kpc)\n";
    report << "AGN Black Hole: M_BH = " << variables["M_BH"] << " kg (" << variables["M_BH"]/1.989e30 << " Msun)\n";
    report << "  Spin: omega_spin = " << variables["omega_spin"] << " rad/s\n";
    report << "  B-field: B = " << variables["B"] << " T\n";
    report << "Star Cluster: M_cluster = " << variables["M_cluster"] << " kg (" << variables["M_cluster"]/1.989e30 << " Msun)\n";
    report << "Dust: rho_dust = " << variables["rho_dust"] << " kg/m^3\n";
    report << "Spatial: r = " << variables["r"] << " m (" << variables["r"]/3.086e19 << " kpc)\n";
    report << "Temporal: t = " << variables["t"] << " s (" << variables["t"]/3.156e7 << " yr)\n";
    report << "  Merger timescale: tau_merge = " << variables["tau_merge"] << " s (" << variables["tau_merge"]/3.156e7 << " yr)\n";
    report << "Redshift: z = " << variables["z"] << "\n";
    report << "Ug components: Ug1=" << variables["Ug1"] << ", Ug2=" << variables["Ug2"] 
           << ", Ug3=" << variables["Ug3"] << ", Ug4=" << variables["Ug4"] << "\n";
    report << "Ui = " << variables["Ui"] << "\n";
    
    double g_test = computeG(variables["t"], variables["r"]);
    report << "Sample g(r=" << variables["r"]/3.086e19 << " kpc, t=" << variables["t"]/3.156e7 
           << " yr) = " << g_test << " m/s^2\n";
    report << "Saved states: " << ngc1316_saved_states.size() << "\n";
    report << "===============================================\n";
    return report.str();
}

bool NGC1316UQFFModule::validateConsistency() {
    bool valid = true;
    
    // Check M = M_visible + M_DM
    double expected_M = variables["M_visible"] + variables["M_DM"];
    if (std::abs(variables["M"] - expected_M) > 1e10) {
        valid = false;
    }
    
    // Check Delta_p = hbar / Delta_x
    double expected_Delta_p = variables["hbar"] / variables["Delta_x"];
    if (std::abs(variables["Delta_p"] - expected_Delta_p) > 1e-40) {
        valid = false;
    }
    
    // Check rho = rho_dust
    if (std::abs(variables["rho"] - variables["rho_dust"]) > 1e-25) {
        valid = false;
    }
    
    // Check M0 = M
    if (std::abs(variables["M0"] - variables["M"]) > 1e10) {
        valid = false;
    }
    
    // Physical bounds
    if (variables["M"] <= 0 || variables["M_BH"] <= 0 || variables["r"] <= 0 || 
        variables["rho_dust"] < 0 || variables["B"] < 0) {
        valid = false;
    }
    
    return valid;
}

void NGC1316UQFFModule::autoCorrectAnomalies() {
    // Correct M
    variables["M"] = variables["M_visible"] + variables["M_DM"];
    variables["M0"] = variables["M"];
    
    // Correct Delta_p
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];
    
    // Correct rho
    variables["rho"] = variables["rho_dust"];
    
    // Enforce physical bounds
    if (variables["M"] <= 0) variables["M"] = 1e11 * 1.989e30;
    if (variables["M_visible"] <= 0) variables["M_visible"] = 3.5e11 * 1.989e30;
    if (variables["M_DM"] <= 0) variables["M_DM"] = 1.5e11 * 1.989e30;
    if (variables["M_BH"] <= 0) variables["M_BH"] = 1e8 * 1.989e30;
    if (variables["r"] <= 0) variables["r"] = 46e3 * 3.086e19;
    if (variables["rho_dust"] < 0) variables["rho_dust"] = 1e-21;
    if (variables["B"] < 0) variables["B"] = 1e-4;
    if (variables["tau_merge"] <= 0) variables["tau_merge"] = 1e9 * 3.156e7;
}

// Example usage
// #include "NGC1316UQFFModule.h"
// int main() {
//     NGC1316UQFFModule mod;
//     double t = 2e9 * 3.156e7;  // 2 Gyr
//     double r = 20e3 * 3.086e19;  // 20 kpc
//     double g = mod.computeG(t, r);
//     std::cout << "g_NGC1316 = " << g << " m/s²\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("M_spiral", 1.5e10 * 1.989e30);
//     mod.printVariables();
//
//     // ===== Demonstrate Dynamic Self-Update & Self-Expansion =====
//     
//     // 1. Variable management
//     mod.createVariable("F_jets", 1e-15);
//     mod.cloneVariable("M_BH", "M_BH_backup");
//     std::cout << "Variables: " << mod.listVariables().size() << " total\n";
//     
//     // 2. Batch operations on NGC 1316 masses
//     std::vector<std::string> mass_group = {"M_visible", "M_DM", "M_spiral"};
//     mod.scaleVariableGroup(mass_group, 1.1);  // 10% mass increase from merger
//     
//     // 3. Self-expansion
//     mod.expandMassScale(1.2);  // Post-merger 20% mass growth
//     mod.expandSpatialScale(1.1);  // Tidal expansion
//     mod.expandTimeScale(1.05);  // Temporal scaling
//     std::cout << "After expansion: M = " << mod.exportState()["M"]/1.989e30 << " Msun\n";
//     
//     // 4. Self-refinement
//     mod.autoRefineParameters(1e-10);
//     std::map<std::string, double> obs = {{"M_BH", 1.2e8*1.989e30}, {"z", 0.0059}};  // Updated observations
//     mod.calibrateToObservations(obs);
//     
//     // 5. Parameter exploration (find optimal AGN configuration)
//     auto agn_objective = [](NGC1316UQFFModule& m) {
//         double g = m.computeG(2e9*3.156e7, 5e3*3.086e19);  // Inner 5 kpc
//         return -std::abs(g - 1e-9);  // Target specific AGN gravity signature
//     };
//     mod.optimizeForMetric(agn_objective);
//     
//     // 6. Generate merger scenario variations
//     auto variations = mod.generateVariations(5);
//     std::cout << "Generated " << variations.size() << " merger scenarios\n";
//     
//     // 7. State management
//     mod.saveState("post_merger_config");
//     mod.expandMassScale(0.9);  // Test different scenario
//     mod.restoreState("post_merger_config");  // Revert to saved
//     std::cout << "Saved states: " << mod.listSavedStates().size() << "\n";
//     
//     // 8. Sensitivity analysis for dust density
//     auto dust_sensitivity = mod.sensitivityAnalysis("rho_dust", 0.1);
//     std::cout << "Dust sensitivity at various radii: " << dust_sensitivity.size() << " samples\n";
//     
//     // 9. System validation
//     bool valid = mod.validateConsistency();
//     std::cout << "System consistency: " << (valid ? "VALID" : "INVALID") << "\n";
//     if (!valid) mod.autoCorrectAnomalies();
//     
//     // 10. Comprehensive report
//     std::cout << mod.generateReport();
//     
//     // 11. Adaptive evolution (optimize for merger dynamics over 10 generations)
//     auto merger_fitness = [](NGC1316UQFFModule& m) {
//         double g1 = m.computeG(1e9*3.156e7, 20e3*3.086e19);  // 1 Gyr, 20 kpc
//         double g2 = m.computeG(3e9*3.156e7, 30e3*3.086e19);  // 3 Gyr, 30 kpc
//         return -(std::abs(g1 - 5e-10) + std::abs(g2 - 3e-10));  // Dual-epoch targets
//     };
//     mod.evolveSystem(10, merger_fitness);
//     std::cout << "Evolved system over 10 generations\n";
//     
//     // 12. Final state export
//     auto final_state = mod.exportState();
//     std::cout << "Final M_spiral = " << final_state["M_spiral"]/1.989e30 << " Msun\n";
//     std::cout << "Final M_BH = " << final_state["M_BH"]/1.989e30 << " Msun\n";
//     std::cout << "Final rho_dust = " << final_state["rho_dust"] << " kg/m^3\n";
//
//     return 0;
// }
// Compile: g++ -o ngc1316_sim base.cpp NGC1316UQFFModule.cpp -lm
// Sample Output: g_NGC1316 ~ 2e37 m/s� (DM/fluid dominant; repulsive terms advance framework).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

NGC1316UQFFModule Evaluation

Strengths :
-Modular, extensible design for modeling NGC 1316 galaxy gravity, including merger history, tidal forces, star cluster disruption, dust lanes, AGN jets, and dark matter.
- Comprehensive physics : gravity, cosmological expansion, magnetic fields, environmental / tidal effects, quantum, fluid(dust), and DM terms.
- Dynamic variable management via std::map enables runtime updates and system adaptation.
- Clear separation of computation functions(e.g., Ug1�Ug4, F_env, quantum, fluid, DM), aiding maintainability.
- NGC 1316 - specific parameters are initialized for realistic simulation; supports easy modification.
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