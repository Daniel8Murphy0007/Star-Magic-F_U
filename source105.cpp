// GalacticBlackHoleModule.h
// Modular C++ implementation of the Mass of the Galactic Black Hole (M_bh) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes M_bh=8.15e36 kg ?4.1e6 M_sun; scales M_bh / d_g in Universal Buoyancy U_bi and Ug4.
// Pluggable: #include "GalacticBlackHoleModule.h"
// GalacticBlackHoleModule mod; mod.computeU_b1(); mod.updateVariable("M_bh", new_value);
// Variables in std::map; example for Sun at t=0, t_n=0.
// Approximations: cos(? t_n)=1; (1 + ?_sw ?_vac,sw)?1; ?=0.001 s^-1; f_feedback=0.1.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef GALACTIC_BLACK_HOLE_MODULE_H
#define GALACTIC_BLACK_HOLE_MODULE_H

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

class GalacticBlackHoleModule {
private:
    std::map<std::string, double> variables;
    double computeM_bhInMsun();
    double computeMbhOverDg();
    double computeU_b1();
    double computeU_g4();

public:
    // Constructor: Initialize with framework defaults
    GalacticBlackHoleModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeM_bh();  // 8.15e36 kg
    double computeM_bhInMsun();  // ?4.1e6 M_sun
    double computeMbhOverDg();  // M_bh / d_g (kg/m)
    double computeU_b1();  // Universal Buoyancy example (J/m^3)
    double computeU_g4();  // Ug4 example (J/m^3)

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
    void expandBlackHoleScale(double M_bh_factor, double feedback_factor);  // M_bh and f_feedback
    void expandGalacticScale(double d_g_factor, double Omega_g_factor);     // Distance and rotation
    void expandBuoyancyScale(double beta_factor, double U_g1_factor);       // Buoyancy params

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

#endif // GALACTIC_BLACK_HOLE_MODULE_H

// GalacticBlackHoleModule.cpp
#include "GalacticBlackHoleModule.h"

// Constructor: Set framework defaults
GalacticBlackHoleModule::GalacticBlackHoleModule() {
    // Universal constants
    variables["M_sun"] = 1.989e30;                  // kg
    variables["M_bh"] = 8.15e36;                    // kg (Sgr A*)

    // Shared params for terms
    variables["beta_1"] = 0.6;                      // Unitless
    variables["U_g1"] = 1.39e26;                    // J/m^3
    variables["Omega_g"] = 7.3e-16;                 // rad/s
    variables["d_g"] = 2.55e20;                     // m
    variables["epsilon_sw"] = 0.001;                // Unitless
    variables["rho_vac_sw"] = 8e-21;                // J/m^3
    variables["U_UA"] = 1.0;                        // Normalized
    variables["t_n"] = 0.0;                         // s
    variables["pi"] = 3.141592653589793;

    // Ug4 params
    variables["k_4"] = 1.0;                         // Unitless
    variables["rho_vac_SCm"] = 7.09e-37;            // J/m^3
    variables["alpha"] = 0.001;                     // s^-1
    variables["f_feedback"] = 0.1;                  // Unitless
}

// Update variable
void GalacticBlackHoleModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void GalacticBlackHoleModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void GalacticBlackHoleModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute M_bh (kg)
double GalacticBlackHoleModule::computeM_bh() {
    return variables["M_bh"];
}

// M_bh in M_sun
double GalacticBlackHoleModule::computeM_bhInMsun() {
    return computeM_bh() / variables["M_sun"];
}

// M_bh / d_g (kg/m)
double GalacticBlackHoleModule::computeMbhOverDg() {
    return computeM_bh() / variables["d_g"];
}

// U_b1 example (J/m^3)
double GalacticBlackHoleModule::computeU_b1() {
    double beta_1 = variables["beta_1"];
    double U_g1 = variables["U_g1"];
    double Omega_g = variables["Omega_g"];
    double mbh_over_dg = computeMbhOverDg();
    double swirl_factor = 1.0 + variables["epsilon_sw"] * variables["rho_vac_sw"];
    double U_UA = variables["U_UA"];
    double cos_term = std::cos(variables["pi"] * variables["t_n"]);
    return -beta_1 * U_g1 * Omega_g * mbh_over_dg * swirl_factor * U_UA * cos_term;
}

// U_g4 example (J/m^3)
double GalacticBlackHoleModule::computeU_g4() {
    double k_4 = variables["k_4"];
    double rho_vac_SCm = variables["rho_vac_SCm"];
    double mbh_over_dg = computeMbhOverDg();
    double exp_term = std::exp( - variables["alpha"] * variables["t_n"] );
    double cos_term = std::cos(variables["pi"] * variables["t_n"]);
    double feedback_factor = 1.0 + variables["f_feedback"];
    return k_4 * (rho_vac_SCm * computeM_bh() / variables["d_g"]) * exp_term * cos_term * feedback_factor;
}

// Equation text
std::string GalacticBlackHoleModule::getEquationText() {
    return "U_bi = -?_i U_gi ?_g (M_bh / d_g) (1 + ?_sw ?_vac,sw) U_UA cos(? t_n)\n"
           "U_g4 = k_4 (?_vac,[SCm] M_bh / d_g) e^{-? t} cos(? t_n) (1 + f_feedback)\n"
           "Where M_bh = 8.15e36 kg ?4.1e6 M_sun (Sgr A*).\n"
           "M_bh / d_g ?3.20e16 kg/m;\n"
           "Example U_b1 ? -1.94e27 J/m�; U_g4 ?2.50e-20 J/m� (t_n=0).\n"
           "Role: Scales SMBH gravity in buoyancy/Ug4; drives galactic dynamics/mergers.\n"
           "UQFF: Central mass for star formation/nebulae; resolves parsec problem.";
}

// Print variables
void GalacticBlackHoleModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// ===== Implementation: 25-Method Dynamic Self-Update & Self-Expansion Capabilities =====

namespace galactic_black_hole_saved_states {
    std::map<std::string, std::map<std::string, double>> saved_states;
}

// Variable Management (5 methods)
void GalacticBlackHoleModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void GalacticBlackHoleModule::removeVariable(const std::string& name) {
    variables.erase(name);
}

void GalacticBlackHoleModule::cloneVariable(const std::string& source, const std::string& target) {
    if (variables.find(source) != variables.end()) {
        variables[target] = variables[source];
    }
}

std::string GalacticBlackHoleModule::listVariables() {
    std::ostringstream oss;
    for (const auto& pair : variables) {
        oss << pair.first << " = " << pair.second << "\n";
    }
    return oss.str();
}

std::string GalacticBlackHoleModule::getSystemName() {
    return "Galactic_Black_Hole_UQFF";
}

// Batch Operations (2 methods)
void GalacticBlackHoleModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = func(variables[name]);
        }
    }
}

void GalacticBlackHoleModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double x) { return x * factor; });
}

// Self-Expansion (4 methods)
void GalacticBlackHoleModule::expandParameterSpace(const std::vector<std::string>& params, double expansion_factor) {
    for (const auto& param : params) {
        if (variables.find(param) != variables.end()) {
            variables[param] *= expansion_factor;
        }
    }
}

void GalacticBlackHoleModule::expandBlackHoleScale(double M_bh_factor, double feedback_factor) {
    // M_bh: black hole mass (kg, ~8.15e36 for Sgr A*)
    // f_feedback: feedback factor (unitless, ~0.1)
    variables["M_bh"] *= M_bh_factor;
    variables["f_feedback"] *= feedback_factor;
}

void GalacticBlackHoleModule::expandGalacticScale(double d_g_factor, double Omega_g_factor) {
    // d_g: galactic distance to Sgr A* (m, ~2.55e20 = 27,000 ly)
    // Omega_g: galactic rotation rate (rad/s, ~7.3e-16)
    variables["d_g"] *= d_g_factor;
    variables["Omega_g"] *= Omega_g_factor;
}

void GalacticBlackHoleModule::expandBuoyancyScale(double beta_factor, double U_g1_factor) {
    // β_1: buoyancy coupling (unitless, ~0.6)
    // U_g1: internal dipole gravity (J/m³, ~1.39e26)
    variables["beta_1"] *= beta_factor;
    variables["U_g1"] *= U_g1_factor;
}

// Self-Refinement (3 methods)
void GalacticBlackHoleModule::autoRefineParameters(const std::string& target_metric, double target_value) {
    // Example: target U_b1 by adjusting M_bh and Omega_g
    if (target_metric == "U_b1") {
        double current = std::abs(computeU_b1());
        if (current > 0) {
            double ratio = std::abs(target_value) / current;
            variables["M_bh"] *= std::sqrt(ratio);
            variables["Omega_g"] *= std::sqrt(ratio);
        }
    } else if (target_metric == "U_g4") {
        double current = computeU_g4();
        if (current > 0) {
            double ratio = target_value / current;
            variables["M_bh"] *= std::sqrt(ratio);
            variables["k_4"] *= std::sqrt(ratio);
        }
    } else if (target_metric == "M_bh_Msun") {
        double target_kg = target_value * variables["M_sun"];
        variables["M_bh"] = target_kg;
    }
}

void GalacticBlackHoleModule::calibrateToObservations(const std::map<std::string, double>& observed) {
    for (const auto& obs : observed) {
        if (variables.find(obs.first) != variables.end()) {
            variables[obs.first] = obs.second;
        }
    }
}

void GalacticBlackHoleModule::optimizeForMetric(const std::string& metric) {
    if (metric == "maximize_U_b1") {
        variables["M_bh"] *= 1.2;
        variables["Omega_g"] *= 1.15;
        variables["beta_1"] *= 1.1;
    } else if (metric == "maximize_U_g4") {
        variables["M_bh"] *= 1.2;
        variables["k_4"] *= 1.15;
        variables["f_feedback"] *= 1.1;
    } else if (metric == "enhance_feedback") {
        variables["f_feedback"] *= 1.5;
    } else if (metric == "increase_mass") {
        variables["M_bh"] *= 1.5; // 1.5 dex mass increase
    }
}

// Parameter Exploration (1 method)
std::vector<std::map<std::string, double>> GalacticBlackHoleModule::generateVariations(int count, double variance) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-variance, variance);

    for (int i = 0; i < count; ++i) {
        std::map<std::string, double> variant = variables;
        for (auto& pair : variant) {
            if (pair.first != "pi" && pair.first != "t_n" && pair.first != "M_sun" && pair.first != "U_UA") {
                double factor = 1.0 + dis(gen);
                pair.second *= factor;
                if (pair.second < 0 && pair.first != "beta_1") {
                    pair.second = std::abs(pair.second);
                }
            }
        }
        variations.push_back(variant);
    }
    return variations;
}

// Adaptive Evolution (2 methods)
void GalacticBlackHoleModule::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-mutation_rate, mutation_rate);

    for (auto& pair : variables) {
        if (pair.first != "pi" && pair.first != "t_n" && pair.first != "M_sun" && pair.first != "U_UA") {
            pair.second *= (1.0 + dis(gen));
            if (pair.second < 0 && pair.first != "beta_1") {
                pair.second = std::abs(pair.second);
            }
        }
    }
}

void GalacticBlackHoleModule::evolveSystem(int generations, std::function<double()> fitness_func) {
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
void GalacticBlackHoleModule::saveState(const std::string& label) {
    galactic_black_hole_saved_states::saved_states[label] = variables;
}

void GalacticBlackHoleModule::restoreState(const std::string& label) {
    if (galactic_black_hole_saved_states::saved_states.find(label) != galactic_black_hole_saved_states.end()) {
        variables = galactic_black_hole_saved_states::saved_states[label];
    }
}

std::vector<std::string> GalacticBlackHoleModule::listSavedStates() {
    std::vector<std::string> labels;
    for (const auto& pair : galactic_black_hole_saved_states::saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string GalacticBlackHoleModule::exportState() {
    std::ostringstream oss;
    oss << "System: " << getSystemName() << "\n";
    oss << listVariables();
    return oss.str();
}

// System Analysis (4 methods)
std::map<std::string, double> GalacticBlackHoleModule::sensitivityAnalysis(const std::vector<std::string>& params) {
    std::map<std::string, double> sensitivity;
    double baseline_ub1 = std::abs(computeU_b1());

    for (const auto& param : params) {
        if (variables.find(param) != variables.end()) {
            double original = variables[param];
            variables[param] *= 1.01; // +1% perturbation
            double perturbed_ub1 = std::abs(computeU_b1());
            variables[param] = original;
            if (baseline_ub1 > 0) {
                sensitivity[param] = std::abs((perturbed_ub1 - baseline_ub1) / baseline_ub1);
            } else {
                sensitivity[param] = 0.0;
            }
        }
    }
    return sensitivity;
}

std::string GalacticBlackHoleModule::generateReport() {
    std::ostringstream oss;
    oss << "========== Galactic Black Hole Module Report ==========\n";
    oss << "System: " << getSystemName() << "\n\n";
    oss << "Key Parameters:\n";
    oss << "  M_bh (Sgr A*): " << std::scientific << variables["M_bh"] << " kg\n";
    oss << "  M_bh (solar masses): " << computeM_bhInMsun() << " M_sun (~4.1 million)\n";
    oss << "  d_g (distance): " << variables["d_g"] << " m (~27,000 ly)\n";
    oss << "  M_bh/d_g: " << computeMbhOverDg() << " kg/m\n";
    oss << "  Ω_g (rotation): " << variables["Omega_g"] << " rad/s\n";
    oss << "  f_feedback: " << variables["f_feedback"] << " (unitless)\n";
    oss << "  β_1 (buoyancy): " << variables["beta_1"] << "\n";
    oss << "  k_4 (Ug4 coupling): " << variables["k_4"] << "\n\n";
    
    oss << "Computed Values at t_n=0:\n";
    double ub1 = computeU_b1();
    double ug4 = computeU_g4();
    oss << "  U_b1 (buoyancy): " << ub1 << " J/m³\n";
    oss << "  U_g4 (star-BH interaction): " << ug4 << " J/m³\n";
    oss << "  |U_b1| / U_g1: " << (std::abs(ub1) / variables["U_g1"]) << "\n\n";
    
    oss << "Mass Comparisons:\n";
    oss << "  M_bh / M_sun: " << computeM_bhInMsun() << " (supermassive)\n";
    oss << "  Typical SMBH: 10⁶ - 10⁹ M_sun\n";
    oss << "  Sgr A*: ~4.1×10⁶ M_sun (confirmed observations)\n\n";
    
    oss << "Physical Context:\n";
    oss << "  Sagittarius A* at galactic center (~8 kpc from Sun)\n";
    oss << "  M_bh/d_g scales SMBH influence on buoyancy and Ug4\n";
    oss << "  U_b1: Opposes gravity via galactic rotation (negative sign)\n";
    oss << "  U_g4: Star-black hole gravitational interactions\n";
    oss << "  f_feedback: AGN feedback from BH accretion/jets\n";
    oss << "  Resolves parsec problem, star formation regulation\n";
    oss << "  Critical for galactic dynamics, mergers, quasar jets\n";
    oss << "========================================================\n";
    return oss.str();
}

bool GalacticBlackHoleModule::validateConsistency() {
    bool valid = true;
    if (variables["M_bh"] <= 0) valid = false;
    if (variables["d_g"] <= 0) valid = false;
    if (variables["Omega_g"] < 0) valid = false;
    if (variables["M_sun"] <= 0) valid = false;
    if (computeM_bhInMsun() < 1e5 || computeM_bhInMsun() > 1e10) {
        // Warning: unusual SMBH mass
        std::cerr << "Warning: M_bh = " << computeM_bhInMsun() << " M_sun is outside typical SMBH range (1e5 - 1e10)\n";
    }
    return valid;
}

void GalacticBlackHoleModule::autoCorrectAnomalies() {
    if (variables["M_bh"] <= 0) variables["M_bh"] = 8.15e36;
    if (variables["d_g"] <= 0) variables["d_g"] = 2.55e20;
    if (variables["Omega_g"] < 0) variables["Omega_g"] = 7.3e-16;
    if (variables["M_sun"] <= 0) variables["M_sun"] = 1.989e30;
    if (variables["beta_1"] < 0) variables["beta_1"] = 0.6;
    if (variables["k_4"] <= 0) variables["k_4"] = 1.0;
    if (variables["f_feedback"] < 0) variables["f_feedback"] = 0.1;
}

// Example usage in base program (snippet)
// #include "GalacticBlackHoleModule.h"
// int main() {
//     GalacticBlackHoleModule mod;
//     double m_sun = mod.computeM_bhInMsun();
//     std::cout << "M_bh ? " << m_sun << " M_sun\n";
//     double u_b1 = mod.computeU_b1();
//     std::cout << "U_b1 = " << u_b1 << " J/m�\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("M_bh", 9e36);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o bh_mass_test bh_mass_test.cpp GalacticBlackHoleModule.cpp -lm
// Sample: M_bh ?4.1e6 M_sun; U_b1 ? -1.94e27 J/m�; scales SMBH influence.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

// ===== Enhanced Example: 18-Step Demonstration of Dynamic Capabilities =====
// int main() {
//     GalacticBlackHoleModule mod;
//     std::cout << "===== Galactic Black Hole Module: Enhanced 18-Step Demo =====\n\n";
//
//     // Step 1: Initial report
//     std::cout << "Step 1: Initial Configuration (Sgr A*)\n";
//     std::cout << mod.generateReport() << "\n";
//
//     // Step 2: Track key SMBH quantities
//     std::cout << "Step 2: Create Tracking Variables\n";
//     double m_msun = mod.computeM_bhInMsun();
//     double mbh_dg = mod.computeMbhOverDg();
//     double ub1_base = mod.computeU_b1();
//     double ug4_base = mod.computeU_g4();
//     mod.createVariable("M_bh_Msun_baseline", m_msun);
//     mod.createVariable("Mbh_over_dg_baseline", mbh_dg);
//     mod.createVariable("U_b1_baseline", ub1_base);
//     mod.createVariable("U_g4_baseline", ug4_base);
//     mod.createVariable("schwarzschild_radius", 2 * 6.674e-11 * mod.variables["M_bh"] / (3e8 * 3e8)); // Rs = 2GM/c²
//     std::cout << "  M_bh: " << std::scientific << m_msun << " M_sun\n";
//     std::cout << "  M_bh/d_g: " << mbh_dg << " kg/m\n";
//     std::cout << "  U_b1: " << ub1_base << " J/m³\n";
//     std::cout << "  U_g4: " << ug4_base << " J/m³\n";
//     std::cout << "  Schwarzschild radius: " << mod.variables["schwarzschild_radius"] << " m\n\n";
//
//     // Step 3: M_bh mass variations
//     std::cout << "Step 3: M_bh Mass Scaling Effects\n";
//     mod.saveState("baseline");
//     std::vector<double> mass_factors = {0.5, 0.8, 1.0, 1.5, 2.0, 3.0};
//     for (double factor : mass_factors) {
//         mod.updateVariable("M_bh", 8.15e36 * factor);
//         double msun = mod.computeM_bhInMsun();
//         double ub1 = mod.computeU_b1();
//         double ug4 = mod.computeU_g4();
//         std::cout << "  M_bh x" << factor << " (" << msun << " M_sun): U_b1=" << ub1 
//                   << " J/m³, U_g4=" << ug4 << " J/m³\n";
//     }
//     mod.restoreState("baseline");
//     std::cout << "\n";
//
//     // Step 4: Distance variations (d_g)
//     std::cout << "Step 4: Galactic Distance Scaling\n";
//     std::vector<double> dist_factors = {0.5, 0.8, 1.0, 1.2, 1.5};
//     for (double factor : dist_factors) {
//         mod.updateVariable("d_g", 2.55e20 * factor);
//         double ly = mod.variables["d_g"] / 9.461e15; // m to light years
//         double mbh_dg_new = mod.computeMbhOverDg();
//         double ub1 = mod.computeU_b1();
//         std::cout << "  d_g x" << factor << " (" << ly << " ly): M_bh/d_g=" << mbh_dg_new 
//                   << " kg/m, U_b1=" << ub1 << " J/m³\n";
//     }
//     mod.restoreState("baseline");
//     std::cout << "\n";
//
//     // Step 5: Feedback factor variations
//     std::cout << "Step 5: AGN Feedback Factor Scaling\n";
//     std::vector<double> feedback_vals = {0.0, 0.05, 0.1, 0.2, 0.5, 1.0};
//     for (double f : feedback_vals) {
//         mod.updateVariable("f_feedback", f);
//         double ug4 = mod.computeU_g4();
//         std::cout << "  f_feedback=" << f << ": U_g4=" << ug4 << " J/m³ ";
//         std::cout << "(" << std::showpos << ((ug4/ug4_base - 1.0)*100) << std::noshowpos << "%)\n";
//     }
//     mod.restoreState("baseline");
//     std::cout << "\n";
//
//     // Step 6: Expand black hole scale (M_bh and feedback)
//     std::cout << "Step 6: Expand Black Hole Scale (M_bh x1.5, f_feedback x1.2)\n";
//     mod.expandBlackHoleScale(1.5, 1.2);
//     std::cout << "  New M_bh: " << mod.variables["M_bh"] << " kg (" << mod.computeM_bhInMsun() << " M_sun)\n";
//     std::cout << "  New f_feedback: " << mod.variables["f_feedback"] << "\n";
//     std::cout << "  New U_b1: " << mod.computeU_b1() << " J/m³\n";
//     std::cout << "  New U_g4: " << mod.computeU_g4() << " J/m³\n\n";
//
//     // Step 7: Expand galactic scale (distance and rotation)
//     std::cout << "Step 7: Expand Galactic Scale (d_g x1.1, Ω_g x1.15)\n";
//     mod.expandGalacticScale(1.1, 1.15);
//     std::cout << "  New d_g: " << mod.variables["d_g"] << " m (" << (mod.variables["d_g"]/9.461e15) << " ly)\n";
//     std::cout << "  New Ω_g: " << mod.variables["Omega_g"] << " rad/s\n";
//     std::cout << "  New M_bh/d_g: " << mod.computeMbhOverDg() << " kg/m\n";
//     std::cout << "  New U_b1: " << mod.computeU_b1() << " J/m³\n\n";
//
//     // Step 8: Expand buoyancy scale
//     std::cout << "Step 8: Expand Buoyancy Scale (β_1 x1.2, U_g1 x1.1)\n";
//     mod.expandBuoyancyScale(1.2, 1.1);
//     std::cout << "  New β_1: " << mod.variables["beta_1"] << "\n";
//     std::cout << "  New U_g1: " << mod.variables["U_g1"] << " J/m³\n";
//     std::cout << "  New U_b1: " << mod.computeU_b1() << " J/m³\n\n";
//
//     // Step 9: Parameter variations
//     std::cout << "Step 9: Generate 10 Parameter Variations (±8%)\n";
//     auto variations = mod.generateVariations(10, 0.08);
//     std::cout << "  Generated " << variations.size() << " configurations\n";
//     std::vector<double> msun_range, ub1_range;
//     for (const auto& var : variations) {
//         double msun_val = var.at("M_bh") / var.at("M_sun");
//         msun_range.push_back(msun_val);
//     }
//     auto msun_minmax = std::minmax_element(msun_range.begin(), msun_range.end());
//     std::cout << "  M_bh range: " << *msun_minmax.first << " to " << *msun_minmax.second << " M_sun\n\n";
//
//     // Step 10: Sensitivity analysis
//     std::cout << "Step 10: Sensitivity Analysis (U_b1 response to ±1% changes)\n";
//     std::vector<std::string> sens_params = {"M_bh", "d_g", "Omega_g", "beta_1", "U_g1", "f_feedback"};
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
//     // Step 11: Auto-refine to target U_b1
//     std::cout << "Step 11: Auto-Refine to Target |U_b1| = 2.5e27 J/m³\n";
//     double target_ub1 = -2.5e27;
//     double before_ub1 = mod.computeU_b1();
//     mod.autoRefineParameters("U_b1", target_ub1);
//     double after_ub1 = mod.computeU_b1();
//     std::cout << "  Before: " << before_ub1 << " J/m³\n";
//     std::cout << "  After: " << after_ub1 << " J/m³\n";
//     std::cout << "  Target: " << target_ub1 << " J/m³\n";
//     std::cout << "  New M_bh: " << mod.computeM_bhInMsun() << " M_sun\n\n";
//
//     // Step 12: Target specific M_bh mass
//     std::cout << "Step 12: Target Specific M_bh = 5.0e6 M_sun\n";
//     double target_msun = 5.0e6;
//     double before_msun = mod.computeM_bhInMsun();
//     mod.autoRefineParameters("M_bh_Msun", target_msun);
//     double after_msun = mod.computeM_bhInMsun();
//     std::cout << "  Before: " << before_msun << " M_sun\n";
//     std::cout << "  After: " << after_msun << " M_sun\n";
//     std::cout << "  New M_bh: " << mod.variables["M_bh"] << " kg\n\n";
//
//     // Step 13: Calibration to observations
//     std::cout << "Step 13: Calibrate to Observational Data\n";
//     std::map<std::string, double> observations = {
//         {"M_bh", 8.5e36},       // Slightly higher mass estimate
//         {"d_g", 2.6e20},        // ~27,500 ly
//         {"f_feedback", 0.12}    // Enhanced feedback
//     };
//     mod.calibrateToObservations(observations);
//     std::cout << "  Calibrated M_bh: " << mod.computeM_bhInMsun() << " M_sun\n";
//     std::cout << "  Calibrated d_g: " << (mod.variables["d_g"]/9.461e15) << " ly\n";
//     std::cout << "  Calibrated f_feedback: " << mod.variables["f_feedback"] << "\n";
//     std::cout << "  New U_b1: " << mod.computeU_b1() << " J/m³\n\n";
//
//     // Step 14: Optimize for maximum buoyancy
//     std::cout << "Step 14: Optimize for Maximum |U_b1| (Buoyancy)\n";
//     double before_opt = std::abs(mod.computeU_b1());
//     mod.optimizeForMetric("maximize_U_b1");
//     double after_opt = std::abs(mod.computeU_b1());
//     std::cout << "  |U_b1| before: " << before_opt << " J/m³\n";
//     std::cout << "  |U_b1| after: " << after_opt << " J/m³\n";
//     std::cout << "  Increase: " << ((after_opt/before_opt - 1.0)*100) << "%\n\n";
//
//     // Step 15: Optimize for enhanced feedback
//     std::cout << "Step 15: Optimize for Enhanced AGN Feedback\n";
//     mod.saveState("before_feedback_opt");
//     double f_before = mod.variables["f_feedback"];
//     mod.optimizeForMetric("enhance_feedback");
//     double f_after = mod.variables["f_feedback"];
//     std::cout << "  f_feedback: " << f_before << " -> " << f_after << "\n";
//     std::cout << "  New U_g4: " << mod.computeU_g4() << " J/m³\n\n";
//
//     // Step 16: Parameter mutation
//     std::cout << "Step 16: Mutate Parameters (±5%)\n";
//     mod.restoreState("before_feedback_opt");
//     mod.saveState("before_mutation");
//     double mbh_before = mod.variables["M_bh"];
//     mod.mutateParameters(0.05);
//     double mbh_after = mod.variables["M_bh"];
//     std::cout << "  M_bh: " << (mbh_before/mod.variables["M_sun"]) << " -> " 
//               << (mbh_after/mod.variables["M_sun"]) << " M_sun\n";
//     std::cout << "  Change: " << ((mbh_after/mbh_before - 1.0)*100) << "%\n\n";
//
//     // Step 17: System evolution (maximize |U_b1|)
//     std::cout << "Step 17: Evolve System (8 generations, maximize |U_b1|)\n";
//     mod.restoreState("before_mutation");
//     double initial_fitness = std::abs(mod.computeU_b1());
//     mod.evolveSystem(8, [&mod]() { return std::abs(mod.computeU_b1()); });
//     double final_fitness = std::abs(mod.computeU_b1());
//     std::cout << "  Initial |U_b1|: " << initial_fitness << " J/m³\n";
//     std::cout << "  Evolved |U_b1|: " << final_fitness << " J/m³\n";
//     std::cout << "  Improvement: " << ((final_fitness / initial_fitness - 1.0) * 100) << "%\n";
//     std::cout << "  Final M_bh: " << mod.computeM_bhInMsun() << " M_sun\n\n";
//
//     // Step 18: Validation and state export
//     std::cout << "Step 18: Validate Consistency and Export\n";
//     bool valid = mod.validateConsistency();
//     std::cout << "  System valid: " << (valid ? "YES" : "NO") << "\n";
//     if (!valid) {
//         std::cout << "  Running auto-correction...\n";
//         mod.autoCorrectAnomalies();
//         std::cout << "  Post-correction valid: " << (mod.validateConsistency() ? "YES" : "NO") << "\n";
//     }
//     std::cout << "  Final M_bh: " << mod.computeM_bhInMsun() << " M_sun\n";
//     std::cout << "  Final U_b1: " << mod.computeU_b1() << " J/m³\n";
//     std::cout << "  Final U_g4: " << mod.computeU_g4() << " J/m³\n\n";
//     mod.saveState("evolved_optimal");
//     auto saved = mod.listSavedStates();
//     std::cout << "  Saved states (" << saved.size() << "): ";
//     for (const auto& s : saved) std::cout << s << " ";
//     std::cout << "\n\n";
//     std::cout << "Final System Export:\n";
//     std::string exported = mod.exportState();
//     std::cout << exported << "\n";
//     std::cout << "Demo complete: 18 steps executed successfully!\n";
//     std::cout << "================================================================\n";
//
//     return 0;
// }

GalacticBlackHoleModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeM_bh, computeM_bhInMsun, computeMbhOverDg, computeU_b1, computeU_g4) are clear, concise, and variable - driven.
- Output and debugging functions(printVariables, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.
- Approximations and physical context are clearly stated, aiding scientific understanding.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid indices, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in galactic black hole mass modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.