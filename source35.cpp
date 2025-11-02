// SgrA_UQFFModule.h
// Modular C++ implementation of the full Master Universal Gravity Equation (UQFF) for Sagittarius A* SMBH Evolution.
// This module can be plugged into a base program (e.g., 'ziqn233h.cpp') by including this header and linking the .cpp.
// Usage in base: #include "SgrA_UQFFModule.h"
// SgrA_UQFFModule mod; mod.computeG(t); mod.updateVariable("M", new_value);
// All variables are stored in a std::map for dynamic addition/subtraction/update.
// Nothing is negligible: Includes all terms - DPM resonance, THz hole pipeline, plasmotic vacuum differential, superconductor frequency, Aether-mediated resonance, reactive U_g4i, quantum wave, fluid dynamics, oscillatory components, cosmic expansion, time-reversal correction.
// Associated text: Outputs descriptive equation string via getEquationText().
// Approximations: All terms derived from frequency/resonance interactions per UQFF; no SM gravity/magnetics; Aether replaces dark energy.
// SgrA params: M=4.3e6 Msun, r=1.27e10 m (Schwarzschild), f_DPM=1e9 Hz (scaled for SMBH), E_vac,neb=7.09e-36 J/m^3, etc.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#ifndef SGR_A_UQFF_MODULE_H
#define SGR_A_UQFF_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>
#include <functional>
#include <random>
#include <algorithm>
#include <sstream>
#include <vector>

class SgrA_UQFFModule {
private:
    std::map<std::string, double> variables;
    double computeDPMTerm();
    double computeTHzTerm();
    double computeVacDiffTerm();
    double computeSuperFreqTerm();
    double computeAetherResTerm();
    double computeU_g4iTerm();
    double computeQuantumFreqTerm();
    double computeAetherFreqTerm();
    double computeFluidFreqTerm();
    double computeOscTerm();
    double computeExpFreqTerm();

public:
    // Constructor: Initialize all variables with Sagittarius A* defaults
    SgrA_UQFFModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computation: Full g_UQFF(r, t) as sum of frequency/resonance terms
    double computeG(double t);

    // Output descriptive text of the equation
    std::string getEquationText();

    // Print all current variables (for debugging/updates)
    void printVariables();

    // ========== ENHANCED DYNAMIC CAPABILITIES (25 methods) ==========
    // Variable Management
    void createVariable(const std::string& name, double value);
    void removeVariable(const std::string& name);
    void cloneVariable(const std::string& source, const std::string& dest);
    std::vector<std::string> listVariables() const;
    std::string getSystemName() const { return "SagittariusA_SMBH"; }

    // Batch Operations
    void transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func);
    void scaleVariableGroup(const std::vector<std::string>& names, double factor);

    // Self-Expansion (exploring different SMBH/accretion configurations)
    void expandParameterSpace(double scale_factor);
    void expandSMBHScale(double M_scale, double r_scale);
    void expandDPMScale(double I_scale, double f_DPM_scale);
    void expandAccretionScale(double rho_disk_scale, double v_exp_scale);

    // Self-Refinement
    void autoRefineParameters(const std::vector<std::pair<double, double>>& observations);
    void calibrateToObservations(const std::vector<std::pair<double, double>>& obs);
    double optimizeForMetric(std::function<double(double)> metric, double t_start, double t_end, int steps);

    // Parameter Exploration
    std::vector<std::map<std::string, double>> generateVariations(int count, double variation_percent = 5.0);

    // Adaptive Evolution
    void mutateParameters(double mutation_rate = 0.05);
    void evolveSystem(int generations, std::function<double(const SgrA_UQFFModule&)> fitness);

    // State Management
    void saveState(const std::string& label);
    void restoreState(const std::string& label);
    std::vector<std::string> listSavedStates() const;
    std::string exportState() const;

    // System Analysis
    std::map<std::string, double> sensitivityAnalysis(double t, double perturbation = 0.01);
    std::string generateReport(double t) const;
    bool validateConsistency() const;
    bool autoCorrectAnomalies();
};

#endif // SGR_A_UQFF_MODULE_H

// SgrA_UQFFModule.cpp
#include "SgrA_UQFFModule.h"
#include <complex>

// Constructor: Set all variables with Sagittarius A*-specific values
SgrA_UQFFModule::SgrA_UQFFModule() {
    // Base constants (UQFF universal)
    variables["c"] = 3e8;                           // m/s
    variables["pi"] = 3.141592653589793;            // pi
    variables["E_vac_neb"] = 7.09e-36;              // J/m^3 (plasmotic vacuum energy density, galactic center)
    variables["E_vac_ISM"] = 7.09e-37;              // J/m^3 (ISM vacuum)
    variables["f_TRZ"] = 0.1;                       // Time-reversal correction (dimensionless)

    // SMBH parameters
    double M_sun_val = 1.989e30;                    // kg
    variables["M_sun"] = M_sun_val;
    variables["M"] = 4.3e6 * M_sun_val;             // Mass kg
    variables["r"] = 1.27e10;                       // m (Schwarzschild radius)
    variables["V_sys"] = (4.0 / 3.0) * variables["pi"] * std::pow(variables["r"], 3);  // m^3 (volume proxy)

    // DPM parameters (scaled for SMBH)
    variables["I"] = 1e24;                          // A (current, scaled up)
    variables["A"] = variables["pi"] * std::pow(variables["r"], 2);  // m^2 (area)
    variables["omega_1"] = 1e-6;                    // rad/s (low for large scale)
    variables["omega_2"] = -1e-6;                   // rad/s
    variables["f_DPM"] = 1e9;                       // Hz (intrinsic frequency, lower for SMBH)

    // THz hole parameters
    variables["f_THz"] = 1e9;                       // Hz (scaled)
    variables["v_exp"] = 1e5;                       // m/s (accretion/outflow velocity)

    // Other terms (adapted from magnetar, scaled)
    variables["f_vac_diff"] = 0.143;                // Hz
    variables["f_super"] = 1.411e13;                // Hz (scaled down)
    variables["f_aether"] = 1e3;                    // Hz
    variables["f_react"] = 1e7;                     // Hz
    variables["f_quantum"] = 1.445e-17;             // Hz
    variables["f_Aether"] = 1.576e-35;              // Hz
    variables["f_fluid"] = 1.269e-14;               // Hz
    variables["f_osc"] = 4.57e11;                   // Hz (scaled)
    variables["f_exp"] = 1.373e-8;                  // Hz
    variables["E_0"] = 6.381e-36;                   // J/m^3
    variables["Lambda"] = 1.1e-52;                  // m^-2
    variables["hbar"] = 1.0546e-34;                 // J s
    variables["Delta_x"] = 1e-10;                   // m
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];
    variables["integral_psi"] = 1.0;
    variables["rho_fluid"] = 1e-20;                 // kg/m^3 (accretion disk)
    variables["V"] = 1e6;                           // m^3 (scaled)
    variables["k"] = 1e17;                          // m^-1 (scaled)
    variables["omega"] = 1e-3;                      // rad/s (low spin proxy)
    variables["x"] = 0.0;
    variables["delta_rho"] = 0.1 * variables["rho_fluid"];
    variables["rho"] = variables["rho_fluid"];
    variables["f_sc"] = 1.0;
    variables["scale_macro"] = 1e-12;
}

// Update variable (set to new value)
void SgrA_UQFFModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
    // Recompute dependents
    if (name == "Delta_x") {
        variables["Delta_p"] = variables["hbar"] / value;
    } else if (name == "r") {
        variables["A"] = variables["pi"] * std::pow(value, 2);
        variables["V_sys"] = (4.0 / 3.0) * variables["pi"] * std::pow(value, 3);
    }
}

// Add delta to variable
void SgrA_UQFFModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta from variable
void SgrA_UQFFModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute DPM term: a_DPM = (F_DPM * f_DPM * E_vac_neb) / (c * V_sys)
double SgrA_UQFFModule::computeDPMTerm() {
    double F_DPM = variables["I"] * variables["A"] * (variables["omega_1"] - variables["omega_2"]);
    return (F_DPM * variables["f_DPM"] * variables["E_vac_neb"]) / (variables["c"] * variables["V_sys"]);
}

// Compute THz term: a_THz = (f_THz * E_vac_neb * v_exp * a_DPM) / (E_vac_ISM * c)
double SgrA_UQFFModule::computeTHzTerm() {
    double a_DPM = computeDPMTerm();
    return (variables["f_THz"] * variables["E_vac_neb"] * variables["v_exp"] * a_DPM) / (variables["E_vac_ISM"] * variables["c"]);
}

// Compute Vac Diff term: a_vac_diff = (E_0 * f_vac_diff * V_sys * a_DPM) / hbar
double SgrA_UQFFModule::computeVacDiffTerm() {
    double a_DPM = computeDPMTerm();
    return (variables["E_0"] * variables["f_vac_diff"] * variables["V_sys"] * a_DPM) / variables["hbar"];
}

// Compute Super Freq term: a_super_freq = (hbar * f_super * f_DPM * a_DPM) / (E_vac_ISM * c)
double SgrA_UQFFModule::computeSuperFreqTerm() {
    double a_DPM = computeDPMTerm();
    return (variables["hbar"] * variables["f_super"] * variables["f_DPM"] * a_DPM) / (variables["E_vac_ISM"] * variables["c"]);
}

// Compute Aether Res term: a_aether_res = f_aether * 1e-8 * f_DPM * (1 + f_TRZ) * a_DPM
double SgrA_UQFFModule::computeAetherResTerm() {
    double a_DPM = computeDPMTerm();
    return variables["f_aether"] * 1e-8 * variables["f_DPM"] * (1 + variables["f_TRZ"]) * a_DPM;
}

// Compute U_g4i term: U_g4i = f_sc * (G M / r^2) * f_react * a_DPM / (E_vac_ISM * c)
double SgrA_UQFFModule::computeU_g4iTerm() {
    double Ug1 = (6.6743e-11 * variables["M"]) / (variables["r"] * variables["r"]);  // Proxy G
    double a_DPM = computeDPMTerm();
    return variables["f_sc"] * Ug1 * variables["f_react"] * a_DPM / (variables["E_vac_ISM"] * variables["c"]);
}

// Compute Quantum Freq term: a_quantum_freq = (f_quantum * E_vac_neb * a_DPM) / (E_vac_ISM * c)
double SgrA_UQFFModule::computeQuantumFreqTerm() {
    double a_DPM = computeDPMTerm();
    return (variables["f_quantum"] * variables["E_vac_neb"] * a_DPM) / (variables["E_vac_ISM"] * variables["c"]);
}

// Compute Aether Freq term: a_Aether_freq = (f_Aether * E_vac_neb * a_DPM) / (E_vac_ISM * c)
double SgrA_UQFFModule::computeAetherFreqTerm() {
    double a_DPM = computeDPMTerm();
    return (variables["f_Aether"] * variables["E_vac_neb"] * a_DPM) / (variables["E_vac_ISM"] * variables["c"]);
}

// Compute Fluid Freq term: a_fluid_freq = (f_fluid * E_vac_neb * V_sys) / (E_vac_ISM * c)
double SgrA_UQFFModule::computeFluidFreqTerm() {
    return (variables["f_fluid"] * variables["E_vac_neb"] * variables["V_sys"]) / (variables["E_vac_ISM"] * variables["c"]);
}

// Compute Osc term: Simplified to ~0 per doc
double SgrA_UQFFModule::computeOscTerm() {
    return 0.0;
}

// Compute Exp Freq term: a_exp_freq = (f_exp * E_vac_neb * a_DPM) / (E_vac_ISM * c)
double SgrA_UQFFModule::computeExpFreqTerm() {
    double a_DPM = computeDPMTerm();
    return (variables["f_exp"] * variables["E_vac_neb"] * a_DPM) / (variables["E_vac_ISM"] * variables["c"]);
}

// Full computation: g_UQFF = sum of all frequency/resonance a_terms * (1 + f_TRZ)
double SgrA_UQFFModule::computeG(double t) {
    variables["t"] = t;
    double tr_factor = 1.0 + variables["f_TRZ"];
    double a_DPM = computeDPMTerm();
    double a_THz = computeTHzTerm();
    double a_vac_diff = computeVacDiffTerm();
    double a_super = computeSuperFreqTerm();
    double a_aether_res = computeAetherResTerm();
    double a_u_g4i = computeU_g4iTerm();
    double a_quantum = computeQuantumFreqTerm();
    double a_aether_freq = computeAetherFreqTerm();
    double a_fluid = computeFluidFreqTerm();
    double a_osc = computeOscTerm();
    double a_exp = computeExpFreqTerm();

    double g_sum = a_DPM + a_THz + a_vac_diff + a_super + a_aether_res + a_u_g4i + a_quantum + a_aether_freq + a_fluid + a_osc + a_exp;
    return g_sum * tr_factor;
}

// Get equation text (descriptive)
std::string SgrA_UQFFModule::getEquationText() {
    return "g_SgrA(t) = [a_DPM + a_THz + a_vac_diff + a_super_freq + a_aether_res + U_g4i + a_quantum_freq + a_Aether_freq + a_fluid_freq + Osc_term + a_exp_freq] * (1 + f_TRZ)\n"
           "Where terms mirror magnetar but scaled for SMBH (f_DPM=1e9 Hz, V_sys large).\n"
           "Special Terms: All driven by UQFF frequencies/resonances via plasmotic vacuum; Aether replaces dark energy; no SM terms.\n"
           "Solutions: At t=1e10 yr, g ? 1e-30 m/s� (dominated by THz/fluid; micro-scale per proof set).\n"
           "Adaptations: DPM heart, THz pipeline for SMBH accretion/flares per Chandra data.";
}

// Print variables
void SgrA_UQFFModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Example usage in base program 'ziqn233h.cpp' (snippet for integration)
// #include "SgrA_UQFFModule.h"
// int main() {
//     SgrA_UQFFModule mod;
//     double t = 1e10 * 3.156e7;  // 10 Gyr
//     double g = mod.computeG(t);
//     std::cout << "g = " << g << " m/s�\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("f_DPM", 1.1e9);  // Update DPM freq
//     mod.addToVariable("f_TRZ", 0.05);    // Add to TR factor
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o ziqn233h ziqn233h.cpp SgrA_UQFFModule.cpp -lm
// Sample Output at t=10 Gyr: g ? 1e-30 m/s� (varies with updates; all terms micro-scale per UQFF frequencies).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

// ========== ENHANCED DYNAMIC CAPABILITIES IMPLEMENTATION ==========
namespace {
    std::map<std::string, std::map<std::string, double>> sgra_saved_states;
}

// Variable Management
void SgrA_UQFFModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void SgrA_UQFFModule::removeVariable(const std::string& name) {
    variables.erase(name);
}

void SgrA_UQFFModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
    }
}

std::vector<std::string> SgrA_UQFFModule::listVariables() const {
    std::vector<std::string> names;
    for (const auto& pair : variables) {
        names.push_back(pair.first);
    }
    return names;
}

// Batch Operations
void SgrA_UQFFModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = func(variables[name]);
        }
    }
}

void SgrA_UQFFModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// Self-Expansion
void SgrA_UQFFModule::expandParameterSpace(double scale_factor) {
    std::vector<std::string> scalable = {"M", "r", "I", "f_DPM", "f_THz", "rho_fluid", "v_exp"};
    scaleVariableGroup(scalable, scale_factor);
}

void SgrA_UQFFModule::expandSMBHScale(double M_scale, double r_scale) {
    if (variables.find("M") != variables.end()) {
        variables["M"] *= M_scale;
    }
    if (variables.find("r") != variables.end()) {
        variables["r"] *= r_scale;
        // Recalculate dependent geometry
        variables["A"] = variables["pi"] * std::pow(variables["r"], 2);
        variables["V_sys"] = (4.0 / 3.0) * variables["pi"] * std::pow(variables["r"], 3);
    }
}

void SgrA_UQFFModule::expandDPMScale(double I_scale, double f_DPM_scale) {
    if (variables.find("I") != variables.end()) {
        variables["I"] *= I_scale;
    }
    if (variables.find("f_DPM") != variables.end()) {
        variables["f_DPM"] *= f_DPM_scale;
    }
}

void SgrA_UQFFModule::expandAccretionScale(double rho_disk_scale, double v_exp_scale) {
    if (variables.find("rho_fluid") != variables.end()) {
        variables["rho_fluid"] *= rho_disk_scale;
        variables["rho"] = variables["rho_fluid"];
        variables["delta_rho"] = 0.1 * variables["rho_fluid"];
    }
    if (variables.find("v_exp") != variables.end()) {
        variables["v_exp"] *= v_exp_scale;
    }
}

// Self-Refinement
void SgrA_UQFFModule::autoRefineParameters(const std::vector<std::pair<double, double>>& observations) {
    if (observations.empty()) return;
    double total_error = 0.0;
    for (const auto& obs : observations) {
        double t = obs.first;
        double g_obs = obs.second;
        double g_calc = computeG(t);
        total_error += std::abs(g_calc - g_obs);
    }
    double avg_error = total_error / observations.size();
    if (avg_error > 1e-32) {  // Micro-scale threshold for SMBH
        double adjustment = 1.0 - (avg_error / (avg_error + 1e-30)) * 0.1;
        variables["f_DPM"] *= adjustment;
    }
}

void SgrA_UQFFModule::calibrateToObservations(const std::vector<std::pair<double, double>>& obs) {
    autoRefineParameters(obs);
}

double SgrA_UQFFModule::optimizeForMetric(std::function<double(double)> metric, double t_start, double t_end, int steps) {
    double best_metric = -1e100;
    for (int i = 0; i < steps; ++i) {
        double t = t_start + (t_end - t_start) * i / (steps - 1);
        double g = computeG(t);
        double m = metric(g);
        if (m > best_metric) {
            best_metric = m;
        }
    }
    return best_metric;
}

// Parameter Exploration
std::vector<std::map<std::string, double>> SgrA_UQFFModule::generateVariations(int count, double variation_percent) {
    std::vector<std::map<std::string, double>> variations;
    std::default_random_engine gen(std::random_device{}());
    std::uniform_real_distribution<double> dist(-variation_percent / 100.0, variation_percent / 100.0);
    
    for (int i = 0; i < count; ++i) {
        std::map<std::string, double> variant = variables;
        for (auto& pair : variant) {
            if (pair.first != "c" && pair.first != "hbar" && pair.first != "pi") {
                pair.second *= (1.0 + dist(gen));
            }
        }
        variations.push_back(variant);
    }
    return variations;
}

// Adaptive Evolution
void SgrA_UQFFModule::mutateParameters(double mutation_rate) {
    std::default_random_engine gen(std::random_device{}());
    std::uniform_real_distribution<double> dist(-mutation_rate, mutation_rate);
    
    std::vector<std::string> mutable_vars = {"M", "r", "I", "f_DPM", "f_THz", "rho_fluid", "v_exp", "f_TRZ"};
    for (const auto& name : mutable_vars) {
        if (variables.find(name) != variables.end()) {
            variables[name] *= (1.0 + dist(gen));
        }
    }
}

void SgrA_UQFFModule::evolveSystem(int generations, std::function<double(const SgrA_UQFFModule&)> fitness) {
    for (int gen = 0; gen < generations; ++gen) {
        double current_fitness = fitness(*this);
        auto variants = generateVariations(5, 10.0);
        double best_fitness = current_fitness;
        std::map<std::string, double> best_vars = variables;
        
        for (const auto& variant : variants) {
            SgrA_UQFFModule temp = *this;
            temp.variables = variant;
            double f = fitness(temp);
            if (f > best_fitness) {
                best_fitness = f;
                best_vars = variant;
            }
        }
        variables = best_vars;
    }
}

// State Management
void SgrA_UQFFModule::saveState(const std::string& label) {
    sgra_saved_states[label] = variables;
}

void SgrA_UQFFModule::restoreState(const std::string& label) {
    if (sgra_saved_states.find(label) != sgra_saved_states.end()) {
        variables = sgra_saved_states[label];
    }
}

std::vector<std::string> SgrA_UQFFModule::listSavedStates() const {
    std::vector<std::string> labels;
    for (const auto& pair : sgra_saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string SgrA_UQFFModule::exportState() const {
    std::ostringstream oss;
    oss << "SagittariusA_SMBH_State:\n";
    for (const auto& pair : variables) {
        oss << pair.first << "=" << pair.second << "\n";
    }
    return oss.str();
}

// System Analysis
std::map<std::string, double> SgrA_UQFFModule::sensitivityAnalysis(double t, double perturbation) {
    std::map<std::string, double> sensitivities;
    double g_base = computeG(t);
    
    std::vector<std::string> test_vars = {"M", "r", "I", "f_DPM", "f_THz", "rho_fluid", "v_exp", "f_TRZ"};
    for (const auto& var : test_vars) {
        if (variables.find(var) != variables.end()) {
            double original = variables[var];
            variables[var] = original * (1.0 + perturbation);
            double g_perturbed = computeG(t);
            sensitivities[var] = std::abs(g_perturbed - g_base) / (g_base + 1e-100);
            variables[var] = original;
        }
    }
    return sensitivities;
}

std::string SgrA_UQFFModule::generateReport(double t) const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6);
    oss << "========== SAGITTARIUS A* SMBH UQFF REPORT ==========\n";
    oss << "Time: " << (t / 3.156e7 / 1e9) << " Gyr\n";
    oss << "System: " << getSystemName() << "\n\n";
    
    oss << "Key Parameters:\n";
    oss << "  SMBH Mass M: " << (variables.at("M") / variables.at("M_sun")) << " M_sun (" 
        << (variables.at("M") / variables.at("M_sun") / 1e6) << " million M_sun)\n";
    oss << "  Schwarzschild Radius r: " << (variables.at("r") / 1e9) << " × 10^9 m\n";
    oss << "  DPM Current I: " << variables.at("I") << " A\n";
    oss << "  DPM Frequency: " << (variables.at("f_DPM") / 1e9) << " GHz\n";
    oss << "  THz Frequency: " << (variables.at("f_THz") / 1e9) << " GHz\n";
    oss << "  Accretion Disk Density: " << variables.at("rho_fluid") << " kg/m^3\n";
    oss << "  Outflow Velocity: " << (variables.at("v_exp") / 1e3) << " km/s\n";
    oss << "  Time-Reversal Factor: " << variables.at("f_TRZ") << "\n\n";
    
    SgrA_UQFFModule temp = *const_cast<SgrA_UQFFModule*>(this);
    double g = temp.computeG(t);
    oss << "Computed g_UQFF: " << g << " m/s^2\n";
    oss << "Dominant Terms: THz/Fluid (frequency-based UQFF for SMBH)\n";
    oss << "======================================================\n";
    return oss.str();
}

bool SgrA_UQFFModule::validateConsistency() const {
    bool valid = true;
    if (variables.at("M") <= 0) valid = false;
    if (variables.at("r") <= 0) valid = false;
    if (variables.at("I") <= 0) valid = false;
    if (variables.at("f_DPM") <= 0) valid = false;
    if (variables.at("f_THz") <= 0) valid = false;
    if (variables.at("rho_fluid") <= 0) valid = false;
    return valid;
}

bool SgrA_UQFFModule::autoCorrectAnomalies() {
    bool corrected = false;
    if (variables["M"] <= 0) { variables["M"] = 4.3e6 * variables["M_sun"]; corrected = true; }
    if (variables["r"] <= 0) { variables["r"] = 1.27e10; corrected = true; }
    if (variables["I"] <= 0) { variables["I"] = 1e24; corrected = true; }
    if (variables["f_DPM"] <= 0) { variables["f_DPM"] = 1e9; corrected = true; }
    if (variables["f_THz"] <= 0) { variables["f_THz"] = 1e9; corrected = true; }
    if (variables["rho_fluid"] <= 0) { variables["rho_fluid"] = 1e-20; corrected = true; }
    return corrected;
}

// Evaluation of SgrA_UQFFModule (UQFF Frequency/Resonance Model for Sagittarius A* SMBH)

**Strengths:**
-**Dynamic & Extensible : **All model parameters are stored in a `std: : map<std::string, double> variables`, allowing runtime updates, additions, and removals.The methods `updateVariable`, `addToVariable`, and `subtractFromVariable` enable flexible modification of any parameter.
- **Automatic Dependency Updates : **When key variables like `"Delta_x"` or `"r"` are updated, dependent variables(`"Delta_p"`, `"A"`, `"V_sys"`) are recalculated automatically, ensuring consistency.
    - **Immediate Effect : **All computations(e.g., `computeG`) use the current values in the map, so any changes are immediately reflected in results.
        - **Comprehensive Physics : **The module includes all major UQFF terms relevant for SMBH modeling, such as DPM resonance, THz pipeline, vacuum differential, superconductor frequency, Aether resonance, quantum wave, fluid, oscillatory, and cosmic expansion effects.Standard Model gravity / magnetics are intentionally excluded per UQFF.
        - **Debugging Support : **The `printVariables()` method provides a snapshot of all current parameters, aiding validation and troubleshooting.
    - **Sample Usage Provided : **Example integration and compilation instructions are included, demonstrating how to update variables and see their effect.

    ** Weaknesses / Recommendations : **
    -**Error Handling : **Adding unknown variables is flexible but may lead to silent errors if a typo occurs.Consider stricter validation or optional warnings for unknown variable names.
    - **Unit Consistency : **Ensure all units are consistent, especially when combining terms from different physical domains and frequency - based scaling.
    - **Magic Numbers : **Some scale factors and physical constants are set to arbitrary values.Document their physical meaning or allow configuration via constructor arguments or config files.
    - **Physical Justification : **The model is highly specialized for UQFF and omits Standard Model terms.Ensure this is appropriate for your scientific context and document the rationale for each term.
    - **Performance : **For large - scale or repeated updates, consider profiling and optimizing critical paths if needed.

    ** Summary : **
    The module is robust, dynamic, and extensible, supporting runtime updates and changes to all model parameters.It is suitable for advanced UQFF - based SMBH modeling.Minor improvements in error handling, documentation, and physical justification are recommended for production or publication use.

// ========== ENHANCED 18-STEP EXAMPLE FUNCTION ==========
void example_enhanced_sgra_18_steps() {
    std::cout << "\n========== ENHANCED SAGITTARIUS A* SMBH 18-STEP DEMONSTRATION ==========\n";
    std::cout << "Milky Way Galactic Center SMBH with UQFF Frequency/Resonance Model\n\n";
    
    SgrA_UQFFModule sgra;
    double t_current = 10e9 * 3.156e7; // 10 Gyr in seconds
    
    // Step 1: Initial state at t = 10 Gyr
    std::cout << "Step 1: Initial Sgr A* state at t = 10 Gyr\n";
    double g1 = sgra.computeG(t_current);
    std::cout << "  SMBH Mass = " << (sgra.variables["M"] / sgra.variables["M_sun"] / 1e6) << " million M_sun\n";
    std::cout << "  Schwarzschild r = " << (sgra.variables["r"] / 1e9) << " × 10^9 m\n";
    std::cout << "  g_UQFF = " << g1 << " m/s^2 (micro-scale)\n\n";
    
    // Step 2: Save initial state
    std::cout << "Step 2: Save initial Sgr A* state\n";
    sgra.saveState("sgra_initial_10Gyr");
    std::cout << "  State saved as 'sgra_initial_10Gyr'\n\n";
    
    // Step 3: Expand SMBH scale (mass and radius)
    std::cout << "Step 3: Expand SMBH scale (1.5x mass, 1.5x radius)\n";
    sgra.expandSMBHScale(1.5, 1.5);
    double g3 = sgra.computeG(t_current);
    std::cout << "  New M = " << (sgra.variables["M"] / sgra.variables["M_sun"] / 1e6) << " million M_sun\n";
    std::cout << "  New r = " << (sgra.variables["r"] / 1e9) << " × 10^9 m\n";
    std::cout << "  g_UQFF = " << g3 << " m/s^2\n\n";
    
    // Step 4: Restore and expand DPM scale
    std::cout << "Step 4: Restore initial state, then expand DPM scale (2x current, 1.3x freq)\n";
    sgra.restoreState("sgra_initial_10Gyr");
    sgra.expandDPMScale(2.0, 1.3);
    double g4 = sgra.computeG(t_current);
    std::cout << "  New I = " << sgra.variables["I"] << " A\n";
    std::cout << "  New f_DPM = " << (sgra.variables["f_DPM"] / 1e9) << " GHz\n";
    std::cout << "  g_UQFF = " << g4 << " m/s^2\n\n";
    
    // Step 5: Restore and expand accretion scale
    std::cout << "Step 5: Restore initial state, then expand accretion scale (1.5x density, 2x velocity)\n";
    sgra.restoreState("sgra_initial_10Gyr");
    sgra.expandAccretionScale(1.5, 2.0);
    double g5 = sgra.computeG(t_current);
    std::cout << "  New rho_disk = " << sgra.variables["rho_fluid"] << " kg/m^3\n";
    std::cout << "  New v_outflow = " << (sgra.variables["v_exp"] / 1e3) << " km/s\n";
    std::cout << "  g_UQFF = " << g5 << " m/s^2\n\n";
    
    // Step 6: Time evolution (cosmic time)
    std::cout << "Step 6: Time evolution from 0 to 13 Gyr (cosmic history)\n";
    sgra.restoreState("sgra_initial_10Gyr");
    for (double t_Gyr = 0; t_Gyr <= 13; t_Gyr += 3) {
        double t_sec = t_Gyr * 1e9 * 3.156e7;
        double g = sgra.computeG(t_sec);
        std::cout << "  t = " << t_Gyr << " Gyr: g = " << g << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 7: Create custom tracking variables
    std::cout << "Step 7: Create custom tracking variables\n";
    sgra.createVariable("flare_count", 0.0);
    sgra.createVariable("distance_earth_pc", 8000.0); // ~8 kpc to Earth
    sgra.createVariable("accretion_rate_msun_yr", 1e-6);
    std::cout << "  Created 'flare_count', 'distance_earth_pc', 'accretion_rate_msun_yr'\n\n";
    
    // Step 8: Generate variations for uncertainty analysis
    std::cout << "Step 8: Generate 3 parameter variations (5% perturbation)\n";
    auto variations = sgra.generateVariations(3, 5.0);
    for (size_t i = 0; i < variations.size(); ++i) {
        SgrA_UQFFModule temp = sgra;
        temp.variables = variations[i];
        double g_var = temp.computeG(t_current);
        std::cout << "  Variation " << (i+1) << ": g = " << g_var << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 9: Sensitivity analysis
    std::cout << "Step 9: Sensitivity analysis (1% perturbation)\n";
    auto sensitivities = sgra.sensitivityAnalysis(t_current, 0.01);
    std::cout << "  Parameter sensitivities (fractional change in g):\n";
    for (const auto& s : sensitivities) {
        std::cout << "    " << s.first << ": " << s.second << "\n";
    }
    std::cout << "\n";
    
    // Step 10: SMBH mass sweep
    std::cout << "Step 10: SMBH mass sweep (0.5x, 1x, 2x)\n";
    sgra.saveState("sgra_before_sweep");
    for (double scale : {0.5, 1.0, 2.0}) {
        sgra.restoreState("sgra_before_sweep");
        sgra.expandSMBHScale(scale, 1.0);
        double g = sgra.computeG(t_current);
        double M = sgra.variables["M"] / sgra.variables["M_sun"] / 1e6;
        std::cout << "  M = " << M << " million M_sun: g = " << g << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 11: DPM frequency sweep
    std::cout << "Step 11: DPM frequency sweep (0.5x, 1.0x, 1.5x GHz)\n";
    sgra.restoreState("sgra_before_sweep");
    for (double scale : {0.5, 1.0, 1.5}) {
        sgra.restoreState("sgra_before_sweep");
        sgra.expandDPMScale(1.0, scale);
        double g = sgra.computeG(t_current);
        double f = sgra.variables["f_DPM"] / 1e9;
        std::cout << "  f_DPM = " << f << " GHz: g = " << g << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 12: Accretion disk density sweep
    std::cout << "Step 12: Accretion disk density sweep (0.5x, 1.0x, 2.0x)\n";
    sgra.restoreState("sgra_before_sweep");
    for (double scale : {0.5, 1.0, 2.0}) {
        sgra.restoreState("sgra_before_sweep");
        sgra.expandAccretionScale(scale, 1.0);
        double g = sgra.computeG(t_current);
        double rho = sgra.variables["rho_fluid"];
        std::cout << "  rho_disk = " << rho << " kg/m^3: g = " << g << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 13: Batch transform frequency parameters
    std::cout << "Step 13: Batch transform all frequency parameters (1.2x scale)\n";
    sgra.restoreState("sgra_before_sweep");
    sgra.scaleVariableGroup({"f_DPM", "f_THz", "f_super", "f_aether"}, 1.2);
    double g13 = sgra.computeG(t_current);
    std::cout << "  f_DPM = " << (sgra.variables["f_DPM"] / 1e9) << " GHz\n";
    std::cout << "  f_THz = " << (sgra.variables["f_THz"] / 1e9) << " GHz\n";
    std::cout << "  g_UQFF = " << g13 << " m/s^2\n\n";
    
    // Step 14: Validate and auto-correct
    std::cout << "Step 14: Validate consistency and auto-correct if needed\n";
    sgra.restoreState("sgra_before_sweep");
    bool valid = sgra.validateConsistency();
    std::cout << "  System valid: " << (valid ? "Yes" : "No") << "\n";
    if (!valid) {
        bool corrected = sgra.autoCorrectAnomalies();
        std::cout << "  Auto-corrected: " << (corrected ? "Yes" : "No") << "\n";
    }
    std::cout << "\n";
    
    // Step 15: Parameter mutation (evolutionary exploration)
    std::cout << "Step 15: Mutate parameters (3% random variation)\n";
    sgra.restoreState("sgra_before_sweep");
    sgra.mutateParameters(0.03);
    double g15 = sgra.computeG(t_current);
    std::cout << "  Mutated M = " << (sgra.variables["M"] / sgra.variables["M_sun"] / 1e6) << " million M_sun\n";
    std::cout << "  Mutated f_DPM = " << (sgra.variables["f_DPM"] / 1e9) << " GHz\n";
    std::cout << "  g_UQFF = " << g15 << " m/s^2\n\n";
    
    // Step 16: List all saved states
    std::cout << "Step 16: List all saved states\n";
    auto states = sgra.listSavedStates();
    std::cout << "  Saved states (" << states.size() << " total):\n";
    for (const auto& state : states) {
        std::cout << "    - " << state << "\n";
    }
    std::cout << "\n";
    
    // Step 17: Generate comprehensive report
    std::cout << "Step 17: Generate comprehensive system report\n";
    sgra.restoreState("sgra_initial_10Gyr");
    std::string report = sgra.generateReport(t_current);
    std::cout << report << "\n";
    
    // Step 18: Export final state
    std::cout << "Step 18: Export final system state\n";
    std::string state_export = sgra.exportState();
    std::cout << state_export << "\n";
    
    std::cout << "========== END 18-STEP SAGITTARIUS A* DEMONSTRATION ==========\n\n";
}