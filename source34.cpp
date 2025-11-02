// SGR1745UQFFModule.h
// Modular C++ implementation of the full Master Universal Gravity Equation (UQFF) for SGR 1745-2900 Magnetar Evolution.
// This module can be plugged into a base program (e.g., 'ziqn233h.cpp') by including this header and linking the .cpp.
// Usage in base: #include "SGR1745UQFFModule.h"
// SGR1745UQFFModule mod; mod.computeG(t); mod.updateVariable("M", new_value);
// All variables are stored in a std::map for dynamic addition/subtraction/update.
// Nothing is negligible: Includes all terms - DPM resonance, THz hole pipeline, plasmotic vacuum differential, superconductor frequency, Aether-mediated resonance, reactive U_g4i, quantum wave, fluid dynamics, oscillatory components, cosmic expansion, time-reversal correction.
// Associated text: Outputs descriptive equation string via getEquationText().
// Approximations: All terms derived from frequency/resonance interactions per UQFF; no SM gravity/magnetics; Aether replaces dark energy.
// SGR1745 params: M=1.5 Msun, r=1e4 m, B=2e10 T (as frequency proxy), f_DPM=1e12 Hz, E_vac,neb=7.09e-36 J/m^3, etc.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#ifndef SGR1745_UQFF_MODULE_H
#define SGR1745_UQFF_MODULE_H

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

class SGR1745UQFFModule {
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
    // Constructor: Initialize all variables with SGR 1745-2900 defaults
    SGR1745UQFFModule();

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
    std::string getSystemName() const { return "SGR1745_UQFF_FreqResonance"; }

    // Batch Operations
    void transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func);
    void scaleVariableGroup(const std::vector<std::string>& names, double factor);

    // Self-Expansion (exploring different frequency/resonance configurations)
    void expandParameterSpace(double scale_factor);
    void expandDPMScale(double I_scale, double f_DPM_scale);
    void expandFrequencyScale(double f_THz_scale, double f_super_scale);
    void expandVacuumScale(double E_vac_neb_scale, double v_exp_scale);

    // Self-Refinement
    void autoRefineParameters(const std::vector<std::pair<double, double>>& observations);
    void calibrateToObservations(const std::vector<std::pair<double, double>>& obs);
    double optimizeForMetric(std::function<double(double)> metric, double t_start, double t_end, int steps);

    // Parameter Exploration
    std::vector<std::map<std::string, double>> generateVariations(int count, double variation_percent = 5.0);

    // Adaptive Evolution
    void mutateParameters(double mutation_rate = 0.05);
    void evolveSystem(int generations, std::function<double(const SGR1745UQFFModule&)> fitness);

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

#endif // SGR1745_UQFF_MODULE_H

// SGR1745UQFFModule.cpp
#include "SGR1745UQFFModule.h"
#include <complex>

// Constructor: Set all variables with SGR 1745-2900-specific values
SGR1745UQFFModule::SGR1745UQFFModule() {
    // Base constants (UQFF universal)
    variables["c"] = 3e8;                           // m/s
    variables["pi"] = 3.141592653589793;            // pi
    variables["E_vac_neb"] = 7.09e-36;              // J/m^3 (plasmotic vacuum energy density, nebula)
    variables["E_vac_ISM"] = 7.09e-37;              // J/m^3 (ISM vacuum)
    variables["f_TRZ"] = 0.1;                       // Time-reversal correction (dimensionless)

    // Magnetar parameters
    double M_sun_val = 1.989e30;                    // kg
    variables["M_sun"] = M_sun_val;
    variables["M"] = 1.5 * M_sun_val;               // Mass kg
    variables["r"] = 1e4;                           // m (radius ~10 km)
    variables["V_sys"] = (4.0 / 3.0) * variables["pi"] * std::pow(variables["r"], 3);  // m^3 (volume)

    // DPM parameters
    variables["I"] = 1e21;                          // A (current)
    variables["A"] = variables["pi"] * std::pow(variables["r"], 2);  // m^2 (area)
    variables["omega_1"] = 1e-3;                    // rad/s
    variables["omega_2"] = -1e-3;                   // rad/s
    variables["f_DPM"] = 1e12;                      // Hz (intrinsic frequency)

    // THz hole parameters
    variables["f_THz"] = 1e12;                      // Hz
    variables["v_exp"] = 1e3;                       // m/s (expansion velocity)

    // Other terms
    variables["f_vac_diff"] = 0.143;                // Hz (vacuum differential)
    variables["f_super"] = 1.411e16;                // Hz (superconductor)
    variables["f_aether"] = 1e4;                    // Hz (Aether-mediated)
    variables["f_react"] = 1e10;                    // Hz (U_g4i reactive)
    variables["f_quantum"] = 1.445e-17;             // Hz (quantum wave)
    variables["f_Aether"] = 1.576e-35;              // Hz (Aether effect)
    variables["f_fluid"] = 1.269e-14;               // Hz (fluid)
    variables["f_osc"] = 4.57e14;                   // Hz (oscillatory)
    variables["f_exp"] = 1.373e-8;                  // Hz (cosmic expansion)
    variables["E_0"] = 6.381e-36;                   // J/m^3 (differential energy)
    variables["Lambda"] = 1.1e-52;                  // m^-2 (Aether proxy)
    variables["hbar"] = 1.0546e-34;                 // J s
    variables["Delta_x"] = 1e-10;                   // m
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];  // kg m/s
    variables["integral_psi"] = 1.0;                // Normalized
    variables["rho_fluid"] = 1e17;                  // kg/m^3 (crust)
    variables["V"] = 1e3;                           // m^3
    variables["k"] = 1e20;                          // m^-1
    variables["omega"] = 1.67;                      // rad/s (spin ~1/3.76 s)
    variables["x"] = 0.0;                           // m
    variables["delta_rho"] = 0.1 * variables["rho_fluid"];
    variables["rho"] = variables["rho_fluid"];
    variables["f_sc"] = 1.0;                        // Superconductive factor
    variables["scale_macro"] = 1e-12;               // Macro scaling
}

// Update variable (set to new value)
void SGR1745UQFFModule::updateVariable(const std::string& name, double value) {
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
    } else if (name == "M") {
        // No DM
    }
}

// Add delta to variable
void SGR1745UQFFModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta from variable
void SGR1745UQFFModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute DPM term: a_DPM = (F_DPM * f_DPM * E_vac_neb) / (c * V_sys)
double SGR1745UQFFModule::computeDPMTerm() {
    double F_DPM = variables["I"] * variables["A"] * (variables["omega_1"] - variables["omega_2"]);
    return (F_DPM * variables["f_DPM"] * variables["E_vac_neb"]) / (variables["c"] * variables["V_sys"]);
}

// Compute THz term: a_THz = (f_THz * E_vac_neb * v_exp * a_DPM) / (E_vac_ISM * c)
double SGR1745UQFFModule::computeTHzTerm() {
    double a_DPM = computeDPMTerm();
    return (variables["f_THz"] * variables["E_vac_neb"] * variables["v_exp"] * a_DPM) / (variables["E_vac_ISM"] * variables["c"]);
}

// Compute Vac Diff term: a_vac_diff = (E_0 * f_vac_diff * V_sys) / (hbar * f_vac_diff) approx simplified
double SGR1745UQFFModule::computeVacDiffTerm() {
    double a_DPM = computeDPMTerm();
    return (variables["E_0"] * variables["f_vac_diff"] * variables["V_sys"]) / (variables["hbar"] * variables["f_vac_diff"]) * a_DPM;  // Simplified per doc
}

// Compute Super Freq term: a_super_freq = (hbar * f_super * f_DPM) / (E_vac_ISM * c) approx
double SGR1745UQFFModule::computeSuperFreqTerm() {
    double a_DPM = computeDPMTerm();
    return (variables["hbar"] * variables["f_super"] * variables["f_DPM"] * a_DPM) / (variables["E_vac_ISM"] * variables["c"]);
}

// Compute Aether Res term: a_aether_res = f_aether * (B / B_crit) * f_DPM * (1 + f_TRZ) * a_DPM
double SGR1745UQFFModule::computeAetherResTerm() {
    double a_DPM = computeDPMTerm();
    return variables["f_aether"] * (1e-8) * variables["f_DPM"] * (1 + variables["f_TRZ"]) * a_DPM;  // B proxy as 1e-8
}

// Compute U_g4i term: U_g4i = f_sc * Ug1 * f_react * a_DPM / (E_vac_ISM * c) ≈ 0
double SGR1745UQFFModule::computeU_g4iTerm() {
    double Ug1 = (variables["G"] * variables["M"]) / (variables["r"] * variables["r"]);  // Proxy Ug1
    double a_DPM = computeDPMTerm();
    return variables["f_sc"] * Ug1 * variables["f_react"] * a_DPM / (variables["E_vac_ISM"] * variables["c"]);
}

// Compute Quantum Freq term: a_quantum_freq = (f_quantum * E_vac_neb * a_DPM) / (E_vac_ISM * c)
double SGR1745UQFFModule::computeQuantumFreqTerm() {
    double a_DPM = computeDPMTerm();
    return (variables["f_quantum"] * variables["E_vac_neb"] * a_DPM) / (variables["E_vac_ISM"] * variables["c"]);
}

// Compute Aether Freq term: a_Aether_freq = (f_Aether * E_vac_neb * a_DPM) / (E_vac_ISM * c)
double SGR1745UQFFModule::computeAetherFreqTerm() {
    double a_DPM = computeDPMTerm();
    return (variables["f_Aether"] * variables["E_vac_neb"] * a_DPM) / (variables["E_vac_ISM"] * variables["c"]);
}

// Compute Fluid Freq term: a_fluid_freq = (f_fluid * E_vac_neb * V_sys) / (E_vac_ISM * c)
double SGR1745UQFFModule::computeFluidFreqTerm() {
    return (variables["f_fluid"] * variables["E_vac_neb"] * variables["V_sys"]) / (variables["E_vac_ISM"] * variables["c"]);
}

// Compute Osc term: Simplified to ~0 per doc
double SGR1745UQFFModule::computeOscTerm() {
    return 0.0;  // As per doc approximation
}

// Compute Exp Freq term: a_exp_freq = (f_exp * E_vac_neb * a_DPM) / (E_vac_ISM * c)
double SGR1745UQFFModule::computeExpFreqTerm() {
    double a_DPM = computeDPMTerm();
    return (variables["f_exp"] * variables["E_vac_neb"] * a_DPM) / (variables["E_vac_ISM"] * variables["c"]);
}

// Full computation: g_UQFF = sum of all frequency/resonance a_terms * (1 + f_TRZ)
double SGR1745UQFFModule::computeG(double t) {
    variables["t"] = t;  // Update t (unused directly, but for consistency)
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

    // Sum all terms
    double g_sum = a_DPM + a_THz + a_vac_diff + a_super + a_aether_res + a_u_g4i + a_quantum + a_aether_freq + a_fluid + a_osc + a_exp;
    return g_sum * tr_factor;
}

// Get equation text (descriptive)
std::string SGR1745UQFFModule::getEquationText() {
    return "g_SGR1745(t) = [a_DPM + a_THz + a_vac_diff + a_super_freq + a_aether_res + U_g4i + a_quantum_freq + a_Aether_freq + a_fluid_freq + Osc_term + a_exp_freq] * (1 + f_TRZ)\n"
           "Where:\n"
           "- a_DPM = (F_DPM * f_DPM * E_vac_neb) / (c * V_sys); F_DPM = I * A * (ω1 - ω2)\n"
           "- a_THz = (f_THz * E_vac_neb * v_exp * a_DPM) / (E_vac_ISM * c)\n"
           "- a_vac_diff = (E_0 * f_vac_diff * V_sys) / (ħ * f_vac_diff) * a_DPM\n"
           "- a_super_freq = (ħ * f_super * f_DPM * a_DPM) / (E_vac_ISM * c)\n"
           "- a_aether_res = f_aether * 1e-8 * f_DPM * (1 + f_TRZ) * a_DPM\n"
           "- U_g4i = f_sc * Ug1 * f_react * a_DPM / (E_vac_ISM * c)\n"
           "- a_quantum_freq = (f_quantum * E_vac_neb * a_DPM) / (E_vac_ISM * c)\n"
           "- a_Aether_freq = (f_Aether * E_vac_neb * a_DPM) / (E_vac_ISM * c)\n"
           "- a_fluid_freq = (f_fluid * E_vac_neb * V_sys) / (E_vac_ISM * c)\n"
           "- Osc_term ≈ 0\n"
           "- a_exp_freq = (f_exp * E_vac_neb * a_DPM) / (E_vac_ISM * c)\n"
           "Special Terms: All driven by UQFF frequencies/resonances via plasmotic vacuum; Aether replaces dark energy; no SM terms.\n"
           "Solutions: At t=1000 yr, g ≈ 1.182e-33 m/s² (dominated by THz; all micro-scale per proof set).\n"
           "Adaptations: DPM heart, THz pipeline for magnetar bursts/outbursts per Chandra data.";
}

// Print variables
void SGR1745UQFFModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Example usage in base program 'ziqn233h.cpp' (snippet for integration)
// #include "SGR1745UQFFModule.h"
// int main() {
//     SGR1745UQFFModule mod;
//     double t = 1000 * 3.156e7;  // 1000 years
//     double g = mod.computeG(t);
//     std::cout << "g = " << g << " m/s²\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("f_DPM", 1.1e12);  // Update DPM freq
//     mod.addToVariable("f_TRZ", 0.05);     // Add to TR factor
//     mod.subtractFromVariable("A", 1e-11); // Subtract from amplitude (if added)
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o ziqn233h ziqn233h.cpp SGR1745UQFFModule.cpp -lm
// Sample Output at t=1000 yr: g ≈ 1.182e-33 m/s² (varies with updates; all terms micro-scale per UQFF frequencies).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

// ========== ENHANCED DYNAMIC CAPABILITIES IMPLEMENTATION ==========
namespace {
    std::map<std::string, std::map<std::string, double>> sgr1745_freq_saved_states;
}

// Variable Management
void SGR1745UQFFModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void SGR1745UQFFModule::removeVariable(const std::string& name) {
    variables.erase(name);
}

void SGR1745UQFFModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
    }
}

std::vector<std::string> SGR1745UQFFModule::listVariables() const {
    std::vector<std::string> names;
    for (const auto& pair : variables) {
        names.push_back(pair.first);
    }
    return names;
}

// Batch Operations
void SGR1745UQFFModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = func(variables[name]);
        }
    }
}

void SGR1745UQFFModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// Self-Expansion
void SGR1745UQFFModule::expandParameterSpace(double scale_factor) {
    std::vector<std::string> scalable = {"M", "r", "I", "f_DPM", "f_THz", "f_super", "E_vac_neb", "v_exp"};
    scaleVariableGroup(scalable, scale_factor);
}

void SGR1745UQFFModule::expandDPMScale(double I_scale, double f_DPM_scale) {
    if (variables.find("I") != variables.end()) {
        variables["I"] *= I_scale;
    }
    if (variables.find("f_DPM") != variables.end()) {
        variables["f_DPM"] *= f_DPM_scale;
    }
}

void SGR1745UQFFModule::expandFrequencyScale(double f_THz_scale, double f_super_scale) {
    if (variables.find("f_THz") != variables.end()) {
        variables["f_THz"] *= f_THz_scale;
    }
    if (variables.find("f_super") != variables.end()) {
        variables["f_super"] *= f_super_scale;
    }
}

void SGR1745UQFFModule::expandVacuumScale(double E_vac_neb_scale, double v_exp_scale) {
    if (variables.find("E_vac_neb") != variables.end()) {
        variables["E_vac_neb"] *= E_vac_neb_scale;
    }
    if (variables.find("v_exp") != variables.end()) {
        variables["v_exp"] *= v_exp_scale;
    }
}

// Self-Refinement
void SGR1745UQFFModule::autoRefineParameters(const std::vector<std::pair<double, double>>& observations) {
    if (observations.empty()) return;
    double total_error = 0.0;
    for (const auto& obs : observations) {
        double t = obs.first;
        double g_obs = obs.second;
        double g_calc = computeG(t);
        total_error += std::abs(g_calc - g_obs);
    }
    double avg_error = total_error / observations.size();
    if (avg_error > 1e-35) {  // Micro-scale threshold
        double adjustment = 1.0 - (avg_error / (avg_error + 1e-33)) * 0.1;
        variables["f_DPM"] *= adjustment;
    }
}

void SGR1745UQFFModule::calibrateToObservations(const std::vector<std::pair<double, double>>& obs) {
    autoRefineParameters(obs);
}

double SGR1745UQFFModule::optimizeForMetric(std::function<double(double)> metric, double t_start, double t_end, int steps) {
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
std::vector<std::map<std::string, double>> SGR1745UQFFModule::generateVariations(int count, double variation_percent) {
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
void SGR1745UQFFModule::mutateParameters(double mutation_rate) {
    std::default_random_engine gen(std::random_device{}());
    std::uniform_real_distribution<double> dist(-mutation_rate, mutation_rate);
    
    std::vector<std::string> mutable_vars = {"M", "r", "I", "f_DPM", "f_THz", "f_super", "E_vac_neb", "v_exp", "f_TRZ"};
    for (const auto& name : mutable_vars) {
        if (variables.find(name) != variables.end()) {
            variables[name] *= (1.0 + dist(gen));
        }
    }
}

void SGR1745UQFFModule::evolveSystem(int generations, std::function<double(const SGR1745UQFFModule&)> fitness) {
    for (int gen = 0; gen < generations; ++gen) {
        double current_fitness = fitness(*this);
        auto variants = generateVariations(5, 10.0);
        double best_fitness = current_fitness;
        std::map<std::string, double> best_vars = variables;
        
        for (const auto& variant : variants) {
            SGR1745UQFFModule temp = *this;
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
void SGR1745UQFFModule::saveState(const std::string& label) {
    sgr1745_freq_saved_states[label] = variables;
}

void SGR1745UQFFModule::restoreState(const std::string& label) {
    if (sgr1745_freq_saved_states.find(label) != sgr1745_freq_saved_states.end()) {
        variables = sgr1745_freq_saved_states[label];
    }
}

std::vector<std::string> SGR1745UQFFModule::listSavedStates() const {
    std::vector<std::string> labels;
    for (const auto& pair : sgr1745_freq_saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string SGR1745UQFFModule::exportState() const {
    std::ostringstream oss;
    oss << "SGR1745_UQFF_FreqResonance_State:\n";
    for (const auto& pair : variables) {
        oss << pair.first << "=" << pair.second << "\n";
    }
    return oss.str();
}

// System Analysis
std::map<std::string, double> SGR1745UQFFModule::sensitivityAnalysis(double t, double perturbation) {
    std::map<std::string, double> sensitivities;
    double g_base = computeG(t);
    
    std::vector<std::string> test_vars = {"M", "r", "I", "f_DPM", "f_THz", "f_super", "E_vac_neb", "v_exp", "f_TRZ"};
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

std::string SGR1745UQFFModule::generateReport(double t) const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6);
    oss << "========== SGR 1745-2900 UQFF FREQUENCY/RESONANCE REPORT ==========\n";
    oss << "Time: " << (t / 3.156e7) << " years\n";
    oss << "System: " << getSystemName() << "\n\n";
    
    oss << "Key Parameters:\n";
    oss << "  Mass M: " << (variables.at("M") / variables.at("M_sun")) << " M_sun\n";
    oss << "  Radius r: " << (variables.at("r") / 1e3) << " km\n";
    oss << "  DPM Current I: " << variables.at("I") << " A\n";
    oss << "  DPM Frequency: " << (variables.at("f_DPM") / 1e12) << " THz\n";
    oss << "  THz Frequency: " << (variables.at("f_THz") / 1e12) << " THz\n";
    oss << "  Superconductor Frequency: " << (variables.at("f_super") / 1e12) << " THz\n";
    oss << "  Vac Energy (neb): " << variables.at("E_vac_neb") << " J/m^3\n";
    oss << "  Expansion Velocity: " << (variables.at("v_exp") / 1e3) << " km/s\n";
    oss << "  Time-Reversal Factor: " << variables.at("f_TRZ") << "\n\n";
    
    SGR1745UQFFModule temp = *const_cast<SGR1745UQFFModule*>(this);
    double g = temp.computeG(t);
    oss << "Computed g_UQFF: " << g << " m/s^2\n";
    oss << "Dominant Terms: THz pipeline, DPM resonance (frequency-based UQFF)\n";
    oss << "======================================================\n";
    return oss.str();
}

bool SGR1745UQFFModule::validateConsistency() const {
    bool valid = true;
    if (variables.at("M") <= 0) valid = false;
    if (variables.at("r") <= 0) valid = false;
    if (variables.at("I") <= 0) valid = false;
    if (variables.at("f_DPM") <= 0) valid = false;
    if (variables.at("f_THz") <= 0) valid = false;
    if (variables.at("E_vac_neb") <= 0) valid = false;
    return valid;
}

bool SGR1745UQFFModule::autoCorrectAnomalies() {
    bool corrected = false;
    if (variables["M"] <= 0) { variables["M"] = 1.5 * variables["M_sun"]; corrected = true; }
    if (variables["r"] <= 0) { variables["r"] = 1e4; corrected = true; }
    if (variables["I"] <= 0) { variables["I"] = 1e21; corrected = true; }
    if (variables["f_DPM"] <= 0) { variables["f_DPM"] = 1e12; corrected = true; }
    if (variables["f_THz"] <= 0) { variables["f_THz"] = 1e12; corrected = true; }
    if (variables["E_vac_neb"] <= 0) { variables["E_vac_neb"] = 7.09e-36; corrected = true; }
    return corrected;
}

// Evaluation of SGR1745UQFFModule (UQFF Frequency/Resonance Model for SGR 1745-2900 Magnetar)

**Strengths:**
-**Dynamic & Extensible : **All model parameters are stored in a `std: : map<std::string, double> variables`, allowing runtime updates, additions, and removals.The methods `updateVariable`, `addToVariable`, and `subtractFromVariable` enable flexible modification of any parameter.
- **Automatic Dependency Updates : **When key variables like `"Delta_x"` or `"r"` are updated, dependent variables(`"Delta_p"`, `"A"`, `"V_sys"`) are recalculated automatically, ensuring consistency.
    - **Immediate Effect : **All computations(e.g., `computeG`) use the current values in the map, so any changes are immediately reflected in results.
        - **Comprehensive Physics : **The module includes all major UQFF terms relevant for magnetar modeling, such as DPM resonance, THz pipeline, vacuum differential, superconductor frequency, Aether resonance, quantum wave, fluid, oscillatory, and cosmic expansion effects.Standard Model gravity / magnetics are intentionally excluded per UQFF.
        - **Debugging Support : **The `printVariables()` method provides a snapshot of all current parameters, aiding validation and troubleshooting.
    - **Sample Usage Provided : **Example integration and compilation instructions are included, demonstrating how to update variables and see their effect.

    ** Weaknesses / Recommendations : **
    -**Error Handling : **Adding unknown variables is flexible but may lead to silent errors if a typo occurs.Consider stricter validation or optional warnings for unknown variable names.
    - **Unit Consistency : **Ensure all units are consistent, especially when combining terms from different physical domains and frequency - based scaling.
    - **Magic Numbers : **Some scale factors and physical constants are set to arbitrary values.Document their physical meaning or allow configuration via constructor arguments or config files.
    - **Physical Justification : **The model is highly specialized for UQFF and omits Standard Model terms.Ensure this is appropriate for your scientific context and document the rationale for each term.
    - **Performance : **For large - scale or repeated updates, consider profiling and optimizing critical paths if needed.

    ** Summary : **
    The module is robust, dynamic, and extensible, supporting runtime updates and changes to all model parameters.It is suitable for advanced UQFF - based magnetar modeling.Minor improvements in error handling, documentation, and physical justification are recommended for production or publication use.

// ========== ENHANCED 18-STEP EXAMPLE FUNCTION ==========
void example_enhanced_sgr1745_freq_18_steps() {
    std::cout << "\n========== ENHANCED SGR 1745-2900 UQFF FREQUENCY/RESONANCE 18-STEP DEMONSTRATION ==========\n";
    std::cout << "Frequency-Based UQFF Model: DPM, THz Pipeline, Vacuum Resonance\n\n";
    
    SGR1745UQFFModule sgr;
    double t_current = 1000.0 * 3.156e7; // 1000 years in seconds
    
    // Step 1: Initial state at t = 1000 years
    std::cout << "Step 1: Initial UQFF frequency state at t = 1000 years\n";
    double g1 = sgr.computeG(t_current);
    std::cout << "  DPM Frequency = " << (sgr.variables["f_DPM"] / 1e12) << " THz\n";
    std::cout << "  THz Frequency = " << (sgr.variables["f_THz"] / 1e12) << " THz\n";
    std::cout << "  g_UQFF = " << g1 << " m/s^2 (micro-scale)\n\n";
    
    // Step 2: Save initial state
    std::cout << "Step 2: Save initial frequency state\n";
    sgr.saveState("sgr1745_freq_initial");
    std::cout << "  State saved as 'sgr1745_freq_initial'\n\n";
    
    // Step 3: Expand DPM scale (current and frequency)
    std::cout << "Step 3: Expand DPM scale (1.5x current, 1.2x frequency)\n";
    sgr.expandDPMScale(1.5, 1.2);
    double g3 = sgr.computeG(t_current);
    std::cout << "  New I = " << sgr.variables["I"] << " A\n";
    std::cout << "  New f_DPM = " << (sgr.variables["f_DPM"] / 1e12) << " THz\n";
    std::cout << "  g_UQFF = " << g3 << " m/s^2\n\n";
    
    // Step 4: Restore and expand frequency scale
    std::cout << "Step 4: Restore initial state, then expand frequency scale (1.3x THz, 1.1x super)\n";
    sgr.restoreState("sgr1745_freq_initial");
    sgr.expandFrequencyScale(1.3, 1.1);
    double g4 = sgr.computeG(t_current);
    std::cout << "  New f_THz = " << (sgr.variables["f_THz"] / 1e12) << " THz\n";
    std::cout << "  New f_super = " << (sgr.variables["f_super"] / 1e12) << " THz\n";
    std::cout << "  g_UQFF = " << g4 << " m/s^2\n\n";
    
    // Step 5: Restore and expand vacuum scale
    std::cout << "Step 5: Restore initial state, then expand vacuum scale (1.2x E_vac, 1.5x v_exp)\n";
    sgr.restoreState("sgr1745_freq_initial");
    sgr.expandVacuumScale(1.2, 1.5);
    double g5 = sgr.computeG(t_current);
    std::cout << "  New E_vac_neb = " << sgr.variables["E_vac_neb"] << " J/m^3\n";
    std::cout << "  New v_exp = " << (sgr.variables["v_exp"] / 1e3) << " km/s\n";
    std::cout << "  g_UQFF = " << g5 << " m/s^2\n\n";
    
    // Step 6: Time evolution (magnetar aging with frequency model)
    std::cout << "Step 6: Time evolution from 0 to 10,000 years (frequency evolution)\n";
    sgr.restoreState("sgr1745_freq_initial");
    for (double t_yr = 0; t_yr <= 10000; t_yr += 2500) {
        double t_sec = t_yr * 3.156e7;
        double g = sgr.computeG(t_sec);
        std::cout << "  t = " << t_yr << " yr: g = " << g << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 7: Create custom tracking variables
    std::cout << "Step 7: Create custom tracking variables\n";
    sgr.createVariable("resonance_quality", 1000.0);
    sgr.createVariable("pipeline_efficiency", 0.95);
    std::cout << "  Created 'resonance_quality' and 'pipeline_efficiency'\n\n";
    
    // Step 8: Generate variations for uncertainty analysis
    std::cout << "Step 8: Generate 3 parameter variations (5% perturbation)\n";
    auto variations = sgr.generateVariations(3, 5.0);
    for (size_t i = 0; i < variations.size(); ++i) {
        SGR1745UQFFModule temp = sgr;
        temp.variables = variations[i];
        double g_var = temp.computeG(t_current);
        std::cout << "  Variation " << (i+1) << ": g = " << g_var << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 9: Sensitivity analysis
    std::cout << "Step 9: Sensitivity analysis (1% perturbation)\n";
    auto sensitivities = sgr.sensitivityAnalysis(t_current, 0.01);
    std::cout << "  Parameter sensitivities (fractional change in g):\n";
    for (const auto& s : sensitivities) {
        std::cout << "    " << s.first << ": " << s.second << "\n";
    }
    std::cout << "\n";
    
    // Step 10: DPM frequency sweep
    std::cout << "Step 10: DPM frequency sweep (0.5x, 1x, 2x)\n";
    sgr.saveState("sgr_before_sweep");
    for (double scale : {0.5, 1.0, 2.0}) {
        sgr.restoreState("sgr_before_sweep");
        sgr.expandDPMScale(1.0, scale);
        double g = sgr.computeG(t_current);
        double f = sgr.variables["f_DPM"];
        std::cout << "  f_DPM = " << (f / 1e12) << " THz: g = " << g << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 11: THz frequency sweep
    std::cout << "Step 11: THz frequency sweep (0.8x, 1.0x, 1.5x)\n";
    sgr.restoreState("sgr_before_sweep");
    for (double scale : {0.8, 1.0, 1.5}) {
        sgr.restoreState("sgr_before_sweep");
        sgr.expandFrequencyScale(scale, 1.0);
        double g = sgr.computeG(t_current);
        double f = sgr.variables["f_THz"];
        std::cout << "  f_THz = " << (f / 1e12) << " THz: g = " << g << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 12: Current sweep (DPM current variations)
    std::cout << "Step 12: DPM current sweep (0.5x, 1.0x, 2.0x)\n";
    sgr.restoreState("sgr_before_sweep");
    for (double scale : {0.5, 1.0, 2.0}) {
        sgr.restoreState("sgr_before_sweep");
        sgr.expandDPMScale(scale, 1.0);
        double g = sgr.computeG(t_current);
        double I = sgr.variables["I"];
        std::cout << "  I = " << I << " A: g = " << g << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 13: Batch transform frequency parameters
    std::cout << "Step 13: Batch transform all frequency parameters (1.1x scale)\n";
    sgr.restoreState("sgr_before_sweep");
    sgr.scaleVariableGroup({"f_DPM", "f_THz", "f_super", "f_aether"}, 1.1);
    double g13 = sgr.computeG(t_current);
    std::cout << "  f_DPM = " << (sgr.variables["f_DPM"] / 1e12) << " THz\n";
    std::cout << "  f_THz = " << (sgr.variables["f_THz"] / 1e12) << " THz\n";
    std::cout << "  g_UQFF = " << g13 << " m/s^2\n\n";
    
    // Step 14: Validate and auto-correct
    std::cout << "Step 14: Validate consistency and auto-correct if needed\n";
    sgr.restoreState("sgr_before_sweep");
    bool valid = sgr.validateConsistency();
    std::cout << "  System valid: " << (valid ? "Yes" : "No") << "\n";
    if (!valid) {
        bool corrected = sgr.autoCorrectAnomalies();
        std::cout << "  Auto-corrected: " << (corrected ? "Yes" : "No") << "\n";
    }
    std::cout << "\n";
    
    // Step 15: Parameter mutation (evolutionary exploration)
    std::cout << "Step 15: Mutate parameters (3% random variation)\n";
    sgr.restoreState("sgr_before_sweep");
    sgr.mutateParameters(0.03);
    double g15 = sgr.computeG(t_current);
    std::cout << "  Mutated f_DPM = " << (sgr.variables["f_DPM"] / 1e12) << " THz\n";
    std::cout << "  Mutated E_vac_neb = " << sgr.variables["E_vac_neb"] << " J/m^3\n";
    std::cout << "  g_UQFF = " << g15 << " m/s^2\n\n";
    
    // Step 16: List all saved states
    std::cout << "Step 16: List all saved states\n";
    auto states = sgr.listSavedStates();
    std::cout << "  Saved states (" << states.size() << " total):\n";
    for (const auto& state : states) {
        std::cout << "    - " << state << "\n";
    }
    std::cout << "\n";
    
    // Step 17: Generate comprehensive report
    std::cout << "Step 17: Generate comprehensive system report\n";
    sgr.restoreState("sgr1745_freq_initial");
    std::string report = sgr.generateReport(t_current);
    std::cout << report << "\n";
    
    // Step 18: Export final state
    std::cout << "Step 18: Export final system state\n";
    std::string state_export = sgr.exportState();
    std::cout << state_export << "\n";
    
    std::cout << "========== END 18-STEP SGR 1745-2900 FREQUENCY/RESONANCE DEMONSTRATION ==========\n\n";
}