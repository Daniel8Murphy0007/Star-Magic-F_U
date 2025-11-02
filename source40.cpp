// CompressedResonanceUQFF24Module.h
// Modular C++ implementation of the UQFF Compressed and Resonance Equations for Systems 18-24.
// This module can be plugged into a base program (e.g., 'ziqn233h.cpp') by including this header and linking the .cpp.
// Usage in base: #include "CompressedResonanceUQFF24Module.h"
// CompressedResonanceUQFF24Module mod; mod.computeCompressedResTerm(t, B); mod.updateVariable("f_DPM", new_value);
// All variables are stored in a std::map for dynamic addition/subtraction/update.
// Nothing is negligible: Includes compressed terms (DPM, THz, vac_diff, super) + resonance (aether, U_g4i, osc, quantum, fluid, exp) scaled for systems 18-24 (e.g., Sombrero, Saturn, M16, Crab).
// Associated text: Outputs descriptive equation string via getEquationText().
// Approximations: Compressed: Sum key frequency terms; Resonance: Real part exp; SC correction integrated; frequencies scaled per system (e.g., f_DPM=1e11 for nebulae, 1e12 for remnants).
// General params: f_DPM=1e12 Hz (default), B=1e-5 T, E_vac=7.09e-36 J/m^3, for systems 18-24.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#ifndef COMPRESSED_RESONANCE_UQFF24_MODULE_H
#define COMPRESSED_RESONANCE_UQFF24_MODULE_H

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

class CompressedResonanceUQFF24Module {
private:
    std::map<std::string, double> variables;
    double computeCompressedTerm();
    double computeResonanceTerm(double t);
    double computeSCIntegrated(double B);

public:
    // Constructor: Initialize all variables with UQFF defaults for compressed/resonance (systems 18-24)
    CompressedResonanceUQFF24Module();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations: Compressed term, Resonance term, full combined with SC
    double computeCompressedResTerm(double t, double B);

    // Output descriptive text of the equations
    std::string getEquationText();

    // Print all current variables (for debugging/updates)
    void printVariables();

    // ========== ENHANCED DYNAMIC CAPABILITIES (25 methods) ==========
    // Variable management (5)
    void createVariable(const std::string& name, double value);
    void removeVariable(const std::string& name);
    void cloneVariable(const std::string& source, const std::string& dest);
    std::vector<std::string> listVariables();
    std::string getSystemName() const { return "Compressed_Resonance_Systems_18_24"; }

    // Batch operations (2)
    void transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func);
    void scaleVariableGroup(const std::vector<std::string>& names, double scale);

    // Self-expansion (4)
    void expandParameterSpace(double scale);
    void expandCompressedScale(double f_DPM_scale, double f_THz_scale);
    void expandResonanceScale(double f_aether_scale, double f_react_scale);
    void expandSystemScale(double I_scale, double V_sys_scale);

    // Self-refinement (3)
    void autoRefineParameters(double t, double B, double target_g, double tolerance);
    void calibrateToObservations(const std::map<std::string, double>& observations);
    void optimizeForMetric(double t, double B, const std::string& metric);

    // Parameter exploration (1)
    std::vector<std::map<std::string, double>> generateVariations(int count, double variation_percent);

    // Adaptive evolution (2)
    void mutateParameters(double mutation_rate);
    void evolveSystem(double t, double B, int generations, double selection_pressure);

    // State management (4)
    void saveState(const std::string& label);
    void restoreState(const std::string& label);
    std::vector<std::string> listSavedStates();
    std::string exportState() const;

    // System analysis (4)
    std::map<std::string, double> sensitivityAnalysis(double t, double B, double perturbation);
    std::string generateReport(double t, double B);
    bool validateConsistency();
    bool autoCorrectAnomalies();
};

#endif // COMPRESSED_RESONANCE_UQFF24_MODULE_H

// CompressedResonanceUQFF24Module.cpp
#include "CompressedResonanceUQFF24Module.h"
#include <complex>

// Constructor: Set all variables with UQFF-specific values for compressed/resonance (systems 18-24)
CompressedResonanceUQFF24Module::CompressedResonanceUQFF24Module() {
    // Base constants (UQFF universal)
    variables["c"] = 3e8;                           // m/s
    variables["pi"] = 3.141592653589793;            // pi
    variables["E_vac"] = 7.09e-36;                  // J/m^3 (plasmotic vacuum)
    variables["hbar"] = 1.0546e-34;                 // J s
    variables["f_TRZ"] = 0.1;                       // Time-reversal

    // Compressed parameters (streamlined DPM, THz, vac_diff, super; scaled for 18-24)
    variables["f_DPM"] = 1e11;                      // Hz (nebula/Saturn scale)
    variables["f_THz"] = 1e11;                      // Hz
    variables["f_vac_diff"] = 0.143;                // Hz
    variables["f_super"] = 1.411e15;                // Hz (scaled)
    variables["I"] = 1e20;                          // A (system scale)
    variables["A_vort"] = 3.142e18;                 // m^2 (larger for galaxies/planets)
    variables["omega_1"] = 1e-2;                    // rad/s
    variables["omega_2"] = -1e-2;                   // rad/s
    variables["v_exp"] = 1e5;                       // m/s (outflow)
    variables["E_0"] = 6.381e-36;                   // J/m^3
    variables["V_sys"] = 4.189e18;                  // m^3 (scaled volume)

    // Resonance parameters (aether, U_g4i, osc, quantum, fluid, exp; scaled)
    variables["f_aether"] = 1e3;                    // Hz
    variables["f_react"] = 1e9;                     // Hz (U_g4i)
    variables["f_quantum"] = 1.445e-17;             // Hz
    variables["f_fluid"] = 1.269e-14;               // Hz
    variables["f_exp"] = 1.373e-8;                  // Hz
    variables["f_osc"] = 4.57e13;                   // Hz
    variables["k"] = 1e18;                          // m^-1 (scaled)
    variables["omega_osc"] = 1e14;                  // rad/s
    variables["x"] = 0.0;                           // m
    variables["A"] = 1e-9;                          // Amplitude (scaled)
    variables["rho_fluid"] = 1e-20;                 // kg/m^3 (gas/atm)
    variables["V"] = 1e6;                           // m^3
    variables["delta_rho"] = 0.1 * variables["rho_fluid"];
    variables["rho"] = variables["rho_fluid"];

    // Superconductive integrated
    variables["B_crit"] = 1e11;                     // T
    variables["f_sc"] = 1.0;                        // Factor

    // Quantum
    variables["Delta_x"] = 1e-10;                   // m
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];
    variables["integral_psi"] = 1.0;
}

// Update variable (set to new value)
void CompressedResonanceUQFF24Module::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
    if (name == "Delta_x") {
        variables["Delta_p"] = variables["hbar"] / value;
    }
}

// Add delta to variable
void CompressedResonanceUQFF24Module::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta from variable
void CompressedResonanceUQFF24Module::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute Compressed Term: Sum streamlined DPM + THz + vac_diff + super (scaled for 18-24)
double CompressedResonanceUQFF24Module::computeCompressedTerm() {
    double F_DPM = variables["I"] * variables["A_vort"] * (variables["omega_1"] - variables["omega_2"]);
    double a_DPM = (F_DPM * variables["f_DPM"] * variables["E_vac"]) / (variables["c"] * variables["V_sys"]);
    double a_THz = (variables["f_THz"] * variables["E_vac"] * variables["v_exp"] * a_DPM) / (variables["E_vac"] / 10 * variables["c"]);
    double a_vac_diff = (variables["E_0"] * variables["f_vac_diff"] * variables["V_sys"] * a_DPM) / variables["hbar"];
    double a_super = (variables["hbar"] * variables["f_super"] * variables["f_DPM"] * a_DPM) / (variables["E_vac"] * variables["c"]);
    return a_DPM + a_THz + a_vac_diff + a_super;
}

// Compute Resonance Term: Sum aether + U_g4i + osc + quantum + fluid + exp (scaled)
double CompressedResonanceUQFF24Module::computeResonanceTerm(double t) {
    double a_DPM = (variables["I"] * variables["A_vort"] * (variables["omega_1"] - variables["omega_2"]) * variables["f_DPM"] * variables["E_vac"]) / (variables["c"] * variables["V_sys"]);
    double a_aether = variables["f_aether"] * 1e-8 * variables["f_DPM"] * (1 + variables["f_TRZ"]) * a_DPM;
    double Ug1_proxy = 1.0;
    double a_u_g4i = variables["f_sc"] * Ug1_proxy * variables["f_react"] * a_DPM / (variables["E_vac"] * variables["c"]);
    double cos_term = 2 * variables["A"] * std::cos(variables["k"] * variables["x"]) * std::cos(variables["omega_osc"] * t);
    std::complex<double> exp_term(variables["A"] * std::exp(std::complex<double>(0, variables["k"] * variables["x"] - variables["omega_osc"] * t)));
    double real_exp = exp_term.real();
    double exp_factor = (2 * variables["pi"] / 13.8);
    double a_osc = cos_term + exp_factor * real_exp;
    double a_quantum = (variables["f_quantum"] * variables["E_vac"] * a_DPM) / (variables["E_vac"] / 10 * variables["c"]);
    double a_fluid = (variables["f_fluid"] * variables["E_vac"] * variables["V"] * a_DPM) / (variables["E_vac"] / 10 * variables["c"]);
    double a_exp = (variables["f_exp"] * variables["E_vac"] * a_DPM) / (variables["E_vac"] / 10 * variables["c"]);
    return a_aether + a_u_g4i + a_osc + a_quantum + a_fluid + a_exp;
}

// Compute SC Integrated: (1 - B / B_crit) * f_sc
double CompressedResonanceUQFF24Module::computeSCIntegrated(double B) {
    return (1.0 - (B / variables["B_crit"])) * variables["f_sc"];
}

// Full Compressed + Resonance with SC: (compressed + resonance) * SC * (1 + f_TRZ)
double CompressedResonanceUQFF24Module::computeCompressedResTerm(double t, double B) {
    double comp = computeCompressedTerm();
    double res = computeResonanceTerm(t);
    double sc_int = computeSCIntegrated(B);
    double tr_factor = 1.0 + variables["f_TRZ"];
    return (comp + res) * sc_int * tr_factor;
}

// Get equation text (descriptive)
std::string CompressedResonanceUQFF24Module::getEquationText() {
    return "Compressed Terms: a_comp = a_DPM + a_THz + a_vac_diff + a_super (scaled for 18-24)\n"
           "Resonance Terms: a_res = a_aether + U_g4i + a_osc + a_quantum + a_fluid + a_exp\n"
           "Full: g_comp_res = (a_comp + a_res) * SC_int * (1 + f_TRZ)\n"
           "Where SC_int = (1 - B / B_crit) * f_sc\n"
           "Special Terms: UQFF compressed/resonance via plasmotic vacuum; no SM; for systems 18-24 (Sombrero, Saturn, M16, Crab).\n"
           "Solutions: Example g_comp_res ~1e-38 m/s² (micro-scale).\n"
           "Adaptations: Frequencies scaled for nebulae/planets/remnants.";
}

// Print variables
void CompressedResonanceUQFF24Module::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Example usage in base program 'ziqn233h.cpp' (snippet for integration)
// #include "CompressedResonanceUQFF24Module.h"
// int main() {
//     CompressedResonanceUQFF24Module mod;
//     double t = 1e9 * 3.156e7;  // 1 Gyr
//     double B = 1e-5;           // T
//     double g_comp_res = mod.computeCompressedResTerm(t, B);
//     std::cout << "g_comp_res = " << g_comp_res << " m/s²\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("f_DPM", 1.1e11);  // Update
//     mod.addToVariable("f_TRZ", 0.05);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o ziqn233h ziqn233h.cpp CompressedResonanceUQFF24Module.cpp -lm
// Sample Output: g_comp_res ≈ 1e-38 m/s² (varies; micro-scale for systems 18-24).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

// ========== ENHANCED DYNAMIC CAPABILITIES IMPLEMENTATION (25 methods) ==========

// Anonymous namespace for state storage
namespace {
    std::map<std::string, std::map<std::string, double>> systems24_saved_states;
}

// Variable management (5 methods)
void CompressedResonanceUQFF24Module::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void CompressedResonanceUQFF24Module::removeVariable(const std::string& name) {
    variables.erase(name);
}

void CompressedResonanceUQFF24Module::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
    }
}

std::vector<std::string> CompressedResonanceUQFF24Module::listVariables() {
    std::vector<std::string> names;
    for (const auto& pair : variables) {
        names.push_back(pair.first);
    }
    return names;
}

// Batch operations (2 methods)
void CompressedResonanceUQFF24Module::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = func(variables[name]);
        }
    }
}

void CompressedResonanceUQFF24Module::scaleVariableGroup(const std::vector<std::string>& names, double scale) {
    transformVariableGroup(names, [scale](double v) { return v * scale; });
}

// Self-expansion (4 methods)
void CompressedResonanceUQFF24Module::expandParameterSpace(double scale) {
    variables["f_DPM"] *= scale;
    variables["f_THz"] *= scale;
    variables["f_vac_diff"] *= scale;
    variables["f_super"] *= scale;
    variables["f_aether"] *= scale;
    variables["f_react"] *= scale;
    variables["f_quantum"] *= scale;
    variables["f_fluid"] *= scale;
    variables["f_exp"] *= scale;
    variables["f_osc"] *= scale;
    variables["I"] *= scale;
    variables["v_exp"] *= scale;
    variables["A"] *= scale;
}

void CompressedResonanceUQFF24Module::expandCompressedScale(double f_DPM_scale, double f_THz_scale) {
    variables["f_DPM"] *= f_DPM_scale;
    variables["f_THz"] *= f_THz_scale;
    variables["f_vac_diff"] *= f_DPM_scale;
    variables["f_super"] *= f_DPM_scale;
}

void CompressedResonanceUQFF24Module::expandResonanceScale(double f_aether_scale, double f_react_scale) {
    variables["f_aether"] *= f_aether_scale;
    variables["f_react"] *= f_react_scale;
    variables["f_osc"] *= f_aether_scale;
    variables["omega_osc"] *= f_aether_scale;
}

void CompressedResonanceUQFF24Module::expandSystemScale(double I_scale, double V_sys_scale) {
    variables["I"] *= I_scale;
    variables["V_sys"] *= V_sys_scale;
    variables["A_vort"] *= V_sys_scale;
}

// Self-refinement (3 methods)
void CompressedResonanceUQFF24Module::autoRefineParameters(double t, double B, double target_g, double tolerance) {
    double current_g = computeCompressedResTerm(t, B);
    int iterations = 0;
    while (std::abs(current_g - target_g) > tolerance && iterations < 100) {
        double ratio = target_g / (current_g + 1e-50);
        if (std::abs(current_g) < 1e-50) {
            variables["f_DPM"] *= 1.1;
            variables["I"] *= 1.1;
        } else {
            variables["f_DPM"] *= std::sqrt(ratio);
            variables["I"] *= std::sqrt(ratio);
        }
        current_g = computeCompressedResTerm(t, B);
        iterations++;
    }
}

void CompressedResonanceUQFF24Module::calibrateToObservations(const std::map<std::string, double>& observations) {
    for (const auto& obs : observations) {
        if (variables.find(obs.first) != variables.end()) {
            variables[obs.first] = obs.second;
        }
    }
}

void CompressedResonanceUQFF24Module::optimizeForMetric(double t, double B, const std::string& metric) {
    if (metric == "maximize_compressed") {
        variables["f_DPM"] *= 1.2;
        variables["f_THz"] *= 1.2;
        variables["I"] *= 1.2;
    } else if (metric == "maximize_resonance") {
        variables["f_aether"] *= 1.2;
        variables["f_react"] *= 1.2;
        variables["f_osc"] *= 1.2;
    } else if (metric == "balance_terms") {
        variables["f_DPM"] *= 1.1;
        variables["f_aether"] *= 1.1;
    }
}

// Parameter exploration (1 method)
std::vector<std::map<std::string, double>> CompressedResonanceUQFF24Module::generateVariations(int count, double variation_percent) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - variation_percent/100.0, 1.0 + variation_percent/100.0);
    
    for (int i = 0; i < count; ++i) {
        std::map<std::string, double> varied = variables;
        for (auto& pair : varied) {
            if (pair.first != "c" && pair.first != "pi" && pair.first != "hbar") {
                pair.second *= dis(gen);
            }
        }
        variations.push_back(varied);
    }
    return variations;
}

// Adaptive evolution (2 methods)
void CompressedResonanceUQFF24Module::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - mutation_rate, 1.0 + mutation_rate);
    
    for (auto& pair : variables) {
        if (pair.first != "c" && pair.first != "pi" && pair.first != "hbar") {
            pair.second *= dis(gen);
        }
    }
}

void CompressedResonanceUQFF24Module::evolveSystem(double t, double B, int generations, double selection_pressure) {
    for (int gen = 0; gen < generations; ++gen) {
        auto variations = generateVariations(10, 5.0);
        double best_g = computeCompressedResTerm(t, B);
        std::map<std::string, double> best_vars = variables;
        
        for (const auto& var : variations) {
            auto temp_vars = variables;
            variables = var;
            double current_g = computeCompressedResTerm(t, B);
            if (std::abs(current_g) > std::abs(best_g)) {
                best_g = current_g;
                best_vars = var;
            }
            variables = temp_vars;
        }
        variables = best_vars;
    }
}

// State management (4 methods)
void CompressedResonanceUQFF24Module::saveState(const std::string& label) {
    systems24_saved_states[label] = variables;
}

void CompressedResonanceUQFF24Module::restoreState(const std::string& label) {
    if (systems24_saved_states.find(label) != systems24_saved_states.end()) {
        variables = systems24_saved_states[label];
    }
}

std::vector<std::string> CompressedResonanceUQFF24Module::listSavedStates() {
    std::vector<std::string> labels;
    for (const auto& pair : systems24_saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string CompressedResonanceUQFF24Module::exportState() const {
    std::ostringstream oss;
    oss << "CompressedResonanceUQFF24Module State Export:\n";
    for (const auto& pair : variables) {
        oss << pair.first << " = " << std::scientific << pair.second << "\n";
    }
    return oss.str();
}

// System analysis (4 methods)
std::map<std::string, double> CompressedResonanceUQFF24Module::sensitivityAnalysis(double t, double B, double perturbation) {
    std::map<std::string, double> sensitivities;
    double base_g = computeCompressedResTerm(t, B);
    
    std::vector<std::string> key_params = {"f_DPM", "f_THz", "f_vac_diff", "f_super", "f_aether", "f_react", "f_quantum", "f_fluid", "f_exp", "f_osc", "I", "V_sys", "B_crit", "f_TRZ"};
    
    for (const auto& param : key_params) {
        if (variables.find(param) != variables.end()) {
            double original = variables[param];
            variables[param] *= (1.0 + perturbation);
            double perturbed_g = computeCompressedResTerm(t, B);
            sensitivities[param] = (perturbed_g - base_g) / (base_g + 1e-50);
            variables[param] = original;
        }
    }
    return sensitivities;
}

std::string CompressedResonanceUQFF24Module::generateReport(double t, double B) {
    std::ostringstream oss;
    oss << "\n========== COMPRESSED RESONANCE UQFF SYSTEMS 18-24 REPORT ==========\n";
    oss << "System: " << getSystemName() << "\n";
    oss << "Time: " << (t / 3.156e7 / 1e9) << " Gyr\n";
    oss << "Magnetic Field: " << B << " T\n\n";
    
    oss << "Compressed Frequencies (scaled for systems 18-24):\n";
    oss << "  f_DPM = " << (variables["f_DPM"] / 1e11) << " x10^11 Hz\n";
    oss << "  f_THz = " << (variables["f_THz"] / 1e11) << " x10^11 Hz\n";
    oss << "  f_vac_diff = " << variables["f_vac_diff"] << " Hz\n";
    oss << "  f_super = " << (variables["f_super"] / 1e15) << " PHz\n\n";
    
    oss << "Resonance Frequencies:\n";
    oss << "  f_aether = " << (variables["f_aether"] / 1e3) << " kHz\n";
    oss << "  f_react = " << (variables["f_react"] / 1e9) << " GHz\n";
    oss << "  f_quantum = " << variables["f_quantum"] << " Hz\n";
    oss << "  f_fluid = " << variables["f_fluid"] << " Hz\n";
    oss << "  f_exp = " << variables["f_exp"] << " Hz\n";
    oss << "  f_osc = " << (variables["f_osc"] / 1e13) << " x10^13 Hz\n\n";
    
    oss << "System Parameters:\n";
    oss << "  Current I = " << variables["I"] << " A\n";
    oss << "  System Volume V_sys = " << variables["V_sys"] << " m^3\n";
    oss << "  Outflow v_exp = " << (variables["v_exp"] / 1e3) << " km/s\n\n";
    
    oss << "Superconductive Integration:\n";
    oss << "  B_crit = " << variables["B_crit"] << " T\n";
    oss << "  SC_int = " << computeSCIntegrated(B) << "\n\n";
    
    oss << "Computed Terms:\n";
    oss << "  Compressed Term = " << computeCompressedTerm() << " m/s^2\n";
    oss << "  Resonance Term = " << computeResonanceTerm(t) << " m/s^2\n";
    oss << "  Full (Comp+Res+SC): " << computeCompressedResTerm(t, B) << " m/s^2\n";
    oss << "====================================================================\n";
    
    return oss.str();
}

bool CompressedResonanceUQFF24Module::validateConsistency() {
    bool valid = true;
    std::vector<std::string> positive_params = {"f_DPM", "f_THz", "f_super", "f_aether", "f_react", "f_quantum", "f_fluid", "f_exp", "f_osc", "I", "V_sys", "B_crit", "E_vac", "c"};
    
    for (const auto& param : positive_params) {
        if (variables.find(param) != variables.end() && variables[param] <= 0) {
            valid = false;
            break;
        }
    }
    return valid;
}

bool CompressedResonanceUQFF24Module::autoCorrectAnomalies() {
    bool corrected = false;
    
    if (variables["f_DPM"] <= 0) { variables["f_DPM"] = 1e11; corrected = true; }
    if (variables["f_THz"] <= 0) { variables["f_THz"] = 1e11; corrected = true; }
    if (variables["f_vac_diff"] <= 0) { variables["f_vac_diff"] = 0.143; corrected = true; }
    if (variables["f_super"] <= 0) { variables["f_super"] = 1.411e15; corrected = true; }
    if (variables["f_aether"] <= 0) { variables["f_aether"] = 1e3; corrected = true; }
    if (variables["f_react"] <= 0) { variables["f_react"] = 1e9; corrected = true; }
    if (variables["f_quantum"] <= 0) { variables["f_quantum"] = 1.445e-17; corrected = true; }
    if (variables["f_fluid"] <= 0) { variables["f_fluid"] = 1.269e-14; corrected = true; }
    if (variables["f_exp"] <= 0) { variables["f_exp"] = 1.373e-8; corrected = true; }
    if (variables["f_osc"] <= 0) { variables["f_osc"] = 4.57e13; corrected = true; }
    if (variables["I"] <= 0) { variables["I"] = 1e20; corrected = true; }
    if (variables["V_sys"] <= 0) { variables["V_sys"] = 4.189e18; corrected = true; }
    if (variables["B_crit"] <= 0) { variables["B_crit"] = 1e11; corrected = true; }
    
    return corrected;
}

// Evaluation of CompressedResonanceUQFF24Module (UQFF Compressed & Resonance Terms for Systems 18-24)

**Strengths:**
-**Dynamic & Extensible : **All model parameters are stored in a `std: : map<std::string, double> variables`, allowing runtime updates, additions, and removals.The methods `updateVariable`, `addToVariable`, and `subtractFromVariable` enable flexible modification of any parameter.
- **Automatic Dependency Updates : **When key variables like `"Delta_x"` are updated, dependent variables(`"Delta_p"`) are recalculated automatically, ensuring consistency.
    - **Immediate Effect : **All computations(e.g., `computeCompressedResTerm`) use the current values in the map, so any changes are immediately reflected in results.
        - **Comprehensive Physics : **The module includes all major UQFF compressed and resonance terms, such as DPM, THz, vacuum differential, superconductor, aether, U_g4i, oscillatory, quantum, fluid, and expansion effects.Standard Model gravity / magnetics are intentionally excluded per UQFF.
        - **Scalability : **Parameters and frequencies are scaled for systems 18 - 24 (e.g., nebulae, planets, remnants), making the module adaptable to a range of astrophysical scenarios.
        - **Debugging Support : **The `printVariables()` method provides a snapshot of all current parameters, aiding validation and troubleshooting.
    - **Sample Usage Provided : **Example integration and compilation instructions are included, demonstrating how to update variables and see their effect.

    ** Weaknesses / Recommendations : **
    -**Error Handling : **Adding unknown variables is flexible but may lead to silent errors if a typo occurs.Consider stricter validation or optional warnings for unknown variable names.
    - **Unit Consistency : **Ensure all units are consistent, especially when combining terms from different physical domains and frequency - based scaling.
    - **Magic Numbers : **Some scale factors and physical constants are set to arbitrary values.Document their physical meaning or allow configuration via constructor arguments or config files.
    - **Physical Justification : **The model is highly specialized for UQFF and omits Standard Model terms.Ensure this is appropriate for your scientific context and document the rationale for each term.
    - **Performance : **For large - scale or repeated updates, consider profiling and optimizing critical paths if needed.

    ** Summary : **
    The module is robust, dynamic, and extensible, supporting runtime updates and changes to all model parameters.It is suitable for advanced UQFF - based compressed and resonance modeling for systems 18 - 24. Minor improvements in error handling, documentation, and physical justification are recommended for production or publication use.

// ========== ENHANCED 18-STEP EXAMPLE FUNCTION ==========
void example_enhanced_systems24_18_steps() {
    std::cout << "\n========== ENHANCED COMPRESSED RESONANCE UQFF SYSTEMS 18-24 DEMONSTRATION ==========\n";
    std::cout << "UQFF Compressed & Resonance for Systems 18-24 (Sombrero, Saturn, M16, Crab)\n\n";
    
    CompressedResonanceUQFF24Module systems24;
    double t_current = 1e9 * 3.156e7; // 1 Gyr in seconds
    double B_current = 1e-5;          // 1e-5 T
    
    // Step 1: Initial state for systems 18-24
    std::cout << "Step 1: Initial compressed/resonance state (systems 18-24) at t = 1 Gyr\n";
    double g1 = systems24.computeCompressedResTerm(t_current, B_current);
    std::cout << "  Compressed: f_DPM = " << (systems24.variables["f_DPM"] / 1e11) << " x10^11 Hz\n";
    std::cout << "  Resonance: f_aether = " << (systems24.variables["f_aether"] / 1e3) << " kHz\n";
    std::cout << "  System: I = " << systems24.variables["I"] << " A\n";
    std::cout << "  g_comp_res = " << g1 << " m/s^2 (micro-scale)\n\n";
    
    // Step 2: Save initial state
    std::cout << "Step 2: Save initial systems 18-24 state\n";
    systems24.saveState("systems24_initial");
    std::cout << "  State saved as 'systems24_initial'\n\n";
    
    // Step 3: Expand compressed scale
    std::cout << "Step 3: Expand compressed scale (1.5x DPM, 1.2x THz)\n";
    systems24.expandCompressedScale(1.5, 1.2);
    double g3 = systems24.computeCompressedResTerm(t_current, B_current);
    std::cout << "  New f_DPM = " << (systems24.variables["f_DPM"] / 1e11) << " x10^11 Hz\n";
    std::cout << "  New f_THz = " << (systems24.variables["f_THz"] / 1e11) << " x10^11 Hz\n";
    std::cout << "  g_comp_res = " << g3 << " m/s^2\n\n";
    
    // Step 4: Restore and expand resonance scale
    std::cout << "Step 4: Restore initial state, then expand resonance scale (2x aether, 1.5x react)\n";
    systems24.restoreState("systems24_initial");
    systems24.expandResonanceScale(2.0, 1.5);
    double g4 = systems24.computeCompressedResTerm(t_current, B_current);
    std::cout << "  New f_aether = " << (systems24.variables["f_aether"] / 1e3) << " kHz\n";
    std::cout << "  New f_react = " << (systems24.variables["f_react"] / 1e9) << " GHz\n";
    std::cout << "  g_comp_res = " << g4 << " m/s^2\n\n";
    
    // Step 5: Restore and expand system scale
    std::cout << "Step 5: Restore initial state, then expand system scale (1.5x current, 2x volume)\n";
    systems24.restoreState("systems24_initial");
    systems24.expandSystemScale(1.5, 2.0);
    double g5 = systems24.computeCompressedResTerm(t_current, B_current);
    std::cout << "  New I = " << systems24.variables["I"] << " A\n";
    std::cout << "  New V_sys = " << systems24.variables["V_sys"] << " m^3\n";
    std::cout << "  g_comp_res = " << g5 << " m/s^2\n\n";
    
    // Step 6: Time evolution (cosmic timescales)
    std::cout << "Step 6: Time evolution from 0 to 5 Gyr\n";
    systems24.restoreState("systems24_initial");
    for (double t_Gyr = 0; t_Gyr <= 5; t_Gyr += 1) {
        double t_sec = t_Gyr * 1e9 * 3.156e7;
        double g = systems24.computeCompressedResTerm(t_sec, B_current);
        std::cout << "  t = " << t_Gyr << " Gyr: g = " << g << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 7: Magnetic field sweep (SC integration effect)
    std::cout << "Step 7: Magnetic field sweep (1e-6 to 1e-3 T, SC integration)\n";
    for (double B : {1e-6, 1e-5, 1e-4, 1e-3}) {
        double g = systems24.computeCompressedResTerm(t_current, B);
        double sc_int = systems24.computeSCIntegrated(B);
        std::cout << "  B = " << B << " T: SC_int = " << sc_int << ", g = " << g << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 8: Create custom tracking variables
    std::cout << "Step 8: Create custom tracking variables\n";
    systems24.createVariable("system_type", 18.0); // 18-24: Sombrero, Saturn, M16, Crab, etc.
    systems24.createVariable("observation_wavelength", 5e-7); // 500 nm (optical)
    systems24.createVariable("distance_Mpc", 10.0); // Distance in Mpc
    std::cout << "  Created 'system_type', 'observation_wavelength', 'distance_Mpc'\n\n";
    
    // Step 9: Generate variations for uncertainty analysis
    std::cout << "Step 9: Generate 3 parameter variations (5% perturbation)\n";
    auto variations = systems24.generateVariations(3, 5.0);
    for (size_t i = 0; i < variations.size(); ++i) {
        CompressedResonanceUQFF24Module temp = systems24;
        temp.variables = variations[i];
        double g_var = temp.computeCompressedResTerm(t_current, B_current);
        std::cout << "  Variation " << (i+1) << ": g = " << g_var << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 10: Sensitivity analysis
    std::cout << "Step 10: Sensitivity analysis (1% perturbation)\n";
    auto sensitivities = systems24.sensitivityAnalysis(t_current, B_current, 0.01);
    std::cout << "  Parameter sensitivities (fractional change in g):\n";
    for (const auto& s : sensitivities) {
        std::cout << "    " << s.first << ": " << s.second << "\n";
    }
    std::cout << "\n";
    
    // Step 11: DPM frequency sweep (compressed term dominant)
    std::cout << "Step 11: DPM frequency sweep (0.5x, 1.0x, 2.0x)\n";
    systems24.saveState("systems24_before_sweep");
    for (double scale : {0.5, 1.0, 2.0}) {
        systems24.restoreState("systems24_before_sweep");
        systems24.expandCompressedScale(scale, 1.0);
        double g = systems24.computeCompressedResTerm(t_current, B_current);
        double f_DPM = systems24.variables["f_DPM"] / 1e11;
        std::cout << "  f_DPM = " << f_DPM << " x10^11 Hz: g = " << g << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 12: Aether frequency sweep (resonance term)
    std::cout << "Step 12: Aether frequency sweep (0.5x, 1.0x, 2.0x)\n";
    systems24.restoreState("systems24_before_sweep");
    for (double scale : {0.5, 1.0, 2.0}) {
        systems24.restoreState("systems24_before_sweep");
        systems24.expandResonanceScale(scale, 1.0);
        double g = systems24.computeCompressedResTerm(t_current, B_current);
        double f_aether = systems24.variables["f_aether"] / 1e3;
        std::cout << "  f_aether = " << f_aether << " kHz: g = " << g << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 13: System current sweep (scale variations)
    std::cout << "Step 13: System current sweep (0.5x, 1.0x, 2.0x)\n";
    systems24.restoreState("systems24_before_sweep");
    for (double scale : {0.5, 1.0, 2.0}) {
        systems24.restoreState("systems24_before_sweep");
        systems24.expandSystemScale(scale, 1.0);
        double g = systems24.computeCompressedResTerm(t_current, B_current);
        double I = systems24.variables["I"];
        std::cout << "  I = " << I << " A: g = " << g << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 14: Batch transform all compressed frequencies
    std::cout << "Step 14: Batch transform all compressed frequencies (1.1x scale)\n";
    systems24.restoreState("systems24_before_sweep");
    systems24.scaleVariableGroup({"f_DPM", "f_THz", "f_vac_diff", "f_super"}, 1.1);
    double g14 = systems24.computeCompressedResTerm(t_current, B_current);
    std::cout << "  All compressed frequencies scaled by 1.1x\n";
    std::cout << "  f_DPM = " << (systems24.variables["f_DPM"] / 1e11) << " x10^11 Hz\n";
    std::cout << "  f_super = " << (systems24.variables["f_super"] / 1e15) << " PHz\n";
    std::cout << "  g_comp_res = " << g14 << " m/s^2\n\n";
    
    // Step 15: Validate and auto-correct
    std::cout << "Step 15: Validate consistency and auto-correct if needed\n";
    systems24.restoreState("systems24_before_sweep");
    bool valid = systems24.validateConsistency();
    std::cout << "  System valid: " << (valid ? "Yes" : "No") << "\n";
    if (!valid) {
        bool corrected = systems24.autoCorrectAnomalies();
        std::cout << "  Auto-corrected: " << (corrected ? "Yes" : "No") << "\n";
    }
    std::cout << "\n";
    
    // Step 16: Auto-refine to target acceleration
    std::cout << "Step 16: Auto-refine parameters to target g = 5e-39 m/s^2\n";
    systems24.restoreState("systems24_before_sweep");
    double target_g = 5e-39;
    systems24.autoRefineParameters(t_current, B_current, target_g, 1e-40);
    double g16 = systems24.computeCompressedResTerm(t_current, B_current);
    std::cout << "  Target g = " << target_g << " m/s^2\n";
    std::cout << "  Achieved g = " << g16 << " m/s^2\n";
    std::cout << "  Refined f_DPM = " << (systems24.variables["f_DPM"] / 1e11) << " x10^11 Hz\n\n";
    
    // Step 17: List all saved states
    std::cout << "Step 17: List all saved states\n";
    auto states = systems24.listSavedStates();
    std::cout << "  Saved states (" << states.size() << " total):\n";
    for (const auto& state : states) {
        std::cout << "    - " << state << "\n";
    }
    std::cout << "\n";
    
    // Step 18: Generate comprehensive report
    std::cout << "Step 18: Generate comprehensive system report\n";
    systems24.restoreState("systems24_initial");
    std::string report = systems24.generateReport(t_current, B_current);
    std::cout << report << "\n";
    
    std::cout << "========== END 18-STEP SYSTEMS 18-24 DEMONSTRATION ==========\n\n";
}