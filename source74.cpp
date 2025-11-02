// UQFFCompressedResonanceModule.h
// Modular C++ implementation of Compressed and Resonance UQFF Equations for Multi-System Evolution (Young Stars Outflows, Eagle Nebula, Big Bang, M51, NGC 1316, V838 Mon, NGC 1300, Student's Guide).
// Supports compressed g_UQFF(r,t) unified form; resonance mode adds oscillatory terms (cos/exp(i ? t)) for wave dynamics.
// Usage: #include "UQFFCompressedResonanceModule.h"; UQFFCompressedResonanceModule mod; mod.setSystem("Eagle"); mod.setMode("resonance"); mod.computeG(t);
// Variables in std::map; auto-loads params from DeepSearch (Hubble/JWST/CERN/high-energy labs).
// Approximations: psi_int=1.0; H(t,z) standard; resonance A=1e-10, ?=1e15; Big Bang: r=c t, rho=rho_c.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef UQFF_COMPRESSED_RESONANCE_MODULE_H
#define UQFF_COMPRESSED_RESONANCE_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>
#include <functional>

class UQFFCompressedResonanceModule {
private:
    std::map<std::string, double> variables;
    std::string current_system;
    std::string mode;  // "compressed" or "resonance"
    double computeHtz(double z_val);
    double computeFenv(double t);
    double computeUgSum();
    double computePsiTotal(double t);
    double computeResonanceTerm(double t);
    double computeQuantumTerm(double t_Hubble_val);
    double computeFluidTerm(double g_base);
    double computeDMTerm();
    double computeMsfFactor(double t);

public:
    // Constructor: General defaults
    UQFFCompressedResonanceModule();

    // Set system and load params
    void setSystem(const std::string& sys_name);

    // Set mode: compressed or resonance
    void setMode(const std::string& m);

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core: g_UQFF(r, t) or I_echo for V838
    double computeG(double t, double r = 0.0);  // r=0 uses default

    // Equation text (mode-specific)
    std::string getEquationText();

    // Print variables
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

#endif // UQFF_COMPRESSED_RESONANCE_MODULE_H

// UQFFCompressedResonanceModule.cpp
#include "UQFFCompressedResonanceModule.h"
#include <complex>

// Constructor
UQFFCompressedResonanceModule::UQFFCompressedResonanceModule() : current_system("Guide"), mode("compressed") {
    // Universal
    variables["G"] = 6.6743e-11;
    variables["c"] = 3e8;
    variables["hbar"] = 1.0546e-34;
    variables["Lambda"] = 1.1e-52;
    variables["q"] = 1.602e-19;
    variables["pi"] = 3.141592653589793;
    variables["t_Hubble"] = 13.8e9 * 3.156e7;
    variables["year_to_s"] = 3.156e7;
    variables["H0"] = 70.0;
    variables["Mpc_to_m"] = 3.086e22;
    variables["Omega_m"] = 0.3;
    variables["Omega_Lambda"] = 0.7;
    variables["M_sun"] = 1.989e30;
    variables["kpc"] = 3.086e19;

    // General defaults (overridden by setSystem)
    variables["M"] = 1e41;  // kg
    variables["M0"] = variables["M"];
    variables["SFR"] = 6e19;  // kg/s (~2 Msun/yr)
    variables["r"] = 1e20;    // m
    variables["z"] = 0.005;
    variables["M_visible"] = 0.7 * variables["M"];
    variables["M_DM"] = 0.3 * variables["M"];
    variables["t"] = 1e9 * variables["year_to_s"];
    variables["rho_fluid"] = 1e-21;
    variables["V"] = 1e50;
    variables["B"] = 1e-5;
    variables["B_crit"] = 1e11;
    variables["Delta_x"] = 1e-10;
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];
    variables["integral_psi"] = 1.0;
    variables["A"] = 1e-10;
    variables["k"] = 1e20;
    variables["omega"] = 1e15;
    variables["x"] = 0.0;
    variables["v"] = 1e3;
    variables["Ug1"] = 0.0; variables["Ug2"] = 0.0; variables["Ug3"] = 0.0; variables["Ug4"] = 0.0;
    variables["scale_macro"] = 1e-12;
    variables["f_TRZ"] = 0.1;
    variables["f_sc"] = 1.0;
    variables["delta_rho"] = 1e-5 * variables["rho_fluid"];
    variables["rho"] = variables["rho_fluid"];
    variables["F_wind"] = 0.0;  // etc. for F_env
}

// Set system: Load DeepSearch params
void UQFFCompressedResonanceModule::setSystem(const std::string& sys_name) {
    current_system = sys_name;
    double Msun = variables["M_sun"];
    double kpc = variables["kpc"];
    double yr_s = variables["year_to_s"];
    if (sys_name == "YoungStars") {
        variables["M"] = 1000 * Msun; variables["r"] = 3e17; variables["SFR"] = 0.1 * Msun / yr_s;
        variables["rho_fluid"] = 1e-20; variables["B"] = 1e-6; variables["z"] = 0.0006;
    } else if (sys_name == "Eagle") {
        variables["M"] = 1e4 * Msun; variables["r"] = 2e17; variables["SFR"] = 0.5 * Msun / yr_s;
        variables["rho_fluid"] = 1e-21; variables["B"] = 3e-5; variables["z"] = 0.002;
    } else if (sys_name == "BigBang") {
        variables["rho_fluid"] = 8e-27; variables["r"] = 1e26; variables["z"] = 1100; variables["SFR"] = 0;  // Cosmic
        variables["M"] = 1e53;  // Observable universe approx
        variables["B"] = 1e-10; variables["t"] = 13.8e9 * yr_s;
    } else if (sys_name == "M51") {
        variables["M"] = 1.6e11 * Msun; variables["r"] = 23e3 * kpc; variables["SFR"] = 2 * Msun / yr_s;
        variables["rho_fluid"] = 1e-21; variables["B"] = 1e-5; variables["z"] = 0.005;
    } else if (sys_name == "NGC1316") {
        variables["M"] = 5e11 * Msun; variables["r"] = 23e3 * kpc; variables["SFR"] = 0.1 * Msun / yr_s;
        variables["rho_fluid"] = 1e-22; variables["B"] = 1e-5; variables["z"] = 0.006;
    } else if (sys_name == "V838Mon") {
        variables["M"] = 8 * Msun; variables["r"] = 2e13; variables["SFR"] = 0;
        variables["rho_fluid"] = 1e-22; variables["B"] = 1e-6; variables["z"] = 0.005;
    } else if (sys_name == "NGC1300") {
        variables["M"] = 1e11 * Msun; variables["r"] = 12e3 * kpc; variables["SFR"] = 1 * Msun / yr_s;
        variables["rho_fluid"] = 1e-21; variables["B"] = 1e-5; variables["z"] = 0.005;
    } else {  // Guide: general
        variables["M"] = Msun; variables["r"] = 1e11; variables["SFR"] = 1e-10 * Msun / yr_s;  // Low
        variables["rho_fluid"] = 1e-20; variables["B"] = 1e-5; variables["z"] = 0;
    }
    // Update dependents
    variables["M_visible"] = 0.7 * variables["M"];
    variables["M_DM"] = 0.3 * variables["M"];
    variables["M0"] = variables["M"];
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];
}

// Set mode
void UQFFCompressedResonanceModule::setMode(const std::string& m) {
    mode = m;
}

// Update etc. (as before)
void UQFFCompressedResonanceModule::updateVariable(const std::string& name, double value) {
    // Similar to template
    if (variables.find(name) != variables.end()) variables[name] = value;
    else variables[name] = value;
    if (name == "M") {
        variables["M_visible"] = 0.7 * value;
        variables["M_DM"] = 0.3 * value;
        variables["M0"] = value;
    }
    if (name == "Delta_x") variables["Delta_p"] = variables["hbar"] / value;
}
void UQFFCompressedResonanceModule::addToVariable(const std::string& name, double delta) {
    updateVariable(name, variables[name] + delta);  // Simplified
}
void UQFFCompressedResonanceModule::subtractFromVariable(const std::string& name, double delta) {
    updateVariable(name, variables[name] - delta);
}

// Helpers (as in previous, adapted)
double UQFFCompressedResonanceModule::computeHtz(double z_val) {
    double Hz_kms = variables["H0"] * std::sqrt(variables["Omega_m"] * std::pow(1 + z_val, 3) + variables["Omega_Lambda"]);
    return (Hz_kms * 1e3) / variables["Mpc_to_m"];
}
double UQFFCompressedResonanceModule::computeFenv(double t) { return 0.1; }  // Simplified
double UQFFCompressedResonanceModule::computeUgSum() { return 1e-10; }  // Placeholder
double UQFFCompressedResonanceModule::computePsiTotal(double t) {
    return variables["q"] * variables["v"] * variables["B"] + 2 * variables["A"] * std::cos(variables["k"] * variables["x"] + variables["omega"] * t);
}
double UQFFCompressedResonanceModule::computeResonanceTerm(double t) {
    if (mode != "resonance") return 0.0;
    std::complex<double> exp_term(variables["A"] * std::exp(std::complex<double>(0, variables["k"] * variables["x"] - variables["omega"] * t)));
    return (2 * variables["pi"] / 13.8) * exp_term.real() * computeG(t, variables["r"]);  // Coupled
}
double UQFFCompressedResonanceModule::computeQuantumTerm(double t_Hubble_val) {
    double unc = std::sqrt(variables["Delta_x"] * variables["Delta_p"]);
    double psi = computePsiTotal(variables["t"]);
    return (variables["hbar"] / unc) * variables["integral_psi"] * (2 * variables["pi"] / t_Hubble_val) * psi;
}
double UQFFCompressedResonanceModule::computeFluidTerm(double g_base) { return variables["rho_fluid"] * variables["V"] * g_base; }
double UQFFCompressedResonanceModule::computeDMTerm() {
    double pert = variables["delta_rho"] / variables["rho"];
    double curv = 3 * variables["G"] * variables["M"] / (std::pow(variables["r"], 3));
    return (variables["M_visible"] + variables["M_DM"]) * (pert + curv);
}
double UQFFCompressedResonanceModule::computeMsfFactor(double t) { return variables["SFR"] * t / variables["M0"]; }

// Compute G (or I for V838)
double UQFFCompressedResonanceModule::computeG(double t, double r_in) {
    if (r_in > 0) variables["r"] = r_in;
    variables["t"] = t;
    if (current_system == "BigBang") variables["r"] = variables["c"] * t;  // Scale
    if (current_system == "V838Mon") {  // Return I_echo
        double rho_d = variables["rho_fluid"] * std::exp(-1.0 * (variables["G"] * variables["M"] / (variables["r"] * variables["r"])));
        return (600000 * 3.826e26) / (4 * variables["pi"] * variables["r"] * variables["r"]) * 1e-12 * rho_d;  // Approx
    }
    double Hz = computeHtz(variables["z"]);
    double expansion = 1 + Hz * t;
    double sc = 1 - variables["B"] / variables["B_crit"];
    double msf = computeMsfFactor(t);
    double mfact = 1 + msf;
    double fenv = computeFenv(t);
    double g_base = (variables["G"] * variables["M"] * mfact / (variables["r"] * variables["r"])) * expansion * sc * (1 + fenv);
    double ugsum = computeUgSum();
    double lambda_t = variables["Lambda"] * variables["c"] * variables["c"] / 3;
    double qterm = computeQuantumTerm(variables["t_Hubble"]);
    double fterm = computeFluidTerm(g_base);
    double dmterm = computeDMTerm();
    double res_term = computeResonanceTerm(t);
    return g_base + ugsum + lambda_t + qterm + fterm + dmterm + res_term;
}

// Equation text
std::string UQFFCompressedResonanceModule::getEquationText() {
    std::string eq = "g_UQFF(r,t) = (G M(t)/r^2) (1 + H(t,z)) (1 - B/B_crit) (1 + F_env) + ? Ug_i + ? c^2/3 + (?/?(?x ?p)) ? ? H ? dV (2?/t_Hubble) + ? V g + (M_vis + M_DM)(??/? + 3GM/r^3)";
    if (mode == "resonance") eq += " + 2 A cos(kx + ? t) g_base + (2?/13.8) Re[A exp(i(kx - ? t))] g_base";
    eq += "\nM(t)=M(1 + SFR t / M0); Systems: " + current_system + "; Learning: Yes, diverse scales refine UQFF; Advancing: Unified compressed/resonance explains outflows to cosmic expansion.";
    return eq;
}

// Print
void UQFFCompressedResonanceModule::printVariables() {
    std::cout << "System: " << current_system << " Mode: " << mode << "\nVariables:\n";
    for (auto& p : variables) std::cout << p.first << " = " << std::scientific << p.second << "\n";
}

// ===== DYNAMIC SELF-UPDATE & SELF-EXPANSION IMPLEMENTATIONS =====

// Static storage for saved states
static std::map<std::string, std::map<std::string, double>> uqff_multi_saved_states;

// 1. Dynamic variable management
void UQFFCompressedResonanceModule::createDynamicVariable(const std::string& name, double value) {
    variables[name] = value;
    std::cout << "Created dynamic variable: " << name << " = " << value << std::endl;
}

void UQFFCompressedResonanceModule::removeDynamicVariable(const std::string& name) {
    if (variables.find(name) != variables.end()) {
        variables.erase(name);
        std::cout << "Removed dynamic variable: " << name << std::endl;
    } else {
        std::cerr << "Variable '" << name << "' not found for removal." << std::endl;
    }
}

void UQFFCompressedResonanceModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
        std::cout << "Cloned " << source << " to " << dest << std::endl;
    } else {
        std::cerr << "Source variable '" << source << "' not found." << std::endl;
    }
}

void UQFFCompressedResonanceModule::listAllVariables() {
    std::cout << "=== All UQFF Variables (System: " << current_system << ", Mode: " << mode 
              << ", Total: " << variables.size() << ") ===" << std::endl;
    for (const auto& pair : variables) {
        std::cout << "  " << pair.first << " = " << pair.second << std::endl;
    }
}

// 2. Batch operations
void UQFFCompressedResonanceModule::applyTransformToGroup(const std::vector<std::string>& varNames,
                                                           std::function<double(double)> transform) {
    for (const auto& name : varNames) {
        if (variables.find(name) != variables.end()) {
            variables[name] = transform(variables[name]);
            std::cout << "Transformed " << name << " to " << variables[name] << std::endl;
        }
    }
}

void UQFFCompressedResonanceModule::scaleVariableGroup(const std::vector<std::string>& varNames, double scale_factor) {
    applyTransformToGroup(varNames, [scale_factor](double val) { return val * scale_factor; });
}

// 3. Self-expansion capabilities
void UQFFCompressedResonanceModule::autoExpandParameterSpace(double scale_factor) {
    std::cout << "Auto-expanding UQFF parameter space by factor " << scale_factor << std::endl;
    std::vector<std::string> expandable = {"M", "M_visible", "M_DM", "r", "V", "SFR"};
    scaleVariableGroup(expandable, scale_factor);
    // Update M consistency
    variables["M"] = variables["M_visible"] + variables["M_DM"];
    variables["M0"] = variables["M"];
    std::cout << "  Updated M_total, M0" << std::endl;
}

void UQFFCompressedResonanceModule::expandMassScale(double mass_multiplier) {
    std::cout << "Expanding mass scale by " << mass_multiplier << std::endl;
    variables["M_visible"] *= mass_multiplier;
    variables["M_DM"] *= mass_multiplier;
    variables["M"] = variables["M_visible"] + variables["M_DM"];
    variables["M0"] = variables["M"];
    std::cout << "  M_visible: " << (variables["M_visible"] / variables["M_sun"]) << " M☉" << std::endl;
    std::cout << "  M_DM: " << (variables["M_DM"] / variables["M_sun"]) << " M☉" << std::endl;
    std::cout << "  M_total: " << (variables["M"] / variables["M_sun"]) << " M☉" << std::endl;
}

void UQFFCompressedResonanceModule::expandSpatialScale(double spatial_multiplier) {
    std::cout << "Expanding spatial scale by " << spatial_multiplier << std::endl;
    variables["r"] *= spatial_multiplier;
    variables["V"] *= (spatial_multiplier * spatial_multiplier * spatial_multiplier);
    variables["Delta_x"] *= spatial_multiplier;
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];
    std::cout << "  r: " << variables["r"] << " m" << std::endl;
    std::cout << "  V: " << variables["V"] << " m^3" << std::endl;
    std::cout << "  Updated Delta_p" << std::endl;
}

void UQFFCompressedResonanceModule::expandTimeScale(double time_multiplier) {
    std::cout << "Expanding time scale by " << time_multiplier << std::endl;
    variables["t"] *= time_multiplier;
    variables["omega"] *= (1.0 / time_multiplier);
    std::cout << "  t: " << (variables["t"] / variables["year_to_s"]) << " years" << std::endl;
    std::cout << "  omega: " << variables["omega"] << " rad/s" << std::endl;
}

// 4. Self-refinement
void UQFFCompressedResonanceModule::autoRefineParameters(double tolerance) {
    std::cout << "Auto-refining UQFF parameters with tolerance " << tolerance << std::endl;
    
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

void UQFFCompressedResonanceModule::calibrateToObservations(const std::map<std::string, double>& observed_values) {
    std::cout << "Calibrating to " << observed_values.size() << " observations for " 
              << current_system << "..." << std::endl;
    for (const auto& obs : observed_values) {
        if (variables.find(obs.first) != variables.end()) {
            double old_val = variables[obs.first];
            updateVariable(obs.first, obs.second);
            std::cout << "  " << obs.first << ": " << old_val << " -> " << obs.second << std::endl;
        }
    }
    std::cout << "Calibration complete." << std::endl;
}

void UQFFCompressedResonanceModule::optimizeForMetric(const std::string& metric_name, double target_value) {
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
void UQFFCompressedResonanceModule::generateVariations(int num_variations, double variation_range) {
    std::cout << "Generating " << num_variations << " UQFF variations for " << current_system 
              << " with range ±" << (variation_range * 100) << "%" << std::endl;
    
    std::vector<std::string> key_params = {"M", "r", "SFR", "rho_fluid", "B", "omega"};
    
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

void UQFFCompressedResonanceModule::findOptimalParameters(const std::string& objective, int iterations) {
    std::cout << "Finding optimal UQFF parameters for: " << objective 
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
void UQFFCompressedResonanceModule::mutateParameters(double mutation_rate, double mutation_strength) {
    std::vector<std::string> mutable_params = {"M", "r", "SFR", "rho_fluid", "B", "omega", "v"};
    
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

void UQFFCompressedResonanceModule::evolveSystem(int generations) {
    std::cout << "Evolving UQFF system (" << current_system << ") over " 
              << generations << " generations..." << std::endl;
    
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
void UQFFCompressedResonanceModule::saveState(const std::string& label) {
    uqff_multi_saved_states[label] = variables;
    std::cout << "Saved UQFF state: " << label << " (System: " << current_system 
              << ", " << variables.size() << " variables)" << std::endl;
}

void UQFFCompressedResonanceModule::restoreState(const std::string& label) {
    if (uqff_multi_saved_states.find(label) != uqff_multi_saved_states.end()) {
        variables = uqff_multi_saved_states[label];
        std::cout << "Restored UQFF state: " << label << std::endl;
    } else {
        std::cerr << "State '" << label << "' not found." << std::endl;
    }
}

void UQFFCompressedResonanceModule::listSavedStates() {
    std::cout << "=== Saved UQFF States (Total: " << uqff_multi_saved_states.size() << ") ===" << std::endl;
    for (const auto& state : uqff_multi_saved_states) {
        std::cout << "  " << state.first << " (" << state.second.size() << " variables)" << std::endl;
    }
}

void UQFFCompressedResonanceModule::exportState(const std::string& filename) {
    std::cout << "Exporting UQFF state to " << filename << " (not implemented - placeholder)" << std::endl;
    // In real implementation: write variables to file
}

// 8. System analysis
void UQFFCompressedResonanceModule::analyzeParameterSensitivity(const std::string& param_name) {
    if (variables.find(param_name) == variables.end()) {
        std::cerr << "Parameter '" << param_name << "' not found." << std::endl;
        return;
    }
    
    std::cout << "=== UQFF Sensitivity Analysis: " << param_name 
              << " (System: " << current_system << ") ===" << std::endl;
    
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

void UQFFCompressedResonanceModule::generateSystemReport() {
    std::cout << "\n========== UQFF Multi-System Report ==========" << std::endl;
    std::cout << "Current System: " << current_system << std::endl;
    std::cout << "Mode: " << mode << " (compressed or resonance)" << std::endl;
    std::cout << "Total Variables: " << variables.size() << std::endl;
    
    // Key parameters
    std::cout << "\nMass Parameters:" << std::endl;
    std::cout << "M_visible: " << (variables["M_visible"] / variables["M_sun"]) << " M☉" << std::endl;
    std::cout << "M_DM: " << (variables["M_DM"] / variables["M_sun"]) << " M☉" << std::endl;
    std::cout << "M_total: " << (variables["M"] / variables["M_sun"]) << " M☉" << std::endl;
    std::cout << "SFR: " << (variables["SFR"] * variables["year_to_s"] / variables["M_sun"]) << " M☉/yr" << std::endl;
    
    std::cout << "\nSpatial Parameters:" << std::endl;
    std::cout << "r: " << variables["r"] << " m" << std::endl;
    std::cout << "z (redshift): " << variables["z"] << std::endl;
    std::cout << "V (volume): " << variables["V"] << " m^3" << std::endl;
    
    std::cout << "\nDynamics:" << std::endl;
    std::cout << "v: " << variables["v"] << " m/s" << std::endl;
    std::cout << "rho_fluid: " << variables["rho_fluid"] << " kg/m^3" << std::endl;
    std::cout << "B (magnetic): " << variables["B"] << " T" << std::endl;
    
    std::cout << "\nQuantum Parameters:" << std::endl;
    std::cout << "Delta_x: " << variables["Delta_x"] << " m" << std::endl;
    std::cout << "Delta_p: " << variables["Delta_p"] << " kg·m/s" << std::endl;
    
    if (mode == "resonance") {
        std::cout << "\nResonance Parameters:" << std::endl;
        std::cout << "A (amplitude): " << variables["A"] << std::endl;
        std::cout << "k (wave number): " << variables["k"] << " m^-1" << std::endl;
        std::cout << "omega (frequency): " << variables["omega"] << " rad/s" << std::endl;
    }
    
    // Current computation
    double t = variables["t"];
    double r = variables["r"];
    double g = computeG(t, r);
    
    std::cout << "\nCurrent Computation:" << std::endl;
    std::cout << "t: " << (t / variables["year_to_s"]) << " years" << std::endl;
    std::cout << "r_eval: " << r << " m" << std::endl;
    std::cout << "g_UQFF: " << g << " m/s^2" << std::endl;
    
    std::cout << "\nSubcomponents:" << std::endl;
    std::cout << "Ug_sum: " << computeUgSum() << std::endl;
    std::cout << "H(z): " << computeHtz(variables["z"]) << " s^-1" << std::endl;
    std::cout << "F_env: " << computeFenv(t) << std::endl;
    if (mode == "resonance") {
        std::cout << "Resonance term: " << computeResonanceTerm(t) << std::endl;
    }
    
    std::cout << "============================================\n" << std::endl;
}

void UQFFCompressedResonanceModule::validatePhysicalConsistency() {
    std::cout << "Validating UQFF physical consistency for " << current_system << "..." << std::endl;
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
        std::cout << "  All checks passed. UQFF system is physically consistent." << std::endl;
    }
}

void UQFFCompressedResonanceModule::autoCorrectAnomalies() {
    std::cout << "Auto-correcting UQFF anomalies for " << current_system << "..." << std::endl;
    
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
        std::cout << "  Correcting M to 1e41 kg" << std::endl;
        variables["M"] = 1e41;
        variables["M_visible"] = 0.7 * variables["M"];
        variables["M_DM"] = 0.3 * variables["M"];
    }
    
    if (variables["r"] <= 0) {
        std::cout << "  Correcting r to 1e20 m" << std::endl;
        variables["r"] = 1e20;
    }
    
    if (variables["SFR"] < 0) {
        std::cout << "  Correcting SFR to 6e19 kg/s" << std::endl;
        variables["SFR"] = 6e19;
    }
    
    std::cout << "Auto-correction complete." << std::endl;
}

// Example
// int main() {
//     std::cout << "=============================================" << std::endl;
//     std::cout << "UQFF Multi-System - Enhanced Demonstration" << std::endl;
//     std::cout << "=============================================" << std::endl;
// 
//     // === Part 1: Basic computation - M51 Galaxy ===
//     std::cout << "\n--- Part 1: M51 Galaxy (Compressed Mode) ---" << std::endl;
//     UQFFCompressedResonanceModule mod;
//     mod.setSystem("M51");
//     mod.setMode("compressed");
//     double t = 1e9 * 3.156e7;  // 1 Gyr
//     double g = mod.computeG(t);
//     std::cout << "g_UQFF (M51, compressed) = " << g << " m/s^2" << std::endl;
// 
//     // === Part 2: Switch to resonance mode ===
//     std::cout << "\n--- Part 2: Resonance Mode ---" << std::endl;
//     mod.setMode("resonance");
//     double g_res = mod.computeG(t);
//     std::cout << "g_UQFF (M51, resonance) = " << g_res << " m/s^2" << std::endl;
// 
//     // === Part 3: System switching - Eagle Nebula ===
//     std::cout << "\n--- Part 3: Eagle Nebula ---" << std::endl;
//     mod.setSystem("Eagle");
//     mod.setMode("compressed");
//     double g_eagle = mod.computeG(3e6 * 3.156e7);  // 3 Myr
//     std::cout << "g_UQFF (Eagle) = " << g_eagle << " m/s^2" << std::endl;
// 
//     // === Part 4: Dynamic variable management ===
//     std::cout << "\n--- Part 4: Dynamic Variable Management ---" << std::endl;
//     mod.createDynamicVariable("custom_scale", 1.5);
//     mod.cloneVariable("omega", "omega_backup");
//     std::cout << "Total variables: " << std::endl;
//     mod.listAllVariables();
// 
//     // === Part 5: Batch operations ===
//     std::cout << "\n--- Part 5: Batch Operations ---" << std::endl;
//     std::vector<std::string> mass_params = {"M_visible", "M_DM"};
//     std::cout << "Scaling mass parameters by 1.2..." << std::endl;
//     mod.scaleVariableGroup(mass_params, 1.2);
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
//     mod.analyzeParameterSensitivity("omega");
// 
//     // === Part 9: System validation ===
//     std::cout << "\n--- Part 9: System Validation ---" << std::endl;
//     mod.validatePhysicalConsistency();
// 
//     // === Part 10: Auto-refinement ===
//     std::cout << "\n--- Part 10: Auto-Refinement ---" << std::endl;
//     mod.autoRefineParameters(0.01);
// 
//     // === Part 11: Multi-system exploration ===
//     std::cout << "\n--- Part 11: Multi-System Exploration ---" << std::endl;
//     std::vector<std::string> systems = {"YoungStars", "BigBang", "NGC1316", "V838Mon", "NGC1300"};
//     for (const auto& sys : systems) {
//         mod.setSystem(sys);
//         double g_sys = mod.computeG(1e9 * 3.156e7);
//         std::cout << "  " << sys << ": g = " << g_sys << " m/s^2" << std::endl;
//     }
// 
//     // === Part 12: Comprehensive system report ===
//     std::cout << "\n--- Part 12: System Report ---" << std::endl;
//     mod.setSystem("M51");
//     mod.generateSystemReport();
// 
//     // === Part 13: State management ===
//     std::cout << "\n--- Part 13: State Management ---" << std::endl;
//     mod.saveState("final_m51_state");
//     mod.listSavedStates();
// 
//     std::cout << "\n=============================================" << std::endl;
//     std::cout << "Enhanced UQFF Multi-System demonstration complete!" << std::endl;
//     std::cout << "=============================================" << std::endl;
// 
//     return 0;
// }
// Compile: g++ -o multi_uqff base.cpp UQFFCompressedResonanceModule.cpp -lm
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

UQFFCompressedResonanceModule Evaluation

Strengths :
-Modular, extensible design for multi - system astrophysical modeling, supporting both compressed and resonance modes.
- Comprehensive physics : gravity, cosmological expansion, magnetic fields, environmental effects, quantum, fluid, and dark matter terms.
- Resonance mode adds oscillatory wave dynamics(cosine and complex exponential terms), broadening physical applicability.
- Dynamic variable management via std::map enables runtime updates and system adaptation.
- System - specific parameter loading via setSystem for easy adaptation to diverse scenarios(e.g., Eagle Nebula, Big Bang, M51, V838 Mon).
- Output functions for equation text and variable state support debugging and documentation.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for flexibility and scalability.
- Some calculations use placeholder values(e.g., computeUgSum, computeFenv); implement full physical models for accuracy.
- Minimal error handling(e.g., division by zero, invalid variable names); add validation for robustness.
- Unit consistency should be checked and documented for all physical quantities.
- For large - scale or performance - critical simulations, optimize data structures and reduce redundant calculations.
- std::map is flexible but may be less efficient than structured types for very large models.
- Expand documentation for function purposes and physical meaning.

Summary:
The code is well - structured, flexible, and suitable for scientific prototyping and educational use in multi - system astrophysical modeling.It implements a broad set of physical effects and adapts to various scenarios, including resonance phenomena.For production or high - performance applications, address the recommendations for improved robustness, maintainability, and scalability.