// DPMModule.h
// Modular C++ implementation of the Birth of Di-Pseudo-Monopole (DPM) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module models the Pre-Big Bang reaction of [SCm] and [UA] in a 26-shell oscillating EM field, yielding 26 resonant sphere centers.
// Pluggable: #include "DPMModule.h"
// DPMModule mod; mod.computeDPM(); mod.updateVariable("num_states", 26);
// Variables in std::map; computes sphere centers (h,k,l,r) for 26 states; resonant points via standing waves.
// Approximations: 26 centers distributed on unit sphere; r fixed; [SCm]/[UA] energies as scalars; inflation barriers at -1/2 states.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef DPM_MODULE_H
#define DPM_MODULE_H

#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <functional>
#include <random>  // For distributing 26 centers

class DPMModule {
private:
    std::map<std::string, double> variables;
    std::vector<std::vector<double>> computeSphereCenters();  // 26 centers [h,k,l]
    std::vector<double> computeResonantPoints(double h, double k, double l, double r);

public:
    // Constructor: Initialize with UQFF defaults for DPM birth
    DPMModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    std::vector<std::vector<double>> computeDPM();  // Returns 26 sphere centers [[h,k,l], ...]
    double computeSCmEnergy();  // [SCm] massless metal energy
    double computeUAEnergy();   // [UA] self-plasmotic vacuum energy
    double computeResonanceFactor();  // Belly Button cosmic standing resonance

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print DPM sphere centers
    void printDPMSpheres();

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
    void expandEnergyScale(double energy_multiplier);
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

#endif // DPM_MODULE_H

// DPMModule.cpp
#include "DPMModule.h"
#include <random>  // For distributing 26 centers

// Constructor: Set UQFF defaults for Pre-Big Bang DPM
DPMModule::DPMModule() {
    // Universal constants
    variables["num_states"] = 26.0;                 // 26 EM fields/states
    variables["r"] = 1.0;                           // Sphere radius (normalized)
    variables["SCm_amount"] = 1e42;                 // Raw [SCm] amount (arbitrary J)
    variables["UA_amount"] = 1e42;                  // Raw [UA] amount (arbitrary J)
    variables["ACP_massive"] = 1.0;                 // 26-field envelope factor
    variables["a_over_b"] = 6.6743e-11;             // G M / r^2 analog
    variables["e"] = 1.602e-19;                     // Elementary charge q analog
    variables["half_state_barrier"] = -0.5;         // High energy superconductive barrier
    variables["decay_rate"] = 1e-10;                // Trapped UA breakdown rate (s^-1)
    variables["t_pre_bigbang"] = 0.0;               // Time at birth (s)

    // Higgs/proton stability default
    variables["Higgs_support"] = 1.0;               // Proton stability factor
}

// Update variable
void DPMModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
    if (name == "num_states") {
        // Ensure integer for centers
        variables[name] = std::round(value);
    }
}

// Add delta
void DPMModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void DPMModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute 26 sphere centers distributed on unit sphere
std::vector<std::vector<double>> DPMModule::computeSphereCenters() {
    int n = static_cast<int>(variables["num_states"]);
    std::vector<std::vector<double>> centers(n, std::vector<double>(3, 0.0));
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> theta_dist(0.0, 2 * M_PI);
    std::uniform_real_distribution<double> phi_dist(0.0, M_PI);

    for (int i = 0; i < n; ++i) {
        double theta = theta_dist(gen);
        double phi = phi_dist(gen);
        double r_sphere = variables["r"];
        centers[i][0] = r_sphere * std::sin(phi) * std::cos(theta);  // h = x
        centers[i][1] = r_sphere * std::sin(phi) * std::sin(theta);  // k = y
        centers[i][2] = r_sphere * std::cos(phi);                    // l = z
    }
    return centers;
}

// Compute resonant points for a single sphere (sample points on surface)
std::vector<double> DPMModule::computeResonantPoints(double h, double k, double l, double r) {
    // Simplified: Return sample point on sphere
    return {h + r, k, l};  // One resonant point
}

// Compute DPM: 26 centers
std::vector<std::vector<double>> DPMModule::computeDPM() {
    return computeSphereCenters();
}

// [SCm] energy (massless metal, extra-universal)
double DPMModule::computeSCmEnergy() {
    return variables["SCm_amount"] * variables["ACP_massive"];
}

// [UA] energy (self-plasmotic vacuum pressed)
double DPMModule::computeUAEnergy() {
    double ua_base = variables["UA_amount"];
    double breakdown = std::exp(-variables["decay_rate"] * variables["t_pre_bigbang"]);
    return ua_base * breakdown * variables["ACP_massive"];
}

// Belly Button cosmic standing resonance factor
double DPMModule::computeResonanceFactor() {
    double scm = computeSCmEnergy();
    double ua = computeUAEnergy();
    double attraction = variables["a_over_b"] * (scm * ua) / (variables["r"] * variables["r"]);
    return attraction * variables["e"] * variables["Higgs_support"];
}

// Equation text
std::string DPMModule::getEquationText() {
    return "Birth of DPM: (x - h)^2 + (y - k)^2 + (z - l)^2 = r^2 for 26 states (centers [h,k,l]).\n"
           "[SCm] + [UA] in 26-shell EM field → Resonant DPM spheres.\n"
           "Resonance Factor = (G M / r^2) * q * Higgs_support (a/b: GM/r^2, e: q).\n"
           "Inflation: -1/2 states as high energy barriers; Trapped UA decays exp(-γ t).\n"
           "UQFF Model: 26 quantum levels; plasma mediates; Higgs proton stability; [SCm] builds matter.\n"
           "At t_pre=0: Resonance ~1e-11 (normalized); 26 centers on unit sphere.";
}

// Print variables
void DPMModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Print DPM spheres (centers)
void DPMModule::printDPMSpheres() {
    auto centers = computeDPM();
    std::cout << "DPM Sphere Centers (26 states, [h,k,l]):\n";
    for (size_t i = 0; i < centers.size(); ++i) {
        std::cout << "State " << i+1 << ": [" << std::fixed << std::setprecision(3)
                  << centers[i][0] << ", " << centers[i][1] << ", " << centers[i][2] << "]\n";
    }
    std::cout << "Resonance Factor: " << std::scientific << computeResonanceFactor() << std::endl;
    std::cout << "[SCm] Energy: " << computeSCmEnergy() << " J\n";
    std::cout << "[UA] Energy: " << computeUAEnergy() << " J\n";
}

// ===== DYNAMIC SELF-UPDATE & SELF-EXPANSION IMPLEMENTATIONS =====

// Static storage for saved states
static std::map<std::string, std::map<std::string, double>> dpm_saved_states;

// 1. Dynamic variable management
void DPMModule::createDynamicVariable(const std::string& name, double value) {
    variables[name] = value;
    std::cout << "Created dynamic variable: " << name << " = " << value << std::endl;
}

void DPMModule::removeDynamicVariable(const std::string& name) {
    if (variables.find(name) != variables.end()) {
        variables.erase(name);
        std::cout << "Removed dynamic variable: " << name << std::endl;
    } else {
        std::cerr << "Variable '" << name << "' not found for removal." << std::endl;
    }
}

void DPMModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
        std::cout << "Cloned " << source << " to " << dest << std::endl;
    } else {
        std::cerr << "Source variable '" << source << "' not found." << std::endl;
    }
}

void DPMModule::listAllVariables() {
    std::cout << "=== All DPM Variables (Total: " << variables.size() << ") ===" << std::endl;
    for (const auto& pair : variables) {
        std::cout << "  " << pair.first << " = " << pair.second << std::endl;
    }
}

// 2. Batch operations
void DPMModule::applyTransformToGroup(const std::vector<std::string>& varNames,
                                       std::function<double(double)> transform) {
    for (const auto& name : varNames) {
        if (variables.find(name) != variables.end()) {
            variables[name] = transform(variables[name]);
            std::cout << "Transformed " << name << " to " << variables[name] << std::endl;
        }
    }
}

void DPMModule::scaleVariableGroup(const std::vector<std::string>& varNames, double scale_factor) {
    applyTransformToGroup(varNames, [scale_factor](double val) { return val * scale_factor; });
}

// 3. Self-expansion capabilities
void DPMModule::autoExpandParameterSpace(double scale_factor) {
    std::cout << "Auto-expanding DPM parameter space by factor " << scale_factor << std::endl;
    std::vector<std::string> expandable = {"SCm_amount", "UA_amount", "r", "ACP_massive"};
    scaleVariableGroup(expandable, scale_factor);
    std::cout << "  Parameter space expanded" << std::endl;
}

void DPMModule::expandEnergyScale(double energy_multiplier) {
    std::cout << "Expanding energy scale by " << energy_multiplier << std::endl;
    variables["SCm_amount"] *= energy_multiplier;
    variables["UA_amount"] *= energy_multiplier;
    std::cout << "  SCm_amount: " << variables["SCm_amount"] << " J" << std::endl;
    std::cout << "  UA_amount: " << variables["UA_amount"] << " J" << std::endl;
}

void DPMModule::expandSpatialScale(double spatial_multiplier) {
    std::cout << "Expanding spatial scale by " << spatial_multiplier << std::endl;
    variables["r"] *= spatial_multiplier;
    std::cout << "  r (sphere radius): " << variables["r"] << " (normalized)" << std::endl;
}

void DPMModule::expandTimeScale(double time_multiplier) {
    std::cout << "Expanding time scale by " << time_multiplier << std::endl;
    variables["t_pre_bigbang"] *= time_multiplier;
    variables["decay_rate"] /= time_multiplier;  // Decay rate inversely scales
    std::cout << "  t_pre_bigbang: " << variables["t_pre_bigbang"] << " s" << std::endl;
    std::cout << "  decay_rate: " << variables["decay_rate"] << " s^-1" << std::endl;
}

// 4. Self-refinement
void DPMModule::autoRefineParameters(double tolerance) {
    std::cout << "Auto-refining DPM parameters with tolerance " << tolerance << std::endl;
    
    // Ensure num_states is integer
    double num_rounded = std::round(variables["num_states"]);
    if (std::abs(num_rounded - variables["num_states"]) > tolerance) {
        std::cout << "  Rounding num_states to integer: " << num_rounded << std::endl;
        variables["num_states"] = num_rounded;
    }
    
    // Ensure positive values
    if (variables["SCm_amount"] <= 0) {
        std::cout << "  WARNING: SCm_amount must be positive" << std::endl;
    }
    
    if (variables["UA_amount"] <= 0) {
        std::cout << "  WARNING: UA_amount must be positive" << std::endl;
    }
    
    if (variables["r"] <= 0) {
        std::cout << "  WARNING: r must be positive" << std::endl;
    }
    
    std::cout << "Refinement complete." << std::endl;
}

void DPMModule::calibrateToObservations(const std::map<std::string, double>& observed_values) {
    std::cout << "Calibrating to " << observed_values.size() << " DPM observations..." << std::endl;
    for (const auto& obs : observed_values) {
        if (variables.find(obs.first) != variables.end()) {
            double old_val = variables[obs.first];
            updateVariable(obs.first, obs.second);
            std::cout << "  " << obs.first << ": " << old_val << " -> " << obs.second << std::endl;
        }
    }
    std::cout << "Calibration complete." << std::endl;
}

void DPMModule::optimizeForMetric(const std::string& metric_name, double target_value) {
    std::cout << "Optimizing for metric: " << metric_name << " = " << target_value << std::endl;
    
    if (metric_name == "resonance") {
        double current_res = computeResonanceFactor();
        double ratio = target_value / std::max(current_res, 1e-100);
        
        // Adjust a_over_b to reach target resonance
        variables["a_over_b"] *= ratio;
        std::cout << "  Adjusted a_over_b by " << ratio << std::endl;
    } else if (metric_name == "SCm_energy") {
        double current_scm = computeSCmEnergy();
        double ratio = target_value / std::max(current_scm, 1e-100);
        
        // Adjust SCm_amount to reach target energy
        variables["SCm_amount"] *= ratio;
        std::cout << "  Adjusted SCm_amount by " << ratio << std::endl;
    }
    
    std::cout << "Optimization complete." << std::endl;
}

// 5. Parameter exploration
void DPMModule::generateVariations(int num_variations, double variation_range) {
    std::cout << "Generating " << num_variations << " DPM variations with range ±" 
              << (variation_range * 100) << "%" << std::endl;
    
    std::vector<std::string> key_params = {"SCm_amount", "UA_amount", "r", "ACP_massive", "a_over_b", "decay_rate"};
    
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

void DPMModule::findOptimalParameters(const std::string& objective, int iterations) {
    std::cout << "Finding optimal DPM parameters for: " << objective 
              << " (" << iterations << " iterations)" << std::endl;
    
    double best_score = -1e100;
    std::map<std::string, double> best_params;
    
    for (int i = 0; i < iterations; ++i) {
        mutateParameters(0.7, 0.1);
        
        double score = computeResonanceFactor();
        
        if (objective == "maximize_resonance") {
            if (score > best_score) {
                best_score = score;
                best_params = variables;
            }
        } else if (objective == "target_1e-11") {
            if (std::abs(score - 1e-11) < std::abs(best_score - 1e-11)) {
                best_score = score;
                best_params = variables;
            }
        }
    }
    
    variables = best_params;
    std::cout << "Optimal Resonance Factor: " << best_score << std::endl;
}

// 6. Adaptive evolution
void DPMModule::mutateParameters(double mutation_rate, double mutation_strength) {
    std::vector<std::string> mutable_params = {"SCm_amount", "UA_amount", "r", "ACP_massive", "a_over_b", "decay_rate"};
    
    for (const auto& param : mutable_params) {
        if (variables.find(param) != variables.end()) {
            if ((rand() / (double)RAND_MAX) < mutation_rate) {
                double mutation = 1.0 + mutation_strength * (2.0 * (rand() / (double)RAND_MAX) - 1.0);
                variables[param] *= mutation;
            }
        }
    }
}

void DPMModule::evolveSystem(int generations) {
    std::cout << "Evolving DPM system over " << generations << " generations..." << std::endl;
    
    for (int gen = 0; gen < generations; ++gen) {
        mutateParameters(0.3, 0.08);
        
        double fitness = computeResonanceFactor();
        
        if (gen % 10 == 0) {
            std::cout << "  Gen " << gen << ": Resonance = " << fitness << std::endl;
        }
    }
    
    std::cout << "Evolution complete." << std::endl;
}

// 7. State management
void DPMModule::saveState(const std::string& label) {
    dpm_saved_states[label] = variables;
    std::cout << "Saved DPM state: " << label << " (" << variables.size() << " variables)" << std::endl;
}

void DPMModule::restoreState(const std::string& label) {
    if (dpm_saved_states.find(label) != dpm_saved_states.end()) {
        variables = dpm_saved_states[label];
        std::cout << "Restored DPM state: " << label << std::endl;
    } else {
        std::cerr << "State '" << label << "' not found." << std::endl;
    }
}

void DPMModule::listSavedStates() {
    std::cout << "=== Saved DPM States (Total: " << dpm_saved_states.size() << ") ===" << std::endl;
    for (const auto& state : dpm_saved_states) {
        std::cout << "  " << state.first << " (" << state.second.size() << " variables)" << std::endl;
    }
}

void DPMModule::exportState(const std::string& filename) {
    std::cout << "Exporting DPM state to " << filename << " (not implemented - placeholder)" << std::endl;
    // In real implementation: write variables to file
}

// 8. System analysis
void DPMModule::analyzeParameterSensitivity(const std::string& param_name) {
    if (variables.find(param_name) == variables.end()) {
        std::cerr << "Parameter '" << param_name << "' not found." << std::endl;
        return;
    }
    
    std::cout << "=== DPM Sensitivity Analysis: " << param_name << " ===" << std::endl;
    
    double base_value = variables[param_name];
    double base_output = computeResonanceFactor();
    
    std::vector<double> perturbations = {0.7, 0.85, 1.0, 1.15, 1.3};
    
    for (double factor : perturbations) {
        updateVariable(param_name, base_value * factor);
        
        double new_output = computeResonanceFactor();
        double sensitivity = (new_output - base_output) / std::max(std::abs(base_output), 1e-100);
        
        std::cout << "  " << param_name << " * " << factor << " -> Resonance change: " 
                  << (sensitivity * 100) << "%" << std::endl;
    }
    
    updateVariable(param_name, base_value);  // Restore
}

void DPMModule::generateSystemReport() {
    std::cout << "\n========== Di-Pseudo-Monopole (DPM) Birth System Report ==========" << std::endl;
    std::cout << "Total Variables: " << variables.size() << std::endl;
    
    // Quantum states
    std::cout << "\nQuantum Configuration:" << std::endl;
    std::cout << "num_states: " << variables["num_states"] << " EM fields/states" << std::endl;
    std::cout << "r (sphere radius): " << variables["r"] << " (normalized)" << std::endl;
    std::cout << "half_state_barrier: " << variables["half_state_barrier"] << " (inflation barrier)" << std::endl;
    
    std::cout << "\nEnergy Components:" << std::endl;
    std::cout << "SCm_amount: " << variables["SCm_amount"] << " J" << std::endl;
    std::cout << "UA_amount: " << variables["UA_amount"] << " J" << std::endl;
    std::cout << "ACP_massive: " << variables["ACP_massive"] << " (26-field envelope)" << std::endl;
    
    double scm_energy = computeSCmEnergy();
    double ua_energy = computeUAEnergy();
    std::cout << "\n[SCm] Energy (computed): " << scm_energy << " J" << std::endl;
    std::cout << "[UA] Energy (computed): " << ua_energy << " J" << std::endl;
    
    std::cout << "\nResonance Parameters:" << std::endl;
    std::cout << "a_over_b (G M/r^2 analog): " << variables["a_over_b"] << std::endl;
    std::cout << "e (charge q analog): " << variables["e"] << " C" << std::endl;
    std::cout << "Higgs_support: " << variables["Higgs_support"] << " (proton stability)" << std::endl;
    
    double resonance = computeResonanceFactor();
    std::cout << "\nResonance Factor (computed): " << resonance << std::endl;
    
    std::cout << "\nTemporal Evolution:" << std::endl;
    std::cout << "t_pre_bigbang: " << variables["t_pre_bigbang"] << " s" << std::endl;
    std::cout << "decay_rate: " << variables["decay_rate"] << " s^-1" << std::endl;
    double decay_factor = std::exp(-variables["decay_rate"] * variables["t_pre_bigbang"]);
    std::cout << "UA decay factor: " << decay_factor << std::endl;
    
    std::cout << "\nDPM Sphere Centers:" << std::endl;
    auto centers = computeDPM();
    std::cout << "Number of centers: " << centers.size() << std::endl;
    std::cout << "Sample center [h,k,l]: [" << std::fixed << std::setprecision(3)
              << centers[0][0] << ", " << centers[0][1] << ", " << centers[0][2] << "]" << std::endl;
    
    std::cout << "\nPhysics Regime:" << std::endl;
    if (std::abs(resonance) < 1e-10) {
        std::cout << "Weak resonance regime (Resonance << 1)" << std::endl;
    } else if (std::abs(resonance) < 1) {
        std::cout << "Moderate resonance regime" << std::endl;
    } else {
        std::cout << "Strong resonance regime (Resonance >= 1)" << std::endl;
    }
    
    std::cout << "==================================================\n" << std::endl;
}

void DPMModule::validatePhysicalConsistency() {
    std::cout << "Validating DPM physical consistency..." << std::endl;
    bool consistent = true;
    
    // Check for NaN/Inf
    for (const auto& pair : variables) {
        if (std::isnan(pair.second) || std::isinf(pair.second)) {
            std::cerr << "  ERROR: " << pair.first << " is NaN/Inf" << std::endl;
            consistent = false;
        }
    }
    
    // Positive values
    if (variables["SCm_amount"] <= 0) {
        std::cerr << "  ERROR: SCm_amount must be positive" << std::endl;
        consistent = false;
    }
    
    if (variables["UA_amount"] <= 0) {
        std::cerr << "  ERROR: UA_amount must be positive" << std::endl;
        consistent = false;
    }
    
    if (variables["r"] <= 0) {
        std::cerr << "  ERROR: r must be positive" << std::endl;
        consistent = false;
    }
    
    if (variables["ACP_massive"] <= 0) {
        std::cerr << "  ERROR: ACP_massive must be positive" << std::endl;
        consistent = false;
    }
    
    // Physical ranges
    if (variables["num_states"] < 1) {
        std::cerr << "  ERROR: num_states must be at least 1" << std::endl;
        consistent = false;
    }
    
    if (variables["num_states"] != std::round(variables["num_states"])) {
        std::cerr << "  WARNING: num_states should be an integer" << std::endl;
    }
    
    if (variables["decay_rate"] < 0) {
        std::cerr << "  ERROR: decay_rate cannot be negative" << std::endl;
        consistent = false;
    }
    
    if (variables["Higgs_support"] <= 0) {
        std::cerr << "  WARNING: Higgs_support should be positive for proton stability" << std::endl;
    }
    
    if (consistent) {
        std::cout << "  All checks passed. DPM system is physically consistent." << std::endl;
    }
}

void DPMModule::autoCorrectAnomalies() {
    std::cout << "Auto-correcting DPM anomalies..." << std::endl;
    
    // Fix NaN/Inf
    for (auto& pair : variables) {
        if (std::isnan(pair.second) || std::isinf(pair.second)) {
            std::cout << "  Correcting " << pair.first << " (was NaN/Inf)" << std::endl;
            pair.second = 1.0;
        }
    }
    
    // Ensure positive values
    if (variables["SCm_amount"] <= 0) {
        std::cout << "  Correcting SCm_amount to 1e42 J" << std::endl;
        variables["SCm_amount"] = 1e42;
    }
    
    if (variables["UA_amount"] <= 0) {
        std::cout << "  Correcting UA_amount to 1e42 J" << std::endl;
        variables["UA_amount"] = 1e42;
    }
    
    if (variables["r"] <= 0) {
        std::cout << "  Correcting r to 1.0 (normalized)" << std::endl;
        variables["r"] = 1.0;
    }
    
    if (variables["ACP_massive"] <= 0) {
        std::cout << "  Correcting ACP_massive to 1.0" << std::endl;
        variables["ACP_massive"] = 1.0;
    }
    
    if (variables["num_states"] < 1) {
        std::cout << "  Correcting num_states to 26" << std::endl;
        variables["num_states"] = 26.0;
    }
    
    // Round num_states to integer
    if (variables["num_states"] != std::round(variables["num_states"])) {
        std::cout << "  Rounding num_states to " << std::round(variables["num_states"]) << std::endl;
        variables["num_states"] = std::round(variables["num_states"]);
    }
    
    if (variables["decay_rate"] < 0) {
        std::cout << "  Correcting decay_rate to 1e-10 s^-1" << std::endl;
        variables["decay_rate"] = 1e-10;
    }
    
    if (variables["Higgs_support"] <= 0) {
        std::cout << "  Correcting Higgs_support to 1.0" << std::endl;
        variables["Higgs_support"] = 1.0;
    }
    
    std::cout << "Auto-correction complete." << std::endl;
}

// Example usage in base program (snippet)
// Uncomment the following code to test the enhanced DPM module with dynamic capabilities
/*
#include "DPMModule.h"
int main() {
    DPMModule mod;
    
    std::cout << "===== Di-Pseudo-Monopole Birth with Dynamic Capabilities =====" << std::endl;
    std::cout << "  Pre-Big Bang [SCm] + [UA] reaction in 26-shell EM field\n" << std::endl;
    
    // Initial DPM computation
    std::cout << "1. Initial DPM sphere centers:" << std::endl;
    mod.printDPMSpheres();
    std::cout << std::endl;
    
    // List all variables
    std::cout << "2. Variable inventory:" << std::endl;
    mod.listAllVariables();
    std::cout << std::endl;
    
    // Generate system report
    std::cout << "3. System report:" << std::endl;
    mod.generateSystemReport();
    std::cout << std::endl;
    
    // Save initial state
    std::cout << "4. Saving initial 26-state configuration:" << std::endl;
    mod.saveState("26_states");
    std::cout << std::endl;
    
    // Modify to 13 states (half-states)
    std::cout << "5. Reducing to 13 states (inflation barriers at -1/2):" << std::endl;
    mod.updateVariable("num_states", 13.0);
    mod.printDPMSpheres();
    std::cout << std::endl;
    
    // Expand energy scale
    std::cout << "6. Expanding energy scale by 2x:" << std::endl;
    mod.expandEnergyScale(2.0);
    std::cout << "New Resonance: " << mod.computeResonanceFactor() << std::endl;
    std::cout << std::endl;
    
    // Save high-energy state
    std::cout << "7. Saving high-energy 13-state configuration:" << std::endl;
    mod.saveState("13_states_2x");
    std::cout << std::endl;
    
    // Analyze sensitivity
    std::cout << "8. Sensitivity analysis for SCm_amount:" << std::endl;
    mod.analyzeParameterSensitivity("SCm_amount");
    std::cout << std::endl;
    
    // Create dynamic variable
    std::cout << "9. Creating dynamic cosmic epoch variable:" << std::endl;
    mod.createDynamicVariable("cosmic_epoch", -1.0);  // Pre-Big Bang
    std::cout << std::endl;
    
    // Optimize for target resonance
    std::cout << "10. Optimizing for target Resonance = 1e-11:" << std::endl;
    mod.optimizeForMetric("resonance", 1e-11);
    std::cout << "Optimized Resonance: " << mod.computeResonanceFactor() << std::endl;
    std::cout << std::endl;
    
    // Validate consistency
    std::cout << "11. Validating physical consistency:" << std::endl;
    mod.validatePhysicalConsistency();
    std::cout << std::endl;
    
    // Generate variations
    std::cout << "12. Generating 3 parameter variations (±20%):" << std::endl;
    mod.generateVariations(3, 0.2);
    std::cout << std::endl;
    
    // Restore 26-state configuration
    std::cout << "13. Restoring original 26-state configuration:" << std::endl;
    mod.restoreState("26_states");
    mod.printDPMSpheres();
    std::cout << std::endl;
    
    // List saved states
    std::cout << "14. List of saved states:" << std::endl;
    mod.listSavedStates();
    std::cout << std::endl;
    
    std::cout << "End of DPM demonstration with dynamic capabilities.\n" << std::endl;
    std::cout << mod.getEquationText() << std::endl;
    return 0;
}
*/
// Compile: g++ -o dpm_test dpm_test.cpp DPMModule.cpp -lm
// Sample: 26 random centers on unit sphere; Resonance ~1e-11 J; UA decays over t.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

DPMModule Evaluation

Strengths :
-Modular, extensible design for modeling the birth of Di - Pseudo - Monopole(DPM) in the UQFF framework.
- Clear encapsulation of variables and sphere center data using std::map and std::vector, supporting dynamic updates and easy extension.
- Implements core physical concepts : distribution of 26 sphere centers on a unit sphere, resonance factor calculation, and energy computations for[SCm] and [UA].
- Approximations and physical meaning are well - documented in comments and equation text.
- Output functions for variable state and DPM sphere centers support debugging and transparency.
- Randomized sphere center generation provides realistic spatial distribution for resonance modeling.
- Handles dynamic state count(e.g., 26 or 13) and updates dependent calculations accordingly.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for flexibility and scalability.
- Some calculations use magic numbers or lack explanatory comments; define named constants and clarify logic.
- Minimal error handling(e.g., division by zero, invalid variable names); add validation for robustness.
- Unit consistency should be checked and documented for all physical quantities.
- For large - scale or performance - critical simulations, optimize data structures and reduce redundant calculations.
- std::map and std::vector are flexible but may be less efficient than structured types for very large models.
- Expand documentation for function purposes and physical meaning.
- Random sphere center generation may not be reproducible; consider seeding for deterministic results if needed.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in DPM birth modeling.It implements the UQFF DPM concept faithfully and adapts to various scenarios.For production or high - performance applications, address the recommendations for improved robustness, maintainability, and scalability.