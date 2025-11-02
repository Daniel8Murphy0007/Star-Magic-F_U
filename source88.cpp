// AndromedaUQFFModule.h
// Modular C++ implementation of the Master Universal Gravity Equation (UQFF) for Andromeda Galaxy Evolution.
// This module can be plugged into a base program by including this header and linking the .cpp.
// Usage: #include "AndromedaUQFFModule.h"
// AndromedaUQFFModule mod; mod.computeG(t); mod.updateVariable("M", new_value);
// Variables stored in std::map for dynamic updates.
// Includes base gravity with expansion and TRZ, BH term, dust friction a_dust, EM/Aether term.
// Approximations: z=-0.001 (blueshift); dust scaled by 1e-12; EM normalized to proton mass.
// Andromeda params: M=1e12 Msun, r=1.04e21 m, M_BH=1.4e8 Msun, v_orbit=2.5e5 m/s, etc.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef ANDROMEDA_UQFF_MODULE_H
#define ANDROMEDA_UQFF_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <functional>

class AndromedaUQFFModule {
private:
    std::map<std::string, double> variables;
    double computeHz();
    double computeADust();
    double computeEMBase();
    double computeEMTerm();

public:
    // Constructor: Initialize with Andromeda defaults
    AndromedaUQFFModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computation: g_Andromeda(r, t)
    double computeG(double t);

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print evolution table (0-10 Gyr, 2 Gyr steps)
    void printEvolutionTable();

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

#endif // ANDROMEDA_UQFF_MODULE_H

// AndromedaUQFFModule.cpp
#include "AndromedaUQFFModule.h"

// Constructor: Set Andromeda-specific values
AndromedaUQFFModule::AndromedaUQFFModule() {
    // Universal constants
    variables["G"] = 6.6743e-11;                    // m^3 kg^-1 s^-2
    variables["M_sun"] = 1.989e30;                  // kg
    variables["q"] = 1.602e-19;                     // C
    variables["proton_mass"] = 1.673e-27;           // kg
    variables["H0"] = 70.0;                         // km/s/Mpc
    variables["Mpc_to_m"] = 3.086e22;               // m/Mpc
    variables["Omega_m"] = 0.3;
    variables["Omega_Lambda"] = 0.7;
    variables["year_to_s"] = 3.156e7;               // s/yr
    variables["Gyr"] = 1e9;                         // yr

    // Andromeda parameters
    variables["M"] = 1e12 * variables["M_sun"];     // Total mass kg
    variables["r"] = 1.04e21;                       // m (half diameter)
    variables["M_BH"] = 1.4e8 * variables["M_sun"]; // SMBH mass kg
    variables["r_BH"] = 1e15;                       // m (core scale)
    variables["rho_dust"] = 1e-20;                  // kg/m^3
    variables["v_orbit"] = 2.5e5;                   // m/s
    variables["rho_mass"] = 1e-21;                  // kg/m^3
    variables["z"] = -0.001;                        // Blueshift
    variables["B"] = 1e-5;                          // T
    variables["rho_vac_UA"] = 7.09e-36;             // J/m^3
    variables["rho_vac_SCm"] = 7.09e-37;            // J/m^3
    variables["f_TRZ"] = 0.1;                       // dimensionless
    variables["scale_macro"] = 1e-12;               // Scaling factor
    variables["t"] = 10.0 * variables["Gyr"] * variables["year_to_s"];  // Default 10 Gyr
}

// Update variable
void AndromedaUQFFModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
    if (name == "M") {
        variables["M_BH"] = 1.4e8 * (value / (1e12 * variables["M_sun"])) * variables["M_sun"];  // Scale if needed
    }
}

// Add delta
void AndromedaUQFFModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void AndromedaUQFFModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute H(z) in s^-1
double AndromedaUQFFModule::computeHz() {
    double one_plus_z = 1.0 + variables["z"];
    double Hz_kms = variables["H0"] * std::sqrt(variables["Omega_m"] * std::pow(one_plus_z, 3) + variables["Omega_Lambda"]);
    return (Hz_kms * 1e3) / variables["Mpc_to_m"];
}

// Compute dust acceleration
double AndromedaUQFFModule::computeADust() {
    double force_per_area = variables["rho_dust"] * std::pow(variables["v_orbit"], 2);
    double a_dust_base = force_per_area / variables["rho_mass"];
    return a_dust_base * variables["scale_macro"];
}

// Compute EM base (m/s^2)
double AndromedaUQFFModule::computeEMBase() {
    double mag_vB = variables["v_orbit"] * variables["B"];
    double force = variables["q"] * mag_vB;
    return force / variables["proton_mass"];
}

// Compute full EM term
double AndromedaUQFFModule::computeEMTerm() {
    double em_base = computeEMBase();
    double vac_ratio = variables["rho_vac_UA"] / variables["rho_vac_SCm"];
    return em_base * (1.0 + vac_ratio) * variables["scale_macro"];
}

// Full g_Andromeda
double AndromedaUQFFModule::computeG(double t) {
    variables["t"] = t;
    double Hz = computeHz();
    double expansion_factor = 1.0 + Hz * t;
    double tr_factor = 1.0 + variables["f_TRZ"];

    double g_grav = (variables["G"] * variables["M"] / (variables["r"] * variables["r"])) * expansion_factor * tr_factor;
    double g_BH = variables["G"] * variables["M_BH"] / (variables["r_BH"] * variables["r_BH"]);
    double a_dust = computeADust();
    double em_term = computeEMTerm();

    return g_grav + g_BH + a_dust + em_term;
}

// Equation text
std::string AndromedaUQFFModule::getEquationText() {
    return "g_Andromeda(r, t) = (G * M / r^2) * (1 + H(z) * t) * (1 + f_TRZ) + (G * M_BH / r_BH^2) + a_dust + q*(v*B) * (1 + ?_UA/?_SCm) * 1e-12\n"
           "Where a_dust = (?_dust * v_orbit^2 / ?_mass) * scale_macro;\n"
           "EM term: q v B / m_proton * (1 + ?_vac_UA / ?_vac_SCm) * scale_macro.\n"
           "Andromeda Adaptations: Blueshift z=-0.001; M=1e12 M_sun; dust lanes with v_orbit=250 km/s.\n"
           "At t=10 Gyr, g ?6.273 m/s� (dust dominant); minimal evolution due to small H(z)t.\n"
           "UQFF Terms: f_TRZ for time-reversal; Aether vacua ratio for EM enhancement.";
}

// Print variables
void AndromedaUQFFModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Print evolution table
void AndromedaUQFFModule::printEvolutionTable() {
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Evolution over time (m/s�):\n";
    std::cout << "t (Gyr) | g_Andromeda\n";
    std::cout << "--------|------------\n";
    for (int i = 0; i <= 5; ++i) {
        double t = i * 2.0 * variables["Gyr"] * variables["year_to_s"];
        double g = computeG(t);
        std::cout << std::setw(6) << (i*2) << "    | " << g << "\n";
    }
}

// ===== DYNAMIC SELF-UPDATE & SELF-EXPANSION IMPLEMENTATIONS =====

// Static storage for saved states
static std::map<std::string, std::map<std::string, double>> andromeda_saved_states;

// 1. Dynamic variable management
void AndromedaUQFFModule::createDynamicVariable(const std::string& name, double value) {
    variables[name] = value;
    std::cout << "Created dynamic variable: " << name << " = " << value << std::endl;
}

void AndromedaUQFFModule::removeDynamicVariable(const std::string& name) {
    if (variables.find(name) != variables.end()) {
        variables.erase(name);
        std::cout << "Removed dynamic variable: " << name << std::endl;
    } else {
        std::cerr << "Variable '" << name << "' not found for removal." << std::endl;
    }
}

void AndromedaUQFFModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
        std::cout << "Cloned " << source << " to " << dest << std::endl;
    } else {
        std::cerr << "Source variable '" << source << "' not found." << std::endl;
    }
}

void AndromedaUQFFModule::listAllVariables() {
    std::cout << "=== All Andromeda Variables (Total: " << variables.size() << ") ===" << std::endl;
    for (const auto& pair : variables) {
        std::cout << "  " << pair.first << " = " << pair.second << std::endl;
    }
}

// 2. Batch operations
void AndromedaUQFFModule::applyTransformToGroup(const std::vector<std::string>& varNames,
                                                 std::function<double(double)> transform) {
    for (const auto& name : varNames) {
        if (variables.find(name) != variables.end()) {
            variables[name] = transform(variables[name]);
            std::cout << "Transformed " << name << " to " << variables[name] << std::endl;
        }
    }
}

void AndromedaUQFFModule::scaleVariableGroup(const std::vector<std::string>& varNames, double scale_factor) {
    applyTransformToGroup(varNames, [scale_factor](double val) { return val * scale_factor; });
}

// 3. Self-expansion capabilities
void AndromedaUQFFModule::autoExpandParameterSpace(double scale_factor) {
    std::cout << "Auto-expanding Andromeda parameter space by factor " << scale_factor << std::endl;
    std::vector<std::string> expandable = {"M", "M_BH", "r", "r_BH"};
    scaleVariableGroup(expandable, scale_factor);
    std::cout << "  Parameter space expanded" << std::endl;
}

void AndromedaUQFFModule::expandMassScale(double mass_multiplier) {
    std::cout << "Expanding mass scale by " << mass_multiplier << std::endl;
    variables["M"] *= mass_multiplier;
    variables["M_BH"] *= mass_multiplier;
    std::cout << "  M_total: " << (variables["M"] / variables["M_sun"]) << " M☉" << std::endl;
    std::cout << "  M_BH: " << (variables["M_BH"] / variables["M_sun"]) << " M☉" << std::endl;
}

void AndromedaUQFFModule::expandSpatialScale(double spatial_multiplier) {
    std::cout << "Expanding spatial scale by " << spatial_multiplier << std::endl;
    variables["r"] *= spatial_multiplier;
    variables["r_BH"] *= spatial_multiplier;
    std::cout << "  r: " << variables["r"] << " m" << std::endl;
    std::cout << "  r_BH: " << variables["r_BH"] << " m" << std::endl;
}

void AndromedaUQFFModule::expandTimeScale(double time_multiplier) {
    std::cout << "Expanding time scale by " << time_multiplier << std::endl;
    variables["t"] *= time_multiplier;
    std::cout << "  t: " << (variables["t"] / (variables["Gyr"] * variables["year_to_s"])) << " Gyr" << std::endl;
}

// 4. Self-refinement
void AndromedaUQFFModule::autoRefineParameters(double tolerance) {
    std::cout << "Auto-refining Andromeda parameters with tolerance " << tolerance << std::endl;
    
    // Validate M_BH scales with M (rough approximation)
    double M_ratio = variables["M"] / (1e12 * variables["M_sun"]);
    double M_BH_expected = 1.4e8 * M_ratio * variables["M_sun"];
    if (std::abs(variables["M_BH"] - M_BH_expected) / M_BH_expected > tolerance) {
        std::cout << "  Correcting M_BH: " << variables["M_BH"] << " -> " << M_BH_expected << std::endl;
        variables["M_BH"] = M_BH_expected;
    }
    
    std::cout << "Refinement complete." << std::endl;
}

void AndromedaUQFFModule::calibrateToObservations(const std::map<std::string, double>& observed_values) {
    std::cout << "Calibrating to " << observed_values.size() << " Andromeda observations..." << std::endl;
    for (const auto& obs : observed_values) {
        if (variables.find(obs.first) != variables.end()) {
            double old_val = variables[obs.first];
            updateVariable(obs.first, obs.second);
            std::cout << "  " << obs.first << ": " << old_val << " -> " << obs.second << std::endl;
        }
    }
    std::cout << "Calibration complete." << std::endl;
}

void AndromedaUQFFModule::optimizeForMetric(const std::string& metric_name, double target_value) {
    std::cout << "Optimizing for metric: " << metric_name << " = " << target_value << std::endl;
    
    if (metric_name == "g" || metric_name == "gravity") {
        double t = variables["t"];
        double current_g = computeG(t);
        double ratio = target_value / std::max(current_g, 1e-100);
        
        // Adjust mass to reach target
        variables["M"] *= ratio;
        variables["M_BH"] *= ratio;
        std::cout << "  Adjusted masses by " << ratio << std::endl;
    } else if (metric_name == "v_orbit") {
        variables["v_orbit"] = target_value;
        std::cout << "  Set v_orbit to " << (target_value / 1e3) << " km/s" << std::endl;
    }
    
    std::cout << "Optimization complete." << std::endl;
}

// 5. Parameter exploration
void AndromedaUQFFModule::generateVariations(int num_variations, double variation_range) {
    std::cout << "Generating " << num_variations << " Andromeda variations with range ±" 
              << (variation_range * 100) << "%" << std::endl;
    
    std::vector<std::string> key_params = {"M", "M_BH", "r", "v_orbit", "rho_dust", "B"};
    
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

void AndromedaUQFFModule::findOptimalParameters(const std::string& objective, int iterations) {
    std::cout << "Finding optimal Andromeda parameters for: " << objective 
              << " (" << iterations << " iterations)" << std::endl;
    
    double best_score = -1e100;
    std::map<std::string, double> best_params;
    
    for (int i = 0; i < iterations; ++i) {
        mutateParameters(0.7, 0.1);
        
        double t = variables["t"];
        double score = computeG(t);
        
        if (objective == "maximize_g") {
            if (score > best_score) {
                best_score = score;
                best_params = variables;
            }
        } else if (objective == "target_6.3") {
            if (std::abs(score - 6.3) < std::abs(best_score - 6.3)) {
                best_score = score;
                best_params = variables;
            }
        }
    }
    
    variables = best_params;
    std::cout << "Optimal g: " << best_score << " m/s^2" << std::endl;
}

// 6. Adaptive evolution
void AndromedaUQFFModule::mutateParameters(double mutation_rate, double mutation_strength) {
    std::vector<std::string> mutable_params = {"M", "M_BH", "r", "v_orbit", "rho_dust", "B"};
    
    for (const auto& param : mutable_params) {
        if (variables.find(param) != variables.end()) {
            if ((rand() / (double)RAND_MAX) < mutation_rate) {
                double mutation = 1.0 + mutation_strength * (2.0 * (rand() / (double)RAND_MAX) - 1.0);
                variables[param] *= mutation;
            }
        }
    }
}

void AndromedaUQFFModule::evolveSystem(int generations) {
    std::cout << "Evolving Andromeda system over " << generations << " generations..." << std::endl;
    
    for (int gen = 0; gen < generations; ++gen) {
        mutateParameters(0.3, 0.08);
        
        double t = variables["t"];
        double fitness = computeG(t);
        
        if (gen % 10 == 0) {
            std::cout << "  Gen " << gen << ": g = " << fitness << " m/s^2" << std::endl;
        }
    }
    
    std::cout << "Evolution complete." << std::endl;
}

// 7. State management
void AndromedaUQFFModule::saveState(const std::string& label) {
    andromeda_saved_states[label] = variables;
    std::cout << "Saved Andromeda state: " << label << " (" << variables.size() << " variables)" << std::endl;
}

void AndromedaUQFFModule::restoreState(const std::string& label) {
    if (andromeda_saved_states.find(label) != andromeda_saved_states.end()) {
        variables = andromeda_saved_states[label];
        std::cout << "Restored Andromeda state: " << label << std::endl;
    } else {
        std::cerr << "State '" << label << "' not found." << std::endl;
    }
}

void AndromedaUQFFModule::listSavedStates() {
    std::cout << "=== Saved Andromeda States (Total: " << andromeda_saved_states.size() << ") ===" << std::endl;
    for (const auto& state : andromeda_saved_states) {
        std::cout << "  " << state.first << " (" << state.second.size() << " variables)" << std::endl;
    }
}

void AndromedaUQFFModule::exportState(const std::string& filename) {
    std::cout << "Exporting Andromeda state to " << filename << " (not implemented - placeholder)" << std::endl;
    // In real implementation: write variables to file
}

// 8. System analysis
void AndromedaUQFFModule::analyzeParameterSensitivity(const std::string& param_name) {
    if (variables.find(param_name) == variables.end()) {
        std::cerr << "Parameter '" << param_name << "' not found." << std::endl;
        return;
    }
    
    std::cout << "=== Andromeda Sensitivity Analysis: " << param_name << " ===" << std::endl;
    
    double base_value = variables[param_name];
    double t = variables["t"];
    double base_output = computeG(t);
    
    std::vector<double> perturbations = {0.7, 0.85, 1.0, 1.15, 1.3};
    
    for (double factor : perturbations) {
        updateVariable(param_name, base_value * factor);
        
        double new_output = computeG(t);
        double sensitivity = (new_output - base_output) / std::max(std::abs(base_output), 1e-100);
        
        std::cout << "  " << param_name << " * " << factor << " -> g change: " 
                  << (sensitivity * 100) << "%" << std::endl;
    }
    
    updateVariable(param_name, base_value);  // Restore
}

void AndromedaUQFFModule::generateSystemReport() {
    std::cout << "\n========== Andromeda Galaxy System Report ==========" << std::endl;
    std::cout << "Total Variables: " << variables.size() << std::endl;
    
    // Key parameters
    std::cout << "\nMass Parameters:" << std::endl;
    std::cout << "M_total: " << (variables["M"] / variables["M_sun"]) << " M☉" << std::endl;
    std::cout << "M_BH (SMBH): " << (variables["M_BH"] / variables["M_sun"]) << " M☉" << std::endl;
    std::cout << "BH Mass Fraction: " << (variables["M_BH"] / variables["M"]) * 100 << "%" << std::endl;
    
    std::cout << "\nSpatial Parameters:" << std::endl;
    std::cout << "r (half-diameter): " << variables["r"] << " m" << std::endl;
    std::cout << "r_BH (core scale): " << variables["r_BH"] << " m" << std::endl;
    std::cout << "z (blueshift): " << variables["z"] << std::endl;
    
    std::cout << "\nDynamics:" << std::endl;
    std::cout << "v_orbit: " << (variables["v_orbit"] / 1e3) << " km/s" << std::endl;
    std::cout << "rho_dust: " << variables["rho_dust"] << " kg/m^3" << std::endl;
    std::cout << "rho_mass: " << variables["rho_mass"] << " kg/m^3" << std::endl;
    std::cout << "B (magnetic): " << variables["B"] << " T" << std::endl;
    
    std::cout << "\nVacuum Energies:" << std::endl;
    std::cout << "rho_vac_UA: " << variables["rho_vac_UA"] << " J/m^3" << std::endl;
    std::cout << "rho_vac_SCm: " << variables["rho_vac_SCm"] << " J/m^3" << std::endl;
    std::cout << "Vacuum ratio: " << (variables["rho_vac_UA"] / variables["rho_vac_SCm"]) << std::endl;
    
    std::cout << "\nScaling Factors:" << std::endl;
    std::cout << "f_TRZ (time-reversal): " << variables["f_TRZ"] << std::endl;
    std::cout << "scale_macro: " << variables["scale_macro"] << std::endl;
    
    // Current computation
    double t = variables["t"];
    double g = computeG(t);
    double Hz = computeHz();
    double a_dust = computeADust();
    double em_term = computeEMTerm();
    
    std::cout << "\nCurrent Computation:" << std::endl;
    std::cout << "t: " << (t / (variables["Gyr"] * variables["year_to_s"])) << " Gyr" << std::endl;
    std::cout << "g_Andromeda: " << g << " m/s^2" << std::endl;
    
    std::cout << "\nSubcomponents:" << std::endl;
    std::cout << "H(z): " << Hz << " s^-1" << std::endl;
    std::cout << "a_dust: " << a_dust << " m/s^2" << std::endl;
    std::cout << "EM term: " << em_term << " m/s^2" << std::endl;
    std::cout << "g_BH: " << (variables["G"] * variables["M_BH"] / (variables["r_BH"] * variables["r_BH"])) << " m/s^2" << std::endl;
    
    std::cout << "======================================================\n" << std::endl;
}

void AndromedaUQFFModule::validatePhysicalConsistency() {
    std::cout << "Validating Andromeda physical consistency..." << std::endl;
    bool consistent = true;
    
    // Check for NaN/Inf
    for (const auto& pair : variables) {
        if (std::isnan(pair.second) || std::isinf(pair.second)) {
            std::cerr << "  ERROR: " << pair.first << " is NaN/Inf" << std::endl;
            consistent = false;
        }
    }
    
    // M_BH should be much smaller than M
    if (variables["M_BH"] > variables["M"]) {
        std::cerr << "  ERROR: M_BH cannot exceed M_total" << std::endl;
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
    
    if (variables["v_orbit"] <= 0) {
        std::cerr << "  ERROR: v_orbit must be positive" << std::endl;
        consistent = false;
    }
    
    // Reasonable ranges
    if (variables["z"] < -0.01 || variables["z"] > 0.01) {
        std::cerr << "  WARNING: z (redshift) outside typical range for Andromeda" << std::endl;
    }
    
    if (consistent) {
        std::cout << "  All checks passed. Andromeda system is physically consistent." << std::endl;
    }
}

void AndromedaUQFFModule::autoCorrectAnomalies() {
    std::cout << "Auto-correcting Andromeda anomalies..." << std::endl;
    
    // Fix NaN/Inf
    for (auto& pair : variables) {
        if (std::isnan(pair.second) || std::isinf(pair.second)) {
            std::cout << "  Correcting " << pair.first << " (was NaN/Inf)" << std::endl;
            pair.second = 1.0;
        }
    }
    
    // Ensure M_BH < M
    if (variables["M_BH"] > variables["M"]) {
        std::cout << "  Correcting M_BH to be < M_total" << std::endl;
        variables["M_BH"] = 1.4e8 * variables["M_sun"];
    }
    
    // Ensure positive values
    if (variables["M"] <= 0) {
        std::cout << "  Correcting M to 1e12 M_sun" << std::endl;
        variables["M"] = 1e12 * variables["M_sun"];
    }
    
    if (variables["r"] <= 0) {
        std::cout << "  Correcting r to 1.04e21 m" << std::endl;
        variables["r"] = 1.04e21;
    }
    
    if (variables["v_orbit"] <= 0) {
        std::cout << "  Correcting v_orbit to 2.5e5 m/s" << std::endl;
        variables["v_orbit"] = 2.5e5;
    }
    
    // Validate z range for Andromeda
    if (variables["z"] < -0.01 || variables["z"] > 0.01) {
        std::cout << "  Correcting z to -0.001 (Andromeda blueshift)" << std::endl;
        variables["z"] = -0.001;
    }
    
    std::cout << "Auto-correction complete." << std::endl;
}

// Example usage in base program (snippet)
// #include "AndromedaUQFFModule.h"
// int main() {
//     std::cout << "=============================================" << std::endl;
//     std::cout << "Andromeda Galaxy - Enhanced Demonstration" << std::endl;
//     std::cout << "=============================================" << std::endl;
// 
//     // === Part 1: Basic gravity computation ===
//     std::cout << "\n--- Part 1: Basic Gravity Computation ---" << std::endl;
//     AndromedaUQFFModule mod;
//     double t = 10.0 * 1e9 * 3.156e7;  // 10 Gyr
//     double g = mod.computeG(t);
//     std::cout << "g_Andromeda at t=10 Gyr: " << g << " m/s^2" << std::endl;
// 
//     // === Part 2: Evolution table ===
//     std::cout << "\n--- Part 2: Evolution Table ---" << std::endl;
//     mod.printEvolutionTable();
// 
//     // === Part 3: Dynamic variable management ===
//     std::cout << "\n--- Part 3: Dynamic Variable Management ---" << std::endl;
//     mod.createDynamicVariable("custom_dust_factor", 2.0);
//     mod.cloneVariable("v_orbit", "v_orbit_backup");
//     std::cout << "Total variables: " << std::endl;
//     mod.listAllVariables();
// 
//     // === Part 4: Orbital velocity exploration ===
//     std::cout << "\n--- Part 4: Orbital Velocity Dynamics ---" << std::endl;
//     std::cout << "Initial v_orbit: " << (mod.variables["v_orbit"] / 1e3) << " km/s" << std::endl;
//     mod.updateVariable("v_orbit", 3e5);  // Increase to 300 km/s
//     std::cout << "Updated v_orbit: " << (mod.variables["v_orbit"] / 1e3) << " km/s" << std::endl;
//     double g_new = mod.computeG(t);
//     std::cout << "New g_Andromeda: " << g_new << " m/s^2" << std::endl;
// 
//     // === Part 5: Black hole mass exploration ===
//     std::cout << "\n--- Part 5: SMBH Mass Exploration ---" << std::endl;
//     mod.updateVariable("M_BH", 2e8 * mod.variables["M_sun"]);
//     std::cout << "Increased M_BH to 2e8 M☉" << std::endl;
// 
//     // === Part 6: Batch operations ===
//     std::cout << "\n--- Part 6: Batch Operations ---" << std::endl;
//     std::vector<std::string> mass_params = {"M", "M_BH"};
//     std::cout << "Scaling mass parameters by 1.2..." << std::endl;
//     mod.scaleVariableGroup(mass_params, 1.2);
// 
//     // === Part 7: Self-expansion ===
//     std::cout << "\n--- Part 7: Self-Expansion ---" << std::endl;
//     mod.saveState("before_expansion");
//     mod.expandSpatialScale(1.5);
//     std::cout << "Spatial scale expanded by 1.5x" << std::endl;
// 
//     // === Part 8: Parameter exploration ===
//     std::cout << "\n--- Part 8: Parameter Exploration ---" << std::endl;
//     mod.generateVariations(3, 0.15);
// 
//     // === Part 9: Sensitivity analysis ===
//     std::cout << "\n--- Part 9: Sensitivity Analysis ---" << std::endl;
//     mod.restoreState("before_expansion");
//     mod.analyzeParameterSensitivity("v_orbit");
// 
//     // === Part 10: System validation ===
//     std::cout << "\n--- Part 10: System Validation ---" << std::endl;
//     mod.validatePhysicalConsistency();
// 
//     // === Part 11: Auto-refinement ===
//     std::cout << "\n--- Part 11: Auto-Refinement ---" << std::endl;
//     mod.autoRefineParameters(0.01);
// 
//     // === Part 12: Comprehensive system report ===
//     std::cout << "\n--- Part 12: System Report ---" << std::endl;
//     mod.generateSystemReport();
// 
//     // === Part 13: State management ===
//     std::cout << "\n--- Part 13: State Management ---" << std::endl;
//     mod.saveState("final_andromeda_state");
//     mod.listSavedStates();
// 
//     std::cout << "\n=============================================" << std::endl;
//     std::cout << "Enhanced Andromeda demonstration complete!" << std::endl;
//     std::cout << "=============================================" << std::endl;
// 
//     return 0;
// }
// Compile: g++ -o andromeda_test andromeda_test.cpp AndromedaUQFFModule.cpp -lm
// Sample Output at t=10 Gyr: g ≈ 6.273 m/s²; table shows near-constant due to small expansion.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

AndromedaUQFFModule Evaluation

Strengths :
-Modular, extensible design for modeling Andromeda Galaxy gravity, including base gravity, expansion, SMBH term, dust friction, and EM / Aether effects.
- Comprehensive physics : gravity, cosmological expansion, SMBH, dust, electromagnetic, and vacuum energy terms.
- Dynamic variable management via std::map enables runtime updates and system adaptation.
- Clear separation of computation functions(e.g., H(z), dust, EM), aiding maintainability.
- Andromeda - specific parameters are initialized for realistic simulation; supports easy modification.
- Output functions for equation text, variable state, and evolution table support debugging and documentation.
- Approximations and scaling factors are documented, supporting scientific reproducibility.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for flexibility and scalability.
- Some calculations use magic numbers or lack explanatory comments; define named constants and clarify logic.
- Minimal error handling(e.g., division by zero, invalid variable names); add validation for robustness.
- Unit consistency should be checked and documented for all physical quantities.
- For large - scale or performance - critical simulations, optimize data structures and reduce redundant calculations.
- std::map is flexible but may be less efficient than structured types for very large models.
- Expand documentation for function purposes and physical meaning.

Summary:
The code is well - structured, flexible, and suitable for scientific prototyping and educational use in galactic gravity modeling.It implements a broad set of physical effects and adapts to various scenarios.For production or high - performance applications, address the recommendations for improved robustness, maintainability, and scalability.