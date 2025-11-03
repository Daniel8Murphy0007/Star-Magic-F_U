// Ug1DefectModule.h
// Modular C++ implementation of the Ug1 Defect Factor (?_def) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes ?_def = 0.01 * sin(0.001 t) (unitless); scales (1 + ?_def) in Universal Gravity U_g1 term.
// Pluggable: #include "Ug1DefectModule.h"
// Ug1DefectModule mod; mod.computeU_g1(0.0, 1.496e11); mod.updateVariable("amplitude", new_value);
// Variables in std::map; example for Sun at t=0 (?_def=0, U_g1?4.51e31 J/m�); t=1570.8 days: +1%.
// Approximations: ?=0.001 day?�; cos(? t_n)=1 at t_n=0; ?_s=3.38e23 T�m�; ?(M_s/r)?M_s/r�=8.89e7 m/s�.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef UG1_DEFECT_MODULE_H
#define UG1_DEFECT_MODULE_H

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

class Ug1DefectModule {
private:
    std::map<std::string, double> variables;
    double computeDelta_def(double t_day);
    double computeU_g1(double t_day, double r);

public:
    // Constructor: Initialize with framework defaults
    Ug1DefectModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeDelta_def(double t_day);  // 0.01 * sin(0.001 t)
    double computeU_g1(double t_day, double r);  // U_g1 with defect (J/m^3)
    double computePeriod_years();  // ~17.22 years

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // ========== ENHANCED SELF-UPDATE & SELF-EXPANSION CAPABILITIES ==========
    
    // Variable Management
    void createVariable(const std::string& name, double value);
    void removeVariable(const std::string& name);
    void cloneVariable(const std::string& source, const std::string& target);
    std::vector<std::string> listVariables() const;
    std::string getSystemName() const;
    
    // Batch Operations
    void transformVariableGroup(const std::vector<std::string>& varNames, 
                                std::function<double(double)> transformFunc);
    void scaleVariableGroup(const std::vector<std::string>& varNames, double factor);
    
    // Self-Expansion Methods (Domain-Specific for Ug1 Defect)
    void expandParameterSpace(double expansion_factor);
    void expandDefectScale(double amplitude_factor, double freq_factor);
    void expandGravityScale(double ug1_factor, double coupling_factor);
    void expandTemporalScale(double period_factor, double decay_factor);
    
    // Self-Refinement
    void autoRefineParameters();
    void calibrateToObservations(const std::map<std::string, double>& observed_data);
    void optimizeForMetric(const std::string& metric);
    
    // Parameter Exploration
    std::vector<std::map<std::string, double>> generateVariations(int count);
    
    // Adaptive Evolution
    void mutateParameters(double mutation_rate);
    void evolveSystem(int generations);
    
    // State Management
    void saveState(const std::string& label);
    void restoreState(const std::string& label);
    std::vector<std::string> listSavedStates() const;
    std::string exportState() const;
    
    // System Analysis
    std::map<std::string, double> sensitivityAnalysis(const std::string& output_var);
    std::string generateReport() const;
    bool validateConsistency() const;
    void autoCorrectAnomalies();
};

#endif // UG1_DEFECT_MODULE_H

// Ug1DefectModule.cpp
#include "Ug1DefectModule.h"

// Constructor: Set framework defaults (Sun)
Ug1DefectModule::Ug1DefectModule() {
    // Universal constants
    variables["amplitude"] = 0.01;                  // Unitless
    variables["freq"] = 0.001;                      // day?�
    variables["k_1"] = 1.5;                         // Coupling
    variables["mu_s"] = 3.38e23;                    // T�m�
    variables["M_s"] = 1.989e30;                    // kg
    variables["alpha"] = 0.001;                     // day?�
    variables["t_n"] = 0.0;                         // days
    variables["pi"] = 3.141592653589793;
    variables["t_day"] = 0.0;                       // days
    variables["r"] = 1.496e11;                      // m (Earth-Sun example)

    // Derived
    variables["period_days"] = 2.0 * M_PI / variables["freq"];
}

// Update variable
void Ug1DefectModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
        if (name == "freq") {
            variables["period_days"] = 2.0 * M_PI / value;
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void Ug1DefectModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        if (name == "freq") {
            variables["period_days"] = 2.0 * M_PI / variables[name];
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void Ug1DefectModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute ?_def = 0.01 * sin(0.001 t) (t in days)
double Ug1DefectModule::computeDelta_def(double t_day) {
    variables["t_day"] = t_day;
    return variables["amplitude"] * std::sin(variables["freq"] * t_day);
}

// Compute U_g1 = k_1 * ?_s * ?(M_s / r) * exp(-? t) * cos(? t_n) * (1 + ?_def)
double Ug1DefectModule::computeU_g1(double t_day, double r) {
    variables["r"] = r;
    double k_1 = variables["k_1"];
    double mu_s = variables["mu_s"];
    double grad_ms_r = variables["M_s"] / (r * r);  // Approx ?(M_s / r) = M_s / r^2
    double exp_term = std::exp( - variables["alpha"] * t_day );
    double cos_tn = std::cos(variables["pi"] * variables["t_n"]);
    double defect_factor = 1.0 + computeDelta_def(t_day);
    return k_1 * mu_s * grad_ms_r * exp_term * cos_tn * defect_factor;
}

// Period in years (365.25 days/year)
double Ug1DefectModule::computePeriod_years() {
    return variables["period_days"] / 365.25;
}

// Equation text
std::string Ug1DefectModule::getEquationText() {
    return "U_g1 = k_1 * ?_s * ?(M_s / r) * e^{-? t} * cos(? t_n) * (1 + ?_def)\n"
           "Where ?_def = 0.01 * sin(0.001 t) (unitless, t days; period ~17.22 yr).\n"
           "Small oscillatory defect (~�1%) in internal dipole gravity.\n"
           "Example t=0, r=1.496e11 m: ?_def=0, U_g1 ?4.51e31 J/m�;\n"
           "t=1570.8 days: ?_def=0.01, U_g1 ?4.56e31 J/m� (+1.1%).\n"
           "Role: Time-dependent perturbations; internal dynamics/[SCm] variations.\n"
           "UQFF: Cyclic defects in stellar gravity; for formation/nebular stability.";
}

// Print variables
void Ug1DefectModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// ========== ENHANCED SELF-UPDATE & SELF-EXPANSION IMPLEMENTATION ==========

// Variable Management
void Ug1DefectModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
    std::cout << "Created variable '" << name << "' = " << value << std::endl;
}

void Ug1DefectModule::removeVariable(const std::string& name) {
    if (variables.erase(name)) {
        std::cout << "Removed variable '" << name << "'" << std::endl;
    } else {
        std::cerr << "Variable '" << name << "' not found for removal" << std::endl;
    }
}

void Ug1DefectModule::cloneVariable(const std::string& source, const std::string& target) {
    if (variables.find(source) != variables.end()) {
        variables[target] = variables[source];
        std::cout << "Cloned '" << source << "' to '" << target << "'" << std::endl;
    } else {
        std::cerr << "Source variable '" << source << "' not found" << std::endl;
    }
}

std::vector<std::string> Ug1DefectModule::listVariables() const {
    std::vector<std::string> varList;
    for (const auto& pair : variables) {
        varList.push_back(pair.first);
    }
    return varList;
}

std::string Ug1DefectModule::getSystemName() const {
    return "Ug1_Defect_Factor_UQFF";
}

// Batch Operations
void Ug1DefectModule::transformVariableGroup(const std::vector<std::string>& varNames, 
                                              std::function<double(double)> transformFunc) {
    for (const auto& name : varNames) {
        if (variables.find(name) != variables.end()) {
            variables[name] = transformFunc(variables[name]);
        }
    }
    std::cout << "Transformed " << varNames.size() << " variables" << std::endl;
}

void Ug1DefectModule::scaleVariableGroup(const std::vector<std::string>& varNames, double factor) {
    transformVariableGroup(varNames, [factor](double val) { return val * factor; });
}

// Self-Expansion Methods
void Ug1DefectModule::expandParameterSpace(double expansion_factor) {
    // Generic expansion - scale all physics parameters
    std::vector<std::string> physics_vars = {"amplitude", "freq", "k_1", "mu_s", "alpha"};
    scaleVariableGroup(physics_vars, expansion_factor);
    
    // Recalculate derived quantities
    variables["period_days"] = 2.0 * M_PI / variables["freq"];
    
    std::cout << "Expanded parameter space by factor " << expansion_factor << std::endl;
}

void Ug1DefectModule::expandDefectScale(double amplitude_factor, double freq_factor) {
    // Defect-specific expansion: amplitude and frequency of δ_def oscillation
    variables["amplitude"] *= amplitude_factor;
    variables["freq"] *= freq_factor;
    variables["period_days"] = 2.0 * M_PI / variables["freq"];
    
    // Period in years for reference
    double period_years = variables["period_days"] / 365.25;
    
    std::cout << "Expanded defect scale: amplitude=" << variables["amplitude"] 
              << ", freq=" << variables["freq"] << " day⁻¹, period=" 
              << period_years << " years\n";
}

void Ug1DefectModule::expandGravityScale(double ug1_factor, double coupling_factor) {
    // U_g1 gravity expansion: k_1 coupling and μ_s magnetic strength
    variables["k_1"] *= coupling_factor;
    variables["mu_s"] *= ug1_factor;
    
    std::cout << "Expanded gravity scale: k_1=" << variables["k_1"] 
              << ", mu_s=" << variables["mu_s"] << " T²m³\n";
}

void Ug1DefectModule::expandTemporalScale(double period_factor, double decay_factor) {
    // Temporal expansion: oscillation period and decay rate
    variables["freq"] /= period_factor;  // Longer period = lower frequency
    variables["alpha"] *= decay_factor;  // Decay rate
    variables["period_days"] = 2.0 * M_PI / variables["freq"];
    
    double period_years = variables["period_days"] / 365.25;
    
    std::cout << "Expanded temporal scale: period=" << period_years 
              << " years, alpha=" << variables["alpha"] << " day⁻¹\n";
}

// Self-Refinement
void Ug1DefectModule::autoRefineParameters() {
    // Auto-refine to match known physics constraints
    // Amplitude typically in range [0.001, 0.05] for 0.1-5% defect
    if (variables["amplitude"] > 0.05) {
        variables["amplitude"] = 0.05;
        std::cout << "Refined amplitude to physical limit (0.05)\n";
    }
    if (variables["amplitude"] < 0.0) {
        variables["amplitude"] = 0.001;
        std::cout << "Refined amplitude to minimum (0.001)\n";
    }
    
    // Frequency should be small (long periods, ~10-20 years typical)
    if (variables["freq"] < 0.0) {
        variables["freq"] = 0.0001;
        std::cout << "Refined freq to minimum (0.0001 day⁻¹)\n";
    }
    
    // Ensure k_1 is positive
    if (variables["k_1"] <= 0.0) {
        variables["k_1"] = 1.5;
        std::cout << "Refined k_1 to default (1.5)\n";
    }
    
    // Recalculate derived quantities
    variables["period_days"] = 2.0 * M_PI / variables["freq"];
    
    std::cout << "Auto-refinement complete\n";
}

void Ug1DefectModule::calibrateToObservations(const std::map<std::string, double>& observed_data) {
    // Calibrate parameters to match observational data
    for (const auto& obs : observed_data) {
        if (variables.find(obs.first) != variables.end()) {
            double old_val = variables[obs.first];
            variables[obs.first] = obs.second;
            std::cout << "Calibrated '" << obs.first << "': " << old_val << " → " << obs.second << "\n";
        }
    }
    
    // Update derived quantities
    if (observed_data.find("freq") != observed_data.end()) {
        variables["period_days"] = 2.0 * M_PI / variables["freq"];
    }
}

void Ug1DefectModule::optimizeForMetric(const std::string& metric) {
    // Optimize parameters for specific physical scenarios
    if (metric == "solar_cycle") {
        // 11-year solar cycle (~4000 days)
        variables["amplitude"] = 0.01;
        variables["freq"] = 2.0 * M_PI / 4015.0;  // ~11 years
        variables["alpha"] = 0.0001;  // Slow decay
        std::cout << "Optimized for solar cycle (11 years)\n";
        
    } else if (metric == "long_period") {
        // ~20-year period for long-term stellar variations
        variables["amplitude"] = 0.02;
        variables["freq"] = 2.0 * M_PI / 7305.0;  // ~20 years
        variables["alpha"] = 0.0001;
        std::cout << "Optimized for long period (20 years)\n";
        
    } else if (metric == "short_period") {
        // ~5-year period for faster variations
        variables["amplitude"] = 0.005;
        variables["freq"] = 2.0 * M_PI / 1826.0;  // ~5 years
        variables["alpha"] = 0.001;
        std::cout << "Optimized for short period (5 years)\n";
        
    } else if (metric == "high_defect") {
        // Large amplitude defects (~5%)
        variables["amplitude"] = 0.05;
        variables["freq"] = 0.001;
        variables["alpha"] = 0.0001;
        std::cout << "Optimized for high defect amplitude (5%)\n";
        
    } else if (metric == "stable") {
        // Minimal defect, very stable
        variables["amplitude"] = 0.001;
        variables["freq"] = 0.0001;
        variables["alpha"] = 0.0001;
        std::cout << "Optimized for stable/minimal defect\n";
        
    } else {
        std::cout << "Unknown metric '" << metric << "', no optimization applied\n";
        return;
    }
    
    // Update derived quantities
    variables["period_days"] = 2.0 * M_PI / variables["freq"];
}

// Parameter Exploration
std::vector<std::map<std::string, double>> Ug1DefectModule::generateVariations(int count) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.8, 1.2);  // ±20% variation
    
    for (int i = 0; i < count; ++i) {
        std::map<std::string, double> variant = variables;
        
        // Vary key parameters
        variant["amplitude"] *= dis(gen);
        variant["freq"] *= dis(gen);
        variant["k_1"] *= dis(gen);
        variant["alpha"] *= dis(gen);
        
        // Ensure physical constraints
        if (variant["amplitude"] < 0.0) variant["amplitude"] = 0.001;
        if (variant["amplitude"] > 0.05) variant["amplitude"] = 0.05;
        if (variant["freq"] < 0.0) variant["freq"] = 0.0001;
        
        // Update period
        variant["period_days"] = 2.0 * M_PI / variant["freq"];
        
        variations.push_back(variant);
    }
    
    std::cout << "Generated " << count << " parameter variations\n";
    return variations;
}

// Adaptive Evolution
void Ug1DefectModule::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-mutation_rate, mutation_rate);
    
    // Mutate key parameters
    variables["amplitude"] *= (1.0 + dis(gen));
    variables["freq"] *= (1.0 + dis(gen));
    variables["k_1"] *= (1.0 + dis(gen));
    variables["alpha"] *= (1.0 + dis(gen));
    
    // Apply constraints
    autoRefineParameters();
    
    std::cout << "Mutated parameters with rate " << mutation_rate << "\n";
}

void Ug1DefectModule::evolveSystem(int generations) {
    double best_fitness = 0.0;
    std::map<std::string, double> best_params = variables;
    
    for (int gen = 0; gen < generations; ++gen) {
        // Generate variations
        auto variations = generateVariations(5);
        
        // Evaluate fitness (example: balance defect amplitude with period length)
        for (const auto& variant : variations) {
            double amp = variant.at("amplitude");
            double period_years = variant.at("period_days") / 365.25;
            
            // Fitness: moderate amplitude (1-2%) with reasonable period (10-20 years)
            double fitness = 1.0 / (1.0 + std::abs(amp - 0.015)) * 
                            1.0 / (1.0 + std::abs(period_years - 15.0) / 10.0);
            
            if (fitness > best_fitness) {
                best_fitness = fitness;
                best_params = variant;
            }
        }
        
        // Mutate current best
        variables = best_params;
        mutateParameters(0.05);
    }
    
    variables = best_params;
    autoRefineParameters();
    
    std::cout << "Evolved system over " << generations << " generations (fitness=" 
              << best_fitness << ")\n";
}

// State Management
static std::map<std::string, std::map<std::string, double>> saved_states_defect;

void Ug1DefectModule::saveState(const std::string& label) {
    saved_states_defect[label] = variables;
    std::cout << "Saved state '" << label << "'\n";
}

void Ug1DefectModule::restoreState(const std::string& label) {
    if (saved_states_defect.find(label) != saved_states_defect.end()) {
        variables = saved_states_defect[label];
        std::cout << "Restored state '" << label << "'\n";
    } else {
        std::cerr << "State '" << label << "' not found\n";
    }
}

std::vector<std::string> Ug1DefectModule::listSavedStates() const {
    std::vector<std::string> state_list;
    for (const auto& pair : saved_states_defect) {
        state_list.push_back(pair.first);
    }
    return state_list;
}

std::string Ug1DefectModule::exportState() const {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(10);
    for (const auto& pair : variables) {
        oss << pair.first << "=" << pair.second << "\n";
    }
    return oss.str();
}

// System Analysis
std::map<std::string, double> Ug1DefectModule::sensitivityAnalysis(const std::string& output_var) {
    std::map<std::string, double> sensitivities;
    double baseline_output = 0.0;
    
    // Compute baseline (example: U_g1 at t=0, r=1.496e11)
    if (output_var == "U_g1") {
        baseline_output = computeU_g1(0.0, variables["r"]);
    } else if (output_var == "delta_def") {
        baseline_output = computeDelta_def(1570.8);  // At peak
    } else {
        std::cerr << "Unknown output variable '" << output_var << "'\n";
        return sensitivities;
    }
    
    // Test sensitivity to each parameter
    double delta = 0.01;  // 1% perturbation
    std::vector<std::string> params = {"amplitude", "freq", "k_1", "mu_s", "alpha"};
    
    for (const auto& param : params) {
        double original_val = variables[param];
        
        // Perturb upward
        variables[param] = original_val * (1.0 + delta);
        if (param == "freq") {
            variables["period_days"] = 2.0 * M_PI / variables["freq"];
        }
        
        double perturbed_output = 0.0;
        if (output_var == "U_g1") {
            perturbed_output = computeU_g1(0.0, variables["r"]);
        } else if (output_var == "delta_def") {
            perturbed_output = computeDelta_def(1570.8);
        }
        
        // Compute sensitivity (normalized)
        double sensitivity = (perturbed_output - baseline_output) / (baseline_output * delta);
        sensitivities[param] = sensitivity;
        
        // Restore original value
        variables[param] = original_val;
        if (param == "freq") {
            variables["period_days"] = 2.0 * M_PI / variables["freq"];
        }
    }
    
    return sensitivities;
}

std::string Ug1DefectModule::generateReport() const {
    std::ostringstream report;
    report << std::scientific << std::setprecision(3);
    
    report << "========== UG1 DEFECT MODULE REPORT ==========\n";
    report << "System: " << getSystemName() << "\n\n";
    
    report << "Defect Parameters:\n";
    report << "  Amplitude = " << variables.at("amplitude") << " (unitless)\n";
    report << "  Frequency = " << variables.at("freq") << " day⁻¹\n";
    report << "  Period = " << variables.at("period_days") << " days (" 
           << (variables.at("period_days") / 365.25) << " years)\n";
    report << "  Decay rate α = " << variables.at("alpha") << " day⁻¹\n\n";
    
    report << "Gravity Parameters:\n";
    report << "  k_1 = " << variables.at("k_1") << " (coupling)\n";
    report << "  μ_s = " << variables.at("mu_s") << " T²m³\n";
    report << "  M_s = " << variables.at("M_s") << " kg\n\n";
    
    // Compute δ_def at key times
    Ug1DefectModule temp_mod = *this;
    double delta_0 = temp_mod.computeDelta_def(0.0);
    double delta_quarter = temp_mod.computeDelta_def(variables.at("period_days") / 4.0);
    double delta_half = temp_mod.computeDelta_def(variables.at("period_days") / 2.0);
    
    report << "Defect Factor δ_def:\n";
    report << "  At t=0: " << std::fixed << std::setprecision(4) << delta_0 << "\n";
    report << "  At quarter period: " << delta_quarter << " (peak)\n";
    report << "  At half period: " << delta_half << "\n";
    report << std::scientific << std::setprecision(3) << "\n";
    
    // Compute U_g1 at Earth distance
    double r_earth = 1.496e11;
    double ug1_t0 = temp_mod.computeU_g1(0.0, r_earth);
    double ug1_peak = temp_mod.computeU_g1(variables.at("period_days") / 4.0, r_earth);
    double enhancement = ((ug1_peak - ug1_t0) / ug1_t0) * 100.0;
    
    report << "U_g1 at Earth Distance (r=" << r_earth << " m):\n";
    report << "  At t=0: " << ug1_t0 << " J/m³\n";
    report << "  At peak defect: " << ug1_peak << " J/m³\n";
    report << "  Variation: " << std::fixed << std::setprecision(1) << enhancement << "%\n";
    report << std::scientific << std::setprecision(3) << "\n";
    
    report << "Physical Context:\n";
    report << "  Oscillatory defect in U_g1 internal dipole gravity\n";
    report << "  Origin: [SCm] density variations, internal dynamics\n";
    report << "  Applications: Stellar irregularities, sunspot cycles\n";
    report << "  UQFF role: Time-dependent U_g1 perturbations\n";
    
    report << "==============================================\n";
    
    return report.str();
}

bool Ug1DefectModule::validateConsistency() const {
    bool consistent = true;
    
    // Check amplitude range [0.001, 0.05]
    if (variables.at("amplitude") < 0.0 || variables.at("amplitude") > 0.05) {
        std::cerr << "Inconsistency: amplitude out of physical range [0.0, 0.05]\n";
        consistent = false;
    }
    
    // Check positive frequency
    if (variables.at("freq") <= 0.0) {
        std::cerr << "Inconsistency: freq must be positive\n";
        consistent = false;
    }
    
    // Check positive k_1
    if (variables.at("k_1") <= 0.0) {
        std::cerr << "Inconsistency: k_1 must be positive\n";
        consistent = false;
    }
    
    // Check positive mu_s
    if (variables.at("mu_s") <= 0.0) {
        std::cerr << "Inconsistency: mu_s must be positive\n";
        consistent = false;
    }
    
    // Check derived quantity consistency
    double expected_period = 2.0 * M_PI / variables.at("freq");
    double actual_period = variables.at("period_days");
    if (std::abs(expected_period - actual_period) / expected_period > 0.01) {
        std::cerr << "Inconsistency: period_days mismatch\n";
        consistent = false;
    }
    
    if (consistent) {
        std::cout << "Consistency validation: PASSED\n";
    }
    
    return consistent;
}

void Ug1DefectModule::autoCorrectAnomalies() {
    bool corrected = false;
    
    // Correct amplitude range
    if (variables["amplitude"] < 0.0) {
        variables["amplitude"] = 0.001;
        corrected = true;
    }
    if (variables["amplitude"] > 0.05) {
        variables["amplitude"] = 0.05;
        corrected = true;
    }
    
    // Correct frequency
    if (variables["freq"] <= 0.0) {
        variables["freq"] = 0.001;
        corrected = true;
    }
    
    // Correct k_1
    if (variables["k_1"] <= 0.0) {
        variables["k_1"] = 1.5;
        corrected = true;
    }
    
    // Correct mu_s
    if (variables["mu_s"] <= 0.0) {
        variables["mu_s"] = 3.38e23;
        corrected = true;
    }
    
    // Recalculate derived quantities
    variables["period_days"] = 2.0 * M_PI / variables["freq"];
    
    if (corrected) {
        std::cout << "Auto-corrected anomalies\n";
    } else {
        std::cout << "No anomalies found\n";
    }
}

// Example usage in base program (snippet)
// #include "Ug1DefectModule.h"
// int main() {
//     Ug1DefectModule mod;
//     double delta = mod.computeDelta_def(0.0);
//     std::cout << "δ_def (t=0) = " << delta << std::endl;
//     double u_g1 = mod.computeU_g1(1570.8, 1.496e11);
//     std::cout << "U_g1 (t=1570.8 days) = " << u_g1 << " J/m³\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("amplitude", 0.02);
//     mod.printVariables();
//     return 0;
// }

// ========== COMPREHENSIVE ENHANCED DEMONSTRATION ==========
/*
int main() {
    std::cout << std::scientific << std::setprecision(3);
    std::cout << "========== UG1 DEFECT MODULE DEMONSTRATION ==========\n\n";
    
    // ===== Step 1: Initialize Module =====
    Ug1DefectModule mod;
    std::cout << "Step 1: Module initialized with defaults:\n";
    std::cout << "  System: " << mod.getSystemName() << "\n";
    std::cout << "  Amplitude = " << mod.variables["amplitude"] << "\n";
    std::cout << "  Frequency = " << mod.variables["freq"] << " day⁻¹\n";
    std::cout << "  Period = " << mod.computePeriod_years() << " years\n\n";
    
    // ===== Step 2: Baseline Computation =====
    std::cout << "Step 2: Compute baseline δ_def and U_g1:\n";
    double r_earth = 1.496e11;  // meters
    
    double delta_t0 = mod.computeDelta_def(0.0);
    double ug1_t0 = mod.computeU_g1(0.0, r_earth);
    std::cout << "  At t=0 days:\n";
    std::cout << "    δ_def = " << std::fixed << std::setprecision(4) << delta_t0 << "\n";
    std::cout << "    U_g1 = " << std::scientific << std::setprecision(3) << ug1_t0 << " J/m³\n";
    
    double t_quarter = mod.variables["period_days"] / 4.0;
    double delta_peak = mod.computeDelta_def(t_quarter);
    double ug1_peak = mod.computeU_g1(t_quarter, r_earth);
    std::cout << "  At quarter period (" << std::fixed << std::setprecision(1) 
              << t_quarter << " days, peak defect):\n";
    std::cout << "    δ_def = " << std::setprecision(4) << delta_peak << "\n";
    std::cout << "    U_g1 = " << std::scientific << std::setprecision(3) << ug1_peak << " J/m³\n";
    std::cout << "    Enhancement: +" << std::fixed << std::setprecision(1) 
              << ((ug1_peak/ug1_t0 - 1.0)*100) << "%\n\n";
    
    // ===== Step 3: Variable Management =====
    std::cout << "Step 3: Variable Management\n";
    std::cout << std::scientific << std::setprecision(3);
    mod.createVariable("max_defect_amplitude", 0.05);
    std::cout << "  Created 'max_defect_amplitude' = 0.05 (5% limit)\n";
    
    std::vector<std::string> all_vars = mod.listVariables();
    std::cout << "  Total variables: " << all_vars.size() << "\n";
    
    mod.cloneVariable("amplitude", "amplitude_backup");
    std::cout << "  Cloned 'amplitude' → 'amplitude_backup'\n\n";
    
    // ===== Step 4: Defect Expansion (Larger Amplitude) =====
    std::cout << "Step 4: Defect Expansion (Increase Amplitude & Frequency)\n";
    mod.expandDefectScale(2.0, 1.5);  // 2x amplitude, 1.5x frequency
    std::cout << "  After expansion:\n";
    std::cout << "    Amplitude = " << mod.variables["amplitude"] << "\n";
    std::cout << "    Frequency = " << mod.variables["freq"] << " day⁻¹\n";
    std::cout << "    Period = " << mod.computePeriod_years() << " years\n";
    double delta_expanded = mod.computeDelta_def(t_quarter);
    std::cout << "    Peak δ_def = " << std::fixed << std::setprecision(4) << delta_expanded << "\n\n";
    
    // ===== Step 5: Gravity Expansion (Enhanced U_g1) =====
    std::cout << "Step 5: Gravity Expansion (Increase k_1 & μ_s)\n";
    std::cout << std::scientific << std::setprecision(3);
    mod.expandGravityScale(1.5, 1.3);  // 1.5x μ_s, 1.3x k_1
    double ug1_expanded = mod.computeU_g1(0.0, r_earth);
    std::cout << "  After expansion: U_g1 = " << ug1_expanded << " J/m³\n\n";
    
    // ===== Step 6: Temporal Expansion (Longer Period) =====
    std::cout << "Step 6: Temporal Expansion (Increase Period)\n";
    mod.expandTemporalScale(2.0, 1.2);  // 2x period, 1.2x decay
    std::cout << "  After temporal expansion:\n";
    std::cout << "    Period = " << mod.computePeriod_years() << " years\n";
    std::cout << "    Frequency = " << mod.variables["freq"] << " day⁻¹\n";
    std::cout << "    Decay α = " << mod.variables["alpha"] << " day⁻¹\n\n";
    
    // ===== Step 7: Batch Operations (Reset Defect Group) =====
    std::cout << "Step 7: Batch Operations (Scale Defect Parameters)\n";
    std::vector<std::string> defect_group = {"amplitude", "freq"};
    mod.scaleVariableGroup(defect_group, 0.5);  // Reduce to moderate
    std::cout << "  Scaled defect group by 0.5:\n";
    std::cout << "    Amplitude = " << mod.variables["amplitude"] << "\n";
    std::cout << "    Frequency = " << mod.variables["freq"] << " day⁻¹\n";
    std::cout << "    Period = " << mod.computePeriod_years() << " years\n\n";
    
    // ===== Step 8-12: Test Different Period Configurations =====
    std::cout << "Steps 8-12: Test Multiple Period Configurations\n";
    
    // Solar cycle (11 years)
    mod.optimizeForMetric("solar_cycle");
    double period_solar = mod.computePeriod_years();
    double amp_solar = mod.variables["amplitude"];
    std::cout << "  Solar Cycle: " << std::fixed << std::setprecision(1) << period_solar 
              << " years, amplitude=" << std::setprecision(3) << amp_solar << "\n";
    
    // Long period (20 years)
    mod.optimizeForMetric("long_period");
    double period_long = mod.computePeriod_years();
    double amp_long = mod.variables["amplitude"];
    std::cout << "  Long Period: " << period_long << " years, amplitude=" 
              << std::setprecision(3) << amp_long << "\n";
    
    // Short period (5 years)
    mod.optimizeForMetric("short_period");
    double period_short = mod.computePeriod_years();
    double amp_short = mod.variables["amplitude"];
    std::cout << "  Short Period: " << period_short << " years, amplitude=" 
              << std::setprecision(3) << amp_short << "\n";
    
    // High defect (5%)
    mod.optimizeForMetric("high_defect");
    double amp_high = mod.variables["amplitude"];
    std::cout << "  High Defect: amplitude=" << amp_high << " (5% variation)\n";
    
    // Stable (minimal defect)
    mod.optimizeForMetric("stable");
    double amp_stable = mod.variables["amplitude"];
    std::cout << "  Stable: amplitude=" << amp_stable << " (0.1% variation)\n\n";
    
    // ===== Step 13: Auto-Refinement =====
    std::cout << "Step 13: Auto-Refinement\n";
    std::cout << std::scientific << std::setprecision(3);
    mod.updateVariable("amplitude", 0.1);  // Set beyond limit
    std::cout << "  Set amplitude = 0.1 (beyond 5% limit)\n";
    mod.autoRefineParameters();
    std::cout << "  After refinement: amplitude = " << mod.variables["amplitude"] 
              << " (clamped to 0.05)\n\n";
    
    // ===== Step 14: Calibration (Observational Data) =====
    std::cout << "Step 14: Calibration to Observations\n";
    std::map<std::string, double> obs_data;
    obs_data["amplitude"] = 0.012;  // Observed 1.2% defect
    obs_data["freq"] = 0.00086;  // ~20 year period
    mod.calibrateToObservations(obs_data);
    std::cout << "  Calibrated: amplitude = " << mod.variables["amplitude"] << "\n";
    std::cout << "  Calibrated: period = " << mod.computePeriod_years() << " years\n\n";
    
    // ===== Step 15: Parameter Variations =====
    std::cout << "Step 15: Generate Parameter Variations\n";
    std::vector<std::map<std::string, double>> variations = mod.generateVariations(5);
    std::cout << "  Generated " << variations.size() << " parameter variations:\n";
    for (size_t i = 0; i < std::min(size_t(3), variations.size()); ++i) {
        double period_var = variations[i]["period_days"] / 365.25;
        std::cout << "    Variant " << (i+1) << ": amplitude=" << variations[i]["amplitude"] 
                  << ", period=" << std::fixed << std::setprecision(1) << period_var << " years\n";
    }
    std::cout << std::scientific << std::setprecision(3) << "\n";
    
    // ===== Step 16: Mutation =====
    std::cout << "Step 16: Mutate Parameters\n";
    mod.updateVariable("amplitude", 0.01);  // Reset to baseline
    mod.updateVariable("freq", 0.001);
    mod.mutateParameters(0.15);  // 15% mutation rate
    std::cout << "  After mutation: amplitude = " << mod.variables["amplitude"] << "\n";
    std::cout << "  After mutation: period = " << mod.computePeriod_years() << " years\n\n";
    
    // ===== Step 17: System Evolution =====
    std::cout << "Step 17: Evolve System (Optimize Defect Parameters)\n";
    mod.evolveSystem(10);  // 10 generations
    std::cout << "  After evolution: amplitude = " << mod.variables["amplitude"] << "\n";
    std::cout << "  After evolution: period = " << mod.computePeriod_years() << " years\n\n";
    
    // ===== Step 18: State Management =====
    std::cout << "Step 18: State Management\n";
    mod.updateVariable("amplitude", 0.01);
    mod.updateVariable("freq", 0.001);
    mod.saveState("baseline_11yr");
    std::cout << "  Saved state 'baseline_11yr'\n";
    
    mod.updateVariable("amplitude", 0.02);
    mod.updateVariable("freq", 0.0005);
    mod.saveState("high_amplitude_20yr");
    std::cout << "  Saved state 'high_amplitude_20yr'\n";
    
    std::vector<std::string> saved = mod.listSavedStates();
    std::cout << "  Total saved states: " << saved.size() << "\n";
    
    mod.restoreState("baseline_11yr");
    std::cout << "  Restored 'baseline_11yr': amplitude = " << mod.variables["amplitude"] << "\n\n";
    
    // ===== Step 19: Export State =====
    std::cout << "Step 19: Export State\n";
    std::string exported = mod.exportState();
    std::cout << "  Exported " << exported.length() << " bytes of state data\n";
    std::cout << "  (Can be saved to file for archival/restoration)\n\n";
    
    // ===== Step 20: Sensitivity Analysis =====
    std::cout << "Step 20: Sensitivity Analysis (U_g1 response)\n";
    std::map<std::string, double> sensitivity = mod.sensitivityAnalysis("U_g1");
    std::cout << "  Sensitivity of U_g1 to parameter changes:\n";
    for (const auto& pair : sensitivity) {
        std::cout << "    " << pair.first << ": " << std::fixed << std::setprecision(2) 
                  << pair.second << "\n";
    }
    std::cout << std::scientific << std::setprecision(3) << "\n";
    
    // ===== Step 21: Validation =====
    std::cout << "Step 21: Consistency Validation\n";
    bool valid = mod.validateConsistency();
    std::cout << "  System consistency: " << (valid ? "VALID" : "INVALID") << "\n";
    if (!valid) {
        mod.autoCorrectAnomalies();
        std::cout << "  Auto-corrected anomalies\n";
    }
    std::cout << "\n";
    
    // ===== Step 22: Generate Full Report =====
    std::cout << "Step 22: Generate Full Report\n";
    std::string report = mod.generateReport();
    std::cout << report << "\n";
    
    // ===== Step 23-26: Defect Oscillation Over Time =====
    std::cout << "Steps 23-26: Defect Oscillation Analysis (One Complete Cycle)\n";
    mod.updateVariable("amplitude", 0.01);
    mod.updateVariable("freq", 0.001);
    double period = mod.variables["period_days"];
    
    std::cout << "  Time (years) | δ_def    | U_g1 (J/m³)  | Variation\n";
    std::cout << "  -----------------------------------------------------------\n";
    
    struct TimePoint {
        double fraction;
        std::string label;
    };
    
    std::vector<TimePoint> time_points = {
        {0.00, "Start"},
        {0.25, "Quarter (peak)"},
        {0.50, "Half"},
        {0.75, "Three-quarter"},
        {1.00, "Full cycle"}
    };
    
    for (const auto& tp : time_points) {
        double t_days = period * tp.fraction;
        double t_years = t_days / 365.25;
        double delta = mod.computeDelta_def(t_days);
        double ug1 = mod.computeU_g1(t_days, r_earth);
        double variation = ((ug1 / ug1_t0) - 1.0) * 100.0;
        
        std::cout << "  " << std::fixed << std::setprecision(2) << std::setw(12) << std::left << t_years
                  << " | " << std::setprecision(4) << std::setw(8) << delta
                  << " | " << std::scientific << std::setprecision(3) << std::setw(12) << ug1
                  << " | " << std::fixed << std::setprecision(2) << std::showpos << variation << "%\n"
                  << std::noshowpos;
    }
    
    std::cout << "\n========== DEMONSTRATION COMPLETE ==========\n";
    std::cout << "Ug1 Defect module validated across complete oscillation cycle.\n";
    std::cout << "Defect mechanism: Small amplitude (~1%) oscillations with long periods (~17 years).\n";
    std::cout << "Physical origin: [SCm] density variations, internal stellar dynamics.\n";
    std::cout << "Applications: Solar cycle modeling, stellar irregularities, gravity perturbations.\n";
    std::cout << "UQFF Integration: Time-dependent U_g1 internal dipole gravity factor.\n";
    
    return 0;
}
*/
// Compile: g++ -o defect_test defect_test.cpp Ug1DefectModule.cpp -lm
// Sample: ?_def=0 at t=0; U_g1?4.56e31 J/m� at peak (+1%); period~17.22 yr.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

Ug1DefectModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeDelta_def, computeU_g1, computePeriod_years) are clear, concise, and variable - driven.
- Automatic recalculation of derived variables(period_days) when frequency changes.
- Output and debugging functions(printVariables, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.
- Models small oscillatory defect in internal dipole gravity, supporting time - dependent perturbation analysis.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in defect factor modeling for universal gravity.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.