// SolarCycleFrequencyModule.h
// Modular C++ implementation of the Solar Cycle Frequency (?_c) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes ?_c = 2? / 3.96e8 s?� (~1.59e-8 rad/s, period ~12.55 years); used in sin(?_c t) for ?_j in U_m.
// Pluggable: #include "SolarCycleFrequencyModule.h"
// SolarCycleFrequencyModule mod; mod.computeMuJExample(0.0); mod.updateVariable("period", new_value);
// Variables in std::map; example for Sun at t=0 (sin=0, ?_j=3.38e23 T�m�); t~1 year: slight increase.
// Approximations: Period=3.96e8 s (~12.55 yr); base B_j=1e3 T.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef SOLAR_CYCLE_FREQUENCY_MODULE_H
#define SOLAR_CYCLE_FREQUENCY_MODULE_H

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

class SolarCycleFrequencyModule {
private:
    std::map<std::string, double> variables;
    double computeOmega_c();
    double computeSinOmegaCT(double t);

public:
    // Constructor: Initialize with framework defaults
    SolarCycleFrequencyModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeOmega_c();  // 2? / period s?�
    double computeSinOmegaCT(double t);  // sin(?_c t)
    double computeMuJExample(double t);  // (10^3 + 0.4 sin(?_c t)) * 3.38e20 T�m�

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // ===== ENHANCED METHODS =====
    // Variable Management
    void createVariable(const std::string& name, double value);
    void removeVariable(const std::string& name);
    void cloneVariable(const std::string& source, const std::string& dest);
    std::vector<std::string> listVariables() const;
    std::string getSystemName() const;
    double getVariable(const std::string& name) const { return variables.at(name); }

    // Batch Operations
    void transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func);
    void scaleVariableGroup(const std::vector<std::string>& names, double factor);

    // Self-Expansion: Domain-Specific Scales
    void expandParameterSpace(double frequency_scale, double magnetic_scale, double time_scale);
    void expandFrequencyScale(double omega_factor, double period_factor);  // ω_c and period
    void expandMagneticScale(double mu_factor, double amplitude_factor);   // μ_j and sin amplitude
    void expandCycleScale(double cycle_factor, double phase_factor);       // Cycle characteristics

    // Self-Refinement
    void autoRefineParameters(const std::string& target, double goal);
    void calibrateToObservations(const std::map<std::string, double>& observations);
    void optimizeForMetric(const std::string& metric);

    // Parameter Exploration
    std::vector<std::map<std::string, double>> generateVariations(int count, double variation_pct);

    // Adaptive Evolution
    void mutateParameters(double mutation_rate);
    void evolveSystem(int generations, std::function<double()> fitness_func);

    // State Management
    void saveState(const std::string& label);
    void restoreState(const std::string& label);
    std::vector<std::string> listSavedStates() const;
    std::string exportState() const;

    // System Analysis
    std::map<std::string, double> sensitivityAnalysis(const std::vector<std::string>& params);
    std::string generateReport() const;
    bool validateConsistency() const;
    void autoCorrectAnomalies();
};

#endif // SOLAR_CYCLE_FREQUENCY_MODULE_H

// SolarCycleFrequencyModule.cpp
#include "SolarCycleFrequencyModule.h"

// Constructor: Set framework defaults
SolarCycleFrequencyModule::SolarCycleFrequencyModule() {
    // Universal constants
    variables["pi"] = 3.141592653589793;
    variables["period"] = 3.96e8;                   // s (~12.55 years)
    variables["base_mu"] = 3.38e20;                 // T�m�
    variables["B_j"] = 1e3;                         // Base T
    variables["t"] = 0.0;                           // s

    // Derived
    variables["omega_c"] = computeOmega_c();
}

// Update variable
void SolarCycleFrequencyModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
        if (name == "period") {
            variables["omega_c"] = computeOmega_c();
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void SolarCycleFrequencyModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        if (name == "period") {
            variables["omega_c"] = computeOmega_c();
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void SolarCycleFrequencyModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute ?_c = 2? / period
double SolarCycleFrequencyModule::computeOmega_c() {
    return 2.0 * variables["pi"] / variables["period"];
}

// Compute sin(?_c t)
double SolarCycleFrequencyModule::computeSinOmegaCT(double t) {
    variables["t"] = t;
    return std::sin(variables["omega_c"] * t);
}

// Example ?_j = (10^3 + 0.4 sin(?_c t)) * 3.38e20
double SolarCycleFrequencyModule::computeMuJExample(double t) {
    double sin_omega = computeSinOmegaCT(t);
    double b_j = variables["B_j"] + 0.4 * sin_omega;
    return b_j * variables["base_mu"];
}

// Equation text
std::string SolarCycleFrequencyModule::getEquationText() {
    return "?_c = 2? / 3.96e8 s?� ?1.59e-8 rad/s (period ~12.55 yr, near 11-yr solar cycle);\n"
           "In U_m: ?_j = (10^3 + 0.4 sin(?_c t)) * 3.38e20 T�m� (cyclic magnetic variation).\n"
           "In U_g3: ... cos(?_s t ?) ... (?_s Sun rotation, but ?_c for cycle).\n"
           "Example t=0: sin=0 ? ?_j=3.38e23 T�m�;\n"
           "t=3.14e7 s (~1 yr): sin?0.477 ? ?_j?3.381e23 T�m� (+0.019%).\n"
           "Role: Models solar cycle periodicity; magnetic activity in strings/fields.\n"
           "UQFF: Cyclic effects in jets/nebulae/formation; near 11-yr Hale cycle.";
}

// Print variables
void SolarCycleFrequencyModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// ===== ENHANCED METHODS IMPLEMENTATION =====

namespace solar_cycle_frequency_saved_states {
    static std::map<std::string, std::map<std::string, double>> saved_states;
}

// Variable Management
void SolarCycleFrequencyModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void SolarCycleFrequencyModule::removeVariable(const std::string& name) {
    variables.erase(name);
}

void SolarCycleFrequencyModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
    }
}

std::vector<std::string> SolarCycleFrequencyModule::listVariables() const {
    std::vector<std::string> names;
    for (const auto& pair : variables) {
        names.push_back(pair.first);
    }
    return names;
}

std::string SolarCycleFrequencyModule::getSystemName() const {
    return "Solar_Cycle_Frequency_UQFF";
}

// Batch Operations
void SolarCycleFrequencyModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = func(variables[name]);
        }
    }
    // Update derived
    variables["omega_c"] = computeOmega_c();
}

void SolarCycleFrequencyModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// Self-Expansion: Domain-Specific Scales
void SolarCycleFrequencyModule::expandParameterSpace(double frequency_scale, double magnetic_scale, double time_scale) {
    variables["omega_c"] *= frequency_scale;
    variables["period"] *= (1.0 / frequency_scale);  // Inverse relationship
    variables["base_mu"] *= magnetic_scale;
    variables["B_j"] *= magnetic_scale;
    variables["t"] *= time_scale;
}

void SolarCycleFrequencyModule::expandFrequencyScale(double omega_factor, double period_factor) {
    variables["omega_c"] *= omega_factor;
    variables["period"] *= period_factor;
    // Recalculate to maintain consistency
    variables["omega_c"] = computeOmega_c();
}

void SolarCycleFrequencyModule::expandMagneticScale(double mu_factor, double amplitude_factor) {
    variables["base_mu"] *= mu_factor;
    variables["B_j"] *= mu_factor;
    // Amplitude factor would affect the 0.4 coefficient in sin term
    // Store as variable for future use
    if (variables.find("sin_amplitude") == variables.end()) {
        variables["sin_amplitude"] = 0.4;
    }
    variables["sin_amplitude"] *= amplitude_factor;
}

void SolarCycleFrequencyModule::expandCycleScale(double cycle_factor, double phase_factor) {
    // Cycle factor affects period (longer/shorter cycles)
    variables["period"] *= cycle_factor;
    variables["omega_c"] = computeOmega_c();
    
    // Phase factor could shift the cycle (stored for future use)
    if (variables.find("phase_offset") == variables.end()) {
        variables["phase_offset"] = 0.0;
    }
    variables["phase_offset"] += phase_factor;
}

// Self-Refinement
void SolarCycleFrequencyModule::autoRefineParameters(const std::string& target, double goal) {
    if (target == "omega_c") {
        variables["omega_c"] = goal;
        // Update period to match: period = 2π / ω_c
        variables["period"] = 2.0 * variables["pi"] / goal;
    } else if (target == "period") {
        variables["period"] = goal;
        variables["omega_c"] = computeOmega_c();
    } else if (target == "period_years") {
        // Convert years to seconds: 1 year ≈ 3.156×10⁷ s
        double seconds_per_year = 3.15576e7;
        variables["period"] = goal * seconds_per_year;
        variables["omega_c"] = computeOmega_c();
    } else if (target == "frequency_hz") {
        // Convert Hz to rad/s: ω = 2π f
        variables["omega_c"] = 2.0 * variables["pi"] * goal;
        variables["period"] = 2.0 * variables["pi"] / variables["omega_c"];
    } else if (target == "base_magnetic_moment") {
        variables["base_mu"] = goal;
    } else if (target == "magnetic_field") {
        variables["B_j"] = goal;
    } else if (target == "sin_amplitude") {
        if (variables.find("sin_amplitude") == variables.end()) {
            variables["sin_amplitude"] = 0.4;
        }
        variables["sin_amplitude"] = goal;
    }
}

void SolarCycleFrequencyModule::calibrateToObservations(const std::map<std::string, double>& observations) {
    for (const auto& obs : observations) {
        if (variables.find(obs.first) != variables.end()) {
            variables[obs.first] = obs.second;
            if (obs.first == "period") {
                variables["omega_c"] = computeOmega_c();
            }
        }
    }
}

void SolarCycleFrequencyModule::optimizeForMetric(const std::string& metric) {
    if (metric == "standard_11yr") {
        // Standard 11-year solar cycle: period = 11 × 3.15576e7 s ≈ 3.47×10⁸ s
        variables["period"] = 11.0 * 3.15576e7;
        variables["omega_c"] = computeOmega_c();
    } else if (metric == "hale_22yr") {
        // Full Hale cycle (magnetic polarity reversal): 22 years
        variables["period"] = 22.0 * 3.15576e7;
        variables["omega_c"] = computeOmega_c();
    } else if (metric == "standard_12_55yr") {
        // Framework default: 12.55 years
        variables["period"] = 3.96e8;
        variables["omega_c"] = computeOmega_c();
    } else if (metric == "fast_cycle") {
        // Faster cycle: 7 years
        variables["period"] = 7.0 * 3.15576e7;
        variables["omega_c"] = computeOmega_c();
    } else if (metric == "slow_cycle") {
        // Slower cycle: 20 years
        variables["period"] = 20.0 * 3.15576e7;
        variables["omega_c"] = computeOmega_c();
    } else if (metric == "high_amplitude") {
        // Increase magnetic variation amplitude
        if (variables.find("sin_amplitude") == variables.end()) {
            variables["sin_amplitude"] = 0.4;
        }
        variables["sin_amplitude"] = 1.0;
    } else if (metric == "low_amplitude") {
        // Decrease magnetic variation amplitude
        if (variables.find("sin_amplitude") == variables.end()) {
            variables["sin_amplitude"] = 0.4;
        }
        variables["sin_amplitude"] = 0.1;
    }
}

// Parameter Exploration
std::vector<std::map<std::string, double>> SolarCycleFrequencyModule::generateVariations(int count, double variation_pct) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - variation_pct, 1.0 + variation_pct);
    
    for (int i = 0; i < count; ++i) {
        std::map<std::string, double> variant = variables;
        for (auto& pair : variant) {
            if (pair.first != "pi" && pair.first != "omega_c") {
                pair.second *= dis(gen);
            }
        }
        // Recalculate derived
        variant["omega_c"] = 2.0 * variant["pi"] / variant["period"];
        variations.push_back(variant);
    }
    return variations;
}

// Adaptive Evolution
void SolarCycleFrequencyModule::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - mutation_rate, 1.0 + mutation_rate);
    
    for (auto& pair : variables) {
        if (pair.first != "pi" && pair.first != "omega_c") {
            pair.second *= dis(gen);
        }
    }
    // Recalculate derived
    variables["omega_c"] = computeOmega_c();
}

void SolarCycleFrequencyModule::evolveSystem(int generations, std::function<double()> fitness_func) {
    double best_fitness = fitness_func();
    std::map<std::string, double> best_state = variables;
    
    for (int gen = 0; gen < generations; ++gen) {
        mutateParameters(0.05);
        double fitness = fitness_func();
        if (fitness > best_fitness) {
            best_fitness = fitness;
            best_state = variables;
        } else {
            variables = best_state;  // Revert if worse
        }
    }
    variables = best_state;
}

// State Management
void SolarCycleFrequencyModule::saveState(const std::string& label) {
    solar_cycle_frequency_saved_states::saved_states[label] = variables;
}

void SolarCycleFrequencyModule::restoreState(const std::string& label) {
    if (solar_cycle_frequency_saved_states::saved_states.find(label) != solar_cycle_frequency_saved_states.end()) {
        variables = solar_cycle_frequency_saved_states::saved_states[label];
    }
}

std::vector<std::string> SolarCycleFrequencyModule::listSavedStates() const {
    std::vector<std::string> labels;
    for (const auto& pair : solar_cycle_frequency_saved_states::saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string SolarCycleFrequencyModule::exportState() const {
    std::ostringstream oss;
    oss << "SolarCycleFrequency_State_Export:\n";
    for (const auto& pair : variables) {
        oss << pair.first << "=" << std::scientific << pair.second << "\n";
    }
    return oss.str();
}

// System Analysis
std::map<std::string, double> SolarCycleFrequencyModule::sensitivityAnalysis(const std::vector<std::string>& params) {
    std::map<std::string, double> sensitivities;
    double t_test = 3.15576e7;  // Test at 1 year
    double baseline = computeMuJExample(t_test);
    
    for (const auto& param : params) {
        if (variables.find(param) != variables.end() && param != "pi") {
            double original = variables[param];
            variables[param] = original * 1.01;
            
            if (param == "period") {
                variables["omega_c"] = computeOmega_c();
            }
            
            double perturbed = computeMuJExample(t_test);
            sensitivities[param] = (perturbed - baseline) / baseline;
            
            // Restore
            variables[param] = original;
            if (param == "period") {
                variables["omega_c"] = computeOmega_c();
            }
        }
    }
    return sensitivities;
}

std::string SolarCycleFrequencyModule::generateReport() const {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(3);
    oss << "===== Solar Cycle Frequency Module Report =====\n";
    oss << "System: " << getSystemName() << "\n\n";
    
    oss << "Frequency Parameters:\n";
    oss << "  Period = " << variables.at("period") << " s (";
    double period_years = variables.at("period") / 3.15576e7;
    oss << std::fixed << std::setprecision(2) << period_years << " years)\n";
    oss << std::scientific;
    oss << "  ω_c = " << variables.at("omega_c") << " rad/s\n";
    double frequency_hz = variables.at("omega_c") / (2.0 * variables.at("pi"));
    oss << "  f = " << frequency_hz << " Hz\n\n";
    
    oss << "Magnetic Parameters:\n";
    oss << "  Base μ = " << variables.at("base_mu") << " T²·m²\n";
    oss << "  B_j (base field) = " << variables.at("B_j") << " T\n";
    double sin_amp = (variables.find("sin_amplitude") != variables.end()) ? 
                     variables.at("sin_amplitude") : 0.4;
    oss << "  Sin amplitude = " << std::fixed << std::setprecision(2) << sin_amp << "\n";
    oss << "  Current t = " << std::scientific << variables.at("t") << " s\n\n";
    
    oss << "Cyclic Magnetic Moment Examples (μ_j):\n";
    std::vector<double> test_times = {0.0, 0.25, 0.5, 0.75, 1.0, 2.0, 5.0, 10.0};
    for (double t_years : test_times) {
        double t_seconds = t_years * 3.15576e7;
        double sin_val = std::sin(variables.at("omega_c") * t_seconds);
        double b_eff = variables.at("B_j") + sin_amp * sin_val;
        double mu_j = b_eff * variables.at("base_mu");
        oss << "  t=" << std::fixed << std::setprecision(2) << t_years << " yr: ";
        oss << "sin(ω_c t)=" << std::fixed << std::setprecision(3) << sin_val << ", ";
        oss << "μ_j=" << std::scientific << mu_j << " T²·m²\n";
    }
    oss << "\n";
    
    // Cycle characteristics
    oss << "Solar Cycle Characteristics:\n";
    if (period_years >= 10.0 && period_years <= 12.0) {
        oss << "  Standard solar cycle range (10-12 years)\n";
    } else if (period_years >= 20.0 && period_years <= 24.0) {
        oss << "  Hale cycle range (full magnetic polarity reversal, ~22 years)\n";
    } else if (period_years < 10.0) {
        oss << "  Fast cycle (<10 years)\n";
    } else {
        oss << "  Extended cycle (>12 years)\n";
    }
    
    oss << "\nPhysical Interpretation:\n";
    oss << "  Cycle period: " << std::fixed << std::setprecision(1) << period_years << " years\n";
    oss << "  Angular frequency: " << std::scientific << variables.at("omega_c") << " rad/s\n";
    oss << "  Magnetic variation: ±" << std::fixed << std::setprecision(1) 
        << (sin_amp * 100.0 / variables.at("B_j")) << "% around baseline\n";
    
    oss << "\n  Applications:\n";
    oss << "    - Solar magnetic activity: 11-year sunspot cycle modulation\n";
    oss << "    - Heliospheric variations: Cyclic magnetic field reversals\n";
    oss << "    - Planetary influences: Solar wind pressure variations\n";
    oss << "    - Stellar jets: Periodic magnetic string activity\n";
    oss << "    - Star formation: Cyclic accretion disk dynamics\n";
    oss << "    - Exoplanet habitability: Long-term stellar activity patterns\n";
    
    return oss.str();
}

bool SolarCycleFrequencyModule::validateConsistency() const {
    bool valid = true;
    
    // Check period is positive
    if (variables.find("period") != variables.end() && variables.at("period") <= 0) {
        std::cerr << "Error: period <= 0 (cycle period must be positive)\n";
        valid = false;
    }
    
    // Check ω_c consistency
    if (variables.find("omega_c") != variables.end()) {
        double expected = 2.0 * variables.at("pi") / variables.at("period");
        double actual = variables.at("omega_c");
        if (std::abs(expected - actual) / expected > 1e-9) {
            std::cerr << "Error: omega_c inconsistent (expected " << expected << ", got " << actual << ")\n";
            valid = false;
        }
    }
    
    // Check period is in reasonable range (1-100 years)
    double period_years = variables.at("period") / 3.15576e7;
    if (period_years < 1.0 || period_years > 100.0) {
        std::cerr << "Warning: Period outside typical range [1, 100] years (current: " 
                  << period_years << " years)\n";
    }
    
    // Check magnetic parameters are positive
    if (variables.find("base_mu") != variables.end() && variables.at("base_mu") <= 0) {
        std::cerr << "Error: base_mu <= 0 (magnetic moment must be positive)\n";
        valid = false;
    }
    
    if (variables.find("B_j") != variables.end() && variables.at("B_j") <= 0) {
        std::cerr << "Error: B_j <= 0 (magnetic field must be positive)\n";
        valid = false;
    }
    
    return valid;
}

void SolarCycleFrequencyModule::autoCorrectAnomalies() {
    // Reset period to typical value if out of range
    double period_years = variables["period"] / 3.15576e7;
    if (variables["period"] <= 0 || period_years > 100.0 || period_years < 1.0) {
        variables["period"] = 3.96e8;  // Standard 12.55 years
        variables["omega_c"] = computeOmega_c();
    }
    
    // Ensure magnetic parameters are positive and reasonable
    if (variables["base_mu"] <= 0 || variables["base_mu"] > 1e30) {
        variables["base_mu"] = 3.38e20;  // Standard value
    }
    
    if (variables["B_j"] <= 0 || variables["B_j"] > 1e6) {
        variables["B_j"] = 1e3;  // Standard value
    }
    
    // Ensure pi is correct
    if (std::abs(variables["pi"] - 3.141592653589793) > 1e-12) {
        variables["pi"] = 3.141592653589793;
        variables["omega_c"] = computeOmega_c();
    }
}

// Example usage in base program (snippet)
int main() {
    SolarCycleFrequencyModule module;
    std::cout << "===== Solar Cycle Frequency Module Enhanced Demonstration =====\n\n";
    
    // Step 1: Report initial state
    std::cout << "STEP 1: Initial Configuration (Period = 3.96e8 s ≈ 12.55 years)\n";
    std::cout << module.generateReport() << "\n";
    
    // Step 2: Compute μ_j at various times
    std::cout << "STEP 2: Magnetic Moment Evolution Over Multiple Cycles\n";
    std::vector<double> test_years = {0.0, 2.5, 5.0, 7.5, 10.0, 12.55, 25.1, 50.0};
    for (double t_yr : test_years) {
        double t_sec = t_yr * 3.15576e7;
        double mu_j = module.computeMuJExample(t_sec);
        double sin_val = module.computeSinOmegaCT(t_sec);
        std::cout << "  t=" << std::fixed << std::setprecision(2) << t_yr << " yr: ";
        std::cout << "sin(ω_c t)=" << std::fixed << std::setprecision(3) << sin_val << ", ";
        std::cout << "μ_j=" << std::scientific << mu_j << " T²·m²\n";
    }
    std::cout << "\n";
    
    // Step 3: Save initial state
    std::cout << "STEP 3: Save Initial State\n";
    module.saveState("standard_12_55yr");
    std::cout << "State saved as 'standard_12_55yr'\n\n";
    
    // Step 4: Test standard 11-year solar cycle
    std::cout << "STEP 4: Test Standard 11-Year Solar Cycle\n";
    module.optimizeForMetric("standard_11yr");
    double omega_11 = module.getVariable("omega_c");
    double period_11 = module.getVariable("period");
    double mu_11 = module.computeMuJExample(3.15576e7);  // At 1 year
    std::cout << "11-year cycle:\n";
    std::cout << "  Period = " << std::scientific << period_11 << " s (";
    std::cout << std::fixed << std::setprecision(2) << (period_11 / 3.15576e7) << " years)\n";
    std::cout << "  ω_c = " << std::scientific << omega_11 << " rad/s\n";
    std::cout << "  μ_j at t=1 yr = " << mu_11 << " T²·m²\n\n";
    
    // Step 5: Test Hale 22-year cycle
    std::cout << "STEP 5: Test Hale 22-Year Magnetic Cycle\n";
    module.restoreState("standard_12_55yr");
    module.optimizeForMetric("hale_22yr");
    double omega_22 = module.getVariable("omega_c");
    double period_22 = module.getVariable("period");
    double mu_22 = module.computeMuJExample(3.15576e7);
    std::cout << "22-year Hale cycle:\n";
    std::cout << "  Period = " << std::scientific << period_22 << " s (";
    std::cout << std::fixed << std::setprecision(2) << (period_22 / 3.15576e7) << " years)\n";
    std::cout << "  ω_c = " << std::scientific << omega_22 << " rad/s\n";
    std::cout << "  μ_j at t=1 yr = " << mu_22 << " T²·m²\n\n";
    
    // Step 6: Test fast 7-year cycle
    std::cout << "STEP 6: Test Fast 7-Year Cycle\n";
    module.restoreState("standard_12_55yr");
    module.optimizeForMetric("fast_cycle");
    double omega_fast = module.getVariable("omega_c");
    double period_fast = module.getVariable("period");
    std::cout << "Fast 7-year cycle:\n";
    std::cout << "  Period = " << std::scientific << period_fast << " s (";
    std::cout << std::fixed << std::setprecision(2) << (period_fast / 3.15576e7) << " years)\n";
    std::cout << "  ω_c = " << std::scientific << omega_fast << " rad/s\n\n";
    
    // Step 7: Test slow 20-year cycle
    std::cout << "STEP 7: Test Slow 20-Year Cycle\n";
    module.restoreState("standard_12_55yr");
    module.optimizeForMetric("slow_cycle");
    double omega_slow = module.getVariable("omega_c");
    double period_slow = module.getVariable("period");
    std::cout << "Slow 20-year cycle:\n";
    std::cout << "  Period = " << std::scientific << period_slow << " s (";
    std::cout << std::fixed << std::setprecision(2) << (period_slow / 3.15576e7) << " years)\n";
    std::cout << "  ω_c = " << std::scientific << omega_slow << " rad/s\n\n";
    
    // Step 8: Restore and expand frequency scale
    std::cout << "STEP 8: Expand Frequency Scale (ω_c x2, faster oscillation)\n";
    module.restoreState("standard_12_55yr");
    module.expandFrequencyScale(2.0, 0.5);  // Double frequency, halve period
    double omega_expanded = module.getVariable("omega_c");
    double period_expanded = module.getVariable("period");
    std::cout << "Expanded frequency:\n";
    std::cout << "  ω_c = " << std::scientific << omega_expanded << " rad/s (2x)\n";
    std::cout << "  Period = " << period_expanded << " s (";
    std::cout << std::fixed << std::setprecision(2) << (period_expanded / 3.15576e7) << " years)\n\n";
    
    // Step 9: Restore and expand magnetic scale
    std::cout << "STEP 9: Expand Magnetic Scale (μ_base x2, B_j x2)\n";
    module.restoreState("standard_12_55yr");
    module.expandMagneticScale(2.0, 1.0);
    double mu_base_new = module.getVariable("base_mu");
    double b_j_new = module.getVariable("B_j");
    double mu_scaled = module.computeMuJExample(0.0);
    std::cout << "Expanded magnetics:\n";
    std::cout << "  Base μ = " << std::scientific << mu_base_new << " T²·m² (2x)\n";
    std::cout << "  B_j = " << b_j_new << " T (2x)\n";
    std::cout << "  μ_j at t=0 = " << mu_scaled << " T²·m²\n\n";
    
    // Step 10: Restore and expand cycle scale
    std::cout << "STEP 10: Expand Cycle Scale (Period x1.5)\n";
    module.restoreState("standard_12_55yr");
    module.expandCycleScale(1.5, 0.0);
    double period_cycle = module.getVariable("period");
    double omega_cycle = module.getVariable("omega_c");
    std::cout << "Expanded cycle:\n";
    std::cout << "  Period = " << std::scientific << period_cycle << " s (";
    std::cout << std::fixed << std::setprecision(2) << (period_cycle / 3.15576e7) << " years)\n";
    std::cout << "  ω_c = " << std::scientific << omega_cycle << " rad/s\n\n";
    
    // Step 11: Sensitivity analysis
    std::cout << "STEP 11: Sensitivity Analysis (at t=1 year)\n";
    module.restoreState("standard_12_55yr");
    std::vector<std::string> params = {"period", "base_mu", "B_j", "t"};
    auto sensitivities = module.sensitivityAnalysis(params);
    for (const auto& pair : sensitivities) {
        std::cout << "  ∂μ_j/∂" << pair.first << " ≈ " << std::scientific << pair.second << " (normalized)\n";
    }
    std::cout << "\n";
    
    // Step 12: Generate variations
    std::cout << "STEP 12: Generate Parameter Variations (5 variants, ±10%)\n";
    auto variations = module.generateVariations(5, 0.1);
    for (int i = 0; i < variations.size(); ++i) {
        double var_period = variations[i]["period"];
        double var_period_yr = var_period / 3.15576e7;
        double var_omega = variations[i]["omega_c"];
        std::cout << "  Variant " << (i+1) << ": Period=" << std::fixed << std::setprecision(2) 
                  << var_period_yr << " yr, ω_c=" << std::scientific << var_omega << " rad/s\n";
    }
    std::cout << "\n";
    
    // Step 13: Auto-refine to target period
    std::cout << "STEP 13: Auto-Refine to Target Period (11.0 years)\n";
    module.restoreState("standard_12_55yr");
    module.autoRefineParameters("period_years", 11.0);
    double refined_period = module.getVariable("period");
    double refined_omega = module.getVariable("omega_c");
    double refined_period_yr = refined_period / 3.15576e7;
    std::cout << "Refined parameters:\n";
    std::cout << "  Period = " << std::fixed << std::setprecision(2) << refined_period_yr 
              << " years (target: 11.0)\n";
    std::cout << "  ω_c = " << std::scientific << refined_omega << " rad/s\n\n";
    
    // Step 14: Auto-refine to target frequency
    std::cout << "STEP 14: Auto-Refine to Target Angular Frequency (2e-8 rad/s)\n";
    module.restoreState("standard_12_55yr");
    module.autoRefineParameters("omega_c", 2.0e-8);
    double omega_target = module.getVariable("omega_c");
    double period_target = module.getVariable("period");
    std::cout << "Refined ω_c = " << std::scientific << omega_target << " rad/s (target: 2e-8)\n";
    std::cout << "Resulting period = " << std::fixed << std::setprecision(2) 
              << (period_target / 3.15576e7) << " years\n\n";
    
    // Step 15: Auto-refine magnetic amplitude
    std::cout << "STEP 15: Auto-Refine Sin Amplitude (1.0 for larger variation)\n";
    module.restoreState("standard_12_55yr");
    module.autoRefineParameters("sin_amplitude", 1.0);
    double sin_amp = module.getVariable("sin_amplitude");
    std::cout << "Refined sin amplitude = " << std::fixed << std::setprecision(2) 
              << sin_amp << " (was 0.4)\n";
    std::cout << "Variation range: ±" << std::fixed << std::setprecision(1) 
              << (sin_amp * 100.0 / module.getVariable("B_j")) << "% of baseline B_j\n\n";
    
    // Step 16: Test high amplitude cycle
    std::cout << "STEP 16: Optimize for High Amplitude Metric\n";
    module.restoreState("standard_12_55yr");
    module.optimizeForMetric("high_amplitude");
    double amp_high = module.getVariable("sin_amplitude");
    double mu_high_max = module.computeMuJExample(module.getVariable("period") / 4.0);
    double mu_high_min = module.computeMuJExample(3.0 * module.getVariable("period") / 4.0);
    std::cout << "High amplitude (sin_amp=" << std::fixed << std::setprecision(1) << amp_high << "):\n";
    std::cout << "  μ_j max ≈ " << std::scientific << mu_high_max << " T²·m²\n";
    std::cout << "  μ_j min ≈ " << mu_high_min << " T²·m²\n";
    std::cout << "  Variation: " << std::fixed << std::setprecision(1) 
              << ((mu_high_max - mu_high_min) / mu_high_max * 100.0) << "%\n\n";
    
    // Step 17: Calibrate to observations
    std::cout << "STEP 17: Calibrate to Observational Data\n";
    module.restoreState("standard_12_55yr");
    std::map<std::string, double> observations;
    observations["period"] = 11.0 * 3.15576e7;  // Observed 11-year cycle
    observations["B_j"] = 1.2e3;                 // Observed higher baseline field
    module.calibrateToObservations(observations);
    std::cout << "Calibrated parameters:\n";
    std::cout << "  Period = " << std::fixed << std::setprecision(2) 
              << (module.getVariable("period") / 3.15576e7) << " years\n";
    std::cout << "  B_j = " << std::scientific << module.getVariable("B_j") << " T\n";
    std::cout << "  ω_c = " << module.getVariable("omega_c") << " rad/s\n\n";
    
    // Step 18: Mutate parameters
    std::cout << "STEP 18: Mutate Parameters (5% random variation)\n";
    module.mutateParameters(0.05);
    double mutated_period = module.getVariable("period") / 3.15576e7;
    double mutated_omega = module.getVariable("omega_c");
    std::cout << "Mutated parameters:\n";
    std::cout << "  Period = " << std::fixed << std::setprecision(2) << mutated_period << " years\n";
    std::cout << "  ω_c = " << std::scientific << mutated_omega << " rad/s\n\n";
    
    // Step 19: Validate consistency
    std::cout << "STEP 19: Validate Consistency\n";
    bool valid = module.validateConsistency();
    std::cout << "Consistency check: " << (valid ? "PASSED" : "FAILED") << "\n\n";
    
    // Step 20: Introduce anomaly and auto-correct
    std::cout << "STEP 20: Introduce Anomaly and Auto-Correct\n";
    module.createVariable("period_anomaly", -100.0);
    module.removeVariable("period");
    module.createVariable("period", -100.0);  // Invalid negative period
    std::cout << "Introduced invalid period = " << module.getVariable("period") << " s\n";
    module.autoCorrectAnomalies();
    std::cout << "Auto-corrected period = " << std::scientific << module.getVariable("period") 
              << " s (";
    std::cout << std::fixed << std::setprecision(2) << (module.getVariable("period") / 3.15576e7) 
              << " years)\n";
    bool valid_after = module.validateConsistency();
    std::cout << "Consistency after correction: " << (valid_after ? "PASSED" : "FAILED") << "\n\n";
    
    // Step 21: List saved states and export
    std::cout << "STEP 21: List Saved States and Export Final State\n";
    auto states = module.listSavedStates();
    std::cout << "Saved states:\n";
    for (const auto& state : states) {
        std::cout << "  - " << state << "\n";
    }
    std::cout << "\n";
    std::string exported = module.exportState();
    std::cout << exported << "\n";
    
    std::cout << "===== Demonstration Complete =====\n";
    return 0;
}

// Original commented example preserved below for reference:
// #include "SolarCycleFrequencyModule.h"
// int main() {
//     SolarCycleFrequencyModule mod;
//     double omega = mod.computeOmega_c();
//     std::cout << "?_c ? " << omega << " rad/s\n";
//     double mu = mod.computeMuJExample(0.0);
//     std::cout << "?_j (t=0) = " << mu << " T�m�\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("period", 3.8e8);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o solar_cycle_test solar_cycle_test.cpp SolarCycleFrequencyModule.cpp -lm
// Sample: ?_c?1.59e-8 rad/s; ?_j?3.38e23 T�m�; periodic variation.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

SolarCycleFrequencyModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeOmega_c, computeSinOmegaCT, computeMuJExample) are clear, concise, and variable - driven.
- Automatic recalculation of derived variables(omega_c) when dependencies change(e.g., period).
- Output and debugging functions(printVariables, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.
- Models solar cycle periodicity and magnetic activity with realistic parameters.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in solar cycle frequency modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.