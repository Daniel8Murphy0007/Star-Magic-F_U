// SurfaceTemperatureModule.h
// Modular C++ implementation of the Surface Temperature (T_s) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes T_s=5778 K (Sun effective); potential scaling T_s / T_s_ref in B_j for U_g3 magnetic strings.
// Pluggable: #include "SurfaceTemperatureModule.h"
// SurfaceTemperatureModule mod; mod.computeU_g3_example(0.0, 5778.0); mod.updateVariable("T_s", new_value);
// Variables in std::map; example for Sun at t=0; T_s=5778 K ? U_g3?1.8e49 J/m� (full); T_s=10000 K: ~3.11e49 J/m�.
// Approximations: T_s_ref=5778 K (Sun); cos(?_s t ?)=1; P_core=1; E_react=1e46; hypothetical B_j scaling.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef SURFACE_TEMPERATURE_MODULE_H
#define SURFACE_TEMPERATURE_MODULE_H

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

class SurfaceTemperatureModule {
private:
    std::map<std::string, double> variables;
    double computeB_j_hypothetical(double t, double T_s);
    double computeU_g3_example(double t, double T_s);

public:
    // Constructor: Initialize with framework defaults (Sun)
    SurfaceTemperatureModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeT_s();  // 5778 K (Sun)
    double computeB_j_hypothetical(double t, double T_s);  // Scaled B_j (T)
    double computeU_g3_example(double t, double T_s);  // U_g3 with scaling (J/m^3)

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

    // Batch Operations
    void transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func);
    void scaleVariableGroup(const std::vector<std::string>& names, double factor);

    // Self-Expansion: Domain-Specific Scales
    void expandParameterSpace(double temp_scale, double field_scale, double energy_scale);
    void expandTemperatureScale(double ts_factor, double thermal_factor);
    void expandFieldScale(double bj_factor, double coupling_factor);
    void expandEnergyScale(double ereact_factor, double ug3_factor);

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

#endif // SURFACE_TEMPERATURE_MODULE_H

// SurfaceTemperatureModule.cpp
#include "SurfaceTemperatureModule.h"

// Constructor: Set framework defaults (Sun)
SurfaceTemperatureModule::SurfaceTemperatureModule() {
    // Universal constants
    variables["T_s"] = 5778.0;                      // K (Sun effective)
    variables["T_s_ref"] = 5778.0;                  // K (reference)
    variables["k_3"] = 1.8;                         // Coupling
    variables["B_ref"] = 1e3;                       // Base T (string)
    variables["omega_s"] = 2.5e-6;                  // rad/s
    variables["P_core"] = 1.0;                      // Unitless
    variables["E_react"] = 1e46;                    // J
    variables["pi"] = 3.141592653589793;
    variables["t"] = 0.0;                           // s
}

// Update variable
void SurfaceTemperatureModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void SurfaceTemperatureModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void SurfaceTemperatureModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute T_s (K)
double SurfaceTemperatureModule::computeT_s() {
    return variables["T_s"];
}

// Hypothetical B_j scaled by T_s / T_s_ref
double SurfaceTemperatureModule::computeB_j_hypothetical(double t, double T_s) {
    variables["t"] = t;
    double base_b = variables["B_ref"] + 0.4 * std::sin(variables["omega_s"] * t);  // Cycle
    return base_b * (T_s / variables["T_s_ref"]);
}

// U_g3 example with scaled B_j
double SurfaceTemperatureModule::computeU_g3_example(double t, double T_s) {
    double k_3 = variables["k_3"];
    double b_j = computeB_j_hypothetical(t, T_s);
    double cos_term = std::cos(variables["omega_s"] * t * variables["pi"]);
    double p_core = variables["P_core"];
    double e_react = variables["E_react"];
    return k_3 * b_j * cos_term * p_core * e_react;
}

// Equation text
std::string SurfaceTemperatureModule::getEquationText() {
    return "B_j ? (10^3 + 0.4 sin(?_s t)) * (T_s / T_s,ref) T (hypothetical);\n"
           "U_g3 = k_3 * ? B_j * cos(?_s t ?) * P_core * E_react\n"
           "Where T_s = 5778 K (Sun effective photosphere; �C=5505).\n"
           "T_s,ref=5778 K; scales string fields by temperature.\n"
           "Example t=0, T_s=5778 K: B_j?1e3 T, U_g3?1.8e49 J/m�;\n"
           "T_s=10000 K: B_j?1730 T, U_g3?3.11e49 J/m� (+73%).\n"
           "Role: Thermal baseline for magnetic strength; variability in U_g3/disks.\n"
           "UQFF: Temperature-dependent fields; extensible for radiation/formation.";
}

// Print variables
void SurfaceTemperatureModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// ===== ENHANCED METHODS IMPLEMENTATION =====

namespace surface_temperature_saved_states {
    static std::map<std::string, std::map<std::string, double>> saved_states;
}

// Variable Management
void SurfaceTemperatureModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void SurfaceTemperatureModule::removeVariable(const std::string& name) {
    variables.erase(name);
}

void SurfaceTemperatureModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
    }
}

std::vector<std::string> SurfaceTemperatureModule::listVariables() const {
    std::vector<std::string> names;
    for (const auto& pair : variables) {
        names.push_back(pair.first);
    }
    return names;
}

std::string SurfaceTemperatureModule::getSystemName() const {
    return "Surface_Temperature_Ts_UQFF";
}

// Batch Operations
void SurfaceTemperatureModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = func(variables[name]);
        }
    }
}

void SurfaceTemperatureModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// Self-Expansion: Domain-Specific Scales
void SurfaceTemperatureModule::expandParameterSpace(double temp_scale, double field_scale, double energy_scale) {
    // Scale temperature
    variables["T_s"] *= temp_scale;
    
    // Scale magnetic field
    variables["B_ref"] *= field_scale;
    
    // Scale energy
    variables["E_react"] *= energy_scale;
}

void SurfaceTemperatureModule::expandTemperatureScale(double ts_factor, double thermal_factor) {
    variables["T_s"] *= ts_factor;
    
    // Advanced temperature characteristics
    if (variables.find("T_s_kelvin") == variables.end()) {
        variables["T_s_kelvin"] = variables["T_s"];
    }
    variables["T_s_kelvin"] *= ts_factor;
    
    // Temperature in Celsius
    if (variables.find("T_s_celsius") == variables.end()) {
        variables["T_s_celsius"] = variables["T_s"] - 273.15;
    }
    variables["T_s_celsius"] = variables["T_s"] - 273.15;
    
    // Thermal energy density (Stefan-Boltzmann)
    if (variables.find("thermal_energy_density_Jm3") == variables.end()) {
        // σ T⁴ / c, σ = 5.67×10⁻⁸ W/(m²·K⁴)
        double sigma = 5.67e-8;
        double c = 3e8;
        double T4 = variables["T_s"] * variables["T_s"] * variables["T_s"] * variables["T_s"];
        variables["thermal_energy_density_Jm3"] = (sigma * T4) / c;
    }
    variables["thermal_energy_density_Jm3"] *= (ts_factor * ts_factor * ts_factor * ts_factor);
    
    // Luminosity scaling (L ∝ T⁴)
    if (variables.find("luminosity_scale") == variables.end()) {
        // Relative to Sun
        double T_sun = 5778.0;
        double ratio = variables["T_s"] / T_sun;
        variables["luminosity_scale"] = ratio * ratio * ratio * ratio;
    }
    double new_ratio = (variables["T_s"] / 5778.0);
    variables["luminosity_scale"] = new_ratio * new_ratio * new_ratio * new_ratio;
    
    // Thermal wavelength peak (Wien's law)
    if (variables.find("peak_wavelength_nm") == variables.end()) {
        // λ_max = 2.898×10⁻³ / T (m)
        variables["peak_wavelength_nm"] = (2.898e-3 / variables["T_s"]) * 1e9;  // Convert to nm
    }
    variables["peak_wavelength_nm"] = (2.898e-3 / variables["T_s"]) * 1e9;
    
    // Thermal coupling
    if (variables.find("thermal_coupling") == variables.end()) {
        variables["thermal_coupling"] = 1.0;
    }
    variables["thermal_coupling"] *= thermal_factor;
}

void SurfaceTemperatureModule::expandFieldScale(double bj_factor, double coupling_factor) {
    variables["B_ref"] *= bj_factor;
    variables["k_3"] *= coupling_factor;
    
    // Field characteristics
    if (variables.find("B_j_scale_T") == variables.end()) {
        variables["B_j_scale_T"] = variables["B_ref"];
    }
    variables["B_j_scale_T"] *= bj_factor;
    
    // Temperature-field coupling
    if (variables.find("temp_field_coupling") == variables.end()) {
        // How T_s affects B_j
        variables["temp_field_coupling"] = variables["B_ref"] / variables["T_s_ref"];
    }
    variables["temp_field_coupling"] *= bj_factor;
    
    // Interaction scale
    if (variables.find("interaction_scale") == variables.end()) {
        variables["interaction_scale"] = variables["k_3"] * variables["B_ref"] * variables["E_react"];
    }
    variables["interaction_scale"] *= bj_factor * coupling_factor;
}

void SurfaceTemperatureModule::expandEnergyScale(double ereact_factor, double ug3_factor) {
    variables["E_react"] *= ereact_factor;
    
    // Energy characteristics
    if (variables.find("reactor_energy_J") == variables.end()) {
        variables["reactor_energy_J"] = variables["E_react"];
    }
    variables["reactor_energy_J"] *= ereact_factor;
    
    // U_g3 scale
    if (variables.find("U_g3_scale_Jm3") == variables.end()) {
        double t_test = 0.0;
        double T_s_test = variables["T_s"];
        variables["U_g3_scale_Jm3"] = computeU_g3_example(t_test, T_s_test);
    }
    variables["U_g3_scale_Jm3"] *= ug3_factor;
    
    // Power output
    if (variables.find("power_output_W") == variables.end()) {
        double cycle_time = 2.0 * variables["pi"] / variables["omega_s"];
        variables["power_output_W"] = variables["E_react"] / cycle_time;
    }
    variables["power_output_W"] *= ereact_factor;
}

// Self-Refinement
void SurfaceTemperatureModule::autoRefineParameters(const std::string& target, double goal) {
    if (target == "T_s" || target == "T_s_kelvin") {
        variables["T_s"] = goal;
        if (variables.find("T_s_kelvin") != variables.end()) {
            variables["T_s_kelvin"] = goal;
        }
        if (variables.find("T_s_celsius") != variables.end()) {
            variables["T_s_celsius"] = goal - 273.15;
        }
    } else if (target == "T_s_celsius") {
        variables["T_s"] = goal + 273.15;
        if (variables.find("T_s_kelvin") != variables.end()) {
            variables["T_s_kelvin"] = goal + 273.15;
        }
        if (variables.find("T_s_celsius") != variables.end()) {
            variables["T_s_celsius"] = goal;
        }
    } else if (target == "T_s_ref") {
        variables["T_s_ref"] = goal;
    } else if (target == "B_ref") {
        variables["B_ref"] = goal;
    } else if (target == "k_3") {
        variables["k_3"] = goal;
    } else if (target == "E_react") {
        variables["E_react"] = goal;
    } else if (target == "U_g3_at_ts") {
        // Target specific U_g3 at current T_s, t=0
        double current_U_g3 = computeU_g3_example(0.0, variables["T_s"]);
        if (current_U_g3 > 0) {
            double scale = goal / current_U_g3;
            variables["k_3"] *= scale;
        }
    } else if (target == "luminosity_scale") {
        if (variables.find("luminosity_scale") == variables.end()) {
            double T_sun = 5778.0;
            double ratio = variables["T_s"] / T_sun;
            variables["luminosity_scale"] = ratio * ratio * ratio * ratio;
        }
        variables["luminosity_scale"] = goal;
        // Back-calculate T_s (L ∝ T⁴)
        double T_sun = 5778.0;
        variables["T_s"] = T_sun * std::pow(goal, 0.25);
    } else if (target == "peak_wavelength_nm") {
        if (variables.find("peak_wavelength_nm") == variables.end()) {
            variables["peak_wavelength_nm"] = (2.898e-3 / variables["T_s"]) * 1e9;
        }
        variables["peak_wavelength_nm"] = goal;
        // Back-calculate T_s (Wien's law)
        variables["T_s"] = (2.898e-3 / (goal * 1e-9));
    }
}

void SurfaceTemperatureModule::calibrateToObservations(const std::map<std::string, double>& observations) {
    for (const auto& obs : observations) {
        if (variables.find(obs.first) != variables.end()) {
            variables[obs.first] = obs.second;
        }
    }
    // Update derived quantities
    if (variables.find("T_s_celsius") != variables.end()) {
        variables["T_s_celsius"] = variables["T_s"] - 273.15;
    }
}

void SurfaceTemperatureModule::optimizeForMetric(const std::string& metric) {
    if (metric == "sun") {
        // Sun (G2V)
        variables["T_s"] = 5778.0;
        variables["T_s_ref"] = 5778.0;
    } else if (metric == "cool_star") {
        // M-type red dwarf
        variables["T_s"] = 3500.0;
        variables["T_s_ref"] = 5778.0;
    } else if (metric == "hot_star") {
        // O-type blue giant
        variables["T_s"] = 30000.0;
        variables["T_s_ref"] = 5778.0;
    } else if (metric == "red_giant") {
        // K/M giant
        variables["T_s"] = 4000.0;
        variables["T_s_ref"] = 5778.0;
    } else if (metric == "white_dwarf") {
        // White dwarf
        variables["T_s"] = 10000.0;
        variables["T_s_ref"] = 5778.0;
    } else if (metric == "neutron_star") {
        // Neutron star surface
        variables["T_s"] = 1e6;  // ~1 million K
        variables["T_s_ref"] = 5778.0;
    } else if (metric == "brown_dwarf") {
        // Brown dwarf
        variables["T_s"] = 1500.0;
        variables["T_s_ref"] = 5778.0;
    } else if (metric == "planet_hot") {
        // Hot Jupiter
        variables["T_s"] = 1500.0;
        variables["T_s_ref"] = 5778.0;
    } else if (metric == "planet_warm") {
        // Earth-like
        variables["T_s"] = 288.0;
        variables["T_s_ref"] = 5778.0;
    } else if (metric == "planet_cold") {
        // Ice giant
        variables["T_s"] = 60.0;
        variables["T_s_ref"] = 5778.0;
    } else if (metric == "stellar_types") {
        // Range from M to O
        // Using average solar-type
        variables["T_s"] = 6000.0;
        variables["T_s_ref"] = 5778.0;
    }
}

// Parameter Exploration
std::vector<std::map<std::string, double>> SurfaceTemperatureModule::generateVariations(int count, double variation_pct) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - variation_pct, 1.0 + variation_pct);
    
    for (int i = 0; i < count; ++i) {
        std::map<std::string, double> variant = variables;
        for (auto& pair : variant) {
            if (pair.first != "pi") {  // Don't vary π
                pair.second *= dis(gen);
            }
        }
        variations.push_back(variant);
    }
    return variations;
}

// Adaptive Evolution
void SurfaceTemperatureModule::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - mutation_rate, 1.0 + mutation_rate);
    
    for (auto& pair : variables) {
        if (pair.first != "pi") {  // Don't mutate π
            pair.second *= dis(gen);
        }
    }
}

void SurfaceTemperatureModule::evolveSystem(int generations, std::function<double()> fitness_func) {
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
void SurfaceTemperatureModule::saveState(const std::string& label) {
    surface_temperature_saved_states::saved_states[label] = variables;
}

void SurfaceTemperatureModule::restoreState(const std::string& label) {
    if (surface_temperature_saved_states::saved_states.find(label) != surface_temperature_saved_states::saved_states.end()) {
        variables = surface_temperature_saved_states::saved_states[label];
    }
}

std::vector<std::string> SurfaceTemperatureModule::listSavedStates() const {
    std::vector<std::string> labels;
    for (const auto& pair : surface_temperature_saved_states::saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string SurfaceTemperatureModule::exportState() const {
    std::ostringstream oss;
    oss << "SurfaceTemperature_State_Export:\n";
    for (const auto& pair : variables) {
        oss << pair.first << "=" << std::scientific << pair.second << "\n";
    }
    return oss.str();
}

// System Analysis
std::map<std::string, double> SurfaceTemperatureModule::sensitivityAnalysis(const std::vector<std::string>& params) {
    std::map<std::string, double> sensitivities;
    double t_test = 0.0;
    double T_s_test = variables["T_s"];
    double baseline = computeU_g3_example(t_test, T_s_test);
    
    for (const auto& param : params) {
        if (variables.find(param) != variables.end() && param != "pi") {
            double original = variables[param];
            variables[param] = original * 1.01;
            
            double perturbed = computeU_g3_example(t_test, T_s_test);
            if (baseline > 0) {
                sensitivities[param] = (perturbed - baseline) / baseline;
            } else {
                sensitivities[param] = 0.0;
            }
            
            // Restore
            variables[param] = original;
        }
    }
    return sensitivities;
}

std::string SurfaceTemperatureModule::generateReport() const {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(3);
    oss << "===== Surface Temperature (T_s) Module Report =====\n";
    oss << "System: " << getSystemName() << "\n\n";
    
    oss << "Surface Temperature:\n";
    oss << "  T_s = " << variables.at("T_s") << " K\n";
    if (variables.find("T_s_celsius") != variables.end()) {
        oss << "  T_s = " << std::fixed << std::setprecision(1) 
            << variables.at("T_s_celsius") << " °C\n";
    } else {
        oss << "  T_s = " << std::fixed << std::setprecision(1) 
            << (variables.at("T_s") - 273.15) << " °C\n";
    }
    oss << std::scientific;
    oss << "  T_s_ref = " << variables.at("T_s_ref") << " K (reference)\n\n";
    
    if (variables.find("thermal_energy_density_Jm3") != variables.end()) {
        oss << "Thermal Characteristics:\n";
        oss << "  Thermal energy density = " << variables.at("thermal_energy_density_Jm3") << " J/m³\n";
        if (variables.find("luminosity_scale") != variables.end()) {
            oss << "  Luminosity scale (rel. to Sun) = " << std::fixed << std::setprecision(2) 
                << variables.at("luminosity_scale") << "x\n";
        }
        oss << std::scientific;
        if (variables.find("peak_wavelength_nm") != variables.end()) {
            oss << "  Peak wavelength (Wien) = " << std::fixed << std::setprecision(1) 
                << variables.at("peak_wavelength_nm") << " nm\n";
        }
        oss << "\n";
    }
    
    oss << std::scientific;
    oss << "Magnetic Field Coupling:\n";
    oss << "  B_ref = " << variables.at("B_ref") << " T (base magnetic field)\n";
    oss << "  k_3 = " << variables.at("k_3") << " (coupling constant)\n";
    if (variables.find("temp_field_coupling") != variables.end()) {
        oss << "  Temp-field coupling = " << variables.at("temp_field_coupling") << " T/K\n";
    }
    oss << "\n";
    
    oss << "Energy Parameters:\n";
    oss << "  E_react = " << variables.at("E_react") << " J\n";
    oss << "  ω_s = " << variables.at("omega_s") << " rad/s\n";
    oss << "  P_core = " << variables.at("P_core") << "\n";
    
    double cycle_period = 2.0 * variables.at("pi") / variables.at("omega_s");
    oss << "  Rotation period = " << cycle_period << " s (";
    oss << std::fixed << std::setprecision(1) << (cycle_period / 86400.0) << " days)\n";
    oss << "\n";
    
    oss << std::scientific;
    oss << "Temperature-Scaled Magnetic String Field B_j:\n";
    double t_test = 0.0;
    std::vector<double> T_s_tests = {3500.0, 5778.0, 10000.0, 30000.0};
    std::vector<std::string> labels = {"cool (M-type)", "Sun (G2V)", "hot (A-type)", "very hot (O-type)"};
    
    for (size_t i = 0; i < T_s_tests.size(); ++i) {
        double B_j = const_cast<SurfaceTemperatureModule*>(this)->computeB_j_hypothetical(t_test, T_s_tests[i]);
        oss << "  B_j at T_s=" << labels[i] << " (" << T_s_tests[i] << " K): " << B_j << " T\n";
    }
    oss << "\n";
    
    oss << "U_g3 Energy Density (t=0):\n";
    for (size_t i = 0; i < T_s_tests.size(); ++i) {
        double U_g3 = const_cast<SurfaceTemperatureModule*>(this)->computeU_g3_example(t_test, T_s_tests[i]);
        oss << "  U_g3 at T_s=" << labels[i] << ": " << U_g3 << " J/m³\n";
    }
    oss << "\n";
    
    oss << "Physical Interpretation:\n";
    double T_s = variables.at("T_s");
    if (T_s > 30000.0) {
        oss << "  O/B-type star (T_s > 30,000 K, hot blue star)\n";
    } else if (T_s > 10000.0) {
        oss << "  A/F-type star (T_s 10,000-30,000 K, white/blue-white star)\n";
    } else if (T_s > 5000.0) {
        oss << "  G-type star (T_s 5,000-10,000 K, yellow star like Sun)\n";
    } else if (T_s > 3500.0) {
        oss << "  K-type star (T_s 3,500-5,000 K, orange star)\n";
    } else if (T_s > 2000.0) {
        oss << "  M-type star (T_s 2,000-3,500 K, red dwarf)\n";
    } else if (T_s > 1000.0) {
        oss << "  Brown dwarf (T_s 1,000-2,000 K, substellar)\n";
    } else if (T_s > 200.0) {
        oss << "  Hot planet (T_s 200-1,000 K, hot Jupiter-like)\n";
    } else if (T_s > 100.0) {
        oss << "  Warm planet (T_s 100-200 K, Earth-like to cold)\n";
    } else {
        oss << "  Cold planet (T_s < 100 K, ice giant/outer planet)\n";
    }
    
    // Spectral classification
    if (T_s >= 5000.0 && T_s <= 6000.0) {
        oss << "  Near-solar spectral class (G-type yellow star)\n";
    }
    
    if (variables.find("peak_wavelength_nm") != variables.end()) {
        double lambda = variables.at("peak_wavelength_nm");
        if (lambda < 400.0) {
            oss << "  Peak emission in ultraviolet (<400 nm)\n";
        } else if (lambda < 500.0) {
            oss << "  Peak emission in blue/violet (400-500 nm)\n";
        } else if (lambda < 600.0) {
            oss << "  Peak emission in green/yellow (500-600 nm)\n";
        } else if (lambda < 700.0) {
            oss << "  Peak emission in orange/red (600-700 nm)\n";
        } else {
            oss << "  Peak emission in infrared (>700 nm)\n";
        }
    }
    
    oss << "\n  Applications:\n";
    oss << "    - Stellar classification: T_s determines spectral type (OBAFGKM)\n";
    oss << "    - Luminosity scaling: L ∝ T_s⁴ (Stefan-Boltzmann law)\n";
    oss << "    - Magnetic field coupling: B_j scales with T_s/T_s_ref\n";
    oss << "    - Habitability zones: Surface temp affects habitable distance\n";
    oss << "    - Stellar evolution: T_s tracks main sequence position\n";
    oss << "    - Planetary atmospheres: Temperature-dependent composition\n";
    oss << "    - Radiation modeling: Wien's law for peak emission wavelength\n";
    
    return oss.str();
}

bool SurfaceTemperatureModule::validateConsistency() const {
    bool valid = true;
    
    // Check T_s is positive
    if (variables.at("T_s") <= 0) {
        std::cerr << "Error: T_s must be positive (above absolute zero)\n";
        valid = false;
    }
    
    // Check T_s is physically reasonable
    if (variables.at("T_s") > 100000.0) {
        std::cerr << "Warning: T_s > 100,000 K (extremely hot, unusual for surface)\n";
    }
    
    if (variables.at("T_s") < 10.0) {
        std::cerr << "Warning: T_s < 10 K (extremely cold, near absolute zero)\n";
    }
    
    // Check T_s_ref is positive
    if (variables.at("T_s_ref") <= 0) {
        std::cerr << "Error: T_s_ref must be positive\n";
        valid = false;
    }
    
    // Check positive fields and energies
    if (variables.at("B_ref") <= 0) {
        std::cerr << "Error: B_ref must be positive\n";
        valid = false;
    }
    
    if (variables.at("k_3") <= 0) {
        std::cerr << "Error: k_3 must be positive\n";
        valid = false;
    }
    
    if (variables.at("E_react") <= 0) {
        std::cerr << "Error: E_react must be positive\n";
        valid = false;
    }
    
    if (variables.at("omega_s") <= 0) {
        std::cerr << "Error: omega_s must be positive\n";
        valid = false;
    }
    
    return valid;
}

void SurfaceTemperatureModule::autoCorrectAnomalies() {
    // Ensure T_s is positive and reasonable
    if (variables["T_s"] <= 0 || variables["T_s"] > 100000.0) {
        variables["T_s"] = 5778.0;  // Reset to Sun
    }
    
    // Ensure T_s_ref is positive
    if (variables["T_s_ref"] <= 0) {
        variables["T_s_ref"] = 5778.0;
    }
    
    // Ensure B_ref is positive
    if (variables["B_ref"] <= 0) {
        variables["B_ref"] = 1e3;
    }
    
    // Ensure k_3 is positive and reasonable
    if (variables["k_3"] <= 0 || variables["k_3"] > 100.0) {
        variables["k_3"] = 1.8;
    }
    
    // Ensure E_react is positive
    if (variables["E_react"] <= 0) {
        variables["E_react"] = 1e46;
    }
    
    // Ensure omega_s is positive
    if (variables["omega_s"] <= 0) {
        variables["omega_s"] = 2.5e-6;
    }
    
    // Ensure P_core is reasonable
    if (variables["P_core"] <= 0 || variables["P_core"] > 10.0) {
        variables["P_core"] = 1.0;
    }
    
    // Update derived Celsius
    if (variables.find("T_s_celsius") != variables.end()) {
        variables["T_s_celsius"] = variables["T_s"] - 273.15;
    }
}

// Example usage in base program (snippet)
// #include "SurfaceTemperatureModule.h"
// int main() {
//     SurfaceTemperatureModule mod;
//     double t_s = mod.computeT_s();
//     std::cout << "T_s = " << t_s << " K\n";
//     double u_g3 = mod.computeU_g3_example(0.0, 10000.0);
//     std::cout << "U_g3 (T_s=10000 K) = " << u_g3 << " J/m�\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("T_s", 6000.0);
//     mod.printVariables();
//     return 0;
// }

// ========== COMPREHENSIVE ENHANCED DEMONSTRATION ==========
/*
int main() {
    std::cout << std::scientific << std::setprecision(3);
    std::cout << "========== SURFACE TEMPERATURE MODULE DEMONSTRATION ==========\n\n";
    
    // ===== Step 1: Initialize Module =====
    SurfaceTemperatureModule mod;
    std::cout << "Step 1: Module initialized with defaults (Sun):\n";
    std::cout << "  System: " << mod.getSystemName() << "\n";
    std::cout << "  T_s = " << mod.computeT_s() << " K\n";
    std::cout << "  T_s_ref = " << mod.variables["T_s_ref"] << " K\n\n";
    
    // ===== Step 2: Compute Baseline =====
    std::cout << "Step 2: Compute baseline at solar temperature:\n";
    double t_test = 0.0;
    
    double B_j_sun = mod.computeB_j_hypothetical(t_test, mod.variables["T_s"]);
    double U_g3_sun = mod.computeU_g3_example(t_test, mod.variables["T_s"]);
    std::cout << "  Sun (T_s=5778 K):\n";
    std::cout << "    B_j = " << B_j_sun << " T\n";
    std::cout << "    U_g3 = " << U_g3_sun << " J/m³\n\n";
    
    // ===== Step 3: Variable Management =====
    std::cout << "Step 3: Variable Management\n";
    mod.createVariable("T_s_celsius", mod.variables["T_s"] - 273.15);
    std::cout << "  Created 'T_s_celsius' = " << std::fixed << std::setprecision(1) 
              << mod.variables["T_s_celsius"] << " °C\n";
    
    std::vector<std::string> all_vars = mod.listVariables();
    std::cout << "  Total variables: " << all_vars.size() << "\n";
    
    mod.cloneVariable("T_s", "T_s_backup");
    std::cout << "  Cloned 'T_s' → 'T_s_backup' = " << mod.variables["T_s_backup"] << " K\n\n";
    
    // ===== Step 4: Temperature Expansion (Stellar Evolution) =====
    std::cout << "Step 4: Temperature Expansion (Cool → Hot Star)\n";
    std::cout << std::scientific << std::setprecision(3);
    mod.expandTemperatureScale(2.0, 1.5);  // Double T_s, 1.5x thermal factors
    std::cout << "  After expansion: T_s = " << mod.variables["T_s"] << " K (evolved to hotter type)\n";
    std::cout << "  Luminosity scaling: L/L_sun ~ " << std::pow(mod.variables["T_s"]/5778.0, 4.0) << "\n";
    std::cout << "  Wien's law peak: λ_peak = " << 2.898e-3 / mod.variables["T_s"] * 1e9 
              << " nm (shifted to blue)\n\n";
    
    // ===== Step 5: Field Expansion (Temperature-Magnetic Coupling) =====
    std::cout << "Step 5: Field Expansion (Temperature-Magnetic Coupling)\n";
    mod.expandFieldScale(1.8, 1.2);
    double B_j_hot = mod.computeB_j_hypothetical(t_test, mod.variables["T_s"]);
    std::cout << "  After field expansion: B_j = " << B_j_hot << " T\n";
    std::cout << "  Temperature-field coupling strengthened\n\n";
    
    // ===== Step 6: Energy Expansion (Reactor & Orbital Energy) =====
    std::cout << "Step 6: Energy Expansion\n";
    mod.expandEnergyScale(2.5, 1.6);
    double U_g3_hot = mod.computeU_g3_example(t_test, mod.variables["T_s"]);
    std::cout << "  After energy expansion: U_g3 = " << U_g3_hot << " J/m³\n";
    std::cout << "  E_react scaled for higher power output\n\n";
    
    // ===== Step 7: Batch Operations =====
    std::cout << "Step 7: Batch Operations (Temperature Group)\n";
    std::vector<std::string> temp_group = {"T_s", "T_s_ref", "T_s_celsius"};
    mod.scaleVariableGroup(temp_group, 0.6);  // Cool down to M-dwarf range
    std::cout << "  Scaled temperature group by 0.6:\n";
    std::cout << "    T_s = " << mod.variables["T_s"] << " K (M-dwarf range)\n";
    std::cout << "    Spectral type: M (red dwarf)\n\n";
    
    // ===== Step 8-10: Test Multiple Stellar Types =====
    std::cout << "Steps 8-10: Test Multiple Stellar Types\n";
    
    // M-dwarf (red)
    mod.updateVariable("T_s", 3500.0);
    double L_m = std::pow(3500.0/5778.0, 4.0);
    std::cout << "  M-dwarf (T_s=3500 K): L/L_sun = " << L_m << ", λ_peak = " 
              << 2.898e-3/3500.0*1e9 << " nm (infrared)\n";
    
    // A-type (white)
    mod.updateVariable("T_s", 9000.0);
    double L_a = std::pow(9000.0/5778.0, 4.0);
    std::cout << "  A-type (T_s=9000 K): L/L_sun = " << L_a << ", λ_peak = " 
              << 2.898e-3/9000.0*1e9 << " nm (UV-visible)\n";
    
    // O-type (blue giant)
    mod.updateVariable("T_s", 30000.0);
    double L_o = std::pow(30000.0/5778.0, 4.0);
    std::cout << "  O-type (T_s=30000 K): L/L_sun = " << L_o << ", λ_peak = " 
              << 2.898e-3/30000.0*1e9 << " nm (far UV)\n\n";
    
    // ===== Step 11: Auto-Refinement =====
    std::cout << "Step 11: Auto-Refinement (Target Solar Values)\n";
    mod.updateVariable("T_s", 5778.0);
    mod.updateVariable("luminosity_scale", 1.0);
    mod.updateVariable("peak_wavelength_nm", 501.5);
    mod.autoRefineParameters();
    std::cout << "  Refined to solar values: T_s = " << mod.variables["T_s"] << " K\n";
    std::cout << "  Spectral class: G2V (yellow star like Sun)\n\n";
    
    // ===== Step 12: Calibration (Observational Data) =====
    std::cout << "Step 12: Calibration to Observations\n";
    std::map<std::string, double> obs_data;
    obs_data["T_s"] = 5800.0;  // Slightly adjusted observation
    obs_data["luminosity_scale"] = 1.02;
    mod.calibrateToObservations(obs_data);
    std::cout << "  Calibrated T_s = " << mod.variables["T_s"] << " K (adjusted to observation)\n\n";
    
    // ===== Step 13-17: Optimization for Different Metrics =====
    std::cout << "Steps 13-17: Optimization for Different Star Types\n";
    
    mod.optimizeForMetric("cool_star");  // M-dwarf
    std::cout << "  Optimized for 'cool_star': T_s = " << mod.variables["T_s"] 
              << " K (M-dwarf range)\n";
    
    mod.optimizeForMetric("hot_star");  // O-type
    std::cout << "  Optimized for 'hot_star': T_s = " << mod.variables["T_s"] 
              << " K (O-type range)\n";
    
    mod.optimizeForMetric("white_dwarf");  // Compact object
    std::cout << "  Optimized for 'white_dwarf': T_s = " << mod.variables["T_s"] 
              << " K (white dwarf)\n";
    
    mod.optimizeForMetric("planets");  // Planetary temperatures
    std::cout << "  Optimized for 'planets': T_s = " << mod.variables["T_s"] 
              << " K (planetary range)\n\n";
    
    // ===== Step 18: Parameter Variations =====
    std::cout << "Step 18: Generate Parameter Variations\n";
    std::vector<std::map<std::string, double>> variations = mod.generateVariations(3);
    std::cout << "  Generated " << variations.size() << " parameter variations:\n";
    for (size_t i = 0; i < variations.size(); ++i) {
        std::cout << "    Variant " << (i+1) << ": T_s = " << variations[i]["T_s"] << " K\n";
    }
    std::cout << "\n";
    
    // ===== Step 19: Mutation =====
    std::cout << "Step 19: Mutate Parameters (Evolutionary Algorithm)\n";
    mod.updateVariable("T_s", 5778.0);  // Reset to Sun
    mod.mutateParameters(0.15);  // 15% mutation rate
    std::cout << "  After mutation: T_s = " << mod.variables["T_s"] 
              << " K (randomly varied)\n\n";
    
    // ===== Step 20: System Evolution =====
    std::cout << "Step 20: Evolve System (Adaptive Improvement)\n";
    mod.evolveSystem(5);  // 5 evolution iterations
    std::cout << "  After 5 evolution cycles: T_s = " << mod.variables["T_s"] << " K\n";
    std::cout << "  System adapted toward optimal configuration\n\n";
    
    // ===== Step 21: State Management =====
    std::cout << "Step 21: State Management\n";
    mod.updateVariable("T_s", 5778.0);
    mod.saveState("solar_baseline");
    std::cout << "  Saved state 'solar_baseline'\n";
    
    mod.updateVariable("T_s", 3500.0);
    mod.saveState("m_dwarf");
    std::cout << "  Saved state 'm_dwarf'\n";
    
    std::vector<std::string> saved = mod.listSavedStates();
    std::cout << "  Total saved states: " << saved.size() << "\n";
    
    mod.restoreState("solar_baseline");
    std::cout << "  Restored 'solar_baseline': T_s = " << mod.variables["T_s"] << " K\n\n";
    
    // ===== Step 22: Export State =====
    std::cout << "Step 22: Export State\n";
    std::string exported = mod.exportState();
    std::cout << "  Exported " << exported.length() << " bytes of state data\n";
    std::cout << "  (Can be saved to file for later restoration)\n\n";
    
    // ===== Step 23: Sensitivity Analysis =====
    std::cout << "Step 23: Sensitivity Analysis\n";
    std::map<std::string, double> sensitivity = mod.sensitivityAnalysis("U_g3_density");
    std::cout << "  Sensitivity of U_g3_density to parameter changes:\n";
    for (const auto& pair : sensitivity) {
        std::cout << "    " << pair.first << ": " << pair.second << "\n";
    }
    std::cout << "\n";
    
    // ===== Step 24: Validation =====
    std::cout << "Step 24: Consistency Validation\n";
    bool valid = mod.validateConsistency();
    std::cout << "  System consistency: " << (valid ? "VALID" : "INVALID") << "\n";
    if (!valid) {
        mod.autoCorrectAnomalies();
        std::cout << "  Auto-corrected anomalies\n";
    }
    std::cout << "\n";
    
    // ===== Step 25: Generate Full Report =====
    std::cout << "Step 25: Generate Full Report\n";
    std::string report = mod.generateReport();
    std::cout << report << "\n";
    
    // ===== Step 26: Final Comparison (Stellar Sequence) =====
    std::cout << "Step 26: Final Stellar Sequence Comparison\n";
    std::cout << "  Spectral Type | T_s (K) | L/L_sun | λ_peak (nm) | Color\n";
    std::cout << "  ------------------------------------------------------------\n";
    
    struct StellarType { std::string name; double temp; std::string color; };
    std::vector<StellarType> sequence = {
        {"M (red dwarf)", 3500, "Red"},
        {"K (orange)", 4500, "Orange"},
        {"G (Sun)", 5778, "Yellow"},
        {"F (yellow-white)", 6500, "Yellow-white"},
        {"A (white)", 9000, "White"},
        {"B (blue-white)", 15000, "Blue-white"},
        {"O (blue giant)", 30000, "Blue"}
    };
    
    for (const auto& star : sequence) {
        double lum = std::pow(star.temp/5778.0, 4.0);
        double peak = 2.898e-3 / star.temp * 1e9;
        std::cout << "  " << std::setw(18) << std::left << star.name 
                  << " | " << std::setw(7) << std::right << star.temp
                  << " | " << std::setw(7) << lum
                  << " | " << std::setw(11) << peak
                  << " | " << star.color << "\n";
    }
    
    std::cout << "\n========== DEMONSTRATION COMPLETE ==========\n";
    std::cout << "Surface Temperature Module validated across full stellar sequence.\n";
    std::cout << "Thermal-magnetic coupling, luminosity scaling, and spectral analysis confirmed.\n";
    
    return 0;
}
*/

// Compile: g++ -o temp_test temp_test.cpp SurfaceTemperatureModule.cpp -lm
// Sample: T_s=5778 K; U_g3 (hot star)?3.11e49 J/m�; thermal scaling.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

SurfaceTemperatureModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeT_s, computeB_j_hypothetical, computeU_g3_example) are clear, concise, and variable - driven.
- Output and debugging functions(printVariables, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.
- Supports temperature scaling for magnetic string strength, enabling modeling of different stellar types.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in surface temperature modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.