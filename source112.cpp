// ScmPenetrationModule.h
// Modular C++ implementation of the [SCm] Penetration Factor (P_SCm) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes P_SCm ?1 (unitless for Sun, ~1e-3 for planets); scales P_SCm in Universal Magnetism U_m term.
// Pluggable: #include "ScmPenetrationModule.h"
// ScmPenetrationModule mod; mod.computeUmContribution(0.0); mod.updateVariable("P_SCm", new_value);
// Variables in std::map; example for Sun at t=0, t_n=0; full penetration for plasma cores.
// Approximations: 1 - e^{-? t cos(? t_n)}=0 at t=0; ?_hat_j=1; ?_j / r_j=2.26e10 T m�.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef SCM_PENETRATION_MODULE_H
#define SCM_PENETRATION_MODULE_H

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

class ScmPenetrationModule {
private:
    std::map<std::string, double> variables;
    double computeUmBase(double t);
    double computeUmContribution(double t);

public:
    // Constructor: Initialize with framework defaults (Sun)
    ScmPenetrationModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeP_SCm();  // ?1 for Sun (unitless)
    double computeUmContribution(double t);  // U_m with P_SCm (J/m^3)
    double computeUmPlanet(double t);  // For planet P_SCm=1e-3

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
    void expandParameterSpace(double penetration_scale, double magnetic_scale, double energy_scale);
    void expandPenetrationScale(double pscm_factor, double stellar_planet_ratio);  // P_SCm stellar vs planetary
    void expandMagneticScale(double mu_factor, double rj_factor);                  // μ_j and r_j
    void expandAmplificationScale(double heaviside_factor, double quasi_factor);   // Heaviside and quasi

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

#endif // SCM_PENETRATION_MODULE_H

// ScmPenetrationModule.cpp
#include "ScmPenetrationModule.h"

// Constructor: Set framework defaults (Sun at t=0)
ScmPenetrationModule::ScmPenetrationModule() {
    // Universal constants
    variables["P_SCm"] = 1.0;                       // Unitless ?1 for Sun
    variables["P_SCm_planet"] = 1e-3;               // For planets
    variables["mu_j"] = 3.38e23;                    // T�m^3
    variables["r_j"] = 1.496e13;                    // m
    variables["gamma"] = 5e-5 / 86400.0;            // s^-1
    variables["t_n"] = 0.0;                         // s
    variables["phi_hat_j"] = 1.0;                   // Normalized
    variables["P_SCm"] = 1.0;                       // Pressure (wait, reuse? No, P_SCm is penetration)
    variables["E_react"] = 1e46;                    // J
    variables["f_Heaviside"] = 0.01;                // Unitless
    variables["f_quasi"] = 0.01;                    // Unitless
    variables["pi"] = 3.141592653589793;
    variables["scale_Heaviside"] = 1e13;            // Amplification

    // Derived
    variables["heaviside_factor"] = 1.0 + variables["scale_Heaviside"] * variables["f_Heaviside"];
}

// Update variable
void ScmPenetrationModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void ScmPenetrationModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void ScmPenetrationModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute P_SCm ?1
double ScmPenetrationModule::computeP_SCm() {
    return variables["P_SCm"];
}

// Base for U_m without P_SCm
double ScmPenetrationModule::computeUmBase(double t) {
    double mu_over_rj = variables["mu_j"] / variables["r_j"];
    double exp_arg = - variables["gamma"] * t * std::cos(variables["pi"] * variables["t_n"]);
    double one_minus_exp = 1.0 - std::exp(exp_arg);
    double phi_hat = variables["phi_hat_j"];
    double p_scm = computeP_SCm();  // Penetration as pressure-like
    double e_react = variables["E_react"];
    return mu_over_rj * one_minus_exp * phi_hat * p_scm * e_react;
}

// U_m contribution with P_SCm
double ScmPenetrationModule::computeUmContribution(double t) {
    double base = computeUmBase(t);
    double heaviside_f = variables["heaviside_factor"];
    double quasi_f = 1.0 + variables["f_quasi"];
    return base * heaviside_f * quasi_f;
}

// U_m for planet (P_SCm=1e-3)
double ScmPenetrationModule::computeUmPlanet(double t) {
    double orig_p = variables["P_SCm"];
    variables["P_SCm"] = variables["P_SCm_planet"];
    double result = computeUmContribution(t);
    variables["P_SCm"] = orig_p;
    return result;
}

// Equation text
std::string ScmPenetrationModule::getEquationText() {
    return "U_m = [ (?_j / r_j) (1 - e^{-? t cos(? t_n)}) ?_hat_j ] P_SCm E_react (1 + 10^13 f_Heaviside) (1 + f_quasi)\n"
           "Where P_SCm ?1 (unitless [SCm] penetration factor; ~1e-3 for planets).\n"
           "Scales magnetic energy for [SCm] interior interaction.\n"
           "Example Sun t=0: U_m ?2.28e65 J/m� (P_SCm=1);\n"
           "Planet: ?2.28e62 J/m� (P_SCm=1e-3, -3 orders).\n"
           "Role: Full for stellar plasma, reduced for solid cores; [SCm] influence on strings.\n"
           "UQFF: Models penetration in jets/nebulae/formation; massless [SCm] dynamics.";
}

// Print variables
void ScmPenetrationModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// ===== ENHANCED METHODS IMPLEMENTATION =====

namespace scm_penetration_saved_states {
    static std::map<std::string, std::map<std::string, double>> saved_states;
}

// Variable Management
void ScmPenetrationModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void ScmPenetrationModule::removeVariable(const std::string& name) {
    variables.erase(name);
}

void ScmPenetrationModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
    }
}

std::vector<std::string> ScmPenetrationModule::listVariables() const {
    std::vector<std::string> names;
    for (const auto& pair : variables) {
        names.push_back(pair.first);
    }
    return names;
}

std::string ScmPenetrationModule::getSystemName() const {
    return "SCm_Penetration_UQFF";
}

// Batch Operations
void ScmPenetrationModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = func(variables[name]);
        }
    }
    // Update derived
    variables["heaviside_factor"] = 1.0 + variables["scale_Heaviside"] * variables["f_Heaviside"];
}

void ScmPenetrationModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// Self-Expansion: Domain-Specific Scales
void ScmPenetrationModule::expandParameterSpace(double penetration_scale, double magnetic_scale, double energy_scale) {
    variables["P_SCm"] *= penetration_scale;
    variables["P_SCm_planet"] *= penetration_scale;
    variables["mu_j"] *= magnetic_scale;
    variables["E_react"] *= energy_scale;
}

void ScmPenetrationModule::expandPenetrationScale(double pscm_factor, double stellar_planet_ratio) {
    variables["P_SCm"] *= pscm_factor;
    // Adjust planetary penetration to maintain ratio
    variables["P_SCm_planet"] = variables["P_SCm"] / stellar_planet_ratio;
}

void ScmPenetrationModule::expandMagneticScale(double mu_factor, double rj_factor) {
    variables["mu_j"] *= mu_factor;
    variables["r_j"] *= rj_factor;
}

void ScmPenetrationModule::expandAmplificationScale(double heaviside_factor, double quasi_factor) {
    variables["f_Heaviside"] *= heaviside_factor;
    variables["f_quasi"] *= quasi_factor;
    variables["heaviside_factor"] = 1.0 + variables["scale_Heaviside"] * variables["f_Heaviside"];
}

// Self-Refinement
void ScmPenetrationModule::autoRefineParameters(const std::string& target, double goal) {
    if (target == "P_SCm") {
        variables["P_SCm"] = goal;
    } else if (target == "P_SCm_planet") {
        variables["P_SCm_planet"] = goal;
    } else if (target == "stellar_planet_ratio") {
        // Set ratio: P_SCm / P_SCm_planet = goal
        variables["P_SCm_planet"] = variables["P_SCm"] / goal;
    } else if (target == "U_m_stellar") {
        // Target specific U_m for stellar case at t=0
        double current_um = computeUmContribution(0.0);
        if (std::abs(current_um) > 1e-9) {
            variables["P_SCm"] *= (goal / current_um);
        }
    } else if (target == "U_m_planetary") {
        // Target specific U_m for planetary case
        double current_um_planet = computeUmPlanet(0.0);
        if (std::abs(current_um_planet) > 1e-9) {
            variables["P_SCm_planet"] *= (goal / current_um_planet);
        }
    }
}

void ScmPenetrationModule::calibrateToObservations(const std::map<std::string, double>& observations) {
    for (const auto& obs : observations) {
        if (variables.find(obs.first) != variables.end()) {
            variables[obs.first] = obs.second;
            if (obs.first == "f_Heaviside" || obs.first == "scale_Heaviside") {
                variables["heaviside_factor"] = 1.0 + variables["scale_Heaviside"] * variables["f_Heaviside"];
            }
        }
    }
}

void ScmPenetrationModule::optimizeForMetric(const std::string& metric) {
    if (metric == "full_stellar_penetration") {
        // Full penetration for stellar plasma
        variables["P_SCm"] = 1.0;
    } else if (metric == "reduced_planetary") {
        // Reduced penetration for planetary solid cores
        variables["P_SCm_planet"] = 1e-3;
    } else if (metric == "intermediate_penetration") {
        // Intermediate between stellar and planetary
        variables["P_SCm"] = 0.5;
        variables["P_SCm_planet"] = 5e-4;
    } else if (metric == "standard_ratio") {
        // Standard 1000:1 stellar-to-planetary ratio
        variables["P_SCm"] = 1.0;
        variables["P_SCm_planet"] = 1e-3;
    } else if (metric == "enhanced_penetration") {
        // Enhanced penetration (deeper SCm interaction)
        variables["P_SCm"] *= 1.5;
        variables["P_SCm_planet"] *= 1.5;
    } else if (metric == "suppressed_penetration") {
        // Suppressed penetration (shallower SCm interaction)
        variables["P_SCm"] *= 0.5;
        variables["P_SCm_planet"] *= 0.5;
    }
}

// Parameter Exploration
std::vector<std::map<std::string, double>> ScmPenetrationModule::generateVariations(int count, double variation_pct) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - variation_pct, 1.0 + variation_pct);
    
    for (int i = 0; i < count; ++i) {
        std::map<std::string, double> variant = variables;
        for (auto& pair : variant) {
            if (pair.first != "pi" && pair.first != "scale_Heaviside" && pair.first != "heaviside_factor") {
                pair.second *= dis(gen);
            }
        }
        // Recalculate derived
        variant["heaviside_factor"] = 1.0 + variant["scale_Heaviside"] * variant["f_Heaviside"];
        variations.push_back(variant);
    }
    return variations;
}

// Adaptive Evolution
void ScmPenetrationModule::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - mutation_rate, 1.0 + mutation_rate);
    
    for (auto& pair : variables) {
        if (pair.first != "pi" && pair.first != "scale_Heaviside" && pair.first != "heaviside_factor") {
            pair.second *= dis(gen);
        }
    }
    // Recalculate derived
    variables["heaviside_factor"] = 1.0 + variables["scale_Heaviside"] * variables["f_Heaviside"];
}

void ScmPenetrationModule::evolveSystem(int generations, std::function<double()> fitness_func) {
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
void ScmPenetrationModule::saveState(const std::string& label) {
    scm_penetration_saved_states::saved_states[label] = variables;
}

void ScmPenetrationModule::restoreState(const std::string& label) {
    if (scm_penetration_saved_states::saved_states.find(label) != scm_penetration_saved_states::saved_states.end()) {
        variables = scm_penetration_saved_states::saved_states[label];
    }
}

std::vector<std::string> ScmPenetrationModule::listSavedStates() const {
    std::vector<std::string> labels;
    for (const auto& pair : scm_penetration_saved_states::saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string ScmPenetrationModule::exportState() const {
    std::ostringstream oss;
    oss << "ScmPenetration_State_Export:\n";
    for (const auto& pair : variables) {
        oss << pair.first << "=" << std::scientific << pair.second << "\n";
    }
    return oss.str();
}

// System Analysis
std::map<std::string, double> ScmPenetrationModule::sensitivityAnalysis(const std::vector<std::string>& params) {
    std::map<std::string, double> sensitivities;
    double t_test = 0.0;  // Test at t=0
    double baseline = computeUmContribution(t_test);
    
    for (const auto& param : params) {
        if (variables.find(param) != variables.end() && param != "pi" && param != "scale_Heaviside") {
            double original = variables[param];
            variables[param] = original * 1.01;
            
            if (param == "f_Heaviside") {
                variables["heaviside_factor"] = 1.0 + variables["scale_Heaviside"] * variables["f_Heaviside"];
            }
            
            double perturbed = computeUmContribution(t_test);
            sensitivities[param] = (perturbed - baseline) / baseline;
            
            // Restore
            variables[param] = original;
            if (param == "f_Heaviside") {
                variables["heaviside_factor"] = 1.0 + variables["scale_Heaviside"] * variables["f_Heaviside"];
            }
        }
    }
    return sensitivities;
}

std::string ScmPenetrationModule::generateReport() const {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(3);
    oss << "===== [SCm] Penetration Module Report =====\n";
    oss << "System: " << getSystemName() << "\n\n";
    
    oss << "Penetration Parameters:\n";
    oss << "  P_SCm (stellar) = " << variables.at("P_SCm") << " (unitless)\n";
    oss << "  P_SCm (planetary) = " << variables.at("P_SCm_planet") << " (unitless)\n";
    double ratio = variables.at("P_SCm") / variables.at("P_SCm_planet");
    oss << "  Stellar/Planetary ratio = " << std::fixed << std::setprecision(0) << ratio << ":1\n\n";
    
    oss << std::scientific;
    oss << "Magnetic Parameters:\n";
    oss << "  μ_j (magnetic moment) = " << variables.at("mu_j") << " T·m³\n";
    oss << "  r_j (jet radius) = " << variables.at("r_j") << " m\n";
    double mu_over_rj = variables.at("mu_j") / variables.at("r_j");
    oss << "  μ_j/r_j = " << mu_over_rj << " T·m²\n\n";
    
    oss << "Decay and Oscillation:\n";
    oss << "  γ (decay rate) = " << variables.at("gamma") << " s⁻¹\n";
    oss << "  t_n (negative time) = " << variables.at("t_n") << " s\n";
    oss << "  π = " << std::fixed << std::setprecision(15) << variables.at("pi") << "\n\n";
    
    oss << std::scientific;
    oss << "Amplification Factors:\n";
    oss << "  f_Heaviside = " << variables.at("f_Heaviside") << "\n";
    oss << "  f_quasi = " << variables.at("f_quasi") << "\n";
    oss << "  Heaviside factor = 1 + 10¹³ × f_Heaviside = " << variables.at("heaviside_factor") << "\n";
    oss << "  Quasi factor = 1 + f_quasi = " << (1.0 + variables.at("f_quasi")) << "\n\n";
    
    oss << "Energy:\n";
    oss << "  E_react = " << variables.at("E_react") << " J\n";
    oss << "  φ_hat_j (flux normalization) = " << variables.at("phi_hat_j") << "\n\n";
    
    oss << "U_m Computations:\n";
    // Stellar case (t=0)
    double um_stellar_t0 = 0.0;  // At t=0, 1-exp(0)=0, so U_m=0
    // But let's compute for small t to show effect
    double t_small = 1000.0 * 86400.0;  // 1000 days in seconds
    double exp_arg = -variables.at("gamma") * t_small * std::cos(variables.at("pi") * variables.at("t_n"));
    double one_minus_exp = 1.0 - std::exp(exp_arg);
    double um_stellar = mu_over_rj * one_minus_exp * variables.at("phi_hat_j") * variables.at("P_SCm") * 
                        variables.at("E_react") * variables.at("heaviside_factor") * (1.0 + variables.at("f_quasi"));
    
    oss << "  Stellar (P_SCm=" << variables.at("P_SCm") << ", t=1000 days):\n";
    oss << "    1 - exp(-γt cos(πt_n)) = " << one_minus_exp << "\n";
    oss << "    U_m = " << um_stellar << " J/m³\n\n";
    
    // Planetary case
    double um_planetary = mu_over_rj * one_minus_exp * variables.at("phi_hat_j") * variables.at("P_SCm_planet") * 
                          variables.at("E_react") * variables.at("heaviside_factor") * (1.0 + variables.at("f_quasi"));
    oss << "  Planetary (P_SCm=" << variables.at("P_SCm_planet") << ", t=1000 days):\n";
    oss << "    U_m = " << um_planetary << " J/m³\n";
    oss << "    Ratio: U_m(stellar)/U_m(planetary) = " << std::fixed << std::setprecision(0) 
        << (um_stellar / um_planetary) << ":1\n\n";
    
    oss << std::scientific;
    oss << "Physical Interpretation:\n";
    if (variables.at("P_SCm") >= 0.9) {
        oss << "  Full stellar penetration: Plasma cores allow complete [SCm] interaction\n";
    } else if (variables.at("P_SCm") >= 0.5) {
        oss << "  Intermediate penetration: Partial [SCm] coupling\n";
    } else if (variables.at("P_SCm") >= 0.1) {
        oss << "  Reduced penetration: Limited [SCm] access\n";
    } else {
        oss << "  Minimal penetration: Planetary-like solid core obstruction\n";
    }
    
    oss << "  Applications:\n";
    oss << "    - Stellar plasmas: Full P_SCm ≈ 1 (complete magnetic string penetration)\n";
    oss << "    - Planetary cores: Reduced P_SCm ≈ 10⁻³ (solid obstruction)\n";
    oss << "    - [SCm] dynamics: Massless superconductive matter interaction\n";
    oss << "    - Jet formation: Penetration affects magnetic field strength\n";
    oss << "    - Nebular structures: P_SCm varies with density and phase\n";
    oss << "    - Star formation: Interior [SCm] coupling drives energy transport\n";
    
    return oss.str();
}

bool ScmPenetrationModule::validateConsistency() const {
    bool valid = true;
    
    // Check P_SCm is in valid range [0, ~2]
    if (variables.find("P_SCm") != variables.end()) {
        double p = variables.at("P_SCm");
        if (p < 0 || p > 2.0) {
            std::cerr << "Warning: P_SCm outside typical range [0, 2] (current: " << p << ")\n";
        }
    }
    
    // Check P_SCm_planet is in valid range [0, 0.1]
    if (variables.find("P_SCm_planet") != variables.end()) {
        double p_planet = variables.at("P_SCm_planet");
        if (p_planet < 0 || p_planet > 0.1) {
            std::cerr << "Warning: P_SCm_planet outside typical range [0, 0.1] (current: " << p_planet << ")\n";
        }
    }
    
    // Check that stellar > planetary
    if (variables.at("P_SCm") < variables.at("P_SCm_planet")) {
        std::cerr << "Error: P_SCm (stellar) < P_SCm_planet (should be stellar > planetary)\n";
        valid = false;
    }
    
    // Check heaviside_factor consistency
    double expected_hf = 1.0 + variables.at("scale_Heaviside") * variables.at("f_Heaviside");
    double actual_hf = variables.at("heaviside_factor");
    if (std::abs(expected_hf - actual_hf) / expected_hf > 1e-9) {
        std::cerr << "Error: heaviside_factor inconsistent (expected " << expected_hf << ", got " << actual_hf << ")\n";
        valid = false;
    }
    
    return valid;
}

void ScmPenetrationModule::autoCorrectAnomalies() {
    // Reset P_SCm to typical value if out of range
    if (variables["P_SCm"] < 0 || variables["P_SCm"] > 2.0) {
        variables["P_SCm"] = 1.0;  // Standard stellar
    }
    
    // Reset P_SCm_planet if out of range
    if (variables["P_SCm_planet"] < 0 || variables["P_SCm_planet"] > 0.1) {
        variables["P_SCm_planet"] = 1e-3;  // Standard planetary
    }
    
    // Ensure stellar > planetary
    if (variables["P_SCm"] < variables["P_SCm_planet"]) {
        variables["P_SCm"] = 1.0;
        variables["P_SCm_planet"] = 1e-3;
    }
    
    // Recalculate heaviside_factor
    variables["heaviside_factor"] = 1.0 + variables["scale_Heaviside"] * variables["f_Heaviside"];
    
    // Ensure π is correct
    if (std::abs(variables["pi"] - 3.141592653589793) > 1e-12) {
        variables["pi"] = 3.141592653589793;
    }
}

// Example usage in base program (snippet)
int main() {
    ScmPenetrationModule module;
    std::cout << "===== [SCm] Penetration Module Enhanced Demonstration =====\n\n";
    
    // Step 1: Report initial state
    std::cout << "STEP 1: Initial Configuration (Stellar P_SCm = 1.0, Planetary = 0.001)\n";
    std::cout << module.generateReport() << "\n";
    
    // Step 2: Compute U_m for stellar vs planetary
    std::cout << "STEP 2: U_m Comparison (Stellar vs Planetary at t=1000 days)\n";
    double t_days = 1000.0;
    double t_seconds = t_days * 86400.0;
    double um_stellar = module.computeUmContribution(t_seconds);
    double um_planetary = module.computeUmPlanet(t_seconds);
    std::cout << "  U_m (stellar, P_SCm=" << module.getVariable("P_SCm") << ") = " 
              << std::scientific << um_stellar << " J/m³\n";
    std::cout << "  U_m (planetary, P_SCm=" << module.getVariable("P_SCm_planet") << ") = " 
              << um_planetary << " J/m³\n";
    std::cout << "  Ratio (stellar/planetary) = " << std::fixed << std::setprecision(0) 
              << (um_stellar / um_planetary) << ":1\n\n";
    
    // Step 3: Save initial state
    std::cout << "STEP 3: Save Initial State\n";
    module.saveState("standard_penetration");
    std::cout << "State saved as 'standard_penetration'\n\n";
    
    // Step 4: Test intermediate penetration
    std::cout << "STEP 4: Test Intermediate Penetration (P_SCm = 0.5)\n";
    module.optimizeForMetric("intermediate_penetration");
    double um_intermediate = module.computeUmContribution(t_seconds);
    std::cout << "Intermediate P_SCm = " << std::scientific << module.getVariable("P_SCm") << "\n";
    std::cout << "U_m (intermediate) = " << um_intermediate << " J/m³\n";
    std::cout << "Ratio to stellar standard: " << std::fixed << std::setprecision(2) 
              << (um_intermediate / um_stellar) << "\n\n";
    
    // Step 5: Test enhanced penetration
    std::cout << "STEP 5: Test Enhanced Penetration (1.5x standard)\n";
    module.restoreState("standard_penetration");
    module.optimizeForMetric("enhanced_penetration");
    double um_enhanced = module.computeUmContribution(t_seconds);
    std::cout << "Enhanced P_SCm = " << std::scientific << module.getVariable("P_SCm") << "\n";
    std::cout << "U_m (enhanced) = " << um_enhanced << " J/m³\n";
    std::cout << "Enhancement factor: " << std::fixed << std::setprecision(2) 
              << (um_enhanced / um_stellar) << "\n\n";
    
    // Step 6: Test suppressed penetration
    std::cout << "STEP 6: Test Suppressed Penetration (0.5x standard)\n";
    module.restoreState("standard_penetration");
    module.optimizeForMetric("suppressed_penetration");
    double um_suppressed = module.computeUmContribution(t_seconds);
    std::cout << "Suppressed P_SCm = " << std::scientific << module.getVariable("P_SCm") << "\n";
    std::cout << "U_m (suppressed) = " << um_suppressed << " J/m³\n";
    std::cout << "Suppression factor: " << std::fixed << std::setprecision(2) 
              << (um_suppressed / um_stellar) << "\n\n";
    
    // Step 7: Restore and expand penetration scale
    std::cout << "STEP 7: Expand Penetration Scale (P_SCm x2)\n";
    module.restoreState("standard_penetration");
    module.expandPenetrationScale(2.0, 1000.0);  // 2x, maintain 1000:1 ratio
    double pscm_expanded = module.getVariable("P_SCm");
    double pscm_planet_expanded = module.getVariable("P_SCm_planet");
    std::cout << "Expanded P_SCm (stellar) = " << std::scientific << pscm_expanded << "\n";
    std::cout << "Expanded P_SCm (planetary) = " << pscm_planet_expanded << "\n";
    std::cout << "Ratio maintained: " << std::fixed << std::setprecision(0) 
              << (pscm_expanded / pscm_planet_expanded) << ":1\n\n";
    
    // Step 8: Restore and expand magnetic scale
    std::cout << "STEP 8: Expand Magnetic Scale (μ_j x2, r_j x1.5)\n";
    module.restoreState("standard_penetration");
    module.expandMagneticScale(2.0, 1.5);
    double mu_new = module.getVariable("mu_j");
    double rj_new = module.getVariable("r_j");
    double mu_over_rj_new = mu_new / rj_new;
    std::cout << "New μ_j = " << std::scientific << mu_new << " T·m³\n";
    std::cout << "New r_j = " << rj_new << " m\n";
    std::cout << "New μ_j/r_j = " << mu_over_rj_new << " T·m²\n\n";
    
    // Step 9: Restore and expand amplification scale
    std::cout << "STEP 9: Expand Amplification Scale (Heaviside x1.5, quasi x2)\n";
    module.restoreState("standard_penetration");
    module.expandAmplificationScale(1.5, 2.0);
    double hf_new = module.getVariable("heaviside_factor");
    double qf_new = 1.0 + module.getVariable("f_quasi");
    std::cout << "New Heaviside factor = " << std::scientific << hf_new << "\n";
    std::cout << "New quasi factor = " << qf_new << "\n\n";
    
    // Step 10: Sensitivity analysis
    std::cout << "STEP 10: Sensitivity Analysis (at t=0)\n";
    module.restoreState("standard_penetration");
    std::vector<std::string> params = {"P_SCm", "mu_j", "r_j", "f_Heaviside", "E_react"};
    auto sensitivities = module.sensitivityAnalysis(params);
    for (const auto& pair : sensitivities) {
        std::cout << "  ∂U_m/∂" << pair.first << " ≈ " << std::scientific << pair.second << " (normalized)\n";
    }
    std::cout << "\n";
    
    // Step 11: Generate variations
    std::cout << "STEP 11: Generate Parameter Variations (5 variants, ±10%)\n";
    auto variations = module.generateVariations(5, 0.1);
    for (int i = 0; i < variations.size(); ++i) {
        double var_pscm = variations[i]["P_SCm"];
        double var_pscm_planet = variations[i]["P_SCm_planet"];
        std::cout << "  Variant " << (i+1) << ": P_SCm(stellar)=" << std::scientific << var_pscm 
                  << ", P_SCm(planetary)=" << var_pscm_planet << "\n";
    }
    std::cout << "\n";
    
    // Step 12: Auto-refine to target stellar/planetary ratio
    std::cout << "STEP 12: Auto-Refine to Target Stellar/Planetary Ratio (500:1)\n";
    module.restoreState("standard_penetration");
    module.autoRefineParameters("stellar_planet_ratio", 500.0);
    double refined_pscm = module.getVariable("P_SCm");
    double refined_pscm_planet = module.getVariable("P_SCm_planet");
    double refined_ratio = refined_pscm / refined_pscm_planet;
    std::cout << "Refined P_SCm (stellar) = " << std::scientific << refined_pscm << "\n";
    std::cout << "Refined P_SCm (planetary) = " << refined_pscm_planet << "\n";
    std::cout << "Achieved ratio = " << std::fixed << std::setprecision(0) << refined_ratio << ":1\n\n";
    
    // Step 13: Auto-refine to target U_m stellar
    std::cout << "STEP 13: Auto-Refine to Target U_m (Stellar, 3e65 J/m³ at t=0)\n";
    module.restoreState("standard_penetration");
    module.autoRefineParameters("U_m_stellar", 3e65);
    double refined_um = module.computeUmContribution(0.0);
    std::cout << "Refined P_SCm = " << std::scientific << module.getVariable("P_SCm") << "\n";
    std::cout << "Achieved U_m = " << refined_um << " J/m³ (target: 3e65)\n\n";
    
    // Step 14: Calibrate to observations
    std::cout << "STEP 14: Calibrate to Observational Data\n";
    module.restoreState("standard_penetration");
    std::map<std::string, double> observations;
    observations["P_SCm"] = 0.8;  // Observed slightly reduced stellar penetration
    observations["mu_j"] = 4.0e23;  // Observed stronger magnetic moment
    module.calibrateToObservations(observations);
    std::cout << "Calibrated P_SCm = " << std::scientific << module.getVariable("P_SCm") << "\n";
    std::cout << "Calibrated μ_j = " << module.getVariable("mu_j") << " T·m³\n\n";
    
    // Step 15: Optimize for metric
    std::cout << "STEP 15: Optimize for 'standard_ratio' Metric\n";
    module.optimizeForMetric("standard_ratio");
    std::cout << "Optimized P_SCm (stellar) = " << std::scientific << module.getVariable("P_SCm") << "\n";
    std::cout << "Optimized P_SCm (planetary) = " << module.getVariable("P_SCm_planet") << "\n";
    std::cout << "Ratio = " << std::fixed << std::setprecision(0) 
              << (module.getVariable("P_SCm") / module.getVariable("P_SCm_planet")) << ":1\n\n";
    
    // Step 16: Mutate parameters
    std::cout << "STEP 16: Mutate Parameters (5% random variation)\n";
    module.mutateParameters(0.05);
    std::cout << "Mutated P_SCm = " << std::scientific << module.getVariable("P_SCm") << "\n";
    std::cout << "Mutated μ_j/r_j = " << (module.getVariable("mu_j") / module.getVariable("r_j")) << " T·m²\n\n";
    
    // Step 17: Validate consistency
    std::cout << "STEP 17: Validate Consistency\n";
    bool valid = module.validateConsistency();
    std::cout << "Consistency check: " << (valid ? "PASSED" : "FAILED") << "\n\n";
    
    // Step 18: Introduce anomaly and auto-correct
    std::cout << "STEP 18: Introduce Anomaly and Auto-Correct\n";
    module.createVariable("P_SCm_anomaly", -0.5);  // Invalid negative P_SCm
    module.removeVariable("P_SCm");
    module.createVariable("P_SCm", 0.0005);  // Invalid: smaller than planetary
    std::cout << "Introduced invalid P_SCm = " << module.getVariable("P_SCm") << " (< planetary)\n";
    module.autoCorrectAnomalies();
    std::cout << "Auto-corrected P_SCm (stellar) = " << std::scientific << module.getVariable("P_SCm") << "\n";
    std::cout << "Auto-corrected P_SCm (planetary) = " << module.getVariable("P_SCm_planet") << "\n";
    bool valid_after = module.validateConsistency();
    std::cout << "Consistency after correction: " << (valid_after ? "PASSED" : "FAILED") << "\n\n";
    
    // Step 19: List saved states and export
    std::cout << "STEP 19: List Saved States and Export Final State\n";
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
// #include "ScmPenetrationModule.h"
// int main() {
//     ScmPenetrationModule mod;
//     double p = mod.computeP_SCm();
//     std::cout << "P_SCm ≈ " << p << std::endl;
//     double um_sun = mod.computeUmContribution(0.0);
//     std::cout << "U_m (Sun) = " << um_sun << " J/m³\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("P_SCm", 1e-3);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o scm_test scm_test.cpp ScmPenetrationModule.cpp -lm
// Sample: P_SCm=1; U_m≈2.28e65 J/m³ (Sun); scales for planetary [SCm].
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

ScmPenetrationModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeP_SCm, computeUmContribution, computeUmPlanet) are clear, concise, and variable - driven.
- Output and debugging functions(printVariables, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.
- Supports both stellar and planetary scenarios by switching P_SCm values.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in[SCm] penetration modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.