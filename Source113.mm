// ScmReactivityDecayModule.h
// Modular C++ implementation of the [SCm] Reactivity Decay Rate (?) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes ?=0.0005 day?� (~5.8e-6 s?�); used in E_react = 10^46 * exp(-? t) for decay in U_m, U_bi, etc.
// Pluggable: #include "ScmReactivityDecayModule.h"
// ScmReactivityDecayModule mod; mod.computeE_react(0.0); mod.updateVariable("kappa_day", new_value);
// Variables in std::map; example for Sun at t=0 (E_react=1e46); t=2000 days: ~3.68e45.
// Approximations: t in days; timescale ~5.5 years; integrates into U_m example.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef SCM_REACTIVITY_DECAY_MODULE_H
#define SCM_REACTIVITY_DECAY_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class ScmReactivityDecayModule {
private:
    std::map<std::string, double> variables;
    double computeKappa_s();  // ? in s?�
    double computeE_react(double t_day);
    double computeUmExample(double t_day);

public:
    // Constructor: Initialize with framework defaults
    ScmReactivityDecayModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeKappa_day();  // 0.0005 day?�
    double computeKappa_s();    // ~5.8e-6 s?�
    double computeE_react(double t_day);  // 1e46 * exp(-? t)
    double computeUmExample(double t_day);  // Simplified U_m with E_react (J/m�)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print decay effects
    void printDecayEffects(double t_day = 2000.0);
};

#endif // SCM_REACTIVITY_DECAY_MODULE_H

// ScmReactivityDecayModule.cpp
#include "ScmReactivityDecayModule.h"

// Constructor: Set framework defaults
ScmReactivityDecayModule::ScmReactivityDecayModule() {
    // Universal constants
    variables["kappa_day"] = 0.0005;                // day?�
    variables["day_to_s"] = 86400.0;                // s/day
    variables["E_react_base"] = 1e46;               // J
    variables["t_day"] = 0.0;                       // days
    variables["mu_over_rj"] = 2.26e10;              // T m� (example)
    variables["P_SCm"] = 1.0;                       // Normalized
    variables["heaviside_f"] = 1e11 + 1.0;          // 1 + 10^13 * 0.01
    variables["quasi_f"] = 1.01;                    // 1 + 0.01
    variables["one_minus_exp"] = 1.0;               // At t=0

    // Derived
    variables["kappa_s"] = computeKappa_s();
}

// Update variable
void ScmReactivityDecayModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
        if (name == "kappa_day") {
            variables["kappa_s"] = computeKappa_s();
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void ScmReactivityDecayModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        if (name == "kappa_day") {
            variables["kappa_s"] = computeKappa_s();
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void ScmReactivityDecayModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute ? (day?�)
double ScmReactivityDecayModule::computeKappa_day() {
    return variables["kappa_day"];
}

// Compute ? in s?�
double ScmReactivityDecayModule::computeKappa_s() {
    return computeKappa_day() / variables["day_to_s"];
}

// Compute E_react = 1e46 * exp(-? t)
double ScmReactivityDecayModule::computeE_react(double t_day) {
    variables["t_day"] = t_day;
    double arg = - computeKappa_day() * t_day;
    return variables["E_react_base"] * std::exp(arg);
}

// Simplified U_m example with E_react
double ScmReactivityDecayModule::computeUmExample(double t_day) {
    double e_react = computeE_react(t_day);
    double one_minus_exp = variables["one_minus_exp"];  // Placeholder; full would compute
    double phi_hat = 1.0;
    double p_scm = variables["P_SCm"];
    double heaviside_f = variables["heaviside_f"];
    double quasi_f = variables["quasi_f"];
    return (variables["mu_over_rj"] * one_minus_exp * phi_hat) * p_scm * e_react * heaviside_f * quasi_f;
}

// Equation text
std::string ScmReactivityDecayModule::getEquationText() {
    return "E_react = 10^46 * exp(-? t) (t days); ?=0.0005 day?� (~5.8e-6 s?�, timescale ~5.5 years).\n"
           "In U_m, U_bi, U_i, U_gi: ... * E_react * ... (decays [SCm] reactivity).\n"
           "Example t=0: E_react=1e46 J; t=2000 days: ~3.68e45 J (~36.8%).\n"
           "U_m (t=0): ?2.28e65 J/m�; t=2000: ?8.39e64 J/m�.\n"
           "Role: Gradual [SCm]-[UA] interaction loss; temporal evolution in jets/nebulae/mergers.\n"
           "UQFF: Models reactivity decay; energy dissipation over cosmic time.";
}

// Print variables
void ScmReactivityDecayModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Print effects
void ScmReactivityDecayModule::printDecayEffects(double t_day) {
    double e_react = computeE_react(t_day);
    double um_ex = computeUmExample(t_day);
    double fraction = e_react / variables["E_react_base"];
    std::cout << "[SCm] Decay Effects at t=" << t_day << " days:\n";
    std::cout << "E_react = " << std::scientific << e_react << " J (" << fraction << " of initial)\n";
    std::cout << "U_m example = " << um_ex << " J/m�\n";
}

// Example usage in base program (snippet)
// #include "ScmReactivityDecayModule.h"
// int main() {
//     ScmReactivityDecayModule mod;
//     double kappa = mod.computeKappa_day();
//     std::cout << "? = " << kappa << " day?�\n";
//     mod.printDecayEffects(2000.0);
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("kappa_day", 0.001);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o scm_decay_test scm_decay_test.cpp ScmReactivityDecayModule.cpp -lm
// Sample: ?=5e-4 day?�; t=2000 days: E_react?3.68e45 J; U_m?8.39e64 J/m�.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

ScmReactivityDecayModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeKappa_day, computeKappa_s, computeE_react, computeUmExample) are clear, concise, and variable - driven.
- Automatic recalculation of derived variables(kappa_s) when dependencies change.
- Output and debugging functions(printVariables, printDecayEffects, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in[SCm] reactivity decay modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.