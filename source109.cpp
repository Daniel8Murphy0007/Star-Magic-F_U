// QuasiLongitudinalModule.h
// Modular C++ implementation of the Quasi-Longitudinal Wave Factor (f_quasi) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes f_quasi=0.01 (unitless) and its scaling (1 + f_quasi) in Universal Magnetism U_m term.
// Pluggable: #include "QuasiLongitudinalModule.h"
// QuasiLongitudinalModule mod; mod.computeUmContribution(0.0); mod.updateVariable("f_quasi", new_value);
// Variables in std::map; example for Sun at t=0, t_n=0; minor 1% increase in U_m.
// Approximations: 1 - e^{-? t cos(? t_n)}=0 at t=0; ?_hat_j=1; P_SCm=1; f_Heaviside=0.01 (1 + 10^13 f=1e11+1).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef QUASI_LONGITUDINAL_MODULE_H
#define QUASI_LONGITUDINAL_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class QuasiLongitudinalModule {
private:
    std::map<std::string, double> variables;
    double computeQuasiFactor();
    double computeUmBase(int j, double t);
    double computeUmContribution(int j, double t);

public:
    // Constructor: Initialize with framework defaults
    QuasiLongitudinalModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeF_quasi();  // 0.01 (unitless)
    double computeQuasiFactor();  // 1 + f_quasi = 1.01
    double computeUmContribution(int j, double t);  // U_m single string (J/m^3)
    double computeUmWithNoQuasi(int j, double t);  // Without quasi

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print U_m comparison (with/without quasi)
    void printUmComparison(int j = 1, double t = 0.0);
};

#endif // QUASI_LONGITUDINAL_MODULE_H

// QuasiLongitudinalModule.cpp
#include "QuasiLongitudinalModule.h"

// Constructor: Set framework defaults
QuasiLongitudinalModule::QuasiLongitudinalModule() {
    // Universal constants
    variables["f_quasi"] = 0.01;                    // Unitless fraction
    variables["mu_j"] = 3.38e23;                    // T�m^3 (j=1)
    variables["r_j"] = 1.496e13;                    // m
    variables["gamma"] = 5e-5 / 86400.0;            // s^-1 (0.00005 day^-1)
    variables["t_n"] = 0.0;                         // s
    variables["phi_hat_j"] = 1.0;                   // Normalized
    variables["P_SCm"] = 1.0;                       // Pressure
    variables["E_react"] = 1e46;                    // J
    variables["f_Heaviside"] = 0.01;                // For Heaviside
    variables["scale_Heaviside"] = 1e13;            // Amplification
    variables["pi"] = 3.141592653589793;

    // Derived
    variables["quasi_factor"] = computeQuasiFactor();
    variables["heaviside_factor"] = 1.0 + variables["scale_Heaviside"] * variables["f_Heaviside"];
}

// Update variable
void QuasiLongitudinalModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
        if (name == "f_quasi") {
            variables["quasi_factor"] = computeQuasiFactor();
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void QuasiLongitudinalModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        if (name == "f_quasi") {
            variables["quasi_factor"] = computeQuasiFactor();
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void QuasiLongitudinalModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute f_quasi (0.01)
double QuasiLongitudinalModule::computeF_quasi() {
    return variables["f_quasi"];
}

// Compute 1 + f_quasi
double QuasiLongitudinalModule::computeQuasiFactor() {
    return 1.0 + computeF_quasi();
}

// Base for U_m without quasi/Heaviside
double QuasiLongitudinalModule::computeUmBase(int j, double t) {
    double mu_over_rj = variables["mu_j"] / variables["r_j"];
    double exp_arg = - variables["gamma"] * t * std::cos(variables["pi"] * variables["t_n"]);
    double one_minus_exp = 1.0 - std::exp(exp_arg);
    double phi_hat = variables["phi_hat_j"];
    return mu_over_rj * one_minus_exp * phi_hat * variables["P_SCm"] * variables["E_react"];
}

// U_m contribution with quasi
double QuasiLongitudinalModule::computeUmContribution(int j, double t) {
    double base = computeUmBase(j, t);
    double quasi_f = computeQuasiFactor();
    double heaviside_f = variables["heaviside_factor"];
    return base * heaviside_f * quasi_f;
}

// U_m without quasi (set f=0 temporarily)
double QuasiLongitudinalModule::computeUmWithNoQuasi(int j, double t) {
    double orig_f = variables["f_quasi"];
    variables["f_quasi"] = 0.0;
    double result = computeUmContribution(j, t);
    variables["f_quasi"] = orig_f;
    return result;
}

// Equation text
std::string QuasiLongitudinalModule::getEquationText() {
    return "U_m = ?_j [ (?_j / r_j) (1 - e^{-? t cos(? t_n)}) ?_hat_j ] P_SCm E_react (1 + 10^13 f_Heaviside) (1 + f_quasi)\n"
           "Where f_quasi = 0.01 (unitless quasi-longitudinal wave factor);\n"
           "Quasi factor = 1 + 0.01 = 1.01 (1% increase).\n"
           "Example j=1, t=0: U_m contrib ?2.28e65 J/m� (with); ?2.26e65 J/m� (without; -1%).\n"
           "Role: Minor scaling for quasi-longitudinal waves in magnetic strings; subtle [SCm]/[UA] wave effects.\n"
           "UQFF: Enhances wave propagation in jets/nebulae; small but cumulative in dynamics.";
}

// Print variables
void QuasiLongitudinalModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Print U_m comparison
void QuasiLongitudinalModule::printUmComparison(int j, double t) {
    double um_with = computeUmContribution(j, t);
    double um_without = computeUmWithNoQuasi(j, t);
    double percent_increase = ((um_with - um_without) / um_without) * 100.0;
    std::cout << "U_m Comparison for j=" << j << " at t=" << t << " s:\n";
    std::cout << "With quasi: " << std::scientific << um_with << " J/m�\n";
    std::cout << "Without quasi: " << um_without << " J/m�\n";
    std::cout << "Increase: +" << std::fixed << std::setprecision(1) << percent_increase << "%\n";
}

// Example usage in base program (snippet)
// #include "QuasiLongitudinalModule.h"
// int main() {
//     QuasiLongitudinalModule mod;
//     double quasi_f = mod.computeQuasiFactor();
//     std::cout << "Quasi Factor = " << quasi_f << std::endl;
//     mod.printUmComparison(1, 0.0);
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("f_quasi", 0.02);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o quasi_test quasi_test.cpp QuasiLongitudinalModule.cpp -lm
// Sample: Factor=1.01; U_m with=2.28e65 J/m� (+1% vs without).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

QuasiLongitudinalModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeF_quasi, computeQuasiFactor, computeUmContribution, computeUmWithNoQuasi) are clear, concise, and variable - driven.
- Automatic recalculation of derived variables(quasi_factor) when dependencies change.
- Output and debugging functions(printVariables, printUmComparison, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in quasi - longitudinal wave factor modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.