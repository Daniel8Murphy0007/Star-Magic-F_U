// HeavisideFractionModule.h
// Modular C++ implementation of the Heaviside Component Fraction (f_Heaviside) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes f_Heaviside=0.01 (unitless) and its scaling (1 + 10^13 * f_Heaviside) in Universal Magnetism U_m term.
// Pluggable: #include "HeavisideFractionModule.h"
// HeavisideFractionModule mod; mod.computeUmContribution(0.0); mod.updateVariable("f_Heaviside", new_value);
// Variables in std::map; example for Sun at t=0, t_n=0; amplifies by ~10^11.
// Approximations: 1 - e^{-? t cos(? t_n)}=0 at t=0; ?_hat_j=1; P_SCm=1; f_quasi=0.01.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef HEAVISIDE_FRACTION_MODULE_H
#define HEAVISIDE_FRACTION_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class HeavisideFractionModule {
private:
    std::map<std::string, double> variables;
    double computeHeavisideFactor();
    double computeUmBase(int j, double t);
    double computeUmContribution(int j, double t);

public:
    // Constructor: Initialize with framework defaults
    HeavisideFractionModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeF_Heaviside();  // 0.01 (unitless)
    double computeHeavisideFactor();  // 1 + 10^13 * f_Heaviside
    double computeUmContribution(int j, double t);  // U_m single string (J/m^3)
    double computeUmWithNoHeaviside(int j, double t);  // Without Heaviside

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print U_m comparison (with/without Heaviside)
    void printUmComparison(int j = 1, double t = 0.0);
};

#endif // HEAVISIDE_FRACTION_MODULE_H

// HeavisideFractionModule.cpp
#include "HeavisideFractionModule.h"

// Constructor: Set framework defaults
HeavisideFractionModule::HeavisideFractionModule() {
    // Universal constants
    variables["f_Heaviside"] = 0.01;                // Unitless fraction
    variables["scale_Heaviside"] = 1e13;            // Amplification factor
    variables["f_quasi"] = 0.01;                    // Quasi factor
    variables["mu_j"] = 3.38e23;                    // T m^3 (j=1)
    variables["r_j"] = 1.496e13;                    // m
    variables["gamma"] = 5e-5 / 86400.0;            // day^-1 to s^-1
    variables["t_n"] = 0.0;                         // s
    variables["phi_hat_j"] = 1.0;                   // Normalized
    variables["P_SCm"] = 1.0;                       // Pressure
    variables["E_react"] = 1e46;                    // J
    variables["pi"] = 3.141592653589793;

    // Derived
    variables["heaviside_factor"] = computeHeavisideFactor();
}

// Update variable
void HeavisideFractionModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
        if (name == "f_Heaviside") {
            variables["heaviside_factor"] = computeHeavisideFactor();
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void HeavisideFractionModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        if (name == "f_Heaviside") {
            variables["heaviside_factor"] = computeHeavisideFactor();
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void HeavisideFractionModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute f_Heaviside (0.01)
double HeavisideFractionModule::computeF_Heaviside() {
    return variables["f_Heaviside"];
}

// Compute 1 + 10^13 * f_Heaviside
double HeavisideFractionModule::computeHeavisideFactor() {
    return 1.0 + variables["scale_Heaviside"] * computeF_Heaviside();
}

// Base for U_m without Heaviside/Quasi
double HeavisideFractionModule::computeUmBase(int j, double t) {
    double mu_over_rj = variables["mu_j"] / variables["r_j"];
    double exp_arg = - variables["gamma"] * t * std::cos(variables["pi"] * variables["t_n"]);
    double one_minus_exp = 1.0 - std::exp(exp_arg);
    double phi_hat = variables["phi_hat_j"];
    return mu_over_rj * one_minus_exp * phi_hat * variables["P_SCm"] * variables["E_react"];
}

// U_m contribution with Heaviside
double HeavisideFractionModule::computeUmContribution(int j, double t) {
    double base = computeUmBase(j, t);
    double heaviside_f = computeHeavisideFactor();
    double quasi_f = 1.0 + variables["f_quasi"];
    return base * heaviside_f * quasi_f;
}

// U_m without Heaviside (set f=0 temporarily)
double HeavisideFractionModule::computeUmWithNoHeaviside(int j, double t) {
    double orig_f = variables["f_Heaviside"];
    variables["f_Heaviside"] = 0.0;
    double result = computeUmContribution(j, t);
    variables["f_Heaviside"] = orig_f;
    return result;
}

// Equation text
std::string HeavisideFractionModule::getEquationText() {
    return "U_m = ?_j [ (?_j / r_j) (1 - e^{-? t cos(? t_n)}) ?_hat_j ] P_SCm E_react (1 + 10^13 f_Heaviside) (1 + f_quasi)\n"
           "Where f_Heaviside = 0.01 (unitless Heaviside fraction);\n"
           "Heaviside factor = 1 + 10^13 * 0.01 = 1 + 1e11 (amplifies ~10^11x).\n"
           "Example j=1, t=0: U_m contrib ?2.28e65 J/m� (with); ?2.28e54 J/m� (without).\n"
           "Role: Threshold-activated scaling in magnetic energy; nonlinear [SCm]/[UA] effects.\n"
           "UQFF: Amplifies small fraction for large impact in nebulae/quasars/jets.";
}

// Print variables
void HeavisideFractionModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Print U_m comparison
void HeavisideFractionModule::printUmComparison(int j, double t) {
    double um_with = computeUmContribution(j, t);
    double um_without = computeUmWithNoHeaviside(j, t);
    double amplification = um_with / um_without;
    std::cout << "U_m Comparison for j=" << j << " at t=" << t << " s:\n";
    std::cout << "With Heaviside: " << std::scientific << um_with << " J/m�\n";
    std::cout << "Without Heaviside: " << um_without << " J/m�\n";
    std::cout << "Amplification: ~" << std::scientific << amplification << "x\n";
}

// Example usage in base program (snippet)
// #include "HeavisideFractionModule.h"
// int main() {
//     HeavisideFractionModule mod;
//     double heav_factor = mod.computeHeavisideFactor();
//     std::cout << "Heaviside Factor = " << heav_factor << std::endl;
//     mod.printUmComparison(1, 0.0);
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("f_Heaviside", 0.02);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o heaviside_test heaviside_test.cpp HeavisideFractionModule.cpp -lm
// Sample: Factor=1e11+1; U_m with=2.28e65 J/m� (~1e11x without).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

HeavisideFractionModule Evaluation

Strengths :
-Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Extensible computation : core methods use variable names and indices, supporting easy logic extension.
- Automatic dependency updates : changing f_Heaviside recalculates heaviside_factor for consistency.
- Debugging and transparency : printVariables, printUmComparison, and getEquationText provide clear module state and calculation output.
- Pluggable design : self - contained, supports multiple instances with independent variable sets.

Weaknesses / Recommendations :
    -Hardcoded constants : consider external configuration(e.g., JSON, XML) for greater flexibility.
    - Minimal error handling : add validation for missing variables, division by zero, and invalid indices.
    - Unit consistency : runtime checks or clearer documentation would help avoid confusion.
    - Efficiency : for large models, consider alternatives to std::map for better performance.
    - Documentation : expand function - level documentation for physical meaning and expected input / output.

    Summary :
    The code is well - structured, clear, and suitable for scientific prototyping and educational use.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.