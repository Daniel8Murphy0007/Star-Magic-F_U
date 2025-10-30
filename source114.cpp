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

// Example usage in base program (snippet)
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