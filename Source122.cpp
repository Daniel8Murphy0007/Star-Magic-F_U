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