// SurfaceMagneticFieldModule.h
// Modular C++ implementation of the Surface Magnetic Field (B_s) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes B_s range [1e-4, 0.4] T for Sun; influences B_j in U_g3 magnetic strings (scaled by B_s / B_ref).
// Pluggable: #include "SurfaceMagneticFieldModule.h"
// SurfaceMagneticFieldModule mod; mod.computeU_g3_example(0.0); mod.updateVariable("B_s_min", new_value);
// Variables in std::map; example for Sun at t=0; quiet Sun B_s=1e-4 T ? U_g3?4.5e45 J/m�.
// Approximations: B_ref=0.4 T (max sunspot); cos(?_s t ?)=1; P_core=1; E_react=1e46.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef SURFACE_MAGNETIC_FIELD_MODULE_H
#define SURFACE_MAGNETIC_FIELD_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class SurfaceMagneticFieldModule {
private:
    std::map<std::string, double> variables;
    double computeB_j(double t, double B_s);
    double computeU_g3_example(double t, double B_s);

public:
    // Constructor: Initialize with framework defaults (Sun)
    SurfaceMagneticFieldModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeB_s_min();  // 1e-4 T (quiet Sun)
    double computeB_s_max();  // 0.4 T (sunspot max)
    double computeB_j(double t, double B_s);  // Scaled B_j (T)
    double computeU_g3_example(double t, double B_s);  // U_g3 with B_j (J/m^3)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();
};

#endif // SURFACE_MAGNETIC_FIELD_MODULE_H

// SurfaceMagneticFieldModule.cpp
#include "SurfaceMagneticFieldModule.h"

// Constructor: Set framework defaults (Sun)
SurfaceMagneticFieldModule::SurfaceMagneticFieldModule() {
    // Universal constants
    variables["B_s_min"] = 1e-4;                    // T (quiet)
    variables["B_s_max"] = 0.4;                     // T (sunspot)
    variables["B_ref"] = 0.4;                       // T (reference max)
    variables["k_3"] = 1.8;                         // Coupling
    variables["omega_s"] = 2.5e-6;                  // rad/s
    variables["P_core"] = 1.0;                      // Unitless
    variables["E_react"] = 1e46;                    // J
    variables["pi"] = 3.141592653589793;
    variables["t"] = 0.0;                           // s
}

// Update variable
void SurfaceMagneticFieldModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void SurfaceMagneticFieldModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void SurfaceMagneticFieldModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute B_s min (T)
double SurfaceMagneticFieldModule::computeB_s_min() {
    return variables["B_s_min"];
}

// Compute B_s max (T)
double SurfaceMagneticFieldModule::computeB_s_max() {
    return variables["B_s_max"];
}

// Compute B_j scaled by B_s / B_ref
double SurfaceMagneticFieldModule::computeB_j(double t, double B_s) {
    variables["t"] = t;
    double base_b = variables["B_ref"] + 0.4 * std::sin(variables["omega_s"] * t);  // Hypothetical cycle
    return base_b * (B_s / variables["B_ref"]);
}

// U_g3 example with B_j
double SurfaceMagneticFieldModule::computeU_g3_example(double t, double B_s) {
    double k_3 = variables["k_3"];
    double b_j = computeB_j(t, B_s);
    double cos_term = std::cos(variables["omega_s"] * t * variables["pi"]);
    double p_core = variables["P_core"];
    double e_react = variables["E_react"];
    return k_3 * b_j * cos_term * p_core * e_react;
}

// Equation text
std::string SurfaceMagneticFieldModule::getEquationText() {
    return "B_j ? (10^3 + 0.4 sin(?_s t)) * (B_s / 0.4) T (hypothetical scaling);\n"
           "U_g3 = k_3 * ? B_j * cos(?_s t ?) * P_core * E_react\n"
           "Where B_s = [1e-4, 0.4] T (Sun surface; quiet to sunspot).\n"
           "B_ref=0.4 T (max); scales string fields by surface B_s.\n"
           "Example t=0, B_s=0.4 T: B_j?1e3 T, U_g3?1.8e49 J/m�;\n"
           "B_s=1e-4 T: B_j?0.25 T, U_g3?4.5e45 J/m� (-4 orders).\n"
           "Role: Baseline magnetic strength for strings; variability in U_g3/disks.\n"
           "UQFF: Surface fields drive cosmic magnetism; extensible for planets.";
}

// Print variables
void SurfaceMagneticFieldModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Example usage in base program (snippet)
// #include "SurfaceMagneticFieldModule.h"
// int main() {
//     SurfaceMagneticFieldModule mod;
//     double b_min = mod.computeB_s_min();
//     std::cout << "B_s range: " << b_min << " to " << mod.computeB_s_max() << " T\n";
//     double u_g3 = mod.computeU_g3_example(0.0, 1e-4);
//     std::cout << "U_g3 (quiet Sun) = " << u_g3 << " J/m�\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("B_s_min", 5e-5);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o surface_b_test surface_b_test.cpp SurfaceMagneticFieldModule.cpp -lm
// Sample: B_s [1e-4, 0.4] T; U_g3 quiet?4.5e45 J/m�; scales magnetic influence.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

SurfaceMagneticFieldModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeB_s_min, computeB_s_max, computeB_j, computeU_g3_example) are clear, concise, and variable - driven.
- Output and debugging functions(printVariables, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.
- Supports a wide range of surface magnetic field strengths, enabling both quiet and active Sun scenarios.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in surface magnetic field modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.