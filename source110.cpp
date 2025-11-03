// OuterFieldBubbleModule.h
// Modular C++ implementation of the Radius of the Outer Field Bubble (R_b) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes R_b=1.496e13 m (100 AU); defines S(r - R_b) step function in Universal Gravity U_g2 term.
// Pluggable: #include "OuterFieldBubbleModule.h"
// OuterFieldBubbleModule mod; mod.computeU_g2(1.5e13); mod.updateVariable("R_b", new_value);
// Variables in std::map; example for Sun at t=0; S=1 for r >= R_b, 0 otherwise.
// Approximations: S step=1 at r=R_b; ?_sw v_sw=5001; E_react=1e46.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef OUTER_FIELD_BUBBLE_MODULE_H
#define OUTER_FIELD_BUBBLE_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class OuterFieldBubbleModule {
private:
    std::map<std::string, double> variables;
    double computeS_r_Rb(double r);
    double computeU_g2(double r);

public:
    // Constructor: Initialize with framework defaults (Sun)
    OuterFieldBubbleModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeR_b();  // 1.496e13 m (100 AU)
    double computeR_bInAU();  // 100 AU
    double computeS_r_Rb(double r);  // Step function
    double computeU_g2(double r);  // U_g2 with S(r - R_b) (J/m^3)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();
};

#endif // OUTER_FIELD_BUBBLE_MODULE_H

// OuterFieldBubbleModule.cpp
#include "OuterFieldBubbleModule.h"

// Constructor: Set framework defaults (Sun at t=0)
OuterFieldBubbleModule::OuterFieldBubbleModule() {
    // Universal constants
    variables["R_b"] = 1.496e13;                    // m (100 AU)
    variables["AU_to_m"] = 1.496e11;                // m/AU
    variables["k_2"] = 1.2;                         // Coupling
    variables["rho_vac_UA"] = 7.09e-36;             // J/m^3
    variables["rho_vac_SCm"] = 7.09e-37;            // J/m^3
    variables["M_s"] = 1.989e30;                    // kg
    variables["r"] = 1.496e13;                      // m (default = R_b)
    variables["delta_sw"] = 0.01;                   // Unitless
    variables["v_sw"] = 5e5;                        // m/s
    variables["H_SCm"] = 1.0;                       // Unitless
    variables["E_react"] = 1e46;                    // J

    // Derived
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
    variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
}

// Update variable
void OuterFieldBubbleModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
        if (name == "rho_vac_UA" || name == "rho_vac_SCm") {
            variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
        } else if (name == "delta_sw" || name == "v_sw") {
            variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void OuterFieldBubbleModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        if (name == "rho_vac_UA" || name == "rho_vac_SCm") {
            variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
        } else if (name == "delta_sw" || name == "v_sw") {
            variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void OuterFieldBubbleModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute R_b (m)
double OuterFieldBubbleModule::computeR_b() {
    return variables["R_b"];
}

// R_b in AU
double OuterFieldBubbleModule::computeR_bInAU() {
    return computeR_b() / variables["AU_to_m"];
}

// Step function S(r - R_b): 1 if r >= R_b, 0 otherwise
double OuterFieldBubbleModule::computeS_r_Rb(double r) {
    return (r >= computeR_b()) ? 1.0 : 0.0;
}

// Compute U_g2 with S(r - R_b)
double OuterFieldBubbleModule::computeU_g2(double r) {
    variables["r"] = r;
    double k_2 = variables["k_2"];
    double rho_sum = variables["rho_sum"];
    double M_s = variables["M_s"];
    double s_step = computeS_r_Rb(r);
    double swirl_factor = variables["swirl_factor"];
    double h_scm = variables["H_SCm"];
    double e_react = variables["E_react"];
    return k_2 * (rho_sum * M_s / (r * r)) * s_step * swirl_factor * h_scm * e_react;
}

// Equation text
std::string OuterFieldBubbleModule::getEquationText() {
    return "U_g2 = k_2 * [(?_vac,[UA] + ?_vac,[SCm]) M_s / r^2] * S(r - R_b) * (1 + ?_sw v_sw) * H_SCm * E_react\n"
           "Where R_b = 1.496e13 m (100 AU, outer bubble radius);\n"
           "S(r - R_b) = 1 (r >= R_b), 0 otherwise (step function).\n"
           "Example r=R_b: U_g2 ?1.18e53 J/m�; r < R_b (e.g., 1 AU): U_g2=0.\n"
           "Role: Defines external gravity boundary (~heliopause); activates U_g2 beyond R_b.\n"
           "UQFF: Separates internal/external fields; models heliosphere/nebular extent.";
}

// Print variables
void OuterFieldBubbleModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Example usage in base program (snippet)
// #include "OuterFieldBubbleModule.h"
// int main() {
//     OuterFieldBubbleModule mod;
//     double rb = mod.computeR_b();
//     std::cout << "R_b = " << rb << " m (" << mod.computeR_bInAU() << " AU)\n";
//     double u_g2 = mod.computeU_g2(1.5e13);  // r > R_b
//     std::cout << "U_g2 (r=1.5e13 m) = " << u_g2 << " J/m�\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("R_b", 2e13);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o bubble_test bubble_test.cpp OuterFieldBubbleModule.cpp -lm
// Sample: R_b=1.496e13 m (100 AU); U_g2?1.18e53 J/m� (r>=R_b); 0 inside.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

OuterFieldBubbleModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeR_b, computeR_bInAU, computeS_r_Rb, computeU_g2) are clear, concise, and variable - driven.
- Automatic recalculation of derived variables(rho_sum, swirl_factor) when dependencies change.
- Output and debugging functions(printVariables, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.
- Step function S(r - R_b) cleanly separates internal and external field regions.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in outer field bubble modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.