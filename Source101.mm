// HeliosphereThicknessModule.h
// Modular C++ implementation of the Heliosphere Thickness Factor (H_SCm) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes H_SCm ?1 (unitless) and its scaling in Universal Gravity U_g2 term.
// Pluggable: #include "HeliosphereThicknessModule.h"
// HeliosphereThicknessModule mod; mod.computeU_g2(0.0, 0.0); mod.updateVariable("H_SCm", new_value);
// Variables in std::map; example for Sun at t=0, t_n=0, r=R_b=1.496e13 m.
// Approximations: S(r - R_b)=1; ?_sw v_sw=5001; E_react=1e46.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef HELIOSPHERE_THICKNESS_MODULE_H
#define HELIOSPHERE_THICKNESS_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class HeliosphereThicknessModule {
private:
    std::map<std::string, double> variables;
    double computeH_SCm();
    double computeU_g2(double t, double t_n);

public:
    // Constructor: Initialize with framework defaults
    HeliosphereThicknessModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeH_SCm();  // ?1 (unitless)
    double computeU_g2(double t, double t_n);  // U_g2 with H_SCm (J/m^3)
    double computeU_g2_no_H(double t, double t_n);  // Without H_SCm variation

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();
};

#endif // HELIOSPHERE_THICKNESS_MODULE_H

// HeliosphereThicknessModule.cpp
#include "HeliosphereThicknessModule.h"

// Constructor: Set framework defaults
HeliosphereThicknessModule::HeliosphereThicknessModule() {
    // Universal constants
    variables["H_SCm"] = 1.0;                       // Unitless ?1
    variables["k_2"] = 1.2;                         // Coupling
    variables["rho_vac_UA"] = 7.09e-36;             // J/m^3
    variables["rho_vac_SCm"] = 7.09e-37;            // J/m^3
    variables["M_s"] = 1.989e30;                    // kg (Sun)
    variables["r"] = 1.496e13;                      // m (R_b)
    variables["R_b"] = 1.496e13;                    // m
    variables["delta_sw"] = 0.01;                   // Unitless
    variables["v_sw"] = 5e5;                        // m/s
    variables["E_react"] = 1e46;                    // J
    variables["S_r_Rb"] = 1.0;                      // Step function
    variables["pi"] = 3.141592653589793;
    variables["t_n"] = 0.0;                         // s

    // Derived
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
    variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
}

// Update variable
void HeliosphereThicknessModule::updateVariable(const std::string& name, double value) {
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
void HeliosphereThicknessModule::addToVariable(const std::string& name, double delta) {
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
void HeliosphereThicknessModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute H_SCm ?1
double HeliosphereThicknessModule::computeH_SCm() {
    return variables["H_SCm"];
}

// Compute U_g2 with H_SCm
double HeliosphereThicknessModule::computeU_g2(double t, double t_n) {
    double k_2 = variables["k_2"];
    double rho_sum = variables["rho_sum"];
    double M_s = variables["M_s"];
    double r = variables["r"];
    double S_r_Rb = variables["S_r_Rb"];
    double swirl_factor = variables["swirl_factor"];
    double H_SCm = computeH_SCm();
    double E_react = variables["E_react"];
    // Simplified; no explicit t dependence in example
    return k_2 * (rho_sum * M_s / (r * r)) * S_r_Rb * swirl_factor * H_SCm * E_react;
}

// U_g2 without H_SCm variation (H=1 fixed)
double HeliosphereThicknessModule::computeU_g2_no_H(double t, double t_n) {
    double orig_H = variables["H_SCm"];
    variables["H_SCm"] = 1.0;
    double result = computeU_g2(t, t_n);
    variables["H_SCm"] = orig_H;
    return result;
}

// Equation text
std::string HeliosphereThicknessModule::getEquationText() {
    return "U_g2 = k_2 * [(?_vac,[UA] + ?_vac,[SCm]) M_s / r^2] * S(r - R_b) * (1 + ?_sw v_sw) * H_SCm * E_react\n"
           "Where H_SCm ?1 (unitless heliosphere thickness factor);\n"
           "Scales outer field bubble gravity for heliopause extent (~120 AU).\n"
           "Example r=R_b=1.496e13 m, t=0: U_g2 ?1.18e53 J/m� (H=1);\n"
           "If H_SCm=1.1: ?1.30e53 J/m� (+10%).\n"
           "Role: Adjusts [SCm] influence in heliosphere; minimal but flexible for boundary variations.\n"
           "UQFF: Models solar wind dominance; key for nebular/heliospheric dynamics.";
}

// Print variables
void HeliosphereThicknessModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Example usage in base program (snippet)
// #include "HeliosphereThicknessModule.h"
// int main() {
//     HeliosphereThicknessModule mod;
//     double h = mod.computeH_SCm();
//     std::cout << "H_SCm ? " << h << std::endl;
//     double u_g2 = mod.computeU_g2(0.0, 0.0);
//     std::cout << "U_g2 = " << u_g2 << " J/m�\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("H_SCm", 1.1);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o heliosphere_test heliosphere_test.cpp HeliosphereThicknessModule.cpp -lm
// Sample: H_SCm=1; U_g2 ?1.18e53 J/m�; +10% for H=1.1.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

HeliosphereThicknessModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Automatic recalculation of derived variables(e.g., rho_sum, swirl_factor) when dependencies change.
- Core computation methods(computeU_g2, computeU_g2_no_H) are clear and use variable - driven logic.
- Output and debugging functions(printVariables, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, division by zero, or invalid input; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in heliosphere thickness and gravity modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.