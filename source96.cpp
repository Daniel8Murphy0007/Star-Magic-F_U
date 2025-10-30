// GalacticDistanceModule.h
// Modular C++ implementation of the Distance from Galactic Center (d_g) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes d_g=2.55e20 m (~27,000 ly) and conversions; scales M_bh / d_g in U_bi and Ug4.
// Pluggable: #include "GalacticDistanceModule.h"
// GalacticDistanceModule mod; mod.computeMbhOverDg(); mod.updateVariable("d_g", new_value);
// Variables in std::map; example for Sun at t=0, t_n=0.
// Approximations: cos(π t_n)=1; ε_sw * ρ_vac,sw ≈0; α=0.001 s^-1; f_feedback=0.1.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef GALACTIC_DISTANCE_MODULE_H
#define GALACTIC_DISTANCE_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class GalacticDistanceModule {
private:
    std::map<std::string, double> variables;
    double computeDgInLy();
    double computeDgInPc();
    double computeMbhOverDg();
    double computeU_b1();
    double computeU_g4();

public:
    // Constructor: Initialize with framework defaults
    GalacticDistanceModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeDg();  // d_g in m (2.55e20)
    double computeDgInLy();
    double computeDgInPc();
    double computeMbhOverDg();  // M_bh / d_g (kg/m)
    double computeU_b1();  // Universal Buoyancy example (J/m^3)
    double computeU_g4();  // Ug4 example (J/m^3)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();
};

#endif // GALACTIC_DISTANCE_MODULE_H

// GalacticDistanceModule.cpp
#include "GalacticDistanceModule.h"

// Constructor: Set framework defaults
GalacticDistanceModule::GalacticDistanceModule() {
    // Universal constants
    variables["c"] = 2.998e8;                       // m/s
    variables["year_to_s"] = 3.156e7;               // s/yr
    variables["ly_to_m"] = variables["c"] * variables["year_to_s"];  // ≈9.461e15 m/ly
    variables["pc_to_ly"] = 3.262;                  // ly/pc
    variables["pi"] = 3.141592653589793;

    // Galactic params
    variables["d_g"] = 2.55e20;                     // m (~27,000 ly)
    variables["M_bh"] = 8.15e36;                    // kg (Sgr A* mass)

    // U_bi params
    variables["beta_1"] = 0.6;                      // Unitless
    variables["U_g1"] = 1.39e26;                    // J/m^3
    variables["Omega_g"] = 7.3e-16;                 // rad/s
    variables["epsilon_sw"] = 0.001;                // Unitless
    variables["rho_vac_sw"] = 8e-21;                // J/m^3
    variables["U_UA"] = 1.0;                        // Normalized
    variables["t_n"] = 0.0;                         // s

    // Ug4 params
    variables["k_4"] = 1.0;                         // Unitless
    variables["rho_vac_SCm"] = 7.09e-37;            // J/m^3
    variables["alpha"] = 0.001;                     // s^-1
    variables["f_feedback"] = 0.1;                  // Unitless
}

// Update variable
void GalacticDistanceModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void GalacticDistanceModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void GalacticDistanceModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute d_g (m)
double GalacticDistanceModule::computeDg() {
    return variables["d_g"];
}

// d_g in ly
double GalacticDistanceModule::computeDgInLy() {
    return computeDg() / variables["ly_to_m"];
}

// d_g in pc
double GalacticDistanceModule::computeDgInPc() {
    return computeDgInLy() / variables["pc_to_ly"];
}

// M_bh / d_g (kg/m)
double GalacticDistanceModule::computeMbhOverDg() {
    return variables["M_bh"] / computeDg();
}

// U_b1 example (J/m^3)
double GalacticDistanceModule::computeU_b1() {
    double beta_1 = variables["beta_1"];
    double U_g1 = variables["U_g1"];
    double Omega_g = variables["Omega_g"];
    double mbh_over_dg = computeMbhOverDg();
    double swirl_factor = 1.0 + variables["epsilon_sw"] * variables["rho_vac_sw"];
    double U_UA = variables["U_UA"];
    double cos_term = std::cos(variables["pi"] * variables["t_n"]);
    return -beta_1 * U_g1 * Omega_g * mbh_over_dg * swirl_factor * U_UA * cos_term;
}

// U_g4 example (J/m^3)
double GalacticDistanceModule::computeU_g4() {
    double k_4 = variables["k_4"];
    double rho_vac_SCm = variables["rho_vac_SCm"];
    double mbh_over_dg = computeMbhOverDg();
    double exp_term = std::exp( - variables["alpha"] * variables["t_n"] );
    double cos_term = std::cos(variables["pi"] * variables["t_n"]);
    double feedback_factor = 1.0 + variables["f_feedback"];
    return k_4 * (rho_vac_SCm * variables["M_bh"]) / computeDg() * exp_term * cos_term * feedback_factor;
}

// Equation text
std::string GalacticDistanceModule::getEquationText() {
    return "U_bi = -β_i U_gi Ω_g (M_bh / d_g) (1 + ε_sw ρ_vac,sw) U_UA cos(π t_n)\n"
           "U_g4 = k_4 (ρ_vac,[SCm] M_bh / d_g) e^{-α t} cos(π t_n) (1 + f_feedback)\n"
           "Where d_g = 2.55e20 m (~27,000 ly / 8,260 pc; Sun to Sgr A*).\n"
           "M_bh / d_g ≈3.20e16 kg/m;\n"
           "Example U_b1 ≈ -1.94e27 J/m³; U_g4 ≈2.50e-20 J/m³ (t_n=0).\n"
           "Role: Scales SMBH influence on buoyancy/Ug4; galactic dynamics in nebulae/disks.\n"
           "UQFF: Enables merger resolution (final parsec); star formation modulation.";
}

// Print variables
void GalacticDistanceModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Example usage in base program (snippet)
// #include "GalacticDistanceModule.h"
// int main() {
//     GalacticDistanceModule mod;
//     double dg_ly = mod.computeDgInLy();
//     std::cout << "d_g ≈ " << dg_ly << " ly\n";
//     double u_b1 = mod.computeU_b1();
//     std::cout << "U_b1 = " << u_b1 << " J/m³\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("d_g", 2.6e20);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o galactic_test galactic_test.cpp GalacticDistanceModule.cpp -lm
// Sample: d_g ≈2.70e4 ly; U_b1 ≈ -1.94e27 J/m³; scales SMBH at galactic distances.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

GalacticDistanceModule Evaluation

Strengths :
-Modular, extensible design for modeling galactic center distance and its role in UQFF calculations.
- Clear encapsulation of variables using std::map, supporting dynamic updates and easy extension.
- Implements core physical concepts : conversion of d_g between units(m, ly, pc), SMBH scaling(M_bh / d_g), and contributions to universal buoyancy(U_b1) and gravity(Ug4).
- Approximations and physical meaning are well - documented in comments and equation text.
- Output functions for variable state and equation text support debugging and transparency.
- Handles dynamic updates to variables and recalculates dependent terms as needed.
- Example calculations and conversion functions provide scientific context and validation.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for flexibility and scalability.
- Some calculations use magic numbers or lack explanatory comments; define named constants and clarify logic.
- Minimal error handling(e.g., division by zero, invalid variable names); add validation for robustness.
- Unit consistency should be checked and documented for all physical quantities.
- For large - scale or performance - critical simulations, optimize data structures and reduce redundant calculations.
- std::map is flexible but may be less efficient than structured types for very large models.
- Expand documentation for function purposes and physical meaning.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in galactic distance modeling.It implements the UQFF distance concept faithfully and adapts to various scenarios.For production or high - performance applications, address the recommendations for improved robustness, maintainability, and scalability.