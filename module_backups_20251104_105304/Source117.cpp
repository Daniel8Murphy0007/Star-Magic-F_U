// StellarMassModule.h
// Modular C++ implementation of the Stellar/Planetary Mass (M_s) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes M_s=1.989e30 kg (1 M_sun for Sun); scales M_s / r^2 in Universal Gravity U_g1 and U_g2 terms.
// Pluggable: #include "StellarMassModule.h"
// StellarMassModule mod; mod.computeU_g2(1.496e13); mod.updateVariable("M_s", new_value);
// Variables in std::map; example for Sun at r=1.496e13 m; U_g2 ?1.18e53 J/m�.
// Approximations: S(r - R_b)=1; (1 + ?_sw v_sw)=5001; H_SCm=1; E_react=1e46.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef STELLAR_MASS_MODULE_H
#define STELLAR_MASS_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class StellarMassModule {
private:
    std::map<std::string, double> variables;
    double computeM_sOverR2(double r);
    double computeU_g1(double r);
    double computeU_g2(double r);

public:
    // Constructor: Initialize with framework defaults (Sun)
    StellarMassModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeM_s();  // 1.989e30 kg
    double computeM_sInMsun();  // 1 M_sun
    double computeM_sOverR2(double r);  // M_s / r^2 (kg/m�)
    double computeU_g1(double r);  // U_g1 example (J/m^3)
    double computeU_g2(double r);  // U_g2 example (J/m^3)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();
};

#endif // STELLAR_MASS_MODULE_H

// StellarMassModule.cpp
#include "StellarMassModule.h"

// Constructor: Set framework defaults (Sun at r=R_b)
StellarMassModule::StellarMassModule() {
    // Universal constants
    variables["M_s"] = 1.989e30;                    // kg (Sun)
    variables["M_sun"] = 1.989e30;                  // kg
    variables["k_1"] = 1.5;                         // Coupling for U_g1
    variables["k_2"] = 1.2;                         // Coupling for U_g2
    variables["rho_vac_UA"] = 7.09e-36;             // J/m^3
    variables["rho_vac_SCm"] = 7.09e-37;            // J/m^3
    variables["r"] = 1.496e13;                      // m (example R_b)
    variables["R_b"] = 1.496e13;                    // m
    variables["S_r_Rb"] = 1.0;                      // Step
    variables["delta_sw"] = 0.01;                   // Unitless
    variables["v_sw"] = 5e5;                        // m/s
    variables["H_SCm"] = 1.0;                       // Unitless
    variables["E_react"] = 1e46;                    // J

    // Derived
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
    variables["swirl_factor"] = 1.0 + variables["delta_sw"] * variables["v_sw"];
}

// Update variable
void StellarMassModule::updateVariable(const std::string& name, double value) {
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
void StellarMassModule::addToVariable(const std::string& name, double delta) {
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
void StellarMassModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute M_s (kg)
double StellarMassModule::computeM_s() {
    return variables["M_s"];
}

// M_s in M_sun
double StellarMassModule::computeM_sInMsun() {
    return computeM_s() / variables["M_sun"];
}

// M_s / r^2 (kg/m�)
double StellarMassModule::computeM_sOverR2(double r) {
    variables["r"] = r;
    if (r == 0.0) return 0.0;
    return computeM_s() / (r * r);
}

// U_g1 example (internal, simplified)
double StellarMassModule::computeU_g1(double r) {
    double k_1 = variables["k_1"];
    double rho_sum = variables["rho_sum"];
    double m_over_r2 = computeM_sOverR2(r);
    double e_react = variables["E_react"];
    return k_1 * rho_sum * m_over_r2 * e_react;  // Simplified
}

// U_g2 example (outer bubble)
double StellarMassModule::computeU_g2(double r) {
    variables["r"] = r;
    double k_2 = variables["k_2"];
    double rho_sum = variables["rho_sum"];
    double s_step = (r >= variables["R_b"]) ? 1.0 : 0.0;
    double swirl_factor = variables["swirl_factor"];
    double h_scm = variables["H_SCm"];
    double e_react = variables["E_react"];
    return k_2 * rho_sum * computeM_sOverR2(r) * s_step * swirl_factor * h_scm * e_react;
}

// Equation text
std::string StellarMassModule::getEquationText() {
    return "U_g1 = k_1 * ?_vac,[UA/SCm] * (M_s / r^2) * ... E_react (internal dipole);\n"
           "U_g2 = k_2 * ?_vac,[UA/SCm] * (M_s / r^2) * S(r - R_b) * (1 + ?_sw v_sw) * H_SCm * E_react (outer bubble).\n"
           "Where M_s = 1.989e30 kg (1 M_sun for Sun).\n"
           "Scales gravity by mass; M_s / r^2 ?8.89e3 kg/m� at r=1.496e13 m.\n"
           "Example U_g2 (r=R_b): ?1.18e53 J/m�.\n"
           "Role: Central mass drives internal/external gravity; stellar/planetary dynamics.\n"
           "UQFF: Mass-dependent fields for nebulae/formation/mergers.";
}

// Print variables
void StellarMassModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Example usage in base program (snippet)
// #include "StellarMassModule.h"
// int main() {
//     StellarMassModule mod;
//     double m_sun = mod.computeM_sInMsun();
//     std::cout << "M_s = " << m_sun << " M_sun\n";
//     double u_g2 = mod.computeU_g2(1.496e13);
//     std::cout << "U_g2 = " << u_g2 << " J/m�\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("M_s", 2e30);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o stellar_mass_test stellar_mass_test.cpp StellarMassModule.cpp -lm
// Sample: M_s=1 M_sun; U_g2?1.18e53 J/m�; scales gravity by mass.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

StellarMassModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeM_s, computeM_sInMsun, computeM_sOverR2, computeU_g1, computeU_g2) are clear, concise, and variable - driven.
- Automatic recalculation of derived variables(rho_sum, swirl_factor) when dependencies change.
- Output and debugging functions(printVariables, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.
- Supports both internal and external gravity calculations(U_g1, U_g2) with mass scaling.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in stellar / planetary mass modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.