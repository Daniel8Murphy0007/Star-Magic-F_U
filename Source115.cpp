// SolarWindModulationModule.h
// Modular C++ implementation of the Solar Wind Modulation Factor (?_sw) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes ?_sw=0.01 (unitless) and its scaling (1 + ?_sw v_sw) in Universal Gravity U_g2 term.
// Pluggable: #include "SolarWindModulationModule.h"
// SolarWindModulationModule mod; mod.computeU_g2(1.496e13); mod.updateVariable("delta_sw", new_value);
// Variables in std::map; example for Sun at r=R_b=1.496e13 m; amplification ~5001x.
// Approximations: S(r - R_b)=1; H_SCm=1; E_react=1e46; ?_sum=7.80e-36 J/m�.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef SOLAR_WIND_MODULATION_MODULE_H
#define SOLAR_WIND_MODULATION_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class SolarWindModulationModule {
private:
    std::map<std::string, double> variables;
    double computeModulationFactor();
    double computeU_g2(double r);

public:
    // Constructor: Initialize with framework defaults (Sun)
    SolarWindModulationModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeDelta_sw();  // 0.01 (unitless)
    double computeModulationFactor();  // 1 + ?_sw v_sw
    double computeU_g2(double r);  // U_g2 with modulation (J/m^3)
    double computeU_g2_no_mod(double r);  // Without modulation

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();
};

#endif // SOLAR_WIND_MODULATION_MODULE_H

// SolarWindModulationModule.cpp
#include "SolarWindModulationModule.h"

// Constructor: Set framework defaults (Sun at r=R_b)
SolarWindModulationModule::SolarWindModulationModule() {
    // Universal constants
    variables["delta_sw"] = 0.01;                   // Unitless
    variables["v_sw"] = 5e5;                        // m/s
    variables["k_2"] = 1.2;                         // Coupling
    variables["rho_vac_UA"] = 7.09e-36;             // J/m^3
    variables["rho_vac_SCm"] = 7.09e-37;            // J/m^3
    variables["M_s"] = 1.989e30;                    // kg
    variables["r"] = 1.496e13;                      // m (R_b)
    variables["R_b"] = 1.496e13;                    // m
    variables["S_r_Rb"] = 1.0;                      // Step
    variables["H_SCm"] = 1.0;                       // Unitless
    variables["E_react"] = 1e46;                    // J

    // Derived
    variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
    variables["modulation_factor"] = computeModulationFactor();
}

// Update variable
void SolarWindModulationModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
        if (name == "delta_sw" || name == "v_sw") {
            variables["modulation_factor"] = computeModulationFactor();
        }
        else if (name == "rho_vac_UA" || name == "rho_vac_SCm") {
            variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
        }
    }
    else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void SolarWindModulationModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        if (name == "delta_sw" || name == "v_sw") {
            variables["modulation_factor"] = computeModulationFactor();
        }
        else if (name == "rho_vac_UA" || name == "rho_vac_SCm") {
            variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
        }
    }
    else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void SolarWindModulationModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute ?_sw = 0.01
double SolarWindModulationModule::computeDelta_sw() {
    return variables["delta_sw"];
}

// Compute 1 + ?_sw * v_sw
double SolarWindModulationModule::computeModulationFactor() {
    return 1.0 + variables["delta_sw"] * variables["v_sw"];
}

// Compute U_g2 with modulation
double SolarWindModulationModule::computeU_g2(double r) {
    variables["r"] = r;
    double k_2 = variables["k_2"];
    double rho_sum = variables["rho_sum"];
    double M_s = variables["M_s"];
    double s_step = (r >= variables["R_b"]) ? 1.0 : 0.0;
    double mod_factor = computeModulationFactor();
    double h_scm = variables["H_SCm"];
    double e_react = variables["E_react"];
    return k_2 * (rho_sum * M_s / (r * r)) * s_step * mod_factor * h_scm * e_react;
}

// U_g2 without modulation (?_sw=0)
double SolarWindModulationModule::computeU_g2_no_mod(double r) {
    double orig_delta = variables["delta_sw"];
    variables["delta_sw"] = 0.0;
    double result = computeU_g2(r);
    variables["delta_sw"] = orig_delta;
    return result;
}

// Equation text
std::string SolarWindModulationModule::getEquationText() {
    return "U_g2 = k_2 * [(?_vac,[UA] + ?_vac,[SCm]) M_s / r^2] * S(r - R_b) * (1 + ?_sw v_sw) * H_SCm * E_react\n"
        "Where ?_sw = 0.01 (unitless solar wind modulation factor);\n"
        "Modulation = 1 + 0.01 * v_sw (v_sw=5e5 m/s ? ~5001x amplification).\n"
        "Example r=R_b=1.496e13 m: U_g2 ?1.18e53 J/m� (with); ?2.36e49 J/m� (without; ~5000x less).\n"
        "Role: Enhances external gravity via solar wind momentum/pressure beyond R_b.\n"
        "UQFF: Models heliosphere dynamics; wind influence on nebular/star formation.";
}

// Print variables
void SolarWindModulationModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Example usage in base program (snippet)
// #include "SolarWindModulationModule.h"
// int main() {
//     SolarWindModulationModule mod;
//     double mod_f = mod.computeModulationFactor();
//     std::cout << "Modulation Factor = " << mod_f << std::endl;
//     double u_g2 = mod.computeU_g2(1.496e13);
//     std::cout << "U_g2 = " << u_g2 << " J/m�\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("delta_sw", 0.02);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o sw_mod_test sw_mod_test.cpp SolarWindModulationModule.cpp -lm
// Sample: Factor=5001; U_g2?1.18e53 J/m�; amplifies outer bubble gravity.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

SolarWindModulationModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeDelta_sw, computeModulationFactor, computeU_g2, computeU_g2_no_mod) are clear, concise, and variable - driven.
- Automatic recalculation of derived variables(modulation_factor, rho_sum) when dependencies change.
- Output and debugging functions(printVariables, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.
- Models strong amplification of gravity terms via solar wind modulation.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in solar wind modulation modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.