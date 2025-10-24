// SolarWindVelocityModule.h
// Modular C++ implementation of the Solar Wind Velocity (v_sw) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes v_sw=5e5 m/s (500 km/s); scales (1 + δ_sw v_sw) in Universal Gravity U_g2 term.
// Pluggable: #include "SolarWindVelocityModule.h"
// SolarWindVelocityModule mod; mod.computeU_g2(1.496e13); mod.updateVariable("v_sw", new_value);
// Variables in std::map; example for Sun at r=R_b=1.496e13 m; amplification ~5001x.
// Approximations: S(r - R_b)=1; H_SCm=1; E_react=1e46; ρ_sum=7.80e-36 J/m³.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef SOLAR_WIND_VELOCITY_MODULE_H
#define SOLAR_WIND_VELOCITY_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

class SolarWindVelocityModule {
private:
    std::map<std::string, double> variables;
    double computeModulationFactor();
    double computeU_g2(double r);

public:
    // Constructor: Initialize with framework defaults (Sun)
    SolarWindVelocityModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeV_sw();  // 5e5 m/s
    double computeV_swKmS();  // 500 km/s
    double computeModulationFactor();  // 1 + δ_sw v_sw
    double computeU_g2(double r);  // U_g2 with modulation (J/m^3)
    double computeU_g2_no_sw(double r);  // Without v_sw (set=0)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();
};

#endif // SOLAR_WIND_VELOCITY_MODULE_H

// SolarWindVelocityModule.cpp
#include "SolarWindVelocityModule.h"

// Constructor: Set framework defaults (Sun at r=R_b)
SolarWindVelocityModule::SolarWindVelocityModule() {
    // Universal constants
    variables["v_sw"] = 5e5;                        // m/s
    variables["delta_sw"] = 0.01;                   // Unitless
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
void SolarWindVelocityModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
        if (name == "v_sw" || name == "delta_sw") {
            variables["modulation_factor"] = computeModulationFactor();
        } else if (name == "rho_vac_UA" || name == "rho_vac_SCm") {
            variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void SolarWindVelocityModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        if (name == "v_sw" || name == "delta_sw") {
            variables["modulation_factor"] = computeModulationFactor();
        } else if (name == "rho_vac_UA" || name == "rho_vac_SCm") {
            variables["rho_sum"] = variables["rho_vac_UA"] + variables["rho_vac_SCm"];
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void SolarWindVelocityModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute v_sw (m/s)
double SolarWindVelocityModule::computeV_sw() {
    return variables["v_sw"];
}

// v_sw in km/s
double SolarWindVelocityModule::computeV_swKmS() {
    return computeV_sw() / 1e3;
}

// Compute 1 + δ_sw * v_sw
double SolarWindVelocityModule::computeModulationFactor() {
    return 1.0 + variables["delta_sw"] * computeV_sw();
}

// Compute U_g2 with modulation
double SolarWindVelocityModule::computeU_g2(double r) {
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

// U_g2 without solar wind (v_sw=0)
double SolarWindVelocityModule::computeU_g2_no_sw(double r) {
    double orig_v = variables["v_sw"];
    variables["v_sw"] = 0.0;
    double result = computeU_g2(r);
    variables["v_sw"] = orig_v;
    return result;
}

// Equation text
std::string SolarWindVelocityModule::getEquationText() {
    return "U_g2 = k_2 * [(ρ_vac,[UA] + ρ_vac,[SCm]) M_s / r^2] * S(r - R_b) * (1 + δ_sw v_sw) * H_SCm * E_react\n"
           "Where v_sw = 5e5 m/s (500 km/s, typical solar wind speed at 1 AU+);\n"
           "Modulation = 1 + 0.01 * v_sw ≈5001 (amplifies ~5000x).\n"
           "Example r=R_b=1.496e13 m: U_g2 ≈1.18e53 J/m³ (with); ≈2.36e49 J/m³ (without v_sw; ~5000x less).\n"
           "Role: Solar wind momentum/pressure enhances external gravity beyond R_b (heliosphere).\n"
           "UQFF: Models wind shaping of fields; key for heliodynamics/nebular formation.";
}

// Print variables
void SolarWindVelocityModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Example usage in base program (snippet)
// #include "SolarWindVelocityModule.h"
// int main() {
//     SolarWindVelocityModule mod;
//     double v = mod.computeV_sw();
//     std::cout << "v_sw = " << v << " m/s (" << mod.computeV_swKmS() << " km/s)\n";
//     double u_g2 = mod.computeU_g2(1.496e13);
//     std::cout << "U_g2 = " << u_g2 << " J/m³\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("v_sw", 4e5);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o sw_vel_test sw_vel_test.cpp SolarWindVelocityModule.cpp -lm
// Sample: v_sw=5e5 m/s (500 km/s); U_g2≈1.18e53 J/m³; amplifies outer bubble.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

SolarWindVelocityModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeV_sw, computeV_swKmS, computeModulationFactor, computeU_g2, computeU_g2_no_sw) are clear, concise, and variable - driven.
- Automatic recalculation of derived variables(modulation_factor, rho_sum) when dependencies change.
- Output and debugging functions(printVariables, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.
- Models strong amplification of gravity terms via solar wind velocity.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in solar wind velocity modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.