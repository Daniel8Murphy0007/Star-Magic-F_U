// BuoyancyCouplingModule.h
// Modular C++ implementation of the Buoyancy Coupling Constants (?_i) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes the Universal Buoyancy terms U_bi = -?_i * U_gi * ?_g * (M_bh / d_g) * E_react for i=1 to 4 (Ug1-Ug4).
// Pluggable: #include "BuoyancyCouplingModule.h"
// BuoyancyCouplingModule mod; mod.computeU_bi(1); mod.updateVariable("beta", new_value);
// Variables in std::map; ?_i=0.6 uniform (unitless); opposes gravity with 60% scaling.
// Approximations: cos(? t_n)=1 at t_n=0; ?_sw * ?_vac,sw ?0; U_UA=1; computes per i or sum.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef BUOYANCY_COUPLING_MODULE_H
#define BUOYANCY_COUPLING_MODULE_H

#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

class BuoyancyCouplingModule {
private:
    std::map<std::string, double> variables;
    std::vector<double> computeAllU_bi();

public:
    // Constructor: Initialize with framework defaults
    BuoyancyCouplingModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeBeta(int i);  // ?_i = 0.6 for all i
    double computeU_bi(int i);  // U_bi for specific i (Ug1-4)
    std::vector<double> computeAllU_bi();  // All four U_bi
    double computeF_U_contribution();  // Sum ?_i terms in F_U

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print all U_bi
    void printU_bi();
};

#endif // BUOYANCY_COUPLING_MODULE_H

// BuoyancyCouplingModule.cpp
#include "BuoyancyCouplingModule.h"

// Constructor: Set framework defaults
BuoyancyCouplingModule::BuoyancyCouplingModule() {
    // Universal constants
    variables["beta"] = 0.6;                        // ?_i uniform (unitless)
    variables["Omega_g"] = 7.3e-16;                 // rad/s (galactic spin)
    variables["M_bh"] = 8.15e36;                    // kg (black hole mass)
    variables["d_g"] = 2.55e20;                     // m (galactic distance)
    variables["E_react"] = 1.0;                     // Reactive energy (normalized)
    variables["epsilon_sw"] = 0.001;                // Swirl factor
    variables["rho_vac_sw"] = 8e-21;                // J/m^3
    variables["U_UA"] = 1.0;                        // Universal Aether factor
    variables["t_n"] = 0.0;                         // Time node (s)
    variables["pi"] = 3.141592653589793;

    // U_gi defaults (example from doc for Ug1; others placeholder)
    variables["U_g1"] = 1.39e26;                    // J/m^3
    variables["U_g2"] = 1e25;                       // Placeholder J/m^3
    variables["U_g3"] = 1e24;                       // Placeholder J/m^3
    variables["U_g4"] = 1e23;                       // Placeholder J/m^3
}

// Update variable
void BuoyancyCouplingModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
    // Dependency: If U_gi changes, U_bi updates in compute
}

// Add delta
void BuoyancyCouplingModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void BuoyancyCouplingModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute ?_i (uniform 0.6)
double BuoyancyCouplingModule::computeBeta(int i) {
    return variables["beta"];  // Uniform for all i
}

// Compute U_bi for specific i
double BuoyancyCouplingModule::computeU_bi(int i) {
    std::string ug_key = "U_g" + std::to_string(i);
    if (variables.find(ug_key) == variables.end()) {
        std::cerr << "U_g" << i << " not found. Using 1e26 default." << std::endl;
        return -0.6 * 1e26 * variables["Omega_g"] * (variables["M_bh"] / variables["d_g"]) * variables["E_react"];
    }
    double U_gi = variables[ug_key];
    double beta_i = computeBeta(i);
    double M_bh_over_d_g = variables["M_bh"] / variables["d_g"];
    double swirl_factor = 1.0 + variables["epsilon_sw"] * variables["rho_vac_sw"];
    double cos_term = std::cos(variables["pi"] * variables["t_n"]);
    return -beta_i * U_gi * variables["Omega_g"] * M_bh_over_d_g * variables["E_react"] * swirl_factor * variables["U_UA"] * cos_term;
}

// Compute all U_bi
std::vector<double> BuoyancyCouplingModule::computeAllU_bi() {
    std::vector<double> u_bi(4);
    for (int i = 1; i <= 4; ++i) {
        u_bi[i-1] = computeU_bi(i);
    }
    return u_bi;
}

// Contribution to F_U (sum of -?_i terms)
double BuoyancyCouplingModule::computeF_U_contribution() {
    auto all_u_bi = computeAllU_bi();
    double sum = 0.0;
    for (double val : all_u_bi) {
        sum += val;
    }
    return sum;
}

// Equation text
std::string BuoyancyCouplingModule::getEquationText() {
    return "U_bi = -?_i * U_gi * ?_g * (M_bh / d_g) * E_react * (1 + ?_sw * ?_vac,sw) * U_UA * cos(? t_n)\n"
           "Where ?_i = 0.6 (unitless, uniform for i=1-4: Ug1-Ug4);\n"
           "Opposes gravity: 60% scaling of gravitational term.\n"
           "In F_U: ? [k_i U_gi - ?_i U_gi ?_g (M_bh/d_g) E_react] + other terms.\n"
           "Role: Stabilizes systems (e.g., molecular clouds, nebulae); counteracts Ug collapse.\n"
           "Example Ug1: U_b1 ? -1.94e27 J/m� (at t_n=0, Sun params)."
           "UQFF: Uniform buoyancy across scales; tunable for refinements.";
}

// Print variables
void BuoyancyCouplingModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Print U_bi
void BuoyancyCouplingModule::printU_bi() {
    auto all_u_bi = computeAllU_bi();
    std::cout << "Universal Buoyancy Terms U_bi (J/m�):\n";
    for (int i = 1; i <= 4; ++i) {
        std::cout << "U_b" << i << " = " << std::scientific << all_u_bi[i-1] << std::endl;
    }
    std::cout << "F_U Buoyancy Contribution (sum): " << std::scientific << computeF_U_contribution() << std::endl;
}

// Example usage in base program (snippet)
// #include "BuoyancyCouplingModule.h"
// int main() {
//     BuoyancyCouplingModule mod;
//     double u_b1 = mod.computeU_bi(1);
//     std::cout << "U_b1 = " << u_b1 << " J/m�\n";
//     mod.printU_bi();
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("beta", 0.7);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o buoyancy_test buoyancy_test.cpp BuoyancyCouplingModule.cpp -lm
// Sample Output: U_b1 ? -1.94e27 J/m�; sum opposes gravity by ~60% scaled.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

BuoyancyCouplingModule Evaluation

Strengths :
-Modular, extensible design for computing buoyancy coupling constants(?_i) and universal buoyancy terms(U_bi) in the UQFF framework.
- Clear encapsulation of variables and buoyancy terms using std::map and std::vector, supporting dynamic updates and easy extension.
- Implements core physical concepts : ?_i scaling, opposition to gravity, and contributions to the unified field(F_U).
- Approximations and physical meaning are well - documented in comments and equation text.
- Output functions for variable state and U_bi terms support debugging and transparency.
- Uniform ?_i(0.6) simplifies analysis and tuning; supports per - term and summed contributions.
- Handles dynamic updates to variables and recalculates dependent terms as needed.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for flexibility and scalability.
- Some calculations use magic numbers or lack explanatory comments; define named constants and clarify logic.
- Minimal error handling(e.g., division by zero, invalid variable names); add validation for robustness.
- Unit consistency should be checked and documented for all physical quantities.
- For large - scale or performance - critical simulations, optimize data structures and reduce redundant calculations.
- std::map and std::vector are flexible but may be less efficient than structured types for very large models.
- Expand documentation for function purposes and physical meaning.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in buoyancy coupling modeling.It implements the UQFF buoyancy concept faithfully and adapts to various scenarios.For production or high - performance applications, address the recommendations for improved robustness, maintainability, and scalability.