// UgCouplingModule.h
// Modular C++ implementation of the Coupling Constants for Ug Ranges (k_i) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes scaled Universal Gravity terms k_i * U_gi for i=1-4 (Ug1-Ug4), with k1=1.5, k2=1.2, k3=1.8, k4=1.0 (unitless).
// Pluggable: #include "UgCouplingModule.h"
// UgCouplingModule mod; mod.computeSumK_Ugi(); mod.updateVariable("U_g1", new_value);
// Variables in std::map; sum contributes to F_U; placeholders for full U_gi equations.
// Approximations: t_n=0, cos(? t_n)=1; ?_def=0, etc.; example values from Sun at t=0.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef UG_COUPLING_MODULE_H
#define UG_COUPLING_MODULE_H

#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

class UgCouplingModule {
private:
    std::map<std::string, double> variables;
    std::vector<double> k_values;  // [k1, k2, k3, k4]
    std::vector<double> computeAllK_Ugi();

public:
    // Constructor: Initialize with framework defaults
    UgCouplingModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeK_i(int i);  // k_i for specific i (1-4)
    double computeU_gi(int i);  // Placeholder U_gi (J/m^3)
    double computeK_Ugi(int i);  // k_i * U_gi
    std::vector<double> computeAllK_Ugi();  // All four k_i * U_gi
    double computeSumK_Ugi();  // Sum for F_U contribution

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print all k_i * U_gi
    void printK_Ugi();
};

#endif // UG_COUPLING_MODULE_H

// UgCouplingModule.cpp
#include "UgCouplingModule.h"

// Constructor: Set framework defaults
UgCouplingModule::UgCouplingModule() {
    // Coupling constants (unitless)
    k_values = {1.5, 1.2, 1.8, 1.0};               // k1=1.5, k2=1.2, k3=1.8, k4=1.0

    // U_gi defaults (example from Sun at t=0, J/m^3)
    variables["U_g1"] = 1.39e26;                    // Internal Dipole
    variables["U_g2"] = 1.18e53;                    // Outer Field Bubble
    variables["U_g3"] = 1.8e49;                     // Magnetic Strings Disk
    variables["U_g4"] = 2.50e-20;                   // Star-Black Hole Interactions

    // Shared params (placeholders)
    variables["mu_s"] = 1.0;                        // Magnetic moment
    variables["M_s"] = 1.989e30;                    // Stellar mass kg
    variables["r"] = 1e11;                          // m
    variables["alpha"] = 1e-10;                     // Decay rate s^-1
    variables["t_n"] = 0.0;                         // s
    variables["pi"] = 3.141592653589793;
    variables["delta_def"] = 0.0;                   // Deformation
    variables["rho_vac_UA"] = 7.09e-36;             // J/m^3
    variables["rho_vac_SCm"] = 7.09e-37;            // J/m^3
    variables["S_r_Rb"] = 1.0;                      // Step function
    variables["delta_sw"] = 0.0;                    // Swirl deformation
    variables["v_sw"] = 0.0;                        // Solar wind velocity
    variables["H_SCm"] = 1.0;                       // Heaviside SCm
    variables["E_react"] = 1.0;                     // Reactive energy
    variables["M_bh"] = 8.15e36;                    // kg
    variables["d_g"] = 2.55e20;                     // m
    variables["f_feedback"] = 0.0;                  // Feedback factor
}

// Update variable
void UgCouplingModule::updateVariable(const std::string& name, double value) {
    variables[name] = value;
}

// Add delta
void UgCouplingModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void UgCouplingModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute k_i (1-based index)
double UgCouplingModule::computeK_i(int i) {
    if (i < 1 || i > 4) {
        std::cerr << "Invalid i: " << i << ". Using k1." << std::endl;
        return k_values[0];
    }
    return k_values[i-1];
}

// Placeholder compute U_gi (simplified; full eqs require more params)
double UgCouplingModule::computeU_gi(int i) {
    std::string key = "U_g" + std::to_string(i);
    if (variables.find(key) != variables.end()) {
        return variables[key];
    }
    std::cerr << "U_g" << i << " not defined. Returning 0." << std::endl;
    return 0.0;
}

// Compute k_i * U_gi
double UgCouplingModule::computeK_Ugi(int i) {
    return computeK_i(i) * computeU_gi(i);
}

// Compute all k_i * U_gi
std::vector<double> UgCouplingModule::computeAllK_Ugi() {
    std::vector<double> k_ugi(4);
    for (int i = 1; i <= 4; ++i) {
        k_ugi[i-1] = computeK_Ugi(i);
    }
    return k_ugi;
}

// Sum k_i * U_gi for F_U
double UgCouplingModule::computeSumK_Ugi() {
    auto all = computeAllK_Ugi();
    double sum = 0.0;
    for (double val : all) {
        sum += val;
    }
    return sum;
}

// Equation text
std::string UgCouplingModule::getEquationText() {
    return "F_U = ? [k_i * U_gi(r,t,M_s,?_s,T_s,B_s,?_vac,[SCm],?_vac,[UA],t_n) - ?_i * ... ] + other terms\n"
           "k_i (unitless): k1=1.5 (Ug1 Internal Dipole), k2=1.2 (Ug2 Outer Bubble),\n"
           "k3=1.8 (Ug3 Magnetic Disk), k4=1.0 (Ug4 Star-BH).\n"
           "U_g1 = k1 * ?_s ?(M_s/r) e^{-? t} cos(? t_n) (1+?_def);\n"
           "U_g2 = k2 * (?_UA + ?_SCm) M_s / r^2 * S(r-R_b) (1+?_sw v_sw) H_SCm E_react;\n"
           "U_g3 = k3 * (?_SCm + ?_UA) ?_g M_s / d_g * cos(? t_n) (1 + f_mag_str);\n"
           "U_g4 = k4 * ?_SCm M_bh / d_g * e^{-? t} cos(? t_n) (1 + f_feedback).\n"
           "Example Sun t=0: ? k_i U_gi ?1.42e53 J/m� (Ug2 dominant).\n"
           "Role: Scales Ug strengths; normalizes contributions in F_U.";
}

// Print variables
void UgCouplingModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
    std::cout << "k_i values: ";
    for (size_t j = 0; j < k_values.size(); ++j) {
        std::cout << "k" << (j+1) << "=" << k_values[j] << " ";
    }
    std::cout << std::endl;
}

// Print k_i * U_gi
void UgCouplingModule::printK_Ugi() {
    auto all = computeAllK_Ugi();
    std::cout << "Scaled Ug Terms k_i * U_gi (J/m�):\n";
    for (int i = 1; i <= 4; ++i) {
        std::cout << "k" << i << " * U_g" << i << " = " << std::scientific << all[i-1] << std::endl;
    }
    std::cout << "Sum ? k_i U_gi = " << std::scientific << computeSumK_Ugi() << std::endl;
}

// Example usage in base program (snippet)
// #include "UgCouplingModule.h"
// int main() {
//     UgCouplingModule mod;
//     double sum = mod.computeSumK_Ugi();
//     std::cout << "? k_i U_gi = " << sum << " J/m�\n";
//     mod.printK_Ugi();
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("U_g3", 2e49);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o ug_coupling_test ug_coupling_test.cpp UgCouplingModule.cpp -lm
// Sample Output: Sum ?1.42e53 J/m�; k3 amplifies Ug3 by 1.8.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

UgCouplingModule Evaluation

Strengths :
-Modular, extensible design for computing scaled universal gravity terms(k_i * U_gi) in the UQFF framework.
- Clear encapsulation of variables and coupling constants using std::map and std::vector, supporting dynamic updates and easy extension.
- Implements core physical concepts : scaling of Ug terms via k_i, contribution to the unified field(F_U), and separation of physical roles for Ug1 - Ug4.
- Approximations and physical meaning are well - documented in comments and equation text.
- Output functions for variable state and scaled Ug terms support debugging and transparency.
- Handles dynamic updates to variables and recalculates dependent terms as needed.
- Example values and equation text provide context for scientific use and validation.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for flexibility and scalability.
- Some calculations use magic numbers or lack explanatory comments; define named constants and clarify logic.
- Minimal error handling(e.g., invalid index for k_i / U_gi, division by zero); add validation for robustness.
- Unit consistency should be checked and documented for all physical quantities.
- For large - scale or performance - critical simulations, optimize data structures and reduce redundant calculations.
- std::map and std::vector are flexible but may be less efficient than structured types for very large models.
- Expand documentation for function purposes and physical meaning.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in gravity coupling modeling.It implements the UQFF coupling concept faithfully and adapts to various scenarios.For production or high - performance applications, address the recommendations for improved robustness, maintainability, and scalability.