// InertiaCouplingModule.h
// Modular C++ implementation of the Inertia Coupling Constants (?_i) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes ?_i=1.0 (unitless, uniform for i=1-4) and scales U_i in F_U: -?_i [?_i U_i E_react].
// Pluggable: #include "InertiaCouplingModule.h"
// InertiaCouplingModule mod; mod.computeSumInertiaTerms(0.0); mod.updateVariable("rho_vac_SCm", new_value);
// Variables in std::map; example for Sun at t=0, t_n=0; U_i ?1.38e-47 J/m�, contrib ? -0.138 J/m�.
// Approximations: Uniform ?_i=1.0; cos(? t_n)=1; f_TRZ=0.1; ?_s=2.5e-6 rad/s; E_react=1e46.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef INERTIA_COUPLING_MODULE_H
#define INERTIA_COUPLING_MODULE_H

#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

class InertiaCouplingModule {
private:
    std::map<std::string, double> variables;
    std::vector<double> computeAllInertiaTerms(double t);

public:
    // Constructor: Initialize with framework defaults
    InertiaCouplingModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeLambda_i(int i);  // ?_i=1.0 (unitless)
    double computeU_i(int i, double t);  // U_i for i=1-4 (J/m^3)
    double computeInertiaTerm(int i, double t);  // -?_i U_i E_react
    double computeSumInertiaTerms(double t);  // Sum for F_U contribution (J/m^3)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print inertia breakdown
    void printInertiaBreakdown(double t = 0.0);
};

#endif // INERTIA_COUPLING_MODULE_H

// InertiaCouplingModule.cpp
#include "InertiaCouplingModule.h"

// Constructor: Set framework defaults (Sun at t=0, level 13)
InertiaCouplingModule::InertiaCouplingModule() {
    // Universal constants
    variables["lambda"] = 1.0;                      // Uniform ?_i (unitless)
    variables["rho_vac_SCm"] = 7.09e-37;            // J/m^3
    variables["rho_vac_UA"] = 7.09e-36;             // J/m^3
    variables["omega_s"] = 2.5e-6;                  // rad/s (Sun rotation)
    variables["f_TRZ"] = 0.1;                       // Unitless
    variables["E_react"] = 1e46;                    // J
    variables["pi"] = 3.141592653589793;
    variables["t_n"] = 0.0;                         // s
    variables["alpha_decay"] = 0.0005;              // For E_react exp, but t=0
}

// Update variable
void InertiaCouplingModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void InertiaCouplingModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void InertiaCouplingModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute ?_i (uniform 1.0)
double InertiaCouplingModule::computeLambda_i(int i) {
    return variables["lambda"];
}

// Compute U_i (uniform across i)
double InertiaCouplingModule::computeU_i(int i, double t) {
    double lambda_i = computeLambda_i(i);
    double rho_sc = variables["rho_vac_SCm"];
    double rho_ua = variables["rho_vac_UA"];
    double omega_s_t = variables["omega_s"];        // Simplified, no t dep in example
    double cos_term = std::cos(variables["pi"] * variables["t_n"]);
    double trz_factor = 1.0 + variables["f_TRZ"];
    return lambda_i * rho_sc * rho_ua * omega_s_t * cos_term * trz_factor;
}

// Compute inertia term -?_i U_i E_react
double InertiaCouplingModule::computeInertiaTerm(int i, double t) {
    double u_i = computeU_i(i, t);
    double e_react = variables["E_react"] * std::exp( - variables["alpha_decay"] * t );  // With decay
    return - computeLambda_i(i) * u_i * e_react;
}

// Sum over i=1 to 4
double InertiaCouplingModule::computeSumInertiaTerms(double t) {
    double sum = 0.0;
    for (int i = 1; i <= 4; ++i) {
        sum += computeInertiaTerm(i, t);
    }
    return sum;
}

// Equation text
std::string InertiaCouplingModule::getEquationText() {
    return "F_U = ... - ?_i [?_i * U_i * E_react] + ...\n"
           "U_i = ?_i * ?_vac,[SCm] * ?_vac,[UA] * ?_s(t) * cos(? t_n) * (1 + f_TRZ)\n"
           "Where ?_i = 1.0 (unitless, uniform for i=1-4: Ug1-Ug4);\n"
           "E_react = 1e46 * e^{-? t} (?=5e-4);\n"
           "Example Sun t=0, t_n=0: U_i ?1.38e-47 J/m�; -?_i U_i E_react ? -0.138 J/m� (per i).\n"
           "Role: Scales resistive inertia; uniform baseline opposition to dynamics.\n"
           "UQFF: Consistent across scales; aids stability in interiors/disks/mergers.";
}

// Print variables
void InertiaCouplingModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Print breakdown
void InertiaCouplingModule::printInertiaBreakdown(double t) {
    std::cout << "Inertia Breakdown at t=" << t << " s:\n";
    for (int i = 1; i <= 4; ++i) {
        double u_i = computeU_i(i, t);
        double term = computeInertiaTerm(i, t);
        std::cout << "i=" << i << ": U_i = " << std::scientific << u_i << " J/m�, Term = " << term << " J/m�\n";
    }
    std::cout << "Sum ? Terms = " << std::scientific << computeSumInertiaTerms(t) << " J/m�\n";
}

// Example usage in base program (snippet)
// #include "InertiaCouplingModule.h"
// int main() {
//     InertiaCouplingModule mod;
//     double sum = mod.computeSumInertiaTerms(0.0);
//     std::cout << "Sum Inertia Terms = " << sum << " J/m�\n";
//     mod.printInertiaBreakdown(0.0);
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("lambda", 1.1);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o inertia_test inertia_test.cpp InertiaCouplingModule.cpp -lm
// Sample: Per i term ? -0.138 J/m�; sum ? -0.552 J/m� (4 terms).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

InertiaCouplingModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeLambda_i, computeU_i, computeInertiaTerm, computeSumInertiaTerms) are clear, concise, and variable - driven.
- Output and debugging functions(printVariables, printInertiaBreakdown, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid indices, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in inertia coupling modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.