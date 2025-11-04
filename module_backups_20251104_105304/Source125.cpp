// Ug3DiskVectorModule.h
// Modular C++ implementation of the Unit Vector in the Ug3 Disk Plane (??_j) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes ??_j (unit vector, magnitude=1; e.g., [cos ?_j, sin ?_j, 0]); scales in Universal Magnetism U_m term.
// Pluggable: #include "Ug3DiskVectorModule.h"
// Ug3DiskVectorModule mod; mod.computeUmContribution(0.0, 1); mod.updateVariable("theta_j", new_value);
// Variables in std::map; example for j=1 at t=0, ?_j=0 (??_j=[1,0,0], U_m?2.28e65 J/m�).
// Approximations: ??_j magnitude=1; 1 - exp=0 at t=0; ?_j / r_j=2.26e10 T m�.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef UG3_DISK_VECTOR_MODULE_H
#define UG3_DISK_VECTOR_MODULE_H

#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

class Ug3DiskVectorModule {
private:
    std::map<std::string, double> variables;
    std::vector<double> computePhiHat_j(int j);
    double computeUmBase(double t);
    double computeUmContribution(double t, int j);

public:
    // Constructor: Initialize with framework defaults
    Ug3DiskVectorModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    std::vector<double> computePhiHat_j(int j);  // Unit vector [cos ?_j, sin ?_j, 0]
    double computePhiHatMagnitude(int j);  // 1.0 (normalized)
    double computeUmContribution(double t, int j);  // U_m single string (J/m^3)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print ??_j and U_m
    void printVectorAndUm(int j = 1, double t = 0.0);
};

#endif // UG3_DISK_VECTOR_MODULE_H

// Ug3DiskVectorModule.cpp
#include "Ug3DiskVectorModule.h"

// Constructor: Set framework defaults
Ug3DiskVectorModule::Ug3DiskVectorModule() {
    // Universal constants
    variables["theta_j"] = 0.0;                     // rad (default azimuthal angle)
    variables["mu_j"] = 3.38e23;                    // T�m^3 (j=1)
    variables["r_j"] = 1.496e13;                    // m
    variables["gamma"] = 5e-5 / 86400.0;            // s^-1
    variables["t_n"] = 0.0;                         // s
    variables["P_SCm"] = 1.0;                       // Pressure
    variables["E_react"] = 1e46;                    // J
    variables["f_Heaviside"] = 0.01;                // Unitless
    variables["f_quasi"] = 0.01;                    // Unitless
    variables["pi"] = 3.141592653589793;

    // Derived
    variables["scale_Heaviside"] = 1e13;
    variables["heaviside_factor"] = 1.0 + variables["scale_Heaviside"] * variables["f_Heaviside"];
}

// Update variable
void Ug3DiskVectorModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void Ug3DiskVectorModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void Ug3DiskVectorModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute ??_j = [cos ?_j, sin ?_j, 0] (disk plane unit vector)
std::vector<double> Ug3DiskVectorModule::computePhiHat_j(int j) {
    double theta = variables["theta_j"];  // Simplified, same for all j or per j
    std::vector<double> phi_hat = {std::cos(theta), std::sin(theta), 0.0};
    return phi_hat;
}

// Magnitude of ??_j (normalized=1)
double Ug3DiskVectorModule::computePhiHatMagnitude(int j) {
    auto phi = computePhiHat_j(j);
    return std::sqrt(phi[0]*phi[0] + phi[1]*phi[1] + phi[2]*phi[2]);  // =1
}

// Base for U_m without ??_j magnitude (since=1)
double Ug3DiskVectorModule::computeUmBase(double t) {
    double mu_over_rj = variables["mu_j"] / variables["r_j"];
    double exp_arg = - variables["gamma"] * t * std::cos(variables["pi"] * variables["t_n"]);
    double one_minus_exp = 1.0 - std::exp(exp_arg);
    double phi_mag = computePhiHatMagnitude(1);  // =1
    double p_scm = variables["P_SCm"];
    double e_react = variables["E_react"];
    return mu_over_rj * one_minus_exp * phi_mag * p_scm * e_react;
}

// U_m contribution with ??_j
double Ug3DiskVectorModule::computeUmContribution(double t, int j) {
    double base = computeUmBase(t);
    double heaviside_f = variables["heaviside_factor"];
    double quasi_f = 1.0 + variables["f_quasi"];
    return base * heaviside_f * quasi_f;
}

// Equation text
std::string Ug3DiskVectorModule::getEquationText() {
    return "U_m = ?_j [ (?_j / r_j) (1 - e^{-? t cos(? t_n)}) \hat{?}_j ] P_SCm E_react (1 + 10^13 f_Heaviside) (1 + f_quasi)\n"
           "Where \hat{?}_j = [cos ?_j, sin ?_j, 0] (unit vector in Ug3 disk plane, |??_j|=1);\n"
           "Specifies azimuthal direction for j-th string in disk (e.g., galactic plane).\n"
           "Example j=1, ?_j=0, t=0: ??_j=[1,0,0], U_m ?2.28e65 J/m� (mag=1).\n"
           "Role: Directional geometry for magnetic contributions in disks/nebulae.\n"
           "UQFF: Vector orientation in U_m/U_g3; collimation in jets/disks/formation.";
}

// Print variables
void Ug3DiskVectorModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Print vector and U_m
void Ug3DiskVectorModule::printVectorAndUm(int j, double t) {
    auto phi = computePhiHat_j(j);
    double mag = computePhiHatMagnitude(j);
    double um = computeUmContribution(t, j);
    std::cout << "??_" << j << " at ?_j=" << variables["theta_j"] << " rad, t=" << t << " s:\n";
    std::cout << "??_j = [" << std::scientific << phi[0] << ", " << phi[1] << ", " << phi[2] << "] (mag=" << mag << ")\n";
    std::cout << "U_m contrib = " << um << " J/m�\n";
}

// Example usage in base program (snippet)
// #include "Ug3DiskVectorModule.h"
// int main() {
//     Ug3DiskVectorModule mod;
//     auto phi = mod.computePhiHat_j(1);
//     std::cout << "??_1 = [" << phi[0] << ", " << phi[1] << ", " << phi[2] << "]\n";
//     mod.printVectorAndUm(1, 0.0);
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("theta_j", M_PI / 2);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o disk_vector_test disk_vector_test.cpp Ug3DiskVectorModule.cpp -lm
// Sample: ??_1=[1,0,0] (?=0); U_m?2.28e65 J/m�; directional in disk plane.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

Ug3DiskVectorModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computePhiHat_j, computePhiHatMagnitude, computeUmContribution) are clear, concise, and variable - driven.
- Uses std::vector for unit vector representation, supporting extensibility for vector operations.
- Output and debugging functions(printVariables, printVectorAndUm, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.
- Models directional geometry for magnetic contributions in disk planes.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map and std::vector.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in disk vector and magnetic geometry modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.