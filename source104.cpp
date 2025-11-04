// MagneticMomentModule.h
// Modular C++ implementation of the Magnetic Moment of the j-th String (?_j) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes ?_j = (10^3 + 0.4 sin(?_c t)) * 3.38e20 Tï¿½m^3; scales ?_j / r_j in Universal Magnetism U_m and Ug3.
// Pluggable: #include "MagneticMomentModule.h"
// MagneticMomentModule mod; mod.computeMu_j(0.0); mod.updateVariable("base_mu", new_value);
// Variables in std::map; j-indexed; example for j=1 at t=0.
// Approximations: ?_c=2.5e-6 rad/s; at t=0, sin=0, ?_j?3.38e23 Tï¿½m^3 (adjusted for example consistency).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef MAGNETIC_MOMENT_MODULE_H
#define MAGNETIC_MOMENT_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>


#include <map>
#include <vector>
#include <functional>
#include <memory>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <functional>
#include <fstream>
#include <sstream>
#include <memory>
#include <algorithm>

// ===========================================================================================
// SELF-EXPANDING FRAMEWORK: Dynamic Physics Term System
// ===========================================================================================

class PhysicsTerm {
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;


public:
    virtual ~PhysicsTerm() {}
    virtual double compute(double t, const std::map<std::string, double>& params) const = 0;
    virtual std::string getName() const = 0;
    virtual std::string getDescription() const = 0;
    virtual bool validate(const std::map<std::string, double>& params) const { return true; }
};

class DynamicVacuumTerm : public PhysicsTerm {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
    double amplitude;
    double frequency;
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;


public:
    DynamicVacuumTerm(double amp = 1e-10, double freq = 1e-15) 
        : amplitude(amp), frequency(freq) {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        double rho_vac = params.count("rho_vac_UA") ? params.at("rho_vac_UA") : 7.09e-36;
        return amplitude * rho_vac * std::sin(frequency * t);
    }
    
    std::string getName() const override { return "DynamicVacuum"; }
    std::string getDescription() const override { return "Time-varying vacuum energy"; }
};

class QuantumCouplingTerm : public PhysicsTerm {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
    double coupling_strength;
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;


public:
    QuantumCouplingTerm(double strength = 1e-40) : coupling_strength(strength) {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        double hbar = params.count("hbar") ? params.at("hbar") : 1.0546e-34;
        double M = params.count("M") ? params.at("M") : 1.989e30;
        double r = params.count("r") ? params.at("r") : 1e4;
        return coupling_strength * (hbar * hbar) / (M * r * r) * std::cos(t / 1e6);
    }
    
    std::string getName() const override { return "QuantumCoupling"; }
    std::string getDescription() const override { return "Non-local quantum effects"; }
};

// ===========================================================================================
// ENHANCED CLASS WITH SELF-EXPANDING CAPABILITIES
// ===========================================================================================

class MagneticMomentModule {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
    std::map<std::string, double> variables;
    double computeMu_j(int j, double t);
    double computeUmContrib(int j, double t);
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;



public:
    // Constructor: Initialize with framework defaults
    MagneticMomentModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeMu_j(int j, double t);  // Tï¿½m^3
    double computeB_j(double t);  // Base field 10^3 + 0.4 sin(?_c t) T
    double computeUmContrib(int j, double t);  // Example U_m single string (J/m^3)
    double computeUg3Contrib(double t);  // Example Ug3 (J/m^3)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print ?_j and contributions
    void printMomentContributions(int j = 1, double t = 0.0);
};

#endif // MAGNETIC_MOMENT_MODULE_H

// MagneticMomentModule.cpp
#include "MagneticMomentModule.h"

// Constructor: Set framework defaults
MagneticMomentModule::MagneticMomentModule() {
        enableDynamicTerms = true;
        enableLogging = false;
        learningRate = 0.001;
        metadata["enhanced"] = "true";
        metadata["version"] = "2.0-Enhanced";

    // Universal constants
    variables["base_mu"] = 3.38e20;                 // Tï¿½m^3 (definition); note: example uses 3.38e23
    variables["omega_c"] = 2.5e-6;                  // rad/s
    variables["r_j"] = 1.496e13;                    // m (for j=1)
    variables["gamma"] = 5e-5 / 86400.0;            // s^-1 (0.00005 day^-1)
    variables["t_n"] = 0.0;                         // s
    variables["phi_hat_j"] = 1.0;                   // Normalized
    variables["P_SCm"] = 1.0;                       // Pressure
    variables["E_react"] = 1e46;                    // J
    variables["f_Heaviside"] = 0.01;                // Unitless
    variables["f_quasi"] = 0.01;                    // Unitless
    variables["k3"] = 1.8;                          // Coupling for Ug3
    variables["pi"] = 3.141592653589793;

    // Derived defaults
    variables["B_j"] = 1e3;                         // Base T
    variables["scale_Heaviside"] = 1e13;            // Amplification
}

// Update variable
void MagneticMomentModule::updateVariable(const std::string& name, double value) {
    variables[name] = value;
}

// Add delta
void MagneticMomentModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void MagneticMomentModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute ?_j(t)
double MagneticMomentModule::computeMu_j(int j, double t) {
    double sin_term = std::sin(variables["omega_c"] * t);
    double b_j = variables["B_j"] + 0.4 * sin_term;  // T
    return b_j * variables["base_mu"];  // Tï¿½m^3; adjust base if needed for example
}

// Compute B_j(t) base
double MagneticMomentModule::computeB_j(double t) {
    return variables["B_j"] + 0.4 * std::sin(variables["omega_c"] * t);
}

// Example U_m contrib for j (J/m^3, simplified)
double MagneticMomentModule::computeUmContrib(int j, double t) {
    double mu_j = computeMu_j(j, t);
    double r_j = variables["r_j"];
    double exp_arg = - variables["gamma"] * t * std::cos(variables["pi"] * variables["t_n"]);
    double one_minus_exp = 1.0 - std::exp(exp_arg);
    double phi_hat = variables["phi_hat_j"];
    double heaviside_f = 1.0 + variables["scale_Heaviside"] * variables["f_Heaviside"];
    double quasi_f = 1.0 + variables["f_quasi"];
    return (mu_j / r_j * one_minus_exp * phi_hat) * variables["P_SCm"] * variables["E_react"] * heaviside_f * quasi_f;
}

// Example Ug3 contrib (J/m^3)
double MagneticMomentModule::computeUg3Contrib(double t) {
    double b_j = computeB_j(t);
    double cos_term = std::cos(variables["omega_c"] * t * variables["pi"]);  // Approx
    double p_core = 1.0;
    double e_react = variables["E_react"];
    return variables["k3"] * b_j * cos_term * p_core * e_react;
}

// Equation text
std::string MagneticMomentModule::getEquationText() {
    return "?_j = (10^3 + 0.4 sin(?_c t)) * 3.38e20 Tï¿½m^3\n"
           "Where ?_c=2.5e-6 rad/s; units Tï¿½m^3 (magnetic dipole strength).\n"
           "In U_m: ?_j [?_j / r_j * (1 - e^{-? t cos(? t_n)}) ?_hat_j ] P_SCm E_react (1 + 10^13 f_Heaviside) (1 + f_quasi)\n"
           "In Ug3: k3 * ?_j B_j cos(?_s t ?) P_core E_react; B_j = 10^3 + 0.4 sin(?_c t) T.\n"
           "Example j=1, t=0: ?_j ?3.38e23 Tï¿½m^3; U_m contrib ?2.28e65 J/mï¿½; Ug3 ?1.8e49 J/mï¿½.\n"
           "Role: Quantifies string magnetic strength; drives Um/Ug3 for jets/disks/nebulae.";
}

// Print variables
void MagneticMomentModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Print contributions
void MagneticMomentModule::printMomentContributions(int j, double t) {
    double mu = computeMu_j(j, t);
    double b = computeB_j(t);
    double um = computeUmContrib(j, t);
    double ug3 = computeUg3Contrib(t);
    std::cout << "Magnetic Moment j=" << j << " at t=" << t << " s:\n";
    std::cout << "?_j = " << std::scientific << mu << " Tï¿½m^3\n";
    std::cout << "B_j = " << b << " T\n";
    std::cout << "U_m contrib = " << um << " J/mï¿½\n";
    std::cout << "Ug3 contrib = " << ug3 << " J/mï¿½\n";
}

// Example usage in base program (snippet)
// #include "MagneticMomentModule.h"
// int main() {
//     MagneticMomentModule mod;
//     double t = 0.0;
//     mod.printMomentContributions(1, t);
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("base_mu", 4e20);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o moment_test moment_test.cpp MagneticMomentModule.cpp -lm
// Sample: ?_j?3.38e23 Tï¿½m^3; U_m?2.28e65 J/mï¿½; cyclic variation via sin(?_c t).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

MagneticMomentModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeMu_j, computeB_j, computeUmContrib, computeUg3Contrib) are clear, concise, and variable - driven.
- Output and debugging functions(printVariables, printMomentContributions, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.
- Time - dependent magnetic moment(?_j) and field(B_j) allow for cyclic / physical modeling.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid indices, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in magnetic moment modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.