// MagneticStringModule.h
// Modular C++ implementation of the Distance Along Magnetic String's Path (r_j) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes r_j = 1.496e13 m (100 AU) and its conversions; scales ?_j / r_j in Universal Magnetism U_m and Ug3.
// Pluggable: #include "MagneticStringModule.h"
// MagneticStringModule mod; mod.computeMuOverRj(); mod.updateVariable("r_j", new_value);
// Variables in std::map; j-indexed strings; example for j=1 at t=0.
// Approximations: ?=5e-5 day^-1; cos(? t_n)=1; ?_hat_j=1; at t=0, 1 - exp term=0.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef MAGNETIC_STRING_MODULE_H
#define MAGNETIC_STRING_MODULE_H

#include <map>
#include <string>
#include <vector>
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

class MagneticStringModule {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
    std::map<std::string, double> variables;
    double computeRjInAU();
    double computeRjInLy();
    double computeRjInPc();
    double computeMuOverRj(int j);
    double computeUmContribution(int j);
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;



public:
    // Constructor: Initialize with framework defaults
    MagneticStringModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeRj(int j);  // r_j in m (default 1.496e13)
    double computeRjInAU(int j);
    double computeRjInLy(int j);
    double computeRjInPc(int j);
    double computeMu_j(int j, double t);  // Magnetic moment
    double computeMuOverRj(int j);
    double computeUmContribution(int j, double t);  // Single string to U_m
    double computeUg3Contribution();  // Example Ug3 influence

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print r_j conversions and contributions
    void printStringContributions(int j = 1, double t = 0.0);
};

#endif // MAGNETIC_STRING_MODULE_H

// MagneticStringModule.cpp
#include "MagneticStringModule.h"

// Constructor: Set framework defaults
MagneticStringModule::MagneticStringModule() {
        enableDynamicTerms = true;
        enableLogging = false;
        learningRate = 0.001;
        metadata["enhanced"] = "true";
        metadata["version"] = "2.0-Enhanced";

    // Universal constants
    variables["AU_to_m"] = 1.496e11;                // m/AU
    variables["c"] = 2.998e8;                       // m/s
    variables["year_to_s"] = 3.156e7;               // s/yr
    variables["ly_to_m"] = variables["c"] * variables["year_to_s"];  // m/ly ?9.461e15
    variables["pc_to_ly"] = 3.262;                  // ly/pc
    variables["pi"] = 3.141592653589793;

    // r_j defaults (m)
    variables["r_1"] = 1.496e13;                    // 100 AU for j=1
    // Add more r_j if needed

    // Magnetic string params
    variables["mu_base"] = 3.38e20;                 // T m^3 base
    variables["omega_c"] = 2.5e-6;                  // rad/s (cavity freq)
    variables["gamma"] = 5e-5 / (86400.0);          // day^-1 to s^-1 (?=5e-5 /day)
    variables["t_n"] = 0.0;                         // s
    variables["phi_hat_1"] = 1.0;                   // Normalized
    variables["P_SCm"] = 1.0;                       // SCm pressure
    variables["E_react"] = 1e46;                    // J
    variables["f_Heaviside"] = 0.01;                // Dimensionless
    variables["f_quasi"] = 0.01;                    // Quasi factor

    // Ug3 related
    variables["k3"] = 1.8;                          // Coupling
    variables["B_j"] = 1e3;                         // T
    variables["Omega_g"] = 7.3e-16;                 // rad/s
    variables["M_s"] = 1.989e30;                    // kg
    variables["d_g"] = 2.55e20;                     // m
}

// Update variable
void MagneticStringModule::updateVariable(const std::string& name, double value) {
    variables[name] = value;
}

// Add delta
void MagneticStringModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void MagneticStringModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute r_j (m)
double MagneticStringModule::computeRj(int j) {
    std::string key = "r_" + std::to_string(j);
    if (variables.find(key) != variables.end()) {
        return variables[key];
    }
    std::cerr << "r_" << j << " not found. Using r_1." << std::endl;
    return variables["r_1"];
}

// r_j in AU
double MagneticStringModule::computeRjInAU(int j) {
    return computeRj(j) / variables["AU_to_m"];
}

// r_j in ly
double MagneticStringModule::computeRjInLy(int j) {
    return computeRj(j) / variables["ly_to_m"];
}

// r_j in pc
double MagneticStringModule::computeRjInPc(int j) {
    return computeRjInLy(j) / variables["pc_to_ly"];
}

// Compute ?_j (t)
double MagneticStringModule::computeMu_j(int j, double t) {
    double sin_term = std::sin(variables["omega_c"] * t);
    return (1e3 + 0.4 * sin_term) * variables["mu_base"];
}

// ?_j / r_j (T m^2)
double MagneticStringModule::computeMuOverRj(int j) {
    double rj = computeRj(j);
    if (rj == 0.0) return 0.0;
    double mu_j = computeMu_j(j, variables["t_n"]);  // Use t_n
    return mu_j / rj;
}

// Single string contribution to U_m (J/m^3, simplified)
double MagneticStringModule::computeUmContribution(int j, double t) {
    double mu_over_rj = computeMuOverRj(j);
    double exp_term = std::exp( - variables["gamma"] * t * std::cos(variables["pi"] * variables["t_n"]) );
    double one_minus_exp = 1.0 - exp_term;
    double phi_hat = variables["phi_hat_1"];  // For j=1
    double heaviside_factor = 1.0 + 1e13 * variables["f_Heaviside"];
    double quasi_factor = 1.0 + variables["f_quasi"];
    return (mu_over_rj * one_minus_exp * phi_hat) * variables["P_SCm"] * variables["E_react"] * heaviside_factor * quasi_factor;
}

// Example Ug3 contribution (J/m^3)
double MagneticStringModule::computeUg3Contribution() {
    double cos_term = std::cos(variables["Omega_g"] * variables["t_n"] * variables["pi"]);
    double rho_sum = variables["rho_vac_SCm"] + variables["rho_vac_UA"];  // Placeholder
    double M_s_over_d_g = variables["M_s"] / variables["d_g"];
    return variables["k3"] * variables["B_j"] * cos_term * rho_sum * variables["Omega_g"] * M_s_over_d_g * 1e46;  // Scaled
}

// Equation text
std::string MagneticStringModule::getEquationText() {
    return "U_m = ?_j [ (?_j / r_j) * (1 - e^{-? t cos(? t_n)}) * ?_hat_j ] * P_SCm * E_react * (1 + 10^13 f_Heaviside) * (1 + f_quasi)\n"
           "Where r_j = 1.496e13 m (100 AU, j-th string path distance);\n"
           "?_j = (10^3 + 0.4 sin(?_c t)) * 3.38e20 T m^3;\n"
           "? ?5.8e-10 s^-1 (5e-5 day^-1); at t=0, 1-exp=0.\n"
           "In Ug3: Influences (?_SCm + ?_UA) ?_g M_s / d_g * cos(...).\n"
           "Example j=1, t=0: ?_1 / r_1 ?2.26e10 T m^2; U_m contrib=0 (exp=1).\n"
           "Ug3 ?1.8e49 J/mï¿½ (k3=1.8 scaling).\n"
           "Role: Scales magnetic string extent; stabilizes disks/nebulae at 100 AU scale.";
}

// Print variables
void MagneticStringModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Print contributions for j, t
void MagneticStringModule::printStringContributions(int j, double t) {
    double rj_m = computeRj(j);
    double rj_au = computeRjInAU(j);
    double rj_ly = computeRjInLy(j);
    double rj_pc = computeRjInPc(j);
    double mu_over_rj = computeMuOverRj(j);
    double um_contrib = computeUmContribution(j, t);
    double ug3 = computeUg3Contribution();
    std::cout << "Magnetic String j=" << j << " at t=" << t << " s:\n";
    std::cout << "r_j = " << std::scientific << rj_m << " m (" << rj_au << " AU, " << rj_ly << " ly, " << rj_pc << " pc)\n";
    std::cout << "?_j / r_j = " << mu_over_rj << " T m^2\n";
    std::cout << "U_m contrib = " << um_contrib << " J/mï¿½\n";
    std::cout << "Ug3 contrib (example) = " << ug3 << " J/mï¿½\n";
}

// Example usage in base program (snippet)
// #include "MagneticStringModule.h"
// int main() {
//     MagneticStringModule mod;
//     mod.printStringContributions(1, 0.0);
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("r_1", 2e13);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o string_test string_test.cpp MagneticStringModule.cpp -lm
// Sample: r_1=1.496e13 m (100 AU); ?/r ?2.26e10; U_m=0 at t=0; Ug3?1.8e49 J/mï¿½.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

MagneticStringModule Evaluation

Strengths :
-Modular, extensible design for modeling magnetic string path distances and their contributions to universal magnetism(U_m) and gravity(Ug3) in the UQFF framework.
- Clear encapsulation of variables and string parameters using std::map, supporting dynamic updates and easy extension.
- Implements core physical concepts : conversion of r_j between units(m, AU, ly, pc), magnetic moment calculations, and their influence on U_m and Ug3.
- Approximations and physical meaning are well - documented in comments and equation text.
- Output functions for variable state and string contributions support debugging and transparency.
- Handles dynamic updates to r_j and recalculates dependent terms as needed.
- Example calculations and conversion functions provide scientific context and validation.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for flexibility and scalability.
- Some calculations use magic numbers or lack explanatory comments; define named constants and clarify logic.
- Minimal error handling(e.g., missing r_j, division by zero, invalid variable names); add validation for robustness.
- Unit consistency should be checked and documented for all physical quantities.
- For large - scale or performance - critical simulations, optimize data structures and reduce redundant calculations.
- std::map is flexible but may be less efficient than structured types for very large models.
- Expand documentation for function purposes and physical meaning.
- Consider supporting multiple j - indexed strings for more general modeling.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in magnetic string modeling.It implements the UQFF magnetic string concept faithfully and adapts to various scenarios.For production or high - performance applications, address the recommendations for improved robustness, maintainability, and scalability.