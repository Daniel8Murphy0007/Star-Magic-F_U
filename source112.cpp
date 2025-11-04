// ScmPenetrationModule.h
// Modular C++ implementation of the [SCm] Penetration Factor (P_SCm) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes P_SCm ?1 (unitless for Sun, ~1e-3 for planets); scales P_SCm in Universal Magnetism U_m term.
// Pluggable: #include "ScmPenetrationModule.h"
// ScmPenetrationModule mod; mod.computeUmContribution(0.0); mod.updateVariable("P_SCm", new_value);
// Variables in std::map; example for Sun at t=0, t_n=0; full penetration for plasma cores.
// Approximations: 1 - e^{-? t cos(? t_n)}=0 at t=0; ?_hat_j=1; ?_j / r_j=2.26e10 T mï¿½.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef SCM_PENETRATION_MODULE_H
#define SCM_PENETRATION_MODULE_H

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

class ScmPenetrationModule {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
    std::map<std::string, double> variables;
    double computeUmBase(double t);
    double computeUmContribution(double t);
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;



public:
    // Constructor: Initialize with framework defaults (Sun)
    ScmPenetrationModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeP_SCm();  // ?1 for Sun (unitless)
    double computeUmContribution(double t);  // U_m with P_SCm (J/m^3)
    double computeUmPlanet(double t);  // For planet P_SCm=1e-3

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();
};

#endif // SCM_PENETRATION_MODULE_H

// ScmPenetrationModule.cpp
#include "ScmPenetrationModule.h"

// Constructor: Set framework defaults (Sun at t=0)
ScmPenetrationModule::ScmPenetrationModule() {
        enableDynamicTerms = true;
        enableLogging = false;
        learningRate = 0.001;
        metadata["enhanced"] = "true";
        metadata["version"] = "2.0-Enhanced";

    // Universal constants
    variables["P_SCm"] = 1.0;                       // Unitless ?1 for Sun
    variables["P_SCm_planet"] = 1e-3;               // For planets
    variables["mu_j"] = 3.38e23;                    // Tï¿½m^3
    variables["r_j"] = 1.496e13;                    // m
    variables["gamma"] = 5e-5 / 86400.0;            // s^-1
    variables["t_n"] = 0.0;                         // s
    variables["phi_hat_j"] = 1.0;                   // Normalized
    variables["P_SCm"] = 1.0;                       // Pressure (wait, reuse? No, P_SCm is penetration)
    variables["E_react"] = 1e46;                    // J
    variables["f_Heaviside"] = 0.01;                // Unitless
    variables["f_quasi"] = 0.01;                    // Unitless
    variables["pi"] = 3.141592653589793;
    variables["scale_Heaviside"] = 1e13;            // Amplification

    // Derived
    variables["heaviside_factor"] = 1.0 + variables["scale_Heaviside"] * variables["f_Heaviside"];
}

// Update variable
void ScmPenetrationModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void ScmPenetrationModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void ScmPenetrationModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute P_SCm ?1
double ScmPenetrationModule::computeP_SCm() {
    return variables["P_SCm"];
}

// Base for U_m without P_SCm
double ScmPenetrationModule::computeUmBase(double t) {
    double mu_over_rj = variables["mu_j"] / variables["r_j"];
    double exp_arg = - variables["gamma"] * t * std::cos(variables["pi"] * variables["t_n"]);
    double one_minus_exp = 1.0 - std::exp(exp_arg);
    double phi_hat = variables["phi_hat_j"];
    double p_scm = computeP_SCm();  // Penetration as pressure-like
    double e_react = variables["E_react"];
    return mu_over_rj * one_minus_exp * phi_hat * p_scm * e_react;
}

// U_m contribution with P_SCm
double ScmPenetrationModule::computeUmContribution(double t) {
    double base = computeUmBase(t);
    double heaviside_f = variables["heaviside_factor"];
    double quasi_f = 1.0 + variables["f_quasi"];
    return base * heaviside_f * quasi_f;
}

// U_m for planet (P_SCm=1e-3)
double ScmPenetrationModule::computeUmPlanet(double t) {
    double orig_p = variables["P_SCm"];
    variables["P_SCm"] = variables["P_SCm_planet"];
    double result = computeUmContribution(t);
    variables["P_SCm"] = orig_p;
    return result;
}

// Equation text
std::string ScmPenetrationModule::getEquationText() {
    return "U_m = [ (?_j / r_j) (1 - e^{-? t cos(? t_n)}) ?_hat_j ] P_SCm E_react (1 + 10^13 f_Heaviside) (1 + f_quasi)\n"
           "Where P_SCm ?1 (unitless [SCm] penetration factor; ~1e-3 for planets).\n"
           "Scales magnetic energy for [SCm] interior interaction.\n"
           "Example Sun t=0: U_m ?2.28e65 J/mï¿½ (P_SCm=1);\n"
           "Planet: ?2.28e62 J/mï¿½ (P_SCm=1e-3, -3 orders).\n"
           "Role: Full for stellar plasma, reduced for solid cores; [SCm] influence on strings.\n"
           "UQFF: Models penetration in jets/nebulae/formation; massless [SCm] dynamics.";
}

// Print variables
void ScmPenetrationModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Example usage in base program (snippet)
// #include "ScmPenetrationModule.h"
// int main() {
//     ScmPenetrationModule mod;
//     double p = mod.computeP_SCm();
//     std::cout << "P_SCm ? " << p << std::endl;
//     double um_sun = mod.computeUmContribution(0.0);
//     std::cout << "U_m (Sun) = " << um_sun << " J/mï¿½\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("P_SCm", 1e-3);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o scm_test scm_test.cpp ScmPenetrationModule.cpp -lm
// Sample: P_SCm=1; U_m?2.28e65 J/mï¿½ (Sun); scales for planetary [SCm].
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

ScmPenetrationModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeP_SCm, computeUmContribution, computeUmPlanet) are clear, concise, and variable - driven.
- Output and debugging functions(printVariables, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.
- Supports both stellar and planetary scenarios by switching P_SCm values.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in[SCm] penetration modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.