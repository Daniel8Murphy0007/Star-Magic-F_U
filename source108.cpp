// CorePenetrationModule.h
// Modular C++ implementation of the Planetary Core Penetration Factor (P_core) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes P_core ?1 (unitless for Sun, ~1e-3 for planets); scales P_core in Universal Gravity U_g3 term.
// Pluggable: #include "CorePenetrationModule.h"
// CorePenetrationModule mod; mod.computeU_g3(0.0); mod.updateVariable("P_core", new_value);
// Variables in std::map; example for Sun at t=0; planet mode with P_core=1e-3.
// Approximations: cos(?_s t ?)=1 at t=0; E_react=1e46; B_j=1e3 T.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef CORE_PENETRATION_MODULE_H
#define CORE_PENETRATION_MODULE_H

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

class CorePenetrationModule {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
    std::map<std::string, double> variables;
    double computeU_g3(double t);
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;



public:
    // Constructor: Initialize with framework defaults (Sun)
    CorePenetrationModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeP_core();  // ?1 for Sun (unitless)
    double computeU_g3(double t);  // U_g3 with P_core (J/m^3)
    double computeU_g3_planet(double t);  // For planet P_core=1e-3

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();
};

#endif // CORE_PENETRATION_MODULE_H

// CorePenetrationModule.cpp
#include "CorePenetrationModule.h"

// Constructor: Set framework defaults (Sun at t=0)
CorePenetrationModule::CorePenetrationModule() {
        enableDynamicTerms = true;
        enableLogging = false;
        learningRate = 0.001;
        metadata["enhanced"] = "true";
        metadata["version"] = "2.0-Enhanced";

    // Universal constants
    variables["P_core"] = 1.0;                      // Unitless ?1 for Sun
    variables["k_3"] = 1.8;                         // Coupling
    variables["B_j"] = 1e3;                         // T (base)
    variables["omega_s"] = 2.5e-6;                  // rad/s
    variables["P_core_planet"] = 1e-3;              // For planets
    variables["E_react"] = 1e46;                    // J
    variables["pi"] = 3.141592653589793;
    variables["t"] = 0.0;                           // s
}

// Update variable
void CorePenetrationModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void CorePenetrationModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void CorePenetrationModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute P_core ?1
double CorePenetrationModule::computeP_core() {
    return variables["P_core"];
}

// Compute U_g3 with P_core
double CorePenetrationModule::computeU_g3(double t) {
    variables["t"] = t;
    double k_3 = variables["k_3"];
    double b_j = variables["B_j"];  // Simplified, no sin term in example
    double cos_term = std::cos(variables["omega_s"] * t * variables["pi"]);
    double p_core = computeP_core();
    double e_react = variables["E_react"];
    return k_3 * b_j * cos_term * p_core * e_react;
}

// U_g3 for planet (P_core=1e-3)
double CorePenetrationModule::computeU_g3_planet(double t) {
    double orig_p = variables["P_core"];
    variables["P_core"] = variables["P_core_planet"];
    double result = computeU_g3(t);
    variables["P_core"] = orig_p;
    return result;
}

// Equation text
std::string CorePenetrationModule::getEquationText() {
    return "U_g3 = k_3 * ?_j B_j(r,?,t,?_vac,[SCm]) * cos(?_s(t) t ?) * P_core * E_react\n"
           "Where P_core ?1 (unitless for Sun, ~1e-3 for planets; core penetration).\n"
           "Scales magnetic disk gravity for core [SCm] influence.\n"
           "Example Sun t=0: U_g3 ?1.8e49 J/mï¿½ (P_core=1);\n"
           "Planet: ?1.8e46 J/mï¿½ (P_core=1e-3, -3 orders).\n"
           "Role: Adjusts core interactions; full for stellar plasma, reduced for solid cores.\n"
           "UQFF: Models penetration in nebulae/star formation/disks.";
}

// Print variables
void CorePenetrationModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Example usage in base program (snippet)
// #include "CorePenetrationModule.h"
// int main() {
//     CorePenetrationModule mod;
//     double p = mod.computeP_core();
//     std::cout << "P_core ? " << p << std::endl;
//     double u_g3_sun = mod.computeU_g3(0.0);
//     std::cout << "U_g3 (Sun) = " << u_g3_sun << " J/mï¿½\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("P_core", 1e-3);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o core_test core_test.cpp CorePenetrationModule.cpp -lm
// Sample: P_core=1; U_g3?1.8e49 J/mï¿½ (Sun); scales for planetary cores.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

CorePenetrationModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeP_core, computeU_g3, computeU_g3_planet) are clear, concise, and variable - driven.
- Output and debugging functions(printVariables, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.
- Supports both stellar and planetary scenarios by switching P_core values.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in core penetration modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.