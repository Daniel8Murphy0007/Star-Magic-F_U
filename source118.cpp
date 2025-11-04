// StellarRotationModule.h
// Modular C++ implementation of the Stellar/Planetary Rotation Rate (?_s) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes ?_s=2.5e-6 rad/s (~29-day Sun period); scales ?_s(t) in U_g3 cos(?_s t ?) and U_i ?_s cos(? t_n).
// Pluggable: #include "StellarRotationModule.h"
// StellarRotationModule mod; mod.computeU_g3(0.0); mod.updateVariable("omega_s", new_value);
// Variables in std::map; example for Sun at t=0, t_n=0; U_g3 ?1.8e49 J/mï¿½, U_i ?1.38e-47 J/mï¿½.
// Approximations: cos(? t_n)=1; f_TRZ=0.1; ?_i=1.0; ?_vac sum=7.80e-36 J/mï¿½.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef STELLAR_ROTATION_MODULE_H
#define STELLAR_ROTATION_MODULE_H

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

class StellarRotationModule {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
    std::map<std::string, double> variables;
    double computeOmega_s_t(double t);  // ?_s(t), simplified constant
    double computeU_g3(double t);
    double computeU_i(double t, double t_n);
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;



public:
    // Constructor: Initialize with framework defaults (Sun)
    StellarRotationModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeOmega_s();  // 2.5e-6 rad/s
    double computeOmega_s_t(double t);  // ?_s(t) (rad/s)
    double computePeriod_days();  // ~29 days
    double computeU_g3(double t);  // U_g3 example (J/m^3)
    double computeU_i(double t, double t_n);  // U_i example (J/m^3)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();
};

#endif // STELLAR_ROTATION_MODULE_H

// StellarRotationModule.cpp
#include "StellarRotationModule.h"

// Constructor: Set framework defaults (Sun at t=0)
StellarRotationModule::StellarRotationModule() {
        enableDynamicTerms = true;
        enableLogging = false;
        learningRate = 0.001;
        metadata["enhanced"] = "true";
        metadata["version"] = "2.0-Enhanced";

    // Universal constants
    variables["omega_s"] = 2.5e-6;                  // rad/s
    variables["k_3"] = 1.8;                         // Coupling U_g3
    variables["B_j"] = 1e3;                         // T
    variables["P_core"] = 1.0;                      // Unitless
    variables["E_react"] = 1e46;                    // J
    variables["lambda_i"] = 1.0;                    // Unitless U_i
    variables["rho_vac_SCm"] = 7.09e-37;            // J/m^3
    variables["rho_vac_UA"] = 7.09e-36;             // J/m^3
    variables["f_TRZ"] = 0.1;                       // Unitless
    variables["pi"] = 3.141592653589793;
    variables["t"] = 0.0;                           // s
    variables["t_n"] = 0.0;                         // s

    // Derived
    variables["rho_sum"] = variables["rho_vac_SCm"] + variables["rho_vac_UA"];
    variables["day_to_s"] = 86400.0;                // s/day
}

// Update variable
void StellarRotationModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
        if (name == "rho_vac_SCm" || name == "rho_vac_UA") {
            variables["rho_sum"] = variables["rho_vac_SCm"] + variables["rho_vac_UA"];
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void StellarRotationModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        if (name == "rho_vac_SCm" || name == "rho_vac_UA") {
            variables["rho_sum"] = variables["rho_vac_SCm"] + variables["rho_vac_UA"];
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void StellarRotationModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute ?_s (rad/s)
double StellarRotationModule::computeOmega_s() {
    return variables["omega_s"];
}

// ?_s(t), simplified as constant (no t dep in example)
double StellarRotationModule::computeOmega_s_t(double t) {
    variables["t"] = t;
    return computeOmega_s();  // Constant for Sun
}

// Rotation period in days
double StellarRotationModule::computePeriod_days() {
    double period_s = 2.0 * M_PI / computeOmega_s();
    return period_s / variables["day_to_s"];
}

// U_g3 example
double StellarRotationModule::computeU_g3(double t) {
    double k_3 = variables["k_3"];
    double b_j = variables["B_j"];
    double cos_term = std::cos(variables["omega_s"] * t * variables["pi"]);
    double p_core = variables["P_core"];
    double e_react = variables["E_react"];
    return k_3 * b_j * cos_term * p_core * e_react;
}

// U_i example
double StellarRotationModule::computeU_i(double t, double t_n) {
    double lambda_i = variables["lambda_i"];
    double rho_sc = variables["rho_vac_SCm"];
    double rho_ua = variables["rho_vac_UA"];
    double omega_s_t = computeOmega_s_t(t);
    double cos_pi_tn = std::cos(variables["pi"] * t_n);
    double trz_factor = 1.0 + variables["f_TRZ"];
    return lambda_i * rho_sc * rho_ua * omega_s_t * cos_pi_tn * trz_factor;
}

// Equation text
std::string StellarRotationModule::getEquationText() {
    return "U_g3 = k_3 * ? B_j * cos(?_s(t) t ?) * P_core * E_react\n"
           "U_i = ?_i * ?_vac,[SCm] * ?_vac,[UA] * ?_s(t) * cos(? t_n) * (1 + f_TRZ)\n"
           "Where ?_s = 2.5e-6 rad/s (~29-day Sun equatorial rotation);\n"
           "Scales rotational oscillations/inertia.\n"
           "Example t=0, t_n=0: U_g3 ?1.8e49 J/mï¿½; U_i ?1.38e-47 J/mï¿½.\n"
           "Role: Introduces spin in disk gravity/inertia; stellar/planetary dynamics.\n"
           "UQFF: Rotational effects in nebulae/disks/formation/mergers.";
}

// Print variables
void StellarRotationModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Example usage in base program (snippet)
// #include "StellarRotationModule.h"
// int main() {
//     StellarRotationModule mod;
//     double omega = mod.computeOmega_s();
//     std::cout << "?_s = " << omega << " rad/s (~" << mod.computePeriod_days() << " days)\n";
//     double u_g3 = mod.computeU_g3(0.0);
//     std::cout << "U_g3 = " << u_g3 << " J/mï¿½\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("omega_s", 3e-6);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o rotation_test rotation_test.cpp StellarRotationModule.cpp -lm
// Sample: ?_s=2.5e-6 rad/s (~29 days); U_g3?1.8e49 J/mï¿½; scales rotation.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

StellarRotationModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeOmega_s, computeOmega_s_t, computePeriod_days, computeU_g3, computeU_i) are clear, concise, and variable - driven.
- Automatic recalculation of derived variables(rho_sum) when dependencies change.
- Output and debugging functions(printVariables, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.
- Models both gravity and inertia effects of stellar / planetary rotation.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in stellar / planetary rotation modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.