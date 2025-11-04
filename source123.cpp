// TimeReversalZoneModule.h
// Modular C++ implementation of the Time-Reversal Zone Factor (f_TRZ) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes f_TRZ=0.1 (unitless); scales (1 + f_TRZ) in Universal Inertia U_i term for TRZ enhancement.
// Pluggable: #include "TimeReversalZoneModule.h"
// TimeReversalZoneModule mod; mod.computeU_i(0.0, 0.0); mod.updateVariable("f_TRZ", new_value);
// Variables in std::map; example for Sun at t=0, t_n=0; U_i ?1.38e-47 J/mï¿½ (with, +10%); without: ?1.25e-47 J/mï¿½.
// Approximations: ?_i=1.0; cos(? t_n)=1; ?_s=2.5e-6 rad/s; ?_sum=7.80e-36 J/mï¿½.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef TIME_REVERSAL_ZONE_MODULE_H
#define TIME_REVERSAL_ZONE_MODULE_H

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

class TimeReversalZoneModule {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
    std::map<std::string, double> variables;
    double computeTRZFactor();
    double computeU_i_base(double t, double t_n);
    double computeU_i(double t, double t_n);
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;



public:
    // Constructor: Initialize with framework defaults
    TimeReversalZoneModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeF_TRZ();  // 0.1 (unitless)
    double computeTRZFactor();  // 1 + f_TRZ = 1.1
    double computeU_i(double t, double t_n);  // U_i with TRZ (J/m^3)
    double computeU_i_no_TRZ(double t, double t_n);  // Without TRZ

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print U_i comparison (with/without TRZ)
    void printUiComparison(double t = 0.0, double t_n = 0.0);
};

#endif // TIME_REVERSAL_ZONE_MODULE_H

// TimeReversalZoneModule.cpp
#include "TimeReversalZoneModule.h"

// Constructor: Set framework defaults (Sun at t=0, level 13)
TimeReversalZoneModule::TimeReversalZoneModule() {
        enableDynamicTerms = true;
        enableLogging = false;
        learningRate = 0.001;
        metadata["enhanced"] = "true";
        metadata["version"] = "2.0-Enhanced";

    // Universal constants
    variables["f_TRZ"] = 0.1;                       // Unitless TRZ factor
    variables["lambda_i"] = 1.0;                    // Coupling
    variables["rho_vac_SCm"] = 7.09e-37;            // J/m^3
    variables["rho_vac_UA"] = 7.09e-36;             // J/m^3
    variables["omega_s"] = 2.5e-6;                  // rad/s
    variables["pi"] = 3.141592653589793;
    variables["t"] = 0.0;                           // s
    variables["t_n"] = 0.0;                         // s

    // Derived
    variables["rho_product"] = variables["rho_vac_SCm"] * variables["rho_vac_UA"];
    variables["trz_factor"] = computeTRZFactor();
}

// Update variable
void TimeReversalZoneModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
        if (name == "f_TRZ") {
            variables["trz_factor"] = computeTRZFactor();
        } else if (name == "rho_vac_SCm" || name == "rho_vac_UA") {
            variables["rho_product"] = variables["rho_vac_SCm"] * variables["rho_vac_UA"];
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void TimeReversalZoneModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        if (name == "f_TRZ") {
            variables["trz_factor"] = computeTRZFactor();
        } else if (name == "rho_vac_SCm" || name == "rho_vac_UA") {
            variables["rho_product"] = variables["rho_vac_SCm"] * variables["rho_vac_UA"];
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void TimeReversalZoneModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute f_TRZ (0.1)
double TimeReversalZoneModule::computeF_TRZ() {
    return variables["f_TRZ"];
}

// Compute 1 + f_TRZ
double TimeReversalZoneModule::computeTRZFactor() {
    return 1.0 + computeF_TRZ();
}

// Base U_i without TRZ factor
double TimeReversalZoneModule::computeU_i_base(double t, double t_n) {
    double lambda_i = variables["lambda_i"];
    double rho_product = variables["rho_product"];
    double omega_s_t = variables["omega_s"];        // Simplified constant
    double cos_pi_tn = std::cos(variables["pi"] * t_n);
    return lambda_i * rho_product * omega_s_t * cos_pi_tn;
}

// U_i with TRZ
double TimeReversalZoneModule::computeU_i(double t, double t_n) {
    variables["t"] = t;
    double base = computeU_i_base(t, t_n);
    double trz_f = computeTRZFactor();
    return base * trz_f;
}

// U_i without TRZ (f=0)
double TimeReversalZoneModule::computeU_i_no_TRZ(double t, double t_n) {
    double orig_f = variables["f_TRZ"];
    variables["f_TRZ"] = 0.0;
    double result = computeU_i(t, t_n);
    variables["f_TRZ"] = orig_f;
    return result;
}

// Equation text
std::string TimeReversalZoneModule::getEquationText() {
    return "U_i = ?_i * ?_vac,[SCm] * ?_vac,[UA] * ?_s(t) * cos(? t_n) * (1 + f_TRZ)\n"
           "Where f_TRZ = 0.1 (unitless time-reversal zone factor; +10% negentropic enhancement);\n"
           "TRZ: Regions for time-reversal/negentropy (COP>1, vacuum extraction).\n"
           "Example Sun t=0, t_n=0: U_i ?1.38e-47 J/mï¿½ (with); ?1.25e-47 J/mï¿½ (without; -9.1%).\n"
           "In F_U: -? ?_i U_i E_react (resistive, TRZ-boosted).\n"
           "Role: Stabilizes via negentropy; TRZ in nebulae/formation/mergers/biology.\n"
           "UQFF: Integrates pondermotive force/time asymmetry; Aether superfluid effects.";
}

// Print variables
void TimeReversalZoneModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Print comparison
void TimeReversalZoneModule::printUiComparison(double t, double t_n) {
    double u_i_with = computeU_i(t, t_n);
    double u_i_without = computeU_i_no_TRZ(t, t_n);
    double percent_increase = ((u_i_with - u_i_without) / u_i_without) * 100.0;
    std::cout << "U_i Comparison at t=" << t << " s, t_n=" << t_n << ":\n";
    std::cout << "With TRZ: " << std::scientific << u_i_with << " J/mï¿½\n";
    std::cout << "Without TRZ: " << u_i_without << " J/mï¿½\n";
    std::cout << "Increase: +" << std::fixed << std::setprecision(1) << percent_increase << "%\n";
}

// Example usage in base program (snippet)
// #include "TimeReversalZoneModule.h"
// int main() {
//     TimeReversalZoneModule mod;
//     double f_trz = mod.computeF_TRZ();
//     std::cout << "f_TRZ = " << f_trz << std::endl;
//     mod.printUiComparison(0.0, 0.0);
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("f_TRZ", 0.2);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o trz_test trz_test.cpp TimeReversalZoneModule.cpp -lm
// Sample: f_TRZ=0.1; U_i with=1.38e-47 J/mï¿½ (+10% vs without).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

TimeReversalZoneModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeF_TRZ, computeTRZFactor, computeU_i, computeU_i_no_TRZ) are clear, concise, and variable - driven.
- Automatic recalculation of derived variables(trz_factor, rho_product) when dependencies change.
- Output and debugging functions(printVariables, printUiComparison, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.
- Models time - reversal zone enhancement and its effect on inertia terms.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in time - reversal zone factor modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.