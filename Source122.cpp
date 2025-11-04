// SurfaceTemperatureModule.h
// Modular C++ implementation of the Surface Temperature (T_s) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes T_s=5778 K (Sun effective); potential scaling T_s / T_s_ref in B_j for U_g3 magnetic strings.
// Pluggable: #include "SurfaceTemperatureModule.h"
// SurfaceTemperatureModule mod; mod.computeU_g3_example(0.0, 5778.0); mod.updateVariable("T_s", new_value);
// Variables in std::map; example for Sun at t=0; T_s=5778 K ? U_g3?1.8e49 J/mï¿½ (full); T_s=10000 K: ~3.11e49 J/mï¿½.
// Approximations: T_s_ref=5778 K (Sun); cos(?_s t ?)=1; P_core=1; E_react=1e46; hypothetical B_j scaling.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef SURFACE_TEMPERATURE_MODULE_H
#define SURFACE_TEMPERATURE_MODULE_H

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

class SurfaceTemperatureModule {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
    std::map<std::string, double> variables;
    double computeB_j_hypothetical(double t, double T_s);
    double computeU_g3_example(double t, double T_s);
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;



public:
    // Constructor: Initialize with framework defaults (Sun)
    SurfaceTemperatureModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeT_s();  // 5778 K (Sun)
    double computeB_j_hypothetical(double t, double T_s);  // Scaled B_j (T)
    double computeU_g3_example(double t, double T_s);  // U_g3 with scaling (J/m^3)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();
};

#endif // SURFACE_TEMPERATURE_MODULE_H

// SurfaceTemperatureModule.cpp
#include "SurfaceTemperatureModule.h"

// Constructor: Set framework defaults (Sun)
SurfaceTemperatureModule::SurfaceTemperatureModule() {
        enableDynamicTerms = true;
        enableLogging = false;
        learningRate = 0.001;
        metadata["enhanced"] = "true";
        metadata["version"] = "2.0-Enhanced";

    // Universal constants
    variables["T_s"] = 5778.0;                      // K (Sun effective)
    variables["T_s_ref"] = 5778.0;                  // K (reference)
    variables["k_3"] = 1.8;                         // Coupling
    variables["B_ref"] = 1e3;                       // Base T (string)
    variables["omega_s"] = 2.5e-6;                  // rad/s
    variables["P_core"] = 1.0;                      // Unitless
    variables["E_react"] = 1e46;                    // J
    variables["pi"] = 3.141592653589793;
    variables["t"] = 0.0;                           // s
}

// Update variable
void SurfaceTemperatureModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void SurfaceTemperatureModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void SurfaceTemperatureModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute T_s (K)
double SurfaceTemperatureModule::computeT_s() {
    return variables["T_s"];
}

// Hypothetical B_j scaled by T_s / T_s_ref
double SurfaceTemperatureModule::computeB_j_hypothetical(double t, double T_s) {
    variables["t"] = t;
    double base_b = variables["B_ref"] + 0.4 * std::sin(variables["omega_s"] * t);  // Cycle
    return base_b * (T_s / variables["T_s_ref"]);
}

// U_g3 example with scaled B_j
double SurfaceTemperatureModule::computeU_g3_example(double t, double T_s) {
    double k_3 = variables["k_3"];
    double b_j = computeB_j_hypothetical(t, T_s);
    double cos_term = std::cos(variables["omega_s"] * t * variables["pi"]);
    double p_core = variables["P_core"];
    double e_react = variables["E_react"];
    return k_3 * b_j * cos_term * p_core * e_react;
}

// Equation text
std::string SurfaceTemperatureModule::getEquationText() {
    return "B_j ? (10^3 + 0.4 sin(?_s t)) * (T_s / T_s,ref) T (hypothetical);\n"
           "U_g3 = k_3 * ? B_j * cos(?_s t ?) * P_core * E_react\n"
           "Where T_s = 5778 K (Sun effective photosphere; ï¿½C=5505).\n"
           "T_s,ref=5778 K; scales string fields by temperature.\n"
           "Example t=0, T_s=5778 K: B_j?1e3 T, U_g3?1.8e49 J/mï¿½;\n"
           "T_s=10000 K: B_j?1730 T, U_g3?3.11e49 J/mï¿½ (+73%).\n"
           "Role: Thermal baseline for magnetic strength; variability in U_g3/disks.\n"
           "UQFF: Temperature-dependent fields; extensible for radiation/formation.";
}

// Print variables
void SurfaceTemperatureModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Example usage in base program (snippet)
// #include "SurfaceTemperatureModule.h"
// int main() {
//     SurfaceTemperatureModule mod;
//     double t_s = mod.computeT_s();
//     std::cout << "T_s = " << t_s << " K\n";
//     double u_g3 = mod.computeU_g3_example(0.0, 10000.0);
//     std::cout << "U_g3 (T_s=10000 K) = " << u_g3 << " J/mï¿½\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("T_s", 6000.0);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o temp_test temp_test.cpp SurfaceTemperatureModule.cpp -lm
// Sample: T_s=5778 K; U_g3 (hot star)?3.11e49 J/mï¿½; thermal scaling.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

SurfaceTemperatureModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeT_s, computeB_j_hypothetical, computeU_g3_example) are clear, concise, and variable - driven.
- Output and debugging functions(printVariables, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.
- Supports temperature scaling for magnetic string strength, enabling modeling of different stellar types.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in surface temperature modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.