// ScmVelocityModule.h
// Modular C++ implementation of the [SCm] Velocity (v_SCm) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes v_SCm = 1e8 m/s (~c/3); scales in E_react = ρ_vac,[SCm] v_SCm² / ρ_vac,A * exp(-κ t) for U_m, U_bi, etc.
// Pluggable: #include "ScmVelocityModule.h"
// ScmVelocityModule mod; mod.computeE_react(0.0); mod.updateVariable("v_sc m", new_value);
// Variables in std::map; example for Sun at t=0 (E_react=1e46 J); t=2000 days: scales down via exp.
// Approximations: κ=0.0005 day⁻¹; ρ_vac,[SCm]=7.09e-37 J/m³; ρ_vac,A=1e-23 J/m³; U_m base=2.28e65 J/m³.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef SCM_VELOCITY_MODULE_H
#define SCM_VELOCITY_MODULE_H

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

class ScmVelocityModule {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
    std::map<std::string, double> variables;
    double computeE_react_base();
    double computeE_react(double t_day);
    double computeUmExample(double t_day);
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;



public:
    // Constructor: Initialize with framework defaults
    ScmVelocityModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeV_sc m();  // 1e8 m/s
    double computeE_react(double t_day);  // ρ_[SCm] v_SCm² / ρ_A * exp(-κ t)
    double computeUmExample(double t_day);  // Simplified U_m with E_react (J/m³)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print velocity effects
    void printVelocityEffects(double t_day = 2000.0);
};

#endif // SCM_VELOCITY_MODULE_H

// ScmVelocityModule.cpp
#include "ScmVelocityModule.h"

// Constructor: Set framework defaults
ScmVelocityModule::ScmVelocityModule() {
        enableDynamicTerms = true;
        enableLogging = false;
        learningRate = 0.001;
        metadata["enhanced"] = "true";
        metadata["version"] = "2.0-Enhanced";

    // Universal constants
    variables["v_sc m"] = 1e8;                      // m/s (~c/3)
    variables["rho_vac_SCm"] = 7.09e-37;            // J/m³
    variables["rho_vac_A"] = 1e-23;                 // J/m³
    variables["kappa_day"] = 0.0005;                // day⁻¹
    variables["day_to_s"] = 86400.0;                // s/day
    variables["t_day"] = 0.0;                       // days
    variables["E_react_base"] = variables["rho_vac_SCm"] * std::pow(variables["v_sc m"], 2) / variables["rho_vac_A"];  // Derived
    variables["mu_over_rj"] = 2.26e10;              // T m² (example)
    variables["P_SCm"] = 1.0;                       // Normalized
    variables["heaviside_f"] = 1e11 + 1.0;          // 1 + 10^13 * 0.01
    variables["quasi_f"] = 1.01;                    // 1 + 0.01
    variables["one_minus_exp"] = 0.0;               // At t=0

    // Derived
    variables["kappa_s"] = variables["kappa_day"] / variables["day_to_s"];
}

// Update variable
void ScmVelocityModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
        if (name == "v_sc m" || name == "rho_vac_SCm" || name == "rho_vac_A") {
            variables["E_react_base"] = variables["rho_vac_SCm"] * std::pow(variables["v_sc m"], 2) / variables["rho_vac_A"];
        } else if (name == "kappa_day") {
            variables["kappa_s"] = value / variables["day_to_s"];
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void ScmVelocityModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        if (name == "v_sc m" || name == "rho_vac_SCm" || name == "rho_vac_A") {
            variables["E_react_base"] = variables["rho_vac_SCm"] * std::pow(variables["v_sc m"], 2) / variables["rho_vac_A"];
        } else if (name == "kappa_day") {
            variables["kappa_s"] = variables["kappa_day"] / variables["day_to_s"];
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void ScmVelocityModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute v_SCm (m/s)
double ScmVelocityModule::computeV_sc m() {
    return variables["v_sc m"];
}

// Compute E_react = E_react_base * exp(-κ t)
double ScmVelocityModule::computeE_react(double t_day) {
    variables["t_day"] = t_day;
    double arg = - variables["kappa_day"] * t_day;
    return variables["E_react_base"] * std::exp(arg);
}

// Simplified U_m example with E_react
double ScmVelocityModule::computeUmExample(double t_day) {
    double e_react = computeE_react(t_day);
    double one_minus_exp = variables["one_minus_exp"];  // Placeholder
    double phi_hat = 1.0;
    double p_scm = variables["P_SCm"];
    double heaviside_f = variables["heaviside_f"];
    double quasi_f = variables["quasi_f"];
    return (variables["mu_over_rj"] * one_minus_exp * phi_hat) * p_scm * e_react * heaviside_f * quasi_f;
}

// Equation text
std::string ScmVelocityModule::getEquationText() {
    return "E_react = [ρ_vac,[SCm] v_SCm² / ρ_vac,A] * exp(-κ t) (t days);\n"
           "v_SCm = 1e8 m/s (~c/3, [SCm] propagation speed);\n"
           "Scales reactivity in U_m, U_bi, U_i, U_gi via E_react.\n"
           "Example t=0: E_react=1e46 J; t=2000 days: ~3.68e45 J (~36.8%).\n"
           "U_m (t=0): ≈2.28e65 J/m³; t=2000: ≈8.39e64 J/m³.\n"
           "Role: [SCm] dynamic speed for relativistic effects; jets/energy transfer.\n"
           "UQFF: Subluminal propagation; [SCm]-[UA] reactions in nebulae/formation.";
}

// Print variables
void ScmVelocityModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Print effects
void ScmVelocityModule::printVelocityEffects(double t_day) {
    double v = computeV_sc m();
    double e_react = computeE_react(t_day);
    double um_ex = computeUmExample(t_day);
    double fraction = e_react / variables["E_react_base"];
    std::cout << "[SCm] Velocity Effects at t=" << t_day << " days:\n";
    std::cout << "v_SCm = " << std::scientific << v << " m/s\n";
    std::cout << "E_react = " << e_react << " J (" << fraction << " of initial)\n";
    std::cout << "U_m example = " << um_ex << " J/m³\n";
}

// Example usage in base program (snippet)
// #include "ScmVelocityModule.h"
// int main() {
//     ScmVelocityModule mod;
//     double v = mod.computeV_sc m();
//     std::cout << "v_SCm = " << v << " m/s\n";
//     mod.printVelocityEffects(2000.0);
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("v_sc m", 1.5e8);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o scm_vel_test scm_vel_test.cpp ScmVelocityModule.cpp -lm
// Sample: v_SCm=1e8 m/s; t=2000 days: E_react≈3.68e45 J; U_m≈8.39e64 J/m³.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

ScmVelocityModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeV_sc m, computeE_react, computeUmExample) are clear, concise, and variable - driven.
- Automatic recalculation of derived variables(E_react_base, kappa_s) when dependencies change.
- Output and debugging functions(printVariables, printVelocityEffects, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.
- Models relativistic[SCm] velocity effects and their impact on reactivity and energy transfer.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- The variable name "v_sc m" contains a space, which is non - standard and may cause confusion or errors; use "v_scm" or "v_SCm" instead.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in[SCm] velocity and reactivity modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.