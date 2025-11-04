// ReciprocationDecayModule.h
// Modular C++ implementation of the Reciprocation Decay Rate (?) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes ?=0.00005 day?ï¿½ (~5.8e-10 s?ï¿½); used in exp(-? t cos(? t_n)) for U_m decay.
// Pluggable: #include "ReciprocationDecayModule.h"
// ReciprocationDecayModule mod; mod.computeOneMinusExp(1000.0, 0.0); mod.updateVariable("gamma_day", new_value);
// Variables in std::map; example for t=1000 days, t_n=0; 1-exp ?0.049.
// Approximations: cos(? t_n)=1; timescale ~55 years; ?_j / r_j=2.26e10 T mï¿½.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef RECIPROCATION_DECAY_MODULE_H
#define RECIPROCATION_DECAY_MODULE_H

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

class ReciprocationDecayModule {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
    std::map<std::string, double> variables;
    double computeGamma_s();  // ? in s?ï¿½
    double computeCosPiTn(double t_n);
    double computeExpTerm(double t_day, double t_n);
    double computeOneMinusExp(double t_day, double t_n);
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;



public:
    // Constructor: Initialize with framework defaults
    ReciprocationDecayModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeGamma_day();  // 0.00005 day?ï¿½
    double computeGamma_s();    // ~5.8e-10 s?ï¿½
    double computeCosPiTn(double t_n);  // cos(? t_n)
    double computeExpTerm(double t_day, double t_n);  // exp(-? t cos(? t_n))
    double computeOneMinusExp(double t_day, double t_n);  // 1 - exp(...)
    double computeUmExample(double t_day, double t_n, double mu_over_rj = 2.26e10);  // Simplified U_m contrib

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print decay effects
    void printDecayEffects(double t_day = 1000.0, double t_n = 0.0);
};

#endif // RECIPROCATION_DECAY_MODULE_H

// ReciprocationDecayModule.cpp
#include "ReciprocationDecayModule.h"

// Constructor: Set framework defaults
ReciprocationDecayModule::ReciprocationDecayModule() {
        enableDynamicTerms = true;
        enableLogging = false;
        learningRate = 0.001;
        metadata["enhanced"] = "true";
        metadata["version"] = "2.0-Enhanced";

    // Universal constants
    variables["gamma_day"] = 0.00005;               // day?ï¿½
    variables["day_to_s"] = 86400.0;                // s/day
    variables["t_n"] = 0.0;                         // days
    variables["t_day"] = 0.0;                       // days
    variables["pi"] = 3.141592653589793;
    variables["mu_over_rj"] = 2.26e10;              // T mï¿½
    variables["P_SCm"] = 1.0;                       // Normalized
    variables["E_react"] = 1e46;                    // J
    variables["heaviside_f"] = 1e11 + 1.0;          // 1 + 10^13 * 0.01
    variables["quasi_f"] = 1.01;                    // 1 + 0.01

    // Derived
    variables["gamma_s"] = computeGamma_s();
}

// Update variable
void ReciprocationDecayModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
        if (name == "gamma_day") {
            variables["gamma_s"] = computeGamma_s();
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void ReciprocationDecayModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        if (name == "gamma_day") {
            variables["gamma_s"] = computeGamma_s();
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void ReciprocationDecayModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute ? (day?ï¿½)
double ReciprocationDecayModule::computeGamma_day() {
    return variables["gamma_day"];
}

// Compute ? in s?ï¿½
double ReciprocationDecayModule::computeGamma_s() {
    return computeGamma_day() / variables["day_to_s"];
}

// Compute cos(? t_n)
double ReciprocationDecayModule::computeCosPiTn(double t_n) {
    variables["t_n"] = t_n;
    return std::cos(variables["pi"] * t_n);
}

// Compute exp(-? t cos(? t_n)) (t in days)
double ReciprocationDecayModule::computeExpTerm(double t_day, double t_n) {
    variables["t_day"] = t_day;
    double cos_pi_tn = computeCosPiTn(t_n);
    double arg = - computeGamma_day() * t_day * cos_pi_tn;
    return std::exp(arg);
}

// Compute 1 - exp(-? t cos(? t_n))
double ReciprocationDecayModule::computeOneMinusExp(double t_day, double t_n) {
    return 1.0 - computeExpTerm(t_day, t_n);
}

// Simplified U_m example (J/mï¿½)
double ReciprocationDecayModule::computeUmExample(double t_day, double t_n, double mu_over_rj) {
    double one_minus_exp = computeOneMinusExp(t_day, t_n);
    double phi_hat = 1.0;
    double p_scm = variables["P_SCm"];
    double e_react = variables["E_react"];
    double heaviside_f = variables["heaviside_f"];
    double quasi_f = variables["quasi_f"];
    return (mu_over_rj * one_minus_exp * phi_hat) * p_scm * e_react * heaviside_f * quasi_f;
}

// Equation text
std::string ReciprocationDecayModule::getEquationText() {
    return "? = 0.00005 day?ï¿½ (~5.8e-10 s?ï¿½; timescale ~55 years);\n"
           "In U_m: ... (1 - exp(-? t cos(? t_n))) ... (t days, reciprocating decay/growth).\n"
           "Negative cos(? t_n): exp(+? t) >1 (growth, negentropic TRZ).\n"
           "Example t=1000 days, t_n=0: 1-exp ?0.049, U_m ?1.12e66 J/mï¿½.\n"
           "UQFF: Slow decay for magnetic strings; cyclic via cos(? t_n) in jets/nebulae/mergers.";
}

// Print variables
void ReciprocationDecayModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Print effects
void ReciprocationDecayModule::printDecayEffects(double t_day, double t_n) {
    double cos_pi = computeCosPiTn(t_n);
    double exp_val = computeExpTerm(t_day, t_n);
    double one_minus = computeOneMinusExp(t_day, t_n);
    double um_ex = computeUmExample(t_day, t_n);
    std::cout << "Decay Effects at t=" << t_day << " days, t_n=" << t_n << ":\n";
    std::cout << "cos(? t_n) = " << cos_pi << "\n";
    std::cout << "exp(-? t cos(? t_n)) = " << exp_val << "\n";
    std::cout << "1 - exp(...) = " << one_minus << "\n";
    std::cout << "U_m example contrib = " << um_ex << " J/mï¿½\n";
}

// Example usage in base program (snippet)
// #include "ReciprocationDecayModule.h"
// int main() {
//     ReciprocationDecayModule mod;
//     double gamma = mod.computeGamma_day();
//     std::cout << "? = " << gamma << " day?ï¿½\n";
//     mod.printDecayEffects(1000.0, 0.0);
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("gamma_day", 0.0001);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o decay_test decay_test.cpp ReciprocationDecayModule.cpp -lm
// Sample: ?=5e-5 day?ï¿½; t=1000 days: 1-exp?0.049; U_m?1.12e66 J/mï¿½.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

ReciprocationDecayModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeGamma_day, computeGamma_s, computeCosPiTn, computeExpTerm, computeOneMinusExp, computeUmExample) are clear, concise, and variable - driven.
- Automatic recalculation of derived variables(gamma_s) when dependencies change.
- Output and debugging functions(printVariables, printDecayEffects, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in reciprocation decay modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.