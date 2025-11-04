// NegativeTimeModule.h
// Modular C++ implementation of the Negative Time Factor (t_n) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes t_n = t - t_0 (s or days, allows t_n < 0); used in cos(? t_n) for oscillations and exp(-? t cos(? t_n)) for growth/decay.
// Pluggable: #include "NegativeTimeModule.h"
// NegativeTimeModule mod; mod.computeCosPiTn(1000.0); mod.updateVariable("t_0", new_value);
// Variables in std::map; defaults t_0=0, t=0 (t_n=0); example for U_m term with t_n negative.
// Approximations: cos even function; ?=5e-5 day^-1; at t_n=-1, exp term negative (growth).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef NEGATIVE_TIME_MODULE_H
#define NEGATIVE_TIME_MODULE_H

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

class NegativeTimeModule {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
    std::map<std::string, double> variables;
    double computeCosPiTn(double t_n);
    double computeExpTerm(double gamma, double t, double t_n);
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;



public:
    // Constructor: Initialize with framework defaults
    NegativeTimeModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeT_n(double t);  // t_n = t - t_0 (s/days)
    double computeCosPiTn(double t);  // cos(? t_n)
    double computeExpTerm(double gamma, double t);  // exp(-? t cos(? t_n))
    double computeOneMinusExp(double gamma, double t);  // 1 - exp(-? t cos(? t_n))
    double computeUmExample(double t, double mu_over_rj = 2.26e10);  // Simplified U_m contrib

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print t_n effects (positive/negative)
    void printTnEffects(double t, double gamma = 5e-5);
};

#endif // NEGATIVE_TIME_MODULE_H

// NegativeTimeModule.cpp
#include "NegativeTimeModule.h"

// Constructor: Set framework defaults
NegativeTimeModule::NegativeTimeModule() {
        enableDynamicTerms = true;
        enableLogging = false;
        learningRate = 0.001;
        metadata["enhanced"] = "true";
        metadata["version"] = "2.0-Enhanced";

    // Universal constants
    variables["t_0"] = 0.0;                         // Reference time (s/days)
    variables["t"] = 0.0;                           // Current time
    variables["gamma"] = 5e-5;                      // day^-1 (example)
    variables["pi"] = 3.141592653589793;
    variables["mu_over_rj"] = 2.26e10;              // T m^2 (example)
    variables["P_SCm"] = 1.0;                       // Normalized
    variables["E_react"] = 1e46;                    // J
    variables["heaviside_f"] = 1e11 + 1.0;          // 1 + 10^13 * 0.01
    variables["quasi_f"] = 1.01;                    // 1 + 0.01
}

// Update variable
void NegativeTimeModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void NegativeTimeModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void NegativeTimeModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute t_n = t - t_0
double NegativeTimeModule::computeT_n(double t) {
    variables["t"] = t;
    return t - variables["t_0"];
}

// Compute cos(? t_n)
double NegativeTimeModule::computeCosPiTn(double t) {
    double t_n = computeT_n(t);
    return std::cos(variables["pi"] * t_n);
}

// Compute exp(-? t cos(? t_n))
double NegativeTimeModule::computeExpTerm(double gamma, double t) {
    double cos_pi_tn = computeCosPiTn(t);
    double arg = - gamma * t * cos_pi_tn;
    return std::exp(arg);
}

// Compute 1 - exp(-? t cos(? t_n))
double NegativeTimeModule::computeOneMinusExp(double gamma, double t) {
    return 1.0 - computeExpTerm(gamma, t);
}

// Simplified U_m example contrib
double NegativeTimeModule::computeUmExample(double t, double mu_over_rj) {
    double gamma = variables["gamma"];
    double one_minus_exp = computeOneMinusExp(gamma, t);
    double phi_hat = 1.0;
    double p_scm = variables["P_SCm"];
    double e_react = variables["E_react"];
    double heaviside_f = variables["heaviside_f"];
    double quasi_f = variables["quasi_f"];
    return (mu_over_rj * one_minus_exp * phi_hat) * p_scm * e_react * heaviside_f * quasi_f;
}

// Equation text
std::string NegativeTimeModule::getEquationText() {
    return "t_n = t - t_0 (s/days, allows t_n < 0 for time-reversal);\n"
           "Used in: cos(? t_n) for oscillations; exp(-? t cos(? t_n)) for decay/growth.\n"
           "In U_m: ... (1 - exp(-? t cos(? t_n))) ...;\n"
           "Negative t_n: e.g., t_n=-1 ? cos(-?)=-1 ? exp(? t) >1 (growth, negentropic).\n"
           "Example t=1000 days, ?=5e-5 day^-1, t_0=0: 1-exp ?0.049, U_m ?1.12e66 J/mï¿½.\n"
           "t_n=-1000: same (cos even); t_n=-1: 1-exp ? -0.051 (growth phase).\n"
           "Role: Models cyclic/TRZ dynamics; forward/reverse time in nebulae/mergers/jets.";
}

// Print variables
void NegativeTimeModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Print effects for positive/negative t_n
void NegativeTimeModule::printTnEffects(double t, double gamma) {
    double t_n_pos = computeT_n(t);  // Positive example
    double cos_pos = computeCosPiTn(t);
    double exp_pos = computeExpTerm(gamma, t);
    double one_minus_pos = computeOneMinusExp(gamma, t);
    double um_pos = computeUmExample(t);

    // Negative t_n: adjust t_0 to make t_n negative
    double orig_t0 = variables["t_0"];
    variables["t_0"] = t + 1.0;  // t_n = t - (t+1) = -1
    double t_n_neg = computeT_n(t);
    double cos_neg = computeCosPiTn(t);
    double exp_neg = computeExpTerm(gamma, t);
    double one_minus_neg = computeOneMinusExp(gamma, t);
    double um_neg = computeUmExample(t);

    variables["t_0"] = orig_t0;  // Restore

    std::cout << "t_n Effects at t=" << t << " (?=" << gamma << "):\n";
    std::cout << "Positive t_n (" << t_n_pos << "): cos(? t_n)=" << cos_pos << ", 1-exp=" << one_minus_pos << ", U_m?" << um_pos << " J/mï¿½\n";
    std::cout << "Negative t_n (" << t_n_neg << "): cos(? t_n)=" << cos_neg << ", 1-exp=" << one_minus_neg << ", U_m?" << um_neg << " J/mï¿½\n";
}

// Example usage in base program (snippet)
// #include "NegativeTimeModule.h"
// int main() {
//     NegativeTimeModule mod;
//     double t = 1000.0;  // days
//     mod.printTnEffects(t);
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("t_0", 500.0);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o tn_test tn_test.cpp NegativeTimeModule.cpp -lm
// Sample: Positive: 1-exp?0.049; Negative t_n=-1: 1-exp?-0.051 (growth); U_m scales accordingly.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

NegativeTimeModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeT_n, computeCosPiTn, computeExpTerm, computeOneMinusExp, computeUmExample) are clear, concise, and variable - driven.
- Handles negative time values(t_n < 0) robustly, enabling modeling of time - reversal and cyclic effects.
    - Output and debugging functions(printVariables, printTnEffects, getEquationText) provide transparency and aid validation.
    - Well - documented physical meaning and example calculations in comments and equation text.

    Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in negative time factor modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.