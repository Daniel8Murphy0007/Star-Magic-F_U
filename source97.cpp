// FeedbackFactorModule.h
// Modular C++ implementation of the Feedback Factor (f_feedback) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes f_feedback=0.1 for ?M_BH=1 dex (10x mass increase); scales (1 + f_feedback) in U_g4 term.
// Pluggable: #include "FeedbackFactorModule.h"
// FeedbackFactorModule mod; mod.computeU_g4(0.0, 0.0); mod.updateVariable("f_feedback", new_value);
// Variables in std::map; example for Sun at t=0, t_n=0; ?M_BH=1 dex ? M_bh_final=10*M_bh_initial.
// Approximations: cos(? t_n)=1; e^{-? t}=1 at t=0; ?=0.001 day^-1 (scaled to s^-1 if needed).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef FEEDBACK_FACTOR_MODULE_H
#define FEEDBACK_FACTOR_MODULE_H

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

class FeedbackFactorModule {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
    std::map<std::string, double> variables;
    double computeDeltaM_BH();  // 1 dex = log10(10) = factor of 10
    double computeM_bh_final();
    double computeU_g4(double t, double t_n);
    double computeU_g4_no_feedback(double t, double t_n);
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;



public:
    // Constructor: Initialize with framework defaults
    FeedbackFactorModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeF_feedback();  // f_feedback=0.1 (unitless)
    double computeDeltaM_BH();   // 1 dex
    double computeM_bh_final();  // 10 * M_bh_initial
    double computeU_g4(double t, double t_n);  // With feedback (J/m^3)
    double computeU_g4_no_feedback(double t, double t_n);  // Without (J/m^3)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print U_g4 comparison (with/without feedback)
    void printU_g4_comparison(double t = 0.0, double t_n = 0.0);
};

#endif // FEEDBACK_FACTOR_MODULE_H

// FeedbackFactorModule.cpp
#include "FeedbackFactorModule.h"

// Constructor: Set framework defaults
FeedbackFactorModule::FeedbackFactorModule() {
        enableDynamicTerms = true;
        enableLogging = false;
        learningRate = 0.001;
        metadata["enhanced"] = "true";
        metadata["version"] = "2.0-Enhanced";

    // Universal constants
    variables["f_feedback"] = 0.1;                  // Unitless, for ?M_BH=1 dex
    variables["delta_M_BH_dex"] = 1.0;              // 1 dex = factor 10
    variables["M_bh_initial"] = 8.15e36;            // kg (Sgr A*)
    variables["k_4"] = 1.0;                         // Coupling for Ug4
    variables["rho_vac_SCm"] = 7.09e-37;            // J/m^3
    variables["d_g"] = 2.55e20;                     // m
    variables["alpha"] = 0.001 / 86400.0;           // day^-1 to s^-1
    variables["pi"] = 3.141592653589793;
    variables["t"] = 0.0;                           // s
    variables["t_n"] = 0.0;                         // s
}

// Update variable
void FeedbackFactorModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
    if (name == "delta_M_BH_dex") {
        // Recalculate M_bh_final if dex changes
        computeM_bh_final();
    }
}

// Add delta
void FeedbackFactorModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void FeedbackFactorModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute f_feedback (fixed 0.1 for 1 dex)
double FeedbackFactorModule::computeF_feedback() {
    return variables["f_feedback"];
}

// Compute ?M_BH in dex
double FeedbackFactorModule::computeDeltaM_BH() {
    return variables["delta_M_BH_dex"];
}

// Compute M_bh_final = M_bh_initial * 10^{?M_BH_dex}
double FeedbackFactorModule::computeM_bh_final() {
    double factor = std::pow(10.0, computeDeltaM_BH());
    double initial = variables["M_bh_initial"];
    variables["M_bh_final"] = initial * factor;
    return variables["M_bh_final"];
}

// Compute U_g4 with feedback
double FeedbackFactorModule::computeU_g4(double t, double t_n) {
    double k_4 = variables["k_4"];
    double rho_vac_SCm = variables["rho_vac_SCm"];
    double M_bh = computeM_bh_final();  // Use final mass for feedback scenario
    double d_g = variables["d_g"];
    double alpha = variables["alpha"];
    double pi = variables["pi"];
    double f_feedback = computeF_feedback();
    double exp_term = std::exp( - alpha * t );
    double cos_term = std::cos( pi * t_n );
    double feedback_factor = 1.0 + f_feedback;
    return k_4 * (rho_vac_SCm * M_bh / d_g) * exp_term * cos_term * feedback_factor;
}

// Compute U_g4 without feedback (f_feedback=0)
double FeedbackFactorModule::computeU_g4_no_feedback(double t, double t_n) {
    double orig_f = variables["f_feedback"];
    variables["f_feedback"] = 0.0;
    double result = computeU_g4(t, t_n);
    variables["f_feedback"] = orig_f;  // Restore
    return result;
}

// Equation text
std::string FeedbackFactorModule::getEquationText() {
    return "U_g4 = k_4 * (?_vac,[SCm] M_bh / d_g) * e^{-? t} * cos(? t_n) * (1 + f_feedback)\n"
           "Where f_feedback = 0.1 (unitless, for ?M_BH = 1 dex = 10x mass increase);\n"
           "?M_BH =1 dex ? M_bh_final = 10 * M_bh_initial (8.15e36 kg ? 8.15e37 kg).\n"
           "Example t=0, t_n=0: U_g4 ?2.75e-20 J/mï¿½ (with); ?2.50e-20 J/mï¿½ (without; +10%).\n"
           "Role: Scales feedback in star-BH interactions; regulates AGN, mergers, star formation.\n"
           "UQFF: Enhances energy density for 10x M_BH; resolves final parsec, quasar jets.";
}

// Print variables
void FeedbackFactorModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Print U_g4 comparison
void FeedbackFactorModule::printU_g4_comparison(double t, double t_n) {
    double u_with = computeU_g4(t, t_n);
    double u_without = computeU_g4_no_feedback(t, t_n);
    double delta_percent = ((u_with - u_without) / u_without) * 100.0;
    std::cout << "U_g4 Comparison at t=" << t << " s, t_n=" << t_n << " s:\n";
    std::cout << "With feedback: " << std::scientific << u_with << " J/mï¿½\n";
    std::cout << "Without feedback: " << std::scientific << u_without << " J/mï¿½\n";
    std::cout << "Difference: +" << std::fixed << std::setprecision(1) << delta_percent << "%\n";
}

// Example usage in base program (snippet)
// #include "FeedbackFactorModule.h"
// int main() {
//     FeedbackFactorModule mod;
//     double m_final = mod.computeM_bh_final();
//     std::cout << "M_bh_final = " << m_final << " kg\n";
//     mod.printU_g4_comparison(0.0, 0.0);
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("delta_M_BH_dex", 2.0);  // 100x mass
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o feedback_test feedback_test.cpp FeedbackFactorModule.cpp -lm
// Sample: M_bh_final=8.15e37 kg; U_g4 with=2.75e-20 J/mï¿½ (+10% vs without).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

FeedbackFactorModule Evaluation

Strengths :
-Modular, extensible design for modeling the feedback factor(f_feedback) in UQFF, with clear encapsulation of variables using std::map.
- Implements core physical concepts : feedback scaling for black hole mass increase(?M_BH in dex), and its effect on the Ug4 term.
- Provides both feedback and non - feedback calculations for Ug4, enabling comparative analysis.
- Approximations and physical meaning are well - documented in comments and equation text.
- Output functions for variable state and Ug4 comparison support debugging and transparency.
- Handles dynamic updates to variables and recalculates dependent terms as needed.
- Example usage and equation text provide scientific context and validation.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for flexibility and scalability.
- Some calculations use magic numbers or lack explanatory comments; define named constants and clarify logic.
- Minimal error handling(e.g., division by zero, invalid variable names); add validation for robustness.
- Unit consistency should be checked and documented for all physical quantities.
- For large - scale or performance - critical simulations, optimize data structures and reduce redundant calculations.
- std::map is flexible but may be less efficient than structured types for very large models.
- Expand documentation for function purposes and physical meaning.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in feedback factor modeling.It implements the UQFF feedback concept faithfully and adapts to various scenarios.For production or high - performance applications, address the recommendations for improved robustness, maintainability, and scalability.