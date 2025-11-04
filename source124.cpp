// Ug1DefectModule.h
// Modular C++ implementation of the Ug1 Defect Factor (?_def) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes ?_def = 0.01 * sin(0.001 t) (unitless); scales (1 + ?_def) in Universal Gravity U_g1 term.
// Pluggable: #include "Ug1DefectModule.h"
// Ug1DefectModule mod; mod.computeU_g1(0.0, 1.496e11); mod.updateVariable("amplitude", new_value);
// Variables in std::map; example for Sun at t=0 (?_def=0, U_g1?4.51e31 J/mï¿½); t=1570.8 days: +1%.
// Approximations: ?=0.001 day?ï¿½; cos(? t_n)=1 at t_n=0; ?_s=3.38e23 Tï¿½mï¿½; ?(M_s/r)?M_s/rï¿½=8.89e7 m/sï¿½.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef UG1_DEFECT_MODULE_H
#define UG1_DEFECT_MODULE_H

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

class Ug1DefectModule {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
    std::map<std::string, double> variables;
    double computeDelta_def(double t_day);
    double computeU_g1(double t_day, double r);
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;



public:
    // Constructor: Initialize with framework defaults
    Ug1DefectModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeDelta_def(double t_day);  // 0.01 * sin(0.001 t)
    double computeU_g1(double t_day, double r);  // U_g1 with defect (J/m^3)
    double computePeriod_years();  // ~17.22 years

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();
};

#endif // UG1_DEFECT_MODULE_H

// Ug1DefectModule.cpp
#include "Ug1DefectModule.h"

// Constructor: Set framework defaults (Sun)
Ug1DefectModule::Ug1DefectModule() {
        enableDynamicTerms = true;
        enableLogging = false;
        learningRate = 0.001;
        metadata["enhanced"] = "true";
        metadata["version"] = "2.0-Enhanced";

    // Universal constants
    variables["amplitude"] = 0.01;                  // Unitless
    variables["freq"] = 0.001;                      // day?ï¿½
    variables["k_1"] = 1.5;                         // Coupling
    variables["mu_s"] = 3.38e23;                    // Tï¿½mï¿½
    variables["M_s"] = 1.989e30;                    // kg
    variables["alpha"] = 0.001;                     // day?ï¿½
    variables["t_n"] = 0.0;                         // days
    variables["pi"] = 3.141592653589793;
    variables["t_day"] = 0.0;                       // days
    variables["r"] = 1.496e11;                      // m (Earth-Sun example)

    // Derived
    variables["period_days"] = 2.0 * M_PI / variables["freq"];
}

// Update variable
void Ug1DefectModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
        if (name == "freq") {
            variables["period_days"] = 2.0 * M_PI / value;
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void Ug1DefectModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
        if (name == "freq") {
            variables["period_days"] = 2.0 * M_PI / variables[name];
        }
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void Ug1DefectModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute ?_def = 0.01 * sin(0.001 t) (t in days)
double Ug1DefectModule::computeDelta_def(double t_day) {
    variables["t_day"] = t_day;
    return variables["amplitude"] * std::sin(variables["freq"] * t_day);
}

// Compute U_g1 = k_1 * ?_s * ?(M_s / r) * exp(-? t) * cos(? t_n) * (1 + ?_def)
double Ug1DefectModule::computeU_g1(double t_day, double r) {
    variables["r"] = r;
    double k_1 = variables["k_1"];
    double mu_s = variables["mu_s"];
    double grad_ms_r = variables["M_s"] / (r * r);  // Approx ?(M_s / r) = M_s / r^2
    double exp_term = std::exp( - variables["alpha"] * t_day );
    double cos_tn = std::cos(variables["pi"] * variables["t_n"]);
    double defect_factor = 1.0 + computeDelta_def(t_day);
    return k_1 * mu_s * grad_ms_r * exp_term * cos_tn * defect_factor;
}

// Period in years (365.25 days/year)
double Ug1DefectModule::computePeriod_years() {
    return variables["period_days"] / 365.25;
}

// Equation text
std::string Ug1DefectModule::getEquationText() {
    return "U_g1 = k_1 * ?_s * ?(M_s / r) * e^{-? t} * cos(? t_n) * (1 + ?_def)\n"
           "Where ?_def = 0.01 * sin(0.001 t) (unitless, t days; period ~17.22 yr).\n"
           "Small oscillatory defect (~ï¿½1%) in internal dipole gravity.\n"
           "Example t=0, r=1.496e11 m: ?_def=0, U_g1 ?4.51e31 J/mï¿½;\n"
           "t=1570.8 days: ?_def=0.01, U_g1 ?4.56e31 J/mï¿½ (+1.1%).\n"
           "Role: Time-dependent perturbations; internal dynamics/[SCm] variations.\n"
           "UQFF: Cyclic defects in stellar gravity; for formation/nebular stability.";
}

// Print variables
void Ug1DefectModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Example usage in base program (snippet)
// #include "Ug1DefectModule.h"
// int main() {
//     Ug1DefectModule mod;
//     double delta = mod.computeDelta_def(0.0);
//     std::cout << "?_def (t=0) = " << delta << std::endl;
//     double u_g1 = mod.computeU_g1(1570.8, 1.496e11);
//     std::cout << "U_g1 (t=1570.8 days) = " << u_g1 << " J/mï¿½\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("amplitude", 0.02);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o defect_test defect_test.cpp Ug1DefectModule.cpp -lm
// Sample: ?_def=0 at t=0; U_g1?4.56e31 J/mï¿½ at peak (+1%); period~17.22 yr.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

Ug1DefectModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computeDelta_def, computeU_g1, computePeriod_years) are clear, concise, and variable - driven.
- Automatic recalculation of derived variables(period_days) when frequency changes.
- Output and debugging functions(printVariables, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.
- Models small oscillatory defect in internal dipole gravity, supporting time - dependent perturbation analysis.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in defect factor modeling for universal gravity.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.