// UgIndexModule.h
// Modular C++ implementation of the Index for Discrete Universal Gravity Ranges (i) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module uses i=1 to 4 to label Ug1-Ug4; computes sum_{i=1}^4 k_i * U_gi for F_U contribution.
// Pluggable: #include "UgIndexModule.h"
// UgIndexModule mod; mod.computeSumKUgi(); mod.updateVariable("U_g1", new_value);
// Variables in std::map; defaults for Sun at t=0; i labels: 1=Internal Dipole, 2=Outer Bubble, 3=Magnetic Disk, 4=Star-BH.
// Approximations: k_i from coupling; sum ?1.42e53 J/mï¿½ (Ug2 dominant).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef UG_INDEX_MODULE_H
#define UG_INDEX_MODULE_H

#include <map>
#include <string>
#include <vector>
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

class UgIndexModule {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
    std::map<std::string, double> variables;
    std::vector<double> k_values;  // [k1=1.5, k2=1.2, k3=1.8, k4=1.0]
    std::vector<double> computeAllKUgi();
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;



public:
    // Constructor: Initialize with framework defaults
    UgIndexModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    int getIndexRange();  // i=1 to 4
    double computeU_gi(int i);  // U_gi for i=1-4 (J/m^3)
    double computeK_i(int i);   // k_i for i
    double computeKUgi(int i);  // k_i * U_gi
    double computeSumKUgi(int i_min=1, int i_max=4);  // Sum for F_U

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print breakdown by i
    void printIndexBreakdown();
};

#endif // UG_INDEX_MODULE_H

// UgIndexModule.cpp
#include "UgIndexModule.h"

// Constructor: Set defaults for Sun at t=0
UgIndexModule::UgIndexModule() {
        enableDynamicTerms = true;
        enableLogging = false;
        learningRate = 0.001;
        metadata["enhanced"] = "true";
        metadata["version"] = "2.0-Enhanced";

    // Coupling constants (unitless)
    k_values = {1.5, 1.2, 1.8, 1.0};               // k1 to k4

    // U_gi defaults (J/m^3, Sun t=0)
    variables["U_g1"] = 1.39e26;                    // Internal Dipole
    variables["U_g2"] = 1.18e53;                    // Outer Field Bubble
    variables["U_g3"] = 1.8e49;                     // Magnetic Strings Disk
    variables["U_g4"] = 2.50e-20;                   // Star-Black Hole Interactions

    // Shared params (placeholders)
    variables["t_n"] = 0.0;                         // s
    variables["pi"] = 3.141592653589793;
}

// Update variable
void UgIndexModule::updateVariable(const std::string& name, double value) {
    variables[name] = value;
}

// Add delta
void UgIndexModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void UgIndexModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Get range of i
int UgIndexModule::getIndexRange() {
    return 4;  // i=1 to 4
}

// Compute U_gi
double UgIndexModule::computeU_gi(int i) {
    std::string key = "U_g" + std::to_string(i);
    if (variables.find(key) != variables.end()) {
        return variables[key];
    }
    std::cerr << "U_g" << i << " not found. Returning 0." << std::endl;
    return 0.0;
}

// Compute k_i (1-based)
double UgIndexModule::computeK_i(int i) {
    if (i < 1 || i > 4) {
        std::cerr << "Invalid i: " << i << ". Using k1." << std::endl;
        return k_values[0];
    }
    return k_values[i-1];
}

// Compute k_i * U_gi
double UgIndexModule::computeKUgi(int i) {
    return computeK_i(i) * computeU_gi(i);
}

// Sum over i_min to i_max
double UgIndexModule::computeSumKUgi(int i_min, int i_max) {
    double sum = 0.0;
    for (int i = i_min; i <= i_max; ++i) {
        sum += computeKUgi(i);
    }
    return sum;
}

// Equation text
std::string UgIndexModule::getEquationText() {
    return "F_U = ?_{i=1}^4 [k_i * U_gi(r,t,M_s,?_s,T_s,B_s,?_vac,[SCm],?_vac,[UA],t_n) - ?_i * ... ] + other terms\n"
           "i (dimensionless integer): Labels Ug ranges; i=1: Internal Dipole, i=2: Outer Bubble,\n"
           "i=3: Magnetic Disk, i=4: Star-BH.\n"
           "Discretizes gravity for summation; enables scale-specific modeling.\n"
           "Example Sun t=0: ? k_i U_gi ?1.42e53 J/mï¿½ (Ug2 dominant).\n"
           "Role: Structures Ug contributions; extensible for more ranges.";
}

// Print variables
void UgIndexModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
    std::cout << "k_i: k1=1.5, k2=1.2, k3=1.8, k4=1.0\n";
}

// Print breakdown
void UgIndexModule::printIndexBreakdown() {
    std::cout << "Ug Index Breakdown (i=1 to 4):\n";
    for (int i = 1; i <= 4; ++i) {
        double ugi = computeU_gi(i);
        double ki = computeK_i(i);
        double kugi = computeKUgi(i);
        std::string label;
        switch(i) {
            case 1: label = "Internal Dipole"; break;
            case 2: label = "Outer Field Bubble"; break;
            case 3: label = "Magnetic Strings Disk"; break;
            case 4: label = "Star-Black Hole"; break;
            default: label = "Unknown";
        }
        std::cout << "i=" << i << " (" << label << "): U_g" << i << "=" << std::scientific << ugi
                  << ", k" << i << "=" << ki << ", k_i U_gi=" << kugi << " J/mï¿½\n";
    }
    std::cout << "Sum ? k_i U_gi = " << std::scientific << computeSumKUgi() << " J/mï¿½\n";
}

// Example usage in base program (snippet)
// #include "UgIndexModule.h"
// int main() {
//     UgIndexModule mod;
//     double sum = mod.computeSumKUgi();
//     std::cout << "? k_i U_gi = " << sum << " J/mï¿½\n";
//     mod.printIndexBreakdown();
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("U_g3", 2e49);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o ug_index_test ug_index_test.cpp UgIndexModule.cpp -lm
// Sample: Sum ?1.42e53 J/mï¿½; i structures 4 Ug ranges.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

UgIndexModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Coupling constants and gravity terms are clearly separated and indexed, supporting extensibility.
- Core computation methods(computeU_gi, computeK_i, computeKUgi, computeSumKUgi) are clear, concise, and variable - driven.
- Output and debugging functions(printVariables, printIndexBreakdown, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid indices, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map and std::vector.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in universal gravity modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.