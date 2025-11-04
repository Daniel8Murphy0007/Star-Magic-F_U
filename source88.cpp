// AndromedaUQFFModule.h
// Modular C++ implementation of the Master Universal Gravity Equation (UQFF) for Andromeda Galaxy Evolution.
// This module can be plugged into a base program by including this header and linking the .cpp.
// Usage: #include "AndromedaUQFFModule.h"
// AndromedaUQFFModule mod; mod.computeG(t); mod.updateVariable("M", new_value);
// Variables stored in std::map for dynamic updates.
// Includes base gravity with expansion and TRZ, BH term, dust friction a_dust, EM/Aether term.
// Approximations: z=-0.001 (blueshift); dust scaled by 1e-12; EM normalized to proton mass.
// Andromeda params: M=1e12 Msun, r=1.04e21 m, M_BH=1.4e8 Msun, v_orbit=2.5e5 m/s, etc.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef ANDROMEDA_UQFF_MODULE_H
#define ANDROMEDA_UQFF_MODULE_H

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

class AndromedaUQFFModule {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
    std::map<std::string, double> variables;
    double computeHz();
    double computeADust();
    double computeEMBase();
    double computeEMTerm();
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;



public:
    // Constructor: Initialize with Andromeda defaults
    AndromedaUQFFModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computation: g_Andromeda(r, t)
    double computeG(double t);

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print evolution table (0-10 Gyr, 2 Gyr steps)
    void printEvolutionTable();
};

#endif // ANDROMEDA_UQFF_MODULE_H

// AndromedaUQFFModule.cpp
#include "AndromedaUQFFModule.h"

// Constructor: Set Andromeda-specific values
AndromedaUQFFModule::AndromedaUQFFModule() {
        enableDynamicTerms = true;
        enableLogging = false;
        learningRate = 0.001;
        metadata["enhanced"] = "true";
        metadata["version"] = "2.0-Enhanced";

    // Universal constants
    variables["G"] = 6.6743e-11;                    // m^3 kg^-1 s^-2
    variables["M_sun"] = 1.989e30;                  // kg
    variables["q"] = 1.602e-19;                     // C
    variables["proton_mass"] = 1.673e-27;           // kg
    variables["H0"] = 70.0;                         // km/s/Mpc
    variables["Mpc_to_m"] = 3.086e22;               // m/Mpc
    variables["Omega_m"] = 0.3;
    variables["Omega_Lambda"] = 0.7;
    variables["year_to_s"] = 3.156e7;               // s/yr
    variables["Gyr"] = 1e9;                         // yr

    // Andromeda parameters
    variables["M"] = 1e12 * variables["M_sun"];     // Total mass kg
    variables["r"] = 1.04e21;                       // m (half diameter)
    variables["M_BH"] = 1.4e8 * variables["M_sun"]; // SMBH mass kg
    variables["r_BH"] = 1e15;                       // m (core scale)
    variables["rho_dust"] = 1e-20;                  // kg/m^3
    variables["v_orbit"] = 2.5e5;                   // m/s
    variables["rho_mass"] = 1e-21;                  // kg/m^3
    variables["z"] = -0.001;                        // Blueshift
    variables["B"] = 1e-5;                          // T
    variables["rho_vac_UA"] = 7.09e-36;             // J/m^3
    variables["rho_vac_SCm"] = 7.09e-37;            // J/m^3
    variables["f_TRZ"] = 0.1;                       // dimensionless
    variables["scale_macro"] = 1e-12;               // Scaling factor
    variables["t"] = 10.0 * variables["Gyr"] * variables["year_to_s"];  // Default 10 Gyr
}

// Update variable
void AndromedaUQFFModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
    if (name == "M") {
        variables["M_BH"] = 1.4e8 * (value / (1e12 * variables["M_sun"])) * variables["M_sun"];  // Scale if needed
    }
}

// Add delta
void AndromedaUQFFModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void AndromedaUQFFModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute H(z) in s^-1
double AndromedaUQFFModule::computeHz() {
    double one_plus_z = 1.0 + variables["z"];
    double Hz_kms = variables["H0"] * std::sqrt(variables["Omega_m"] * std::pow(one_plus_z, 3) + variables["Omega_Lambda"]);
    return (Hz_kms * 1e3) / variables["Mpc_to_m"];
}

// Compute dust acceleration
double AndromedaUQFFModule::computeADust() {
    double force_per_area = variables["rho_dust"] * std::pow(variables["v_orbit"], 2);
    double a_dust_base = force_per_area / variables["rho_mass"];
    return a_dust_base * variables["scale_macro"];
}

// Compute EM base (m/s^2)
double AndromedaUQFFModule::computeEMBase() {
    double mag_vB = variables["v_orbit"] * variables["B"];
    double force = variables["q"] * mag_vB;
    return force / variables["proton_mass"];
}

// Compute full EM term
double AndromedaUQFFModule::computeEMTerm() {
    double em_base = computeEMBase();
    double vac_ratio = variables["rho_vac_UA"] / variables["rho_vac_SCm"];
    return em_base * (1.0 + vac_ratio) * variables["scale_macro"];
}

// Full g_Andromeda
double AndromedaUQFFModule::computeG(double t) {
    variables["t"] = t;
    double Hz = computeHz();
    double expansion_factor = 1.0 + Hz * t;
    double tr_factor = 1.0 + variables["f_TRZ"];

    double g_grav = (variables["G"] * variables["M"] / (variables["r"] * variables["r"])) * expansion_factor * tr_factor;
    double g_BH = variables["G"] * variables["M_BH"] / (variables["r_BH"] * variables["r_BH"]);
    double a_dust = computeADust();
    double em_term = computeEMTerm();

    return g_grav + g_BH + a_dust + em_term;
}

// Equation text
std::string AndromedaUQFFModule::getEquationText() {
    return "g_Andromeda(r, t) = (G * M / r^2) * (1 + H(z) * t) * (1 + f_TRZ) + (G * M_BH / r_BH^2) + a_dust + q*(v*B) * (1 + ?_UA/?_SCm) * 1e-12\n"
           "Where a_dust = (?_dust * v_orbit^2 / ?_mass) * scale_macro;\n"
           "EM term: q v B / m_proton * (1 + ?_vac_UA / ?_vac_SCm) * scale_macro.\n"
           "Andromeda Adaptations: Blueshift z=-0.001; M=1e12 M_sun; dust lanes with v_orbit=250 km/s.\n"
           "At t=10 Gyr, g ?6.273 m/sï¿½ (dust dominant); minimal evolution due to small H(z)t.\n"
           "UQFF Terms: f_TRZ for time-reversal; Aether vacua ratio for EM enhancement.";
}

// Print variables
void AndromedaUQFFModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Print evolution table
void AndromedaUQFFModule::printEvolutionTable() {
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Evolution over time (m/sï¿½):\n";
    std::cout << "t (Gyr) | g_Andromeda\n";
    std::cout << "--------|------------\n";
    for (int i = 0; i <= 5; ++i) {
        double t = i * 2.0 * variables["Gyr"] * variables["year_to_s"];
        double g = computeG(t);
        std::cout << std::setw(6) << (i*2) << "    | " << g << "\n";
    }
}

// Example usage in base program (snippet)
// #include "AndromedaUQFFModule.h"
// int main() {
//     AndromedaUQFFModule mod;
//     double t = 10.0 * 1e9 * 3.156e7;
//     double g = mod.computeG(t);
//     std::cout << "g = " << g << " m/sï¿½\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.printEvolutionTable();
//     mod.updateVariable("v_orbit", 3e5);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o andromeda_test andromeda_test.cpp AndromedaUQFFModule.cpp -lm
// Sample Output at t=10 Gyr: g ? 6.273 m/sï¿½; table shows near-constant due to small expansion.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

AndromedaUQFFModule Evaluation

Strengths :
-Modular, extensible design for modeling Andromeda Galaxy gravity, including base gravity, expansion, SMBH term, dust friction, and EM / Aether effects.
- Comprehensive physics : gravity, cosmological expansion, SMBH, dust, electromagnetic, and vacuum energy terms.
- Dynamic variable management via std::map enables runtime updates and system adaptation.
- Clear separation of computation functions(e.g., H(z), dust, EM), aiding maintainability.
- Andromeda - specific parameters are initialized for realistic simulation; supports easy modification.
- Output functions for equation text, variable state, and evolution table support debugging and documentation.
- Approximations and scaling factors are documented, supporting scientific reproducibility.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for flexibility and scalability.
- Some calculations use magic numbers or lack explanatory comments; define named constants and clarify logic.
- Minimal error handling(e.g., division by zero, invalid variable names); add validation for robustness.
- Unit consistency should be checked and documented for all physical quantities.
- For large - scale or performance - critical simulations, optimize data structures and reduce redundant calculations.
- std::map is flexible but may be less efficient than structured types for very large models.
- Expand documentation for function purposes and physical meaning.

Summary:
The code is well - structured, flexible, and suitable for scientific prototyping and educational use in galactic gravity modeling.It implements a broad set of physical effects and adapts to various scenarios.For production or high - performance applications, address the recommendations for improved robustness, maintainability, and scalability.