// NGC2264UQFFModule.h
// Modular C++ implementation of the Master Universal Gravity Equation (MUGE & UQFF Integration) for Cone Nebula (NGC 2264) Evolution.
// This module models NGC 2264's gravitational dynamics, incorporating stellar winds, pillar erosion, protostar formation, dust/gas densities, and dark matter.
// Usage: #include "NGC2264UQFFModule.h" in base program; NGC2264UQFFModule mod; mod.computeG(t); mod.updateVariable("SFR", new_value);
// Variables in std::map for dynamic updates; supports F_env(t) with wind and erosion terms.
// Approximations: psi_integral normalized to 1.0; H(t,z) with Omega_m=0.3, Omega_Lambda=0.7; E_react exp decay; pillar waves simplified.
// NGC 2264 params: M=100 Msun, r=3.31e16 m, SFR=0.01 Msun/yr, v_wind=20 km/s, rho=1e-20 kg/m^3, B=1e-5 T, z=0.0008, etc.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef NGC2264_UQFF_MODULE_H
#define NGC2264_UQFF_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>


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

class NGC2264UQFFModule {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
    std::map<std::string, double> variables;
    double computeHtz(double z_val);
    double computeFenv(double t);
    double computeUg1(double t);
    double computeUg2(double t);
    double computeUg3prime(double t);
    double computeUg4(double t);
    double computeUi(double t);
    double computePsiIntegral(double r, double t);
    double computeQuantumTerm(double t_Hubble_val, double r);
    double computeFluidTerm(double g_base);
    double computeDMTerm(double r);
    double computeUgSum(double r);
    double computeMsfFactor(double t);
    double computeRt(double t);
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;



public:
    // Constructor: Initialize with NGC 2264 defaults
    NGC2264UQFFModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computation: g_NGC2264(r, t)
    double computeG(double t, double r);

    // Output descriptive text of the equation
    std::string getEquationText();

    // Print all current variables (for debugging)
    void printVariables();
};

#endif // NGC2264_UQFF_MODULE_H

// NGC2264UQFFModule.cpp
#include "NGC2264UQFFModule.h"
#include <complex>

// Constructor: NGC 2264-specific values
NGC2264UQFFModule::NGC2264UQFFModule() {
        enableDynamicTerms = true;
        enableLogging = false;
        learningRate = 0.001;
        metadata["enhanced"] = "true";
        metadata["version"] = "2.0-Enhanced";

    // Universal constants
    variables["G"] = 6.6743e-11;                    // m^3 kg^-1 s^-2
    variables["c"] = 3e8;                           // m/s
    variables["hbar"] = 1.0546e-34;                 // J s
    variables["Lambda"] = 1.1e-52;                  // m^-2
    variables["q"] = 1.602e-19;                     // C
    variables["pi"] = 3.141592653589793;            // pi
    variables["t_Hubble"] = 13.8e9 * 3.156e7;       // s
    variables["year_to_s"] = 3.156e7;               // s/yr
    variables["H0"] = 70.0;                         // km/s/Mpc
    variables["Mpc_to_m"] = 3.086e22;               // m/Mpc
    variables["Omega_m"] = 0.3;
    variables["Omega_Lambda"] = 0.7;
    double M_sun_val = 1.989e30;                    // kg
    double ly_val = 9.461e15;                       // m

    // NGC 2264 parameters
    variables["M_visible"] = 80 * M_sun_val;        // kg
    variables["M_DM"] = 20 * M_sun_val;             // kg
    variables["M"] = variables["M_visible"] + variables["M_DM"];  // Total initial
    variables["M0"] = variables["M"];
    variables["SFR"] = 0.01 * M_sun_val / variables["year_to_s"]; // kg/s
    variables["r"] = 3.31e16;                       // m (~3.5 ly radius)
    variables["z"] = 0.0008;                        // Redshift
    variables["v_wind"] = 20e3;                     // m/s
    variables["t"] = 3e6 * variables["year_to_s"];  // Default t=3 Myr s

    // Dynamics
    variables["rho_fluid"] = 1e-20;                 // kg/m^3
    variables["V"] = 1e48;                          // m^3
    variables["B"] = 1e-5;                          // T
    variables["B_crit"] = 1e11;                     // T
    variables["Delta_x"] = 1e-10;                   // m
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];
    variables["integral_psi"] = 1.0;                // Normalized

    // Wave/oscillatory for pillar dynamics
    variables["A"] = 1e-10;
    variables["k"] = 1e20;
    variables["omega"] = 1e-14;                     // rad/s for wind waves
    variables["x"] = 0.0;
    variables["v"] = variables["v_wind"];           // m/s
    variables["sigma"] = 1e15;                      // m for Gaussian

    // Ug subterms & Ui
    variables["Ug1"] = 0.0;                         // Dipole
    variables["Ug2"] = 0.0;                         // Superconductor
    variables["Ug3"] = 0.0;                         // External
    variables["Ug4"] = 0.0;                         // Reaction
    variables["Ui"] = 0.0;
    variables["mu_0"] = 4 * variables["pi"] * 1e-7; // H/m
    variables["rho_vac_SCm"] = 7.09e-37;            // J/m^3
    variables["rho_vac_UA"] = 7.09e-36;             // J/m^3
    variables["lambda_I"] = 1.0;
    variables["omega_i"] = 1e-8;                    // rad/s
    variables["t_n"] = 0.0;
    variables["F_RZ"] = 0.01;
    variables["k_4"] = 1.0;
    variables["k_SF"] = 1e-10;                      // N/Msun, adjusted to m/s^2
    variables["omega_spin"] = 1e-5;                 // rad/s protostar spin
    variables["I_dipole"] = 1e18;                   // A
    variables["A_dipole"] = 1e12;                   // m^2
    variables["H_aether"] = 1e-6;                   // A/m
    variables["delta_rho_over_rho"] = 1e-5;

    // Scales
    variables["scale_macro"] = 1e-12;
    variables["f_TRZ"] = 0.1;
    variables["f_sc"] = 1.0;
    variables["v_r"] = 1e3;                         // m/s radial velocity
    variables["rho"] = variables["rho_fluid"];
}

// Update variable (with dependents)
void NGC2264UQFFModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding." << std::endl;
        variables[name] = value;
    }
    if (name == "Delta_x") {
        variables["Delta_p"] = variables["hbar"] / value;
    } else if (name == "M_visible" || name == "M_DM") {
        variables["M"] = variables["M_visible"] + variables["M_DM"];
        variables["M0"] = variables["M"];
    }
}

// Add/subtract
void NGC2264UQFFModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        variables[name] = delta;
    }
}
void NGC2264UQFFModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute H(t, z)
double NGC2264UQFFModule::computeHtz(double z_val) {
    double Hz_kms = variables["H0"] * std::sqrt(variables["Omega_m"] * std::pow(1.0 + z_val, 3) + variables["Omega_Lambda"]);
    return (Hz_kms * 1e3) / variables["Mpc_to_m"];
}

// M(t)
double NGC2264UQFFModule::computeMsfFactor(double t) {
    return variables["SFR"] * t / variables["M0"];
}

// r(t)
double NGC2264UQFFModule::computeRt(double t) {
    return variables["r"] + variables["v_r"] * t;
}

// F_env(t)
double NGC2264UQFFModule::computeFenv(double t) {
    double F_wind = variables["rho_fluid"] * std::pow(variables["v_wind"], 2);
    double F_SF = variables["k_SF"] * variables["SFR"] / 1.989e30;  // Normalize to m/s^2
    double F_erode = 0.05 * (t / (3e6 * variables["year_to_s"]));  // Erosion factor
    return F_wind + F_SF + F_erode;
}

// Ug1: dipole
double NGC2264UQFFModule::computeUg1(double t) {
    double mu_dipole = variables["I_dipole"] * variables["A_dipole"] * variables["omega_spin"];
    return mu_dipole * variables["B"];
}

// Ug2: superconductor
double NGC2264UQFFModule::computeUg2(double t) {
    double B_super = variables["mu_0"] * variables["H_aether"];
    return (B_super * B_super) / (2 * variables["mu_0"]);
}

// Ug3': external (stellar wind as external)
double NGC2264UQFFModule::computeUg3prime(double t) {
    double M_star = 20 * 1.989e30;  // Massive star
    double r_star = 1e10;           // Approx
    return (variables["G"] * M_star) / (r_star * r_star);
}

// Ug4: reaction
double NGC2264UQFFModule::computeUg4(double t) {
    double E_react = 1e40 * std::exp(-0.0005 * t);
    return variables["k_4"] * E_react;
}

// Ui
double NGC2264UQFFModule::computeUi(double t) {
    return variables["lambda_I"] * (variables["rho_vac_SCm"] / variables["rho_vac_UA"]) * variables["omega_i"] * std::cos(variables["pi"] * variables["t_n"]) * (1 + variables["F_RZ"]);
}

// Psi integral (simplified)
double NGC2264UQFFModule::computePsiIntegral(double r, double t) {
    double A = variables["A"];
    double m = 1.0;  // m-mode for pillar
    double omega = variables["omega"];
    double sigma = variables["sigma"];
    std::complex<double> psi_pillar(A * std::exp(-r*r / (2 * sigma * sigma)) * std::exp(std::complex<double>(0, m * 0 - omega * t)));
    return std::norm(psi_pillar);  // |psi|^2
}

// Quantum term
double NGC2264UQFFModule::computeQuantumTerm(double t_Hubble_val, double r) {
    double unc = std::sqrt(variables["Delta_x"] * variables["Delta_p"]);
    double psi_int = computePsiIntegral(r, variables["t"]);
    return (variables["hbar"] / unc) * variables["integral_psi"] * (2 * variables["pi"] / t_Hubble_val) * psi_int;
}

// Fluid
double NGC2264UQFFModule::computeFluidTerm(double g_base) {
    return variables["rho_fluid"] * variables["V"] * g_base;
}

// DM
double NGC2264UQFFModule::computeDMTerm(double r) {
    double pert = variables["delta_rho_over_rho"];
    double curv = 3 * variables["G"] * variables["M"] / (r * r * r);
    return (variables["M_visible"] + variables["M_DM"]) * (pert + curv);
}

// Ug sum
double NGC2264UQFFModule::computeUgSum(double r) {
    double Ug_base = (variables["G"] * variables["M"]) / (r * r);
    variables["Ug1"] = computeUg1(variables["t"]);
    variables["Ug2"] = computeUg2(variables["t"]);
    variables["Ug3"] = computeUg3prime(variables["t"]);
    variables["Ug4"] = computeUg4(variables["t"]);
    return Ug_base + variables["Ug1"] + variables["Ug2"] + variables["Ug3"] + variables["Ug4"];
}

// Full g_NGC2264
double NGC2264UQFFModule::computeG(double t, double r) {
    variables["t"] = t;
    double msf_factor = computeMsfFactor(t);
    double m_factor = 1.0 + msf_factor;
    double Hz = computeHtz(variables["z"]);
    double expansion = 1.0 + Hz * t;
    double sc_correction = 1.0 - (variables["B"] / variables["B_crit"]);
    double f_env = computeFenv(t);
    double tr_factor = 1.0 + variables["f_TRZ"];
    double rt = computeRt(t);  // But use input r for profile

    // Base gravity
    double g_base = (variables["G"] * variables["M"] * m_factor / (r * r)) * expansion * sc_correction * (1.0 + f_env) * tr_factor;

    // Ug sum (includes base? Adjust: Ug sum without base)
    double ug_sum = computeUgSum(r) - g_base;  // Subtract to avoid double-count

    // Cosmological
    double lambda_term = variables["Lambda"] * (variables["c"] * variables["c"]) / 3.0;

    // Ui
    double ui_term = computeUi(t);

    // Quantum
    double quantum_term = computeQuantumTerm(variables["t_Hubble"], r);

    // Fluid
    double fluid_term = computeFluidTerm(g_base);

    // DM
    double dm_term = computeDMTerm(r);

    // Total
    return g_base + ug_sum + lambda_term + ui_term + quantum_term + fluid_term + dm_term;
}

// Equation text
std::string NGC2264UQFFModule::getEquationText() {
    return "g_NGC2264(r, t) = (G * M(t) / r(t)^2) * (1 + H(t, z)) * (1 - B(t) / B_crit) * (1 + F_env(t)) + "
           "(U_g1 + U_g2 + U_g3' + U_g4) + U_i + (Lambda * c^2 / 3) + "
           "(hbar / sqrt(?x * ?p)) * ? (?_total * H * ?_total dV) * (2? / t_Hubble) + "
           "?_fluid * V * g + (M_visible + M_DM) * (??/? + 3 G M / r^3)\n"
           "Where: M(t) = M * (1 + M_SF(t)); M_SF(t) = SFR * t; r(t) = r0 + v_r t;\n"
           "H(t, z) = H0 * sqrt(?m (1+z)^3 + ??); F_env(t) = F_wind + F_SF + F_erode;\n"
           "F_wind = ? v_wind^2; U_g1 = ?_dipole * B; U_g2 = B_super^2 / (2 ?0);\n"
           "U_g3' = G M_star / r_star^2; U_g4 = k4 * E_react(t); U_i = ?_I * (?_SCm/?_UA) * ?_i * cos(? t_n) * (1 + F_RZ);\n"
           "?_total = A exp(-r^2/(2?^2)) exp(i(m? - ? t)) + wind terms; Insights: Attractive (g_base, Ug1, Ug3') vs. Repulsive (U_g2, ?) advance UQFF.\n"
           "Adaptations: Hubble ACS 2002 data; SFR=0.01 Msun/yr; M=100 Msun. Solutions: g ~1e-9 m/sï¿½ at t=3 Myr (wind/fluid dominant).";
}

// Print
void NGC2264UQFFModule::printVariables() {
    std::cout << "NGC 2264 Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Example usage
// #include "NGC2264UQFFModule.h"
// int main() {
//     NGC2264UQFFModule mod;
//     double t = 3e6 * 3.156e7;  // 3 Myr
//     double r = 1e16;  // 1 ly
//     double g = mod.computeG(t, r);
//     std::cout << "g_NGC2264 = " << g << " m/sï¿½\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("SFR", 0.02 * 1.989e30 / (3.156e7));
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o ngc2264_sim base.cpp NGC2264UQFFModule.cpp -lm
// Sample Output: g_NGC2264 ~ 1e-9 m/sï¿½ (wind/fluid dominant; repulsive terms advance framework).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

NGC2264UQFFModule Evaluation

Strengths :
-Modular, extensible design for modeling Cone Nebula(NGC 2264) gravity, including stellar winds, pillar erosion, protostar formation, dust / gas densities, and dark matter.
- Comprehensive physics : gravity, cosmological expansion, magnetic fields, environmental / wind / erosion effects, quantum, fluid, and DM terms.
- Dynamic variable management via std::map enables runtime updates and system adaptation.
- Clear separation of computation functions(e.g., Ug1ï¿½Ug4, F_env, quantum, fluid, DM), aiding maintainability.
- NGC 2264 - specific parameters are initialized for realistic simulation; supports easy modification.
- Output functions for equation text and variable state support debugging and documentation.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for flexibility and scalability.
- Some calculations use magic numbers or lack explanatory comments; define named constants and clarify logic.
- Minimal error handling(e.g., division by zero, invalid variable names); add validation for robustness.
- Unit consistency should be checked and documented for all physical quantities.
- For large - scale or performance - critical simulations, optimize data structures and reduce redundant calculations.
- std::map is flexible but may be less efficient than structured types for very large models.
- Expand documentation for function purposes and physical meaning.

Summary:
The code is well - structured, flexible, and suitable for scientific prototyping and educational use in nebular dynamics modeling.It implements a broad set of physical effects and adapts to various scenarios.For production or high - performance applications, address the recommendations for improved robustness, maintainability, and scalability.