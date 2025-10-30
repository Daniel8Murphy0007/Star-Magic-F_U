// M51UQFFModule.h
// Modular C++ implementation of the Master Universal Gravity Equation (MUGE & UQFF Integration) for Whirlpool Galaxy (M51) Evolution.
// This module models M51's gravitational dynamics, incorporating interaction with NGC 5195, star formation, black hole torus/jets, spiral arm density waves, and dark matter.
// Usage: #include "M51UQFFModule.h" in base program; M51UQFFModule mod; mod.computeG(t); mod.updateVariable("SFR", new_value);
// Variables in std::map for dynamic updates; supports F_env(t) with tidal and SF terms.
// Approximations: psi_integral normalized to 1.0; H(t,z) with Omega_m=0.3, Omega_Lambda=0.7; E_react exp decay; density waves simplified.
// M51 params: M=1.6e11 Msun, r=23.58 kpc, SFR=1 Msun/yr, M_BH=1e6 Msun, M_NGC5195=1e10 Msun, d=50 kpc, B=1e-5 T, z=0.002, etc.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef M51_UQFF_MODULE_H
#define M51_UQFF_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>

class M51UQFFModule {
private:
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

public:
    // Constructor: Initialize with M51 defaults
    M51UQFFModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computation: g_M51(r, t)
    double computeG(double t, double r);

    // Output descriptive text of the equation
    std::string getEquationText();

    // Print all current variables (for debugging)
    void printVariables();
};

#endif // M51_UQFF_MODULE_H

// M51UQFFModule.cpp
#include "M51UQFFModule.h"
#include <complex>

// Constructor: M51-specific values
M51UQFFModule::M51UQFFModule() {
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
    double kpc_val = 3.086e19;                      // m

    // M51 parameters
    variables["M_visible"] = 1.2e11 * M_sun_val;    // kg
    variables["M_DM"] = 4e10 * M_sun_val;           // kg
    variables["M"] = variables["M_visible"] + variables["M_DM"];  // Total initial
    variables["M0"] = variables["M"];
    variables["SFR"] = 1 * M_sun_val / variables["year_to_s"];    // kg/s
    variables["r"] = 23.58e3 * kpc_val;             // m
    variables["z"] = 0.002;                         // Redshift
    variables["M_NGC5195"] = 1e10 * M_sun_val;      // kg
    variables["d_NGC5195"] = 50e3 * kpc_val;        // m
    variables["M_BH"] = 1e6 * M_sun_val;            // kg (central BH)
    variables["t"] = 5e8 * variables["year_to_s"];  // Default t=500 Myr s

    // Dynamics
    variables["rho_fluid"] = 1e-20;                 // kg/m^3
    variables["V"] = 1e50;                          // m^3
    variables["B"] = 1e-5;                          // T
    variables["B_crit"] = 1e11;                     // T
    variables["Delta_x"] = 1e-10;                   // m
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];
    variables["integral_psi"] = 1.0;                // Normalized

    // Wave/oscillatory for spiral arms
    variables["A"] = 1e-10;
    variables["k"] = 1e20;
    variables["omega"] = 1e-15;                     // rad/s for density waves
    variables["x"] = 0.0;
    variables["v"] = 1e3;                           // m/s
    variables["sigma"] = 1e3 * kpc_val;             // m for Gaussian

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
    variables["omega_spin"] = 1e-4;                 // rad/s BH spin
    variables["I_dipole"] = 1e20;                   // A
    variables["A_dipole"] = 1e15;                   // m^2
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
void M51UQFFModule::updateVariable(const std::string& name, double value) {
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
    } else if (name == "SFR") {
        // Adjust units if needed
    }
}

// Add/subtract
void M51UQFFModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        variables[name] = delta;
    }
}
void M51UQFFModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute H(t, z)
double M51UQFFModule::computeHtz(double z_val) {
    double Hz_kms = variables["H0"] * std::sqrt(variables["Omega_m"] * std::pow(1.0 + z_val, 3) + variables["Omega_Lambda"]);
    return (Hz_kms * 1e3) / variables["Mpc_to_m"];
}

// M(t)
double M51UQFFModule::computeMsfFactor(double t) {
    return variables["SFR"] * t / variables["M0"];
}

// r(t)
double M51UQFFModule::computeRt(double t) {
    return variables["r"] + variables["v_r"] * t;
}

// F_env(t)
double M51UQFFModule::computeFenv(double t) {
    double F_tidal = (variables["G"] * variables["M_NGC5195"]) / (variables["d_NGC5195"] * variables["d_NGC5195"]);
    double F_SF = variables["k_SF"] * variables["SFR"] / 1.989e30;  // Normalize to m/s^2
    return F_tidal + F_SF;
}

// Ug1: dipole
double M51UQFFModule::computeUg1(double t) {
    double mu_dipole = variables["I_dipole"] * variables["A_dipole"] * variables["omega_spin"];
    return mu_dipole * variables["B"];
}

// Ug2: superconductor
double M51UQFFModule::computeUg2(double t) {
    double B_super = variables["mu_0"] * variables["H_aether"];
    return (B_super * B_super) / (2 * variables["mu_0"]);
}

// Ug3': external
double M51UQFFModule::computeUg3prime(double t) {
    return (variables["G"] * variables["M_NGC5195"]) / (variables["d_NGC5195"] * variables["d_NGC5195"]);
}

// Ug4: reaction
double M51UQFFModule::computeUg4(double t) {
    double E_react = 1e46 * std::exp(-0.0005 * t);
    return variables["k_4"] * E_react;
}

// Ui
double M51UQFFModule::computeUi(double t) {
    return variables["lambda_I"] * (variables["rho_vac_SCm"] / variables["rho_vac_UA"]) * variables["omega_i"] * std::cos(variables["pi"] * variables["t_n"]) * (1 + variables["F_RZ"]);
}

// Psi integral (simplified)
double M51UQFFModule::computePsiIntegral(double r, double t) {
    double A = variables["A"];
    double m = 2.0;
    double omega = variables["omega"];
    double sigma = variables["sigma"];
    std::complex<double> psi_spiral(A * std::exp(-r*r / (2 * sigma * sigma)) * std::exp(std::complex<double>(0, m * 0 - omega * t)));
    return std::norm(psi_spiral);  // |psi|^2
}

// Quantum term
double M51UQFFModule::computeQuantumTerm(double t_Hubble_val, double r) {
    double unc = std::sqrt(variables["Delta_x"] * variables["Delta_p"]);
    double psi_int = computePsiIntegral(r, variables["t"]);
    return (variables["hbar"] / unc) * variables["integral_psi"] * (2 * variables["pi"] / t_Hubble_val) * psi_int;
}

// Fluid
double M51UQFFModule::computeFluidTerm(double g_base) {
    return variables["rho_fluid"] * variables["V"] * g_base;
}

// DM
double M51UQFFModule::computeDMTerm(double r) {
    double pert = variables["delta_rho_over_rho"];
    double curv = 3 * variables["G"] * variables["M"] / (r * r * r);
    return (variables["M_visible"] + variables["M_DM"]) * (pert + curv);
}

// Ug sum
double M51UQFFModule::computeUgSum(double r) {
    double Ug_base = (variables["G"] * variables["M"]) / (r * r);
    variables["Ug1"] = computeUg1(variables["t"]);
    variables["Ug2"] = computeUg2(variables["t"]);
    variables["Ug3"] = computeUg3prime(variables["t"]);
    variables["Ug4"] = computeUg4(variables["t"]);
    return Ug_base + variables["Ug1"] + variables["Ug2"] + variables["Ug3"] + variables["Ug4"];
}

// Full g_M51
double M51UQFFModule::computeG(double t, double r) {
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
std::string M51UQFFModule::getEquationText() {
    return "g_M51(r, t) = (G * M(t) / r(t)^2) * (1 + H(t, z)) * (1 - B(t) / B_crit) * (1 + F_env(t)) + "
           "(U_g1 + U_g2 + U_g3' + U_g4) + U_i + (Lambda * c^2 / 3) + "
           "(hbar / sqrt(?x * ?p)) * ? (?_total * H * ?_total dV) * (2? / t_Hubble) + "
           "?_fluid * V * g + (M_visible + M_DM) * (??/? + 3 G M / r^3)\n"
           "Where: M(t) = M * (1 + M_SF(t)); M_SF(t) = SFR * t; r(t) = r0 + v_r t;\n"
           "H(t, z) = H0 * sqrt(?m (1+z)^3 + ??); F_env(t) = F_tidal + F_SF;\n"
           "F_tidal = G M_NGC5195 / d^2; U_g1 = ?_dipole * B; U_g2 = B_super^2 / (2 ?0);\n"
           "U_g3' = G M_NGC5195 / d^2; U_g4 = k4 * E_react(t); U_i = ?_I * (?_SCm/?_UA) * ?_i * cos(? t_n) * (1 + F_RZ);\n"
           "?_total = A exp(-r^2/(2?^2)) exp(i(m? - ? t)) + BH terms; Insights: Attractive (g_base, Ug1, Ug3') vs. Repulsive (U_g1, U_g2, ?) components advance UQFF.\n"
           "Adaptations: Hubble ACS/WFPC2 data; SFR=1 Msun/yr; M=1.6e11 Msun. Solutions: g ~3e36 m/s� at t=500 Myr (DM dominant).";
}

// Print
void M51UQFFModule::printVariables() {
    std::cout << "M51 Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Example usage
// #include "M51UQFFModule.h"
// int main() {
//     M51UQFFModule mod;
//     double t = 5e8 * 3.156e7;  // 500 Myr
//     double r = 10e3 * 3.086e19;  // 10 kpc
//     double g = mod.computeG(t, r);
//     std::cout << "g_M51 = " << g << " m/s�\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("SFR", 2 * mod.variables["SFR"]);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o m51_sim base.cpp M51UQFFModule.cpp -lm
// Sample Output: g_M51 ~ 3e36 m/s� (DM/fluid dominant; repulsive terms advance framework).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

M51UQFFModule Evaluation

Strengths :
-Modular, extensible design for modeling M51 galaxy gravity, including interaction with NGC 5195, star formation, black hole, spiral arms, and dark matter.
- Comprehensive physics : gravity, cosmological expansion, magnetic fields, environmental / tidal effects, quantum, fluid, and DM terms.
- Dynamic variable management via std::map enables runtime updates and system adaptation.
- Clear separation of computation functions(e.g., Ug1�Ug4, F_env, quantum, fluid, DM), aiding maintainability.
- M51 - specific parameters are initialized for realistic simulation; supports easy modification.
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
The code is well - structured, flexible, and suitable for scientific prototyping and educational use in galactic dynamics modeling.It implements a broad set of physical effects and adapts to various scenarios.For production or high - performance applications, address the recommendations for improved robustness, maintainability, and scalability.