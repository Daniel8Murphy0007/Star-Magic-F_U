
// UQFFBuoyancyModule.h
// Modular C++ implementation of the full Master Unified Field Equation (F_U_Bi_i & UQFF Integration) for Buoyancy Equations across M104, NGC 4839, Chandra and Webb, NGC 346, NGC 1672.
// This module can be plugged into a base program (e.g., 'uqff_buoyancy_sim.cpp') by including this header and linking the .cpp.
// Usage in base: #include "UQFFBuoyancyModule.h"
// UQFFBuoyancyModule mod; mod.computeFBi(system); mod.updateVariable("F_rel", {new_real, new_imag});
// All variables are stored in a std::map for dynamic addition/subtraction/update, using complex<double> for real/imaginary components.
// Nothing is negligible: Includes all terms - base force, momentum, gravity, vacuum stability, LENR resonance, activation, directed energy, magnetic resonance, neutron, relativistic, neutrino, Sweet vac, Kozima drop.
// Associated text: Outputs descriptive equation string via getEquationText().
// Approximations: Integral approximated as integrand * x2 (quadratic root); imag parts small and not fully scaled; LENR dominant due to low ω_0; x2 from quadratic solver approx; F_rel from 1998 LEP.
// Multi-system params: M104 M=7.17e41 kg r=9.46e20 m; NGC 4839 M=1.46e45 kg r=3.09e22 m; Chandra/Webb (composite M=1e30 kg r=1e10 m); NGC 346 M=1e36 kg r=2.36e17 m; NGC 1672 M=3.978e40 kg r=4.40e20 m.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 22, 2025.

#ifndef UQFF_BUOYANCY_MODULE_H
#define UQFF_BUOYANCY_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>

using cdouble = std::complex<double>;

class UQFFBuoyancyModule {
private:
    std::map<std::string, cdouble> variables;
    cdouble computeIntegrand(double t, const std::string& system);
    cdouble computeDPM_resonance(const std::string& system);
    cdouble computeX2(const std::string& system);
    cdouble computeQuadraticRoot(cdouble a, cdouble b, cdouble c);
    cdouble computeLENRTerm(const std::string& system);
    double computeG(double t, const std::string& system);
    cdouble computeQ_wave(double t, const std::string& system);
    cdouble computeUb1(const std::string& system);
    cdouble computeUi(double t, const std::string& system);
    void setSystemParams(const std::string& system);

public:
    // Constructor: Initialize all variables with multi-system defaults
    UQFFBuoyancyModule();

    // Dynamic variable operations (complex)
    void updateVariable(const std::string& name, cdouble value);
    void addToVariable(const std::string& name, cdouble delta);
    void subtractFromVariable(const std::string& name, cdouble delta);

    // Core computation: Full F_U_Bi_i(r, t) for system (approx integral)
    cdouble computeFBi(const std::string& system, double t);

    // Sub-equations
    cdouble computeCompressed(const std::string& system, double t);  // Integrand
    cdouble computeResonant(const std::string& system);
    cdouble computeBuoyancy(const std::string& system);
    cdouble computeSuperconductive(const std::string& system, double t);
    double computeCompressedG(const std::string& system, double t);  // g(r,t)

    // Output descriptive text of the equation
    std::string getEquationText(const std::string& system);

    // Print all current variables (for debugging/updates)
    void printVariables();
};

#endif // UQFF_BUOYANCY_MODULE_H
// SurfaceMagneticFieldModule.cpp
#include "SurfaceMagneticFieldModule.h"
// Compute minimum surface magnetic field (quiet Sun)
double SurfaceMagneticFieldModule::computeB_s_min() {
    return variables["B_s_min"];
}
// Compute maximum surface magnetic field (sunspot max)
double SurfaceMagneticFieldModule::computeB_s_max() {
    return variables["B_s_max"];
}
// Compute scaled B_j (T)
double SurfaceMagneticFieldModule::computeB_j(double t, double B_s) {
    variables["t"] = t;
    double base_b = variables["B_ref"] + 0.4 * std::sin(variables["omega_s"] * t);  // Hypothetical cycle
    return base_b * (B_s / variables["B_ref"]);
}
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    } else {
        variables[name] += delta;
    }
}
// Subtract delta
void SurfaceMagneticFieldModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}
// Compute B_s min (T)
double SurfaceMagneticFieldModule::computeB_s_min() {
    return variables["B_s_min"];
}
// Compute B_s max (T)
double SurfaceMagneticFieldModule::computeB_s_max() {
    return variables["B_s_max"];
}
// Compute B_j scaled by B_s / B_ref
double SurfaceMagneticFieldModule::computeB_j(double t, double B_s) {
    variables["t"] = t;
    double base_b = variables["B_ref"] + 0.4 * std::sin(variables["omega_s"] * t);  // Hypothetical cycle
    return base_b * (B_s / variables["B_ref"]);
}
// U_g3 example with B_j
double SurfaceMagneticFieldModule::computeU_g3_example(double t, double B_s) {
    double k_3 = variables["k_3"];
    double b_j = computeB_j(t, B_s);
    double cos_term = std::cos(variables["omega_s"] * t * variables["pi"]);
    double p_core = variables["P_core"];
    double e_react = variables["E_react"];
    return k_3 * b_j * cos_term * p_core * e_react;
}
// Equation text
std::string SurfaceMagneticFieldModule::getEquationText() {
    return "B_j ? (10^3 + 0.4 sin(?_s t)) * (B_s / 0.4) T (hypothetical scaling);\n"
           "U_g3 = k_3 * ? B_j * cos(?_s t ?) * P_core * E_react\n"
           "Where B_s = [1e-4, 0.4] T (Sun surface; quiet to sunspot).\n";
}
// Update variable
void SurfaceMagneticFieldModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}
// Add delta
void SurfaceMagneticFieldModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << delta << std::endl;
        variables[name] = delta;
    }
}
// Constructor: Initialize all variables
SurfaceMagneticFieldModule::SurfaceMagneticFieldModule() {
    variables["B_s_min"] = 1e-4;                    // T
    variables["B_s_max"] = 0.4;                     // T
    variables["B_ref"] = 0.4;                       // T
    variables["omega_s"] = 2.7e-6;                  // rad/s (approx solar rotation rate)
    variables["k_3"] = 1.0e3;                       // hypothetical constant
    variables["P_core"] = 3.8e26;                   // W (solar core power)
    variables["E_react"] = 1.6e-13;                 // J (hypothetical reaction energy)
    variables["pi"] = 3.141592653589793;
    variables["t"] = 0.0;                           // s
}
// Compute B_j scaled by B_s / B_ref
double SurfaceMagneticFieldModule::computeB_j(double t, double B_s) {
    variables["t"] = t;
    double base_b = variables["B_ref"] + 0.4 * std::sin(variables["omega_s"] * t);  // Hypothetical cycle
    return base_b * (B_s / variables["B_ref"]);
}
// U_g3 example with B_j
double SurfaceMagneticFieldModule::computeU_g3_example(double t, double B_s) {
    double k_3 = variables["k_3"];
    double b_j = computeB_j(t, B_s);
    double cos_term = std::cos(variables["omega_s"] * t * variables["pi"]);
    double p_core = variables["P_core"];
    double e_react = variables["E_react"];
    return k_3 * b_j * cos_term * p_core * e_react;
}
// Equation text
std::string SurfaceMagneticFieldModule::getEquationText() {
    return "B_j ? (0.4 + 0.4 sin(?_s t)) * (B_s / 0.4) T (hypothetical scaling);\n"
           "U_g3 = k_3 * ? B_j * cos(?_s t ?) * P_core * E_react\n"
           "Where B_s = [1e-4, 0.4] T (Sun surface; quiet to sunspot).\n";
}
// Compute B_j scaled by B_s / B_ref
double SurfaceMagneticFieldModule::computeB_j(double t, double B_s) {
    variables["t"] = t;
    double base_b = variables["B_ref"] + 0.4 * std::sin(variables["omega_s"] * t);  // Hypothetical cycle
    return base_b * (B_s / variables["B_ref"]);
}
}
    } else {
        variables[name] += delta;
    }
}
}
// Subtract delta
void SurfaceMagneticFieldModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}
// Constructor: Initialize all variables
SurfaceMagneticFieldModule::SurfaceMagneticFieldModule() {
    variables["B_s_min"] = 1e-4;                    // T (quiet)
    variables["B_s_max"] = 0.4;                     // T (sunspot)
    variables["B_ref"] = 0.4;                       // T (reference max)
    variables["k_3"] = 1.8;                         // Coupling
    variables["omega_s"] = 2.5e-6;                  // rad/s
    variables["P_core"] = 1.0;                      // Unitless
    variables["E_react"] = 1e46;                    // J
    variables["pi"] = 3.141592653589793;
    variables["t"] = 0.0;                           // s
}
// Compute B_j scaled by B_s / B_ref
double SurfaceMagneticFieldModule::computeB_j(double t, double B_s) {
    variables["t"] = t;
    double base_b = variables["B_ref"] + 0.4 * std::sin(variables["omega_s"] * t);  // Hypothetical cycle
    return base_b * (B_s / variables["B_ref"]);
}
// U_g3 example with B_j
double SurfaceMagneticFieldModule::computeU_g3_example(double t, double B_s) {
    double k_3 = variables["k_3"];
    double b_j = computeB_j(t, B_s);
    double cos_term = std::cos(variables["omega_s"] * t * variables["pi"]);
    double p_core = variables["P_core"];
    double e_react = variables["E_react"];
    return k_3 * b_j * cos_term * p_core * e_react;
}
// Equation text
std::string SurfaceMagneticFieldModule::getEquationText() {
    return "B_j ? (0.4 + 0.4 sin(?_s t)) * (B_s / 0.4) T (hypothetical scaling);\n"
           "U_g3 = k_3 * ? B_j * cos(?_s t ?) * P_core * E_react\n"
           "Where B_s = [1e-4, 0.4] T (Sun surface; quiet to sunspot).\n";
}
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}
// Subtract delta
void SurfaceMagneticFieldModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}
}
// Constructor: Initialize all variables
SurfaceMagneticFieldModule::SurfaceMagneticFieldModule() {
    variables["B_s_min"] = 1e-4;                    // T (quiet)
    variables["B_s_max"] = 0.4;                     // T (sunspot)
    variables["B_ref"] = 0.4;                       // T (reference max)
    variables["k_3"] = 1.8;                         // Coupling
    variables["omega_s"] = 2.5e-6;                  // rad/s
    variables["P_core"] = 1.0;                      // Unitless
    variables["E_react"] = 1e46;                    // J
    variables["pi"] = 3.141592653589793;
    variables["t"] = 0.0;                           // s
}
}// Compute B_j scaled by B_s / B_ref
double SurfaceMagneticFieldModule::computeB_j(double t, double B_s) {
    variables["t"] = t;
    double base_b = variables["B_ref"] + 0.4 * std::sin(variables["omega_s"] * t);  // Hypothetical cycle
    return base_b * (B_s / variables["B_ref"]);
}
// U_g3 example with B_j
double SurfaceMagneticFieldModule::computeU_g3_example(double t, double B_s) {
    double k_3 = variables["k_3"];
    double b_j = computeB_j(t, B_s);
    double cos_term = std::cos(variables["omega_s"] * t * variables["pi"]);
    double p_core = variables["P_core"];
    double e_react = variables["E_react"];
    return k_3 * b_j * cos_term * p_core * e_react;
}
// Equation text
std::string SurfaceMagneticFieldModule::getEquationText() {
    return "B_j ? (0.4 + 0.4 sin(?_s t)) * (B_s / 0.4) T (hypothetical scaling);\n"
           "U_g3 = k_3 * ? B_j * cos(?_s t ?) * P_core * E_react\n"
           "Where B_s = [1e-4, 0.4] T (Sun surface; quiet to sunspot).\n";
}// Update variable
void SurfaceMagneticFieldModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}
// Add delta
void SurfaceMagneticFieldModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}
// Subtract delta
void SurfaceMagneticFieldModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}
}
// Constructor: Initialize all variables
SurfaceMagneticFieldModule::SurfaceMagneticFieldModule() {
    variables["B_s_min"] = 1e-4;                    // T (quiet)
    variables["B_s_max"] = 0.4;                     // T (sunspot)
    variables["B_ref"] = 0.4;                       // T (reference max)
    variables["k_3"] = 1.8;                         // Coupling
    variables["omega_s"] = 2.5e-6;                  // rad/s
    variables["P_core"] = 1.0;                      // Unitless
    variables["E_react"] = 1e46;                    // J
    variables["pi"] = 3.141592653589793;
    variables["t"] = 0.0;                           // s
}
// Compute B_j scaled by B_s / B_ref
double SurfaceMagneticFieldModule::computeB_j(double t, double B_s) {
    variables["t"] = t;
    double base_b = variables["B_ref"] + 0.4 * std::sin(variables["omega_s"] * t);  // Hypothetical cycle
    return base_b * (B_s / variables["B_ref"]);
}
// U_g3 example with B_j
double SurfaceMagneticFieldModule::computeU_g3_example(double t, double B_s) {
    double k_3 = variables["k_3"];
    double b_j = computeB_j(t, B_s);
    double cos_term = std::cos(variables["omega_s"] * t * variables["pi"]);
    double p_core = variables["P_core"];
    double e_react = variables["E_react"];
    return k_3 * b_j * cos_term * p_core * e_react;
}
// Equation text
std::string SurfaceMagneticFieldModule::getEquationText() {
    return "B_j ? (0.4 + 0.4 sin(?_s t)) * (B_s / 0.4) T (hypothetical scaling);\n"
           "U_g3 = k_3 * ? B_j * cos(?_s t ?) * P_core * E_react\n"
           "Where B_s = [1e-4, 0.4] T (Sun surface; quiet to sunspot).\n";
}// Update variable
void SurfaceMagneticFieldModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}
// Add delta
void SurfaceMagneticFieldModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void SurfaceMagneticFieldModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}
}

// Constructor: Initialize all variables
SurfaceMagneticFieldModule::SurfaceMagneticFieldModule() {
    variables["B_s_min"] = 1e-4;                    // T (quiet)
    variables["B_s_max"] = 0.4;                     // T (sunspot)
    variables["B_ref"] = 0.4;                       // T (reference max)
    variables["k_3"] = 1.8;                         // Coupling
    variables["omega_s"] = 2.5e-6;                  // rad/s
    variables["P_core"] = 1.0;                      // Unitless
    variables["E_react"] = 1e46;                    // J
    variables["pi"] = 3.141592653589793;
    variables["t"] = 0.0;                           // s
}
// Compute B_j scaled by B_s / B_ref
double SurfaceMagneticFieldModule::computeB_j(double t, double B_s) {
    variables["t"] = t;
    double base_b = variables["B_ref"] + 0.4 * std::sin(variables["omega_s"] * t);  // Hypothetical cycle
    return base_b * (B_s / variables["B_ref"]);
}
// U_g3 example with B_j
double SurfaceMagneticFieldModule::computeU_g3_example(double t, double B_s) {
    double k_3 = variables["k_3"];
    double b_j = computeB_j(t, B_s);
    double cos_term = std::cos(variables["omega_s"] * t * variables["pi"]);
    double p_core = variables["P_core"];
    double e_react = variables["E_react"];
    return k_3 * b_j * cos_term * p_core * e_react;
}
// Equation text
std::string SurfaceMagneticFieldModule::getEquationText() {
    return "B_j ? (0.4 + 0.4 sin(?_s t)) * (B_s / 0.4) T (hypothetical scaling);\n"
           "U_g3 = k_3 * ? B_j * cos(?_s t ?) * P_core * E_react\n"
           "Where B_s = [1e-4, 0.4] T (Sun surface; quiet to sunspot).\n";
}// Compute B_j scaled by B_s / B_ref
double SurfaceMagneticFieldModule::computeB_j(double t, double B_s) {
    variables["t"] = t;
    double base_b = variables["B_ref"] + 0.4 * std::sin(variables["omega_s"] * t);  // Hypothetical cycle
    return base_b * (B_s / variables["B_ref"]);
}
// U_g3 example with B_j
double SurfaceMagneticFieldModule::computeU_g3_example(double t, double B_s) {
    double k_3 = variables["k_3"];
    double b_j = computeB_j(t, B_s);
    double cos_term = std::cos(variables["omega_s"] * t * variables["pi"]);
    double p_core = variables["P_core"];
    double e_react = variables["E_react"];
    return k_3 * b_j * cos_term * p_core * e_react;
}
// Equation text
std::string SurfaceMagneticFieldModule::getEquationText() {
    return "B_j ? (0.4 + 0.4 sin(?_s t)) * (B_s / 0.4) T (hypothetical scaling);\n"
           "U_g3 = k_3 * ? B_j * cos(?_s t ?) * P_core * E_react\n"
           "Where B_s = [1e-4, 0.4] T (Sun surface; quiet to sunspot).\n";
}// Update variable
void SurfaceMagneticFieldModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}
}
// Add delta
void SurfaceMagneticFieldModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}
}
// Subtract delta
void SurfaceMagneticFieldModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}
}
}
// Constructor: Initialize all variables
SurfaceMagneticFieldModule::SurfaceMagneticFieldModule() {
    variables["B_s_min"] = 1e-4;                    // T (quiet)
    variables["B_s_max"] = 0.4;                     // T (sunspot)
    variables["B_ref"] = 0.4;                       // T (reference max)
    variables["k_3"] = 1.8;                         // Coupling
    variables["omega_s"] = 2.5e-6;                  // rad/s
    variables["P_core"] = 1.0;                      // Unitless
    variables["E_react"] = 1e46;                    // J
    variables["pi"] = 3.141592653589793;
    variables["t"] = 0.0;                           // s
}
}// Compute B_j scaled by B_s / B_ref
double SurfaceMagneticFieldModule::computeB_j(double t, double B_s) {
    variables["t"] = t;
    double base_b = variables["B_ref"] + 0.4 * std::sin(variables["omega_s"] * t);  // Hypothetical cycle
    return base_b * (B_s / variables["B_ref"]);
}
}
// U_g3 example with B_j
double SurfaceMagneticFieldModule::computeU_g3_example(double t, double B_s) {
    double k_3 = variables["k_3"];
    double b_j = computeB_j(t, B_s);
    double cos_term = std::cos(variables["omega_s"] * t * variables["pi"]);
    double p_core = variables["P_core"];
    double e_react = variables["E_react"];
    return k_3 * b_j * cos_term * p_core * e_react;
}
}
}// Equation text
std::string SurfaceMagneticFieldModule::getEquationText() {
    return "B_j ? (0.4 + 0.4 sin(?_s t)) * (B_s / 0.4) T (hypothetical scaling);\n"
           "U_g3 = k_3 * ? B_j * cos(?_s t ?) * P_core * E_react\n"
           "Where B_s = [1e-4, 0.4] T (Sun surface; quiet to sunspot).\n";
}// Update variable
void SurfaceMagneticFieldModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}
}
// UQFFBuoyancyModule.cpp
#include "UQFFBuoyancyModule.h"
#include <complex>

// Constructor: Set all variables with multi-system defaults
UQFFBuoyancyModule::UQFFBuoyancyModule() {
    double pi_val = 3.141592653589793;
    cdouble zero = {0.0, 0.0};
    cdouble i_small = {0.0, 1e-37};

    // Base constants (universal)
    variables["G"] = {6.6743e-11, 0.0};
    variables["c"] = {3e8, 0.0};
    variables["hbar"] = {1.0546e-34, 0.0};
    variables["q"] = {1.6e-19, 0.0};
    variables["pi"] = {pi_val, 0.0};
    variables["m_e"] = {9.11e-31, 0.0};
    variables["mu_B"] = {9.274e-24, 0.0};
    variables["g_Lande"] = {2.0, 0.0};
    variables["k_B"] = {1.38e-23, 0.0};
    variables["mu0"] = {4 * pi_val * 1e-7, 0.0};

    // Shared params
    variables["F_rel"] = {4.30e33, 0.0};  // Relativistic coherence from LEP 1998
    variables["F0"] = {1.83e71, 0.0};
    variables["V"] = {1e-3, 0.0};  // Default particle velocity
    variables["theta"] = {pi_val / 4, 0.0};  // 45 deg
    variables["phi"] = {pi_val / 4, 0.0};
    variables["omega_act"] = {2 * pi_val * 300, 0.0};
    variables["k_act"] = {1e-6, 0.0};
    variables["k_DE"] = {1e-30, 0.0};
    variables["k_neutron"] = {1e10, 0.0};
    variables["sigma_n"] = {1e-4, 0.0};
    variables["k_rel"] = {1e-10, 0.0};
    variables["E_cm_astro"] = {1.24e24, 0.0};
    variables["E_cm"] = {3.0264e-8, 0.0};  // 189 GeV in J
    variables["F_neutrino"] = {9.07e-42, 1e-43};
    variables["k_LENR"] = {1e-10, 0.0};
    variables["omega_LENR"] = {2 * pi_val * 1.25e12, 0.0};
    variables["rho_vac_UA"] = {7.09e-36, 1e-37};
    variables["DPM_momentum"] = {0.93, 0.05};
    variables["DPM_gravity"] = {1.0, 0.1};
    variables["DPM_stability"] = {0.01, 0.001};
    variables["beta_i"] = {0.6, 0.0};
    variables["V_infl_UA"] = {1e-6, 1e-7};
    variables["rho_vac_A"] = {1e-30, 1e-31};
    variables["a_universal"] = {1e12, 1e11};
    variables["lambda_i"] = {1.0, 0.0};
    variables["rho_vac_SCm"] = {7.09e-37, 1e-38};
    variables["omega_s"] = {2.5e-6, 0.0};
    variables["k_s"] = {1.8, 0.0};
    variables["P_core"] = {1.0, 0.0};
    variables["E_react"] = {1e46, 0.0};
    variables["t"] = {0.0, 0.0};
    // System-specific params (M104 as default)
    setSystemParams("M104");
}
// Set system-specific parameters
void UQFFBuoyancyModule::setSystemParams(const std::string& system) {
    if (system == "M104") {
        // M104 specific parameters
        variables["k_s"] = {1.8, 0.0};
        variables["P_core"] = {1.0, 0.0};
        variables["E_react"] = {1e46, 0.0};
        variables["t"] = {0.0, 0.0};
    }
}
// Update variable
void UQFFBuoyancyModule::updateVariable(const std::string& name, cdouble value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}
}
// Add delta
void UQFFBuoyancyModule::addToVariable(const std::string& name, cdouble delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}
}
// Add delta
void UQFFBuoyancyModule::addToVariable(const std::string& name, cdouble delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}
}
// Add delta
void UQFFBuoyancyModule::addToVariable(const std::string& name, cdouble delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}
}
}
// Subtract delta
void UQFFBuoyancyModule::subtractFromVariable(const std::string& name, cdouble delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] -= delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with negative delta " << delta << std::endl;
        variables[name] = -delta;
    }
}
}
}
// Subtract delta
void UQFFBuoyancyModule::subtractFromVariable(const std::string& name, cdouble delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] -= delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with negative delta " << delta << std::endl;
        variables[name] = -delta;
    }
}
}
}
// Subtract delta
void UQFFBuoyancyModule::subtractFromVariable(const std::string& name, cdouble delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] -= delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with negative delta " << delta << std::endl;
        variables[name] = -delta;
    }
}
}
}
// Compute integrand for F_U_Bi_i
cdouble UQFFBuoyancyModule::computeIntegrand(double t, const std::
string& system) {
    cdouble F0 = variables["F0"];
    cdouble V = variables["V"];
    cdouble theta = variables["theta"];
    cdouble phi = variables["phi"];
    cdouble omega_act = variables["omega_act"];
    cdouble k_act = variables["k_act"];
    cdouble k_DE = variables["k_DE"];
    cdouble k_neutron = variables["k_neutron"];
    cdouble sigma_n = variables["sigma_n"];
    cdouble k_rel = variables["k_rel"];
    cdouble E_cm_astro = variables["E_cm_astro"];
    cdouble E_cm = variables["E_cm"];
    cdouble F_neutrino = variables["F_neutrino"];
    cdouble k_LENR = variables["k_LENR"];
    cdouble omega_LENR = variables["omega_LENR"];
    cdouble rho_vac_UA = variables["rho_vac_UA"];
    cdouble DPM_momentum = variables["DPM_momentum"];
    cdouble DPM_gravity = variables["DPM_gravity"];
    cdouble DPM_stability = variables["DPM_stability"];
    cdouble beta_i = variables["beta_i"];
    cdouble V_infl_UA = variables["V_infl_UA"];
    cdouble rho_vac_A = variables["rho_vac_A"];
    cdouble a_universal = variables["a_universal"];
    cdouble lambda_i = variables["lambda_i"];
    cdouble rho_vac_SCm = variables["rho_vac_SCm"];

    // Compute terms
    cdouble term_activation = k_act * std::sin(omega_act * t) * std::exp(-k_act * t);
    cdouble term_directed_energy = k_DE * V * std::sin(theta) * std::cos(phi);
    cdouble term_neutron = k_neutron * sigma_n * E_cm / E_cm_astro;
    cdouble term_relativistic = k_rel * (E_cm / (variables["m_e"] * std::pow(variables["c"], 2.0)));
    cdouble term_neutrino = F_neutrino;
    cdouble term_LENR = computeLENRTerm(system);
    cdouble term_vacuum_stability = rho_vac_UA * DPM_stability;
    
    return F0 * (term_activation + term_directed_energy + term_neutron + term_relativistic +
                 term_neutrino + term_LENR + term_vacuum_stability);
}
}
// Subtract delta
void UQFFBuoyancyModule::subtractFromVariable(const std::string& name, cdouble delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] -= delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with negative delta " << delta << std::endl;
        variables[name] = -delta;
    }
}
}
// Set system-specific parameters
void UQFFBuoyancyModule::setSystemParams(const std::string& system) {
    if (system == "M104") {
        // M104 specific parameters
        variables["k_s"] = {1.8, 0.0};
        variables["P_core"] = {1.0, 0.0};
        variables["E_react"] = {1e46, 0.0};
        variables["t"] = {0.0, 0.0};
    }
}
// Compute LENR term
cdouble UQFFBuoyancyModule::computeLENRTerm(const std::string& system) {
    cdouble k_LENR = variables["k_LENR"];
    cdouble omega_LENR = variables["omega_LENR"];
    cdouble t = variables["t"];
    return k_LENR * std::sin(omega_LENR * t) * std::exp(-k_LENR * t);
}
// Compute integrand for F_U_Bi_i
cdouble UQFFBuoyancyModule::computeIntegrand(double t, const std::
string& system) {
    cdouble F0 = variables["F0"];
    cdouble V = variables["V"];
    cdouble theta = variables["theta"];
    cdouble phi = variables["phi"];
    cdouble omega_act = variables["omega_act"];
    cdouble k_act = variables["k_act"];
    cdouble k_DE = variables["k_DE"];
    cdouble k_neutron = variables["k_neutron"];
    cdouble sigma_n = variables["sigma_n"];
    cdouble k_rel = variables["k_rel"];
    cdouble E_cm_astro = variables["E_cm_astro"];
    cdouble E_cm = variables["E_cm"];
    cdouble F_neutrino = variables["F_neutrino"];
    cdouble k_LENR = variables["k_LENR"];
    cdouble omega_LENR = variables["omega_LENR"];
    cdouble rho_vac_UA = variables["rho_vac_UA"];
    cdouble DPM_momentum = variables["DPM_momentum"];
    cdouble DPM_gravity = variables["DPM_gravity"];
    cdouble DPM_stability = variables["DPM_stability"];
    cdouble beta_i = variables["beta_i"];
    cdouble V_infl_UA = variables["V_infl_UA"];
    cdouble rho_vac_A = variables["rho_vac_A"];
    cdouble a_universal = variables["a_universal"];
    cdouble lambda_i = variables["lambda_i"];
    cdouble rho_vac_SCm = variables["rho_vac_SCm"];

    // Compute terms
    cdouble term_activation = k_act * std::sin(omega_act * t) * std::exp(-k_act * t);
    cdouble term_directed_energy = k_DE * V * std::sin(theta) * std::cos(phi);
    cdouble term_neutron = k_neutron * sigma_n * E_cm / E_cm_astro;
    cdouble term_relativistic = k_rel * (E_cm / (variables["m_e"] * std::pow(variables["c"], 2.0)));
    cdouble term_neutrino = F_neutrino;
    cdouble term_LENR = computeLENRTerm(system);
    cdouble term_vacuum_stability = rho_vac_UA * DPM_stability;
    
    return F0 * (term_activation + term_directed_energy + term_neutron + term_relativistic +
                 term_neutrino + term_LENR + term_vacuum_stability);
}
}
// U_g3 example with B_j
double SurfaceMagneticFieldModule::computeU_g3_example(double t, double B_s) {
    double k_3 = variables["k_3"];
    double b_j = computeB_j(t, B_s);
    double cos_term = std::cos(variables["omega_s"] * t * variables["pi"]);
    double p_core = variables["P_core"];
    double e_react = variables["E_react"];
    return k_3 * b_j * cos_term * p_core * e_react;
}
}