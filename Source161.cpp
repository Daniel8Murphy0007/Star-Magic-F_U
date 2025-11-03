
// UQFFBuoyancyAstroModule.h
// Modular C++ implementation of the full Master Unified Field Equation (F_U_Bi_i & UQFF Integration) for Buoyancy Equations across J1610+1811, PLCK G287.0+32.9, PSZ2 G181.06+48.47, ASKAP J1832-0911, Sonification Collection.
// This module can be plugged into a base program (e.g., 'uqff_buoyancy_sim.cpp') by including this header and linking the .cpp.
// Usage in base: #include "UQFFBuoyancyAstroModule.h"
// UQFFBuoyancyAstroModule mod; mod.computeFBi(system); mod.updateVariable("F_rel", {new_real, new_imag});
// All variables are stored in a std::map for dynamic addition/subtraction/update, using complex<double> for real/imaginary components.
// Nothing is negligible: Includes all terms - base force, momentum, gravity, vacuum stability, LENR resonance, activation, directed energy, magnetic resonance, neutron, relativistic, neutrino, Sweet vac, Kozima drop.
// Associated text: Outputs descriptive equation string via getEquationText().
// Approximations: Integral approximated as integrand * x2 (quadratic root); imag parts small and not fully scaled; LENR dominant due to low ω_0; x2 from quadratic solver approx; F_rel from 1998 LEP.
// Multi-system params: J1610+1811 M=2.785e30 kg r=3.09e15 m; PLCK G287.0+32.9 M=1.989e44 kg r=3.09e22 m; PSZ2 G181.06+48.47 M=1.989e44 kg r=3.09e22 m; ASKAP J1832-0911 M=2.785e30 kg r=4.63e16 m; Sonification Collection M=1.989e31 kg r=6.17e16 m.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 22, 2025.

#ifndef UQFF_BUOYANCY_ASTRO_MODULE_H
#define UQFF_BUOYANCY_ASTRO_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>
#include <sstream>

using cdouble = std::complex<double>;

class UQFFBuoyancyAstroModule {
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
    UQFFBuoyancyAstroModule();

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

#endif // UQFF_BUOYANCY_ASTRO_MODULE_H

// UQFFBuoyancyAstroModule.cpp
#include "UQFFBuoyancyAstroModule.h"

// Constructor: Set all variables with multi-system defaults
UQFFBuoyancyAstroModule::UQFFBuoyancyAstroModule() {
    double pi_val = 3.141592653589793;

    // Base constants (universal)
    variables["G"] = {6.6743e-11, 0.0}; // Gravitational constant
    variables["c"] = {3e8, 0.0};
    variables["hbar"] = {1.0546e-34, 0.0};
    variables["q"] = {1.6e-19, 0.0};
    variables["pi"] = {pi_val, 0.0};
    variables["m_e"] = {9.11e-31, 0.0};
    variables["mu_B"] = {9.274e-24, 0.0};
    variables["g_Lande"] = {2.0, 0.0};
    variables["k_B"] = {1.38e-23, 0.0};
    variables["mu0"] = {4 * pi_val * 1e-7, 0.0};
    variables["epsilon0"] = {8.85e-12, 0.0}; // For quadratic terms

    // Shared params from document
    variables["F0"] = {1.83e71, 0.0};
    variables["V"] = {1e-3, 0.0};  // Default particle velocity
    variables["theta"] = {pi_val / 4, 0.0};  // 45 deg
    variables["phi"] = {pi_val / 4, 0.0};
    variables["omega_act"] = {2 * pi_val * 300, 0.0};
    variables["k_act"] = {1e-6, 0.0};
    variables["k_DE"] = {1e-30, 0.0};
    variables["k_neutron"] = {1e10, 0.0};
    variables["k_rel"] = {1e-10, 0.0};
    variables["k_relativistic"] = {1e-20, 0.0};
    variables["k_neutrino"] = {1e-15, 0.0};
    variables["k_Sweet"] = {1e-25, 0.0};
    variables["k_Kozima"] = {1e-18, 0.0};
    variables["omega_0_LENR"] = {2 * pi_val * 1.25e12, 0.0};  // LENR resonance freq 1.2-1.3 THz
    variables["k_LENR"] = {1e-10, 0.0};
    variables["rho_vac_UA"] = {7.09e-36, 0.0};  // Vacuum energy density
    variables["sigma_n"] = {1e-4, 0.0};  // Scaled for astrophysical densities
    variables["E_cm"] = {189.0, 0.0};  // GeV
    variables["E_cm_astro_local_adj_eff_enhanced"] = {1.24e24, 0.0};  // events/m3
    variables["DPM_stability"] = {0.01, 0.0};
    variables["DPM_momentum"] = {0.93, 0.0};
    variables["DPM_gravity"] = {1.0, 0.0};
    variables["DPM_light"] = {0.01, 0.0};
    variables["DPM_resonance"] = {1.0, 0.0};  // Default
    variables["phase"] = {2.36e-3, 0.0};  // s^-1
    variables["curvature"] = {1e-22, 0.0};
    variables["k_10_13"] = {1e-13, 0.0};  // For light term in quadratic
    variables["k_b_term"] = {2.51e-5, 0.0};  // Constant in b
    variables["c_constant"] = {-3.06e175, 0.0};  // Constant in c
    variables["c_inv_r2_coeff"] = {1e-29, 0.0};  // 10^{-29}/r^2 in c
    variables["a_eps_coeff"] = {1.38e-41, 0.0};  // Coefficient for first term in a (as per document)

    // System-specific params will be set in setSystemParams()
}

// Set system-specific parameters
void UQFFBuoyancyAstroModule::setSystemParams(const std::string& system)
{
    double pi_val = variables["pi"].real();
    if (system == "J1610+1811") {
        this->variables["M"] = {2.785e30, 0.0};
        this->variables["r"] = {3.09e15, 0.0};
        this->variables["T"] = {1e4, 0.0};
        this->variables["L_X"] = {1e31, 0.0};
        this->variables["B0"] = {1e-4, 0.0};
        this->variables["omega0"] = {1e-12, 0.0};
        this->variables["Mach"] = {1.0, 0.0};  // ℳ
        this->variables["C"] = {1.0, 0.0};
        this->variables["theta"] = {pi_val / 4, 0.0};
        this->variables["t"] = {3.156e10, 0.0};
    } else if (system == "PLCK_G287.0+32.9") {
        this->variables["M"] = {1.989e44, 0.0};
        this->variables["r"] = {3.09e22, 0.0};
        this->variables["T"] = {1e7, 0.0};
        this->variables["L_X"] = {1e38, 0.0};
        this->variables["B0"] = {1e-4, 0.0};
        this->variables["omega0"] = {1e-15, 0.0};
        this->variables["Mach"] = {1.5, 0.0};
        this->variables["C"] = {1.2, 0.0};
        this->variables["theta"] = {pi_val / 4, 0.0};
        this->variables["t"] = {1.42e17, 0.0};
    } else if (system == "PSZ2_G181.06+48.47") {
        this->variables["M"] = {1.989e44, 0.0};
        this->variables["r"] = {3.09e22, 0.0};
        this->variables["T"] = {1e7, 0.0};
        this->variables["L_X"] = {1e39, 0.0};
        this->variables["B0"] = {1e-4, 0.0};
        this->variables["omega0"] = {1e-15, 0.0};
        this->variables["Mach"] = {1.5, 0.0};
        this->variables["C"] = {1.2, 0.0};
        this->variables["theta"] = {pi_val / 4, 0.0};
        this->variables["t"] = {2.36e17, 0.0};
    } else if (system == "ASKAP_J1832-0911") {
        this->variables["M"] = {2.785e30, 0.0};
        this->variables["r"] = {4.63e16, 0.0};
        this->variables["T"] = {1e4, 0.0};
        this->variables["L_X"] = {1e31, 0.0};
        this->variables["B0"] = {1e-4, 0.0};
        this->variables["omega0"] = {1e-12, 0.0};  // From 44-min period approx
        this->variables["Mach"] = {1.0, 0.0};
        this->variables["C"] = {1.0, 0.0};
        this->variables["theta"] = {pi_val / 4, 0.0};
        this->variables["t"] = {3.156e10, 0.0};
    } else if (system == "SonificationCollection") {
        this->variables["M"] = {1.989e31, 0.0};
        this->variables["r"] = {6.17e16, 0.0};
        this->variables["T"] = {1e5, 0.0};
        this->variables["L_X"] = {1e33, 0.0};
        this->variables["B0"] = {1e-5, 0.0};
        this->variables["omega0"] = {1e-12, 0.0};
        this->variables["Mach"] = {1.0, 0.0};
        this->variables["C"] = {1.0, 0.0};
        this->variables["theta"] = {pi_val / 4, 0.0};
        this->variables["t"] = {3.156e14, 0.0};
    }
}

// Dynamic variable operations (complex)
void UQFFBuoyancyAstroModule::updateVariable(const std::string& name, cdouble value) {
    this->variables[name] = value;
}
void UQFFBuoyancyAstroModule::addToVariable(const std::string& name, cdouble delta) {
    this->variables[name] += delta;
}
void UQFFBuoyancyAstroModule::subtractFromVariable(const std::string& name, cdouble delta) {
    this->variables[name] -= delta;
}

// Core computation: Full F_U_Bi_i(r, t) for system (approx integral)
cdouble UQFFBuoyancyAstroModule::computeFBi(const std::string& system, double t) {
    setSystemParams(system);
    cdouble integrand = computeIntegrand(t, system);
    cdouble x2 = computeX2(system);
    cdouble f_bi_i = integrand * x2;
    double cos_theta = cos(variables["theta"].real());
    cdouble momentum_term = (variables["m_e"] * variables["c"] * variables["c"] / (variables["r"] * variables["r"])) * variables["DPM_momentum"] * cos_theta;
    cdouble gravity_term = (variables["G"] * variables["M"] / (variables["r"] * variables["r"])) * variables["DPM_gravity"];
    cdouble f_bi = -variables["F0"] + momentum_term + gravity_term + f_bi_i;
    return f_bi;
}

// Sub-equations
cdouble UQFFBuoyancyAstroModule::computeCompressed(const std::string& system, double t) {
    setSystemParams(system);
    return computeIntegrand(t, system);
}
cdouble UQFFBuoyancyAstroModule::computeResonant(const std::string& system) {
    setSystemParams(system);
    return computeDPM_resonance(system);
}
cdouble UQFFBuoyancyAstroModule::computeBuoyancy(const std::string& system) {
    setSystemParams(system);
    return computeUb1(system);
}
cdouble UQFFBuoyancyAstroModule::computeSuperconductive(const std::string& system, double t) {
    setSystemParams(system);
    return computeUi(t, system);
}
double UQFFBuoyancyAstroModule::computeCompressedG(const std::string& system, double t) {
    setSystemParams(system);
    return computeG(t, system);
}

// Output descriptive text of the equation
std::string UQFFBuoyancyAstroModule::getEquationText(const std::string& system) {
    setSystemParams(system);
    std::ostringstream oss;
    oss << "F_U_Bi_i(r, t) = Integral[Integrand(r, t) dt] approximated as Integrand * x2\n";
    oss << "Where Integrand includes terms for base force, momentum, gravity, vacuum stability, LENR resonance, activation, directed energy, magnetic resonance, neutron, relativistic, neutrino, Sweet vac, Kozima drop.\n";
    oss << "LENR Resonance: F_LENR = k_LENR * (ω_LENR / ω_0)^2\n";
    oss << "Activation: F_act = k_act * cos(ω_act t)\n";
    oss << "Directed Energy: F_DE = k_DE * L_X\n";
    oss << "Magnetic Resonance: F_res = 2 q B_0 V sinθ * DPM_resonance\n";
    oss << "Neutron Drop: F_neutron = k_neutron * σ_n\n";
    oss << "Relativistic: F_rel = k_rel * (E_cm_astro_local_adj_eff_enhanced / E_cm)^2\n";
    oss << "Neutrino: F_neutrino = k_neutrino * L_X\n";
    oss << "Sweet Vac: F_sweet = k_Sweet * ρ_vac_UA\n";
    oss << "Kozima Drop: F_kozima = k_Kozima * σ_n\n";
    oss << "Relativistic Correction: F_relativ = k_relativistic * (V / c)^2 * F0\n";
    oss << "System: " << system << "\n";
    return oss.str();
}

// Print all current variables (for debugging/updates)
void UQFFBuoyancyAstroModule::printVariables() {
    for (const auto& pair : variables) {
        std::cout << std::setw(15) << pair.first << " : " << pair.second << std::endl;
    }
}

// Compute integrand for F_U_Bi_i
cdouble UQFFBuoyancyAstroModule::computeIntegrand(double t, const std::string& system) {
    setSystemParams(system);
    double pi_val = variables["pi"].real();
    double sin_theta = sin(variables["theta"].real());
    double cos_theta = cos(variables["theta"].real());
    cdouble dpm_res = computeDPM_resonance(system);
    cdouble f_lenr = computeLENRTerm(system);
    cdouble f_act = variables["k_act"] * cos(variables["omega_act"].real() * t);
    cdouble f_de = variables["k_DE"] * variables["L_X"];
    cdouble f_neutron = variables["k_neutron"] * variables["sigma_n"];
    cdouble f_rel = variables["k_rel"] * pow(variables["E_cm_astro_local_adj_eff_enhanced"].real() / variables["E_cm"].real(), 2.0);
    cdouble f_res = 2.0 * variables["q"].real() * variables["B0"].real() * variables["V"].real() * sin_theta * dpm_res;
    cdouble f_relativ = variables["k_relativistic"] * pow(variables["V"].real() / variables["c"].real(), 2.0) * variables["F0"];
    cdouble f_neutrino = variables["k_neutrino"] * variables["L_X"];
    cdouble f_sweet = variables["k_Sweet"] * variables["rho_vac_UA"];
    cdouble f_kozima = variables["k_Kozima"] * variables["sigma_n"];
    cdouble momentum_term = (variables["m_e"] * variables["c"] * variables["c"] / (variables["r"] * variables["r"])) * variables["DPM_momentum"] * cos_theta;
    cdouble gravity_term = (variables["G"] * variables["M"] / (variables["r"] * variables["r"])) * variables["DPM_gravity"];
    cdouble vac_term = variables["rho_vac_UA"] * variables["DPM_stability"];
    cdouble integrand = -variables["F0"] + momentum_term + gravity_term + vac_term + f_lenr + f_act + f_de + f_res + f_neutron + f_rel + f_relativ + f_neutrino + f_sweet + f_kozima;
    return integrand;
}

// Compute DPM resonance term
cdouble UQFFBuoyancyAstroModule::computeDPM_resonance(const std::string& system) {
    setSystemParams(system);
    double g_lande = variables["g_Lande"].real();
    double mu_b = variables["mu_B"].real();
    double b0 = variables["B0"].real();
    double hbar_omega0 = variables["hbar"].real() * variables["omega0"].real();
    if (hbar_omega0 == 0.0) return {0.0, 0.0};
    return {g_lande * mu_b * b0 / hbar_omega0, 0.0};
}

// Compute x2 from quadratic root approximation (negative root as per doc)
cdouble UQFFBuoyancyAstroModule::computeX2(const std::string& system) {
    setSystemParams(system);
    double r_real = variables["r"].real();
    double r2 = r_real * r_real;
    double t_val = variables["T"].real();
    double m = variables["M"].real();
    double pi_val = variables["pi"].real();
    // a terms from parsed document formula
    double term1_num = variables["a_eps_coeff"].real() * variables["q"].real();
    double term1_denom = 4.0 * pi_val * variables["epsilon0"].real() * r2 * t_val;
    double term1 = term1_num / term1_denom;
    double term2 = variables["G"].real() * m / r2;
    double term3 = pow(variables["c"].real(), 4.0) * variables["k_10_13"].real() / r2 * variables["DPM_light"].real();
    cdouble a = {term1 + term2 + term3, 0.0};
    // b from parsed
    double term_b1 = variables["k_b_term"].real();
    double term_b2 = t_val / r2;
    double term_b3 = 2.0 * variables["phase"].real();
    cdouble b = {term_b1 + term_b2 + term_b3, 0.0};
    // c from parsed
    double term_c1 = variables["c_constant"].real();
    double term_c2 = variables["c_inv_r2_coeff"].real() / r2;
    double term_c3 = variables["curvature"].real();
    cdouble c = {term_c1 + term_c2 + term_c3, 0.0};
    return computeQuadraticRoot(a, b, c);
}

// Compute quadratic root (negative branch: [-b - sqrt(b^2 - 4ac)] / 2a)
cdouble UQFFBuoyancyAstroModule::computeQuadraticRoot(cdouble a, cdouble b, cdouble c) {
    cdouble disc = b * b - 4.0 * a * c;
    double disc_real = disc.real();
    if (disc_real < 0.0) disc_real = 0.0;  // Force real for approximation
    cdouble sqrt_disc = {sqrt(disc_real), 0.0};
    return (-b - sqrt_disc) / (2.0 * a);
}

// Compute LENR term
cdouble UQFFBuoyancyAstroModule::computeLENRTerm(const std::string& system) {
    setSystemParams(system);
    double omega0_real = variables["omega0"].real();
    if (omega0_real == 0.0) return {0.0, 0.0};
    return variables["k_LENR"] * pow(variables["omega_0_LENR"].real() / omega0_real, 2.0);
}

// Compute gravitational acceleration g(r,t) - as per document
double UQFFBuoyancyAstroModule::computeG(double t, const std::string& system) {
    return -1.07e16;
}

// Compute Q_wave term - as per document
cdouble UQFFBuoyancyAstroModule::computeQ_wave(double t, const std::string& system) {
    return {3.11e5, 0.0};
}

// Compute Ub1 buoyancy term
cdouble UQFFBuoyancyAstroModule::computeUb1(const std::string& system) {
    return computeIntegrand(0.0, system);
}

// Compute Ui superconductive term
cdouble UQFFBuoyancyAstroModule::computeUi(double t, const std::string& system) {
    return computeQ_wave(t, system);
}

// Dynamic method test
void UQFFBuoyancyAstroModule::enableSelfLearning(bool enable) {
    self_learning_enabled = enable;
}
