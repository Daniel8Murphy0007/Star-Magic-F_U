
// UQFFBuoyancyCNBModule.h
// Modular C++ implementation of the full Master Unified Field Equation (F_U_Bi_i & UQFF Integration) for Buoyancy Equations across J1610+1811, PLCK G287.0+32.9, PSZ2 G181.06+48.47, ASKAP J1832-0911, Sonification Collection, Centaurus A with CNB integration.
// This module can be plugged into a base program (e.g., 'uqff_buoyancy_sim.cpp') by including this header and linking the .cpp.
// Usage in base: #include "UQFFBuoyancyCNBModule.h"
// UQFFBuoyancyCNBModule mod; mod.computeFBi(system); mod.updateVariable("F_rel", {new_real, new_imag});
// All variables are stored in a std::map for dynamic addition/subtraction/update, using complex<double> for real/imaginary components.
// Nothing is negligible: Includes all terms - base force, momentum, gravity, vacuum stability, LENR resonance, activation, directed energy, magnetic resonance, neutron, relativistic, neutrino (CNB), Sweet vac, Kozima drop.
// Associated text: Outputs descriptive equation string via getEquationText().
// Approximations: Integral approximated as integrand * x2 (quadratic root); imag parts small and not fully scaled; LENR dominant due to low omega_0; x2 from quadratic solver approx; F_rel from 1998 LEP.
// Multi-system params: J1610+1811 M=2.785e30 kg r=3.09e15 m; PLCK G287.0+32.9 M=1.989e44 kg r=3.09e22 m; PSZ2 G181.06+48.47 M=1.989e44 kg r=3.09e22 m; ASKAP J1832-0911 M=2.785e30 kg r=4.63e16 m; Sonification Collection M=1.989e31 kg r=6.17e16 m; Centaurus A M=1.094e38 kg r=6.17e17 m.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 22, 2025.

#ifndef UQFF_BUOYANCY_CNB_MODULE_H
#define UQFF_BUOYANCY_CNB_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>
#include <sstream>

using cdouble = std::complex<double>;

class UQFFBuoyancyCNBModule {
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
    UQFFBuoyancyCNBModule();

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

#endif // UQFF_BUOYANCY_CNB_MODULE_H

// UQFFBuoyancyCNBModule.cpp
#include "UQFFBuoyancyCNBModule.h"
#include <complex>

// Constructor: Set all variables with multi-system defaults
UQFFBuoyancyCNBModule::UQFFBuoyancyCNBModule() {
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
    variables["F_rel"] = {4.30e33, 0.0};  // Relativistic coherence from LEP 1998
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
    variables["k_neutrino"] = {1e-10, 0.0};  // For CNB
    variables["k_Sweet"] = {1e-25, 0.0};
    variables["k_Kozima"] = {1e-18, 0.0};
    variables["omega_0_LENR"] = {2 * pi_val * 1.25e12, 0.0};  // LENR resonance freq 1.2-1.3 THz
    variables["k_LENR"] = {1e-10, 0.0};
    variables["rho_vac_UA"] = {7.09e-36, 0.0};  // Vacuum energy density
    variables["sigma_n"] = {1e-4, 0.0};  // Scaled for astrophysical densities
    variables["sigma_CNB"] = {1e-49, 0.0};  // CNB cross-section m^2
    variables["n_CNB"] = {3.36e8, 0.0};  // CNB number density m^-3
    variables["E_CNB"] = {2.69e-23, 0.0};  // CNB average energy J
    variables["E_cm"] = {189.0, 0.0};  // GeV
    variables["E_cm_astro_local_adj_eff_enhanced"] = {1.24e24, 0.0};  // events/m3
    variables["DPM_stability"] = {0.01, 0.0};
    variables["DPM_momentum"] = {0.93, 0.0};
    variables["DPM_gravity"] = {1.0, 0.0};
    variables["DPM_light"] = {0.01, 0.0};
    variables["phase"] = {2.36e-3, 0.0};  // s^-1
    variables["curvature"] = {1e-22, 0.0};
    variables["k_10_13"] = {1e-13, 0.0};  // For light term in quadratic
    variables["k_b_term"] = {2.51e-5, 0.0};  // Constant in b
    variables["c_constant"] = {-3.06e175, 0.0};  // Constant in c
    variables["c_inv_r2_coeff"] = {1e-29, 0.0};  // 10^{-29}/r^2 in c
    variables["a_eps_coeff"] = {1.38e-41, 0.0};  // Coefficient for first term in a

    // System-specific params will be set in setSystemParams()
}

// Set system-specific parameters
void UQFFBuoyancyCNBModule::setSystemParams(const std::string& system)
{
    double pi_val = variables["pi"].real();
    if (system == "J1610+1811") {
        this->variables["M"] = {2.785e30, 0.0};
        this->variables["r"] = {3.09e15, 0.0};
        this->variables["T"] = {1e4, 0.0};
        this->variables["L_X"] = {1e31, 0.0};
        this->variables["B0"] = {1e-4, 0.0};
        this->variables["omega0"] = {1e-12, 0.0};
        this->variables["Mach"] = {1.0, 0.0};
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
        this->variables["omega0"] = {1e-12, 0.0};
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
    } else if (system == "CentaurusA") {
        this->variables["M"] = {1.094e38, 0.0};
        this->variables["r"] = {6.17e17, 0.0};
        this->variables["T"] = {1e4, 0.0};
        this->variables["L_X"] = {1e36, 0.0};
        this->variables["B0"] = {1e-4, 0.0};
        this->variables["omega0"] = {1e-15, 0.0};
        this->variables["Mach"] = {1.5, 0.0};
        this->variables["C"] = {1.2, 0.0};
        this->variables["theta"] = {pi_val / 4, 0.0};
        this->variables["t"] = {3.472e14, 0.0};
    }
}

// Dynamic variable operations (complex)
void UQFFBuoyancyCNBModule::updateVariable(const std::string& name, cdouble value) {
    this->variables[name] = value;
}
void UQFFBuoyancyCNBModule::addToVariable(const std::string& name, cdouble delta) {
    this->variables[name] += delta;
}
void UQFFBuoyancyCNBModule::subtractFromVariable(const std::string& name, cdouble delta) {
    this->variables[name] -= delta;
}

// Core computation: Full F_U_Bi_i(r, t) for system (approx integral)
cdouble UQFFBuoyancyCNBModule::computeFBi(const std::string& system, double t) {
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
cdouble UQFFBuoyancyCNBModule::computeCompressed(const std::string& system, double t) {
    setSystemParams(system);
    return computeIntegrand(t, system);
}
cdouble UQFFBuoyancyCNBModule::computeResonant(const std::string& system) {
    setSystemParams(system);
    return computeDPM_resonance(system);
}
cdouble UQFFBuoyancyCNBModule::computeBuoyancy(const std::string& system) {
    setSystemParams(system);
    return computeUb1(system);
}
cdouble UQFFBuoyancyCNBModule::computeSuperconductive(const std::string& system, double t) {
    setSystemParams(system);
    return computeUi(t, system);
}
double UQFFBuoyancyCNBModule::computeCompressedG(const std::string& system, double t) {
    setSystemParams(system);
    return computeG(t, system);
}

// Output descriptive text of the equation
std::string UQFFBuoyancyCNBModule::getEquationText(const std::string& system) {
    setSystemParams(system);
    std::ostringstream oss;
    oss << "F_U_Bi_i(r, t) = Integral[Integrand(r, t) dt] approximated as Integrand * x2\n";
    oss << "Where Integrand includes terms for base force, momentum, gravity, vacuum stability, LENR resonance, activation, directed energy, magnetic resonance, neutron, relativistic, neutrino (CNB), Sweet vac, Kozima drop.\n";
    oss << "LENR Resonance: F_LENR = k_LENR * (omega_LENR / omega_0)^2\n";
    oss << "Activation: F_act = k_act * cos(omega_act t)\n";
    oss << "Directed Energy: F_DE = k_DE * L_X\n";
    oss << "Magnetic Resonance: F_res = 2 q B_0 V sin theta * DPM_resonance\n";
    oss << "Neutron Drop: F_neutron = k_neutron * sigma_n\n";
    oss << "Relativistic: F_rel = k_rel * (E_cm_astro_local_adj_eff_enhanced / E_cm)^2 = 4.30e33 N\n";
    oss << "CNB Neutrino: F_neutrino = k_neutrino * sigma_CNB * n_CNB * E_CNB approx 9.07e-42 N\n";
    oss << "Sweet Vac: F_sweet = k_Sweet * rho_vac_UA\n";
    oss << "Kozima Drop: F_kozima = k_Kozima * sigma_n\n";
    oss << "Relativistic Correction: F_relativ = k_relativistic * (V / c)^2 * F0\n";
    oss << "System: " << system << "\n";
    return oss.str();
}

// Print all current variables (for debugging/updates)
void UQFFBuoyancyCNBModule::printVariables() {
    for (const auto& pair : variables) {
        std::cout << std::setw(15) << pair.first << " : " << pair.second << std::endl;
    }
}

// Compute integrand for F_U_Bi_i
cdouble UQFFBuoyancyCNBModule::computeIntegrand(double t, const std::string& system) {
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
    cdouble f_neutrino = variables["k_neutrino"] * variables["sigma_CNB"] * variables["n_CNB"] * variables["E_CNB"];
    cdouble f_res = 2.0 * variables["q"].real() * variables["B0"].real() * variables["V"].real() * sin_theta * dpm_res;
    cdouble f_relativ = variables["k_relativistic"] * pow(variables["V"].real() / variables["c"].real(), 2.0) * variables["F0"];
    cdouble f_sweet = variables["k_Sweet"] * variables["rho_vac_UA"];
    cdouble f_kozima = variables["k_Kozima"] * variables["sigma_n"];
    cdouble momentum_term = (variables["m_e"] * variables["c"] * variables["c"] / (variables["r"] * variables["r"])) * variables["DPM_momentum"] * cos_theta;
    cdouble gravity_term = (variables["G"] * variables["M"] / (variables["r"] * variables["r"])) * variables["DPM_gravity"];
    cdouble vac_term = variables["rho_vac_UA"] * variables["DPM_stability"];
    cdouble integrand = -variables["F0"] + momentum_term + gravity_term + vac_term + f_lenr + f_act + f_de + f_res + f_neutron + f_rel + f_neutrino + f_relativ + f_sweet + f_kozima;
    return integrand;
}

// Compute DPM resonance term
cdouble UQFFBuoyancyCNBModule::computeDPM_resonance(const std::string& system) {
    setSystemParams(system);
    double g_lande = variables["g_Lande"].real();
    double mu_b = variables["mu_B"].real();
    double b0 = variables["B0"].real();
    double hbar_omega0 = variables["hbar"].real() * variables["omega0"].real();
    if (hbar_omega0 == 0.0) return {0.0, 0.0};
    return {g_lande * mu_b * b0 / hbar_omega0, 0.0};
}

// Compute x2 from quadratic root approximation (negative root as per doc)
cdouble UQFFBuoyancyCNBModule::computeX2(const std::string& system) {
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
cdouble UQFFBuoyancyCNBModule::computeQuadraticRoot(cdouble a, cdouble b, cdouble c) {
    cdouble disc = b * b - 4.0 * a * c;
    double disc_real = disc.real();
    if (disc_real < 0.0) disc_real = 0.0;  // Force real for approximation
    cdouble sqrt_disc = {sqrt(disc_real), 0.0};
    return (-b - sqrt_disc) / (2.0 * a);
}

// Compute LENR term
cdouble UQFFBuoyancyCNBModule::computeLENRTerm(const std::string& system) {
    setSystemParams(system);
    double omega0_real = variables["omega0"].real();
    if (omega0_real == 0.0) return {0.0, 0.0};
    return variables["k_LENR"] * pow(variables["omega_0_LENR"].real() / omega0_real, 2.0);
}

// Compute gravitational acceleration g(r,t) - as per document
double UQFFBuoyancyCNBModule::computeG(double t, const std::string& system) {
    return -1.07e16;
}

// Compute Q_wave term - as per document
cdouble UQFFBuoyancyCNBModule::computeQ_wave(double t, const std::string& system) {
    return {3.11e5, 0.0};
}

// Compute Ub1 buoyancy term
cdouble UQFFBuoyancyCNBModule::computeUb1(const std::string& system) {
    return computeIntegrand(0.0, system);
}

// Compute Ui superconductive term
cdouble UQFFBuoyancyCNBModule::computeUi(double t, const std::string& system) {
    return computeQ_wave(t, system);
}

// ===== CNB-ENHANCED BUOYANCY FUNCTIONS =====

// Enhanced Ub1 with CNB coupling
cdouble UQFFBuoyancyCNBModule::computeUb1_CNB(const std::string& system) {
    setSystemParams(system);
    
    // CNB parameters
    cdouble sigma_CNB = variables["sigma_CNB"];
    cdouble n_CNB = variables["n_CNB"];
    cdouble E_CNB = variables["E_CNB"];
    cdouble rho_vac = variables["rho_vac_UA"];
    cdouble c = variables["c"];
    
    // CNB enhancement factor
    cdouble CNB_coupling = sigma_CNB * n_CNB * E_CNB;
    cdouble base_buoyancy = rho_vac * c * c / 3.0;
    cdouble enhancement = 1.0 + CNB_coupling / (rho_vac * c * c);
    
    // System scaling
    cdouble M = variables["M"];
    cdouble r = variables["r"];
    cdouble scaling = std::sqrt(M / (r * r * r));
    
    return base_buoyancy * enhancement * scaling;
}

// Enhanced Ui with CNB integration  
cdouble UQFFBuoyancyCNBModule::computeUi_CNB(double t, const std::string& system) {
    setSystemParams(system);
    
    // CNB superconductive interaction
    cdouble sigma_CNB = variables["sigma_CNB"];
    cdouble n_CNB = variables["n_CNB"];
    cdouble E_CNB = variables["E_CNB"];
    cdouble k_neutrino = variables["k_neutrino"];
    cdouble c = variables["c"];
    
    // CNB flux and neutrino coupling
    cdouble CNB_flux = n_CNB * c;
    cdouble coupling = k_neutrino * sigma_CNB * E_CNB;
    
    // Quantum coherence
    cdouble hbar = variables["hbar"];
    cdouble omega = 2.0 * variables["pi"] * E_CNB / hbar;
    cdouble coherence = std::exp(cdouble(0.0, omega.real() * t));
    
    // Mass scaling
    cdouble M = variables["M"];
    cdouble r = variables["r"];
    cdouble mass_scale = M / (4.0 * variables["pi"] * r * r * r);
    
    return CNB_flux * coupling * coherence * mass_scale;
}

// ===== SURFACE MAGNETIC FIELD MODULE FOR CNB SYSTEMS =====

#ifndef SURFACE_MAGNETIC_FIELD_MODULE_CNB_H
#define SURFACE_MAGNETIC_FIELD_MODULE_CNB_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

class SurfaceMagneticFieldModule {
private:
    std::map<std::string, double> variables;
    std::map<std::string, std::vector<double>> variable_history;
    std::map<std::string, std::string> variable_dependencies;
    bool self_learning_enabled;
    double learning_rate;
    int update_counter;
    
    // Dynamic helper functions
    void updateDependencies(const std::string& changed_var);
    double computeGradient(const std::string& var, const std::string& target);
    void recordHistory(const std::string& name, double value);

public:
    // Constructor with CNB-enhanced dynamic capabilities
    SurfaceMagneticFieldModule();
    
    // Core magnetic field computations for CNB systems
    double computeB_j(double t, double B_s);
    double computeB_s_min();
    double computeB_s_max();
    double computeU_g3_example(double t, double B_s);
    double computeCNB_MagneticCoupling(double B_field, double CNB_flux);
    
    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    
    // Advanced dynamic capabilities for CNB systems
    void autoCalibrate(const std::string& observable, double target_value, double tolerance = 0.01);
    void adaptiveUpdate(double dt, const std::string& feedback_param = "");
    void scaleToCNBData(const std::map<std::string, double>& cnb_data);
    void addCustomVariable(const std::string& name, double value, const std::string& dependency = "");
    std::map<std::string, double> getVariableHistory(const std::string& name, int steps = 10);
    void enableSelfLearning(bool enable);
    void exportState(const std::string& filename);
    void importState(const std::string& filename);
    
    // Enhanced magnetic field equations for CNB systems
    std::string getEquationText();
};

#endif // SURFACE_MAGNETIC_FIELD_MODULE_CNB_H

// ===== SURFACE MAGNETIC FIELD MODULE IMPLEMENTATION FOR CNB SYSTEMS =====

// Enhanced SurfaceMagneticFieldModule constructor with CNB capabilities
SurfaceMagneticFieldModule::SurfaceMagneticFieldModule() {
    // Initialize dynamic capabilities
    self_learning_enabled = false;
    learning_rate = 0.08;  // Higher learning rate for CNB magnetic field variability
    update_counter = 0;
    
    // Universal constants
    variables["B_s_min"] = 1e-6;                    // T (quiet CNB systems)
    variables["B_s_max"] = 1.2;                     // T (active CNB systems)
    variables["B_ref"] = 1.2;                       // T (reference max for CNB)
    variables["k_3"] = 2.5;                         // Enhanced coupling for CNB
    variables["omega_s"] = 1.8e-5;                  // rad/s (CNB system frequency)
    variables["P_core"] = 1.0;                      // Unitless
    variables["E_react"] = 5e46;                    // J (enhanced for CNB systems)
    variables["pi"] = 3.141592653589793;
    variables["t"] = 0.0;                           // s
    
    // CNB-specific parameters
    variables["CNB_temperature"] = 1.95;            // K
    variables["CNB_density"] = 4e8;                 // m^-3
    variables["CNB_coupling"] = 2.3e-12;           // CNB-magnetic coupling
    variables["magnetic_diffusion"] = 5e-11;       // m/s (enhanced for CNB)
    variables["convection_velocity"] = 2.5e3;      // m/s (higher for active systems)
    variables["neutrino_flux"] = 3.3e10;           // m^-2 s^-1
    
    // System evolution parameters
    variables["evolution_timescale"] = 1e15;       // s (CNB evolution timescale)
    variables["thermal_coupling"] = 1.5e-8;        // Thermal-magnetic coupling
}

// Compute minimum surface magnetic field for CNB systems
double SurfaceMagneticFieldModule::computeB_s_min() {
    return variables["B_s_min"];
}

// Compute maximum surface magnetic field for CNB systems
double SurfaceMagneticFieldModule::computeB_s_max() {
    return variables["B_s_max"];
}

// Compute scaled B_j based on time t and surface field B_s for CNB systems
double SurfaceMagneticFieldModule::computeB_j(double t, double B_s) {
    variables["t"] = t;
    
    // Enhanced magnetic field evolution with CNB coupling
    double cnb_oscillation = 0.6 * std::sin(variables["omega_s"] * t);
    double thermal_factor = 1.0 + variables["thermal_coupling"] * variables["CNB_temperature"];
    double base_b = variables["B_ref"] + cnb_oscillation * thermal_factor;
    
    return base_b * (B_s / variables["B_ref"]);
}

// Compute CNB-magnetic field coupling
double SurfaceMagneticFieldModule::computeCNB_MagneticCoupling(double B_field, double CNB_flux) {
    double coupling_strength = variables["CNB_coupling"];
    double density_factor = variables["CNB_density"] / 4e8;  // Normalized to standard CNB density
    double temperature_factor = variables["CNB_temperature"] / 1.95;  // Normalized to CNB temperature
    
    return coupling_strength * B_field * CNB_flux * density_factor * temperature_factor;
}

// Compute U_g3 example with CNB enhancement
double SurfaceMagneticFieldModule::computeU_g3_example(double t, double B_s) {
    double k_3 = variables["k_3"];
    double b_j = computeB_j(t, B_s);
    double cos_term = std::cos(variables["omega_s"] * t * variables["pi"]);
    double p_core = variables["P_core"];
    double e_react = variables["E_react"];
    
    // CNB enhancement factor
    double cnb_enhancement = 1.0 + computeCNB_MagneticCoupling(b_j, variables["neutrino_flux"]);
    
    return k_3 * b_j * cos_term * p_core * e_react * cnb_enhancement;
}

// Enhanced updateVariable with CNB dynamic capabilities
void SurfaceMagneticFieldModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
    
    // Record history for dynamic capabilities
    recordHistory(name, value);
    
    // Update dependencies
    updateDependencies(name);
    
    // Increment update counter
    update_counter++;
    
    // Trigger self-learning if enabled
    if (self_learning_enabled && update_counter % 4 == 0) {  // More frequent for CNB systems
        adaptiveUpdate(1.0, name);
    }
}

// Auto-calibrate magnetic field parameters for CNB systems
void SurfaceMagneticFieldModule::autoCalibrate(const std::string& observable, double target_value, double tolerance) {
    if (variables.find(observable) == variables.end()) {
        std::cerr << "Observable '" << observable << "' not found for calibration." << std::endl;
        return;
    }
    
    double current_value = variables[observable];
    double error = std::abs(current_value - target_value) / target_value;
    
    if (error > tolerance) {
        // CNB-specific parameter adjustment
        std::vector<std::string> tunable_params = {"B_ref", "k_3", "omega_s", "P_core", "CNB_coupling", "thermal_coupling"};
        
        for (const auto& param : tunable_params) {
            double gradient = computeGradient(param, observable);
            if (std::abs(gradient) > 1e-25) {
                double adjustment = learning_rate * (target_value - current_value) / gradient;
                variables[param] += adjustment;
                recordHistory(param, variables[param]);
            }
        }
        
        std::cout << "Auto-calibrated CNB magnetic " << observable << " from " << current_value 
                  << " to target " << target_value << " (error: " << error << ")" << std::endl;
    }
}

// Adaptive magnetic field evolution for CNB systems
void SurfaceMagneticFieldModule::adaptiveUpdate(double dt, const std::string& feedback_param) {
    if (!self_learning_enabled) return;
    
    // CNB evolution timescale
    double evolution_factor = std::exp(-dt / variables["evolution_timescale"]);
    
    // CNB temperature evolution
    variables["CNB_temperature"] *= (1.0 + 0.0001 * std::sin(dt / 1e12));
    
    // Adaptive magnetic field reference with CNB coupling
    double cnb_factor = variables["CNB_density"] / 4e8;
    variables["B_ref"] = variables["B_s_max"] * (0.2 + 0.8 * cnb_factor);
    
    // Enhanced magnetic diffusion effects for CNB systems
    double diffusion_decay = std::exp(-dt * variables["magnetic_diffusion"] / 1e4);
    variables["k_3"] *= diffusion_decay;
    
    // CNB-driven convection enhancement
    double convection_enhancement = 1.0 + 0.15 * variables["convection_velocity"] / 2.5e3;
    variables["omega_s"] *= convection_enhancement;
    
    // Neutrino flux variability
    variables["neutrino_flux"] *= (1.0 + 0.001 * std::cos(dt / 1e10));
    
    recordHistory("adaptive_time", dt);
    std::cout << "CNB magnetic adaptive update: B_ref=" << variables["B_ref"] 
              << ", CNB_T=" << variables["CNB_temperature"] << std::endl;
}

// Scale to CNB observational data
void SurfaceMagneticFieldModule::scaleToCNBData(const std::map<std::string, double>& cnb_data) {
    for (const auto& data : cnb_data) {
        if (data.first == "cnb_temperature" && variables.find("CNB_temperature") != variables.end()) {
            double scaling = data.second / variables["CNB_temperature"];
            variables["CNB_temperature"] = data.second;
            variables["thermal_coupling"] *= scaling;
        }
        
        if (data.first == "neutrino_flux") {
            variables["neutrino_flux"] = data.second;
            variables["CNB_coupling"] *= std::sqrt(data.second / 3.3e10);
        }
        
        if (data.first == "magnetic_field_strength") {
            double scaling = data.second / variables["B_s_max"];
            variables["B_s_max"] = data.second;
            variables["B_ref"] *= scaling;
        }
    }
    std::cout << "Scaled CNB magnetic module to " << cnb_data.size() << " CNB observations." << std::endl;
}

// Add custom magnetic variables for CNB systems
void SurfaceMagneticFieldModule::addCustomVariable(const std::string& name, double value, const std::string& dependency) {
    variables[name] = value;
    if (!dependency.empty()) {
        variable_dependencies[name] = dependency;
    }
    recordHistory(name, value);
    std::cout << "Added custom CNB magnetic variable: " << name << " = " << value << std::endl;
}

// Get magnetic parameter history for CNB systems
std::map<std::string, double> SurfaceMagneticFieldModule::getVariableHistory(const std::string& name, int steps) {
    std::map<std::string, double> history;
    if (variable_history.find(name) != variable_history.end()) {
        auto& hist = variable_history[name];
        int start = std::max(0, (int)hist.size() - steps);
        for (int i = start; i < (int)hist.size(); i++) {
            history["step_" + std::to_string(i)] = hist[i];
        }
    }
    return history;
}

// Enable CNB magnetic self-learning
void SurfaceMagneticFieldModule::enableSelfLearning(bool enable) {
    self_learning_enabled = enable;
    if (enable) {
        std::cout << "CNB magnetic self-learning enabled with rate: " << learning_rate << std::endl;
    } else {
        std::cout << "CNB magnetic self-learning disabled." << std::endl;
    }
}

// Export CNB magnetic state
void SurfaceMagneticFieldModule::exportState(const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        file << "# SurfaceMagneticFieldModule CNB State Export" << std::endl;
        file << "update_counter=" << update_counter << std::endl;
        file << "learning_rate=" << learning_rate << std::endl;
        file << "self_learning_enabled=" << (self_learning_enabled ? 1 : 0) << std::endl;
        
        for (const auto& var : variables) {
            file << var.first << "=" << var.second << std::endl;
        }
        file.close();
        std::cout << "CNB magnetic state exported to: " << filename << std::endl;
    }
}

// Import CNB magnetic state
void SurfaceMagneticFieldModule::importState(const std::string& filename) {
    std::ifstream file(filename);
    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            if (line[0] == '#') continue;
            
            size_t eq_pos = line.find('=');
            if (eq_pos != std::string::npos) {
                std::string key = line.substr(0, eq_pos);
                std::string value_str = line.substr(eq_pos + 1);
                
                if (key == "update_counter") {
                    update_counter = std::stoi(value_str);
                } else if (key == "learning_rate") {
                    learning_rate = std::stod(value_str);
                } else if (key == "self_learning_enabled") {
                    self_learning_enabled = (std::stoi(value_str) == 1);
                } else {
                    variables[key] = std::stod(value_str);
                }
            }
        }
        file.close();
        std::cout << "CNB magnetic state imported from: " << filename << std::endl;
    }
}

// Enhanced equation text for CNB systems
std::string SurfaceMagneticFieldModule::getEquationText() {
    return "CNB-Enhanced Magnetic Field Equations:\n"
           "B_j  (B_ref + 0.6 sin(?_s t) * T_thermal) * (B_s / B_ref) T\n"
           "U_g3 = k_3 *  B_j * cos(?_s t p) * P_core * E_react * (1 + CNB_coupling)\n"
           "CNB_coupling = ?_CNB * B_field * ?_flux * (?_CNB/?_0) * (T_CNB/T_0)\n"
           "Where:\n"
           "- B_s = [1e-6, 1.2] T (CNB system range)\n"
           "- T_CNB = " + std::to_string(variables["CNB_temperature"]) + " K (cosmic neutrino background temperature)\n"
           "- ?_flux = " + std::to_string(variables["neutrino_flux"]) + " m s (neutrino flux)\n"
           "- ?_CNB = " + std::to_string(variables["CNB_coupling"]) + " (CNB-magnetic coupling)\n"
           "CNB Systems: J1610+1811, PLCK G287.0+32.9, PSZ2 G181.06+48.47,\n"
           "ASKAP J1832-0911, Sonification Collection, Centaurus A\n"
           "Enhanced Features: Thermal coupling, neutrino flux modulation,\n"
           "adaptive CNB temperature evolution, enhanced diffusion rates.";
}

// Helper functions for CNB magnetic field module
void SurfaceMagneticFieldModule::updateDependencies(const std::string& changed_var) {
    if (changed_var == "CNB_temperature") {
        // Update thermal coupling based on CNB temperature
        variables["thermal_coupling"] = 1.5e-8 * (variables["CNB_temperature"] / 1.95);
    }
    
    if (changed_var == "B_s_max") {
        variables["B_ref"] = variables["B_s_max"];
    }
    
    if (changed_var == "neutrino_flux") {
        // Update CNB coupling based on neutrino flux
        variables["CNB_coupling"] = 2.3e-12 * std::sqrt(variables["neutrino_flux"] / 3.3e10);
    }
    
    if (changed_var == "omega_s") {
        // Update evolution timescale
        variables["evolution_timescale"] = 2 * M_PI / variables["omega_s"] * 1e9;
    }
}

double SurfaceMagneticFieldModule::computeGradient(const std::string& var, const std::string& target) {
    if (variables.find(var) == variables.end() || variables.find(target) == variables.end()) {
        return 0.0;
    }
    
    double original_value = variables[var];
    double original_target = variables[target];
    
    // Small perturbation
    double delta = original_value * 1e-6;
    variables[var] += delta;
    
    // Recompute target with CNB coupling
    double new_target = computeB_j(variables["t"], variables["B_ref"]);
    
    // Restore original value
    variables[var] = original_value;
    
    return (new_target - original_target) / delta;
}

void SurfaceMagneticFieldModule::recordHistory(const std::string& name, double value) {
    variable_history[name].push_back(value);
    
    // Keep only last 150 values for CNB systems (more history)
    if (variable_history[name].size() > 150) {
        variable_history[name].erase(variable_history[name].begin());
    }
}
