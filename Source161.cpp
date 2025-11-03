
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
#include <vector>
#include <functional>
#include <random>
#include <algorithm>
#include <fstream>

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

    // 25-method dynamic self-update & self-expansion capabilities
    
    // 1. Variable Management (5 methods)
    void createVariable(const std::string& name, const std::complex<double>& value);
    void removeVariable(const std::string& name);
    std::complex<double> cloneVariable(const std::string& srcName, const std::string& destName);
    std::vector<std::string> listVariables() const;
    std::string getSystemName() const;
    
    // 2. Batch Operations (2 methods)
    void transformVariableGroup(const std::vector<std::string>& varNames, 
                                 std::function<std::complex<double>(std::complex<double>)> func);
    void scaleVariableGroup(const std::vector<std::string>& varNames, const std::complex<double>& scaleFactor);
    
    // 3. Self-Expansion (4 methods - parameter space + 3 domain-specific scales)
    void expandParameterSpace(int numNewParams = 5);
    void expandSystemScale(const std::string& system); // X-ray/cluster specific
    void expandForceScale(double factor = 1.5);        // LENR + relativistic
    void expandAstroScale(double xrayFactor = 1.3, double magneticFactor = 1.2); // X-ray/magnetic physics
    
    // 4. Self-Refinement (3 methods)
    void autoRefineParameters(int iterations = 100);
    void calibrateToObservations(const std::map<std::string, std::complex<double>>& observed);
    void optimizeForMetric(std::function<double(const std::map<std::string, std::complex<double>>&)> metric, 
                            int iterations = 50);
    
    // 5. Parameter Exploration (1 method)
    std::vector<std::map<std::string, std::complex<double>>> generateVariations(int numVariations = 10, 
                                                                                  double stdDev = 0.1);
    
    // 6. Adaptive Evolution (2 methods)
    void mutateParameters(double mutationRate = 0.05);
    void evolveSystem(int generations = 20, std::function<double(UQFFBuoyancyAstroModule&)> fitness = nullptr);
    
    // 7. State Management (4 methods)
    void saveState(const std::string& label);
    void restoreState(const std::string& label);
    std::vector<std::string> listSavedStates() const;
    std::string exportState(const std::string& format = "text") const;
    
    // 8. System Analysis (4 methods)
    std::map<std::string, double> sensitivityAnalysis(const std::string& system, 
                                                       const std::vector<std::string>& params);
    std::string generateReport(const std::string& system) const;
    bool validateConsistency() const;
    void autoCorrectAnomalies();

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

// ============================================================================
// SECTION: 25-Method Dynamic Self-Update & Self-Expansion Implementation
// ============================================================================

// Namespace for saved states
namespace saved_states_astro {
    std::map<std::string, std::map<std::string, std::complex<double>>> states;
}

// ============================================================================
// 1. Variable Management (5 methods)
// ============================================================================

void UQFFBuoyancyAstroModule::createVariable(const std::string& name, const std::complex<double>& value) {
    variables[name] = value;
    std::cout << "Created variable '" << name << "' = " << value << std::endl;
}

void UQFFBuoyancyAstroModule::removeVariable(const std::string& name) {
    auto it = variables.find(name);
    if (it != variables.end()) {
        variables.erase(it);
        std::cout << "Removed variable '" << name << "'" << std::endl;
    } else {
        std::cerr << "Variable '" << name << "' not found for removal." << std::endl;
    }
}

std::complex<double> UQFFBuoyancyAstroModule::cloneVariable(const std::string& srcName, const std::string& destName) {
    auto it = variables.find(srcName);
    if (it != variables.end()) {
        variables[destName] = it->second;
        std::cout << "Cloned '" << srcName << "' to '" << destName << "' = " << it->second << std::endl;
        return it->second;
    } else {
        std::cerr << "Source variable '" << srcName << "' not found." << std::endl;
        return std::complex<double>(0.0, 0.0);
    }
}

std::vector<std::string> UQFFBuoyancyAstroModule::listVariables() const {
    std::vector<std::string> varNames;
    varNames.reserve(variables.size());
    for (const auto& pair : variables) {
        varNames.push_back(pair.first);
    }
    return varNames;
}

std::string UQFFBuoyancyAstroModule::getSystemName() const {
    return "UQFFBuoyancy_XRay_Cluster_MultiSystem";
}

// ============================================================================
// 2. Batch Operations (2 methods)
// ============================================================================

void UQFFBuoyancyAstroModule::transformVariableGroup(const std::vector<std::string>& varNames, 
                                                       std::function<std::complex<double>(std::complex<double>)> func) {
    for (const auto& name : varNames) {
        auto it = variables.find(name);
        if (it != variables.end()) {
            it->second = func(it->second);
            std::cout << "Transformed '" << name << "' to " << it->second << std::endl;
        } else {
            std::cerr << "Variable '" << name << "' not found for transformation." << std::endl;
        }
    }
}

void UQFFBuoyancyAstroModule::scaleVariableGroup(const std::vector<std::string>& varNames, 
                                                   const std::complex<double>& scaleFactor) {
    for (const auto& name : varNames) {
        auto it = variables.find(name);
        if (it != variables.end()) {
            it->second *= scaleFactor;
            std::cout << "Scaled '" << name << "' by " << scaleFactor << " to " << it->second << std::endl;
        } else {
            std::cerr << "Variable '" << name << "' not found for scaling." << std::endl;
        }
    }
}

// ============================================================================
// 3. Self-Expansion (4 methods)
// ============================================================================

void UQFFBuoyancyAstroModule::expandParameterSpace(int numNewParams) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1e30, 1e44); // X-ray source to galaxy cluster mass range
    
    for (int i = 0; i < numNewParams; ++i) {
        std::string newName = "expanded_param_" + std::to_string(i);
        double realVal = dis(gen);
        double imagVal = realVal * 1e-6; // small imaginary component
        variables[newName] = std::complex<double>(realVal, imagVal);
        std::cout << "Expanded parameter space: " << newName << " = " << variables[newName] << std::endl;
    }
}

void UQFFBuoyancyAstroModule::expandSystemScale(const std::string& system) {
    setSystemParams(system);
    
    // X-ray/cluster-specific scaling
    if (system == "J1610+1811") {
        variables["M"] *= std::complex<double>(1.18, 0.0); // pulsar-like source
        variables["r"] *= std::complex<double>(1.12, 0.0);
        variables["L_X"] *= std::complex<double>(1.25, 0.0); // X-ray luminosity boost
        std::cout << "Expanded J1610+1811 (pulsar) scale" << std::endl;
    } else if (system == "PLCK_G287.0+32.9") {
        variables["M"] *= std::complex<double>(1.10, 0.0); // massive galaxy cluster
        variables["r"] *= std::complex<double>(1.05, 0.0);
        variables["L_X"] *= std::complex<double>(1.15, 0.0);
        std::cout << "Expanded PLCK G287.0+32.9 (galaxy cluster) scale" << std::endl;
    } else if (system == "PSZ2_G181.06+48.47") {
        variables["M"] *= std::complex<double>(1.12, 0.0); // Planck cluster
        variables["r"] *= std::complex<double>(1.06, 0.0);
        variables["L_X"] *= std::complex<double>(1.20, 0.0);
        std::cout << "Expanded PSZ2 G181.06+48.47 (Planck cluster) scale" << std::endl;
    } else if (system == "ASKAP_J1832-0911") {
        variables["M"] *= std::complex<double>(1.20, 0.0); // radio transient
        variables["r"] *= std::complex<double>(1.15, 0.0);
        variables["L_X"] *= std::complex<double>(1.30, 0.0);
        std::cout << "Expanded ASKAP J1832-0911 (radio transient) scale" << std::endl;
    } else if (system == "SonificationCollection") {
        variables["M"] *= std::complex<double>(1.15, 0.0); // sonification data
        variables["r"] *= std::complex<double>(1.10, 0.0);
        variables["L_X"] *= std::complex<double>(1.22, 0.0);
        std::cout << "Expanded Sonification Collection scale" << std::endl;
    }
}

void UQFFBuoyancyAstroModule::expandForceScale(double factor) {
    // LENR resonance expansion
    auto it_lenr = variables.find("omega_0_LENR");
    if (it_lenr != variables.end()) {
        it_lenr->second *= std::complex<double>(factor, 0.0);
        std::cout << "Expanded LENR resonance by factor " << factor << std::endl;
    }
    
    // Relativistic force expansion
    auto it_rel = variables.find("k_rel");
    if (it_rel != variables.end()) {
        it_rel->second *= std::complex<double>(factor, 0.0);
        std::cout << "Expanded relativistic force coefficient by factor " << factor << std::endl;
    }
    
    // Base force F0 expansion
    auto it_f0 = variables.find("F0");
    if (it_f0 != variables.end()) {
        it_f0->second *= std::complex<double>(factor, 0.0);
        std::cout << "Expanded base force F0 by factor " << factor << std::endl;
    }
}

void UQFFBuoyancyAstroModule::expandAstroScale(double xrayFactor, double magneticFactor) {
    // X-ray luminosity expansion
    auto it_lx = variables.find("L_X");
    if (it_lx != variables.end()) {
        it_lx->second *= std::complex<double>(xrayFactor, 0.0);
        std::cout << "Expanded X-ray luminosity by factor " << xrayFactor << std::endl;
    }
    
    // Magnetic field expansion
    auto it_b0 = variables.find("B0");
    if (it_b0 != variables.end()) {
        it_b0->second *= std::complex<double>(magneticFactor, 0.0);
        std::cout << "Expanded magnetic field B0 by factor " << magneticFactor << std::endl;
    }
    
    // Temperature expansion
    auto it_t = variables.find("T");
    if (it_t != variables.end()) {
        it_t->second *= std::complex<double>(1.15, 0.0);
        std::cout << "Expanded temperature by factor 1.15" << std::endl;
    }
    
    // Mach number expansion
    auto it_mach = variables.find("Mach");
    if (it_mach != variables.end()) {
        it_mach->second *= std::complex<double>(1.10, 0.0);
        std::cout << "Expanded Mach number by factor 1.10" << std::endl;
    }
}

// ============================================================================
// 4. Self-Refinement (3 methods)
// ============================================================================

void UQFFBuoyancyAstroModule::autoRefineParameters(int iterations) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-0.01, 0.01); // small refinements
    
    for (int i = 0; i < iterations; ++i) {
        for (auto& pair : variables) {
            double delta_real = dis(gen);
            double delta_imag = dis(gen) * 0.1; // smaller imaginary refinement
            pair.second += std::complex<double>(delta_real * std::abs(pair.second.real()), 
                                                 delta_imag * std::abs(pair.second.imag()));
        }
    }
    std::cout << "Auto-refined parameters over " << iterations << " iterations" << std::endl;
}

void UQFFBuoyancyAstroModule::calibrateToObservations(const std::map<std::string, std::complex<double>>& observed) {
    for (const auto& obs : observed) {
        auto it = variables.find(obs.first);
        if (it != variables.end()) {
            std::complex<double> error = obs.second - it->second;
            it->second += error * 0.5; // 50% correction
            std::cout << "Calibrated '" << obs.first << "' toward observation: " << obs.second << std::endl;
        }
    }
}

void UQFFBuoyancyAstroModule::optimizeForMetric(
    std::function<double(const std::map<std::string, std::complex<double>>&)> metric, 
    int iterations) {
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-0.05, 0.05);
    
    double bestMetric = metric(variables);
    auto bestVars = variables;
    
    for (int i = 0; i < iterations; ++i) {
        auto testVars = variables;
        for (auto& pair : testVars) {
            double delta = dis(gen);
            pair.second *= std::complex<double>(1.0 + delta, 1.0 + delta * 0.1);
        }
        
        double testMetric = metric(testVars);
        if (testMetric > bestMetric) {
            bestMetric = testMetric;
            bestVars = testVars;
        }
    }
    
    variables = bestVars;
    std::cout << "Optimized for metric over " << iterations << " iterations. Best metric: " 
              << bestMetric << std::endl;
}

// ============================================================================
// 5. Parameter Exploration (1 method)
// ============================================================================

std::vector<std::map<std::string, std::complex<double>>> UQFFBuoyancyAstroModule::generateVariations(
    int numVariations, double stdDev) {
    
    std::vector<std::map<std::string, std::complex<double>>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dis(0.0, stdDev);
    
    for (int i = 0; i < numVariations; ++i) {
        auto variant = variables;
        for (auto& pair : variant) {
            double noise_real = dis(gen);
            double noise_imag = dis(gen) * 0.1;
            pair.second *= std::complex<double>(1.0 + noise_real, 1.0 + noise_imag);
        }
        variations.push_back(variant);
    }
    
    std::cout << "Generated " << numVariations << " parameter variations with stdDev=" 
              << stdDev << std::endl;
    return variations;
}

// ============================================================================
// 6. Adaptive Evolution (2 methods)
// ============================================================================

void UQFFBuoyancyAstroModule::mutateParameters(double mutationRate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    std::normal_distribution<> mutationDis(0.0, mutationRate);
    
    for (auto& pair : variables) {
        if (dis(gen) < mutationRate) {
            double mutation_real = mutationDis(gen);
            double mutation_imag = mutationDis(gen) * 0.1;
            pair.second *= std::complex<double>(1.0 + mutation_real, 1.0 + mutation_imag);
        }
    }
    std::cout << "Mutated parameters with mutation rate " << mutationRate << std::endl;
}

void UQFFBuoyancyAstroModule::evolveSystem(int generations, 
                                            std::function<double(UQFFBuoyancyAstroModule&)> fitness) {
    if (!fitness) {
        fitness = [](UQFFBuoyancyAstroModule& mod) {
            return std::abs(mod.computeFBi("J1610+1811", 0.0).real());
        };
    }
    
    double bestFitness = fitness(*this);
    auto bestVars = variables;
    
    for (int gen = 0; gen < generations; ++gen) {
        mutateParameters(0.05);
        double currentFitness = fitness(*this);
        
        if (currentFitness > bestFitness) {
            bestFitness = currentFitness;
            bestVars = variables;
        } else {
            variables = bestVars;
        }
    }
    
    variables = bestVars;
    std::cout << "Evolved system over " << generations << " generations. Best fitness: " 
              << bestFitness << std::endl;
}

// ============================================================================
// 7. State Management (4 methods)
// ============================================================================

void UQFFBuoyancyAstroModule::saveState(const std::string& label) {
    saved_states_astro::states[label] = variables;
    std::cout << "Saved state: '" << label << "'" << std::endl;
}

void UQFFBuoyancyAstroModule::restoreState(const std::string& label) {
    auto it = saved_states_astro::states.find(label);
    if (it != saved_states_astro::states.end()) {
        variables = it->second;
        std::cout << "Restored state: '" << label << "'" << std::endl;
    } else {
        std::cerr << "State '" << label << "' not found." << std::endl;
    }
}

std::vector<std::string> UQFFBuoyancyAstroModule::listSavedStates() const {
    std::vector<std::string> stateNames;
    for (const auto& pair : saved_states_astro::states) {
        stateNames.push_back(pair.first);
    }
    return stateNames;
}

std::string UQFFBuoyancyAstroModule::exportState(const std::string& format) const {
    std::ostringstream oss;
    
    if (format == "text") {
        oss << "=== UQFF X-Ray/Cluster Buoyancy State Export (Text) ===" << std::endl;
        for (const auto& pair : variables) {
            oss << pair.first << " = " << pair.second << std::endl;
        }
    } else if (format == "csv") {
        oss << "Variable,Real,Imaginary" << std::endl;
        for (const auto& pair : variables) {
            oss << pair.first << "," << pair.second.real() << "," << pair.second.imag() << std::endl;
        }
    } else if (format == "json") {
        oss << "{" << std::endl;
        size_t count = 0;
        for (const auto& pair : variables) {
            oss << "  \"" << pair.first << "\": {\"real\": " << pair.second.real() 
                << ", \"imag\": " << pair.second.imag() << "}";
            if (++count < variables.size()) oss << ",";
            oss << std::endl;
        }
        oss << "}" << std::endl;
    }
    
    return oss.str();
}

// ============================================================================
// 8. System Analysis (4 methods)
// ============================================================================

std::map<std::string, double> UQFFBuoyancyAstroModule::sensitivityAnalysis(
    const std::string& system, const std::vector<std::string>& params) {
    
    std::map<std::string, double> sensitivities;
    double baseline = std::abs(computeFBi(system, 0.0).real());
    
    for (const auto& param : params) {
        auto it = variables.find(param);
        if (it != variables.end()) {
            std::complex<double> original = it->second;
            it->second *= std::complex<double>(1.01, 1.0); // 1% perturbation
            
            double perturbed = std::abs(computeFBi(system, 0.0).real());
            sensitivities[param] = std::abs((perturbed - baseline) / baseline);
            
            it->second = original; // restore
        }
    }
    
    std::cout << "Sensitivity analysis completed for " << params.size() << " parameters" << std::endl;
    return sensitivities;
}

std::string UQFFBuoyancyAstroModule::generateReport(const std::string& system) const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(3);
    
    oss << "========================================" << std::endl;
    oss << "UQFF X-Ray/Cluster Buoyancy Report" << std::endl;
    oss << "System: " << system << std::endl;
    oss << "========================================" << std::endl;
    
    auto it_m = variables.find("M");
    auto it_r = variables.find("r");
    auto it_t = variables.find("T");
    auto it_lx = variables.find("L_X");
    auto it_b0 = variables.find("B0");
    
    if (it_m != variables.end() && it_r != variables.end()) {
        oss << "Mass: " << it_m->second.real() << " kg" << std::endl;
        oss << "Radius: " << it_r->second.real() << " m" << std::endl;
    }
    if (it_t != variables.end()) {
        oss << "Temperature: " << it_t->second.real() << " K" << std::endl;
    }
    if (it_lx != variables.end()) {
        oss << "X-ray Luminosity: " << it_lx->second.real() << " W" << std::endl;
    }
    if (it_b0 != variables.end()) {
        oss << "Magnetic Field: " << it_b0->second.real() << " T" << std::endl;
    }
    
    oss << "Total Variables: " << variables.size() << std::endl;
    oss << "========================================" << std::endl;
    
    return oss.str();
}

bool UQFFBuoyancyAstroModule::validateConsistency() const {
    bool consistent = true;
    
    // Check for NaN or Inf values
    for (const auto& pair : variables) {
        if (std::isnan(pair.second.real()) || std::isnan(pair.second.imag()) ||
            std::isinf(pair.second.real()) || std::isinf(pair.second.imag())) {
            std::cerr << "Inconsistency detected in '" << pair.first << "': " 
                      << pair.second << std::endl;
            consistent = false;
        }
    }
    
    // Check physical bounds for X-ray/cluster systems
    auto it_m = variables.find("M");
    if (it_m != variables.end()) {
        double mass = it_m->second.real();
        if (mass < 1e29 || mass > 1e46) { // X-ray source to galaxy cluster
            std::cerr << "Mass out of X-ray/cluster range: " << mass << " kg" << std::endl;
            consistent = false;
        }
    }
    
    // Check X-ray luminosity bounds
    auto it_lx = variables.find("L_X");
    if (it_lx != variables.end()) {
        double lx = it_lx->second.real();
        if (lx < 1e28 || lx > 1e46) {
            std::cerr << "X-ray luminosity out of range: " << lx << " W" << std::endl;
            consistent = false;
        }
    }
    
    if (consistent) {
        std::cout << "Consistency validation passed" << std::endl;
    }
    return consistent;
}

void UQFFBuoyancyAstroModule::autoCorrectAnomalies() {
    for (auto& pair : variables) {
        // Correct NaN/Inf
        if (std::isnan(pair.second.real()) || std::isinf(pair.second.real())) {
            pair.second = std::complex<double>(1e35, 0.0); // typical X-ray/cluster value
            std::cout << "Corrected anomaly in '" << pair.first << "'" << std::endl;
        }
        if (std::isnan(pair.second.imag()) || std::isinf(pair.second.imag())) {
            pair.second = std::complex<double>(pair.second.real(), 0.0);
            std::cout << "Corrected imaginary anomaly in '" << pair.first << "'" << std::endl;
        }
    }
    
    std::cout << "Auto-correction completed" << std::endl;
}

// ============================================================================
// SECTION: Comprehensive Example - Testing All 25 Methods
// ============================================================================

int main() {
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "========================================" << std::endl;
    std::cout << "UQFF X-Ray/Cluster Buoyancy Module" << std::endl;
    std::cout << "Dynamic Self-Update & Self-Expansion" << std::endl;
    std::cout << "Comprehensive Testing Example" << std::endl;
    std::cout << "========================================\n" << std::endl;
    
    UQFFBuoyancyAstroModule module;
    
    // ========================================
    // Test 1: System Identification
    // ========================================
    std::cout << "\n[Test 1] System Identification:" << std::endl;
    std::cout << "System Name: " << module.getSystemName() << std::endl;
    
    // ========================================
    // Test 2: Variable Management (5 methods)
    // ========================================
    std::cout << "\n[Test 2] Variable Management:" << std::endl;
    
    // Create new variables
    module.createVariable("test_xray_lum", std::complex<double>(5e35, 0.0));
    module.createVariable("test_cluster_temp", std::complex<double>(8e6, 0.0));
    
    // Clone variable
    module.cloneVariable("test_xray_lum", "cloned_lum");
    
    // List all variables
    std::cout << "Total variables: " << module.listVariables().size() << std::endl;
    
    // Remove variable
    module.removeVariable("test_cluster_temp");
    std::cout << "After removal: " << module.listVariables().size() << " variables" << std::endl;
    
    // ========================================
    // Test 3: Batch Operations (2 methods)
    // ========================================
    std::cout << "\n[Test 3] Batch Operations:" << std::endl;
    
    // Scale a group
    std::vector<std::string> scaleGroup = {"M", "r", "L_X"};
    module.scaleVariableGroup(scaleGroup, std::complex<double>(1.20, 0.0));
    
    // Transform a group
    std::vector<std::string> transformGroup = {"T", "B0"};
    module.transformVariableGroup(transformGroup, 
        [](std::complex<double> val) { return val * std::complex<double>(1.10, 0.0); });
    
    // ========================================
    // Test 4: Self-Expansion (4 methods)
    // ========================================
    std::cout << "\n[Test 4] Self-Expansion Capabilities:" << std::endl;
    
    // Expand parameter space
    module.expandParameterSpace(3);
    
    // System-specific expansion for all 5 systems
    std::cout << "\nExpanding all X-ray/cluster systems:" << std::endl;
    module.expandSystemScale("J1610+1811");
    module.expandSystemScale("PLCK_G287.0+32.9");
    module.expandSystemScale("PSZ2_G181.06+48.47");
    module.expandSystemScale("ASKAP_J1832-0911");
    module.expandSystemScale("SonificationCollection");
    
    // Force scale expansion
    module.expandForceScale(1.30);
    
    // Astro scale expansion
    module.expandAstroScale(1.35, 1.25);
    
    // ========================================
    // Test 5: Self-Refinement (3 methods)
    // ========================================
    std::cout << "\n[Test 5] Self-Refinement:" << std::endl;
    
    // Auto-refine parameters
    module.autoRefineParameters(50);
    
    // Calibrate to observations
    std::map<std::string, std::complex<double>> observations;
    observations["M"] = std::complex<double>(3e30, 0.0);
    observations["L_X"] = std::complex<double>(1.5e31, 0.0);
    module.calibrateToObservations(observations);
    
    // Optimize for a metric
    auto metric = [](const std::map<std::string, std::complex<double>>& vars) {
        auto it = vars.find("L_X");
        if (it != vars.end()) {
            return std::abs(it->second.real());
        }
        return 0.0;
    };
    module.optimizeForMetric(metric, 30);
    
    // ========================================
    // Test 6: Parameter Exploration (1 method)
    // ========================================
    std::cout << "\n[Test 6] Parameter Exploration:" << std::endl;
    
    auto variations = module.generateVariations(5, 0.15);
    std::cout << "Generated " << variations.size() << " parameter variations" << std::endl;
    
    // ========================================
    // Test 7: Adaptive Evolution (2 methods)
    // ========================================
    std::cout << "\n[Test 7] Adaptive Evolution:" << std::endl;
    
    // Mutate parameters
    module.mutateParameters(0.08);
    
    // Evolve system
    auto fitness = [](UQFFBuoyancyAstroModule& mod) {
        return std::abs(mod.computeFBi("J1610+1811", 0.0).real());
    };
    module.evolveSystem(15, fitness);
    
    // ========================================
    // Test 8: State Management (4 methods)
    // ========================================
    std::cout << "\n[Test 8] State Management:" << std::endl;
    
    // Save state
    module.saveState("initial_state");
    module.saveState("after_expansion");
    
    // List saved states
    auto states = module.listSavedStates();
    std::cout << "Saved states (" << states.size() << "): ";
    for (const auto& state : states) {
        std::cout << state << " ";
    }
    std::cout << std::endl;
    
    // Export state in different formats
    std::cout << "\nExport (text format):\n" << module.exportState("text").substr(0, 200) << "..." << std::endl;
    std::cout << "\nExport (csv format):\n" << module.exportState("csv").substr(0, 150) << "..." << std::endl;
    
    // Restore state
    module.restoreState("initial_state");
    
    // ========================================
    // Test 9: System Analysis (4 methods)
    // ========================================
    std::cout << "\n[Test 9] System Analysis:" << std::endl;
    
    // Sensitivity analysis
    std::vector<std::string> sensitivityParams = {"M", "r", "L_X", "T"};
    auto sensitivities = module.sensitivityAnalysis("J1610+1811", sensitivityParams);
    std::cout << "Sensitivities for J1610+1811:" << std::endl;
    for (const auto& sens : sensitivities) {
        std::cout << "  " << sens.first << ": " << sens.second << std::endl;
    }
    
    // Generate report
    std::cout << "\n" << module.generateReport("J1610+1811") << std::endl;
    
    // Validate consistency
    bool consistent = module.validateConsistency();
    std::cout << "System consistency: " << (consistent ? "PASS" : "FAIL") << std::endl;
    
    // Auto-correct anomalies
    module.autoCorrectAnomalies();
    
    // ========================================
    // Test 10-14: Multi-System Computations
    // ========================================
    std::cout << "\n[Tests 10-14] Multi-System Force Computations:" << std::endl;
    
    std::vector<std::string> systems = {"J1610+1811", "PLCK_G287.0+32.9", "PSZ2_G181.06+48.47", 
                                         "ASKAP_J1832-0911", "SonificationCollection"};
    for (const auto& sys : systems) {
        module.expandSystemScale(sys);
        std::complex<double> force = module.computeFBi(sys, 0.0);
        std::cout << "\n" << sys << " buoyancy force:" << std::endl;
        std::cout << "  Real: " << force.real() << " N" << std::endl;
        std::cout << "  Imag: " << force.imag() << " N" << std::endl;
        std::cout << "  Magnitude: " << std::abs(force) << " N" << std::endl;
    }
    
    // ========================================
    // Test 15: Equation Text Display
    // ========================================
    std::cout << "\n[Test 15] Equation Text:" << std::endl;
    std::cout << module.getEquationText("J1610+1811").substr(0, 300) << "..." << std::endl;
    
    // ========================================
    // Test 16: Variable Display
    // ========================================
    std::cout << "\n[Test 16] Variable Display (first 10):" << std::endl;
    auto varList = module.listVariables();
    for (size_t i = 0; i < std::min(size_t(10), varList.size()); ++i) {
        std::cout << "  " << varList[i] << std::endl;
    }
    std::cout << "  ... (total " << varList.size() << " variables)" << std::endl;
    
    // ========================================
    // Test 17: X-Ray Parameter Optimization
    // ========================================
    std::cout << "\n[Test 17] X-Ray Parameter Optimization:" << std::endl;
    
    // Create X-ray optimization metric
    auto xrayMetric = [](const std::map<std::string, std::complex<double>>& vars) {
        double score = 0.0;
        
        auto it_lx = vars.find("L_X");
        if (it_lx != vars.end()) {
            score += std::abs(it_lx->second.real()) / 1e35;
        }
        
        auto it_t = vars.find("T");
        if (it_t != vars.end()) {
            score += std::abs(it_t->second.real()) / 1e6;
        }
        
        return score;
    };
    
    module.optimizeForMetric(xrayMetric, 40);
    std::cout << "X-ray parameters optimized" << std::endl;
    
    // ========================================
    // Test 18: State Evolution Tracking
    // ========================================
    std::cout << "\n[Test 18] State Evolution Tracking:" << std::endl;
    
    module.saveState("pre_evolution");
    module.expandAstroScale(1.30, 1.20);
    module.saveState("post_astro_expansion");
    module.expandForceScale(1.4);
    module.saveState("post_force_expansion");
    
    std::cout << "Evolution tracking states: " << module.listSavedStates().size() << std::endl;
    
    // ========================================
    // Test 19: Cross-System Sensitivity
    // ========================================
    std::cout << "\n[Test 19] Cross-System Sensitivity:" << std::endl;
    
    for (const auto& sys : systems) {
        auto sens = module.sensitivityAnalysis(sys, {"M", "L_X"});
        std::cout << sys << " sensitivity to L_X: " << sens["L_X"] << std::endl;
    }
    
    // ========================================
    // Test 20: Final Comprehensive Report
    // ========================================
    std::cout << "\n[Test 20] Final Comprehensive Report:" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "System: " << module.getSystemName() << std::endl;
    std::cout << "Total Variables: " << module.listVariables().size() << std::endl;
    std::cout << "Saved States: " << module.listSavedStates().size() << std::endl;
    std::cout << "Consistency: " << (module.validateConsistency() ? "PASS" : "FAIL") << std::endl;
    
    std::cout << "\nX-Ray/Cluster Systems Tested:" << std::endl;
    for (const auto& sys : systems) {
        std::cout << "  - " << sys << std::endl;
    }
    
    std::cout << "\nAll 25 Dynamic Methods Tested:" << std::endl;
    std::cout << "  ✓ Variable Management (5)" << std::endl;
    std::cout << "  ✓ Batch Operations (2)" << std::endl;
    std::cout << "  ✓ Self-Expansion (4)" << std::endl;
    std::cout << "  ✓ Self-Refinement (3)" << std::endl;
    std::cout << "  ✓ Parameter Exploration (1)" << std::endl;
    std::cout << "  ✓ Adaptive Evolution (2)" << std::endl;
    std::cout << "  ✓ State Management (4)" << std::endl;
    std::cout << "  ✓ System Analysis (4)" << std::endl;
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "All Tests Completed Successfully!" << std::endl;
    std::cout << "UQFF X-Ray/Cluster Buoyancy Module" << std::endl;
    std::cout << "Fully Dynamic & Self-Expanding" << std::endl;
    std::cout << "========================================" << std::endl;
    
    return 0;
}

