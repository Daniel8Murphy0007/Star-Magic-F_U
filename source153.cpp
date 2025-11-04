// Abell2256UQFFModule.h
// Modular C++ implementation of the full Master Unified Field Equation (F_U_Bi_i & UQFF Integration) for Abell 2256 Galaxy Cluster Evolution.
// This module can be plugged into a base program (e.g., 'abell_sim.cpp') by including this header and linking the .cpp.
// Usage in base: #include "Abell2256UQFFModule.h"
// Abell2256UQFFModule mod; mod.computeF(t); mod.updateVariable("M", {new_real, new_imag});
// All variables are stored in a std::map for dynamic addition/subtraction/update, using complex<double> for real/imaginary components.
// Nothing is negligible: Includes all terms - base force, momentum, gravity, vacuum stability, LENR resonance, activation, directed energy, magnetic resonance, neutron, relativistic, neutrino.
// Associated text: Outputs descriptive equation string via getEquationText().
// Approximations: Integral approximated as integrand * x2 (quadratic root); imag parts small and not fully scaled; LENR dominant due to low Ï‰_0; x2 from quadratic solver approx.
// Abell 2256 params: M=1.23e45 kg, r=3.93e22 m, L_X=3.7e37 W, B0=1e-9 T, t=6.31e15 s, Ï‰_0=1e-15 s^-1, etc.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 11, 2025.

#ifndef ABELL2256_UQFF_MODULE_H
#define ABELL2256_UQFF_MODULE_H

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
using cdouble = std::complex<double>;

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

class Abell2256UQFFModule {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
    std::map<std::string, cdouble> variables;
    cdouble computeIntegrand(double t);
    cdouble computeDPM_resonance();
    cdouble computeX2();
    cdouble computeQuadraticRoot(cdouble a, cdouble b, cdouble c);
    cdouble computeLENRTerm();
    double computeG(double t);
    cdouble computeQ_wave(double t);
    cdouble computeUb1();
    cdouble computeUi(double t);
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;



public:
    // Constructor: Initialize all variables with Abell 2256 defaults
    Abell2256UQFFModule();

    // Dynamic variable operations (complex)
    void updateVariable(const std::string& name, cdouble value);
    void addToVariable(const std::string& name, cdouble delta);
    void subtractFromVariable(const std::string& name, cdouble delta);

    // Core computation: Full F_U_Bi_i(r, t) for Abell 2256 (approx integral)
    cdouble computeF(double t);

    // Sub-equations
    cdouble computeCompressed(double t);  // Integrand
    cdouble computeResonant();
    cdouble computeBuoyancy();
    cdouble computeSuperconductive(double t);
    double computeCompressedG(double t);  // g(r,t)

    // Output descriptive text of the equation
    std::string getEquationText();

    // Print all current variables (for debugging/updates)
    void printVariables();
};

#endif // ABELL2256_UQFF_MODULE_H

// Abell2256UQFFModule.cpp
#include "Abell2256UQFFModule.h"
#include <complex>

// Constructor: Set all variables with Abell 2256-specific values
Abell2256UQFFModule::Abell2256UQFFModule() {
        enableDynamicTerms = true;
        enableLogging = false;
        learningRate = 0.001;
        metadata["enhanced"] = "true";
        metadata["version"] = "2.0-Enhanced";

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

    // Abell 2256 parameters
    variables["M"] = {1.23e45, 0.0};
    variables["r"] = {3.93e22, 0.0};
    variables["L_X"] = {3.7e37, 0.0};
    variables["B0"] = {1e-9, 0.0};
    variables["omega0"] = {1e-15, 0.0};
    variables["theta"] = {pi_val / 4, 0.0};  // 45 deg
    variables["t"] = {6.31e15, 0.0};  // Default t
    variables["rho_gas"] = {5e-24, 0.0};
    variables["V"] = {1e-3, 0.0};  // Particle velocity
    variables["F0"] = {1.83e71, 0.0};

    // Vacuum and DPM
    variables["rho_vac_UA"] = {7.09e-36, 1e-37};
    variables["DPM_momentum"] = {0.93, 0.05};
    variables["DPM_gravity"] = {1.0, 0.1};
    variables["DPM_stability"] = {0.01, 0.001};

    // LENR and activation
    variables["k_LENR"] = {1e-10, 0.0};
    variables["omega_LENR"] = {2 * pi_val * 1.25e12, 0.0};
    variables["k_act"] = {1e-6, 0.0};
    variables["omega_act"] = {2 * pi_val * 300, 0.0};
    variables["phi"] = {pi_val / 4, 0.0};

    // Other couplings
    variables["k_DE"] = {1e-30, 0.0};
    variables["k_neutron"] = {1e10, 0.0};
    variables["sigma_n"] = {1e-4, 0.0};
    variables["k_rel"] = {1e-10, 0.0};
    variables["E_cm_astro"] = {1.24e24, 0.0};  // Refined, imag 0 for simplicity
    variables["E_cm"] = {2.18e-6, 0.0};  // 13.6 TeV in J
    variables["F_neutrino"] = {9.07e-42, 1e-43};

    // Quadratic approx
    variables["x2"] = {-1.35e172, 0.0};  // Refined approx root

    // Buoyancy
    variables["beta_i"] = {0.6, 0.0};
    variables["V_infl_UA"] = {1e-6, 1e-7};
    variables["rho_vac_A"] = {1e-30, 1e-31};
    variables["a_universal"] = {1e12, 1e11};

    // Superconductive
    variables["lambda_i"] = {1.0, 0.0};
    variables["rho_vac_SCm"] = {7.09e-37, 1e-38};
    variables["omega_s"] = {2.5e-6, 1e-7};
    variables["f_TRZ"] = {0.1, 0.0};
    variables["t_scale"] = {1e16, 0.0};  // For t_n = t / t_scale
}

// Update variable (set to new complex value)
void Abell2256UQFFModule::updateVariable(const std::string& name, cdouble value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
    // Dependencies: e.g., if "B0" updated, but computed on fly
}

// Add delta (complex) to variable
void Abell2256UQFFModule::addToVariable(const std::string& name, cdouble delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta (complex)
void Abell2256UQFFModule::subtractFromVariable(const std::string& name, cdouble delta) {
    addToVariable(name, -delta);
}

// Compute DPM_resonance
cdouble Abell2256UQFFModule::computeDPM_resonance() {
    cdouble g = variables["g_Lande"];
    cdouble muB = variables["mu_B"];
    cdouble B = variables["B0"];
    cdouble hbar = variables["hbar"];
    cdouble omega0 = variables["omega0"];
    // Use refined real form
    return (g * muB * B / (hbar * omega0)).real();  // Return as complex with imag 0
}

// Compute LENR term
cdouble Abell2256UQFFModule::computeLENRTerm() {
    cdouble k = variables["k_LENR"];
    cdouble omegaL = variables["omega_LENR"];
    cdouble omega0 = variables["omega0"];
    return k * pow(omegaL / omega0, 2.0);
}

// Compute integrand for F_U_Bi_i
cdouble Abell2256UQFFModule::computeIntegrand(double t_user) {
    variables["t"] = {t_user, 0.0};
    double cos_theta = cos(variables["theta"].real());
    double sin_theta = sin(variables["theta"].real());
    double cos_act = cos(variables["omega_act"].real() * t_user + variables["phi"].real());

    cdouble term_base = -variables["F0"];
    cdouble term_mom = (variables["m_e"] * pow(variables["c"], 2) / pow(variables["r"], 2)) * variables["DPM_momentum"] * cos_theta;
    cdouble term_grav = (variables["G"] * variables["M"] / pow(variables["r"], 2)) * variables["DPM_gravity"];
    cdouble term_vac = variables["rho_vac_UA"] * variables["DPM_stability"];
    cdouble term_LENR = computeLENRTerm();
    cdouble term_act = variables["k_act"] * cos_act;
    cdouble term_DE = variables["k_DE"] * variables["L_X"];
    cdouble term_res = 2 * variables["q"] * variables["B0"] * variables["V"] * sin_theta * computeDPM_resonance();
    cdouble term_neut = variables["k_neutron"] * variables["sigma_n"];
    cdouble term_rel = variables["k_rel"] * pow(variables["E_cm_astro"] / variables["E_cm"], 2);
    cdouble term_neutrino = variables["F_neutrino"];

    return term_base + term_mom + term_grav + term_vac + term_LENR + term_act + term_DE + term_res + term_neut + term_rel + term_neutrino;
}

// Approx x2 (hardcoded refined for stability; dynamic via var)
cdouble Abell2256UQFFModule::computeX2() {
    return variables["x2"];
}

// Quadratic root helper (for future refinement)
cdouble Abell2256UQFFModule::computeQuadraticRoot(cdouble a, cdouble b, cdouble c) {
    cdouble disc = sqrt(b*b - 4*a*c);
    return (-b - disc) / (2*a);  // Negative root approx
}

// Full F_U_Bi_i approx as integrand * x2
cdouble Abell2256UQFFModule::computeF(double t) {
    cdouble integ = computeIntegrand(t);
    cdouble x2_val = computeX2();
    return integ * x2_val;
}

// Compressed (integrand)
cdouble Abell2256UQFFModule::computeCompressed(double t) {
    return computeIntegrand(t);
}

// Resonant DPM
cdouble Abell2256UQFFModule::computeResonant() {
    return computeDPM_resonance();
}

// Buoyancy Ub1
cdouble Abell2256UQFFModule::computeBuoyancy() {
    cdouble beta = variables["beta_i"];
    cdouble V = variables["V_infl_UA"];
    cdouble rho = variables["rho_vac_A"];
    cdouble a = variables["a_universal"];
    return beta * V * rho * a;
}

// Superconductive Ui
cdouble Abell2256UQFFModule::computeSuperconductive(double t) {
    double tn = t / variables["t_scale"].real();
    cdouble lambda = variables["lambda_i"];
    cdouble rho_sc = variables["rho_vac_SCm"];
    cdouble rho_ua = variables["rho_vac_UA"];
    cdouble omega_s = variables["omega_s"];
    double cos_term = cos(pi_val * tn);
    cdouble f_trz = variables["f_TRZ"];
    return lambda * (rho_sc / rho_ua * omega_s * cos_term * (1 + f_trz.real()));
}

// Compressed g(r,t)
double Abell2256UQFFModule::computeCompressedG(double t) {
    double G_val = variables["G"].real();
    double M_val = variables["M"].real();
    double rho = variables["rho_gas"].real();
    double r_val = variables["r"].real();
    double kB = variables["k_B"].real();
    double T_val = 8e7;  // Fixed for calc
    double m_e_val = variables["m_e"].real();
    double c_val = variables["c"].real();
    double dpm_curv = 1e-22;  // From list

    double term1 = - (G_val * M_val * rho) / r_val;
    double term2 = - (kB * T_val * rho) / (m_e_val * c_val * c_val);
    double term3 = dpm_curv * pow(c_val, 4) / (G_val * r_val * r_val);

    return term1 + term2 + term3;
}

// Resonant Q_wave
cdouble Abell2256UQFFModule::computeQ_wave(double t) {
    double mu0_val = variables["mu0"].real();
    double B_val = variables["B0"].real();
    cdouble dpm_res = computeDPM_resonance();
    double rho = variables["rho_gas"].real();
    double v = 1.7e6;  // Velocity dispersion
    double dpm_phase = 2.36e-3;
    double t_val = t;

    cdouble term1 = 0.5 * mu0_val * B_val * B_val * dpm_res;
    cdouble term2 = 0.5 * rho * v * v * dpm_phase * t_val;

    return term1 + term2;
}

// Get equation text (descriptive)
std::string Abell2256UQFFModule::getEquationText() {
    return "F_U_{Bi_i} = \\int_0^{x_2} \\left[ -F_0 + \\left( \\frac{m_e c^2}{r^2} \\right) DPM_{momentum} \\cos\\theta + \\left( \\frac{G M}{r^2} \\right) DPM_{gravity} + \\rho_{vac,[UA]} DPM_{stability} + k_{LENR} \\left( \\frac{\\omega_{LENR}}{\\omega_0} \\right)^2 + k_{act} \\cos(\\omega_{act} t + \\phi) + k_{DE} L_X + 2 q B_0 V \\sin\\theta DPM_{resonance} + k_{neutron} \\sigma_n + k_{rel} \\left( \\frac{E_{cm,astro}}{E_{cm}} \\right)^2 + F_{neutrino} \\right] dx \\approx -8.32 \\times 10^{217} + i \\cdot (-6.75 \\times 10^{160}) N (approx; imag scaled separately in framework)\n"
           "Compressed: F_U_{Bi_i,integrand} = sum of terms \\approx 6.16 \\times 10^{45} N\n"
           "Resonant: DPM_{resonance} = g \\mu_B B_0 / (\\hbar \\omega_0) \\approx 1.76 \\times 10^{17}\n"
           "Buoyancy: Ub1 = \\beta_i \\cdot V_{infl,[UA]} \\cdot \\rho_{vac,A} \\cdot a_{universal} \\approx 6 \\times 10^{-19} + i \\cdot 6.6 \\times 10^{-20} N\n"
           "Superconductive: Ui = \\lambda_i \\left( \\frac{\\rho_{vac,[SCm]}}{\\rho_{vac,[UA]}} \\cdot \\omega_s(t) \\cdot \\cos(\\pi t_n) \\cdot (1 + f_{TRZ}) \\right) \\approx 1.38 \\times 10^{-47} + i \\cdot 7.80 \\times 10^{-51} J/m^3\n"
           "Compressed g(r,t) = - (G M \\rho_{gas}) / r - (k_B T \\rho_{gas}) / (m_e c^2) + DPM_{curvature} (c^4 / (G r^2)) \\approx -1.05 \\times 10^{-11} m/s^2\n"
           "Q_wave \\approx (1/2) \\mu_0 B_0^2 DPM_{resonance} + (1/2) \\rho_{gas} v^2 DPM_{phase} t \\approx 1.07 \\times 10^{-4} J/m^3\n"
           "Adaptations for Abell 2256: Merger shocks, radio halo/relics, ICM gas; z=0.058; M500=1.23e45 kg; validated with spectral index -1.56, velocity ~1700 km/s.";
}

// Print variables (complex)
void Abell2256UQFFModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << std::setprecision(10)
                  << pair.second.real() << " + i " << pair.second.imag() << std::endl;
    }
}

// Example usage in base program 'abell_sim.cpp' (snippet for integration)
// #include "Abell2256UQFFModule.h"
// #include <complex>
// int main() {
//     Abell2256UQFFModule mod;
//     double t = 6.31e15;  // 0.2 Gyr
//     auto F = mod.computeF(t);
//     std::cout << "F = " << F.real() << " + i " << F.imag() << " N\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("M", {1.5e45, 0.0});  // Update mass
//     mod.addToVariable("f_TRZ", {0.05, 0.0});  // Add to TR factor
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o abell_sim abell_sim.cpp Abell2256UQFFModule.cpp -lm
// Sample Output at t=0.2 Gyr: F â‰ˆ -8.32e217 + i (large; approx per framework; dominant real from LENR * x2).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 11, 2025.
// ===== MISSING ABELL 2256 BUOYANCY FUNCTION IMPLEMENTATIONS =====

// Compute Ub1 buoyancy term for galaxy cluster merger dynamics
cdouble Abell2256UQFFModule::computeUb1() {
    // Enhanced buoyancy calculation for Abell 2256 galaxy cluster
    cdouble beta = variables["beta_i"];
    cdouble V = variables["V_infl_UA"];
    cdouble rho_vac = variables["rho_vac_UA"];
    cdouble rho_A = variables["rho_vac_A"];
    cdouble a = variables["a_universal"];
    cdouble G = variables["G"];
    cdouble M = variables["M"];
    cdouble r = variables["r"];
    cdouble c = variables["c"];
    
    // Base buoyancy force
    cdouble base_buoyancy = beta * V * rho_A * a;
    
    // Galaxy cluster gravitational enhancement
    cdouble cluster_grav = G * M / (r * r);
    
    // Intracluster medium (ICM) buoyancy effects
    cdouble rho_ICM = variables["rho_gas"];  // ICM gas density
    cdouble T_ICM = 8e7;  // ICM temperature ~80 keV
    cdouble k_B = variables["k_B"];
    cdouble pressure_term = rho_ICM * k_B * T_ICM / (variables["m_e"] * c * c);
    
    // Radio halo and relic contributions (Abell 2256 specific)
    cdouble L_X = variables["L_X"];
    cdouble B0 = variables["B0"];
    cdouble radio_enhancement = std::sqrt(L_X / 1e37) * std::sqrt(B0 / 1e-9);
    
    // Merger shock dynamics (major-minor merger at z=0.058)
    double velocity_dispersion = 1700e3;  // m/s, characteristic velocity
    cdouble merger_factor = velocity_dispersion * velocity_dispersion / (c * c);
    
    // Dark matter halo effects
    double M500 = M.real();  // M500 mass
    cdouble dm_factor = std::log10(M500 / 1e45);  // Mass scaling
    
    // Cosmic ray pressure and relativistic effects
    cdouble relativistic_pressure = variables["k_rel"] * pressure_term * merger_factor;
    
    return base_buoyancy + cluster_grav * pressure_term * radio_enhancement * (1.0 + merger_factor) * dm_factor + relativistic_pressure;
}

// Compute Ui superconductive term for galaxy cluster magnetism and dynamics
cdouble Abell2256UQFFModule::computeUi(double t) {
    double pi_val = variables["pi"].real();
    double tn = t / variables["t_scale"].real();
    cdouble lambda = variables["lambda_i"];
    cdouble rho_sc = variables["rho_vac_SCm"];
    cdouble rho_ua = variables["rho_vac_UA"];
    cdouble omega_s = variables["omega_s"];
    cdouble f_trz = variables["f_TRZ"];
    
    // Time-dependent oscillations
    double cos_term = cos(pi_val * tn);
    
    // Abell 2256 specific magnetic field enhancements
    cdouble B0 = variables["B0"];
    cdouble mu0 = variables["mu0"];
    cdouble magnetic_energy = B0 * B0 / (2.0 * mu0);
    
    // ICM turbulence and magnetic field amplification
    double turbulent_velocity = 500e3;  // m/s, turbulent velocity in ICM
    cdouble turbulence_factor = turbulent_velocity / variables["c"];
    
    // Merger-induced magnetic field enhancement
    double merger_age = 6.31e15;  // s, ~0.2 Gyr since merger
    cdouble time_evolution = std::exp(-t / (2.0 * merger_age));  // Decay factor
    
    // Radio halo magnetic field coupling
    cdouble L_X = variables["L_X"];
    cdouble synchrotron_scaling = std::sqrt(L_X * magnetic_energy);
    
    // Relativistic particle acceleration in merger shocks
    cdouble shock_acceleration = variables["k_rel"] * turbulence_factor * turbulence_factor;
    
    // Dark matter interaction effects
    cdouble M = variables["M"];
    cdouble r = variables["r"];
    cdouble dm_coupling = M / (4.0 * pi_val * r * r * r);  // Mass density profile
    
    // Enhanced superconductive term with galaxy cluster physics
    cdouble base_ui = lambda * (rho_sc / rho_ua) * omega_s * cos_term * (1.0 + f_trz);
    cdouble cluster_enhancement = magnetic_energy * turbulence_factor * time_evolution * shock_acceleration;
    cdouble radio_coupling = synchrotron_scaling * dm_coupling;
    
    return base_ui * cluster_enhancement + radio_coupling;
}

// ===== ENHANCED DYNAMIC CAPABILITIES =====

// Auto-calibrate parameters to match observational targets
void Abell2256UQFFModule::autoCalibrate(const std::string& observable, double target_value, double tolerance) {
    if (variables.find(observable) == variables.end()) {
        std::cerr << "Observable '" << observable << "' not found for calibration." << std::endl;
        return;
    }
    
    double current_value = variables[observable].real();
    double error = std::abs(current_value - target_value) / target_value;
    
    if (error > tolerance) {
        // Gradient-based parameter adjustment
        std::vector<std::string> tunable_params = {"M", "r", "B0", "rho_gas", "L_X"};
        
        for (const auto& param : tunable_params) {
            cdouble gradient = computeGradient(param, observable);
            if (std::abs(gradient) > 1e-20) {
                cdouble adjustment = learning_rate * (target_value - current_value) / gradient;
                variables[param] += adjustment;
                recordHistory(param, variables[param]);
            }
        }
        
        std::cout << "Auto-calibrated " << observable << " from " << current_value 
                  << " to target " << target_value << " (error: " << error << ")" << std::endl;
    }
}

// Adaptive parameter updates based on temporal evolution
void Abell2256UQFFModule::adaptiveUpdate(double dt, const std::string& feedback_param) {
    if (!self_learning_enabled) return;
    
    // Merger evolution timescale
    double merger_timescale = 6.31e15;  // ~0.2 Gyr
    double evolution_factor = std::exp(-dt / merger_timescale);
    
    // Adaptive magnetic field evolution
    cdouble B0_old = variables["B0"];
    variables["B0"] *= evolution_factor;
    
    // Adaptive velocity dispersion decay
    double velocity_decay = 0.99;  // 1% decay per update
    if (variables.find("velocity_dispersion") != variables.end()) {
        variables["velocity_dispersion"] *= velocity_decay;
    } else {
        variables["velocity_dispersion"] = {1700e3 * velocity_decay, 0.0};
    }
    
    // ICM cooling and heating balance
    cdouble T_ICM_factor = 1.0 + 0.01 * std::sin(2 * M_PI * dt / merger_timescale);
    if (variables.find("T_ICM") != variables.end()) {
        variables["T_ICM"] *= T_ICM_factor;
    } else {
        variables["T_ICM"] = {8e7 * T_ICM_factor.real(), 0.0};
    }
    
    recordHistory("adaptive_time", {dt, 0.0});
    std::cout << "Adaptive update: B0=" << variables["B0"].real() 
 * Enhanced: November 04, 2025 - Added self-expanding capabilities
              << ", v_disp=" << variables["velocity_dispersion"].real() << std::endl;
}

// Scale parameters to match observational data
void Abell2256UQFFModule::scaleToObservations(const std::map<std::string, double>& observations) {
    for (const auto& obs : observations) {
        if (variables.find(obs.first) != variables.end()) {
            double scaling = obs.second / variables[obs.first].real();
            variables[obs.first] *= scaling;
            
            // Scale related parameters
            if (obs.first == "L_X") {
                variables["B0"] *= std::sqrt(scaling);  // B scales with sqrt(L_X)
                variables["rho_gas"] *= scaling;        // Gas density scales with L_X
            }
            if (obs.first == "M") {
                variables["r"] *= std::pow(scaling, 1.0/3.0);  // r scales with M^(1/3)
            }
        }
    }
    std::cout << "Scaled to " << observations.size() << " observational constraints." << std::endl;
}

// Add custom variables with dependency tracking
void Abell2256UQFFModule::addCustomVariable(const std::string& name, cdouble value, const std::string& dependency) {
    variables[name] = value;
    if (!dependency.empty()) {
        variable_dependencies[name] = dependency;
    }
    recordHistory(name, value);
    std::cout << "Added custom variable: " << name << " = " << value << std::endl;
}

// Get variable evolution history
std::map<std::string, cdouble> Abell2256UQFFModule::getVariableHistory(const std::string& name, int steps) {
    std::map<std::string, cdouble> history;
    if (variable_history.find(name) != variable_history.end()) {
        auto& hist = variable_history[name];
        int start = std::max(0, (int)hist.size() - steps);
        for (int i = start; i < (int)hist.size(); i++) {
            history["step_" + std::to_string(i)] = hist[i];
        }
    }
    return history;
}

// Enable/disable self-learning capabilities
void Abell2256UQFFModule::enableSelfLearning(bool enable) {
    self_learning_enabled = enable;
    if (enable) {
        std::cout << "Self-learning enabled with rate: " << learning_rate << std::endl;
    } else {
        std::cout << "Self-learning disabled." << std::endl;
    }
}

// Export current state for persistence
void Abell2256UQFFModule::exportState(const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        file << "# Abell2256UQFFModule State Export" << std::endl;
        file << "update_counter=" << update_counter << std::endl;
        file << "learning_rate=" << learning_rate << std::endl;
        file << "self_learning_enabled=" << (self_learning_enabled ? 1 : 0) << std::endl;
        
        for (const auto& var : variables) {
            file << var.first << "=" << var.second.real() << "," << var.second.imag() << std::endl;
        }
        file.close();
        std::cout << "State exported to: " << filename << std::endl;
    }
}

// Import state from file
void Abell2256UQFFModule::importState(const std::string& filename) {
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
                    size_t comma_pos = value_str.find(',');
                    if (comma_pos != std::string::npos) {
                        double real_part = std::stod(value_str.substr(0, comma_pos));
                        double imag_part = std::stod(value_str.substr(comma_pos + 1));
                        variables[key] = {real_part, imag_part};
                    }
                }
            }
        }
        file.close();
        std::cout << "State imported from: " << filename << std::endl;
    }
}

// Helper function: Update dependent variables
void Abell2256UQFFModule::updateDependencies(const std::string& changed_var) {
    // Automatic dependency updates
    if (changed_var == "M") {
        // Update Schwarzschild radius
        cdouble M = variables["M"];
        cdouble G = variables["G"];
        cdouble c = variables["c"];
        variables["r_schwarzschild"] = 2.0 * G * M / (c * c);
    }
    
    if (changed_var == "B0") {
        // Update magnetic energy density
        cdouble B0 = variables["B0"];
        cdouble mu0 = variables["mu0"];
        variables["u_magnetic"] = B0 * B0 / (2.0 * mu0);
    }
    
    if (changed_var == "L_X") {
        // Update implied temperature
        cdouble L_X = variables["L_X"];
        variables["T_implied"] = std::pow(L_X / 1e37, 0.25) * 8e7;  // Rough scaling
    }
}

// Helper function: Compute numerical gradient
cdouble Abell2256UQFFModule::computeGradient(const std::string& var, const std::string& target) {
    if (variables.find(var) == variables.end() || variables.find(target) == variables.end()) {
        return {0.0, 0.0};
    }
    
    cdouble original_value = variables[var];
    cdouble original_target = variables[target];
    
    // Small perturbation
    cdouble delta = original_value * 1e-6;
    variables[var] += delta;
    
    // Recompute target (simplified - would need full recalculation in practice)
    cdouble new_target = computeF(variables["t"].real());  // Use main computation
    
    // Restore original value
    variables[var] = original_value;
    
    return (new_target - original_target) / delta;
}

// Helper function: Record variable history
void Abell2256UQFFModule::recordHistory(const std::string& name, cdouble value) {
    variable_history[name].push_back(value);
    
    // Keep only last 100 values to prevent memory bloat
    if (variable_history[name].size() > 100) {
        variable_history[name].erase(variable_history[name].begin());
    }
}
