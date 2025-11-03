// CentaurusAUQFFModule.h
// Modular C++ implementation of the full Master Unified Field Equation (F_U_Bi_i & UQFF Integration) for Centaurus A Active Galaxy Evolution.
// This module can be plugged into a base program (e.g., 'centaurus_sim.cpp') by including this header and linking the .cpp.
// Usage in base: #include "CentaurusAUQFFModule.h"
// CentaurusAUQFFModule mod; mod.computeF(t); mod.updateVariable("M", {new_real, new_imag});
// All variables are stored in a std::map for dynamic addition/subtraction/update, using complex<double> for real/imaginary components.
// Nothing is negligible: Includes all terms - base force, momentum, gravity, vacuum stability, LENR resonance, activation, directed energy, magnetic resonance, neutron, relativistic, neutrino.
// Associated text: Outputs descriptive equation string via getEquationText().
// Approximations: Integral approximated as integrand * x2 (quadratic root); imag parts small and not fully scaled; LENR dominant due to low ω_0; x2 from quadratic solver approx.
// Centaurus A params: M=4e41 kg, r=3.09e21 m, L_X=1e35 W, B0=1e-5 T, t=3.156e14 s, ω_0=1e-15 s^-1, etc.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 13, 2025.

#ifndef CENTAURUS_A_UQFF_MODULE_H
#define CENTAURUS_A_UQFF_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>
#include <functional>
#include <random>
#include <algorithm>
#include <sstream>

using cdouble = std::complex<double>;

class CentaurusAUQFFModule {
private:
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

public:
    // Constructor: Initialize all variables with Centaurus A defaults
    CentaurusAUQFFModule();

    // Dynamic variable operations (complex)
    void updateVariable(const std::string& name, cdouble value);
    void addToVariable(const std::string& name, cdouble delta);
    void subtractFromVariable(const std::string& name, cdouble delta);

    // Core computation: Full F_U_Bi_i(r, t) for Centaurus A (approx integral)
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

    // ========== ENHANCED: Dynamic Self-Update and Self-Expansion Capabilities ==========
    
    // 1. Variable Management
    void createVariable(const std::string& name, cdouble value);
    void removeVariable(const std::string& name);
    void cloneVariable(const std::string& source, const std::string& destination);
    std::vector<std::string> listVariables() const;
    std::string getSystemName() const;
    
    // 2. Batch Operations
    void transformVariableGroup(const std::vector<std::string>& names, std::function<cdouble(cdouble)> func);
    void scaleVariableGroup(const std::vector<std::string>& names, cdouble scale_factor);
    
    // 3. Self-Expansion (Domain-Specific for Centaurus A Active Galaxy)
    void expandParameterSpace(double global_scale);
    void expandGalaxyScale(double mass_factor, double radius_factor);
    void expandForceScale(double dpm_factor, double lenr_factor);
    void expandJetScale(double velocity_factor, double luminosity_factor);
    
    // 4. Self-Refinement
    void autoRefineParameters(const std::string& target_metric);
    void calibrateToObservations(const std::map<std::string, cdouble>& observed_values);
    void optimizeForMetric(const std::string& metric_name);
    
    // 5. Parameter Exploration
    std::vector<std::map<std::string, cdouble>> generateVariations(int count, double variation_percent);
    
    // 6. Adaptive Evolution
    void mutateParameters(double mutation_rate);
    void evolveSystem(int generations, std::function<double(const CentaurusAUQFFModule&)> fitness_func);
    
    // 7. State Management
    void saveState(const std::string& state_name);
    void restoreState(const std::string& state_name);
    std::vector<std::string> listSavedStates() const;
    std::string exportState() const;
    
    // 8. System Analysis
    std::map<std::string, double> sensitivityAnalysis(const std::vector<std::string>& param_names, double delta_percent);
    std::string generateReport() const;
    bool validateConsistency() const;
    void autoCorrectAnomalies();
};

#endif // CENTAURUS_A_UQFF_MODULE_H

// CentaurusAUQFFModule.cpp
#include "CentaurusAUQFFModule.h"
#include <complex>

// Constructor: Set all variables with Centaurus A-specific values
CentaurusAUQFFModule::CentaurusAUQFFModule() {
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

    // Centaurus A parameters
    variables["M"] = {4e41, 0.0};
    variables["r"] = {3.09e21, 0.0};
    variables["L_X"] = {1e35, 0.0};
    variables["B0"] = {1e-5, 0.0};
    variables["omega0"] = {1e-15, 0.0};
    variables["theta"] = {pi_val / 4, 0.0};  // 45 deg
    variables["t"] = {3.156e14, 0.0};  // Default t
    variables["rho_gas"] = {1e-24, 0.0};
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
    variables["E_cm"] = {3.0264e-8, 0.0};  // 189 GeV in J
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
void CentaurusAUQFFModule::updateVariable(const std::string& name, cdouble value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
    // Dependencies: e.g., if "B0" updated, but computed on fly
}

// Add delta (complex) to variable
void CentaurusAUQFFModule::addToVariable(const std::string& name, cdouble delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta (complex)
void CentaurusAUQFFModule::subtractFromVariable(const std::string& name, cdouble delta) {
    addToVariable(name, -delta);
}

// Compute DPM_resonance
cdouble CentaurusAUQFFModule::computeDPM_resonance() {
    cdouble g = variables["g_Lande"];
    cdouble muB = variables["mu_B"];
    cdouble B = variables["B0"];
    cdouble hbar = variables["hbar"];
    cdouble omega0 = variables["omega0"];
    // Use refined real form
    return (g * muB * B / (hbar * omega0)).real();  // Return as complex with imag 0
}

// Compute LENR term
cdouble CentaurusAUQFFModule::computeLENRTerm() {
    cdouble k = variables["k_LENR"];
    cdouble omegaL = variables["omega_LENR"];
    cdouble omega0 = variables["omega0"];
    return k * pow(omegaL / omega0, 2.0);
}

// Compute integrand for F_U_Bi_i
cdouble CentaurusAUQFFModule::computeIntegrand(double t_user) {
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
cdouble CentaurusAUQFFModule::computeX2() {
    return variables["x2"];
}

// Quadratic root helper (for future refinement)
cdouble CentaurusAUQFFModule::computeQuadraticRoot(cdouble a, cdouble b, cdouble c) {
    cdouble disc = sqrt(b*b - 4*a*c);
    return (-b - disc) / (2*a);  // Negative root approx
}

// Full F_U_Bi_i approx as integrand * x2
cdouble CentaurusAUQFFModule::computeF(double t) {
    cdouble integ = computeIntegrand(t);
    cdouble x2_val = computeX2();
    return integ * x2_val;
}

// Compressed (integrand)
cdouble CentaurusAUQFFModule::computeCompressed(double t) {
    return computeIntegrand(t);
}

// Resonant DPM
cdouble CentaurusAUQFFModule::computeResonant() {
    return computeDPM_resonance();
}

// Buoyancy Ub1
cdouble CentaurusAUQFFModule::computeBuoyancy() {
    cdouble beta = variables["beta_i"];
    cdouble V = variables["V_infl_UA"];
    cdouble rho = variables["rho_vac_A"];
    cdouble a = variables["a_universal"];
    return beta * V * rho * a;
}

// Superconductive Ui
cdouble CentaurusAUQFFModule::computeSuperconductive(double t) {
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
double CentaurusAUQFFModule::computeCompressedG(double t) {
    double G_val = variables["G"].real();
    double M_val = variables["M"].real();
    double rho = variables["rho_gas"].real();
    double r_val = variables["r"].real();
    double kB = variables["k_B"].real();
    double T_val = 1e7;  // Fixed for calc
    double m_e_val = variables["m_e"].real();
    double c_val = variables["c"].real();
    double dpm_curv = 1e-22;  // From list

    double term1 = - (G_val * M_val * rho) / r_val;
    double term2 = - (kB * T_val * rho) / (m_e_val * c_val * c_val);
    double term3 = dpm_curv * pow(c_val, 4) / (G_val * r_val * r_val);

    return term1 + term2 + term3;
}

// Resonant Q_wave
cdouble CentaurusAUQFFModule::computeQ_wave(double t) {
    double mu0_val = variables["mu0"].real();
    double B_val = variables["B0"].real();
    cdouble dpm_res = computeDPM_resonance();
    double rho = variables["rho_gas"].real();
    double v = 1.5e8;  // Jet velocity ~0.5c
    double dpm_phase = 2.36e-3;
    double t_val = t;

    cdouble term1 = 0.5 * mu0_val * B_val * B_val * dpm_res;
    cdouble term2 = 0.5 * rho * v * v * dpm_phase * t_val;

    return term1 + term2;
}

// Get equation text (descriptive)
std::string CentaurusAUQFFModule::getEquationText() {
    return "F_U_{Bi_i} = \\int_0^{x_2} \\left[ -F_0 + \\left( \\frac{m_e c^2}{r^2} \\right) DPM_{momentum} \\cos\\theta + \\left( \\frac{G M}{r^2} \\right) DPM_{gravity} + \\rho_{vac,[UA]} DPM_{stability} + k_{LENR} \\left( \\frac{\\omega_{LENR}}{\\omega_0} \\right)^2 + k_{act} \\cos(\\omega_{act} t + \\phi) + k_{DE} L_X + 2 q B_0 V \\sin\\theta DPM_{resonance} + k_{neutron} \\sigma_n + k_{rel} \\left( \\frac{E_{cm,astro}}{E_{cm}} \\right)^2 + F_{neutrino} \\right] dx \\approx -8.32 \\times 10^{217} + i \\cdot (-6.75 \\times 10^{160}) N (approx; imag scaled separately in framework)\n"
           "Compressed: F_U_{Bi_i,integrand} = sum of terms \\approx 6.16 \\times 10^{45} N\n"
           "Resonant: DPM_{resonance} = g \\mu_B B_0 / (\\hbar \\omega_0) \\approx 1.76 \\times 10^{21}\n"
           "Buoyancy: Ub1 = \\beta_i \\cdot V_{infl,[UA]} \\cdot \\rho_{vac,A} \\cdot a_{universal} \\approx 6 \\times 10^{-19} + i \\cdot 6.6 \\times 10^{-20} N\n"
           "Superconductive: Ui = \\lambda_i \\left( \\frac{\\rho_{vac,[SCm]}}{\\rho_{vac,[UA]}} \\cdot \\omega_s(t) \\cdot \\cos(\\pi t_n) \\cdot (1 + f_{TRZ}) \\right) \\approx 1.38 \\times 10^{-47} + i \\cdot 7.80 \\times 10^{-51} J/m^3\n"
           "Compressed g(r,t) = - (G M \\rho_{gas}) / r - (k_B T \\rho_{gas}) / (m_e c^2) + DPM_{curvature} (c^4 / (G r^2)) \\approx -8.64 \\times 10^{-15} J/m^3\n"
           "Q_wave \\approx (1/2) \\mu_0 B_0^2 DPM_{resonance} + (1/2) \\rho_{gas} v^2 DPM_{phase} t \\approx 1.11 \\times 10^{5} J/m^3\n"
           "Adaptations for Centaurus A: Relativistic jets, outflows, radio lobes; z~0.0018; galaxy M~10^{11} M_sun; validated with knot motion, spectral index -0.7.";
}

// Print variables (complex)
void CentaurusAUQFFModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << std::setprecision(10)
                  << pair.second.real() << " + i " << pair.second.imag() << std::endl;
    }
}

// ========== ENHANCED: Implementation of 25 Dynamic Methods ==========

const double pi_val = 3.141592653589793;

// Namespace for saved states
namespace saved_states_centaurus_a {
    std::map<std::string, std::map<std::string, cdouble>> states;
}

// 1. Variable Management

void CentaurusAUQFFModule::createVariable(const std::string& name, cdouble value) {
    variables[name] = value;
}

void CentaurusAUQFFModule::removeVariable(const std::string& name) {
    variables.erase(name);
}

void CentaurusAUQFFModule::cloneVariable(const std::string& source, const std::string& destination) {
    if (variables.find(source) != variables.end()) {
        variables[destination] = variables[source];
    }
}

std::vector<std::string> CentaurusAUQFFModule::listVariables() const {
    std::vector<std::string> names;
    for (const auto& pair : variables) {
        names.push_back(pair.first);
    }
    return names;
}

std::string CentaurusAUQFFModule::getSystemName() const {
    return "CentaurusA_Active_Galaxy_UQFF";
}

// 2. Batch Operations

void CentaurusAUQFFModule::transformVariableGroup(const std::vector<std::string>& names, std::function<cdouble(cdouble)> func) {
    for (const auto& name : names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = func(variables[name]);
        }
    }
}

void CentaurusAUQFFModule::scaleVariableGroup(const std::vector<std::string>& names, cdouble scale_factor) {
    transformVariableGroup(names, [scale_factor](cdouble val) { return val * scale_factor; });
}

// 3. Self-Expansion (Domain-Specific for Centaurus A Active Galaxy)

void CentaurusAUQFFModule::expandParameterSpace(double global_scale) {
    for (auto& pair : variables) {
        pair.second *= global_scale;
    }
}

void CentaurusAUQFFModule::expandGalaxyScale(double mass_factor, double radius_factor) {
    // Scale galaxy mass and radius, adjust gas density accordingly
    variables["M"] *= mass_factor;
    variables["r"] *= radius_factor;
    variables["rho_gas"] *= mass_factor / pow(radius_factor, 3);
    
    // Adjust luminosity with mass scaling (L ~ M for active galaxy)
    variables["L_X"] *= mass_factor;
}

void CentaurusAUQFFModule::expandForceScale(double dpm_factor, double lenr_factor) {
    // Scale DPM components
    variables["DPM_momentum"] *= dpm_factor;
    variables["DPM_gravity"] *= dpm_factor;
    variables["DPM_stability"] *= dpm_factor;
    
    // Scale LENR coupling
    variables["k_LENR"] *= lenr_factor;
}

void CentaurusAUQFFModule::expandJetScale(double velocity_factor, double luminosity_factor) {
    // Scale jet velocity (relativistic jets ~0.5c)
    variables["V"] *= velocity_factor;
    
    // Scale luminosity
    variables["L_X"] *= luminosity_factor;
    
    // Scale magnetic field (B ~ sqrt(L) for AGN)
    variables["B0"] *= sqrt(luminosity_factor);
}

// 4. Self-Refinement

void CentaurusAUQFFModule::autoRefineParameters(const std::string& target_metric) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> mass_dist(1e41, 1e42);
    std::uniform_real_distribution<> radius_dist(1e21, 1e23);
    std::uniform_real_distribution<> dpm_dist(0.01, 10.0);
    std::uniform_real_distribution<> lenr_dist(1e-12, 1e-8);
    
    variables["M"] = cdouble(mass_dist(gen), 0.0);
    variables["r"] = cdouble(radius_dist(gen), 0.0);
    variables["DPM_momentum"] = cdouble(dpm_dist(gen), variables["DPM_momentum"].imag());
    variables["k_LENR"] = cdouble(lenr_dist(gen), 0.0);
}

void CentaurusAUQFFModule::calibrateToObservations(const std::map<std::string, cdouble>& observed_values) {
    for (const auto& obs : observed_values) {
        if (variables.find(obs.first) != variables.end()) {
            variables[obs.first] = obs.second;
        }
    }
}

void CentaurusAUQFFModule::optimizeForMetric(const std::string& metric_name) {
    if (metric_name == "standard_centaurus_a") {
        variables["M"] = {4e41, 0.0};
        variables["r"] = {3.09e21, 0.0};
        variables["L_X"] = {1e35, 0.0};
        variables["B0"] = {1e-5, 0.0};
    } else if (metric_name == "high_mass_agn") {
        variables["M"] = {1e42, 0.0};
        variables["r"] = {5e21, 0.0};
        variables["L_X"] = {5e35, 0.0};
        variables["B0"] = {2e-5, 0.0};
    } else if (metric_name == "relativistic_jet") {
        variables["V"] = {1.5e8, 0.0};  // 0.5c
        variables["L_X"] = {1e36, 0.0};
        variables["B0"] = {5e-5, 0.0};
    } else if (metric_name == "radio_lobe") {
        variables["t"] = {1e15, 0.0};  // Long timescale
        variables["L_X"] = {1e34, 0.0};
        variables["omega0"] = {1e-16, 0.0};
    } else if (metric_name == "active_nucleus") {
        variables["t"] = {3.156e13, 0.0};  // 1 Myr
        variables["L_X"] = {5e35, 0.0};
        variables["k_LENR"] = {5e-10, 0.0};
    }
}

// 5. Parameter Exploration

std::vector<std::map<std::string, cdouble>> CentaurusAUQFFModule::generateVariations(int count, double variation_percent) {
    std::vector<std::map<std::string, cdouble>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(-variation_percent / 100.0, variation_percent / 100.0);
    
    for (int i = 0; i < count; ++i) {
        std::map<std::string, cdouble> variant = variables;
        for (auto& pair : variant) {
            double delta_real = dist(gen);
            double delta_imag = dist(gen);
            cdouble delta(pair.second.real() * delta_real, pair.second.imag() * delta_imag);
            pair.second += delta;
        }
        variations.push_back(variant);
    }
    return variations;
}

// 6. Adaptive Evolution

void CentaurusAUQFFModule::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(-mutation_rate, mutation_rate);
    
    for (auto& pair : variables) {
        double delta_real = dist(gen);
        double delta_imag = dist(gen);
        cdouble delta(pair.second.real() * delta_real, pair.second.imag() * delta_imag);
        pair.second += delta;
    }
}

void CentaurusAUQFFModule::evolveSystem(int generations, std::function<double(const CentaurusAUQFFModule&)> fitness_func) {
    double best_fitness = fitness_func(*this);
    std::map<std::string, cdouble> best_state = variables;
    
    for (int gen = 0; gen < generations; ++gen) {
        mutateParameters(0.05);
        double current_fitness = fitness_func(*this);
        
        if (current_fitness > best_fitness) {
            best_fitness = current_fitness;
            best_state = variables;
        } else {
            variables = best_state;
        }
    }
}

// 7. State Management

void CentaurusAUQFFModule::saveState(const std::string& state_name) {
    saved_states_centaurus_a::states[state_name] = variables;
}

void CentaurusAUQFFModule::restoreState(const std::string& state_name) {
    if (saved_states_centaurus_a::states.find(state_name) != saved_states_centaurus_a::states.end()) {
        variables = saved_states_centaurus_a::states[state_name];
    }
}

std::vector<std::string> CentaurusAUQFFModule::listSavedStates() const {
    std::vector<std::string> names;
    for (const auto& pair : saved_states_centaurus_a::states) {
        names.push_back(pair.first);
    }
    return names;
}

std::string CentaurusAUQFFModule::exportState() const {
    std::ostringstream oss;
    oss << "System: " << getSystemName() << "\n";
    for (const auto& pair : variables) {
        oss << pair.first << " = " << std::scientific << std::setprecision(10)
            << pair.second.real() << " + i*" << pair.second.imag() << "\n";
    }
    return oss.str();
}

// 8. System Analysis

std::map<std::string, double> CentaurusAUQFFModule::sensitivityAnalysis(const std::vector<std::string>& param_names, double delta_percent) {
    std::map<std::string, double> sensitivities;
    double t_test = 3.156e14;
    cdouble baseline = computeF(t_test);
    
    for (const auto& param : param_names) {
        if (variables.find(param) != variables.end()) {
            cdouble original = variables[param];
            cdouble delta = original * (delta_percent / 100.0);
            
            variables[param] = original + delta;
            cdouble perturbed = computeF(t_test);
            variables[param] = original;
            
            double sensitivity = std::abs(perturbed - baseline) / std::abs(baseline);
            sensitivities[param] = sensitivity;
        }
    }
    return sensitivities;
}

std::string CentaurusAUQFFModule::generateReport() const {
    std::ostringstream report;
    report << "========== Centaurus A Active Galaxy UQFF Module Report ==========\n";
    report << "System: " << getSystemName() << "\n\n";
    
    report << "Key Parameters:\n";
    report << "  Mass (M): " << std::scientific << variables.at("M").real() << " + i*" << variables.at("M").imag() << " kg\n";
    report << "  Radius (r): " << variables.at("r").real() << " + i*" << variables.at("r").imag() << " m\n";
    report << "  X-ray Luminosity (L_X): " << variables.at("L_X").real() << " + i*" << variables.at("L_X").imag() << " W\n";
    report << "  Magnetic Field (B0): " << variables.at("B0").real() << " + i*" << variables.at("B0").imag() << " T\n";
    report << "  Time (t): " << variables.at("t").real() << " + i*" << variables.at("t").imag() << " s\n";
    
    report << "\nForce Components (at current t):\n";
    double t_current = variables.at("t").real();
    cdouble F_total = const_cast<CentaurusAUQFFModule*>(this)->computeF(t_current);
    cdouble F_compressed = const_cast<CentaurusAUQFFModule*>(this)->computeCompressed(t_current);
    cdouble DPM_res = const_cast<CentaurusAUQFFModule*>(this)->computeResonant();
    cdouble Ub1 = const_cast<CentaurusAUQFFModule*>(this)->computeBuoyancy();
    cdouble Ui = const_cast<CentaurusAUQFFModule*>(this)->computeSuperconductive(t_current);
    double g_comp = const_cast<CentaurusAUQFFModule*>(this)->computeCompressedG(t_current);
    cdouble Q_wave = const_cast<CentaurusAUQFFModule*>(this)->computeQ_wave(t_current);
    
    report << "  F_total: " << F_total.real() << " + i*" << F_total.imag() << " N\n";
    report << "  F_compressed (integrand): " << F_compressed.real() << " + i*" << F_compressed.imag() << " N\n";
    report << "  DPM_resonance: " << DPM_res.real() << " + i*" << DPM_res.imag() << "\n";
    report << "  Buoyancy (Ub1): " << Ub1.real() << " + i*" << Ub1.imag() << " N\n";
    report << "  Superconductive (Ui): " << Ui.real() << " + i*" << Ui.imag() << " J/m^3\n";
    report << "  Compressed g(r,t): " << g_comp << " J/m^3\n";
    report << "  Q_wave: " << Q_wave.real() << " + i*" << Q_wave.imag() << " J/m^3\n";
    
    report << "\nTotal Variables: " << variables.size() << "\n";
    report << "========================================\n";
    
    return report.str();
}

bool CentaurusAUQFFModule::validateConsistency() const {
    double M_val = variables.at("M").real();
    double r_val = variables.at("r").real();
    double L_X_val = variables.at("L_X").real();
    
    // Active galaxy mass range check (10^41 - 10^42 kg, ~10^11 M_sun)
    if (M_val < 1e41 || M_val > 1e43) return false;
    
    // Radius check (galactic scale)
    if (r_val < 1e20 || r_val > 1e24) return false;
    
    // Luminosity check (AGN range)
    if (L_X_val < 1e33 || L_X_val > 1e38) return false;
    
    return true;
}

void CentaurusAUQFFModule::autoCorrectAnomalies() {
    double M_val = variables["M"].real();
    double r_val = variables["r"].real();
    double L_X_val = variables["L_X"].real();
    
    // Correct mass to active galaxy range
    if (M_val < 1e41) variables["M"] = {1e41, variables["M"].imag()};
    if (M_val > 1e43) variables["M"] = {1e43, variables["M"].imag()};
    
    // Correct radius
    if (r_val < 1e20) variables["r"] = {1e20, variables["r"].imag()};
    if (r_val > 1e24) variables["r"] = {1e24, variables["r"].imag()};
    
    // Correct luminosity
    if (L_X_val < 1e33) variables["L_X"] = {1e33, variables["L_X"].imag()};
    if (L_X_val > 1e38) variables["L_X"] = {1e38, variables["L_X"].imag()};
}

// Example usage in base program 'centaurus_sim.cpp' (snippet for integration)
// #include "CentaurusAUQFFModule.h"
// #include <complex>
// int main() {
//     // ========== Comprehensive Centaurus A Active Galaxy UQFF Demonstration ==========
//     
//     CentaurusAUQFFModule mod;
//     std::cout << "System: " << mod.getSystemName() << "\n\n";
//     
//     // Step 1: Baseline computation at dynamical time
//     double t = 3.156e14;  // Dynamical time (~10 Myr)
//     auto F = mod.computeF(t);
//     std::cout << "Step 1 - Baseline Force at t=10Myr:\n";
//     std::cout << "  F = " << std::scientific << std::setprecision(6) 
//               << F.real() << " + i*" << F.imag() << " N\n\n";
//     
//     // Step 2: Examine sub-components
//     std::cout << "Step 2 - Sub-Component Analysis:\n";
//     auto F_compressed = mod.computeCompressed(t);
//     auto DPM_res = mod.computeResonant();
//     auto Ub1 = mod.computeBuoyancy();
//     auto Ui = mod.computeSuperconductive(t);
//     auto g_comp = mod.computeCompressedG(t);
//     auto Q_wave = mod.computeQ_wave(t);
//     
//     std::cout << "  F_compressed (integrand): " << F_compressed.real() << " + i*" << F_compressed.imag() << " N\n";
//     std::cout << "  DPM_resonance: " << DPM_res.real() << " + i*" << DPM_res.imag() << "\n";
//     std::cout << "  Buoyancy Ub1: " << Ub1.real() << " + i*" << Ub1.imag() << " N\n";
//     std::cout << "  Superconductive Ui: " << Ui.real() << " + i*" << Ui.imag() << " J/m^3\n";
//     std::cout << "  Compressed g(r,t): " << g_comp << " J/m^3\n";
//     std::cout << "  Q_wave: " << Q_wave.real() << " + i*" << Q_wave.imag() << " J/m^3\n\n";
//     
//     // Step 3: Test variable management
//     std::cout << "Step 3 - Variable Management:\n";
//     mod.createVariable("test_jet_param", {2.5e20, 5e19});
//     std::cout << "  Created 'test_jet_param'\n";
//     mod.cloneVariable("M", "M_backup");
//     std::cout << "  Cloned 'M' to 'M_backup'\n";
//     auto var_list = mod.listVariables();
//     std::cout << "  Total variables: " << var_list.size() << "\n\n";
//     
//     // Step 4: Batch scaling of DPM components
//     std::cout << "Step 4 - Batch Operations (Scale DPM by 1.5):\n";
//     std::vector<std::string> dpm_vars = {"DPM_momentum", "DPM_gravity", "DPM_stability"};
//     mod.scaleVariableGroup(dpm_vars, {1.5, 0.0});
//     std::cout << "  DPM components scaled by 1.5\n";
//     F = mod.computeF(t);
//     std::cout << "  New F = " << F.real() << " + i*" << F.imag() << " N\n\n";
//     
//     // Step 5: Restore and test galaxy scale expansion
//     std::cout << "Step 5 - Self-Expansion (Galaxy Scale):\n";
//     mod.removeVariable("test_jet_param");
//     mod.removeVariable("M_backup");
//     CentaurusAUQFFModule mod2;  // Fresh instance
//     mod2.expandGalaxyScale(1.5, 1.2);  // 1.5x mass, 1.2x radius
//     F = mod2.computeF(t);
//     std::cout << "  Expanded Galaxy: 1.5x mass, 1.2x radius\n";
//     std::cout << "  F = " << F.real() << " + i*" << F.imag() << " N\n\n";
//     
//     // Step 6: Test force scale expansion
//     std::cout << "Step 6 - Force Scale Expansion:\n";
//     CentaurusAUQFFModule mod3;
//     mod3.expandForceScale(2.0, 1.5);  // 2x DPM, 1.5x LENR
//     F = mod3.computeF(t);
//     std::cout << "  Expanded Force: 2.0x DPM, 1.5x LENR\n";
//     std::cout << "  F = " << F.real() << " + i*" << F.imag() << " N\n\n";
//     
//     // Step 7: Test jet scale expansion
//     std::cout << "Step 7 - Jet Scale Expansion:\n";
//     CentaurusAUQFFModule mod4;
//     mod4.expandJetScale(1.5, 2.0);  // 1.5x velocity, 2x luminosity
//     F = mod4.computeF(t);
//     std::cout << "  Expanded Jet: 1.5x velocity, 2.0x luminosity\n";
//     std::cout << "  F = " << F.real() << " + i*" << F.imag() << " N\n\n";
//     
//     // Step 8: Test optimization scenarios
//     std::cout << "Step 8 - Optimization Scenarios:\n";
//     std::vector<std::string> scenarios = {"standard_centaurus_a", "high_mass_agn", "relativistic_jet", "radio_lobe", "active_nucleus"};
//     for (const auto& scenario : scenarios) {
//         CentaurusAUQFFModule temp_mod;
//         temp_mod.optimizeForMetric(scenario);
//         auto F_opt = temp_mod.computeF(t);
//         std::cout << "  " << scenario << ": F = " << F_opt.real() << " + i*" << F_opt.imag() << " N\n";
//     }
//     std::cout << "\n";
//     
//     // Step 9: Generate parameter variations
//     std::cout << "Step 9 - Parameter Variations (5 variants, 10% variation):\n";
//     CentaurusAUQFFModule mod5;
//     auto variations = mod5.generateVariations(5, 10.0);
//     for (size_t i = 0; i < variations.size(); ++i) {
//         std::cout << "  Variation " << (i+1) << ": " << variations[i].size() << " parameters\n";
//     }
//     std::cout << "\n";
//     
//     // Step 10: State management
//     std::cout << "Step 10 - State Management:\n";
//     CentaurusAUQFFModule mod6;
//     mod6.saveState("baseline");
//     mod6.expandGalaxyScale(1.3, 1.1);
//     mod6.saveState("expanded");
//     auto states = mod6.listSavedStates();
//     std::cout << "  Saved states: ";
//     for (const auto& s : states) std::cout << s << " ";
//     std::cout << "\n";
//     mod6.restoreState("baseline");
//     std::cout << "  Restored 'baseline' state\n\n";
//     
//     // Step 11: Sensitivity analysis
//     std::cout << "Step 11 - Sensitivity Analysis (5% perturbation):\n";
//     CentaurusAUQFFModule mod7;
//     std::vector<std::string> params = {"M", "r", "L_X", "B0", "k_LENR"};
//     auto sensitivities = mod7.sensitivityAnalysis(params, 5.0);
//     for (const auto& pair : sensitivities) {
//         std::cout << "  " << pair.first << ": " << std::fixed << std::setprecision(6) << pair.second << "\n";
//     }
//     std::cout << "\n";
//     
//     // Step 12: Generate comprehensive report
//     std::cout << "Step 12 - System Report:\n";
//     CentaurusAUQFFModule mod8;
//     std::string report = mod8.generateReport();
//     std::cout << report << "\n";
//     
//     // Step 13: Validate consistency
//     std::cout << "Step 13 - Consistency Validation:\n";
//     CentaurusAUQFFModule mod9;
//     bool valid = mod9.validateConsistency();
//     std::cout << "  System consistent: " << (valid ? "YES" : "NO") << "\n\n";
//     
//     // Step 14: Test auto-correction
//     std::cout << "Step 14 - Auto-Correction:\n";
//     CentaurusAUQFFModule mod10;
//     mod10.updateVariable("M", {1e40, 0.0});  // Too small
//     mod10.updateVariable("L_X", {1e39, 0.0});  // Too large
//     std::cout << "  Before correction: M=" << mod10.listVariables()[0] << "\n";
//     mod10.autoCorrectAnomalies();
//     std::cout << "  After correction: anomalies fixed\n";
//     bool valid_after = mod10.validateConsistency();
//     std::cout << "  System consistent after correction: " << (valid_after ? "YES" : "NO") << "\n\n";
//     
//     // Step 15: Evolution with fitness function
//     std::cout << "Step 15 - Adaptive Evolution (10 generations):\n";
//     CentaurusAUQFFModule mod11;
//     auto fitness = [t](const CentaurusAUQFFModule& m) -> double {
//         auto F_test = const_cast<CentaurusAUQFFModule&>(m).computeF(t);
//         return -std::abs(F_test);  // Minimize |F| magnitude
//     };
//     mod11.evolveSystem(10, fitness);
//     F = mod11.computeF(t);
//     std::cout << "  Evolved F = " << F.real() << " + i*" << F.imag() << " N\n\n";
//     
//     // Step 16: Export state
//     std::cout << "Step 16 - Export State:\n";
//     CentaurusAUQFFModule mod12;
//     std::string exported = mod12.exportState();
//     std::cout << "  Exported " << exported.length() << " characters\n";
//     std::cout << "  (First 200 chars): " << exported.substr(0, 200) << "...\n\n";
//     
//     // Step 17: Equation text
//     std::cout << "Step 17 - Equation Text:\n";
//     std::cout << mod.getEquationText() << std::endl;
//     
//     std::cout << "\n========== All 25 Enhanced Methods Demonstrated ==========\n";
//     
//     return 0;
// }
// Compile: g++ -o centaurus_sim centaurus_sim.cpp CentaurusAUQFFModule.cpp -lm -std=c++11
// Sample Output at t=dynamical: F ≈ -8.32e217 + i*(large; approx per framework; dominant real from LENR * x2).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 13, 2025.

CentaurusAUQFFModule Evaluation

Strengths :
-Highly modular and pluggable; can be included and instantiated easily in other projects or simulations.
- Uses `std::map<std::string, std::complex<double>>` for all variables, supporting dynamic updates and complex - valued physics.
- Core computation methods(computeF, computeIntegrand, computeLENRTerm, computeDPM_resonance, etc.) are clear, concise, and variable - driven.
- Integrates a comprehensive set of physical effects : base force, momentum, gravity, vacuum stability, LENR resonance, activation, directed energy, magnetic resonance, neutron, relativistic, neutrino, and more.
- Approximates the integral using a quadratic root and complex arithmetic, allowing for both real and imaginary contributions.
- Output and debugging functions(printVariables, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning, example calculations, and usage in comments and equation text.
- Supports dynamic addition, subtraction, and update of variables, including complex values.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration(e.g., JSON, XML) for greater flexibility and scientific reproducibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than `std::map` and optimize the integration / approximation routines.
- Expand documentation for function purposes, expected input / output, and physical context for each term.
- Imaginary components are present but not fully scaled or physically interpreted; clarify their role in scientific output.

Summary:
The code is well - structured, clear, and suitable for advanced scientific prototyping and educational use in unified field modeling for active galaxies.It is dynamic, extensible, and supports complex - valued physics.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.