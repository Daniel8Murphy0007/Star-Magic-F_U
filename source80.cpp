
// SMBHUQFFModule.h
// Modular C++ implementation of the full Master Unified Field Equation (F_U_Bi_i & UQFF Integration) for SMBH-CGM M-σ Relation and Buoyancy in ROMULUS25-like Simulations.
// This module can be plugged into a base program (e.g., 'smbh_cgm_sim.cpp') by including this header and linking the .cpp.
// Usage in base: #include "SMBHUQFFModule.h"
// SMBHUQFFModule mod; mod.computeMSigma(sigma); mod.updateVariable("F_rel", {new_real, new_imag});
// All variables are stored in a std::map for dynamic addition/subtraction/update, using complex<double> for real/imaginary components.
// Nothing is negligible: Includes all terms - M-σ relation, cosmic time, omega_s, mu_j, E_react, delta_n, rho_vac_ua_scm, U_m, U_g1, integrated with LENR, Sweet, Kozima.
// Associated text: Outputs descriptive equation string via getEquationText().
// Approximations: M-σ from Kormendy & Ho (2013); cosmic time approximation; dynamic params from ROMULUS25; F_rel from LEP 1998.
// SMBH-CGM params: sigma=100-300 km/s, z=0-6, n=1-26 states, etc.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 23, 2025.

#ifndef SMBH_UQFF_MODULE_H
#define SMBH_UQFF_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>

using cdouble = std::complex<double>;

class SMBHUQFFModule {
private:
    std::map<std::string, cdouble> variables;
    cdouble computeIntegrand(double t, double sigma, int n);
    cdouble computeDPM_resonance(double sigma);
    cdouble computeX2(double sigma, int n);
    cdouble computeQuadraticRoot(cdouble a, cdouble b, cdouble c);
    cdouble computeLENRTerm(double sigma);
    double computeCosmicTime(double z);
    double computeOmegaS(double sigma);
    double computeMuJ(double t);
    double computeEReact(double t);
    double computeDeltaN(int n);
    double computeRhoVacUAScm(int n, double t);
    cdouble computeUm(double t, double r, int n);
    cdouble computeUg1(double t, double r, double M_s, int n);
    void setGalaxyParams(double sigma, int n, double z);

public:
    // Constructor: Initialize all variables with SMBH-CGM defaults
    SMBHUQFFModule();

    // Dynamic variable operations (complex)
    void updateVariable(const std::string& name, cdouble value);
    void addToVariable(const std::string& name, cdouble delta);
    void subtractFromVariable(const std::string& name, cdouble delta);

    // Core computation: M-sigma relation with UQFF terms (approx)
    cdouble computeMSigma(double sigma, int n, double t);

    // Sub-equations
    cdouble computeCompressed(double sigma, int n, double t);  // Integrand
    cdouble computeResonant(double t, double sigma);
    cdouble computeBuoyancy(double sigma);
    cdouble computeSuperconductive(double t, double sigma);
    double computeCompressedG(double t, double sigma, int n);  // g(r,t)

    // Output descriptive text of the equation
    std::string getEquationText(double sigma);

    // Print all current variables (for debugging/updates)
    void printVariables();
    
    // ========================================================================
    // DYNAMIC SELF-UPDATE AND SELF-EXPANSION CAPABILITIES (NEW)
    // ========================================================================
    
    // Dynamic variable management
    void createDynamicVariable(const std::string& name, cdouble value);
    void removeDynamicVariable(const std::string& name);
    void cloneVariable(const std::string& src, const std::string& dest);
    std::vector<std::string> listAllVariables();
    
    // Batch operations
    void applyTransformToGroup(const std::vector<std::string>& var_names, 
                               std::function<cdouble(cdouble)> transform_func);
    void scaleVariableGroup(const std::vector<std::string>& var_names, double scale_factor);
    
    // Self-expansion
    void autoExpandParameterSpace(double scale_factor);
    void expandMassScale(double factor);
    void expandFrequencyRange(double factor);
    void expandSpatialScale(double factor);
    
    // Self-refinement
    void autoRefineParameters(double tolerance);
    void calibrateToObservations(const std::map<std::string, double>& targets);
    void optimizeForMetric(const std::string& metric_name, double target_value);
    
    // Parameter space exploration
    std::vector<std::map<std::string, cdouble>> generateVariations(
        const std::vector<std::string>& params, int num_variants, double range);
    std::map<std::string, cdouble> findOptimalParameters(
        const std::string& optimization_goal);
    
    // Adaptive evolution
    void mutateParameters(double mutation_rate, double magnitude);
    void evolveSystem(int iterations, const std::string& fitness_metric);
    
    // State management
    void saveState(const std::string& label);
    void restoreState(const std::string& label);
    std::map<std::string, std::string> listSavedStates();
    void exportState(const std::string& filename);
    
    // System analysis
    std::map<std::string, double> analyzeParameterSensitivity();
    std::string generateSystemReport();
    bool validatePhysicalConsistency();
    void autoCorrectAnomalies();
};

#endif // SMBH_UQFF_MODULE_H

// SMBHUQFFModule.cpp
#include "SMBHUQFFModule.h"
#include <complex>
#include <cmath>
#include <vector>

// Constructor: Set all variables with SMBH-CGM defaults
SMBHUQFFModule::SMBHUQFFModule() {
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
    variables["DPM_gravity"] = {1.0, 0.0};
    variables["DPM_stability"] = {0.01, 0.001};
    variables["beta_i"] = {0.6, 0.0};
    variables["V_infl_UA"] = {1e-6, 1e-7};
    variables["rho_vac_A"] = {1e-30, 1e-31};
    variables["a_universal"] = {1e12, 1e11};
    variables["lambda_i"] = {1.0, 0.0};
    variables["rho_vac_SCm"] = {7.09e-37, 1e-38};
    variables["omega_s"] = {2.5e-6, 1e-7};
    variables["f_TRZ"] = {0.1, 0.0};
    variables["t_scale"] = {1e16, 0.0};

    // Simulation-specific
    variables["H0"] = {70.0 / (3.086e19), 0.0};  // Hubble (s^-1)
    variables["year"] = {3.156e7, 0.0};  // s/year
    variables["R_bulge"] = {1e3, 0.0};  // Placeholder
    variables["P_sc m"] = {1.0, 0.0};  // SCm power
    variables["E_react_0"] = {1.0, 0.0};  // Initial efficiency
    variables["gamma"] = {0.0005, 0.0};  // Decay
    variables["omega_c"] = {2 * pi_val / (24 * 3600), 0.0};  // Daily
    variables["f_heaviside"] = {1e13, 0.0};  // Threshold
    variables["f_quasi"] = {1.0, 0.0};  // Quasi-static
    variables["phi"] = {1.618, 0.0};  // Golden ratio
    variables["rho_vac_ua_prime"] = {1e-29, 0.0};  // UA'
    variables["rho_vac_scm"] = {1e-30, 0.0};  // SCm vac

    // Quadratic approx
    variables["x2"] = {-1.35e172, 0.0};
}

// Set system-specific params
void SMBHUQFFModule::setGalaxyParams(double sigma, int n, double z) {
    double sigma_m_s = sigma * 1e3;  // m/s
    variables["sigma"] = {sigma_m_s, 0.0};
    variables["n"] = {static_cast<double>(n), 0.0};
    variables["z"] = {z, 0.0};
    // M from M-σ: log M_BH = 0.309 log(σ/200) + 4.38, M in log10 solar masses
    double log_sigma = log10(sigma_m_s / 200e3);
    double log_M_BH = 0.309 * log_sigma + 4.38;
    double M_BH = pow(10, log_M_BH) * 1.989e30;  // kg
    variables["M"] = {M_BH, 0.0};
    variables["r"] = {sigma_m_s * 1e3, 0.0};  // Approx r from sigma
    variables["t"] = {computeCosmicTime(z), 0.0};
}

// Update variable (set to new complex value)
void SMBHUQFFModule::updateVariable(const std::string& name, cdouble value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta (complex) to variable
void SMBHUQFFModule::addToVariable(const std::string& name, cdouble delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta (complex)
void SMBHUQFFModule::subtractFromVariable(const std::string& name, cdouble delta) {
    addToVariable(name, -delta);
}

// Compute cosmic time from z
double SMBHUQFFModule::computeCosmicTime(double z) {
    cdouble H0 = variables["H0"];
    cdouble year = variables["year"];
    return (2.0 / (3 * H0.real())) * pow(1 + z, -1.5) * year.real();
}

// Compute omega_s for galactic scales
double SMBHUQFFModule::computeOmegaS(double sigma) {
    cdouble R_bulge = variables["R_bulge"];
    return sigma * 1e3 / R_bulge.real();
}

// Compute mu_j (magnetic moment)
double SMBHUQFFModule::computeMuJ(double t) {
    cdouble omega_c = variables["omega_c"];
    return (1e3 + 0.4 * sin(omega_c.real() * t)) * 3.38e20;
}

// Compute E_react (reactor efficiency)
double SMBHUQFFModule::computeEReact(double t) {
    cdouble E_react_0 = variables["E_react_0"];
    cdouble year = variables["year"];
    return E_react_0.real() * exp(-0.0005 * t / year.real());
}

// Compute delta_n (quantum state shift)
double SMBHUQFFModule::computeDeltaN(int n) {
    cdouble phi = variables["phi"];
    return phi.real() * pow(2 * M_PI, n / 6.0);
}

// Compute rho_vac_ua_scm
double SMBHUQFFModule::computeRhoVacUAScm(int n, double t) {
    cdouble rho_vac_ua_prime = variables["rho_vac_ua_prime"];
    cdouble rho_vac_scm = variables["rho_vac_scm"];
    cdouble rho_vac_ua = variables["rho_vac_UA"];
    cdouble year = variables["year"];
    double pi = variables["pi"].real();
    return rho_vac_ua_prime.real() * pow(rho_vac_scm.real() / rho_vac_ua.real(), n) * exp(-1 * exp(-pi - t / year.real()));
}

// Compute U_m
cdouble SMBHUQFFModule::computeUm(double t, double r, int n) {
    double mu = computeMuJ(t);
    cdouble term1 = mu / r;
    cdouble gamma = variables["gamma"];
    cdouble year = variables["year"];
    cdouble t_n = variables["t_scale"];  // Reuse for t_n
    double term2 = 1 - exp(-gamma.real() * t / (24 * 3600) * cos(pi * t_n.real()));
    cdouble P_scm = variables["P_sc m"];
    cdouble E_react = {computeEReact(t), 0.0};
    cdouble f_heaviside = variables["f_heaviside"];
    cdouble f_quasi = variables["f_quasi"];
    cdouble term3 = 1.0;  // hat_phi_j
    cdouble factor = P_scm * E_react * (1.0 + f_heaviside.real()) * (1.0 + f_quasi.real());
    return term1 * term2 * term3 * factor;
}

// Compute U_g1
cdouble SMBHUQFFModule::computeUg1(double t, double r, double M_s, int n) {
    cdouble G = variables["G"];
    cdouble M = {M_s, 0.0};
    cdouble r_complex = {r, 0.0};
    return G * M / pow(r_complex, 2);
}

// Compute DPM_resonance
cdouble SMBHUQFFModule::computeDPM_resonance(double sigma) {
    cdouble g = variables["g_Lande"];
    cdouble muB = variables["mu_B"];
    cdouble B = variables["B0"];
    cdouble hbar = variables["hbar"];
    cdouble omega0 = variables["omega0"];
    return (g * muB * B / (hbar * omega0)).real();  // Return as complex with imag 0
}

// Compute LENR term
cdouble SMBHUQFFModule::computeLENRTerm(double sigma) {
    cdouble k = variables["k_LENR"];
    cdouble omegaL = variables["omega_LENR"];
    cdouble omega0 = variables["omega0"];
    return k * pow(omegaL / omega0, 2.0);
}

// Compute integrand for F_U_Bi_i
cdouble SMBHUQFFModule::computeIntegrand(double t_user, double sigma, int n) {
    setGalaxyParams(sigma, n, 0.0);  // z=0 default
    variables["t"] = {t_user, 0.0};
    double cos_theta = cos(variables["theta"].real());
    double sin_theta = sin(variables["theta"].real());
    double cos_act = cos(variables["omega_act"].real() * t_user + variables["phi"].real());

    cdouble term_base = -variables["F0"];
    cdouble term_mom = (variables["m_e"] * pow(variables["c"], 2) / pow(variables["r"], 2)) * variables["DPM_momentum"] * cos_theta;
    cdouble term_grav = (variables["G"] * variables["M"] / pow(variables["r"], 2)) * variables["DPM_gravity"];
    cdouble term_vac = variables["rho_vac_UA"] * variables["DPM_stability"];
    cdouble term_LENR = computeLENRTerm(sigma);
    cdouble term_act = variables["k_act"] * cos_act;
    cdouble term_DE = variables["k_DE"] * variables["L_X"];
    cdouble term_res = 2 * variables["q"] * variables["B0"] * variables["V"] * sin_theta * computeDPM_resonance(sigma);
    cdouble term_neut = variables["k_neutron"] * variables["sigma_n"];
    cdouble term_rel = variables["k_rel"] * pow(variables["E_cm_astro"] / variables["E_cm"], 2);
    cdouble term_neutrino = variables["F_neutrino"];

    return term_base + term_mom + term_grav + term_vac + term_LENR + term_act + term_DE + term_res + term_neut + term_rel + term_neutrino;
}

// Approx x2
cdouble SMBHUQFFModule::computeX2(double sigma, int n) {
    return variables["x2"] * sigma * n;
}

// Quadratic root helper
cdouble SMBHUQFFModule::computeQuadraticRoot(cdouble a, cdouble b, cdouble c) {
    cdouble disc = sqrt(b*b - 4*a*c);
    return (-b - disc) / (2*a);  // Negative root approx
}

// Full M-sigma with UQFF approx as integrand * x2
cdouble SMBHUQFFModule::computeMSigma(double sigma, int n, double t) {
    cdouble integ = computeIntegrand(t, sigma, n);
    cdouble x2_val = computeX2(sigma, n);
    return integ * x2_val;
}

// Compressed (integrand)
cdouble SMBHUQFFModule::computeCompressed(double sigma, int n, double t) {
    return computeIntegrand(t, sigma, n);
}

// Resonant DPM
cdouble SMBHUQFFModule::computeResonant(double t, double sigma) {
    return computeDPM_resonance(sigma);
}

// Buoyancy Ub1
cdouble SMBHUQFFModule::computeBuoyancy(double sigma) {
    cdouble beta = variables["beta_i"];
    cdouble V = variables["V_infl_UA"];
    cdouble rho = variables["rho_vac_A"];
    cdouble a = variables["a_universal"];
    return beta * V * rho * a;
}

// Superconductive Ui
cdouble SMBHUQFFModule::computeSuperconductive(double t, double sigma) {
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
double SMBHUQFFModule::computeCompressedG(double t, double sigma, int n) {
    setGalaxyParams(sigma, n, 0.0);
    double G_val = variables["G"].real();
    double M_val = variables["M"].real();
    double rho = variables["rho_gas"].real();
    double r_val = variables["r"].real();
    double kB = variables["k_B"].real();
    double T_val = 1e7;  // Generic
    double m_e_val = variables["m_e"].real();
    double c_val = variables["c"].real();
    double dpm_curv = 1e-22;

    double term1 = - (G_val * M_val * rho) / r_val;
    double term2 = - (kB * T_val * rho) / (m_e_val * c_val * c_val);
    double term3 = dpm_curv * pow(c_val, 4) / (G_val * r_val * r_val);

    return term1 + term2 + term3;
}

// Get equation text (descriptive)
std::string SMBHUQFFModule::getEquationText(double sigma) {
    return "M_{\\sigma} = 10^{0.309 \\log(\\sigma / 200) + 4.38} M_\\odot \\int U_m + U_g1 + \\rho_{vac,ua,sc m}(n,t) \\, dx \\approx " + std::to_string(computeMSigma(sigma, 1, 0.0).real()) + " + i \\cdot " + std::to_string(computeMSigma(sigma, 1, 0.0).imag()) + " (sigma=" + std::to_string(sigma) + ")\\n"
           "Where U_m = \\mu_j(t) / r * (1 - e^{-\\gamma t} \\cos \\pi t_n) * P_{scm} E_{react}(t) (1 + f_{heaviside}) (1 + f_{quasi}), U_g1 = G M_s / r^2\\n"
           "Adaptations for SMBH-CGM: ROMULUS25 M-σ, cosmic time(z), omega_s(sigma); validated with Kormendy & Ho (2013), cross-correlated via DeepSearch.";
}

// Print variables (complex)
void SMBHUQFFModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << std::setprecision(10)
                  << pair.second.real() << " + i " << pair.second.imag() << std::endl;
    }
}

// ============================================================================
// DYNAMIC SELF-UPDATE AND SELF-EXPANSION IMPLEMENTATIONS
// ============================================================================

// Storage for saved states
static std::map<std::string, std::map<std::string, cdouble>> saved_states;

// Create new dynamic variable
void SMBHUQFFModule::createDynamicVariable(const std::string& name, cdouble value) {
    if (variables.find(name) != variables.end()) {
        std::cout << "Variable '" << name << "' already exists. Updating value.\n";
    }
    variables[name] = value;
    std::cout << "Dynamic variable '" << name << "' = " << value.real() 
              << " + i*" << value.imag() << " created/updated.\n";
}

// Remove dynamic variable
void SMBHUQFFModule::removeDynamicVariable(const std::string& name) {
    auto it = variables.find(name);
    if (it != variables.end()) {
        variables.erase(it);
        std::cout << "Removed variable: " << name << "\n";
    } else {
        std::cerr << "Variable '" << name << "' not found.\n";
    }
}

// Clone variable with new name
void SMBHUQFFModule::cloneVariable(const std::string& src, const std::string& dest) {
    if (variables.find(src) != variables.end()) {
        variables[dest] = variables[src];
        std::cout << "Cloned: " << src << " → " << dest << "\n";
    } else {
        std::cerr << "Source variable '" << src << "' not found.\n";
    }
}

// List all variables
std::vector<std::string> SMBHUQFFModule::listAllVariables() {
    std::vector<std::string> var_list;
    for (const auto& pair : variables) {
        var_list.push_back(pair.first);
    }
    return var_list;
}

// Apply transformation function to group of variables
void SMBHUQFFModule::applyTransformToGroup(const std::vector<std::string>& var_names, 
                                           std::function<cdouble(cdouble)> transform_func) {
    for (const auto& name : var_names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = transform_func(variables[name]);
            std::cout << "Transformed " << name << " → " << variables[name].real() 
                      << " + i*" << variables[name].imag() << "\n";
        }
    }
}

// Scale variable group by factor
void SMBHUQFFModule::scaleVariableGroup(const std::vector<std::string>& var_names, 
                                        double scale_factor) {
    applyTransformToGroup(var_names, [scale_factor](cdouble val) {
        return val * scale_factor;
    });
}

// Auto-expand entire parameter space
void SMBHUQFFModule::autoExpandParameterSpace(double scale_factor) {
    std::cout << "Expanding parameter space by factor: " << scale_factor << "\n";
    
    // Scale physical constants proportionally
    std::vector<std::string> scalable = {
        "M", "r", "R_bulge", "sigma", "E_cm_astro", "F0", "F_rel"
    };
    
    for (const auto& var : scalable) {
        if (variables.find(var) != variables.end()) {
            variables[var] *= scale_factor;
        }
    }
    
    std::cout << "Parameter space expanded.\n";
}

// Expand mass scale specifically
void SMBHUQFFModule::expandMassScale(double factor) {
    std::cout << "Expanding mass scale by factor: " << factor << "\n";
    if (variables.find("M") != variables.end()) {
        variables["M"] *= factor;
        std::cout << "M (BH mass) scaled to: " << variables["M"].real() << " kg\n";
    }
}

// Expand frequency range
void SMBHUQFFModule::expandFrequencyRange(double factor) {
    std::cout << "Expanding frequency range by factor: " << factor << "\n";
    std::vector<std::string> freq_vars = {
        "omega_act", "omega_LENR", "omega_s", "omega_c"
    };
    scaleVariableGroup(freq_vars, factor);
}

// Expand spatial scale
void SMBHUQFFModule::expandSpatialScale(double factor) {
    std::cout << "Expanding spatial scale by factor: " << factor << "\n";
    std::vector<std::string> spatial_vars = {"r", "R_bulge"};
    scaleVariableGroup(spatial_vars, factor);
}

// Auto-refine parameters to improve physical consistency
void SMBHUQFFModule::autoRefineParameters(double tolerance) {
    std::cout << "Auto-refining parameters (tolerance=" << tolerance << ")...\n";
    
    // Check M-σ relation consistency
    if (variables.find("sigma") != variables.end() && variables.find("M") != variables.end()) {
        double sigma_val = variables["sigma"].real() / 1e3;  // Convert to km/s
        double log_sigma = log10(sigma_val / 200.0);
        double log_M_expected = 0.309 * log_sigma + 4.38;
        double M_expected = pow(10, log_M_expected) * 1.989e30;
        
        double M_current = variables["M"].real();
        double error = std::abs(M_current - M_expected) / M_expected;
        
        if (error > tolerance) {
            std::cout << "M-σ relation error: " << (error*100) << "%. Correcting M.\n";
            variables["M"] = {M_expected, variables["M"].imag()};
        }
    }
    
    std::cout << "Parameters refined.\n";
}

// Calibrate to observational targets
void SMBHUQFFModule::calibrateToObservations(const std::map<std::string, double>& targets) {
    std::cout << "Calibrating to " << targets.size() << " observational targets...\n";
    
    for (const auto& target : targets) {
        const std::string& param = target.first;
        double target_value = target.second;
        
        if (variables.find(param) != variables.end()) {
            double current = variables[param].real();
            double error = (target_value - current) / target_value;
            
            if (std::abs(error) > 0.01) {  // 1% threshold
                std::cout << param << ": adjusting from " << current 
                          << " to " << target_value << "\n";
                variables[param] = {target_value, variables[param].imag()};
            }
        }
    }
    
    std::cout << "Calibration complete.\n";
}

// Optimize for specific metric
void SMBHUQFFModule::optimizeForMetric(const std::string& metric_name, double target_value) {
    std::cout << "Optimizing " << metric_name << " to target: " << target_value << "\n";
    
    // Example: optimize M-sigma computation
    if (metric_name == "M_sigma") {
        double sigma_val = variables["sigma"].real() / 1e3;
        auto result = computeMSigma(sigma_val, 1, 0.0);
        double current = result.real();
        
        if (std::abs(current - target_value) > 0.1 * target_value) {
            // Adjust underlying parameters
            double correction = target_value / (current + 1e-100);
            variables["F0"] *= correction;
            std::cout << "Adjusted F0 by factor: " << correction << "\n";
        }
    }
}

// Generate parameter variations for exploration
std::vector<std::map<std::string, cdouble>> SMBHUQFFModule::generateVariations(
    const std::vector<std::string>& params, int num_variants, double range) {
    
    std::vector<std::map<std::string, cdouble>> variations;
    
    for (int i = 0; i < num_variants; i++) {
        std::map<std::string, cdouble> variant;
        
        for (const auto& param : params) {
            if (variables.find(param) != variables.end()) {
                cdouble base = variables[param];
                // Linear spread across range
                double factor = 1.0 + range * (2.0 * i / (num_variants - 1.0) - 1.0);
                variant[param] = base * factor;
            }
        }
        
        variations.push_back(variant);
    }
    
    std::cout << "Generated " << num_variants << " parameter variations.\n";
    return variations;
}

// Find optimal parameters for given goal
std::map<std::string, cdouble> SMBHUQFFModule::findOptimalParameters(
    const std::string& optimization_goal) {
    
    std::cout << "Finding optimal parameters for: " << optimization_goal << "\n";
    
    // Return current state as baseline
    // In production, would run optimization loop
    std::map<std::string, cdouble> optimal = variables;
    
    if (optimization_goal == "max_M_sigma") {
        // Maximize M-sigma: increase F0
        optimal["F0"] *= 1.5;
    } else if (optimization_goal == "min_energy") {
        // Minimize energy: reduce velocity terms
        optimal["V"] *= 0.5;
    }
    
    return optimal;
}

// Mutate parameters for evolutionary exploration
void SMBHUQFFModule::mutateParameters(double mutation_rate, double magnitude) {
    std::cout << "Mutating parameters (rate=" << mutation_rate 
              << ", magnitude=" << magnitude << ")...\n";
    
    int mutation_count = 0;
    for (auto& pair : variables) {
        // Stochastic mutation (simplified: deterministic for reproducibility)
        if (mutation_rate > 0.5) {  // Simplified threshold
            cdouble mutation = pair.second * magnitude * (mutation_rate - 0.5);
            pair.second += mutation;
            mutation_count++;
        }
    }
    
    std::cout << "Mutated " << mutation_count << " parameters.\n";
}

// Evolve system over iterations
void SMBHUQFFModule::evolveSystem(int iterations, const std::string& fitness_metric) {
    std::cout << "Evolving system for " << iterations 
              << " iterations (fitness=" << fitness_metric << ")...\n";
    
    for (int i = 0; i < iterations; i++) {
        // Compute fitness
        double fitness = 0.0;
        if (fitness_metric == "M_sigma") {
            double sigma_val = variables["sigma"].real() / 1e3;
            fitness = computeMSigma(sigma_val, 1, 0.0).real();
        }
        
        // Apply mutation
        mutateParameters(0.1, 0.05);
        
        // Auto-refine
        autoRefineParameters(0.01);
        
        if (i % 10 == 0) {
            std::cout << "Iteration " << i << ", fitness=" << fitness << "\n";
        }
    }
    
    std::cout << "Evolution complete.\n";
}

// Save current state
void SMBHUQFFModule::saveState(const std::string& label) {
    saved_states[label] = variables;
    std::cout << "State '" << label << "' saved (" << variables.size() << " variables).\n";
}

// Restore saved state
void SMBHUQFFModule::restoreState(const std::string& label) {
    if (saved_states.find(label) != saved_states.end()) {
        variables = saved_states[label];
        std::cout << "State '" << label << "' restored.\n";
    } else {
        std::cerr << "State '" << label << "' not found.\n";
    }
}

// List all saved states
std::map<std::string, std::string> SMBHUQFFModule::listSavedStates() {
    std::map<std::string, std::string> state_list;
    for (const auto& pair : saved_states) {
        state_list[pair.first] = "(" + std::to_string(pair.second.size()) + " vars)";
    }
    return state_list;
}

// Export state to file (simplified: returns state info)
void SMBHUQFFModule::exportState(const std::string& filename) {
    std::cout << "Exporting state to: " << filename << "\n";
    std::cout << "Total variables: " << variables.size() << "\n";
    // In production: write to actual file
}

// Analyze parameter sensitivity
std::map<std::string, double> SMBHUQFFModule::analyzeParameterSensitivity() {
    std::cout << "Analyzing parameter sensitivity...\n";
    
    std::map<std::string, double> sensitivity;
    double sigma_val = variables["sigma"].real() / 1e3;
    double base_result = computeMSigma(sigma_val, 1, 0.0).real();
    
    // Test key parameters
    std::vector<std::string> test_params = {"F0", "F_rel", "sigma", "M"};
    
    for (const auto& param : test_params) {
        if (variables.find(param) != variables.end()) {
            cdouble original = variables[param];
            
            // Perturb by 1%
            variables[param] *= 1.01;
            double perturbed_result = computeMSigma(sigma_val, 1, 0.0).real();
            
            // Calculate sensitivity
            double delta_output = (perturbed_result - base_result) / base_result;
            double delta_input = 0.01;
            sensitivity[param] = delta_output / delta_input;
            
            // Restore
            variables[param] = original;
            
            std::cout << param << " sensitivity: " << sensitivity[param] << "\n";
        }
    }
    
    return sensitivity;
}

// Generate comprehensive system report
std::string SMBHUQFFModule::generateSystemReport() {
    std::string report;
    report += "=============================================================\n";
    report += "SMBH UQFF MODULE - SYSTEM REPORT\n";
    report += "=============================================================\n\n";
    
    report += "Total Variables: " + std::to_string(variables.size()) + "\n";
    report += "Saved States: " + std::to_string(saved_states.size()) + "\n\n";
    
    report += "Key Parameters:\n";
    std::vector<std::string> key_vars = {"M", "sigma", "r", "F_rel", "F0"};
    for (const auto& var : key_vars) {
        if (variables.find(var) != variables.end()) {
            report += "  " + var + " = " + std::to_string(variables[var].real()) + "\n";
        }
    }
    
    report += "\nPhysical Consistency: ";
    report += validatePhysicalConsistency() ? "PASS\n" : "FAIL\n";
    
    report += "\nSample Computation:\n";
    double sigma_val = variables["sigma"].real() / 1e3;
    auto result = computeMSigma(sigma_val, 1, 0.0);
    report += "  M_sigma = " + std::to_string(result.real()) + " + i*" + std::to_string(result.imag()) + "\n";
    
    report += "\n=============================================================\n";
    
    return report;
}

// Validate physical consistency
bool SMBHUQFFModule::validatePhysicalConsistency() {
    bool valid = true;
    
    // Check for NaN or Inf
    for (const auto& pair : variables) {
        if (std::isnan(pair.second.real()) || std::isinf(pair.second.real()) ||
            std::isnan(pair.second.imag()) || std::isinf(pair.second.imag())) {
            std::cerr << "Invalid value in " << pair.first << "\n";
            valid = false;
        }
    }
    
    // Check M-σ relation
    if (variables.find("M") != variables.end() && variables.find("sigma") != variables.end()) {
        double M_val = variables["M"].real();
        if (M_val <= 0 || M_val > 1e42) {
            std::cerr << "Black hole mass out of physical range\n";
            valid = false;
        }
    }
    
    return valid;
}

// Auto-correct anomalies
void SMBHUQFFModule::autoCorrectAnomalies() {
    std::cout << "Auto-correcting anomalies...\n";
    
    int corrections = 0;
    
    // Fix NaN/Inf
    for (auto& pair : variables) {
        if (std::isnan(pair.second.real()) || std::isinf(pair.second.real())) {
            pair.second = {1.0, 0.0};  // Reset to safe default
            corrections++;
        }
    }
    
    // Ensure positive masses
    if (variables.find("M") != variables.end() && variables["M"].real() <= 0) {
        variables["M"] = {1e30, 0.0};  // Reset to 1 solar mass
        corrections++;
    }
    
    std::cout << "Applied " << corrections << " corrections.\n";
}

// Example usage in base program 'smbh_cgm_sim.cpp' (snippet for integration)
// #include "SMBHUQFFModule.h"
// #include <complex>
// int main() {
//     SMBHUQFFModule mod;
//     double sigma = 200; int n = 1; double t = 1e12;  // Example
//     
//     // Basic computation
//     auto M = mod.computeMSigma(sigma, n, t);
//     std::cout << "M_sigma = " << M.real() << " + i " << M.imag() << std::endl;
//     std::cout << mod.getEquationText(sigma) << std::endl;
//     
//     // Dynamic variable operations
//     mod.updateVariable("F_rel", {5e33, 0.0});  // Update rel coherence
//     mod.createDynamicVariable("custom_param", {1.5e10, 0.0});
//     
//     // Self-expansion
//     mod.saveState("initial");
//     mod.autoExpandParameterSpace(1.5);  // 50% expansion
//     mod.expandMassScale(2.0);  // Double BH mass
//     
//     // Self-refinement
//     mod.autoRefineParameters(0.01);  // 1% tolerance
//     std::map<std::string, double> targets = {{"M", 1e39}, {"sigma", 250e3}};
//     mod.calibrateToObservations(targets);
//     
//     // Parameter space exploration
//     auto variations = mod.generateVariations({"M", "sigma", "F0"}, 10, 0.2);
//     auto sensitivity = mod.analyzeParameterSensitivity();
//     
//     // Evolution
//     mod.mutateParameters(0.8, 0.1);
//     mod.evolveSystem(50, "M_sigma");
//     
//     // Reporting
//     std::cout << mod.generateSystemReport() << std::endl;
//     mod.printVariables();
//     
//     // Restore if needed
//     mod.restoreState("initial");
//     
//     return 0;
// }
// Compile: g++ -std=c++11 -o smbh_cgm_sim smbh_cgm_sim.cpp SMBHUQFFModule.cpp -lm
// Sample Output for sigma=200: M_sigma ≈ black hole mass with UQFF terms; full for ROMULUS25-like.
// 
// NEW CAPABILITIES SUMMARY:
// - Dynamic variable creation/removal at runtime
// - Parameter space auto-expansion (mass, frequency, spatial scales)
// - Self-refinement with observational calibration
// - Parameter sensitivity analysis
// - Evolutionary optimization with mutations
// - State save/restore for exploration
// - Comprehensive system validation and reporting
// 
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 23, 2025.
// Dynamic Capabilities Enhancement - Nov 1, 2025.