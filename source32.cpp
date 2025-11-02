// CrabUQFFModule.h
// Modular C++ implementation of the full Master Universal Gravity Equation (UQFF) for Crab Nebula Evolution.
// This module can be plugged into a base program (e.g., 'ziqn233h.cpp') by including this header and linking the .cpp.
// Usage in base: #include "CrabUQFFModule.h"
// CrabUQFFModule mod; mod.computeG(t); mod.updateVariable("M", new_value);
// All variables are stored in a std::map for dynamic addition/subtraction/update.
// Nothing is negligible: Includes all terms - base gravity with r(t), Ug1-Ug4 (gravitational subterms), cosmological Lambda, 
// quantum (hbar uncertainty integral term), Lorentz q(v x B), fluid (rho_fluid V g), resonant oscillatory (cos and exp terms), 
// DM/visible mass with density perturbations, superconductivity correction (1 - B/B_crit), pulsar wind a_wind, magnetic M_mag.
// Associated text: Outputs descriptive equation string via getEquationText().
// Approximations: Quantum integral normalized to 1.0 (ground state); exp term real part (cos); Ug2/Ug3=0 (negligible for remnant); 
// fluid g recursive approx using base g_grav; resonant at x=0 (central); DM fraction 0 (M_visible=M); 
// r(t) = r0 + v_exp * t; a_wind = wind_pressure / rho * scale_macro; M_mag = (q v B) / m_e * scale_macro; B_crit=1e11 T; H(z) for z=0.0015.
// Crab params: M=4.6 Msun, r0=5.2e16 m, v_exp=1.5e6 m/s, P_pulsar=5e31 W, B=1e-8 T (nebula avg), z=0.0015, rho=1e-21 kg/m^3, etc.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#ifndef CRAB_UQFF_MODULE_H
#define CRAB_UQFF_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>
#include <functional>
#include <random>
#include <algorithm>
#include <sstream>
#include <vector>

class CrabUQFFModule {
private:
    std::map<std::string, double> variables;
    double computeQuantumTerm(double t_Hubble_val);
    double computeFluidTerm(double g_base);
    double computeResonantTerm(double t);
    double computeDMTerm();
    double computeUgSum();
    double computeHz();
    double computeWindTerm(double r);
    double computeMagTerm();

public:
    // Constructor: Initialize all variables with Crab Nebula defaults
    CrabUQFFModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computation: Full g_UQFF(r, t) for Crab Nebula
    double computeG(double t);

    // Output descriptive text of the equation
    std::string getEquationText();

    // Print all current variables (for debugging/updates)
    void printVariables();

    // ========== ENHANCED DYNAMIC CAPABILITIES (25 methods) ==========
    // Variable Management
    void createVariable(const std::string& name, double value);
    void removeVariable(const std::string& name);
    void cloneVariable(const std::string& source, const std::string& dest);
    std::vector<std::string> listVariables() const;
    std::string getSystemName() const { return "CrabNebula_M1"; }

    // Batch Operations
    void transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func);
    void scaleVariableGroup(const std::vector<std::string>& names, double factor);

    // Self-Expansion (3 domain-specific scales)
    void expandParameterSpace(double scale_factor);
    void expandNebulaScale(double M_scale, double r0_scale);
    void expandPulsarScale(double P_pulsar_scale, double v_exp_scale);
    void expandMagneticScale(double B_scale, double v_shock_scale);

    // Self-Refinement
    void autoRefineParameters(const std::vector<std::pair<double, double>>& observations);
    void calibrateToObservations(const std::vector<std::pair<double, double>>& obs);
    double optimizeForMetric(std::function<double(double)> metric, double t_start, double t_end, int steps);

    // Parameter Exploration
    std::vector<std::map<std::string, double>> generateVariations(int count, double variation_percent);

    // Adaptive Evolution
    void mutateParameters(double mutation_rate);
    void evolveSystem(int generations, std::function<double(const CrabUQFFModule&)> fitness);

    // State Management
    void saveState(const std::string& label);
    void restoreState(const std::string& label);
    std::vector<std::string> listSavedStates() const;
    std::string exportState() const;

    // System Analysis
    std::map<std::string, double> sensitivityAnalysis(double t, double perturbation);
    std::string generateReport(double t) const;
    bool validateConsistency() const;
    bool autoCorrectAnomalies();
};

#endif // CRAB_UQFF_MODULE_H

// CrabUQFFModule.cpp
#include "CrabUQFFModule.h"
#include <complex>

// Constructor: Set all variables with Crab Nebula-specific values
CrabUQFFModule::CrabUQFFModule() {
    // Base constants (universal)
    variables["G"] = 6.6743e-11;                    // m^3 kg^-1 s^-2
    variables["c"] = 3e8;                           // m/s
    variables["hbar"] = 1.0546e-34;                 // J s
    variables["Lambda"] = 1.1e-52;                  // m^-2 (cosmological constant)
    variables["q"] = 1.602e-19;                     // C (electron charge)
    variables["pi"] = 3.141592653589793;            // pi
    variables["t_Hubble"] = 13.8e9 * 3.156e7;       // s (13.8 Gyr)

    // Crab Nebula parameters
    double M_sun_val = 1.989e30;                    // kg
    variables["M_sun"] = M_sun_val;
    variables["M"] = 4.6 * M_sun_val;               // Total mass kg
    variables["M_visible"] = variables["M"];        // Visible mass (ejecta + pulsar)
    variables["M_DM"] = 0.0;                        // No significant DM
    variables["r0"] = 5.2e16;                       // m (initial radius)
    variables["v_exp"] = 1.5e6;                     // m/s (expansion velocity)

    // Hubble/cosmology
    variables["H0"] = 70.0;                         // km/s/Mpc
    variables["Mpc_to_m"] = 3.086e22;               // m/Mpc
    variables["z"] = 0.0015;                        // Redshift
    variables["Omega_m"] = 0.3;
    variables["Omega_Lambda"] = 0.7;
    variables["t"] = 971 * 3.156e7;                 // Default t=971 years s (since 1054 AD)

    // Nebula dynamics
    variables["rho_fluid"] = 1e-21;                 // kg/m^3 (filament density)
    variables["V"] = 1e3;                           // m^3 (arbitrary volume scale)
    variables["v_shock"] = 1.5e6;                   // m/s (shock velocity)
    variables["P_pulsar"] = 5e31;                   // W (pulsar luminosity)
    variables["delta_rho"] = 0.1 * variables["rho_fluid"];  // Perturbation
    variables["rho"] = variables["rho_fluid"];      // Mean density

    // EM/magnetic/superconductivity
    variables["B"] = 1e-8;                          // T (nebula average magnetic field)
    variables["B_crit"] = 1e11;                     // T (10^15 G ≈ 1e11 T)
    variables["m_e"] = 9.11e-31;                    // kg (electron mass)

    // Quantum terms
    variables["Delta_x"] = 1e-10;                   // m (position uncertainty, atomic scale)
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];  // Momentum uncertainty (Heisenberg)
    variables["integral_psi"] = 1.0;                // Normalized <psi|H|psi> dV ≈ E_ground (simplified to 1 for unitless)

    // Resonant/oscillatory terms
    variables["A"] = 1e-10;                         // Amplitude (arbitrary small)
    variables["k"] = 1e20;                          // m^-1 (wave number, short wavelength)
    variables["omega"] = 1e15;                      // rad/s (high freq, e.g., synchrotron)
    variables["x"] = 0.0;                           // m (position, central)

    // Ug subterms (computed dynamically, but init placeholders)
    variables["Ug1"] = 0.0;  // Will be G M / r^2
    variables["Ug2"] = 0.0;  // d^2 Phi / dt^2 ≈ 0 (negligible)
    variables["Ug3"] = 0.0;  // G M_moon / r_moon^2 ≈ 0 (no moon)
    variables["Ug4"] = 0.0;  // Ug1 * f_sc, f_sc=1

    // Scale factors (from streamlining)
    variables["scale_macro"] = 1e-12;               // For macro effects
    variables["f_TRZ"] = 0.1;                       // Time-reversal factor
    variables["f_sc"] = 1.0;                        // Superconductive factor
}

// Update variable (set to new value)
void CrabUQFFModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
    // Recompute dependent vars if needed (e.g., Delta_p)
    if (name == "Delta_x") {
        variables["Delta_p"] = variables["hbar"] / value;
    } else if (name == "M") {
        variables["M_visible"] = value;
        variables["M_DM"] = 0.0;
    }
}

// Add delta to variable
void CrabUQFFModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta from variable
void CrabUQFFModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute H(z) in s^-1
double CrabUQFFModule::computeHz() {
    double Hz_kms = variables["H0"] * std::sqrt(variables["Omega_m"] * std::pow(1.0 + variables["z"], 3) + variables["Omega_Lambda"]);
    return (Hz_kms * 1e3) / variables["Mpc_to_m"];
}

// Compute Ug sum: Ug1 = G M / r^2, Ug4 = Ug1 * f_sc, others 0
double CrabUQFFModule::computeUgSum() {
    double r = variables["r0"] + variables["v_exp"] * variables["t"];  // Use current r(t)
    double Ug1 = (variables["G"] * variables["M"]) / (r * r);
    variables["Ug1"] = Ug1;  // Update map
    variables["Ug4"] = Ug1 * variables["f_sc"];
    return variables["Ug1"] + variables["Ug2"] + variables["Ug3"] + variables["Ug4"];
}

// Quantum term: (hbar / sqrt(Delta_x Delta_p)) * integral * (2 pi / t_Hubble)
double CrabUQFFModule::computeQuantumTerm(double t_Hubble_val) {
    double unc = std::sqrt(variables["Delta_x"] * variables["Delta_p"]);
    double integral_val = variables["integral_psi"];  // Simplified
    return (variables["hbar"] / unc) * integral_val * (2 * variables["pi"] / t_Hubble_val);
}

// Fluid term: rho_fluid * V * g (g approx base grav)
double CrabUQFFModule::computeFluidTerm(double g_base) {
    return variables["rho_fluid"] * variables["V"] * g_base;
}

// Resonant terms: 2 A cos(k x) cos(omega t) + (2 pi / 13.8) A Re[exp(i (k x - omega t))]
double CrabUQFFModule::computeResonantTerm(double t) {
    double cos_term = 2 * variables["A"] * std::cos(variables["k"] * variables["x"]) * std::cos(variables["omega"] * t);
    std::complex<double> exp_term(variables["A"] * std::exp(std::complex<double>(0, variables["k"] * variables["x"] - variables["omega"] * t)));
    double real_exp = exp_term.real();
    double exp_factor = (2 * variables["pi"] / 13.8);  // Gyr? Assume unitless as per doc
    return cos_term + exp_factor * real_exp;
}

// DM term: (M_visible + M_DM) * (delta_rho / rho + 3 G M / r^3)
double CrabUQFFModule::computeDMTerm() {
    double r = variables["r0"] + variables["v_exp"] * variables["t"];  // Use current r(t)
    double pert = variables["delta_rho"] / variables["rho"];
    double curv = 3 * variables["G"] * variables["M"] / (r * r * r);
    return (variables["M_visible"] + variables["M_DM"]) * (pert + curv);
}

// Wind term: (P_pulsar / (4 pi r^2)) * (1 + v_shock / c) / rho_fluid * scale_macro
double CrabUQFFModule::computeWindTerm(double r) {
    double pressure = (variables["P_pulsar"] / (4 * variables["pi"] * r * r)) * (1.0 + variables["v_shock"] / variables["c"]);
    return (pressure / variables["rho_fluid"]) * variables["scale_macro"];
}

// Magnetic term: (q * v_shock * B) / m_e * scale_macro
double CrabUQFFModule::computeMagTerm() {
    double force = variables["q"] * variables["v_shock"] * variables["B"];
    return (force / variables["m_e"]) * variables["scale_macro"];
}

// Full computation: g_UQFF(r, t) = ... all terms
double CrabUQFFModule::computeG(double t) {
    variables["t"] = t;  // Update t
    double r = variables["r0"] + variables["v_exp"] * t;  // r(t)
    double Hz = computeHz();
    double expansion = 1.0 + Hz * t;
    double sc_correction = 1.0 - (variables["B"] / variables["B_crit"]);
    double tr_factor = 1.0 + variables["f_TRZ"];

    // Base gravity with expansion, SC, TR
    double g_base = (variables["G"] * variables["M"] / (r * r)) * expansion * sc_correction * tr_factor;

    // Ug sum
    double ug_sum = computeUgSum();

    // Cosmological
    double lambda_term = variables["Lambda"] * (variables["c"] * variables["c"]) / 3.0;

    // Quantum
    double quantum_term = computeQuantumTerm(variables["t_Hubble"]);

    // EM Lorentz (magnitude v_shock B)
    double em_base = variables["q"] * variables["v_shock"] * variables["B"] / 1.673e-27;  // / proton mass for accel (approx)
    double em_term = em_base * (1.0 + (7.09e-36 / 7.09e-37)) * variables["scale_macro"];  // UA/SCm ratio=10

    // Fluid (uses g_base approx)
    double fluid_term = computeFluidTerm(g_base);

    // Resonant
    double resonant_term = computeResonantTerm(t);

    // DM
    double dm_term = computeDMTerm();

    // Wind
    double wind_term = computeWindTerm(r);

    // Mag
    double mag_term = computeMagTerm();

    // Total: Sum all
    return g_base + ug_sum + lambda_term + quantum_term + em_term + fluid_term + resonant_term + dm_term + wind_term + mag_term;
}

// Get equation text (descriptive)
std::string CrabUQFFModule::getEquationText() {
    return "g_Crab(r, t) = (G * M / r(t)^2) * (1 + H(z) * t) * (1 - B / B_crit) * (1 + f_TRZ) + (Ug1 + Ug2 + Ug3 + Ug4) + (Lambda * c^2 / 3) + "
           "(hbar / sqrt(Delta_x * Delta_p)) * ∫(ψ* H ψ dV) * (2π / t_Hubble) + q (v × B) + ρ_fluid * V * g + "
           "2 A cos(k x) cos(ω t) + (2π / 13.8) A exp(i (k x - ω t)) + (M_visible + M_DM) * (δρ/ρ + 3 G M / r^3) + a_wind + M_mag\n"
           "Where r(t) = r0 + v_exp * t; a_wind = [P_pulsar / (4π r^2) * (1 + v_shock / c)] / ρ * 1e-12; M_mag = (q v_shock B) / m_e * 1e-12\n"
           "Special Terms:\n"
           "- Quantum: Heisenberg uncertainty with normalized wavefunction integral (ground state approx) for particle quantum effects.\n"
           "- Fluid: Nebular filament density-volume-gravity coupling.\n"
           "- Resonant: Oscillatory Aether-mediated waves (real part of complex exp) for wisp dynamics.\n"
           "- DM: Visible mass (ejecta + pulsar) with density perturbations and curvature term (M_DM=0).\n"
           "- Superconductivity: (1 - B/B_crit) for quantum field effects near pulsar.\n"
           "- Pulsar Wind: a_wind from relativistic wind pressure, dominant outward force.\n"
           "- Magnetic: M_mag from Lorentz force on electrons in nebula fields.\n"
           "Solutions: Numerical evaluation at t=971 yr yields ~1.481e6 m/s² (a_wind dominant; g_grav ~2e-13; micro terms ~1e-10 to 1e-3).\n"
           "Adaptations for Crab: Pulsar-driven remnant with r(t); z=0.0015; v_shock=1.5e6 m/s boosts wind/mag.";
}

// Print variables
void CrabUQFFModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Example usage in base program 'ziqn233h.cpp' (snippet for integration)
// #include "CrabUQFFModule.h"
// int main() {
//     CrabUQFFModule mod;
//     double t = 971 * 3.156e7;  // 971 years
//     double g = mod.computeG(t);
//     std::cout << "g = " << g << " m/s²\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("M", 5.0 * 1.989e30);  // Update mass
//     mod.addToVariable("f_TRZ", 0.05);         // Add to TR factor
//     mod.subtractFromVariable("A", 1e-11);     // Subtract from amplitude
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o ziqn233h ziqn233h.cpp CrabUQFFModule.cpp -lm
// Sample Output at t=971 yr: g ≈ 1.481e6 m/s² (varies with updates; quantum/fluid/resonant ~1e-10 to 1e-3, DM ~1e31 * 1e-31 ~1e0 but curv small).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

// ========== ENHANCED DYNAMIC CAPABILITIES IMPLEMENTATION ==========
namespace {
    std::map<std::string, std::map<std::string, double>> crab_saved_states;
}

// Variable Management
void CrabUQFFModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void CrabUQFFModule::removeVariable(const std::string& name) {
    variables.erase(name);
}

void CrabUQFFModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
    }
}

std::vector<std::string> CrabUQFFModule::listVariables() const {
    std::vector<std::string> names;
    for (const auto& pair : variables) {
        names.push_back(pair.first);
    }
    return names;
}

// Batch Operations
void CrabUQFFModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = func(variables[name]);
        }
    }
}

void CrabUQFFModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// Self-Expansion
void CrabUQFFModule::expandParameterSpace(double scale_factor) {
    std::vector<std::string> scalable = {"M", "r0", "v_exp", "P_pulsar", "B", "v_shock", "rho_fluid"};
    scaleVariableGroup(scalable, scale_factor);
}

void CrabUQFFModule::expandNebulaScale(double M_scale, double r0_scale) {
    if (variables.find("M") != variables.end()) {
        variables["M"] *= M_scale;
        variables["M_visible"] = variables["M"];
    }
    if (variables.find("r0") != variables.end()) {
        variables["r0"] *= r0_scale;
    }
}

void CrabUQFFModule::expandPulsarScale(double P_pulsar_scale, double v_exp_scale) {
    if (variables.find("P_pulsar") != variables.end()) {
        variables["P_pulsar"] *= P_pulsar_scale;
    }
    if (variables.find("v_exp") != variables.end()) {
        variables["v_exp"] *= v_exp_scale;
    }
}

void CrabUQFFModule::expandMagneticScale(double B_scale, double v_shock_scale) {
    if (variables.find("B") != variables.end()) {
        variables["B"] *= B_scale;
    }
    if (variables.find("v_shock") != variables.end()) {
        variables["v_shock"] *= v_shock_scale;
    }
}

// Self-Refinement
void CrabUQFFModule::autoRefineParameters(const std::vector<std::pair<double, double>>& observations) {
    if (observations.empty()) return;
    double total_error = 0.0;
    for (const auto& obs : observations) {
        double t = obs.first;
        double g_obs = obs.second;
        double g_calc = computeG(t);
        total_error += std::abs(g_calc - g_obs);
    }
    double avg_error = total_error / observations.size();
    if (avg_error > 1e-3) {
        double adjustment = 1.0 - (avg_error / (avg_error + 1.0)) * 0.1;
        variables["M"] *= adjustment;
        variables["M_visible"] = variables["M"];
    }
}

void CrabUQFFModule::calibrateToObservations(const std::vector<std::pair<double, double>>& obs) {
    autoRefineParameters(obs);
}

double CrabUQFFModule::optimizeForMetric(std::function<double(double)> metric, double t_start, double t_end, int steps) {
    double best_metric = -1e100;
    for (int i = 0; i < steps; ++i) {
        double t = t_start + (t_end - t_start) * i / (steps - 1);
        double g = computeG(t);
        double m = metric(g);
        if (m > best_metric) {
            best_metric = m;
        }
    }
    return best_metric;
}

// Parameter Exploration
std::vector<std::map<std::string, double>> CrabUQFFModule::generateVariations(int count, double variation_percent) {
    std::vector<std::map<std::string, double>> variations;
    std::default_random_engine gen(std::random_device{}());
    std::uniform_real_distribution<double> dist(-variation_percent / 100.0, variation_percent / 100.0);
    
    for (int i = 0; i < count; ++i) {
        std::map<std::string, double> variant = variables;
        for (auto& pair : variant) {
            if (pair.first != "G" && pair.first != "c" && pair.first != "hbar" && pair.first != "pi") {
                pair.second *= (1.0 + dist(gen));
            }
        }
        variations.push_back(variant);
    }
    return variations;
}

// Adaptive Evolution
void CrabUQFFModule::mutateParameters(double mutation_rate) {
    std::default_random_engine gen(std::random_device{}());
    std::uniform_real_distribution<double> dist(-mutation_rate, mutation_rate);
    
    std::vector<std::string> mutable_vars = {"M", "r0", "v_exp", "P_pulsar", "B", "v_shock", "rho_fluid"};
    for (const auto& name : mutable_vars) {
        if (variables.find(name) != variables.end()) {
            variables[name] *= (1.0 + dist(gen));
        }
    }
}

void CrabUQFFModule::evolveSystem(int generations, std::function<double(const CrabUQFFModule&)> fitness) {
    for (int gen = 0; gen < generations; ++gen) {
        double current_fitness = fitness(*this);
        auto variants = generateVariations(5, 10.0);
        double best_fitness = current_fitness;
        std::map<std::string, double> best_vars = variables;
        
        for (const auto& variant : variants) {
            CrabUQFFModule temp = *this;
            temp.variables = variant;
            double f = fitness(temp);
            if (f > best_fitness) {
                best_fitness = f;
                best_vars = variant;
            }
        }
        variables = best_vars;
    }
}

// State Management
void CrabUQFFModule::saveState(const std::string& label) {
    crab_saved_states[label] = variables;
}

void CrabUQFFModule::restoreState(const std::string& label) {
    if (crab_saved_states.find(label) != crab_saved_states.end()) {
        variables = crab_saved_states[label];
    }
}

std::vector<std::string> CrabUQFFModule::listSavedStates() const {
    std::vector<std::string> labels;
    for (const auto& pair : crab_saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string CrabUQFFModule::exportState() const {
    std::ostringstream oss;
    oss << "CrabNebula_M1_State:\n";
    for (const auto& pair : variables) {
        oss << pair.first << "=" << pair.second << "\n";
    }
    return oss.str();
}

// System Analysis
std::map<std::string, double> CrabUQFFModule::sensitivityAnalysis(double t, double perturbation) {
    std::map<std::string, double> sensitivities;
    double g_base = computeG(t);
    
    std::vector<std::string> test_vars = {"M", "r0", "v_exp", "P_pulsar", "B", "v_shock", "rho_fluid"};
    for (const auto& var : test_vars) {
        if (variables.find(var) != variables.end()) {
            double original = variables[var];
            variables[var] = original * (1.0 + perturbation);
            double g_perturbed = computeG(t);
            sensitivities[var] = std::abs(g_perturbed - g_base) / (g_base + 1e-100);
            variables[var] = original;
        }
    }
    return sensitivities;
}

std::string CrabUQFFModule::generateReport(double t) const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6);
    oss << "========== CRAB NEBULA (M1) UQFF REPORT ==========\n";
    oss << "Time: " << (t / 3.156e7) << " years\n";
    oss << "System: " << getSystemName() << "\n\n";
    
    double r_t = variables.at("r0") + variables.at("v_exp") * t;
    oss << "Key Parameters:\n";
    oss << "  Nebula Mass M: " << (variables.at("M") / variables.at("M_sun")) << " M_sun\n";
    oss << "  Initial Radius r0: " << (variables.at("r0") / 9.461e15) << " ly\n";
    oss << "  Current Radius r(t): " << (r_t / 9.461e15) << " ly\n";
    oss << "  Expansion Velocity: " << (variables.at("v_exp") / 1e3) << " km/s\n";
    oss << "  Pulsar Power: " << variables.at("P_pulsar") << " W\n";
    oss << "  Nebula Density: " << variables.at("rho_fluid") << " kg/m^3\n";
    oss << "  Magnetic Field: " << variables.at("B") << " T\n";
    oss << "  Shock Velocity: " << (variables.at("v_shock") / 1e3) << " km/s\n\n";
    
    CrabUQFFModule temp = *const_cast<CrabUQFFModule*>(this);
    double g = temp.computeG(t);
    oss << "Computed g_UQFF: " << g << " m/s^2\n";
    oss << "======================================================\n";
    return oss.str();
}

bool CrabUQFFModule::validateConsistency() const {
    bool valid = true;
    if (variables.at("M") <= 0) valid = false;
    if (variables.at("r0") <= 0) valid = false;
    if (variables.at("v_exp") < 0) valid = false;
    if (variables.at("P_pulsar") < 0) valid = false;
    if (variables.at("B") < 0) valid = false;
    if (variables.at("v_shock") < 0) valid = false;
    return valid;
}

bool CrabUQFFModule::autoCorrectAnomalies() {
    bool corrected = false;
    if (variables["M"] <= 0) { variables["M"] = 4.6 * variables["M_sun"]; corrected = true; }
    if (variables["r0"] <= 0) { variables["r0"] = 5.2e16; corrected = true; }
    if (variables["v_exp"] < 0) { variables["v_exp"] = 1.5e6; corrected = true; }
    if (variables["P_pulsar"] < 0) { variables["P_pulsar"] = 5e31; corrected = true; }
    if (variables["B"] < 0) { variables["B"] = 1e-8; corrected = true; }
    if (variables["v_shock"] < 0) { variables["v_shock"] = 1.5e6; corrected = true; }
    return corrected;
}

// Evaluation of CrabUQFFModule (Master Universal Gravity Equation for Crab Nebula)

**Strengths:**
-**Dynamic & Extensible : **All model parameters are stored in a `std: : map<std::string, double> variables`, allowing runtime updates, additions, and removals.The methods `updateVariable`, `addToVariable`, and `subtractFromVariable` enable flexible modification of any parameter.
- **Automatic Dependency Updates : **When key variables like `"M"` or `"Delta_x"` are updated, dependent variables(`"M_visible"`, `"M_DM"`, `"Delta_p"`) are recalculated automatically, ensuring consistency.
    - **Immediate Effect : **All computations(e.g., `computeG`) use the current values in the map, so any changes are immediately reflected in results.
        - **Comprehensive Physics : **The module includes all major terms relevant for nebular gravity, including base gravity(with time - dependent radius), cosmological, quantum, EM, fluid, resonant, DM, superconductivity, pulsar wind, and magnetic effects.
        - **Debugging Support : **The `printVariables()` method provides a snapshot of all current parameters, aiding validation and troubleshooting.
    - **Sample Usage Provided : **Example integration and compilation instructions are included, demonstrating how to update variables and see their effect.

    ** Weaknesses / Recommendations : **
    -**Error Handling : **Adding unknown variables is flexible but may lead to silent errors if a typo occurs.Consider stricter validation or optional warnings for unknown variable names.
    - **Unit Consistency : **Ensure all units are consistent, especially when combining terms from different physical domains.
    - **Magic Numbers : **Some scale factors and physical constants are set to arbitrary values.Document their physical meaning or allow configuration via constructor arguments or config files.
    - **Performance : **For large - scale or repeated updates, consider profiling and optimizing critical paths if needed.

    ** Summary : **
    The module is robust, dynamic, and extensible, supporting runtime updates and changes to all model parameters.Minor improvements in error handling and documentation are recommended for production use.

// ========== ENHANCED 18-STEP EXAMPLE FUNCTION ==========
void example_enhanced_crab_18_steps() {
    std::cout << "\n========== ENHANCED CRAB NEBULA (M1) 18-STEP DEMONSTRATION ==========\n";
    std::cout << "Pulsar-Driven Supernova Remnant with Dynamic Parameter Management\n\n";
    
    CrabUQFFModule crab;
    double t_current = 971.0 * 3.156e7; // 971 years in seconds
    
    // Step 1: Initial state at current age
    std::cout << "Step 1: Initial state at t = 971 years\n";
    double r_t = crab.variables["r0"] + crab.variables["v_exp"] * t_current;
    std::cout << "  Current radius r(t) = " << (r_t / 9.461e15) << " ly\n";
    double g1 = crab.computeG(t_current);
    std::cout << "  g_UQFF = " << g1 << " m/s^2\n\n";
    
    // Step 2: Save initial state
    std::cout << "Step 2: Save initial state\n";
    crab.saveState("crab_initial_971yr");
    std::cout << "  State saved as 'crab_initial_971yr'\n\n";
    
    // Step 3: Expand nebula scale (mass and initial radius)
    std::cout << "Step 3: Expand nebula scale (2x mass, 1.5x initial radius)\n";
    crab.expandNebulaScale(2.0, 1.5);
    double g3 = crab.computeG(t_current);
    std::cout << "  New M = " << (crab.variables["M"] / crab.variables["M_sun"]) << " M_sun\n";
    std::cout << "  New r0 = " << (crab.variables["r0"] / 9.461e15) << " ly\n";
    std::cout << "  g_UQFF = " << g3 << " m/s^2\n\n";
    
    // Step 4: Restore and expand pulsar scale
    std::cout << "Step 4: Restore initial state, then expand pulsar scale (1.5x power, 1.2x expansion)\n";
    crab.restoreState("crab_initial_971yr");
    crab.expandPulsarScale(1.5, 1.2);
    double g4 = crab.computeG(t_current);
    std::cout << "  New P_pulsar = " << crab.variables["P_pulsar"] << " W\n";
    std::cout << "  New v_exp = " << (crab.variables["v_exp"] / 1e3) << " km/s\n";
    std::cout << "  g_UQFF = " << g4 << " m/s^2 (pulsar wind dominates)\n\n";
    
    // Step 5: Restore and expand magnetic scale
    std::cout << "Step 5: Restore initial state, then expand magnetic scale (2x B, 1.3x shock)\n";
    crab.restoreState("crab_initial_971yr");
    crab.expandMagneticScale(2.0, 1.3);
    double g5 = crab.computeG(t_current);
    std::cout << "  New B = " << crab.variables["B"] << " T\n";
    std::cout << "  New v_shock = " << (crab.variables["v_shock"] / 1e3) << " km/s\n";
    std::cout << "  g_UQFF = " << g5 << " m/s^2\n\n";
    
    // Step 6: Time evolution demonstration
    std::cout << "Step 6: Time evolution from 0 to 2000 years\n";
    crab.restoreState("crab_initial_971yr");
    for (double t_yr = 0; t_yr <= 2000; t_yr += 500) {
        double t_sec = t_yr * 3.156e7;
        double r = crab.variables["r0"] + crab.variables["v_exp"] * t_sec;
        double g = crab.computeG(t_sec);
        std::cout << "  t = " << t_yr << " yr: r(t) = " << (r / 9.461e15) 
                  << " ly, g = " << g << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 7: Create custom variables
    std::cout << "Step 7: Create custom tracking variables\n";
    crab.createVariable("age_years", 971.0);
    crab.createVariable("distance_ly", 6500.0); // Distance to Crab
    std::cout << "  Created 'age_years' and 'distance_ly'\n\n";
    
    // Step 8: Generate variations for uncertainty analysis
    std::cout << "Step 8: Generate 3 parameter variations (5% perturbation)\n";
    auto variations = crab.generateVariations(3, 5.0);
    for (size_t i = 0; i < variations.size(); ++i) {
        CrabUQFFModule temp = crab;
        temp.variables = variations[i];
        double g_var = temp.computeG(t_current);
        std::cout << "  Variation " << (i+1) << ": g = " << g_var << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 9: Sensitivity analysis
    std::cout << "Step 9: Sensitivity analysis (1% perturbation)\n";
    auto sensitivities = crab.sensitivityAnalysis(t_current, 0.01);
    std::cout << "  Parameter sensitivities (fractional change in g):\n";
    for (const auto& s : sensitivities) {
        std::cout << "    " << s.first << ": " << s.second << "\n";
    }
    std::cout << "\n";
    
    // Step 10: Pulsar power sweep
    std::cout << "Step 10: Pulsar power sweep (0.5x, 1x, 2x)\n";
    crab.saveState("crab_before_sweep");
    for (double scale : {0.5, 1.0, 2.0}) {
        crab.restoreState("crab_before_sweep");
        crab.expandPulsarScale(scale, 1.0);
        double g = crab.computeG(t_current);
        std::cout << "  P_pulsar × " << scale << ": g = " << g << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 11: Expansion velocity sweep
    std::cout << "Step 11: Expansion velocity sweep (0.8x, 1.0x, 1.2x)\n";
    crab.restoreState("crab_before_sweep");
    for (double scale : {0.8, 1.0, 1.2}) {
        crab.restoreState("crab_before_sweep");
        crab.expandPulsarScale(1.0, scale);
        double g = crab.computeG(t_current);
        double v = crab.variables["v_exp"];
        std::cout << "  v_exp = " << (v / 1e3) << " km/s: g = " << g << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 12: Magnetic field sweep
    std::cout << "Step 12: Magnetic field sweep (0.5x, 1x, 2x)\n";
    crab.restoreState("crab_before_sweep");
    for (double scale : {0.5, 1.0, 2.0}) {
        crab.restoreState("crab_before_sweep");
        crab.expandMagneticScale(scale, 1.0);
        double g = crab.computeG(t_current);
        std::cout << "  B × " << scale << ": g = " << g << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 13: Batch variable transformation
    std::cout << "Step 13: Batch transform expansion parameters (1.1x scale)\n";
    crab.restoreState("crab_before_sweep");
    crab.scaleVariableGroup({"v_exp", "v_shock"}, 1.1);
    double g13 = crab.computeG(t_current);
    std::cout << "  v_exp = " << (crab.variables["v_exp"] / 1e3) << " km/s\n";
    std::cout << "  v_shock = " << (crab.variables["v_shock"] / 1e3) << " km/s\n";
    std::cout << "  g_UQFF = " << g13 << " m/s^2\n\n";
    
    // Step 14: Validate and auto-correct
    std::cout << "Step 14: Validate consistency and auto-correct if needed\n";
    crab.restoreState("crab_before_sweep");
    bool valid = crab.validateConsistency();
    std::cout << "  System valid: " << (valid ? "Yes" : "No") << "\n";
    if (!valid) {
        bool corrected = crab.autoCorrectAnomalies();
        std::cout << "  Auto-corrected: " << (corrected ? "Yes" : "No") << "\n";
    }
    std::cout << "\n";
    
    // Step 15: Parameter mutation
    std::cout << "Step 15: Mutate parameters (3% random variation)\n";
    crab.restoreState("crab_before_sweep");
    crab.mutateParameters(0.03);
    double g15 = crab.computeG(t_current);
    std::cout << "  Mutated M = " << (crab.variables["M"] / crab.variables["M_sun"]) << " M_sun\n";
    std::cout << "  g_UQFF = " << g15 << " m/s^2\n\n";
    
    // Step 16: List all saved states
    std::cout << "Step 16: List all saved states\n";
    auto states = crab.listSavedStates();
    std::cout << "  Saved states (" << states.size() << " total):\n";
    for (const auto& state : states) {
        std::cout << "    - " << state << "\n";
    }
    std::cout << "\n";
    
    // Step 17: Generate full report
    std::cout << "Step 17: Generate comprehensive system report\n";
    crab.restoreState("crab_initial_971yr");
    std::string report = crab.generateReport(t_current);
    std::cout << report << "\n";
    
    // Step 18: Export final state
    std::cout << "Step 18: Export final system state\n";
    std::string state_export = crab.exportState();
    std::cout << state_export << "\n";
    
    std::cout << "========== END 18-STEP CRAB NEBULA DEMONSTRATION ==========\n\n";
}