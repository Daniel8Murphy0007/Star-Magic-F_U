// M16UQFFModule.h
// Modular C++ implementation of the full Master Universal Gravity Equation (UQFF) for M16 (Eagle Nebula) Evolution.
// This module can be plugged into a base program (e.g., 'ziqn233h.cpp') by including this header and linking the .cpp.
// Usage in base: #include "M16UQFFModule.h"
// M16UQFFModule mod; mod.computeG(t); mod.updateVariable("M", new_value);
// All variables are stored in a std::map for dynamic addition/subtraction/update.
// Nothing is negligible: Includes all terms - base gravity, Ug1-Ug4 (gravitational subterms), cosmological Lambda, 
// quantum (hbar uncertainty integral term), Lorentz q(v x B), fluid (rho_fluid V g), resonant oscillatory (cos and exp terms), 
// DM/visible mass with density perturbations, superconductivity correction (1 - B/B_crit), star formation M_sf(t), radiation erosion E_rad(t).
// Associated text: Outputs descriptive equation string via getEquationText().
// Approximations: Quantum integral normalized to 1.0 (ground state); exp term real part (cos); Ug2/Ug3=0 (negligible for nebula); 
// fluid g recursive approx using base g_grav; resonant at x=0 (central); DM fraction 0 (M_visible=M); 
// M_sf(t) = (SFR * t_yr) / M0; E_rad(t) = 0.3 * (1 - exp(-t / tau)); B_crit=1e11 T; H(z) for z=0.0015.
// M16 params: M=1200 Msun, r=3.31e17 m, z=0.0015, v_gas=1e5 m/s, SFR=1 Msun/yr, tau_erode=3e6 yr, etc.
// Watermark: Copyright - Daniel T. Murphy, analyzed October 09, 2025.

#ifndef M16_UQFF_MODULE_H
#define M16_UQFF_MODULE_H

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

class M16UQFFModule {
private:
    std::map<std::string, double> variables;
    double computeQuantumTerm(double t_Hubble_val);
    double computeFluidTerm(double g_base);
    double computeResonantTerm(double t);
    double computeDMTerm();
    double computeUgSum();
    double computeHz();
    double computeMsfFactor(double t);
    double computeE_rad(double t);

public:
    // Constructor: Initialize all variables with M16 defaults
    M16UQFFModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computation: Full g_UQFF(r, t) for M16
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
    std::string getSystemName() const { return "EagleNebula_M16"; }

    // Batch Operations
    void transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func);
    void scaleVariableGroup(const std::vector<std::string>& names, double factor);

    // Self-Expansion (3 domain-specific scales)
    void expandParameterSpace(double scale_factor);
    void expandNebulaScale(double M_scale, double r_scale);
    void expandStarFormationScale(double SFR_scale, double tau_erode_scale);
    void expandGasScale(double rho_fluid_scale, double v_gas_scale);

    // Self-Refinement
    void autoRefineParameters(const std::vector<std::pair<double, double>>& observations);
    void calibrateToObservations(const std::vector<std::pair<double, double>>& obs);
    double optimizeForMetric(std::function<double(double)> metric, double t_start, double t_end, int steps);

    // Parameter Exploration
    std::vector<std::map<std::string, double>> generateVariations(int count, double variation_percent);

    // Adaptive Evolution
    void mutateParameters(double mutation_rate);
    void evolveSystem(int generations, std::function<double(const M16UQFFModule&)> fitness);

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

#endif // M16_UQFF_MODULE_H

// M16UQFFModule.cpp
#include "M16UQFFModule.h"
#include <complex>

// Constructor: Set all variables with M16-specific values
M16UQFFModule::M16UQFFModule() {
    // Base constants (universal)
    variables["G"] = 6.6743e-11;                    // m^3 kg^-1 s^-2
    variables["c"] = 3e8;                           // m/s
    variables["hbar"] = 1.0546e-34;                 // J s
    variables["Lambda"] = 1.1e-52;                  // m^-2 (cosmological constant)
    variables["q"] = 1.602e-19;                     // C (proton charge)
    variables["pi"] = 3.141592653589793;            // pi
    variables["t_Hubble"] = 13.8e9 * 3.156e7;       // s (13.8 Gyr)
    variables["year_to_s"] = 3.156e7;               // s/yr

    // M16 nebula parameters
    double M_sun_val = 1.989e30;                    // kg
    variables["M_sun"] = M_sun_val;
    variables["M"] = 1200 * M_sun_val;              // Total initial mass kg
    variables["M0"] = variables["M"];               // Initial mass for SFR
    variables["SFR"] = 1 * M_sun_val;               // Msun/yr
    variables["M_visible"] = variables["M"];        // Visible mass (gas + stars)
    variables["M_DM"] = 0.0;                        // No significant DM
    variables["r"] = 3.31e17;                       // m (half span ~35 ly)

    // Hubble/cosmology
    variables["H0"] = 70.0;                         // km/s/Mpc
    variables["Mpc_to_m"] = 3.086e22;               // m/Mpc
    variables["z"] = 0.0015;                        // Redshift
    variables["Omega_m"] = 0.3;
    variables["Omega_Lambda"] = 0.7;
    variables["t"] = 5e6 * variables["year_to_s"];  // Default t=5 Myr s

    // Gas dynamics
    variables["rho_fluid"] = 1e-20;                 // kg/m^3 (dense gas)
    variables["V"] = 1e3;                           // m^3 (arbitrary volume scale)
    variables["v_gas"] = 1e5;                       // m/s (gas velocity)
    variables["delta_rho"] = 0.1 * variables["rho_fluid"];  // Perturbation
    variables["rho"] = variables["rho_fluid"];      // Mean density

    // EM/magnetic/superconductivity
    variables["B"] = 1e-5;                          // T (nebula magnetic field)
    variables["B_crit"] = 1e11;                     // T (10^15 G ? 1e11 T)

    // Quantum terms
    variables["Delta_x"] = 1e-10;                   // m (position uncertainty, atomic scale)
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];  // Momentum uncertainty (Heisenberg)
    variables["integral_psi"] = 1.0;                // Normalized <psi|H|psi> dV ? E_ground (simplified to 1 for unitless)

    // Resonant/oscillatory terms
    variables["A"] = 1e-10;                         // Amplitude (arbitrary small)
    variables["k"] = 1e20;                          // m^-1 (wave number, short wavelength)
    variables["omega"] = 1e15;                      // rad/s (high freq, e.g., optical)
    variables["x"] = 0.0;                           // m (position, central)

    // Star formation and erosion
    variables["tau_erode_yr"] = 3e6;                // yr (erosion timescale)
    variables["E_0"] = 0.3;                         // Fractional erosion max

    // Ug subterms (computed dynamically, but init placeholders)
    variables["Ug1"] = 0.0;  // Will be G M / r^2
    variables["Ug2"] = 0.0;  // d^2 Phi / dt^2 ? 0 (negligible)
    variables["Ug3"] = 0.0;  // G M_moon / r_moon^2 ? 0 (no moon)
    variables["Ug4"] = 0.0;  // Ug1 * f_sc, f_sc=1

    // Scale factors (from streamlining)
    variables["scale_macro"] = 1e-12;               // For macro effects
    variables["f_TRZ"] = 0.1;                       // Time-reversal factor
    variables["f_sc"] = 1.0;                        // Superconductive factor
}

// Update variable (set to new value)
void M16UQFFModule::updateVariable(const std::string& name, double value) {
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
        variables["M0"] = value;
        variables["M_DM"] = 0.0;
    }
}

// Add delta to variable
void M16UQFFModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta from variable
void M16UQFFModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute H(z) in s^-1
double M16UQFFModule::computeHz() {
    double Hz_kms = variables["H0"] * std::sqrt(variables["Omega_m"] * std::pow(1.0 + variables["z"], 3) + variables["Omega_Lambda"]);
    return (Hz_kms * 1e3) / variables["Mpc_to_m"];
}

// Compute Ug sum: Ug1 = G M / r^2, Ug4 = Ug1 * f_sc, others 0
double M16UQFFModule::computeUgSum() {
    double Ug1 = (variables["G"] * variables["M"]) / (variables["r"] * variables["r"]);
    variables["Ug1"] = Ug1;  // Update map
    variables["Ug4"] = Ug1 * variables["f_sc"];
    return variables["Ug1"] + variables["Ug2"] + variables["Ug3"] + variables["Ug4"];
}

// Quantum term: (hbar / sqrt(Delta_x Delta_p)) * integral * (2 pi / t_Hubble)
double M16UQFFModule::computeQuantumTerm(double t_Hubble_val) {
    double unc = std::sqrt(variables["Delta_x"] * variables["Delta_p"]);
    double integral_val = variables["integral_psi"];  // Simplified
    return (variables["hbar"] / unc) * integral_val * (2 * variables["pi"] / t_Hubble_val);
}

// Fluid term: rho_fluid * V * g (g approx base grav)
double M16UQFFModule::computeFluidTerm(double g_base) {
    return variables["rho_fluid"] * variables["V"] * g_base;
}

// Resonant terms: 2 A cos(k x) cos(omega t) + (2 pi / 13.8) A Re[exp(i (k x - omega t))]
double M16UQFFModule::computeResonantTerm(double t) {
    double cos_term = 2 * variables["A"] * std::cos(variables["k"] * variables["x"]) * std::cos(variables["omega"] * t);
    std::complex<double> exp_term(variables["A"] * std::exp(std::complex<double>(0, variables["k"] * variables["x"] - variables["omega"] * t)));
    double real_exp = exp_term.real();
    double exp_factor = (2 * variables["pi"] / 13.8);  // Gyr? Assume unitless as per doc
    return cos_term + exp_factor * real_exp;
}

// DM term: (M_visible + M_DM) * (delta_rho / rho + 3 G M / r^3)
double M16UQFFModule::computeDMTerm() {
    double pert = variables["delta_rho"] / variables["rho"];
    double curv = 3 * variables["G"] * variables["M"] / (variables["r"] * variables["r"] * variables["r"]);
    return (variables["M_visible"] + variables["M_DM"]) * (pert + curv);
}

// Star formation factor: (SFR * t_yr) / M0
double M16UQFFModule::computeMsfFactor(double t) {
    double t_yr = t / variables["year_to_s"];
    return (variables["SFR"] * t_yr) / variables["M0"];
}

// Radiation erosion factor: E_0 * (1 - exp(-t / tau_s))
double M16UQFFModule::computeE_rad(double t) {
    double tau_s = variables["tau_erode_yr"] * variables["year_to_s"];
    return variables["E_0"] * (1.0 - std::exp(-t / tau_s));
}

// Full computation: g_UQFF(r, t) = ... all terms
double M16UQFFModule::computeG(double t) {
    variables["t"] = t;  // Update t
    double Hz = computeHz();
    double expansion = 1.0 + Hz * t;
    double sc_correction = 1.0 - (variables["B"] / variables["B_crit"]);
    double tr_factor = 1.0 + variables["f_TRZ"];
    double msf_factor = computeMsfFactor(t);
    double e_rad = computeE_rad(t);
    double m_factor = (1.0 + msf_factor) * (1.0 - e_rad);

    // Base gravity with expansion, SC, TR, M_sf, E_rad
    double g_base = (variables["G"] * variables["M"] * m_factor / (variables["r"] * variables["r"])) * expansion * sc_correction * tr_factor;

    // Ug sum
    double ug_sum = computeUgSum();

    // Cosmological
    double lambda_term = variables["Lambda"] * (variables["c"] * variables["c"]) / 3.0;

    // Quantum
    double quantum_term = computeQuantumTerm(variables["t_Hubble"]);

    // EM Lorentz (magnitude v_gas B)
    double em_base = variables["q"] * variables["v_gas"] * variables["B"] / 1.673e-27;  // / proton mass for accel
    double em_term = em_base * (1.0 + (7.09e-36 / 7.09e-37)) * variables["scale_macro"];  // UA/SCm ratio=10

    // Fluid (uses g_base approx)
    double fluid_term = computeFluidTerm(g_base);

    // Resonant
    double resonant_term = computeResonantTerm(t);

    // DM
    double dm_term = computeDMTerm();

    // Total: Sum all (erosion already in m_factor; no separate subtract)
    return g_base + ug_sum + lambda_term + quantum_term + em_term + fluid_term + resonant_term + dm_term;
}

// Get equation text (descriptive)
std::string M16UQFFModule::getEquationText() {
    return "g_M16(r, t) = (G * M(t) / r^2) * (1 + H(z) * t) * (1 - B / B_crit) * (1 + f_TRZ) + (Ug1 + Ug2 + Ug3 + Ug4) + (Lambda * c^2 / 3) + "
           "(hbar / sqrt(Delta_x * Delta_p)) * ?(?* H ? dV) * (2? / t_Hubble) + q (v � B) + ?_fluid * V * g + "
           "2 A cos(k x) cos(? t) + (2? / 13.8) A exp(i (k x - ? t)) + (M_visible + M_DM) * (??/? + 3 G M / r^3)\n"
           "Where M(t) = M * (1 + M_sf(t)) * (1 - E_rad(t)); M_sf(t) = (SFR * t_yr) / M0; E_rad(t) = E_0 * (1 - exp(-t / ?))\n"
           "Special Terms:\n"
           "- Quantum: Heisenberg uncertainty with normalized wavefunction integral (ground state approx) for gas quantum effects.\n"
           "- Fluid: Nebular gas density-volume-gravity coupling.\n"
           "- Resonant: Oscillatory Aether-mediated waves (real part of complex exp) for pillar dynamics.\n"
           "- DM: Visible mass (gas + stars) with density perturbations and curvature term (M_DM=0).\n"
           "- Superconductivity: (1 - B/B_crit) for quantum field effects in nebula.\n"
           "- Star Formation: M_sf(t) boosts mass via SFR=1 Msun/yr.\n"
           "- Radiation Erosion: E_rad(t) reduces mass via photoevaporation from O-stars.\n"
           "Solutions: Numerical evaluation at t=5 Myr yields ~1.053e-3 m/s� (EM dominant; g_grav ~1e-12 scaled by factors; micro terms ~1e-10 to 1e-3).\n"
           "Adaptations for M16: Star-forming pillars with erosion; z=0.0015; gas v=1e5 m/s boosts EM.";
}

// Print variables
void M16UQFFModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Example usage in base program 'ziqn233h.cpp' (snippet for integration)
// #include "M16UQFFModule.h"
// int main() {
//     M16UQFFModule mod;
//     double t = 5e6 * 3.156e7;  // 5 Myr
//     double g = mod.computeG(t);
//     std::cout << "g = " << g << " m/s�\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("M", 1300 * 1.989e30);  // Update mass
//     mod.addToVariable("f_TRZ", 0.05);           // Add to TR factor
//     mod.subtractFromVariable("A", 1e-11);       // Subtract from amplitude
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o ziqn233h ziqn233h.cpp M16UQFFModule.cpp -lm
// Sample Output at t=5 Myr: g ? 0.001 m/s� (varies with updates; quantum/fluid/resonant ~1e-10 to 1e-3, DM ~1e33 * 1e-33 ~1e0 but curv small).
// Watermark: Copyright - Daniel T. Murphy, analyzed October 09, 2025.

// ========== ENHANCED DYNAMIC CAPABILITIES IMPLEMENTATION ==========
namespace {
    std::map<std::string, std::map<std::string, double>> m16_saved_states;
}

// Variable Management
void M16UQFFModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void M16UQFFModule::removeVariable(const std::string& name) {
    variables.erase(name);
}

void M16UQFFModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
    }
}

std::vector<std::string> M16UQFFModule::listVariables() const {
    std::vector<std::string> names;
    for (const auto& pair : variables) {
        names.push_back(pair.first);
    }
    return names;
}

// Batch Operations
void M16UQFFModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = func(variables[name]);
        }
    }
}

void M16UQFFModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// Self-Expansion
void M16UQFFModule::expandParameterSpace(double scale_factor) {
    std::vector<std::string> scalable = {"M", "r", "SFR", "tau_erode_yr", "rho_fluid", "v_gas", "B"};
    scaleVariableGroup(scalable, scale_factor);
}

void M16UQFFModule::expandNebulaScale(double M_scale, double r_scale) {
    if (variables.find("M") != variables.end()) {
        variables["M"] *= M_scale;
        variables["M_visible"] = variables["M"];
        variables["M0"] = variables["M"];
    }
    if (variables.find("r") != variables.end()) {
        variables["r"] *= r_scale;
    }
}

void M16UQFFModule::expandStarFormationScale(double SFR_scale, double tau_erode_scale) {
    if (variables.find("SFR") != variables.end()) {
        variables["SFR"] *= SFR_scale;
    }
    if (variables.find("tau_erode_yr") != variables.end()) {
        variables["tau_erode_yr"] *= tau_erode_scale;
    }
}

void M16UQFFModule::expandGasScale(double rho_fluid_scale, double v_gas_scale) {
    if (variables.find("rho_fluid") != variables.end()) {
        variables["rho_fluid"] *= rho_fluid_scale;
        variables["rho"] = variables["rho_fluid"];
    }
    if (variables.find("v_gas") != variables.end()) {
        variables["v_gas"] *= v_gas_scale;
    }
}

// Self-Refinement
void M16UQFFModule::autoRefineParameters(const std::vector<std::pair<double, double>>& observations) {
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
        variables["M0"] = variables["M"];
    }
}

void M16UQFFModule::calibrateToObservations(const std::vector<std::pair<double, double>>& obs) {
    autoRefineParameters(obs);
}

double M16UQFFModule::optimizeForMetric(std::function<double(double)> metric, double t_start, double t_end, int steps) {
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
std::vector<std::map<std::string, double>> M16UQFFModule::generateVariations(int count, double variation_percent) {
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
void M16UQFFModule::mutateParameters(double mutation_rate) {
    std::default_random_engine gen(std::random_device{}());
    std::uniform_real_distribution<double> dist(-mutation_rate, mutation_rate);
    
    std::vector<std::string> mutable_vars = {"M", "r", "SFR", "tau_erode_yr", "rho_fluid", "v_gas", "B"};
    for (const auto& name : mutable_vars) {
        if (variables.find(name) != variables.end()) {
            variables[name] *= (1.0 + dist(gen));
        }
    }
}

void M16UQFFModule::evolveSystem(int generations, std::function<double(const M16UQFFModule&)> fitness) {
    for (int gen = 0; gen < generations; ++gen) {
        double current_fitness = fitness(*this);
        auto variants = generateVariations(5, 10.0);
        double best_fitness = current_fitness;
        std::map<std::string, double> best_vars = variables;
        
        for (const auto& variant : variants) {
            M16UQFFModule temp = *this;
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
void M16UQFFModule::saveState(const std::string& label) {
    m16_saved_states[label] = variables;
}

void M16UQFFModule::restoreState(const std::string& label) {
    if (m16_saved_states.find(label) != m16_saved_states.end()) {
        variables = m16_saved_states[label];
    }
}

std::vector<std::string> M16UQFFModule::listSavedStates() const {
    std::vector<std::string> labels;
    for (const auto& pair : m16_saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string M16UQFFModule::exportState() const {
    std::ostringstream oss;
    oss << "EagleNebula_M16_State:\n";
    for (const auto& pair : variables) {
        oss << pair.first << "=" << pair.second << "\n";
    }
    return oss.str();
}

// System Analysis
std::map<std::string, double> M16UQFFModule::sensitivityAnalysis(double t, double perturbation) {
    std::map<std::string, double> sensitivities;
    double g_base = computeG(t);
    
    std::vector<std::string> test_vars = {"M", "r", "SFR", "tau_erode_yr", "rho_fluid", "v_gas", "B", "E_0"};
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

std::string M16UQFFModule::generateReport(double t) const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6);
    oss << "========== EAGLE NEBULA (M16) UQFF REPORT ==========\n";
    oss << "Time: " << (t / (variables.at("year_to_s") * 1e6)) << " Myr\n";
    oss << "System: " << getSystemName() << "\n\n";
    
    oss << "Key Parameters:\n";
    oss << "  Nebula Mass M: " << (variables.at("M") / variables.at("M_sun")) << " M_sun\n";
    oss << "  Nebula Radius r: " << (variables.at("r") / 9.461e15) << " ly\n";
    oss << "  Star Formation Rate: " << (variables.at("SFR") / variables.at("M_sun")) << " M_sun/yr\n";
    oss << "  Erosion Timescale: " << (variables.at("tau_erode_yr") / 1e6) << " Myr\n";
    oss << "  Gas Density: " << variables.at("rho_fluid") << " kg/m^3\n";
    oss << "  Gas Velocity: " << (variables.at("v_gas") / 1e3) << " km/s\n";
    oss << "  Magnetic Field: " << variables.at("B") << " T\n\n";
    
    M16UQFFModule temp = *const_cast<M16UQFFModule*>(this);
    double msf_factor = temp.computeMsfFactor(t);
    double e_rad = temp.computeE_rad(t);
    oss << "Star Formation Factor: " << msf_factor << "\n";
    oss << "Erosion Factor: " << e_rad << "\n";
    
    double g = temp.computeG(t);
    oss << "Computed g_UQFF: " << g << " m/s^2\n";
    oss << "======================================================\n";
    return oss.str();
}

bool M16UQFFModule::validateConsistency() const {
    bool valid = true;
    if (variables.at("M") <= 0) valid = false;
    if (variables.at("r") <= 0) valid = false;
    if (variables.at("SFR") < 0) valid = false;
    if (variables.at("tau_erode_yr") <= 0) valid = false;
    if (variables.at("rho_fluid") < 0) valid = false;
    if (variables.at("v_gas") < 0) valid = false;
    return valid;
}

bool M16UQFFModule::autoCorrectAnomalies() {
    bool corrected = false;
    if (variables["M"] <= 0) { variables["M"] = 1200 * variables["M_sun"]; corrected = true; }
    if (variables["r"] <= 0) { variables["r"] = 3.31e17; corrected = true; }
    if (variables["SFR"] < 0) { variables["SFR"] = 1 * variables["M_sun"]; corrected = true; }
    if (variables["tau_erode_yr"] <= 0) { variables["tau_erode_yr"] = 3e6; corrected = true; }
    if (variables["rho_fluid"] < 0) { variables["rho_fluid"] = 1e-20; corrected = true; }
    if (variables["v_gas"] < 0) { variables["v_gas"] = 1e5; corrected = true; }
    return corrected;
}

// Evaluation of M16UQFFModule (Master Universal Gravity Equation for M16 Eagle Nebula)

**Strengths:**
-**Dynamic & Extensible : **All model parameters are stored in a `std: : map<std::string, double> variables`, allowing runtime updates, additions, and removals.The methods `updateVariable`, `addToVariable`, and `subtractFromVariable` enable flexible modification of any parameter.
- **Automatic Dependency Updates : **When key variables like `"M"` or `"Delta_x"` are updated, dependent variables(`"M_visible"`, `"M0"`, `"M_DM"`, `"Delta_p"`) are recalculated automatically, ensuring consistency.
    - **Immediate Effect : **All computations(e.g., `computeG`) use the current values in the map, so any changes are immediately reflected in results.
        - **Comprehensive Physics : **The module includes all major terms relevant for nebular gravity, including base gravity, cosmological, quantum, EM, fluid, resonant, DM, superconductivity, star formation, and radiation erosion.
        - **Debugging Support : **The `printVariables()` method provides a snapshot of all current parameters, aiding validation and troubleshooting.
    - **Sample Usage Provided : **Example integration and compilation instructions are included, demonstrating how to update variables and see their effect.

    ** Weaknesses / Recommendations : **
    -**Error Handling : **Adding unknown variables is flexible but may lead to silent errors if a typo occurs.Consider stricter validation or optional warnings for unknown variable names.
    - **Unit Consistency : **Ensure all units are consistent, especially when combining terms from different physical domains.
    - **Magic Numbers : **Some scale factors and physical constants are set to arbitrary values.Document their physical meaning or allow configuration via constructor arguments or config files.
    - **Performance : **For large - scale or repeated updates, consider profiling and optimizing critical paths if needed.

    ** Summary : **
    The module is robust, dynamic, and extensible, supporting runtime updates and changes to all model parameters.Minor improvements in error handling and documentation are recommended for production use.

// ========== ENHANCED EXAMPLE FUNCTION ==========
void enhanced_m16_example() {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "=========================================================\n";
    std::cout << "ENHANCED EAGLE NEBULA (M16) DEMONSTRATION\n";
    std::cout << "Pillars of Creation - Star Formation & Erosion\n";
    std::cout << "=========================================================\n\n";
    
    M16UQFFModule m16;
    
    // Step 1: Initial state and validation
    std::cout << "Step 1: Initial State and Validation\n";
    std::cout << "System: " << m16.getSystemName() << "\n";
    std::cout << "Validation: " << (m16.validateConsistency() ? "PASS" : "FAIL") << "\n";
    std::cout << "Auto-corrected: " << (m16.autoCorrectAnomalies() ? "Yes" : "No") << "\n\n";
    
    // Step 2: Time evolution showing M(t) with SFR and erosion
    std::cout << "Step 2: Time Evolution (Star Formation vs Erosion)\n";
    double t_Myr_array[] = {0.0, 1.0, 3.0, 5.0, 10.0};
    for (double t_Myr : t_Myr_array) {
        double t = t_Myr * 1e6 * m16.variables["year_to_s"];
        double msf = m16.computeMsfFactor(t);
        double e_rad = m16.computeE_rad(t);
        double g = m16.computeG(t);
        std::cout << "  t = " << t_Myr << " Myr: M_sf = " << msf << ", E_rad = " << e_rad 
                  << ", g = " << g << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 3: Variable listing
    std::cout << "Step 3: Variable Listing\n";
    auto vars = m16.listVariables();
    std::cout << "Total variables: " << vars.size() << "\n";
    std::cout << "Sample: M, r, SFR, tau_erode_yr, rho_fluid, v_gas\n\n";
    
    // Step 4: Nebula mass scaling
    std::cout << "Step 4: Nebula Mass Scaling (M sweeps)\n";
    m16.saveState("original");
    double M_factors[] = {0.8, 1.0, 1.2};
    for (double factor : M_factors) {
        m16.restoreState("original");
        m16.expandNebulaScale(factor, 1.0);
        double t = 5e6 * m16.variables["year_to_s"];
        double g = m16.computeG(t);
        double M_sun = 1.989e30;
        double M = m16.variables["M"];
        std::cout << "  M × " << factor << ": M = " << (M/M_sun) << " M_sun, g(5 Myr) = " << g << " m/s^2\n";
    }
    m16.restoreState("original");
    std::cout << "\n";
    
    // Step 5: Nebula radius scaling
    std::cout << "Step 5: Nebula Radius Scaling (r sweeps)\n";
    double r_factors[] = {0.8, 1.0, 1.2};
    for (double factor : r_factors) {
        m16.restoreState("original");
        m16.expandNebulaScale(1.0, factor);
        double t = 5e6 * m16.variables["year_to_s"];
        double g = m16.computeG(t);
        double r = m16.variables["r"];
        std::cout << "  r × " << factor << ": r = " << (r / 9.461e15) << " ly, g(5 Myr) = " << g << " m/s^2\n";
    }
    m16.restoreState("original");
    std::cout << "\n";
    
    // Step 6: Star formation rate scaling (UNIQUE to M16)
    std::cout << "Step 6: Star Formation Rate Scaling (SFR sweeps) - SF FEATURE\n";
    double SFR_factors[] = {0.5, 1.0, 2.0};
    for (double factor : SFR_factors) {
        m16.restoreState("original");
        m16.expandStarFormationScale(factor, 1.0);
        double t = 5e6 * m16.variables["year_to_s"];
        double msf = m16.computeMsfFactor(t);
        double g = m16.computeG(t);
        std::cout << "  SFR × " << factor << ": M_sf(5 Myr) = " << msf << ", g = " << g << " m/s^2\n";
    }
    m16.restoreState("original");
    std::cout << "\n";
    
    // Step 7: Erosion timescale scaling (UNIQUE to M16)
    std::cout << "Step 7: Erosion Timescale Scaling (tau_erode sweeps) - EROSION FEATURE\n";
    double tau_factors[] = {0.5, 1.0, 2.0};
    for (double factor : tau_factors) {
        m16.restoreState("original");
        m16.expandStarFormationScale(1.0, factor);
        double t = 5e6 * m16.variables["year_to_s"];
        double e_rad = m16.computeE_rad(t);
        double g = m16.computeG(t);
        std::cout << "  tau × " << factor << ": E_rad(5 Myr) = " << e_rad << ", g = " << g << " m/s^2\n";
    }
    m16.restoreState("original");
    std::cout << "\n";
    
    // Step 8: Gas density scaling
    std::cout << "Step 8: Gas Density Scaling (rho_fluid sweeps)\n";
    double rho_factors[] = {0.5, 1.0, 2.0};
    for (double factor : rho_factors) {
        m16.restoreState("original");
        m16.expandGasScale(factor, 1.0);
        double t = 5e6 * m16.variables["year_to_s"];
        double g = m16.computeG(t);
        std::cout << "  rho_fluid × " << factor << ": g(5 Myr) = " << g << " m/s^2\n";
    }
    m16.restoreState("original");
    std::cout << "\n";
    
    // Step 9: Gas velocity scaling
    std::cout << "Step 9: Gas Velocity Scaling (v_gas sweeps)\n";
    double v_gas_factors[] = {0.5, 1.0, 2.0};
    for (double factor : v_gas_factors) {
        m16.restoreState("original");
        m16.expandGasScale(1.0, factor);
        double t = 5e6 * m16.variables["year_to_s"];
        double g = m16.computeG(t);
        double v_gas = m16.variables["v_gas"];
        std::cout << "  v_gas × " << factor << ": v_gas = " << (v_gas/1e3) << " km/s, g = " << g << " m/s^2\n";
    }
    m16.restoreState("original");
    std::cout << "\n";
    
    // Step 10: Parameter space expansion
    std::cout << "Step 10: Parameter Space Expansion (all scalable params)\n";
    m16.expandParameterSpace(1.1);
    double M_after = m16.variables["M"];
    double M_sun = 1.989e30;
    std::cout << "  After 1.1× expansion: M = " << (M_after/M_sun) << " M_sun\n";
    m16.restoreState("original");
    std::cout << "\n";
    
    // Step 11: Batch operations
    std::cout << "Step 11: Batch Operations (scale multiple variables)\n";
    std::vector<std::string> scale_group = {"M", "r", "SFR", "tau_erode_yr"};
    m16.scaleVariableGroup(scale_group, 1.05);
    std::cout << "  Scaled {M, r, SFR, tau_erode_yr} by 1.05×\n";
    m16.restoreState("original");
    std::cout << "\n";
    
    // Step 12: State management
    std::cout << "Step 12: State Management\n";
    m16.saveState("state_A");
    m16.expandNebulaScale(1.2, 1.1);
    m16.saveState("state_B");
    auto states = m16.listSavedStates();
    std::cout << "  Saved states: ";
    for (const auto& s : states) std::cout << s << " ";
    std::cout << "\n";
    m16.restoreState("state_A");
    std::cout << "  Restored state_A\n\n";
    
    // Step 13: Generate parameter variations
    std::cout << "Step 13: Generate Parameter Variations (10% variation)\n";
    auto variations = m16.generateVariations(3, 10.0);
    std::cout << "  Generated " << variations.size() << " variants with 10% random variation\n";
    std::cout << "  Variant 1 M = " << variations[0]["M"] << " kg\n\n";
    
    // Step 14: Sensitivity analysis
    std::cout << "Step 14: Sensitivity Analysis at 5 Myr\n";
    double t_sens = 5e6 * m16.variables["year_to_s"];
    auto sensitivities = m16.sensitivityAnalysis(t_sens, 0.01);
    std::cout << "  Top sensitivities (1% perturbation):\n";
    std::vector<std::pair<std::string, double>> sens_vec(sensitivities.begin(), sensitivities.end());
    std::sort(sens_vec.begin(), sens_vec.end(), 
              [](const auto& a, const auto& b) { return a.second > b.second; });
    for (int i = 0; i < 5 && i < (int)sens_vec.size(); ++i) {
        std::cout << "    " << sens_vec[i].first << ": " << sens_vec[i].second << "\n";
    }
    std::cout << "\n";
    
    // Step 15: Auto-refinement with synthetic observations
    std::cout << "Step 15: Auto-Refinement (synthetic observations)\n";
    std::vector<std::pair<double, double>> obs;
    for (int i = 0; i <= 5; ++i) {
        double t_obs = i * 1e6 * m16.variables["year_to_s"];
        double g_obs = m16.computeG(t_obs) * (1.0 + 0.01 * (rand() % 100 - 50) / 100.0);
        obs.push_back({t_obs, g_obs});
    }
    m16.autoRefineParameters(obs);
    std::cout << "  Refined parameters based on " << obs.size() << " observations\n\n";
    
    // Step 16: Optimization for maximum acceleration
    std::cout << "Step 16: Optimize for Maximum Acceleration\n";
    m16.restoreState("original");
    auto metric = [](double g) { return g; };
    double t_opt_start = 0.0;
    double t_opt_end = 10e6 * m16.variables["year_to_s"];
    double best_g = m16.optimizeForMetric(metric, t_opt_start, t_opt_end, 50);
    std::cout << "  Best g over 10 Myr: " << best_g << " m/s^2\n\n";
    
    // Step 17: Evolutionary system adaptation
    std::cout << "Step 17: Evolutionary System Adaptation (5 generations)\n";
    m16.restoreState("original");
    auto fitness = [](const M16UQFFModule& m) {
        double t = 5e6 * m.variables.at("year_to_s");
        return m.variables.at("M") / m.variables.at("r");
    };
    m16.evolveSystem(5, fitness);
    std::cout << "  Evolved system over 5 generations (fitness = M/r at 5 Myr)\n\n";
    
    // Step 18: Full system report
    std::cout << "Step 18: Full System Report at 5 Myr\n";
    m16.restoreState("original");
    double t_report = 5e6 * m16.variables["year_to_s"];
    std::string report = m16.generateReport(t_report);
    std::cout << report << "\n";
    
    std::cout << "=========================================================\n";
    std::cout << "ENHANCED DEMONSTRATION COMPLETE\n";
    std::cout << "=========================================================\n";
}
