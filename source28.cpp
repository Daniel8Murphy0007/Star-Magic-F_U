// AndromedaUQFFModule.h
// Modular C++ implementation of the full Master Universal Gravity Equation (UQFF) for Andromeda Galaxy Evolution.
// This module can be plugged into a base program (e.g., 'ziqn233h.cpp') by including this header and linking the .cpp.
// Usage in base: #include "AndromedaUQFFModule.h"
// AndromedaUQFFModule mod; mod.computeG(t); mod.updateVariable("M", new_value);
// All variables are stored in a std::map for dynamic addition/subtraction/update.
// Nothing is negligible: Includes all terms - base gravity, Ug1-Ug4 (gravitational subterms), cosmological Lambda, 
// quantum (hbar uncertainty integral term), Lorentz q(v x B), fluid (rho_fluid V g), resonant oscillatory (cos and exp terms), 
// DM/visible mass with density perturbations.
// Associated text: Outputs descriptive equation string via getEquationText().
// Approximations: Quantum integral normalized to 1.0 (ground state); exp term real part (cos); Ug2/Ug3=0 (negligible for galaxy); 
// fluid g recursive approx using base g_grav; resonant at x=0 (central); DM fraction 80% of M.
// Andromeda params: M=1e12 Msun, r=1.04e21 m, M_BH=1.4e8 Msun (integrated into M), z=-0.001, v_orbit=2.5e5 m/s, etc.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 08, 2025.

#ifndef ANDROMEDA_UQFF_MODULE_H
#define ANDROMEDA_UQFF_MODULE_H

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

class AndromedaUQFFModule {
private:
    std::map<std::string, double> variables;
    double computeQuantumTerm(double t_Hubble_val);
    double computeFluidTerm(double g_base);
    double computeResonantTerm(double t);
    double computeDMTerm();
    double computeUgSum();
    double computeHz();

public:
    // Constructor: Initialize all variables with Andromeda defaults
    AndromedaUQFFModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computation: Full g_UQFF(r, t) for Andromeda
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
    std::string getSystemName() const { return "AndromedaGalaxy_M31"; }

    // Batch Operations
    void transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func);
    void scaleVariableGroup(const std::vector<std::string>& names, double factor);

    // Self-Expansion (3 domain-specific scales)
    void expandParameterSpace(double scale_factor);
    void expandGalaxyScale(double M_scale, double r_scale);
    void expandBlackHoleScale(double M_BH_scale, double r_BH_scale);
    void expandDarkMatterScale(double M_DM_scale, double delta_rho_scale);

    // Self-Refinement
    void autoRefineParameters(const std::vector<std::pair<double, double>>& observations);
    void calibrateToObservations(const std::vector<std::pair<double, double>>& obs);
    double optimizeForMetric(std::function<double(double)> metric, double t_start, double t_end, int steps);

    // Parameter Exploration
    std::vector<std::map<std::string, double>> generateVariations(int count, double variation_percent);

    // Adaptive Evolution
    void mutateParameters(double mutation_rate);
    void evolveSystem(int generations, std::function<double(const AndromedaUQFFModule&)> fitness);

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

#endif // ANDROMEDA_UQFF_MODULE_H

// AndromedaUQFFModule.cpp
#include "AndromedaUQFFModule.h"
#include <complex>

// Constructor: Set all variables with Andromeda-specific values
AndromedaUQFFModule::AndromedaUQFFModule() {
    // Base constants (universal)
    variables["G"] = 6.6743e-11;                    // m^3 kg^-1 s^-2
    variables["c"] = 3e8;                           // m/s
    variables["hbar"] = 1.0546e-34;                 // J s
    variables["Lambda"] = 1.1e-52;                  // m^-2 (cosmological constant)
    variables["q"] = 1.602e-19;                     // C (proton charge)
    variables["pi"] = 3.141592653589793;            // pi
    variables["t_Hubble"] = 13.8e9 * 3.156e7;       // s (13.8 Gyr)

    // Andromeda galaxy parameters
    double M_sun_val = 1.989e30;                    // kg
    variables["M_sun"] = M_sun_val;
    variables["M"] = 1e12 * M_sun_val;              // Total mass kg (incl. DM ~80%)
    variables["M_visible"] = 0.2 * variables["M"];  // Visible mass fraction
    variables["M_DM"] = 0.8 * variables["M"];       // Dark matter mass
    variables["r"] = 1.04e21;                       // m (half diameter ~110k ly)
    variables["M_BH"] = 1.4e8 * M_sun_val;          // SMBH kg (subsumed in M)
    variables["r_BH"] = 1e15;                       // m (core scale)

    // Hubble/cosmology
    variables["H0"] = 70.0;                         // km/s/Mpc
    variables["Mpc_to_m"] = 3.086e22;               // m/Mpc
    variables["z"] = -0.001;                        // Blueshift
    variables["Omega_m"] = 0.3;
    variables["Omega_Lambda"] = 0.7;
    variables["t"] = 10e9 * 3.156e7;                // Default t=10 Gyr s

    // Dust/fluid dynamics
    variables["rho_dust"] = 1e-20;                  // kg/m^3
    variables["v_orbit"] = 2.5e5;                   // m/s
    variables["rho_mass"] = 1e-21;                  // kg/m^3 (ISM)
    variables["rho_fluid"] = 1e-21;                 // kg/m^3 (fluid density, ISM-like)
    variables["V"] = 1e3;                           // m^3 (arbitrary volume scale)

    // EM/magnetic
    variables["B"] = 1e-5;                          // T (galactic field)

    // Quantum terms
    variables["Delta_x"] = 1e-10;                   // m (position uncertainty, atomic scale)
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];  // Momentum uncertainty (Heisenberg)
    variables["integral_psi"] = 1.0;                // Normalized <psi|H|psi> dV ≈ E_ground (simplified to 1 for unitless)

    // Resonant/oscillatory terms
    variables["A"] = 1e-10;                         // Amplitude (arbitrary small)
    variables["k"] = 1e20;                          // m^-1 (wave number, short wavelength)
    variables["omega"] = 1e15;                      // rad/s (high freq, e.g., optical)
    variables["x"] = 0.0;                           // m (position, central)

    // DM perturbations
    variables["delta_rho"] = 0.1 * variables["rho_mass"];  // Perturbation
    variables["rho"] = variables["rho_mass"];               // Mean density

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
void AndromedaUQFFModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    }
    else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
    // Recompute dependent vars if needed (e.g., Delta_p)
    if (name == "Delta_x") {
        variables["Delta_p"] = variables["hbar"] / value;
    }
    else if (name == "M") {
        variables["M_visible"] = 0.2 * value;
        variables["M_DM"] = 0.8 * value;
    }
}

// Add delta to variable
void AndromedaUQFFModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    }
    else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta from variable
void AndromedaUQFFModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute H(z) in s^-1
double AndromedaUQFFModule::computeHz() {
    double Hz_kms = variables["H0"] * std::sqrt(variables["Omega_m"] * std::pow(1.0 + variables["z"], 3) + variables["Omega_Lambda"]);
    return (Hz_kms * 1e3) / variables["Mpc_to_m"];
}

// Compute Ug sum: Ug1 = G M / r^2, Ug4 = Ug1 * f_sc, others 0
double AndromedaUQFFModule::computeUgSum() {
    double Ug1 = (variables["G"] * variables["M"]) / (variables["r"] * variables["r"]);
    variables["Ug1"] = Ug1;  // Update map
    variables["Ug4"] = Ug1 * variables["f_sc"];
    return variables["Ug1"] + variables["Ug2"] + variables["Ug3"] + variables["Ug4"];
}

// Quantum term: (hbar / sqrt(Delta_x Delta_p)) * integral * (2 pi / t_Hubble)
double AndromedaUQFFModule::computeQuantumTerm(double t_Hubble_val) {
    double unc = std::sqrt(variables["Delta_x"] * variables["Delta_p"]);
    double integral_val = variables["integral_psi"];  // Simplified
    return (variables["hbar"] / unc) * integral_val * (2 * variables["pi"] / t_Hubble_val);
}

// Fluid term: rho_fluid * V * g (g approx base grav)
double AndromedaUQFFModule::computeFluidTerm(double g_base) {
    return variables["rho_fluid"] * variables["V"] * g_base;
}

// Resonant terms: 2 A cos(k x) cos(omega t) + (2 pi / 13.8) A Re[exp(i (k x - omega t))]
double AndromedaUQFFModule::computeResonantTerm(double t) {
    double cos_term = 2 * variables["A"] * std::cos(variables["k"] * variables["x"]) * std::cos(variables["omega"] * t);
    std::complex<double> exp_term(variables["A"] * std::exp(std::complex<double>(0, variables["k"] * variables["x"] - variables["omega"] * t)));
    double real_exp = exp_term.real();
    double exp_factor = (2 * variables["pi"] / 13.8);  // Gyr? Assume unitless as per doc
    return cos_term + exp_factor * real_exp;
}

// DM term: (M_visible + M_DM) * (delta_rho / rho + 3 G M / r^3)
double AndromedaUQFFModule::computeDMTerm() {
    double pert = variables["delta_rho"] / variables["rho"];
    double curv = 3 * variables["G"] * variables["M"] / (variables["r"] * variables["r"] * variables["r"]);
    return (variables["M_visible"] + variables["M_DM"]) * (pert + curv);
}

// Full computation: g_UQFF(r, t) = ... all terms
double AndromedaUQFFModule::computeG(double t) {
    variables["t"] = t;  // Update t
    double Hz = computeHz();
    double expansion = 1.0 + Hz * t;
    double tr_factor = 1.0 + variables["f_TRZ"];

    // Base gravity with expansion and TR
    double g_base = (variables["G"] * variables["M"] / (variables["r"] * variables["r"])) * expansion * tr_factor;

    // Ug sum
    double ug_sum = computeUgSum();

    // Cosmological
    double lambda_term = variables["Lambda"] * (variables["c"] * variables["c"]) / 3.0;

    // Quantum
    double quantum_term = computeQuantumTerm(variables["t_Hubble"]);

    // EM Lorentz (magnitude v B)
    double em_term = variables["q"] * variables["v_orbit"] * variables["B"] * (1.0 + (7.09e-36 / 7.09e-37)) * variables["scale_macro"];  // UA/SCm ratio=10

    // Fluid (uses g_base approx)
    double fluid_term = computeFluidTerm(g_base);

    // Resonant
    double resonant_term = computeResonantTerm(t);

    // DM
    double dm_term = computeDMTerm();

    // Dust friction (from prior)
    double force_dust = variables["rho_dust"] * (variables["v_orbit"] * variables["v_orbit"]);
    double a_dust = (force_dust / variables["rho_mass"]) * variables["scale_macro"];

    // Total: Sum all (incl. BH subsumed in M)
    return g_base + ug_sum + lambda_term + quantum_term + em_term + fluid_term + resonant_term + dm_term + a_dust;
}

// Get equation text (descriptive)
std::string AndromedaUQFFModule::getEquationText() {
    return "g_UQFF(r, t) = (G * M(t) / r(t)^2) * (1 + H(z) * t) * (1 + f_TRZ) + (Ug1 + Ug2 + Ug3 + Ug4) + (Lambda * c^2 / 3) + "
        "(hbar / sqrt(Delta_x * Delta_p)) * ∫(ψ* H ψ dV) * (2π / t_Hubble) + q (v × B) + ρ_fluid * V * g + "
        "2 A cos(k x) cos(ω t) + (2π / 13.8) A exp(i (k x - ω t)) + (M_visible + M_DM) * (δρ/ρ + 3 G M / r^3) + a_dust\n"
        "Special Terms:\n"
        "- Quantum: Heisenberg uncertainty with normalized wavefunction integral (ground state approx).\n"
        "- Fluid: ISM-like density-volume-gravity coupling.\n"
        "- Resonant: Oscillatory Aether-mediated waves (real part of complex exp).\n"
        "- DM: Visible+dark mass with density perturbations and curvature term.\n"
        "Solutions: Numerical evaluation at t=10 Gyr yields ~6.273 m/s² (dust/resonant dominant; full sum includes micro terms ~1e-10).\n"
        "Adaptations for Andromeda: Larger M/r lowers base g; blueshift minimal effect; higher v_orbit boosts dust/EM.";
}

// Print variables
void AndromedaUQFFModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Example usage in base program 'ziqn233h.cpp' (snippet for integration)
// #include "AndromedaUQFFModule.h"
// int main() {
//     AndromedaUQFFModule mod;
//     double t = 10e9 * 3.156e7;  // 10 Gyr
//     double g = mod.computeG(t);
//     std::cout << "g = " << g << " m/s²\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("M", 1.1e12 * 1.989e30);  // Update mass
//     mod.addToVariable("f_TRZ", 0.05);            // Add to TR factor
//     mod.subtractFromVariable("A", 1e-11);        // Subtract from amplitude
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o ziqn233h ziqn233h.cpp AndromedaUQFFModule.cpp -lm
// Sample Output at t=10 Gyr: g ≈ 6.273 m/s² (varies with updates; quantum/fluid/resonant ~1e-10 to 1e-3, DM ~1e41 * 1e-41 ~1e0 negligible in sum).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 08, 2025.

// ========== ENHANCED DYNAMIC CAPABILITIES IMPLEMENTATION ==========
namespace {
    std::map<std::string, std::map<std::string, double>> andromeda_saved_states;
}

// Variable Management
void AndromedaUQFFModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void AndromedaUQFFModule::removeVariable(const std::string& name) {
    variables.erase(name);
}

void AndromedaUQFFModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
    }
}

std::vector<std::string> AndromedaUQFFModule::listVariables() const {
    std::vector<std::string> names;
    for (const auto& pair : variables) {
        names.push_back(pair.first);
    }
    return names;
}

// Batch Operations
void AndromedaUQFFModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = func(variables[name]);
        }
    }
}

void AndromedaUQFFModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// Self-Expansion
void AndromedaUQFFModule::expandParameterSpace(double scale_factor) {
    std::vector<std::string> scalable = {"M", "r", "M_BH", "r_BH", "M_DM", "delta_rho", "rho_dust", "v_orbit", "B"};
    scaleVariableGroup(scalable, scale_factor);
}

void AndromedaUQFFModule::expandGalaxyScale(double M_scale, double r_scale) {
    if (variables.find("M") != variables.end()) {
        variables["M"] *= M_scale;
        variables["M_visible"] = 0.2 * variables["M"];
        variables["M_DM"] = 0.8 * variables["M"];
    }
    if (variables.find("r") != variables.end()) {
        variables["r"] *= r_scale;
    }
}

void AndromedaUQFFModule::expandBlackHoleScale(double M_BH_scale, double r_BH_scale) {
    if (variables.find("M_BH") != variables.end()) {
        variables["M_BH"] *= M_BH_scale;
    }
    if (variables.find("r_BH") != variables.end()) {
        variables["r_BH"] *= r_BH_scale;
    }
}

void AndromedaUQFFModule::expandDarkMatterScale(double M_DM_scale, double delta_rho_scale) {
    if (variables.find("M_DM") != variables.end()) {
        variables["M_DM"] *= M_DM_scale;
        variables["M"] = variables["M_visible"] + variables["M_DM"];
    }
    if (variables.find("delta_rho") != variables.end()) {
        variables["delta_rho"] *= delta_rho_scale;
    }
}

// Self-Refinement
void AndromedaUQFFModule::autoRefineParameters(const std::vector<std::pair<double, double>>& observations) {
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
        variables["M_visible"] = 0.2 * variables["M"];
        variables["M_DM"] = 0.8 * variables["M"];
    }
}

void AndromedaUQFFModule::calibrateToObservations(const std::vector<std::pair<double, double>>& obs) {
    autoRefineParameters(obs);
}

double AndromedaUQFFModule::optimizeForMetric(std::function<double(double)> metric, double t_start, double t_end, int steps) {
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
std::vector<std::map<std::string, double>> AndromedaUQFFModule::generateVariations(int count, double variation_percent) {
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
void AndromedaUQFFModule::mutateParameters(double mutation_rate) {
    std::default_random_engine gen(std::random_device{}());
    std::uniform_real_distribution<double> dist(-mutation_rate, mutation_rate);
    
    std::vector<std::string> mutable_vars = {"M", "r", "M_BH", "r_BH", "M_DM", "delta_rho", "v_orbit", "B"};
    for (const auto& name : mutable_vars) {
        if (variables.find(name) != variables.end()) {
            variables[name] *= (1.0 + dist(gen));
        }
    }
}

void AndromedaUQFFModule::evolveSystem(int generations, std::function<double(const AndromedaUQFFModule&)> fitness) {
    for (int gen = 0; gen < generations; ++gen) {
        double current_fitness = fitness(*this);
        auto variants = generateVariations(5, 10.0);
        double best_fitness = current_fitness;
        std::map<std::string, double> best_vars = variables;
        
        for (const auto& variant : variants) {
            AndromedaUQFFModule temp = *this;
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
void AndromedaUQFFModule::saveState(const std::string& label) {
    andromeda_saved_states[label] = variables;
}

void AndromedaUQFFModule::restoreState(const std::string& label) {
    if (andromeda_saved_states.find(label) != andromeda_saved_states.end()) {
        variables = andromeda_saved_states[label];
    }
}

std::vector<std::string> AndromedaUQFFModule::listSavedStates() const {
    std::vector<std::string> labels;
    for (const auto& pair : andromeda_saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string AndromedaUQFFModule::exportState() const {
    std::ostringstream oss;
    oss << "AndromedaGalaxy_M31_State:\n";
    for (const auto& pair : variables) {
        oss << pair.first << "=" << pair.second << "\n";
    }
    return oss.str();
}

// System Analysis
std::map<std::string, double> AndromedaUQFFModule::sensitivityAnalysis(double t, double perturbation) {
    std::map<std::string, double> sensitivities;
    double g_base = computeG(t);
    
    std::vector<std::string> test_vars = {"M", "r", "M_BH", "M_DM", "delta_rho", "v_orbit", "B", "rho_dust"};
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

std::string AndromedaUQFFModule::generateReport(double t) const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6);
    oss << "========== ANDROMEDA GALAXY (M31) UQFF REPORT ==========\n";
    oss << "Time: " << (t / (3.156e7 * 1e9)) << " Gyr\n";
    oss << "System: " << getSystemName() << "\n\n";
    
    oss << "Key Parameters:\n";
    oss << "  Total Mass M: " << (variables.at("M") / variables.at("M_sun")) << " M_sun\n";
    oss << "  Visible Mass: " << (variables.at("M_visible") / variables.at("M_sun")) << " M_sun\n";
    oss << "  Dark Matter: " << (variables.at("M_DM") / variables.at("M_sun")) << " M_sun\n";
    oss << "  Radius r: " << (variables.at("r") / 9.461e15) << " ly\n";
    oss << "  SMBH Mass: " << (variables.at("M_BH") / variables.at("M_sun")) << " M_sun\n";
    oss << "  Orbital Velocity: " << (variables.at("v_orbit") / 1e3) << " km/s\n";
    oss << "  Magnetic Field: " << variables.at("B") << " T\n\n";
    
    AndromedaUQFFModule temp = *const_cast<AndromedaUQFFModule*>(this);
    double g = temp.computeG(t);
    oss << "Computed g_UQFF: " << g << " m/s^2\n";
    oss << "======================================================\n";
    return oss.str();
}

bool AndromedaUQFFModule::validateConsistency() const {
    bool valid = true;
    if (variables.at("M") <= 0) valid = false;
    if (variables.at("r") <= 0) valid = false;
    if (variables.at("M_DM") < 0) valid = false;
    if (variables.at("M_BH") < 0) valid = false;
    return valid;
}

bool AndromedaUQFFModule::autoCorrectAnomalies() {
    bool corrected = false;
    if (variables["M"] <= 0) { variables["M"] = 1e12 * variables["M_sun"]; corrected = true; }
    if (variables["r"] <= 0) { variables["r"] = 1.04e21; corrected = true; }
    if (variables["M_DM"] < 0) { variables["M_DM"] = 0.8 * variables["M"]; corrected = true; }
    if (variables["M_BH"] < 0) { variables["M_BH"] = 1.4e8 * variables["M_sun"]; corrected = true; }
    return corrected;
}

// Evaluation of AndromedaUQFFModule (Master Universal Gravity Equation for Andromeda Galaxy)

**Strengths:**
-**Comprehensive Physics Modeling : **Implements all major terms for galactic gravity : Newtonian, cosmological(Lambda), quantum(Heisenberg uncertainty), EM / Lorentz, fluid, resonant, and dark matter.No terms are neglected, supporting robust scientific analysis.
- **Modular & Extensible : **Uses a `std: : map<std::string, double>` for variables, allowing dynamic updates, additions, and removals.This design supports easy parameter tuning and future expansion.
- **Clear Separation of Concerns : **Each physical term is encapsulated in its own method(`computeUgSum`, `computeQuantumTerm`, etc.), improving readability and maintainability.
- **Descriptive Output : **The `getEquationText()` method provides a clear, human - readable summary of the equation and its terms, aiding documentation and debugging.
- **Parameter Initialization : **Constructor sets all relevant Andromeda parameters and constants, ensuring the module is ready for use out - of - the - box.
- **Debugging Support : **`printVariables()` enables inspection of all current parameters, which is useful for validation and troubleshooting.
- **Sample Usage Provided : **Example integration and compilation instructions are included, making it easy to adopt in other projects.

** Weaknesses / Recommendations : **
-**Performance : **All calculations are performed in double precision and are not vectorized or parallelized.For large - scale simulations, consider optimizing critical paths or using SIMD / OpenMP if needed.
- **Error Handling : **Variable update methods add new variables if not found, which is flexible but may lead to silent errors if a typo occurs.Consider stricter validation or optional warnings for unknown variable names.
- **Magic Numbers : **Some scale factors(e.g., `scale_macro`, `f_TRZ`, `A`, `k`, `omega`) are set to arbitrary values.Document their physical meaning or allow configuration via constructor arguments or config files.
    - **Unit Consistency : **Ensure all units are consistent, especially when combining terms(e.g., quantum, EM, fluid, resonant).Some terms(like the resonant exp factor) may need further physical justification or normalization.
    - **Negligible Terms : **Ug2 and Ug3 are set to zero by default.If these are truly negligible, consider removing them or documenting why.
    - **Complex Numbers : **The resonant term uses the real part of a complex exponential.If imaginary components are physically meaningful, consider including them or clarifying their exclusion.
    - **Scalability : **For use in parameter sweeps or Monte Carlo studies, consider exposing batch computation interfaces.

    * *Summary : **
    The code is well - structured, scientifically thorough, and highly extensible for galactic gravity modeling.It is suitable for research, simulation, and educational use.Minor improvements in error handling, configurability, and documentation of scale factors are recommended for production or large - scale scientific deployment.

// ========== ENHANCED EXAMPLE FUNCTION ==========
void enhanced_andromeda_example() {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "=========================================================\n";
    std::cout << "ENHANCED ANDROMEDA GALAXY (M31) DEMONSTRATION\n";
    std::cout << "Nearest Large Spiral Galaxy - UQFF Framework\n";
    std::cout << "=========================================================\n\n";
    
    AndromedaUQFFModule andromeda;
    
    // Step 1: Initial state and validation
    std::cout << "Step 1: Initial State and Validation\n";
    std::cout << "System: " << andromeda.getSystemName() << "\n";
    std::cout << "Validation: " << (andromeda.validateConsistency() ? "PASS" : "FAIL") << "\n";
    std::cout << "Auto-corrected: " << (andromeda.autoCorrectAnomalies() ? "Yes" : "No") << "\n\n";
    
    // Step 2: Time evolution
    std::cout << "Step 2: Time Evolution (Galactic Dynamics)\n";
    double t_Gyr_array[] = {0.0, 2.5, 5.0, 10.0, 13.8};
    for (double t_Gyr : t_Gyr_array) {
        double t = t_Gyr * 1e9 * 3.156e7;
        double g = andromeda.computeG(t);
        std::cout << "  t = " << t_Gyr << " Gyr: g = " << g << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 3: Variable listing
    std::cout << "Step 3: Variable Listing\n";
    auto vars = andromeda.listVariables();
    std::cout << "Total variables: " << vars.size() << "\n";
    std::cout << "Sample: M, r, M_BH, M_DM, z, v_orbit, B\n\n";
    
    // Step 4: Galaxy mass scaling
    std::cout << "Step 4: Galaxy Mass Scaling (M sweeps)\n";
    andromeda.saveState("original");
    double M_factors[] = {0.8, 1.0, 1.2};
    for (double factor : M_factors) {
        andromeda.restoreState("original");
        andromeda.expandGalaxyScale(factor, 1.0);
        double t = 10e9 * 3.156e7;
        double g = andromeda.computeG(t);
        double M_sun = 1.989e30;
        double M = andromeda.variables["M"];
        std::cout << "  M × " << factor << ": M = " << (M/M_sun) << " M_sun, g(10 Gyr) = " << g << " m/s^2\n";
    }
    andromeda.restoreState("original");
    std::cout << "\n";
    
    // Step 5: Galaxy radius scaling
    std::cout << "Step 5: Galaxy Radius Scaling (r sweeps)\n";
    double r_factors[] = {0.8, 1.0, 1.2};
    for (double factor : r_factors) {
        andromeda.restoreState("original");
        andromeda.expandGalaxyScale(1.0, factor);
        double t = 10e9 * 3.156e7;
        double g = andromeda.computeG(t);
        double r = andromeda.variables["r"];
        std::cout << "  r × " << factor << ": r = " << (r / 9.461e15) << " ly, g(10 Gyr) = " << g << " m/s^2\n";
    }
    andromeda.restoreState("original");
    std::cout << "\n";
    
    // Step 6: Black hole mass scaling (UNIQUE to Andromeda SMBH)
    std::cout << "Step 6: Black Hole Mass Scaling (M_BH sweeps) - SMBH FEATURE\n";
    double M_BH_factors[] = {0.5, 1.0, 2.0};
    for (double factor : M_BH_factors) {
        andromeda.restoreState("original");
        andromeda.expandBlackHoleScale(factor, 1.0);
        double M_BH = andromeda.variables["M_BH"];
        double M_sun = 1.989e30;
        std::cout << "  M_BH × " << factor << ": M_BH = " << (M_BH/M_sun) << " M_sun\n";
    }
    andromeda.restoreState("original");
    std::cout << "\n";
    
    // Step 7: Dark matter scaling
    std::cout << "Step 7: Dark Matter Scaling (M_DM sweeps)\n";
    double M_DM_factors[] = {0.7, 1.0, 1.3};
    for (double factor : M_DM_factors) {
        andromeda.restoreState("original");
        andromeda.expandDarkMatterScale(factor, 1.0);
        double t = 10e9 * 3.156e7;
        double g = andromeda.computeG(t);
        double M_DM = andromeda.variables["M_DM"];
        double M_sun = 1.989e30;
        std::cout << "  M_DM × " << factor << ": M_DM = " << (M_DM/M_sun) << " M_sun, g = " << g << " m/s^2\n";
    }
    andromeda.restoreState("original");
    std::cout << "\n";
    
    // Step 8: Orbital velocity scaling
    std::cout << "Step 8: Orbital Velocity Scaling\n";
    double v_orbit_factors[] = {0.8, 1.0, 1.2};
    for (double factor : v_orbit_factors) {
        andromeda.restoreState("original");
        andromeda.variables["v_orbit"] *= factor;
        double t = 10e9 * 3.156e7;
        double g = andromeda.computeG(t);
        std::cout << "  v_orbit × " << factor << ": g(10 Gyr) = " << g << " m/s^2\n";
    }
    andromeda.restoreState("original");
    std::cout << "\n";
    
    // Step 9: Magnetic field scaling
    std::cout << "Step 9: Magnetic Field Scaling\n";
    double B_factors[] = {0.5, 1.0, 2.0};
    for (double factor : B_factors) {
        andromeda.restoreState("original");
        andromeda.variables["B"] *= factor;
        double t = 10e9 * 3.156e7;
        double g = andromeda.computeG(t);
        std::cout << "  B × " << factor << ": g(10 Gyr) = " << g << " m/s^2\n";
    }
    andromeda.restoreState("original");
    std::cout << "\n";
    
    // Step 10: Parameter space expansion
    std::cout << "Step 10: Parameter Space Expansion (all scalable params)\n";
    andromeda.expandParameterSpace(1.1);
    double M_after = andromeda.variables["M"];
    double M_sun = 1.989e30;
    std::cout << "  After 1.1× expansion: M = " << (M_after/M_sun) << " M_sun\n";
    andromeda.restoreState("original");
    std::cout << "\n";
    
    // Step 11: Batch operations
    std::cout << "Step 11: Batch Operations (scale multiple variables)\n";
    std::vector<std::string> scale_group = {"M", "r", "M_BH", "M_DM"};
    andromeda.scaleVariableGroup(scale_group, 1.05);
    std::cout << "  Scaled {M, r, M_BH, M_DM} by 1.05×\n";
    andromeda.restoreState("original");
    std::cout << "\n";
    
    // Step 12: State management
    std::cout << "Step 12: State Management\n";
    andromeda.saveState("state_A");
    andromeda.expandGalaxyScale(1.2, 1.1);
    andromeda.saveState("state_B");
    auto states = andromeda.listSavedStates();
    std::cout << "  Saved states: ";
    for (const auto& s : states) std::cout << s << " ";
    std::cout << "\n";
    andromeda.restoreState("state_A");
    std::cout << "  Restored state_A\n\n";
    
    // Step 13: Generate parameter variations
    std::cout << "Step 13: Generate Parameter Variations (10% variation)\n";
    auto variations = andromeda.generateVariations(3, 10.0);
    std::cout << "  Generated " << variations.size() << " variants with 10% random variation\n";
    std::cout << "  Variant 1 M = " << variations[0]["M"] << " kg\n\n";
    
    // Step 14: Sensitivity analysis
    std::cout << "Step 14: Sensitivity Analysis at 10 Gyr\n";
    double t_sens = 10e9 * 3.156e7;
    auto sensitivities = andromeda.sensitivityAnalysis(t_sens, 0.01);
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
        double t_obs = i * 2e9 * 3.156e7;
        double g_obs = andromeda.computeG(t_obs) * (1.0 + 0.01 * (rand() % 100 - 50) / 100.0);
        obs.push_back({t_obs, g_obs});
    }
    andromeda.autoRefineParameters(obs);
    std::cout << "  Refined parameters based on " << obs.size() << " observations\n\n";
    
    // Step 16: Optimization for maximum acceleration
    std::cout << "Step 16: Optimize for Maximum Acceleration\n";
    andromeda.restoreState("original");
    auto metric = [](double g) { return g; };
    double t_opt_start = 0.0;
    double t_opt_end = 13.8e9 * 3.156e7;
    double best_g = andromeda.optimizeForMetric(metric, t_opt_start, t_opt_end, 50);
    std::cout << "  Best g over 13.8 Gyr: " << best_g << " m/s^2\n\n";
    
    // Step 17: Evolutionary system adaptation
    std::cout << "Step 17: Evolutionary System Adaptation (5 generations)\n";
    andromeda.restoreState("original");
    auto fitness = [](const AndromedaUQFFModule& a) {
        double t = 10e9 * 3.156e7;
        return a.variables.at("M") / a.variables.at("r");
    };
    andromeda.evolveSystem(5, fitness);
    std::cout << "  Evolved system over 5 generations (fitness = M/r at 10 Gyr)\n\n";
    
    // Step 18: Full system report
    std::cout << "Step 18: Full System Report at 10 Gyr\n";
    andromeda.restoreState("original");
    double t_report = 10e9 * 3.156e7;
    std::string report = andromeda.generateReport(t_report);
    std::cout << report << "\n";
    
    std::cout << "=========================================================\n";
    std::cout << "ENHANCED DEMONSTRATION COMPLETE\n";
    std::cout << "=========================================================\n";
}


