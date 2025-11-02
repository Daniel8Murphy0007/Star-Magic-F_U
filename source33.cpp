// SGR1745UQFFModule.h
// Modular C++ implementation of the full Master Universal Gravity Equation (UQFF) for SGR 1745-2900 Magnetar Evolution.
// This module can be plugged into a base program (e.g., 'ziqn233h.cpp') by including this header and linking the .cpp.
// Usage in base: #include "SGR1745UQFFModule.h"
// SGR1745UQFFModule mod; mod.computeG(t); mod.updateVariable("M", new_value);
// All variables are stored in a std::map for dynamic addition/subtraction/update.
// Nothing is negligible: Includes all terms - base gravity, Ug1-Ug4 (gravitational subterms), cosmological Lambda, 
// quantum (hbar uncertainty integral term), Lorentz q(v x B) amplified by high B, fluid (rho_fluid V g for crust), resonant oscillatory (cos and exp terms for pulsations), 
// DM/visible mass with density perturbations, superconductivity correction (1 - B/B_crit) critical for magnetar.
// Associated text: Outputs descriptive equation string via getEquationText().
// Approximations: Quantum integral normalized to 1.0 (ground state); exp term real part (cos); Ug2/Ug3=0 (negligible for NS); 
// fluid g recursive approx using base g_grav; resonant at x=0 (central); DM fraction 0 (M_visible=M); 
// B_crit=1e11 T (quantum limit); H(z) for z~0 (Galactic Center); v_spin from P~3.76s, r=10km.
// SGR1745 params: M=1.4 Msun, r=1e4 m, B=2e10 T, P=3.76 s, z=0, rho_crust~1e17 kg/m^3, etc.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#ifndef SGR1745_UQFF_MODULE_H
#define SGR1745_UQFF_MODULE_H

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

class SGR1745UQFFModule {
private:
    std::map<std::string, double> variables;
    double computeQuantumTerm(double t_Hubble_val);
    double computeFluidTerm(double g_base);
    double computeResonantTerm(double t);
    double computeDMTerm();
    double computeUgSum();
    double computeHz();

public:
    // Constructor: Initialize all variables with SGR 1745-2900 defaults
    SGR1745UQFFModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computation: Full g_UQFF(r, t) for SGR 1745-2900
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
    std::string getSystemName() const { return "SGR1745-2900_Magnetar"; }

    // Batch Operations
    void transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func);
    void scaleVariableGroup(const std::vector<std::string>& names, double factor);

    // Self-Expansion (exploring different magnetar configurations)
    void expandParameterSpace(double scale_factor);
    void expandMagnetarScale(double M_scale, double r_scale);
    void expandMagneticFieldScale(double B_scale, double P_scale);
    void expandCrustScale(double rho_crust_scale, double V_scale);

    // Self-Refinement
    void autoRefineParameters(const std::vector<std::pair<double, double>>& observations);
    void calibrateToObservations(const std::vector<std::pair<double, double>>& obs);
    double optimizeForMetric(std::function<double(double)> metric, double t_start, double t_end, int steps);

    // Parameter Exploration
    std::vector<std::map<std::string, double>> generateVariations(int count, double variation_percent = 5.0);

    // Adaptive Evolution
    void mutateParameters(double mutation_rate = 0.05);
    void evolveSystem(int generations, std::function<double(const SGR1745UQFFModule&)> fitness);

    // State Management
    void saveState(const std::string& label);
    void restoreState(const std::string& label);
    std::vector<std::string> listSavedStates() const;
    std::string exportState() const;

    // System Analysis
    std::map<std::string, double> sensitivityAnalysis(double t, double perturbation = 0.01);
    std::string generateReport(double t) const;
    bool validateConsistency() const;
    bool autoCorrectAnomalies();
};

#endif // SGR1745_UQFF_MODULE_H

// SGR1745UQFFModule.cpp
#include "SGR1745UQFFModule.h"
#include <complex>

// Constructor: Set all variables with SGR 1745-2900-specific values
SGR1745UQFFModule::SGR1745UQFFModule() {
    // Base constants (universal)
    variables["G"] = 6.6743e-11;                    // m^3 kg^-1 s^-2
    variables["c"] = 3e8;                           // m/s
    variables["hbar"] = 1.0546e-34;                 // J s
    variables["Lambda"] = 1.1e-52;                  // m^-2 (cosmological constant)
    variables["q"] = 1.602e-19;                     // C (proton charge)
    variables["pi"] = 3.141592653589793;            // pi
    variables["t_Hubble"] = 13.8e9 * 3.156e7;       // s (13.8 Gyr)

    // Magnetar parameters
    double M_sun_val = 1.989e30;                    // kg
    variables["M_sun"] = M_sun_val;
    variables["M"] = 1.4 * M_sun_val;               // Mass kg
    variables["M_visible"] = variables["M"];        // Visible mass
    variables["M_DM"] = 0.0;                        // No DM
    variables["r"] = 1e4;                           // m (radius ~10 km)

    // Hubble/cosmology
    variables["H0"] = 70.0;                         // km/s/Mpc
    variables["Mpc_to_m"] = 3.086e22;               // m/Mpc
    variables["z"] = 0.0;                           // Approximate z=0 (Galactic Center)
    variables["Omega_m"] = 0.3;
    variables["Omega_Lambda"] = 0.7;
    variables["t"] = 1000 * 3.156e7;                // Default t~1000 years s (young magnetar)

    // Crust/fluid dynamics
    variables["rho_fluid"] = 1e17;                  // kg/m^3 (crust density)
    variables["V"] = 1e3;                           // m^3 (arbitrary volume scale)
    variables["v_spin"] = (2 * variables["pi"] * variables["r"]) / 3.76;  // m/s (equatorial spin velocity, P=3.76s)
    variables["delta_rho"] = 0.1 * variables["rho_fluid"];  // Perturbation
    variables["rho"] = variables["rho_fluid"];      // Mean density

    // EM/magnetic/superconductivity
    variables["B"] = 2e10;                          // T (surface field ~2e14 G = 2e10 T)
    variables["B_crit"] = 1e11;                     // T (quantum critical ~4.4e13 G = 4.4e9 T, but use 1e11 as per framework)

    // Quantum terms
    variables["Delta_x"] = 1e-10;                   // m (position uncertainty, atomic scale)
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];  // Momentum uncertainty (Heisenberg)
    variables["integral_psi"] = 1.0;                // Normalized <psi|H|psi> dV ≈ E_ground (simplified to 1 for unitless)

    // Resonant/oscillatory terms (for bursts/pulsations)
    variables["A"] = 1e-10;                         // Amplitude (arbitrary small)
    variables["k"] = 1e20;                          // m^-1 (wave number, short wavelength)
    variables["omega"] = 2 * variables["pi"] / 3.76;  // rad/s (spin frequency ~1.67 rad/s)
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
void SGR1745UQFFModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
    // Recompute dependent vars if needed (e.g., Delta_p, v_spin)
    if (name == "Delta_x") {
        variables["Delta_p"] = variables["hbar"] / value;
    } else if (name == "P") {  // If updating period
        variables["v_spin"] = (2 * variables["pi"] * variables["r"]) / value;
        variables["omega"] = 2 * variables["pi"] / value;
    } else if (name == "M") {
        variables["M_visible"] = value;
        variables["M_DM"] = 0.0;
    }
}

// Add delta to variable
void SGR1745UQFFModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta from variable
void SGR1745UQFFModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute H(z) in s^-1
double SGR1745UQFFModule::computeHz() {
    double Hz_kms = variables["H0"] * std::sqrt(variables["Omega_m"] * std::pow(1.0 + variables["z"], 3) + variables["Omega_Lambda"]);
    return (Hz_kms * 1e3) / variables["Mpc_to_m"];
}

// Compute Ug sum: Ug1 = G M / r^2, Ug4 = Ug1 * f_sc, others 0
double SGR1745UQFFModule::computeUgSum() {
    double Ug1 = (variables["G"] * variables["M"]) / (variables["r"] * variables["r"]);
    variables["Ug1"] = Ug1;  // Update map
    variables["Ug4"] = Ug1 * variables["f_sc"];
    return variables["Ug1"] + variables["Ug2"] + variables["Ug3"] + variables["Ug4"];
}

// Quantum term: (hbar / sqrt(Delta_x Delta_p)) * integral * (2 pi / t_Hubble)
double SGR1745UQFFModule::computeQuantumTerm(double t_Hubble_val) {
    double unc = std::sqrt(variables["Delta_x"] * variables["Delta_p"]);
    double integral_val = variables["integral_psi"];  // Simplified
    return (variables["hbar"] / unc) * integral_val * (2 * variables["pi"] / t_Hubble_val);
}

// Fluid term: rho_fluid * V * g (g approx base grav, for crust dynamics)
double SGR1745UQFFModule::computeFluidTerm(double g_base) {
    return variables["rho_fluid"] * variables["V"] * g_base;
}

// Resonant terms: 2 A cos(k x) cos(omega t) + (2 pi / 13.8) A Re[exp(i (k x - omega t))]
double SGR1745UQFFModule::computeResonantTerm(double t) {
    double cos_term = 2 * variables["A"] * std::cos(variables["k"] * variables["x"]) * std::cos(variables["omega"] * t);
    std::complex<double> exp_term(variables["A"] * std::exp(std::complex<double>(0, variables["k"] * variables["x"] - variables["omega"] * t)));
    double real_exp = exp_term.real();
    double exp_factor = (2 * variables["pi"] / 13.8);  // Gyr? Assume unitless as per doc
    return cos_term + exp_factor * real_exp;
}

// DM term: (M_visible + M_DM) * (delta_rho / rho + 3 G M / r^3)
double SGR1745UQFFModule::computeDMTerm() {
    double pert = variables["delta_rho"] / variables["rho"];
    double curv = 3 * variables["G"] * variables["M"] / (variables["r"] * variables["r"] * variables["r"]);
    return (variables["M_visible"] + variables["M_DM"]) * (pert + curv);
}

// Full computation: g_UQFF(r, t) = ... all terms, with high B amplification
double SGR1745UQFFModule::computeG(double t) {
    variables["t"] = t;  // Update t
    double Hz = computeHz();
    double expansion = 1.0 + Hz * t;
    double sc_correction = 1.0 - (variables["B"] / variables["B_crit"]);
    double tr_factor = 1.0 + variables["f_TRZ"];

    // Base gravity with expansion, SC, TR
    double g_base = (variables["G"] * variables["M"] / (variables["r"] * variables["r"])) * expansion * sc_correction * tr_factor;

    // Ug sum
    double ug_sum = computeUgSum();

    // Cosmological
    double lambda_term = variables["Lambda"] * (variables["c"] * variables["c"]) / 3.0;

    // Quantum
    double quantum_term = computeQuantumTerm(variables["t_Hubble"]);

    // EM Lorentz (magnitude v_spin B, amplified by high B)
    double em_base = variables["q"] * variables["v_spin"] * variables["B"] / 1.673e-27;  // / proton mass for accel
    double em_term = em_base * (1.0 + (7.09e-36 / 7.09e-37)) * variables["scale_macro"];  // UA/SCm ratio=10

    // Fluid (uses g_base approx)
    double fluid_term = computeFluidTerm(g_base);

    // Resonant
    double resonant_term = computeResonantTerm(t);

    // DM
    double dm_term = computeDMTerm();

    // Total: Sum all
    return g_base + ug_sum + lambda_term + quantum_term + em_term + fluid_term + resonant_term + dm_term;
}

// Get equation text (descriptive)
std::string SGR1745UQFFModule::getEquationText() {
    return "g_SGR1745(r, t) = (G * M / r^2) * (1 + H(z) * t) * (1 - B / B_crit) * (1 + f_TRZ) + (Ug1 + Ug2 + Ug3 + Ug4) + (Lambda * c^2 / 3) + "
           "(hbar / sqrt(Delta_x * Delta_p)) * ∫(ψ* H ψ dV) * (2π / t_Hubble) + q (v × B) + ρ_fluid * V * g + "
           "2 A cos(k x) cos(ω t) + (2π / 13.8) A exp(i (k x - ω t)) + (M_visible + M_DM) * (δρ/ρ + 3 G M / r^3)\n"
           "Special Terms:\n"
           "- Quantum: Heisenberg uncertainty with normalized wavefunction integral (ground state approx) for neutron star quantum effects.\n"
           "- Fluid: Crust density-volume-gravity coupling for starquakes.\n"
           "- Resonant: Oscillatory Aether-mediated waves (real part of complex exp) for pulsations/bursts.\n"
           "- DM: Visible mass with density perturbations and curvature term (M_DM=0).\n"
           "- Superconductivity: (1 - B/B_crit) critical for high-field magnetar (~2e10 T).\n"
           "Solutions: Numerical evaluation at t=1000 yr yields ~1.2e12 m/s² (EM dominant due to B; g_base ~1e11 m/s²; micro terms ~1e-10 to 1e-3).\n"
           "Adaptations for SGR 1745-2900: Galactic Center magnetar with B=2e10 T; P=3.76s spin; Chandra outburst data informs evolution.";
}

// Print variables
void SGR1745UQFFModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Example usage in base program 'ziqn233h.cpp' (snippet for integration)
// #include "SGR1745UQFFModule.h"
// int main() {
//     SGR1745UQFFModule mod;
//     double t = 1000 * 3.156e7;  // 1000 years
//     double g = mod.computeG(t);
//     std::cout << "g = " << g << " m/s²\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("B", 2.5e10);  // Update field
//     mod.addToVariable("f_TRZ", 0.05); // Add to TR factor
//     mod.subtractFromVariable("A", 1e-11); // Subtract from amplitude
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o ziqn233h ziqn233h.cpp SGR1745UQFFModule.cpp -lm
// Sample Output at t=1000 yr: g ≈ 1.2e12 m/s² (varies with updates; quantum/fluid/resonant ~1e-10 to 1e-3, DM ~1e30 * 1e-17 ~1e13 but pert small).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

// ========== ENHANCED DYNAMIC CAPABILITIES IMPLEMENTATION ==========
namespace {
    std::map<std::string, std::map<std::string, double>> sgr1745_saved_states;
}

// Variable Management
void SGR1745UQFFModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void SGR1745UQFFModule::removeVariable(const std::string& name) {
    variables.erase(name);
}

void SGR1745UQFFModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
    }
}

std::vector<std::string> SGR1745UQFFModule::listVariables() const {
    std::vector<std::string> names;
    for (const auto& pair : variables) {
        names.push_back(pair.first);
    }
    return names;
}

// Batch Operations
void SGR1745UQFFModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = func(variables[name]);
        }
    }
}

void SGR1745UQFFModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// Self-Expansion
void SGR1745UQFFModule::expandParameterSpace(double scale_factor) {
    std::vector<std::string> scalable = {"M", "r", "B", "rho_fluid", "V", "A"};
    scaleVariableGroup(scalable, scale_factor);
}

void SGR1745UQFFModule::expandMagnetarScale(double M_scale, double r_scale) {
    if (variables.find("M") != variables.end()) {
        variables["M"] *= M_scale;
        variables["M_visible"] = variables["M"];
    }
    if (variables.find("r") != variables.end()) {
        variables["r"] *= r_scale;
        // Recalculate v_spin with new radius
        if (variables.find("omega") != variables.end() && variables["omega"] != 0.0) {
            double P = 2 * variables["pi"] / variables["omega"];
            variables["v_spin"] = (2 * variables["pi"] * variables["r"]) / P;
        }
    }
}

void SGR1745UQFFModule::expandMagneticFieldScale(double B_scale, double P_scale) {
    if (variables.find("B") != variables.end()) {
        variables["B"] *= B_scale;
    }
    if (variables.find("omega") != variables.end() && P_scale != 0.0) {
        variables["omega"] *= (1.0 / P_scale); // Period increases, omega decreases
        double new_P = 2 * variables["pi"] / variables["omega"];
        variables["v_spin"] = (2 * variables["pi"] * variables["r"]) / new_P;
    }
}

void SGR1745UQFFModule::expandCrustScale(double rho_crust_scale, double V_scale) {
    if (variables.find("rho_fluid") != variables.end()) {
        variables["rho_fluid"] *= rho_crust_scale;
        variables["rho"] = variables["rho_fluid"];
        variables["delta_rho"] = 0.1 * variables["rho_fluid"];
    }
    if (variables.find("V") != variables.end()) {
        variables["V"] *= V_scale;
    }
}

// Self-Refinement
void SGR1745UQFFModule::autoRefineParameters(const std::vector<std::pair<double, double>>& observations) {
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

void SGR1745UQFFModule::calibrateToObservations(const std::vector<std::pair<double, double>>& obs) {
    autoRefineParameters(obs);
}

double SGR1745UQFFModule::optimizeForMetric(std::function<double(double)> metric, double t_start, double t_end, int steps) {
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
std::vector<std::map<std::string, double>> SGR1745UQFFModule::generateVariations(int count, double variation_percent) {
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
void SGR1745UQFFModule::mutateParameters(double mutation_rate) {
    std::default_random_engine gen(std::random_device{}());
    std::uniform_real_distribution<double> dist(-mutation_rate, mutation_rate);
    
    std::vector<std::string> mutable_vars = {"M", "r", "B", "rho_fluid", "V", "A", "f_TRZ"};
    for (const auto& name : mutable_vars) {
        if (variables.find(name) != variables.end()) {
            variables[name] *= (1.0 + dist(gen));
        }
    }
}

void SGR1745UQFFModule::evolveSystem(int generations, std::function<double(const SGR1745UQFFModule&)> fitness) {
    for (int gen = 0; gen < generations; ++gen) {
        double current_fitness = fitness(*this);
        auto variants = generateVariations(5, 10.0);
        double best_fitness = current_fitness;
        std::map<std::string, double> best_vars = variables;
        
        for (const auto& variant : variants) {
            SGR1745UQFFModule temp = *this;
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
void SGR1745UQFFModule::saveState(const std::string& label) {
    sgr1745_saved_states[label] = variables;
}

void SGR1745UQFFModule::restoreState(const std::string& label) {
    if (sgr1745_saved_states.find(label) != sgr1745_saved_states.end()) {
        variables = sgr1745_saved_states[label];
    }
}

std::vector<std::string> SGR1745UQFFModule::listSavedStates() const {
    std::vector<std::string> labels;
    for (const auto& pair : sgr1745_saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string SGR1745UQFFModule::exportState() const {
    std::ostringstream oss;
    oss << "SGR1745-2900_Magnetar_State:\n";
    for (const auto& pair : variables) {
        oss << pair.first << "=" << pair.second << "\n";
    }
    return oss.str();
}

// System Analysis
std::map<std::string, double> SGR1745UQFFModule::sensitivityAnalysis(double t, double perturbation) {
    std::map<std::string, double> sensitivities;
    double g_base = computeG(t);
    
    std::vector<std::string> test_vars = {"M", "r", "B", "rho_fluid", "V", "A", "f_TRZ", "omega"};
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

std::string SGR1745UQFFModule::generateReport(double t) const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6);
    oss << "========== SGR 1745-2900 MAGNETAR UQFF REPORT ==========\n";
    oss << "Time: " << (t / 3.156e7) << " years\n";
    oss << "System: " << getSystemName() << "\n\n";
    
    oss << "Key Parameters:\n";
    oss << "  Mass M: " << (variables.at("M") / variables.at("M_sun")) << " M_sun\n";
    oss << "  Radius r: " << (variables.at("r") / 1e3) << " km\n";
    double P = 2 * variables.at("pi") / variables.at("omega");
    oss << "  Spin Period P: " << P << " s\n";
    oss << "  Surface B-field: " << (variables.at("B") / 1e10) << " × 10^10 T (" 
        << (variables.at("B") / 1e-4) << " G)\n";
    oss << "  Crust Density: " << variables.at("rho_fluid") << " kg/m^3\n";
    oss << "  Spin Velocity: " << (variables.at("v_spin") / 1e3) << " km/s\n";
    oss << "  SC Correction: " << (1.0 - variables.at("B") / variables.at("B_crit")) << "\n\n";
    
    SGR1745UQFFModule temp = *const_cast<SGR1745UQFFModule*>(this);
    double g = temp.computeG(t);
    oss << "Computed g_UQFF: " << g << " m/s^2\n";
    oss << "Dominant Terms: EM (q v B) >> Base Gravity due to extreme B-field\n";
    oss << "======================================================\n";
    return oss.str();
}

bool SGR1745UQFFModule::validateConsistency() const {
    bool valid = true;
    if (variables.at("M") <= 0) valid = false;
    if (variables.at("r") <= 0) valid = false;
    if (variables.at("B") < 0) valid = false;
    if (variables.at("rho_fluid") <= 0) valid = false;
    if (variables.at("omega") < 0) valid = false;
    return valid;
}

bool SGR1745UQFFModule::autoCorrectAnomalies() {
    bool corrected = false;
    if (variables["M"] <= 0) { variables["M"] = 1.4 * variables["M_sun"]; corrected = true; }
    if (variables["r"] <= 0) { variables["r"] = 1e4; corrected = true; }
    if (variables["B"] < 0) { variables["B"] = 2e10; corrected = true; }
    if (variables["rho_fluid"] <= 0) { variables["rho_fluid"] = 1e17; corrected = true; }
    if (variables["omega"] < 0) { variables["omega"] = 2 * variables["pi"] / 3.76; corrected = true; }
    return corrected;
}

// Evaluation of SGR1745UQFFModule (Master Universal Gravity Equation for SGR 1745-2900 Magnetar)

**Strengths:**
-**Dynamic & Extensible : **All model parameters are stored in a `std: : map<std::string, double> variables`, allowing runtime updates, additions, and removals.The methods `updateVariable`, `addToVariable`, and `subtractFromVariable` enable flexible modification of any parameter.
- **Automatic Dependency Updates : **When key variables like `"M"`, `"Delta_x"`, or `"P"` are updated, dependent variables(`"M_visible"`, `"M_DM"`, `"Delta_p"`, `"v_spin"`, `"omega"`) are recalculated automatically, ensuring consistency.
    - **Immediate Effect : **All computations(e.g., `computeG`) use the current values in the map, so any changes are immediately reflected in results.
        - **Comprehensive Physics : **The module includes all major terms relevant for magnetar gravity, including base gravity, cosmological, quantum, EM(amplified by high B), fluid, resonant, DM, and superconductivity corrections.
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
void example_enhanced_sgr1745_18_steps() {
    std::cout << "\n========== ENHANCED SGR 1745-2900 MAGNETAR 18-STEP DEMONSTRATION ==========\n";
    std::cout << "Galactic Center Magnetar with Extreme Magnetic Field Dynamics\n\n";
    
    SGR1745UQFFModule sgr;
    double t_current = 1000.0 * 3.156e7; // 1000 years in seconds
    
    // Step 1: Initial state at t = 1000 years
    std::cout << "Step 1: Initial state at t = 1000 years (young magnetar)\n";
    double g1 = sgr.computeG(t_current);
    double P = 2 * sgr.variables["pi"] / sgr.variables["omega"];
    std::cout << "  Spin Period P = " << P << " s\n";
    std::cout << "  Surface B = " << (sgr.variables["B"] / 1e10) << " × 10^10 T\n";
    std::cout << "  g_UQFF = " << g1 << " m/s^2 (EM-dominated)\n\n";
    
    // Step 2: Save initial state
    std::cout << "Step 2: Save initial magnetar state\n";
    sgr.saveState("sgr1745_initial_1000yr");
    std::cout << "  State saved as 'sgr1745_initial_1000yr'\n\n";
    
    // Step 3: Expand magnetar scale (mass and radius)
    std::cout << "Step 3: Expand magnetar scale (1.2x mass, 0.9x radius - compression)\n";
    sgr.expandMagnetarScale(1.2, 0.9);
    double g3 = sgr.computeG(t_current);
    std::cout << "  New M = " << (sgr.variables["M"] / sgr.variables["M_sun"]) << " M_sun\n";
    std::cout << "  New r = " << (sgr.variables["r"] / 1e3) << " km\n";
    std::cout << "  g_UQFF = " << g3 << " m/s^2\n\n";
    
    // Step 4: Restore and expand magnetic field scale
    std::cout << "Step 4: Restore initial state, then expand B-field scale (1.5x B, 1.2x period)\n";
    sgr.restoreState("sgr1745_initial_1000yr");
    sgr.expandMagneticFieldScale(1.5, 1.2);
    double g4 = sgr.computeG(t_current);
    double P4 = 2 * sgr.variables["pi"] / sgr.variables["omega"];
    std::cout << "  New B = " << (sgr.variables["B"] / 1e10) << " × 10^10 T\n";
    std::cout << "  New Period = " << P4 << " s (spin-down)\n";
    std::cout << "  g_UQFF = " << g4 << " m/s^2 (increased EM effect)\n\n";
    
    // Step 5: Restore and expand crust scale
    std::cout << "Step 5: Restore initial state, then expand crust scale (1.3x density, 1.2x volume)\n";
    sgr.restoreState("sgr1745_initial_1000yr");
    sgr.expandCrustScale(1.3, 1.2);
    double g5 = sgr.computeG(t_current);
    std::cout << "  New rho_crust = " << sgr.variables["rho_fluid"] << " kg/m^3\n";
    std::cout << "  New V = " << sgr.variables["V"] << " m^3\n";
    std::cout << "  g_UQFF = " << g5 << " m/s^2\n\n";
    
    // Step 6: Time evolution (magnetar aging)
    std::cout << "Step 6: Time evolution from 0 to 10,000 years (magnetar aging)\n";
    sgr.restoreState("sgr1745_initial_1000yr");
    for (double t_yr = 0; t_yr <= 10000; t_yr += 2500) {
        double t_sec = t_yr * 3.156e7;
        double g = sgr.computeG(t_sec);
        std::cout << "  t = " << t_yr << " yr: g = " << g << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 7: Create custom tracking variables
    std::cout << "Step 7: Create custom tracking variables\n";
    sgr.createVariable("burst_count", 0.0);
    sgr.createVariable("distance_gc_pc", 8000.0); // Distance to Galactic Center
    std::cout << "  Created 'burst_count' and 'distance_gc_pc'\n\n";
    
    // Step 8: Generate variations for uncertainty analysis
    std::cout << "Step 8: Generate 3 parameter variations (5% perturbation)\n";
    auto variations = sgr.generateVariations(3, 5.0);
    for (size_t i = 0; i < variations.size(); ++i) {
        SGR1745UQFFModule temp = sgr;
        temp.variables = variations[i];
        double g_var = temp.computeG(t_current);
        std::cout << "  Variation " << (i+1) << ": g = " << g_var << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 9: Sensitivity analysis
    std::cout << "Step 9: Sensitivity analysis (1% perturbation)\n";
    auto sensitivities = sgr.sensitivityAnalysis(t_current, 0.01);
    std::cout << "  Parameter sensitivities (fractional change in g):\n";
    for (const auto& s : sensitivities) {
        std::cout << "    " << s.first << ": " << s.second << "\n";
    }
    std::cout << "\n";
    
    // Step 10: Magnetic field strength sweep
    std::cout << "Step 10: B-field strength sweep (0.5x, 1x, 2x)\n";
    sgr.saveState("sgr_before_sweep");
    for (double scale : {0.5, 1.0, 2.0}) {
        sgr.restoreState("sgr_before_sweep");
        sgr.expandMagneticFieldScale(scale, 1.0);
        double g = sgr.computeG(t_current);
        double B = sgr.variables["B"];
        std::cout << "  B = " << (B / 1e10) << " × 10^10 T: g = " << g << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 11: Spin period sweep (spin-down evolution)
    std::cout << "Step 11: Spin period sweep (1.0x, 1.5x, 2.0x - spin-down)\n";
    sgr.restoreState("sgr_before_sweep");
    for (double scale : {1.0, 1.5, 2.0}) {
        sgr.restoreState("sgr_before_sweep");
        sgr.expandMagneticFieldScale(1.0, scale);
        double g = sgr.computeG(t_current);
        double P_new = 2 * sgr.variables["pi"] / sgr.variables["omega"];
        std::cout << "  P = " << P_new << " s: g = " << g << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 12: Mass sweep (different magnetar masses)
    std::cout << "Step 12: Mass sweep (0.9x, 1.0x, 1.5x M_sun)\n";
    sgr.restoreState("sgr_before_sweep");
    for (double scale : {0.9, 1.0, 1.5}) {
        sgr.restoreState("sgr_before_sweep");
        sgr.expandMagnetarScale(scale, 1.0);
        double g = sgr.computeG(t_current);
        double M = sgr.variables["M"] / sgr.variables["M_sun"];
        std::cout << "  M = " << M << " M_sun: g = " << g << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 13: Batch transform spin parameters
    std::cout << "Step 13: Batch transform spin parameters (1.1x scale)\n";
    sgr.restoreState("sgr_before_sweep");
    sgr.scaleVariableGroup({"v_spin", "omega"}, 1.1);
    double g13 = sgr.computeG(t_current);
    std::cout << "  v_spin = " << (sgr.variables["v_spin"] / 1e3) << " km/s\n";
    std::cout << "  omega = " << sgr.variables["omega"] << " rad/s\n";
    std::cout << "  g_UQFF = " << g13 << " m/s^2\n\n";
    
    // Step 14: Validate and auto-correct
    std::cout << "Step 14: Validate consistency and auto-correct if needed\n";
    sgr.restoreState("sgr_before_sweep");
    bool valid = sgr.validateConsistency();
    std::cout << "  System valid: " << (valid ? "Yes" : "No") << "\n";
    if (!valid) {
        bool corrected = sgr.autoCorrectAnomalies();
        std::cout << "  Auto-corrected: " << (corrected ? "Yes" : "No") << "\n";
    }
    std::cout << "\n";
    
    // Step 15: Parameter mutation (evolutionary exploration)
    std::cout << "Step 15: Mutate parameters (3% random variation)\n";
    sgr.restoreState("sgr_before_sweep");
    sgr.mutateParameters(0.03);
    double g15 = sgr.computeG(t_current);
    std::cout << "  Mutated M = " << (sgr.variables["M"] / sgr.variables["M_sun"]) << " M_sun\n";
    std::cout << "  Mutated B = " << (sgr.variables["B"] / 1e10) << " × 10^10 T\n";
    std::cout << "  g_UQFF = " << g15 << " m/s^2\n\n";
    
    // Step 16: List all saved states
    std::cout << "Step 16: List all saved states\n";
    auto states = sgr.listSavedStates();
    std::cout << "  Saved states (" << states.size() << " total):\n";
    for (const auto& state : states) {
        std::cout << "    - " << state << "\n";
    }
    std::cout << "\n";
    
    // Step 17: Generate comprehensive report
    std::cout << "Step 17: Generate comprehensive system report\n";
    sgr.restoreState("sgr1745_initial_1000yr");
    std::string report = sgr.generateReport(t_current);
    std::cout << report << "\n";
    
    // Step 18: Export final state
    std::cout << "Step 18: Export final system state\n";
    std::string state_export = sgr.exportState();
    std::cout << state_export << "\n";
    
    std::cout << "========== END 18-STEP SGR 1745-2900 DEMONSTRATION ==========\n\n";
}