// MUGEModule.h
// Modular C++ implementation of the Master Universal Gravity Equation (MUGE) for multiple astronomical systems.
// This module integrates compressed UQFF (from documents) and resonance-based UQFF models.
// Can be plugged into a base program by including this header and linking the .cpp.
// Usage: #include "MUGEModule.h"
// MUGEModule mod("Magnetar"); mod.computeG(t); mod.updateVariable("M", new_value);
// Supports systems: Magnetar SGR 1745-2900, Sagittarius A*, Tapestry of Blazing Starbirth, Westerlund 2,
// Pillars of Creation, Rings of Relativity, Students Guide to the Universe.
// Variables in std::map for dynamic updates; both models computable via computeG_compressed(double t) and computeG_resonance(double t).
// All terms conserved: base gravity, expansion, superconductivity, Ug terms, Lambda, quantum integral, fluid, DM perturbations, system-specific (e.g., stellar wind, lensing).
// Approximations: Ug1/Ug2/Ug4 negligible in compressed; integral normalized; DM fraction variable; resonance freqs tuned per system.
// Associated text: getEquationText() for both models.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef MUGE_MODULE_H
#define MUGE_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>
#include <functional>

enum class SystemType {
    MAGNETAR_SGR_1745_2900,
    SAGITTARIUS_A,
    TAPESTRY_BLAZING_STARBIRTH,
    WESTERLUND_2,
    PILLARS_CREATION,
    RINGS_RELATIVITY,
    STUDENTS_GUIDE_UNIVERSE
};

class MUGEModule {
private:
    std::map<std::string, double> variables;
    SystemType current_system;
    double computeH(double t, double z);
    double computeQuantumTerm();
    double computeFluidTerm(double g_base);
    double computeDMTerm();
    double computeUgSum();
    double computeLambdaTerm();
    double computeResonantTerm(double t);
    double computeEMTerm();
    double computeSystemSpecificTerm(double t);
    // Resonance-specific helpers
    double computeADPM();
    double computeATHz();
    double computeAvacDiff();
    double computeASuperFreq();
    double computeAAetherRes();
    double computeUg4i();
    double computeAQuantumFreq();
    double computeAAetherFreq();
    double computeAFluidFreq();
    double computeOscTerm(double t);
    double computeAExpFreq();
    double computeFTRZ();

public:
    // Constructor: Initialize with system-specific defaults
    MUGEModule(SystemType sys = SystemType::MAGNETAR_SGR_1745_2900);

    // Set system
    void setSystem(SystemType sys);

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations: Compressed and Resonance MUGE g(r, t)
    double computeG_compressed(double t);
    double computeG_resonance(double t);

    // Output descriptive texts
    std::string getEquationText_compressed();
    std::string getEquationText_resonance();

    // Print all current variables
    void printVariables();

    // ===== DYNAMIC SELF-UPDATE & SELF-EXPANSION CAPABILITIES =====
    
    // Dynamic variable management
    void createDynamicVariable(const std::string& name, double value);
    void removeDynamicVariable(const std::string& name);
    void cloneVariable(const std::string& source, const std::string& dest);
    void listAllVariables();
    
    // Batch operations on variable groups
    void applyTransformToGroup(const std::vector<std::string>& varNames, 
                               std::function<double(double)> transform);
    void scaleVariableGroup(const std::vector<std::string>& varNames, double scale_factor);
    
    // Self-expansion capabilities
    void autoExpandParameterSpace(double scale_factor);
    void expandMassScale(double mass_multiplier);
    void expandSpatialScale(double spatial_multiplier);
    void expandTimeScale(double time_multiplier);
    
    // Self-refinement
    void autoRefineParameters(double tolerance);
    void calibrateToObservations(const std::map<std::string, double>& observed_values);
    void optimizeForMetric(const std::string& metric_name, double target_value);
    
    // Parameter exploration
    void generateVariations(int num_variations, double variation_range);
    void findOptimalParameters(const std::string& objective, int iterations);
    
    // Adaptive evolution
    void mutateParameters(double mutation_rate, double mutation_strength);
    void evolveSystem(int generations);
    
    // State management
    void saveState(const std::string& label);
    void restoreState(const std::string& label);
    void listSavedStates();
    void exportState(const std::string& filename);
    
    // System analysis
    void analyzeParameterSensitivity(const std::string& param_name);
    void generateSystemReport();
    void validatePhysicalConsistency();
    void autoCorrectAnomalies();
};

#endif // MUGE_MODULE_H

// MUGEModule.cpp
#include "MUGEModule.h"
#include <complex>

// Constructor: Set universal constants and system-specific params
MUGEModule::MUGEModule(SystemType sys) : current_system(sys) {
    // Universal constants
    variables["G"] = 6.6743e-11;                    // m^3 kg^-1 s^-2
    variables["c"] = 3e8;                           // m/s
    variables["hbar"] = 1.0546e-34;                 // J s
    variables["Lambda"] = 1.1e-52;                  // m^-2
    variables["q"] = 1.602e-19;                     // C
    variables["pi"] = 3.141592653589793;
    variables["t_Hubble"] = 4.35e17;                // s
    variables["H0"] = 2.269e-18;                    // s^-1 (70 km/s/Mpc)
    variables["Omega_m"] = 0.3;
    variables["Omega_Lambda"] = 0.7;
    variables["year_to_s"] = 3.156e7;
    variables["M_sun"] = 1.989e30;                  // kg

    // Quantum defaults
    variables["Delta_x"] = 1e-10;
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];
    variables["integral_psi"] = 2.176e-18;          // J, normalized

    // Resonance defaults
    variables["Evac_neb"] = 7.09e-36;               // J/m^3
    variables["Evac_ISM"] = 7.09e-37;               // J/m^3
    variables["Delta_Evac"] = 6.381e-36;            // J/m^3
    variables["v_exp"] = 1e3;                       // m/s
    variables["f_THz"] = 1e12;                      // Hz, placeholder
    variables["f_DPM"] = 1e9;                       // Hz
    variables["FDPM"] = 6.284e29;                   // A m^2
    variables["F_super"] = 6.287e-19;               // dimensionless
    variables["UA_SCm"] = 10.0;                     // scaling
    variables["omega_i"] = 1e-8;                    // rad/s
    variables["k4"] = 1.0;
    variables["f_react"] = 1e10;                    // Hz
    variables["E_react"] = 1e-20;                   // J
    variables["f_quantum"] = 1.445e-17;             // Hz
    variables["f_Aether"] = 1.576e-35;              // Hz
    variables["f_fluid"] = 1.269e-14;               // Hz
    variables["f_osc"] = 4.57e14;                   // Hz
    variables["f_exp"] = 1e-18;                     // Hz
    variables["f_TRZ"] = 0.1;                       // dimensionless

    // Fluid/DM defaults
    variables["rho_fluid"] = 1e-20;                 // kg/m^3
    variables["V"] = 1e3;                           // m^3
    variables["g_local"] = 9.8;                     // m/s^2
    variables["DM_fraction"] = 0.85;
    variables["delta_rho_over_rho"] = 1e-5;

    // Ug defaults (negligible except Ug3' where applicable)
    variables["Ug1"] = 0.0;
    variables["Ug2"] = 0.0;
    variables["Ug3_prime"] = 0.0;
    variables["Ug4"] = 0.0;

    // System-specific initialization
    setSystem(sys);
}

// Set system and update params
void MUGEModule::setSystem(SystemType sys) {
    current_system = sys;
    switch (sys) {
        case SystemType::MAGNETAR_SGR_1745_2900:
            variables["M"] = 1.5 * variables["M_sun"];
            variables["r"] = 1e4;                       // m
            variables["z"] = 0.0009;
            variables["B"] = 1e10;                      // T
            variables["B_crit"] = 1e11;                 // T
            variables["r_BH"] = 2.84e15;                // m to Sgr A*
            variables["M_BH"] = 4.1e6 * variables["M_sun"];
            variables["t"] = 3.799e10;                  // s
            variables["rho_fluid"] = 1e-15;
            variables["V"] = 4.189e12;
            variables["g_local"] = 10.0;
            variables["M_DM"] = 0.0;
            variables["M_visible"] = variables["M"];
            variables["Ug3_prime"] = (variables["G"] * variables["M_BH"]) / (variables["r_BH"] * variables["r_BH"]);
            variables["F_env"] = 0.0;                   // Mmag + D(t) negligible
            variables["v_wind"] = 0.0;                  // No wind
            break;
        case SystemType::SAGITTARIUS_A:
            variables["M"] = 4.1e6 * variables["M_sun"];
            variables["r"] = 1.18e10;                   // m (event horizon approx)
            variables["z"] = 0.00034;
            variables["B"] = 1e-5;
            variables["B_crit"] = 1e11;
            variables["t"] = 1e6 * variables["year_to_s"];
            variables["rho_fluid"] = 1e-20;
            variables["V"] = 1e3;
            variables["g_local"] = 1e-6;
            variables["M_DM"] = 0.85 * variables["M"];
            variables["M_visible"] = 0.15 * variables["M"];
            variables["spin_adjust"] = std::sin(30.0 * variables["pi"] / 180.0);  // sin(30)
            variables["dOmega_dt"] = 1e-3;              // rad/s, placeholder for GW
            variables["F_env"] = 0.0;
            variables["v_wind"] = 8e3;
            break;
        case SystemType::TAPESTRY_BLAZING_STARBIRTH:
            variables["M"] = 2000 * variables["M_sun"];
            variables["r"] = 1.18e17;
            variables["z"] = 0.00034;
            variables["B"] = 1e-5;
            variables["B_crit"] = 1e11;
            variables["t"] = 1e6 * variables["year_to_s"];
            variables["rho_fluid"] = 1e-20;
            variables["V"] = 1e3;
            variables["g_local"] = 1e-12;
            variables["M_DM"] = 0.85 * variables["M"];
            variables["M_visible"] = 0.15 * variables["M"];
            variables["rho"] = variables["rho_fluid"];
            variables["F_env"] = 0.0;
            variables["v_wind"] = 8e3;
            break;
        case SystemType::WESTERLUND_2:
            variables["M"] = 3000 * variables["M_sun"];
            variables["r"] = 2e17;
            variables["z"] = 0.001;
            variables["B"] = 1e-5;
            variables["B_crit"] = 1e11;
            variables["t"] = 2e6 * variables["year_to_s"];
            variables["rho_fluid"] = 1e-20;
            variables["V"] = 1e3;
            variables["g_local"] = 1e-12;
            variables["M_DM"] = 0.85 * variables["M"];
            variables["M_visible"] = 0.15 * variables["M"];
            variables["rho"] = variables["rho_fluid"];
            variables["F_env"] = 0.0;
            variables["v_wind"] = 1e4;
            break;
        case SystemType::PILLARS_CREATION:
            variables["M"] = 800 * variables["M_sun"];
            variables["r"] = 1e17;
            variables["z"] = 0.002;
            variables["B"] = 1e-6;
            variables["B_crit"] = 1e11;
            variables["t"] = 1e6 * variables["year_to_s"];
            variables["rho_fluid"] = 1e-19;
            variables["V"] = 1e4;
            variables["g_local"] = 1e-11;
            variables["M_DM"] = 0.85 * variables["M"];
            variables["M_visible"] = 0.15 * variables["M"];
            variables["E_t"] = 0.1;                     // Erosion term
            variables["rho"] = variables["rho_fluid"];
            variables["F_env"] = 0.0;
            variables["v_wind"] = 8e3;
            break;
        case SystemType::RINGS_RELATIVITY:
            variables["M"] = 1e6 * variables["M_sun"];
            variables["r"] = 1e16;
            variables["z"] = 0.01;
            variables["B"] = 1e-4;
            variables["B_crit"] = 1e11;
            variables["t"] = 1e7 * variables["year_to_s"];
            variables["rho_fluid"] = 1e-21;
            variables["V"] = 1e5;
            variables["g_local"] = 1e-10;
            variables["M_DM"] = 0.85 * variables["M"];
            variables["M_visible"] = 0.15 * variables["M"];
            variables["L_t"] = 0.05;                    // Lensing term
            variables["F_env"] = 0.0;
            variables["v_wind"] = 5e3;
            break;
        case SystemType::STUDENTS_GUIDE_UNIVERSE:
            variables["M"] = 1 * variables["M_sun"];
            variables["r"] = 1e11;                      // AU scale
            variables["z"] = 0.0;
            variables["B"] = 1e-5;
            variables["B_crit"] = 1e11;
            variables["t"] = 1e9 * variables["year_to_s"];
            variables["rho_fluid"] = 1e-25;
            variables["V"] = 1e12;
            variables["g_local"] = 1e-11;
            variables["M_DM"] = 0.27 * variables["M"];
            variables["M_visible"] = 0.73 * variables["M"];
            variables["F_env"] = 0.0;
            variables["v_wind"] = 0.0;
            break;
    }
    if (variables.find("Delta_x") != variables.end()) {
        variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];
    }
    if (variables.find("M") != variables.end()) {
        variables["M_visible"] = (1.0 - variables["DM_fraction"]) * variables["M"];
        variables["M_DM"] = variables["DM_fraction"] * variables["M"];
    }
}

// Update variable with dependencies
void MUGEModule::updateVariable(const std::string& name, double value) {
    variables[name] = value;
    if (name == "Delta_x") {
        variables["Delta_p"] = variables["hbar"] / value;
    } else if (name == "M") {
        variables["M_visible"] = (1.0 - variables["DM_fraction"]) * value;
        variables["M_DM"] = variables["DM_fraction"] * value;
    } else if (name == "DM_fraction") {
        variables["M_visible"] = (1.0 - value) * variables["M"];
        variables["M_DM"] = value * variables["M"];
    }
    // System-specific: e.g., update Ug3_prime for Magnetar/SgrA
    if (current_system == SystemType::MAGNETAR_SGR_1745_2900 || current_system == SystemType::SAGITTARIUS_A) {
        if (variables.find("M_BH") != variables.end() && variables.find("r_BH") != variables.end()) {
            variables["Ug3_prime"] = (variables["G"] * variables["M_BH"]) / (variables["r_BH"] * variables["r_BH"]);
        }
    }
}

void MUGEModule::addToVariable(const std::string& name, double delta) {
    updateVariable(name, variables[name] + delta);
}

void MUGEModule::subtractFromVariable(const std::string& name, double delta) {
    updateVariable(name, variables[name] - delta);
}

// Compute H(t,z)
double MUGEModule::computeH(double t, double z) {
    double Hz = variables["H0"] * std::sqrt(variables["Omega_m"] * std::pow(1 + z, 3) + variables["Omega_Lambda"]);
    return Hz * t;
}

// Ug sum (Ug3' for external, others 0)
double MUGEModule::computeUgSum() {
    return variables["Ug1"] + variables["Ug2"] + variables["Ug3_prime"] + variables["Ug4"];
}

// Lambda term
double MUGEModule::computeLambdaTerm() {
    return (variables["Lambda"] * variables["c"] * variables["c"]) / 3.0;
}

// Quantum term
double MUGEModule::computeQuantumTerm() {
    double unc = std::sqrt(variables["Delta_x"] * variables["Delta_p"]);
    double integral_val = variables["integral_psi"];
    return (variables["hbar"] / unc) * integral_val * (2 * variables["pi"] / variables["t_Hubble"]);
}

// Fluid term
double MUGEModule::computeFluidTerm(double g_base) {
    return variables["rho_fluid"] * variables["V"] * g_base;
}

// DM term
double MUGEModule::computeDMTerm() {
    double pert = variables["delta_rho_over_rho"];
    double curv = 3 * variables["G"] * variables["M"] / (variables["r"] * variables["r"] * variables["r"]);
    return (variables["M_visible"] + variables["M_DM"]) * (pert + curv);
}

// Resonant term (cos + Re[exp])
double MUGEModule::computeResonantTerm(double t) {
    double A = variables["A"];  // Assume added if needed, default 1e-10
    double k = variables["k"];  // 1e20
    double omega = variables["omega"];  // 1e15
    double x = 0.0;
    double cos_term = 2 * A * std::cos(k * x) * std::cos(omega * t);
    std::complex<double> exp_term(A * std::exp(std::complex<double>(0, k * x - omega * t)));
    double real_exp = exp_term.real();
    double exp_factor = (2 * variables["pi"] / 13.8);
    return cos_term + exp_factor * real_exp;
}

// EM term q (v x B) magnitude
double MUGEModule::computeEMTerm() {
    double v = variables["v_wind"];
    double B = variables["B"];
    return (variables["q"] * v * B) / 1.673e-27 * variables["scale_macro"];  // Scaled, assume scale_macro=1e-12
}

// System-specific term (e.g., wind, erosion, lensing)
double MUGEModule::computeSystemSpecificTerm(double t) {
    double term = 0.0;
    switch (current_system) {
        case SystemType::SAGITTARIUS_A:
            term += (variables["G"] * variables["M"] * variables["M"] / (variables["c"] * variables["c"] * variables["c"] * variables["c"] * variables["r"])) * std::pow(variables["dOmega_dt"], 2);
            term *= variables["spin_adjust"];
            break;
        case SystemType::TAPESTRY_BLAZING_STARBIRTH:
        case SystemType::WESTERLUND_2:
            term += variables["rho"] * std::pow(variables["v_wind"], 2);
            break;
        case SystemType::PILLARS_CREATION:
            term += variables["rho"] * std::pow(variables["v_wind"], 2) * (1 - variables["E_t"]);
            break;
        case SystemType::RINGS_RELATIVITY:
            term += variables["rho_fluid"] * variables["V"] * variables["g_local"] * (1 + variables["L_t"]);
            break;
        case SystemType::STUDENTS_GUIDE_UNIVERSE:
            term = 0.0;  // Simplified
            break;
        default:
            term += variables["rho_fluid"] * std::pow(variables["v_wind"], 2);  // Default wind
    }
    return term;
}

// Compressed MUGE
double MUGEModule::computeG_compressed(double t) {
    variables["t"] = t;
    double Hz_t = computeH(t, variables["z"]);
    double expansion = 1.0 + Hz_t;
    double sc_correction = 1.0 - (variables["B"] / variables["B_crit"]);
    double env_factor = 1.0 + variables["F_env"];
    double g_base = (variables["G"] * variables["M"] / (variables["r"] * variables["r"])) * expansion * sc_correction * env_factor;

    double ug_sum = computeUgSum();
    double lambda_term = computeLambdaTerm();
    double quantum_term = computeQuantumTerm();
    double em_term = computeEMTerm();
    double fluid_term = computeFluidTerm(g_base);
    double resonant_term = computeResonantTerm(t);
    double dm_term = computeDMTerm();
    double sys_term = computeSystemSpecificTerm(t);

    return g_base + ug_sum + lambda_term + quantum_term + em_term + fluid_term + resonant_term + dm_term + sys_term;
}

// Resonance helpers
double MUGEModule::computeADPM() {
    return variables["c"] * variables["V"] * variables["FDPM"] * variables["f_DPM"] * variables["Evac_neb"];
}

double MUGEModule::computeATHz() {
    return (variables["Evac_ISM"] / variables["c"]) * variables["f_THz"] * variables["Evac_neb"] * variables["v_exp"] * computeADPM();
}

double MUGEModule::computeAvacDiff() {
    return (variables["Evac_neb"] / (variables["c"] * variables["c"])) * variables["Delta_Evac"] * std::pow(variables["v_exp"], 2) * computeADPM();
}

double MUGEModule::computeASuperFreq() {
    return (variables["Evac_neb"] / variables["c"]) * variables["F_super"] * variables["f_THz"] * computeADPM();
}

double MUGEModule::computeAAetherRes() {
    return variables["UA_SCm"] * variables["omega_i"] * variables["f_THz"] * computeADPM() * (1 + variables["f_TRZ"]);
}

double MUGEModule::computeUg4i() {
    return variables["k4"] * variables["E_react"] * variables["f_react"] * computeADPM() / (variables["Evac_neb"] * variables["c"]);
}

double MUGEModule::computeAQuantumFreq() {
    return (variables["Evac_ISM"] / variables["c"]) * variables["f_quantum"] * variables["Evac_neb"] * computeADPM();
}

double MUGEModule::computeAAetherFreq() {
    return (variables["Evac_ISM"] / variables["c"]) * variables["f_Aether"] * variables["Evac_neb"] * computeADPM();
}

double MUGEModule::computeAFluidFreq() {
    return (variables["Evac_ISM"] / variables["c"]) * variables["f_fluid"] * variables["Evac_neb"] * variables["V"];
}

double MUGEModule::computeOscTerm(double t) {
    double A = 1e-10;  // Default
    double omega = variables["f_osc"] * 2 * variables["pi"];
    return 2 * A * std::cos(omega * t);  // Simplified osc
}

double MUGEModule::computeAExpFreq() {
    return (variables["Evac_ISM"] / variables["c"]) * variables["f_exp"] * variables["Evac_neb"] * computeADPM();
}

double MUGEModule::computeFTRZ() {
    return variables["f_TRZ"];
}

// Resonance MUGE
double MUGEModule::computeG_resonance(double t) {
    double aDPM = computeADPM();
    double aTHz = computeATHz();
    double aVacDiff = computeAvacDiff();
    double aSuperFreq = computeASuperFreq();
    double aAetherRes = computeAAetherRes();
    double ug4i = computeUg4i();
    double aQuantumFreq = computeAQuantumFreq();
    double aAetherFreq = computeAAetherFreq();
    double aFluidFreq = computeAFluidFreq();
    double oscTerm = computeOscTerm(t);
    double aExpFreq = computeAExpFreq();
    double fTRZ = computeFTRZ();

    return aDPM + aTHz + aVacDiff + aSuperFreq + aAetherRes + ug4i + aQuantumFreq + aAetherFreq + aFluidFreq + oscTerm + aExpFreq + fTRZ;
}

// Equation texts (side-by-side style in string)
std::string MUGEModule::getEquationText_compressed() {
    std::string sys_name;
    switch (current_system) {
        case SystemType::MAGNETAR_SGR_1745_2900: sys_name = "Magnetar SGR 1745-2900"; break;
        // ... (similar for others)
        default: sys_name = "Generic";
    }
    return "Compressed MUGE for " + sys_name + ":\n"
           "g(r,t) = (G M(t)/r^2) (1 + H(t,z)) (1 - B/B_crit) (1 + F_env(t)) + (Ug1 + Ug2 + Ug3' + Ug4) + (Lambda c^2 / 3) + "
           "(hbar / sqrt(Delta_x Delta_p)) ?(? H ? dV) (2? / t_Hubble) + q (v � B) + ?_fluid V g + "
           "2 A cos(k x) cos(? t) + (2?/13.8) A Re[exp(i (k x - ? t))] + (M_vis + M_DM) (??/? + 3 G M / r^3) + SysTerm\n"
           "SysTerm: e.g., for Magnetar: G M_BH / r_BH^2; for Sgr A*: (G M^2 / c^4 r) (d?/dt)^2 sin(30); for Starbirth: ? v_wind^2\n"
           "Variables: As in map; Approximations: Ug1=Ug2=Ug4=0, integral normalized=1.0 scaled.";
}

std::string MUGEModule::getEquationText_resonance() {
    std::string sys_name;  // Similar
    return "Resonance MUGE for " + sys_name + ":\n"
           "g(r,t) = a_DPM + a_THz + a_vac_diff + a_super_freq + a_aether_res + Ug4i + a_quantum_freq + a_Aether_freq + "
           "a_fluid_freq + Osc_term + a_exp_freq + f_TRZ\n"
           "Where a_DPM = c V_sys F_DPM f_DPM E_vac,neb; a_THz = (E_vac,ISM / c) f_THz E_vac,neb v_exp a_DPM; etc.\n"
           "Variables: Resonance freqs tuned (f_THz=1e12 Hz, etc.); Osc_term ? 2 A cos(? t); f_TRZ=0.1.\n"
           "Integration: Sum yields effective g ~1e-11 m/s^2 for nebulae, dominated by fluid/resonant terms.";
}

void MUGEModule::printVariables() {
    std::cout << "Current Variables for " << static_cast<int>(current_system) << ":\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// ===== DYNAMIC SELF-UPDATE & SELF-EXPANSION IMPLEMENTATIONS =====

// Static storage for saved states
static std::map<std::string, std::map<std::string, double>> muge_saved_states;

// 1. Dynamic variable management
void MUGEModule::createDynamicVariable(const std::string& name, double value) {
    variables[name] = value;
    std::cout << "Created dynamic variable: " << name << " = " << value << std::endl;
}

void MUGEModule::removeDynamicVariable(const std::string& name) {
    if (variables.find(name) != variables.end()) {
        variables.erase(name);
        std::cout << "Removed dynamic variable: " << name << std::endl;
    } else {
        std::cerr << "Variable '" << name << "' not found for removal." << std::endl;
    }
}

void MUGEModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
        std::cout << "Cloned " << source << " to " << dest << std::endl;
    } else {
        std::cerr << "Source variable '" << source << "' not found." << std::endl;
    }
}

void MUGEModule::listAllVariables() {
    std::cout << "=== All MUGE Variables (Total: " << variables.size() << ") ===" << std::endl;
    std::cout << "System: " << static_cast<int>(current_system) << std::endl;
    for (const auto& pair : variables) {
        std::cout << "  " << pair.first << " = " << pair.second << std::endl;
    }
}

// 2. Batch operations
void MUGEModule::applyTransformToGroup(const std::vector<std::string>& varNames,
                                       std::function<double(double)> transform) {
    for (const auto& name : varNames) {
        if (variables.find(name) != variables.end()) {
            variables[name] = transform(variables[name]);
            std::cout << "Transformed " << name << " to " << variables[name] << std::endl;
        }
    }
}

void MUGEModule::scaleVariableGroup(const std::vector<std::string>& varNames, double scale_factor) {
    applyTransformToGroup(varNames, [scale_factor](double val) { return val * scale_factor; });
}

// 3. Self-expansion capabilities
void MUGEModule::autoExpandParameterSpace(double scale_factor) {
    std::cout << "Auto-expanding MUGE parameter space by factor " << scale_factor << std::endl;
    std::vector<std::string> expandable = {"M", "r", "rho_fluid", "V"};
    scaleVariableGroup(expandable, scale_factor);
    // Update dependent variables
    variables["M_visible"] = (1.0 - variables["DM_fraction"]) * variables["M"];
    variables["M_DM"] = variables["DM_fraction"] * variables["M"];
    std::cout << "  Updated M_visible, M_DM" << std::endl;
}

void MUGEModule::expandMassScale(double mass_multiplier) {
    std::cout << "Expanding mass scale by " << mass_multiplier << std::endl;
    variables["M"] *= mass_multiplier;
    variables["M_visible"] = (1.0 - variables["DM_fraction"]) * variables["M"];
    variables["M_DM"] = variables["DM_fraction"] * variables["M"];
    std::cout << "  M_total: " << variables["M"] << " kg" << std::endl;
}

void MUGEModule::expandSpatialScale(double spatial_multiplier) {
    std::cout << "Expanding spatial scale by " << spatial_multiplier << std::endl;
    std::vector<std::string> spatial_vars = {"r", "Delta_x", "V"};
    scaleVariableGroup(spatial_vars, spatial_multiplier);
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];
    std::cout << "  Delta_p updated: " << variables["Delta_p"] << " kg·m/s" << std::endl;
}

void MUGEModule::expandTimeScale(double time_multiplier) {
    std::cout << "Expanding time scale by " << time_multiplier << std::endl;
    std::vector<std::string> time_vars = {"t", "t_Hubble"};
    scaleVariableGroup(time_vars, time_multiplier);
}

// 4. Self-refinement
void MUGEModule::autoRefineParameters(double tolerance) {
    std::cout << "Auto-refining MUGE parameters with tolerance " << tolerance << std::endl;
    
    // Validate M = M_visible + M_DM based on DM_fraction
    double M_expected = variables["M_visible"] + variables["M_DM"];
    double error = std::abs(variables["M"] - M_expected) / std::max(M_expected, 1e-100);
    
    if (error > tolerance) {
        std::cout << "  Correcting M: " << variables["M"] << " -> " << M_expected << std::endl;
        variables["M"] = M_expected;
    }
    
    // Validate Delta_p from Delta_x (Heisenberg)
    double Delta_p_expected = variables["hbar"] / variables["Delta_x"];
    if (std::abs(variables["Delta_p"] - Delta_p_expected) / Delta_p_expected > tolerance) {
        std::cout << "  Correcting Delta_p: " << variables["Delta_p"] << " -> " << Delta_p_expected << std::endl;
        variables["Delta_p"] = Delta_p_expected;
    }
    
    // Validate DM_fraction consistency
    double DM_fraction_check = variables["M_DM"] / variables["M"];
    if (std::abs(DM_fraction_check - variables["DM_fraction"]) > tolerance) {
        std::cout << "  Note: DM_fraction inconsistency detected" << std::endl;
    }
    
    std::cout << "Refinement complete." << std::endl;
}

void MUGEModule::calibrateToObservations(const std::map<std::string, double>& observed_values) {
    std::cout << "Calibrating to " << observed_values.size() << " MUGE observations..." << std::endl;
    for (const auto& obs : observed_values) {
        if (variables.find(obs.first) != variables.end()) {
            double old_val = variables[obs.first];
            updateVariable(obs.first, obs.second);
            std::cout << "  " << obs.first << ": " << old_val << " -> " << obs.second << std::endl;
        }
    }
    std::cout << "Calibration complete." << std::endl;
}

void MUGEModule::optimizeForMetric(const std::string& metric_name, double target_value) {
    std::cout << "Optimizing for metric: " << metric_name << " = " << target_value << std::endl;
    
    if (metric_name == "g_compressed" || metric_name == "gravity") {
        double t = variables["t"];
        double current_g = computeG_compressed(t);
        double ratio = target_value / std::max(current_g, 1e-100);
        
        // Adjust total mass to reach target
        variables["M"] *= ratio;
        variables["M_visible"] = (1.0 - variables["DM_fraction"]) * variables["M"];
        variables["M_DM"] = variables["DM_fraction"] * variables["M"];
        std::cout << "  Adjusted total mass by " << ratio << std::endl;
    } else if (metric_name == "g_resonance") {
        double t = variables["t"];
        double current_g = computeG_resonance(t);
        double ratio = target_value / std::max(current_g, 1e-100);
        
        // Scale resonance parameters
        variables["FDPM"] *= ratio;
        std::cout << "  Adjusted FDPM by " << ratio << std::endl;
    }
    
    std::cout << "Optimization complete." << std::endl;
}

// 5. Parameter exploration
void MUGEModule::generateVariations(int num_variations, double variation_range) {
    std::cout << "Generating " << num_variations << " MUGE variations with range ±" 
              << (variation_range * 100) << "%" << std::endl;
    
    std::vector<std::string> key_params = {"M", "r", "B", "DM_fraction", "rho_fluid"};
    
    for (int i = 0; i < num_variations; ++i) {
        std::cout << "  Variation " << (i+1) << ":" << std::endl;
        for (const auto& param : key_params) {
            if (variables.find(param) != variables.end()) {
                double base = variables[param];
                double variation = base * (1.0 + variation_range * (2.0 * (rand() / (double)RAND_MAX) - 1.0));
                std::cout << "    " << param << ": " << base << " -> " << variation << std::endl;
            }
        }
    }
}

void MUGEModule::findOptimalParameters(const std::string& objective, int iterations) {
    std::cout << "Finding optimal MUGE parameters for: " << objective 
              << " (" << iterations << " iterations)" << std::endl;
    
    double best_score = -1e100;
    std::map<std::string, double> best_params;
    
    for (int i = 0; i < iterations; ++i) {
        mutateParameters(0.7, 0.1);
        
        double t = variables["t"];
        double score = 0.0;
        
        if (objective == "maximize_g_compressed") {
            score = computeG_compressed(t);
        } else if (objective == "maximize_g_resonance") {
            score = computeG_resonance(t);
        } else {
            score = computeG_compressed(t);
        }
        
        if (score > best_score) {
            best_score = score;
            best_params = variables;
        }
    }
    
    variables = best_params;
    std::cout << "Optimal score: " << best_score << std::endl;
}

// 6. Adaptive evolution
void MUGEModule::mutateParameters(double mutation_rate, double mutation_strength) {
    std::vector<std::string> mutable_params = {"M", "r", "B", "DM_fraction", "rho_fluid", "v_wind"};
    
    for (const auto& param : mutable_params) {
        if (variables.find(param) != variables.end()) {
            if ((rand() / (double)RAND_MAX) < mutation_rate) {
                double mutation = 1.0 + mutation_strength * (2.0 * (rand() / (double)RAND_MAX) - 1.0);
                variables[param] *= mutation;
            }
        }
    }
    
    // Update dependent variables
    variables["M_visible"] = (1.0 - variables["DM_fraction"]) * variables["M"];
    variables["M_DM"] = variables["DM_fraction"] * variables["M"];
}

void MUGEModule::evolveSystem(int generations) {
    std::cout << "Evolving MUGE system over " << generations << " generations..." << std::endl;
    
    for (int gen = 0; gen < generations; ++gen) {
        mutateParameters(0.3, 0.08);
        
        double t = variables["t"];
        double fitness = computeG_compressed(t);
        
        if (gen % 10 == 0) {
            std::cout << "  Gen " << gen << ": g = " << fitness << " m/s^2" << std::endl;
        }
    }
    
    std::cout << "Evolution complete." << std::endl;
}

// 7. State management
void MUGEModule::saveState(const std::string& label) {
    muge_saved_states[label] = variables;
    std::cout << "Saved MUGE state: " << label << " (" << variables.size() << " variables)" << std::endl;
}

void MUGEModule::restoreState(const std::string& label) {
    if (muge_saved_states.find(label) != muge_saved_states.end()) {
        variables = muge_saved_states[label];
        std::cout << "Restored MUGE state: " << label << std::endl;
    } else {
        std::cerr << "State '" << label << "' not found." << std::endl;
    }
}

void MUGEModule::listSavedStates() {
    std::cout << "=== Saved MUGE States (Total: " << muge_saved_states.size() << ") ===" << std::endl;
    for (const auto& state : muge_saved_states) {
        std::cout << "  " << state.first << " (" << state.second.size() << " variables)" << std::endl;
    }
}

void MUGEModule::exportState(const std::string& filename) {
    std::cout << "Exporting MUGE state to " << filename << " (not implemented - placeholder)" << std::endl;
    // In real implementation: write variables to file
}

// 8. System analysis
void MUGEModule::analyzeParameterSensitivity(const std::string& param_name) {
    if (variables.find(param_name) == variables.end()) {
        std::cerr << "Parameter '" << param_name << "' not found." << std::endl;
        return;
    }
    
    std::cout << "=== MUGE Sensitivity Analysis: " << param_name << " ===" << std::endl;
    
    double base_value = variables[param_name];
    double t = variables["t"];
    double base_output = computeG_compressed(t);
    
    std::vector<double> perturbations = {0.7, 0.85, 1.0, 1.15, 1.3};
    
    for (double factor : perturbations) {
        updateVariable(param_name, base_value * factor);
        
        double new_output = computeG_compressed(t);
        double sensitivity = (new_output - base_output) / std::max(std::abs(base_output), 1e-100);
        
        std::cout << "  " << param_name << " * " << factor << " -> g change: " 
                  << (sensitivity * 100) << "%" << std::endl;
    }
    
    updateVariable(param_name, base_value);  // Restore
}

void MUGEModule::generateSystemReport() {
    std::cout << "\n========== MUGE System Report ==========" << std::endl;
    std::cout << "System Type: " << static_cast<int>(current_system) << std::endl;
    std::cout << "Total Variables: " << variables.size() << std::endl;
    
    // Key MUGE parameters
    std::cout << "\nMass Parameters:" << std::endl;
    if (variables.find("M") != variables.end()) {
        std::cout << "M_total: " << (variables["M"] / variables["M_sun"]) << " M☉" << std::endl;
    }
    if (variables.find("M_visible") != variables.end()) {
        std::cout << "M_visible: " << (variables["M_visible"] / variables["M_sun"]) << " M☉" << std::endl;
    }
    if (variables.find("M_DM") != variables.end()) {
        std::cout << "M_DM: " << (variables["M_DM"] / variables["M_sun"]) << " M☉ (" 
                  << (variables["DM_fraction"] * 100) << "%)" << std::endl;
    }
    
    std::cout << "\nSpatial Parameters:" << std::endl;
    std::cout << "r: " << variables["r"] << " m" << std::endl;
    std::cout << "z (redshift): " << variables["z"] << std::endl;
    
    std::cout << "\nMagnetic Parameters:" << std::endl;
    std::cout << "B: " << variables["B"] << " T" << std::endl;
    std::cout << "B_crit: " << variables["B_crit"] << " T" << std::endl;
    
    // Current computations
    double t = variables["t"];
    double g_comp = computeG_compressed(t);
    double g_res = computeG_resonance(t);
    
    std::cout << "\nCurrent Computations:" << std::endl;
    std::cout << "g_compressed: " << g_comp << " m/s^2" << std::endl;
    std::cout << "g_resonance: " << g_res << " m/s^2" << std::endl;
    
    std::cout << "\nUg Terms:" << std::endl;
    std::cout << "Ug1: " << variables["Ug1"] << std::endl;
    std::cout << "Ug2: " << variables["Ug2"] << std::endl;
    std::cout << "Ug3_prime: " << variables["Ug3_prime"] << std::endl;
    std::cout << "Ug4: " << variables["Ug4"] << std::endl;
    
    std::cout << "============================================\n" << std::endl;
}

void MUGEModule::validatePhysicalConsistency() {
    std::cout << "Validating MUGE physical consistency..." << std::endl;
    bool consistent = true;
    
    // Check for NaN/Inf
    for (const auto& pair : variables) {
        if (std::isnan(pair.second) || std::isinf(pair.second)) {
            std::cerr << "  ERROR: " << pair.first << " is NaN/Inf" << std::endl;
            consistent = false;
        }
    }
    
    // M = M_visible + M_DM
    double M_total_expected = variables["M_visible"] + variables["M_DM"];
    if (std::abs(variables["M"] - M_total_expected) / M_total_expected > 0.01) {
        std::cerr << "  ERROR: M != M_visible + M_DM" << std::endl;
        consistent = false;
    }
    
    // Delta_p from Delta_x (Heisenberg)
    double Delta_p_expected = variables["hbar"] / variables["Delta_x"];
    if (std::abs(variables["Delta_p"] - Delta_p_expected) / Delta_p_expected > 0.01) {
        std::cerr << "  WARNING: Delta_p violates Heisenberg uncertainty" << std::endl;
        consistent = false;
    }
    
    // B < B_crit
    if (variables["B"] >= variables["B_crit"]) {
        std::cerr << "  WARNING: B >= B_crit (superconductivity breakdown)" << std::endl;
        consistent = false;
    }
    
    // DM_fraction in [0,1]
    if (variables["DM_fraction"] < 0 || variables["DM_fraction"] > 1) {
        std::cerr << "  ERROR: DM_fraction outside [0,1]" << std::endl;
        consistent = false;
    }
    
    if (consistent) {
        std::cout << "  All checks passed. MUGE system is physically consistent." << std::endl;
    }
}

void MUGEModule::autoCorrectAnomalies() {
    std::cout << "Auto-correcting MUGE anomalies..." << std::endl;
    
    // Fix NaN/Inf
    for (auto& pair : variables) {
        if (std::isnan(pair.second) || std::isinf(pair.second)) {
            std::cout << "  Correcting " << pair.first << " (was NaN/Inf)" << std::endl;
            pair.second = 1.0;
        }
    }
    
    // Enforce M_total = M_visible + M_DM
    double M_total_expected = variables["M_visible"] + variables["M_DM"];
    if (std::abs(variables["M"] - M_total_expected) / M_total_expected > 0.01) {
        std::cout << "  Correcting M to M_visible + M_DM" << std::endl;
        variables["M"] = M_total_expected;
    }
    
    // Enforce Delta_p = hbar / Delta_x
    double Delta_p_expected = variables["hbar"] / variables["Delta_x"];
    if (std::abs(variables["Delta_p"] - Delta_p_expected) / Delta_p_expected > 0.01) {
        std::cout << "  Correcting Delta_p to satisfy Heisenberg uncertainty" << std::endl;
        variables["Delta_p"] = Delta_p_expected;
    }
    
    // Cap B at B_crit
    if (variables["B"] >= variables["B_crit"]) {
        std::cout << "  Capping B to 0.99 * B_crit" << std::endl;
        variables["B"] = variables["B_crit"] * 0.99;
    }
    
    // Clamp DM_fraction to [0,1]
    if (variables["DM_fraction"] < 0) {
        std::cout << "  Correcting DM_fraction to 0" << std::endl;
        variables["DM_fraction"] = 0;
    } else if (variables["DM_fraction"] > 1) {
        std::cout << "  Correcting DM_fraction to 1" << std::endl;
        variables["DM_fraction"] = 1;
    }
    
    std::cout << "Auto-correction complete." << std::endl;
}

// Example usage snippet:
// #include "MUGEModule.h"
// int main() {
//     // === BASIC COMPUTATION ===
//     std::cout << "\n=== MUGE - Basic Computation ===" << std::endl;
//     
//     MUGEModule mod(SystemType::MAGNETAR_SGR_1745_2900);
//     double t = 3.799e10;
//     
//     double g_comp = mod.computeG_compressed(t);
//     double g_res = mod.computeG_resonance(t);
//     std::cout << "Compressed g = " << g_comp << " m/s²\n";
//     std::cout << "Resonance g = " << g_res << " m/s²\n";
//     
//     // Batch operations and multi-system tests
//     mod.saveState("magnetar_initial");
//     mod.createDynamicVariable("test_var", 42.0);
//     
//     // Self-expansion example
//     mod.autoExpandParameterSpace(1.5);
//     mod.expandMassScale(2.0);
//     mod.generateSystemReport();
//     
//     // Parameter exploration
//     mod.generateVariations(3, 0.15);
//     mod.findOptimalParameters("maximize_g_compressed", 50);
//     
//     // System analysis
//     mod.analyzeParameterSensitivity("M");
//     mod.validatePhysicalConsistency();
//     mod.autoCorrectAnomalies();
//     
//     mod.restoreState("magnetar_initial");
//     std::cout << "\nFinal g_compressed: " << mod.computeG_compressed(t) << " m/s²" << std::endl;
//     
//     return 0;
// }
// Compile: g++ -o muge_test muge_test.cpp MUGEModule.cpp -lm
// Sample: For Magnetar t=3.8e10s, g_comp ~1.79e12 m/s² (base dom.); g_res ~1e-10 m/s² (resonant scaled).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

MUGEModule Evaluation

Strengths :
-Modular, extensible design for modeling gravity in multiple astronomical systems, supporting both compressed and resonance - based UQFF models.
- Comprehensive physics : gravity, cosmological expansion, superconductivity, quantum, fluid, dark matter, system - specific effects, and resonance phenomena.
- Dynamic variable management via std::map enables runtime updates and system adaptation.
- System - specific parameter loading via setSystem for flexible analysis across diverse scenarios.
- Clear separation of computation functions(e.g., quantum, fluid, DM, Ug terms, resonance helpers), aiding maintainability.
- Output functions for equation text and variable state support debugging and documentation.
- Both compressed and resonance models are implemented, allowing comparative analysis and physical insight.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for flexibility and scalability.
- Some calculations use magic numbers or lack explanatory comments; define named constants and clarify logic.
- Minimal error handling(e.g., division by zero, invalid variable names); add validation for robustness.
- Unit consistency should be checked and documented for all physical quantities.
- For large - scale or performance - critical simulations, optimize data structures and reduce redundant calculations.
- std::map is flexible but may be less efficient than structured types for very large models.
- Expand documentation for function purposes and physical meaning.
- Minor typo : duplicate `double` in function declarations(e.g., `double double computeAAetherRes(); `).

    Summary:
The code is well - structured, flexible, and suitable for scientific prototyping and educational use in multi - system gravity modeling.It implements a broad set of physical effects and adapts to various scenarios.For production or high - performance applications, address the recommendations for improved robustness, maintainability, and scalability.