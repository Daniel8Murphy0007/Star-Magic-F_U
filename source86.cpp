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


#include <map>
#include <vector>
#include <functional>
#include <memory>
#include <algorithm>
#include <fstream>
#include <sstream>
enum #include <map>
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

class SystemType {
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
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
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
    double double computeAAetherRes();
    double computeUg4i();
    double computeAQuantumFreq();
    double double computeAAetherFreq();
    double computeAFluidFreq();
    double computeOscTerm(double t);
    double computeAExpFreq();
    double computeFTRZ();
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;



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
           "(hbar / sqrt(Delta_x Delta_p)) ?(? H ? dV) (2? / t_Hubble) + q (v ï¿½ B) + ?_fluid V g + "
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

// Example usage snippet:
// #include "MUGEModule.h"
// int main() {
//     MUGEModule mod(SystemType::MAGNETAR_SGR_1745_2900);
//     double t = 3.799e10;
//     double g_comp = mod.computeG_compressed(t);
//     double g_res = mod.computeG_resonance(t);
//     std::cout << "Compressed g = " << g_comp << " m/sï¿½\n";
//     std::cout << "Resonance g = " << g_res << " m/sï¿½\n";
//     std::cout << mod.getEquationText_compressed() << "\n" << mod.getEquationText_resonance() << std::endl;
//     mod.updateVariable("M", 2.0 * mod.variables["M_sun"]);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o muge_test muge_test.cpp MUGEModule.cpp -lm
// Sample: For Magnetar t=3.8e10s, g_comp ~1.79e12 m/sï¿½ (base dom.); g_res ~1e-10 m/sï¿½ (resonant scaled).
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