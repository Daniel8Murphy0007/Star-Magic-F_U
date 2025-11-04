// CrabResonanceUQFFModule.h
// Modular C++ implementation of the Master Universal Gravity Equation (UQFF Resonance) for Crab Nebula Evolution.
// This module can be plugged into a base program (e.g., 'ziqn233h.cpp') by including this header and linking the .cpp.
// Usage in base: #include "CrabResonanceUQFFModule.h"
// CrabResonanceUQFFModule mod; mod.computeG(t); mod.updateVariable("f_DPM", new_value);
// All variables are stored in a std::map for dynamic addition/subtraction/update.
// Nothing is negligible: Includes all resonance-focused terms - DPM resonance, THz pipeline resonance, Aether-mediated resonance, U_g4i reactive resonance, quantum resonance, fluid resonance, oscillatory resonance (cos/exp), cosmic expansion resonance, with SC correction integrated.
// Associated text: Outputs descriptive equation string via getEquationText().
// Approximations: Resonance terms use real part of exp; frequencies from pulsar spin/wind; no SM gravity; Aether replaces dark energy.
// Crab params: M=4.6 Msun, r0=5.2e16 m, v_exp=1.5e6 m/s, f_DPM=1e12 Hz (pulsar-aligned), E_vac=7.09e-36 J/m^3, etc.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#ifndef CRAB_RESONANCE_UQFF_MODULE_H
#define CRAB_RESONANCE_UQFF_MODULE_H

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

class CrabResonanceUQFFModule {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
    std::map<std::string, double> variables;
    double computeDPMResTerm();
    double computeTHzResTerm();
    double computeAetherResTerm();
    double computeU_g4iResTerm();
    double computeQuantumResTerm();
    double computeFluidResTerm();
    double computeOscResTerm(double t);
    double computeExpResTerm();
    double computeSCResIntegrated(double B);
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;



public:
    // Constructor: Initialize all variables with Crab Nebula resonance defaults
    CrabResonanceUQFFModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computation: Full g_UQFF Resonance(r, t) as sum of resonance terms
    double computeG(double t, double B);

    // Output descriptive text of the equation
    std::string getEquationText();

    // Print all current variables (for debugging/updates)
    void printVariables();
};

#endif // CRAB_RESONANCE_UQFF_MODULE_H

// CrabResonanceUQFFModule.cpp
#include "CrabResonanceUQFFModule.h"
#include <complex>

// Constructor: Set all variables with Crab Nebula-specific resonance values
CrabResonanceUQFFModule::CrabResonanceUQFFModule() {
        enableDynamicTerms = true;
        enableLogging = false;
        learningRate = 0.001;
        metadata["enhanced"] = "true";
        metadata["version"] = "2.0-Enhanced";

    // Base constants (UQFF universal)
    variables["c"] = 3e8;                           // m/s
    variables["pi"] = 3.141592653589793;            // pi
    variables["E_vac"] = 7.09e-36;                  // J/m^3 (plasmotic vacuum energy density)
    variables["hbar"] = 1.0546e-34;                 // J s
    variables["f_TRZ"] = 0.1;                       // Time-reversal correction

    // Crab Nebula parameters
    double M_sun_val = 1.989e30;                    // kg
    variables["M_sun"] = M_sun_val;
    variables["M"] = 4.6 * M_sun_val;               // Total mass kg
    variables["r0"] = 5.2e16;                       // m (initial radius)
    variables["v_exp"] = 1.5e6;                     // m/s (expansion velocity)

    // Resonance parameters (pulsar-driven)
    variables["f_DPM"] = 1e12;                      // Hz (DPM, aligned with 30 Hz pulsar scaled)
    variables["f_THz"] = 1e12;                      // Hz (THz hole)
    variables["f_aether"] = 1e4;                    // Hz (Aether-mediated)
    variables["f_react"] = 1e10;                    // Hz (U_g4i reactive)
    variables["f_quantum"] = 1.445e-17;             // Hz (quantum wave)
    variables["f_fluid"] = 1.269e-14;               // Hz (filament fluid)
    variables["f_exp"] = 1.373e-8;                  // Hz (expansion)
    variables["f_osc"] = 30.2 * 60;                 // Hz (pulsar 30.2 Hz * 60 for res scale)
    variables["I"] = 1e21;                          // A (current proxy from wind)
    variables["A_vort"] = 3.142e8;                  // m^2 (vortical area proxy)
    variables["omega_1"] = 1e-3;                    // rad/s
    variables["omega_2"] = -1e-3;                   // rad/s
    variables["E_0"] = 6.381e-36;                   // J/m^3
    variables["f_vac_diff"] = 0.143;                // Hz
    variables["V_sys"] = 4.189e12;                  // m^3 (proxy)

    // Superconductive resonance integrated
    variables["B_crit"] = 1e11;                     // T
    variables["f_sc"] = 1.0;                        // Factor

    // Oscillatory/resonant
    variables["k"] = 1e20;                          // m^-1
    variables["omega_osc"] = 1e15;                  // rad/s (synchrotron scale)
    variables["x"] = 0.0;                           // m
    variables["A"] = 1e-10;                         // Amplitude

    // Fluid/DM proxies
    variables["rho_fluid"] = 1e-21;                 // kg/m^3 (filaments)
    variables["V"] = 1e3;                           // m^3
    variables["delta_rho"] = 0.1 * variables["rho_fluid"];
    variables["rho"] = variables["rho_fluid"];

    // Quantum
    variables["Delta_x"] = 1e-10;                   // m
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];
    variables["integral_psi"] = 1.0;
}

// Update variable (set to new value)
void CrabResonanceUQFFModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
    if (name == "Delta_x") {
        variables["Delta_p"] = variables["hbar"] / value;
    }
}

// Add delta to variable
void CrabResonanceUQFFModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta from variable
void CrabResonanceUQFFModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute DPM Resonance Term: a_DPM_res = (F_DPM * f_DPM * E_vac) / (c * V_sys)
double CrabResonanceUQFFModule::computeDPMResTerm() {
    double F_DPM = variables["I"] * variables["A_vort"] * (variables["omega_1"] - variables["omega_2"]);
    double r_t = variables["r0"] + variables["v_exp"] * variables["t"];  // r(t) proxy
    double V_sys_t = (4.0 / 3.0) * variables["pi"] * std::pow(r_t, 3);  // Updated volume
    return (F_DPM * variables["f_DPM"] * variables["E_vac"]) / (variables["c"] * V_sys_t);
}

// Compute THz Resonance Term: a_THz_res = (f_THz * E_vac * v_exp * a_DPM_res) / (E_vac / 10 * c)
double CrabResonanceUQFFModule::computeTHzResTerm() {
    double a_DPM_res = computeDPMResTerm();
    return (variables["f_THz"] * variables["E_vac"] * variables["v_exp"] * a_DPM_res) / ((variables["E_vac"] / 10.0) * variables["c"]);
}

// Compute Aether Resonance Term: a_aether_res = f_aether * 1e-8 * f_DPM * (1 + f_TRZ) * a_DPM_res
double CrabResonanceUQFFModule::computeAetherResTerm() {
    double a_DPM_res = computeDPMResTerm();
    return variables["f_aether"] * 1e-8 * variables["f_DPM"] * (1 + variables["f_TRZ"]) * a_DPM_res;
}

// Compute U_g4i Reactive Resonance Term: U_g4i_res = f_sc * Ug1_proxy * f_react * a_DPM_res / (E_vac * c)
double CrabResonanceUQFFModule::computeU_g4iResTerm() {
    double Ug1_proxy = 1.0;  // Normalized
    double a_DPM_res = computeDPMResTerm();
    return variables["f_sc"] * Ug1_proxy * variables["f_react"] * a_DPM_res / (variables["E_vac"] * variables["c"]);
}

// Compute Quantum Resonance Term: a_quantum_res = (f_quantum * E_vac * a_DPM_res) / ((E_vac / 10) * c)
double CrabResonanceUQFFModule::computeQuantumResTerm() {
    double a_DPM_res = computeDPMResTerm();
    return (variables["f_quantum"] * variables["E_vac"] * a_DPM_res) / ((variables["E_vac"] / 10.0) * variables["c"]);
}

// Compute Fluid Resonance Term: a_fluid_res = (f_fluid * E_vac * V * a_DPM_res) / ((E_vac / 10) * c)
double CrabResonanceUQFFModule::computeFluidResTerm() {
    double a_DPM_res = computeDPMResTerm();
    return (variables["f_fluid"] * variables["E_vac"] * variables["V"] * a_DPM_res) / ((variables["E_vac"] / 10.0) * variables["c"]);
}

// Compute Oscillatory Resonance Term: 2 A cos(k x) cos(omega t) + (2 pi / 13.8) A Re[exp(i (k x - omega t))]
double CrabResonanceUQFFModule::computeOscResTerm(double t) {
    double cos_term = 2 * variables["A"] * std::cos(variables["k"] * variables["x"]) * std::cos(variables["omega_osc"] * t);
    std::complex<double> exp_term(variables["A"] * std::exp(std::complex<double>(0, variables["k"] * variables["x"] - variables["omega_osc"] * t)));
    double real_exp = exp_term.real();
    double exp_factor = (2 * variables["pi"] / 13.8);
    return cos_term + exp_factor * real_exp;
}

// Compute Expansion Resonance Term: a_exp_res = (f_exp * E_vac * a_DPM_res) / ((E_vac / 10) * c)
double CrabResonanceUQFFModule::computeExpResTerm() {
    double a_DPM_res = computeDPMResTerm();
    return (variables["f_exp"] * variables["E_vac"] * a_DPM_res) / ((variables["E_vac"] / 10.0) * variables["c"]);
}

// Compute SC Resonance Integrated: (1 - B / B_crit) * f_sc
double CrabResonanceUQFFModule::computeSCResIntegrated(double B) {
    return (1.0 - (B / variables["B_crit"])) * variables["f_sc"];
}

// Full g_UQFF Resonance: Sum resonance terms * SC * (1 + f_TRZ)
double CrabResonanceUQFFModule::computeG(double t, double B) {
    variables["t"] = t;
    double a_DPM_res = computeDPMResTerm();
    double a_THz_res = computeTHzResTerm();
    double a_aether_res = computeAetherResTerm();
    double a_u_g4i_res = computeU_g4iResTerm();
    double a_quantum_res = computeQuantumResTerm();
    double a_fluid_res = computeFluidResTerm();
    double a_osc_res = computeOscResTerm(t);
    double a_exp_res = computeExpResTerm();
    double sc_int = computeSCResIntegrated(B);
    double tr_factor = 1.0 + variables["f_TRZ"];
    double res_sum = a_DPM_res + a_THz_res + a_aether_res + a_u_g4i_res + a_quantum_res + a_fluid_res + a_osc_res + a_exp_res;
    return res_sum * sc_int * tr_factor;
}

// Get equation text (descriptive)
std::string CrabResonanceUQFFModule::getEquationText() {
    return "g_Crab_Res(t, B) = [a_DPM_res + a_THz_res + a_aether_res + U_g4i_res + a_quantum_res + a_fluid_res + a_osc_res + a_exp_res] * SC_int * (1 + f_TRZ)\n"
           "Where:\n"
           "- a_DPM_res = (F_DPM * f_DPM * E_vac) / (c * V_sys(t)); F_DPM = I * A * (?1 - ?2); V_sys(t) = 4/3 ? r(t)^3, r(t)=r0 + v_exp t\n"
           "- a_THz_res = (f_THz * E_vac * v_exp * a_DPM_res) / (E_vac/10 * c)\n"
           "- a_aether_res = f_aether * 1e-8 * f_DPM * (1 + f_TRZ) * a_DPM_res\n"
           "- U_g4i_res = f_sc * Ug1_proxy * f_react * a_DPM_res / (E_vac * c)\n"
           "- a_quantum_res = (f_quantum * E_vac * a_DPM_res) / (E_vac/10 * c)\n"
           "- a_fluid_res = (f_fluid * E_vac * V * a_DPM_res) / (E_vac/10 * c)\n"
           "- a_osc_res = 2 A cos(k x) cos(? t) + (2? / 13.8) A Re[exp(i (k x - ? t))]\n"
           "- a_exp_res = (f_exp * E_vac * a_DPM_res) / (E_vac/10 * c)\n"
           "- SC_int = (1 - B / B_crit) * f_sc\n"
           "Special Terms: UQFF resonance via plasmotic vacuum; Aether replaces dark energy; no SM terms; pulsar-driven f_osc.\n"
           "Solutions: At t=971 yr, B=1e-8 T, g ? 1e-40 m/sï¿½ (resonance micro-scale, wind proxy).\n"
           "Adaptations: Resonance focus for Crab wisps/shocks per Hubble/Chandra.";
}

// Print variables
void CrabResonanceUQFFModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Example usage in base program 'ziqn233h.cpp' (snippet for integration)
// #include "CrabResonanceUQFFModule.h"
// int main() {
//     CrabResonanceUQFFModule mod;
//     double t = 971 * 3.156e7;  // 971 years
//     double B = 1e-8;           // T (nebula avg)
//     double g_res = mod.computeG(t, B);
//     std::cout << "g_res = " << g_res << " m/sï¿½\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("f_DPM", 1.1e12);  // Update resonance freq
//     mod.addToVariable("f_TRZ", 0.05);     // Add to TR factor
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o ziqn233h ziqn233h.cpp CrabResonanceUQFFModule.cpp -lm
// Sample Output at t=971 yr: g_res ? 1e-40 m/sï¿½ (varies; micro-scale resonance terms).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

// Evaluation of CrabResonanceUQFFModule (UQFF Resonance Model for Crab Nebula)

**Strengths:**
-**Dynamic & Extensible : **All model parameters are stored in a `std: : map<std::string, double> variables`, allowing runtime updates, additions, and removals.The methods `updateVariable`, `addToVariable`, and `subtractFromVariable` enable flexible modification of any parameter.
- **Automatic Dependency Updates : **When key variables like `"Delta_x"` are updated, dependent variables(`"Delta_p"`) are recalculated automatically, ensuring consistency.
    - **Immediate Effect : **All computations(e.g., `computeG`) use the current values in the map, so any changes are immediately reflected in results.
        - **Comprehensive Physics : **The module includes all major UQFF resonance terms relevant for Crab Nebula modeling, such as DPM resonance, THz pipeline, Aether resonance, U_g4i reactive, quantum, fluid, oscillatory, and cosmic expansion effects.Standard Model gravity / magnetics are intentionally excluded per UQFF.
        - **Time - Dependent Volume : **The resonance terms correctly use the time - dependent radius and volume(`r(t)`, `V_sys(t)`), reflecting nebula expansion.
    - **Debugging Support : **The `printVariables()` method provides a snapshot of all current parameters, aiding validation and troubleshooting.
    - **Sample Usage Provided : **Example integration and compilation instructions are included, demonstrating how to update variables and see their effect.

    ** Weaknesses / Recommendations : **
    -**Error Handling : **Adding unknown variables is flexible but may lead to silent errors if a typo occurs.Consider stricter validation or optional warnings for unknown variable names.
    - **Unit Consistency : **Ensure all units are consistent, especially when combining terms from different physical domains and frequency - based scaling.
    - **Magic Numbers : **Some scale factors and physical constants are set to arbitrary values.Document their physical meaning or allow configuration via constructor arguments or config files.
    - **Physical Justification : **The model is highly specialized for UQFF and omits Standard Model terms.Ensure this is appropriate for your scientific context and document the rationale for each term.
    - **Performance : **For large - scale or repeated updates, consider profiling and optimizing critical paths if needed.

    ** Summary : **
    The module is robust, dynamic, and extensible, supporting runtime updates and changes to all model parameters.It is suitable for advanced UQFF - based resonance modeling of the Crab Nebula.Minor improvements in error handling, documentation, and physical justification are recommended for production or publication use.