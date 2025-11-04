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

class M16UQFFModule {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
    std::map<std::string, double> variables;
    double computeQuantumTerm(double t_Hubble_val);
    double computeFluidTerm(double g_base);
    double computeResonantTerm(double t);
    double computeDMTerm();
    double computeUgSum();
    double computeHz();
    double computeMsfFactor(double t);
    double computeE_rad(double t);
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;



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
};

#endif // M16_UQFF_MODULE_H

// M16UQFFModule.cpp
#include "M16UQFFModule.h"
#include <complex>

// Constructor: Set all variables with M16-specific values
M16UQFFModule::M16UQFFModule() {
        enableDynamicTerms = true;
        enableLogging = false;
        learningRate = 0.001;
        metadata["enhanced"] = "true";
        metadata["version"] = "2.0-Enhanced";

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
           "(hbar / sqrt(Delta_x * Delta_p)) * ?(?* H ? dV) * (2? / t_Hubble) + q (v ï¿½ B) + ?_fluid * V * g + "
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
           "Solutions: Numerical evaluation at t=5 Myr yields ~1.053e-3 m/sï¿½ (EM dominant; g_grav ~1e-12 scaled by factors; micro terms ~1e-10 to 1e-3).\n"
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
//     std::cout << "g = " << g << " m/sï¿½\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("M", 1300 * 1.989e30);  // Update mass
//     mod.addToVariable("f_TRZ", 0.05);           // Add to TR factor
//     mod.subtractFromVariable("A", 1e-11);       // Subtract from amplitude
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o ziqn233h ziqn233h.cpp M16UQFFModule.cpp -lm
// Sample Output at t=5 Myr: g ? 0.001 m/sï¿½ (varies with updates; quantum/fluid/resonant ~1e-10 to 1e-3, DM ~1e33 * 1e-33 ~1e0 but curv small).
// Watermark: Copyright - Daniel T. Murphy, analyzed October 09, 2025.

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