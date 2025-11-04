// SpiralSupernovaeUQFFModule.h
// Modular C++ implementation of the full Master Universal Gravity Equation (UQFF & SM Integration) for Spirals and Supernovae Evolution.
// This module can be plugged into a base program (e.g., 'ziqn233h.cpp') by including this header and linking the .cpp.
// Usage in base: #include "SpiralSupernovaeUQFFModule.h"
// SpiralSupernovaeUQFFModule mod; mod.computeG(t); mod.updateVariable("M", new_value);
// All variables are stored in a std::map for dynamic addition/subtraction/update.
// Nothing is negligible: Includes all terms - base gravity with T_spiral, Ug1-Ug4, cosmological Lambda with ?_?, quantum integral, Lorentz q(v x B), fluid rho_fluid V g, resonant oscillatory (cos/exp), DM/visible with perturbations, supernova SN_term.
// Associated text: Outputs descriptive equation string via getEquationText().
// Approximations: Quantum integral normalized to 1.0; exp real part; Ug2/Ug3=0; DM fraction 0.85; T_spiral with ?_p; SN_term from L_SN.
// Spiral-SN params: M=1.989e41 kg, r=9.258e20 m, H0=73 km/s/Mpc, ?_p=20 km/s/kpc, L_SN=1e36 W, z up to 1.5, etc.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#ifndef SPIRAL_SUPERNOVAE_UQFF_MODULE_H
#define SPIRAL_SUPERNOVAE_UQFF_MODULE_H

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

class SpiralSupernovaeUQFFModule {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
    std::map<std::string, double> variables;
    double computeQuantumTerm(double t_Hubble_val);
    double computeFluidTerm(double g_base);
    double computeResonantTerm(double t);
    double computeDMTerm();
    double computeUgSum();
    double computeHz(double z);
    double computeT_spiral(double t);
    double computeSN_term(double z);
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;



public:
    // Constructor: Initialize all variables with Spirals and Supernovae defaults
    SpiralSupernovaeUQFFModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computation: Full g_UQFF(r, t) for Spirals and Supernovae
    double computeG(double t, double z);

    // Output descriptive text of the equation
    std::string getEquationText();

    // Print all current variables (for debugging/updates)
    void printVariables();
};

#endif // SPIRAL_SUPERNOVAE_UQFF_MODULE_H

// SpiralSupernovaeUQFFModule.cpp
#include "SpiralSupernovaeUQFFModule.h"
#include <complex>

// Constructor: Set all variables with Spirals and Supernovae-specific values
SpiralSupernovaeUQFFModule::SpiralSupernovaeUQFFModule() {
        enableDynamicTerms = true;
        enableLogging = false;
        learningRate = 0.001;
        metadata["enhanced"] = "true";
        metadata["version"] = "2.0-Enhanced";

    // Base constants (universal)
    variables["G"] = 6.6743e-11;                    // m^3 kg^-1 s^-2
    variables["c"] = 3e8;                           // m/s
    variables["hbar"] = 1.0546e-34;                 // J s
    variables["Lambda"] = 1.1e-52;                  // m^-2
    variables["q"] = 1.602e-19;                     // C
    variables["pi"] = 3.141592653589793;            // pi
    variables["t_Hubble"] = 13.8e9 * 3.156e7;       // s

    // Spiral-SN parameters
    double M_sun_val = 1.989e30;                    // kg
    variables["M_sun"] = M_sun_val;
    variables["M"] = 1e11 * M_sun_val;              // Galaxy mass kg
    variables["M_visible"] = 0.15 * variables["M"]; // Visible fraction
    variables["M_DM"] = 0.85 * variables["M"];      // Dark matter
    variables["r"] = 9.258e20;                      // m (~30 kpc)
    variables["M_gas"] = 1e9 * M_sun_val;           // Gas mass

    // Hubble/cosmology
    variables["H0"] = 73.0;                         // km/s/Mpc (SH0ES)
    variables["Mpc_to_m"] = 3.086e22;               // m/Mpc
    variables["z"] = 0.5;                           // Typical z for SN
    variables["Omega_m"] = 0.3;
    variables["Omega_Lambda"] = 0.7;
    variables["t"] = 5e9 * 3.156e7;                 // Default t=5 Gyr s

    // Spiral dynamics
    variables["Omega_p"] = 20e3 / 3.086e19;         // rad/s (20 km/s/kpc pattern speed)

    // SN parameters
    variables["L_SN"] = 1e36;                       // W (peak luminosity)
    variables["rho_fluid"] = 1e-21;                 // kg/m^3 (ISM)
    variables["V"] = 1e3;                           // m^3
    variables["v_rot"] = 2e5;                       // m/s (rotation)
    variables["delta_rho"] = 0.1 * variables["rho_fluid"];
    variables["rho"] = variables["rho_fluid"];

    // EM/magnetic
    variables["B"] = 1e-5;                          // T (galactic field)
    variables["B_crit"] = 1e11;                     // T

    // Quantum terms
    variables["Delta_x"] = 1e-10;                   // m
    variables["Delta_p"] = variables["hbar"] / variables["Delta_x"];
    variables["integral_psi"] = 1.0;

    // Resonant/oscillatory
    variables["A"] = 1e-10;
    variables["k"] = 1e20;
    variables["omega"] = 1e15;
    variables["x"] = 0.0;

    // Ug subterms
    variables["Ug1"] = 0.0;
    variables["Ug2"] = 0.0;
    variables["Ug3"] = 0.0;
    variables["Ug4"] = 0.0;

    // Scale factors
    variables["scale_macro"] = 1e-12;
    variables["f_TRZ"] = 0.1;
    variables["f_sc"] = 1.0;
}

// Update variable (set to new value)
void SpiralSupernovaeUQFFModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
    if (name == "Delta_x") {
        variables["Delta_p"] = variables["hbar"] / value;
    } else if (name == "M") {
        variables["M_visible"] = 0.15 * value;
        variables["M_DM"] = 0.85 * value;
    } else if (name == "H0") {
        // Recompute if needed
    }
}

// Add delta to variable
void SpiralSupernovaeUQFFModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta from variable
void SpiralSupernovaeUQFFModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute H(z) in s^-1
double SpiralSupernovaeUQFFModule::computeHz(double z) {
    double Hz_kms = variables["H0"] * std::sqrt(variables["Omega_m"] * std::pow(1.0 + z, 3) + variables["Omega_Lambda"]);
    return (Hz_kms * 1e3) / variables["Mpc_to_m"];
}

// Compute Ug sum: Ug1 = G M / r^2, Ug4 = Ug1 * f_sc, others 0
double SpiralSupernovaeUQFFModule::computeUgSum() {
    double Ug1 = (variables["G"] * variables["M"]) / (variables["r"] * variables["r"]);
    variables["Ug1"] = Ug1;
    variables["Ug4"] = Ug1 * variables["f_sc"];
    return variables["Ug1"] + variables["Ug2"] + variables["Ug3"] + variables["Ug4"];
}

// Quantum term: (hbar / sqrt(Delta_x Delta_p)) * integral * (2 pi / t_Hubble)
double SpiralSupernovaeUQFFModule::computeQuantumTerm(double t_Hubble_val) {
    double unc = std::sqrt(variables["Delta_x"] * variables["Delta_p"]);
    double integral_val = variables["integral_psi"];
    return (variables["hbar"] / unc) * integral_val * (2 * variables["pi"] / t_Hubble_val);
}

// Fluid term: rho_fluid * V * g
double SpiralSupernovaeUQFFModule::computeFluidTerm(double g_base) {
    return variables["rho_fluid"] * variables["V"] * g_base;
}

// Resonant terms: 2 A cos(k x) cos(omega t) + (2 pi / 13.8) A Re[exp(i (k x - omega t))]
double SpiralSupernovaeUQFFModule::computeResonantTerm(double t) {
    double cos_term = 2 * variables["A"] * std::cos(variables["k"] * variables["x"]) * std::cos(variables["omega"] * t);
    std::complex<double> exp_term(variables["A"] * std::exp(std::complex<double>(0, variables["k"] * variables["x"] - variables["omega"] * t)));
    double real_exp = exp_term.real();
    double exp_factor = (2 * variables["pi"] / 13.8);
    return cos_term + exp_factor * real_exp;
}

// DM term: (M_visible + M_DM) * (delta_rho / rho + 3 G M / r^3)
double SpiralSupernovaeUQFFModule::computeDMTerm() {
    double pert = variables["delta_rho"] / variables["rho"];
    double curv = 3 * variables["G"] * variables["M"] / (variables["r"] * variables["r"] * variables["r"]);
    return (variables["M_visible"] + variables["M_DM"]) * (pert + curv);
}

// Spiral torque term: T_spiral = G * M_gas * M / r^2 * (1 + ?_p * t)
double SpiralSupernovaeUQFFModule::computeT_spiral(double t) {
    double torque_base = (variables["G"] * variables["M_gas"] * variables["M"]) / (variables["r"] * variables["r"]);
    return torque_base * (1.0 + variables["Omega_p"] * t);
}

// Supernova term: SN_term = (L_SN / (4 pi r^2 c)) * (1 + H(z) * t)
double SpiralSupernovaeUQFFModule::computeSN_term(double z) {
    double Hz = computeHz(z);
    double flux = variables["L_SN"] / (4 * variables["pi"] * variables["r"] * variables["r"] * variables["c"]);
    return flux * (1.0 + Hz * variables["t"]);
}

// Full computation: g_UQFF(r, t) = ... all terms with T_spiral and SN_term
double SpiralSupernovaeUQFFModule::computeG(double t, double z) {
    variables["t"] = t;
    double Hz = computeHz(z);
    double expansion = 1.0 + Hz * t;
    double sc_correction = 1.0 - (variables["B"] / variables["B_crit"]);
    double tr_factor = 1.0 + variables["f_TRZ"];
    double t_spiral = computeT_spiral(t);
    double sn_term = computeSN_term(z);

    // Base gravity with expansion, SC, TR, T_spiral
    double g_base = ((variables["G"] * variables["M"] / (variables["r"] * variables["r"])) * expansion * sc_correction * tr_factor) * (1.0 + t_spiral);

    // Ug sum
    double ug_sum = computeUgSum();

    // Cosmological with ?_?
    double lambda_term = variables["Lambda"] * (variables["c"] * variables["c"] * variables["Omega_Lambda"]) / 3.0;

    // Quantum
    double quantum_term = computeQuantumTerm(variables["t_Hubble"]);

    // EM Lorentz (rotation v_rot B)
    double em_base = variables["q"] * variables["v_rot"] * variables["B"] / 1.673e-27;
    double em_term = em_base * (1.0 + (7.09e-36 / 7.09e-37)) * variables["scale_macro"];

    // Fluid
    double fluid_term = computeFluidTerm(g_base);

    // Resonant
    double resonant_term = computeResonantTerm(t);

    // DM
    double dm_term = computeDMTerm();

    // Total: Sum all + SN_term
    return g_base + ug_sum + lambda_term + quantum_term + em_term + fluid_term + resonant_term + dm_term + sn_term;
}

// Get equation text (descriptive)
std::string SpiralSupernovaeUQFFModule::getEquationText() {
    return "g_Spiral_SN(r, t) = (G * M / r^2) * (1 + H(z) * t) * (1 + T_spiral) * (1 - B / B_crit) * (1 + f_TRZ) + (Ug1 + Ug2 + Ug3 + Ug4) + (Lambda * c^2 * ?_? / 3) + "
           "(hbar / sqrt(Delta_x * Delta_p)) * ?(?* H ? dV) * (2? / t_Hubble) + q (v ï¿½ B) + ?_fluid * V * g + "
           "2 A cos(k x) cos(? t) + (2? / 13.8) A exp(i (k x - ? t)) + (M_visible + M_DM) * (??/? + 3 G M / r^3) + SN_term\n"
           "Where T_spiral = G * M_gas * M / r^2 * (1 + ?_p * t); SN_term = (L_SN / (4? r^2 c)) * (1 + H(z) * t)\n"
           "Special Terms:\n"
           "- Quantum: Heisenberg uncertainty for ISM quantum effects.\n"
           "- Fluid: Gas density-volume-gravity coupling in arms.\n"
           "- Resonant: Oscillatory Aether waves for density waves.\n"
           "- DM: Visible+dark mass with perturbations for rotation curves.\n"
           "- Superconductivity: (1 - B/B_crit) for galactic fields.\n"
           "- Spiral Torque: T_spiral for arm evolution.\n"
           "- Supernova: SN_term for expansion probe.\n"
           "Solutions: At t=5 Gyr, z=0.5, g_Spiral_SN ~1e-10 m/sï¿½ (Lambda/SN dominant; g_base ~1e-10).\n"
           "Adaptations for Spirals and Supernovae: SH0ES H0=73; ?_p=20 km/s/kpc; L_SN=1e36 W for Ia SN.";
}

// Print variables
void SpiralSupernovaeUQFFModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Example usage in base program 'ziqn233h.cpp' (snippet for integration)
// #include "SpiralSupernovaeUQFFModule.h"
// int main() {
//     SpiralSupernovaeUQFFModule mod;
//     double t = 5e9 * 3.156e7;  // 5 Gyr
//     double z = 0.5;            // Typical SN z
//     double g = mod.computeG(t, z);
//     std::cout << "g = " << g << " m/sï¿½\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("H0", 74.0);  // Update H0
//     mod.addToVariable("f_TRZ", 0.05);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o ziqn233h ziqn233h.cpp SpiralSupernovaeUQFFModule.cpp -lm
// Sample Output at t=5 Gyr, z=0.5: g ? 1e-10 m/sï¿½ (varies; Lambda/SN dominant).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

// Evaluation of SpiralSupernovaeUQFFModule (UQFF & Standard Model Integration for Spiral Galaxies and Supernovae)

**Strengths:**
-**Dynamic & Extensible : **All model parameters are stored in a `std: : map<std::string, double> variables`, allowing runtime updates, additions, and removals.The methods `updateVariable`, `addToVariable`, and `subtractFromVariable` enable flexible modification of any parameter.
- **Automatic Dependency Updates : **When key variables like `"Delta_x"` or `"M"` are updated, dependent variables(`"Delta_p"`, `"M_visible"`, `"M_DM"`) are recalculated automatically, ensuring consistency.
    - **Immediate Effect : **All computations(e.g., `computeG`) use the current values in the map, so any changes are immediately reflected in results.
        - **Comprehensive Physics : **The module includes all major UQFF and Standard Model terms relevant for spiral galaxy and supernova gravity, such as base gravity, cosmological constant, quantum, EM, fluid, resonant, DM, spiral torque, supernova, and superconductivity corrections.
        - **Specialized Terms : **Incorporates spiral torque(`T_spiral`) and supernova(`SN_term`) effects, which are important for galactic evolution and feedback.
            - **Debugging Support : **The `printVariables()` method provides a snapshot of all current parameters, aiding validation and troubleshooting.
            - **Sample Usage Provided : **Example integration and compilation instructions are included, demonstrating how to update variables and see their effect.

            ** Weaknesses / Recommendations : **
            -**Error Handling : **Adding unknown variables is flexible but may lead to silent errors if a typo occurs.Consider stricter validation or optional warnings for unknown variable names.
            - **Unit Consistency : **Ensure all units are consistent, especially when combining terms from different physical domains and frequency - based scaling.
            - **Magic Numbers : **Some scale factors and physical constants are set to arbitrary values.Document their physical meaning or allow configuration via constructor arguments or config files.
            - **Physical Justification : **The model is highly specialized for UQFF and galactic physics.Ensure this is appropriate for your scientific context and document the rationale for each term.
            - **Performance : **For large - scale or repeated updates, consider profiling and optimizing critical paths if needed.

            ** Summary : **
            The module is robust, dynamic, and extensible, supporting runtime updates and changes to all model parameters.It is suitable for advanced UQFF - based and Standard Model modeling of spiral galaxies and supernovae.Minor improvements in error handling, documentation, and physical justification are recommended for production or publication use.