
// UQFFBuoyancyModule.h
// Modular C++ implementation of the full Master Unified Field Equation (F_U_Bi_i & UQFF Integration) for Buoyancy Equations across M74, Eagle Nebula (M16), M84, Centaurus A, Supernova Survey.
// This module can be plugged into a base program (e.g., 'uqff_buoyancy_sim.cpp') by including this header and linking the .cpp.
// Usage in base: #include "UQFFBuoyancyModule.h"
// UQFFBuoyancyModule mod; mod.computeFBi(system); mod.updateVariable("F_rel", {new_real, new_imag});
// All variables are stored in a std::map for dynamic addition/subtraction/update, using complex<double> for real/imaginary components.
// Nothing is negligible: Includes all terms - base force, momentum, gravity, vacuum stability, LENR resonance, activation, directed energy, magnetic resonance, neutron, relativistic, neutrino, Sweet vac, Kozima drop.
// Associated text: Outputs descriptive equation string via getEquationText().
// Approximations: Integral approximated as integrand * x2 (quadratic root); imag parts small and not fully scaled; LENR dominant due to low ω_0; x2 from quadratic solver approx; F_rel from 1998 LEP; g(r,t) and Q_wave dynamic per system to address repetition.
// Multi-system params: M74 M=7.17e41 kg r=9.46e20 m; M16 M=1e36 kg r=2.36e17 m; M84 M=1.46e45 kg r=3.09e22 m; Centaurus A M=4e41 kg r=3.09e21 m; Supernova Survey (generic M=1e30 kg r=1e10 m).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 22, 2025.

#ifndef UQFF_BUOYANCY_MODULE_H
#define UQFF_BUOYANCY_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>

using cdouble = std::complex<double>;

class UQFFBuoyancyModule {
private:
    std::map<std::string, cdouble> variables;
    cdouble computeIntegrand(double t, const std::string& system);
    cdouble computeDPM_resonance(const std::string& system);
    cdouble computeX2(const std::string& system);
    cdouble computeQuadraticRoot(cdouble a, cdouble b, cdouble c);
    cdouble computeLENRTerm(const std::string& system);
    double computeG(double t, const std::string& system);
    cdouble computeQ_wave(double t, const std::string& system);
    cdouble computeUb1(const std::string& system);
    cdouble computeUi(double t, const std::string& system);
    void setSystemParams(const std::string& system);

public:
    // Constructor: Initialize all variables with multi-system defaults
    UQFFBuoyancyModule();

    // Dynamic variable operations (complex)
    void updateVariable(const std::string& name, cdouble value);
    void addToVariable(const std::string& name, cdouble delta);
    void subtractFromVariable(const std::string& name, cdouble delta);

    // Core computation: Full F_U_Bi_i(r, t) for system (approx integral)
    cdouble computeFBi(const std::string& system, double t);

    // Sub-equations
    cdouble computeCompressed(const std::string& system, double t);  // Integrand
    cdouble computeResonant(const std::string& system);
    cdouble computeBuoyancy(const std::string& system);
    cdouble computeSuperconductive(const std::string& system, double t);
    double computeCompressedG(const std::string& system, double t);  // g(r,t)

    // Output descriptive text of the equation
    std::string getEquationText(const std::string& system);

    // Print all current variables (for debugging/updates)
    void printVariables();
};

#endif // UQFF_BUOYANCY_MODULE_H
// UQFFBuoyancyModule.cpp
#include "UQFFBuoyancyModule.h"
// Constructor: Initialize all variables with multi-system defaults
UQFFBuoyancyModule::UQFFBuoyancyModule() {
    // Universal constants
    variables["G"] = cdouble(6.67430e-11, 0.0);          // m^3 kg^-1 s^-2
    variables["c"] = cdouble(299792458.0, 0.0);          // m/s
    variables["h_bar"] = cdouble(1.054571817e-34, 0.0);   // J·s
    variables["e_charge"] = cdouble(1.602176634e-19, 0.0); // C
    variables["epsilon_0"] = cdouble(8.854187817e-12, 0.0); // F/m
    variables["mu_0"] = cdouble(1.25663706212e-6, 0.0);   // N/A^2
    variables["pi"] = cdouble(3.141592653589793, 0.0);

    // System-specific defaults (M74 as default)
    setSystemParams("M74");

    // Dynamic variables initialized to zero
    variables["F_rel"] = cdouble(0.0, 0.0);
    variables["F_U_Bi_i"] = cdouble(0.0, 0.0);
}
// Set system-specific parameters
void UQFFBuoyancyModule::setSystemParams(const std::string& system)
{
    if (system == "M74") {
        variables["M_system"] = cdouble(7.17e41, 0.0);   // kg
        variables["r_system"] = cdouble(9.46e20, 0.0);   // m
    } else if (system == "M16") {
        variables["M_system"] = cdouble(1e36, 0.0);      // kg
        variables["r_system"] = cdouble(2.36e17, 0.0);   // m
    } else if (system == "M84") {
        variables["M_system"] = cdouble(1.46e45, 0.0);   // kg
        variables["r_system"] = cdouble(3.09e22, 0.0);   // m
    } else if (system == "CentaurusA") {
        variables["M_system"] = cdouble(4e41, 0.0);      // kg
        variables["r_system"] = cdouble(3.09e21, 0.0);   // m
    } else if (system == "SupernovaSurvey") {
        variables["M_system"] = cdouble(1e30, 0.0);      // kg
        variables["r_system"] = cdouble(1e10, 0.0);      // m
    } else {
        std::cerr << "Unknown system: " << system << ". Using M74 defaults." << std::endl;
        variables["M_system"] = cdouble(7.17e41, 0.0);   // kg
        variables["r_system"] = cdouble(9.46e20, 0.0);   // m
    }
}
// UQFFBuoyancyModule.h
// Modular C++ implementation of the Surface Magnetic Field Module for stellar magnetic field modeling (Sun).
// This module can be plugged into a base program (e.g., 'surface_magnetic_field_sim.cpp') by including this header and linking the .cpp.
// Usage in base: #include "SurfaceMagneticFieldModule.h"
// SurfaceMagneticFieldModule mod; mod.computeU_g3_example(0.0); mod.updateVariable("B_s_min", new_value);
// All variables are stored in a std::map for dynamic addition/subtraction/update.
// Associated text: Outputs descriptive equation string via getEquationText().
// Approximations: B_ref=0.4 T (max sunspot); cos(?_s t ?)=1; P_core=1; E_react=1e46.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.
#ifndef SURFACE_MAGNETIC_FIELD_MODULE_H
#define SURFACE_MAGNETIC_FIELD_MODULE_H
#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
class SurfaceMagneticFieldModule {
private:
    std::map<std::string, double> variables;
    double computeB_j(double t, double B_s);
    double computeU_g3_example(double t, double B_s);
public:
    // Constructor: Initialize with framework defaults (Sun)
    SurfaceMagneticFieldModule();
    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeB_s_min();  // 1e-4 T (quiet Sun)
    double computeB_s_max();  // 0.4 T (sunspot max)
    double computeB_j(double t, double B_s);  // Scaled B_j (T)
    double computeU_g3_example(double t, double B_s);  // U_g3 with B_j (J/m^3)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();
};

#endif // SURFACE_MAGNETIC_FIELD_MODULE_H
// SurfaceMagneticFieldModule.cpp
#include "SurfaceMagneticFieldModule.h"

SurfaceMagneticFieldModule::SurfaceMagneticFieldModule() {
    // Initialize default variables
    variables["B_s_min"] = 1e-4;
    variables["B_s_max"] = 0.4;
    variables["U_g3"] = 0.0;
}

void SurfaceMagneticFieldModule::updateVariable(const std::string& name, double value) {
    variables[name] = value;
}

void SurfaceMagneticFieldModule::addToVariable(const std::string& name, double delta) {
    variables[name] += delta;
}

void SurfaceMagneticFieldModule::subtractFromVariable(const std::string& name, double delta) {
    variables[name] -= delta;
}

double SurfaceMagneticFieldModule::computeB_s_min() {
    return variables["B_s_min"];
}

double SurfaceMagneticFieldModule::computeB_s_max() {
    return variables["B_s_max"];
}

double SurfaceMagneticFieldModule::computeB_j(double t, double B_s) {
    // Placeholder implementation
    return B_s * std::sin(t);
}

double SurfaceMagneticFieldModule::computeU_g3_example(double t, double B_s) {
    // Placeholder implementation
    return B_s * std::cos(t);
}

std::string SurfaceMagneticFieldModule::getEquationText() {
    std::ostringstream oss;
    oss << "B_s_min = " << computeB_s_min() << " T\n";
    oss << "B_s_max = " << computeB_s_max() << " T\n";
    return oss.str();
}

void SurfaceMagneticFieldModule::printVariables() {
    for (const auto& var : variables) {
        std::cout << var.first << " = " << var.second << "\n";
    }
}
// UQFFBuoyancyModule.cpp
#include "UQFFBuoyancyModule.h"
#include <complex>

// Constructor: Set all variables with multi-system defaults
UQFFBuoyancyModule::UQFFBuoyancyModule() {
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
    variables["omega_LENR"] = {2 * pi_val * 1e9, 0.0};
    variables["k_LENR"] = {1e-9, 0.0};
    variables["x2"] = {1e3, 0.0};  // Quadratic root approx
    variables["DPM_resonance"] = {1e5, 0.0};
    variables["LENR_term"] = {1e10, 0.0};
    variables["M_system"] = {7.17e41, 0.0};  // M74 default
    variables["r_system"] = {9.46e20, 0.0};
    variables["F_U_Bi_i"] = zero;
}
// Core computation: Full F_U_Bi_i(r, t) for system (approx integral)
cdouble UQFFBuoyancyModule::computeFBi(const std::string& system, double t) {
    setSystemParams(system);
    cdouble integrand = computeIntegrand(t, system);
    cdouble x2 = computeX2(system);
    cdouble result = integrand * x2;
    return result;
}
// Dynamic variable operations (complex)
void UQFFBuoyancyModule::updateVariable(const std::string& name, cdouble value) {
    variables[name] = value;
}
void UQFFBuoyancyModule::addToVariable(const std::string& name, cdouble delta) {
    variables[name] += delta;
}
void UQFFBuoyancyModule::subtractFromVariable(const std::string& name, cdouble delta) {
    variables[name] -= delta;
}
// Core computation: Full F_U_Bi_i(r, t) for system (approx integral)
cdouble UQFFBuoyancyModule::computeFBi(const std::string& system, double t) {
    setSystemParams(system);
    cdouble integrand = computeIntegrand(t, system);
    cdouble x2 = computeX2(system);
    cdouble result = integrand * x2;
    return result;
}

// ===== SURFACE MAGNETIC FIELD MODULE FOR STELLAR SYSTEMS =====

#ifndef SURFACE_MAGNETIC_FIELD_MODULE_STELLAR_H
#define SURFACE_MAGNETIC_FIELD_MODULE_STELLAR_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

class SurfaceMagneticFieldModule {
private:
    std::map<std::string, double> variables;
    std::map<std::string, std::vector<double>> variable_history;
    std::map<std::string, std::string> variable_dependencies;
    bool self_learning_enabled;
    double learning_rate;
    int update_counter;
    
    // Dynamic helper functions
    void updateDependencies(const std::string& changed_var);
    double computeGradient(const std::string& var, const std::string& target);
    void recordHistory(const std::string& name, double value);

public:
    // Constructor with stellar-enhanced dynamic capabilities
    SurfaceMagneticFieldModule();
    
    // Core magnetic field computations for stellar systems
    double computeB_j(double t, double B_s);
    double computeB_s_min();
    double computeB_s_max();
    double computeU_g3_example(double t, double B_s);
    double computeStellarCoupling(double B_field, double stellar_luminosity);
    
    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    
    // Advanced dynamic capabilities for stellar systems
    void autoCalibrate(const std::string& observable, double target_value, double tolerance = 0.01);
    void adaptiveUpdate(double dt, const std::string& feedback_param = "");
    void scaleToStellarData(const std::map<std::string, double>& stellar_data);
    void addCustomVariable(const std::string& name, double value, const std::string& dependency = "");
    std::map<std::string, double> getVariableHistory(const std::string& name, int steps = 10);
    void enableSelfLearning(bool enable);
    void exportState(const std::string& filename);
    void importState(const std::string& filename);
    
    // Enhanced magnetic field equations for stellar systems
    std::string getEquationText();
};

#endif // SURFACE_MAGNETIC_FIELD_MODULE_STELLAR_H

// ===== SURFACE MAGNETIC FIELD MODULE IMPLEMENTATION FOR STELLAR SYSTEMS =====

// Enhanced SurfaceMagneticFieldModule constructor with stellar capabilities
SurfaceMagneticFieldModule::SurfaceMagneticFieldModule() {
    // Initialize dynamic capabilities
    self_learning_enabled = false;
    learning_rate = 0.09;  // Higher learning rate for stellar variability
    update_counter = 0;
    
    // Universal constants
    variables["B_s_min"] = 5e-8;                    // T (quiet stellar systems)
    variables["B_s_max"] = 3.0;                     // T (active stellar systems)
    variables["B_ref"] = 3.0;                       // T (reference max for stars)
    variables["k_3"] = 3.2;                         // Enhanced coupling for stellar systems
    variables["omega_s"] = 2.1e-5;                  // rad/s (stellar system frequency)
    variables["P_core"] = 1.0;                      // Unitless
    variables["E_react"] = 7e46;                    // J (enhanced for stellar systems)
    variables["pi"] = 3.141592653589793;
    variables["t"] = 0.0;                           // s
    
    // Stellar-specific parameters
    variables["luminosity_coupling"] = 6.8e-17;    // Luminosity-magnetic coupling
    variables["stellar_rotation"] = 3.2e-6;        // rad/s (typical stellar rotation)
    variables["magnetic_diffusion"] = 6e-11;       // m/s (stellar diffusion)
    variables["convection_velocity"] = 5e3;        // m/s (stellar convection)
    variables["nuclear_burning_rate"] = 3.8e26;    // W (typical stellar power)
    
    // System evolution parameters for stellar systems
    variables["evolution_timescale"] = 3e14;       // s (stellar evolution timescale)
    variables["thermal_coupling"] = 1.8e-8;        // Thermal-magnetic coupling
    variables["flare_frequency"] = 0.1;            // day (stellar flare frequency)
    variables["chromospheric_heating"] = 2.5e4;    // K (chromospheric temperature)
    variables["coronal_heating"] = 1.5e6;          // K (coronal temperature)
    variables["stellar_wind_velocity"] = 4e5;      // m/s (stellar wind speed)
}

// Compute minimum surface magnetic field for stellar systems
double SurfaceMagneticFieldModule::computeB_s_min() {
    return variables["B_s_min"];
}

// Compute maximum surface magnetic field for stellar systems
double SurfaceMagneticFieldModule::computeB_s_max() {
    return variables["B_s_max"];
}

// Compute scaled B_j based on time t and surface field B_s for stellar systems
double SurfaceMagneticFieldModule::computeB_j(double t, double B_s) {
    variables["t"] = t;
    
    // Enhanced magnetic field evolution with stellar coupling
    double stellar_oscillation = 0.8 * std::sin(variables["omega_s"] * t);
    double rotation_factor = 1.0 + 0.3 * std::sin(variables["stellar_rotation"] * t);
    double nuclear_factor = 1.0 + 0.1 * std::sin(t / (365.25 * 86400));  // Annual nuclear cycle
    double flare_factor = 1.0 + variables["flare_frequency"] * std::exp(-std::fmod(t, 86400) / 7200);  // Daily flare cycle
    double thermal_factor = 1.0 + variables["thermal_coupling"] * variables["chromospheric_heating"] / 2.5e4;
    double wind_factor = 1.0 + 0.05 * variables["stellar_wind_velocity"] / 4e5;
    
    double base_b = variables["B_ref"] + stellar_oscillation * rotation_factor * nuclear_factor * flare_factor * thermal_factor * wind_factor;
    
    return base_b * (B_s / variables["B_ref"]);
}

// Compute stellar-magnetic field coupling
double SurfaceMagneticFieldModule::computeStellarCoupling(double B_field, double stellar_luminosity) {
    double coupling_strength = variables["luminosity_coupling"];
    double rotation_factor = variables["stellar_rotation"] / 3.2e-6;  // Normalized to typical stellar rotation
    double nuclear_factor = variables["nuclear_burning_rate"] / 3.8e26;  // Normalized to solar luminosity
    double wind_factor = variables["stellar_wind_velocity"] / 4e5;  // Normalized wind speed
    double thermal_gradient = variables["coronal_heating"] / variables["chromospheric_heating"];
    
    return coupling_strength * B_field * stellar_luminosity * rotation_factor * nuclear_factor * wind_factor * thermal_gradient;
}

// Compute U_g3 example with stellar enhancement
double SurfaceMagneticFieldModule::computeU_g3_example(double t, double B_s) {
    double k_3 = variables["k_3"];
    double b_j = computeB_j(t, B_s);
    double cos_term = std::cos(variables["omega_s"] * t * variables["pi"]);
    double p_core = variables["P_core"];
    double e_react = variables["E_react"];
    
    // Stellar enhancement factor
    double stellar_enhancement = 1.0 + computeStellarCoupling(b_j, variables["nuclear_burning_rate"]);
    
    // Stellar flare activity
    double flare_factor = 1.0 + variables["flare_frequency"] * std::sin(t / 3600);  // Hourly flare modulation
    
    // Coronal mass ejection effects
    double cme_factor = 1.0 + 0.02 * std::sin(t / (7 * 86400));  // Weekly CME cycle
    
    return k_3 * b_j * cos_term * p_core * e_react * stellar_enhancement * flare_factor * cme_factor;
}

// Enhanced updateVariable with stellar dynamic capabilities
void SurfaceMagneticFieldModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
    
    // Record history for dynamic capabilities
    recordHistory(name, value);
    
    // Update dependencies
    updateDependencies(name);
    
    // Increment update counter
    update_counter++;
    
    // Trigger self-learning if enabled
    if (self_learning_enabled && update_counter % 3 == 0) {  // High frequency for stellar variability
        adaptiveUpdate(1.0, name);
    }
}

// Auto-calibrate magnetic field parameters for stellar systems
void SurfaceMagneticFieldModule::autoCalibrate(const std::string& observable, double target_value, double tolerance) {
    if (variables.find(observable) == variables.end()) {
        std::cerr << "Observable '" << observable << "' not found for calibration." << std::endl;
        return;
    }
    
    double current_value = variables[observable];
    double error = std::abs(current_value - target_value) / target_value;
    
    if (error > tolerance) {
        // Stellar-specific parameter adjustment
        std::vector<std::string> tunable_params = {"B_ref", "k_3", "omega_s", "P_core", "luminosity_coupling", "thermal_coupling", "stellar_rotation", "flare_frequency"};
        
        for (const auto& param : tunable_params) {
            double gradient = computeGradient(param, observable);
            if (std::abs(gradient) > 1e-25) {
                double adjustment = learning_rate * (target_value - current_value) / gradient;
                variables[param] += adjustment;
                recordHistory(param, variables[param]);
            }
        }
        
        std::cout << "Auto-calibrated stellar magnetic " << observable << " from " << current_value 
                  << " to target " << target_value << " (error: " << error << ")" << std::endl;
    }
}

// Adaptive magnetic field evolution for stellar systems
void SurfaceMagneticFieldModule::adaptiveUpdate(double dt, const std::string& feedback_param) {
    if (!self_learning_enabled) return;
    
    // Stellar evolution timescale
    double evolution_factor = std::exp(-dt / variables["evolution_timescale"]);
    
    // Nuclear burning rate evolution
    variables["nuclear_burning_rate"] *= (1.0 + 0.0005 * std::sin(dt / 1e11));
    
    // Adaptive magnetic field reference with stellar coupling
    double stellar_factor = variables["nuclear_burning_rate"] / 3.8e26;
    variables["B_ref"] = variables["B_s_max"] * (0.5 + 0.5 * stellar_factor);
    
    // Enhanced magnetic diffusion effects for stellar systems
    double diffusion_decay = std::exp(-dt * variables["magnetic_diffusion"] / 1e4);
    variables["k_3"] *= diffusion_decay;
    
    // Stellar rotation effects (magnetic braking)
    double rotation_decay = std::exp(-dt / 1e13);  // Magnetic braking timescale
    variables["stellar_rotation"] *= rotation_decay;
    variables["omega_s"] *= (1.0 + 0.15 * variables["convection_velocity"] / 5e3);
    
    // Flare frequency evolution
    variables["flare_frequency"] *= (1.0 + 0.003 * std::cos(dt / 1e8));
    
    // Stellar wind evolution
    variables["stellar_wind_velocity"] *= (1.0 + 0.001 * std::sin(dt / 1e10));
    
    // Chromospheric and coronal heating evolution
    variables["chromospheric_heating"] *= (1.0 + 0.0002 * std::sin(dt / 1e9));
    variables["coronal_heating"] *= (1.0 + 0.0003 * std::cos(dt / 1e9));
    
    recordHistory("adaptive_time", dt);
    std::cout << "Stellar magnetic adaptive update: B_ref=" << variables["B_ref"] 
              << ", nuclear_rate=" << variables["nuclear_burning_rate"] << std::endl;
}

// Scale to stellar observational data
void SurfaceMagneticFieldModule::scaleToStellarData(const std::map<std::string, double>& stellar_data) {
    for (const auto& data : stellar_data) {
        if (data.first == "nuclear_burning_rate") {
            double scaling = data.second / variables["nuclear_burning_rate"];
            variables["nuclear_burning_rate"] = data.second;
            variables["thermal_coupling"] *= scaling;
            variables["luminosity_coupling"] *= std::sqrt(scaling);
        }
        
        if (data.first == "stellar_rotation") {
            variables["stellar_rotation"] = data.second;
            variables["omega_s"] *= data.second / 3.2e-6;
            variables["flare_frequency"] *= std::sqrt(data.second / 3.2e-6);  // Rotation-activity relation
        }
        
        if (data.first == "magnetic_field_strength") {
            double scaling = data.second / variables["B_s_max"];
            variables["B_s_max"] = data.second;
            variables["B_ref"] *= scaling;
        }
        
        if (data.first == "stellar_wind_velocity") {
            variables["stellar_wind_velocity"] = data.second;
            variables["coronal_heating"] *= std::pow(data.second / 4e5, 0.5);  // Wind-temperature relation
        }
        
        if (data.first == "flare_frequency") {
            variables["flare_frequency"] = data.second;
            variables["chromospheric_heating"] *= (1.0 + 0.1 * data.second / 0.1);
        }
    }
    std::cout << "Scaled stellar magnetic module to " << stellar_data.size() << " stellar observations." << std::endl;
}

// Add custom magnetic variables for stellar systems
void SurfaceMagneticFieldModule::addCustomVariable(const std::string& name, double value, const std::string& dependency) {
    variables[name] = value;
    if (!dependency.empty()) {
        variable_dependencies[name] = dependency;
    }
    recordHistory(name, value);
    std::cout << "Added custom stellar magnetic variable: " << name << " = " << value << std::endl;
}

// Get magnetic parameter history for stellar systems
std::map<std::string, double> SurfaceMagneticFieldModule::getVariableHistory(const std::string& name, int steps) {
    std::map<std::string, double> history;
    if (variable_history.find(name) != variable_history.end()) {
        auto& hist = variable_history[name];
        int start = std::max(0, (int)hist.size() - steps);
        for (int i = start; i < (int)hist.size(); i++) {
            history["step_" + std::to_string(i)] = hist[i];
        }
    }
    return history;
}

// Enable stellar magnetic self-learning
void SurfaceMagneticFieldModule::enableSelfLearning(bool enable) {
    self_learning_enabled = enable;
    if (enable) {
        std::cout << "Stellar magnetic self-learning enabled with rate: " << learning_rate << std::endl;
    } else {
        std::cout << "Stellar magnetic self-learning disabled." << std::endl;
    }
}

// Export stellar magnetic state
void SurfaceMagneticFieldModule::exportState(const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        file << "# SurfaceMagneticFieldModule Stellar State Export" << std::endl;
        file << "update_counter=" << update_counter << std::endl;
        file << "learning_rate=" << learning_rate << std::endl;
        file << "self_learning_enabled=" << (self_learning_enabled ? 1 : 0) << std::endl;
        
        for (const auto& var : variables) {
            file << var.first << "=" << var.second << std::endl;
        }
        file.close();
        std::cout << "Stellar magnetic state exported to: " << filename << std::endl;
    }
}

// Import stellar magnetic state
void SurfaceMagneticFieldModule::importState(const std::string& filename) {
    std::ifstream file(filename);
    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            if (line[0] == '#') continue;
            
            size_t eq_pos = line.find('=');
            if (eq_pos != std::string::npos) {
                std::string key = line.substr(0, eq_pos);
                std::string value_str = line.substr(eq_pos + 1);
                
                if (key == "update_counter") {
                    update_counter = std::stoi(value_str);
                } else if (key == "learning_rate") {
                    learning_rate = std::stod(value_str);
                } else if (key == "self_learning_enabled") {
                    self_learning_enabled = (std::stoi(value_str) == 1);
                } else {
                    variables[key] = std::stod(value_str);
                }
            }
        }
        file.close();
        std::cout << "Stellar magnetic state imported from: " << filename << std::endl;
    }
}

// Enhanced equation text for stellar systems
std::string SurfaceMagneticFieldModule::getEquationText() {
    return "Stellar-Enhanced Magnetic Field Equations:\n"
           "B_j  (B_ref + 0.8 sin(?_s t) * O_rot * N_nuclear * F_flare * T_thermal * W_wind) * (B_s / B_ref) T\n"
           "U_g3 = k_3 *  B_j * cos(?_s t p) * P_core * E_react * (1 + Stellar_coupling) * F_flare * CME_factor\n"
           "Stellar_coupling = ?_stellar * B_field * L_nuclear * (O/O_0) * (N/N_0) * (v_wind/v_0) * (T_cor/T_chr)\n"
           "Where:\n"
           "- B_s = [5e-8, 3.0] T (stellar system range)\n"
           "- L_nuclear = " + std::to_string(variables["nuclear_burning_rate"]) + " W (nuclear burning rate)\n"
           "- O_rot = " + std::to_string(variables["stellar_rotation"]) + " rad/s (stellar rotation)\n"
           "- ?_stellar = " + std::to_string(variables["luminosity_coupling"]) + " (stellar coupling)\n"
           "- F_flare = " + std::to_string(variables["flare_frequency"]) + " day (flare frequency)\n"
           "- v_wind = " + std::to_string(variables["stellar_wind_velocity"]) + " m/s (stellar wind velocity)\n"
           "Stellar Systems: M74, Eagle Nebula (M16), M84, Centaurus A,\n"
           "Supernova Survey\n"
           "Enhanced Features: Nuclear burning coupling, stellar rotation effects,\n"
           "flare activity modulation, coronal heating, stellar wind coupling.";
}

// Helper functions for stellar magnetic field module
void SurfaceMagneticFieldModule::updateDependencies(const std::string& changed_var) {
    if (changed_var == "nuclear_burning_rate") {
        // Update thermal coupling based on nuclear burning rate
        variables["thermal_coupling"] = 1.8e-8 * (variables["nuclear_burning_rate"] / 3.8e26);
        // Update coronal heating from nuclear rate
        variables["coronal_heating"] = 1.5e6 * std::sqrt(variables["nuclear_burning_rate"] / 3.8e26);
    }
    
    if (changed_var == "B_s_max") {
        variables["B_ref"] = variables["B_s_max"];
    }
    
    if (changed_var == "stellar_rotation") {
        // Update omega_s based on stellar rotation
        variables["omega_s"] = 2.1e-5 * (variables["stellar_rotation"] / 3.2e-6);
        // Update flare frequency from rotation-activity relation
        variables["flare_frequency"] = 0.1 * std::sqrt(variables["stellar_rotation"] / 3.2e-6);
        // Update stellar wind velocity
        variables["stellar_wind_velocity"] = 4e5 * std::pow(variables["stellar_rotation"] / 3.2e-6, 0.3);
    }
    
    if (changed_var == "omega_s") {
        // Update convection velocity
        variables["convection_velocity"] = 5e3 * (variables["omega_s"] / 2.1e-5);
    }
    
    if (changed_var == "flare_frequency") {
        // Update chromospheric heating from flare activity
        variables["chromospheric_heating"] = 2.5e4 * (1.0 + variables["flare_frequency"] / 0.1);
    }
}

double SurfaceMagneticFieldModule::computeGradient(const std::string& var, const std::string& target) {
    if (variables.find(var) == variables.end() || variables.find(target) == variables.end()) {
        return 0.0;
    }
    
    double original_value = variables[var];
    double original_target = variables[target];
    
    // Small perturbation
    double delta = original_value * 1e-6;
    variables[var] += delta;
    
    // Recompute target with stellar coupling
    double new_target = computeB_j(variables["t"], variables["B_ref"]);
    
    // Restore original value
    variables[var] = original_value;
    
    return (new_target - original_target) / delta;
}

void SurfaceMagneticFieldModule::recordHistory(const std::string& name, double value) {
    variable_history[name].push_back(value);
    
    // Keep only last 80 values for stellar systems (shorter history for fast variability)
    if (variable_history[name].size() > 80) {
        variable_history[name].erase(variable_history[name].begin());
    }
}
