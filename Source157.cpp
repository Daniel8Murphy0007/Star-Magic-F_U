// UQFFBuoyancyModule157.h
// UQFF Buoyancy for Observational Systems: M104, NGC 4839, Chandra/Webb, NGC 346, NGC 1672
// Copyright - Daniel T. Murphy, analyzed Oct 22, 2025.

#ifndef UQFF_BUOYANCY_MODULE_H
#define UQFF_BUOYANCY_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>


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

class UQFFBuoyancyModule157 {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
    std::map<std::string, double> variables;
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;



public:
    UQFFBuoyancyModule157();
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);
    double computeUb1(const std::string& system);
    double computeUi(const std::string& system);
    double computeFBi(const std::string& system);
    double computeDPM_resonance(const std::string& system);
    double computeX2(const std::string& system);
    double computeG(const std::string& system);
    double computeQ_wave(const std::string& system);
    std::string getEquationText(const std::string& system);
    void printVariables();
};

#endif // UQFF_BUOYANCY_MODULE_H

// Implementation
UQFFBuoyancyModule157::UQFFBuoyancyModule157() {
        enableDynamicTerms = true;
        enableLogging = false;
        learningRate = 0.001;
        metadata["enhanced"] = "true";
        metadata["version"] = "2.0-Enhanced";

    variables["c"] = 299792458.0;
    variables["h"] = 6.62607015e-34;
    variables["k_B"] = 1.380649e-23;
    variables["G"] = 6.67430e-11;
    variables["M_s"] = 1.989e30;
    variables["R_s"] = 6.96e8;
    variables["M104_mass"] = 8e11 * 1.989e30;
    variables["M104_dist"] = 9.55 * 3.086e22;
    variables["M104_factor"] = 1.25;
    variables["NGC4839_mass"] = 2.5e12 * 1.989e30;
    variables["NGC4839_dist"] = 97.7 * 3.086e22;
    variables["NGC4839_factor"] = 1.45;
    variables["CW_factor"] = 1.85;
    variables["CW_efficiency"] = 0.95;
    variables["NGC346_mass"] = 1.8e5 * 1.989e30;
    variables["NGC346_dist"] = 60 * 1000 * 3.086e22;
    variables["NGC346_factor"] = 2.1;
    variables["NGC1672_mass"] = 1.2e11 * 1.989e30;
    variables["NGC1672_dist"] = 18.6 * 3.086e22;
    variables["NGC1672_factor"] = 1.65;
    variables["t_obs"] = 0.0;
    variables["omega_obs"] = 2 * M_PI / (365.25 * 24 * 3600);
    variables["rho_vac"] = 9.47e-27;
    variables["beta_1"] = 1.2e-3;
    variables["beta_i"] = 0.85e-3;
    variables["alpha"] = 0.0073;
    variables["tau_decay"] = 887.7;
}

double UQFFBuoyancyModule157::computeUb1(const std::string& system) {
    double base_Ub1 = variables["beta_1"] * variables["rho_vac"] * std::pow(variables["c"], 2);
    double system_factor = 1.0;
    double time_modulation = 1.0 + 0.15 * std::sin(variables["omega_obs"] * variables["t_obs"]);
    
    if (system == "M104") {
        system_factor = variables["M104_factor"];
        double distance_factor = std::pow(variables["M104_dist"] / (10 * 3.086e22), -0.3);
        base_Ub1 *= distance_factor;
    } else if (system == "NGC4839") {
        system_factor = variables["NGC4839_factor"];
        double cluster_factor = 1.0 + 0.25 * std::log(variables["NGC4839_mass"] / (1e12 * variables["M_s"]));
        base_Ub1 *= cluster_factor;
    } else if (system == "Chandra_Webb") {
        system_factor = variables["CW_factor"];
        base_Ub1 *= variables["CW_efficiency"];
    } else if (system == "NGC346") {
        system_factor = variables["NGC346_factor"];
        double sf_enhancement = 1.0 + 0.8 * std::tanh(variables["NGC346_mass"] / (1e5 * variables["M_s"]));
        base_Ub1 *= sf_enhancement;
    } else if (system == "NGC1672") {
        system_factor = variables["NGC1672_factor"];
        double spiral_factor = 1.0 + 0.35 * std::cos(2 * variables["omega_obs"] * variables["t_obs"]);
        base_Ub1 *= spiral_factor;
    }
    
    return base_Ub1 * system_factor * time_modulation;
}

double UQFFBuoyancyModule157::computeUi(const std::string& system) {
    double base_Ui = variables["beta_i"] * variables["rho_vac"] * std::pow(variables["c"], 2);
    double system_factor = 1.0;
    if (system == "M104") system_factor = variables["M104_factor"] * 0.95;
    else if (system == "NGC4839") system_factor = variables["NGC4839_factor"] * 1.1;
    else if (system == "Chandra_Webb") system_factor = variables["CW_factor"] * 1.25;
    else if (system == "NGC346") system_factor = variables["NGC346_factor"] * 0.85;
    else if (system == "NGC1672") system_factor = variables["NGC1672_factor"] * 1.05;
    return base_Ui * system_factor;
}

double UQFFBuoyancyModule157::computeFBi(const std::string& system) {
    double Ub1 = computeUb1(system);
    double Ui = computeUi(system);
    double coupling_strength = variables["alpha"] * std::sqrt(Ub1 * Ui);
    double system_coupling = 1.0;
    if (system == "M104") system_coupling = 1.0 + 0.15 * std::log(variables["M104_mass"] / variables["M_s"]);
    else if (system == "NGC4839") system_coupling = 1.0 + 0.25 * std::log(variables["NGC4839_mass"] / variables["M_s"]);
    else if (system == "Chandra_Webb") system_coupling = 1.0 + 0.4 * variables["CW_efficiency"];
    else if (system == "NGC346") system_coupling = 1.0 + 0.3 * std::log(variables["NGC346_mass"] / variables["M_s"]);
    else if (system == "NGC1672") system_coupling = 1.0 + 0.2 * std::log(variables["NGC1672_mass"] / variables["M_s"]);
    return coupling_strength * system_coupling;
}

double UQFFBuoyancyModule157::computeDPM_resonance(const std::string& system) {
    double base_frequency = variables["c"] / (2 * M_PI * variables["R_s"]);
    double system_resonance = 1.0;
    if (system == "M104") system_resonance = 1.0 + 0.2 * std::sqrt(variables["M104_mass"] / variables["M_s"]);
    else if (system == "NGC4839") system_resonance = 1.0 + 0.35 * std::log(variables["NGC4839_mass"] / (1e12 * variables["M_s"]));
    else if (system == "Chandra_Webb") system_resonance = 1.0 + 0.5 * variables["CW_efficiency"];
    else if (system == "NGC346") system_resonance = 1.0 + 0.8 * std::tanh(variables["NGC346_mass"] / (1e5 * variables["M_s"]));
    else if (system == "NGC1672") system_resonance = 1.0 + 0.3 * std::sin(4 * variables["omega_obs"] * variables["t_obs"]);
    return base_frequency * system_resonance;
}

double UQFFBuoyancyModule157::computeX2(const std::string& system) {
    double base_X2 = std::pow(variables["h"] * variables["c"] / (variables["k_B"] * 2.7), 2);
    double system_enhancement = 1.0;
    if (system == "M104") system_enhancement = variables["M104_factor"] * 0.9;
    else if (system == "NGC4839") system_enhancement = variables["NGC4839_factor"] * 1.2;
    else if (system == "Chandra_Webb") system_enhancement = variables["CW_factor"] * variables["CW_efficiency"];
    else if (system == "NGC346") system_enhancement = variables["NGC346_factor"] * 1.1;
    else if (system == "NGC1672") system_enhancement = variables["NGC1672_factor"] * 1.05;
    return base_X2 * system_enhancement;
}

double UQFFBuoyancyModule157::computeG(const std::string& system) {
    double modified_G = variables["G"];
    double system_modification = 1.0;
    if (system == "M104") system_modification = 1.0 + 1e-6 * std::log(variables["M104_mass"] / (1e11 * variables["M_s"]));
    else if (system == "NGC4839") system_modification = 1.0 + 2e-6 * std::log(variables["NGC4839_mass"] / (1e12 * variables["M_s"]));
    else if (system == "Chandra_Webb") system_modification = 1.0 + 1.5e-6 * variables["CW_efficiency"];
    else if (system == "NGC346") system_modification = 1.0 + 3e-6 * std::tanh(variables["NGC346_mass"] / (1e5 * variables["M_s"]));
    else if (system == "NGC1672") system_modification = 1.0 + 1.2e-6 * std::log(variables["NGC1672_mass"] / (1e11 * variables["M_s"]));
    return modified_G * system_modification;
}

double UQFFBuoyancyModule157::computeQ_wave(const std::string& system) {
    double base_Q = variables["rho_vac"] * std::pow(variables["c"], 3) / variables["h"];
    double system_wave_factor = 1.0;
    if (system == "M104") {
        double distance_modulation = std::sin(variables["c"] * variables["t_obs"] / variables["M104_dist"]);
        system_wave_factor = variables["M104_factor"] * (1.0 + 0.1 * distance_modulation);
    } else if (system == "NGC4839") {
        double cluster_wave = std::cos(variables["c"] * variables["t_obs"] / variables["NGC4839_dist"]);
        system_wave_factor = variables["NGC4839_factor"] * (1.0 + 0.15 * cluster_wave);
    } else if (system == "Chandra_Webb") {
        system_wave_factor = variables["CW_factor"] * (1.0 + 0.2 * variables["CW_efficiency"]);
    } else if (system == "NGC346") {
        system_wave_factor = variables["NGC346_factor"] * (1.0 + 0.3);
    } else if (system == "NGC1672") {
        system_wave_factor = variables["NGC1672_factor"] * (1.0 + 0.25);
    }
    return base_Q * system_wave_factor;
}

void UQFFBuoyancyModule157::updateVariable(const std::string& name, double value) { variables[name] = value; }
void UQFFBuoyancyModule157::addToVariable(const std::string& name, double delta) { variables[name] += delta; }
void UQFFBuoyancyModule157::subtractFromVariable(const std::string& name, double delta) { variables[name] -= delta; }

std::string UQFFBuoyancyModule157::getEquationText(const std::string& system) {
    return "UQFF Buoyancy Module 157 - System: " + system + " - Ub1, Ui, FBi with observational enhancements";
}

void UQFFBuoyancyModule157::printVariables() {
    std::cout << "=== UQFF Buoyancy Module 157 Variables ===" << std::endl;
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << pair.second << std::endl;
    }
    std::cout << "==========================================" << std::endl;
}

// ===== SURFACE MAGNETIC FIELD MODULE FOR OBSERVATIONAL SYSTEMS =====

#ifndef SURFACE_MAGNETIC_FIELD_MODULE_OBS_H
#define SURFACE_MAGNETIC_FIELD_MODULE_OBS_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

class SurfaceMagneticFieldModule {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
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
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;



public:
    // Constructor with observational-enhanced dynamic capabilities
    SurfaceMagneticFieldModule();
    
    // Core magnetic field computations for observational systems
    double computeB_j(double t, double B_s);
    double computeB_s_min();
    double computeB_s_max();
    double computeU_g3_example(double t, double B_s);
    double computeObservationalCoupling(double B_field, double luminosity);
    
    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    
    // Advanced dynamic capabilities for observational systems
    void autoCalibrate(const std::string& observable, double target_value, double tolerance = 0.01);
    void adaptiveUpdate(double dt, const std::string& feedback_param = "");
    void scaleToObservationalData(const std::map<std::string, double>& obs_data);
    void addCustomVariable(const std::string& name, double value, const std::string& dependency = "");
    std::map<std::string, double> getVariableHistory(const std::string& name, int steps = 10);
    void enableSelfLearning(bool enable);
    void exportState(const std::string& filename);
    void importState(const std::string& filename);
    
    // Enhanced magnetic field equations for observational systems
    std::string getEquationText();
};

#endif // SURFACE_MAGNETIC_FIELD_MODULE_OBS_H

// ===== SURFACE MAGNETIC FIELD MODULE IMPLEMENTATION FOR OBSERVATIONAL SYSTEMS =====

// Enhanced SurfaceMagneticFieldModule constructor with observational capabilities
SurfaceMagneticFieldModule::SurfaceMagneticFieldModule() {
    // Initialize dynamic capabilities
    self_learning_enabled = false;
    learning_rate = 0.06;  // Moderate learning rate for observational stability
    update_counter = 0;
    
    // Universal constants
    variables["B_s_min"] = 5e-7;                    // T (quiet observational systems)
    variables["B_s_max"] = 1.5;                     // T (active observational systems)
    variables["B_ref"] = 1.5;                       // T (reference max for observations)
    variables["k_3"] = 2.2;                         // Standard coupling for observations
    variables["omega_s"] = 1.6e-5;                  // rad/s (observational system frequency)
    variables["P_core"] = 1.0;                      // Unitless
    variables["E_react"] = 4.5e46;                  // J (standard for observational systems)
    variables["pi"] = 3.141592653589793;
    variables["t"] = 0.0;                           // s
    
    // Observational-specific parameters
    variables["luminosity_coupling"] = 3.5e-15;    // Luminosity-magnetic coupling
    variables["spectral_index"] = -2.1;            // Typical radio spectral index
    variables["magnetic_diffusion"] = 3e-11;       // m/s (standard diffusion)
    variables["convection_velocity"] = 2e3;        // m/s (moderate convection)
    variables["observation_frequency"] = 1.4e9;    // Hz (typical radio frequency)
    
    // System evolution parameters for observations
    variables["evolution_timescale"] = 8e14;       // s (observational evolution timescale)
    variables["thermal_coupling"] = 1.2e-8;        // Thermal-magnetic coupling
    variables["flux_density"] = 1e-26;             // W m Hz (typical flux density)
    variables["angular_resolution"] = 1e-3;        // arcsec (typical resolution)
}

// Compute minimum surface magnetic field for observational systems
double SurfaceMagneticFieldModule::computeB_s_min() {
    return variables["B_s_min"];
}

// Compute maximum surface magnetic field for observational systems
double SurfaceMagneticFieldModule::computeB_s_max() {
    return variables["B_s_max"];
}

// Compute scaled B_j based on time t and surface field B_s for observational systems
double SurfaceMagneticFieldModule::computeB_j(double t, double B_s) {
    variables["t"] = t;
    
    // Enhanced magnetic field evolution with observational coupling
    double obs_oscillation = 0.5 * std::sin(variables["omega_s"] * t);
    double thermal_factor = 1.0 + variables["thermal_coupling"] * std::pow(variables["flux_density"] / 1e-26, 0.3);
    double base_b = variables["B_ref"] + obs_oscillation * thermal_factor;
    
    return base_b * (B_s / variables["B_ref"]);
}

// Compute observational-magnetic field coupling
double SurfaceMagneticFieldModule::computeObservationalCoupling(double B_field, double luminosity) {
    double coupling_strength = variables["luminosity_coupling"];
    double spectral_factor = std::pow(variables["observation_frequency"] / 1.4e9, variables["spectral_index"]);
    double flux_factor = variables["flux_density"] / 1e-26;
    
    return coupling_strength * B_field * luminosity * spectral_factor * flux_factor;
}

// Compute U_g3 example with observational enhancement
double SurfaceMagneticFieldModule::computeU_g3_example(double t, double B_s) {
    double k_3 = variables["k_3"];
    double b_j = computeB_j(t, B_s);
    double cos_term = std::cos(variables["omega_s"] * t * variables["pi"]);
    double p_core = variables["P_core"];
    double e_react = variables["E_react"];
    
    // Observational enhancement factor
    double obs_enhancement = 1.0 + computeObservationalCoupling(b_j, variables["flux_density"] * 1e26);
    
    return k_3 * b_j * cos_term * p_core * e_react * obs_enhancement;
}

// Enhanced updateVariable with observational dynamic capabilities
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
    if (self_learning_enabled && update_counter % 5 == 0) {  // Standard frequency for observations
        adaptiveUpdate(1.0, name);
    }
}

// Auto-calibrate magnetic field parameters for observational systems
void SurfaceMagneticFieldModule::autoCalibrate(const std::string& observable, double target_value, double tolerance) {
    if (variables.find(observable) == variables.end()) {
        std::cerr << "Observable '" << observable << "' not found for calibration." << std::endl;
        return;
    }
    
    double current_value = variables[observable];
    double error = std::abs(current_value - target_value) / target_value;
    
    if (error > tolerance) {
        // Observational-specific parameter adjustment
        std::vector<std::string> tunable_params = {"B_ref", "k_3", "omega_s", "P_core", "luminosity_coupling", "thermal_coupling"};
        
        for (const auto& param : tunable_params) {
            double gradient = computeGradient(param, observable);
            if (std::abs(gradient) > 1e-25) {
                double adjustment = learning_rate * (target_value - current_value) / gradient;
                variables[param] += adjustment;
                recordHistory(param, variables[param]);
            }
        }
        
        std::cout << "Auto-calibrated observational magnetic " << observable << " from " << current_value 
                  << " to target " << target_value << " (error: " << error << ")" << std::endl;
    }
}

// Adaptive magnetic field evolution for observational systems
void SurfaceMagneticFieldModule::adaptiveUpdate(double dt, const std::string& feedback_param) {
    if (!self_learning_enabled) return;
    
    // Observational evolution timescale
    double evolution_factor = std::exp(-dt / variables["evolution_timescale"]);
    
    // Flux density evolution
    variables["flux_density"] *= (1.0 + 0.0002 * std::sin(dt / 1e11));
    
    // Adaptive magnetic field reference with observational coupling
    double obs_factor = variables["flux_density"] / 1e-26;
    variables["B_ref"] = variables["B_s_max"] * (0.3 + 0.7 * std::sqrt(obs_factor));
    
    // Standard magnetic diffusion effects for observational systems
    double diffusion_decay = std::exp(-dt * variables["magnetic_diffusion"] / 1e4);
    variables["k_3"] *= diffusion_decay;
    
    // Observational-driven frequency evolution
    double frequency_enhancement = 1.0 + 0.1 * variables["convection_velocity"] / 2e3;
    variables["omega_s"] *= frequency_enhancement;
    
    // Spectral index variability
    variables["spectral_index"] *= (1.0 + 0.0005 * std::cos(dt / 1e10));
    
    recordHistory("adaptive_time", dt);
    std::cout << "Observational magnetic adaptive update: B_ref=" << variables["B_ref"] 
 * Enhanced: November 04, 2025 - Added self-expanding capabilities
              << ", flux_density=" << variables["flux_density"] << std::endl;
}

// Scale to observational data
void SurfaceMagneticFieldModule::scaleToObservationalData(const std::map<std::string, double>& obs_data) {
    for (const auto& data : obs_data) {
        if (data.first == "flux_density") {
            double scaling = data.second / variables["flux_density"];
            variables["flux_density"] = data.second;
            variables["luminosity_coupling"] *= std::sqrt(scaling);
        }
        
        if (data.first == "observation_frequency") {
            variables["observation_frequency"] = data.second;
            double freq_factor = data.second / 1.4e9;
            variables["spectral_index"] *= std::pow(freq_factor, 0.1);
        }
        
        if (data.first == "magnetic_field_strength") {
            double scaling = data.second / variables["B_s_max"];
            variables["B_s_max"] = data.second;
            variables["B_ref"] *= scaling;
        }
        
        if (data.first == "angular_resolution") {
            variables["angular_resolution"] = data.second;
            variables["thermal_coupling"] *= std::sqrt(1e-3 / data.second);
        }
    }
    std::cout << "Scaled observational magnetic module to " << obs_data.size() << " observations." << std::endl;
}

// Add custom magnetic variables for observational systems
void SurfaceMagneticFieldModule::addCustomVariable(const std::string& name, double value, const std::string& dependency) {
    variables[name] = value;
    if (!dependency.empty()) {
        variable_dependencies[name] = dependency;
    }
    recordHistory(name, value);
    std::cout << "Added custom observational magnetic variable: " << name << " = " << value << std::endl;
}

// Get magnetic parameter history for observational systems
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

// Enable observational magnetic self-learning
void SurfaceMagneticFieldModule::enableSelfLearning(bool enable) {
    self_learning_enabled = enable;
    if (enable) {
        std::cout << "Observational magnetic self-learning enabled with rate: " << learning_rate << std::endl;
    } else {
        std::cout << "Observational magnetic self-learning disabled." << std::endl;
    }
}

// Export observational magnetic state
void SurfaceMagneticFieldModule::exportState(const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        file << "# SurfaceMagneticFieldModule Observational State Export" << std::endl;
        file << "update_counter=" << update_counter << std::endl;
        file << "learning_rate=" << learning_rate << std::endl;
        file << "self_learning_enabled=" << (self_learning_enabled ? 1 : 0) << std::endl;
        
        for (const auto& var : variables) {
            file << var.first << "=" << var.second << std::endl;
        }
        file.close();
        std::cout << "Observational magnetic state exported to: " << filename << std::endl;
    }
}

// Import observational magnetic state
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
        std::cout << "Observational magnetic state imported from: " << filename << std::endl;
    }
}

// Enhanced equation text for observational systems
std::string SurfaceMagneticFieldModule::getEquationText() {
    return "Observational-Enhanced Magnetic Field Equations:\n"
           "B_j  (B_ref + 0.5 sin(ω_s t) * T_thermal) * (B_s / B_ref) T\n"
           "U_g3 = k_3 *  B_j * cos(ω_s t π) * P_core * E_react * (1 + Obs_coupling)\n"
           "Obs_coupling = κ_obs * B_field * L * (ν/ν_0)^α * (S/S_0)\n"
           "Where:\n"
           "- B_s = [5e-7, 1.5] T (observational system range)\n"
           "- S = " + std::to_string(variables["flux_density"]) + " W m Hz (flux density)\n"
           "- ν = " + std::to_string(variables["observation_frequency"]) + " Hz (observation frequency)\n"
           "- α = " + std::to_string(variables["spectral_index"]) + " (spectral index)\n"
           "- κ_obs = " + std::to_string(variables["luminosity_coupling"]) + " (observational coupling)\n"
           "Observational Systems: M104 (Sombrero), NGC 4839, Chandra/Webb,\n"
           "NGC 346, NGC 1672\n"
           "Enhanced Features: Flux density coupling, spectral index evolution,\n"
           "frequency-dependent scaling, angular resolution effects.";
}

// Helper functions for observational magnetic field module
void SurfaceMagneticFieldModule::updateDependencies(const std::string& changed_var) {
    if (changed_var == "flux_density") {
        // Update thermal coupling based on flux density
        double flux_factor = variables["flux_density"] / 1e-26;
        variables["thermal_coupling"] = 1.2e-8 * std::sqrt(flux_factor);
    }
    
    if (changed_var == "B_s_max") {
        variables["B_ref"] = variables["B_s_max"];
    }
    
    if (changed_var == "observation_frequency") {
        // Update luminosity coupling based on frequency
        double freq_factor = variables["observation_frequency"] / 1.4e9;
        variables["luminosity_coupling"] = 3.5e-15 * std::pow(freq_factor, 0.2);
    }
    
    if (changed_var == "omega_s") {
        // Update evolution timescale
        variables["evolution_timescale"] = 2 * M_PI / variables["omega_s"] * 1e9;
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
    
    // Recompute target with observational coupling
    double new_target = computeB_j(variables["t"], variables["B_ref"]);
    
    // Restore original value
    variables[var] = original_value;
    
    return (new_target - original_target) / delta;
}

void SurfaceMagneticFieldModule::recordHistory(const std::string& name, double value) {
    variable_history[name].push_back(value);
    
    // Keep only last 100 values for observational systems (standard history)
    if (variable_history[name].size() > 100) {
        variable_history[name].erase(variable_history[name].begin());
    }
}
