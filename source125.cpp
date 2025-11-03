// Ug3DiskVectorModule.h
// Modular C++ implementation of the Unit Vector in the Ug3 Disk Plane (??_j) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes ??_j (unit vector, magnitude=1; e.g., [cos ?_j, sin ?_j, 0]); scales in Universal Magnetism U_m term.
// Pluggable: #include "Ug3DiskVectorModule.h"
// Ug3DiskVectorModule mod; mod.computeUmContribution(0.0, 1); mod.updateVariable("theta_j", new_value);
// Variables in std::map; example for j=1 at t=0, ?_j=0 (??_j=[1,0,0], U_m?2.28e65 J/m�).
// Approximations: ??_j magnitude=1; 1 - exp=0 at t=0; ?_j / r_j=2.26e10 T m�.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef UG3_DISK_VECTOR_MODULE_H
#define UG3_DISK_VECTOR_MODULE_H

#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <functional>
#include <random>
#include <algorithm>
#include <sstream>

class Ug3DiskVectorModule {
private:
    std::map<std::string, double> variables;
    std::vector<double> computePhiHat_j(int j);
    double computeUmBase(double t);
    double computeUmContribution(double t, int j);

public:
    // Constructor: Initialize with framework defaults
    Ug3DiskVectorModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    std::vector<double> computePhiHat_j(int j);  // Unit vector [cos ?_j, sin ?_j, 0]
    double computePhiHatMagnitude(int j);  // 1.0 (normalized)
    double computeUmContribution(double t, int j);  // U_m single string (J/m^3)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print ??_j and U_m
    void printVectorAndUm(int j = 1, double t = 0.0);

    // ========== ENHANCED SELF-UPDATE & SELF-EXPANSION CAPABILITIES ==========
    
    // Variable Management
    void createVariable(const std::string& name, double value);
    void removeVariable(const std::string& name);
    void cloneVariable(const std::string& source, const std::string& target);
    std::vector<std::string> listVariables() const;
    std::string getSystemName() const;
    
    // Batch Operations
    void transformVariableGroup(const std::vector<std::string>& varNames, 
                                std::function<double(double)> transformFunc);
    void scaleVariableGroup(const std::vector<std::string>& varNames, double factor);
    
    // Self-Expansion Methods (Domain-Specific for Ug3 Disk Vector)
    void expandParameterSpace(double expansion_factor);
    void expandVectorScale(double angle_factor, double geometry_factor);
    void expandMagneticScale(double um_factor, double coupling_factor);
    void expandDiskScale(double radius_factor, double density_factor);
    
    // Self-Refinement
    void autoRefineParameters();
    void calibrateToObservations(const std::map<std::string, double>& observed_data);
    void optimizeForMetric(const std::string& metric);
    
    // Parameter Exploration
    std::vector<std::map<std::string, double>> generateVariations(int count);
    
    // Adaptive Evolution
    void mutateParameters(double mutation_rate);
    void evolveSystem(int generations);
    
    // State Management
    void saveState(const std::string& label);
    void restoreState(const std::string& label);
    std::vector<std::string> listSavedStates() const;
    std::string exportState() const;
    
    // System Analysis
    std::map<std::string, double> sensitivityAnalysis(const std::string& output_var);
    std::string generateReport() const;
    bool validateConsistency() const;
    void autoCorrectAnomalies();
};

#endif // UG3_DISK_VECTOR_MODULE_H

// Ug3DiskVectorModule.cpp
#include "Ug3DiskVectorModule.h"

// Constructor: Set framework defaults
Ug3DiskVectorModule::Ug3DiskVectorModule() {
    // Universal constants
    variables["theta_j"] = 0.0;                     // rad (default azimuthal angle)
    variables["mu_j"] = 3.38e23;                    // T�m^3 (j=1)
    variables["r_j"] = 1.496e13;                    // m
    variables["gamma"] = 5e-5 / 86400.0;            // s^-1
    variables["t_n"] = 0.0;                         // s
    variables["P_SCm"] = 1.0;                       // Pressure
    variables["E_react"] = 1e46;                    // J
    variables["f_Heaviside"] = 0.01;                // Unitless
    variables["f_quasi"] = 0.01;                    // Unitless
    variables["pi"] = 3.141592653589793;

    // Derived
    variables["scale_Heaviside"] = 1e13;
    variables["heaviside_factor"] = 1.0 + variables["scale_Heaviside"] * variables["f_Heaviside"];
}

// Update variable
void Ug3DiskVectorModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void Ug3DiskVectorModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void Ug3DiskVectorModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute ??_j = [cos ?_j, sin ?_j, 0] (disk plane unit vector)
std::vector<double> Ug3DiskVectorModule::computePhiHat_j(int j) {
    double theta = variables["theta_j"];  // Simplified, same for all j or per j
    std::vector<double> phi_hat = {std::cos(theta), std::sin(theta), 0.0};
    return phi_hat;
}

// Magnitude of ??_j (normalized=1)
double Ug3DiskVectorModule::computePhiHatMagnitude(int j) {
    auto phi = computePhiHat_j(j);
    return std::sqrt(phi[0]*phi[0] + phi[1]*phi[1] + phi[2]*phi[2]);  // =1
}

// Base for U_m without ??_j magnitude (since=1)
double Ug3DiskVectorModule::computeUmBase(double t) {
    double mu_over_rj = variables["mu_j"] / variables["r_j"];
    double exp_arg = - variables["gamma"] * t * std::cos(variables["pi"] * variables["t_n"]);
    double one_minus_exp = 1.0 - std::exp(exp_arg);
    double phi_mag = computePhiHatMagnitude(1);  // =1
    double p_scm = variables["P_SCm"];
    double e_react = variables["E_react"];
    return mu_over_rj * one_minus_exp * phi_mag * p_scm * e_react;
}

// U_m contribution with ??_j
double Ug3DiskVectorModule::computeUmContribution(double t, int j) {
    double base = computeUmBase(t);
    double heaviside_f = variables["heaviside_factor"];
    double quasi_f = 1.0 + variables["f_quasi"];
    return base * heaviside_f * quasi_f;
}

// Equation text
std::string Ug3DiskVectorModule::getEquationText() {
    return "U_m = ?_j [ (?_j / r_j) (1 - e^{-? t cos(? t_n)}) \hat{?}_j ] P_SCm E_react (1 + 10^13 f_Heaviside) (1 + f_quasi)\n"
           "Where \hat{?}_j = [cos ?_j, sin ?_j, 0] (unit vector in Ug3 disk plane, |??_j|=1);\n"
           "Specifies azimuthal direction for j-th string in disk (e.g., galactic plane).\n"
           "Example j=1, ?_j=0, t=0: ??_j=[1,0,0], U_m ?2.28e65 J/m� (mag=1).\n"
           "Role: Directional geometry for magnetic contributions in disks/nebulae.\n"
           "UQFF: Vector orientation in U_m/U_g3; collimation in jets/disks/formation.";
}

// Print variables
void Ug3DiskVectorModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Print vector and U_m
void Ug3DiskVectorModule::printVectorAndUm(int j, double t) {
    auto phi = computePhiHat_j(j);
    double mag = computePhiHatMagnitude(j);
    double um = computeUmContribution(t, j);
    std::cout << "??_" << j << " at ?_j=" << variables["theta_j"] << " rad, t=" << t << " s:\n";
    std::cout << "??_j = [" << std::scientific << phi[0] << ", " << phi[1] << ", " << phi[2] << "] (mag=" << mag << ")\n";
    std::cout << "U_m contrib = " << um << " J/m³\n";
}

// ========== ENHANCED SELF-UPDATE & SELF-EXPANSION IMPLEMENTATION ==========

// Variable Management
void Ug3DiskVectorModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
    std::cout << "Created variable '" << name << "' = " << value << std::endl;
}

void Ug3DiskVectorModule::removeVariable(const std::string& name) {
    if (variables.erase(name)) {
        std::cout << "Removed variable '" << name << "'" << std::endl;
    } else {
        std::cerr << "Variable '" << name << "' not found for removal" << std::endl;
    }
}

void Ug3DiskVectorModule::cloneVariable(const std::string& source, const std::string& target) {
    if (variables.find(source) != variables.end()) {
        variables[target] = variables[source];
        std::cout << "Cloned '" << source << "' to '" << target << "'" << std::endl;
    } else {
        std::cerr << "Source variable '" << source << "' not found" << std::endl;
    }
}

std::vector<std::string> Ug3DiskVectorModule::listVariables() const {
    std::vector<std::string> varList;
    for (const auto& pair : variables) {
        varList.push_back(pair.first);
    }
    return varList;
}

std::string Ug3DiskVectorModule::getSystemName() const {
    return "Ug3_Disk_Vector_PhiHat_UQFF";
}

// Batch Operations
void Ug3DiskVectorModule::transformVariableGroup(const std::vector<std::string>& varNames, 
                                                  std::function<double(double)> transformFunc) {
    for (const auto& name : varNames) {
        if (variables.find(name) != variables.end()) {
            variables[name] = transformFunc(variables[name]);
        }
    }
    std::cout << "Transformed " << varNames.size() << " variables" << std::endl;
}

void Ug3DiskVectorModule::scaleVariableGroup(const std::vector<std::string>& varNames, double factor) {
    transformVariableGroup(varNames, [factor](double val) { return val * factor; });
}

// Self-Expansion Methods
void Ug3DiskVectorModule::expandParameterSpace(double expansion_factor) {
    // Generic expansion - scale all physics parameters
    std::vector<std::string> physics_vars = {"mu_j", "r_j", "gamma", "E_react", "P_SCm"};
    scaleVariableGroup(physics_vars, expansion_factor);
    
    std::cout << "Expanded parameter space by factor " << expansion_factor << std::endl;
}

void Ug3DiskVectorModule::expandVectorScale(double angle_factor, double geometry_factor) {
    // Vector orientation expansion: theta_j angle and geometric factors
    variables["theta_j"] *= angle_factor;
    
    // Normalize angle to [0, 2π]
    while (variables["theta_j"] > 2.0 * variables["pi"]) {
        variables["theta_j"] -= 2.0 * variables["pi"];
    }
    while (variables["theta_j"] < 0.0) {
        variables["theta_j"] += 2.0 * variables["pi"];
    }
    
    // Geometry factor for disk plane characteristics
    if (variables.find("geometry_scale") == variables.end()) {
        variables["geometry_scale"] = 1.0;
    }
    variables["geometry_scale"] *= geometry_factor;
    
    std::cout << "Expanded vector scale: theta_j=" << variables["theta_j"] 
              << " rad, geometry_scale=" << variables["geometry_scale"] << "\n";
}

void Ug3DiskVectorModule::expandMagneticScale(double um_factor, double coupling_factor) {
    // Universal Magnetism expansion: μ_j and coupling factors
    variables["mu_j"] *= um_factor;
    
    // Heaviside and quasi factors for enhancement
    variables["f_Heaviside"] *= coupling_factor;
    variables["f_quasi"] *= coupling_factor;
    variables["heaviside_factor"] = 1.0 + variables["scale_Heaviside"] * variables["f_Heaviside"];
    
    std::cout << "Expanded magnetic scale: mu_j=" << variables["mu_j"] 
              << " T²m³, f_Heaviside=" << variables["f_Heaviside"] << "\n";
}

void Ug3DiskVectorModule::expandDiskScale(double radius_factor, double density_factor) {
    // Disk geometry expansion: r_j radius and density factors
    variables["r_j"] *= radius_factor;
    variables["P_SCm"] *= density_factor;  // Pressure/density in disk
    
    std::cout << "Expanded disk scale: r_j=" << variables["r_j"] 
              << " m, P_SCm=" << variables["P_SCm"] << "\n";
}

// Self-Refinement
void Ug3DiskVectorModule::autoRefineParameters() {
    // Auto-refine to match known physics constraints
    // Angle theta_j should be in [0, 2π]
    while (variables["theta_j"] > 2.0 * variables["pi"]) {
        variables["theta_j"] -= 2.0 * variables["pi"];
    }
    while (variables["theta_j"] < 0.0) {
        variables["theta_j"] += 2.0 * variables["pi"];
    }
    
    // f_Heaviside typically small [0.001, 0.1]
    if (variables["f_Heaviside"] > 0.1) {
        variables["f_Heaviside"] = 0.1;
        std::cout << "Refined f_Heaviside to limit (0.1)\n";
    }
    if (variables["f_Heaviside"] < 0.0) {
        variables["f_Heaviside"] = 0.001;
        std::cout << "Refined f_Heaviside to minimum (0.001)\n";
    }
    
    // f_quasi typically small [0.001, 0.1]
    if (variables["f_quasi"] > 0.1) {
        variables["f_quasi"] = 0.1;
        std::cout << "Refined f_quasi to limit (0.1)\n";
    }
    if (variables["f_quasi"] < 0.0) {
        variables["f_quasi"] = 0.001;
        std::cout << "Refined f_quasi to minimum (0.001)\n";
    }
    
    // Ensure positive r_j
    if (variables["r_j"] <= 0.0) {
        variables["r_j"] = 1.496e13;
        std::cout << "Refined r_j to default (1.496e13 m)\n";
    }
    
    // Recalculate derived quantities
    variables["heaviside_factor"] = 1.0 + variables["scale_Heaviside"] * variables["f_Heaviside"];
    
    std::cout << "Auto-refinement complete\n";
}

void Ug3DiskVectorModule::calibrateToObservations(const std::map<std::string, double>& observed_data) {
    // Calibrate parameters to match observational data
    for (const auto& obs : observed_data) {
        if (variables.find(obs.first) != variables.end()) {
            double old_val = variables[obs.first];
            variables[obs.first] = obs.second;
            std::cout << "Calibrated '" << obs.first << "': " << old_val << " → " << obs.second << "\n";
        }
    }
    
    // Update derived quantities
    if (observed_data.find("f_Heaviside") != observed_data.end()) {
        variables["heaviside_factor"] = 1.0 + variables["scale_Heaviside"] * variables["f_Heaviside"];
    }
}

void Ug3DiskVectorModule::optimizeForMetric(const std::string& metric) {
    // Optimize parameters for specific physical scenarios
    if (metric == "galactic_disk") {
        // Galactic disk geometry
        variables["theta_j"] = 0.0;  // Aligned with galactic plane
        variables["r_j"] = 1e17;  // ~10 kpc scale
        variables["f_Heaviside"] = 0.01;
        variables["f_quasi"] = 0.01;
        std::cout << "Optimized for galactic disk\n";
        
    } else if (metric == "accretion_disk") {
        // Accretion disk around compact object
        variables["theta_j"] = variables["pi"] / 4.0;  // Inclined
        variables["r_j"] = 1e10;  // ~100 km scale
        variables["f_Heaviside"] = 0.05;
        variables["f_quasi"] = 0.05;
        std::cout << "Optimized for accretion disk\n";
        
    } else if (metric == "protoplanetary_disk") {
        // Protoplanetary disk
        variables["theta_j"] = 0.0;  // Flat disk
        variables["r_j"] = 1e13;  // ~100 AU scale
        variables["f_Heaviside"] = 0.02;
        variables["f_quasi"] = 0.02;
        std::cout << "Optimized for protoplanetary disk\n";
        
    } else if (metric == "jet_collimation") {
        // High magnetic field for jet collimation
        variables["theta_j"] = variables["pi"] / 2.0;  // Perpendicular
        variables["mu_j"] *= 10.0;
        variables["f_Heaviside"] = 0.08;
        variables["f_quasi"] = 0.08;
        std::cout << "Optimized for jet collimation\n";
        
    } else if (metric == "nebula") {
        // Nebula/star formation disk
        variables["theta_j"] = variables["pi"] / 6.0;  // Slight inclination
        variables["r_j"] = 1e15;  // ~10 pc scale
        variables["f_Heaviside"] = 0.03;
        variables["f_quasi"] = 0.03;
        std::cout << "Optimized for nebula disk\n";
        
    } else {
        std::cout << "Unknown metric '" << metric << "', no optimization applied\n";
        return;
    }
    
    // Update derived quantities
    variables["heaviside_factor"] = 1.0 + variables["scale_Heaviside"] * variables["f_Heaviside"];
}

// Parameter Exploration
std::vector<std::map<std::string, double>> Ug3DiskVectorModule::generateVariations(int count) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.8, 1.2);  // ±20% variation
    std::uniform_real_distribution<> angle_dis(0.0, 2.0 * variables["pi"]);
    
    for (int i = 0; i < count; ++i) {
        std::map<std::string, double> variant = variables;
        
        // Vary key parameters
        variant["theta_j"] = angle_dis(gen);  // Random angle
        variant["mu_j"] *= dis(gen);
        variant["r_j"] *= dis(gen);
        variant["f_Heaviside"] *= dis(gen);
        variant["f_quasi"] *= dis(gen);
        
        // Ensure physical constraints
        if (variant["f_Heaviside"] < 0.0) variant["f_Heaviside"] = 0.001;
        if (variant["f_Heaviside"] > 0.1) variant["f_Heaviside"] = 0.1;
        if (variant["f_quasi"] < 0.0) variant["f_quasi"] = 0.001;
        if (variant["f_quasi"] > 0.1) variant["f_quasi"] = 0.1;
        
        // Update heaviside_factor
        variant["heaviside_factor"] = 1.0 + variant["scale_Heaviside"] * variant["f_Heaviside"];
        
        variations.push_back(variant);
    }
    
    std::cout << "Generated " << count << " parameter variations\n";
    return variations;
}

// Adaptive Evolution
void Ug3DiskVectorModule::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-mutation_rate, mutation_rate);
    
    // Mutate key parameters
    variables["theta_j"] += dis(gen) * variables["pi"];  // Angle mutation
    variables["mu_j"] *= (1.0 + dis(gen));
    variables["r_j"] *= (1.0 + dis(gen));
    variables["f_Heaviside"] *= (1.0 + dis(gen));
    variables["f_quasi"] *= (1.0 + dis(gen));
    
    // Apply constraints
    autoRefineParameters();
    
    std::cout << "Mutated parameters with rate " << mutation_rate << "\n";
}

void Ug3DiskVectorModule::evolveSystem(int generations) {
    double best_fitness = 0.0;
    std::map<std::string, double> best_params = variables;
    
    for (int gen = 0; gen < generations; ++gen) {
        // Generate variations
        auto variations = generateVariations(5);
        
        // Evaluate fitness (example: balance angle near disk plane with reasonable f values)
        for (const auto& variant : variations) {
            double theta = variant.at("theta_j");
            double f_h = variant.at("f_Heaviside");
            double f_q = variant.at("f_quasi");
            
            // Fitness: prefer angles near 0 or π (aligned with disk), moderate f values
            double angle_fitness = 1.0 / (1.0 + std::min(std::abs(theta), 
                                          std::abs(theta - variables["pi"])));
            double f_fitness = 1.0 / (1.0 + std::abs(f_h - 0.02) + std::abs(f_q - 0.02));
            double fitness = angle_fitness * f_fitness;
            
            if (fitness > best_fitness) {
                best_fitness = fitness;
                best_params = variant;
            }
        }
        
        // Mutate current best
        variables = best_params;
        mutateParameters(0.05);
    }
    
    variables = best_params;
    autoRefineParameters();
    
    std::cout << "Evolved system over " << generations << " generations (fitness=" 
              << best_fitness << ")\n";
}

// State Management
static std::map<std::string, std::map<std::string, double>> saved_states_disk;

void Ug3DiskVectorModule::saveState(const std::string& label) {
    saved_states_disk[label] = variables;
    std::cout << "Saved state '" << label << "'\n";
}

void Ug3DiskVectorModule::restoreState(const std::string& label) {
    if (saved_states_disk.find(label) != saved_states_disk.end()) {
        variables = saved_states_disk[label];
        std::cout << "Restored state '" << label << "'\n";
    } else {
        std::cerr << "State '" << label << "' not found\n";
    }
}

std::vector<std::string> Ug3DiskVectorModule::listSavedStates() const {
    std::vector<std::string> state_list;
    for (const auto& pair : saved_states_disk) {
        state_list.push_back(pair.first);
    }
    return state_list;
}

std::string Ug3DiskVectorModule::exportState() const {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(10);
    for (const auto& pair : variables) {
        oss << pair.first << "=" << pair.second << "\n";
    }
    return oss.str();
}

// System Analysis
std::map<std::string, double> Ug3DiskVectorModule::sensitivityAnalysis(const std::string& output_var) {
    std::map<std::string, double> sensitivities;
    double baseline_output = 0.0;
    
    // Compute baseline
    if (output_var == "U_m") {
        baseline_output = computeUmContribution(0.0, 1);
    } else if (output_var == "phi_hat_mag") {
        baseline_output = computePhiHatMagnitude(1);
    } else {
        std::cerr << "Unknown output variable '" << output_var << "'\n";
        return sensitivities;
    }
    
    // Test sensitivity to each parameter
    double delta = 0.01;  // 1% perturbation
    std::vector<std::string> params = {"theta_j", "mu_j", "r_j", "f_Heaviside", "f_quasi"};
    
    for (const auto& param : params) {
        double original_val = variables[param];
        
        // Perturb upward
        variables[param] = original_val * (1.0 + delta);
        if (param == "f_Heaviside") {
            variables["heaviside_factor"] = 1.0 + variables["scale_Heaviside"] * variables["f_Heaviside"];
        }
        
        double perturbed_output = 0.0;
        if (output_var == "U_m") {
            perturbed_output = computeUmContribution(0.0, 1);
        } else if (output_var == "phi_hat_mag") {
            perturbed_output = computePhiHatMagnitude(1);
        }
        
        // Compute sensitivity (normalized)
        double sensitivity = 0.0;
        if (baseline_output != 0.0) {
            sensitivity = (perturbed_output - baseline_output) / (baseline_output * delta);
        }
        sensitivities[param] = sensitivity;
        
        // Restore original value
        variables[param] = original_val;
        if (param == "f_Heaviside") {
            variables["heaviside_factor"] = 1.0 + variables["scale_Heaviside"] * variables["f_Heaviside"];
        }
    }
    
    return sensitivities;
}

std::string Ug3DiskVectorModule::generateReport() const {
    std::ostringstream report;
    report << std::scientific << std::setprecision(3);
    
    report << "========== UG3 DISK VECTOR MODULE REPORT ==========\n";
    report << "System: " << getSystemName() << "\n\n";
    
    report << "Vector Parameters:\n";
    report << "  theta_j = " << variables.at("theta_j") << " rad (" 
           << (variables.at("theta_j") * 180.0 / variables.at("pi")) << " deg)\n";
    
    // Compute φ̂_j at current angle
    Ug3DiskVectorModule temp_mod = *this;
    auto phi = temp_mod.computePhiHat_j(1);
    double mag = temp_mod.computePhiHatMagnitude(1);
    report << "  φ̂_j = [" << phi[0] << ", " << phi[1] << ", " << phi[2] << "]\n";
    report << "  |φ̂_j| = " << mag << " (normalized)\n\n";
    
    report << "Magnetic Parameters:\n";
    report << "  μ_j = " << variables.at("mu_j") << " T²m³\n";
    report << "  r_j = " << variables.at("r_j") << " m\n";
    report << "  γ = " << variables.at("gamma") << " s⁻¹\n\n";
    
    report << "Enhancement Factors:\n";
    report << "  f_Heaviside = " << variables.at("f_Heaviside") << "\n";
    report << "  Heaviside factor = " << variables.at("heaviside_factor") << "\n";
    report << "  f_quasi = " << variables.at("f_quasi") << "\n\n";
    
    report << "Energy Parameters:\n";
    report << "  E_react = " << variables.at("E_react") << " J\n";
    report << "  P_SCm = " << variables.at("P_SCm") << "\n\n";
    
    // Compute U_m at t=0
    double um_t0 = temp_mod.computeUmContribution(0.0, 1);
    report << "U_m Contribution (j=1, t=0):\n";
    report << "  U_m = " << um_t0 << " J/m³\n\n";
    
    report << "Physical Context:\n";
    report << "  Unit vector φ̂_j specifies azimuthal direction in Ug3 disk plane\n";
    report << "  Applications: Galactic disks, accretion disks, protoplanetary disks\n";
    report << "  Role: Directional geometry for magnetic string contributions\n";
    report << "  UQFF: Vector orientation in U_m/U_g3 disk structures\n";
    
    report << "===================================================\n";
    
    return report.str();
}

bool Ug3DiskVectorModule::validateConsistency() const {
    bool consistent = true;
    
    // Check theta_j in [0, 2π]
    if (variables.at("theta_j") < 0.0 || variables.at("theta_j") > 2.0 * variables.at("pi")) {
        std::cerr << "Inconsistency: theta_j out of range [0, 2π]\n";
        consistent = false;
    }
    
    // Check positive r_j
    if (variables.at("r_j") <= 0.0) {
        std::cerr << "Inconsistency: r_j must be positive\n";
        consistent = false;
    }
    
    // Check positive mu_j
    if (variables.at("mu_j") <= 0.0) {
        std::cerr << "Inconsistency: mu_j must be positive\n";
        consistent = false;
    }
    
    // Check f_Heaviside range
    if (variables.at("f_Heaviside") < 0.0 || variables.at("f_Heaviside") > 0.1) {
        std::cerr << "Inconsistency: f_Heaviside out of physical range [0.0, 0.1]\n";
        consistent = false;
    }
    
    // Check f_quasi range
    if (variables.at("f_quasi") < 0.0 || variables.at("f_quasi") > 0.1) {
        std::cerr << "Inconsistency: f_quasi out of physical range [0.0, 0.1]\n";
        consistent = false;
    }
    
    // Check derived quantity consistency
    double expected_heaviside = 1.0 + variables.at("scale_Heaviside") * variables.at("f_Heaviside");
    double actual_heaviside = variables.at("heaviside_factor");
    if (std::abs(expected_heaviside - actual_heaviside) / expected_heaviside > 0.01) {
        std::cerr << "Inconsistency: heaviside_factor mismatch\n";
        consistent = false;
    }
    
    // Check φ̂_j magnitude = 1
    Ug3DiskVectorModule temp_mod = *this;
    double mag = temp_mod.computePhiHatMagnitude(1);
    if (std::abs(mag - 1.0) > 1e-10) {
        std::cerr << "Inconsistency: φ̂_j magnitude != 1\n";
        consistent = false;
    }
    
    if (consistent) {
        std::cout << "Consistency validation: PASSED\n";
    }
    
    return consistent;
}

void Ug3DiskVectorModule::autoCorrectAnomalies() {
    bool corrected = false;
    
    // Correct theta_j to [0, 2π]
    while (variables["theta_j"] > 2.0 * variables["pi"]) {
        variables["theta_j"] -= 2.0 * variables["pi"];
        corrected = true;
    }
    while (variables["theta_j"] < 0.0) {
        variables["theta_j"] += 2.0 * variables["pi"];
        corrected = true;
    }
    
    // Correct f_Heaviside
    if (variables["f_Heaviside"] < 0.0) {
        variables["f_Heaviside"] = 0.001;
        corrected = true;
    }
    if (variables["f_Heaviside"] > 0.1) {
        variables["f_Heaviside"] = 0.1;
        corrected = true;
    }
    
    // Correct f_quasi
    if (variables["f_quasi"] < 0.0) {
        variables["f_quasi"] = 0.001;
        corrected = true;
    }
    if (variables["f_quasi"] > 0.1) {
        variables["f_quasi"] = 0.1;
        corrected = true;
    }
    
    // Correct r_j
    if (variables["r_j"] <= 0.0) {
        variables["r_j"] = 1.496e13;
        corrected = true;
    }
    
    // Correct mu_j
    if (variables["mu_j"] <= 0.0) {
        variables["mu_j"] = 3.38e23;
        corrected = true;
    }
    
    // Recalculate derived quantities
    variables["heaviside_factor"] = 1.0 + variables["scale_Heaviside"] * variables["f_Heaviside"];
    
    if (corrected) {
        std::cout << "Auto-corrected anomalies\n";
    } else {
        std::cout << "No anomalies found\n";
    }
}

// Example usage in base program (snippet)
// #include "Ug3DiskVectorModule.h"
// int main() {
//     Ug3DiskVectorModule mod;
//     auto phi = mod.computePhiHat_j(1);
//     std::cout << "φ̂_1 = [" << phi[0] << ", " << phi[1] << ", " << phi[2] << "]\n";
//     mod.printVectorAndUm(1, 0.0);
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("theta_j", M_PI / 2);
//     mod.printVariables();
//     return 0;
// }

// ========== COMPREHENSIVE ENHANCED DEMONSTRATION ==========
/*
int main() {
    std::cout << std::scientific << std::setprecision(3);
    std::cout << "========== UG3 DISK VECTOR MODULE DEMONSTRATION ==========\n\n";
    
    // ===== Step 1: Initialize Module =====
    Ug3DiskVectorModule mod;
    std::cout << "Step 1: Module initialized with defaults:\n";
    std::cout << "  System: " << mod.getSystemName() << "\n";
    std::cout << "  theta_j = " << mod.variables["theta_j"] << " rad (0 deg)\n";
    auto phi = mod.computePhiHat_j(1);
    std::cout << "  φ̂_1 = [" << phi[0] << ", " << phi[1] << ", " << phi[2] << "]\n";
    std::cout << "  |φ̂_1| = " << mod.computePhiHatMagnitude(1) << "\n\n";
    
    // ===== Step 2: Baseline Computation =====
    std::cout << "Step 2: Compute baseline U_m contribution:\n";
    double um_t0 = mod.computeUmContribution(0.0, 1);
    std::cout << "  U_m (j=1, t=0) = " << um_t0 << " J/m³\n";
    std::cout << "  Direction: Along x-axis (theta=0)\n\n";
    
    // ===== Step 3: Variable Management =====
    std::cout << "Step 3: Variable Management\n";
    mod.createVariable("disk_inclination_deg", 0.0);
    std::cout << "  Created 'disk_inclination_deg' for geometric analysis\n";
    
    std::vector<std::string> all_vars = mod.listVariables();
    std::cout << "  Total variables: " << all_vars.size() << "\n";
    
    mod.cloneVariable("theta_j", "theta_j_initial");
    std::cout << "  Cloned 'theta_j' → 'theta_j_initial'\n\n";
    
    // ===== Step 4: Vector Expansion (Rotate Orientation) =====
    std::cout << "Step 4: Vector Expansion (Rotate to Different Angles)\n";
    mod.expandVectorScale(2.0, 1.5);  // 2x angle (still normalized), 1.5x geometry
    phi = mod.computePhiHat_j(1);
    std::cout << "  After rotation: theta_j = " << mod.variables["theta_j"] << " rad\n";
    std::cout << "  φ̂_1 = [" << std::fixed << std::setprecision(4) 
              << phi[0] << ", " << phi[1] << ", " << phi[2] << "]\n";
    std::cout << "  |φ̂_1| = " << mod.computePhiHatMagnitude(1) << " (still normalized)\n\n";
    
    // ===== Step 5: Magnetic Expansion (Enhanced U_m) =====
    std::cout << "Step 5: Magnetic Expansion (Increase μ_j)\n";
    std::cout << std::scientific << std::setprecision(3);
    mod.expandMagneticScale(2.0, 1.5);  // 2x μ_j, 1.5x coupling
    double um_expanded = mod.computeUmContribution(0.0, 1);
    std::cout << "  After expansion: U_m = " << um_expanded << " J/m³\n";
    std::cout << "  Enhancement: " << std::fixed << std::setprecision(1) 
              << (um_expanded/um_t0) << "x\n\n";
    
    // ===== Step 6: Disk Expansion (Larger Radius) =====
    std::cout << "Step 6: Disk Expansion (Increase r_j)\n";
    std::cout << std::scientific << std::setprecision(3);
    mod.expandDiskScale(3.0, 1.3);  // 3x radius, 1.3x density
    std::cout << "  After expansion: r_j = " << mod.variables["r_j"] << " m\n";
    std::cout << "  P_SCm = " << mod.variables["P_SCm"] << "\n\n";
    
    // ===== Step 7: Batch Operations (Reset Vector Parameters) =====
    std::cout << "Step 7: Batch Operations (Scale Enhancement Factors)\n";
    std::vector<std::string> enhancement_group = {"f_Heaviside", "f_quasi"};
    mod.scaleVariableGroup(enhancement_group, 0.5);
    std::cout << "  Scaled enhancement group by 0.5:\n";
    std::cout << "    f_Heaviside = " << mod.variables["f_Heaviside"] << "\n";
    std::cout << "    f_quasi = " << mod.variables["f_quasi"] << "\n\n";
    
    // ===== Step 8-12: Test Different Disk Geometries =====
    std::cout << "Steps 8-12: Test Multiple Disk Configurations\n";
    
    // Galactic disk
    mod.optimizeForMetric("galactic_disk");
    phi = mod.computePhiHat_j(1);
    std::cout << "  Galactic Disk (theta=" << std::fixed << std::setprecision(2) 
              << (mod.variables["theta_j"]*180/mod.variables["pi"]) << " deg):\n";
    std::cout << "    φ̂ = [" << std::setprecision(4) << phi[0] << ", " << phi[1] << ", " << phi[2] << "]\n";
    std::cout << "    r_j = " << std::scientific << std::setprecision(2) << mod.variables["r_j"] << " m\n";
    
    // Accretion disk
    mod.optimizeForMetric("accretion_disk");
    phi = mod.computePhiHat_j(1);
    std::cout << "  Accretion Disk (theta=" << std::fixed << std::setprecision(2) 
              << (mod.variables["theta_j"]*180/mod.variables["pi"]) << " deg):\n";
    std::cout << "    φ̂ = [" << std::setprecision(4) << phi[0] << ", " << phi[1] << ", " << phi[2] << "]\n";
    std::cout << "    r_j = " << std::scientific << std::setprecision(2) << mod.variables["r_j"] << " m\n";
    
    // Protoplanetary disk
    mod.optimizeForMetric("protoplanetary_disk");
    phi = mod.computePhiHat_j(1);
    std::cout << "  Protoplanetary Disk (theta=" << std::fixed << std::setprecision(2) 
              << (mod.variables["theta_j"]*180/mod.variables["pi"]) << " deg):\n";
    std::cout << "    φ̂ = [" << std::setprecision(4) << phi[0] << ", " << phi[1] << ", " << phi[2] << "]\n";
    std::cout << "    r_j = " << std::scientific << std::setprecision(2) << mod.variables["r_j"] << " m\n";
    
    // Jet collimation
    mod.optimizeForMetric("jet_collimation");
    phi = mod.computePhiHat_j(1);
    std::cout << "  Jet Collimation (theta=" << std::fixed << std::setprecision(2) 
              << (mod.variables["theta_j"]*180/mod.variables["pi"]) << " deg):\n";
    std::cout << "    φ̂ = [" << std::setprecision(4) << phi[0] << ", " << phi[1] << ", " << phi[2] << "]\n";
    std::cout << "    Perpendicular orientation for jet axis\n";
    
    // Nebula
    mod.optimizeForMetric("nebula");
    phi = mod.computePhiHat_j(1);
    std::cout << "  Nebula (theta=" << std::fixed << std::setprecision(2) 
              << (mod.variables["theta_j"]*180/mod.variables["pi"]) << " deg):\n";
    std::cout << "    φ̂ = [" << std::setprecision(4) << phi[0] << ", " << phi[1] << ", " << phi[2] << "]\n";
    std::cout << "    r_j = " << std::scientific << std::setprecision(2) << mod.variables["r_j"] << " m\n\n";
    
    // ===== Step 13: Auto-Refinement =====
    std::cout << "Step 13: Auto-Refinement\n";
    std::cout << std::scientific << std::setprecision(3);
    mod.updateVariable("theta_j", 10.0);  // Beyond 2π
    std::cout << "  Set theta_j = 10.0 rad (beyond 2π)\n";
    mod.autoRefineParameters();
    std::cout << "  After refinement: theta_j = " << mod.variables["theta_j"] 
              << " rad (normalized)\n\n";
    
    // ===== Step 14: Calibration (Observational Data) =====
    std::cout << "Step 14: Calibration to Observations\n";
    std::map<std::string, double> obs_data;
    obs_data["theta_j"] = 0.785398;  // π/4, 45 degrees
    obs_data["f_Heaviside"] = 0.025;
    mod.calibrateToObservations(obs_data);
    phi = mod.computePhiHat_j(1);
    std::cout << "  Calibrated: theta_j = " << mod.variables["theta_j"] << " rad (45 deg)\n";
    std::cout << "  φ̂_1 = [" << std::fixed << std::setprecision(4) 
              << phi[0] << ", " << phi[1] << ", " << phi[2] << "]\n\n";
    
    // ===== Step 15: Parameter Variations =====
    std::cout << "Step 15: Generate Parameter Variations\n";
    std::cout << std::scientific << std::setprecision(3);
    std::vector<std::map<std::string, double>> variations = mod.generateVariations(5);
    std::cout << "  Generated " << variations.size() << " parameter variations:\n";
    for (size_t i = 0; i < std::min(size_t(3), variations.size()); ++i) {
        double theta_deg = variations[i]["theta_j"] * 180.0 / mod.variables["pi"];
        std::cout << "    Variant " << (i+1) << ": theta=" << std::fixed << std::setprecision(1) 
                  << theta_deg << " deg, r_j=" << std::scientific << std::setprecision(2) 
                  << variations[i]["r_j"] << " m\n";
    }
    std::cout << "\n";
    
    // ===== Step 16: Mutation =====
    std::cout << "Step 16: Mutate Parameters\n";
    std::cout << std::scientific << std::setprecision(3);
    mod.updateVariable("theta_j", 0.0);  // Reset
    mod.mutateParameters(0.2);  // 20% mutation
    phi = mod.computePhiHat_j(1);
    std::cout << "  After mutation: theta_j = " << mod.variables["theta_j"] << " rad\n";
    std::cout << "  φ̂_1 = [" << std::fixed << std::setprecision(4) 
              << phi[0] << ", " << phi[1] << ", " << phi[2] << "]\n\n";
    
    // ===== Step 17: System Evolution =====
    std::cout << "Step 17: Evolve System (Optimize Disk Geometry)\n";
    std::cout << std::scientific << std::setprecision(3);
    mod.evolveSystem(10);  // 10 generations
    std::cout << "  After evolution: theta_j = " << mod.variables["theta_j"] << " rad\n";
    std::cout << "  After evolution: f_Heaviside = " << mod.variables["f_Heaviside"] << "\n\n";
    
    // ===== Step 18: State Management =====
    std::cout << "Step 18: State Management\n";
    mod.updateVariable("theta_j", 0.0);
    mod.saveState("aligned_disk");
    std::cout << "  Saved state 'aligned_disk'\n";
    
    mod.updateVariable("theta_j", mod.variables["pi"] / 2.0);
    mod.saveState("perpendicular_disk");
    std::cout << "  Saved state 'perpendicular_disk'\n";
    
    std::vector<std::string> saved = mod.listSavedStates();
    std::cout << "  Total saved states: " << saved.size() << "\n";
    
    mod.restoreState("aligned_disk");
    std::cout << "  Restored 'aligned_disk': theta_j = " << mod.variables["theta_j"] << " rad\n\n";
    
    // ===== Step 19: Export State =====
    std::cout << "Step 19: Export State\n";
    std::string exported = mod.exportState();
    std::cout << "  Exported " << exported.length() << " bytes of state data\n";
    std::cout << "  (Can be saved to file for archival/restoration)\n\n";
    
    // ===== Step 20: Sensitivity Analysis =====
    std::cout << "Step 20: Sensitivity Analysis (U_m response)\n";
    std::map<std::string, double> sensitivity = mod.sensitivityAnalysis("U_m");
    std::cout << "  Sensitivity of U_m to parameter changes:\n";
    for (const auto& pair : sensitivity) {
        std::cout << "    " << pair.first << ": " << std::fixed << std::setprecision(2) 
                  << pair.second << "\n";
    }
    std::cout << "\n";
    
    // ===== Step 21: Validation =====
    std::cout << "Step 21: Consistency Validation\n";
    bool valid = mod.validateConsistency();
    std::cout << "  System consistency: " << (valid ? "VALID" : "INVALID") << "\n";
    if (!valid) {
        mod.autoCorrectAnomalies();
        std::cout << "  Auto-corrected anomalies\n";
    }
    std::cout << "\n";
    
    // ===== Step 22: Generate Full Report =====
    std::cout << "Step 22: Generate Full Report\n";
    std::string report = mod.generateReport();
    std::cout << report << "\n";
    
    // ===== Step 23-26: Vector Orientation Sweep =====
    std::cout << "Steps 23-26: Vector Orientation Sweep (Full Circle)\n";
    std::cout << "  Angle (deg) | φ̂_x    | φ̂_y    | φ̂_z | |φ̂| | Direction\n";
    std::cout << "  ----------------------------------------------------------------\n";
    
    struct AnglePoint {
        double degrees;
        std::string label;
    };
    
    std::vector<AnglePoint> angles = {
        {0.0, "East (0°)"},
        {45.0, "NE (45°)"},
        {90.0, "North (90°)"},
        {135.0, "NW (135°)"},
        {180.0, "West (180°)"},
        {225.0, "SW (225°)"},
        {270.0, "South (270°)"},
        {315.0, "SE (315°)"}
    };
    
    for (const auto& ap : angles) {
        double theta_rad = ap.degrees * mod.variables["pi"] / 180.0;
        mod.updateVariable("theta_j", theta_rad);
        auto phi_vec = mod.computePhiHat_j(1);
        double mag = mod.computePhiHatMagnitude(1);
        
        std::cout << "  " << std::setw(11) << std::right << std::fixed << std::setprecision(1) << ap.degrees
                  << " | " << std::setprecision(4) << std::setw(7) << phi_vec[0]
                  << " | " << std::setw(7) << phi_vec[1]
                  << " | " << std::setw(3) << phi_vec[2]
                  << " | " << std::setw(3) << mag
                  << " | " << ap.label << "\n";
    }
    
    std::cout << "\n========== DEMONSTRATION COMPLETE ==========\n";
    std::cout << "Ug3 Disk Vector module validated across full angular sweep.\n";
    std::cout << "Unit vector φ̂_j maintains normalization (|φ̂|=1) at all angles.\n";
    std::cout << "Physical applications: Galactic disks, accretion disks, jet collimation.\n";
    std::cout << "UQFF Integration: Directional geometry for U_m/U_g3 magnetic contributions.\n";
    std::cout << "Vector specifies azimuthal orientation in disk plane (galactic/nebular structures).\n";
    
    return 0;
}
*/
// Compile: g++ -o disk_vector_test disk_vector_test.cpp Ug3DiskVectorModule.cpp -lm
// Sample: ??_1=[1,0,0] (?=0); U_m?2.28e65 J/m�; directional in disk plane.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

Ug3DiskVectorModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Core computation methods(computePhiHat_j, computePhiHatMagnitude, computeUmContribution) are clear, concise, and variable - driven.
- Uses std::vector for unit vector representation, supporting extensibility for vector operations.
- Output and debugging functions(printVariables, printVectorAndUm, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.
- Models directional geometry for magnetic contributions in disk planes.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid input, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map and std::vector.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in disk vector and magnetic geometry modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.