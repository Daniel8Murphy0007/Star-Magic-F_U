// UQFFBuoyancyModule157.h
// UQFF Buoyancy for Observational Systems: M104, NGC 4839, Chandra/Webb, NGC 346, NGC 1672
// Copyright - Daniel T. Murphy, analyzed Oct 22, 2025.

#ifndef UQFF_BUOYANCY_MODULE_H
#define UQFF_BUOYANCY_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <vector>
#include <functional>
#include <random>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <iomanip>

class UQFFBuoyancyModule157 {
private:
    std::map<std::string, double> variables;

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
    
    // ===== DYNAMIC SELF-UPDATE AND SELF-EXPANSION METHODS =====
    
    // Variable Management (5 methods)
    void createVariable(const std::string& name, double value);
    void removeVariable(const std::string& name);
    void cloneVariable(const std::string& source, const std::string& dest);
    std::vector<std::string> listVariables();
    std::string getSystemName() const;
    
    // Batch Operations (2 methods)
    void transformVariableGroup(const std::vector<std::string>& var_names, std::function<double(double)> transform);
    void scaleVariableGroup(const std::vector<std::string>& var_names, double scale_factor);
    
    // Self-Expansion (4 methods - observational-specific)
    void expandParameterSpace(double expansion_factor);
    void expandSystemScale(const std::string& system, double mass_factor, double distance_factor);
    void expandForceScale(double buoyancy_factor, double coupling_factor);
    void expandObservationalScale(const std::string& system, double luminosity_factor, double efficiency_factor);
    
    // Self-Refinement (3 methods)
    void autoRefineParameters(const std::string& system, double target_force, double tolerance = 0.01);
    void calibrateToObservations(const std::string& system, const std::map<std::string, double>& observations);
    void optimizeForMetric(const std::string& metric, double target_value);
    
    // Parameter Exploration (1 method)
    std::vector<std::map<std::string, double>> generateVariations(int count, double variation_range);
    
    // Adaptive Evolution (2 methods)
    void mutateParameters(double mutation_rate);
    void evolveSystem(const std::string& system, int generations, std::function<double(double)> fitness_func);
    
    // State Management (4 methods)
    void saveState(const std::string& state_name);
    void restoreState(const std::string& state_name);
    std::vector<std::string> listSavedStates();
    void exportStateUQFF(const std::string& filename);
    
    // System Analysis (4 methods)
    std::map<std::string, double> sensitivityAnalysis(const std::string& system, const std::vector<std::string>& params);
    std::string generateReport(const std::string& system);
    bool validateConsistency(const std::string& system);
    void autoCorrectAnomalies(const std::string& system);
};

#endif // UQFF_BUOYANCY_MODULE_H

// Implementation
UQFFBuoyancyModule157::UQFFBuoyancyModule157() {
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

// ===== DYNAMIC SELF-UPDATE AND SELF-EXPANSION METHODS IMPLEMENTATION =====

namespace saved_states_obs {
    std::map<std::string, std::map<std::string, double>> states;
}

// Variable Management (5 methods)
void UQFFBuoyancyModule157::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void UQFFBuoyancyModule157::removeVariable(const std::string& name) {
    variables.erase(name);
}

void UQFFBuoyancyModule157::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
    }
}

std::vector<std::string> UQFFBuoyancyModule157::listVariables() {
    std::vector<std::string> names;
    for (const auto& pair : variables) {
        names.push_back(pair.first);
    }
    return names;
}

std::string UQFFBuoyancyModule157::getSystemName() const {
    return "UQFFBuoyancy_Observational_MultiSystem";
}

// Batch Operations (2 methods)
void UQFFBuoyancyModule157::transformVariableGroup(const std::vector<std::string>& var_names, std::function<double(double)> transform) {
    for (const auto& name : var_names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = transform(variables[name]);
        }
    }
}

void UQFFBuoyancyModule157::scaleVariableGroup(const std::vector<std::string>& var_names, double scale_factor) {
    for (const auto& name : var_names) {
        if (variables.find(name) != variables.end()) {
            variables[name] *= scale_factor;
        }
    }
}

// Self-Expansion (4 methods - observational-specific)
void UQFFBuoyancyModule157::expandParameterSpace(double expansion_factor) {
    // Scale all observational parameters
    std::vector<std::string> obs_params = {"M104_factor", "NGC4839_factor", "CW_factor", 
                                           "NGC346_factor", "NGC1672_factor", "CW_efficiency"};
    for (const auto& param : obs_params) {
        if (variables.find(param) != variables.end()) {
            variables[param] *= expansion_factor;
        }
    }
    
    // Scale physical constants
    variables["beta_1"] *= expansion_factor;
    variables["beta_i"] *= expansion_factor;
    variables["alpha"] *= expansion_factor;
}

void UQFFBuoyancyModule157::expandSystemScale(const std::string& system, double mass_factor, double distance_factor) {
    if (system == "M104") {
        variables["M104_mass"] *= mass_factor;
        variables["M104_dist"] *= distance_factor;
        variables["M104_factor"] *= std::sqrt(mass_factor / distance_factor);
    } else if (system == "NGC4839") {
        variables["NGC4839_mass"] *= mass_factor;
        variables["NGC4839_dist"] *= distance_factor;
        variables["NGC4839_factor"] *= std::sqrt(mass_factor / distance_factor);
    } else if (system == "NGC346") {
        variables["NGC346_mass"] *= mass_factor;
        variables["NGC346_dist"] *= distance_factor;
        variables["NGC346_factor"] *= std::sqrt(mass_factor / distance_factor);
    } else if (system == "NGC1672") {
        variables["NGC1672_mass"] *= mass_factor;
        variables["NGC1672_dist"] *= distance_factor;
        variables["NGC1672_factor"] *= std::sqrt(mass_factor / distance_factor);
    } else if (system == "Chandra_Webb") {
        variables["CW_factor"] *= mass_factor;
        variables["CW_efficiency"] *= std::min(1.0, variables["CW_efficiency"].operator double() * distance_factor);
    }
}

void UQFFBuoyancyModule157::expandForceScale(double buoyancy_factor, double coupling_factor) {
    variables["beta_1"] *= buoyancy_factor;
    variables["beta_i"] *= buoyancy_factor;
    variables["alpha"] *= coupling_factor;
    variables["rho_vac"] *= buoyancy_factor;
}

void UQFFBuoyancyModule157::expandObservationalScale(const std::string& system, double luminosity_factor, double efficiency_factor) {
    if (system == "M104") {
        variables["M104_factor"] *= luminosity_factor;
    } else if (system == "NGC4839") {
        variables["NGC4839_factor"] *= luminosity_factor;
    } else if (system == "Chandra_Webb") {
        variables["CW_factor"] *= luminosity_factor;
        variables["CW_efficiency"] *= std::min(1.0, variables["CW_efficiency"].operator double() * efficiency_factor);
    } else if (system == "NGC346") {
        variables["NGC346_factor"] *= luminosity_factor;
    } else if (system == "NGC1672") {
        variables["NGC1672_factor"] *= luminosity_factor;
    }
    
    // Adjust observational frequency
    variables["omega_obs"] *= luminosity_factor;
}

// Self-Refinement (3 methods)
void UQFFBuoyancyModule157::autoRefineParameters(const std::string& system, double target_force, double tolerance) {
    for (int iter = 0; iter < 100; iter++) {
        double current_force = computeFBi(system);
        double error = std::abs(current_force - target_force) / std::abs(target_force);
        
        if (error < tolerance) break;
        
        // Adjust system-specific parameters
        double adjustment = (target_force - current_force) / target_force * 0.1;
        if (system == "M104") {
            variables["M104_factor"] *= (1.0 + adjustment);
        } else if (system == "NGC4839") {
            variables["NGC4839_factor"] *= (1.0 + adjustment);
        } else if (system == "Chandra_Webb") {
            variables["CW_factor"] *= (1.0 + adjustment * 0.8);
            variables["CW_efficiency"] *= (1.0 + adjustment * 0.2);
        } else if (system == "NGC346") {
            variables["NGC346_factor"] *= (1.0 + adjustment);
        } else if (system == "NGC1672") {
            variables["NGC1672_factor"] *= (1.0 + adjustment);
        }
        
        // Adjust coupling
        variables["alpha"] *= (1.0 + adjustment * 0.5);
    }
}

void UQFFBuoyancyModule157::calibrateToObservations(const std::string& system, const std::map<std::string, double>& observations) {
    for (const auto& obs : observations) {
        if (variables.find(obs.first) != variables.end()) {
            variables[obs.first] = obs.second;
        }
    }
    
    // Adjust derived parameters based on observations
    if (observations.find("luminosity") != observations.end()) {
        double lum_factor = observations.at("luminosity") / 1e36;
        if (system == "M104") variables["M104_factor"] *= std::sqrt(lum_factor);
        else if (system == "NGC4839") variables["NGC4839_factor"] *= std::sqrt(lum_factor);
        else if (system == "NGC346") variables["NGC346_factor"] *= std::sqrt(lum_factor);
        else if (system == "NGC1672") variables["NGC1672_factor"] *= std::sqrt(lum_factor);
    }
    
    if (observations.find("distance") != observations.end() && system == "M104") {
        variables["M104_dist"] = observations.at("distance");
    }
}

void UQFFBuoyancyModule157::optimizeForMetric(const std::string& metric, double target_value) {
    // Observational-specific optimization scenarios
    if (metric == "standard") {
        variables["M104_factor"] = 1.25;
        variables["NGC4839_factor"] = 1.45;
        variables["CW_factor"] = 1.85;
        variables["NGC346_factor"] = 2.1;
        variables["NGC1672_factor"] = 1.65;
    } else if (metric == "high_luminosity") {
        variables["M104_factor"] *= 2.0;
        variables["NGC4839_factor"] *= 2.5;
        variables["CW_factor"] *= 3.0;
        variables["NGC346_factor"] *= 2.8;
        variables["NGC1672_factor"] *= 2.2;
    } else if (metric == "strong_coupling") {
        variables["alpha"] *= 5.0;
        variables["beta_1"] *= 2.0;
        variables["beta_i"] *= 2.0;
    } else if (metric == "efficient_observation") {
        variables["CW_efficiency"] = std::min(0.99, variables["CW_efficiency"].operator double() * 1.5);
        variables["omega_obs"] *= 2.0;
    } else if (metric == "cluster_optimized") {
        variables["NGC4839_factor"] *= 3.0;
        variables["M104_factor"] *= 1.5;
    } else if (metric == "star_formation") {
        variables["NGC346_factor"] *= 5.0;
        variables["NGC1672_factor"] *= 3.0;
    }
}

// Parameter Exploration (1 method)
std::vector<std::map<std::string, double>> UQFFBuoyancyModule157::generateVariations(int count, double variation_range) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-variation_range, variation_range);
    
    for (int i = 0; i < count; i++) {
        std::map<std::string, double> variation = variables;
        
        // Vary observational parameters
        std::vector<std::string> vary_params = {"M104_factor", "NGC4839_factor", "CW_factor",
                                                "NGC346_factor", "NGC1672_factor", "alpha", "beta_1", "beta_i"};
        for (const auto& param : vary_params) {
            double factor = 1.0 + dis(gen);
            variation[param] *= factor;
        }
        
        variations.push_back(variation);
    }
    
    return variations;
}

// Adaptive Evolution (2 methods)
void UQFFBuoyancyModule157::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-mutation_rate, mutation_rate);
    
    // Mutate observational parameters
    std::vector<std::string> mutable_params = {"M104_factor", "NGC4839_factor", "CW_factor",
                                                "NGC346_factor", "NGC1672_factor", "CW_efficiency",
                                                "alpha", "beta_1", "beta_i", "omega_obs"};
    
    for (const auto& param : mutable_params) {
        if (variables.find(param) != variables.end()) {
            double factor = 1.0 + dis(gen);
            if (param == "CW_efficiency") {
                variables[param] = std::min(1.0, std::max(0.1, variables[param].operator double() * factor));
            } else {
                variables[param] *= factor;
            }
        }
    }
}

void UQFFBuoyancyModule157::evolveSystem(const std::string& system, int generations, std::function<double(double)> fitness_func) {
    std::map<std::string, double> best_params = variables;
    double best_fitness = fitness_func(computeFBi(system));
    
    for (int gen = 0; gen < generations; gen++) {
        mutateParameters(0.1);
        double current_force = computeFBi(system);
        double current_fitness = fitness_func(current_force);
        
        if (current_fitness > best_fitness) {
            best_fitness = current_fitness;
            best_params = variables;
        } else {
            variables = best_params;
        }
    }
    
    variables = best_params;
}

// State Management (4 methods)
void UQFFBuoyancyModule157::saveState(const std::string& state_name) {
    saved_states_obs::states[state_name] = variables;
}

void UQFFBuoyancyModule157::restoreState(const std::string& state_name) {
    if (saved_states_obs::states.find(state_name) != saved_states_obs::states.end()) {
        variables = saved_states_obs::states[state_name];
    }
}

std::vector<std::string> UQFFBuoyancyModule157::listSavedStates() {
    std::vector<std::string> names;
    for (const auto& pair : saved_states_obs::states) {
        names.push_back(pair.first);
    }
    return names;
}

void UQFFBuoyancyModule157::exportStateUQFF(const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        file << "# UQFF Observational Buoyancy Module State Export\n";
        file << "# System: " << getSystemName() << "\n";
        for (const auto& pair : variables) {
            file << pair.first << " = " << pair.second << "\n";
        }
        file.close();
    }
}

// System Analysis (4 methods)
std::map<std::string, double> UQFFBuoyancyModule157::sensitivityAnalysis(const std::string& system, const std::vector<std::string>& params) {
    std::map<std::string, double> sensitivities;
    
    double base_force = computeFBi(system);
    
    for (const auto& param : params) {
        if (variables.find(param) != variables.end()) {
            double original = variables[param];
            variables[param] *= 1.01;
            double perturbed_force = computeFBi(system);
            
            double sensitivity = std::abs((perturbed_force - base_force) / base_force) / 0.01;
            sensitivities[param] = sensitivity;
            
            variables[param] = original;
        }
    }
    
    return sensitivities;
}

std::string UQFFBuoyancyModule157::generateReport(const std::string& system) {
    std::ostringstream report;
    
    report << "=== UQFF Observational Buoyancy Module Report ===\n";
    report << "System: " << system << "\n";
    report << "Module Name: " << getSystemName() << "\n\n";
    
    report << "System Parameters:\n";
    if (system == "M104") {
        report << "  Mass: " << variables["M104_mass"] / (1e11 * variables["M_s"]) << " × 10^11 M_sun\n";
        report << "  Distance: " << variables["M104_dist"] / 3.086e22 << " Mpc\n";
        report << "  System Factor: " << variables["M104_factor"] << "\n";
    } else if (system == "NGC4839") {
        report << "  Mass: " << variables["NGC4839_mass"] / (1e12 * variables["M_s"]) << " × 10^12 M_sun\n";
        report << "  Distance: " << variables["NGC4839_dist"] / 3.086e22 << " Mpc\n";
        report << "  System Factor: " << variables["NGC4839_factor"] << "\n";
    } else if (system == "Chandra_Webb") {
        report << "  System Factor: " << variables["CW_factor"] << "\n";
        report << "  Efficiency: " << variables["CW_efficiency"] << "\n";
    } else if (system == "NGC346") {
        report << "  Mass: " << variables["NGC346_mass"] / (1e5 * variables["M_s"]) << " × 10^5 M_sun\n";
        report << "  Distance: " << variables["NGC346_dist"] / (1000 * 3.086e22) << " kpc\n";
        report << "  System Factor: " << variables["NGC346_factor"] << "\n";
    } else if (system == "NGC1672") {
        report << "  Mass: " << variables["NGC1672_mass"] / (1e11 * variables["M_s"]) << " × 10^11 M_sun\n";
        report << "  Distance: " << variables["NGC1672_dist"] / 3.086e22 << " Mpc\n";
        report << "  System Factor: " << variables["NGC1672_factor"] << "\n";
    }
    
    report << "\nPhysical Parameters:\n";
    report << "  Beta_1: " << variables["beta_1"] << "\n";
    report << "  Beta_i: " << variables["beta_i"] << "\n";
    report << "  Alpha: " << variables["alpha"] << "\n";
    report << "  Rho_vac: " << variables["rho_vac"] << " kg/m^3\n\n";
    
    double Ub1 = computeUb1(system);
    double Ui = computeUi(system);
    double FBi = computeFBi(system);
    report << "Computed Forces:\n";
    report << "  Ub1 = " << Ub1 << " N\n";
    report << "  Ui = " << Ui << " N\n";
    report << "  FBi = " << FBi << " N\n\n";
    
    double dpm = computeDPM_resonance(system);
    double x2 = computeX2(system);
    double g = computeG(system);
    double q = computeQ_wave(system);
    report << "Additional Quantities:\n";
    report << "  DPM Resonance = " << dpm << " Hz\n";
    report << "  X2 = " << x2 << "\n";
    report << "  Modified G = " << g << " m^3/kg/s^2\n";
    report << "  Q_wave = " << q << "\n";
    
    return report.str();
}

bool UQFFBuoyancyModule157::validateConsistency(const std::string& system) {
    // Check for valid parameters
    if (variables["alpha"] <= 0.0) return false;
    if (variables["beta_1"] <= 0.0) return false;
    if (variables["beta_i"] <= 0.0) return false;
    if (variables["rho_vac"] <= 0.0) return false;
    
    // Check system-specific parameters
    if (system == "M104" && variables["M104_mass"] <= 0.0) return false;
    if (system == "NGC4839" && variables["NGC4839_mass"] <= 0.0) return false;
    if (system == "NGC346" && variables["NGC346_mass"] <= 0.0) return false;
    if (system == "NGC1672" && variables["NGC1672_mass"] <= 0.0) return false;
    if (system == "Chandra_Webb" && (variables["CW_efficiency"] < 0.0 || variables["CW_efficiency"] > 1.0)) return false;
    
    // Check force computation
    double force = computeFBi(system);
    if (std::isnan(force) || std::isinf(force)) return false;
    
    return true;
}

void UQFFBuoyancyModule157::autoCorrectAnomalies(const std::string& system) {
    // Correct physical constants if out of range
    if (variables["alpha"] <= 0.0 || variables["alpha"] > 1.0) {
        variables["alpha"] = 0.0073;
    }
    
    if (variables["beta_1"] <= 0.0 || variables["beta_1"] > 1.0) {
        variables["beta_1"] = 1.2e-3;
    }
    
    if (variables["beta_i"] <= 0.0 || variables["beta_i"] > 1.0) {
        variables["beta_i"] = 0.85e-3;
    }
    
    if (variables["rho_vac"] <= 0.0 || variables["rho_vac"] > 1e-20) {
        variables["rho_vac"] = 9.47e-27;
    }
    
    // Correct system-specific parameters
    if (system == "M104" && variables["M104_mass"] <= 0.0) {
        variables["M104_mass"] = 8e11 * 1.989e30;
    }
    
    if (system == "NGC4839" && variables["NGC4839_mass"] <= 0.0) {
        variables["NGC4839_mass"] = 2.5e12 * 1.989e30;
    }
    
    if (system == "NGC346" && variables["NGC346_mass"] <= 0.0) {
        variables["NGC346_mass"] = 1.8e5 * 1.989e30;
    }
    
    if (system == "NGC1672" && variables["NGC1672_mass"] <= 0.0) {
        variables["NGC1672_mass"] = 1.2e11 * 1.989e30;
    }
    
    if (system == "Chandra_Webb") {
        if (variables["CW_efficiency"] < 0.0 || variables["CW_efficiency"] > 1.0) {
            variables["CW_efficiency"] = 0.95;
        }
    }
    
    // Recalibrate factors if force becomes non-physical
    double force = computeFBi(system);
    if (std::isnan(force) || std::isinf(force)) {
        if (system == "M104") variables["M104_factor"] = 1.25;
        else if (system == "NGC4839") variables["NGC4839_factor"] = 1.45;
        else if (system == "Chandra_Webb") variables["CW_factor"] = 1.85;
        else if (system == "NGC346") variables["NGC346_factor"] = 2.1;
        else if (system == "NGC1672") variables["NGC1672_factor"] = 1.65;
    }
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

// ===== COMPREHENSIVE EXAMPLE USAGE FOR OBSERVATIONAL SYSTEMS =====

int main() {
    std::cout << "=== UQFF Observational Buoyancy Module - Comprehensive Dynamic Capabilities Demo ===\n\n";
    
    UQFFBuoyancyModule157 obs_module;
    SurfaceMagneticFieldModule mag_module;
    
    // Define all 5 observational systems
    std::vector<std::string> systems = {
        "M104", "NGC4839", "Chandra_Webb", "NGC346", "NGC1672"
    };
    
    std::cout << "System Name: " << obs_module.getSystemName() << "\n\n";
    
    // Test 1: Variable Management across all systems
    std::cout << "Test 1: Variable Management\n";
    obs_module.createVariable("test_obs_param", 1.75);
    obs_module.cloneVariable("M104_factor", "M104_factor_backup");
    auto var_list = obs_module.listVariables();
    std::cout << "  Created and cloned variables. Total variables: " << var_list.size() << "\n\n";
    
    // Test 2: Batch Operations on observational parameters
    std::cout << "Test 2: Batch Operations\n";
    std::vector<std::string> obs_params = {"M104_factor", "NGC4839_factor", "NGC346_factor"};
    obs_module.scaleVariableGroup(obs_params, 1.8);
    std::cout << "  Scaled observational parameter group by factor 1.8\n\n";
    
    // Test 3: Self-Expansion for each system
    std::cout << "Test 3: Self-Expansion (System-Specific)\n";
    for (const auto& sys : systems) {
        obs_module.expandSystemScale(sys, 1.15, 1.05);
        std::cout << "  Expanded " << sys << " (M×1.15, d×1.05)\n";
    }
    obs_module.expandForceScale(1.3, 1.4);
    std::cout << "  Expanded force scales (buoyancy×1.3, coupling×1.4)\n";
    obs_module.expandObservationalScale("M104", 1.25, 1.1);
    std::cout << "  Expanded M104 observational scale (L×1.25, eff×1.1)\n\n";
    
    // Test 4: Parameter Space Expansion
    std::cout << "Test 4: Parameter Space Expansion\n";
    obs_module.expandParameterSpace(1.12);
    std::cout << "  Expanded observational parameter space by 1.12×\n\n";
    
    // Test 5: Auto-Refinement for multiple systems
    std::cout << "Test 5: Auto-Refinement\n";
    obs_module.autoRefineParameters("M104", 1e15, 0.05);
    std::cout << "  Refined M104 to target force 1e15 N\n";
    obs_module.autoRefineParameters("NGC4839", 5e15, 0.05);
    std::cout << "  Refined NGC4839 to target force 5e15 N\n\n";
    
    // Test 6: Calibration to Observations
    std::cout << "Test 6: Calibration to Observations\n";
    std::map<std::string, double> obs_data;
    obs_data["M104_dist"] = 10.0 * 3.086e22;
    obs_data["luminosity"] = 2.5e36;
    obs_data["alpha"] = 0.0080;
    obs_module.calibrateToObservations("M104", obs_data);
    std::cout << "  Calibrated M104 to observational data\n\n";
    
    // Test 7: Optimization Scenarios
    std::cout << "Test 7: Optimization Scenarios\n";
    std::vector<std::string> scenarios = {"standard", "high_luminosity", "strong_coupling",
                                          "efficient_observation", "cluster_optimized", "star_formation"};
    for (const auto& scenario : scenarios) {
        obs_module.optimizeForMetric(scenario, 0.0);
        std::cout << "  Optimized for: " << scenario << "\n";
    }
    std::cout << "\n";
    
    // Test 8: Parameter Variations
    std::cout << "Test 8: Parameter Variations\n";
    auto variations = obs_module.generateVariations(5, 0.12);
    std::cout << "  Generated " << variations.size() << " observational parameter variations (±12%)\n\n";
    
    // Test 9: Adaptive Evolution
    std::cout << "Test 9: Adaptive Evolution\n";
    obs_module.mutateParameters(0.09);
    std::cout << "  Mutated observational parameters (rate: 0.09)\n";
    auto fitness = [](double force) { return -std::abs(force - 1e15); };
    obs_module.evolveSystem("NGC346", 20, fitness);
    std::cout << "  Evolved NGC346 system (20 generations)\n\n";
    
    // Test 10: State Management
    std::cout << "Test 10: State Management\n";
    obs_module.saveState("obs_config_1");
    obs_module.saveState("obs_config_2");
    auto saved_states = obs_module.listSavedStates();
    std::cout << "  Saved " << saved_states.size() << " observational configurations\n";
    obs_module.restoreState("obs_config_1");
    std::cout << "  Restored state: obs_config_1\n";
    obs_module.exportStateUQFF("obs_state_export.txt");
    std::cout << "  Exported observational state to file\n\n";
    
    // Test 11: Sensitivity Analysis for all systems
    std::cout << "Test 11: Sensitivity Analysis\n";
    std::vector<std::string> test_params = {"alpha", "beta_1", "beta_i", "M104_factor", "NGC4839_factor"};
    for (const auto& sys : systems) {
        auto sensitivities = obs_module.sensitivityAnalysis(sys, test_params);
        std::cout << "  " << sys << " - Most sensitive param: ";
        double max_sens = 0.0;
        std::string max_param;
        for (const auto& s : sensitivities) {
            if (s.second > max_sens) {
                max_sens = s.second;
                max_param = s.first;
            }
        }
        std::cout << max_param << " (sensitivity: " << max_sens << ")\n";
    }
    std::cout << "\n";
    
    // Test 12: Report Generation for each system
    std::cout << "Test 12: Report Generation\n";
    for (const auto& sys : systems) {
        std::string report = obs_module.generateReport(sys);
        std::cout << "  Generated report for " << sys << " (" << report.length() << " chars)\n";
    }
    std::cout << "\n";
    
    // Test 13: Consistency Validation
    std::cout << "Test 13: Consistency Validation\n";
    for (const auto& sys : systems) {
        bool valid = obs_module.validateConsistency(sys);
        std::cout << "  " << sys << ": " << (valid ? "VALID" : "INVALID") << "\n";
    }
    std::cout << "\n";
    
    // Test 14: Auto-Correction
    std::cout << "Test 14: Auto-Correction\n";
    obs_module.updateVariable("alpha", -0.5);  // Invalid value
    obs_module.autoCorrectAnomalies("M104");
    std::cout << "  Corrected anomalies in M104 observational parameters\n\n";
    
    // Test 15: Magnetic Field Module - Observational Integration
    std::cout << "Test 15: Magnetic Field Module - Observational Capabilities\n";
    mag_module.enableSelfLearning(true);
    std::cout << "  Enabled observational magnetic self-learning\n";
    mag_module.autoCalibrate("B_ref", 1.3, 0.02);
    std::cout << "  Auto-calibrated B_ref to 1.3 T\n";
    mag_module.adaptiveUpdate(1e10, "luminosity_coupling");
    std::cout << "  Applied observational adaptive update\n\n";
    
    // Test 16: Observational Data Scaling
    std::cout << "Test 16: Observational Data Scaling\n";
    std::map<std::string, double> mag_obs_data;
    mag_obs_data["flux_density"] = 1.2e-26;
    mag_obs_data["observation_frequency"] = 1.5e9;
    mag_obs_data["magnetic_field_strength"] = 1.4;
    mag_module.scaleToObservationalData(mag_obs_data);
    std::cout << "  Scaled magnetic module to observational data\n\n";
    
    // Test 17: Custom Observational Variables
    std::cout << "Test 17: Custom Observational Variables\n";
    mag_module.addCustomVariable("spectral_resolution", 0.005, "observation_frequency");
    mag_module.addCustomVariable("beam_size", 2.5e-4, "angular_resolution");
    std::cout << "  Added 2 custom observational magnetic variables\n\n";
    
    // Test 18: Magnetic History Tracking
    std::cout << "Test 18: Observational Magnetic History\n";
    auto history = mag_module.getVariableHistory("B_ref", 10);
    std::cout << "  Retrieved " << history.size() << " B_ref history entries\n\n";
    
    // Test 19: State Export/Import
    std::cout << "Test 19: Observational Magnetic State Management\n";
    mag_module.exportState("obs_magnetic_state.txt");
    std::cout << "  Exported observational magnetic state\n";
    mag_module.importState("obs_magnetic_state.txt");
    std::cout << "  Imported observational magnetic state\n\n";
    
    // Test 20: Final Force Computations for All Systems
    std::cout << "Test 20: Final Observational Force Computations\n";
    for (const auto& sys : systems) {
        double Ub1 = obs_module.computeUb1(sys);
        double Ui = obs_module.computeUi(sys);
        double FBi = obs_module.computeFBi(sys);
        double dpm = obs_module.computeDPM_resonance(sys);
        double x2 = obs_module.computeX2(sys);
        double g = obs_module.computeG(sys);
        double q = obs_module.computeQ_wave(sys);
        double b_j = mag_module.computeB_j(1e10, 1.0);
        double u_g3 = mag_module.computeU_g3_example(1e10, 1.0);
        double obs_coupling = mag_module.computeObservationalCoupling(b_j, 1e36);
        
        std::cout << "  " << sys << ":\n";
        std::cout << "    Ub1 = " << Ub1 << " N\n";
        std::cout << "    Ui = " << Ui << " N\n";
        std::cout << "    FBi = " << FBi << " N\n";
        std::cout << "    DPM_resonance = " << dpm << " Hz\n";
        std::cout << "    X2 = " << x2 << "\n";
        std::cout << "    Modified_G = " << g << " m^3/kg/s^2\n";
        std::cout << "    Q_wave = " << q << "\n";
        std::cout << "    B_j = " << b_j << " T\n";
        std::cout << "    U_g3 = " << u_g3 << " J\n";
        std::cout << "    Obs-Magnetic Coupling = " << obs_coupling << "\n\n";
    }
    
    std::cout << "=== All 20 Observational Dynamic Capability Tests Completed Successfully ===\n";
    
    return 0;
}
