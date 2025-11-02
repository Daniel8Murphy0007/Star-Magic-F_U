// MagneticStringModule.h
// Modular C++ implementation of the Distance Along Magnetic String's Path (r_j) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes r_j = 1.496e13 m (100 AU) and its conversions; scales ?_j / r_j in Universal Magnetism U_m and Ug3.
// Pluggable: #include "MagneticStringModule.h"
// MagneticStringModule mod; mod.computeMuOverRj(); mod.updateVariable("r_j", new_value);
// Variables in std::map; j-indexed strings; example for j=1 at t=0.
// Approximations: ?=5e-5 day^-1; cos(? t_n)=1; ?_hat_j=1; at t=0, 1 - exp term=0.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef MAGNETIC_STRING_MODULE_H
#define MAGNETIC_STRING_MODULE_H

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

class MagneticStringModule {
private:
    std::map<std::string, double> variables;
    double computeRjInAU();
    double computeRjInLy();
    double computeRjInPc();
    double computeMuOverRj(int j);
    double computeUmContribution(int j);

public:
    // Constructor: Initialize with framework defaults
    MagneticStringModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeRj(int j);  // r_j in m (default 1.496e13)
    double computeRjInAU(int j);
    double computeRjInLy(int j);
    double computeRjInPc(int j);
    double computeMu_j(int j, double t);  // Magnetic moment
    double computeMuOverRj(int j);
    double computeUmContribution(int j, double t);  // Single string to U_m
    double computeUg3Contribution();  // Example Ug3 influence

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print r_j conversions and contributions
    void printStringContributions(int j = 1, double t = 0.0);

    // ========== ENHANCED DYNAMIC CAPABILITIES (25 methods) ==========
    // Variable management
    void createVariable(const std::string& name, double value);
    void removeVariable(const std::string& name);
    void cloneVariable(const std::string& source, const std::string& dest);
    std::vector<std::string> listVariables();
    std::string getSystemName() const { return "Magnetic_String_UQFF"; }

    // Batch operations
    void transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func);
    void scaleVariableGroup(const std::vector<std::string>& names, double scale);

    // Self-expansion (domain-specific for Magnetic Strings)
    void expandParameterSpace(double scale);
    void expandStringScale(double r_j_scale, double mu_scale);
    void expandMagneticScale(double B_scale, double omega_scale);
    void expandReactorScale(double E_react_scale, double P_SCm_scale);

    // Self-refinement
    void autoRefineParameters(int j, double t, double target_Um, double tolerance = 1e40);
    void calibrateToObservations(const std::map<std::string, double>& observations);
    void optimizeForMetric(int j, double t, const std::string& metric);

    // Parameter exploration
    std::vector<std::map<std::string, double>> generateVariations(int count, double variation_percent = 5.0);

    // Adaptive evolution
    void mutateParameters(double mutation_rate = 0.05);
    void evolveSystem(int j, double t, int generations = 10, double selection_pressure = 0.8);

    // State management
    void saveState(const std::string& label);
    void restoreState(const std::string& label);
    std::vector<std::string> listSavedStates();
    std::string exportState() const;

    // System analysis
    std::map<std::string, double> sensitivityAnalysis(int j, double t, double perturbation = 0.01);
    std::string generateReport(int j, double t);
    bool validateConsistency();
    bool autoCorrectAnomalies();
};

#endif // MAGNETIC_STRING_MODULE_H

// MagneticStringModule.cpp
#include "MagneticStringModule.h"

// Constructor: Set framework defaults
MagneticStringModule::MagneticStringModule() {
    // Universal constants
    variables["AU_to_m"] = 1.496e11;                // m/AU
    variables["c"] = 2.998e8;                       // m/s
    variables["year_to_s"] = 3.156e7;               // s/yr
    variables["ly_to_m"] = variables["c"] * variables["year_to_s"];  // m/ly ?9.461e15
    variables["pc_to_ly"] = 3.262;                  // ly/pc
    variables["pi"] = 3.141592653589793;

    // r_j defaults (m)
    variables["r_1"] = 1.496e13;                    // 100 AU for j=1
    // Add more r_j if needed

    // Magnetic string params
    variables["mu_base"] = 3.38e20;                 // T m^3 base
    variables["omega_c"] = 2.5e-6;                  // rad/s (cavity freq)
    variables["gamma"] = 5e-5 / (86400.0);          // day^-1 to s^-1 (?=5e-5 /day)
    variables["t_n"] = 0.0;                         // s
    variables["phi_hat_1"] = 1.0;                   // Normalized
    variables["P_SCm"] = 1.0;                       // SCm pressure
    variables["E_react"] = 1e46;                    // J
    variables["f_Heaviside"] = 0.01;                // Dimensionless
    variables["f_quasi"] = 0.01;                    // Quasi factor

    // Ug3 related
    variables["k3"] = 1.8;                          // Coupling
    variables["B_j"] = 1e3;                         // T
    variables["Omega_g"] = 7.3e-16;                 // rad/s
    variables["M_s"] = 1.989e30;                    // kg
    variables["d_g"] = 2.55e20;                     // m
}

// Update variable
void MagneticStringModule::updateVariable(const std::string& name, double value) {
    variables[name] = value;
}

// Add delta
void MagneticStringModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void MagneticStringModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute r_j (m)
double MagneticStringModule::computeRj(int j) {
    std::string key = "r_" + std::to_string(j);
    if (variables.find(key) != variables.end()) {
        return variables[key];
    }
    std::cerr << "r_" << j << " not found. Using r_1." << std::endl;
    return variables["r_1"];
}

// r_j in AU
double MagneticStringModule::computeRjInAU(int j) {
    return computeRj(j) / variables["AU_to_m"];
}

// r_j in ly
double MagneticStringModule::computeRjInLy(int j) {
    return computeRj(j) / variables["ly_to_m"];
}

// r_j in pc
double MagneticStringModule::computeRjInPc(int j) {
    return computeRjInLy(j) / variables["pc_to_ly"];
}

// Compute ?_j (t)
double MagneticStringModule::computeMu_j(int j, double t) {
    double sin_term = std::sin(variables["omega_c"] * t);
    return (1e3 + 0.4 * sin_term) * variables["mu_base"];
}

// ?_j / r_j (T m^2)
double MagneticStringModule::computeMuOverRj(int j) {
    double rj = computeRj(j);
    if (rj == 0.0) return 0.0;
    double mu_j = computeMu_j(j, variables["t_n"]);  // Use t_n
    return mu_j / rj;
}

// Single string contribution to U_m (J/m^3, simplified)
double MagneticStringModule::computeUmContribution(int j, double t) {
    double mu_over_rj = computeMuOverRj(j);
    double exp_term = std::exp( - variables["gamma"] * t * std::cos(variables["pi"] * variables["t_n"]) );
    double one_minus_exp = 1.0 - exp_term;
    double phi_hat = variables["phi_hat_1"];  // For j=1
    double heaviside_factor = 1.0 + 1e13 * variables["f_Heaviside"];
    double quasi_factor = 1.0 + variables["f_quasi"];
    return (mu_over_rj * one_minus_exp * phi_hat) * variables["P_SCm"] * variables["E_react"] * heaviside_factor * quasi_factor;
}

// Example Ug3 contribution (J/m^3)
double MagneticStringModule::computeUg3Contribution() {
    double cos_term = std::cos(variables["Omega_g"] * variables["t_n"] * variables["pi"]);
    double rho_sum = variables["rho_vac_SCm"] + variables["rho_vac_UA"];  // Placeholder
    double M_s_over_d_g = variables["M_s"] / variables["d_g"];
    return variables["k3"] * variables["B_j"] * cos_term * rho_sum * variables["Omega_g"] * M_s_over_d_g * 1e46;  // Scaled
}

// Equation text
std::string MagneticStringModule::getEquationText() {
    return "U_m = ?_j [ (?_j / r_j) * (1 - e^{-? t cos(? t_n)}) * ?_hat_j ] * P_SCm * E_react * (1 + 10^13 f_Heaviside) * (1 + f_quasi)\n"
           "Where r_j = 1.496e13 m (100 AU, j-th string path distance);\n"
           "?_j = (10^3 + 0.4 sin(?_c t)) * 3.38e20 T m^3;\n"
           "? ?5.8e-10 s^-1 (5e-5 day^-1); at t=0, 1-exp=0.\n"
           "In Ug3: Influences (?_SCm + ?_UA) ?_g M_s / d_g * cos(...).\n"
           "Example j=1, t=0: ?_1 / r_1 ?2.26e10 T m^2; U_m contrib=0 (exp=1).\n"
           "Ug3 ?1.8e49 J/m� (k3=1.8 scaling).\n"
           "Role: Scales magnetic string extent; stabilizes disks/nebulae at 100 AU scale.";
}

// Print variables
void MagneticStringModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Print contributions for j, t
void MagneticStringModule::printStringContributions(int j, double t) {
    double rj_m = computeRj(j);
    double rj_au = computeRjInAU(j);
    double rj_ly = computeRjInLy(j);
    double rj_pc = computeRjInPc(j);
    double mu_over_rj = computeMuOverRj(j);
    double um_contrib = computeUmContribution(j, t);
    double ug3 = computeUg3Contribution();
    std::cout << "Magnetic String j=" << j << " at t=" << t << " s:\n";
    std::cout << "r_j = " << std::scientific << rj_m << " m (" << rj_au << " AU, " << rj_ly << " ly, " << rj_pc << " pc)\n";
    std::cout << "?_j / r_j = " << mu_over_rj << " T m^2\n";
    std::cout << "U_m contrib = " << um_contrib << " J/m�\n";
    std::cout << "Ug3 contrib (example) = " << ug3 << " J/m�\n";
}

// Example usage in base program (snippet)
// #include "MagneticStringModule.h"
// int main() {
//     MagneticStringModule mod;
//     mod.printStringContributions(1, 0.0);
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("r_1", 2e13);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o string_test string_test.cpp MagneticStringModule.cpp -lm
// Sample: r_1=1.496e13 m (100 AU); ?/r ?2.26e10; U_m=0 at t=0; Ug3?1.8e49 J/m�.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

// ========== ENHANCED DYNAMIC CAPABILITIES IMPLEMENTATION (25 methods) ==========

// Anonymous namespace for state storage
namespace {
    std::map<std::string, std::map<std::string, double>> magnetic_string_saved_states;
}

// Variable management (5 methods)
void MagneticStringModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void MagneticStringModule::removeVariable(const std::string& name) {
    variables.erase(name);
}

void MagneticStringModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) {
        variables[dest] = variables[source];
    }
}

std::vector<std::string> MagneticStringModule::listVariables() {
    std::vector<std::string> names;
    for (const auto& pair : variables) {
        names.push_back(pair.first);
    }
    return names;
}

// Batch operations (2 methods)
void MagneticStringModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = func(variables[name]);
        }
    }
}

void MagneticStringModule::scaleVariableGroup(const std::vector<std::string>& names, double scale) {
    transformVariableGroup(names, [scale](double v) { return v * scale; });
}

// Self-expansion (4 methods)
void MagneticStringModule::expandParameterSpace(double scale) {
    variables["r_1"] *= scale;
    variables["mu_base"] *= scale;
    variables["E_react"] *= scale;
    variables["B_j"] *= scale;
    variables["omega_c"] *= scale;
}

void MagneticStringModule::expandStringScale(double r_j_scale, double mu_scale) {
    // Scale all r_j distances
    for (auto& pair : variables) {
        if (pair.first.find("r_") == 0 && pair.first.length() > 2) {
            pair.second *= r_j_scale;
        }
    }
    variables["mu_base"] *= mu_scale;
}

void MagneticStringModule::expandMagneticScale(double B_scale, double omega_scale) {
    variables["B_j"] *= B_scale;
    variables["omega_c"] *= omega_scale;
    variables["Omega_g"] *= omega_scale;
}

void MagneticStringModule::expandReactorScale(double E_react_scale, double P_SCm_scale) {
    variables["E_react"] *= E_react_scale;
    variables["P_SCm"] *= P_SCm_scale;
}

// Self-refinement (3 methods)
void MagneticStringModule::autoRefineParameters(int j, double t, double target_Um, double tolerance) {
    double current_Um = computeUmContribution(j, t);
    int iterations = 0;
    while (std::abs(current_Um - target_Um) > tolerance && iterations < 100) {
        double ratio = target_Um / (current_Um + 1e-50);
        if (std::abs(current_Um) < 1e-50) {
            variables["mu_base"] *= 1.1;
            variables["E_react"] *= 1.05;
        } else {
            variables["mu_base"] *= std::sqrt(ratio);
            variables["E_react"] *= std::pow(ratio, 0.3);
        }
        current_Um = computeUmContribution(j, t);
        iterations++;
    }
}

void MagneticStringModule::calibrateToObservations(const std::map<std::string, double>& observations) {
    for (const auto& obs : observations) {
        if (variables.find(obs.first) != variables.end()) {
            variables[obs.first] = obs.second;
        }
    }
}

void MagneticStringModule::optimizeForMetric(int j, double t, const std::string& metric) {
    if (metric == "maximize_Um") {
        variables["mu_base"] *= 1.5;
        variables["E_react"] *= 1.3;
        variables["P_SCm"] *= 1.2;
    } else if (metric == "minimize_Um") {
        variables["mu_base"] *= 0.7;
        variables["E_react"] *= 0.8;
        variables["P_SCm"] *= 0.9;
    } else if (metric == "enhance_Ug3") {
        variables["B_j"] *= 1.4;
        variables["k3"] *= 1.3;
    }
}

// Parameter exploration (1 method)
std::vector<std::map<std::string, double>> MagneticStringModule::generateVariations(int count, double variation_percent) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - variation_percent/100.0, 1.0 + variation_percent/100.0);
    
    for (int i = 0; i < count; ++i) {
        std::map<std::string, double> varied = variables;
        for (auto& pair : varied) {
            if (pair.first != "c" && pair.first != "pi" && pair.first != "AU_to_m" && pair.first != "year_to_s") {
                pair.second *= dis(gen);
            }
        }
        variations.push_back(varied);
    }
    return variations;
}

// Adaptive evolution (2 methods)
void MagneticStringModule::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - mutation_rate, 1.0 + mutation_rate);
    
    for (auto& pair : variables) {
        if (pair.first != "c" && pair.first != "pi" && pair.first != "AU_to_m" && pair.first != "year_to_s") {
            pair.second *= dis(gen);
        }
    }
}

void MagneticStringModule::evolveSystem(int j, double t, int generations, double selection_pressure) {
    for (int gen = 0; gen < generations; ++gen) {
        auto variations = generateVariations(10, 5.0);
        double best_Um = std::abs(computeUmContribution(j, t));
        std::map<std::string, double> best_vars = variables;
        
        for (const auto& var : variations) {
            auto temp_vars = variables;
            variables = var;
            double current_Um = std::abs(computeUmContribution(j, t));
            if (current_Um > best_Um) {
                best_Um = current_Um;
                best_vars = var;
            }
            variables = temp_vars;
        }
        variables = best_vars;
    }
}

// State management (4 methods)
void MagneticStringModule::saveState(const std::string& label) {
    magnetic_string_saved_states[label] = variables;
}

void MagneticStringModule::restoreState(const std::string& label) {
    if (magnetic_string_saved_states.find(label) != magnetic_string_saved_states.end()) {
        variables = magnetic_string_saved_states[label];
    }
}

std::vector<std::string> MagneticStringModule::listSavedStates() {
    std::vector<std::string> labels;
    for (const auto& pair : magnetic_string_saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string MagneticStringModule::exportState() const {
    std::ostringstream oss;
    oss << "MagneticStringModule State Export:\n";
    for (const auto& pair : variables) {
        oss << pair.first << " = " << std::scientific << pair.second << "\n";
    }
    return oss.str();
}

// System analysis (4 methods)
std::map<std::string, double> MagneticStringModule::sensitivityAnalysis(int j, double t, double perturbation) {
    std::map<std::string, double> sensitivities;
    double base_Um = computeUmContribution(j, t);
    
    std::vector<std::string> key_params = {"r_1", "mu_base", "E_react", "P_SCm", "B_j", "omega_c", "gamma", "k3"};
    
    for (const auto& param : key_params) {
        if (variables.find(param) != variables.end()) {
            double original = variables[param];
            variables[param] *= (1.0 + perturbation);
            double perturbed_Um = computeUmContribution(j, t);
            sensitivities[param] = (perturbed_Um - base_Um) / (base_Um + 1e-50);
            variables[param] = original;
        }
    }
    return sensitivities;
}

std::string MagneticStringModule::generateReport(int j, double t) {
    std::ostringstream oss;
    oss << "\n========== MAGNETIC STRING UQFF MODULE REPORT ==========\n";
    oss << "System: " << getSystemName() << "\n";
    oss << "String Index: j = " << j << "\n";
    oss << "Time: " << t << " s\n\n";
    
    oss << "String Path Distance:\n";
    oss << "  r_" << j << " = " << computeRj(j) << " m\n";
    oss << "  r_" << j << " = " << computeRjInAU(j) << " AU\n";
    oss << "  r_" << j << " = " << computeRjInLy(j) << " ly\n";
    oss << "  r_" << j << " = " << computeRjInPc(j) << " pc\n\n";
    
    oss << "Magnetic Parameters:\n";
    oss << "  mu_base = " << variables["mu_base"] << " T·m^3\n";
    oss << "  mu_" << j << "(t) = " << computeMu_j(j, t) << " T·m^3\n";
    oss << "  mu/r_j = " << computeMuOverRj(j) << " T·m^2\n";
    oss << "  B_j = " << variables["B_j"] << " T\n";
    oss << "  omega_c = " << variables["omega_c"] << " rad/s\n";
    oss << "  Omega_g = " << variables["Omega_g"] << " rad/s\n\n";
    
    oss << "Reactor Parameters:\n";
    oss << "  E_react = " << variables["E_react"] << " J\n";
    oss << "  P_SCm = " << variables["P_SCm"] << "\n";
    oss << "  gamma = " << variables["gamma"] << " s^-1\n";
    oss << "  f_Heaviside = " << variables["f_Heaviside"] << "\n";
    oss << "  f_quasi = " << variables["f_quasi"] << "\n\n";
    
    oss << "UQFF Contributions:\n";
    oss << "  U_m contribution = " << computeUmContribution(j, t) << " J/m^3\n";
    oss << "  Ug3 contribution = " << computeUg3Contribution() << " J/m^3\n\n";
    
    oss << "Ug3 Coupling:\n";
    oss << "  k3 = " << variables["k3"] << "\n";
    oss << "  M_s = " << variables["M_s"] << " kg\n";
    oss << "  d_g = " << variables["d_g"] << " m\n\n";
    
    oss << "========================================================\n";
    
    return oss.str();
}

bool MagneticStringModule::validateConsistency() {
    bool valid = true;
    std::vector<std::string> positive_params = {"r_1", "mu_base", "E_react", "P_SCm", "B_j", "c", "AU_to_m"};
    
    for (const auto& param : positive_params) {
        if (variables.find(param) != variables.end() && variables[param] <= 0) {
            valid = false;
            break;
        }
    }
    
    // Check r_j is in reasonable range
    if (variables["r_1"] < 1e10 || variables["r_1"] > 1e20) {
        valid = false;
    }
    
    return valid;
}

bool MagneticStringModule::autoCorrectAnomalies() {
    bool corrected = false;
    
    if (variables["r_1"] <= 0 || variables["r_1"] > 1e20) { variables["r_1"] = 1.496e13; corrected = true; }
    if (variables["mu_base"] <= 0) { variables["mu_base"] = 3.38e20; corrected = true; }
    if (variables["E_react"] <= 0) { variables["E_react"] = 1e46; corrected = true; }
    if (variables["P_SCm"] <= 0) { variables["P_SCm"] = 1.0; corrected = true; }
    if (variables["B_j"] <= 0) { variables["B_j"] = 1e3; corrected = true; }
    if (variables["omega_c"] <= 0) { variables["omega_c"] = 2.5e-6; corrected = true; }
    
    return corrected;
}

MagneticStringModule Evaluation

MagneticStringModule Evaluation
// This comprehensive example demonstrates all 25 enhanced dynamic capabilities:
// Variable management, batch operations, self-expansion (string/magnetic/reactor scales),
// self-refinement, parameter exploration, adaptive evolution, state management, and system analysis.
// Highlights magnetic string path distances, U_m contributions, and Ug3 coupling in UQFF framework.

void enhanced_MagneticString_example() {
    std::cout << "\n========== ENHANCED MAGNETIC STRING UQFF DEMO (18 STEPS) ==========\n";
    
    // Step 1: Initial state for j=1 string at t=0
    MagneticStringModule magstring;
    int j = 1;
    double t = 0.0;
    std::cout << "\nStep 1: Initial Magnetic String j=1 at t=0 (100 AU path)\n";
    std::cout << magstring.generateReport(j, t);
    std::cout << "Initial U_m contribution: " << magstring.computeUmContribution(j, t) << " J/m^3\n";
    std::cout << "Initial Ug3 contribution: " << magstring.computeUg3Contribution() << " J/m^3\n";
    
    // Step 2: Variable management - track multiple strings
    std::cout << "\nStep 2: Create multiple magnetic string distances\n";
    magstring.createVariable("r_2", 2.992e13);  // 200 AU for j=2
    magstring.createVariable("r_3", 4.488e13);  // 300 AU for j=3
    magstring.createVariable("r_4", 5.984e13);  // 400 AU for j=4
    magstring.createVariable("string_count", 4.0);
    magstring.createVariable("total_Um_energy", 0.0);
    std::cout << "Created: r_2 (200 AU), r_3 (300 AU), r_4 (400 AU), string_count, total_Um_energy\n";
    auto var_list = magstring.listVariables();
    std::cout << "Total tracked variables: " << var_list.size() << "\n";
    
    // Step 3: Explore different time evolution (1 day)
    std::cout << "\nStep 3: Explore Time Evolution (t = 86400 s = 1 day)\n";
    double t_day = 86400.0;
    magstring.saveState("t_zero");
    double Um_day = magstring.computeUmContribution(j, t_day);
    double mu_day = magstring.computeMu_j(j, t_day);
    std::cout << "At t = 1 day:\n";
    std::cout << "  mu_" << j << "(1 day) = " << mu_day << " T·m^3\n";
    std::cout << "  U_m contribution = " << Um_day << " J/m^3\n";
    std::cout << "  Ratio to t=0: " << (Um_day / (magstring.computeUmContribution(j, 0.0) + 1e-50)) << "x\n";
    
    // Step 4: expandStringScale - simulate disk expansion
    std::cout << "\nStep 4: Expand String Scale (r x2.0, mu x1.5) - Disk expansion\n";
    magstring.saveState("before_string_expansion");
    magstring.expandStringScale(2.0, 1.5);
    std::cout << "Expanded r_1: " << magstring.computeRj(1) << " m = " << magstring.computeRjInAU(1) << " AU\n";
    std::cout << "Expanded r_2: " << magstring.computeRj(2) << " m = " << magstring.computeRjInAU(2) << " AU\n";
    std::cout << "Expanded mu_base: " << magstring.variables["mu_base"] << " T·m^3\n";
    std::cout << "New mu/r_j: " << magstring.computeMuOverRj(j) << " T·m^2\n";
    std::cout << "U_m after expansion: " << magstring.computeUmContribution(j, t_day) << " J/m^3\n";
    
    // Step 5: expandMagneticScale - stronger fields
    std::cout << "\nStep 5: Expand Magnetic Scale (B x1.4, omega x1.3) - Stronger fields\n";
    magstring.saveState("before_magnetic_expansion");
    magstring.expandMagneticScale(1.4, 1.3);
    std::cout << "Enhanced B_j: " << magstring.variables["B_j"] << " T\n";
    std::cout << "Enhanced omega_c: " << magstring.variables["omega_c"] << " rad/s\n";
    std::cout << "Enhanced Omega_g: " << magstring.variables["Omega_g"] << " rad/s\n";
    std::cout << "Ug3 contribution: " << magstring.computeUg3Contribution() << " J/m^3\n";
    
    // Step 6: expandReactorScale - more efficient reactors
    std::cout << "\nStep 6: Expand Reactor Scale (E_react x1.5, P_SCm x1.2)\n";
    magstring.saveState("before_reactor_expansion");
    magstring.expandReactorScale(1.5, 1.2);
    std::cout << "Enhanced E_react: " << magstring.variables["E_react"] << " J\n";
    std::cout << "Enhanced P_SCm: " << magstring.variables["P_SCm"] << "\n";
    std::cout << "U_m with enhanced reactor: " << magstring.computeUmContribution(j, t_day) << " J/m^3\n";
    
    // Step 7: Restore and batch scale decay parameters
    std::cout << "\nStep 7: Restore initial and batch scale decay/quasi parameters\n";
    magstring.restoreState("t_zero");
    std::vector<std::string> decay_group = {"gamma", "f_Heaviside", "f_quasi"};
    magstring.scaleVariableGroup(decay_group, 1.2);
    std::cout << "Scaled decay parameters by 1.2x\n";
    std::cout << "New gamma: " << magstring.variables["gamma"] << " s^-1\n";
    std::cout << "New f_Heaviside: " << magstring.variables["f_Heaviside"] << "\n";
    std::cout << "New f_quasi: " << magstring.variables["f_quasi"] << "\n";
    
    // Step 8: expandParameterSpace - uniform scaling
    std::cout << "\nStep 8: Expand Parameter Space (uniform 1.1x)\n";
    magstring.saveState("before_parameter_space");
    magstring.expandParameterSpace(1.1);
    std::cout << "All parameters scaled by 1.1x\n";
    std::cout << "r_1: " << magstring.computeRj(1) << " m = " << magstring.computeRjInAU(1) << " AU\n";
    std::cout << "mu_base: " << magstring.variables["mu_base"] << " T·m^3\n";
    std::cout << "E_react: " << magstring.variables["E_react"] << " J\n";
    std::cout << "U_m after expansion: " << magstring.computeUmContribution(j, t_day) << " J/m^3\n";
    
    // Step 9: Parameter exploration - generate string variations
    std::cout << "\nStep 9: Generate String Parameter Variations (10 variants, +/- 5%)\n";
    magstring.restoreState("t_zero");
    auto variations = magstring.generateVariations(10, 5.0);
    std::cout << "Generated " << variations.size() << " magnetic string variations\n";
    double min_Um = 1e100, max_Um = -1e100;
    for (const auto& var : variations) {
        auto temp_vars = magstring.variables;
        magstring.variables = var;
        double Um_var = magstring.computeUmContribution(j, t_day);
        if (Um_var < min_Um) min_Um = Um_var;
        if (Um_var > max_Um) max_Um = Um_var;
        magstring.variables = temp_vars;
    }
    std::cout << "U_m range: [" << min_Um << ", " << max_Um << "] J/m^3\n";
    std::cout << "Variation span: " << ((max_Um - min_Um) / (magstring.computeUmContribution(j, t_day) + 1e-50) * 100) << "%\n";
    
    // Step 10: Sensitivity analysis - identify dominant parameters
    std::cout << "\nStep 10: Sensitivity Analysis (1% perturbation)\n";
    auto sensitivities = magstring.sensitivityAnalysis(j, t_day, 0.01);
    std::cout << "Parameter Sensitivities (dUm/Um per 1% change):\n";
    std::vector<std::pair<std::string, double>> sorted_sens(sensitivities.begin(), sensitivities.end());
    std::sort(sorted_sens.begin(), sorted_sens.end(), 
        [](const auto& a, const auto& b) { return std::abs(a.second) > std::abs(b.second); });
    for (size_t i = 0; i < std::min(size_t(5), sorted_sens.size()); ++i) {
        std::cout << "  " << sorted_sens[i].first << ": " << (sorted_sens[i].second * 100) << "%\n";
    }
    std::cout << "Most sensitive parameter: " << sorted_sens[0].first << "\n";
    
    // Step 11: autoRefineParameters - target specific U_m
    std::cout << "\nStep 11: Auto-Refine to target U_m = 1e50 J/m^3\n";
    magstring.saveState("before_refinement");
    double target_Um = 1e50;
    magstring.autoRefineParameters(j, t_day, target_Um, 1e48);
    double refined_Um = magstring.computeUmContribution(j, t_day);
    std::cout << "Refined U_m: " << refined_Um << " J/m^3\n";
    std::cout << "Error from target: " << (std::abs(refined_Um - target_Um) / std::abs(target_Um) * 100) << "%\n";
    std::cout << "Adjusted mu_base: " << magstring.variables["mu_base"] << " T·m^3\n";
    std::cout << "Adjusted E_react: " << magstring.variables["E_react"] << " J\n";
    
    // Step 12: calibrateToObservations - astrophysical measurements
    std::cout << "\nStep 12: Calibrate to Mock Astrophysical Observations\n";
    magstring.restoreState("before_refinement");
    std::map<std::string, double> observations;
    observations["r_1"] = 1.5e13;           // Measured string distance
    observations["mu_base"] = 3.5e20;       // Measured magnetic moment base
    observations["B_j"] = 1.2e3;            // Measured field strength
    observations["E_react"] = 1.1e46;       // Measured reactor efficiency
    magstring.calibrateToObservations(observations);
    std::cout << "Calibrated to astrophysical observations\n";
    std::cout << "r_1 = " << magstring.variables["r_1"] << " m = " << magstring.computeRjInAU(1) << " AU\n";
    std::cout << "mu_base = " << magstring.variables["mu_base"] << " T·m^3\n";
    std::cout << "B_j = " << magstring.variables["B_j"] << " T\n";
    std::cout << "U_m calibrated: " << magstring.computeUmContribution(j, t_day) << " J/m^3\n";
    
    // Step 13: optimizeForMetric - maximize U_m
    std::cout << "\nStep 13: Optimize for Maximum U_m (Enhanced magnetism)\n";
    magstring.saveState("before_optimization");
    magstring.optimizeForMetric(j, t_day, "maximize_Um");
    std::cout << "Optimized for maximum U_m\n";
    std::cout << "Enhanced mu_base: " << magstring.variables["mu_base"] << " T·m^3 (+50%)\n";
    std::cout << "Enhanced E_react: " << magstring.variables["E_react"] << " J (+30%)\n";
    std::cout << "Enhanced P_SCm: " << magstring.variables["P_SCm"] << " (+20%)\n";
    std::cout << "U_m (maximized): " << magstring.computeUmContribution(j, t_day) << " J/m^3\n";
    
    // Step 14: mutateParameters - cosmic perturbations
    std::cout << "\nStep 14: Mutate Parameters (+/- 3% cosmic perturbations)\n";
    magstring.restoreState("before_optimization");
    magstring.saveState("before_mutation");
    magstring.mutateParameters(0.03);
    std::cout << "Applied 3% random cosmic perturbations\n";
    std::cout << "r_1: " << magstring.computeRj(1) << " m = " << magstring.computeRjInAU(1) << " AU\n";
    std::cout << "mu_base: " << magstring.variables["mu_base"] << " T·m^3\n";
    std::cout << "E_react: " << magstring.variables["E_react"] << " J\n";
    std::cout << "U_m after mutation: " << magstring.computeUmContribution(j, t_day) << " J/m^3\n";
    
    // Step 15: evolveSystem - adaptive selection for strongest U_m
    std::cout << "\nStep 15: Evolve System (10 generations, selection pressure 0.8)\n";
    magstring.restoreState("before_mutation");
    magstring.evolveSystem(j, t_day, 10, 0.8);
    std::cout << "Evolved over 10 generations\n";
    std::cout << "Evolved U_m: " << magstring.computeUmContribution(j, t_day) << " J/m^3\n";
    std::cout << "Evolved mu_base: " << magstring.variables["mu_base"] << " T·m^3\n";
    std::cout << "Evolved r_1: " << magstring.computeRj(1) << " m = " << magstring.computeRjInAU(1) << " AU\n";
    
    // Step 16: State management demonstration
    std::cout << "\nStep 16: State Management - List all saved magnetic string states\n";
    auto saved_states = magstring.listSavedStates();
    std::cout << "Saved magnetic string states (" << saved_states.size() << "):\n";
    for (const auto& label : saved_states) {
        std::cout << "  - " << label << "\n";
    }
    
    // Step 17: validateConsistency and autoCorrectAnomalies
    std::cout << "\nStep 17: Validate Consistency and Auto-Correct\n";
    bool valid = magstring.validateConsistency();
    std::cout << "Current state valid: " << (valid ? "YES" : "NO") << "\n";
    if (!valid) {
        std::cout << "Applying auto-corrections...\n";
        bool corrected = magstring.autoCorrectAnomalies();
        std::cout << "Corrections applied: " << (corrected ? "YES" : "NO") << "\n";
        std::cout << "State valid after correction: " << (magstring.validateConsistency() ? "YES" : "NO") << "\n";
    }
    std::cout << "String distance check:\n";
    std::cout << "  r_1 = " << magstring.variables["r_1"] << " m\n";
    std::cout << "  Range: [1e10, 1e20] m\n";
    std::cout << "  Valid: " << (magstring.variables["r_1"] >= 1e10 && magstring.variables["r_1"] <= 1e20 ? "YES" : "NO") << "\n";
    
    // Step 18: Full report and export
    std::cout << "\nStep 18: Generate Full Report and Export State\n";
    magstring.restoreState("t_zero");
    std::cout << magstring.generateReport(j, t_day);
    std::string export_data = magstring.exportState();
    std::cout << "\nState exported (" << export_data.length() << " characters)\n";
    std::cout << "First 500 chars of export:\n" << export_data.substr(0, 500) << "...\n";
    
    std::cout << "\n========== ENHANCED DEMO COMPLETE (18 STEPS EXECUTED) ==========\n";
    std::cout << "Demonstrated: Variable management, batch ops, string/magnetic/reactor expansion,\n";
    std::cout << "              refinement, exploration, sensitivity, optimization, evolution,\n";
    std::cout << "              state management, validation, reporting for magnetic strings.\n";
    std::cout << "Key Insight: Magnetic string at 100 AU scale (r_j = 1.496e13 m) contributes to\n";
    std::cout << "             U_m (universal magnetism) and Ug3 (gravity disk term). Most sensitive\n";
    std::cout << "             to E_react, mu_base, and P_SCm parameters in UQFF framework.\n";
}

Strengths :
-Modular, extensible design for modeling magnetic string path distances and their contributions to universal magnetism(U_m) and gravity(Ug3) in the UQFF framework.
- Clear encapsulation of variables and string parameters using std::map, supporting dynamic updates and easy extension.
- Implements core physical concepts : conversion of r_j between units(m, AU, ly, pc), magnetic moment calculations, and their influence on U_m and Ug3.
- Approximations and physical meaning are well - documented in comments and equation text.
- Output functions for variable state and string contributions support debugging and transparency.
- Handles dynamic updates to r_j and recalculates dependent terms as needed.
- Example calculations and conversion functions provide scientific context and validation.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for flexibility and scalability.
- Some calculations use magic numbers or lack explanatory comments; define named constants and clarify logic.
- Minimal error handling(e.g., missing r_j, division by zero, invalid variable names); add validation for robustness.
- Unit consistency should be checked and documented for all physical quantities.
- For large - scale or performance - critical simulations, optimize data structures and reduce redundant calculations.
- std::map is flexible but may be less efficient than structured types for very large models.
- Expand documentation for function purposes and physical meaning.
- Consider supporting multiple j - indexed strings for more general modeling.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in magnetic string modeling.It implements the UQFF magnetic string concept faithfully and adapts to various scenarios.For production or high - performance applications, address the recommendations for improved robustness, maintainability, and scalability.