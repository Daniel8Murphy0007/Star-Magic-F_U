// UgIndexModule.h
// Modular C++ implementation of the Index for Discrete Universal Gravity Ranges (i) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module uses i=1 to 4 to label Ug1-Ug4; computes sum_{i=1}^4 k_i * U_gi for F_U contribution.
// Pluggable: #include "UgIndexModule.h"
// UgIndexModule mod; mod.computeSumKUgi(); mod.updateVariable("U_g1", new_value);
// Variables in std::map; defaults for Sun at t=0; i labels: 1=Internal Dipole, 2=Outer Bubble, 3=Magnetic Disk, 4=Star-BH.
// Approximations: k_i from coupling; sum ?1.42e53 J/m� (Ug2 dominant).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef UG_INDEX_MODULE_H
#define UG_INDEX_MODULE_H

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

class UgIndexModule {
private:
    std::map<std::string, double> variables;
    std::vector<double> k_values;  // [k1=1.5, k2=1.2, k3=1.8, k4=1.0]
    std::vector<double> computeAllKUgi();

public:
    // Constructor: Initialize with framework defaults
    UgIndexModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    int getIndexRange();  // i=1 to 4
    double computeU_gi(int i);  // U_gi for i=1-4 (J/m^3)
    double computeK_i(int i);   // k_i for i
    double computeKUgi(int i);  // k_i * U_gi
    double computeSumKUgi(int i_min=1, int i_max=4);  // Sum for F_U

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print breakdown by i
    void printIndexBreakdown();

    // ===== Enhanced: 25-Method Dynamic Self-Update & Self-Expansion Capabilities =====

    // Variable Management (5 methods)
    void createVariable(const std::string& name, double value);
    void removeVariable(const std::string& name);
    void cloneVariable(const std::string& source, const std::string& target);
    std::string listVariables();
    std::string getSystemName();

    // Batch Operations (2 methods)
    void transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func);
    void scaleVariableGroup(const std::vector<std::string>& names, double factor);

    // Self-Expansion (4 methods: 1 generic + 3 domain-specific)
    void expandParameterSpace(const std::vector<std::string>& params, double expansion_factor);
    void expandGravityScale(double U_g_factor);                      // Scale all U_g1-U_g4
    void expandCouplingScale(double k_factor);                        // Scale all k_i
    void expandRangeScale(int range_i, double U_factor, double k_factor); // Scale specific i

    // Self-Refinement (3 methods)
    void autoRefineParameters(const std::string& target_metric, double target_value);
    void calibrateToObservations(const std::map<std::string, double>& observed);
    void optimizeForMetric(const std::string& metric);

    // Parameter Exploration (1 method)
    std::vector<std::map<std::string, double>> generateVariations(int count, double variance);

    // Adaptive Evolution (2 methods)
    void mutateParameters(double mutation_rate);
    void evolveSystem(int generations, std::function<double()> fitness_func);

    // State Management (4 methods)
    void saveState(const std::string& label);
    void restoreState(const std::string& label);
    std::vector<std::string> listSavedStates();
    std::string exportState();

    // System Analysis (4 methods)
    std::map<std::string, double> sensitivityAnalysis(const std::vector<std::string>& params);
    std::string generateReport();
    bool validateConsistency();
    void autoCorrectAnomalies();
};

#endif // UG_INDEX_MODULE_H

// UgIndexModule.cpp
#include "UgIndexModule.h"

// Constructor: Set defaults for Sun at t=0
UgIndexModule::UgIndexModule() {
    // Coupling constants (unitless)
    k_values = {1.5, 1.2, 1.8, 1.0};               // k1 to k4

    // U_gi defaults (J/m^3, Sun t=0)
    variables["U_g1"] = 1.39e26;                    // Internal Dipole
    variables["U_g2"] = 1.18e53;                    // Outer Field Bubble
    variables["U_g3"] = 1.8e49;                     // Magnetic Strings Disk
    variables["U_g4"] = 2.50e-20;                   // Star-Black Hole Interactions

    // Shared params (placeholders)
    variables["t_n"] = 0.0;                         // s
    variables["pi"] = 3.141592653589793;
}

// Update variable
void UgIndexModule::updateVariable(const std::string& name, double value) {
    variables[name] = value;
}

// Add delta
void UgIndexModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void UgIndexModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Get range of i
int UgIndexModule::getIndexRange() {
    return 4;  // i=1 to 4
}

// Compute U_gi
double UgIndexModule::computeU_gi(int i) {
    std::string key = "U_g" + std::to_string(i);
    if (variables.find(key) != variables.end()) {
        return variables[key];
    }
    std::cerr << "U_g" << i << " not found. Returning 0." << std::endl;
    return 0.0;
}

// Compute k_i (1-based)
double UgIndexModule::computeK_i(int i) {
    if (i < 1 || i > 4) {
        std::cerr << "Invalid i: " << i << ". Using k1." << std::endl;
        return k_values[0];
    }
    return k_values[i-1];
}

// Compute k_i * U_gi
double UgIndexModule::computeKUgi(int i) {
    return computeK_i(i) * computeU_gi(i);
}

// Sum over i_min to i_max
double UgIndexModule::computeSumKUgi(int i_min, int i_max) {
    double sum = 0.0;
    for (int i = i_min; i <= i_max; ++i) {
        sum += computeKUgi(i);
    }
    return sum;
}

// Equation text
std::string UgIndexModule::getEquationText() {
    return "F_U = ?_{i=1}^4 [k_i * U_gi(r,t,M_s,?_s,T_s,B_s,?_vac,[SCm],?_vac,[UA],t_n) - ?_i * ... ] + other terms\n"
           "i (dimensionless integer): Labels Ug ranges; i=1: Internal Dipole, i=2: Outer Bubble,\n"
           "i=3: Magnetic Disk, i=4: Star-BH.\n"
           "Discretizes gravity for summation; enables scale-specific modeling.\n"
           "Example Sun t=0: ? k_i U_gi ?1.42e53 J/m� (Ug2 dominant).\n"
           "Role: Structures Ug contributions; extensible for more ranges.";
}

// Print variables
void UgIndexModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
    std::cout << "k_i: k1=1.5, k2=1.2, k3=1.8, k4=1.0\n";
}

// Print breakdown
void UgIndexModule::printIndexBreakdown() {
    std::cout << "Ug Index Breakdown (i=1 to 4):\n";
    for (int i = 1; i <= 4; ++i) {
        double ugi = computeU_gi(i);
        double ki = computeK_i(i);
        double kugi = computeKUgi(i);
        std::string label;
        switch(i) {
            case 1: label = "Internal Dipole"; break;
            case 2: label = "Outer Field Bubble"; break;
            case 3: label = "Magnetic Strings Disk"; break;
            case 4: label = "Star-Black Hole"; break;
            default: label = "Unknown";
        }
        std::cout << "i=" << i << " (" << label << "): U_g" << i << "=" << std::scientific << ugi
                  << ", k" << i << "=" << ki << ", k_i U_gi=" << kugi << " J/m�\n";
    }
    std::cout << "Sum ? k_i U_gi = " << std::scientific << computeSumKUgi() << " J/m�\n";
}

// ===== Implementation: 25-Method Dynamic Self-Update & Self-Expansion Capabilities =====

namespace ug_index_saved_states {
    std::map<std::string, std::map<std::string, double>> saved_states;
    std::map<std::string, std::vector<double>> saved_k_values;
}

// Variable Management (5 methods)
void UgIndexModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void UgIndexModule::removeVariable(const std::string& name) {
    variables.erase(name);
}

void UgIndexModule::cloneVariable(const std::string& source, const std::string& target) {
    if (variables.find(source) != variables.end()) {
        variables[target] = variables[source];
    }
}

std::string UgIndexModule::listVariables() {
    std::ostringstream oss;
    for (const auto& pair : variables) {
        oss << pair.first << " = " << pair.second << "\n";
    }
    oss << "k_values: [";
    for (size_t i = 0; i < k_values.size(); ++i) {
        oss << k_values[i];
        if (i < k_values.size() - 1) oss << ", ";
    }
    oss << "]\n";
    return oss.str();
}

std::string UgIndexModule::getSystemName() {
    return "Ug_Index_UQFF";
}

// Batch Operations (2 methods)
void UgIndexModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        if (variables.find(name) != variables.end()) {
            variables[name] = func(variables[name]);
        }
    }
}

void UgIndexModule::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    transformVariableGroup(names, [factor](double x) { return x * factor; });
}

// Self-Expansion (4 methods)
void UgIndexModule::expandParameterSpace(const std::vector<std::string>& params, double expansion_factor) {
    for (const auto& param : params) {
        if (variables.find(param) != variables.end()) {
            variables[param] *= expansion_factor;
        }
    }
}

void UgIndexModule::expandGravityScale(double U_g_factor) {
    // Scale all U_g1 through U_g4
    for (int i = 1; i <= 4; ++i) {
        std::string key = "U_g" + std::to_string(i);
        if (variables.find(key) != variables.end()) {
            variables[key] *= U_g_factor;
        }
    }
}

void UgIndexModule::expandCouplingScale(double k_factor) {
    // Scale all k_i coupling constants
    for (size_t i = 0; i < k_values.size(); ++i) {
        k_values[i] *= k_factor;
    }
}

void UgIndexModule::expandRangeScale(int range_i, double U_factor, double k_factor) {
    // Scale specific U_gi and k_i
    if (range_i >= 1 && range_i <= 4) {
        std::string key = "U_g" + std::to_string(range_i);
        if (variables.find(key) != variables.end()) {
            variables[key] *= U_factor;
        }
        k_values[range_i - 1] *= k_factor;
    }
}

// Self-Refinement (3 methods)
void UgIndexModule::autoRefineParameters(const std::string& target_metric, double target_value) {
    // Example: target sum of k_i * U_gi by adjusting all U_gi proportionally
    if (target_metric == "sum_kUgi") {
        double current = computeSumKUgi();
        if (current > 0) {
            double ratio = target_value / current;
            for (int i = 1; i <= 4; ++i) {
                std::string key = "U_g" + std::to_string(i);
                if (variables.find(key) != variables.end()) {
                    variables[key] *= ratio;
                }
            }
        }
    } else if (target_metric.substr(0, 4) == "kUgi") {
        // Target specific k_i * U_gi
        int i = std::stoi(target_metric.substr(4, 1));
        double current = computeKUgi(i);
        if (current > 0) {
            double ratio = target_value / current;
            std::string key = "U_g" + std::to_string(i);
            variables[key] *= ratio;
        }
    }
}

void UgIndexModule::calibrateToObservations(const std::map<std::string, double>& observed) {
    for (const auto& obs : observed) {
        if (variables.find(obs.first) != variables.end()) {
            variables[obs.first] = obs.second;
        }
        // Handle k_i calibration
        if (obs.first.substr(0, 2) == "k_" && obs.first.length() == 3) {
            int idx = obs.first[2] - '1';
            if (idx >= 0 && idx < 4) {
                k_values[idx] = obs.second;
            }
        }
    }
}

void UgIndexModule::optimizeForMetric(const std::string& metric) {
    if (metric == "maximize_sum") {
        expandGravityScale(1.1);
        expandCouplingScale(1.05);
    } else if (metric == "minimize_sum") {
        expandGravityScale(0.9);
        expandCouplingScale(0.95);
    } else if (metric == "balance_ranges") {
        // Equalize contribution from each i
        double avg = computeSumKUgi() / 4.0;
        for (int i = 1; i <= 4; ++i) {
            double current = computeKUgi(i);
            if (current > 0) {
                double ratio = avg / current;
                std::string key = "U_g" + std::to_string(i);
                variables[key] *= ratio;
            }
        }
    } else if (metric == "enhance_Ug2") {
        // Boost Ug2 (dominant outer bubble)
        variables["U_g2"] *= 1.2;
    }
}

// Parameter Exploration (1 method)
std::vector<std::map<std::string, double>> UgIndexModule::generateVariations(int count, double variance) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-variance, variance);

    for (int i = 0; i < count; ++i) {
        std::map<std::string, double> variant = variables;
        for (auto& pair : variant) {
            if (pair.first != "pi" && pair.first != "t_n") {
                double factor = 1.0 + dis(gen);
                pair.second *= factor;
            }
        }
        variations.push_back(variant);
    }
    return variations;
}

// Adaptive Evolution (2 methods)
void UgIndexModule::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-mutation_rate, mutation_rate);

    // Mutate U_gi values
    for (auto& pair : variables) {
        if (pair.first != "pi" && pair.first != "t_n") {
            pair.second *= (1.0 + dis(gen));
            if (pair.second < 0) pair.second = std::abs(pair.second);
        }
    }

    // Mutate k_i values
    for (size_t i = 0; i < k_values.size(); ++i) {
        k_values[i] *= (1.0 + dis(gen));
        if (k_values[i] < 0) k_values[i] = std::abs(k_values[i]);
    }
}

void UgIndexModule::evolveSystem(int generations, std::function<double()> fitness_func) {
    double best_fitness = fitness_func();
    std::map<std::string, double> best_state = variables;
    std::vector<double> best_k = k_values;

    for (int gen = 0; gen < generations; ++gen) {
        mutateParameters(0.05);
        double current_fitness = fitness_func();
        if (current_fitness > best_fitness) {
            best_fitness = current_fitness;
            best_state = variables;
            best_k = k_values;
        } else {
            variables = best_state;
            k_values = best_k;
        }
    }
    variables = best_state;
    k_values = best_k;
}

// State Management (4 methods)
void UgIndexModule::saveState(const std::string& label) {
    ug_index_saved_states::saved_states[label] = variables;
    ug_index_saved_states::saved_k_values[label] = k_values;
}

void UgIndexModule::restoreState(const std::string& label) {
    if (ug_index_saved_states::saved_states.find(label) != ug_index_saved_states::saved_states.end()) {
        variables = ug_index_saved_states::saved_states[label];
        k_values = ug_index_saved_states::saved_k_values[label];
    }
}

std::vector<std::string> UgIndexModule::listSavedStates() {
    std::vector<std::string> labels;
    for (const auto& pair : ug_index_saved_states::saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string UgIndexModule::exportState() {
    std::ostringstream oss;
    oss << "System: " << getSystemName() << "\n";
    oss << listVariables();
    return oss.str();
}

// System Analysis (4 methods)
std::map<std::string, double> UgIndexModule::sensitivityAnalysis(const std::vector<std::string>& params) {
    std::map<std::string, double> sensitivity;
    double baseline = computeSumKUgi();

    for (const auto& param : params) {
        if (variables.find(param) != variables.end()) {
            double original = variables[param];
            variables[param] *= 1.01; // +1% perturbation
            double perturbed = computeSumKUgi();
            variables[param] = original;
            sensitivity[param] = std::abs((perturbed - baseline) / baseline);
        }
    }

    // Sensitivity for k_i
    for (int i = 1; i <= 4; ++i) {
        std::string k_name = "k_" + std::to_string(i);
        double original_k = k_values[i - 1];
        k_values[i - 1] *= 1.01;
        double perturbed = computeSumKUgi();
        k_values[i - 1] = original_k;
        sensitivity[k_name] = std::abs((perturbed - baseline) / baseline);
    }

    return sensitivity;
}

std::string UgIndexModule::generateReport() {
    std::ostringstream oss;
    oss << "========== Ug Index Module Report ==========\n";
    oss << "System: " << getSystemName() << "\n\n";
    oss << "Index Structure (i=1 to 4):\n";
    oss << "  i=1: Internal Dipole (stellar core dynamics)\n";
    oss << "  i=2: Outer Field Bubble (heliosphere ~100 AU)\n";
    oss << "  i=3: Magnetic Strings Disk (Kuiper belt scale)\n";
    oss << "  i=4: Star-Black Hole Interactions (galactic)\n\n";
    
    oss << "Coupling Constants:\n";
    for (int i = 1; i <= 4; ++i) {
        oss << "  k_" << i << " = " << k_values[i-1] << " (unitless)\n";
    }
    oss << "\n";

    oss << "Universal Gravity Components:\n";
    double total = 0.0;
    for (int i = 1; i <= 4; ++i) {
        double ugi = computeU_gi(i);
        double ki = computeK_i(i);
        double kugi = computeKUgi(i);
        total += kugi;
        oss << "  U_g" << i << " = " << std::scientific << ugi << " J/m^3\n";
        oss << "    k_" << i << " * U_g" << i << " = " << kugi << " J/m^3";
        if (total > 0) {
            oss << " (" << (kugi / total * 100) << "% of total)\n";
        } else {
            oss << "\n";
        }
    }
    oss << "\n";

    oss << "Total Gravity Contribution:\n";
    oss << "  Sum(k_i * U_gi) = " << std::scientific << total << " J/m^3\n";
    oss << "  Dominant range: ";
    int dominant = 2; // Default Ug2
    double max_val = computeKUgi(2);
    for (int i = 1; i <= 4; ++i) {
        if (computeKUgi(i) > max_val) {
            max_val = computeKUgi(i);
            dominant = i;
        }
    }
    oss << "i=" << dominant << " (contributes " << (max_val / total * 100) << "%)\n\n";

    oss << "Physical Context:\n";
    oss << "  Discrete indexing enables scale-specific modeling\n";
    oss << "  Each i represents distinct spatial/energy scale\n";
    oss << "  Extensible to additional ranges (i=5, 6, ...)\n";
    oss << "  Sum appears in F_U unified field equation\n";
    oss << "============================================\n";
    return oss.str();
}

bool UgIndexModule::validateConsistency() {
    bool valid = true;
    // Check U_gi are positive
    for (int i = 1; i <= 4; ++i) {
        std::string key = "U_g" + std::to_string(i);
        if (variables.find(key) != variables.end() && variables[key] < 0) {
            valid = false;
        }
    }
    // Check k_i are positive
    for (size_t i = 0; i < k_values.size(); ++i) {
        if (k_values[i] <= 0) valid = false;
    }
    return valid;
}

void UgIndexModule::autoCorrectAnomalies() {
    // Reset U_gi to defaults if negative or invalid
    if (variables["U_g1"] <= 0) variables["U_g1"] = 1.39e26;
    if (variables["U_g2"] <= 0) variables["U_g2"] = 1.18e53;
    if (variables["U_g3"] <= 0) variables["U_g3"] = 1.8e49;
    if (variables["U_g4"] <= 0) variables["U_g4"] = 2.50e-20;
    
    // Reset k_i to defaults if invalid
    for (size_t i = 0; i < k_values.size(); ++i) {
        if (k_values[i] <= 0) {
            k_values = {1.5, 1.2, 1.8, 1.0};
            break;
        }
    }
}

// Example usage in base program (snippet)
// #include "UgIndexModule.h"
// int main() {
//     UgIndexModule mod;
//     double sum = mod.computeSumKUgi();
//     std::cout << "? k_i U_gi = " << sum << " J/m�\n";
//     mod.printIndexBreakdown();
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("U_g3", 2e49);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o ug_index_test ug_index_test.cpp UgIndexModule.cpp -lm
// Sample: Sum ?1.42e53 J/m�; i structures 4 Ug ranges.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

// ===== Enhanced Example: 18-Step Demonstration of Dynamic Capabilities =====
// int main() {
//     UgIndexModule mod;
//     std::cout << "===== Ug Index Module: Enhanced 18-Step Demo =====\n\n";
//
//     // Step 1: Initial report and breakdown
//     std::cout << "Step 1: Initial Configuration\n";
//     std::cout << mod.generateReport() << "\n";
//     mod.printIndexBreakdown();
//     std::cout << "\n";
//
//     // Step 2: Variable tracking for each range
//     std::cout << "Step 2: Create Range-Specific Tracking Variables\n";
//     for (int i = 1; i <= 4; ++i) {
//         std::string contrib_name = "contribution_i" + std::to_string(i);
//         mod.createVariable(contrib_name, mod.computeKUgi(i));
//         std::cout << "  " << contrib_name << " = " << std::scientific << mod.variables[contrib_name] << " J/m^3\n";
//     }
//     double total_contrib = mod.computeSumKUgi();
//     mod.createVariable("total_gravity_contribution", total_contrib);
//     std::cout << "  total_gravity_contribution = " << total_contrib << " J/m^3\n\n";
//
//     // Step 3: Individual range variations
//     std::cout << "Step 3: Explore Individual Range Scaling\n";
//     mod.saveState("baseline");
//     std::vector<double> scale_factors = {0.8, 1.0, 1.2, 1.5, 2.0};
//     for (int i = 1; i <= 4; ++i) {
//         std::cout << "  Range i=" << i << " scaling:\n";
//         for (double factor : scale_factors) {
//             mod.restoreState("baseline");
//             std::string key = "U_g" + std::to_string(i);
//             double original = mod.variables[key];
//             mod.updateVariable(key, original * factor);
//             double new_sum = mod.computeSumKUgi();
//             std::cout << "    U_g" << i << " x" << factor << " -> Sum = " << new_sum 
//                       << " (" << std::showpos << ((new_sum/total_contrib - 1.0)*100) << std::noshowpos << "%)\n";
//         }
//     }
//     mod.restoreState("baseline");
//     std::cout << "\n";
//
//     // Step 4: Coupling constant variations
//     std::cout << "Step 4: Explore Coupling Constant Scaling\n";
//     for (int i = 1; i <= 4; ++i) {
//         std::cout << "  k_" << i << " variations:\n";
//         for (double k_mult : {0.5, 1.0, 1.5, 2.0}) {
//             mod.restoreState("baseline");
//             mod.expandRangeScale(i, 1.0, k_mult);
//             double new_sum = mod.computeSumKUgi();
//             std::cout << "    k_" << i << " x" << k_mult << " -> Sum = " << new_sum 
//                       << " (" << std::showpos << ((new_sum/total_contrib - 1.0)*100) << std::noshowpos << "%)\n";
//         }
//     }
//     mod.restoreState("baseline");
//     std::cout << "\n";
//
//     // Step 5: Expand all gravity scales uniformly
//     std::cout << "Step 5: Expand All Gravity Scales (U_g x1.3)\n";
//     mod.expandGravityScale(1.3);
//     std::cout << "  New Sum = " << mod.computeSumKUgi() << " J/m^3\n";
//     mod.printIndexBreakdown();
//     std::cout << "\n";
//
//     // Step 6: Expand all coupling constants
//     std::cout << "Step 6: Expand All Coupling Constants (k_i x1.1)\n";
//     mod.expandCouplingScale(1.1);
//     std::cout << "  New Sum = " << mod.computeSumKUgi() << " J/m^3\n";
//     mod.printIndexBreakdown();
//     std::cout << "\n";
//
//     // Step 7: Expand specific range (boost Ug2 - outer bubble)
//     std::cout << "Step 7: Expand Specific Range i=2 (U_g2 x1.5, k_2 x1.2)\n";
//     mod.expandRangeScale(2, 1.5, 1.2);
//     std::cout << "  New U_g2 = " << mod.variables["U_g2"] << " J/m^3\n";
//     std::cout << "  New k_2 = " << mod.k_values[1] << "\n";
//     std::cout << "  New Sum = " << mod.computeSumKUgi() << " J/m^3\n\n";
//
//     // Step 8: Parameter variations
//     std::cout << "Step 8: Generate 10 Parameter Variations (±8%)\n";
//     auto variations = mod.generateVariations(10, 0.08);
//     std::cout << "  Generated " << variations.size() << " configurations\n";
//     std::vector<double> var_sums;
//     for (size_t v = 0; v < variations.size(); ++v) {
//         double var_sum = 0.0;
//         for (int i = 1; i <= 4; ++i) {
//             std::string key = "U_g" + std::to_string(i);
//             var_sum += mod.k_values[i-1] * variations[v][key];
//         }
//         var_sums.push_back(var_sum);
//     }
//     auto minmax = std::minmax_element(var_sums.begin(), var_sums.end());
//     std::cout << "  Sum range: " << *minmax.first << " to " << *minmax.second << " J/m^3\n\n";
//
//     // Step 9: Sensitivity analysis
//     std::cout << "Step 9: Sensitivity Analysis (Sum response to ±1% changes)\n";
//     std::vector<std::string> sens_params = {"U_g1", "U_g2", "U_g3", "U_g4"};
//     auto sensitivity = mod.sensitivityAnalysis(sens_params);
//     std::cout << "  Parameter sensitivities:\n";
//     std::vector<std::pair<std::string, double>> sorted_sens(sensitivity.begin(), sensitivity.end());
//     std::sort(sorted_sens.begin(), sorted_sens.end(), 
//               [](const auto& a, const auto& b) { return a.second > b.second; });
//     for (const auto& s : sorted_sens) {
//         std::cout << "    " << s.first << ": " << (s.second * 100) << "% sensitivity\n";
//     }
//     std::cout << "\n";
//
//     // Step 10: Auto-refinement to target sum
//     std::cout << "Step 10: Auto-Refine to Target Sum = 2.0e53 J/m^3\n";
//     double target_sum = 2.0e53;
//     double before_refine = mod.computeSumKUgi();
//     mod.autoRefineParameters("sum_kUgi", target_sum);
//     double after_refine = mod.computeSumKUgi();
//     std::cout << "  Before: " << before_refine << " J/m^3\n";
//     std::cout << "  After: " << after_refine << " J/m^3\n";
//     std::cout << "  Target: " << target_sum << " J/m^3\n";
//     std::cout << "  Error: " << std::abs(after_refine - target_sum) / target_sum * 100 << "%\n\n";
//
//     // Step 11: Target specific range contribution
//     std::cout << "Step 11: Target Specific Range (k_2 * U_g2 = 1.5e53 J/m^3)\n";
//     double target_kg2 = 1.5e53;
//     double before_kg2 = mod.computeKUgi(2);
//     mod.autoRefineParameters("kUgi2", target_kg2);
//     double after_kg2 = mod.computeKUgi(2);
//     std::cout << "  Before: " << before_kg2 << " J/m^3\n";
//     std::cout << "  After: " << after_kg2 << " J/m^3\n";
//     std::cout << "  Target: " << target_kg2 << " J/m^3\n\n";
//
//     // Step 12: Calibration to observations
//     std::cout << "Step 12: Calibrate to Observational Data\n";
//     std::map<std::string, double> observations = {
//         {"U_g1", 1.5e26},
//         {"U_g2", 1.3e53},
//         {"U_g3", 2.0e49},
//         {"k_1", 1.6},
//         {"k_2", 1.3}
//     };
//     mod.calibrateToObservations(observations);
//     std::cout << "  Calibrated U_g1: " << mod.variables["U_g1"] << " J/m^3\n";
//     std::cout << "  Calibrated U_g2: " << mod.variables["U_g2"] << " J/m^3\n";
//     std::cout << "  Calibrated U_g3: " << mod.variables["U_g3"] << " J/m^3\n";
//     std::cout << "  Calibrated k_1: " << mod.k_values[0] << "\n";
//     std::cout << "  Calibrated k_2: " << mod.k_values[1] << "\n";
//     std::cout << "  New Sum: " << mod.computeSumKUgi() << " J/m^3\n\n";
//
//     // Step 13: Optimization for balanced ranges
//     std::cout << "Step 13: Optimize for Balanced Range Contributions\n";
//     mod.saveState("before_balance");
//     mod.optimizeForMetric("balance_ranges");
//     std::cout << "  After balancing:\n";
//     for (int i = 1; i <= 4; ++i) {
//         double contrib = mod.computeKUgi(i);
//         std::cout << "    k_" << i << " * U_g" << i << " = " << contrib << " J/m^3\n";
//     }
//     std::cout << "  New Sum: " << mod.computeSumKUgi() << " J/m^3\n\n";
//
//     // Step 14: Restore and enhance Ug2
//     std::cout << "Step 14: Restore and Enhance Ug2 (Outer Bubble Dominance)\n";
//     mod.restoreState("before_balance");
//     mod.optimizeForMetric("enhance_Ug2");
//     std::cout << "  Enhanced U_g2: " << mod.variables["U_g2"] << " J/m^3\n";
//     std::cout << "  k_2 * U_g2: " << mod.computeKUgi(2) << " J/m^3\n";
//     std::cout << "  New Sum: " << mod.computeSumKUgi() << " J/m^3\n\n";
//
//     // Step 15: Parameter mutation
//     std::cout << "Step 15: Mutate Parameters (±4%)\n";
//     mod.saveState("before_mutation");
//     double ug2_before = mod.variables["U_g2"];
//     double k2_before = mod.k_values[1];
//     mod.mutateParameters(0.04);
//     std::cout << "  U_g2: " << ug2_before << " -> " << mod.variables["U_g2"] << " J/m^3\n";
//     std::cout << "  k_2: " << k2_before << " -> " << mod.k_values[1] << "\n";
//     std::cout << "  New Sum: " << mod.computeSumKUgi() << " J/m^3\n\n";
//
//     // Step 16: System evolution (maximize sum)
//     std::cout << "Step 16: Evolve System (8 generations, maximize Sum)\n";
//     mod.restoreState("before_mutation");
//     double initial_fitness = mod.computeSumKUgi();
//     mod.evolveSystem(8, [&mod]() { return mod.computeSumKUgi(); });
//     double final_fitness = mod.computeSumKUgi();
//     std::cout << "  Initial Sum: " << initial_fitness << " J/m^3\n";
//     std::cout << "  Evolved Sum: " << final_fitness << " J/m^3\n";
//     std::cout << "  Improvement: " << ((final_fitness / initial_fitness - 1.0) * 100) << "%\n\n";
//
//     // Step 17: Validation and auto-correction
//     std::cout << "Step 17: Validate Consistency\n";
//     bool valid = mod.validateConsistency();
//     std::cout << "  System valid: " << (valid ? "YES" : "NO") << "\n";
//     if (!valid) {
//         std::cout << "  Running auto-correction...\n";
//         mod.autoCorrectAnomalies();
//         std::cout << "  Post-correction valid: " << (mod.validateConsistency() ? "YES" : "NO") << "\n";
//     }
//     std::cout << "\n";
//
//     // Step 18: State management and export
//     std::cout << "Step 18: State Management and Export\n";
//     mod.saveState("evolved_optimal");
//     auto saved = mod.listSavedStates();
//     std::cout << "  Saved states (" << saved.size() << "): ";
//     for (const auto& s : saved) std::cout << s << " ";
//     std::cout << "\n\n";
//     std::cout << "Final System Export:\n";
//     std::string exported = mod.exportState();
//     std::cout << exported << "\n";
//     std::cout << "Demo complete: 18 steps executed successfully!\n";
//     std::cout << "================================================\n";
//
//     return 0;
// }

UgIndexModule Evaluation

Strengths :
-Modular and pluggable design; can be included and instantiated easily in other projects.
- Dynamic variable management using std::map allows runtime updates, additions, and removals.
- Coupling constants and gravity terms are clearly separated and indexed, supporting extensibility.
- Core computation methods(computeU_gi, computeK_i, computeKUgi, computeSumKUgi) are clear, concise, and variable - driven.
- Output and debugging functions(printVariables, printIndexBreakdown, getEquationText) provide transparency and aid validation.
- Well - documented physical meaning and example calculations in comments and equation text.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for greater flexibility.
- Minimal error handling for missing variables, invalid indices, or division by zero; add validation for robustness.
- Unit consistency is described in comments but not enforced; runtime checks or clearer documentation would help.
- For large - scale or performance - critical simulations, consider more efficient data structures than std::map and std::vector.
- Expand documentation for function purposes and expected input / output.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in universal gravity modeling.It is dynamic and can be updated or expanded easily.For production or high - performance applications, address the recommendations above for improved robustness, maintainability, and scalability.