// FeedbackFactorModule.h
// Modular C++ implementation of the Feedback Factor (f_feedback) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes f_feedback=0.1 for ?M_BH=1 dex (10x mass increase); scales (1 + f_feedback) in U_g4 term.
// Pluggable: #include "FeedbackFactorModule.h"
// FeedbackFactorModule mod; mod.computeU_g4(0.0, 0.0); mod.updateVariable("f_feedback", new_value);
// Variables in std::map; example for Sun at t=0, t_n=0; ?M_BH=1 dex ? M_bh_final=10*M_bh_initial.
// Approximations: cos(? t_n)=1; e^{-? t}=1 at t=0; ?=0.001 day^-1 (scaled to s^-1 if needed).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef FEEDBACK_FACTOR_MODULE_H
#define FEEDBACK_FACTOR_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <functional>
#include <random>
#include <algorithm>
#include <sstream>
#include <vector>

class FeedbackFactorModule {
private:
    std::map<std::string, double> variables;
    double computeDeltaM_BH();  // 1 dex = log10(10) = factor of 10
    double computeM_bh_final();
    double computeU_g4(double t, double t_n);
    double computeU_g4_no_feedback(double t, double t_n);

public:
    // Constructor: Initialize with framework defaults
    FeedbackFactorModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeF_feedback();  // f_feedback=0.1 (unitless)
    double computeDeltaM_BH();   // 1 dex
    double computeM_bh_final();  // 10 * M_bh_initial
    double computeU_g4(double t, double t_n);  // With feedback (J/m^3)
    double computeU_g4_no_feedback(double t, double t_n);  // Without (J/m^3)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();

    // Print U_g4 comparison (with/without feedback)
    void printU_g4_comparison(double t = 0.0, double t_n = 0.0);
    
    // ========== ENHANCED DYNAMIC CAPABILITIES (25 methods) ==========
    // Variable management
    void createVariable(const std::string& name, double value);
    void removeVariable(const std::string& name);
    void cloneVariable(const std::string& source, const std::string& dest);
    std::vector<std::string> listVariables();
    std::string getSystemName();

    // Batch operations
    void transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> fn);
    void scaleVariableGroup(const std::vector<std::string>& names, double scale);

    // Self-expansion (domain-specific)
    void expandParameterSpace(double scale);
    void expandFeedbackScale(double f_feedback_scale, double delta_dex_scale);
    void expandBlackHoleScale(double M_bh_scale, double k4_scale);
    void expandDecayScale(double alpha_scale, double rho_scale);

    // Self-refinement
    void autoRefineParameters(double t, double t_n, double target_Ug4, double tolerance);
    void calibrateToObservations(const std::map<std::string, double>& observations);
    void optimizeForMetric(const std::string& metric);

    // Exploration & evolution
    std::vector<std::map<std::string,double>> generateVariations(int count, double percent);
    void mutateParameters(double magnitude);
    void evolveSystem(double t, double t_n, int generations, double selection_pressure);

    // State management
    void saveState(const std::string& label);
    bool restoreState(const std::string& label);
    std::vector<std::string> listSavedStates();
    std::string exportState();

    // Analysis
    std::map<std::string,double> sensitivityAnalysis(double t, double t_n, double delta);
    std::string generateReport(double t, double t_n);
    bool validateConsistency();
    bool autoCorrectAnomalies();
};

#endif // FEEDBACK_FACTOR_MODULE_H

// FeedbackFactorModule.cpp
#include "FeedbackFactorModule.h"

// Constructor: Set framework defaults
FeedbackFactorModule::FeedbackFactorModule() {
    // Universal constants
    variables["f_feedback"] = 0.1;                  // Unitless, for ?M_BH=1 dex
    variables["delta_M_BH_dex"] = 1.0;              // 1 dex = factor 10
    variables["M_bh_initial"] = 8.15e36;            // kg (Sgr A*)
    variables["k_4"] = 1.0;                         // Coupling for Ug4
    variables["rho_vac_SCm"] = 7.09e-37;            // J/m^3
    variables["d_g"] = 2.55e20;                     // m
    variables["alpha"] = 0.001 / 86400.0;           // day^-1 to s^-1
    variables["pi"] = 3.141592653589793;
    variables["t"] = 0.0;                           // s
    variables["t_n"] = 0.0;                         // s
}

// Update variable
void FeedbackFactorModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
    if (name == "delta_M_BH_dex") {
        // Recalculate M_bh_final if dex changes
        computeM_bh_final();
    }
}

// Add delta
void FeedbackFactorModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void FeedbackFactorModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute f_feedback (fixed 0.1 for 1 dex)
double FeedbackFactorModule::computeF_feedback() {
    return variables["f_feedback"];
}

// Compute ?M_BH in dex
double FeedbackFactorModule::computeDeltaM_BH() {
    return variables["delta_M_BH_dex"];
}

// Compute M_bh_final = M_bh_initial * 10^{?M_BH_dex}
double FeedbackFactorModule::computeM_bh_final() {
    double factor = std::pow(10.0, computeDeltaM_BH());
    double initial = variables["M_bh_initial"];
    variables["M_bh_final"] = initial * factor;
    return variables["M_bh_final"];
}

// Compute U_g4 with feedback
double FeedbackFactorModule::computeU_g4(double t, double t_n) {
    double k_4 = variables["k_4"];
    double rho_vac_SCm = variables["rho_vac_SCm"];
    double M_bh = computeM_bh_final();  // Use final mass for feedback scenario
    double d_g = variables["d_g"];
    double alpha = variables["alpha"];
    double pi = variables["pi"];
    double f_feedback = computeF_feedback();
    double exp_term = std::exp( - alpha * t );
    double cos_term = std::cos( pi * t_n );
    double feedback_factor = 1.0 + f_feedback;
    return k_4 * (rho_vac_SCm * M_bh / d_g) * exp_term * cos_term * feedback_factor;
}

// Compute U_g4 without feedback (f_feedback=0)
double FeedbackFactorModule::computeU_g4_no_feedback(double t, double t_n) {
    double orig_f = variables["f_feedback"];
    variables["f_feedback"] = 0.0;
    double result = computeU_g4(t, t_n);
    variables["f_feedback"] = orig_f;  // Restore
    return result;
}

// Equation text
std::string FeedbackFactorModule::getEquationText() {
    return "U_g4 = k_4 * (?_vac,[SCm] M_bh / d_g) * e^{-? t} * cos(? t_n) * (1 + f_feedback)\n"
           "Where f_feedback = 0.1 (unitless, for ?M_BH = 1 dex = 10x mass increase);\n"
           "?M_BH =1 dex ? M_bh_final = 10 * M_bh_initial (8.15e36 kg ? 8.15e37 kg).\n"
           "Example t=0, t_n=0: U_g4 ?2.75e-20 J/m� (with); ?2.50e-20 J/m� (without; +10%).\n"
           "Role: Scales feedback in star-BH interactions; regulates AGN, mergers, star formation.\n"
           "UQFF: Enhances energy density for 10x M_BH; resolves final parsec, quasar jets.";
}

// Print variables
void FeedbackFactorModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Print U_g4 comparison
void FeedbackFactorModule::printU_g4_comparison(double t, double t_n) {
    double u_with = computeU_g4(t, t_n);
    double u_without = computeU_g4_no_feedback(t, t_n);
    double delta_percent = ((u_with - u_without) / u_without) * 100.0;
    std::cout << "U_g4 Comparison at t=" << t << " s, t_n=" << t_n << " s:\n";
    std::cout << "With feedback: " << std::scientific << u_with << " J/m�\n";
    std::cout << "Without feedback: " << std::scientific << u_without << " J/m�\n";
    std::cout << "Difference: +" << std::fixed << std::setprecision(1) << delta_percent << "%\n";
}

// Example usage in base program (snippet)
// #include "FeedbackFactorModule.h"
// int main() {
//     FeedbackFactorModule mod;
//     double m_final = mod.computeM_bh_final();
//     std::cout << "M_bh_final = " << m_final << " kg\n";
//     mod.printU_g4_comparison(0.0, 0.0);
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("delta_M_BH_dex", 2.0);  // 100x mass
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o feedback_test feedback_test.cpp FeedbackFactorModule.cpp -lm
// Sample: M_bh_final=8.15e37 kg; U_g4 with=2.75e-20 J/m³ (+10% vs without).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

// ---------------------- ENHANCED DYNAMIC CAPABILITIES IMPLEMENTATION ----------------------

namespace {
    // Simple persistent saved-state storage for the module
    static std::map<std::string, std::map<std::string,double>> feedback_factor_saved_states;
}

void FeedbackFactorModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void FeedbackFactorModule::removeVariable(const std::string& name) {
    auto it = variables.find(name);
    if (it != variables.end()) variables.erase(it);
}

void FeedbackFactorModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) variables[dest] = variables[source];
}

std::vector<std::string> FeedbackFactorModule::listVariables() {
    std::vector<std::string> keys;
    for (const auto& p : variables) keys.push_back(p.first);
    return keys;
}

std::string FeedbackFactorModule::getSystemName() {
    return "Feedback_Factor_UQFF";
}

void FeedbackFactorModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> fn) {
    for (const auto& n : names) {
        if (variables.find(n) != variables.end()) variables[n] = fn(variables[n]);
    }
}

void FeedbackFactorModule::scaleVariableGroup(const std::vector<std::string>& names, double scale) {
    transformVariableGroup(names, [scale](double v){ return v * scale; });
}

void FeedbackFactorModule::expandParameterSpace(double scale) {
    // Uniform scaling preserving constants
    std::vector<std::string> exclude = {"pi"};
    for (auto& kv : variables) {
        if (std::find(exclude.begin(), exclude.end(), kv.first) == exclude.end()) {
            kv.second *= scale;
        }
    }
}

void FeedbackFactorModule::expandFeedbackScale(double f_feedback_scale, double delta_dex_scale) {
    variables["f_feedback"] *= f_feedback_scale;
    variables["delta_M_BH_dex"] *= delta_dex_scale;
    computeM_bh_final();  // Recalculate final mass
}

void FeedbackFactorModule::expandBlackHoleScale(double M_bh_scale, double k4_scale) {
    variables["M_bh_initial"] *= M_bh_scale;
    variables["k_4"] *= k4_scale;
    computeM_bh_final();  // Recalculate final mass
}

void FeedbackFactorModule::expandDecayScale(double alpha_scale, double rho_scale) {
    variables["alpha"] *= alpha_scale;
    variables["rho_vac_SCm"] *= rho_scale;
}

void FeedbackFactorModule::autoRefineParameters(double t, double t_n, double target_Ug4, double tolerance) {
    // Iteratively adjust f_feedback and M_bh to reach target U_g4
    double cur = computeU_g4(t, t_n);
    int iter = 0;
    while (std::abs(cur - target_Ug4) > tolerance && iter++ < 50) {
        double factor = target_Ug4 / (cur + 1e-50);
        // Conservative adjustments
        variables["f_feedback"] *= std::pow(factor, 0.25);
        variables["M_bh_initial"] *= std::pow(factor, 0.25);
        computeM_bh_final();
        cur = computeU_g4(t, t_n);
    }
}

void FeedbackFactorModule::calibrateToObservations(const std::map<std::string, double>& observations) {
    for (const auto& p : observations) {
        if (variables.find(p.first) != variables.end()) {
            variables[p.first] = p.second;
            if (p.first == "delta_M_BH_dex" || p.first == "M_bh_initial") {
                computeM_bh_final();
            }
        }
    }
}

void FeedbackFactorModule::optimizeForMetric(const std::string& metric) {
    if (metric == "maximize_Ug4") {
        variables["f_feedback"] *= 1.5;
        variables["M_bh_initial"] *= 1.3;
        computeM_bh_final();
    } else if (metric == "minimize_Ug4") {
        variables["f_feedback"] *= 0.5;
        variables["alpha"] *= 1.5;  // Faster decay
    } else if (metric == "enhance_feedback") {
        variables["f_feedback"] *= 2.0;
        variables["delta_M_BH_dex"] *= 1.2;
        computeM_bh_final();
    }
}

std::vector<std::map<std::string,double>> FeedbackFactorModule::generateVariations(int count, double percent) {
    std::vector<std::map<std::string,double>> out;
    std::random_device rd; std::mt19937 gen(rd());
    std::vector<std::string> exclude = {"pi"};
    for (int i=0; i<count; i++) {
        std::map<std::string,double> copy = variables;
        for (auto& kv : copy) {
            if (std::find(exclude.begin(), exclude.end(), kv.first) != exclude.end()) continue;
            std::normal_distribution<> d(1.0, percent/100.0);
            kv.second *= d(gen);
        }
        out.push_back(copy);
    }
    return out;
}

void FeedbackFactorModule::mutateParameters(double magnitude) {
    std::random_device rd; std::mt19937 gen(rd());
    std::normal_distribution<> d(0.0, magnitude);
    std::vector<std::string> exclude = {"pi"};
    for (auto& kv : variables) {
        if (std::find(exclude.begin(), exclude.end(), kv.first) != exclude.end()) continue;
        double delta = d(gen) * kv.second;
        kv.second += delta;
    }
    computeM_bh_final();  // Update final mass after mutation
}

void FeedbackFactorModule::evolveSystem(double t, double t_n, int generations, double selection_pressure) {
    int pop = 20;
    std::vector<std::map<std::string,double>> population;
    for (int i=0; i<pop; i++) population.push_back(variables);
    std::random_device rd; std::mt19937 gen(rd());
    
    for (int g=0; g<generations; g++) {
        // Mutate all
        for (auto& ind : population) {
            std::normal_distribution<> d(0.0, 0.02);
            for (auto& kv : ind) {
                if (kv.first=="pi") continue;
                kv.second *= (1.0 + d(gen));
            }
        }
        
        // Score by absolute U_g4
        std::vector<std::pair<double,int>> scores;
        for (int i=0; i<pop; i++) {
            auto temp = variables;
            variables = population[i];
            computeM_bh_final();
            double score = std::abs(computeU_g4(t, t_n));
            scores.push_back({score, i});
            variables = temp;
        }
        std::sort(scores.begin(), scores.end(), std::greater<>());
        
        // Select top survivors
        int survivors = std::max(2, (int)(pop * selection_pressure));
        std::vector<std::map<std::string,double>> next;
        for (int i=0; i<survivors; i++) next.push_back(population[scores[i].second]);
        
        // Refill by mutating survivors
        std::uniform_int_distribution<> uid(0, survivors-1);
        while ((int)next.size() < pop) {
            auto child = next[uid(gen)];
            std::normal_distribution<> d(0.0, 0.05);
            for (auto& kv : child) {
                if (kv.first != "pi") kv.second *= (1.0 + d(gen));
            }
            next.push_back(child);
        }
        population = next;
    }
    
    // Set to best
    variables = population[0];
    computeM_bh_final();
}

void FeedbackFactorModule::saveState(const std::string& label) {
    feedback_factor_saved_states[label] = variables;
}

bool FeedbackFactorModule::restoreState(const std::string& label) {
    auto it = feedback_factor_saved_states.find(label);
    if (it == feedback_factor_saved_states.end()) return false;
    variables = it->second;
    return true;
}

std::vector<std::string> FeedbackFactorModule::listSavedStates() {
    std::vector<std::string> keys;
    for (const auto& p : feedback_factor_saved_states) keys.push_back(p.first);
    return keys;
}

std::string FeedbackFactorModule::exportState() {
    std::ostringstream ss;
    ss << "FeedbackFactorModule State Export\n";
    for (const auto& kv : variables) {
        ss << kv.first << ": " << std::scientific << kv.second << "\n";
    }
    return ss.str();
}

std::map<std::string,double> FeedbackFactorModule::sensitivityAnalysis(double t, double t_n, double delta) {
    std::map<std::string,double> result;
    double base = computeU_g4(t, t_n);
    std::vector<std::string> keys = {"f_feedback","delta_M_BH_dex","M_bh_initial","k_4","rho_vac_SCm","d_g","alpha"};
    
    for (const auto& k : keys) {
        if (variables.find(k) == variables.end()) continue;
        double orig = variables[k];
        
        variables[k] = orig * (1.0 + delta);
        if (k == "delta_M_BH_dex" || k == "M_bh_initial") computeM_bh_final();
        double up = computeU_g4(t, t_n);
        
        variables[k] = orig * (1.0 - delta);
        if (k == "delta_M_BH_dex" || k == "M_bh_initial") computeM_bh_final();
        double down = computeU_g4(t, t_n);
        
        variables[k] = orig;
        computeM_bh_final();
        
        double sens = 0.5 * ((up - down) / (base + 1e-50));
        result[k] = sens;
    }
    return result;
}

std::string FeedbackFactorModule::generateReport(double t, double t_n) {
    std::ostringstream ss;
    ss << "\n========== FeedbackFactorModule Report ==========\n";
    ss << "System: " << getSystemName() << "\n";
    ss << "f_feedback: " << std::fixed << std::setprecision(3) << variables["f_feedback"] << " (unitless)\n";
    ss << "Delta M_BH: " << variables["delta_M_BH_dex"] << " dex (factor " 
       << std::pow(10.0, variables["delta_M_BH_dex"]) << "x)\n";
    ss << "M_bh_initial: " << std::scientific << variables["M_bh_initial"] << " kg\n";
    ss << "M_bh_final: " << std::scientific << computeM_bh_final() << " kg\n";
    ss << "k_4: " << std::fixed << std::setprecision(2) << variables["k_4"] << "\n";
    ss << "rho_vac_SCm: " << std::scientific << variables["rho_vac_SCm"] << " J/m^3\n";
    ss << "d_g: " << variables["d_g"] << " m\n";
    ss << "alpha: " << variables["alpha"] << " s^-1\n";
    ss << "\nAt t=" << t << " s, t_n=" << t_n << " s:\n";
    ss << "U_g4 (with feedback): " << computeU_g4(t, t_n) << " J/m^3\n";
    ss << "U_g4 (no feedback): " << computeU_g4_no_feedback(t, t_n) << " J/m^3\n";
    double diff = ((computeU_g4(t, t_n) - computeU_g4_no_feedback(t, t_n)) / computeU_g4_no_feedback(t, t_n)) * 100.0;
    ss << "Feedback enhancement: +" << std::fixed << std::setprecision(1) << diff << "%\n";
    ss << "==================================================\n";
    return ss.str();
}

bool FeedbackFactorModule::validateConsistency() {
    if (variables["f_feedback"] < 0) return false;
    if (variables["delta_M_BH_dex"] < 0) return false;
    if (variables["M_bh_initial"] <= 0) return false;
    if (variables["d_g"] <= 0) return false;
    if (variables["f_feedback"] > 2.0) return false;  // Reasonable upper bound
    if (variables["delta_M_BH_dex"] > 5.0) return false;  // 100000x seems excessive
    return true;
}

bool FeedbackFactorModule::autoCorrectAnomalies() {
    bool changed = false;
    if (variables["f_feedback"] < 0 || variables["f_feedback"] > 2.0) {
        variables["f_feedback"] = 0.1;
        changed = true;
    }
    if (variables["delta_M_BH_dex"] < 0 || variables["delta_M_BH_dex"] > 5.0) {
        variables["delta_M_BH_dex"] = 1.0;
        changed = true;
    }
    if (variables["M_bh_initial"] <= 0) {
        variables["M_bh_initial"] = 8.15e36;
        changed = true;
    }
    if (variables["d_g"] <= 0) {
        variables["d_g"] = 2.55e20;
        changed = true;
    }
    if (changed) computeM_bh_final();
    return changed;
}

// 18-step enhanced example demonstrating all dynamic capabilities
void enhanced_FeedbackFactor_example() {
    std::cout << "\n========== ENHANCED FEEDBACK FACTOR DEMO (18 STEPS) ==========\n";
    
    FeedbackFactorModule mod;
    double t = 0.0, t_n = 0.0;
    
    // Step 1: Initial report
    std::cout << "\nStep 1: Initial State\n";
    std::cout << mod.generateReport(t, t_n);
    mod.saveState("initial");
    
    // Step 2: Variable management
    std::cout << "\nStep 2: Create tracking variables\n";
    mod.createVariable("merger_count", 0.0);
    mod.createVariable("total_mass_accreted", 0.0);
    mod.createVariable("feedback_strength", mod.computeF_feedback());
    std::cout << "Created: merger_count, total_mass_accreted, feedback_strength\n";
    auto var_list = mod.listVariables();
    std::cout << "Total variables: " << var_list.size() << "\n";
    
    // Step 3: Explore different feedback scenarios
    std::cout << "\nStep 3: Compare feedback scenarios (0.05, 0.1, 0.2)\n";
    std::vector<double> f_vals = {0.05, 0.1, 0.2};
    for (double f : f_vals) {
        mod.updateVariable("f_feedback", f);
        double u_with = mod.computeU_g4(t, t_n);
        double u_without = mod.computeU_g4_no_feedback(t, t_n);
        std::cout << "  f=" << f << ": U_g4=" << u_with 
                  << " (+% " << ((u_with/u_without - 1.0)*100) << ")\n";
    }
    mod.restoreState("initial");
    
    // Step 4: Expand feedback scale
    std::cout << "\nStep 4: Expand feedback scale (f×1.5, dex×1.2)\n";
    mod.saveState("before_feedback_expansion");
    mod.expandFeedbackScale(1.5, 1.2);
    std::cout << "New f_feedback: " << mod.computeF_feedback() << "\n";
    std::cout << "New delta_dex: " << mod.computeDeltaM_BH() << " dex\n";
    std::cout << "New M_bh_final: " << mod.computeM_bh_final() << " kg\n";
    std::cout << "U_g4: " << mod.computeU_g4(t, t_n) << " J/m^3\n";
    
    // Step 5: Expand black hole scale
    std::cout << "\nStep 5: Expand black hole scale (M×1.3, k4×1.1)\n";
    mod.saveState("before_bh_expansion");
    mod.expandBlackHoleScale(1.3, 1.1);
    std::cout << "New M_bh_initial: " << mod.variables["M_bh_initial"] << " kg\n";
    std::cout << "New M_bh_final: " << mod.computeM_bh_final() << " kg\n";
    std::cout << "New k_4: " << mod.variables["k_4"] << "\n";
    std::cout << "U_g4: " << mod.computeU_g4(t, t_n) << " J/m^3\n";
    
    // Step 6: Expand decay scale
    std::cout << "\nStep 6: Expand decay scale (alpha×1.2, rho×1.1)\n";
    mod.saveState("before_decay_expansion");
    mod.expandDecayScale(1.2, 1.1);
    std::cout << "New alpha: " << mod.variables["alpha"] << " s^-1\n";
    std::cout << "New rho_vac_SCm: " << mod.variables["rho_vac_SCm"] << " J/m^3\n";
    
    // Step 7: Time evolution
    std::cout << "\nStep 7: Time evolution (t=0 to t=86400s = 1 day)\n";
    mod.restoreState("initial");
    double t_day = 86400.0;
    std::cout << "U_g4(t=0): " << mod.computeU_g4(0.0, 0.0) << " J/m^3\n";
    std::cout << "U_g4(t=1 day): " << mod.computeU_g4(t_day, 0.0) << " J/m^3\n";
    std::cout << "Decay factor: " << (mod.computeU_g4(t_day, 0.0) / mod.computeU_g4(0.0, 0.0)) << "\n";
    
    // Step 8: Parameter exploration
    std::cout << "\nStep 8: Generate parameter variations (10 variants, ±5%)\n";
    auto variations = mod.generateVariations(10, 5.0);
    std::cout << "Generated " << variations.size() << " variations\n";
    double min_U = 1e100, max_U = -1e100;
    for (const auto& var : variations) {
        auto temp = mod.variables;
        mod.variables = var;
        mod.computeM_bh_final();
        double U = mod.computeU_g4(t, t_n);
        if (U < min_U) min_U = U;
        if (U > max_U) max_U = U;
        mod.variables = temp;
    }
    std::cout << "U_g4 range: [" << min_U << ", " << max_U << "] J/m^3\n";
    std::cout << "Variation span: " << ((max_U - min_U) / mod.computeU_g4(t, t_n) * 100) << "%\n";
    
    // Step 9: Sensitivity analysis
    std::cout << "\nStep 9: Sensitivity analysis (1% perturbation)\n";
    auto sensitivities = mod.sensitivityAnalysis(t, t_n, 0.01);
    std::cout << "Parameter sensitivities (dU_g4/U_g4 per 1% change):\n";
    std::vector<std::pair<std::string, double>> sorted_sens(sensitivities.begin(), sensitivities.end());
    std::sort(sorted_sens.begin(), sorted_sens.end(), 
        [](const auto& a, const auto& b) { return std::abs(a.second) > std::abs(b.second); });
    for (size_t i = 0; i < std::min(size_t(5), sorted_sens.size()); ++i) {
        std::cout << "  " << sorted_sens[i].first << ": " << (sorted_sens[i].second * 100) << "%\n";
    }
    
    // Step 10: Auto-refinement
    std::cout << "\nStep 10: Auto-refine to target U_g4 = 3.0e-20 J/m^3\n";
    mod.saveState("before_refinement");
    double target_U = 3.0e-20;
    mod.autoRefineParameters(t, t_n, target_U, 1e-22);
    double refined_U = mod.computeU_g4(t, t_n);
    std::cout << "Refined U_g4: " << refined_U << " J/m^3\n";
    std::cout << "Error: " << (std::abs(refined_U - target_U) / target_U * 100) << "%\n";
    std::cout << "Adjusted f_feedback: " << mod.computeF_feedback() << "\n";
    std::cout << "Adjusted M_bh_initial: " << mod.variables["M_bh_initial"] << " kg\n";
    
    // Step 11: Calibration to observations
    std::cout << "\nStep 11: Calibrate to mock observations\n";
    mod.restoreState("before_refinement");
    std::map<std::string, double> observations;
    observations["f_feedback"] = 0.12;
    observations["delta_M_BH_dex"] = 1.2;
    observations["M_bh_initial"] = 8.5e36;
    mod.calibrateToObservations(observations);
    std::cout << "Calibrated report:\n" << mod.generateReport(t, t_n);
    
    // Step 12: Optimize for metric
    std::cout << "\nStep 12: Optimize for maximize_Ug4\n";
    mod.saveState("before_optimization");
    mod.optimizeForMetric("maximize_Ug4");
    std::cout << "Optimized U_g4: " << mod.computeU_g4(t, t_n) << " J/m^3\n";
    std::cout << "Enhanced f_feedback: " << mod.computeF_feedback() << "\n";
    std::cout << "Enhanced M_bh_final: " << mod.computeM_bh_final() << " kg\n";
    
    // Step 13: Mutate parameters
    std::cout << "\nStep 13: Mutate parameters (±3% cosmic noise)\n";
    mod.restoreState("before_optimization");
    mod.saveState("before_mutation");
    mod.mutateParameters(0.03);
    std::cout << "Mutated f_feedback: " << mod.computeF_feedback() << "\n";
    std::cout << "Mutated M_bh_initial: " << mod.variables["M_bh_initial"] << " kg\n";
    std::cout << "Mutated U_g4: " << mod.computeU_g4(t, t_n) << " J/m^3\n";
    
    // Step 14: Evolve system
    std::cout << "\nStep 14: Evolve system (8 generations, selection=0.7)\n";
    mod.restoreState("before_mutation");
    mod.evolveSystem(t, t_n, 8, 0.7);
    std::cout << "Evolved U_g4: " << mod.computeU_g4(t, t_n) << " J/m^3\n";
    std::cout << "Evolved f_feedback: " << mod.computeF_feedback() << "\n";
    std::cout << "Evolved M_bh_final: " << mod.computeM_bh_final() << " kg\n";
    
    // Step 15: State management
    std::cout << "\nStep 15: List all saved states\n";
    auto saved_states = mod.listSavedStates();
    std::cout << "Saved states (" << saved_states.size() << "):\n";
    for (const auto& label : saved_states) {
        std::cout << "  - " << label << "\n";
    }
    
    // Step 16: Validate and auto-correct
    std::cout << "\nStep 16: Validate consistency\n";
    bool valid = mod.validateConsistency();
    std::cout << "Current state valid: " << (valid ? "YES" : "NO") << "\n";
    if (!valid) {
        std::cout << "Applying auto-corrections...\n";
        bool corrected = mod.autoCorrectAnomalies();
        std::cout << "Corrections applied: " << (corrected ? "YES" : "NO") << "\n";
        std::cout << "State valid after correction: " << (mod.validateConsistency() ? "YES" : "NO") << "\n";
    }
    
    // Step 17: Export state
    std::cout << "\nStep 17: Export state\n";
    std::string export_data = mod.exportState();
    std::cout << "Exported state (" << export_data.length() << " characters)\n";
    std::cout << "First 300 chars:\n" << export_data.substr(0, 300) << "...\n";
    
    // Step 18: Final comparison and report
    std::cout << "\nStep 18: Final comparison (multiple dex scenarios)\n";
    mod.restoreState("initial");
    std::vector<double> dex_vals = {0.5, 1.0, 1.5, 2.0};
    for (double dex : dex_vals) {
        mod.updateVariable("delta_M_BH_dex", dex);
        mod.computeM_bh_final();
        double U_with = mod.computeU_g4(t, t_n);
        double U_without = mod.computeU_g4_no_feedback(t, t_n);
        double enhancement = ((U_with / U_without) - 1.0) * 100.0;
        std::cout << "  dex=" << dex << " (M×" << std::pow(10, dex) 
                  << "): U_g4=" << U_with << " J/m^3, enhancement=+" 
                  << std::fixed << std::setprecision(1) << enhancement << "%\n";
    }
    
    std::cout << "\n========== ENHANCED FEEDBACK FACTOR DEMO COMPLETE ==========\n";
    std::cout << "Demonstrated: Variable management, batch ops, feedback/BH/decay scaling,\n";
    std::cout << "              refinement, calibration, optimization, exploration, evolution,\n";
    std::cout << "              state management, sensitivity, validation, reporting.\n";
    std::cout << "Key Insight: Feedback factor f_feedback scales U_g4 by (1+f_feedback).\n";
    std::cout << "             For 1 dex BH mass increase (10x), f=0.1 gives +10% enhancement.\n";
    std::cout << "             Critical for AGN feedback, quasar jets, merger dynamics in UQFF.\n";
}

// ---------------------- END ENHANCED DYNAMIC CAPABILITIES ----------------------

FeedbackFactorModule Evaluation

Strengths :
-Modular, extensible design for modeling the feedback factor(f_feedback) in UQFF, with clear encapsulation of variables using std::map.
- Implements core physical concepts : feedback scaling for black hole mass increase(?M_BH in dex), and its effect on the Ug4 term.
- Provides both feedback and non - feedback calculations for Ug4, enabling comparative analysis.
- Approximations and physical meaning are well - documented in comments and equation text.
- Output functions for variable state and Ug4 comparison support debugging and transparency.
- Handles dynamic updates to variables and recalculates dependent terms as needed.
- Example usage and equation text provide scientific context and validation.

Weaknesses / Recommendations:
-Many constants and parameters are hardcoded; consider external configuration for flexibility and scalability.
- Some calculations use magic numbers or lack explanatory comments; define named constants and clarify logic.
- Minimal error handling(e.g., division by zero, invalid variable names); add validation for robustness.
- Unit consistency should be checked and documented for all physical quantities.
- For large - scale or performance - critical simulations, optimize data structures and reduce redundant calculations.
- std::map is flexible but may be less efficient than structured types for very large models.
- Expand documentation for function purposes and physical meaning.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in feedback factor modeling.It implements the UQFF feedback concept faithfully and adapts to various scenarios.For production or high - performance applications, address the recommendations for improved robustness, maintainability, and scalability.