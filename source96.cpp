// GalacticDistanceModule.h
// Modular C++ implementation of the Distance from Galactic Center (d_g) in the Universal Quantum Field Superconductive Framework (UQFF).
// This module computes d_g=2.55e20 m (~27,000 ly) and conversions; scales M_bh / d_g in U_bi and Ug4.
// Pluggable: #include "GalacticDistanceModule.h"
// GalacticDistanceModule mod; mod.computeMbhOverDg(); mod.updateVariable("d_g", new_value);
// Variables in std::map; example for Sun at t=0, t_n=0.
// Approximations: cos(π t_n)=1; ε_sw * ρ_vac,sw ≈0; α=0.001 s^-1; f_feedback=0.1.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

#ifndef GALACTIC_DISTANCE_MODULE_H
#define GALACTIC_DISTANCE_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <functional>
#include <random>
#include <algorithm>
#include <sstream>

class GalacticDistanceModule {
private:
    std::map<std::string, double> variables;
    double computeDgInLy();
    double computeDgInPc();
    double computeMbhOverDg();
    double computeU_b1();
    double computeU_g4();

public:
    // Constructor: Initialize with framework defaults
    GalacticDistanceModule();

    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    void addToVariable(const std::string& name, double delta);
    void subtractFromVariable(const std::string& name, double delta);

    // Core computations
    double computeDg();  // d_g in m (2.55e20)
    double computeDgInLy();
    double computeDgInPc();
    double computeMbhOverDg();  // M_bh / d_g (kg/m)
    double computeU_b1();  // Universal Buoyancy example (J/m^3)
    double computeU_g4();  // Ug4 example (J/m^3)

    // Output descriptive text
    std::string getEquationText();

    // Print all current variables
    void printVariables();
    
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
    void expandGalacticScale(double d_g_scale, double M_bh_scale);
    void expandBlackHoleScale(double M_bh_scale, double k4_scale);
    void expandEnvironmentScale(double rho_scale, double epsilon_scale);

    // Self-refinement
    void autoRefineParameters(double target_Ub1, double tolerance);
    void calibrateToObservations(const std::map<std::string, double>& observations);
    void optimizeForMetric(const std::string& metric);

    // Exploration & evolution
    std::vector<std::map<std::string,double>> generateVariations(int count, double percent);
    void mutateParameters(double magnitude);
    void evolveSystem(int generations, double selection_pressure);

    // State management
    void saveState(const std::string& label);
    bool restoreState(const std::string& label);
    std::vector<std::string> listSavedStates();
    std::string exportState();

    // Analysis
    std::map<std::string,double> sensitivityAnalysis(double delta);
    std::string generateReport();
    bool validateConsistency();
    bool autoCorrectAnomalies();
};

#endif // GALACTIC_DISTANCE_MODULE_H

// GalacticDistanceModule.cpp
#include "GalacticDistanceModule.h"

// Constructor: Set framework defaults
GalacticDistanceModule::GalacticDistanceModule() {
    // Universal constants
    variables["c"] = 2.998e8;                       // m/s
    variables["year_to_s"] = 3.156e7;               // s/yr
    variables["ly_to_m"] = variables["c"] * variables["year_to_s"];  // ≈9.461e15 m/ly
    variables["pc_to_ly"] = 3.262;                  // ly/pc
    variables["pi"] = 3.141592653589793;

    // Galactic params
    variables["d_g"] = 2.55e20;                     // m (~27,000 ly)
    variables["M_bh"] = 8.15e36;                    // kg (Sgr A* mass)

    // U_bi params
    variables["beta_1"] = 0.6;                      // Unitless
    variables["U_g1"] = 1.39e26;                    // J/m^3
    variables["Omega_g"] = 7.3e-16;                 // rad/s
    variables["epsilon_sw"] = 0.001;                // Unitless
    variables["rho_vac_sw"] = 8e-21;                // J/m^3
    variables["U_UA"] = 1.0;                        // Normalized
    variables["t_n"] = 0.0;                         // s

    // Ug4 params
    variables["k_4"] = 1.0;                         // Unitless
    variables["rho_vac_SCm"] = 7.09e-37;            // J/m^3
    variables["alpha"] = 0.001;                     // s^-1
    variables["f_feedback"] = 0.1;                  // Unitless
}

// Update variable
void GalacticDistanceModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}

// Add delta
void GalacticDistanceModule::addToVariable(const std::string& name, double delta) {
    if (variables.find(name) != variables.end()) {
        variables[name] += delta;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with delta " << delta << std::endl;
        variables[name] = delta;
    }
}

// Subtract delta
void GalacticDistanceModule::subtractFromVariable(const std::string& name, double delta) {
    addToVariable(name, -delta);
}

// Compute d_g (m)
double GalacticDistanceModule::computeDg() {
    return variables["d_g"];
}

// d_g in ly
double GalacticDistanceModule::computeDgInLy() {
    return computeDg() / variables["ly_to_m"];
}

// d_g in pc
double GalacticDistanceModule::computeDgInPc() {
    return computeDgInLy() / variables["pc_to_ly"];
}

// M_bh / d_g (kg/m)
double GalacticDistanceModule::computeMbhOverDg() {
    return variables["M_bh"] / computeDg();
}

// U_b1 example (J/m^3)
double GalacticDistanceModule::computeU_b1() {
    double beta_1 = variables["beta_1"];
    double U_g1 = variables["U_g1"];
    double Omega_g = variables["Omega_g"];
    double mbh_over_dg = computeMbhOverDg();
    double swirl_factor = 1.0 + variables["epsilon_sw"] * variables["rho_vac_sw"];
    double U_UA = variables["U_UA"];
    double cos_term = std::cos(variables["pi"] * variables["t_n"]);
    return -beta_1 * U_g1 * Omega_g * mbh_over_dg * swirl_factor * U_UA * cos_term;
}

// U_g4 example (J/m^3)
double GalacticDistanceModule::computeU_g4() {
    double k_4 = variables["k_4"];
    double rho_vac_SCm = variables["rho_vac_SCm"];
    double mbh_over_dg = computeMbhOverDg();
    double exp_term = std::exp( - variables["alpha"] * variables["t_n"] );
    double cos_term = std::cos(variables["pi"] * variables["t_n"]);
    double feedback_factor = 1.0 + variables["f_feedback"];
    return k_4 * (rho_vac_SCm * variables["M_bh"]) / computeDg() * exp_term * cos_term * feedback_factor;
}

// Equation text
std::string GalacticDistanceModule::getEquationText() {
    return "U_bi = -β_i U_gi Ω_g (M_bh / d_g) (1 + ε_sw ρ_vac,sw) U_UA cos(π t_n)\n"
           "U_g4 = k_4 (ρ_vac,[SCm] M_bh / d_g) e^{-α t} cos(π t_n) (1 + f_feedback)\n"
           "Where d_g = 2.55e20 m (~27,000 ly / 8,260 pc; Sun to Sgr A*).\n"
           "M_bh / d_g ≈3.20e16 kg/m;\n"
           "Example U_b1 ≈ -1.94e27 J/m³; U_g4 ≈2.50e-20 J/m³ (t_n=0).\n"
           "Role: Scales SMBH influence on buoyancy/Ug4; galactic dynamics in nebulae/disks.\n"
           "UQFF: Enables merger resolution (final parsec); star formation modulation.";
}

// Print variables
void GalacticDistanceModule::printVariables() {
    std::cout << "Current Variables:\n";
    for (const auto& pair : variables) {
        std::cout << pair.first << " = " << std::scientific << pair.second << std::endl;
    }
}

// Example usage in base program (snippet)
// #include "GalacticDistanceModule.h"
// int main() {
//     GalacticDistanceModule mod;
//     double dg_ly = mod.computeDgInLy();
//     std::cout << "d_g ≈ " << dg_ly << " ly\n";
//     double u_b1 = mod.computeU_b1();
//     std::cout << "U_b1 = " << u_b1 << " J/m³\n";
//     std::cout << mod.getEquationText() << std::endl;
//     mod.updateVariable("d_g", 2.6e20);
//     mod.printVariables();
//     return 0;
// }
// Compile: g++ -o galactic_test galactic_test.cpp GalacticDistanceModule.cpp -lm
// Sample: d_g ≈2.70e4 ly; U_b1 ≈ -1.94e27 J/m³; scales SMBH at galactic distances.
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.

// ---------------------- ENHANCED DYNAMIC CAPABILITIES IMPLEMENTATION ----------------------

namespace {
    // Simple persistent saved-state storage for the module
    static std::map<std::string, std::map<std::string,double>> galactic_distance_saved_states;
}

void GalacticDistanceModule::createVariable(const std::string& name, double value) {
    variables[name] = value;
}

void GalacticDistanceModule::removeVariable(const std::string& name) {
    auto it = variables.find(name);
    if (it != variables.end()) variables.erase(it);
}

void GalacticDistanceModule::cloneVariable(const std::string& source, const std::string& dest) {
    if (variables.find(source) != variables.end()) variables[dest] = variables[source];
}

std::vector<std::string> GalacticDistanceModule::listVariables() {
    std::vector<std::string> keys;
    for (const auto& p : variables) keys.push_back(p.first);
    return keys;
}

std::string GalacticDistanceModule::getSystemName() {
    return "Galactic_Distance_UQFF";
}

void GalacticDistanceModule::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> fn) {
    for (const auto& n : names) {
        if (variables.find(n) != variables.end()) variables[n] = fn(variables[n]);
    }
}

void GalacticDistanceModule::scaleVariableGroup(const std::vector<std::string>& names, double scale) {
    transformVariableGroup(names, [scale](double v){ return v * scale; });
}

void GalacticDistanceModule::expandParameterSpace(double scale) {
    // Uniform scaling that preserves constants
    std::vector<std::string> exclude = {"c","pi","year_to_s","ly_to_m","pc_to_ly"};
    for (auto& kv : variables) {
        if (std::find(exclude.begin(), exclude.end(), kv.first) == exclude.end()) {
            kv.second *= scale;
        }
    }
}

void GalacticDistanceModule::expandGalacticScale(double d_g_scale, double M_bh_scale) {
    variables["d_g"] *= d_g_scale;
    variables["M_bh"] *= M_bh_scale;
}

void GalacticDistanceModule::expandBlackHoleScale(double M_bh_scale, double k4_scale) {
    variables["M_bh"] *= M_bh_scale;
    variables["k_4"] *= k4_scale;
}

void GalacticDistanceModule::expandEnvironmentScale(double rho_scale, double epsilon_scale) {
    variables["rho_vac_SCm"] *= rho_scale;
    variables["epsilon_sw"] *= epsilon_scale;
}

void GalacticDistanceModule::autoRefineParameters(double target_Ub1, double tolerance) {
    // Very small iterative tuner adjusting beta_1 and M_bh
    double cur = computeU_b1();
    int iter = 0;
    while (std::abs(cur - target_Ub1) > tolerance && iter++ < 50) {
        double factor = target_Ub1 / (cur + 1e-50);
        // conservative adjustments
        variables["beta_1"] *= std::pow(factor, 0.33);
        variables["M_bh"] *= std::pow(factor, 0.33);
        cur = computeU_b1();
    }
}

void GalacticDistanceModule::calibrateToObservations(const std::map<std::string, double>& observations) {
    for (const auto& p : observations) {
        if (variables.find(p.first) != variables.end()) variables[p.first] = p.second;
    }
}

void GalacticDistanceModule::optimizeForMetric(const std::string& metric) {
    if (metric == "maximize_Ub1") {
        variables["beta_1"] *= 1.5;
        variables["M_bh"] *= 1.3;
    } else if (metric == "minimize_Ug4") {
        variables["k_4"] *= 0.5;
        variables["alpha"] *= 1.2; // faster decay
    }
}

std::vector<std::map<std::string,double>> GalacticDistanceModule::generateVariations(int count, double percent) {
    std::vector<std::map<std::string,double>> out;
    std::random_device rd; std::mt19937 gen(rd());
    std::vector<std::string> exclude = {"c","pi","year_to_s","ly_to_m","pc_to_ly"};
    for (int i=0;i<count;i++) {
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

void GalacticDistanceModule::mutateParameters(double magnitude) {
    std::random_device rd; std::mt19937 gen(rd());
    std::normal_distribution<> d(0.0, magnitude);
    std::vector<std::string> exclude = {"c","pi","year_to_s","ly_to_m","pc_to_ly"};
    for (auto& kv : variables) {
        if (std::find(exclude.begin(), exclude.end(), kv.first) != exclude.end()) continue;
        double delta = d(gen) * kv.second;
        kv.second += delta;
    }
}

void GalacticDistanceModule::evolveSystem(int generations, double selection_pressure) {
    int pop = 20;
    std::vector<std::map<std::string,double>> population;
    for (int i=0;i<pop;i++) population.push_back(variables);
    std::random_device rd; std::mt19937 gen(rd());
    for (int g=0; g<generations; g++) {
        // mutate all
        for (auto& ind : population) {
            std::normal_distribution<> d(0.0, 0.02);
            for (auto& kv : ind) {
                if (kv.first=="c"||kv.first=="pi") continue;
                kv.second *= (1.0 + d(gen));
            }
        }
        // score by absolute U_b1
        std::vector<std::pair<double,int>> scores;
        for (int i=0;i<pop;i++) {
            // temporary evaluate
            auto temp = variables;
            variables = population[i];
            double score = std::abs(computeU_b1());
            scores.push_back({score,i});
            variables = temp;
        }
        std::sort(scores.begin(), scores.end(), std::greater<>());
        // select top based on selection pressure
        int survivors = std::max(2, (int)(pop * selection_pressure));
        std::vector<std::map<std::string,double>> next;
        for (int i=0;i<survivors;i++) next.push_back(population[scores[i].second]);
        // refill by mutating survivors
        std::uniform_int_distribution<> uid(0,survivors-1);
        while ((int)next.size() < pop) {
            auto child = next[uid(gen)];
            std::normal_distribution<> d(0.0, 0.05);
            for (auto& kv : child) kv.second *= (1.0 + d(gen));
            next.push_back(child);
        }
        population = next;
    }
    // set variables to best
    auto best = population[0];
    variables = best;
}

void GalacticDistanceModule::saveState(const std::string& label) {
    galactic_distance_saved_states[label] = variables;
}

bool GalacticDistanceModule::restoreState(const std::string& label) {
    auto it = galactic_distance_saved_states.find(label);
    if (it == galactic_distance_saved_states.end()) return false;
    variables = it->second;
    return true;
}

std::vector<std::string> GalacticDistanceModule::listSavedStates() {
    std::vector<std::string> keys;
    for (const auto& p : galactic_distance_saved_states) keys.push_back(p.first);
    return keys;
}

std::string GalacticDistanceModule::exportState() {
    std::ostringstream ss;
    ss << "GalacticDistanceModule State Export\n";
    for (const auto& kv : variables) ss << kv.first << ": " << kv.second << "\n";
    return ss.str();
}

std::map<std::string,double> GalacticDistanceModule::sensitivityAnalysis(double delta) {
    std::map<std::string,double> result;
    double base = computeU_b1();
    std::vector<std::string> keys = {"d_g","M_bh","beta_1","U_g1","Omega_g","epsilon_sw","rho_vac_SCm","alpha","k_4"};
    for (const auto& k : keys) {
        if (variables.find(k) == variables.end()) continue;
        double orig = variables[k];
        variables[k] = orig * (1.0 + delta);
        double up = computeU_b1();
        variables[k] = orig * (1.0 - delta);
        double down = computeU_b1();
        variables[k] = orig;
        double sens = 0.5 * ((up - down) / (base + 1e-50));
        result[k] = sens;
    }
    return result;
}

std::string GalacticDistanceModule::generateReport() {
    std::ostringstream ss;
    ss << "GalacticDistanceModule Report\n";
    ss << "d_g (m): " << std::scientific << variables["d_g"] << "\n";
    ss << "d_g (ly): " << std::fixed << computeDgInLy() << " ly\n";
    ss << "d_g (pc): " << std::fixed << computeDgInPc() << " pc\n";
    ss << "M_bh: " << std::scientific << variables["M_bh"] << " kg\n";
    ss << "M_bh/d_g: " << std::scientific << computeMbhOverDg() << " kg/m\n";
    ss << "U_b1: " << std::scientific << computeU_b1() << " J/m^3\n";
    ss << "U_g4: " << std::scientific << computeU_g4() << " J/m^3\n";
    return ss.str();
}

bool GalacticDistanceModule::validateConsistency() {
    if (variables["d_g"] <= 0) return false;
    if (variables["M_bh"] <= 0) return false;
    if (variables["d_g"] < 1e18 || variables["d_g"] > 1e23) return false;
    return true;
}

bool GalacticDistanceModule::autoCorrectAnomalies() {
    bool changed = false;
    if (variables["d_g"] < 1e18 || variables["d_g"] > 1e23) { variables["d_g"] = 2.55e20; changed = true; }
    if (variables["M_bh"] <= 0) { variables["M_bh"] = 8.15e36; changed = true; }
    if (variables["rho_vac_SCm"] <= 0) { variables["rho_vac_SCm"] = 7.09e-37; changed = true; }
    return changed;
}

// 18-step enhanced example demonstrating dynamic capabilities
void enhanced_GalacticDistance_example() {
    std::cout << "\n========== ENHANCED GALACTIC DISTANCE DEMO (18 STEPS) ==========\n";
    GalacticDistanceModule mod;
    mod.saveState("initial");

    // Step 1: report initial state
    std::cout << "Step 1: Initial report:\n" << mod.generateReport() << std::endl;

    // Step 2: create observational overrides
    std::cout << "Step 2: Add observation variables and small perturbs\n";
    mod.createVariable("observed_d_g", mod.computeDg()*1.02);
    mod.createVariable("observed_M_bh", mod.computeDg()*0.0 + mod.variables["M_bh"]*1.01);

    // Step 3: convert units
    std::cout << "Step 3: d_g in ly/pc: " << mod.computeDgInLy() << " ly, " << mod.computeDgInPc() << " pc\n";

    // Step 4: expand galactic scale (simulate larger halo)
    std::cout << "Step 4: Expand galactic scale (d_g x1.1, M_bh x1.05)\n";
    mod.expandGalacticScale(1.1, 1.05);
    std::cout << "d_g now: " << mod.computeDg() << " m\n";

    // Step 5: expand black hole and environment
    std::cout << "Step 5: Expand black hole scale and environment\n";
    mod.expandBlackHoleScale(1.2, 1.1);
    mod.expandEnvironmentScale(1.05, 1.02);

    // Step 6: run sensitivity analysis
    std::cout << "Step 6: Sensitivity analysis (1% delta)\n";
    auto sens = mod.sensitivityAnalysis(0.01);
    for (const auto& s : sens) std::cout << "  " << s.first << ": " << s.second << "\n";

    // Step 7: generate variations
    std::cout << "Step 7: Generate variations (5 variants, 3%)\n";
    auto vars = mod.generateVariations(5, 3.0);
    std::cout << "  Variants generated: " << vars.size() << "\n";

    // Step 8: auto refine towards a target U_b1
    std::cout << "Step 8: Auto-refine parameters toward target U_b1 = -1e27 J/m^3\n";
    mod.autoRefineParameters(-1e27, 1e24);
    std::cout << "  Refined U_b1: " << mod.computeU_b1() << "\n";

    // Step 9: calibrate to mock observations
    std::cout << "Step 9: Calibrate to mock observations\n";
    std::map<std::string,double> obs; obs["d_g"] = 2.6e20; obs["M_bh"] = 8.3e36;
    mod.calibrateToObservations(obs);
    std::cout << "  Calibrated report:\n" << mod.generateReport() << std::endl;

    // Step 10: save a state for comparison
    mod.saveState("post_calibration");

    // Step 11: optimize metric (maximize U_b1)
    std::cout << "Step 11: Optimize for maximize_Ub1\n";
    mod.optimizeForMetric("maximize_Ub1");
    std::cout << "  Optimized U_b1: " << mod.computeU_b1() << "\n";

    // Step 12: mutate parameters (cosmic noise)
    std::cout << "Step 12: Mutate parameters (2%)\n";
    mod.mutateParameters(0.02);
    std::cout << "  Mutated U_b1: " << mod.computeU_b1() << "\n";

    // Step 13: evolve system (5 generations)
    std::cout << "Step 13: Evolve system (5 gens)\n";
    mod.evolveSystem(5, 0.6);
    std::cout << "  Evolved U_b1: " << mod.computeU_b1() << "\n";

    // Step 14: list saved states
    std::cout << "Step 14: Saved states:\n";
    for (auto s : mod.listSavedStates()) std::cout << "  - " << s << "\n";

    // Step 15: restore initial and compare
    std::cout << "Step 15: Restore initial state and compare\n";
    mod.restoreState("initial");
    std::cout << "  After restore U_b1: " << mod.computeU_b1() << "\n";

    // Step 16: validate and auto-correct anomalies
    std::cout << "Step 16: Validate consistency: " << (mod.validateConsistency() ? "OK" : "BAD") << "\n";
    if (!mod.validateConsistency()) {
        mod.autoCorrectAnomalies();
        std::cout << "  Corrections applied. New U_b1: " << mod.computeU_b1() << "\n";
    }

    // Step 17: export state and show a snippet
    std::cout << "Step 17: Export state (snippet):\n" << mod.exportState().substr(0,400) << "\n";

    // Step 18: final report
    std::cout << "Step 18: Final report:\n" << mod.generateReport() << std::endl;

    std::cout << "\n========== ENHANCED GALACTIC DISTANCE DEMO COMPLETE ==========\n";
}

// ---------------------- END ENHANCED DYNAMIC CAPABILITIES ----------------------

GalacticDistanceModule Evaluation

Strengths :
 - Modular, extensible design for modeling galactic center distance and its role in UQFF calculations.
 - Clear encapsulation of variables using std::map, supporting dynamic updates and easy extension.
 - Implements core physical concepts : conversion of d_g between units(m, ly, pc), SMBH scaling(M_bh / d_g), and contributions to universal buoyancy(U_b1) and gravity(Ug4).
 - Approximations and physical meaning are well - documented in comments and equation text.
 - Output functions for variable state and equation text support debugging and transparency.
 - Handles dynamic updates to variables and recalculates dependent terms as needed.
 - Example calculations and conversion functions provide scientific context and validation.

Weaknesses / Recommendations:
 - Many constants and parameters are hardcoded; consider external configuration for flexibility and scalability.
 - Some calculations use magic numbers or lack explanatory comments; define named constants and clarify logic.
 - Minimal error handling(e.g., division by zero, invalid variable names); add validation for robustness.
 - Unit consistency should be checked and documented for all physical quantities.
 - For large - scale or performance - critical simulations, optimize data structures and reduce redundant calculations.
 - std::map is flexible but may be less efficient than structured types for very large models.
 - Expand documentation for function purposes and physical meaning.

Summary:
The code is well - structured, clear, and suitable for scientific prototyping and educational use in galactic distance modeling. It implements the UQFF distance concept faithfully and adapts to various scenarios. For production or high - performance applications, address the recommendations for improved robustness, maintainability, and scalability.