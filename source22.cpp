/**
 * ================================================================================================
 * Header: BubbleNebula.h
 *
 * Description: C++ Module for Bubble Nebula (NGC 7635) Class
 *              This is the twelfth module in a series of 500+ code files for the Universal Quantum
 *              Field Framework (UQFF) simulations, focusing on emission nebula evolution and
 *              gravity equations derived from Hubble datasets, high-energy lab simulations, and
 *              UQFF refinements (dated May 09, 2025, updated for full term inclusion on October 08, 2025).
 *
 * Purpose: Encapsulates the Master Universal Gravity Equation (MUGE) for Bubble Nebula evolution.
 *          Includes ALL terms: base gravity (static M), cosmic expansion (H_0), magnetic correction (static B),
 *          expansion E(t), UQFF Ug components with f_TRZ, Lambda, quantum uncertainty, scaled EM with [UA],
 *          fluid dynamics, oscillatory waves, DM/density perturbations, and stellar wind feedback (pressure / density for acc).
 *          Supports dynamic variable updates for all parameters.
 *
 * Integration: Designed for inclusion in base program 'ziqn233h.cpp' (not present here).
 *              Instantiate class in main: BubbleNebula bubble;
 *              Compute: double g = bubble.compute_g_Bubble(t);
 *
 * Key Features:
 *   - Default values from UQFF document: M = 46 * Msun, r = 4.731e16 m (5 ly), B = 1e-6 T,
 *     E_0 = 0.1, tau_exp = 4 Myr, rho_wind = 1e-21 kg/m^3, v_wind = 1.8e6 m/s.
 *   - Units handled: Msun to kg, ly to m; wind term as (rho * v_wind^2) / rho_fluid for acceleration.
 *   - Setter methods for updates: setVar(double new_val) or addToVar(double delta)/subtractFromVar(double delta).
 *   - Computes g_Bubble(r, t) with every term explicitly included.
 *
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: October 08, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

#ifndef BUBBLE_NEBULA_H
#define BUBBLE_NEBULA_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <map>
#include <string>
#include <functional>
#include <random>
#include <algorithm>
#include <sstream>
#include <vector>

class BubbleNebula {
private:
    // Core parameters (mutable for updates)
    double G;               // Gravitational constant
    double M;               // Total mass (kg)
    double r;               // Radius (m)
    double H0;              // Hubble constant (s^-1)
    double B;               // Static magnetic field (T)
    double B_crit;          // Critical B field (T)
    double Lambda;          // Cosmological constant
    double c_light;         // Speed of light
    double q_charge;        // Charge (proton)
    double gas_v;           // Gas velocity for EM (m/s)
    double f_TRZ;           // Time-reversal factor
    double E_0;             // Initial expansion factor
    double tau_exp;         // Expansion timescale (s)
    double rho_wind;        // Wind density (kg/m^3)
    double v_wind;          // Wind velocity (m/s)
    double rho_fluid;       // Fluid density (kg/m^3)
    double rho_vac_UA;      // UA vacuum density (J/m^3)
    double rho_vac_SCm;     // SCm vacuum density (J/m^3)
    double scale_EM;        // EM scaling factor
    double proton_mass;     // Proton mass for EM acceleration

    // Additional parameters for full inclusion of terms
    double hbar;            // Reduced Planck's constant
    double t_Hubble;        // Hubble time (s)
    double delta_x;         // Position uncertainty (m)
    double delta_p;         // Momentum uncertainty (kg m/s)
    double integral_psi;    // Wavefunction integral approximation
    double A_osc;           // Oscillatory amplitude (m/s^2)
    double k_osc;           // Wave number (1/m)
    double omega_osc;       // Angular frequency (rad/s)
    double x_pos;           // Position for oscillation (m)
    double t_Hubble_gyr;    // Hubble time in Gyr
    double M_DM_factor;     // Dark matter mass fraction
    double delta_rho_over_rho; // Density perturbation fraction

    // Computed caches (updated on demand)
    double ug1_base;        // Cached Ug1 = G*M/r^2

public:
    // Constructor with default UQFF values
    BubbleNebula() {
        initializeDefaults();
    }

    // Destructor (empty)
    ~BubbleNebula() {}

    // Initialization method (called in constructor)
    void initializeDefaults() {
        G = 6.6743e-11;
        double M_sun = 1.989e30;
        M = 46.0 * M_sun;
        double ly_to_m = 9.461e15;
        r = 5.0 * ly_to_m;
        H0 = 2.184e-18;
        B = 1e-6;
        B_crit = 1e11;
        Lambda = 1.1e-52;
        c_light = 3e8;
        q_charge = 1.602e-19;
        gas_v = 1e5;
        f_TRZ = 0.1;
        E_0 = 0.1;
        tau_exp = 4e6 * 3.156e7;
        rho_wind = 1e-21;
        v_wind = 1.8e6;
        rho_fluid = 1e-21;
        rho_vac_UA = 7.09e-36;
        rho_vac_SCm = 7.09e-37;
        scale_EM = 1e-12;
        proton_mass = 1.673e-27;

        // Full terms defaults
        hbar = 1.0546e-34;
        t_Hubble = 13.8e9 * 3.156e7;
        t_Hubble_gyr = 13.8;
        delta_x = 1e-10;
        delta_p = hbar / delta_x;
        integral_psi = 1.0;
        A_osc = 1e-10;
        k_osc = 1.0 / r;
        omega_osc = 2 * M_PI / (r / c_light);
        x_pos = r;
        M_DM_factor = 0.1;
        delta_rho_over_rho = 1e-5;

        updateCache();
    }

    // Cache update for efficiency (call after parameter changes)
    void updateCache() {
        ug1_base = (G * M) / (r * r);
    }

    // Universal setter for any variable (by name, for flexibility)
    bool setVariable(const std::string& varName, double newValue) {
        if (varName == "G") { G = newValue; }
        else if (varName == "M") { M = newValue; }
        else if (varName == "r") { r = newValue; }
        else if (varName == "H0") { H0 = newValue; }
        else if (varName == "B") { B = newValue; }
        else if (varName == "B_crit") { B_crit = newValue; }
        else if (varName == "Lambda") { Lambda = newValue; }
        else if (varName == "c_light") { c_light = newValue; }
        else if (varName == "q_charge") { q_charge = newValue; }
        else if (varName == "gas_v") { gas_v = newValue; }
        else if (varName == "f_TRZ") { f_TRZ = newValue; }
        else if (varName == "E_0") { E_0 = newValue; }
        else if (varName == "tau_exp") { tau_exp = newValue; }
        else if (varName == "rho_wind") { rho_wind = newValue; }
        else if (varName == "v_wind") { v_wind = newValue; }
        else if (varName == "rho_fluid") { rho_fluid = newValue; }
        else if (varName == "rho_vac_UA") { rho_vac_UA = newValue; }
        else if (varName == "rho_vac_SCm") { rho_vac_SCm = newValue; }
        else if (varName == "scale_EM") { scale_EM = newValue; }
        else if (varName == "proton_mass") { proton_mass = newValue; }
        // Full terms
        else if (varName == "hbar") { hbar = newValue; }
        else if (varName == "t_Hubble") { t_Hubble = newValue; }
        else if (varName == "t_Hubble_gyr") { t_Hubble_gyr = newValue; }
        else if (varName == "delta_x") { delta_x = newValue; }
        else if (varName == "delta_p") { delta_p = newValue; }
        else if (varName == "integral_psi") { integral_psi = newValue; }
        else if (varName == "A_osc") { A_osc = newValue; }
        else if (varName == "k_osc") { k_osc = newValue; }
        else if (varName == "omega_osc") { omega_osc = newValue; }
        else if (varName == "x_pos") { x_pos = newValue; }
        else if (varName == "M_DM_factor") { M_DM_factor = newValue; }
        else if (varName == "delta_rho_over_rho") { delta_rho_over_rho = newValue; }
        else {
            std::cerr << "Error: Unknown variable '" << varName << "'." << std::endl;
            return false;
        }
        updateCache();
        return true;
    }

    // Addition method for variables
    bool addToVariable(const std::string& varName, double delta) {
        return setVariable(varName, getVariable(varName) + delta);
    }

    // Subtraction method for variables
    bool subtractFromVariable(const std::string& varName, double delta) {
        return addToVariable(varName, -delta);
    }

    // Getter for any variable (helper for add/subtract)
    double getVariable(const std::string& varName) const {
        if (varName == "G") return G;
        else if (varName == "M") return M;
        else if (varName == "r") return r;
        else if (varName == "H0") return H0;
        else if (varName == "B") return B;
        else if (varName == "B_crit") return B_crit;
        else if (varName == "Lambda") return Lambda;
        else if (varName == "c_light") return c_light;
        else if (varName == "q_charge") return q_charge;
        else if (varName == "gas_v") return gas_v;
        else if (varName == "f_TRZ") return f_TRZ;
        else if (varName == "E_0") return E_0;
        else if (varName == "tau_exp") return tau_exp;
        else if (varName == "rho_wind") return rho_wind;
        else if (varName == "v_wind") return v_wind;
        else if (varName == "rho_fluid") return rho_fluid;
        else if (varName == "rho_vac_UA") return rho_vac_UA;
        else if (varName == "rho_vac_SCm") return rho_vac_SCm;
        else if (varName == "scale_EM") return scale_EM;
        else if (varName == "proton_mass") return proton_mass;
        // Full terms
        else if (varName == "hbar") return hbar;
        else if (varName == "t_Hubble") return t_Hubble;
        else if (varName == "t_Hubble_gyr") return t_Hubble_gyr;
        else if (varName == "delta_x") return delta_x;
        else if (varName == "delta_p") return delta_p;
        else if (varName == "integral_psi") return integral_psi;
        else if (varName == "A_osc") return A_osc;
        else if (varName == "k_osc") return k_osc;
        else if (varName == "omega_osc") return omega_osc;
        else if (varName == "x_pos") return x_pos;
        else if (varName == "M_DM_factor") return M_DM_factor;
        else if (varName == "delta_rho_over_rho") return delta_rho_over_rho;
        else {
            std::cerr << "Error: Unknown variable '" << varName << "'." << std::endl;
            return 0.0;
        }
    }

    // E(t) computation
    double E_t(double t) const {
        return E_0 * (1 - exp(-t / tau_exp));
    }

    // Ug terms computation
    double compute_Ug(double Et) const {
        double Ug1 = ug1_base;
        double Ug2 = 0.0;
        double Ug3 = 0.0;
        double corr_B = 1 - B / B_crit;
        double Ug4 = Ug1 * corr_B;
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1 + f_TRZ) * (1 - Et);
    }

    // Volume computation for fluid
    double compute_V() const {
        return (4.0 / 3.0) * M_PI * r * r * r;
    }

    // Main MUGE computation (includes ALL terms)
    double compute_g_Bubble(double t) const {
        if (t < 0) {
            std::cerr << "Error: Time t must be non-negative." << std::endl;
            return 0.0;
        }

        double Et = E_t(t);

        // Term 1: Base + H0 + B + E corrections
        double corr_H = 1 + H0 * t;
        double corr_B = 1 - B / B_crit;
        double corr_E = 1 - Et;
        double term1 = ug1_base * corr_H * corr_B * corr_E;

        // Term 2: UQFF Ug with f_TRZ and E
        double term2 = compute_Ug(Et);

        // Term 3: Lambda
        double term3 = (Lambda * c_light * c_light) / 3.0;

        // Term 4: Scaled EM with UA
        double cross_vB = gas_v * B;  // Magnitude, assuming perpendicular
        double em_base = (q_charge * cross_vB) / proton_mass;
        double corr_UA = 1 + (rho_vac_UA / rho_vac_SCm);
        double term4 = (em_base * corr_UA) * scale_EM;

        // Quantum uncertainty term
        double sqrt_unc = sqrt(delta_x * delta_p);
        double term_q = (hbar / sqrt_unc) * integral_psi * (2 * M_PI / t_Hubble);

        // Fluid term (effective acceleration)
        double V = compute_V();
        double term_fluid = (rho_fluid * V * ug1_base) / M;

        // Oscillatory terms (real parts)
        double term_osc1 = 2 * A_osc * cos(k_osc * x_pos) * cos(omega_osc * t);
        double arg = k_osc * x_pos - omega_osc * t;
        double term_osc2 = (2 * M_PI / t_Hubble_gyr) * A_osc * cos(arg);
        double term_osc = term_osc1 + term_osc2;

        // DM and density perturbation term (converted to acceleration)
        double M_dm = M * M_DM_factor;
        double pert1 = delta_rho_over_rho;
        double pert2 = 3 * G * M / (r * r * r);
        double term_dm_force_like = (M + M_dm) * (pert1 + pert2);
        double term_DM = term_dm_force_like / M;

        // Stellar wind feedback term (pressure / density for acceleration)
        double wind_pressure = rho_wind * v_wind * v_wind;
        double term_wind = wind_pressure / rho_fluid;

        // Total g_Bubble (all terms summed)
        return term1 + term2 + term3 + term4 + term_q + term_fluid + term_osc + term_DM + term_wind;
    }

    // Debug/Output method (for transparency in base program)
    void printParameters(std::ostream& os = std::cout) const {
        os << std::fixed << std::setprecision(3);
        os << "Bubble Nebula Parameters:" << std::endl;
        os << "G: " << G << ", M: " << M << ", r: " << r << std::endl;
        os << "H0: " << H0 << ", B: " << B << ", B_crit: " << B_crit << std::endl;
        os << "f_TRZ: " << f_TRZ << ", E_0: " << E_0 << ", tau_exp: " << tau_exp << std::endl;
        os << "rho_fluid: " << rho_fluid << ", rho_wind: " << rho_wind << ", v_wind: " << v_wind << std::endl;
        os << "gas_v: " << gas_v << ", M_DM_factor: " << M_DM_factor << std::endl;
        os << "A_osc: " << A_osc << ", delta_rho_over_rho: " << delta_rho_over_rho << std::endl;
        os << "ug1_base: " << ug1_base << std::endl;
    }

    // Example computation at t=2 Myr (for testing)
    double exampleAt2Myr() const {
        double t_example = 2e6 * 3.156e7;
        return compute_g_Bubble(t_example);
    }

    // ========== ENHANCED DYNAMIC CAPABILITIES (25 methods) ==========

    // --- Variable Management (5 methods) ---
    bool createVariable(const std::string& name, double value);
    bool removeVariable(const std::string& name);
    bool cloneVariable(const std::string& src, const std::string& dest);
    std::vector<std::string> listVariables() const;
    std::string getSystemName() const;

    // --- Batch Operations (2 methods) ---
    bool transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func);
    bool scaleVariableGroup(const std::vector<std::string>& names, double factor);

    // --- Self-Expansion (4 methods) ---
    void expandParameterSpace(double factor);
    void expandNebulaScale(double M_scale, double r_scale);
    void expandExpansionScale(double E_0_scale, double tau_exp_scale);
    void expandWindMagneticScale(double rho_wind_scale, double v_wind_scale, double B_scale);

    // --- Self-Refinement (3 methods) ---
    void autoRefineParameters(const std::vector<std::pair<double, double>>& observations);
    void calibrateToObservations(const std::vector<double>& times, const std::vector<double>& g_obs);
    double optimizeForMetric(std::function<double(double)> metric, double t_start, double t_end, int steps);

    // --- Parameter Exploration (1 method) ---
    std::vector<std::map<std::string, double>> generateVariations(int count, double variation_pct);

    // --- Adaptive Evolution (2 methods) ---
    void mutateParameters(double mutation_rate, std::mt19937& rng);
    void evolveSystem(int generations, std::function<double(const BubbleNebula&)> fitness);

    // --- State Management (4 methods) ---
    bool saveState(const std::string& stateName);
    bool restoreState(const std::string& stateName);
    std::vector<std::string> listSavedStates() const;
    std::string exportState() const;

    // --- System Analysis (4 methods) ---
    std::map<std::string, double> sensitivityAnalysis(double t, double delta_pct);
    std::string generateReport(double t) const;
    bool validateConsistency() const;
    bool autoCorrectAnomalies();
};

#endif // BUBBLE_NEBULA_H

// ========== IMPLEMENTATION OF ENHANCED METHODS (Outside class) ==========

// Anonymous namespace for state storage
namespace {
    std::map<std::string, std::map<std::string, double>> bubble_saved_states;
}

// --- Variable Management (5 methods) ---
bool BubbleNebula::createVariable(const std::string& name, double value) {
    return setVariable(name, value);
}

bool BubbleNebula::removeVariable(const std::string& name) {
    std::cerr << "Warning: Cannot remove built-in variable '" << name << "' in BubbleNebula class." << std::endl;
    return false;
}

bool BubbleNebula::cloneVariable(const std::string& src, const std::string& dest) {
    double val = getVariable(src);
    return setVariable(dest, val);
}

std::vector<std::string> BubbleNebula::listVariables() const {
    return {"G", "M", "r", "H0", "B", "B_crit", "Lambda", "c_light", "q_charge",
            "gas_v", "f_TRZ", "E_0", "tau_exp", "rho_wind", "v_wind", "rho_fluid",
            "rho_vac_UA", "rho_vac_SCm", "scale_EM", "proton_mass", "hbar", "t_Hubble",
            "t_Hubble_gyr", "delta_x", "delta_p", "integral_psi", "A_osc", "k_osc",
            "omega_osc", "x_pos", "M_DM_factor", "delta_rho_over_rho"};
}

std::string BubbleNebula::getSystemName() const {
    return "BubbleNebula";
}

// --- Batch Operations (2 methods) ---
bool BubbleNebula::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        double val = getVariable(name);
        if (!setVariable(name, func(val))) return false;
    }
    return true;
}

bool BubbleNebula::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    return transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// --- Self-Expansion (4 methods) ---
void BubbleNebula::expandParameterSpace(double factor) {
    std::vector<std::string> expandable = {"M", "r", "B", "rho_fluid", "rho_wind", "A_osc", "M_DM_factor"};
    scaleVariableGroup(expandable, factor);
}

void BubbleNebula::expandNebulaScale(double M_scale, double r_scale) {
    setVariable("M", getVariable("M") * M_scale);
    setVariable("r", getVariable("r") * r_scale);
}

void BubbleNebula::expandExpansionScale(double E_0_scale, double tau_exp_scale) {
    setVariable("E_0", getVariable("E_0") * E_0_scale);
    setVariable("tau_exp", getVariable("tau_exp") * tau_exp_scale);
}

void BubbleNebula::expandWindMagneticScale(double rho_wind_scale, double v_wind_scale, double B_scale) {
    setVariable("rho_wind", getVariable("rho_wind") * rho_wind_scale);
    setVariable("v_wind", getVariable("v_wind") * v_wind_scale);
    setVariable("B", getVariable("B") * B_scale);
    setVariable("B_crit", getVariable("B_crit") * B_scale);
}

// --- Self-Refinement (3 methods) ---
void BubbleNebula::autoRefineParameters(const std::vector<std::pair<double, double>>& observations) {
    if (observations.empty()) return;
    
    double sum_error = 0.0;
    for (const auto& obs : observations) {
        double t = obs.first;
        double g_obs = obs.second;
        double g_calc = compute_g_Bubble(t);
        sum_error += std::abs(g_calc - g_obs);
    }
    double avg_error = sum_error / observations.size();
    
    if (avg_error > 1e-6) {
        double adj_factor = 1.0 - std::min(0.1, avg_error / 1e6);
        setVariable("E_0", getVariable("E_0") * adj_factor);
        setVariable("tau_exp", getVariable("tau_exp") * (2.0 - adj_factor));
    }
}

void BubbleNebula::calibrateToObservations(const std::vector<double>& times, const std::vector<double>& g_obs) {
    if (times.size() != g_obs.size() || times.empty()) return;
    
    std::vector<std::pair<double, double>> obs;
    for (size_t i = 0; i < times.size(); ++i) {
        obs.push_back({times[i], g_obs[i]});
    }
    
    for (int iter = 0; iter < 5; ++iter) {
        autoRefineParameters(obs);
    }
}

double BubbleNebula::optimizeForMetric(std::function<double(double)> metric, double t_start, double t_end, int steps) {
    double best_score = -1e100;
    double dt = (t_end - t_start) / steps;
    
    for (int i = 0; i <= steps; ++i) {
        double t = t_start + i * dt;
        double g = compute_g_Bubble(t);
        double score = metric(g);
        if (score > best_score) best_score = score;
    }
    return best_score;
}

// --- Parameter Exploration (1 method) ---
std::vector<std::map<std::string, double>> BubbleNebula::generateVariations(int count, double variation_pct) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-variation_pct/100.0, variation_pct/100.0);
    
    auto vars = listVariables();
    for (int i = 0; i < count; ++i) {
        std::map<std::string, double> variant;
        for (const auto& v : vars) {
            double val = getVariable(v);
            double variation = val * (1.0 + dis(gen));
            variant[v] = variation;
        }
        variations.push_back(variant);
    }
    return variations;
}

// --- Adaptive Evolution (2 methods) ---
void BubbleNebula::mutateParameters(double mutation_rate, std::mt19937& rng) {
    std::uniform_real_distribution<> dis(-mutation_rate, mutation_rate);
    auto vars = listVariables();
    
    for (const auto& v : vars) {
        if (v == "c_light" || v == "G" || v == "hbar") continue;
        double val = getVariable(v);
        double delta = val * dis(rng);
        setVariable(v, val + delta);
    }
}

void BubbleNebula::evolveSystem(int generations, std::function<double(const BubbleNebula&)> fitness) {
    std::random_device rd;
    std::mt19937 rng(rd());
    
    double best_fitness = fitness(*this);
    saveState("evolution_best");
    
    for (int gen = 0; gen < generations; ++gen) {
        saveState("evolution_temp");
        mutateParameters(0.05, rng);
        
        double new_fitness = fitness(*this);
        if (new_fitness > best_fitness) {
            best_fitness = new_fitness;
            saveState("evolution_best");
        } else {
            restoreState("evolution_temp");
        }
    }
    
    restoreState("evolution_best");
}

// --- State Management (4 methods) ---
bool BubbleNebula::saveState(const std::string& stateName) {
    std::map<std::string, double> state;
    auto vars = listVariables();
    for (const auto& v : vars) {
        state[v] = getVariable(v);
    }
    bubble_saved_states[stateName] = state;
    return true;
}

bool BubbleNebula::restoreState(const std::string& stateName) {
    auto it = bubble_saved_states.find(stateName);
    if (it == bubble_saved_states.end()) return false;
    
    for (const auto& pair : it->second) {
        setVariable(pair.first, pair.second);
    }
    return true;
}

std::vector<std::string> BubbleNebula::listSavedStates() const {
    std::vector<std::string> names;
    for (const auto& pair : bubble_saved_states) {
        names.push_back(pair.first);
    }
    return names;
}

std::string BubbleNebula::exportState() const {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(6);
    oss << "BubbleNebula State Export:\n";
    auto vars = listVariables();
    for (const auto& v : vars) {
        oss << v << " = " << getVariable(v) << "\n";
    }
    return oss.str();
}

// --- System Analysis (4 methods) ---
std::map<std::string, double> BubbleNebula::sensitivityAnalysis(double t, double delta_pct) {
    std::map<std::string, double> sensitivities;
    double g_base = compute_g_Bubble(t);
    
    auto vars = listVariables();
    for (const auto& v : vars) {
        if (v == "c_light" || v == "G" || v == "hbar") continue;
        
        double original = getVariable(v);
        double delta = original * delta_pct / 100.0;
        
        setVariable(v, original + delta);
        double g_plus = compute_g_Bubble(t);
        setVariable(v, original);
        
        double sensitivity = (g_base != 0.0) ? std::abs((g_plus - g_base) / g_base) : 0.0;
        sensitivities[v] = sensitivity;
    }
    
    return sensitivities;
}

std::string BubbleNebula::generateReport(double t) const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6);
    oss << "============================================\n";
    oss << "BUBBLE NEBULA REPORT\n";
    oss << "NGC 7635 Emission Nebula\n";
    oss << "============================================\n";
    oss << "Time: t = " << t << " s (" << (t/3.156e7/1e6) << " Myr)\n\n";
    
    oss << "Physical Parameters:\n";
    double M_sun = 1.989e30;
    oss << "  Mass M = " << M << " kg (" << (M/M_sun) << " M_sun)\n";
    double ly_to_m = 9.461e15;
    oss << "  Radius r = " << r << " m (" << (r/ly_to_m) << " ly)\n";
    oss << "  Hubble constant H0 = " << H0 << " s^-1\n";
    oss << "  Magnetic field B = " << B << " T (B_crit = " << B_crit << " T)\n";
    oss << "  f_TRZ = " << f_TRZ << "\n";
    oss << "  Expansion E_0 = " << E_0 << ", tau_exp = " << tau_exp << " s\n";
    oss << "  E(t) = " << E_t(t) << "\n";
    oss << "  Wind density rho_wind = " << rho_wind << " kg/m^3, v_wind = " << v_wind << " m/s\n";
    oss << "  Fluid density rho_fluid = " << rho_fluid << " kg/m^3\n";
    oss << "  Gas velocity = " << gas_v << " m/s\n";
    oss << "  DM factor = " << M_DM_factor << "\n\n";
    
    oss << "Computed Acceleration:\n";
    oss << "  g_Bubble(t) = " << compute_g_Bubble(t) << " m/s^2\n\n";
    
    oss << "UQFF Terms:\n";
    double Et = E_t(t);
    double corr_H = 1 + H0 * t;
    double corr_B = 1 - B / B_crit;
    double corr_E = 1 - Et;
    oss << "  Base (with H0, B, E(t)): " << (ug1_base * corr_H * corr_B * corr_E) << " m/s^2\n";
    oss << "  Ug total: " << compute_Ug(Et) << " m/s^2\n";
    oss << "  Lambda: " << ((Lambda * c_light * c_light) / 3.0) << " m/s^2\n";
    
    double cross_vB = gas_v * B;
    double em_base = (q_charge * cross_vB) / proton_mass;
    double corr_UA = 1 + (rho_vac_UA / rho_vac_SCm);
    oss << "  EM (scaled with UA): " << (em_base * corr_UA * scale_EM) << " m/s^2\n";
    
    double sqrt_unc = sqrt(delta_x * delta_p);
    oss << "  Quantum: " << ((hbar / sqrt_unc) * integral_psi * (2 * M_PI / t_Hubble)) << " m/s^2\n";
    
    double V = compute_V();
    oss << "  Fluid: " << ((rho_fluid * V * ug1_base) / M) << " m/s^2\n";
    
    oss << "  Oscillatory: (combined real parts)\n";
    
    double M_dm = M * M_DM_factor;
    double pert1 = delta_rho_over_rho;
    double pert2 = 3 * G * M / (r * r * r);
    double term_dm_force_like = (M + M_dm) * (pert1 + pert2);
    oss << "  DM: " << (term_dm_force_like / M) << " m/s^2\n";
    
    double wind_pressure = rho_wind * v_wind * v_wind;
    oss << "  Stellar Wind Feedback: " << (wind_pressure / rho_fluid) << " m/s^2\n";
    
    oss << "============================================\n";
    return oss.str();
}

bool BubbleNebula::validateConsistency() const {
    bool valid = true;
    
    if (M <= 0 || r <= 0) { std::cerr << "Error: M and r must be positive.\n"; valid = false; }
    if (B < 0 || B_crit <= 0) { std::cerr << "Error: B, B_crit must be non-negative/positive.\n"; valid = false; }
    if (H0 < 0) { std::cerr << "Error: H0 must be non-negative.\n"; valid = false; }
    if (tau_exp <= 0) { std::cerr << "Error: Expansion timescale must be positive.\n"; valid = false; }
    if (E_0 < 0 || E_0 > 1.0) { std::cerr << "Warning: Expansion factor E_0 outside [0,1].\n"; }
    if (rho_fluid <= 0 || rho_wind < 0) { std::cerr << "Error: Fluid/wind densities must be positive/non-negative.\n"; valid = false; }
    if (v_wind < 0 || gas_v < 0) { std::cerr << "Error: Velocities must be non-negative.\n"; valid = false; }
    if (M_DM_factor < 0 || M_DM_factor > 1.0) { std::cerr << "Warning: DM factor outside [0,1].\n"; }
    
    return valid;
}

bool BubbleNebula::autoCorrectAnomalies() {
    bool corrected = false;
    
    double M_sun = 1.989e30;
    double ly_to_m = 9.461e15;
    
    if (M <= 0) { M = 46.0 * M_sun; corrected = true; }
    if (r <= 0) { r = 5.0 * ly_to_m; corrected = true; }
    if (B < 0) { B = 1e-6; corrected = true; }
    if (B_crit <= 0) { B_crit = 1e11; corrected = true; }
    if (H0 < 0) { H0 = 2.184e-18; corrected = true; }
    if (tau_exp <= 0) { tau_exp = 4e6 * 3.156e7; corrected = true; }
    if (E_0 < 0) { E_0 = 0.0; corrected = true; }
    if (E_0 > 1.0) { E_0 = 1.0; corrected = true; }
    if (rho_fluid <= 0) { rho_fluid = 1e-21; corrected = true; }
    if (rho_wind < 0) { rho_wind = 1e-21; corrected = true; }
    if (v_wind < 0) { v_wind = 1.8e6; corrected = true; }
    if (gas_v < 0) { gas_v = 1e5; corrected = true; }
    if (M_DM_factor < 0) { M_DM_factor = 0.1; corrected = true; }
    if (M_DM_factor > 1.0) { M_DM_factor = 1.0; corrected = true; }
    
    if (corrected) updateCache();
    return corrected;
}

// ========== ENHANCED EXAMPLE FUNCTION ==========
void enhanced_bubble_example() {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "=========================================================\n";
    std::cout << "ENHANCED BUBBLE NEBULA DEMONSTRATION\n";
    std::cout << "NGC 7635 Emission Nebula with Expansion\n";
    std::cout << "=========================================================\n\n";
    
    BubbleNebula bubble;
    
    // Step 1: Initial state and validation
    std::cout << "Step 1: Initial State and Validation\n";
    std::cout << "System: " << bubble.getSystemName() << "\n";
    std::cout << "Validation: " << (bubble.validateConsistency() ? "PASS" : "FAIL") << "\n";
    std::cout << "Auto-corrected: " << (bubble.autoCorrectAnomalies() ? "Yes" : "No") << "\n\n";
    
    // Step 2: Time evolution showing expansion E(t)
    std::cout << "Step 2: Time Evolution (Nebula Expansion E(t))\n";
    double t_Myr_array[] = {0.0, 1.0, 2.0, 4.0, 8.0};
    for (double t_Myr : t_Myr_array) {
        double t = t_Myr * 1e6 * 3.156e7;
        double Et = bubble.E_t(t);
        double g = bubble.compute_g_Bubble(t);
        std::cout << "  t = " << t_Myr << " Myr: E(t) = " << Et << ", g = " << g << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 3: Variable listing
    std::cout << "Step 3: Variable Listing\n";
    auto vars = bubble.listVariables();
    std::cout << "Total variables: " << vars.size() << "\n";
    std::cout << "Sample: " << vars[0] << ", " << vars[1] << ", " << vars[11] << " (E_0), " 
              << vars[12] << " (tau_exp)\n\n";
    
    // Step 4: Nebula mass scaling
    std::cout << "Step 4: Nebula Mass Scaling (M sweeps)\n";
    bubble.saveState("original");
    double M_factors[] = {0.5, 1.0, 2.0};
    for (double factor : M_factors) {
        bubble.restoreState("original");
        bubble.expandNebulaScale(factor, 1.0);
        double t = 2e6 * 3.156e7;
        double g = bubble.compute_g_Bubble(t);
        double M_sun = 1.989e30;
        double M = bubble.getVariable("M");
        std::cout << "  M × " << factor << ": M = " << (M/M_sun) << " M_sun, g(2 Myr) = " << g << " m/s^2\n";
    }
    bubble.restoreState("original");
    std::cout << "\n";
    
    // Step 5: Expansion factor scaling (UNIQUE to Bubble Nebula)
    std::cout << "Step 5: Expansion Factor Scaling (E_0 sweeps) - BUBBLE EXPANSION FEATURE\n";
    double E_0_factors[] = {0.5, 1.0, 2.0};
    for (double factor : E_0_factors) {
        bubble.restoreState("original");
        bubble.expandExpansionScale(factor, 1.0);
        double t = 2e6 * 3.156e7;
        double Et = bubble.E_t(t);
        double g = bubble.compute_g_Bubble(t);
        std::cout << "  E_0 × " << factor << ": E(2 Myr) = " << Et << ", g = " << g << " m/s^2\n";
    }
    bubble.restoreState("original");
    std::cout << "\n";
    
    // Step 6: Expansion timescale sweeps
    std::cout << "Step 6: Expansion Timescale Sweeps (tau_exp)\n";
    double tau_exp_factors[] = {0.5, 1.0, 2.0};
    for (double factor : tau_exp_factors) {
        bubble.restoreState("original");
        bubble.expandExpansionScale(1.0, factor);
        double t = 2e6 * 3.156e7;
        double Et = bubble.E_t(t);
        std::cout << "  tau_exp × " << factor << ": E(2 Myr) = " << Et << "\n";
    }
    bubble.restoreState("original");
    std::cout << "\n";
    
    // Step 7: Wind velocity scaling
    std::cout << "Step 7: Wind Velocity Scaling\n";
    double v_wind_factors[] = {0.5, 1.0, 2.0};
    for (double factor : v_wind_factors) {
        bubble.restoreState("original");
        bubble.expandWindMagneticScale(1.0, factor, 1.0);
        double t = 2e6 * 3.156e7;
        double g = bubble.compute_g_Bubble(t);
        std::cout << "  v_wind × " << factor << ": g(2 Myr) = " << g << " m/s^2\n";
    }
    bubble.restoreState("original");
    std::cout << "\n";
    
    // Step 8: Magnetic field scaling
    std::cout << "Step 8: Magnetic Field Scaling\n";
    double B_factors[] = {0.5, 1.0, 2.0};
    for (double factor : B_factors) {
        bubble.restoreState("original");
        bubble.expandWindMagneticScale(1.0, 1.0, factor);
        double t = 2e6 * 3.156e7;
        double g = bubble.compute_g_Bubble(t);
        std::cout << "  B × " << factor << ": g(2 Myr) = " << g << " m/s^2\n";
    }
    bubble.restoreState("original");
    std::cout << "\n";
    
    // Step 9: Parameter space expansion
    std::cout << "Step 9: Parameter Space Expansion (all scalable params)\n";
    bubble.expandParameterSpace(1.2);
    double M_after = bubble.getVariable("M");
    double M_sun = 1.989e30;
    std::cout << "  After 1.2× expansion: M = " << (M_after/M_sun) << " M_sun\n";
    bubble.restoreState("original");
    std::cout << "\n";
    
    // Step 10: Batch operations
    std::cout << "Step 10: Batch Operations (scale multiple variables)\n";
    std::vector<std::string> scale_group = {"M", "r", "B"};
    bubble.scaleVariableGroup(scale_group, 1.1);
    std::cout << "  Scaled {M, r, B} by 1.1×\n";
    bubble.restoreState("original");
    std::cout << "\n";
    
    // Step 11: State management
    std::cout << "Step 11: State Management\n";
    bubble.saveState("state_A");
    bubble.expandExpansionScale(1.5, 1.2);
    bubble.saveState("state_B");
    auto states = bubble.listSavedStates();
    std::cout << "  Saved states: ";
    for (const auto& s : states) std::cout << s << " ";
    std::cout << "\n";
    bubble.restoreState("state_A");
    std::cout << "  Restored state_A\n\n";
    
    // Step 12: Generate parameter variations
    std::cout << "Step 12: Generate Parameter Variations (5% variation)\n";
    auto variations = bubble.generateVariations(3, 5.0);
    std::cout << "  Generated " << variations.size() << " variants with 5% random variation\n";
    std::cout << "  Variant 1 M = " << variations[0]["M"] << " kg\n\n";
    
    // Step 13: Sensitivity analysis
    std::cout << "Step 13: Sensitivity Analysis at 2 Myr\n";
    double t_sens = 2e6 * 3.156e7;
    auto sensitivities = bubble.sensitivityAnalysis(t_sens, 1.0);
    std::cout << "  Top sensitivities (1% perturbation):\n";
    std::vector<std::pair<std::string, double>> sens_vec(sensitivities.begin(), sensitivities.end());
    std::sort(sens_vec.begin(), sens_vec.end(), 
              [](const auto& a, const auto& b) { return a.second > b.second; });
    for (int i = 0; i < 5 && i < (int)sens_vec.size(); ++i) {
        std::cout << "    " << sens_vec[i].first << ": " << sens_vec[i].second << "\n";
    }
    std::cout << "\n";
    
    // Step 14: Auto-refinement with synthetic observations
    std::cout << "Step 14: Auto-Refinement (synthetic observations)\n";
    std::vector<std::pair<double, double>> obs;
    for (int i = 0; i <= 5; ++i) {
        double t_obs = i * 1e6 * 3.156e7;
        double g_obs = bubble.compute_g_Bubble(t_obs) * (1.0 + 0.01 * (rand() % 100 - 50) / 100.0);
        obs.push_back({t_obs, g_obs});
    }
    bubble.autoRefineParameters(obs);
    std::cout << "  Refined parameters based on " << obs.size() << " observations\n\n";
    
    // Step 15: Optimization for maximum acceleration
    std::cout << "Step 15: Optimize for Maximum Acceleration\n";
    bubble.restoreState("original");
    auto metric = [](double g) { return g; };
    double t_opt_start = 0.0;
    double t_opt_end = 8e6 * 3.156e7;
    double best_g = bubble.optimizeForMetric(metric, t_opt_start, t_opt_end, 50);
    std::cout << "  Best g over 8 Myr: " << best_g << " m/s^2\n\n";
    
    // Step 16: Evolutionary system adaptation
    std::cout << "Step 16: Evolutionary System Adaptation (5 generations)\n";
    bubble.restoreState("original");
    auto fitness = [](const BubbleNebula& b) {
        double t = 2e6 * 3.156e7;
        return b.compute_g_Bubble(t);
    };
    bubble.evolveSystem(5, fitness);
    std::cout << "  Evolved system over 5 generations (fitness = g at 2 Myr)\n\n";
    
    // Step 17: Full system report
    std::cout << "Step 17: Full System Report at 2 Myr\n";
    bubble.restoreState("original");
    double t_report = 2e6 * 3.156e7;
    std::string report = bubble.generateReport(t_report);
    std::cout << report << "\n";
    
    // Step 18: Full state export
    std::cout << "Step 18: Full State Export\n";
    std::string exported = bubble.exportState();
    std::cout << "Exported state (first 500 chars):\n";
    std::cout << exported.substr(0, 500) << "...\n\n";
    
    std::cout << "=========================================================\n";
    std::cout << "ENHANCED DEMONSTRATION COMPLETE\n";
    std::cout << "=========================================================\n";
}