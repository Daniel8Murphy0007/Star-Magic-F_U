/**
 * ================================================================================================
 * Header: StarbirthTapestry.h
 *
 * Description: C++ Module for "Tapestry of Blazing Starbirth" (NGC 2014 & NGC 2020) Class
 *              This is the fourth module in a series of 500+ code files for the Universal Quantum
 *              Field Framework (UQFF) simulations, focusing on star-forming region evolution and
 *              gravity equations derived from Hubble datasets, high-energy lab simulations, and
 *              UQFF refinements (dated May 09, 2025, updated for full term inclusion on October 08, 2025).
 *
 * Purpose: Encapsulates the Master Universal Gravity Equation (MUGE) for "Tapestry of Blazing
 *          Starbirth" evolution in the LMC. Includes ALL terms: base gravity with star formation
 *          growth M(t), cosmic expansion (H_0), magnetic correction (static B), UQFF Ug components
 *          with f_TRZ, Lambda, quantum uncertainty, scaled EM with [UA], fluid dynamics, oscillatory
 *          waves, DM/density perturbations, and stellar wind feedback (pressure / density for acc).
 *          Supports dynamic variable updates for all parameters.
 *
 * Integration: Designed for inclusion in base program 'ziqn233h.cpp' (not present here).
 *              Instantiate class in main: StarbirthTapestry tapestry;
 *              Compute: double g = tapestry.compute_g_Starbirth(t);
 *
 * Key Features:
 *   - Default values from UQFF document: M_initial = 240 Msun, r = 10 ly, B = 1e-6 T, etc.
 *   - M(t) = M_initial * (1 + M_dot_factor * exp(-t / tau_SF)), with M_dot_factor = 10000 / 240.
 *   - Units handled: ly to m, Msun to kg; wind term as (rho * v_wind^2) / rho_fluid for acceleration.
 *   - Setter methods for updates: setVar(double new_val) or addToVar(double delta)/subtractFromVar(double delta).
 *   - Computes g_Starbirth(r, t) with every term explicitly included.
 *
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: October 08, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

#ifndef STARBIRTH_TAPESTRY_H
#define STARBIRTH_TAPESTRY_H

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

class StarbirthTapestry {
private:
    // Core parameters (mutable for updates)
    double G;               // Gravitational constant
    double M_initial;       // Initial mass (kg)
    double r;               // Radius (m)
    double H0;              // Hubble constant (s^-1)
    double B;               // Static magnetic field (T)
    double B_crit;          // Critical B field (T)
    double Lambda;          // Cosmological constant
    double c_light;         // Speed of light
    double q_charge;        // Charge (proton)
    double gas_v;           // Gas velocity for EM (m/s)
    double f_TRZ;           // Time-reversal factor
    double M_dot_factor;    // Star formation factor (dimensionless)
    double tau_SF;          // Star formation timescale (s)
    double rho_wind;        // Wind density (kg/m^3)
    double v_wind;          // Wind velocity (m/s)
    double rho_fluid;       // Fluid density (kg/m^3)
    double rho_vac_UA;      // UA vacuum density (J/m^3)
    double rho_vac_SCm;     // SCm vacuum density (J/m^3)
    double scale_EM;        // EM scaling factor

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
    double proton_mass;     // Proton mass for EM acceleration

    // Computed caches (updated on demand)
    double ug1_base;        // Cached Ug1 for initial M

public:
    // Constructor with default UQFF values
    StarbirthTapestry() {
        initializeDefaults();
    }

    // Destructor (empty)
    ~StarbirthTapestry() {}

    // Initialization method (called in constructor)
    void initializeDefaults() {
        G = 6.6743e-11;
        double M_sun = 1.989e30;
        double M_initial_sun = 240.0;
        M_initial = M_initial_sun * M_sun;
        double ly_to_m = 9.461e15;
        r = 10.0 * ly_to_m;
        H0 = 2.184e-18;
        B = 1e-6;
        B_crit = 1e11;
        Lambda = 1.1e-52;
        c_light = 3e8;
        q_charge = 1.602e-19;
        gas_v = 1e5;
        f_TRZ = 0.1;
        double gas_mass_sun = 10000.0;
        M_dot_factor = gas_mass_sun / M_initial_sun;
        tau_SF = 5e6 * 3.156e7;
        rho_wind = 1e-21;
        v_wind = 2e6;
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
        A_osc = 1e-10;  // Small for nebula scale
        k_osc = 1.0 / r;
        omega_osc = 2 * M_PI / (r / c_light);
        x_pos = r;
        M_DM_factor = 0.1;
        delta_rho_over_rho = 1e-5;

        updateCache();
    }

    // Cache update for efficiency (call after parameter changes)
    void updateCache() {
        ug1_base = (G * M_initial) / (r * r);
    }

    // Universal setter for any variable (by name, for flexibility)
    bool setVariable(const std::string& varName, double newValue) {
        if (varName == "G") { G = newValue; }
        else if (varName == "M_initial") { M_initial = newValue; }
        else if (varName == "r") { r = newValue; }
        else if (varName == "H0") { H0 = newValue; }
        else if (varName == "B") { B = newValue; }
        else if (varName == "B_crit") { B_crit = newValue; }
        else if (varName == "Lambda") { Lambda = newValue; }
        else if (varName == "c_light") { c_light = newValue; }
        else if (varName == "q_charge") { q_charge = newValue; }
        else if (varName == "gas_v") { gas_v = newValue; }
        else if (varName == "f_TRZ") { f_TRZ = newValue; }
        else if (varName == "M_dot_factor") { M_dot_factor = newValue; }
        else if (varName == "tau_SF") { tau_SF = newValue; }
        else if (varName == "rho_wind") { rho_wind = newValue; }
        else if (varName == "v_wind") { v_wind = newValue; }
        else if (varName == "rho_fluid") { rho_fluid = newValue; }
        else if (varName == "rho_vac_UA") { rho_vac_UA = newValue; }
        else if (varName == "rho_vac_SCm") { rho_vac_SCm = newValue; }
        else if (varName == "scale_EM") { scale_EM = newValue; }
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
        else if (varName == "proton_mass") { proton_mass = newValue; }
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
        else if (varName == "M_initial") return M_initial;
        else if (varName == "r") return r;
        else if (varName == "H0") return H0;
        else if (varName == "B") return B;
        else if (varName == "B_crit") return B_crit;
        else if (varName == "Lambda") return Lambda;
        else if (varName == "c_light") return c_light;
        else if (varName == "q_charge") return q_charge;
        else if (varName == "gas_v") return gas_v;
        else if (varName == "f_TRZ") return f_TRZ;
        else if (varName == "M_dot_factor") return M_dot_factor;
        else if (varName == "tau_SF") return tau_SF;
        else if (varName == "rho_wind") return rho_wind;
        else if (varName == "v_wind") return v_wind;
        else if (varName == "rho_fluid") return rho_fluid;
        else if (varName == "rho_vac_UA") return rho_vac_UA;
        else if (varName == "rho_vac_SCm") return rho_vac_SCm;
        else if (varName == "scale_EM") return scale_EM;
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
        else if (varName == "proton_mass") return proton_mass;
        else {
            std::cerr << "Error: Unknown variable '" << varName << "'." << std::endl;
            return 0.0;
        }
    }

    // M(t) computation
    double M_t(double t) const {
        double M_dot = M_dot_factor * exp(-t / tau_SF);
        return M_initial * (1 + M_dot);
    }

    // Ug terms computation
    double compute_Ug(double Mt) const {
        double Ug1 = (G * Mt) / (r * r);
        double Ug2 = 0.0;
        double Ug3 = 0.0;
        double corr_B = 1 - B / B_crit;
        double Ug4 = Ug1 * corr_B;
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1 + f_TRZ);
    }

    // Volume computation for fluid
    double compute_V() const {
        return (4.0 / 3.0) * M_PI * r * r * r;
    }

    // Main MUGE computation (includes ALL terms)
    double compute_g_Starbirth(double t) const {
        if (t < 0) {
            std::cerr << "Error: Time t must be non-negative." << std::endl;
            return 0.0;
        }

        double Mt = M_t(t);
        double ug1_t = (G * Mt) / (r * r);

        // Term 1: Base + H0 + B corrections
        double corr_H = 1 + H0 * t;
        double corr_B = 1 - B / B_crit;
        double term1 = ug1_t * corr_H * corr_B;

        // Term 2: UQFF Ug with f_TRZ
        double term2 = compute_Ug(Mt);

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
        double term_fluid = (rho_fluid * V * ug1_t) / Mt;

        // Oscillatory terms (real parts)
        double term_osc1 = 2 * A_osc * cos(k_osc * x_pos) * cos(omega_osc * t);
        double arg = k_osc * x_pos - omega_osc * t;
        double term_osc2 = (2 * M_PI / t_Hubble_gyr) * A_osc * cos(arg);
        double term_osc = term_osc1 + term_osc2;

        // DM and density perturbation term (converted to acceleration)
        double M_dm = Mt * M_DM_factor;
        double pert1 = delta_rho_over_rho;
        double pert2 = 3 * G * Mt / (r * r * r);
        double term_dm_force_like = (Mt + M_dm) * (pert1 + pert2);
        double term_DM = term_dm_force_like / Mt;

        // Stellar wind feedback term (pressure / density for acceleration)
        double wind_pressure = rho_wind * v_wind * v_wind;
        double term_wind = wind_pressure / rho_fluid;

        // Total g_Starbirth (all terms summed)
        return term1 + term2 + term3 + term4 + term_q + term_fluid + term_osc + term_DM + term_wind;
    }

    // Debug/Output method (for transparency in base program)
    void printParameters(std::ostream& os = std::cout) const {
        os << std::fixed << std::setprecision(3);
        os << "Tapestry of Blazing Starbirth Parameters:" << std::endl;
        os << "G: " << G << ", M_initial: " << M_initial << ", r: " << r << std::endl;
        os << "H0: " << H0 << ", B: " << B << ", B_crit: " << B_crit << std::endl;
        os << "f_TRZ: " << f_TRZ << ", M_dot_factor: " << M_dot_factor << ", tau_SF: " << tau_SF << std::endl;
        os << "rho_fluid: " << rho_fluid << ", rho_wind: " << rho_wind << ", v_wind: " << v_wind << std::endl;
        os << "gas_v: " << gas_v << ", M_DM_factor: " << M_DM_factor << std::endl;
        os << "A_osc: " << A_osc << ", delta_rho_over_rho: " << delta_rho_over_rho << std::endl;
        os << "ug1_base: " << ug1_base << std::endl;
    }

    // Example computation at t=2.5 Myr (for testing)
    double exampleAt2_5Myr() const {
        double t_example = 2.5e6 * 3.156e7;
        return compute_g_Starbirth(t_example);
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
    void expandStarFormationScale(double M_dot_factor_scale, double tau_SF_scale);
    void expandWindFeedbackScale(double rho_wind_scale, double v_wind_scale);
    void expandMagneticFluidScale(double B_scale, double rho_fluid_scale);

    // --- Self-Refinement (3 methods) ---
    void autoRefineParameters(const std::vector<std::pair<double, double>>& observations);
    void calibrateToObservations(const std::vector<double>& times, const std::vector<double>& g_obs);
    double optimizeForMetric(std::function<double(double)> metric, double t_start, double t_end, int steps);

    // --- Parameter Exploration (1 method) ---
    std::vector<std::map<std::string, double>> generateVariations(int count, double variation_pct);

    // --- Adaptive Evolution (2 methods) ---
    void mutateParameters(double mutation_rate, std::mt19937& rng);
    void evolveSystem(int generations, std::function<double(const StarbirthTapestry&)> fitness);

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

#endif // STARBIRTH_TAPESTRY_H

// ========== IMPLEMENTATION OF ENHANCED METHODS (Outside class) ==========

// Anonymous namespace for state storage (external to class since class uses member variables)
namespace {
    std::map<std::string, std::map<std::string, double>> starbirth_saved_states;
}

// --- Variable Management (5 methods) ---
bool StarbirthTapestry::createVariable(const std::string& name, double value) {
    // For starbirth class with member variables, we sync to member if it exists
    return setVariable(name, value);
}

bool StarbirthTapestry::removeVariable(const std::string& name) {
    // Cannot remove built-in member variables, return false
    std::cerr << "Warning: Cannot remove built-in variable '" << name << "' in starbirth class." << std::endl;
    return false;
}

bool StarbirthTapestry::cloneVariable(const std::string& src, const std::string& dest) {
    double val = getVariable(src);
    return setVariable(dest, val);
}

std::vector<std::string> StarbirthTapestry::listVariables() const {
    return {"G", "M_initial", "r", "H0", "B", "B_crit", "Lambda", "c_light", "q_charge",
            "gas_v", "f_TRZ", "M_dot_factor", "tau_SF", "rho_wind", "v_wind", "rho_fluid",
            "rho_vac_UA", "rho_vac_SCm", "scale_EM", "hbar", "t_Hubble", "t_Hubble_gyr",
            "delta_x", "delta_p", "integral_psi", "A_osc", "k_osc", "omega_osc", "x_pos",
            "M_DM_factor", "delta_rho_over_rho", "proton_mass"};
}

std::string StarbirthTapestry::getSystemName() const {
    return "StarbirthTapestry";
}

// --- Batch Operations (2 methods) ---
bool StarbirthTapestry::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        double val = getVariable(name);
        if (!setVariable(name, func(val))) return false;
    }
    return true;
}

bool StarbirthTapestry::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    return transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// --- Self-Expansion (4 methods) ---
void StarbirthTapestry::expandParameterSpace(double factor) {
    std::vector<std::string> expandable = {"M_initial", "r", "B", "rho_fluid", "rho_wind", "A_osc", "M_DM_factor"};
    scaleVariableGroup(expandable, factor);
}

void StarbirthTapestry::expandStarFormationScale(double M_dot_factor_scale, double tau_SF_scale) {
    setVariable("M_dot_factor", getVariable("M_dot_factor") * M_dot_factor_scale);
    setVariable("tau_SF", getVariable("tau_SF") * tau_SF_scale);
}

void StarbirthTapestry::expandWindFeedbackScale(double rho_wind_scale, double v_wind_scale) {
    setVariable("rho_wind", getVariable("rho_wind") * rho_wind_scale);
    setVariable("v_wind", getVariable("v_wind") * v_wind_scale);
}

void StarbirthTapestry::expandMagneticFluidScale(double B_scale, double rho_fluid_scale) {
    setVariable("B", getVariable("B") * B_scale);
    setVariable("B_crit", getVariable("B_crit") * B_scale);
    setVariable("rho_fluid", getVariable("rho_fluid") * rho_fluid_scale);
}

// --- Self-Refinement (3 methods) ---
void StarbirthTapestry::autoRefineParameters(const std::vector<std::pair<double, double>>& observations) {
    if (observations.empty()) return;
    
    double sum_error = 0.0;
    for (const auto& obs : observations) {
        double t = obs.first;
        double g_obs = obs.second;
        double g_calc = compute_g_Starbirth(t);
        sum_error += std::abs(g_calc - g_obs);
    }
    double avg_error = sum_error / observations.size();
    
    // Simple refinement: adjust M_dot_factor and tau_SF based on error
    if (avg_error > 1e-6) {
        double adj_factor = 1.0 - std::min(0.1, avg_error / 1e6);
        setVariable("M_dot_factor", getVariable("M_dot_factor") * adj_factor);
        setVariable("tau_SF", getVariable("tau_SF") * (2.0 - adj_factor));
    }
}

void StarbirthTapestry::calibrateToObservations(const std::vector<double>& times, const std::vector<double>& g_obs) {
    if (times.size() != g_obs.size() || times.empty()) return;
    
    std::vector<std::pair<double, double>> obs;
    for (size_t i = 0; i < times.size(); ++i) {
        obs.push_back({times[i], g_obs[i]});
    }
    
    // Iterative refinement (5 passes)
    for (int iter = 0; iter < 5; ++iter) {
        autoRefineParameters(obs);
    }
}

double StarbirthTapestry::optimizeForMetric(std::function<double(double)> metric, double t_start, double t_end, int steps) {
    double best_score = -1e100;
    double dt = (t_end - t_start) / steps;
    
    for (int i = 0; i <= steps; ++i) {
        double t = t_start + i * dt;
        double g = compute_g_Starbirth(t);
        double score = metric(g);
        if (score > best_score) best_score = score;
    }
    return best_score;
}

// --- Parameter Exploration (1 method) ---
std::vector<std::map<std::string, double>> StarbirthTapestry::generateVariations(int count, double variation_pct) {
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
void StarbirthTapestry::mutateParameters(double mutation_rate, std::mt19937& rng) {
    std::uniform_real_distribution<> dis(-mutation_rate, mutation_rate);
    auto vars = listVariables();
    
    for (const auto& v : vars) {
        if (v == "c_light" || v == "G" || v == "hbar") continue; // Skip constants
        double val = getVariable(v);
        double delta = val * dis(rng);
        setVariable(v, val + delta);
    }
}

void StarbirthTapestry::evolveSystem(int generations, std::function<double(const StarbirthTapestry&)> fitness) {
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
bool StarbirthTapestry::saveState(const std::string& stateName) {
    std::map<std::string, double> state;
    auto vars = listVariables();
    for (const auto& v : vars) {
        state[v] = getVariable(v);
    }
    starbirth_saved_states[stateName] = state;
    return true;
}

bool StarbirthTapestry::restoreState(const std::string& stateName) {
    auto it = starbirth_saved_states.find(stateName);
    if (it == starbirth_saved_states.end()) return false;
    
    for (const auto& pair : it->second) {
        setVariable(pair.first, pair.second);
    }
    return true;
}

std::vector<std::string> StarbirthTapestry::listSavedStates() const {
    std::vector<std::string> names;
    for (const auto& pair : starbirth_saved_states) {
        names.push_back(pair.first);
    }
    return names;
}

std::string StarbirthTapestry::exportState() const {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(6);
    oss << "StarbirthTapestry State Export:\n";
    auto vars = listVariables();
    for (const auto& v : vars) {
        oss << v << " = " << getVariable(v) << "\n";
    }
    return oss.str();
}

// --- System Analysis (4 methods) ---
std::map<std::string, double> StarbirthTapestry::sensitivityAnalysis(double t, double delta_pct) {
    std::map<std::string, double> sensitivities;
    double g_base = compute_g_Starbirth(t);
    
    auto vars = listVariables();
    for (const auto& v : vars) {
        if (v == "c_light" || v == "G" || v == "hbar") continue; // Skip constants
        
        double original = getVariable(v);
        double delta = original * delta_pct / 100.0;
        
        setVariable(v, original + delta);
        double g_plus = compute_g_Starbirth(t);
        setVariable(v, original);
        
        double sensitivity = (g_base != 0.0) ? std::abs((g_plus - g_base) / g_base) : 0.0;
        sensitivities[v] = sensitivity;
    }
    
    return sensitivities;
}

std::string StarbirthTapestry::generateReport(double t) const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6);
    oss << "============================================\n";
    oss << "TAPESTRY OF BLAZING STARBIRTH REPORT\n";
    oss << "NGC 2014 & NGC 2020 (LMC)\n";
    oss << "============================================\n";
    oss << "Time: t = " << t << " s (" << (t/3.156e7/1e6) << " Myr)\n\n";
    
    oss << "Physical Parameters:\n";
    double M_sun = 1.989e30;
    oss << "  Initial Mass M_initial = " << M_initial << " kg (" << (M_initial/M_sun) << " M_sun)\n";
    oss << "  M(t) = " << M_t(t) << " kg (" << (M_t(t)/M_sun) << " M_sun)\n";
    double ly_to_m = 9.461e15;
    oss << "  Radius r = " << r << " m (" << (r/ly_to_m) << " ly)\n";
    oss << "  Magnetic field B = " << B << " T (B_crit = " << B_crit << " T)\n";
    oss << "  Star formation M_dot_factor = " << M_dot_factor << ", tau_SF = " << tau_SF << " s\n";
    oss << "  Wind density rho_wind = " << rho_wind << " kg/m^3, v_wind = " << v_wind << " m/s\n";
    oss << "  Fluid density rho_fluid = " << rho_fluid << " kg/m^3\n";
    oss << "  Gas velocity = " << gas_v << " m/s\n";
    oss << "  DM factor = " << M_DM_factor << "\n\n";
    
    oss << "Computed Acceleration:\n";
    oss << "  g_Starbirth(t) = " << compute_g_Starbirth(t) << " m/s^2\n\n";
    
    oss << "UQFF Terms:\n";
    double Mt = M_t(t);
    double ug1_t = (G * Mt) / (r * r);
    double corr_H = 1 + H0 * t;
    double corr_B = 1 - B / B_crit;
    oss << "  Base (with H0, B, M(t)): " << (ug1_t * corr_H * corr_B) << " m/s^2\n";
    oss << "  Ug total: " << compute_Ug(Mt) << " m/s^2\n";
    oss << "  Lambda: " << ((Lambda * c_light * c_light) / 3.0) << " m/s^2\n";
    
    double cross_vB = gas_v * B;
    double em_base = (q_charge * cross_vB) / proton_mass;
    double corr_UA = 1 + (rho_vac_UA / rho_vac_SCm);
    oss << "  EM (scaled with UA): " << (em_base * corr_UA * scale_EM) << " m/s^2\n";
    
    double sqrt_unc = sqrt(delta_x * delta_p);
    oss << "  Quantum: " << ((hbar / sqrt_unc) * integral_psi * (2 * M_PI / t_Hubble)) << " m/s^2\n";
    
    double V = compute_V();
    oss << "  Fluid: " << ((rho_fluid * V * ug1_t) / Mt) << " m/s^2\n";
    
    oss << "  Oscillatory: (combined real parts)\n";
    
    double M_dm = Mt * M_DM_factor;
    double pert1 = delta_rho_over_rho;
    double pert2 = 3 * G * Mt / (r * r * r);
    double term_dm_force_like = (Mt + M_dm) * (pert1 + pert2);
    oss << "  DM: " << (term_dm_force_like / Mt) << " m/s^2\n";
    
    double wind_pressure = rho_wind * v_wind * v_wind;
    oss << "  Stellar Wind Feedback: " << (wind_pressure / rho_fluid) << " m/s^2\n";
    
    oss << "============================================\n";
    return oss.str();
}

bool StarbirthTapestry::validateConsistency() const {
    bool valid = true;
    
    if (M_initial <= 0 || r <= 0) { std::cerr << "Error: M_initial and r must be positive.\n"; valid = false; }
    if (B < 0 || B_crit <= 0) { std::cerr << "Error: B, B_crit must be non-negative/positive.\n"; valid = false; }
    if (tau_SF <= 0) { std::cerr << "Error: Star formation timescale must be positive.\n"; valid = false; }
    if (rho_fluid <= 0 || rho_wind < 0) { std::cerr << "Error: Fluid/wind densities must be positive/non-negative.\n"; valid = false; }
    if (v_wind < 0 || gas_v < 0) { std::cerr << "Error: Velocities must be non-negative.\n"; valid = false; }
    if (M_DM_factor < 0 || M_DM_factor > 1.0) { std::cerr << "Warning: DM factor outside [0,1].\n"; }
    
    return valid;
}

bool StarbirthTapestry::autoCorrectAnomalies() {
    bool corrected = false;
    
    double M_sun = 1.989e30;
    double ly_to_m = 9.461e15;
    
    if (M_initial <= 0) { M_initial = 240.0 * M_sun; corrected = true; }
    if (r <= 0) { r = 10.0 * ly_to_m; corrected = true; }
    if (B < 0) { B = 1e-6; corrected = true; }
    if (B_crit <= 0) { B_crit = 1e11; corrected = true; }
    if (tau_SF <= 0) { tau_SF = 5e6 * 3.156e7; corrected = true; }
    if (rho_fluid <= 0) { rho_fluid = 1e-21; corrected = true; }
    if (rho_wind < 0) { rho_wind = 1e-21; corrected = true; }
    if (v_wind < 0) { v_wind = 2e6; corrected = true; }
    if (gas_v < 0) { gas_v = 1e5; corrected = true; }
    if (M_DM_factor < 0) { M_DM_factor = 0.1; corrected = true; }
    if (M_DM_factor > 1.0) { M_DM_factor = 1.0; corrected = true; }
    
    if (corrected) updateCache();
    return corrected;
}

// ========== ENHANCED EXAMPLE FUNCTION ==========
void enhancedStarbirthTapestryExample() {
    std::cout << "\n========== ENHANCED STARBIRTH TAPESTRY UQFF EXAMPLE ==========\n\n";
    
    StarbirthTapestry tapestry;
    
    // Step 1: Initial state
    std::cout << "Step 1: Initial Configuration\n";
    tapestry.printParameters();
    double t0 = 0.0;
    std::cout << "g_Starbirth(t=0) = " << tapestry.compute_g_Starbirth(t0) << " m/s^2\n\n";
    
    // Step 2: Time evolution (0, 1, 2.5, 5, 10 Myr)
    std::cout << "Step 2: Time Evolution (0, 1, 2.5, 5, 10 Myr)\n";
    double M_sun = 1.989e30;
    for (double t_myr : {0.0, 1.0, 2.5, 5.0, 10.0}) {
        double t = t_myr * 1e6 * 3.156e7;
        std::cout << "  t = " << t_myr << " Myr: g = " << tapestry.compute_g_Starbirth(t) 
                  << " m/s^2, M(t) = " << (tapestry.M_t(t)/M_sun) << " M_sun\n";
    }
    std::cout << "\n";
    
    // Step 3: Star formation scaling
    std::cout << "Step 3: Star Formation Scaling (M_dot_factor x1.3, tau_SF x0.8)\n";
    tapestry.expandStarFormationScale(1.3, 0.8);
    double t_test = 2.5e6 * 3.156e7;
    std::cout << "After expansion: M_dot_factor = " << tapestry.getVariable("M_dot_factor") 
              << ", tau_SF = " << tapestry.getVariable("tau_SF") << " s\n";
    std::cout << "g_Starbirth(t=2.5 Myr) = " << tapestry.compute_g_Starbirth(t_test) << " m/s^2\n\n";
    
    // Step 4: Wind feedback scaling
    std::cout << "Step 4: Wind Feedback Scaling (rho_wind x1.5, v_wind x1.2)\n";
    tapestry.expandWindFeedbackScale(1.5, 1.2);
    std::cout << "After expansion: rho_wind = " << tapestry.getVariable("rho_wind") 
              << " kg/m^3, v_wind = " << tapestry.getVariable("v_wind") << " m/s\n";
    std::cout << "g_Starbirth(t=2.5 Myr) = " << tapestry.compute_g_Starbirth(t_test) << " m/s^2\n\n";
    
    // Step 5: Magnetic & fluid scaling
    std::cout << "Step 5: Magnetic & Fluid Scaling (B x1.4, rho_fluid x1.3)\n";
    tapestry.expandMagneticFluidScale(1.4, 1.3);
    std::cout << "After expansion: B = " << tapestry.getVariable("B") 
              << " T, rho_fluid = " << tapestry.getVariable("rho_fluid") << " kg/m^3\n";
    std::cout << "g_Starbirth(t=2.5 Myr) = " << tapestry.compute_g_Starbirth(t_test) << " m/s^2\n\n";
    
    // Step 6: State save/restore
    std::cout << "Step 6: State Management\n";
    tapestry.saveState("expanded_state");
    tapestry.setVariable("M_initial", 500.0 * M_sun);
    std::cout << "Modified M_initial to 500 M_sun: g = " << tapestry.compute_g_Starbirth(t_test) << " m/s^2\n";
    tapestry.restoreState("expanded_state");
    std::cout << "Restored state: M_initial = " << (tapestry.getVariable("M_initial")/M_sun) 
              << " M_sun, g = " << tapestry.compute_g_Starbirth(t_test) << " m/s^2\n";
    std::cout << "Saved states: ";
    for (const auto& s : tapestry.listSavedStates()) std::cout << s << " ";
    std::cout << "\n\n";
    
    // Step 7: Sensitivity analysis
    std::cout << "Step 7: Sensitivity Analysis at t=2.5 Myr (top 5 parameters)\n";
    auto sens = tapestry.sensitivityAnalysis(t_test, 1.0);
    std::vector<std::pair<std::string, double>> sens_vec(sens.begin(), sens.end());
    std::sort(sens_vec.begin(), sens_vec.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
    for (int i = 0; i < std::min(5, (int)sens_vec.size()); ++i) {
        std::cout << "  " << sens_vec[i].first << ": " << sens_vec[i].second << "\n";
    }
    std::cout << "\n";
    
    // Step 8: Generate variations
    std::cout << "Step 8: Generate Parameter Variations (3 variants, 10% variation)\n";
    auto variations = tapestry.generateVariations(3, 10.0);
    for (size_t i = 0; i < variations.size(); ++i) {
        std::cout << "  Variant " << (i+1) << ": M_initial = " << (variations[i]["M_initial"]/M_sun) 
                  << " M_sun, v_wind = " << variations[i]["v_wind"] << " m/s\n";
    }
    std::cout << "\n";
    
    // Step 9: Batch transformation
    std::cout << "Step 9: Batch Transform (scale density parameters by 1.2)\n";
    tapestry.transformVariableGroup({"rho_fluid", "rho_wind"}, [](double v) { return v * 1.2; });
    std::cout << "After transform: rho_fluid = " << tapestry.getVariable("rho_fluid") 
              << " kg/m^3, rho_wind = " << tapestry.getVariable("rho_wind") << " kg/m^3\n";
    std::cout << "g_Starbirth(t=2.5 Myr) = " << tapestry.compute_g_Starbirth(t_test) << " m/s^2\n\n";
    
    // Step 10: Consistency validation
    std::cout << "Step 10: Consistency Validation\n";
    bool valid = tapestry.validateConsistency();
    std::cout << "System is " << (valid ? "VALID" : "INVALID") << "\n\n";
    
    // Step 11: Metric optimization
    std::cout << "Step 11: Optimize for Maximum g (t=0 to 10 Myr, 100 steps)\n";
    double max_g = tapestry.optimizeForMetric([](double g) { return g; }, 0.0, 10e6 * 3.156e7, 100);
    std::cout << "Maximum g found: " << max_g << " m/s^2\n\n";
    
    // Step 12: Full system report
    std::cout << "Step 12: Full System Report at t=3 Myr\n";
    double t_report = 3e6 * 3.156e7;
    std::cout << tapestry.generateReport(t_report) << "\n";
    
    // Step 13: Initial mass sweep
    std::cout << "Step 13: Initial Mass Sweep (M_initial = 100, 240, 400, 600 M_sun)\n";
    tapestry.saveState("before_sweep");
    for (double M_val : {100.0, 240.0, 400.0, 600.0}) {
        tapestry.setVariable("M_initial", M_val * M_sun);
        std::cout << "  M_initial = " << M_val << " M_sun: g(t=2.5Myr) = " << tapestry.compute_g_Starbirth(t_test) << " m/s^2\n";
    }
    tapestry.restoreState("before_sweep");
    std::cout << "\n";
    
    // Step 14: Star formation rate sweep
    std::cout << "Step 14: Star Formation Factor Sweep (M_dot_factor = 20, 41.7, 60, 80)\n";
    for (double M_dot : {20.0, 41.7, 60.0, 80.0}) {
        tapestry.setVariable("M_dot_factor", M_dot);
        std::cout << "  M_dot_factor = " << M_dot << ": M(t=2.5Myr) = " << (tapestry.M_t(t_test)/M_sun) 
                  << " M_sun, g = " << tapestry.compute_g_Starbirth(t_test) << " m/s^2\n";
    }
    tapestry.restoreState("before_sweep");
    std::cout << "\n";
    
    // Step 15: Wind velocity sweep
    std::cout << "Step 15: Wind Velocity Sweep (v_wind = 1, 2, 3, 4 x 10^6 m/s)\n";
    for (double v_wind : {1e6, 2e6, 3e6, 4e6}) {
        tapestry.setVariable("v_wind", v_wind);
        std::cout << "  v_wind = " << (v_wind/1e6) << " x 10^6 m/s: g(t=2.5Myr) = " << tapestry.compute_g_Starbirth(t_test) << " m/s^2\n";
    }
    tapestry.restoreState("before_sweep");
    std::cout << "\n";
    
    // Step 16: Magnetic field sweep
    std::cout << "Step 16: Magnetic Field Sweep (B = 0.5, 1.0, 1.5, 2.0 x 10^-6 T)\n";
    for (double B_val : {0.5e-6, 1.0e-6, 1.5e-6, 2.0e-6}) {
        tapestry.setVariable("B", B_val);
        std::cout << "  B = " << (B_val*1e6) << " x 10^-6 T: g(t=2.5Myr) = " << tapestry.compute_g_Starbirth(t_test) << " m/s^2\n";
    }
    tapestry.restoreState("before_sweep");
    std::cout << "\n";
    
    // Step 17: DM factor sweep
    std::cout << "Step 17: Dark Matter Factor Sweep (M_DM_factor = 0.0, 0.1, 0.2, 0.3)\n";
    for (double dm : {0.0, 0.1, 0.2, 0.3}) {
        tapestry.setVariable("M_DM_factor", dm);
        std::cout << "  M_DM_factor = " << dm << ": g(t=2.5Myr) = " << tapestry.compute_g_Starbirth(t_test) << " m/s^2\n";
    }
    tapestry.restoreState("before_sweep");
    std::cout << "\n";
    
    // Step 18: Export final state
    std::cout << "Step 18: Export Final State\n";
    std::cout << tapestry.exportState() << "\n";
    
    std::cout << "========== ENHANCED EXAMPLE COMPLETE ==========\n\n";
}