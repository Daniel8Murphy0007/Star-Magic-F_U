/**
 * ================================================================================================
 * Header: NGC3603.h
 *
 * Description: C++ Module for NGC 3603 Extreme Star Cluster Class
 *              This is the eleventh module in a series of 500+ code files for the Universal Quantum
 *              Field Framework (UQFF) simulations, focusing on young massive star cluster evolution
 *              and gravity equations derived from Hubble datasets, high-energy lab simulations, and
 *              UQFF refinements (dated May 09, 2025, updated for full term inclusion on October 08, 2025).
 *
 * Purpose: Encapsulates the Master Universal Gravity Equation (MUGE) for NGC 3603 evolution.
 *          Includes ALL terms: base gravity with mass growth M(t), cosmic expansion (H_0), magnetic
 *          correction (static B), UQFF Ug components with f_TRZ, Lambda, quantum uncertainty, scaled
 *          EM with [UA], fluid dynamics, oscillatory waves, DM/density perturbations, stellar wind
 *          feedback (pressure / density for acc), and cavity pressure P(t) / rho_fluid.
 *          Supports dynamic variable updates for all parameters.
 *
 * Integration: Designed for inclusion in base program 'ziqn233h.cpp' (not present here).
 *              Instantiate class in main: NGC3603 ngc3603;
 *              Compute: double g = ngc3603.compute_g_NGC3603(t);
 *
 * Key Features:
 *   - Default values from UQFF document: M0 = 400,000 Msun, r = 8.998e15 m (9.5 ly), tau_SF = 1 Myr,
 *     rho_wind = 1e-20 kg/m^3, v_wind = 2e6 m/s, P0 = 4e-8 (Pa), tau_exp = 1 Myr, B = 1e-5 T.
 *   - Units handled: Msun to kg, ly to m; pressure terms as P / rho_fluid for acceleration.
 *   - Setter methods for updates: setVar(double new_val) or addToVar(double delta)/subtractFromVar(double delta).
 *   - Computes g_NGC3603(r, t) with every term explicitly included.
 *
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: October 08, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

#ifndef NGC_3603_H
#define NGC_3603_H

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

class NGC3603 {
private:
    // Core parameters (mutable for updates)
    double G;               // Gravitational constant
    double M0;              // Initial mass (kg)
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
    double P0;              // Initial pressure (Pa)
    double tau_exp;         // Expansion timescale (s)
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
    double ug1_base;        // Cached Ug1 for initial M0

public:
    // Constructor with default UQFF values
    NGC3603() {
        initializeDefaults();
    }

    // Destructor (empty)
    ~NGC3603() {}

    // Initialization method (called in constructor)
    void initializeDefaults() {
        G = 6.6743e-11;
        double M_sun = 1.989e30;
        M0 = 400000.0 * M_sun;
        double ly_to_m = 9.461e15;
        r = 9.5 * ly_to_m;
        H0 = 2.184e-18;
        B = 1e-5;
        B_crit = 1e11;
        Lambda = 1.1e-52;
        c_light = 3e8;
        q_charge = 1.602e-19;
        gas_v = 1e5;
        f_TRZ = 0.1;
        M_dot_factor = 1.0;
        tau_SF = 1e6 * 3.156e7;
        rho_wind = 1e-20;
        v_wind = 2e6;
        rho_fluid = 1e-20;
        P0 = 4e-8;
        tau_exp = 1e6 * 3.156e7;
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
        ug1_base = (G * M0) / (r * r);
    }

    // Universal setter for any variable (by name, for flexibility)
    bool setVariable(const std::string& varName, double newValue) {
        if (varName == "G") { G = newValue; }
        else if (varName == "M0") { M0 = newValue; }
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
        else if (varName == "P0") { P0 = newValue; }
        else if (varName == "tau_exp") { tau_exp = newValue; }
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
        else if (varName == "M0") return M0;
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
        else if (varName == "P0") return P0;
        else if (varName == "tau_exp") return tau_exp;
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

    // M(t) computation
    double M_t(double t) const {
        double M_dot = M_dot_factor * exp(-t / tau_SF);
        return M0 * (1 + M_dot);
    }

    // P(t) computation
    double P_t(double t) const {
        return P0 * exp(-t / tau_exp);
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
    double compute_g_NGC3603(double t) const {
        if (t < 0) {
            std::cerr << "Error: Time t must be non-negative." << std::endl;
            return 0.0;
        }

        double Mt = M_t(t);
        double Pt = P_t(t);
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

        // Cavity pressure term (P(t) / rho_fluid for acceleration)
        double term_pressure = Pt / rho_fluid;

        // Total g_NGC3603 (all terms summed)
        return term1 + term2 + term3 + term4 + term_q + term_fluid + term_osc + term_DM + term_wind + term_pressure;
    }

    // Debug/Output method (for transparency in base program)
    void printParameters(std::ostream& os = std::cout) const {
        os << std::fixed << std::setprecision(3);
        os << "NGC 3603 Parameters:" << std::endl;
        os << "G: " << G << ", M0: " << M0 << ", r: " << r << std::endl;
        os << "H0: " << H0 << ", B: " << B << ", B_crit: " << B_crit << std::endl;
        os << "f_TRZ: " << f_TRZ << ", M_dot_factor: " << M_dot_factor << ", tau_SF: " << tau_SF << std::endl;
        os << "rho_fluid: " << rho_fluid << ", rho_wind: " << rho_wind << ", v_wind: " << v_wind << std::endl;
        os << "P0: " << P0 << ", tau_exp: " << tau_exp << std::endl;
        os << "gas_v: " << gas_v << ", M_DM_factor: " << M_DM_factor << std::endl;
        os << "A_osc: " << A_osc << ", delta_rho_over_rho: " << delta_rho_over_rho << std::endl;
        os << "ug1_base: " << ug1_base << std::endl;
    }

    // Example computation at t=500k years (for testing)
    double exampleAt500kYears() const {
        double t_example = 5e5 * 3.156e7;
        return compute_g_NGC3603(t_example);
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
    void expandCavityPressureScale(double P0_scale, double tau_exp_scale);
    void expandWindMagneticScale(double rho_wind_scale, double v_wind_scale, double B_scale);

    // --- Self-Refinement (3 methods) ---
    void autoRefineParameters(const std::vector<std::pair<double, double>>& observations);
    void calibrateToObservations(const std::vector<double>& times, const std::vector<double>& g_obs);
    double optimizeForMetric(std::function<double(double)> metric, double t_start, double t_end, int steps);

    // --- Parameter Exploration (1 method) ---
    std::vector<std::map<std::string, double>> generateVariations(int count, double variation_pct);

    // --- Adaptive Evolution (2 methods) ---
    void mutateParameters(double mutation_rate, std::mt19937& rng);
    void evolveSystem(int generations, std::function<double(const NGC3603&)> fitness);

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

#endif // NGC_3603_H

// ========== IMPLEMENTATION OF ENHANCED METHODS (Outside class) ==========

// Anonymous namespace for state storage
namespace {
    std::map<std::string, std::map<std::string, double>> ngc3603_saved_states;
}

// --- Variable Management (5 methods) ---
bool NGC3603::createVariable(const std::string& name, double value) {
    return setVariable(name, value);
}

bool NGC3603::removeVariable(const std::string& name) {
    std::cerr << "Warning: Cannot remove built-in variable '" << name << "' in NGC3603 class." << std::endl;
    return false;
}

bool NGC3603::cloneVariable(const std::string& src, const std::string& dest) {
    double val = getVariable(src);
    return setVariable(dest, val);
}

std::vector<std::string> NGC3603::listVariables() const {
    return {"G", "M0", "r", "H0", "B", "B_crit", "Lambda", "c_light", "q_charge",
            "gas_v", "f_TRZ", "M_dot_factor", "tau_SF", "rho_wind", "v_wind",
            "rho_fluid", "P0", "tau_exp", "rho_vac_UA", "rho_vac_SCm", "scale_EM",
            "proton_mass", "hbar", "t_Hubble", "t_Hubble_gyr", "delta_x", "delta_p",
            "integral_psi", "A_osc", "k_osc", "omega_osc", "x_pos", "M_DM_factor",
            "delta_rho_over_rho"};
}

std::string NGC3603::getSystemName() const {
    return "NGC3603";
}

// --- Batch Operations (2 methods) ---
bool NGC3603::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        double val = getVariable(name);
        if (!setVariable(name, func(val))) return false;
    }
    return true;
}

bool NGC3603::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    return transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// --- Self-Expansion (4 methods) ---
void NGC3603::expandParameterSpace(double factor) {
    std::vector<std::string> expandable = {"M0", "r", "B", "rho_fluid", "rho_wind", "A_osc", "M_DM_factor"};
    scaleVariableGroup(expandable, factor);
}

void NGC3603::expandStarFormationScale(double M_dot_factor_scale, double tau_SF_scale) {
    setVariable("M_dot_factor", getVariable("M_dot_factor") * M_dot_factor_scale);
    setVariable("tau_SF", getVariable("tau_SF") * tau_SF_scale);
}

void NGC3603::expandCavityPressureScale(double P0_scale, double tau_exp_scale) {
    setVariable("P0", getVariable("P0") * P0_scale);
    setVariable("tau_exp", getVariable("tau_exp") * tau_exp_scale);
}

void NGC3603::expandWindMagneticScale(double rho_wind_scale, double v_wind_scale, double B_scale) {
    setVariable("rho_wind", getVariable("rho_wind") * rho_wind_scale);
    setVariable("v_wind", getVariable("v_wind") * v_wind_scale);
    setVariable("B", getVariable("B") * B_scale);
    setVariable("B_crit", getVariable("B_crit") * B_scale);
}

// --- Self-Refinement (3 methods) ---
void NGC3603::autoRefineParameters(const std::vector<std::pair<double, double>>& observations) {
    if (observations.empty()) return;
    
    double sum_error = 0.0;
    for (const auto& obs : observations) {
        double t = obs.first;
        double g_obs = obs.second;
        double g_calc = compute_g_NGC3603(t);
        sum_error += std::abs(g_calc - g_obs);
    }
    double avg_error = sum_error / observations.size();
    
    if (avg_error > 1e-6) {
        double adj_factor = 1.0 - std::min(0.1, avg_error / 1e6);
        setVariable("M_dot_factor", getVariable("M_dot_factor") * adj_factor);
        setVariable("tau_SF", getVariable("tau_SF") * (2.0 - adj_factor));
    }
}

void NGC3603::calibrateToObservations(const std::vector<double>& times, const std::vector<double>& g_obs) {
    if (times.size() != g_obs.size() || times.empty()) return;
    
    std::vector<std::pair<double, double>> obs;
    for (size_t i = 0; i < times.size(); ++i) {
        obs.push_back({times[i], g_obs[i]});
    }
    
    for (int iter = 0; iter < 5; ++iter) {
        autoRefineParameters(obs);
    }
}

double NGC3603::optimizeForMetric(std::function<double(double)> metric, double t_start, double t_end, int steps) {
    double best_score = -1e100;
    double dt = (t_end - t_start) / steps;
    
    for (int i = 0; i <= steps; ++i) {
        double t = t_start + i * dt;
        double g = compute_g_NGC3603(t);
        double score = metric(g);
        if (score > best_score) best_score = score;
    }
    return best_score;
}

// --- Parameter Exploration (1 method) ---
std::vector<std::map<std::string, double>> NGC3603::generateVariations(int count, double variation_pct) {
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
void NGC3603::mutateParameters(double mutation_rate, std::mt19937& rng) {
    std::uniform_real_distribution<> dis(-mutation_rate, mutation_rate);
    auto vars = listVariables();
    
    for (const auto& v : vars) {
        if (v == "c_light" || v == "G" || v == "hbar") continue;
        double val = getVariable(v);
        double delta = val * dis(rng);
        setVariable(v, val + delta);
    }
}

void NGC3603::evolveSystem(int generations, std::function<double(const NGC3603&)> fitness) {
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
bool NGC3603::saveState(const std::string& stateName) {
    std::map<std::string, double> state;
    auto vars = listVariables();
    for (const auto& v : vars) {
        state[v] = getVariable(v);
    }
    ngc3603_saved_states[stateName] = state;
    return true;
}

bool NGC3603::restoreState(const std::string& stateName) {
    auto it = ngc3603_saved_states.find(stateName);
    if (it == ngc3603_saved_states.end()) return false;
    
    for (const auto& pair : it->second) {
        setVariable(pair.first, pair.second);
    }
    return true;
}

std::vector<std::string> NGC3603::listSavedStates() const {
    std::vector<std::string> names;
    for (const auto& pair : ngc3603_saved_states) {
        names.push_back(pair.first);
    }
    return names;
}

std::string NGC3603::exportState() const {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(6);
    oss << "NGC3603 State Export:\n";
    auto vars = listVariables();
    for (const auto& v : vars) {
        oss << v << " = " << getVariable(v) << "\n";
    }
    return oss.str();
}

// --- System Analysis (4 methods) ---
std::map<std::string, double> NGC3603::sensitivityAnalysis(double t, double delta_pct) {
    std::map<std::string, double> sensitivities;
    double g_base = compute_g_NGC3603(t);
    
    auto vars = listVariables();
    for (const auto& v : vars) {
        if (v == "c_light" || v == "G" || v == "hbar") continue;
        
        double original = getVariable(v);
        double delta = original * delta_pct / 100.0;
        
        setVariable(v, original + delta);
        double g_plus = compute_g_NGC3603(t);
        setVariable(v, original);
        
        double sensitivity = (g_base != 0.0) ? std::abs((g_plus - g_base) / g_base) : 0.0;
        sensitivities[v] = sensitivity;
    }
    
    return sensitivities;
}

std::string NGC3603::generateReport(double t) const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6);
    oss << "============================================\n";
    oss << "NGC 3603 REPORT\n";
    oss << "Extreme Young Massive Star Cluster\n";
    oss << "============================================\n";
    oss << "Time: t = " << t << " s (" << (t/3.156e7/1e6) << " Myr)\n\n";
    
    oss << "Physical Parameters:\n";
    double M_sun = 1.989e30;
    oss << "  Initial Mass M0 = " << M0 << " kg (" << (M0/M_sun) << " M_sun)\n";
    oss << "  M(t) = " << M_t(t) << " kg (" << (M_t(t)/M_sun) << " M_sun)\n";
    double ly_to_m = 9.461e15;
    oss << "  Radius r = " << r << " m (" << (r/ly_to_m) << " ly)\n";
    oss << "  Hubble constant H0 = " << H0 << " s^-1\n";
    oss << "  Magnetic field B = " << B << " T (B_crit = " << B_crit << " T)\n";
    oss << "  f_TRZ = " << f_TRZ << "\n";
    oss << "  Star formation M_dot_factor = " << M_dot_factor << ", tau_SF = " << tau_SF << " s\n";
    oss << "  Cavity pressure P0 = " << P0 << " Pa, tau_exp = " << tau_exp << " s\n";
    oss << "  P(t) = " << P_t(t) << " Pa\n";
    oss << "  Wind density rho_wind = " << rho_wind << " kg/m^3, v_wind = " << v_wind << " m/s\n";
    oss << "  Fluid density rho_fluid = " << rho_fluid << " kg/m^3\n";
    oss << "  Gas velocity = " << gas_v << " m/s\n";
    oss << "  DM factor = " << M_DM_factor << "\n\n";
    
    oss << "Computed Acceleration:\n";
    oss << "  g_NGC3603(t) = " << compute_g_NGC3603(t) << " m/s^2\n\n";
    
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
    
    double Pt = P_t(t);
    oss << "  Cavity Pressure: " << (Pt / rho_fluid) << " m/s^2\n";
    
    oss << "============================================\n";
    return oss.str();
}

bool NGC3603::validateConsistency() const {
    bool valid = true;
    
    if (M0 <= 0 || r <= 0) { std::cerr << "Error: M0 and r must be positive.\n"; valid = false; }
    if (B < 0 || B_crit <= 0) { std::cerr << "Error: B, B_crit must be non-negative/positive.\n"; valid = false; }
    if (H0 < 0) { std::cerr << "Error: H0 must be non-negative.\n"; valid = false; }
    if (tau_SF <= 0 || tau_exp <= 0) { std::cerr << "Error: Timescales must be positive.\n"; valid = false; }
    if (P0 < 0) { std::cerr << "Warning: Initial pressure P0 is negative.\n"; }
    if (rho_fluid <= 0 || rho_wind < 0) { std::cerr << "Error: Fluid/wind densities must be positive/non-negative.\n"; valid = false; }
    if (v_wind < 0 || gas_v < 0) { std::cerr << "Error: Velocities must be non-negative.\n"; valid = false; }
    if (M_DM_factor < 0 || M_DM_factor > 1.0) { std::cerr << "Warning: DM factor outside [0,1].\n"; }
    
    return valid;
}

bool NGC3603::autoCorrectAnomalies() {
    bool corrected = false;
    
    double M_sun = 1.989e30;
    double ly_to_m = 9.461e15;
    
    if (M0 <= 0) { M0 = 400000.0 * M_sun; corrected = true; }
    if (r <= 0) { r = 9.5 * ly_to_m; corrected = true; }
    if (B < 0) { B = 1e-5; corrected = true; }
    if (B_crit <= 0) { B_crit = 1e11; corrected = true; }
    if (H0 < 0) { H0 = 2.184e-18; corrected = true; }
    if (tau_SF <= 0) { tau_SF = 1e6 * 3.156e7; corrected = true; }
    if (tau_exp <= 0) { tau_exp = 1e6 * 3.156e7; corrected = true; }
    if (P0 < 0) { P0 = 4e-8; corrected = true; }
    if (rho_fluid <= 0) { rho_fluid = 1e-20; corrected = true; }
    if (rho_wind < 0) { rho_wind = 1e-20; corrected = true; }
    if (v_wind < 0) { v_wind = 2e6; corrected = true; }
    if (gas_v < 0) { gas_v = 1e5; corrected = true; }
    if (M_DM_factor < 0) { M_DM_factor = 0.1; corrected = true; }
    if (M_DM_factor > 1.0) { M_DM_factor = 1.0; corrected = true; }
    
    if (corrected) updateCache();
    return corrected;
}

// ========== ENHANCED EXAMPLE FUNCTION ==========
void enhanced_ngc3603_example() {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "=========================================================\n";
    std::cout << "ENHANCED NGC 3603 DEMONSTRATION\n";
    std::cout << "Extreme Young Massive Star Cluster with Cavity Pressure\n";
    std::cout << "=========================================================\n\n";
    
    NGC3603 ngc3603;
    
    // Step 1: Initial state and validation
    std::cout << "Step 1: Initial State and Validation\n";
    std::cout << "System: " << ngc3603.getSystemName() << "\n";
    std::cout << "Validation: " << (ngc3603.validateConsistency() ? "PASS" : "FAIL") << "\n";
    std::cout << "Auto-corrected: " << (ngc3603.autoCorrectAnomalies() ? "Yes" : "No") << "\n\n";
    
    // Step 2: Time evolution showing M(t) and P(t)
    std::cout << "Step 2: Time Evolution (Mass Growth M(t) and Cavity Pressure P(t))\n";
    double t_Myr_array[] = {0.0, 0.5, 1.0, 2.0, 5.0};
    for (double t_Myr : t_Myr_array) {
        double t = t_Myr * 1e6 * 3.156e7;
        double M_sun = 1.989e30;
        double Mt = ngc3603.M_t(t);
        double Pt = ngc3603.P_t(t);
        double g = ngc3603.compute_g_NGC3603(t);
        std::cout << "  t = " << t_Myr << " Myr: M(t) = " << (Mt/M_sun) << " M_sun, P(t) = " << Pt 
                  << " Pa, g = " << g << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 3: Variable listing
    std::cout << "Step 3: Variable Listing\n";
    auto vars = ngc3603.listVariables();
    std::cout << "Total variables: " << vars.size() << "\n";
    std::cout << "Sample: " << vars[0] << ", " << vars[1] << ", " << vars[11] << " (M_dot_factor), " 
              << vars[16] << " (P0)\n\n";
    
    // Step 4: Star formation scaling
    std::cout << "Step 4: Star Formation Scaling (M_dot_factor sweeps)\n";
    ngc3603.saveState("original");
    double M_dot_factors[] = {0.5, 1.0, 2.0};
    for (double factor : M_dot_factors) {
        ngc3603.restoreState("original");
        ngc3603.expandStarFormationScale(factor, 1.0);
        double t = 1e6 * 3.156e7;
        double M_sun = 1.989e30;
        double Mt = ngc3603.M_t(t);
        std::cout << "  M_dot_factor × " << factor << ": M(1 Myr) = " << (Mt/M_sun) << " M_sun\n";
    }
    ngc3603.restoreState("original");
    std::cout << "\n";
    
    // Step 5: Cavity pressure scaling (UNIQUE to NGC 3603)
    std::cout << "Step 5: Cavity Pressure Scaling (P0 sweeps) - CAVITY FEATURE\n";
    double P0_factors[] = {0.5, 1.0, 2.0};
    for (double factor : P0_factors) {
        ngc3603.restoreState("original");
        ngc3603.expandCavityPressureScale(factor, 1.0);
        double t = 1e6 * 3.156e7;
        double Pt = ngc3603.P_t(t);
        double g = ngc3603.compute_g_NGC3603(t);
        std::cout << "  P0 × " << factor << ": P(1 Myr) = " << Pt << " Pa, g = " << g << " m/s^2\n";
    }
    ngc3603.restoreState("original");
    std::cout << "\n";
    
    // Step 6: Cavity expansion timescale sweeps
    std::cout << "Step 6: Cavity Expansion Timescale Sweeps (tau_exp)\n";
    double tau_exp_factors[] = {0.5, 1.0, 2.0};
    for (double factor : tau_exp_factors) {
        ngc3603.restoreState("original");
        ngc3603.expandCavityPressureScale(1.0, factor);
        double t = 1e6 * 3.156e7;
        double Pt = ngc3603.P_t(t);
        std::cout << "  tau_exp × " << factor << ": P(1 Myr) = " << Pt << " Pa\n";
    }
    ngc3603.restoreState("original");
    std::cout << "\n";
    
    // Step 7: Wind velocity scaling
    std::cout << "Step 7: Wind Velocity Scaling\n";
    double v_wind_factors[] = {0.5, 1.0, 2.0};
    for (double factor : v_wind_factors) {
        ngc3603.restoreState("original");
        ngc3603.expandWindMagneticScale(1.0, factor, 1.0);
        double t = 1e6 * 3.156e7;
        double g = ngc3603.compute_g_NGC3603(t);
        std::cout << "  v_wind × " << factor << ": g(1 Myr) = " << g << " m/s^2\n";
    }
    ngc3603.restoreState("original");
    std::cout << "\n";
    
    // Step 8: Magnetic field scaling
    std::cout << "Step 8: Magnetic Field Scaling\n";
    double B_factors[] = {0.5, 1.0, 2.0};
    for (double factor : B_factors) {
        ngc3603.restoreState("original");
        ngc3603.expandWindMagneticScale(1.0, 1.0, factor);
        double t = 1e6 * 3.156e7;
        double g = ngc3603.compute_g_NGC3603(t);
        std::cout << "  B × " << factor << ": g(1 Myr) = " << g << " m/s^2\n";
    }
    ngc3603.restoreState("original");
    std::cout << "\n";
    
    // Step 9: Parameter space expansion
    std::cout << "Step 9: Parameter Space Expansion (all scalable params)\n";
    ngc3603.expandParameterSpace(1.2);
    double M0_after = ngc3603.getVariable("M0");
    double M_sun = 1.989e30;
    std::cout << "  After 1.2× expansion: M0 = " << (M0_after/M_sun) << " M_sun\n";
    ngc3603.restoreState("original");
    std::cout << "\n";
    
    // Step 10: Batch operations
    std::cout << "Step 10: Batch Operations (scale multiple variables)\n";
    std::vector<std::string> scale_group = {"M0", "r", "B"};
    ngc3603.scaleVariableGroup(scale_group, 1.1);
    std::cout << "  Scaled {M0, r, B} by 1.1×\n";
    ngc3603.restoreState("original");
    std::cout << "\n";
    
    // Step 11: State management
    std::cout << "Step 11: State Management\n";
    ngc3603.saveState("state_A");
    ngc3603.expandStarFormationScale(1.5, 1.2);
    ngc3603.saveState("state_B");
    auto states = ngc3603.listSavedStates();
    std::cout << "  Saved states: ";
    for (const auto& s : states) std::cout << s << " ";
    std::cout << "\n";
    ngc3603.restoreState("state_A");
    std::cout << "  Restored state_A\n\n";
    
    // Step 12: Generate parameter variations
    std::cout << "Step 12: Generate Parameter Variations (5% variation)\n";
    auto variations = ngc3603.generateVariations(3, 5.0);
    std::cout << "  Generated " << variations.size() << " variants with 5% random variation\n";
    std::cout << "  Variant 1 M0 = " << variations[0]["M0"] << " kg\n\n";
    
    // Step 13: Sensitivity analysis
    std::cout << "Step 13: Sensitivity Analysis at 500k years\n";
    double t_sens = 5e5 * 3.156e7;
    auto sensitivities = ngc3603.sensitivityAnalysis(t_sens, 1.0);
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
        double g_obs = ngc3603.compute_g_NGC3603(t_obs) * (1.0 + 0.01 * (rand() % 100 - 50) / 100.0);
        obs.push_back({t_obs, g_obs});
    }
    ngc3603.autoRefineParameters(obs);
    std::cout << "  Refined parameters based on " << obs.size() << " observations\n\n";
    
    // Step 15: Optimization for maximum acceleration
    std::cout << "Step 15: Optimize for Maximum Acceleration\n";
    ngc3603.restoreState("original");
    auto metric = [](double g) { return g; };
    double t_opt_start = 0.0;
    double t_opt_end = 5e6 * 3.156e7;
    double best_g = ngc3603.optimizeForMetric(metric, t_opt_start, t_opt_end, 50);
    std::cout << "  Best g over 5 Myr: " << best_g << " m/s^2\n\n";
    
    // Step 16: Evolutionary system adaptation
    std::cout << "Step 16: Evolutionary System Adaptation (5 generations)\n";
    ngc3603.restoreState("original");
    auto fitness = [](const NGC3603& n) {
        double t = 5e5 * 3.156e7;
        return n.compute_g_NGC3603(t);
    };
    ngc3603.evolveSystem(5, fitness);
    std::cout << "  Evolved system over 5 generations (fitness = g at 500k years)\n\n";
    
    // Step 17: Full system report
    std::cout << "Step 17: Full System Report at 500k years\n";
    ngc3603.restoreState("original");
    double t_report = 5e5 * 3.156e7;
    std::string report = ngc3603.generateReport(t_report);
    std::cout << report << "\n";
    
    // Step 18: Full state export
    std::cout << "Step 18: Full State Export\n";
    std::string exported = ngc3603.exportState();
    std::cout << "Exported state (first 500 chars):\n";
    std::cout << exported.substr(0, 500) << "...\n\n";
    
    std::cout << "=========================================================\n";
    std::cout << "ENHANCED DEMONSTRATION COMPLETE\n";
    std::cout << "=========================================================\n";
}