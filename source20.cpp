/**
 * ================================================================================================
 * Header: GalaxyNGC2525.h
 *
 * Description: C++ Module for Galaxy NGC 2525 Class
 *              This is the tenth module in a series of 500+ code files for the Universal Quantum
 *              Field Framework (UQFF) simulations, focusing on barred spiral galaxy evolution and
 *              gravity equations derived from Hubble datasets, high-energy lab simulations, and
 *              UQFF refinements (dated May 09, 2025, updated for full term inclusion on October 08, 2025).
 *
 * Purpose: Encapsulates the Master Universal Gravity Equation (MUGE) for NGC 2525 evolution.
 *          Includes ALL terms: base gravity (static M), cosmic expansion (H(z)), magnetic correction (static B),
 *          black hole influence, UQFF Ug components with f_TRZ, Lambda, quantum uncertainty, scaled EM with [UA],
 *          fluid dynamics, oscillatory waves, DM/density perturbations, and supernova mass loss - (G * M_SN(t)) / r^2.
 *          Supports dynamic variable updates for all parameters.
 *
 * Integration: Designed for inclusion in base program 'ziqn233h.cpp' (not present here).
 *              Instantiate class in main: GalaxyNGC2525 ngc2525;
 *              Compute: double g = ngc2525.compute_g_NGC2525(t);
 *
 * Key Features:
 *   - Default values from UQFF document: M ? 1.0000225e10 Msun, r = 2.836e20 m, M_BH = 2.25e7 Msun,
 *     r_BH = 1.496e11 m, z ? 0.016, Hz ? 2.19e-18 s^-1, M_SN0 = 1.4 Msun, tau_SN = 1 yr, B = 1e-5 T.
 *   - Units handled: Msun to kg, ly to m; SN term as mass loss acceleration.
 *   - Setter methods for updates: setVar(double new_val) or addToVar(double delta)/subtractFromVar(double delta).
 *   - Computes g_NGC2525(r, t) with every term explicitly included.
 *
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: October 08, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

#ifndef GALAXY_NGC_2525_H
#define GALAXY_NGC_2525_H

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

class GalaxyNGC2525 {
private:
    // Core parameters (mutable for updates)
    double G;               // Gravitational constant
    double M;               // Total galaxy mass (kg)
    double r;               // Galaxy radius (m)
    double Hz;              // Hubble parameter at z (s^-1)
    double B;               // Static magnetic field (T)
    double B_crit;          // Critical B field (T)
    double Lambda;          // Cosmological constant
    double c_light;         // Speed of light
    double q_charge;        // Charge (proton)
    double gas_v;           // Gas velocity for EM (m/s)
    double f_TRZ;           // Time-reversal factor
    double M_BH;            // Black hole mass (kg)
    double r_BH;            // Black hole influence radius (m)
    double M_SN0;           // Initial SN mass (kg)
    double tau_SN;          // SN decay timescale (s)
    double rho_vac_UA;      // UA vacuum density (J/m^3)
    double rho_vac_SCm;     // SCm vacuum density (J/m^3)
    double scale_EM;        // EM scaling factor
    double proton_mass;     // Proton mass for EM acceleration
    double z_gal;           // Galaxy redshift

    // Additional parameters for full inclusion of terms
    double hbar;            // Reduced Planck's constant
    double t_Hubble;        // Hubble time (s)
    double delta_x;         // Position uncertainty (m)
    double delta_p;         // Momentum uncertainty (kg m/s)
    double integral_psi;    // Wavefunction integral approximation
    double rho_fluid;       // Fluid density (kg/m^3)
    double A_osc;           // Oscillatory amplitude (m/s^2)
    double k_osc;           // Wave number (1/m)
    double omega_osc;       // Angular frequency (rad/s)
    double x_pos;           // Position for oscillation (m)
    double t_Hubble_gyr;    // Hubble time in Gyr
    double M_DM_factor;     // Dark matter mass fraction
    double delta_rho_over_rho; // Density perturbation fraction

    // Computed caches (updated on demand)
    double ug1_base;        // Cached Ug1 = G*M/r^2
    double g_BH;            // Cached BH acceleration

public:
    // Constructor with default UQFF values
    GalaxyNGC2525() {
        initializeDefaults();
    }

    // Destructor (empty)
    ~GalaxyNGC2525() {}

    // Initialization method (called in constructor)
    void initializeDefaults() {
        G = 6.6743e-11;
        double M_sun = 1.989e30;
        M = (1e10 + 2.25e7) * M_sun;
        r = 2.836e20;
        z_gal = 0.016;
        double Hz_kms = 70 * sqrt(0.3 * pow(1 + z_gal, 3) + 0.7);  // km/s/Mpc
        Hz = (Hz_kms * 1000 / 3.086e19);  // s^-1
        B = 1e-5;
        B_crit = 1e11;
        Lambda = 1.1e-52;
        c_light = 3e8;
        q_charge = 1.602e-19;
        gas_v = 1e5;
        f_TRZ = 0.1;
        M_BH = 2.25e7 * M_sun;
        r_BH = 1.496e11;
        M_SN0 = 1.4 * M_sun;
        tau_SN = 1 * 3.156e7;
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
        rho_fluid = 1e-21;
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
        g_BH = (G * M_BH) / (r_BH * r_BH);
    }

    // Universal setter for any variable (by name, for flexibility)
    bool setVariable(const std::string& varName, double newValue) {
        if (varName == "G") { G = newValue; }
        else if (varName == "M") { M = newValue; }
        else if (varName == "r") { r = newValue; }
        else if (varName == "Hz") { Hz = newValue; }
        else if (varName == "B") { B = newValue; }
        else if (varName == "B_crit") { B_crit = newValue; }
        else if (varName == "Lambda") { Lambda = newValue; }
        else if (varName == "c_light") { c_light = newValue; }
        else if (varName == "q_charge") { q_charge = newValue; }
        else if (varName == "gas_v") { gas_v = newValue; }
        else if (varName == "f_TRZ") { f_TRZ = newValue; }
        else if (varName == "M_BH") { M_BH = newValue; }
        else if (varName == "r_BH") { r_BH = newValue; }
        else if (varName == "M_SN0") { M_SN0 = newValue; }
        else if (varName == "tau_SN") { tau_SN = newValue; }
        else if (varName == "rho_vac_UA") { rho_vac_UA = newValue; }
        else if (varName == "rho_vac_SCm") { rho_vac_SCm = newValue; }
        else if (varName == "scale_EM") { scale_EM = newValue; }
        else if (varName == "proton_mass") { proton_mass = newValue; }
        else if (varName == "z_gal") { z_gal = newValue; }
        // Full terms
        else if (varName == "hbar") { hbar = newValue; }
        else if (varName == "t_Hubble") { t_Hubble = newValue; }
        else if (varName == "t_Hubble_gyr") { t_Hubble_gyr = newValue; }
        else if (varName == "delta_x") { delta_x = newValue; }
        else if (varName == "delta_p") { delta_p = newValue; }
        else if (varName == "integral_psi") { integral_psi = newValue; }
        else if (varName == "rho_fluid") { rho_fluid = newValue; }
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
        else if (varName == "Hz") return Hz;
        else if (varName == "B") return B;
        else if (varName == "B_crit") return B_crit;
        else if (varName == "Lambda") return Lambda;
        else if (varName == "c_light") return c_light;
        else if (varName == "q_charge") return q_charge;
        else if (varName == "gas_v") return gas_v;
        else if (varName == "f_TRZ") return f_TRZ;
        else if (varName == "M_BH") return M_BH;
        else if (varName == "r_BH") return r_BH;
        else if (varName == "M_SN0") return M_SN0;
        else if (varName == "tau_SN") return tau_SN;
        else if (varName == "rho_vac_UA") return rho_vac_UA;
        else if (varName == "rho_vac_SCm") return rho_vac_SCm;
        else if (varName == "scale_EM") return scale_EM;
        else if (varName == "proton_mass") return proton_mass;
        else if (varName == "z_gal") return z_gal;
        // Full terms
        else if (varName == "hbar") return hbar;
        else if (varName == "t_Hubble") return t_Hubble;
        else if (varName == "t_Hubble_gyr") return t_Hubble_gyr;
        else if (varName == "delta_x") return delta_x;
        else if (varName == "delta_p") return delta_p;
        else if (varName == "integral_psi") return integral_psi;
        else if (varName == "rho_fluid") return rho_fluid;
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

    // M_SN(t) computation
    double M_SN_t(double t) const {
        return M_SN0 * exp(-t / tau_SN);
    }

    // Ug terms computation
    double compute_Ug() const {
        double Ug1 = ug1_base;
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
    double compute_g_NGC2525(double t) const {
        if (t < 0) {
            std::cerr << "Error: Time t must be non-negative." << std::endl;
            return 0.0;
        }

        double MSNt = M_SN_t(t);

        // Term 1: Base + Hz + B corrections
        double corr_H = 1 + Hz * t;
        double corr_B = 1 - B / B_crit;
        double term1 = ug1_base * corr_H * corr_B;

        // BH term
        double term_BH = g_BH;

        // Term 2: UQFF Ug with f_TRZ
        double term2 = compute_Ug();

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

        // SN mass loss term (negative acceleration)
        double term_SN = -(G * MSNt) / (r * r);

        // Total g_NGC2525 (all terms summed)
        return term1 + term_BH + term2 + term3 + term4 + term_q + term_fluid + term_osc + term_DM + term_SN;
    }

    // Debug/Output method (for transparency in base program)
    void printParameters(std::ostream& os = std::cout) const {
        os << std::fixed << std::setprecision(3);
        os << "NGC 2525 Parameters:" << std::endl;
        os << "G: " << G << ", M: " << M << ", r: " << r << std::endl;
        os << "Hz: " << Hz << ", B: " << B << ", B_crit: " << B_crit << std::endl;
        os << "f_TRZ: " << f_TRZ << ", M_BH: " << M_BH << ", r_BH: " << r_BH << std::endl;
        os << "M_SN0: " << M_SN0 << ", tau_SN: " << tau_SN << std::endl;
        os << "rho_fluid: " << rho_fluid << ", M_DM_factor: " << M_DM_factor << std::endl;
        os << "A_osc: " << A_osc << ", delta_rho_over_rho: " << delta_rho_over_rho << std::endl;
        os << "ug1_base: " << ug1_base << ", g_BH: " << g_BH << std::endl;
    }

    // Example computation at t=7 years (for testing)
    double exampleAt7Years() const {
        double t_example = 7 * 3.156e7;
        return compute_g_NGC2525(t_example);
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
    void expandGalacticScale(double M_scale, double r_scale);
    void expandBlackHoleScale(double M_BH_scale, double r_BH_scale);
    void expandSupernovaScale(double M_SN0_scale, double tau_SN_scale);

    // --- Self-Refinement (3 methods) ---
    void autoRefineParameters(const std::vector<std::pair<double, double>>& observations);
    void calibrateToObservations(const std::vector<double>& times, const std::vector<double>& g_obs);
    double optimizeForMetric(std::function<double(double)> metric, double t_start, double t_end, int steps);

    // --- Parameter Exploration (1 method) ---
    std::vector<std::map<std::string, double>> generateVariations(int count, double variation_pct);

    // --- Adaptive Evolution (2 methods) ---
    void mutateParameters(double mutation_rate, std::mt19937& rng);
    void evolveSystem(int generations, std::function<double(const GalaxyNGC2525&)> fitness);

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

#endif // GALAXY_NGC_2525_H

// ========== IMPLEMENTATION OF ENHANCED METHODS (Outside class) ==========

// Anonymous namespace for state storage
namespace {
    std::map<std::string, std::map<std::string, double>> ngc2525_saved_states;
}

// --- Variable Management (5 methods) ---
bool GalaxyNGC2525::createVariable(const std::string& name, double value) {
    return setVariable(name, value);
}

bool GalaxyNGC2525::removeVariable(const std::string& name) {
    std::cerr << "Warning: Cannot remove built-in variable '" << name << "' in GalaxyNGC2525 class." << std::endl;
    return false;
}

bool GalaxyNGC2525::cloneVariable(const std::string& src, const std::string& dest) {
    double val = getVariable(src);
    return setVariable(dest, val);
}

std::vector<std::string> GalaxyNGC2525::listVariables() const {
    return {"G", "M", "r", "Hz", "B", "B_crit", "Lambda", "c_light", "q_charge",
            "gas_v", "f_TRZ", "M_BH", "r_BH", "M_SN0", "tau_SN", "rho_vac_UA",
            "rho_vac_SCm", "scale_EM", "proton_mass", "z_gal", "hbar", "t_Hubble",
            "t_Hubble_gyr", "delta_x", "delta_p", "integral_psi", "rho_fluid",
            "A_osc", "k_osc", "omega_osc", "x_pos", "M_DM_factor", "delta_rho_over_rho"};
}

std::string GalaxyNGC2525::getSystemName() const {
    return "GalaxyNGC2525";
}

// --- Batch Operations (2 methods) ---
bool GalaxyNGC2525::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        double val = getVariable(name);
        if (!setVariable(name, func(val))) return false;
    }
    return true;
}

bool GalaxyNGC2525::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    return transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// --- Self-Expansion (4 methods) ---
void GalaxyNGC2525::expandParameterSpace(double factor) {
    std::vector<std::string> expandable = {"M", "r", "B", "rho_fluid", "A_osc", "M_DM_factor"};
    scaleVariableGroup(expandable, factor);
}

void GalaxyNGC2525::expandGalacticScale(double M_scale, double r_scale) {
    setVariable("M", getVariable("M") * M_scale);
    setVariable("r", getVariable("r") * r_scale);
}

void GalaxyNGC2525::expandBlackHoleScale(double M_BH_scale, double r_BH_scale) {
    setVariable("M_BH", getVariable("M_BH") * M_BH_scale);
    setVariable("r_BH", getVariable("r_BH") * r_BH_scale);
}

void GalaxyNGC2525::expandSupernovaScale(double M_SN0_scale, double tau_SN_scale) {
    setVariable("M_SN0", getVariable("M_SN0") * M_SN0_scale);
    setVariable("tau_SN", getVariable("tau_SN") * tau_SN_scale);
}

// --- Self-Refinement (3 methods) ---
void GalaxyNGC2525::autoRefineParameters(const std::vector<std::pair<double, double>>& observations) {
    if (observations.empty()) return;
    
    double sum_error = 0.0;
    for (const auto& obs : observations) {
        double t = obs.first;
        double g_obs = obs.second;
        double g_calc = compute_g_NGC2525(t);
        sum_error += std::abs(g_calc - g_obs);
    }
    double avg_error = sum_error / observations.size();
    
    if (avg_error > 1e-6) {
        double adj_factor = 1.0 - std::min(0.1, avg_error / 1e6);
        setVariable("M", getVariable("M") * adj_factor);
        setVariable("f_TRZ", getVariable("f_TRZ") * (2.0 - adj_factor));
    }
}

void GalaxyNGC2525::calibrateToObservations(const std::vector<double>& times, const std::vector<double>& g_obs) {
    if (times.size() != g_obs.size() || times.empty()) return;
    
    std::vector<std::pair<double, double>> obs;
    for (size_t i = 0; i < times.size(); ++i) {
        obs.push_back({times[i], g_obs[i]});
    }
    
    for (int iter = 0; iter < 5; ++iter) {
        autoRefineParameters(obs);
    }
}

double GalaxyNGC2525::optimizeForMetric(std::function<double(double)> metric, double t_start, double t_end, int steps) {
    double best_score = -1e100;
    double dt = (t_end - t_start) / steps;
    
    for (int i = 0; i <= steps; ++i) {
        double t = t_start + i * dt;
        double g = compute_g_NGC2525(t);
        double score = metric(g);
        if (score > best_score) best_score = score;
    }
    return best_score;
}

// --- Parameter Exploration (1 method) ---
std::vector<std::map<std::string, double>> GalaxyNGC2525::generateVariations(int count, double variation_pct) {
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
void GalaxyNGC2525::mutateParameters(double mutation_rate, std::mt19937& rng) {
    std::uniform_real_distribution<> dis(-mutation_rate, mutation_rate);
    auto vars = listVariables();
    
    for (const auto& v : vars) {
        if (v == "c_light" || v == "G" || v == "hbar") continue;
        double val = getVariable(v);
        double delta = val * dis(rng);
        setVariable(v, val + delta);
    }
}

void GalaxyNGC2525::evolveSystem(int generations, std::function<double(const GalaxyNGC2525&)> fitness) {
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
bool GalaxyNGC2525::saveState(const std::string& stateName) {
    std::map<std::string, double> state;
    auto vars = listVariables();
    for (const auto& v : vars) {
        state[v] = getVariable(v);
    }
    ngc2525_saved_states[stateName] = state;
    return true;
}

bool GalaxyNGC2525::restoreState(const std::string& stateName) {
    auto it = ngc2525_saved_states.find(stateName);
    if (it == ngc2525_saved_states.end()) return false;
    
    for (const auto& pair : it->second) {
        setVariable(pair.first, pair.second);
    }
    return true;
}

std::vector<std::string> GalaxyNGC2525::listSavedStates() const {
    std::vector<std::string> names;
    for (const auto& pair : ngc2525_saved_states) {
        names.push_back(pair.first);
    }
    return names;
}

std::string GalaxyNGC2525::exportState() const {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(6);
    oss << "GalaxyNGC2525 State Export:\n";
    auto vars = listVariables();
    for (const auto& v : vars) {
        oss << v << " = " << getVariable(v) << "\n";
    }
    return oss.str();
}

// --- System Analysis (4 methods) ---
std::map<std::string, double> GalaxyNGC2525::sensitivityAnalysis(double t, double delta_pct) {
    std::map<std::string, double> sensitivities;
    double g_base = compute_g_NGC2525(t);
    
    auto vars = listVariables();
    for (const auto& v : vars) {
        if (v == "c_light" || v == "G" || v == "hbar") continue;
        
        double original = getVariable(v);
        double delta = original * delta_pct / 100.0;
        
        setVariable(v, original + delta);
        double g_plus = compute_g_NGC2525(t);
        setVariable(v, original);
        
        double sensitivity = (g_base != 0.0) ? std::abs((g_plus - g_base) / g_base) : 0.0;
        sensitivities[v] = sensitivity;
    }
    
    return sensitivities;
}

std::string GalaxyNGC2525::generateReport(double t) const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6);
    oss << "============================================\n";
    oss << "GALAXY NGC 2525 REPORT\n";
    oss << "Barred Spiral Galaxy with Supernova\n";
    oss << "============================================\n";
    oss << "Time: t = " << t << " s (" << (t/3.156e7) << " years)\n\n";
    
    oss << "Physical Parameters:\n";
    double M_sun = 1.989e30;
    oss << "  Galaxy Mass M = " << M << " kg (" << (M/M_sun) << " M_sun)\n";
    double ly_to_m = 9.461e15;
    oss << "  Galaxy Radius r = " << r << " m (" << (r/ly_to_m) << " ly)\n";
    oss << "  Galaxy Redshift z_gal = " << z_gal << "\n";
    oss << "  Hubble Parameter Hz = " << Hz << " s^-1\n";
    oss << "  Magnetic field B = " << B << " T (B_crit = " << B_crit << " T)\n";
    oss << "  Black Hole M_BH = " << M_BH << " kg (" << (M_BH/M_sun) << " M_sun)\n";
    oss << "  Black Hole Radius r_BH = " << r_BH << " m (" << (r_BH/1.496e11) << " AU)\n";
    oss << "  f_TRZ = " << f_TRZ << "\n";
    oss << "  Supernova M_SN0 = " << M_SN0 << " kg (" << (M_SN0/M_sun) << " M_sun)\n";
    oss << "  SN Timescale tau_SN = " << tau_SN << " s (" << (tau_SN/3.156e7) << " years)\n";
    oss << "  M_SN(t) = " << M_SN_t(t) << " kg (" << (M_SN_t(t)/M_sun) << " M_sun)\n";
    oss << "  Fluid density rho_fluid = " << rho_fluid << " kg/m^3\n";
    oss << "  Gas velocity = " << gas_v << " m/s\n";
    oss << "  DM factor = " << M_DM_factor << "\n\n";
    
    oss << "Computed Acceleration:\n";
    oss << "  g_NGC2525(t) = " << compute_g_NGC2525(t) << " m/s^2\n\n";
    
    oss << "UQFF Terms:\n";
    double corr_H = 1 + Hz * t;
    double corr_B = 1 - B / B_crit;
    oss << "  Base (with Hz, B): " << (ug1_base * corr_H * corr_B) << " m/s^2\n";
    oss << "  Black Hole: " << g_BH << " m/s^2\n";
    oss << "  Ug total: " << compute_Ug() << " m/s^2\n";
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
    
    double MSNt = M_SN_t(t);
    oss << "  Supernova Mass Loss: " << (-(G * MSNt) / (r * r)) << " m/s^2\n";
    
    oss << "============================================\n";
    return oss.str();
}

bool GalaxyNGC2525::validateConsistency() const {
    bool valid = true;
    
    if (M <= 0 || r <= 0) { std::cerr << "Error: M and r must be positive.\n"; valid = false; }
    if (B < 0 || B_crit <= 0) { std::cerr << "Error: B, B_crit must be non-negative/positive.\n"; valid = false; }
    if (Hz < 0) { std::cerr << "Error: Hz must be non-negative.\n"; valid = false; }
    if (z_gal < 0) { std::cerr << "Warning: Galaxy redshift z_gal is negative.\n"; }
    if (M_BH <= 0 || r_BH <= 0) { std::cerr << "Error: M_BH and r_BH must be positive.\n"; valid = false; }
    if (M_SN0 < 0 || tau_SN <= 0) { std::cerr << "Error: M_SN0 must be non-negative, tau_SN positive.\n"; valid = false; }
    if (rho_fluid <= 0) { std::cerr << "Error: Fluid density must be positive.\n"; valid = false; }
    if (gas_v < 0) { std::cerr << "Error: Gas velocity must be non-negative.\n"; valid = false; }
    if (M_DM_factor < 0 || M_DM_factor > 1.0) { std::cerr << "Warning: DM factor outside [0,1].\n"; }
    
    return valid;
}

bool GalaxyNGC2525::autoCorrectAnomalies() {
    bool corrected = false;
    
    double M_sun = 1.989e30;
    
    if (M <= 0) { M = (1e10 + 2.25e7) * M_sun; corrected = true; }
    if (r <= 0) { r = 2.836e20; corrected = true; }
    if (B < 0) { B = 1e-5; corrected = true; }
    if (B_crit <= 0) { B_crit = 1e11; corrected = true; }
    if (Hz < 0) { Hz = 2.19e-18; corrected = true; }
    if (z_gal < 0) { z_gal = 0.016; corrected = true; }
    if (M_BH <= 0) { M_BH = 2.25e7 * M_sun; corrected = true; }
    if (r_BH <= 0) { r_BH = 1.496e11; corrected = true; }
    if (M_SN0 < 0) { M_SN0 = 1.4 * M_sun; corrected = true; }
    if (tau_SN <= 0) { tau_SN = 1 * 3.156e7; corrected = true; }
    if (rho_fluid <= 0) { rho_fluid = 1e-21; corrected = true; }
    if (gas_v < 0) { gas_v = 1e5; corrected = true; }
    if (M_DM_factor < 0) { M_DM_factor = 0.1; corrected = true; }
    if (M_DM_factor > 1.0) { M_DM_factor = 1.0; corrected = true; }
    
    if (corrected) updateCache();
    return corrected;
}

// ========== ENHANCED EXAMPLE FUNCTION ==========
void enhanced_ngc2525_example() {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "=========================================================\n";
    std::cout << "ENHANCED GALAXY NGC 2525 DEMONSTRATION\n";
    std::cout << "Barred Spiral Galaxy with Supernova Mass Loss\n";
    std::cout << "=========================================================\n\n";
    
    GalaxyNGC2525 ngc2525;
    
    // Step 1: Initial state and validation
    std::cout << "Step 1: Initial State and Validation\n";
    std::cout << "System: " << ngc2525.getSystemName() << "\n";
    std::cout << "Validation: " << (ngc2525.validateConsistency() ? "PASS" : "FAIL") << "\n";
    std::cout << "Auto-corrected: " << (ngc2525.autoCorrectAnomalies() ? "Yes" : "No") << "\n\n";
    
    // Step 2: Time evolution showing supernova mass loss M_SN(t)
    std::cout << "Step 2: Time Evolution (Supernova Mass Loss)\n";
    double t_years_array[] = {0.0, 1.0, 5.0, 10.0, 20.0};
    for (double t_years : t_years_array) {
        double t = t_years * 3.156e7;
        double M_sun = 1.989e30;
        double MSNt = ngc2525.M_SN_t(t);
        double g = ngc2525.compute_g_NGC2525(t);
        std::cout << "  t = " << t_years << " years: M_SN(t) = " << (MSNt/M_sun) << " M_sun, g = " << g << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 3: Variable listing
    std::cout << "Step 3: Variable Listing\n";
    auto vars = ngc2525.listVariables();
    std::cout << "Total variables: " << vars.size() << "\n";
    std::cout << "Sample: " << vars[0] << ", " << vars[1] << ", " << vars[11] << " (M_BH), " 
              << vars[13] << " (M_SN0)\n\n";
    
    // Step 4: Galactic mass scaling
    std::cout << "Step 4: Galactic Mass Scaling (M sweeps)\n";
    ngc2525.saveState("original");
    double M_factors[] = {0.5, 1.0, 2.0};
    for (double factor : M_factors) {
        ngc2525.restoreState("original");
        ngc2525.expandGalacticScale(factor, 1.0);
        double t = 7 * 3.156e7;
        double g = ngc2525.compute_g_NGC2525(t);
        double M_sun = 1.989e30;
        double M = ngc2525.getVariable("M");
        std::cout << "  M × " << factor << ": M = " << (M/M_sun) << " M_sun, g(7 years) = " << g << " m/s^2\n";
    }
    ngc2525.restoreState("original");
    std::cout << "\n";
    
    // Step 5: Black hole scaling (UNIQUE to galaxy with central BH)
    std::cout << "Step 5: Black Hole Scaling (M_BH sweeps) - CENTRAL BH FEATURE\n";
    double M_BH_factors[] = {0.5, 1.0, 2.0};
    for (double factor : M_BH_factors) {
        ngc2525.restoreState("original");
        ngc2525.expandBlackHoleScale(factor, 1.0);
        double t = 7 * 3.156e7;
        double g = ngc2525.compute_g_NGC2525(t);
        double M_sun = 1.989e30;
        double M_BH = ngc2525.getVariable("M_BH");
        std::cout << "  M_BH × " << factor << ": M_BH = " << (M_BH/M_sun) << " M_sun, g(7 years) = " << g << " m/s^2\n";
    }
    ngc2525.restoreState("original");
    std::cout << "\n";
    
    // Step 6: Supernova mass loss scaling (UNIQUE to SN galaxy)
    std::cout << "Step 6: Supernova Initial Mass Scaling (M_SN0 sweeps) - SN FEATURE\n";
    double M_SN0_factors[] = {0.5, 1.0, 2.0};
    for (double factor : M_SN0_factors) {
        ngc2525.restoreState("original");
        ngc2525.expandSupernovaScale(factor, 1.0);
        double t = 7 * 3.156e7;
        double M_sun = 1.989e30;
        double MSNt = ngc2525.M_SN_t(t);
        double g = ngc2525.compute_g_NGC2525(t);
        std::cout << "  M_SN0 × " << factor << ": M_SN(7yr) = " << (MSNt/M_sun) << " M_sun, g = " << g << " m/s^2\n";
    }
    ngc2525.restoreState("original");
    std::cout << "\n";
    
    // Step 7: Supernova decay timescale sweeps
    std::cout << "Step 7: Supernova Timescale Sweeps (tau_SN)\n";
    double tau_SN_factors[] = {0.5, 1.0, 2.0};
    for (double factor : tau_SN_factors) {
        ngc2525.restoreState("original");
        ngc2525.expandSupernovaScale(1.0, factor);
        double t = 7 * 3.156e7;
        double M_sun = 1.989e30;
        double MSNt = ngc2525.M_SN_t(t);
        std::cout << "  tau_SN × " << factor << ": M_SN(7yr) = " << (MSNt/M_sun) << " M_sun\n";
    }
    ngc2525.restoreState("original");
    std::cout << "\n";
    
    // Step 8: Black hole radius scaling
    std::cout << "Step 8: Black Hole Influence Radius Scaling\n";
    double r_BH_factors[] = {0.5, 1.0, 2.0};
    for (double factor : r_BH_factors) {
        ngc2525.restoreState("original");
        ngc2525.expandBlackHoleScale(1.0, factor);
        double t = 7 * 3.156e7;
        double g = ngc2525.compute_g_NGC2525(t);
        std::cout << "  r_BH × " << factor << ": g(7 years) = " << g << " m/s^2\n";
    }
    ngc2525.restoreState("original");
    std::cout << "\n";
    
    // Step 9: Parameter space expansion
    std::cout << "Step 9: Parameter Space Expansion (all scalable params)\n";
    ngc2525.expandParameterSpace(1.2);
    double M_after = ngc2525.getVariable("M");
    double M_sun = 1.989e30;
    std::cout << "  After 1.2× expansion: M = " << (M_after/M_sun) << " M_sun\n";
    ngc2525.restoreState("original");
    std::cout << "\n";
    
    // Step 10: Batch operations
    std::cout << "Step 10: Batch Operations (scale multiple variables)\n";
    std::vector<std::string> scale_group = {"M", "r", "M_BH"};
    ngc2525.scaleVariableGroup(scale_group, 1.1);
    std::cout << "  Scaled {M, r, M_BH} by 1.1×\n";
    ngc2525.restoreState("original");
    std::cout << "\n";
    
    // Step 11: State management
    std::cout << "Step 11: State Management\n";
    ngc2525.saveState("state_A");
    ngc2525.expandGalacticScale(1.5, 1.2);
    ngc2525.saveState("state_B");
    auto states = ngc2525.listSavedStates();
    std::cout << "  Saved states: ";
    for (const auto& s : states) std::cout << s << " ";
    std::cout << "\n";
    ngc2525.restoreState("state_A");
    std::cout << "  Restored state_A\n\n";
    
    // Step 12: Generate parameter variations
    std::cout << "Step 12: Generate Parameter Variations (5% variation)\n";
    auto variations = ngc2525.generateVariations(3, 5.0);
    std::cout << "  Generated " << variations.size() << " variants with 5% random variation\n";
    std::cout << "  Variant 1 M = " << variations[0]["M"] << " kg\n\n";
    
    // Step 13: Sensitivity analysis
    std::cout << "Step 13: Sensitivity Analysis at 7 years\n";
    double t_sens = 7 * 3.156e7;
    auto sensitivities = ngc2525.sensitivityAnalysis(t_sens, 1.0);
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
        double t_obs = i * 5 * 3.156e7;
        double g_obs = ngc2525.compute_g_NGC2525(t_obs) * (1.0 + 0.01 * (rand() % 100 - 50) / 100.0);
        obs.push_back({t_obs, g_obs});
    }
    ngc2525.autoRefineParameters(obs);
    std::cout << "  Refined parameters based on " << obs.size() << " observations\n\n";
    
    // Step 15: Optimization for maximum acceleration
    std::cout << "Step 15: Optimize for Maximum Acceleration\n";
    ngc2525.restoreState("original");
    auto metric = [](double g) { return g; };
    double t_opt_start = 0.0;
    double t_opt_end = 20 * 3.156e7;
    double best_g = ngc2525.optimizeForMetric(metric, t_opt_start, t_opt_end, 50);
    std::cout << "  Best g over 20 years: " << best_g << " m/s^2\n\n";
    
    // Step 16: Evolutionary system adaptation
    std::cout << "Step 16: Evolutionary System Adaptation (5 generations)\n";
    ngc2525.restoreState("original");
    auto fitness = [](const GalaxyNGC2525& n) {
        double t = 7 * 3.156e7;
        return n.compute_g_NGC2525(t);
    };
    ngc2525.evolveSystem(5, fitness);
    std::cout << "  Evolved system over 5 generations (fitness = g at 7 years)\n\n";
    
    // Step 17: Full system report
    std::cout << "Step 17: Full System Report at 7 years\n";
    ngc2525.restoreState("original");
    double t_report = 7 * 3.156e7;
    std::string report = ngc2525.generateReport(t_report);
    std::cout << report << "\n";
    
    // Step 18: Full state export
    std::cout << "Step 18: Full State Export\n";
    std::string exported = ngc2525.exportState();
    std::cout << "Exported state (first 500 chars):\n";
    std::cout << exported.substr(0, 500) << "...\n\n";
    
    std::cout << "=========================================================\n";
    std::cout << "ENHANCED DEMONSTRATION COMPLETE\n";
    std::cout << "=========================================================\n";
}