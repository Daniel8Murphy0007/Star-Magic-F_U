/**
 * ================================================================================================
 * Header: RingsOfRelativity.h
 *
 * Description: C++ Module for "Rings of Relativity" (GAL-CLUS-022058s Einstein Ring) Class
 *              This is the eighth module in a series of 500+ code files for the Universal Quantum
 *              Field Framework (UQFF) simulations, focusing on gravitational lensing ring evolution
 *              and gravity equations derived from Hubble datasets, high-energy lab simulations, and
 *              UQFF refinements (dated May 09, 2025, updated for full term inclusion on October 08, 2025).
 *
 * Purpose: Encapsulates the Master Universal Gravity Equation (MUGE) for "Rings of Relativity" evolution.
 *          Includes ALL terms: base gravity (static M), cosmic expansion (H(z)), magnetic correction (static B),
 *          lensing amplification L(t) (static), UQFF Ug components with f_TRZ, Lambda, quantum uncertainty,
 *          scaled EM with [UA], fluid dynamics, oscillatory waves, DM/density perturbations, and stellar
 *          wind feedback (pressure / density for acc, added for completeness). Supports dynamic variable updates.
 *
 * Integration: Designed for inclusion in base program 'ziqn233h.cpp' (not present here).
 *              Instantiate class in main: RingsOfRelativity rings;
 *              Compute: double g = rings.compute_g_Rings(t);
 *
 * Key Features:
 *   - Default values from UQFF document: M = 1e14 Msun, r = 3.086e20 m (~10 kpc), B = 1e-5 T,
 *     H_z for z=0.5 ? 2.42e-18 s^-1, L = (G*M)/(c^2*r) * 0.67, etc.
 *   - Units handled: Msun to kg, kpc to m; wind term as (rho * v_wind^2) / rho_fluid for acceleration.
 *   - Setter methods for updates: setVar(double new_val) or addToVar(double delta)/subtractFromVar(double delta).
 *   - Computes g_Rings(r, t) with every term explicitly included.
 *
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: October 08, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

#ifndef RINGS_OF_RELATIVITY_H
#define RINGS_OF_RELATIVITY_H

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

class RingsOfRelativity {
private:
    // Core parameters (mutable for updates)
    double G;               // Gravitational constant
    double M;               // Lensing mass (kg)
    double r;               // Einstein radius (m)
    double Hz;              // Hubble parameter at z (s^-1)
    double B;               // Static magnetic field (T)
    double B_crit;          // Critical B field (T)
    double Lambda;          // Cosmological constant
    double c_light;         // Speed of light
    double q_charge;        // Charge (proton)
    double gas_v;           // Gas velocity for EM (m/s)
    double f_TRZ;           // Time-reversal factor
    double L_factor;        // Lensing factor (D_LS / D_S ? 0.67)
    double rho_vac_UA;      // UA vacuum density (J/m^3)
    double rho_vac_SCm;     // SCm vacuum density (J/m^3)
    double scale_EM;        // EM scaling factor
    double proton_mass;     // Proton mass for EM acceleration
    double z_lens;          // Lens redshift

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
    double rho_wind;        // Wind density (kg/m^3)
    double v_wind;          // Wind velocity (m/s)

    // Computed caches (updated on demand)
    double ug1_base;        // Cached Ug1 = G*M/r^2
    double L_t;             // Cached lensing term

public:
    // Constructor with default UQFF values
    RingsOfRelativity() {
        initializeDefaults();
    }

    // Destructor (empty)
    ~RingsOfRelativity() {}

    // Initialization method (called in constructor)
    void initializeDefaults() {
        G = 6.6743e-11;
        double M_sun = 1.989e30;
        M = 1e14 * M_sun;
        r = 3.086e20;
        z_lens = 0.5;
        double Hz_kms = 70 * sqrt(0.3 * pow(1 + z_lens, 3) + 0.7);  // km/s/Mpc
        Hz = (Hz_kms * 1000 / 3.086e19);  // s^-1
        B = 1e-5;
        B_crit = 1e11;
        Lambda = 1.1e-52;
        c_light = 3e8;
        q_charge = 1.602e-19;
        gas_v = 1e5;
        f_TRZ = 0.1;
        L_factor = 0.67;
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
        A_osc = 1e-12;  // Small for lensing scale
        k_osc = 1.0 / r;
        omega_osc = 2 * M_PI / (r / c_light);
        x_pos = r;
        M_DM_factor = 0.1;
        delta_rho_over_rho = 1e-5;
        rho_wind = 1e-21;
        v_wind = 2e6;

        updateCache();
    }

    // Cache update for efficiency (call after parameter changes)
    void updateCache() {
        ug1_base = (G * M) / (r * r);
        L_t = ((G * M) / (pow(c_light, 2) * r)) * L_factor;
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
        else if (varName == "L_factor") { L_factor = newValue; }
        else if (varName == "rho_vac_UA") { rho_vac_UA = newValue; }
        else if (varName == "rho_vac_SCm") { rho_vac_SCm = newValue; }
        else if (varName == "scale_EM") { scale_EM = newValue; }
        else if (varName == "proton_mass") { proton_mass = newValue; }
        else if (varName == "z_lens") { z_lens = newValue; }
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
        else if (varName == "rho_wind") { rho_wind = newValue; }
        else if (varName == "v_wind") { v_wind = newValue; }
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
        else if (varName == "L_factor") return L_factor;
        else if (varName == "rho_vac_UA") return rho_vac_UA;
        else if (varName == "rho_vac_SCm") return rho_vac_SCm;
        else if (varName == "scale_EM") return scale_EM;
        else if (varName == "proton_mass") return proton_mass;
        else if (varName == "z_lens") return z_lens;
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
        else if (varName == "rho_wind") return rho_wind;
        else if (varName == "v_wind") return v_wind;
        else {
            std::cerr << "Error: Unknown variable '" << varName << "'." << std::endl;
            return 0.0;
        }
    }

    // Ug terms computation
    double compute_Ug(double /*Mt*/) const {  // Mt static as M
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
    double compute_g_Rings(double t) const {
        if (t < 0) {
            std::cerr << "Error: Time t must be non-negative." << std::endl;
            return 0.0;
        }

        // Term 1: Base + Hz + B + L corrections
        double corr_H = 1 + Hz * t;
        double corr_B = 1 - B / B_crit;
        double corr_L = 1 + L_t;
        double term1 = ug1_base * corr_H * corr_B * corr_L;

        // Term 2: UQFF Ug with f_TRZ
        double term2 = compute_Ug(0);  // No Mt variation

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

        // Total g_Rings (all terms summed)
        return term1 + term2 + term3 + term4 + term_q + term_fluid + term_osc + term_DM + term_wind;
    }

    // Debug/Output method (for transparency in base program)
    void printParameters(std::ostream& os = std::cout) const {
        os << std::fixed << std::setprecision(3);
        os << "Rings of Relativity Parameters:" << std::endl;
        os << "G: " << G << ", M: " << M << ", r: " << r << std::endl;
        os << "Hz: " << Hz << ", B: " << B << ", B_crit: " << B_crit << std::endl;
        os << "f_TRZ: " << f_TRZ << ", L_t: " << L_t << ", L_factor: " << L_factor << std::endl;
        os << "rho_fluid: " << rho_fluid << ", rho_wind: " << rho_wind << ", v_wind: " << v_wind << std::endl;
        os << "gas_v: " << gas_v << ", M_DM_factor: " << M_DM_factor << std::endl;
        os << "A_osc: " << A_osc << ", delta_rho_over_rho: " << delta_rho_over_rho << std::endl;
        os << "ug1_base: " << ug1_base << std::endl;
    }

    // Example computation at t=5 Gyr (for testing)
    double exampleAt5Gyr() const {
        double t_example = 5e9 * 3.156e7;
        return compute_g_Rings(t_example);
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
    void expandLensingScale(double M_scale, double L_factor_scale);
    void expandRedshiftScale(double z_scale, double Hz_scale);
    void expandMagneticWindScale(double B_scale, double rho_wind_scale, double v_wind_scale);

    // --- Self-Refinement (3 methods) ---
    void autoRefineParameters(const std::vector<std::pair<double, double>>& observations);
    void calibrateToObservations(const std::vector<double>& times, const std::vector<double>& g_obs);
    double optimizeForMetric(std::function<double(double)> metric, double t_start, double t_end, int steps);

    // --- Parameter Exploration (1 method) ---
    std::vector<std::map<std::string, double>> generateVariations(int count, double variation_pct);

    // --- Adaptive Evolution (2 methods) ---
    void mutateParameters(double mutation_rate, std::mt19937& rng);
    void evolveSystem(int generations, std::function<double(const RingsOfRelativity&)> fitness);

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

#endif // RINGS_OF_RELATIVITY_H

// ========== IMPLEMENTATION OF ENHANCED METHODS (Outside class) ==========

// Anonymous namespace for state storage
namespace {
    std::map<std::string, std::map<std::string, double>> rings_saved_states;
}

// --- Variable Management (5 methods) ---
bool RingsOfRelativity::createVariable(const std::string& name, double value) {
    return setVariable(name, value);
}

bool RingsOfRelativity::removeVariable(const std::string& name) {
    std::cerr << "Warning: Cannot remove built-in variable '" << name << "' in RingsOfRelativity class." << std::endl;
    return false;
}

bool RingsOfRelativity::cloneVariable(const std::string& src, const std::string& dest) {
    double val = getVariable(src);
    return setVariable(dest, val);
}

std::vector<std::string> RingsOfRelativity::listVariables() const {
    return {"G", "M", "r", "Hz", "B", "B_crit", "Lambda", "c_light", "q_charge",
            "gas_v", "f_TRZ", "L_factor", "rho_vac_UA", "rho_vac_SCm", "scale_EM",
            "proton_mass", "z_lens", "hbar", "t_Hubble", "t_Hubble_gyr", "delta_x",
            "delta_p", "integral_psi", "rho_fluid", "A_osc", "k_osc", "omega_osc",
            "x_pos", "M_DM_factor", "delta_rho_over_rho", "rho_wind", "v_wind"};
}

std::string RingsOfRelativity::getSystemName() const {
    return "RingsOfRelativity";
}

// --- Batch Operations (2 methods) ---
bool RingsOfRelativity::transformVariableGroup(const std::vector<std::string>& names, std::function<double(double)> func) {
    for (const auto& name : names) {
        double val = getVariable(name);
        if (!setVariable(name, func(val))) return false;
    }
    return true;
}

bool RingsOfRelativity::scaleVariableGroup(const std::vector<std::string>& names, double factor) {
    return transformVariableGroup(names, [factor](double v) { return v * factor; });
}

// --- Self-Expansion (4 methods) ---
void RingsOfRelativity::expandParameterSpace(double factor) {
    std::vector<std::string> expandable = {"M", "r", "B", "rho_fluid", "rho_wind", "A_osc", "M_DM_factor"};
    scaleVariableGroup(expandable, factor);
}

void RingsOfRelativity::expandLensingScale(double M_scale, double L_factor_scale) {
    setVariable("M", getVariable("M") * M_scale);
    setVariable("L_factor", getVariable("L_factor") * L_factor_scale);
}

void RingsOfRelativity::expandRedshiftScale(double z_scale, double Hz_scale) {
    setVariable("z_lens", getVariable("z_lens") * z_scale);
    setVariable("Hz", getVariable("Hz") * Hz_scale);
}

void RingsOfRelativity::expandMagneticWindScale(double B_scale, double rho_wind_scale, double v_wind_scale) {
    setVariable("B", getVariable("B") * B_scale);
    setVariable("B_crit", getVariable("B_crit") * B_scale);
    setVariable("rho_wind", getVariable("rho_wind") * rho_wind_scale);
    setVariable("v_wind", getVariable("v_wind") * v_wind_scale);
}

// --- Self-Refinement (3 methods) ---
void RingsOfRelativity::autoRefineParameters(const std::vector<std::pair<double, double>>& observations) {
    if (observations.empty()) return;
    
    double sum_error = 0.0;
    for (const auto& obs : observations) {
        double t = obs.first;
        double g_obs = obs.second;
        double g_calc = compute_g_Rings(t);
        sum_error += std::abs(g_calc - g_obs);
    }
    double avg_error = sum_error / observations.size();
    
    if (avg_error > 1e-6) {
        double adj_factor = 1.0 - std::min(0.1, avg_error / 1e6);
        setVariable("L_factor", getVariable("L_factor") * adj_factor);
        setVariable("f_TRZ", getVariable("f_TRZ") * (2.0 - adj_factor));
    }
}

void RingsOfRelativity::calibrateToObservations(const std::vector<double>& times, const std::vector<double>& g_obs) {
    if (times.size() != g_obs.size() || times.empty()) return;
    
    std::vector<std::pair<double, double>> obs;
    for (size_t i = 0; i < times.size(); ++i) {
        obs.push_back({times[i], g_obs[i]});
    }
    
    for (int iter = 0; iter < 5; ++iter) {
        autoRefineParameters(obs);
    }
}

double RingsOfRelativity::optimizeForMetric(std::function<double(double)> metric, double t_start, double t_end, int steps) {
    double best_score = -1e100;
    double dt = (t_end - t_start) / steps;
    
    for (int i = 0; i <= steps; ++i) {
        double t = t_start + i * dt;
        double g = compute_g_Rings(t);
        double score = metric(g);
        if (score > best_score) best_score = score;
    }
    return best_score;
}

// --- Parameter Exploration (1 method) ---
std::vector<std::map<std::string, double>> RingsOfRelativity::generateVariations(int count, double variation_pct) {
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
void RingsOfRelativity::mutateParameters(double mutation_rate, std::mt19937& rng) {
    std::uniform_real_distribution<> dis(-mutation_rate, mutation_rate);
    auto vars = listVariables();
    
    for (const auto& v : vars) {
        if (v == "c_light" || v == "G" || v == "hbar") continue;
        double val = getVariable(v);
        double delta = val * dis(rng);
        setVariable(v, val + delta);
    }
}

void RingsOfRelativity::evolveSystem(int generations, std::function<double(const RingsOfRelativity&)> fitness) {
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
bool RingsOfRelativity::saveState(const std::string& stateName) {
    std::map<std::string, double> state;
    auto vars = listVariables();
    for (const auto& v : vars) {
        state[v] = getVariable(v);
    }
    rings_saved_states[stateName] = state;
    return true;
}

bool RingsOfRelativity::restoreState(const std::string& stateName) {
    auto it = rings_saved_states.find(stateName);
    if (it == rings_saved_states.end()) return false;
    
    for (const auto& pair : it->second) {
        setVariable(pair.first, pair.second);
    }
    return true;
}

std::vector<std::string> RingsOfRelativity::listSavedStates() const {
    std::vector<std::string> names;
    for (const auto& pair : rings_saved_states) {
        names.push_back(pair.first);
    }
    return names;
}

std::string RingsOfRelativity::exportState() const {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(6);
    oss << "RingsOfRelativity State Export:\n";
    auto vars = listVariables();
    for (const auto& v : vars) {
        oss << v << " = " << getVariable(v) << "\n";
    }
    return oss.str();
}

// --- System Analysis (4 methods) ---
std::map<std::string, double> RingsOfRelativity::sensitivityAnalysis(double t, double delta_pct) {
    std::map<std::string, double> sensitivities;
    double g_base = compute_g_Rings(t);
    
    auto vars = listVariables();
    for (const auto& v : vars) {
        if (v == "c_light" || v == "G" || v == "hbar") continue;
        
        double original = getVariable(v);
        double delta = original * delta_pct / 100.0;
        
        setVariable(v, original + delta);
        double g_plus = compute_g_Rings(t);
        setVariable(v, original);
        
        double sensitivity = (g_base != 0.0) ? std::abs((g_plus - g_base) / g_base) : 0.0;
        sensitivities[v] = sensitivity;
    }
    
    return sensitivities;
}

std::string RingsOfRelativity::generateReport(double t) const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6);
    oss << "============================================\n";
    oss << "RINGS OF RELATIVITY REPORT\n";
    oss << "GAL-CLUS-022058s Einstein Ring\n";
    oss << "============================================\n";
    oss << "Time: t = " << t << " s (" << (t/3.156e7/1e9) << " Gyr)\n\n";
    
    oss << "Physical Parameters:\n";
    double M_sun = 1.989e30;
    oss << "  Lensing Mass M = " << M << " kg (" << (M/M_sun) << " M_sun)\n";
    double kpc_to_m = 3.086e19;
    oss << "  Einstein Radius r = " << r << " m (" << (r/kpc_to_m) << " kpc)\n";
    oss << "  Lens Redshift z_lens = " << z_lens << "\n";
    oss << "  Hubble Parameter Hz = " << Hz << " s^-1\n";
    oss << "  Magnetic field B = " << B << " T (B_crit = " << B_crit << " T)\n";
    oss << "  Lensing Factor L_factor = " << L_factor << ", L_t = " << L_t << "\n";
    oss << "  f_TRZ = " << f_TRZ << "\n";
    oss << "  Wind density rho_wind = " << rho_wind << " kg/m^3, v_wind = " << v_wind << " m/s\n";
    oss << "  Fluid density rho_fluid = " << rho_fluid << " kg/m^3\n";
    oss << "  Gas velocity = " << gas_v << " m/s\n";
    oss << "  DM factor = " << M_DM_factor << "\n\n";
    
    oss << "Computed Acceleration:\n";
    oss << "  g_Rings(t) = " << compute_g_Rings(t) << " m/s^2\n\n";
    
    oss << "UQFF Terms:\n";
    double corr_H = 1 + Hz * t;
    double corr_B = 1 - B / B_crit;
    double corr_L = 1 + L_t;
    oss << "  Base (with Hz, B, L): " << (ug1_base * corr_H * corr_B * corr_L) << " m/s^2\n";
    oss << "  Ug total: " << compute_Ug(0) << " m/s^2\n";
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

bool RingsOfRelativity::validateConsistency() const {
    bool valid = true;
    
    if (M <= 0 || r <= 0) { std::cerr << "Error: M and r must be positive.\n"; valid = false; }
    if (B < 0 || B_crit <= 0) { std::cerr << "Error: B, B_crit must be non-negative/positive.\n"; valid = false; }
    if (Hz < 0) { std::cerr << "Error: Hz must be non-negative.\n"; valid = false; }
    if (z_lens < 0) { std::cerr << "Warning: Lens redshift z_lens is negative.\n"; }
    if (L_factor < 0 || L_factor > 1.0) { std::cerr << "Warning: Lensing factor L_factor outside [0,1].\n"; }
    if (rho_fluid <= 0 || rho_wind < 0) { std::cerr << "Error: Fluid/wind densities must be positive/non-negative.\n"; valid = false; }
    if (v_wind < 0 || gas_v < 0) { std::cerr << "Error: Velocities must be non-negative.\n"; valid = false; }
    if (M_DM_factor < 0 || M_DM_factor > 1.0) { std::cerr << "Warning: DM factor outside [0,1].\n"; }
    
    return valid;
}

bool RingsOfRelativity::autoCorrectAnomalies() {
    bool corrected = false;
    
    double M_sun = 1.989e30;
    
    if (M <= 0) { M = 1e14 * M_sun; corrected = true; }
    if (r <= 0) { r = 3.086e20; corrected = true; }
    if (B < 0) { B = 1e-5; corrected = true; }
    if (B_crit <= 0) { B_crit = 1e11; corrected = true; }
    if (Hz < 0) { Hz = 2.42e-18; corrected = true; }
    if (z_lens < 0) { z_lens = 0.5; corrected = true; }
    if (L_factor < 0) { L_factor = 0.0; corrected = true; }
    if (L_factor > 1.0) { L_factor = 1.0; corrected = true; }
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
void enhanced_rings_example() {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "=========================================================\n";
    std::cout << "ENHANCED RINGS OF RELATIVITY DEMONSTRATION\n";
    std::cout << "GAL-CLUS-022058s Einstein Ring (z=0.5)\n";
    std::cout << "=========================================================\n\n";
    
    RingsOfRelativity rings;
    
    // Step 1: Initial state and validation
    std::cout << "Step 1: Initial State and Validation\n";
    std::cout << "System: " << rings.getSystemName() << "\n";
    std::cout << "Validation: " << (rings.validateConsistency() ? "PASS" : "FAIL") << "\n";
    std::cout << "Auto-corrected: " << (rings.autoCorrectAnomalies() ? "Yes" : "No") << "\n\n";
    
    // Step 2: Time evolution showing gravitational lensing effects
    std::cout << "Step 2: Time Evolution (Einstein Ring gravitational lensing)\n";
    double t_Gyr_array[] = {0.0, 1.0, 2.0, 5.0, 10.0};
    for (double t_Gyr : t_Gyr_array) {
        double t = t_Gyr * 1e9 * 3.156e7;
        double g = rings.compute_g_Rings(t);
        std::cout << "  t = " << t_Gyr << " Gyr: g = " << g << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 3: Variable listing
    std::cout << "Step 3: Variable Listing\n";
    auto vars = rings.listVariables();
    std::cout << "Total variables: " << vars.size() << "\n";
    std::cout << "Sample: " << vars[0] << ", " << vars[1] << ", " << vars[11] << " (L_factor), " 
              << vars[16] << " (z_lens)\n\n";
    
    // Step 4: Lensing mass scaling
    std::cout << "Step 4: Lensing Mass Scaling (M sweeps)\n";
    rings.saveState("original");
    double M_factors[] = {0.5, 1.0, 2.0};
    for (double factor : M_factors) {
        rings.restoreState("original");
        rings.expandLensingScale(factor, 1.0);
        double t = 5e9 * 3.156e7;
        double g = rings.compute_g_Rings(t);
        double M_sun = 1.989e30;
        double M = rings.getVariable("M");
        std::cout << "  M × " << factor << ": M = " << (M/M_sun) << " M_sun, g(5 Gyr) = " << g << " m/s^2\n";
    }
    rings.restoreState("original");
    std::cout << "\n";
    
    // Step 5: Lensing factor scaling (UNIQUE to Einstein Ring)
    std::cout << "Step 5: Lensing Factor Scaling (L_factor sweeps) - EINSTEIN RING FEATURE\n";
    double L_factors[] = {0.3, 0.67, 1.0};
    for (double factor : L_factors) {
        rings.restoreState("original");
        rings.setVariable("L_factor", factor);
        double t = 5e9 * 3.156e7;
        double g = rings.compute_g_Rings(t);
        std::cout << "  L_factor = " << factor << ": g(5 Gyr) = " << g << " m/s^2\n";
    }
    rings.restoreState("original");
    std::cout << "\n";
    
    // Step 6: Redshift scaling (cosmological distance effects)
    std::cout << "Step 6: Redshift Scaling (z_lens and Hz sweeps)\n";
    double z_values[] = {0.25, 0.5, 1.0};
    for (double z : z_values) {
        rings.restoreState("original");
        double Hz_kms = 70 * sqrt(0.3 * pow(1 + z, 3) + 0.7);
        double Hz_new = (Hz_kms * 1000 / 3.086e19);
        rings.setVariable("z_lens", z);
        rings.setVariable("Hz", Hz_new);
        double t = 5e9 * 3.156e7;
        double g = rings.compute_g_Rings(t);
        std::cout << "  z_lens = " << z << ", Hz = " << Hz_new << ": g(5 Gyr) = " << g << " m/s^2\n";
    }
    rings.restoreState("original");
    std::cout << "\n";
    
    // Step 7: Magnetic field scaling
    std::cout << "Step 7: Magnetic Field Scaling\n";
    double B_factors[] = {0.5, 1.0, 2.0};
    for (double factor : B_factors) {
        rings.restoreState("original");
        rings.expandMagneticWindScale(factor, 1.0, 1.0);
        double t = 5e9 * 3.156e7;
        double g = rings.compute_g_Rings(t);
        std::cout << "  B × " << factor << ": g(5 Gyr) = " << g << " m/s^2\n";
    }
    rings.restoreState("original");
    std::cout << "\n";
    
    // Step 8: Wind velocity scaling
    std::cout << "Step 8: Wind Velocity Scaling\n";
    double v_wind_factors[] = {0.5, 1.0, 2.0};
    for (double factor : v_wind_factors) {
        rings.restoreState("original");
        rings.expandMagneticWindScale(1.0, 1.0, factor);
        double t = 5e9 * 3.156e7;
        double g = rings.compute_g_Rings(t);
        std::cout << "  v_wind × " << factor << ": g(5 Gyr) = " << g << " m/s^2\n";
    }
    rings.restoreState("original");
    std::cout << "\n";
    
    // Step 9: Parameter space expansion
    std::cout << "Step 9: Parameter Space Expansion (all scalable params)\n";
    rings.expandParameterSpace(1.2);
    double M_after = rings.getVariable("M");
    double M_sun = 1.989e30;
    std::cout << "  After 1.2× expansion: M = " << (M_after/M_sun) << " M_sun\n";
    rings.restoreState("original");
    std::cout << "\n";
    
    // Step 10: Batch operations
    std::cout << "Step 10: Batch Operations (scale multiple variables)\n";
    std::vector<std::string> scale_group = {"M", "r", "B"};
    rings.scaleVariableGroup(scale_group, 1.1);
    std::cout << "  Scaled {M, r, B} by 1.1×\n";
    rings.restoreState("original");
    std::cout << "\n";
    
    // Step 11: State management
    std::cout << "Step 11: State Management\n";
    rings.saveState("state_A");
    rings.expandLensingScale(1.5, 1.2);
    rings.saveState("state_B");
    auto states = rings.listSavedStates();
    std::cout << "  Saved states: ";
    for (const auto& s : states) std::cout << s << " ";
    std::cout << "\n";
    rings.restoreState("state_A");
    std::cout << "  Restored state_A\n\n";
    
    // Step 12: Generate parameter variations
    std::cout << "Step 12: Generate Parameter Variations (5% variation)\n";
    auto variations = rings.generateVariations(3, 5.0);
    std::cout << "  Generated " << variations.size() << " variants with 5% random variation\n";
    std::cout << "  Variant 1 M = " << variations[0]["M"] << " kg\n\n";
    
    // Step 13: Sensitivity analysis
    std::cout << "Step 13: Sensitivity Analysis at 5 Gyr\n";
    double t_sens = 5e9 * 3.156e7;
    auto sensitivities = rings.sensitivityAnalysis(t_sens, 1.0);
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
        double t_obs = i * 2e9 * 3.156e7;
        double g_obs = rings.compute_g_Rings(t_obs) * (1.0 + 0.01 * (rand() % 100 - 50) / 100.0);
        obs.push_back({t_obs, g_obs});
    }
    rings.autoRefineParameters(obs);
    std::cout << "  Refined parameters based on " << obs.size() << " observations\n\n";
    
    // Step 15: Optimization for maximum acceleration
    std::cout << "Step 15: Optimize for Maximum Acceleration\n";
    rings.restoreState("original");
    auto metric = [](double g) { return g; };
    double t_opt_start = 0.0;
    double t_opt_end = 10e9 * 3.156e7;
    double best_g = rings.optimizeForMetric(metric, t_opt_start, t_opt_end, 50);
    std::cout << "  Best g over 10 Gyr: " << best_g << " m/s^2\n\n";
    
    // Step 16: Evolutionary system adaptation
    std::cout << "Step 16: Evolutionary System Adaptation (5 generations)\n";
    rings.restoreState("original");
    auto fitness = [](const RingsOfRelativity& r) {
        double t = 5e9 * 3.156e7;
        return r.compute_g_Rings(t);
    };
    rings.evolveSystem(5, fitness);
    std::cout << "  Evolved system over 5 generations (fitness = g at 5 Gyr)\n\n";
    
    // Step 17: Full system report
    std::cout << "Step 17: Full System Report at 5 Gyr\n";
    rings.restoreState("original");
    double t_report = 5e9 * 3.156e7;
    std::string report = rings.generateReport(t_report);
    std::cout << report << "\n";
    
    // Step 18: Full state export
    std::cout << "Step 18: Full State Export\n";
    std::string exported = rings.exportState();
    std::cout << "Exported state (first 500 chars):\n";
    std::cout << exported.substr(0, 500) << "...\n\n";
    
    std::cout << "=========================================================\n";
    std::cout << "ENHANCED DEMONSTRATION COMPLETE\n";
    std::cout << "=========================================================\n";
}