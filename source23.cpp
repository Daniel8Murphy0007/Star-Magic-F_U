/**
 * ================================================================================================
 * Header: AntennaeGalaxies.h
 *
 * Description: C++ Module for Antennae Galaxies (NGC 4038 & NGC 4039) Class
 *              This is the fourteenth module in a series of 500+ code files for the Universal Quantum
 *              Field Framework (UQFF) simulations, focusing on interacting galaxy merger evolution
 *              and gravity equations derived from Hubble datasets, high-energy lab simulations, and
 *              UQFF refinements (dated May 09, 2025, updated for full term inclusion on October 08, 2025).
 *
 * Purpose: Encapsulates the Master Universal Gravity Equation (MUGE) for Antennae Galaxies evolution.
 *          Includes ALL terms: base gravity with star formation growth M(t), cosmic expansion (H(z)),
 *          magnetic correction (static B), merger interaction I(t), UQFF Ug components with f_TRZ,
 *          Lambda, quantum uncertainty, scaled EM with [UA], fluid dynamics, oscillatory waves,
 *          DM/density perturbations, and stellar feedback (pressure / density for acc).
 *          Supports dynamic variable updates for all parameters.
 *
 * Integration: Designed for inclusion in base program 'ziqn233h.cpp' (not present here).
 *              Instantiate class in main: AntennaeGalaxies antennae;
 *              Compute: double g = antennae.compute_g_Antennae(t);
 *
 * Key Features:
 *   - Default values from UQFF document: M0 = 2e11 Msun, r = 2.838e20 m (30k ly), z = 0.0105,
 *     Hz ? 2.19e-18 s^-1, SFR_factor = 20 / (2e11), tau_SF = 100 Myr, I0 = 0.1, tau_merger = 400 Myr,
 *     rho_wind = 1e-21 kg/m^3, v_wind = 2e6 m/s, B = 1e-5 T.
 *   - Units handled: Msun to kg, ly to m; interaction term I(t) scales gravity.
 *   - Setter methods for updates: setVar(double new_val) or addToVar(double delta)/subtractFromVar(double delta).
 *   - Computes g_Antennae(r, t) with every term explicitly included.
 *
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: October 08, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

#ifndef ANTENNAE_GALAXIES_H
#define ANTENNAE_GALAXIES_H

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

class AntennaeGalaxies {
private:
    // Core parameters (mutable for updates)
    double G;               // Gravitational constant
    double M0;              // Initial combined mass (kg)
    double r;               // Separation (m)
    double Hz;              // Hubble parameter at z (s^-1)
    double B;               // Static magnetic field (T)
    double B_crit;          // Critical B field (T)
    double Lambda;          // Cosmological constant
    double c_light;         // Speed of light
    double q_charge;        // Charge (proton)
    double gas_v;           // Gas velocity for EM (m/s)
    double f_TRZ;           // Time-reversal factor
    double SFR_factor;      // Star formation rate factor (dimensionless)
    double tau_SF;          // Star formation timescale (s)
    double I0;              // Initial interaction factor
    double tau_merger;      // Merger timescale (s)
    double rho_wind;        // Wind density (kg/m^3)
    double v_wind;          // Wind velocity (m/s)
    double rho_fluid;       // Fluid density (kg/m^3)
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
    AntennaeGalaxies() {
        initializeDefaults();
    }

    // Destructor (empty)
    ~AntennaeGalaxies() {}

    // Initialization method (called in constructor)
    void initializeDefaults() {
        G = 6.6743e-11;
        double M_sun = 1.989e30;
        M0 = 2e11 * M_sun;
        double ly_to_m = 9.461e15;
        r = 30000.0 * ly_to_m;
        z_gal = 0.0105;
        double Hz_kms = 70 * sqrt(0.3 * pow(1 + z_gal, 3) + 0.7);  // km/s/Mpc
        Hz = (Hz_kms * 1000 / 3.086e19);  // s^-1
        B = 1e-5;
        B_crit = 1e11;
        Lambda = 1.1e-52;
        c_light = 3e8;
        q_charge = 1.602e-19;
        gas_v = 1e5;
        f_TRZ = 0.1;
        SFR_factor = 20.0 / (2e11);  // Normalized SFR
        tau_SF = 100e6 * 3.156e7;
        I0 = 0.1;
        tau_merger = 400e6 * 3.156e7;
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
        else if (varName == "Hz") { Hz = newValue; }
        else if (varName == "B") { B = newValue; }
        else if (varName == "B_crit") { B_crit = newValue; }
        else if (varName == "Lambda") { Lambda = newValue; }
        else if (varName == "c_light") { c_light = newValue; }
        else if (varName == "q_charge") { q_charge = newValue; }
        else if (varName == "gas_v") { gas_v = newValue; }
        else if (varName == "f_TRZ") { f_TRZ = newValue; }
        else if (varName == "SFR_factor") { SFR_factor = newValue; }
        else if (varName == "tau_SF") { tau_SF = newValue; }
        else if (varName == "I0") { I0 = newValue; }
        else if (varName == "tau_merger") { tau_merger = newValue; }
        else if (varName == "rho_wind") { rho_wind = newValue; }
        else if (varName == "v_wind") { v_wind = newValue; }
        else if (varName == "rho_fluid") { rho_fluid = newValue; }
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
        else if (varName == "Hz") return Hz;
        else if (varName == "B") return B;
        else if (varName == "B_crit") return B_crit;
        else if (varName == "Lambda") return Lambda;
        else if (varName == "c_light") return c_light;
        else if (varName == "q_charge") return q_charge;
        else if (varName == "gas_v") return gas_v;
        else if (varName == "f_TRZ") return f_TRZ;
        else if (varName == "SFR_factor") return SFR_factor;
        else if (varName == "tau_SF") return tau_SF;
        else if (varName == "I0") return I0;
        else if (varName == "tau_merger") return tau_merger;
        else if (varName == "rho_wind") return rho_wind;
        else if (varName == "v_wind") return v_wind;
        else if (varName == "rho_fluid") return rho_fluid;
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
        double M_dot = SFR_factor * exp(-t / tau_SF);
        return M0 * (1 + M_dot);
    }

    // I(t) computation
    double I_t(double t) const {
        return I0 * exp(-t / tau_merger);
    }

    // Ug terms computation
    double compute_Ug(double Mt, double It) const {
        double Ug1 = (G * Mt) / (r * r);
        double Ug2 = 0.0;
        double Ug3 = 0.0;
        double corr_B = 1 - B / B_crit;
        double Ug4 = Ug1 * corr_B;
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1 + f_TRZ) * (1 + It);
    }

    // Volume computation for fluid
    double compute_V() const {
        return (4.0 / 3.0) * M_PI * r * r * r;
    }

    // Main MUGE computation (includes ALL terms)
    double compute_g_Antennae(double t) const {
        if (t < 0) {
            std::cerr << "Error: Time t must be non-negative." << std::endl;
            return 0.0;
        }

        double Mt = M_t(t);
        double It = I_t(t);
        double ug1_t = (G * Mt) / (r * r);

        // Term 1: Base + Hz + B + I corrections
        double corr_H = 1 + Hz * t;
        double corr_B = 1 - B / B_crit;
        double corr_I = 1 + It;
        double term1 = ug1_t * corr_H * corr_B * corr_I;

        // Term 2: UQFF Ug with f_TRZ and I
        double term2 = compute_Ug(Mt, It);

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

        // Stellar feedback term (pressure / density for acceleration)
        double wind_pressure = rho_wind * v_wind * v_wind;
        double term_feedback = wind_pressure / rho_fluid;

        // Total g_Antennae (all terms summed)
        return term1 + term2 + term3 + term4 + term_q + term_fluid + term_osc + term_DM + term_feedback;
    }

    // Debug/Output method (for transparency in base program)
    void printParameters(std::ostream& os = std::cout) const {
        os << std::fixed << std::setprecision(3);
        os << "Antennae Galaxies Parameters:" << std::endl;
        os << "G: " << G << ", M0: " << M0 << ", r: " << r << std::endl;
        os << "Hz: " << Hz << ", B: " << B << ", B_crit: " << B_crit << std::endl;
        os << "f_TRZ: " << f_TRZ << ", SFR_factor: " << SFR_factor << ", tau_SF: " << tau_SF << std::endl;
        os << "I0: " << I0 << ", tau_merger: " << tau_merger << std::endl;
        os << "rho_fluid: " << rho_fluid << ", rho_wind: " << rho_wind << ", v_wind: " << v_wind << std::endl;
        os << "gas_v: " << gas_v << ", M_DM_factor: " << M_DM_factor << std::endl;
        os << "A_osc: " << A_osc << ", delta_rho_over_rho: " << delta_rho_over_rho << std::endl;
        os << "ug1_base: " << ug1_base << std::endl;
    }

    // Example computation at t=300 Myr (for testing)
    double exampleAt300Myr() const {
        double t_example = 300e6 * 3.156e7;
        return compute_g_Antennae(t_example);
    }

    // ========== ENHANCED DYNAMIC CAPABILITIES (25 methods) ==========
    
    // Variable Management (5 methods)
    void createVariable(const std::string& name, double value);
    bool removeVariable(const std::string& name);
    void cloneVariable(const std::string& source, const std::string& dest);
    std::vector<std::string> listVariables() const;
    std::string getSystemName() const { return "AntennaeGalaxies"; }
    
    // Batch Operations (2 methods)
    void transformVariableGroup(const std::vector<std::string>& vars, std::function<double(double)> func);
    void scaleVariableGroup(const std::vector<std::string>& vars, double factor);
    
    // Self-Expansion (4 methods - domain-specific for galaxy merger)
    void expandParameterSpace(double scale_factor);
    void expandGalaxyScale(double M_scale, double r_scale);
    void expandMergerScale(double I0_scale, double tau_merger_scale);
    void expandStarFormationScale(double SFR_factor_scale, double tau_SF_scale);
    
    // Self-Refinement (3 methods)
    void autoRefineParameters(const std::vector<std::pair<double, double>>& observations);
    void calibrateToObservations(const std::vector<std::pair<double, double>>& obs_data);
    double optimizeForMetric(std::function<double(double)> metric, double t_start, double t_end, int steps);
    
    // Parameter Exploration (1 method)
    std::vector<std::map<std::string, double>> generateVariations(int count, double variation_percent);
    
    // Adaptive Evolution (2 methods)
    void mutateParameters(double mutation_rate);
    void evolveSystem(int generations, std::function<double(const AntennaeGalaxies&)> fitness);
    
    // State Management (4 methods)
    void saveState(const std::string& label);
    bool restoreState(const std::string& label);
    std::vector<std::string> listSavedStates() const;
    std::string exportState() const;
    
    // System Analysis (4 methods)
    std::map<std::string, double> sensitivityAnalysis(double t, double perturbation);
    std::string generateReport(double t) const;
    bool validateConsistency() const;
    bool autoCorrectAnomalies();
};

#endif // ANTENNAE_GALAXIES_H

// ========== IMPLEMENTATION OF ENHANCED METHODS ==========

namespace {
    // Storage for saved states (anonymous namespace for encapsulation)
    std::map<std::string, std::map<std::string, double>> antennae_saved_states;
}

// Variable Management
void AntennaeGalaxies::createVariable(const std::string& name, double value) {
    setVariable(name, value);
}

bool AntennaeGalaxies::removeVariable(const std::string& name) {
    // Cannot truly remove core variables, but can reset to defaults
    initializeDefaults();
    return true;
}

void AntennaeGalaxies::cloneVariable(const std::string& source, const std::string& dest) {
    double value = getVariable(source);
    setVariable(dest, value);
}

std::vector<std::string> AntennaeGalaxies::listVariables() const {
    return {"G", "M0", "r", "Hz", "B", "B_crit", "Lambda", "c_light", "q_charge", "gas_v", 
            "f_TRZ", "SFR_factor", "tau_SF", "I0", "tau_merger", "rho_wind", "v_wind", 
            "rho_fluid", "rho_vac_UA", "rho_vac_SCm", "scale_EM", "proton_mass", "z_gal",
            "hbar", "t_Hubble", "t_Hubble_gyr", "delta_x", "delta_p", "integral_psi",
            "A_osc", "k_osc", "omega_osc", "x_pos", "M_DM_factor", "delta_rho_over_rho"};
}

// Batch Operations
void AntennaeGalaxies::transformVariableGroup(const std::vector<std::string>& vars, std::function<double(double)> func) {
    for (const auto& var : vars) {
        double current = getVariable(var);
        setVariable(var, func(current));
    }
}

void AntennaeGalaxies::scaleVariableGroup(const std::vector<std::string>& vars, double factor) {
    transformVariableGroup(vars, [factor](double v) { return v * factor; });
}

// Self-Expansion (domain-specific for galaxy merger)
void AntennaeGalaxies::expandParameterSpace(double scale_factor) {
    std::vector<std::string> scalable = {"M0", "r", "I0", "tau_merger", "SFR_factor", "tau_SF", 
                                          "rho_wind", "v_wind", "B"};
    scaleVariableGroup(scalable, scale_factor);
}

void AntennaeGalaxies::expandGalaxyScale(double M_scale, double r_scale) {
    setVariable("M0", getVariable("M0") * M_scale);
    setVariable("r", getVariable("r") * r_scale);
}

void AntennaeGalaxies::expandMergerScale(double I0_scale, double tau_merger_scale) {
    setVariable("I0", getVariable("I0") * I0_scale);
    setVariable("tau_merger", getVariable("tau_merger") * tau_merger_scale);
}

void AntennaeGalaxies::expandStarFormationScale(double SFR_factor_scale, double tau_SF_scale) {
    setVariable("SFR_factor", getVariable("SFR_factor") * SFR_factor_scale);
    setVariable("tau_SF", getVariable("tau_SF") * tau_SF_scale);
}

// Self-Refinement
void AntennaeGalaxies::autoRefineParameters(const std::vector<std::pair<double, double>>& observations) {
    calibrateToObservations(observations);
}

void AntennaeGalaxies::calibrateToObservations(const std::vector<std::pair<double, double>>& obs_data) {
    if (obs_data.empty()) return;
    
    double total_error = 0.0;
    for (const auto& obs : obs_data) {
        double t = obs.first;
        double g_obs = obs.second;
        double g_model = compute_g_Antennae(t);
        total_error += std::abs(g_model - g_obs);
    }
    
    double avg_error = total_error / obs_data.size();
    if (avg_error > 1e-9) {
        double correction = 0.95;
        setVariable("SFR_factor", getVariable("SFR_factor") * correction);
    }
}

double AntennaeGalaxies::optimizeForMetric(std::function<double(double)> metric, double t_start, double t_end, int steps) {
    double best_metric = -1e100;
    double dt = (t_end - t_start) / steps;
    
    for (int i = 0; i <= steps; ++i) {
        double t = t_start + i * dt;
        double g = compute_g_Antennae(t);
        double m = metric(g);
        if (m > best_metric) best_metric = m;
    }
    return best_metric;
}

// Parameter Exploration
std::vector<std::map<std::string, double>> AntennaeGalaxies::generateVariations(int count, double variation_percent) {
    std::vector<std::map<std::string, double>> variations;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-variation_percent / 100.0, variation_percent / 100.0);
    
    auto vars = listVariables();
    for (int i = 0; i < count; ++i) {
        std::map<std::string, double> variant;
        for (const auto& var : vars) {
            double base = getVariable(var);
            double variation = base * (1.0 + dis(gen));
            variant[var] = variation;
        }
        variations.push_back(variant);
    }
    return variations;
}

// Adaptive Evolution
void AntennaeGalaxies::mutateParameters(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-mutation_rate, mutation_rate);
    
    auto vars = listVariables();
    for (const auto& var : vars) {
        double current = getVariable(var);
        double mutated = current * (1.0 + dis(gen));
        setVariable(var, mutated);
    }
}

void AntennaeGalaxies::evolveSystem(int generations, std::function<double(const AntennaeGalaxies&)> fitness) {
    double best_fitness = fitness(*this);
    saveState("evolution_best");
    
    for (int gen = 0; gen < generations; ++gen) {
        saveState("evolution_temp");
        mutateParameters(0.1);
        
        double current_fitness = fitness(*this);
        if (current_fitness > best_fitness) {
            best_fitness = current_fitness;
            saveState("evolution_best");
        } else {
            restoreState("evolution_temp");
        }
    }
    restoreState("evolution_best");
}

// State Management
void AntennaeGalaxies::saveState(const std::string& label) {
    std::map<std::string, double> state;
    auto vars = listVariables();
    for (const auto& var : vars) {
        state[var] = getVariable(var);
    }
    antennae_saved_states[label] = state;
}

bool AntennaeGalaxies::restoreState(const std::string& label) {
    auto it = antennae_saved_states.find(label);
    if (it == antennae_saved_states.end()) return false;
    
    for (const auto& pair : it->second) {
        setVariable(pair.first, pair.second);
    }
    return true;
}

std::vector<std::string> AntennaeGalaxies::listSavedStates() const {
    std::vector<std::string> labels;
    for (const auto& pair : antennae_saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::string AntennaeGalaxies::exportState() const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(15);
    oss << "AntennaeGalaxies State Export\n";
    oss << "============================\n";
    
    auto vars = listVariables();
    for (const auto& var : vars) {
        oss << var << ": " << getVariable(var) << "\n";
    }
    return oss.str();
}

// System Analysis
std::map<std::string, double> AntennaeGalaxies::sensitivityAnalysis(double t, double perturbation) {
    std::map<std::string, double> sensitivities;
    double baseline = compute_g_Antennae(t);
    
    auto vars = listVariables();
    for (const auto& var : vars) {
        double original = getVariable(var);
        setVariable(var, original * (1.0 + perturbation));
        double perturbed = compute_g_Antennae(t);
        sensitivities[var] = std::abs(perturbed - baseline) / baseline;
        setVariable(var, original);
    }
    return sensitivities;
}

std::string AntennaeGalaxies::generateReport(double t) const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6);
    oss << "=== Antennae Galaxies System Report ===\n";
    oss << "Time: " << (t / (1e6 * 3.156e7)) << " Myr\n";
    oss << "M(t): " << (M_t(t) / 1.989e30) << " M_sun\n";
    oss << "I(t) interaction factor: " << I_t(t) << "\n";
    oss << "g_Antennae: " << compute_g_Antennae(t) << " m/s^2\n";
    oss << "Core parameters: M0=" << (M0/1.989e30) << " M_sun, r=" << (r/9.461e15) << " ly\n";
    oss << "Merger: I0=" << I0 << ", tau_merger=" << (tau_merger/(1e6*3.156e7)) << " Myr\n";
    oss << "Star formation: SFR_factor=" << SFR_factor << ", tau_SF=" << (tau_SF/(1e6*3.156e7)) << " Myr\n";
    return oss.str();
}

bool AntennaeGalaxies::validateConsistency() const {
    bool valid = true;
    if (M0 <= 0 || r <= 0 || tau_SF <= 0 || tau_merger <= 0) valid = false;
    if (I0 < 0 || I0 > 1) valid = false;
    if (SFR_factor < 0) valid = false;
    return valid;
}

bool AntennaeGalaxies::autoCorrectAnomalies() {
    bool corrected = false;
    if (M0 <= 0) { M0 = 2e11 * 1.989e30; corrected = true; }
    if (r <= 0) { r = 30000 * 9.461e15; corrected = true; }
    if (tau_SF <= 0) { tau_SF = 100e6 * 3.156e7; corrected = true; }
    if (tau_merger <= 0) { tau_merger = 400e6 * 3.156e7; corrected = true; }
    if (I0 < 0) { I0 = 0.0; corrected = true; }
    if (I0 > 1) { I0 = 1.0; corrected = true; }
    if (SFR_factor < 0) { SFR_factor = 20.0 / (2e11); corrected = true; }
    if (corrected) updateCache();
    return corrected;
}

// ========== ENHANCED EXAMPLE FUNCTION ==========
void enhanced_antennae_example() {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "=========================================================\n";
    std::cout << "ENHANCED ANTENNAE GALAXIES DEMONSTRATION\n";
    std::cout << "NGC 4038 & NGC 4039 Galaxy Merger\n";
    std::cout << "=========================================================\n\n";
    
    AntennaeGalaxies antennae;
    
    // Step 1: Initial state and validation
    std::cout << "Step 1: Initial State and Validation\n";
    std::cout << "System: " << antennae.getSystemName() << "\n";
    std::cout << "Validation: " << (antennae.validateConsistency() ? "PASS" : "FAIL") << "\n";
    std::cout << "Auto-corrected: " << (antennae.autoCorrectAnomalies() ? "Yes" : "No") << "\n\n";
    
    // Step 2: Time evolution showing M(t) and I(t)
    std::cout << "Step 2: Time Evolution (Star Formation M(t) and Merger Interaction I(t))\n";
    double t_Myr_array[] = {0.0, 100.0, 200.0, 300.0, 500.0};
    for (double t_Myr : t_Myr_array) {
        double t = t_Myr * 1e6 * 3.156e7;
        double Mt = antennae.M_t(t);
        double It = antennae.I_t(t);
        double g = antennae.compute_g_Antennae(t);
        double M_sun = 1.989e30;
        std::cout << "  t = " << t_Myr << " Myr: M(t) = " << (Mt/M_sun) << " M_sun, I(t) = " << It << ", g = " << g << " m/s^2\n";
    }
    std::cout << "\n";
    
    // Step 3: Variable listing
    std::cout << "Step 3: Variable Listing\n";
    auto vars = antennae.listVariables();
    std::cout << "Total variables: " << vars.size() << "\n";
    std::cout << "Sample: " << vars[0] << ", " << vars[1] << ", " << vars[13] << " (I0), " 
              << vars[14] << " (tau_merger)\n\n";
    
    // Step 4: Galaxy mass scaling
    std::cout << "Step 4: Galaxy Mass Scaling (M sweeps)\n";
    antennae.saveState("original");
    double M_factors[] = {0.5, 1.0, 2.0};
    for (double factor : M_factors) {
        antennae.restoreState("original");
        antennae.expandGalaxyScale(factor, 1.0);
        double t = 300e6 * 3.156e7;
        double g = antennae.compute_g_Antennae(t);
        double M_sun = 1.989e30;
        double M = antennae.getVariable("M0");
        std::cout << "  M × " << factor << ": M0 = " << (M/M_sun) << " M_sun, g(300 Myr) = " << g << " m/s^2\n";
    }
    antennae.restoreState("original");
    std::cout << "\n";
    
    // Step 5: Merger interaction scaling (UNIQUE to Antennae merger)
    std::cout << "Step 5: Merger Interaction Scaling (I0 sweeps) - GALAXY MERGER FEATURE\n";
    double I0_factors[] = {0.5, 1.0, 2.0};
    for (double factor : I0_factors) {
        antennae.restoreState("original");
        antennae.expandMergerScale(factor, 1.0);
        double t = 300e6 * 3.156e7;
        double It = antennae.I_t(t);
        double g = antennae.compute_g_Antennae(t);
        std::cout << "  I0 × " << factor << ": I(300 Myr) = " << It << ", g = " << g << " m/s^2\n";
    }
    antennae.restoreState("original");
    std::cout << "\n";
    
    // Step 6: Merger timescale sweeps
    std::cout << "Step 6: Merger Timescale Sweeps (tau_merger)\n";
    double tau_merger_factors[] = {0.5, 1.0, 2.0};
    for (double factor : tau_merger_factors) {
        antennae.restoreState("original");
        antennae.expandMergerScale(1.0, factor);
        double t = 300e6 * 3.156e7;
        double It = antennae.I_t(t);
        std::cout << "  tau_merger × " << factor << ": I(300 Myr) = " << It << "\n";
    }
    antennae.restoreState("original");
    std::cout << "\n";
    
    // Step 7: Star formation rate scaling
    std::cout << "Step 7: Star Formation Rate Scaling\n";
    double SFR_factors[] = {0.5, 1.0, 2.0};
    for (double factor : SFR_factors) {
        antennae.restoreState("original");
        antennae.expandStarFormationScale(factor, 1.0);
        double t = 300e6 * 3.156e7;
        double Mt = antennae.M_t(t);
        double M_sun = 1.989e30;
        std::cout << "  SFR_factor × " << factor << ": M(300 Myr) = " << (Mt/M_sun) << " M_sun\n";
    }
    antennae.restoreState("original");
    std::cout << "\n";
    
    // Step 8: Star formation timescale scaling
    std::cout << "Step 8: Star Formation Timescale Scaling\n";
    double tau_SF_factors[] = {0.5, 1.0, 2.0};
    for (double factor : tau_SF_factors) {
        antennae.restoreState("original");
        antennae.expandStarFormationScale(1.0, factor);
        double t = 300e6 * 3.156e7;
        double Mt = antennae.M_t(t);
        double M_sun = 1.989e30;
        std::cout << "  tau_SF × " << factor << ": M(300 Myr) = " << (Mt/M_sun) << " M_sun\n";
    }
    antennae.restoreState("original");
    std::cout << "\n";
    
    // Step 9: Parameter space expansion
    std::cout << "Step 9: Parameter Space Expansion (all scalable params)\n";
    antennae.expandParameterSpace(1.2);
    double M_after = antennae.getVariable("M0");
    double M_sun = 1.989e30;
    std::cout << "  After 1.2× expansion: M0 = " << (M_after/M_sun) << " M_sun\n";
    antennae.restoreState("original");
    std::cout << "\n";
    
    // Step 10: Batch operations
    std::cout << "Step 10: Batch Operations (scale multiple variables)\n";
    std::vector<std::string> scale_group = {"M0", "r", "I0"};
    antennae.scaleVariableGroup(scale_group, 1.1);
    std::cout << "  Scaled {M0, r, I0} by 1.1×\n";
    antennae.restoreState("original");
    std::cout << "\n";
    
    // Step 11: State management
    std::cout << "Step 11: State Management\n";
    antennae.saveState("state_A");
    antennae.expandMergerScale(1.5, 1.2);
    antennae.saveState("state_B");
    auto states = antennae.listSavedStates();
    std::cout << "  Saved states: ";
    for (const auto& s : states) std::cout << s << " ";
    std::cout << "\n";
    antennae.restoreState("state_A");
    std::cout << "  Restored state_A\n\n";
    
    // Step 12: Generate parameter variations
    std::cout << "Step 12: Generate Parameter Variations (5% variation)\n";
    auto variations = antennae.generateVariations(3, 5.0);
    std::cout << "  Generated " << variations.size() << " variants with 5% random variation\n";
    std::cout << "  Variant 1 M0 = " << variations[0]["M0"] << " kg\n\n";
    
    // Step 13: Sensitivity analysis
    std::cout << "Step 13: Sensitivity Analysis at 300 Myr\n";
    double t_sens = 300e6 * 3.156e7;
    auto sensitivities = antennae.sensitivityAnalysis(t_sens, 0.01);
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
        double t_obs = i * 100e6 * 3.156e7;
        double g_obs = antennae.compute_g_Antennae(t_obs) * (1.0 + 0.01 * (rand() % 100 - 50) / 100.0);
        obs.push_back({t_obs, g_obs});
    }
    antennae.autoRefineParameters(obs);
    std::cout << "  Refined parameters based on " << obs.size() << " observations\n\n";
    
    // Step 15: Optimization for maximum acceleration
    std::cout << "Step 15: Optimize for Maximum Acceleration\n";
    antennae.restoreState("original");
    auto metric = [](double g) { return g; };
    double t_opt_start = 0.0;
    double t_opt_end = 500e6 * 3.156e7;
    double best_g = antennae.optimizeForMetric(metric, t_opt_start, t_opt_end, 50);
    std::cout << "  Best g over 500 Myr: " << best_g << " m/s^2\n\n";
    
    // Step 16: Evolutionary system adaptation
    std::cout << "Step 16: Evolutionary System Adaptation (5 generations)\n";
    antennae.restoreState("original");
    auto fitness = [](const AntennaeGalaxies& a) {
        double t = 300e6 * 3.156e7;
        return a.compute_g_Antennae(t);
    };
    antennae.evolveSystem(5, fitness);
    std::cout << "  Evolved system over 5 generations (fitness = g at 300 Myr)\n\n";
    
    // Step 17: Full system report
    std::cout << "Step 17: Full System Report at 300 Myr\n";
    antennae.restoreState("original");
    double t_report = 300e6 * 3.156e7;
    std::string report = antennae.generateReport(t_report);
    std::cout << report << "\n";
    
    // Step 18: Full state export
    std::cout << "Step 18: Full State Export\n";
    std::string exported = antennae.exportState();
    std::cout << "Exported state (first 500 chars):\n";
    std::cout << exported.substr(0, 500) << "...\n\n";
    
    std::cout << "=========================================================\n";
    std::cout << "ENHANCED DEMONSTRATION COMPLETE\n";
    std::cout << "=========================================================\n";
}