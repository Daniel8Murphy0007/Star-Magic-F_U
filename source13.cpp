/**
 * ================================================================================================
 * Header: MagnetarSGR1745_2900.h
 *
 * Description: C++ Module for SGR 1745-2900 Magnetar Class
 *              This is the second module in a series of 500+ code files for the Universal Quantum
 *              Field Framework (UQFF) simulations, focusing on magnetar evolution and gravity
 *              equations derived from Chandra X-ray Observatory datasets, high-energy lab simulations,
 *              and UQFF refinements (dated May 11, 2025, updated for full term inclusion on October 08, 2025).
 *
 * Purpose: Encapsulates the Master Universal Gravity Equation (MUGE) for SGR 1745-2900 magnetar
 *          evolution, including black hole proximity (Sgr A*), magnetic energy, and outburst decay.
 *          Includes ALL terms: base gravity, cosmic expansion (H(z)), magnetic decay, BH influence,
 *          UQFF Ug components with f_sc, Lambda, quantum uncertainty, EM, fluid, oscillatory waves,
 *          DM/density perturbations, magnetic energy (effective g), and decay power (cumulative energy effective g).
 *          Supports dynamic variable updates for all parameters.
 *
 * Integration: Designed for inclusion in base program 'ziqn233h.cpp' (not present here).
 *              Instantiate class in main: MagnetarSGR1745_2900 mag;
 *              Compute: double g = mag.compute_g_Magnetar(t);
 *
 * Key Features:
 *   - Default values from UQFF document, with numerical computations for H(z), M_mag, etc.
 *   - Units handled: Energy terms converted to effective acceleration via / (M * r); decay as cumulative energy / (M * r).
 *   - Setter methods for updates: setVar(double new_val) or addToVar(double delta)/subtractFromVar(double delta).
 *   - Computes g_Magnetar(r, t) with every term explicitly included.
 *
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: October 08, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

#ifndef MAGNETAR_SGR1745_2900_H
#define MAGNETAR_SGR1745_2900_H

#include <iostream>
#include <cmath>
#include <iomanip>

class MagnetarSGR1745_2900 {
private:
    // Core parameters (mutable for updates)
    double G;               // Gravitational constant
    double M;               // Magnetar mass
    double r;               // Radius
    double Hz;              // Hubble parameter at z (s^-1)
    double B0;              // Initial magnetic field
    double tau_B;           // B decay timescale (s) - not used in this eq, but for consistency
    double B_crit;          // Critical B field
    double Lambda;          // Cosmological constant
    double c_light;         // Speed of light
    double q_charge;        // Charge (proton)
    double v_surf;          // Surface velocity
    double f_sc;            // Superconductive factor (computed as 1 - B/B_crit)
    double rho_vac_UA;      // UA vacuum density - not used here
    double rho_vac_SCm;     // SCm vacuum density - not used here
    double P_init;          // Initial rotation period (s)
    double tau_Omega;       // Omega decay timescale (s)
    double scale_EM;        // EM scaling factor
    double proton_mass;     // Proton mass for EM acceleration
    double M_BH;            // Black hole mass
    double r_BH;            // Distance to black hole
    double mu0;             // Vacuum permeability
    double L0_W;            // Initial luminosity (W)
    double tau_decay;       // Decay timescale (s)

    // Additional parameters for full inclusion of terms
    double hbar;            // Reduced Planck's constant
    double t_Hubble;        // Hubble time (s)
    double t_Hubble_gyr;    // Hubble time in Gyr
    double delta_x;         // Position uncertainty (m)
    double delta_p;         // Momentum uncertainty (kg m/s)
    double integral_psi;    // Wavefunction integral approximation
    double rho_fluid;       // Fluid density (kg/m^3)
    double A_osc;           // Oscillatory amplitude (m/s^2)
    double k_osc;           // Wave number (1/m)
    double omega_osc;       // Angular frequency (rad/s)
    double x_pos;           // Position for oscillation (m)
    double M_DM_factor;     // Dark matter mass fraction
    double delta_rho_over_rho; // Density perturbation fraction

    // Computed caches (updated on demand)
    double ug1_base;        // Cached Ug1 = G*M/r^2
    double B;               // Current B (static in this model)

public:
    // Constructor with default UQFF values
    MagnetarSGR1745_2900() {
        initializeDefaults();
    }

    // Destructor (empty)
    ~MagnetarSGR1745_2900() {}

    // Initialization method (called in constructor)
    void initializeDefaults() {
        G = 6.6743e-11;
        M = 1.4 * 1.989e30;
        r = 1e4;
        Hz = 2.269e-18;  // Computed H(z)
        B0 = 2e10;
        B = B0;  // Static for this model
        tau_B = 4000 * 3.15576e7;  // Default, not used
        B_crit = 1e11;
        Lambda = 1.1e-52;
        c_light = 3e8;
        q_charge = 1.602e-19;
        v_surf = 1e6;
        f_sc = 1 - (B / B_crit);  // Initial
        rho_vac_UA = 7.09e-36;
        rho_vac_SCm = 7.09e-37;
        P_init = 3.76;  // Pulse period
        tau_Omega = 10000 * 3.15576e7;
        scale_EM = 1e-12;
        proton_mass = 1.673e-27;
        M_BH = 4e6 * 1.989e30;
        r_BH = 2.83e16;
        mu0 = 4 * M_PI * 1e-7;
        L0_W = 5e28;  // 5e35 erg/s = 5e28 W
        tau_decay = 3.5 * 365.25 * 24 * 3600;  // 3.5 years in s

        // Full terms defaults
        hbar = 1.0546e-34;
        t_Hubble = 13.8e9 * 3.15576e7;
        t_Hubble_gyr = 13.8;
        delta_x = 1e-10;
        delta_p = hbar / delta_x;
        integral_psi = 1.0;
        rho_fluid = 1e17;
        A_osc = 1e10;
        k_osc = 1.0 / r;
        omega_osc = 2 * M_PI / P_init;
        x_pos = r;
        M_DM_factor = 0.1;
        delta_rho_over_rho = 1e-5;

        updateCache();
    }

    // Cache update for efficiency (call after parameter changes)
    void updateCache() {
        ug1_base = (G * M) / (r * r);
        f_sc = 1 - (B / B_crit);
    }

    // Universal setter for any variable (by name, for flexibility)
    bool setVariable(const std::string& varName, double newValue) {
        if (varName == "G") { G = newValue; }
        else if (varName == "M") { M = newValue; }
        else if (varName == "r") { r = newValue; }
        else if (varName == "Hz") { Hz = newValue; }
        else if (varName == "B0") { B0 = newValue; B = newValue; }
        else if (varName == "tau_B") { tau_B = newValue; }
        else if (varName == "B_crit") { B_crit = newValue; }
        else if (varName == "Lambda") { Lambda = newValue; }
        else if (varName == "c_light") { c_light = newValue; }
        else if (varName == "q_charge") { q_charge = newValue; }
        else if (varName == "v_surf") { v_surf = newValue; }
        else if (varName == "f_sc") { f_sc = newValue; }
        else if (varName == "rho_vac_UA") { rho_vac_UA = newValue; }
        else if (varName == "rho_vac_SCm") { rho_vac_SCm = newValue; }
        else if (varName == "P_init") { P_init = newValue; }
        else if (varName == "tau_Omega") { tau_Omega = newValue; }
        else if (varName == "scale_EM") { scale_EM = newValue; }
        else if (varName == "proton_mass") { proton_mass = newValue; }
        else if (varName == "M_BH") { M_BH = newValue; }
        else if (varName == "r_BH") { r_BH = newValue; }
        else if (varName == "mu0") { mu0 = newValue; }
        else if (varName == "L0_W") { L0_W = newValue; }
        else if (varName == "tau_decay") { tau_decay = newValue; }
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
        else if (varName == "B0") return B0;
        else if (varName == "tau_B") return tau_B;
        else if (varName == "B_crit") return B_crit;
        else if (varName == "Lambda") return Lambda;
        else if (varName == "c_light") return c_light;
        else if (varName == "q_charge") return q_charge;
        else if (varName == "v_surf") return v_surf;
        else if (varName == "f_sc") return f_sc;
        else if (varName == "rho_vac_UA") return rho_vac_UA;
        else if (varName == "rho_vac_SCm") return rho_vac_SCm;
        else if (varName == "P_init") return P_init;
        else if (varName == "tau_Omega") return tau_Omega;
        else if (varName == "scale_EM") return scale_EM;
        else if (varName == "proton_mass") return proton_mass;
        else if (varName == "M_BH") return M_BH;
        else if (varName == "r_BH") return r_BH;
        else if (varName == "mu0") return mu0;
        else if (varName == "L0_W") return L0_W;
        else if (varName == "tau_decay") return tau_decay;
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

    // B(t) - static B for this model
    double B_t(double /*t*/) const {
        return B;
    }

    // Omega(t) computation
    double Omega_t(double t) const {
        return (2 * M_PI / P_init) * exp(-t / tau_Omega);
    }

    // dOmega/dt computation
    double dOmega_dt(double t) const {
        double omega0 = 2 * M_PI / P_init;
        return omega0 * (-1.0 / tau_Omega) * exp(-t / tau_Omega);
    }

    // Ug terms computation
    double compute_Ug() const {
        double Ug1 = ug1_base;
        double Ug2 = 0.0;
        double Ug3 = 0.0;
        double Ug4 = Ug1 * f_sc;
        return Ug1 + Ug2 + Ug3 + Ug4;
    }

    // Volume computation
    double compute_V() const {
        return (4.0 / 3.0) * M_PI * r * r * r;
    }

    // Magnetic energy M_mag (J)
    double compute_M_mag() const {
        double V = compute_V();
        return (B_t(0) * B_t(0) / (2 * mu0)) * V;
    }

    // Cumulative decay energy up to t (J)
    double compute_cumulative_D(double t) const {
        double exp_term = exp(-t / tau_decay);
        return L0_W * tau_decay * (1 - exp_term);
    }

    // Main MUGE computation (includes ALL terms)
    double compute_g_Magnetar(double t) const {
        if (t < 0) {
            std::cerr << "Error: Time t must be non-negative." << std::endl;
            return 0.0;
        }

        double Bt = B_t(t);
        double dOdt = dOmega_dt(t);

        // f_sc update
        double current_f_sc = 1 - (Bt / B_crit);

        // Term 1: Base + H(z) + B corrections
        double corr_H = 1 + Hz * t;
        double corr_B = current_f_sc;
        double term1 = ug1_base * corr_H * corr_B;

        // BH term
        double term_BH = (G * M_BH) / (r_BH * r_BH);

        // Term 2: UQFF Ug 
        double term2 = compute_Ug();

        // Term 3: Lambda
        double term3 = (Lambda * c_light * c_light) / 3.0;

        // Term 4: Scaled EM (v x B magnitude)
        double cross_vB = v_surf * Bt;
        double em_base = (q_charge * cross_vB) / proton_mass;
        double term4 = em_base * scale_EM;  // UA not used here

        // Term 5: GW (assumed same as previous)
        double gw_prefactor = (G * M * M) / (pow(c_light, 4) * r);
        double term5 = gw_prefactor * (dOdt * dOdt);

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

        // Magnetic energy term (effective g)
        double M_mag = compute_M_mag();
        double term_mag = M_mag / (M * r);

        // Decay term (cumulative energy effective g)
        double cum_D = compute_cumulative_D(t);
        double term_decay = cum_D / (M * r);

        // Total g_Magnetar (all terms summed)
        return term1 + term_BH + term2 + term3 + term4 + term5 + term_q + term_fluid + term_osc + term_DM + term_mag + term_decay;
    }

    // Debug/Output method (for transparency in base program)
    void printParameters(std::ostream& os = std::cout) const {
        os << std::fixed << std::setprecision(3);
        os << "SGR 1745-2900 Parameters:" << std::endl;
        os << "G: " << G << ", M: " << M << ", r: " << r << std::endl;
        os << "Hz: " << Hz << ", B: " << B << ", M_BH: " << M_BH << ", r_BH: " << r_BH << std::endl;
        os << "L0_W: " << L0_W << ", tau_decay: " << tau_decay << std::endl;
        os << "f_sc: " << f_sc << ", rho_fluid: " << rho_fluid << ", M_DM_factor: " << M_DM_factor << std::endl;
        os << "A_osc: " << A_osc << ", delta_rho_over_rho: " << delta_rho_over_rho << std::endl;
        double M_mag = compute_M_mag();
        os << "M_mag (J): " << M_mag << ", ug1_base: " << ug1_base << std::endl;
    }

    // Example computation at t=1 year (for testing)
    double exampleAtOneYear() const {
        double t_example = 1.0 * 365.25 * 24 * 3600;
        return compute_g_Magnetar(t_example);
    }
};

#endif // MAGNETAR_SGR1745_2900_H