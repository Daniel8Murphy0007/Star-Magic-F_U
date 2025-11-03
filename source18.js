/**
 * ================================================================================================
 * Module: source18.js (PillarsOfCreation)
 *
 * Description: JavaScript Module for Pillars of Creation (Eagle Nebula) Class
 *              This is the seventh module in a series of 500+ code files for the Universal Quantum
 *              Field Framework (UQFF) simulations, focusing on star-forming pillars evolution and
 *              gravity equations derived from Hubble datasets, high-energy lab simulations, and
 *              UQFF refinements (dated May 09, 2025, updated for full term inclusion on October 08, 2025).
 *
 * Purpose: Encapsulates the Master Universal Gravity Equation (MUGE) for Pillars of Creation evolution.
 *          Includes ALL terms: base gravity with mass growth M(t), cosmic expansion (H_0), magnetic
 *          correction (static B), erosion E(t), UQFF Ug components with f_TRZ, Lambda, quantum uncertainty,
 *          scaled EM with [UA], fluid dynamics, oscillatory waves, DM/density perturbations, and stellar
 *          wind feedback (pressure / density for acc). Supports dynamic variable updates for all parameters.
 *
 * Integration: Designed for Node.js or browser environments.
 *              Usage: const pillars = new PillarsOfCreation();
 *                     const g = pillars.compute_g_Pillars(t);
 *
 * Key Features:
 *   - Default values from UQFF document: M_initial = 10,100 Msun, r = 4.731e16 m (5 ly),
 *     B = 1e-6 T, M_dot_factor ≈ 0.9901, tau_SF = 1 Myr, E_0 = 0.1, tau_erosion = 1 Myr,
 *     rho_wind = 1e-21 kg/m³, v_wind = 2e6 m/s.
 *   - Units handled: Msun to kg; wind term as (rho * v_wind²) / rho_fluid for acceleration.
 *   - Setter methods for updates: setVariable(name, new_val) or addToVariable(name, delta)/subtractFromVariable(name, delta).
 *   - Computes g_Pillars(r, t) with every term explicitly included.
 *   - UNIQUE EROSION FEATURE: E(t) = E_0 * exp(-t/tau_erosion) models photoevaporation by nearby massive stars
 *
 * Author: Converted by GitHub Copilot from C++ to JavaScript
 * Date: November 3, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

class PillarsOfCreation {
    // Static property for state storage (replaces anonymous namespace)
    static savedStates = {};

    constructor() {
        // Core parameters (all mutable for updates)
        this.G = 0;               // Gravitational constant
        this.M_initial = 0;       // Initial mass (kg)
        this.r = 0;               // Radius (m)
        this.H0 = 0;              // Hubble constant (s^-1)
        this.B = 0;               // Static magnetic field (T)
        this.B_crit = 0;          // Critical B field (T)
        this.Lambda = 0;          // Cosmological constant
        this.c_light = 0;         // Speed of light
        this.q_charge = 0;        // Charge (proton)
        this.gas_v = 0;           // Gas velocity for EM (m/s)
        this.f_TRZ = 0;           // Time-reversal factor
        this.M_dot_factor = 0;    // Star formation factor (dimensionless)
        this.tau_SF = 0;          // Star formation timescale (s)
        this.E_0 = 0;             // Initial erosion factor
        this.tau_erosion = 0;     // Erosion timescale (s)
        this.rho_wind = 0;        // Wind density (kg/m³)
        this.v_wind = 0;          // Wind velocity (m/s)
        this.rho_fluid = 0;       // Fluid density (kg/m³)
        this.rho_vac_UA = 0;      // UA vacuum density (J/m³)
        this.rho_vac_SCm = 0;     // SCm vacuum density (J/m³)
        this.scale_EM = 0;        // EM scaling factor
        this.proton_mass = 0;     // Proton mass for EM acceleration

        // Additional parameters for full inclusion of terms
        this.hbar = 0;            // Reduced Planck's constant
        this.t_Hubble = 0;        // Hubble time (s)
        this.delta_x = 0;         // Position uncertainty (m)
        this.delta_p = 0;         // Momentum uncertainty (kg m/s)
        this.integral_psi = 0;    // Wavefunction integral approximation
        this.A_osc = 0;           // Oscillatory amplitude (m/s²)
        this.k_osc = 0;           // Wave number (1/m)
        this.omega_osc = 0;       // Angular frequency (rad/s)
        this.x_pos = 0;           // Position for oscillation (m)
        this.t_Hubble_gyr = 0;    // Hubble time in Gyr
        this.M_DM_factor = 0;     // Dark matter mass fraction
        this.delta_rho_over_rho = 0; // Density perturbation fraction

        // Computed cache (updated on demand)
        this.ug1_base = 0;        // Cached Ug1 for initial M

        this.initializeDefaults();
    }

    // Initialization method (called in constructor)
    initializeDefaults() {
        this.G = 6.6743e-11;
        const M_sun = 1.989e30;
        const M_initial_sun = 10100.0;
        this.M_initial = M_initial_sun * M_sun;
        this.r = 4.731e16;
        this.H0 = 2.184e-18;
        this.B = 1e-6;
        this.B_crit = 1e11;
        this.Lambda = 1.1e-52;
        this.c_light = 3e8;
        this.q_charge = 1.602e-19;
        this.gas_v = 1e5;
        this.f_TRZ = 0.1;
        this.M_dot_factor = 1e4 / M_initial_sun;
        this.tau_SF = 1e6 * 3.156e7;
        this.E_0 = 0.1;
        this.tau_erosion = 1e6 * 3.156e7;
        this.rho_wind = 1e-21;
        this.v_wind = 2e6;
        this.rho_fluid = 1e-21;
        this.rho_vac_UA = 7.09e-36;
        this.rho_vac_SCm = 7.09e-37;
        this.scale_EM = 1e-12;
        this.proton_mass = 1.673e-27;

        // Full terms defaults
        this.hbar = 1.0546e-34;
        this.t_Hubble = 13.8e9 * 3.156e7;
        this.t_Hubble_gyr = 13.8;
        this.delta_x = 1e-10;
        this.delta_p = this.hbar / this.delta_x;
        this.integral_psi = 1.0;
        this.A_osc = 1e-10;  // Small for pillar scale
        this.k_osc = 1.0 / this.r;
        this.omega_osc = 2 * Math.PI / (this.r / this.c_light);
        this.x_pos = this.r;
        this.M_DM_factor = 0.1;
        this.delta_rho_over_rho = 1e-5;

        this.updateCache();
    }

    // Cache update for efficiency (call after parameter changes)
    updateCache() {
        this.ug1_base = (this.G * this.M_initial) / (this.r * this.r);
    }

    // Universal setter for any variable (by name, for flexibility)
    setVariable(varName, newValue) {
        if (varName === "G") { this.G = newValue; }
        else if (varName === "M_initial") { this.M_initial = newValue; }
        else if (varName === "r") { this.r = newValue; }
        else if (varName === "H0") { this.H0 = newValue; }
        else if (varName === "B") { this.B = newValue; }
        else if (varName === "B_crit") { this.B_crit = newValue; }
        else if (varName === "Lambda") { this.Lambda = newValue; }
        else if (varName === "c_light") { this.c_light = newValue; }
        else if (varName === "q_charge") { this.q_charge = newValue; }
        else if (varName === "gas_v") { this.gas_v = newValue; }
        else if (varName === "f_TRZ") { this.f_TRZ = newValue; }
        else if (varName === "M_dot_factor") { this.M_dot_factor = newValue; }
        else if (varName === "tau_SF") { this.tau_SF = newValue; }
        else if (varName === "E_0") { this.E_0 = newValue; }
        else if (varName === "tau_erosion") { this.tau_erosion = newValue; }
        else if (varName === "rho_wind") { this.rho_wind = newValue; }
        else if (varName === "v_wind") { this.v_wind = newValue; }
        else if (varName === "rho_fluid") { this.rho_fluid = newValue; }
        else if (varName === "rho_vac_UA") { this.rho_vac_UA = newValue; }
        else if (varName === "rho_vac_SCm") { this.rho_vac_SCm = newValue; }
        else if (varName === "scale_EM") { this.scale_EM = newValue; }
        else if (varName === "proton_mass") { this.proton_mass = newValue; }
        // Full terms
        else if (varName === "hbar") { this.hbar = newValue; }
        else if (varName === "t_Hubble") { this.t_Hubble = newValue; }
        else if (varName === "t_Hubble_gyr") { this.t_Hubble_gyr = newValue; }
        else if (varName === "delta_x") { this.delta_x = newValue; }
        else if (varName === "delta_p") { this.delta_p = newValue; }
        else if (varName === "integral_psi") { this.integral_psi = newValue; }
        else if (varName === "A_osc") { this.A_osc = newValue; }
        else if (varName === "k_osc") { this.k_osc = newValue; }
        else if (varName === "omega_osc") { this.omega_osc = newValue; }
        else if (varName === "x_pos") { this.x_pos = newValue; }
        else if (varName === "M_DM_factor") { this.M_DM_factor = newValue; }
        else if (varName === "delta_rho_over_rho") { this.delta_rho_over_rho = newValue; }
        else {
            console.error(`Error: Unknown variable '${varName}'.`);
            return false;
        }
        this.updateCache();
        return true;
    }

    // Addition method for variables
    addToVariable(varName, delta) {
        return this.setVariable(varName, this.getVariable(varName) + delta);
    }

    // Subtraction method for variables
    subtractFromVariable(varName, delta) {
        return this.addToVariable(varName, -delta);
    }

    // Getter for any variable (helper for add/subtract)
    getVariable(varName) {
        if (varName === "G") return this.G;
        else if (varName === "M_initial") return this.M_initial;
        else if (varName === "r") return this.r;
        else if (varName === "H0") return this.H0;
        else if (varName === "B") return this.B;
        else if (varName === "B_crit") return this.B_crit;
        else if (varName === "Lambda") return this.Lambda;
        else if (varName === "c_light") return this.c_light;
        else if (varName === "q_charge") return this.q_charge;
        else if (varName === "gas_v") return this.gas_v;
        else if (varName === "f_TRZ") return this.f_TRZ;
        else if (varName === "M_dot_factor") return this.M_dot_factor;
        else if (varName === "tau_SF") return this.tau_SF;
        else if (varName === "E_0") return this.E_0;
        else if (varName === "tau_erosion") return this.tau_erosion;
        else if (varName === "rho_wind") return this.rho_wind;
        else if (varName === "v_wind") return this.v_wind;
        else if (varName === "rho_fluid") return this.rho_fluid;
        else if (varName === "rho_vac_UA") return this.rho_vac_UA;
        else if (varName === "rho_vac_SCm") return this.rho_vac_SCm;
        else if (varName === "scale_EM") return this.scale_EM;
        else if (varName === "proton_mass") return this.proton_mass;
        // Full terms
        else if (varName === "hbar") return this.hbar;
        else if (varName === "t_Hubble") return this.t_Hubble;
        else if (varName === "t_Hubble_gyr") return this.t_Hubble_gyr;
        else if (varName === "delta_x") return this.delta_x;
        else if (varName === "delta_p") return this.delta_p;
        else if (varName === "integral_psi") return this.integral_psi;
        else if (varName === "A_osc") return this.A_osc;
        else if (varName === "k_osc") return this.k_osc;
        else if (varName === "omega_osc") return this.omega_osc;
        else if (varName === "x_pos") return this.x_pos;
        else if (varName === "M_DM_factor") return this.M_DM_factor;
        else if (varName === "delta_rho_over_rho") return this.delta_rho_over_rho;
        else {
            console.error(`Error: Unknown variable '${varName}'.`);
            return 0.0;
        }
    }

    // M(t) computation (star formation mass growth)
    M_t(t) {
        const M_dot = this.M_dot_factor * Math.exp(-t / this.tau_SF);
        return this.M_initial * (1 + M_dot);
    }

    // E(t) computation (erosion factor - UNIQUE to Pillars of Creation)
    E_t(t) {
        return this.E_0 * Math.exp(-t / this.tau_erosion);
    }

    // Ug terms computation
    compute_Ug(Mt) {
        const Ug1 = (this.G * Mt) / (this.r * this.r);
        const Ug2 = 0.0;
        const Ug3 = 0.0;
        const corr_B = 1 - this.B / this.B_crit;
        const Ug4 = Ug1 * corr_B;
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1 + this.f_TRZ);
    }

    // Volume computation for fluid
    compute_V() {
        return (4.0 / 3.0) * Math.PI * this.r * this.r * this.r;
    }

    // Main MUGE computation (includes ALL terms)
    compute_g_Pillars(t) {
        if (t < 0) {
            console.error("Error: Time t must be non-negative.");
            return 0.0;
        }

        const Mt = this.M_t(t);
        const Et = this.E_t(t);
        const ug1_t = (this.G * Mt) / (this.r * this.r);

        // Term 1: Base + H0 + B + E corrections
        const corr_H = 1 + this.H0 * t;
        const corr_B = 1 - this.B / this.B_crit;
        const corr_E = 1 - Et;
        const term1 = ug1_t * corr_H * corr_B * corr_E;

        // Term 2: UQFF Ug with f_TRZ
        const term2 = this.compute_Ug(Mt);

        // Term 3: Lambda
        const term3 = (this.Lambda * this.c_light * this.c_light) / 3.0;

        // Term 4: Scaled EM with UA
        const cross_vB = this.gas_v * this.B;  // Magnitude, assuming perpendicular
        const em_base = (this.q_charge * cross_vB) / this.proton_mass;
        const corr_UA = 1 + (this.rho_vac_UA / this.rho_vac_SCm);
        const term4 = (em_base * corr_UA) * this.scale_EM;

        // Quantum uncertainty term
        const sqrt_unc = Math.sqrt(this.delta_x * this.delta_p);
        const term_q = (this.hbar / sqrt_unc) * this.integral_psi * (2 * Math.PI / this.t_Hubble);

        // Fluid term (effective acceleration)
        const V = this.compute_V();
        const term_fluid = (this.rho_fluid * V * ug1_t) / Mt;

        // Oscillatory terms (real parts)
        const term_osc1 = 2 * this.A_osc * Math.cos(this.k_osc * this.x_pos) * Math.cos(this.omega_osc * t);
        const arg = this.k_osc * this.x_pos - this.omega_osc * t;
        const term_osc2 = (2 * Math.PI / this.t_Hubble_gyr) * this.A_osc * Math.cos(arg);
        const term_osc = term_osc1 + term_osc2;

        // DM and density perturbation term (converted to acceleration)
        const M_dm = Mt * this.M_DM_factor;
        const pert1 = this.delta_rho_over_rho;
        const pert2 = 3 * this.G * Mt / (this.r * this.r * this.r);
        const term_dm_force_like = (Mt + M_dm) * (pert1 + pert2);
        const term_DM = term_dm_force_like / Mt;

        // Stellar wind feedback term (pressure / density for acceleration)
        const wind_pressure = this.rho_wind * this.v_wind * this.v_wind;
        const term_wind = wind_pressure / this.rho_fluid;

        // Total g_Pillars (all terms summed)
        return term1 + term2 + term3 + term4 + term_q + term_fluid + term_osc + term_DM + term_wind;
    }

    // Debug/Output method (for transparency)
    printParameters() {
        console.log("Pillars of Creation Parameters:");
        console.log(`G: ${this.G.toFixed(3)}, M_initial: ${this.M_initial.toFixed(3)}, r: ${this.r.toFixed(3)}`);
        console.log(`H0: ${this.H0.toFixed(3)}, B: ${this.B.toFixed(3)}, B_crit: ${this.B_crit.toFixed(3)}`);
        console.log(`f_TRZ: ${this.f_TRZ.toFixed(3)}, M_dot_factor: ${this.M_dot_factor.toFixed(3)}, tau_SF: ${this.tau_SF.toFixed(3)}`);
        console.log(`E_0: ${this.E_0.toFixed(3)}, tau_erosion: ${this.tau_erosion.toFixed(3)}`);
        console.log(`rho_fluid: ${this.rho_fluid.toFixed(3)}, rho_wind: ${this.rho_wind.toFixed(3)}, v_wind: ${this.v_wind.toFixed(3)}`);
        console.log(`gas_v: ${this.gas_v.toFixed(3)}, M_DM_factor: ${this.M_DM_factor.toFixed(3)}`);
        console.log(`A_osc: ${this.A_osc.toFixed(3)}, delta_rho_over_rho: ${this.delta_rho_over_rho.toFixed(3)}`);
        console.log(`ug1_base: ${this.ug1_base.toFixed(3)}`);
    }

    // Example computation at t=500k years (for testing)
    exampleAt500kYears() {
        const t_example = 5e5 * 3.156e7;
        return this.compute_g_Pillars(t_example);
    }

    // ========== ENHANCED DYNAMIC CAPABILITIES (25 methods) ==========

    // --- Variable Management (5 methods) ---
    createVariable(name, value) {
        return this.setVariable(name, value);
    }

    removeVariable(name) {
        console.warn(`Warning: Cannot remove built-in variable '${name}' in PillarsOfCreation class.`);
        return false;
    }

    cloneVariable(src, dest) {
        const val = this.getVariable(src);
        return this.setVariable(dest, val);
    }

    listVariables() {
        return ["G", "M_initial", "r", "H0", "B", "B_crit", "Lambda", "c_light", "q_charge",
                "gas_v", "f_TRZ", "M_dot_factor", "tau_SF", "E_0", "tau_erosion", "rho_wind", "v_wind",
                "rho_fluid", "rho_vac_UA", "rho_vac_SCm", "scale_EM", "proton_mass", "hbar", "t_Hubble",
                "t_Hubble_gyr", "delta_x", "delta_p", "integral_psi", "A_osc", "k_osc", "omega_osc",
                "x_pos", "M_DM_factor", "delta_rho_over_rho"];
    }

    getSystemName() {
        return "PillarsOfCreation";
    }

    // --- Batch Operations (2 methods) ---
    transformVariableGroup(names, func) {
        for (const name of names) {
            const val = this.getVariable(name);
            if (!this.setVariable(name, func(val))) return false;
        }
        return true;
    }

    scaleVariableGroup(names, factor) {
        return this.transformVariableGroup(names, v => v * factor);
    }

    // --- Self-Expansion (4 methods) ---
    expandParameterSpace(factor) {
        const expandable = ["M_initial", "r", "B", "rho_fluid", "rho_wind", "A_osc", "M_DM_factor"];
        return this.scaleVariableGroup(expandable, factor);
    }

    expandStarFormationScale(M_dot_factor_scale, tau_SF_scale) {
        this.setVariable("M_dot_factor", this.getVariable("M_dot_factor") * M_dot_factor_scale);
        this.setVariable("tau_SF", this.getVariable("tau_SF") * tau_SF_scale);
    }

    expandErosionScale(E_0_scale, tau_erosion_scale) {
        this.setVariable("E_0", this.getVariable("E_0") * E_0_scale);
        this.setVariable("tau_erosion", this.getVariable("tau_erosion") * tau_erosion_scale);
    }

    expandWindMagneticScale(rho_wind_scale, v_wind_scale, B_scale) {
        this.setVariable("rho_wind", this.getVariable("rho_wind") * rho_wind_scale);
        this.setVariable("v_wind", this.getVariable("v_wind") * v_wind_scale);
        this.setVariable("B", this.getVariable("B") * B_scale);
        this.setVariable("B_crit", this.getVariable("B_crit") * B_scale);
    }

    // --- Self-Refinement (3 methods) ---
    autoRefineParameters(observations) {
        if (observations.length === 0) return;
        
        let sum_error = 0.0;
        for (const obs of observations) {
            const t = obs[0];
            const g_obs = obs[1];
            const g_calc = this.compute_g_Pillars(t);
            sum_error += Math.abs(g_calc - g_obs);
        }
        const avg_error = sum_error / observations.length;
        
        if (avg_error > 1e-6) {
            const adj_factor = 1.0 - Math.min(0.1, avg_error / 1e6);
            this.setVariable("M_dot_factor", this.getVariable("M_dot_factor") * adj_factor);
            this.setVariable("tau_SF", this.getVariable("tau_SF") * (2.0 - adj_factor));
        }
    }

    calibrateToObservations(times, g_obs) {
        if (times.length !== g_obs.length || times.length === 0) return;
        
        const obs = [];
        for (let i = 0; i < times.length; i++) {
            obs.push([times[i], g_obs[i]]);
        }
        
        for (let iter = 0; iter < 5; iter++) {
            this.autoRefineParameters(obs);
        }
    }

    optimizeForMetric(metric, t_start, t_end, steps) {
        let best_score = -1e100;
        const dt = (t_end - t_start) / steps;
        
        for (let i = 0; i <= steps; i++) {
            const t = t_start + i * dt;
            const g = this.compute_g_Pillars(t);
            const score = metric(g);
            if (score > best_score) best_score = score;
        }
        return best_score;
    }

    // --- Parameter Exploration (1 method) ---
    generateVariations(count, variation_pct) {
        const variations = [];
        const vars = this.listVariables();
        
        for (let i = 0; i < count; i++) {
            const variant = {};
            for (const v of vars) {
                const val = this.getVariable(v);
                const variation = val * (1.0 + (Math.random() * 2 - 1) * (variation_pct / 100.0));
                variant[v] = variation;
            }
            variations.push(variant);
        }
        return variations;
    }

    // --- Adaptive Evolution (2 methods) ---
    mutateParameters(mutation_rate) {
        const vars = this.listVariables();
        
        for (const v of vars) {
            if (v === "c_light" || v === "G" || v === "hbar") continue;
            const val = this.getVariable(v);
            const delta = val * ((Math.random() * 2 - 1) * mutation_rate);
            this.setVariable(v, val + delta);
        }
    }

    evolveSystem(generations, fitness) {
        let best_fitness = fitness(this);
        this.saveState("evolution_best");
        
        for (let gen = 0; gen < generations; gen++) {
            this.saveState("evolution_temp");
            this.mutateParameters(0.05);
            
            const new_fitness = fitness(this);
            if (new_fitness > best_fitness) {
                best_fitness = new_fitness;
                this.saveState("evolution_best");
            } else {
                this.restoreState("evolution_temp");
            }
        }
        
        this.restoreState("evolution_best");
    }

    // --- State Management (4 methods) ---
    saveState(stateName) {
        const state = {};
        const vars = this.listVariables();
        for (const v of vars) {
            state[v] = this.getVariable(v);
        }
        PillarsOfCreation.savedStates[stateName] = state;
        return true;
    }

    restoreState(stateName) {
        const state = PillarsOfCreation.savedStates[stateName];
        if (!state) return false;
        
        for (const [key, value] of Object.entries(state)) {
            this.setVariable(key, value);
        }
        return true;
    }

    listSavedStates() {
        return Object.keys(PillarsOfCreation.savedStates);
    }

    exportState() {
        let output = "PillarsOfCreation State Export:\n";
        const vars = this.listVariables();
        for (const v of vars) {
            output += `${v} = ${this.getVariable(v).toExponential(6)}\n`;
        }
        return output;
    }

    // --- System Analysis (4 methods) ---
    sensitivityAnalysis(t, delta_pct) {
        const sensitivities = {};
        const g_base = this.compute_g_Pillars(t);
        
        const vars = this.listVariables();
        for (const v of vars) {
            if (v === "c_light" || v === "G" || v === "hbar") continue;
            
            const original = this.getVariable(v);
            const delta = original * delta_pct / 100.0;
            
            this.setVariable(v, original + delta);
            const g_plus = this.compute_g_Pillars(t);
            this.setVariable(v, original);
            
            const sensitivity = (g_base !== 0.0) ? Math.abs((g_plus - g_base) / g_base) : 0.0;
            sensitivities[v] = sensitivity;
        }
        
        return sensitivities;
    }

    generateReport(t) {
        let output = "============================================\n";
        output += "PILLARS OF CREATION REPORT\n";
        output += "Eagle Nebula (M16)\n";
        output += "============================================\n";
        output += `Time: t = ${t.toFixed(6)} s (${(t/3.156e7/1e6).toFixed(6)} Myr)\n\n`;
        
        output += "Physical Parameters:\n";
        const M_sun = 1.989e30;
        output += `  Initial Mass M_initial = ${this.M_initial.toFixed(6)} kg (${(this.M_initial/M_sun).toFixed(6)} M_sun)\n`;
        output += `  M(t) = ${this.M_t(t).toFixed(6)} kg (${(this.M_t(t)/M_sun).toFixed(6)} M_sun)\n`;
        const ly_to_m = 9.461e15;
        output += `  Radius r = ${this.r.toFixed(6)} m (${(this.r/ly_to_m).toFixed(6)} ly)\n`;
        output += `  Magnetic field B = ${this.B.toFixed(6)} T (B_crit = ${this.B_crit.toFixed(6)} T)\n`;
        output += `  Star formation M_dot_factor = ${this.M_dot_factor.toFixed(6)}, tau_SF = ${this.tau_SF.toFixed(6)} s\n`;
        output += `  Erosion E_0 = ${this.E_0.toFixed(6)}, tau_erosion = ${this.tau_erosion.toFixed(6)} s\n`;
        output += `  E(t) = ${this.E_t(t).toFixed(6)}\n`;
        output += `  Wind density rho_wind = ${this.rho_wind.toFixed(6)} kg/m³, v_wind = ${this.v_wind.toFixed(6)} m/s\n`;
        output += `  Fluid density rho_fluid = ${this.rho_fluid.toFixed(6)} kg/m³\n`;
        output += `  Gas velocity = ${this.gas_v.toFixed(6)} m/s\n`;
        output += `  DM factor = ${this.M_DM_factor.toFixed(6)}\n\n`;
        
        output += "Computed Acceleration:\n";
        output += `  g_Pillars(t) = ${this.compute_g_Pillars(t).toFixed(6)} m/s²\n\n`;
        
        output += "UQFF Terms:\n";
        const Mt = this.M_t(t);
        const Et = this.E_t(t);
        const ug1_t = (this.G * Mt) / (this.r * this.r);
        const corr_H = 1 + this.H0 * t;
        const corr_B = 1 - this.B / this.B_crit;
        const corr_E = 1 - Et;
        output += `  Base (with H0, B, E, M(t)): ${(ug1_t * corr_H * corr_B * corr_E).toFixed(6)} m/s²\n`;
        output += `  Ug total: ${this.compute_Ug(Mt).toFixed(6)} m/s²\n`;
        output += `  Lambda: ${((this.Lambda * this.c_light * this.c_light) / 3.0).toExponential(6)} m/s²\n`;
        
        const cross_vB = this.gas_v * this.B;
        const em_base = (this.q_charge * cross_vB) / this.proton_mass;
        const corr_UA = 1 + (this.rho_vac_UA / this.rho_vac_SCm);
        output += `  EM (scaled with UA): ${(em_base * corr_UA * this.scale_EM).toExponential(6)} m/s²\n`;
        
        const sqrt_unc = Math.sqrt(this.delta_x * this.delta_p);
        output += `  Quantum: ${((this.hbar / sqrt_unc) * this.integral_psi * (2 * Math.PI / this.t_Hubble)).toExponential(6)} m/s²\n`;
        
        const V = this.compute_V();
        output += `  Fluid: ${((this.rho_fluid * V * ug1_t) / Mt).toExponential(6)} m/s²\n`;
        
        output += "  Oscillatory: (combined real parts)\n";
        
        const M_dm = Mt * this.M_DM_factor;
        const pert1 = this.delta_rho_over_rho;
        const pert2 = 3 * this.G * Mt / (this.r * this.r * this.r);
        const term_dm_force_like = (Mt + M_dm) * (pert1 + pert2);
        output += `  DM: ${(term_dm_force_like / Mt).toExponential(6)} m/s²\n`;
        
        const wind_pressure = this.rho_wind * this.v_wind * this.v_wind;
        output += `  Stellar Wind Feedback: ${(wind_pressure / this.rho_fluid).toExponential(6)} m/s²\n`;
        
        output += "============================================\n";
        return output;
    }

    validateConsistency() {
        let valid = true;
        
        if (this.M_initial <= 0 || this.r <= 0) { console.error("Error: M_initial and r must be positive."); valid = false; }
        if (this.B < 0 || this.B_crit <= 0) { console.error("Error: B, B_crit must be non-negative/positive."); valid = false; }
        if (this.tau_SF <= 0 || this.tau_erosion <= 0) { console.error("Error: Timescales must be positive."); valid = false; }
        if (this.E_0 < 0 || this.E_0 > 1.0) { console.warn("Warning: Erosion factor E_0 outside [0,1]."); }
        if (this.rho_fluid <= 0 || this.rho_wind < 0) { console.error("Error: Fluid/wind densities must be positive/non-negative."); valid = false; }
        if (this.v_wind < 0 || this.gas_v < 0) { console.error("Error: Velocities must be non-negative."); valid = false; }
        if (this.M_DM_factor < 0 || this.M_DM_factor > 1.0) { console.warn("Warning: DM factor outside [0,1]."); }
        
        return valid;
    }

    autoCorrectAnomalies() {
        let corrected = false;
        
        const M_sun = 1.989e30;
        
        if (this.M_initial <= 0) { this.M_initial = 10100.0 * M_sun; corrected = true; }
        if (this.r <= 0) { this.r = 4.731e16; corrected = true; }
        if (this.B < 0) { this.B = 1e-6; corrected = true; }
        if (this.B_crit <= 0) { this.B_crit = 1e11; corrected = true; }
        if (this.tau_SF <= 0) { this.tau_SF = 1e6 * 3.156e7; corrected = true; }
        if (this.tau_erosion <= 0) { this.tau_erosion = 1e6 * 3.156e7; corrected = true; }
        if (this.E_0 < 0) { this.E_0 = 0.0; corrected = true; }
        if (this.E_0 > 1.0) { this.E_0 = 1.0; corrected = true; }
        if (this.rho_fluid <= 0) { this.rho_fluid = 1e-21; corrected = true; }
        if (this.rho_wind < 0) { this.rho_wind = 1e-21; corrected = true; }
        if (this.v_wind < 0) { this.v_wind = 2e6; corrected = true; }
        if (this.gas_v < 0) { this.gas_v = 1e5; corrected = true; }
        if (this.M_DM_factor < 0) { this.M_DM_factor = 0.1; corrected = true; }
        if (this.M_DM_factor > 1.0) { this.M_DM_factor = 1.0; corrected = true; }
        
        if (corrected) this.updateCache();
        return corrected;
    }
}

// ========== ENHANCED EXAMPLE FUNCTION ==========
function enhancedPillarsExample() {
    console.log("=========================================================");
    console.log("ENHANCED PILLARS OF CREATION DEMONSTRATION");
    console.log("Eagle Nebula (M16) - Star Formation with Erosion");
    console.log("=========================================================\n");
    
    const pillars = new PillarsOfCreation();
    
    // Step 1: Initial state and validation
    console.log("Step 1: Initial State and Validation");
    console.log(`System: ${pillars.getSystemName()}`);
    console.log(`Validation: ${pillars.validateConsistency() ? "PASS" : "FAIL"}`);
    console.log(`Auto-corrected: ${pillars.autoCorrectAnomalies() ? "Yes" : "No"}\n`);
    
    // Step 2: Time evolution showing M(t) and E(t)
    console.log("Step 2: Time Evolution (M(t) mass growth and E(t) erosion)");
    const t_Myr_array = [0.0, 0.5, 1.0, 2.0, 5.0];
    for (const t_Myr of t_Myr_array) {
        const t = t_Myr * 1e6 * 3.156e7;
        const M_sun = 1.989e30;
        const Mt = pillars.M_t(t);
        const Et = pillars.E_t(t);
        const g = pillars.compute_g_Pillars(t);
        console.log(`  t = ${t_Myr.toFixed(1)} Myr: M(t) = ${(Mt/M_sun).toFixed(2)} M_sun, E(t) = ${Et.toFixed(6)}, g = ${g.toExponential(6)} m/s²`);
    }
    console.log();
    
    // Step 3: Variable listing
    console.log("Step 3: Variable Listing");
    const vars = pillars.listVariables();
    console.log(`Total variables: ${vars.length}`);
    console.log(`Sample: ${vars[0]}, ${vars[1]}, ${vars[13]} (E_0), ${vars[14]} (tau_erosion)\n`);
    
    // Step 4: Star formation scaling
    console.log("Step 4: Star Formation Scaling (M_dot_factor sweeps)");
    pillars.saveState("original");
    const M_dot_factors = [0.5, 1.0, 2.0];
    for (const factor of M_dot_factors) {
        pillars.restoreState("original");
        pillars.expandStarFormationScale(factor, 1.0);
        const t = 1e6 * 3.156e7;
        const Mt = pillars.M_t(t);
        const M_sun = 1.989e30;
        console.log(`  M_dot_factor × ${factor.toFixed(1)}: M(1 Myr) = ${(Mt/M_sun).toFixed(2)} M_sun`);
    }
    pillars.restoreState("original");
    console.log();
    
    // Step 5: Erosion scaling (UNIQUE to this module)
    console.log("Step 5: Erosion Scaling (E_0 sweeps) - UNIQUE FEATURE");
    const E_0_factors = [0.5, 1.0, 2.0];
    for (const factor of E_0_factors) {
        pillars.restoreState("original");
        pillars.expandErosionScale(factor, 1.0);
        const t = 1e6 * 3.156e7;
        const Et = pillars.E_t(t);
        const g = pillars.compute_g_Pillars(t);
        console.log(`  E_0 × ${factor.toFixed(1)}: E(1 Myr) = ${Et.toFixed(6)}, g = ${g.toExponential(6)} m/s²`);
    }
    pillars.restoreState("original");
    console.log();
    
    // Step 6: Erosion timescale sweeps
    console.log("Step 6: Erosion Timescale Sweeps (tau_erosion)");
    const tau_erosion_factors = [0.5, 1.0, 2.0];
    for (const factor of tau_erosion_factors) {
        pillars.restoreState("original");
        pillars.expandErosionScale(1.0, factor);
        const t = 1e6 * 3.156e7;
        const Et = pillars.E_t(t);
        console.log(`  tau_erosion × ${factor.toFixed(1)}: E(1 Myr) = ${Et.toFixed(6)}`);
    }
    pillars.restoreState("original");
    console.log();
    
    // Step 7: Wind velocity scaling
    console.log("Step 7: Wind Velocity Scaling");
    const v_wind_factors = [0.5, 1.0, 2.0];
    for (const factor of v_wind_factors) {
        pillars.restoreState("original");
        pillars.expandWindMagneticScale(1.0, factor, 1.0);
        const t = 1e6 * 3.156e7;
        const g = pillars.compute_g_Pillars(t);
        console.log(`  v_wind × ${factor.toFixed(1)}: g(1 Myr) = ${g.toExponential(6)} m/s²`);
    }
    pillars.restoreState("original");
    console.log();
    
    // Step 8: Magnetic field scaling
    console.log("Step 8: Magnetic Field Scaling");
    const B_factors = [0.5, 1.0, 2.0];
    for (const factor of B_factors) {
        pillars.restoreState("original");
        pillars.expandWindMagneticScale(1.0, 1.0, factor);
        const t = 1e6 * 3.156e7;
        const g = pillars.compute_g_Pillars(t);
        console.log(`  B × ${factor.toFixed(1)}: g(1 Myr) = ${g.toExponential(6)} m/s²`);
    }
    pillars.restoreState("original");
    console.log();
    
    // Step 9: Parameter space expansion
    console.log("Step 9: Parameter Space Expansion (all scalable params)");
    pillars.expandParameterSpace(1.2);
    const M_after = pillars.getVariable("M_initial");
    const M_sun = 1.989e30;
    console.log(`  After 1.2× expansion: M_initial = ${(M_after/M_sun).toFixed(2)} M_sun`);
    pillars.restoreState("original");
    console.log();
    
    // Step 10: Batch operations
    console.log("Step 10: Batch Operations (scale multiple variables)");
    const scale_group = ["M_initial", "r", "B"];
    pillars.scaleVariableGroup(scale_group, 1.1);
    console.log("  Scaled {M_initial, r, B} by 1.1×");
    pillars.restoreState("original");
    console.log();
    
    // Step 11: State management
    console.log("Step 11: State Management");
    pillars.saveState("state_A");
    pillars.expandStarFormationScale(1.5, 1.5);
    pillars.saveState("state_B");
    const states = pillars.listSavedStates();
    console.log(`  Saved states: ${states.join(" ")}`);
    pillars.restoreState("state_A");
    console.log("  Restored state_A\n");
    
    // Step 12: Generate parameter variations
    console.log("Step 12: Generate Parameter Variations (5% variation)");
    const variations = pillars.generateVariations(3, 5.0);
    console.log(`  Generated ${variations.length} variants with 5% random variation`);
    console.log(`  Variant 1 M_initial = ${variations[0].M_initial.toExponential(6)} kg\n`);
    
    // Step 13: Sensitivity analysis
    console.log("Step 13: Sensitivity Analysis at 1 Myr");
    const t_sens = 1e6 * 3.156e7;
    const sensitivities = pillars.sensitivityAnalysis(t_sens, 1.0);
    console.log("  Top sensitivities (1% perturbation):");
    const sens_vec = Object.entries(sensitivities).sort((a, b) => b[1] - a[1]);
    for (let i = 0; i < 5 && i < sens_vec.length; i++) {
        console.log(`    ${sens_vec[i][0]}: ${sens_vec[i][1].toExponential(3)}`);
    }
    console.log();
    
    // Step 14: Auto-refinement with synthetic observations
    console.log("Step 14: Auto-Refinement (synthetic observations)");
    const obs = [];
    for (let i = 0; i <= 5; i++) {
        const t_obs = i * 1e6 * 3.156e7;
        const g_obs = pillars.compute_g_Pillars(t_obs) * (1.0 + 0.01 * (Math.random() * 100 - 50) / 100.0);
        obs.push([t_obs, g_obs]);
    }
    pillars.autoRefineParameters(obs);
    console.log(`  Refined parameters based on ${obs.length} observations\n`);
    
    // Step 15: Optimization for maximum acceleration
    console.log("Step 15: Optimize for Maximum Acceleration");
    pillars.restoreState("original");
    const metric = g => g;
    const t_opt_start = 0.0;
    const t_opt_end = 5e6 * 3.156e7;
    const best_g = pillars.optimizeForMetric(metric, t_opt_start, t_opt_end, 50);
    console.log(`  Best g over 5 Myr: ${best_g.toExponential(6)} m/s²\n`);
    
    // Step 16: Evolutionary system adaptation
    console.log("Step 16: Evolutionary System Adaptation (5 generations)");
    pillars.restoreState("original");
    const fitness = p => {
        const t = 1e6 * 3.156e7;
        return p.compute_g_Pillars(t);
    };
    pillars.evolveSystem(5, fitness);
    console.log("  Evolved system over 5 generations (fitness = g at 1 Myr)\n");
    
    // Step 17: Full system report
    console.log("Step 17: Full System Report at 1 Myr");
    pillars.restoreState("original");
    const t_report = 1e6 * 3.156e7;
    const report = pillars.generateReport(t_report);
    console.log(report);
    
    // Step 18: Full state export
    console.log("Step 18: Full State Export");
    const exported = pillars.exportState();
    console.log("Exported state (first 500 chars):");
    console.log(exported.substring(0, 500) + "...\n");
    
    console.log("=========================================================");
    console.log("ENHANCED DEMONSTRATION COMPLETE");
    console.log("=========================================================");
}

// Module exports for Node.js and browser
if (typeof module !== 'undefined' && module.exports) {
    module.exports = PillarsOfCreation;
    module.exports.enhancedPillarsExample = enhancedPillarsExample;
}

// Default export for ES6 modules
export default PillarsOfCreation;
export { enhancedPillarsExample };
