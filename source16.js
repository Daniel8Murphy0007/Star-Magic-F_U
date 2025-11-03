/**
 * ================================================================================================
 * Module: source16.js (StarbirthTapestry)
 *
 * Description: JavaScript Module for "Tapestry of Blazing Starbirth" (NGC 2014 & NGC 2020) Class
 *              This is the fourth module in a series of 500+ code files for the Universal Quantum
 *              Field Framework (UQFF) simulations, focusing on star-forming region evolution and
 *              gravity equations derived from Hubble datasets, high-energy lab simulations, and
 *              UQFF refinements (dated May 09, 2025, updated for full term inclusion on October 08, 2025).
 *
 * Purpose: Encapsulates the Master Universal Gravity Equation (MUGE) for "Tapestry of Blazing
 *          Starbirth" evolution in the Large Magellanic Cloud (LMC). Includes ALL terms: base 
 *          gravity with star formation growth M(t), cosmic expansion (H_0), magnetic correction 
 *          (static B), UQFF Ug components with f_TRZ, Lambda, quantum uncertainty, scaled EM with 
 *          [UA], fluid dynamics, oscillatory waves, DM/density perturbations, and stellar wind 
 *          feedback (pressure / density for acc). Supports dynamic variable updates for all parameters.
 *
 * Integration: Designed for Node.js and browser environments with ES6 module syntax.
 *              Usage: import StarbirthTapestry from './source16.js';
 *                     const tapestry = new StarbirthTapestry();
 *                     const g = tapestry.compute_g_Starbirth(t);
 *
 * Key Features:
 *   - Default values from UQFF document: M_initial = 240 Msun, r = 10 ly, B = 1e-6 T, etc.
 *   - M(t) = M_initial * (1 + M_dot_factor * exp(-t / tau_SF)), with M_dot_factor = 10000 / 240.
 *   - Units handled: ly to m, Msun to kg; wind term as (rho * v_wind^2) / rho_fluid for acceleration.
 *   - Setter methods for updates: setVariable(newValue) or addToVariable(delta)/subtractFromVariable(delta).
 *   - Computes g_Starbirth(r, t) with every term explicitly included.
 *   - 25 enhanced dynamic capabilities for analysis and adaptation.
 *
 * Converted from: source16.cpp (C++ module)
 * Original Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Conversion Date: November 03, 2025
 * JavaScript Conversion: GitHub Copilot
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

class StarbirthTapestry {
    // Static property for state storage (replaces anonymous namespace)
    static savedStates = {};

    constructor() {
        // Core parameters (all mutable for updates)
        this.G = 0.0;               // Gravitational constant
        this.M_initial = 0.0;       // Initial mass (kg)
        this.r = 0.0;               // Radius (m)
        this.H0 = 0.0;              // Hubble constant (s^-1)
        this.B = 0.0;               // Static magnetic field (T)
        this.B_crit = 0.0;          // Critical B field (T)
        this.Lambda = 0.0;          // Cosmological constant
        this.c_light = 0.0;         // Speed of light
        this.q_charge = 0.0;        // Charge (proton)
        this.gas_v = 0.0;           // Gas velocity for EM (m/s)
        this.f_TRZ = 0.0;           // Time-reversal factor
        this.M_dot_factor = 0.0;    // Star formation factor (dimensionless)
        this.tau_SF = 0.0;          // Star formation timescale (s)
        this.rho_wind = 0.0;        // Wind density (kg/m^3)
        this.v_wind = 0.0;          // Wind velocity (m/s)
        this.rho_fluid = 0.0;       // Fluid density (kg/m^3)
        this.rho_vac_UA = 0.0;      // UA vacuum density (J/m^3)
        this.rho_vac_SCm = 0.0;     // SCm vacuum density (J/m^3)
        this.scale_EM = 0.0;        // EM scaling factor

        // Additional parameters for full inclusion of terms
        this.hbar = 0.0;            // Reduced Planck's constant
        this.t_Hubble = 0.0;        // Hubble time (s)
        this.delta_x = 0.0;         // Position uncertainty (m)
        this.delta_p = 0.0;         // Momentum uncertainty (kg m/s)
        this.integral_psi = 0.0;    // Wavefunction integral approximation
        this.A_osc = 0.0;           // Oscillatory amplitude (m/s^2)
        this.k_osc = 0.0;           // Wave number (1/m)
        this.omega_osc = 0.0;       // Angular frequency (rad/s)
        this.x_pos = 0.0;           // Position for oscillation (m)
        this.t_Hubble_gyr = 0.0;    // Hubble time in Gyr
        this.M_DM_factor = 0.0;     // Dark matter mass fraction
        this.delta_rho_over_rho = 0.0; // Density perturbation fraction
        this.proton_mass = 0.0;     // Proton mass for EM acceleration

        // Computed caches (updated on demand)
        this.ug1_base = 0.0;        // Cached Ug1 for initial M

        // Initialize with default UQFF values
        this.initializeDefaults();
    }

    /**
     * Initialize all parameters with default UQFF values for NGC 2014 & NGC 2020.
     */
    initializeDefaults() {
        this.G = 6.6743e-11;
        const M_sun = 1.989e30;
        const M_initial_sun = 240.0;
        this.M_initial = M_initial_sun * M_sun;
        const ly_to_m = 9.461e15;
        this.r = 10.0 * ly_to_m;
        this.H0 = 2.184e-18;
        this.B = 1e-6;
        this.B_crit = 1e11;
        this.Lambda = 1.1e-52;
        this.c_light = 3e8;
        this.q_charge = 1.602e-19;
        this.gas_v = 1e5;
        this.f_TRZ = 0.1;
        const gas_mass_sun = 10000.0;
        this.M_dot_factor = gas_mass_sun / M_initial_sun;
        this.tau_SF = 5e6 * 3.156e7;
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
        this.A_osc = 1e-10;  // Small for nebula scale
        this.k_osc = 1.0 / this.r;
        this.omega_osc = 2 * Math.PI / (this.r / this.c_light);
        this.x_pos = this.r;
        this.M_DM_factor = 0.1;
        this.delta_rho_over_rho = 1e-5;

        this.updateCache();
    }

    /**
     * Update cached values for efficiency (call after parameter changes).
     */
    updateCache() {
        this.ug1_base = (this.G * this.M_initial) / (this.r * this.r);
    }

    /**
     * Universal setter for any variable (by name, for flexibility).
     * @param {string} varName - Name of the variable to set.
     * @param {number} newValue - New value to assign.
     * @returns {boolean} - True if successful, false if variable not found.
     */
    setVariable(varName, newValue) {
        const validVars = [
            'G', 'M_initial', 'r', 'H0', 'B', 'B_crit', 'Lambda', 'c_light', 'q_charge',
            'gas_v', 'f_TRZ', 'M_dot_factor', 'tau_SF', 'rho_wind', 'v_wind', 'rho_fluid',
            'rho_vac_UA', 'rho_vac_SCm', 'scale_EM', 'hbar', 't_Hubble', 't_Hubble_gyr',
            'delta_x', 'delta_p', 'integral_psi', 'A_osc', 'k_osc', 'omega_osc', 'x_pos',
            'M_DM_factor', 'delta_rho_over_rho', 'proton_mass'
        ];

        if (validVars.includes(varName)) {
            this[varName] = newValue;
            this.updateCache();
            return true;
        } else {
            console.error(`Error: Unknown variable '${varName}'.`);
            return false;
        }
    }

    /**
     * Add a delta value to an existing variable.
     * @param {string} varName - Name of the variable.
     * @param {number} delta - Value to add.
     * @returns {boolean} - True if successful.
     */
    addToVariable(varName, delta) {
        return this.setVariable(varName, this.getVariable(varName) + delta);
    }

    /**
     * Subtract a delta value from an existing variable.
     * @param {string} varName - Name of the variable.
     * @param {number} delta - Value to subtract.
     * @returns {boolean} - True if successful.
     */
    subtractFromVariable(varName, delta) {
        return this.addToVariable(varName, -delta);
    }

    /**
     * Getter for any variable (helper for add/subtract).
     * @param {string} varName - Name of the variable.
     * @returns {number} - Current value of the variable, or 0.0 if not found.
     */
    getVariable(varName) {
        if (this.hasOwnProperty(varName)) {
            return this[varName];
        } else {
            console.error(`Error: Unknown variable '${varName}'.`);
            return 0.0;
        }
    }

    /**
     * Compute time-dependent mass M(t) with star formation.
     * M(t) = M_initial * (1 + M_dot_factor * exp(-t / tau_SF))
     * @param {number} t - Time in seconds.
     * @returns {number} - Mass at time t in kg.
     */
    M_t(t) {
        const M_dot = this.M_dot_factor * Math.exp(-t / this.tau_SF);
        return this.M_initial * (1 + M_dot);
    }

    /**
     * Compute UQFF Ug terms (Ug1 + Ug2 + Ug3 + Ug4) * (1 + f_TRZ).
     * @param {number} Mt - Current mass at time t.
     * @returns {number} - Total Ug acceleration (m/s^2).
     */
    compute_Ug(Mt) {
        const Ug1 = (this.G * Mt) / (this.r * this.r);
        const Ug2 = 0.0;  // No separate Ug2 component for this system
        const Ug3 = 0.0;  // No separate Ug3 component for this system
        const corr_B = 1 - this.B / this.B_crit;
        const Ug4 = Ug1 * corr_B;
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1 + this.f_TRZ);
    }

    /**
     * Compute volume of the star-forming region (spherical approximation).
     * @returns {number} - Volume in m^3.
     */
    compute_V() {
        return (4.0 / 3.0) * Math.PI * this.r * this.r * this.r;
    }

    /**
     * Main MUGE computation for Tapestry of Blazing Starbirth (includes ALL terms).
     * Returns total gravitational acceleration at distance r and time t.
     * @param {number} t - Time in seconds.
     * @returns {number} - Total acceleration g_Starbirth (m/s^2).
     */
    compute_g_Starbirth(t) {
        if (t < 0) {
            console.error('Error: Time t must be non-negative.');
            return 0.0;
        }

        const Mt = this.M_t(t);
        const ug1_t = (this.G * Mt) / (this.r * this.r);

        // Term 1: Base + H0 + B corrections
        const corr_H = 1 + this.H0 * t;
        const corr_B = 1 - this.B / this.B_crit;
        const term1 = ug1_t * corr_H * corr_B;

        // Term 2: UQFF Ug with f_TRZ
        const term2 = this.compute_Ug(Mt);

        // Term 3: Lambda (cosmological constant)
        const term3 = (this.Lambda * this.c_light * this.c_light) / 3.0;

        // Term 4: Scaled EM with UA
        const cross_vB = this.gas_v * this.B;  // Magnitude, assuming perpendicular
        const em_base = (this.q_charge * cross_vB) / this.proton_mass;
        const corr_UA = 1 + (this.rho_vac_UA / this.rho_vac_SCm);
        const term4 = (em_base * corr_UA) * this.scale_EM;

        // Term 5: Quantum uncertainty term
        const sqrt_unc = Math.sqrt(this.delta_x * this.delta_p);
        const term_q = (this.hbar / sqrt_unc) * this.integral_psi * (2 * Math.PI / this.t_Hubble);

        // Term 6: Fluid term (effective acceleration)
        const V = this.compute_V();
        const term_fluid = (this.rho_fluid * V * ug1_t) / Mt;

        // Term 7: Oscillatory terms (real parts, standing + traveling wave)
        const term_osc1 = 2 * this.A_osc * Math.cos(this.k_osc * this.x_pos) * Math.cos(this.omega_osc * t);
        const arg = this.k_osc * this.x_pos - this.omega_osc * t;
        const term_osc2 = (2 * Math.PI / this.t_Hubble_gyr) * this.A_osc * Math.cos(arg);
        const term_osc = term_osc1 + term_osc2;

        // Term 8: DM and density perturbation term (converted to acceleration)
        const M_dm = Mt * this.M_DM_factor;
        const pert1 = this.delta_rho_over_rho;
        const pert2 = 3 * this.G * Mt / (this.r * this.r * this.r);
        const term_dm_force_like = (Mt + M_dm) * (pert1 + pert2);
        const term_DM = term_dm_force_like / Mt;

        // Term 9: Stellar wind feedback term (pressure / density for acceleration)
        const wind_pressure = this.rho_wind * this.v_wind * this.v_wind;
        const term_wind = wind_pressure / this.rho_fluid;

        // Total g_Starbirth (all terms summed)
        return term1 + term2 + term3 + term4 + term_q + term_fluid + term_osc + term_DM + term_wind;
    }

    /**
     * Debug/Output method (for transparency in console).
     * @param {object} stream - Output stream (default: console).
     */
    printParameters(stream = console) {
        stream.log('Tapestry of Blazing Starbirth Parameters:');
        stream.log(`G: ${this.G.toFixed(3)}, M_initial: ${this.M_initial.toExponential(3)}, r: ${this.r.toExponential(3)}`);
        stream.log(`H0: ${this.H0.toExponential(3)}, B: ${this.B.toExponential(3)}, B_crit: ${this.B_crit.toExponential(3)}`);
        stream.log(`f_TRZ: ${this.f_TRZ.toFixed(3)}, M_dot_factor: ${this.M_dot_factor.toFixed(3)}, tau_SF: ${this.tau_SF.toExponential(3)}`);
        stream.log(`rho_fluid: ${this.rho_fluid.toExponential(3)}, rho_wind: ${this.rho_wind.toExponential(3)}, v_wind: ${this.v_wind.toExponential(3)}`);
        stream.log(`gas_v: ${this.gas_v.toExponential(3)}, M_DM_factor: ${this.M_DM_factor.toFixed(3)}`);
        stream.log(`A_osc: ${this.A_osc.toExponential(3)}, delta_rho_over_rho: ${this.delta_rho_over_rho.toExponential(3)}`);
        stream.log(`ug1_base: ${this.ug1_base.toExponential(3)}`);
    }

    /**
     * Example computation at t=2.5 Myr (for quick testing).
     * @returns {number} - Acceleration at t=2.5 Myr (m/s^2).
     */
    exampleAt2_5Myr() {
        const t_example = 2.5e6 * 3.156e7;
        return this.compute_g_Starbirth(t_example);
    }

    // ========== ENHANCED DYNAMIC CAPABILITIES (25 methods) ==========

    // --- Variable Management (5 methods) ---

    /**
     * Create a new variable (for this class, syncs to existing member variable).
     * @param {string} name - Variable name.
     * @param {number} value - Initial value.
     * @returns {boolean} - True if successful.
     */
    createVariable(name, value) {
        return this.setVariable(name, value);
    }

    /**
     * Remove a variable (not allowed for built-in member variables).
     * @param {string} name - Variable name.
     * @returns {boolean} - Always false for built-in variables.
     */
    removeVariable(name) {
        console.error(`Warning: Cannot remove built-in variable '${name}' in StarbirthTapestry class.`);
        return false;
    }

    /**
     * Clone a variable (copy value from src to dest).
     * @param {string} src - Source variable name.
     * @param {string} dest - Destination variable name.
     * @returns {boolean} - True if successful.
     */
    cloneVariable(src, dest) {
        const val = this.getVariable(src);
        return this.setVariable(dest, val);
    }

    /**
     * List all available variables in the system.
     * @returns {Array<string>} - Array of variable names.
     */
    listVariables() {
        return [
            'G', 'M_initial', 'r', 'H0', 'B', 'B_crit', 'Lambda', 'c_light', 'q_charge',
            'gas_v', 'f_TRZ', 'M_dot_factor', 'tau_SF', 'rho_wind', 'v_wind', 'rho_fluid',
            'rho_vac_UA', 'rho_vac_SCm', 'scale_EM', 'hbar', 't_Hubble', 't_Hubble_gyr',
            'delta_x', 'delta_p', 'integral_psi', 'A_osc', 'k_osc', 'omega_osc', 'x_pos',
            'M_DM_factor', 'delta_rho_over_rho', 'proton_mass'
        ];
    }

    /**
     * Get the system name identifier.
     * @returns {string} - System name.
     */
    getSystemName() {
        return 'StarbirthTapestry';
    }

    // --- Batch Operations (2 methods) ---

    /**
     * Transform a group of variables using a function.
     * @param {Array<string>} names - Variable names to transform.
     * @param {function} func - Transformation function (value => new_value).
     * @returns {boolean} - True if all successful.
     */
    transformVariableGroup(names, func) {
        for (const name of names) {
            const val = this.getVariable(name);
            if (!this.setVariable(name, func(val))) return false;
        }
        return true;
    }

    /**
     * Scale a group of variables by a constant factor.
     * @param {Array<string>} names - Variable names to scale.
     * @param {number} factor - Scaling factor.
     * @returns {boolean} - True if successful.
     */
    scaleVariableGroup(names, factor) {
        return this.transformVariableGroup(names, v => v * factor);
    }

    // --- Self-Expansion (4 methods) ---

    /**
     * Expand parameter space by scaling key physical parameters.
     * @param {number} factor - Scaling factor for expansion.
     */
    expandParameterSpace(factor) {
        const expandable = ['M_initial', 'r', 'B', 'rho_fluid', 'rho_wind', 'A_osc', 'M_DM_factor'];
        this.scaleVariableGroup(expandable, factor);
    }

    /**
     * Expand star formation scale (M_dot_factor and tau_SF).
     * @param {number} M_dot_factor_scale - Scale factor for M_dot_factor.
     * @param {number} tau_SF_scale - Scale factor for tau_SF.
     */
    expandStarFormationScale(M_dot_factor_scale, tau_SF_scale) {
        this.setVariable('M_dot_factor', this.getVariable('M_dot_factor') * M_dot_factor_scale);
        this.setVariable('tau_SF', this.getVariable('tau_SF') * tau_SF_scale);
    }

    /**
     * Expand wind feedback scale (rho_wind and v_wind).
     * @param {number} rho_wind_scale - Scale factor for rho_wind.
     * @param {number} v_wind_scale - Scale factor for v_wind.
     */
    expandWindFeedbackScale(rho_wind_scale, v_wind_scale) {
        this.setVariable('rho_wind', this.getVariable('rho_wind') * rho_wind_scale);
        this.setVariable('v_wind', this.getVariable('v_wind') * v_wind_scale);
    }

    /**
     * Expand magnetic and fluid scale (B, B_crit, rho_fluid).
     * @param {number} B_scale - Scale factor for B and B_crit.
     * @param {number} rho_fluid_scale - Scale factor for rho_fluid.
     */
    expandMagneticFluidScale(B_scale, rho_fluid_scale) {
        this.setVariable('B', this.getVariable('B') * B_scale);
        this.setVariable('B_crit', this.getVariable('B_crit') * B_scale);
        this.setVariable('rho_fluid', this.getVariable('rho_fluid') * rho_fluid_scale);
    }

    // --- Self-Refinement (3 methods) ---

    /**
     * Auto-refine parameters based on observations (time, g_obs pairs).
     * @param {Array<Array<number>>} observations - Array of [time, g_obs] pairs.
     */
    autoRefineParameters(observations) {
        if (observations.length === 0) return;

        let sum_error = 0.0;
        for (const obs of observations) {
            const t = obs[0];
            const g_obs = obs[1];
            const g_calc = this.compute_g_Starbirth(t);
            sum_error += Math.abs(g_calc - g_obs);
        }
        const avg_error = sum_error / observations.length;

        // Simple refinement: adjust M_dot_factor and tau_SF based on error
        if (avg_error > 1e-6) {
            const adj_factor = 1.0 - Math.min(0.1, avg_error / 1e6);
            this.setVariable('M_dot_factor', this.getVariable('M_dot_factor') * adj_factor);
            this.setVariable('tau_SF', this.getVariable('tau_SF') * (2.0 - adj_factor));
        }
    }

    /**
     * Calibrate to observations using iterative refinement.
     * @param {Array<number>} times - Array of time values.
     * @param {Array<number>} g_obs - Array of observed g values.
     */
    calibrateToObservations(times, g_obs) {
        if (times.length !== g_obs.length || times.length === 0) return;

        const obs = times.map((t, i) => [t, g_obs[i]]);

        // Iterative refinement (5 passes)
        for (let iter = 0; iter < 5; iter++) {
            this.autoRefineParameters(obs);
        }
    }

    /**
     * Optimize system for a given metric function over time range.
     * @param {function} metric - Metric function to optimize (g => score).
     * @param {number} t_start - Start time.
     * @param {number} t_end - End time.
     * @param {number} steps - Number of time steps.
     * @returns {number} - Best score found.
     */
    optimizeForMetric(metric, t_start, t_end, steps) {
        let best_score = -1e100;
        const dt = (t_end - t_start) / steps;

        for (let i = 0; i <= steps; i++) {
            const t = t_start + i * dt;
            const g = this.compute_g_Starbirth(t);
            const score = metric(g);
            if (score > best_score) best_score = score;
        }
        return best_score;
    }

    // --- Parameter Exploration (1 method) ---

    /**
     * Generate parameter variations using Monte Carlo sampling.
     * @param {number} count - Number of variations to generate.
     * @param {number} variation_pct - Percentage variation (e.g., 10 for ±10%).
     * @returns {Array<Object>} - Array of parameter variation objects.
     */
    generateVariations(count, variation_pct) {
        const variations = [];
        const vars = this.listVariables();

        for (let i = 0; i < count; i++) {
            const variant = {};
            for (const v of vars) {
                const val = this.getVariable(v);
                const variation = val * (1.0 + (Math.random() * 2 - 1) * variation_pct / 100.0);
                variant[v] = variation;
            }
            variations.push(variant);
        }
        return variations;
    }

    // --- Adaptive Evolution (2 methods) ---

    /**
     * Mutate parameters randomly for evolutionary algorithms.
     * @param {number} mutation_rate - Mutation rate (fractional, e.g., 0.05 for 5%).
     */
    mutateParameters(mutation_rate) {
        const vars = this.listVariables();

        for (const v of vars) {
            if (v === 'c_light' || v === 'G' || v === 'hbar') continue; // Skip constants

            const val = this.getVariable(v);
            const delta = val * (Math.random() * 2 - 1) * mutation_rate;
            this.setVariable(v, val + delta);
        }
    }

    /**
     * Evolve system using genetic algorithm approach.
     * @param {number} generations - Number of generations to evolve.
     * @param {function} fitness - Fitness function (system => score).
     */
    evolveSystem(generations, fitness) {
        let best_fitness = fitness(this);
        this.saveState('evolution_best');

        for (let gen = 0; gen < generations; gen++) {
            this.saveState('evolution_temp');
            this.mutateParameters(0.05);

            const new_fitness = fitness(this);
            if (new_fitness > best_fitness) {
                best_fitness = new_fitness;
                this.saveState('evolution_best');
            } else {
                this.restoreState('evolution_temp');
            }
        }

        this.restoreState('evolution_best');
    }

    // --- State Management (4 methods) ---

    /**
     * Save current state with a given name.
     * @param {string} stateName - Name for the saved state.
     * @returns {boolean} - True if successful.
     */
    saveState(stateName) {
        const state = {};
        const vars = this.listVariables();
        for (const v of vars) {
            state[v] = this.getVariable(v);
        }
        StarbirthTapestry.savedStates[stateName] = state;
        return true;
    }

    /**
     * Restore previously saved state.
     * @param {string} stateName - Name of the state to restore.
     * @returns {boolean} - True if successful, false if state not found.
     */
    restoreState(stateName) {
        const state = StarbirthTapestry.savedStates[stateName];
        if (!state) return false;

        for (const [key, value] of Object.entries(state)) {
            this.setVariable(key, value);
        }
        return true;
    }

    /**
     * List all saved state names.
     * @returns {Array<string>} - Array of saved state names.
     */
    listSavedStates() {
        return Object.keys(StarbirthTapestry.savedStates);
    }

    /**
     * Export current state as formatted string.
     * @returns {string} - Formatted state export.
     */
    exportState() {
        let output = 'StarbirthTapestry State Export:\n';
        const vars = this.listVariables();
        for (const v of vars) {
            output += `${v} = ${this.getVariable(v).toExponential(6)}\n`;
        }
        return output;
    }

    // --- System Analysis (4 methods) ---

    /**
     * Perform sensitivity analysis on all parameters at given time.
     * @param {number} t - Time for analysis.
     * @param {number} delta_pct - Perturbation percentage (e.g., 1 for 1%).
     * @returns {Object} - Object mapping variable names to sensitivity values.
     */
    sensitivityAnalysis(t, delta_pct) {
        const sensitivities = {};
        const g_base = this.compute_g_Starbirth(t);

        const vars = this.listVariables();
        for (const v of vars) {
            if (v === 'c_light' || v === 'G' || v === 'hbar') continue; // Skip constants

            const original = this.getVariable(v);
            const delta = original * delta_pct / 100.0;

            this.setVariable(v, original + delta);
            const g_plus = this.compute_g_Starbirth(t);
            this.setVariable(v, original);

            const sensitivity = (g_base !== 0.0) ? Math.abs((g_plus - g_base) / g_base) : 0.0;
            sensitivities[v] = sensitivity;
        }

        return sensitivities;
    }

    /**
     * Generate comprehensive system report at given time.
     * @param {number} t - Time for report.
     * @returns {string} - Formatted report string.
     */
    generateReport(t) {
        let output = '';
        output += '============================================\n';
        output += 'TAPESTRY OF BLAZING STARBIRTH REPORT\n';
        output += 'NGC 2014 & NGC 2020 (LMC)\n';
        output += '============================================\n';
        output += `Time: t = ${t.toFixed(6)} s (${(t/3.156e7/1e6).toFixed(3)} Myr)\n\n`;

        output += 'Physical Parameters:\n';
        const M_sun = 1.989e30;
        output += `  Initial Mass M_initial = ${this.M_initial.toExponential(6)} kg (${(this.M_initial/M_sun).toFixed(2)} M_sun)\n`;
        output += `  M(t) = ${this.M_t(t).toExponential(6)} kg (${(this.M_t(t)/M_sun).toFixed(2)} M_sun)\n`;
        const ly_to_m = 9.461e15;
        output += `  Radius r = ${this.r.toExponential(6)} m (${(this.r/ly_to_m).toFixed(2)} ly)\n`;
        output += `  Magnetic field B = ${this.B.toExponential(6)} T (B_crit = ${this.B_crit.toExponential(3)} T)\n`;
        output += `  Star formation M_dot_factor = ${this.M_dot_factor.toFixed(3)}, tau_SF = ${this.tau_SF.toExponential(3)} s\n`;
        output += `  Wind density rho_wind = ${this.rho_wind.toExponential(6)} kg/m^3, v_wind = ${this.v_wind.toExponential(3)} m/s\n`;
        output += `  Fluid density rho_fluid = ${this.rho_fluid.toExponential(6)} kg/m^3\n`;
        output += `  Gas velocity = ${this.gas_v.toExponential(3)} m/s\n`;
        output += `  DM factor = ${this.M_DM_factor.toFixed(3)}\n\n`;

        output += 'Computed Acceleration:\n';
        output += `  g_Starbirth(t) = ${this.compute_g_Starbirth(t).toExponential(6)} m/s^2\n\n`;

        output += 'UQFF Terms:\n';
        const Mt = this.M_t(t);
        const ug1_t = (this.G * Mt) / (this.r * this.r);
        const corr_H = 1 + this.H0 * t;
        const corr_B = 1 - this.B / this.B_crit;
        output += `  Base (with H0, B, M(t)): ${(ug1_t * corr_H * corr_B).toExponential(6)} m/s^2\n`;
        output += `  Ug total: ${this.compute_Ug(Mt).toExponential(6)} m/s^2\n`;
        output += `  Lambda: ${((this.Lambda * this.c_light * this.c_light) / 3.0).toExponential(6)} m/s^2\n`;

        const cross_vB = this.gas_v * this.B;
        const em_base = (this.q_charge * cross_vB) / this.proton_mass;
        const corr_UA = 1 + (this.rho_vac_UA / this.rho_vac_SCm);
        output += `  EM (scaled with UA): ${(em_base * corr_UA * this.scale_EM).toExponential(6)} m/s^2\n`;

        const sqrt_unc = Math.sqrt(this.delta_x * this.delta_p);
        output += `  Quantum: ${((this.hbar / sqrt_unc) * this.integral_psi * (2 * Math.PI / this.t_Hubble)).toExponential(6)} m/s^2\n`;

        const V = this.compute_V();
        output += `  Fluid: ${((this.rho_fluid * V * ug1_t) / Mt).toExponential(6)} m/s^2\n`;

        output += `  Oscillatory: (combined real parts)\n`;

        const M_dm = Mt * this.M_DM_factor;
        const pert1 = this.delta_rho_over_rho;
        const pert2 = 3 * this.G * Mt / (this.r * this.r * this.r);
        const term_dm_force_like = (Mt + M_dm) * (pert1 + pert2);
        output += `  DM: ${(term_dm_force_like / Mt).toExponential(6)} m/s^2\n`;

        const wind_pressure = this.rho_wind * this.v_wind * this.v_wind;
        output += `  Stellar Wind Feedback: ${(wind_pressure / this.rho_fluid).toExponential(6)} m/s^2\n`;

        output += '============================================\n';
        return output;
    }

    /**
     * Validate consistency of all parameters.
     * @returns {boolean} - True if all parameters are physically valid.
     */
    validateConsistency() {
        let valid = true;

        if (this.M_initial <= 0 || this.r <= 0) {
            console.error('Error: M_initial and r must be positive.');
            valid = false;
        }
        if (this.B < 0 || this.B_crit <= 0) {
            console.error('Error: B, B_crit must be non-negative/positive.');
            valid = false;
        }
        if (this.tau_SF <= 0) {
            console.error('Error: Star formation timescale must be positive.');
            valid = false;
        }
        if (this.rho_fluid <= 0 || this.rho_wind < 0) {
            console.error('Error: Fluid/wind densities must be positive/non-negative.');
            valid = false;
        }
        if (this.v_wind < 0 || this.gas_v < 0) {
            console.error('Error: Velocities must be non-negative.');
            valid = false;
        }
        if (this.M_DM_factor < 0 || this.M_DM_factor > 1.0) {
            console.warn('Warning: DM factor outside [0,1].');
        }

        return valid;
    }

    /**
     * Auto-correct anomalies in parameters (reset to defaults if invalid).
     * @returns {boolean} - True if any corrections were made.
     */
    autoCorrectAnomalies() {
        let corrected = false;

        const M_sun = 1.989e30;
        const ly_to_m = 9.461e15;

        if (this.M_initial <= 0) { this.M_initial = 240.0 * M_sun; corrected = true; }
        if (this.r <= 0) { this.r = 10.0 * ly_to_m; corrected = true; }
        if (this.B < 0) { this.B = 1e-6; corrected = true; }
        if (this.B_crit <= 0) { this.B_crit = 1e11; corrected = true; }
        if (this.tau_SF <= 0) { this.tau_SF = 5e6 * 3.156e7; corrected = true; }
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

/**
 * Enhanced example demonstrating all StarbirthTapestry capabilities.
 */
function enhancedStarbirthTapestryExample() {
    console.log('\n========== ENHANCED STARBIRTH TAPESTRY UQFF EXAMPLE ==========\n');

    const tapestry = new StarbirthTapestry();

    // Step 1: Initial state
    console.log('Step 1: Initial Configuration');
    tapestry.printParameters();
    const t0 = 0.0;
    console.log(`g_Starbirth(t=0) = ${tapestry.compute_g_Starbirth(t0).toExponential(6)} m/s^2\n`);

    // Step 2: Time evolution (0, 1, 2.5, 5, 10 Myr)
    console.log('Step 2: Time Evolution (0, 1, 2.5, 5, 10 Myr)');
    const M_sun = 1.989e30;
    for (const t_myr of [0.0, 1.0, 2.5, 5.0, 10.0]) {
        const t = t_myr * 1e6 * 3.156e7;
        console.log(`  t = ${t_myr} Myr: g = ${tapestry.compute_g_Starbirth(t).toExponential(6)} m/s^2, M(t) = ${(tapestry.M_t(t)/M_sun).toFixed(2)} M_sun`);
    }
    console.log();

    // Step 3: Star formation scaling
    console.log('Step 3: Star Formation Scaling (M_dot_factor x1.3, tau_SF x0.8)');
    tapestry.expandStarFormationScale(1.3, 0.8);
    const t_test = 2.5e6 * 3.156e7;
    console.log(`After expansion: M_dot_factor = ${tapestry.getVariable('M_dot_factor').toFixed(3)}, tau_SF = ${tapestry.getVariable('tau_SF').toExponential(3)} s`);
    console.log(`g_Starbirth(t=2.5 Myr) = ${tapestry.compute_g_Starbirth(t_test).toExponential(6)} m/s^2\n`);

    // Step 4: Wind feedback scaling
    console.log('Step 4: Wind Feedback Scaling (rho_wind x1.5, v_wind x1.2)');
    tapestry.expandWindFeedbackScale(1.5, 1.2);
    console.log(`After expansion: rho_wind = ${tapestry.getVariable('rho_wind').toExponential(6)} kg/m^3, v_wind = ${tapestry.getVariable('v_wind').toExponential(3)} m/s`);
    console.log(`g_Starbirth(t=2.5 Myr) = ${tapestry.compute_g_Starbirth(t_test).toExponential(6)} m/s^2\n`);

    // Step 5: Magnetic & fluid scaling
    console.log('Step 5: Magnetic & Fluid Scaling (B x1.4, rho_fluid x1.3)');
    tapestry.expandMagneticFluidScale(1.4, 1.3);
    console.log(`After expansion: B = ${tapestry.getVariable('B').toExponential(6)} T, rho_fluid = ${tapestry.getVariable('rho_fluid').toExponential(6)} kg/m^3`);
    console.log(`g_Starbirth(t=2.5 Myr) = ${tapestry.compute_g_Starbirth(t_test).toExponential(6)} m/s^2\n`);

    // Step 6: State save/restore
    console.log('Step 6: State Management');
    tapestry.saveState('expanded_state');
    tapestry.setVariable('M_initial', 500.0 * M_sun);
    console.log(`Modified M_initial to 500 M_sun: g = ${tapestry.compute_g_Starbirth(t_test).toExponential(6)} m/s^2`);
    tapestry.restoreState('expanded_state');
    console.log(`Restored state: M_initial = ${(tapestry.getVariable('M_initial')/M_sun).toFixed(2)} M_sun, g = ${tapestry.compute_g_Starbirth(t_test).toExponential(6)} m/s^2`);
    console.log(`Saved states: ${tapestry.listSavedStates().join(' ')}\n`);

    // Step 7: Sensitivity analysis
    console.log('Step 7: Sensitivity Analysis at t=2.5 Myr (top 5 parameters)');
    const sens = tapestry.sensitivityAnalysis(t_test, 1.0);
    const sens_array = Object.entries(sens).sort((a, b) => b[1] - a[1]);
    for (let i = 0; i < Math.min(5, sens_array.length); i++) {
        console.log(`  ${sens_array[i][0]}: ${sens_array[i][1].toExponential(3)}`);
    }
    console.log();

    // Step 8: Generate variations
    console.log('Step 8: Generate Parameter Variations (3 variants, 10% variation)');
    const variations = tapestry.generateVariations(3, 10.0);
    for (let i = 0; i < variations.length; i++) {
        console.log(`  Variant ${i+1}: M_initial = ${(variations[i].M_initial/M_sun).toFixed(2)} M_sun, v_wind = ${variations[i].v_wind.toExponential(3)} m/s`);
    }
    console.log();

    // Step 9: Batch transformation
    console.log('Step 9: Batch Transform (scale density parameters by 1.2)');
    tapestry.transformVariableGroup(['rho_fluid', 'rho_wind'], v => v * 1.2);
    console.log(`After transform: rho_fluid = ${tapestry.getVariable('rho_fluid').toExponential(6)} kg/m^3, rho_wind = ${tapestry.getVariable('rho_wind').toExponential(6)} kg/m^3`);
    console.log(`g_Starbirth(t=2.5 Myr) = ${tapestry.compute_g_Starbirth(t_test).toExponential(6)} m/s^2\n`);

    // Step 10: Consistency validation
    console.log('Step 10: Consistency Validation');
    const valid = tapestry.validateConsistency();
    console.log(`System is ${valid ? 'VALID' : 'INVALID'}\n`);

    // Step 11: Metric optimization
    console.log('Step 11: Optimize for Maximum g (t=0 to 10 Myr, 100 steps)');
    const max_g = tapestry.optimizeForMetric(g => g, 0.0, 10e6 * 3.156e7, 100);
    console.log(`Maximum g found: ${max_g.toExponential(6)} m/s^2\n`);

    // Step 12: Full system report
    console.log('Step 12: Full System Report at t=3 Myr');
    const t_report = 3e6 * 3.156e7;
    console.log(tapestry.generateReport(t_report));

    // Step 13: Initial mass sweep
    console.log('Step 13: Initial Mass Sweep (M_initial = 100, 240, 400, 600 M_sun)');
    tapestry.saveState('before_sweep');
    for (const M_val of [100.0, 240.0, 400.0, 600.0]) {
        tapestry.setVariable('M_initial', M_val * M_sun);
        console.log(`  M_initial = ${M_val} M_sun: g(t=2.5Myr) = ${tapestry.compute_g_Starbirth(t_test).toExponential(6)} m/s^2`);
    }
    tapestry.restoreState('before_sweep');
    console.log();

    // Step 14: Star formation rate sweep
    console.log('Step 14: Star Formation Factor Sweep (M_dot_factor = 20, 41.7, 60, 80)');
    for (const M_dot of [20.0, 41.7, 60.0, 80.0]) {
        tapestry.setVariable('M_dot_factor', M_dot);
        console.log(`  M_dot_factor = ${M_dot}: M(t=2.5Myr) = ${(tapestry.M_t(t_test)/M_sun).toFixed(2)} M_sun, g = ${tapestry.compute_g_Starbirth(t_test).toExponential(6)} m/s^2`);
    }
    tapestry.restoreState('before_sweep');
    console.log();

    // Step 15: Wind velocity sweep
    console.log('Step 15: Wind Velocity Sweep (v_wind = 1, 2, 3, 4 x 10^6 m/s)');
    for (const v_wind of [1e6, 2e6, 3e6, 4e6]) {
        tapestry.setVariable('v_wind', v_wind);
        console.log(`  v_wind = ${(v_wind/1e6).toFixed(1)} x 10^6 m/s: g(t=2.5Myr) = ${tapestry.compute_g_Starbirth(t_test).toExponential(6)} m/s^2`);
    }
    tapestry.restoreState('before_sweep');
    console.log();

    // Step 16: Magnetic field sweep
    console.log('Step 16: Magnetic Field Sweep (B = 0.5, 1.0, 1.5, 2.0 x 10^-6 T)');
    for (const B_val of [0.5e-6, 1.0e-6, 1.5e-6, 2.0e-6]) {
        tapestry.setVariable('B', B_val);
        console.log(`  B = ${(B_val*1e6).toFixed(1)} x 10^-6 T: g(t=2.5Myr) = ${tapestry.compute_g_Starbirth(t_test).toExponential(6)} m/s^2`);
    }
    tapestry.restoreState('before_sweep');
    console.log();

    // Step 17: DM factor sweep
    console.log('Step 17: Dark Matter Factor Sweep (M_DM_factor = 0.0, 0.1, 0.2, 0.3)');
    for (const dm of [0.0, 0.1, 0.2, 0.3]) {
        tapestry.setVariable('M_DM_factor', dm);
        console.log(`  M_DM_factor = ${dm.toFixed(1)}: g(t=2.5Myr) = ${tapestry.compute_g_Starbirth(t_test).toExponential(6)} m/s^2`);
    }
    tapestry.restoreState('before_sweep');
    console.log();

    // Step 18: Export final state
    console.log('Step 18: Export Final State');
    console.log(tapestry.exportState());

    console.log('========== ENHANCED EXAMPLE COMPLETE ==========\n');
}

// Export for module usage
export default StarbirthTapestry;
export { enhancedStarbirthTapestryExample };
