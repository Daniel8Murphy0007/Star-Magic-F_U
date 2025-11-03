/**
 * ================================================================================================
 * Module: HorseheadNebula.js
 *
 * Description: JavaScript ES6 Module for Horsehead Nebula (Barnard 33) Class
 *              This is the fifteenth module in a series of 500+ code files for the Universal
 *              Quantum Field Framework (UQFF) simulations, focusing on dark nebula evolution
 *              and gravity equations derived from Hubble datasets, high-energy lab simulations,
 *              and UQFF refinements (dated May 09, 2025, updated October 08, 2025).
 *
 * Purpose: Encapsulates the Master Universal Gravity Equation (MUGE) for Horsehead Nebula evolution.
 *          Includes ALL terms: base gravity (static M), cosmic expansion (H_0), magnetic correction
 *          (static B), erosion E(t), UQFF Ug components with f_TRZ, Lambda, quantum uncertainty,
 *          scaled EM with [UA], fluid dynamics, oscillatory waves, DM/density perturbations, and
 *          stellar wind feedback (pressure / density for acc).
 *          Supports dynamic variable updates for all parameters.
 *
 * Key Features:
 *   - Default values from UQFF document: M = 1000 M☉, r = 2.5 ly, B = 1e-6 T, E_0 = 0.1,
 *     tau_erosion = 5 Myr, rho_wind = 1e-21 kg/m³, v_wind = 2e6 m/s.
 *   - E(t) = E_0 × (1 - exp(-t/tau_erosion)) - erosion factor (UNIQUE to dark nebula)
 *   - Units handled: M☉ to kg, ly to m; wind term as (rho × v_wind²) / rho_fluid
 *   - Setter methods: setVariable(varName, newValue), addToVariable, subtractFromVariable
 *   - Computes g_Horsehead(r, t) with every term explicitly included (10 terms total)
 *   - Enhanced dynamic capabilities: 38 methods total (13 core + 25 enhanced)
 *
 * Conversion: Ported from C++ source24.cpp (HorseheadNebula.h) to JavaScript ES6
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript
 * Date: November 03, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

class HorseheadNebula {
    // Static storage for saved states (class-level)
    static savedStates = {};

    constructor() {
        this.initializeDefaults();
    }

    /**
     * Initialize all default parameters from UQFF specifications
     */
    initializeDefaults() {
        // Physical constants
        this.G = 6.6743e-11;                  // Gravitational constant (m³/kg/s²)
        const M_sun = 1.989e30;               // Solar mass (kg)
        this.M = 1000.0 * M_sun;              // Total mass (kg) - 1000 M☉
        const ly_to_m = 9.461e15;             // Light year to meters
        this.r = 2.5 * ly_to_m;               // Radius (m) - 2.5 ly
        this.H0 = 2.184e-18;                  // Hubble constant (s^-1)
        
        // Magnetic and cosmological parameters
        this.B = 1e-6;                        // Static magnetic field (T)
        this.B_crit = 1e11;                   // Critical B field (T)
        this.Lambda = 1.1e-52;                // Cosmological constant (m^-2)
        this.c_light = 3e8;                   // Speed of light (m/s)
        
        // EM parameters
        this.q_charge = 1.602e-19;            // Proton charge (C)
        this.gas_v = 1e5;                     // Gas velocity for EM (m/s)
        this.f_TRZ = 0.1;                     // Time-reversal factor
        
        // Erosion parameters (UNIQUE to dark nebula)
        this.E_0 = 0.1;                       // Initial erosion factor
        this.tau_erosion = 5e6 * 3.156e7;     // Erosion timescale (s) - 5 Myr
        
        // Wind and fluid parameters
        this.rho_wind = 1e-21;                // Wind density (kg/m³)
        this.v_wind = 2e6;                    // Wind velocity (m/s)
        this.rho_fluid = 1e-21;               // Fluid density (kg/m³)
        
        // Vacuum energy densities
        this.rho_vac_UA = 7.09e-36;           // UA vacuum density (J/m³)
        this.rho_vac_SCm = 7.09e-37;          // SCm vacuum density (J/m³)
        
        // Scaling parameters
        this.scale_EM = 1e-12;                // EM scaling factor
        this.proton_mass = 1.673e-27;         // Proton mass for EM acceleration (kg)
        
        // Quantum uncertainty parameters
        this.hbar = 1.0546e-34;               // Reduced Planck's constant (J·s)
        this.t_Hubble = 13.8e9 * 3.156e7;     // Hubble time (s)
        this.t_Hubble_gyr = 13.8;             // Hubble time (Gyr)
        this.delta_x = 1e-10;                 // Position uncertainty (m)
        this.delta_p = this.hbar / this.delta_x;  // Momentum uncertainty (kg·m/s)
        this.integral_psi = 1.0;              // Wavefunction integral approximation
        
        // Oscillatory wave parameters
        this.A_osc = 1e-10;                   // Oscillatory amplitude (m/s²)
        this.k_osc = 1.0 / this.r;            // Wave number (1/m)
        this.omega_osc = 2 * Math.PI / (this.r / this.c_light);  // Angular frequency (rad/s)
        this.x_pos = this.r;                  // Position for oscillation (m)
        
        // Dark matter and density perturbation parameters
        this.M_DM_factor = 0.1;               // Dark matter mass fraction
        this.delta_rho_over_rho = 1e-5;       // Density perturbation fraction
        
        this.updateCache();
    }

    /**
     * Update cached values (call after parameter changes)
     */
    updateCache() {
        this.ug1_base = (this.G * this.M) / (this.r * this.r);
    }

    // ========== CORE PHYSICS METHODS (13 methods) ==========

    /**
     * Compute erosion factor E(t)
     * E(t) = E_0 × (1 - exp(-t/tau_erosion)) - UNIQUE to dark nebula
     * Asymptotic approach to E_0 over time
     * @param {number} t - Time (s)
     * @returns {number} Erosion factor (dimensionless)
     */
    E_t(t) {
        return this.E_0 * (1 - Math.exp(-t / this.tau_erosion));
    }

    /**
     * Compute UQFF Ug components with erosion correction
     * @param {number} Et - Erosion factor at time t
     * @returns {number} Total Ug (m/s²)
     */
    compute_Ug(Et) {
        const Ug1 = this.ug1_base;
        const Ug2 = 0.0;
        const Ug3 = 0.0;
        const corr_B = 1 - this.B / this.B_crit;
        const Ug4 = Ug1 * corr_B;
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1 + this.f_TRZ) * (1 - Et);
    }

    /**
     * Compute volume for fluid dynamics
     * @returns {number} Volume (m³)
     */
    compute_V() {
        return (4.0 / 3.0) * Math.PI * this.r * this.r * this.r;
    }

    /**
     * Main MUGE computation - includes ALL 10 terms
     * @param {number} t - Time (s)
     * @returns {number} Total gravitational acceleration g_Horsehead (m/s²)
     */
    compute_g_Horsehead(t) {
        if (t < 0) {
            console.error("Error: Time t must be non-negative.");
            return 0.0;
        }

        const Et = this.E_t(t);

        // Term 1: Base + H0 + B + E corrections
        const corr_H = 1 + this.H0 * t;
        const corr_B = 1 - this.B / this.B_crit;
        const corr_E = 1 - Et;
        const term1 = this.ug1_base * corr_H * corr_B * corr_E;

        // Term 2: UQFF Ug with f_TRZ and E
        const term2 = this.compute_Ug(Et);

        // Term 3: Lambda (cosmological constant)
        const term3 = (this.Lambda * this.c_light * this.c_light) / 3.0;

        // Term 4: Scaled EM with UA
        const cross_vB = this.gas_v * this.B;  // Magnitude, assuming perpendicular
        const em_base = (this.q_charge * cross_vB) / this.proton_mass;
        const corr_UA = 1 + (this.rho_vac_UA / this.rho_vac_SCm);
        const term4 = (em_base * corr_UA) * this.scale_EM;

        // Term 5: Quantum uncertainty
        const sqrt_unc = Math.sqrt(this.delta_x * this.delta_p);
        const term_q = (this.hbar / sqrt_unc) * this.integral_psi * (2 * Math.PI / this.t_Hubble);

        // Term 6: Fluid dynamics (effective acceleration)
        const V = this.compute_V();
        const term_fluid = (this.rho_fluid * V * this.ug1_base) / this.M;

        // Term 7: Oscillatory waves (real parts)
        const term_osc1 = 2 * this.A_osc * Math.cos(this.k_osc * this.x_pos) * Math.cos(this.omega_osc * t);
        const arg = this.k_osc * this.x_pos - this.omega_osc * t;
        const term_osc2 = (2 * Math.PI / this.t_Hubble_gyr) * this.A_osc * Math.cos(arg);
        const term_osc = term_osc1 + term_osc2;

        // Term 8: Dark matter and density perturbation (converted to acceleration)
        const M_dm = this.M * this.M_DM_factor;
        const pert1 = this.delta_rho_over_rho;
        const pert2 = 3 * this.G * this.M / (this.r * this.r * this.r);
        const term_dm_force_like = (this.M + M_dm) * (pert1 + pert2);
        const term_DM = term_dm_force_like / this.M;

        // Term 9: Stellar wind feedback (pressure / density for acceleration)
        const wind_pressure = this.rho_wind * this.v_wind * this.v_wind;
        const term_wind = wind_pressure / this.rho_fluid;

        // Total g_Horsehead (all 9 terms summed - note: term_osc counts as 1)
        return term1 + term2 + term3 + term4 + term_q + term_fluid + term_osc + term_DM + term_wind;
    }

    /**
     * Universal getter for any variable by name
     * @param {string} varName - Variable name
     * @returns {number} Variable value
     */
    getVariable(varName) {
        const validVars = [
            'G', 'M', 'r', 'H0', 'B', 'B_crit', 'Lambda', 'c_light', 'q_charge', 'gas_v',
            'f_TRZ', 'E_0', 'tau_erosion', 'rho_wind', 'v_wind', 'rho_fluid',
            'rho_vac_UA', 'rho_vac_SCm', 'scale_EM', 'proton_mass',
            'hbar', 't_Hubble', 't_Hubble_gyr', 'delta_x', 'delta_p', 'integral_psi',
            'A_osc', 'k_osc', 'omega_osc', 'x_pos', 'M_DM_factor', 'delta_rho_over_rho'
        ];
        
        if (!validVars.includes(varName)) {
            console.error(`Error: Unknown variable '${varName}'.`);
            return 0.0;
        }
        return this[varName];
    }

    /**
     * Universal setter for any variable by name
     * @param {string} varName - Variable name
     * @param {number} newValue - New value
     * @returns {boolean} Success status
     */
    setVariable(varName, newValue) {
        const validVars = [
            'G', 'M', 'r', 'H0', 'B', 'B_crit', 'Lambda', 'c_light', 'q_charge', 'gas_v',
            'f_TRZ', 'E_0', 'tau_erosion', 'rho_wind', 'v_wind', 'rho_fluid',
            'rho_vac_UA', 'rho_vac_SCm', 'scale_EM', 'proton_mass',
            'hbar', 't_Hubble', 't_Hubble_gyr', 'delta_x', 'delta_p', 'integral_psi',
            'A_osc', 'k_osc', 'omega_osc', 'x_pos', 'M_DM_factor', 'delta_rho_over_rho'
        ];
        
        if (!validVars.includes(varName)) {
            console.error(`Error: Unknown variable '${varName}'.`);
            return false;
        }
        this[varName] = newValue;
        this.updateCache();
        return true;
    }

    /**
     * Add to variable value
     * @param {string} varName - Variable name
     * @param {number} delta - Amount to add
     * @returns {boolean} Success status
     */
    addToVariable(varName, delta) {
        return this.setVariable(varName, this.getVariable(varName) + delta);
    }

    /**
     * Subtract from variable value
     * @param {string} varName - Variable name
     * @param {number} delta - Amount to subtract
     * @returns {boolean} Success status
     */
    subtractFromVariable(varName, delta) {
        return this.addToVariable(varName, -delta);
    }

    /**
     * Example computation at t=3 Myr (for testing)
     * @returns {number} g_Horsehead at 3 Myr (m/s²)
     */
    exampleAt3Myr() {
        const t_example = 3e6 * 3.156e7;
        return this.compute_g_Horsehead(t_example);
    }

    /**
     * Print parameters to console
     */
    printParameters() {
        console.log("Horsehead Nebula Parameters:");
        console.log(`G: ${this.G.toExponential(3)}, M: ${this.M.toExponential(3)}, r: ${this.r.toExponential(3)}`);
        console.log(`H0: ${this.H0.toExponential(3)}, B: ${this.B.toExponential(3)}, B_crit: ${this.B_crit.toExponential(3)}`);
        console.log(`f_TRZ: ${this.f_TRZ.toFixed(3)}, E_0: ${this.E_0.toFixed(3)}, tau_erosion: ${this.tau_erosion.toExponential(3)}`);
        console.log(`rho_fluid: ${this.rho_fluid.toExponential(3)}, rho_wind: ${this.rho_wind.toExponential(3)}, v_wind: ${this.v_wind.toExponential(3)}`);
        console.log(`gas_v: ${this.gas_v.toExponential(3)}, M_DM_factor: ${this.M_DM_factor.toFixed(3)}`);
        console.log(`A_osc: ${this.A_osc.toExponential(3)}, delta_rho_over_rho: ${this.delta_rho_over_rho.toExponential(3)}`);
        console.log(`ug1_base: ${this.ug1_base.toExponential(3)}`);
    }

    // ========== ENHANCED DYNAMIC CAPABILITIES (25 methods) ==========

    // ----- Variable Management (5 methods) -----

    /**
     * Create or update a variable
     * @param {string} name - Variable name
     * @param {number} value - Variable value
     */
    createVariable(name, value) {
        this.setVariable(name, value);
    }

    /**
     * Remove a variable (resets to defaults)
     * @param {string} name - Variable name
     * @returns {boolean} Success status
     */
    removeVariable(name) {
        this.initializeDefaults();
        return true;
    }

    /**
     * Clone a variable value
     * @param {string} source - Source variable name
     * @param {string} dest - Destination variable name
     */
    cloneVariable(source, dest) {
        const value = this.getVariable(source);
        this.setVariable(dest, value);
    }

    /**
     * List all available variables
     * @returns {Array<string>} Array of variable names
     */
    listVariables() {
        return [
            'G', 'M', 'r', 'H0', 'B', 'B_crit', 'Lambda', 'c_light', 'q_charge', 'gas_v',
            'f_TRZ', 'E_0', 'tau_erosion', 'rho_wind', 'v_wind', 'rho_fluid',
            'rho_vac_UA', 'rho_vac_SCm', 'scale_EM', 'proton_mass',
            'hbar', 't_Hubble', 't_Hubble_gyr', 'delta_x', 'delta_p', 'integral_psi',
            'A_osc', 'k_osc', 'omega_osc', 'x_pos', 'M_DM_factor', 'delta_rho_over_rho'
        ];
    }

    /**
     * Get system name
     * @returns {string} System name
     */
    getSystemName() {
        return "HorseheadNebula";
    }

    // ----- Batch Operations (2 methods) -----

    /**
     * Transform a group of variables using a function
     * @param {Array<string>} vars - Variable names
     * @param {Function} func - Transformation function
     */
    transformVariableGroup(vars, func) {
        vars.forEach(varName => {
            const current = this.getVariable(varName);
            this.setVariable(varName, func(current));
        });
    }

    /**
     * Scale a group of variables by a factor
     * @param {Array<string>} vars - Variable names
     * @param {number} factor - Scaling factor
     */
    scaleVariableGroup(vars, factor) {
        this.transformVariableGroup(vars, v => v * factor);
    }

    // ----- Self-Expansion (4 methods - domain-specific for dark nebula) -----

    /**
     * Expand entire parameter space by a scale factor
     * @param {number} scale_factor - Scaling factor for all scalable parameters
     */
    expandParameterSpace(scale_factor) {
        const scalable = ['M', 'r', 'E_0', 'tau_erosion', 'rho_wind', 'v_wind', 'B'];
        this.scaleVariableGroup(scalable, scale_factor);
    }

    /**
     * Expand nebula mass and radius scales
     * @param {number} M_scale - Mass scaling factor
     * @param {number} r_scale - Radius scaling factor
     */
    expandNebulaScale(M_scale, r_scale) {
        this.setVariable('M', this.getVariable('M') * M_scale);
        this.setVariable('r', this.getVariable('r') * r_scale);
    }

    /**
     * Expand erosion scales (UNIQUE to dark nebula)
     * @param {number} E_0_scale - Erosion factor scaling
     * @param {number} tau_erosion_scale - Erosion timescale scaling
     */
    expandErosionScale(E_0_scale, tau_erosion_scale) {
        this.setVariable('E_0', this.getVariable('E_0') * E_0_scale);
        this.setVariable('tau_erosion', this.getVariable('tau_erosion') * tau_erosion_scale);
    }

    /**
     * Expand wind and magnetic field scales
     * @param {number} rho_wind_scale - Wind density scaling
     * @param {number} v_wind_scale - Wind velocity scaling
     * @param {number} B_scale - Magnetic field scaling
     */
    expandWindMagneticScale(rho_wind_scale, v_wind_scale, B_scale) {
        this.setVariable('rho_wind', this.getVariable('rho_wind') * rho_wind_scale);
        this.setVariable('v_wind', this.getVariable('v_wind') * v_wind_scale);
        this.setVariable('B', this.getVariable('B') * B_scale);
    }

    // ----- Self-Refinement (3 methods) -----

    /**
     * Auto-refine parameters based on observations
     * @param {Array<Array<number>>} observations - Array of [time, g] pairs
     */
    autoRefineParameters(observations) {
        this.calibrateToObservations(observations);
    }

    /**
     * Calibrate to observational data
     * @param {Array<Array<number>>} obs_data - Array of [time, g] pairs
     */
    calibrateToObservations(obs_data) {
        if (obs_data.length === 0) return;

        let total_error = 0.0;
        obs_data.forEach(([t, g_obs]) => {
            const g_model = this.compute_g_Horsehead(t);
            total_error += Math.abs(g_model - g_obs);
        });

        const avg_error = total_error / obs_data.length;
        if (avg_error > 1e-9) {
            const correction = 0.95;
            this.setVariable('E_0', this.getVariable('E_0') * correction);
        }
    }

    /**
     * Optimize for a given metric over time range
     * @param {Function} metric - Metric function (takes g, returns score)
     * @param {number} t_start - Start time (s)
     * @param {number} t_end - End time (s)
     * @param {number} steps - Number of time steps
     * @returns {number} Best metric value
     */
    optimizeForMetric(metric, t_start, t_end, steps) {
        let best_metric = -1e100;
        const dt = (t_end - t_start) / steps;

        for (let i = 0; i <= steps; i++) {
            const t = t_start + i * dt;
            const g = this.compute_g_Horsehead(t);
            const m = metric(g);
            if (m > best_metric) best_metric = m;
        }
        return best_metric;
    }

    // ----- Parameter Exploration (1 method) -----

    /**
     * Generate parameter variations for exploration
     * @param {number} count - Number of variations
     * @param {number} variation_percent - Variation percentage
     * @returns {Array<Object>} Array of parameter variation objects
     */
    generateVariations(count, variation_percent) {
        const variations = [];
        const vars = this.listVariables();

        for (let i = 0; i < count; i++) {
            const variant = {};
            vars.forEach(varName => {
                const base = this.getVariable(varName);
                const variation = base * (1.0 + (Math.random() * 2 - 1) * variation_percent / 100.0);
                variant[varName] = variation;
            });
            variations.push(variant);
        }
        return variations;
    }

    // ----- Adaptive Evolution (2 methods) -----

    /**
     * Mutate parameters randomly
     * @param {number} mutation_rate - Mutation rate (fractional)
     */
    mutateParameters(mutation_rate) {
        const vars = this.listVariables();
        vars.forEach(varName => {
            const current = this.getVariable(varName);
            const mutated = current * (1.0 + (Math.random() * 2 - 1) * mutation_rate);
            this.setVariable(varName, mutated);
        });
    }

    /**
     * Evolve system over generations using fitness function
     * @param {number} generations - Number of generations
     * @param {Function} fitness - Fitness function (takes HorseheadNebula instance)
     */
    evolveSystem(generations, fitness) {
        let best_fitness = fitness(this);
        this.saveState('evolution_best');

        for (let gen = 0; gen < generations; gen++) {
            this.saveState('evolution_temp');
            this.mutateParameters(0.1);

            const current_fitness = fitness(this);
            if (current_fitness > best_fitness) {
                best_fitness = current_fitness;
                this.saveState('evolution_best');
            } else {
                this.restoreState('evolution_temp');
            }
        }
        this.restoreState('evolution_best');
    }

    // ----- State Management (4 methods) -----

    /**
     * Save current state with a label
     * @param {string} label - State label
     */
    saveState(label) {
        const state = {};
        this.listVariables().forEach(varName => {
            state[varName] = this.getVariable(varName);
        });
        HorseheadNebula.savedStates[label] = state;
    }

    /**
     * Restore state from a label
     * @param {string} label - State label
     * @returns {boolean} Success status
     */
    restoreState(label) {
        const state = HorseheadNebula.savedStates[label];
        if (!state) return false;

        Object.entries(state).forEach(([varName, value]) => {
            this.setVariable(varName, value);
        });
        return true;
    }

    /**
     * List all saved states
     * @returns {Array<string>} Array of state labels
     */
    listSavedStates() {
        return Object.keys(HorseheadNebula.savedStates);
    }

    /**
     * Export current state as string
     * @returns {string} State export string
     */
    exportState() {
        let output = "HorseheadNebula State Export\n";
        output += "============================\n";

        this.listVariables().forEach(varName => {
            const value = this.getVariable(varName);
            output += `${varName}: ${value.toExponential(15)}\n`;
        });
        return output;
    }

    // ----- System Analysis (4 methods) -----

    /**
     * Perform sensitivity analysis at given time
     * @param {number} t - Time (s)
     * @param {number} perturbation - Perturbation fraction (e.g., 0.01 for 1%)
     * @returns {Object} Sensitivity map {varName: sensitivity}
     */
    sensitivityAnalysis(t, perturbation) {
        const sensitivities = {};
        const baseline = this.compute_g_Horsehead(t);

        this.listVariables().forEach(varName => {
            const original = this.getVariable(varName);
            this.setVariable(varName, original * (1.0 + perturbation));
            const perturbed = this.compute_g_Horsehead(t);
            sensitivities[varName] = Math.abs(perturbed - baseline) / baseline;
            this.setVariable(varName, original);
        });
        return sensitivities;
    }

    /**
     * Generate comprehensive system report
     * @param {number} t - Time (s)
     * @returns {string} Report string
     */
    generateReport(t) {
        const M_sun = 1.989e30;
        const ly_to_m = 9.461e15;
        const Myr_to_s = 1e6 * 3.156e7;

        let report = "=== Horsehead Nebula System Report ===\n";
        report += `Time: ${(t / Myr_to_s).toFixed(6)} Myr\n`;
        report += `E(t) erosion factor: ${this.E_t(t).toFixed(6)}\n`;
        report += `g_Horsehead: ${this.compute_g_Horsehead(t).toExponential(6)} m/s^2\n`;
        report += `Core parameters: M=${(this.M/M_sun).toFixed(0)} M_sun, r=${(this.r/ly_to_m).toFixed(1)} ly\n`;
        report += `Erosion: E_0=${this.E_0.toFixed(3)}, tau_erosion=${(this.tau_erosion/Myr_to_s).toFixed(0)} Myr\n`;
        report += `Wind: rho_wind=${this.rho_wind.toExponential(3)} kg/m^3, v_wind=${this.v_wind.toExponential(3)} m/s\n`;
        return report;
    }

    /**
     * Validate parameter consistency
     * @returns {boolean} Validation status
     */
    validateConsistency() {
        let valid = true;
        if (this.M <= 0 || this.r <= 0 || this.tau_erosion <= 0) valid = false;
        if (this.E_0 < 0 || this.E_0 > 1) valid = false;
        if (this.rho_wind < 0 || this.v_wind < 0) valid = false;
        return valid;
    }

    /**
     * Auto-correct anomalous parameters
     * @returns {boolean} Whether corrections were made
     */
    autoCorrectAnomalies() {
        let corrected = false;
        const M_sun = 1.989e30;
        const ly_to_m = 9.461e15;
        const Myr_to_s = 1e6 * 3.156e7;

        if (this.M <= 0) { this.M = 1000.0 * M_sun; corrected = true; }
        if (this.r <= 0) { this.r = 2.5 * ly_to_m; corrected = true; }
        if (this.tau_erosion <= 0) { this.tau_erosion = 5 * Myr_to_s; corrected = true; }
        if (this.E_0 < 0) { this.E_0 = 0.0; corrected = true; }
        if (this.E_0 > 1) { this.E_0 = 1.0; corrected = true; }
        if (this.rho_wind < 0) { this.rho_wind = 1e-21; corrected = true; }
        if (this.v_wind < 0) { this.v_wind = 2e6; corrected = true; }
        if (corrected) this.updateCache();
        return corrected;
    }
}

// ========== ENHANCED EXAMPLE FUNCTION ==========

function enhanced_horsehead_example() {
    console.log("=========================================================");
    console.log("ENHANCED HORSEHEAD NEBULA DEMONSTRATION");
    console.log("Barnard 33 Dark Nebula with Erosion");
    console.log("=========================================================\n");

    const horsehead = new HorseheadNebula();
    const M_sun = 1.989e30;
    const Myr_to_s = 1e6 * 3.156e7;

    // Step 1: Initial state and validation
    console.log("Step 1: Initial State and Validation");
    console.log(`System: ${horsehead.getSystemName()}`);
    console.log(`Validation: ${horsehead.validateConsistency() ? 'PASS' : 'FAIL'}`);
    console.log(`Auto-corrected: ${horsehead.autoCorrectAnomalies() ? 'Yes' : 'No'}\n`);

    // Step 2: Time evolution showing erosion E(t)
    console.log("Step 2: Time Evolution (Nebula Erosion E(t))");
    const t_Myr_array = [0.0, 1.0, 2.0, 3.0, 5.0];
    t_Myr_array.forEach(t_Myr => {
        const t = t_Myr * Myr_to_s;
        const Et = horsehead.E_t(t);
        const g = horsehead.compute_g_Horsehead(t);
        console.log(`  t = ${t_Myr.toFixed(1)} Myr: E(t) = ${Et.toFixed(6)}, g = ${g.toExponential(6)} m/s^2`);
    });
    console.log();

    // Step 3: Variable listing
    console.log("Step 3: Variable Listing");
    const vars = horsehead.listVariables();
    console.log(`Total variables: ${vars.length}`);
    console.log(`Sample: ${vars[0]}, ${vars[1]}, ${vars[11]} (E_0), ${vars[12]} (tau_erosion)\n`);

    // Step 4: Nebula mass scaling
    console.log("Step 4: Nebula Mass Scaling (M sweeps)");
    horsehead.saveState('original');
    const M_factors = [0.5, 1.0, 2.0];
    M_factors.forEach(factor => {
        horsehead.restoreState('original');
        horsehead.expandNebulaScale(factor, 1.0);
        const t = 3 * Myr_to_s;
        const g = horsehead.compute_g_Horsehead(t);
        const M = horsehead.getVariable('M');
        console.log(`  M × ${factor}: M = ${(M/M_sun).toFixed(0)} M_sun, g(3 Myr) = ${g.toExponential(6)} m/s^2`);
    });
    horsehead.restoreState('original');
    console.log();

    // Step 5: Erosion factor scaling (UNIQUE to dark nebula erosion)
    console.log("Step 5: Erosion Factor Scaling (E_0 sweeps) - DARK NEBULA EROSION FEATURE");
    const E_0_factors = [0.5, 1.0, 2.0];
    E_0_factors.forEach(factor => {
        horsehead.restoreState('original');
        horsehead.expandErosionScale(factor, 1.0);
        const t = 3 * Myr_to_s;
        const Et = horsehead.E_t(t);
        const g = horsehead.compute_g_Horsehead(t);
        console.log(`  E_0 × ${factor}: E(3 Myr) = ${Et.toFixed(6)}, g = ${g.toExponential(6)} m/s^2`);
    });
    horsehead.restoreState('original');
    console.log();

    // Step 6: Erosion timescale sweeps
    console.log("Step 6: Erosion Timescale Sweeps (tau_erosion)");
    const tau_erosion_factors = [0.5, 1.0, 2.0];
    tau_erosion_factors.forEach(factor => {
        horsehead.restoreState('original');
        horsehead.expandErosionScale(1.0, factor);
        const t = 3 * Myr_to_s;
        const Et = horsehead.E_t(t);
        console.log(`  tau_erosion × ${factor}: E(3 Myr) = ${Et.toFixed(6)}`);
    });
    horsehead.restoreState('original');
    console.log();

    // Step 7: Generate report
    console.log("Step 7: Full System Report at 3 Myr");
    const t_report = 3 * Myr_to_s;
    const report = horsehead.generateReport(t_report);
    console.log(report);

    console.log("=========================================================");
    console.log("ENHANCED DEMONSTRATION COMPLETE");
    console.log("=========================================================\n");
}

// Run inline test if executed directly
if (typeof require !== 'undefined' && require.main === module) {
    console.log("Running inline test...\n");
    enhanced_horsehead_example();
}

// Export for use as module
export default HorseheadNebula;
