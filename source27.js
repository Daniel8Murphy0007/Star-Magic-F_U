/**
 * ================================================================================================
 * Module: GalaxyNGC1792.js
 *
 * Description: JavaScript ES6 Module for NGC 1792 ("The Stellar Forge") Starburst Galaxy Class
 *              This is the nineteenth module in a series of 500+ code files for the Universal
 *              Quantum Field Framework (UQFF) simulations, focusing on starburst galaxy evolution
 *              and gravity equations derived from Hubble datasets, high-energy lab simulations,
 *              and UQFF refinements (dated May 09, 2025, updated October 08, 2025).
 *
 * Purpose: Encapsulates the Master Universal Gravity Equation (MUGE) for NGC 1792 evolution.
 *          Includes ALL terms: base gravity with star formation growth M(t), cosmic expansion (H(z)),
 *          magnetic correction (static B), UQFF Ug components with f_TRZ, Lambda, quantum uncertainty,
 *          scaled EM with [UA], fluid dynamics, oscillatory waves, DM/density perturbations, and
 *          supernova feedback (pressure / density for acc).
 *
 * Key Features:
 *   - Default values from UQFF document: M0 = 1×10¹⁰ M☉, r = 80,000 ly, z = 0.0095,
 *     Hz ≈ 2.19×10⁻¹⁸ s⁻¹, SFR_factor = 10/1×10¹⁰ (normalized), tau_SF = 100 Myr,
 *     B = 1×10⁻⁵ T, rho_wind = 1×10⁻²¹ kg/m³, v_wind = 2×10⁶ m/s.
 *   - M(t) = M0 × (1 + SFR_factor × exp(-t/tau_SF)) - mass growth with star formation
 *   - Units handled: M☉ to kg, ly to m; feedback term as (rho × v_wind²) / rho_fluid
 *   - Setter methods: setVariable, addToVariable, subtractFromVariable
 *   - Computes g_NGC1792(r, t) with every term explicitly included (10 terms total)
 *   - Enhanced dynamic capabilities: 38 methods total (13 core + 25 enhanced)
 *
 * Conversion: Ported from C++ source27.cpp (GalaxyNGC1792.h) to JavaScript ES6
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript
 * Date: November 03, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

class GalaxyNGC1792 {
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
        this.M0 = 1e10 * M_sun;               // Initial mass (kg)
        const ly_to_m = 9.461e15;             // Light year to meters
        this.r = 80000.0 * ly_to_m;           // Radius (m) - 80,000 ly
        this.z_gal = 0.0095;                  // Galaxy redshift
        
        // Hubble parameter at z
        const Hz_kms = 70 * Math.sqrt(0.3 * Math.pow(1 + this.z_gal, 3) + 0.7);  // km/s/Mpc
        this.Hz = (Hz_kms * 1000 / 3.086e22); // Convert to s^-1 (correct Mpc to m: 1 Mpc = 3.086e22 m)
        
        // Magnetic field (static)
        this.B = 1e-5;                        // Static magnetic field (T)
        this.B_crit = 1e11;                   // Critical B field (T)
        
        // Cosmological parameters
        this.Lambda = 1.1e-52;                // Cosmological constant (m^-2)
        this.c_light = 3e8;                   // Speed of light (m/s)
        
        // EM parameters
        this.q_charge = 1.602e-19;            // Proton charge (C)
        this.gas_v = 1e5;                     // Gas velocity for EM (m/s)
        this.f_TRZ = 0.1;                     // Time-reversal factor
        
        // Star formation parameters (UNIQUE to starburst)
        this.SFR_factor = 10.0 / 1e10;        // Normalized star formation rate factor
        this.tau_SF = 100e6 * 3.156e7;        // Star formation timescale (s) - 100 Myr
        
        // Wind/feedback parameters (UNIQUE to starburst with supernovae)
        this.rho_wind = 1e-21;                // Wind density (kg/m³)
        this.v_wind = 2e6;                    // Wind velocity (m/s) - supernova-driven
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
        this.ug1_base = (this.G * this.M0) / (this.r * this.r);
    }

    // ========== CORE PHYSICS METHODS (13 methods) ==========

    /**
     * Compute mass M(t) with star formation
     * M(t) = M0 × (1 + SFR_factor × exp(-t/tau_SF)) - mass growth
     * @param {number} t - Time (s)
     * @returns {number} Mass at time t (kg)
     */
    M_t(t) {
        const M_dot = this.SFR_factor * Math.exp(-t / this.tau_SF);
        return this.M0 * (1 + M_dot);
    }

    /**
     * Compute UQFF Ug components with B correction
     * @param {number} Mt - Mass at time t (kg)
     * @returns {number} Total Ug (m/s²)
     */
    compute_Ug(Mt) {
        const Ug1 = (this.G * Mt) / (this.r * this.r);
        const Ug2 = 0.0;
        const Ug3 = 0.0;
        const corr_B = 1 - this.B / this.B_crit;
        const Ug4 = Ug1 * corr_B;
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1 + this.f_TRZ);
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
     * @returns {number} Total gravitational acceleration g_NGC1792 (m/s²)
     */
    compute_g_NGC1792(t) {
        if (t < 0) {
            console.error("Error: Time t must be non-negative.");
            return 0.0;
        }

        const Mt = this.M_t(t);
        const ug1_t = (this.G * Mt) / (this.r * this.r);

        // Term 1: Base + Hz + B corrections
        const corr_H = 1 + this.Hz * t;
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

        // Term 5: Quantum uncertainty
        const sqrt_unc = Math.sqrt(this.delta_x * this.delta_p);
        const term_q = (this.hbar / sqrt_unc) * this.integral_psi * (2 * Math.PI / this.t_Hubble);

        // Term 6: Fluid dynamics (effective acceleration)
        const V = this.compute_V();
        const term_fluid = (this.rho_fluid * V * ug1_t) / Mt;

        // Term 7: Oscillatory waves (real parts)
        const term_osc1 = 2 * this.A_osc * Math.cos(this.k_osc * this.x_pos) * Math.cos(this.omega_osc * t);
        const arg = this.k_osc * this.x_pos - this.omega_osc * t;
        const term_osc2 = (2 * Math.PI / this.t_Hubble_gyr) * this.A_osc * Math.cos(arg);
        const term_osc = term_osc1 + term_osc2;

        // Term 8: Dark matter and density perturbation (converted to acceleration)
        const M_dm = Mt * this.M_DM_factor;
        const pert1 = this.delta_rho_over_rho;
        const pert2 = 3 * this.G * Mt / (this.r * this.r * this.r);
        const term_dm_force_like = (Mt + M_dm) * (pert1 + pert2);
        const term_DM = term_dm_force_like / Mt;

        // Term 9: Supernova feedback (pressure / density for acceleration)
        const wind_pressure = this.rho_wind * this.v_wind * this.v_wind;
        const term_feedback = wind_pressure / this.rho_fluid;

        // Total g_NGC1792 (all 9 terms summed - note: term_osc counts as 1)
        return term1 + term2 + term3 + term4 + term_q + term_fluid + term_osc + term_DM + term_feedback;
    }

    /**
     * Universal getter for any variable by name
     * @param {string} varName - Variable name
     * @returns {number} Variable value
     */
    getVariable(varName) {
        const validVars = [
            'G', 'M0', 'r', 'Hz', 'B', 'B_crit', 'Lambda', 'c_light', 'q_charge',
            'gas_v', 'f_TRZ', 'SFR_factor', 'tau_SF', 'rho_wind', 'v_wind',
            'rho_fluid', 'rho_vac_UA', 'rho_vac_SCm', 'scale_EM', 'proton_mass',
            'z_gal', 'hbar', 't_Hubble', 't_Hubble_gyr', 'delta_x', 'delta_p',
            'integral_psi', 'A_osc', 'k_osc', 'omega_osc', 'x_pos', 'M_DM_factor',
            'delta_rho_over_rho'
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
            'G', 'M0', 'r', 'Hz', 'B', 'B_crit', 'Lambda', 'c_light', 'q_charge',
            'gas_v', 'f_TRZ', 'SFR_factor', 'tau_SF', 'rho_wind', 'v_wind',
            'rho_fluid', 'rho_vac_UA', 'rho_vac_SCm', 'scale_EM', 'proton_mass',
            'z_gal', 'hbar', 't_Hubble', 't_Hubble_gyr', 'delta_x', 'delta_p',
            'integral_psi', 'A_osc', 'k_osc', 'omega_osc', 'x_pos', 'M_DM_factor',
            'delta_rho_over_rho'
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
     * Example computation at t=50 Myr (for testing)
     * @returns {number} g_NGC1792 at 50 Myr (m/s²)
     */
    exampleAt50Myr() {
        const t_example = 50e6 * 3.156e7;
        return this.compute_g_NGC1792(t_example);
    }

    /**
     * Print parameters to console
     */
    printParameters() {
        console.log("NGC 1792 Parameters:");
        console.log(`G: ${this.G.toExponential(3)}, M0: ${this.M0.toExponential(3)}, r: ${this.r.toExponential(3)}`);
        console.log(`Hz: ${this.Hz.toExponential(3)}, B: ${this.B.toExponential(3)}, B_crit: ${this.B_crit.toExponential(3)}`);
        console.log(`f_TRZ: ${this.f_TRZ.toFixed(3)}, SFR_factor: ${this.SFR_factor.toExponential(3)}, tau_SF: ${this.tau_SF.toExponential(3)}`);
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
            'G', 'M0', 'r', 'Hz', 'B', 'B_crit', 'Lambda', 'c_light', 'q_charge',
            'gas_v', 'f_TRZ', 'SFR_factor', 'tau_SF', 'rho_wind', 'v_wind',
            'rho_fluid', 'rho_vac_UA', 'rho_vac_SCm', 'scale_EM', 'proton_mass',
            'z_gal', 'hbar', 't_Hubble', 't_Hubble_gyr', 'delta_x', 'delta_p',
            'integral_psi', 'A_osc', 'k_osc', 'omega_osc', 'x_pos', 'M_DM_factor',
            'delta_rho_over_rho'
        ];
    }

    /**
     * Get system name
     * @returns {string} System name
     */
    getSystemName() {
        return "GalaxyNGC1792";
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

    // ----- Self-Expansion (4 methods - domain-specific for starburst galaxy) -----

    /**
     * Expand entire parameter space by a scale factor
     * @param {number} scale_factor - Scaling factor for all scalable parameters
     */
    expandParameterSpace(scale_factor) {
        const scalable = ['M0', 'r', 'SFR_factor', 'tau_SF', 'rho_wind', 'v_wind', 'B'];
        this.scaleVariableGroup(scalable, scale_factor);
    }

    /**
     * Expand galaxy mass and radius scales
     * @param {number} M0_scale - Mass scaling factor
     * @param {number} r_scale - Radius scaling factor
     */
    expandGalaxyScale(M0_scale, r_scale) {
        this.setVariable('M0', this.getVariable('M0') * M0_scale);
        this.setVariable('r', this.getVariable('r') * r_scale);
    }

    /**
     * Expand starburst scales (UNIQUE to starburst galaxy)
     * @param {number} SFR_factor_scale - SFR factor scaling
     * @param {number} tau_SF_scale - SF timescale scaling
     */
    expandStarburstScale(SFR_factor_scale, tau_SF_scale) {
        this.setVariable('SFR_factor', this.getVariable('SFR_factor') * SFR_factor_scale);
        this.setVariable('tau_SF', this.getVariable('tau_SF') * tau_SF_scale);
    }

    /**
     * Expand wind and magnetic scales (UNIQUE to supernova-driven starburst)
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
            const g_model = this.compute_g_NGC1792(t);
            total_error += Math.abs(g_model - g_obs);
        });

        const avg_error = total_error / obs_data.length;
        if (avg_error > 1e-9) {
            const correction = 0.95;
            this.setVariable('SFR_factor', this.getVariable('SFR_factor') * correction);
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
            const g = this.compute_g_NGC1792(t);
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
     * @param {Function} fitness - Fitness function (takes GalaxyNGC1792 instance)
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
        GalaxyNGC1792.savedStates[label] = state;
    }

    /**
     * Restore state from a label
     * @param {string} label - State label
     * @returns {boolean} Success status
     */
    restoreState(label) {
        const state = GalaxyNGC1792.savedStates[label];
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
        return Object.keys(GalaxyNGC1792.savedStates);
    }

    /**
     * Export current state as string
     * @returns {string} State export string
     */
    exportState() {
        let output = "GalaxyNGC1792 State Export\n";
        output += "===========================\n";

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
        const baseline = this.compute_g_NGC1792(t);

        this.listVariables().forEach(varName => {
            const original = this.getVariable(varName);
            this.setVariable(varName, original * (1.0 + perturbation));
            const perturbed = this.compute_g_NGC1792(t);
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

        let report = "=== NGC 1792 System Report ===\n";
        report += `Time: ${(t / Myr_to_s).toFixed(6)} Myr\n`;
        report += `M(t): ${(this.M_t(t) / M_sun).toExponential(6)} M_sun\n`;
        report += `g_NGC1792: ${this.compute_g_NGC1792(t).toExponential(6)} m/s^2\n`;
        report += `Core parameters: M0=${(this.M0/M_sun).toExponential(3)} M_sun, r=${(this.r/ly_to_m).toFixed(0)} ly\n`;
        report += `Starburst: SFR_factor=${this.SFR_factor.toExponential(3)}, tau_SF=${(this.tau_SF/Myr_to_s).toFixed(0)} Myr\n`;
        report += `Wind: rho_wind=${this.rho_wind.toExponential(3)} kg/m^3, v_wind=${this.v_wind.toExponential(3)} m/s\n`;
        report += `Magnetic: B=${this.B.toExponential(3)} T\n`;
        return report;
    }

    /**
     * Validate parameter consistency
     * @returns {boolean} Validation status
     */
    validateConsistency() {
        let valid = true;
        if (this.M0 <= 0 || this.r <= 0 || this.tau_SF <= 0) valid = false;
        if (this.SFR_factor < 0) valid = false;
        if (this.rho_wind < 0 || this.v_wind < 0 || this.B < 0) valid = false;
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

        if (this.M0 <= 0) { this.M0 = 1e10 * M_sun; corrected = true; }
        if (this.r <= 0) { this.r = 80000 * ly_to_m; corrected = true; }
        if (this.tau_SF <= 0) { this.tau_SF = 100 * Myr_to_s; corrected = true; }
        if (this.SFR_factor < 0) { this.SFR_factor = 10.0 / 1e10; corrected = true; }
        if (this.rho_wind < 0) { this.rho_wind = 1e-21; corrected = true; }
        if (this.v_wind < 0) { this.v_wind = 2e6; corrected = true; }
        if (this.B < 0) { this.B = 1e-5; corrected = true; }
        if (corrected) this.updateCache();
        return corrected;
    }
}

// ========== ENHANCED EXAMPLE FUNCTION ==========

function enhanced_ngc1792_example() {
    console.log("=========================================================");
    console.log("ENHANCED NGC 1792 DEMONSTRATION");
    console.log("The Stellar Forge - Starburst Galaxy");
    console.log("=========================================================\n");

    const ngc = new GalaxyNGC1792();
    const M_sun = 1.989e30;
    const Myr_to_s = 1e6 * 3.156e7;

    // Step 1: Initial state and validation
    console.log("Step 1: Initial State and Validation");
    console.log(`System: ${ngc.getSystemName()}`);
    console.log(`Validation: ${ngc.validateConsistency() ? 'PASS' : 'FAIL'}`);
    console.log(`Auto-corrected: ${ngc.autoCorrectAnomalies() ? 'Yes' : 'No'}\n`);

    // Step 2: Time evolution showing M(t)
    console.log("Step 2: Time Evolution (Starburst Mass Growth M(t))");
    const t_Myr_array = [0.0, 25.0, 50.0, 100.0, 200.0];
    t_Myr_array.forEach(t_Myr => {
        const t = t_Myr * Myr_to_s;
        const Mt = ngc.M_t(t);
        const g = ngc.compute_g_NGC1792(t);
        console.log(`  t = ${t_Myr.toFixed(1)} Myr: M(t) = ${(Mt/M_sun).toExponential(4)} M_sun, g = ${g.toExponential(4)} m/s^2`);
    });
    console.log();

    // Step 3: Starburst rate scaling (UNIQUE to starburst galaxy)
    console.log("Step 3: Starburst Rate Scaling (SFR_factor sweeps) - STARBURST");
    ngc.saveState('original');
    const SFR_factors = [0.5, 1.0, 2.0];
    SFR_factors.forEach(factor => {
        ngc.restoreState('original');
        ngc.expandStarburstScale(factor, 1.0);
        const t = 50 * Myr_to_s;
        const Mt = ngc.M_t(t);
        console.log(`  SFR_factor × ${factor}: M(50 Myr) = ${(Mt/M_sun).toExponential(4)} M_sun`);
    });
    ngc.restoreState('original');
    console.log();

    // Step 4: Supernova wind velocity scaling
    console.log("Step 4: Supernova Wind Velocity Scaling (v_wind sweeps)");
    const v_wind_factors = [0.5, 1.0, 2.0];
    v_wind_factors.forEach(factor => {
        ngc.restoreState('original');
        ngc.expandWindMagneticScale(1.0, factor, 1.0);
        const t = 50 * Myr_to_s;
        const g = ngc.compute_g_NGC1792(t);
        console.log(`  v_wind × ${factor}: g(50 Myr) = ${g.toExponential(6)} m/s^2`);
    });
    ngc.restoreState('original');
    console.log();

    // Step 5: Generate report
    console.log("Step 5: Full System Report at 50 Myr");
    const t_report = 50 * Myr_to_s;
    const report = ngc.generateReport(t_report);
    console.log(report);

    console.log("=========================================================");
    console.log("ENHANCED DEMONSTRATION COMPLETE");
    console.log("=========================================================\n");
}

// Run inline test if executed directly
if (typeof require !== 'undefined' && require.main === module) {
    console.log("Running inline test...\n");
    enhanced_ngc1792_example();
}

// Export for use as module
export default GalaxyNGC1792;
