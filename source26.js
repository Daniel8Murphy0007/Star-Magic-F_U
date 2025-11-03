/**
 * ================================================================================================
 * Module: HUDFGalaxies.js
 *
 * Description: JavaScript ES6 Module for Hubble Ultra Deep Field (HUDF) "Galaxies Galore" Class
 *              This is the eighteenth module in a series of 500+ code files for the Universal
 *              Quantum Field Framework (UQFF) simulations, focusing on cosmic field of galaxies
 *              evolution and gravity equations derived from Hubble datasets, high-energy lab
 *              simulations, and UQFF refinements (dated May 09, 2025, updated October 08, 2025).
 *
 * Purpose: Encapsulates the Master Universal Gravity Equation (MUGE) for HUDF galaxies evolution.
 *          Includes ALL terms: base gravity with formation growth M(t), cosmic expansion (H(z)),
 *          magnetic correction (static B), interactions I(t), UQFF Ug components with f_TRZ,
 *          Lambda, quantum uncertainty, scaled EM with [UA], fluid dynamics, oscillatory waves,
 *          DM/density perturbations, and merger feedback.
 *
 * Key Features:
 *   - Default values from UQFF document: M0 = 1×10¹² M☉ (representative field mass),
 *     r = 1.3×10¹¹ ly (cosmic scale), z_avg = 3.5, Hz ≈ 2.5×10⁻¹⁸ s⁻¹, SFR_factor = 1.0,
 *     tau_SF = 1 Gyr, I0 = 0.05, tau_inter = 1 Gyr, rho_wind = 1×10⁻²² kg/m³,
 *     v_wind = 1×10⁶ m/s, B = 1×10⁻¹⁰ T.
 *   - M(t) = M0 × (1 + SFR_factor × exp(-t/tau_SF)) - mass growth with star formation
 *   - I(t) = I0 × exp(-t/tau_inter) - interaction exponential decay
 *   - Units handled: M☉ to kg, ly to m; interaction term I(t) scales gravity
 *   - Setter methods: setVariable, addToVariable, subtractFromVariable
 *   - Computes g_HUDF(r, t) with every term explicitly included (10 terms total)
 *   - Enhanced dynamic capabilities: 38 methods total (13 core + 25 enhanced)
 *
 * Conversion: Ported from C++ source26.cpp (HUDFGalaxies.h) to JavaScript ES6
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript
 * Date: November 03, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

class HUDFGalaxies {
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
        this.M0 = 1e12 * M_sun;               // Initial field mass (kg)
        const ly_to_m = 9.461e15;             // Light year to meters
        this.r = 1.3e11 * ly_to_m;            // Cosmic scale radius (m) - 130 billion ly
        this.z_avg = 3.5;                     // Average redshift
        
        // Hubble parameter at z_avg
        const Hz_kms = 70 * Math.sqrt(0.3 * Math.pow(1 + this.z_avg, 3) + 0.7);  // km/s/Mpc
        this.Hz = (Hz_kms * 1000 / 3.086e19); // Convert to s^-1
        
        // Magnetic field (static)
        this.B = 1e-10;                       // Static magnetic field (T)
        this.B_crit = 1e11;                   // Critical B field (T)
        
        // Cosmological parameters
        this.Lambda = 1.1e-52;                // Cosmological constant (m^-2)
        this.c_light = 3e8;                   // Speed of light (m/s)
        
        // EM parameters
        this.q_charge = 1.602e-19;            // Proton charge (C)
        this.gas_v = 1e5;                     // Gas velocity for EM (m/s)
        this.f_TRZ = 0.1;                     // Time-reversal factor
        
        // Star formation parameters (UNIQUE to cosmic field)
        this.SFR_factor = 1.0;                // Star formation rate factor
        this.tau_SF = 1e9 * 3.156e7;          // Star formation timescale (s) - 1 Gyr
        
        // Interaction parameters (UNIQUE to galaxy field)
        this.I0 = 0.05;                       // Initial interaction factor
        this.tau_inter = 1e9 * 3.156e7;       // Interaction timescale (s) - 1 Gyr
        
        // Wind/feedback parameters
        this.rho_wind = 1e-22;                // Wind density (kg/m³)
        this.v_wind = 1e6;                    // Wind velocity (m/s)
        this.rho_fluid = 1e-22;               // Fluid density (kg/m³)
        
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
        this.A_osc = 1e-12;                   // Oscillatory amplitude (m/s²)
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
     * Compute interaction factor I(t)
     * I(t) = I0 × exp(-t/tau_inter) - exponential decay
     * @param {number} t - Time (s)
     * @returns {number} Interaction factor (dimensionless)
     */
    I_t(t) {
        return this.I0 * Math.exp(-t / this.tau_inter);
    }

    /**
     * Compute UQFF Ug components with B and I corrections
     * @param {number} Mt - Mass at time t (kg)
     * @param {number} It - Interaction factor at time t
     * @returns {number} Total Ug (m/s²)
     */
    compute_Ug(Mt, It) {
        const Ug1 = (this.G * Mt) / (this.r * this.r);
        const Ug2 = 0.0;
        const Ug3 = 0.0;
        const corr_B = 1 - this.B / this.B_crit;
        const Ug4 = Ug1 * corr_B;
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1 + this.f_TRZ) * (1 + It);
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
     * @returns {number} Total gravitational acceleration g_HUDF (m/s²)
     */
    compute_g_HUDF(t) {
        if (t < 0) {
            console.error("Error: Time t must be non-negative.");
            return 0.0;
        }

        const Mt = this.M_t(t);
        const It = this.I_t(t);
        const ug1_t = (this.G * Mt) / (this.r * this.r);

        // Term 1: Base + Hz + B + I corrections
        const corr_H = 1 + this.Hz * t;
        const corr_B = 1 - this.B / this.B_crit;
        const corr_I = 1 + It;
        const term1 = ug1_t * corr_H * corr_B * corr_I;

        // Term 2: UQFF Ug with f_TRZ and I
        const term2 = this.compute_Ug(Mt, It);

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

        // Term 9: Merger feedback (pressure / density for acceleration)
        const wind_pressure = this.rho_wind * this.v_wind * this.v_wind;
        const term_feedback = wind_pressure / this.rho_fluid;

        // Total g_HUDF (all 9 terms summed - note: term_osc counts as 1)
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
            'gas_v', 'f_TRZ', 'SFR_factor', 'tau_SF', 'I0', 'tau_inter', 'rho_wind',
            'v_wind', 'rho_fluid', 'rho_vac_UA', 'rho_vac_SCm', 'scale_EM', 'proton_mass',
            'z_avg', 'hbar', 't_Hubble', 't_Hubble_gyr', 'delta_x', 'delta_p',
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
            'gas_v', 'f_TRZ', 'SFR_factor', 'tau_SF', 'I0', 'tau_inter', 'rho_wind',
            'v_wind', 'rho_fluid', 'rho_vac_UA', 'rho_vac_SCm', 'scale_EM', 'proton_mass',
            'z_avg', 'hbar', 't_Hubble', 't_Hubble_gyr', 'delta_x', 'delta_p',
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
     * Example computation at t=5 Gyr (for testing)
     * @returns {number} g_HUDF at 5 Gyr (m/s²)
     */
    exampleAt5Gyr() {
        const t_example = 5e9 * 3.156e7;
        return this.compute_g_HUDF(t_example);
    }

    /**
     * Print parameters to console
     */
    printParameters() {
        console.log("HUDF Galaxies Parameters:");
        console.log(`G: ${this.G.toExponential(3)}, M0: ${this.M0.toExponential(3)}, r: ${this.r.toExponential(3)}`);
        console.log(`Hz: ${this.Hz.toExponential(3)}, B: ${this.B.toExponential(3)}, B_crit: ${this.B_crit.toExponential(3)}`);
        console.log(`f_TRZ: ${this.f_TRZ.toFixed(3)}, SFR_factor: ${this.SFR_factor.toFixed(3)}, tau_SF: ${this.tau_SF.toExponential(3)}`);
        console.log(`I0: ${this.I0.toFixed(3)}, tau_inter: ${this.tau_inter.toExponential(3)}`);
        console.log(`rho_fluid: ${this.rho_fluid.toExponential(3)}, rho_wind: ${this.rho_wind.toExponential(3)}, v_wind: ${this.v_wind.toExponential(3)}`);
        console.log(`gas_v: ${this.gas_v.toExponential(3)}, M_DM_factor: ${this.M_DM_factor.toFixed(3)}`);
        console.log(`A_osc: ${this.A_osc.toExponential(3)}, delta_rho_over_rho: ${this.delta_rho_over_rho.toExponential(3)}`);
        console.log(`z_avg: ${this.z_avg.toFixed(2)}, ug1_base: ${this.ug1_base.toExponential(3)}`);
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
            'gas_v', 'f_TRZ', 'SFR_factor', 'tau_SF', 'I0', 'tau_inter', 'rho_wind',
            'v_wind', 'rho_fluid', 'rho_vac_UA', 'rho_vac_SCm', 'scale_EM', 'proton_mass',
            'z_avg', 'hbar', 't_Hubble', 't_Hubble_gyr', 'delta_x', 'delta_p',
            'integral_psi', 'A_osc', 'k_osc', 'omega_osc', 'x_pos', 'M_DM_factor',
            'delta_rho_over_rho'
        ];
    }

    /**
     * Get system name
     * @returns {string} System name
     */
    getSystemName() {
        return "HUDFGalaxies";
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

    // ----- Self-Expansion (4 methods - domain-specific for cosmic field) -----

    /**
     * Expand entire parameter space by a scale factor
     * @param {number} scale_factor - Scaling factor for all scalable parameters
     */
    expandParameterSpace(scale_factor) {
        const scalable = ['M0', 'r', 'I0', 'tau_inter', 'SFR_factor', 'tau_SF',
                          'rho_wind', 'v_wind', 'B'];
        this.scaleVariableGroup(scalable, scale_factor);
    }

    /**
     * Expand cosmic field mass and radius scales
     * @param {number} M0_scale - Mass scaling factor
     * @param {number} r_scale - Radius scaling factor
     */
    expandCosmicFieldScale(M0_scale, r_scale) {
        this.setVariable('M0', this.getVariable('M0') * M0_scale);
        this.setVariable('r', this.getVariable('r') * r_scale);
    }

    /**
     * Expand star formation scales (UNIQUE to cosmic field)
     * @param {number} SFR_factor_scale - SFR factor scaling
     * @param {number} tau_SF_scale - SF timescale scaling
     */
    expandStarFormationScale(SFR_factor_scale, tau_SF_scale) {
        this.setVariable('SFR_factor', this.getVariable('SFR_factor') * SFR_factor_scale);
        this.setVariable('tau_SF', this.getVariable('tau_SF') * tau_SF_scale);
    }

    /**
     * Expand interaction scales (UNIQUE to galaxy field)
     * @param {number} I0_scale - Interaction factor scaling
     * @param {number} tau_inter_scale - Interaction timescale scaling
     */
    expandInteractionScale(I0_scale, tau_inter_scale) {
        this.setVariable('I0', this.getVariable('I0') * I0_scale);
        this.setVariable('tau_inter', this.getVariable('tau_inter') * tau_inter_scale);
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
            const g_model = this.compute_g_HUDF(t);
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
            const g = this.compute_g_HUDF(t);
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
     * @param {Function} fitness - Fitness function (takes HUDFGalaxies instance)
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
        HUDFGalaxies.savedStates[label] = state;
    }

    /**
     * Restore state from a label
     * @param {string} label - State label
     * @returns {boolean} Success status
     */
    restoreState(label) {
        const state = HUDFGalaxies.savedStates[label];
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
        return Object.keys(HUDFGalaxies.savedStates);
    }

    /**
     * Export current state as string
     * @returns {string} State export string
     */
    exportState() {
        let output = "HUDFGalaxies State Export\n";
        output += "=========================\n";

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
        const baseline = this.compute_g_HUDF(t);

        this.listVariables().forEach(varName => {
            const original = this.getVariable(varName);
            this.setVariable(varName, original * (1.0 + perturbation));
            const perturbed = this.compute_g_HUDF(t);
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
        const Gyr_to_s = 1e9 * 3.156e7;

        let report = "=== HUDF Galaxies System Report ===\n";
        report += `Time: ${(t / Gyr_to_s).toFixed(6)} Gyr\n`;
        report += `M(t): ${(this.M_t(t) / M_sun).toExponential(6)} M_sun\n`;
        report += `I(t) interaction factor: ${this.I_t(t).toFixed(6)}\n`;
        report += `g_HUDF: ${this.compute_g_HUDF(t).toExponential(6)} m/s^2\n`;
        report += `Core parameters: M0=${(this.M0/M_sun).toExponential(3)} M_sun, r=${(this.r/ly_to_m).toExponential(3)} ly\n`;
        report += `Cosmic field: z_avg=${this.z_avg.toFixed(2)}, Hz=${this.Hz.toExponential(3)} s^-1\n`;
        report += `Star formation: SFR_factor=${this.SFR_factor.toFixed(3)}, tau_SF=${(this.tau_SF/Gyr_to_s).toFixed(2)} Gyr\n`;
        report += `Interactions: I0=${this.I0.toFixed(3)}, tau_inter=${(this.tau_inter/Gyr_to_s).toFixed(2)} Gyr\n`;
        return report;
    }

    /**
     * Validate parameter consistency
     * @returns {boolean} Validation status
     */
    validateConsistency() {
        let valid = true;
        if (this.M0 <= 0 || this.r <= 0 || this.tau_SF <= 0 || this.tau_inter <= 0) valid = false;
        if (this.I0 < 0 || this.I0 > 1) valid = false;
        if (this.SFR_factor < 0) valid = false;
        if (this.z_avg < 0) valid = false;
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
        const Gyr_to_s = 1e9 * 3.156e7;

        if (this.M0 <= 0) { this.M0 = 1e12 * M_sun; corrected = true; }
        if (this.r <= 0) { this.r = 1.3e11 * ly_to_m; corrected = true; }
        if (this.tau_SF <= 0) { this.tau_SF = Gyr_to_s; corrected = true; }
        if (this.tau_inter <= 0) { this.tau_inter = Gyr_to_s; corrected = true; }
        if (this.I0 < 0) { this.I0 = 0.0; corrected = true; }
        if (this.I0 > 1) { this.I0 = 1.0; corrected = true; }
        if (this.SFR_factor < 0) { this.SFR_factor = 1.0; corrected = true; }
        if (this.z_avg < 0) { this.z_avg = 3.5; corrected = true; }
        if (corrected) this.updateCache();
        return corrected;
    }
}

// ========== ENHANCED EXAMPLE FUNCTION ==========

function enhanced_hudf_example() {
    console.log("=========================================================");
    console.log("ENHANCED HUDF GALAXIES DEMONSTRATION");
    console.log("Hubble Ultra Deep Field - Cosmic Field of Galaxies");
    console.log("=========================================================\n");

    const hudf = new HUDFGalaxies();
    const M_sun = 1.989e30;
    const Gyr_to_s = 1e9 * 3.156e7;
    const ly_to_m = 9.461e15;

    // Step 1: Initial state and validation
    console.log("Step 1: Initial State and Validation");
    console.log(`System: ${hudf.getSystemName()}`);
    console.log(`Validation: ${hudf.validateConsistency() ? 'PASS' : 'FAIL'}`);
    console.log(`Auto-corrected: ${hudf.autoCorrectAnomalies() ? 'Yes' : 'No'}\n`);

    // Step 2: Time evolution showing M(t) and I(t)
    console.log("Step 2: Time Evolution (Mass M(t) and Interaction I(t))");
    const t_Gyr_array = [0.0, 1.0, 2.0, 5.0, 10.0];
    t_Gyr_array.forEach(t_Gyr => {
        const t = t_Gyr * Gyr_to_s;
        const Mt = hudf.M_t(t);
        const It = hudf.I_t(t);
        const g = hudf.compute_g_HUDF(t);
        console.log(`  t = ${t_Gyr.toFixed(1)} Gyr: M(t) = ${(Mt/M_sun).toExponential(4)} M_sun, I(t) = ${It.toFixed(6)}, g = ${g.toExponential(4)} m/s^2`);
    });
    console.log();

    // Step 3: Variable listing
    console.log("Step 3: Variable Listing");
    const vars = hudf.listVariables();
    console.log(`Total variables: ${vars.length}`);
    console.log(`Sample: ${vars[0]}, ${vars[1]}, ${vars[13]} (I0), ${vars[22]} (z_avg)\n`);

    // Step 4: Cosmic field mass scaling
    console.log("Step 4: Cosmic Field Mass Scaling (M0 sweeps)");
    hudf.saveState('original');
    const M_factors = [0.5, 1.0, 2.0];
    M_factors.forEach(factor => {
        hudf.restoreState('original');
        hudf.expandCosmicFieldScale(factor, 1.0);
        const t = 5 * Gyr_to_s;
        const g = hudf.compute_g_HUDF(t);
        const M = hudf.getVariable('M0');
        console.log(`  M0 × ${factor}: M0 = ${(M/M_sun).toExponential(3)} M_sun, g(5 Gyr) = ${g.toExponential(6)} m/s^2`);
    });
    hudf.restoreState('original');
    console.log();

    // Step 5: Cosmic field radius scaling (UNIQUE to HUDF)
    console.log("Step 5: Cosmic Field Radius Scaling (r sweeps) - COSMIC SCALE");
    const r_factors = [0.5, 1.0, 2.0];
    r_factors.forEach(factor => {
        hudf.restoreState('original');
        hudf.expandCosmicFieldScale(1.0, factor);
        const t = 5 * Gyr_to_s;
        const g = hudf.compute_g_HUDF(t);
        const r = hudf.getVariable('r');
        console.log(`  r × ${factor}: r = ${(r/ly_to_m).toExponential(3)} ly, g = ${g.toExponential(6)} m/s^2`);
    });
    hudf.restoreState('original');
    console.log();

    // Step 6: Star formation scaling (UNIQUE to cosmic field)
    console.log("Step 6: Star Formation Rate Scaling (SFR_factor sweeps)");
    const SFR_factors = [0.5, 1.0, 2.0];
    SFR_factors.forEach(factor => {
        hudf.restoreState('original');
        hudf.expandStarFormationScale(factor, 1.0);
        const t = 5 * Gyr_to_s;
        const Mt = hudf.M_t(t);
        console.log(`  SFR_factor × ${factor}: M(5 Gyr) = ${(Mt/M_sun).toExponential(4)} M_sun`);
    });
    hudf.restoreState('original');
    console.log();

    // Step 7: Interaction scaling (UNIQUE to galaxy field)
    console.log("Step 7: Galaxy Interaction Scaling (I0 sweeps)");
    const I0_factors = [0.5, 1.0, 2.0];
    I0_factors.forEach(factor => {
        hudf.restoreState('original');
        hudf.expandInteractionScale(factor, 1.0);
        const t = 5 * Gyr_to_s;
        const It = hudf.I_t(t);
        const g = hudf.compute_g_HUDF(t);
        console.log(`  I0 × ${factor}: I(5 Gyr) = ${It.toFixed(6)}, g = ${g.toExponential(6)} m/s^2`);
    });
    hudf.restoreState('original');
    console.log();

    // Step 8: Generate report
    console.log("Step 8: Full System Report at 5 Gyr");
    const t_report = 5 * Gyr_to_s;
    const report = hudf.generateReport(t_report);
    console.log(report);

    console.log("=========================================================");
    console.log("ENHANCED DEMONSTRATION COMPLETE");
    console.log("=========================================================\n");
}

// Run inline test if executed directly
if (typeof require !== 'undefined' && require.main === module) {
    console.log("Running inline test...\n");
    enhanced_hudf_example();
}

// Export for use as module
export default HUDFGalaxies;
