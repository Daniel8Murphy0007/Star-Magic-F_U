/**
 * ================================================================================================
 * Module: AntennaeGalaxies.js
 *
 * Description: JavaScript ES6 Module for Antennae Galaxies (NGC 4038 & NGC 4039) Class
 *              This is the fourteenth module in a series of 500+ code files for the Universal
 *              Quantum Field Framework (UQFF) simulations, focusing on interacting galaxy merger
 *              evolution and gravity equations derived from Hubble datasets, high-energy lab
 *              simulations, and UQFF refinements (dated May 09, 2025, updated October 08, 2025).
 *
 * Purpose: Encapsulates the Master Universal Gravity Equation (MUGE) for Antennae Galaxies merger.
 *          Includes ALL terms: base gravity with star formation growth M(t), cosmic expansion H(z),
 *          magnetic correction (static B), merger interaction I(t), UQFF Ug components with f_TRZ,
 *          Lambda, quantum uncertainty, scaled EM with [UA], fluid dynamics, oscillatory waves,
 *          DM/density perturbations, and stellar feedback (pressure / density for acc).
 *          Supports dynamic variable updates for all parameters.
 *
 * Key Features:
 *   - Default values from UQFF document: M0 = 2e11 M☉, r = 30,000 ly (2.838e20 m), z = 0.0105,
 *     Hz ≈ 2.19e-18 s^-1, SFR_factor = 20/(2e11), tau_SF = 100 Myr, I0 = 0.1, tau_merger = 400 Myr,
 *     rho_wind = 1e-21 kg/m³, v_wind = 2e6 m/s, B = 1e-5 T.
 *   - M(t) = M0 × (1 + SFR_factor × exp(-t/tau_SF)) - star formation growth
 *   - I(t) = I0 × exp(-t/tau_merger) - merger interaction exponential decay
 *   - Units handled: M☉ to kg, ly to m; interaction term I(t) scales gravity
 *   - Setter methods for updates: setVariable(varName, newValue), addToVariable, subtractFromVariable
 *   - Computes g_Antennae(r, t) with every term explicitly included (10 terms total)
 *   - Enhanced dynamic capabilities: 38 methods total (13 core + 25 enhanced)
 *
 * Conversion: Ported from C++ source23.cpp (AntennaeGalaxies.h) to JavaScript ES6
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript
 * Date: November 03, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

class AntennaeGalaxies {
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
        this.M0 = 2e11 * M_sun;               // Initial combined mass (kg)
        const ly_to_m = 9.461e15;             // Light year to meters
        this.r = 30000.0 * ly_to_m;           // Separation (m) - 30,000 ly
        this.z_gal = 0.0105;                  // Galaxy redshift
        
        // Hubble parameter at z
        const Hz_kms = 70 * Math.sqrt(0.3 * Math.pow(1 + this.z_gal, 3) + 0.7);  // km/s/Mpc
        this.Hz = (Hz_kms * 1000 / 3.086e19); // Convert to s^-1
        
        // Magnetic and cosmological parameters
        this.B = 1e-5;                        // Static magnetic field (T)
        this.B_crit = 1e11;                   // Critical B field (T)
        this.Lambda = 1.1e-52;                // Cosmological constant (m^-2)
        this.c_light = 3e8;                   // Speed of light (m/s)
        
        // EM parameters
        this.q_charge = 1.602e-19;            // Proton charge (C)
        this.gas_v = 1e5;                     // Gas velocity for EM (m/s)
        this.f_TRZ = 0.1;                     // Time-reversal factor
        
        // Star formation parameters
        this.SFR_factor = 20.0 / (2e11);      // Star formation rate factor (dimensionless)
        this.tau_SF = 100e6 * 3.156e7;        // Star formation timescale (s) - 100 Myr
        
        // Merger interaction parameters (UNIQUE to Antennae merger)
        this.I0 = 0.1;                        // Initial interaction factor
        this.tau_merger = 400e6 * 3.156e7;    // Merger timescale (s) - 400 Myr
        
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
        this.ug1_base = (this.G * this.M0) / (this.r * this.r);
    }

    // ========== CORE PHYSICS METHODS (13 methods) ==========

    /**
     * Compute mass M(t) with star formation growth
     * M(t) = M0 × (1 + SFR_factor × exp(-t/tau_SF))
     * @param {number} t - Time (s)
     * @returns {number} Mass at time t (kg)
     */
    M_t(t) {
        const M_dot = this.SFR_factor * Math.exp(-t / this.tau_SF);
        return this.M0 * (1 + M_dot);
    }

    /**
     * Compute merger interaction factor I(t)
     * I(t) = I0 × exp(-t/tau_merger) - UNIQUE to Antennae merger
     * @param {number} t - Time (s)
     * @returns {number} Interaction factor (dimensionless)
     */
    I_t(t) {
        return this.I0 * Math.exp(-t / this.tau_merger);
    }

    /**
     * Compute UQFF Ug components with merger interaction
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
     * @returns {number} Total gravitational acceleration g_Antennae (m/s²)
     */
    compute_g_Antennae(t) {
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

        // Term 9: Stellar feedback (pressure / density for acceleration)
        const wind_pressure = this.rho_wind * this.v_wind * this.v_wind;
        const term_feedback = wind_pressure / this.rho_fluid;

        // Total g_Antennae (all 9 terms summed - note: term_osc counts as 1)
        return term1 + term2 + term3 + term4 + term_q + term_fluid + term_osc + term_DM + term_feedback;
    }

    /**
     * Universal getter for any variable by name
     * @param {string} varName - Variable name
     * @returns {number} Variable value
     */
    getVariable(varName) {
        const validVars = [
            'G', 'M0', 'r', 'Hz', 'B', 'B_crit', 'Lambda', 'c_light', 'q_charge', 'gas_v',
            'f_TRZ', 'SFR_factor', 'tau_SF', 'I0', 'tau_merger', 'rho_wind', 'v_wind',
            'rho_fluid', 'rho_vac_UA', 'rho_vac_SCm', 'scale_EM', 'proton_mass', 'z_gal',
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
            'G', 'M0', 'r', 'Hz', 'B', 'B_crit', 'Lambda', 'c_light', 'q_charge', 'gas_v',
            'f_TRZ', 'SFR_factor', 'tau_SF', 'I0', 'tau_merger', 'rho_wind', 'v_wind',
            'rho_fluid', 'rho_vac_UA', 'rho_vac_SCm', 'scale_EM', 'proton_mass', 'z_gal',
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
     * Example computation at t=300 Myr (for testing)
     * @returns {number} g_Antennae at 300 Myr (m/s²)
     */
    exampleAt300Myr() {
        const t_example = 300e6 * 3.156e7;
        return this.compute_g_Antennae(t_example);
    }

    /**
     * Print parameters to console
     */
    printParameters() {
        console.log("Antennae Galaxies Parameters:");
        console.log(`G: ${this.G.toExponential(3)}, M0: ${this.M0.toExponential(3)}, r: ${this.r.toExponential(3)}`);
        console.log(`Hz: ${this.Hz.toExponential(3)}, B: ${this.B.toExponential(3)}, B_crit: ${this.B_crit.toExponential(3)}`);
        console.log(`f_TRZ: ${this.f_TRZ.toFixed(3)}, SFR_factor: ${this.SFR_factor.toExponential(3)}, tau_SF: ${this.tau_SF.toExponential(3)}`);
        console.log(`I0: ${this.I0.toFixed(3)}, tau_merger: ${this.tau_merger.toExponential(3)}`);
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
            'G', 'M0', 'r', 'Hz', 'B', 'B_crit', 'Lambda', 'c_light', 'q_charge', 'gas_v',
            'f_TRZ', 'SFR_factor', 'tau_SF', 'I0', 'tau_merger', 'rho_wind', 'v_wind',
            'rho_fluid', 'rho_vac_UA', 'rho_vac_SCm', 'scale_EM', 'proton_mass', 'z_gal',
            'hbar', 't_Hubble', 't_Hubble_gyr', 'delta_x', 'delta_p', 'integral_psi',
            'A_osc', 'k_osc', 'omega_osc', 'x_pos', 'M_DM_factor', 'delta_rho_over_rho'
        ];
    }

    /**
     * Get system name
     * @returns {string} System name
     */
    getSystemName() {
        return "AntennaeGalaxies";
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

    // ----- Self-Expansion (4 methods - domain-specific for galaxy merger) -----

    /**
     * Expand entire parameter space by a scale factor
     * @param {number} scale_factor - Scaling factor for all scalable parameters
     */
    expandParameterSpace(scale_factor) {
        const scalable = ['M0', 'r', 'I0', 'tau_merger', 'SFR_factor', 'tau_SF',
                          'rho_wind', 'v_wind', 'B'];
        this.scaleVariableGroup(scalable, scale_factor);
    }

    /**
     * Expand galaxy mass and radius scales
     * @param {number} M_scale - Mass scaling factor
     * @param {number} r_scale - Radius scaling factor
     */
    expandGalaxyScale(M_scale, r_scale) {
        this.setVariable('M0', this.getVariable('M0') * M_scale);
        this.setVariable('r', this.getVariable('r') * r_scale);
    }

    /**
     * Expand merger interaction scales (UNIQUE to Antennae merger)
     * @param {number} I0_scale - Interaction strength scaling factor
     * @param {number} tau_merger_scale - Merger timescale scaling factor
     */
    expandMergerScale(I0_scale, tau_merger_scale) {
        this.setVariable('I0', this.getVariable('I0') * I0_scale);
        this.setVariable('tau_merger', this.getVariable('tau_merger') * tau_merger_scale);
    }

    /**
     * Expand star formation scales
     * @param {number} SFR_factor_scale - SFR factor scaling
     * @param {number} tau_SF_scale - SF timescale scaling
     */
    expandStarFormationScale(SFR_factor_scale, tau_SF_scale) {
        this.setVariable('SFR_factor', this.getVariable('SFR_factor') * SFR_factor_scale);
        this.setVariable('tau_SF', this.getVariable('tau_SF') * tau_SF_scale);
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
            const g_model = this.compute_g_Antennae(t);
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
            const g = this.compute_g_Antennae(t);
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
     * @param {Function} fitness - Fitness function (takes AntennaeGalaxies instance)
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
        AntennaeGalaxies.savedStates[label] = state;
    }

    /**
     * Restore state from a label
     * @param {string} label - State label
     * @returns {boolean} Success status
     */
    restoreState(label) {
        const state = AntennaeGalaxies.savedStates[label];
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
        return Object.keys(AntennaeGalaxies.savedStates);
    }

    /**
     * Export current state as string
     * @returns {string} State export string
     */
    exportState() {
        let output = "AntennaeGalaxies State Export\n";
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
        const baseline = this.compute_g_Antennae(t);

        this.listVariables().forEach(varName => {
            const original = this.getVariable(varName);
            this.setVariable(varName, original * (1.0 + perturbation));
            const perturbed = this.compute_g_Antennae(t);
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

        let report = "=== Antennae Galaxies System Report ===\n";
        report += `Time: ${(t / Myr_to_s).toFixed(6)} Myr\n`;
        report += `M(t): ${(this.M_t(t) / M_sun).toExponential(6)} M_sun\n`;
        report += `I(t) interaction factor: ${this.I_t(t).toFixed(6)}\n`;
        report += `g_Antennae: ${this.compute_g_Antennae(t).toExponential(6)} m/s^2\n`;
        report += `Core parameters: M0=${(this.M0/M_sun).toExponential(3)} M_sun, r=${(this.r/ly_to_m).toFixed(0)} ly\n`;
        report += `Merger: I0=${this.I0.toFixed(3)}, tau_merger=${(this.tau_merger/Myr_to_s).toFixed(0)} Myr\n`;
        report += `Star formation: SFR_factor=${this.SFR_factor.toExponential(3)}, tau_SF=${(this.tau_SF/Myr_to_s).toFixed(0)} Myr\n`;
        return report;
    }

    /**
     * Validate parameter consistency
     * @returns {boolean} Validation status
     */
    validateConsistency() {
        let valid = true;
        if (this.M0 <= 0 || this.r <= 0 || this.tau_SF <= 0 || this.tau_merger <= 0) valid = false;
        if (this.I0 < 0 || this.I0 > 1) valid = false;
        if (this.SFR_factor < 0) valid = false;
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

        if (this.M0 <= 0) { this.M0 = 2e11 * M_sun; corrected = true; }
        if (this.r <= 0) { this.r = 30000 * ly_to_m; corrected = true; }
        if (this.tau_SF <= 0) { this.tau_SF = 100 * Myr_to_s; corrected = true; }
        if (this.tau_merger <= 0) { this.tau_merger = 400 * Myr_to_s; corrected = true; }
        if (this.I0 < 0) { this.I0 = 0.0; corrected = true; }
        if (this.I0 > 1) { this.I0 = 1.0; corrected = true; }
        if (this.SFR_factor < 0) { this.SFR_factor = 20.0 / (2e11); corrected = true; }
        if (corrected) this.updateCache();
        return corrected;
    }
}

// ========== ENHANCED EXAMPLE FUNCTION ==========

function enhanced_antennae_example() {
    console.log("=========================================================");
    console.log("ENHANCED ANTENNAE GALAXIES DEMONSTRATION");
    console.log("NGC 4038 & NGC 4039 Galaxy Merger");
    console.log("=========================================================\n");

    const antennae = new AntennaeGalaxies();
    const M_sun = 1.989e30;
    const Myr_to_s = 1e6 * 3.156e7;

    // Step 1: Initial state and validation
    console.log("Step 1: Initial State and Validation");
    console.log(`System: ${antennae.getSystemName()}`);
    console.log(`Validation: ${antennae.validateConsistency() ? 'PASS' : 'FAIL'}`);
    console.log(`Auto-corrected: ${antennae.autoCorrectAnomalies() ? 'Yes' : 'No'}\n`);

    // Step 2: Time evolution showing M(t) and I(t)
    console.log("Step 2: Time Evolution (Star Formation M(t) and Merger Interaction I(t))");
    const t_Myr_array = [0.0, 100.0, 200.0, 300.0, 500.0];
    t_Myr_array.forEach(t_Myr => {
        const t = t_Myr * Myr_to_s;
        const Mt = antennae.M_t(t);
        const It = antennae.I_t(t);
        const g = antennae.compute_g_Antennae(t);
        console.log(`  t = ${t_Myr.toFixed(1)} Myr: M(t) = ${(Mt/M_sun).toExponential(6)} M_sun, I(t) = ${It.toFixed(6)}, g = ${g.toExponential(6)} m/s^2`);
    });
    console.log();

    // Step 3: Variable listing
    console.log("Step 3: Variable Listing");
    const vars = antennae.listVariables();
    console.log(`Total variables: ${vars.length}`);
    console.log(`Sample: ${vars[0]}, ${vars[1]}, ${vars[13]} (I0), ${vars[14]} (tau_merger)\n`);

    // Step 4: Galaxy mass scaling
    console.log("Step 4: Galaxy Mass Scaling (M sweeps)");
    antennae.saveState('original');
    const M_factors = [0.5, 1.0, 2.0];
    M_factors.forEach(factor => {
        antennae.restoreState('original');
        antennae.expandGalaxyScale(factor, 1.0);
        const t = 300 * Myr_to_s;
        const g = antennae.compute_g_Antennae(t);
        const M = antennae.getVariable('M0');
        console.log(`  M × ${factor}: M0 = ${(M/M_sun).toExponential(3)} M_sun, g(300 Myr) = ${g.toExponential(6)} m/s^2`);
    });
    antennae.restoreState('original');
    console.log();

    // Step 5: Merger interaction scaling (UNIQUE to Antennae merger)
    console.log("Step 5: Merger Interaction Scaling (I0 sweeps) - GALAXY MERGER FEATURE");
    const I0_factors = [0.5, 1.0, 2.0];
    I0_factors.forEach(factor => {
        antennae.restoreState('original');
        antennae.expandMergerScale(factor, 1.0);
        const t = 300 * Myr_to_s;
        const It = antennae.I_t(t);
        const g = antennae.compute_g_Antennae(t);
        console.log(`  I0 × ${factor}: I(300 Myr) = ${It.toFixed(6)}, g = ${g.toExponential(6)} m/s^2`);
    });
    antennae.restoreState('original');
    console.log();

    // Step 6: Merger timescale sweeps
    console.log("Step 6: Merger Timescale Sweeps (tau_merger)");
    const tau_merger_factors = [0.5, 1.0, 2.0];
    tau_merger_factors.forEach(factor => {
        antennae.restoreState('original');
        antennae.expandMergerScale(1.0, factor);
        const t = 300 * Myr_to_s;
        const It = antennae.I_t(t);
        console.log(`  tau_merger × ${factor}: I(300 Myr) = ${It.toFixed(6)}`);
    });
    antennae.restoreState('original');
    console.log();

    // Step 7: Star formation rate scaling
    console.log("Step 7: Star Formation Rate Scaling");
    const SFR_factors = [0.5, 1.0, 2.0];
    SFR_factors.forEach(factor => {
        antennae.restoreState('original');
        antennae.expandStarFormationScale(factor, 1.0);
        const t = 300 * Myr_to_s;
        const Mt = antennae.M_t(t);
        console.log(`  SFR_factor × ${factor}: M(300 Myr) = ${(Mt/M_sun).toExponential(6)} M_sun`);
    });
    antennae.restoreState('original');
    console.log();

    // Step 8: Generate report
    console.log("Step 8: Full System Report at 300 Myr");
    const t_report = 300 * Myr_to_s;
    const report = antennae.generateReport(t_report);
    console.log(report);

    console.log("=========================================================");
    console.log("ENHANCED DEMONSTRATION COMPLETE");
    console.log("=========================================================\n");
}

// Run inline test if executed directly
if (typeof require !== 'undefined' && require.main === module) {
    console.log("Running inline test...\n");
    enhanced_antennae_example();
}

// Export for use as module
export default AntennaeGalaxies;
