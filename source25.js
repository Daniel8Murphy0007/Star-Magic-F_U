/**
 * ================================================================================================
 * Module: NGC1275.js
 *
 * Description: JavaScript ES6 Module for NGC 1275 (Magnetic Monster Perseus A) Class
 *              This is the sixteenth module in a series of 500+ code files for the Universal
 *              Quantum Field Framework (UQFF) simulations, focusing on active galactic nucleus
 *              evolution and gravity equations derived from Hubble datasets, high-energy lab
 *              simulations, and UQFF refinements (dated May 09, 2025, updated October 08, 2025).
 *
 * Purpose: Encapsulates the Master Universal Gravity Equation (MUGE) for NGC 1275 evolution.
 *          Includes ALL terms: base gravity (static M), cosmic expansion (H(z)), magnetic field
 *          B(t), filament support F(t), black hole influence, UQFF Ug components with f_TRZ,
 *          Lambda, quantum uncertainty, scaled EM with [UA], fluid dynamics, oscillatory waves,
 *          DM/density perturbations, and cooling flow term.
 *
 * Key Features:
 *   - Default values from UQFF document: M = 1e11 M☉, r = 200,000 ly, M_BH = 8e8 M☉,
 *     z = 0.0176, Hz ≈ 2.20e-18 s^-1, B0 = 5e-9 T, tau_B = 100 Myr, F0 = 0.1,
 *     tau_fil = 100 Myr, rho_cool = 1e-20 kg/m³, v_cool = 3e3 m/s.
 *   - B(t) = B0 × exp(-t/tau_B) - magnetic field exponential decay
 *   - F(t) = F0 × exp(-t/tau_fil) - filament support exponential decay
 *   - Supermassive black hole contribution: g_BH = G×M_BH/r_BH²
 *   - Units handled: M☉ to kg, ly to m; cooling term as (rho × v_cool²) / rho_fluid
 *   - Setter methods: setVariable, addToVariable, subtractFromVariable
 *   - Computes g_NGC1275(r, t) with every term explicitly included (11 terms total)
 *   - Enhanced dynamic capabilities: 38 methods total (13 core + 25 enhanced)
 *
 * Conversion: Ported from C++ source25.cpp (NGC1275.h) to JavaScript ES6
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript
 * Date: November 03, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

class NGC1275 {
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
        this.M = 1e11 * M_sun;                // Total galaxy mass (kg)
        const ly_to_m = 9.461e15;             // Light year to meters
        this.r = 200000.0 * ly_to_m;          // Radius (m) - 200,000 ly
        this.z_gal = 0.0176;                  // Galaxy redshift
        
        // Hubble parameter at z
        const Hz_kms = 70 * Math.sqrt(0.3 * Math.pow(1 + this.z_gal, 3) + 0.7);  // km/s/Mpc
        this.Hz = (Hz_kms * 1000 / 3.086e19); // Convert to s^-1
        
        // Magnetic field parameters (UNIQUE to magnetic monster)
        this.B0 = 5e-9;                       // Initial magnetic field (T)
        this.tau_B = 100e6 * 3.156e7;         // B decay timescale (s) - 100 Myr
        this.B_crit = 1e11;                   // Critical B field (T)
        
        // Cosmological parameters
        this.Lambda = 1.1e-52;                // Cosmological constant (m^-2)
        this.c_light = 3e8;                   // Speed of light (m/s)
        
        // EM parameters
        this.q_charge = 1.602e-19;            // Proton charge (C)
        this.gas_v = 1e5;                     // Gas velocity for EM (m/s)
        this.f_TRZ = 0.1;                     // Time-reversal factor
        
        // Black hole parameters (UNIQUE to AGN)
        this.M_BH = 8e8 * M_sun;              // Black hole mass (kg) - 800 million M☉
        this.r_BH = 1e18;                     // Black hole radius (m) - influence radius
        
        // Filament parameters (UNIQUE to NGC 1275)
        this.F0 = 0.1;                        // Initial filament factor
        this.tau_fil = 100e6 * 3.156e7;       // Filament timescale (s) - 100 Myr
        
        // Cooling flow parameters
        this.rho_cool = 1e-20;                // Cooling flow density (kg/m³)
        this.v_cool = 3e3;                    // Cooling flow velocity (m/s)
        this.rho_fluid = 1e-20;               // Fluid density (kg/m³)
        
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
        this.g_BH = (this.G * this.M_BH) / (this.r_BH * this.r_BH);
    }

    // ========== CORE PHYSICS METHODS (13 methods) ==========

    /**
     * Compute magnetic field B(t)
     * B(t) = B0 × exp(-t/tau_B) - exponential decay
     * @param {number} t - Time (s)
     * @returns {number} Magnetic field at time t (T)
     */
    B_t(t) {
        return this.B0 * Math.exp(-t / this.tau_B);
    }

    /**
     * Compute filament support factor F(t)
     * F(t) = F0 × exp(-t/tau_fil) - exponential decay
     * @param {number} t - Time (s)
     * @returns {number} Filament factor (dimensionless)
     */
    F_t(t) {
        return this.F0 * Math.exp(-t / this.tau_fil);
    }

    /**
     * Compute UQFF Ug components with B and F corrections
     * @param {number} Bt - Magnetic field at time t (T)
     * @param {number} Ft - Filament factor at time t
     * @returns {number} Total Ug (m/s²)
     */
    compute_Ug(Bt, Ft) {
        const Ug1 = this.ug1_base;
        const Ug2 = 0.0;
        const Ug3 = 0.0;
        const corr_B = 1 - Bt / this.B_crit;
        const Ug4 = Ug1 * corr_B;
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1 + this.f_TRZ) * (1 + Ft);
    }

    /**
     * Compute volume for fluid dynamics
     * @returns {number} Volume (m³)
     */
    compute_V() {
        return (4.0 / 3.0) * Math.PI * this.r * this.r * this.r;
    }

    /**
     * Main MUGE computation - includes ALL 11 terms
     * @param {number} t - Time (s)
     * @returns {number} Total gravitational acceleration g_NGC1275 (m/s²)
     */
    compute_g_NGC1275(t) {
        if (t < 0) {
            console.error("Error: Time t must be non-negative.");
            return 0.0;
        }

        const Bt = this.B_t(t);
        const Ft = this.F_t(t);

        // Term 1: Base + Hz + B + F corrections
        const corr_H = 1 + this.Hz * t;
        const corr_B = 1 - Bt / this.B_crit;
        const corr_F = 1 + Ft;
        const term1 = this.ug1_base * corr_H * corr_B * corr_F;

        // Term BH: Black hole contribution (UNIQUE to AGN)
        const term_BH = this.g_BH;

        // Term 2: UQFF Ug with f_TRZ, B, F
        const term2 = this.compute_Ug(Bt, Ft);

        // Term 3: Lambda (cosmological constant)
        const term3 = (this.Lambda * this.c_light * this.c_light) / 3.0;

        // Term 4: Scaled EM with UA (using B(t))
        const cross_vB = this.gas_v * Bt;  // Magnitude, assuming perpendicular
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

        // Term 9: Cooling flow (pressure / density for acceleration)
        const cool_pressure = this.rho_cool * this.v_cool * this.v_cool;
        const term_cool = cool_pressure / this.rho_fluid;

        // Total g_NGC1275 (all 10 terms summed - note: term_osc counts as 1, term_BH is separate)
        return term1 + term_BH + term2 + term3 + term4 + term_q + term_fluid + term_osc + term_DM + term_cool;
    }

    /**
     * Universal getter for any variable by name
     * @param {string} varName - Variable name
     * @returns {number} Variable value
     */
    getVariable(varName) {
        const validVars = [
            'G', 'M', 'r', 'Hz', 'B0', 'tau_B', 'B_crit', 'Lambda', 'c_light', 'q_charge',
            'gas_v', 'f_TRZ', 'M_BH', 'r_BH', 'F0', 'tau_fil', 'rho_cool', 'v_cool',
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
            'G', 'M', 'r', 'Hz', 'B0', 'tau_B', 'B_crit', 'Lambda', 'c_light', 'q_charge',
            'gas_v', 'f_TRZ', 'M_BH', 'r_BH', 'F0', 'tau_fil', 'rho_cool', 'v_cool',
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
     * Example computation at t=50 Myr (for testing)
     * @returns {number} g_NGC1275 at 50 Myr (m/s²)
     */
    exampleAt50Myr() {
        const t_example = 50e6 * 3.156e7;
        return this.compute_g_NGC1275(t_example);
    }

    /**
     * Print parameters to console
     */
    printParameters() {
        console.log("NGC 1275 Parameters:");
        console.log(`G: ${this.G.toExponential(3)}, M: ${this.M.toExponential(3)}, r: ${this.r.toExponential(3)}`);
        console.log(`Hz: ${this.Hz.toExponential(3)}, B0: ${this.B0.toExponential(3)}, tau_B: ${this.tau_B.toExponential(3)}`);
        console.log(`f_TRZ: ${this.f_TRZ.toFixed(3)}, M_BH: ${this.M_BH.toExponential(3)}, r_BH: ${this.r_BH.toExponential(3)}`);
        console.log(`F0: ${this.F0.toFixed(3)}, tau_fil: ${this.tau_fil.toExponential(3)}`);
        console.log(`rho_fluid: ${this.rho_fluid.toExponential(3)}, rho_cool: ${this.rho_cool.toExponential(3)}, v_cool: ${this.v_cool.toExponential(3)}`);
        console.log(`gas_v: ${this.gas_v.toExponential(3)}, M_DM_factor: ${this.M_DM_factor.toFixed(3)}`);
        console.log(`A_osc: ${this.A_osc.toExponential(3)}, delta_rho_over_rho: ${this.delta_rho_over_rho.toExponential(3)}`);
        console.log(`ug1_base: ${this.ug1_base.toExponential(3)}, g_BH: ${this.g_BH.toExponential(3)}`);
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
            'G', 'M', 'r', 'Hz', 'B0', 'tau_B', 'B_crit', 'Lambda', 'c_light', 'q_charge',
            'gas_v', 'f_TRZ', 'M_BH', 'r_BH', 'F0', 'tau_fil', 'rho_cool', 'v_cool',
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
        return "NGC1275";
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

    // ----- Self-Expansion (4 methods - domain-specific for AGN) -----

    /**
     * Expand entire parameter space by a scale factor
     * @param {number} scale_factor - Scaling factor for all scalable parameters
     */
    expandParameterSpace(scale_factor) {
        const scalable = ['M', 'r', 'M_BH', 'r_BH', 'B0', 'F0', 'tau_B', 'tau_fil',
                          'rho_cool', 'v_cool'];
        this.scaleVariableGroup(scalable, scale_factor);
    }

    /**
     * Expand galaxy mass and radius scales
     * @param {number} M_scale - Mass scaling factor
     * @param {number} r_scale - Radius scaling factor
     */
    expandGalaxyScale(M_scale, r_scale) {
        this.setVariable('M', this.getVariable('M') * M_scale);
        this.setVariable('r', this.getVariable('r') * r_scale);
    }

    /**
     * Expand black hole scales (UNIQUE to AGN)
     * @param {number} M_BH_scale - Black hole mass scaling
     * @param {number} r_BH_scale - Black hole radius scaling
     */
    expandBlackHoleScale(M_BH_scale, r_BH_scale) {
        this.setVariable('M_BH', this.getVariable('M_BH') * M_BH_scale);
        this.setVariable('r_BH', this.getVariable('r_BH') * r_BH_scale);
    }

    /**
     * Expand magnetic and filament scales (UNIQUE to NGC 1275)
     * @param {number} B0_scale - Magnetic field scaling
     * @param {number} F0_scale - Filament factor scaling
     */
    expandMagneticFilamentScale(B0_scale, F0_scale) {
        this.setVariable('B0', this.getVariable('B0') * B0_scale);
        this.setVariable('F0', this.getVariable('F0') * F0_scale);
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
            const g_model = this.compute_g_NGC1275(t);
            total_error += Math.abs(g_model - g_obs);
        });

        const avg_error = total_error / obs_data.length;
        if (avg_error > 1e-9) {
            const correction = 0.95;
            this.setVariable('B0', this.getVariable('B0') * correction);
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
            const g = this.compute_g_NGC1275(t);
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
     * @param {Function} fitness - Fitness function (takes NGC1275 instance)
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
        NGC1275.savedStates[label] = state;
    }

    /**
     * Restore state from a label
     * @param {string} label - State label
     * @returns {boolean} Success status
     */
    restoreState(label) {
        const state = NGC1275.savedStates[label];
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
        return Object.keys(NGC1275.savedStates);
    }

    /**
     * Export current state as string
     * @returns {string} State export string
     */
    exportState() {
        let output = "NGC1275 State Export\n";
        output += "====================\n";

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
        const baseline = this.compute_g_NGC1275(t);

        this.listVariables().forEach(varName => {
            const original = this.getVariable(varName);
            this.setVariable(varName, original * (1.0 + perturbation));
            const perturbed = this.compute_g_NGC1275(t);
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

        let report = "=== NGC 1275 System Report ===\n";
        report += `Time: ${(t / Myr_to_s).toFixed(6)} Myr\n`;
        report += `B(t) magnetic field: ${this.B_t(t).toExponential(6)} T\n`;
        report += `F(t) filament factor: ${this.F_t(t).toFixed(6)}\n`;
        report += `g_NGC1275: ${this.compute_g_NGC1275(t).toExponential(6)} m/s^2\n`;
        report += `Core parameters: M=${(this.M/M_sun).toExponential(3)} M_sun, r=${(this.r/ly_to_m).toFixed(0)} ly\n`;
        report += `Black hole: M_BH=${(this.M_BH/M_sun).toExponential(3)} M_sun, r_BH=${this.r_BH.toExponential(3)} m\n`;
        report += `Magnetic: B0=${this.B0.toExponential(3)} T, tau_B=${(this.tau_B/Myr_to_s).toFixed(0)} Myr\n`;
        report += `Filament: F0=${this.F0.toFixed(3)}, tau_fil=${(this.tau_fil/Myr_to_s).toFixed(0)} Myr\n`;
        report += `Cooling: rho_cool=${this.rho_cool.toExponential(3)} kg/m^3, v_cool=${this.v_cool.toExponential(3)} m/s\n`;
        return report;
    }

    /**
     * Validate parameter consistency
     * @returns {boolean} Validation status
     */
    validateConsistency() {
        let valid = true;
        if (this.M <= 0 || this.r <= 0 || this.M_BH <= 0 || this.r_BH <= 0) valid = false;
        if (this.tau_B <= 0 || this.tau_fil <= 0) valid = false;
        if (this.B0 < 0 || this.F0 < 0 || this.F0 > 1) valid = false;
        if (this.rho_cool < 0 || this.v_cool < 0) valid = false;
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

        if (this.M <= 0) { this.M = 1e11 * M_sun; corrected = true; }
        if (this.r <= 0) { this.r = 200000 * ly_to_m; corrected = true; }
        if (this.M_BH <= 0) { this.M_BH = 8e8 * M_sun; corrected = true; }
        if (this.r_BH <= 0) { this.r_BH = 1e18; corrected = true; }
        if (this.tau_B <= 0) { this.tau_B = 100 * Myr_to_s; corrected = true; }
        if (this.tau_fil <= 0) { this.tau_fil = 100 * Myr_to_s; corrected = true; }
        if (this.B0 < 0) { this.B0 = 5e-9; corrected = true; }
        if (this.F0 < 0) { this.F0 = 0.0; corrected = true; }
        if (this.F0 > 1) { this.F0 = 1.0; corrected = true; }
        if (this.rho_cool < 0) { this.rho_cool = 1e-20; corrected = true; }
        if (this.v_cool < 0) { this.v_cool = 3e3; corrected = true; }
        if (corrected) this.updateCache();
        return corrected;
    }
}

// ========== ENHANCED EXAMPLE FUNCTION ==========

function enhanced_ngc1275_example() {
    console.log("=========================================================");
    console.log("ENHANCED NGC 1275 DEMONSTRATION");
    console.log("Perseus A - Magnetic Monster AGN");
    console.log("=========================================================\n");

    const ngc = new NGC1275();
    const M_sun = 1.989e30;
    const Myr_to_s = 1e6 * 3.156e7;

    // Step 1: Initial state and validation
    console.log("Step 1: Initial State and Validation");
    console.log(`System: ${ngc.getSystemName()}`);
    console.log(`Validation: ${ngc.validateConsistency() ? 'PASS' : 'FAIL'}`);
    console.log(`Auto-corrected: ${ngc.autoCorrectAnomalies() ? 'Yes' : 'No'}\n`);

    // Step 2: Time evolution showing B(t) and F(t)
    console.log("Step 2: Time Evolution (Magnetic Field B(t) and Filament Support F(t))");
    const t_Myr_array = [0.0, 25.0, 50.0, 100.0, 200.0];
    t_Myr_array.forEach(t_Myr => {
        const t = t_Myr * Myr_to_s;
        const Bt = ngc.B_t(t);
        const Ft = ngc.F_t(t);
        const g = ngc.compute_g_NGC1275(t);
        console.log(`  t = ${t_Myr.toFixed(1)} Myr: B(t) = ${Bt.toExponential(6)} T, F(t) = ${Ft.toFixed(6)}, g = ${g.toExponential(6)} m/s^2`);
    });
    console.log();

    // Step 3: Variable listing
    console.log("Step 3: Variable Listing");
    const vars = ngc.listVariables();
    console.log(`Total variables: ${vars.length}`);
    console.log(`Sample: ${vars[0]}, ${vars[1]}, ${vars[12]} (M_BH), ${vars[14]} (F0)\n`);

    // Step 4: Black hole scaling (UNIQUE to AGN with SMBH)
    console.log("Step 4: Black Hole Scaling (M_BH sweeps) - AGN FEATURE");
    ngc.saveState('original');
    const M_BH_factors = [0.5, 1.0, 2.0];
    M_BH_factors.forEach(factor => {
        ngc.restoreState('original');
        ngc.expandBlackHoleScale(factor, 1.0);
        const t = 50 * Myr_to_s;
        const g = ngc.compute_g_NGC1275(t);
        const M_BH = ngc.getVariable('M_BH');
        console.log(`  M_BH × ${factor}: M_BH = ${(M_BH/M_sun).toExponential(3)} M_sun, g(50 Myr) = ${g.toExponential(6)} m/s^2`);
    });
    ngc.restoreState('original');
    console.log();

    // Step 5: Magnetic field scaling (UNIQUE to magnetic monster)
    console.log("Step 5: Magnetic Field Scaling (B0 sweeps) - MAGNETIC MONSTER FEATURE");
    const B0_factors = [0.5, 1.0, 2.0];
    B0_factors.forEach(factor => {
        ngc.restoreState('original');
        ngc.expandMagneticFilamentScale(factor, 1.0);
        const t = 50 * Myr_to_s;
        const Bt = ngc.B_t(t);
        const g = ngc.compute_g_NGC1275(t);
        console.log(`  B0 × ${factor}: B(50 Myr) = ${Bt.toExponential(6)} T, g = ${g.toExponential(6)} m/s^2`);
    });
    ngc.restoreState('original');
    console.log();

    // Step 6: Generate report
    console.log("Step 6: Full System Report at 50 Myr");
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
    enhanced_ngc1275_example();
}

// Export for use as module
export default NGC1275;
