/**
 * ================================================================================================
 * Module: NebulaCrab.js
 *
 * Description: JavaScript ES6 Module for Crab Nebula (M1) UQFF Class
 *              This is the twenty-fourth module in a series of 500+ code files for the Universal
 *              Quantum Field Framework (UQFF) simulations, focusing on Crab Nebula evolution
 *              and gravity equations derived from Hubble datasets, high-energy lab simulations,
 *              and UQFF refinements (dated May 09, 2025, updated October 09, 2025).
 *
 * Purpose: Encapsulates the Master Universal Gravity Equation (MUGE) for Crab Nebula.
 *          Includes ALL terms: base gravity with r(t), UQFF Ug components, Lambda,
 *          quantum uncertainty, EM/Lorentz, fluid dynamics, resonant waves, DM with density
 *          perturbations, superconductivity correction, pulsar wind a_wind, magnetic M_mag.
 *
 * Key Features:
 *   - Default values from UQFF document: M = 4.6 M☉, r0 = 5.2×10¹⁶ m (~5.5 ly),
 *     v_exp = 1.5×10⁶ m/s, P_pulsar = 5×10³¹ W, t = 971 years, z = 0.0015.
 *   - Time-dependent radius: r(t) = r0 + v_exp × t (expanding remnant)
 *   - Pulsar Wind: a_wind = [P_pulsar / (4π r²) × (1 + v_shock/c)] / ρ (dominant outward force)
 *   - Magnetic Force: M_mag = (q × v_shock × B) / m_e (Lorentz acceleration on electrons)
 *   - Units handled: All SI (kg, m, s, T, W)
 *   - Setter methods: setVariable, addToVariable, subtractFromVariable
 *   - Computes g_Crab(r, t) with every term explicitly included (10 major terms)
 *   - Enhanced dynamic capabilities: 38 methods total (13 core + 25 enhanced)
 *
 * Conversion: Ported from C++ source32.cpp (CrabUQFFModule.h) to JavaScript ES6
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript
 * Date: November 03, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

class NebulaCrab {
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
        this.G = 6.6743e-11;                      // Gravitational constant (m³/kg/s²)
        this.c = 3e8;                             // Speed of light (m/s)
        this.hbar = 1.0546e-34;                   // Reduced Planck's constant (J·s)
        this.Lambda = 1.1e-52;                    // Cosmological constant (m⁻²)
        this.q = 1.602e-19;                       // Electron charge (C)
        this.pi = Math.PI;                        // π
        this.t_Hubble = 13.8e9 * 3.156e7;         // Hubble time (s) - 13.8 Gyr
        this.year_to_s = 3.156e7;                 // Seconds per year
        
        // Crab Nebula parameters
        this.M_sun = 1.989e30;                    // Solar mass (kg)
        this.M = 4.6 * this.M_sun;                // Total mass (kg)
        this.M_visible = this.M;                  // Visible mass (ejecta + pulsar)
        this.M_DM = 0.0;                          // No significant DM
        this.r0 = 5.2e16;                         // Initial radius (m) ~5.5 ly
        this.v_exp = 1.5e6;                       // Expansion velocity (m/s)
        
        // Hubble/cosmology
        this.H0 = 70.0;                           // Hubble constant (km/s/Mpc)
        this.Mpc_to_m = 3.086e22;                 // Megaparsec to meters
        this.z = 0.0015;                          // Redshift
        this.Omega_m = 0.3;                       // Matter density parameter
        this.Omega_Lambda = 0.7;                  // Dark energy density parameter
        this.t = 971 * this.year_to_s;            // Default time = 971 years (s) since 1054 AD
        
        // Nebula dynamics
        this.rho_fluid = 1e-21;                   // Filament density (kg/m³)
        this.V = 1e3;                             // Volume scale (m³)
        this.v_shock = 1.5e6;                     // Shock velocity (m/s)
        this.P_pulsar = 5e31;                     // Pulsar luminosity (W)
        this.delta_rho = 0.1 * this.rho_fluid;    // Density perturbation (kg/m³)
        this.rho = this.rho_fluid;                // Mean density (kg/m³)
        
        // EM/magnetic/superconductivity
        this.B = 1e-8;                            // Nebula average magnetic field (T)
        this.B_crit = 1e11;                       // Critical field (10¹⁵ G ≈ 1e11 T)
        this.m_e = 9.11e-31;                      // Electron mass (kg)
        
        // Quantum terms
        this.Delta_x = 1e-10;                     // Position uncertainty (m)
        this.Delta_p = this.hbar / this.Delta_x;  // Momentum uncertainty (kg·m/s)
        this.integral_psi = 1.0;                  // Normalized wavefunction integral
        
        // Resonant/oscillatory terms
        this.A = 1e-10;                           // Oscillatory amplitude (m/s²)
        this.k = 1e20;                            // Wave number (m⁻¹)
        this.omega = 1e15;                        // Angular frequency (rad/s)
        this.x = 0.0;                             // Position (m) - central
        
        // Ug subterms (computed dynamically)
        this.Ug1 = 0.0;                           // G M / r²
        this.Ug2 = 0.0;                           // d²Φ/dt² ≈ 0 (negligible)
        this.Ug3 = 0.0;                           // Moon term ≈ 0 (no moon)
        this.Ug4 = 0.0;                           // Ug1 × f_sc
        
        // Scale factors
        this.scale_macro = 1e-12;                 // Macro effects scaling
        this.f_TRZ = 0.1;                         // Time-reversal factor
        this.f_sc = 1.0;                          // Superconductive factor
        
        // UA/SCm vacuum densities for EM correction
        this.rho_vac_UA = 7.09e-36;               // UA vacuum density (J/m³)
        this.rho_vac_SCm = 7.09e-37;              // SCm vacuum density (J/m³)
        
        this.updateCache();
    }

    /**
     * Update cached values (call after parameter changes)
     */
    updateCache() {
        const r_current = this.r0 + this.v_exp * this.t;
        this.Ug1 = (this.G * this.M) / (r_current * r_current);
        this.Ug4 = this.Ug1 * this.f_sc;
        this.Delta_p = this.hbar / this.Delta_x;
    }

    // ========== CORE PHYSICS METHODS (15 methods) ==========

    /**
     * Compute current radius r(t) = r0 + v_exp × t
     * @param {number} t - Time (s)
     * @returns {number} Current radius (m)
     */
    computeRadius(t) {
        return this.r0 + this.v_exp * t;
    }

    /**
     * Compute Hubble parameter H(z) in s⁻¹
     * @returns {number} Hz (s⁻¹)
     */
    computeHz() {
        const Hz_kms = this.H0 * Math.sqrt(this.Omega_m * Math.pow(1.0 + this.z, 3) + this.Omega_Lambda);
        return (Hz_kms * 1e3) / this.Mpc_to_m;
    }

    /**
     * Compute UQFF Ug sum: Ug1 + Ug2 + Ug3 + Ug4
     * Uses current r(t) for Ug1 = G M / r²
     * @param {number} r - Current radius (m)
     * @returns {number} Total Ug (m/s²)
     */
    computeUgSum(r) {
        this.Ug1 = (this.G * this.M) / (r * r);
        this.Ug4 = this.Ug1 * this.f_sc;
        return this.Ug1 + this.Ug2 + this.Ug3 + this.Ug4;
    }

    /**
     * Quantum term: (hbar / sqrt(Δx Δp)) × ∫ψ*Hψ dV × (2π / t_Hubble)
     * @param {number} t_Hubble_val - Hubble time (s)
     * @returns {number} Quantum contribution (m/s²)
     */
    computeQuantumTerm(t_Hubble_val) {
        const unc = Math.sqrt(this.Delta_x * this.Delta_p);
        const integral_val = this.integral_psi;  // Normalized to 1
        return (this.hbar / unc) * integral_val * (2 * this.pi / t_Hubble_val);
    }

    /**
     * Fluid term: ρ_fluid × V × g (g approximated with base gravity)
     * @param {number} g_base - Base gravitational acceleration (m/s²)
     * @returns {number} Fluid contribution (m/s²)
     */
    computeFluidTerm(g_base) {
        return this.rho_fluid * this.V * g_base;
    }

    /**
     * Resonant terms: 2A cos(kx) cos(ωt) + (2π/13.8) A Re[exp(i(kx - ωt))]
     * @param {number} t - Time (s)
     * @returns {number} Resonant contribution (m/s²)
     */
    computeResonantTerm(t) {
        const cos_term = 2 * this.A * Math.cos(this.k * this.x) * Math.cos(this.omega * t);
        // Complex exponential real part: A cos(kx - ωt)
        const phase = this.k * this.x - this.omega * t;
        const real_exp = this.A * Math.cos(phase);
        const exp_factor = (2 * this.pi) / 13.8;  // Unitless factor
        return cos_term + exp_factor * real_exp;
    }

    /**
     * DM term: (M_visible + M_DM) × (δρ/ρ + 3GM/r³)
     * @param {number} r - Current radius (m)
     * @returns {number} DM contribution (force-like, needs normalization by mass)
     */
    computeDMTerm(r) {
        const pert = this.delta_rho / this.rho;
        const curv = (3 * this.G * this.M) / (r * r * r);
        return (this.M_visible + this.M_DM) * (pert + curv);
    }

    /**
     * Pulsar wind term: [P_pulsar / (4π r²) × (1 + v_shock/c)] / ρ_fluid × scale_macro
     * @param {number} r - Current radius (m)
     * @returns {number} Wind contribution (m/s²)
     */
    computeWindTerm(r) {
        const pressure = (this.P_pulsar / (4 * this.pi * r * r)) * (1.0 + this.v_shock / this.c);
        return (pressure / this.rho_fluid) * this.scale_macro;
    }

    /**
     * Magnetic term: (q × v_shock × B) / m_e × scale_macro
     * @returns {number} Magnetic contribution (m/s²)
     */
    computeMagTerm() {
        const force = this.q * this.v_shock * this.B;
        return (force / this.m_e) * this.scale_macro;
    }

    /**
     * Main MUGE computation - includes ALL 10 major terms with r(t)
     * @param {number} t - Time (s)
     * @returns {number} Total gravitational acceleration g_Crab (m/s²)
     */
    compute_g_Crab(t) {
        if (t < 0) {
            console.error("Error: Time t must be non-negative.");
            return 0.0;
        }

        this.t = t;  // Update time
        const r = this.computeRadius(t);
        const Hz = this.computeHz();
        const expansion = 1.0 + Hz * t;
        const sc_correction = 1.0 - (this.B / this.B_crit);
        const tr_factor = 1.0 + this.f_TRZ;

        // Term 1: Base gravity with expansion, SC, TR
        const g_base = ((this.G * this.M) / (r * r)) * expansion * sc_correction * tr_factor;

        // Term 2: UQFF Ug sum (Ug1 + Ug2 + Ug3 + Ug4)
        const ug_sum = this.computeUgSum(r);

        // Term 3: Cosmological constant
        const lambda_term = (this.Lambda * this.c * this.c) / 3.0;

        // Term 4: Quantum uncertainty
        const quantum_term = this.computeQuantumTerm(this.t_Hubble);

        // Term 5: EM Lorentz (q v_shock × B, with UA/SCm correction)
        const em_base = (this.q * this.v_shock * this.B) / 1.673e-27;  // / proton mass for acceleration
        const corr_UA = 1.0 + (this.rho_vac_UA / this.rho_vac_SCm);  // = 1 + 10 = 11
        const em_term = em_base * corr_UA * this.scale_macro;

        // Term 6: Fluid dynamics
        const fluid_term = this.computeFluidTerm(g_base);

        // Term 7: Resonant oscillatory
        const resonant_term = this.computeResonantTerm(t);

        // Term 8: DM with density perturbations (normalized by M for acceleration)
        const dm_force_like = this.computeDMTerm(r);
        const dm_term = dm_force_like / this.M;  // Convert force-like to acceleration

        // Term 9: Pulsar wind (DOMINANT in Crab Nebula)
        const wind_term = this.computeWindTerm(r);

        // Term 10: Magnetic acceleration
        const mag_term = this.computeMagTerm();

        // Total g_Crab: sum all 10 terms
        return g_base + ug_sum + lambda_term + quantum_term + em_term + fluid_term + resonant_term + dm_term + wind_term + mag_term;
    }

    /**
     * Universal getter for any variable by name
     * @param {string} varName - Variable name
     * @returns {number} Variable value
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
     * Universal setter for any variable by name
     * @param {string} varName - Variable name
     * @param {number} newValue - New value
     * @returns {boolean} Success status
     */
    setVariable(varName, newValue) {
        if (this.hasOwnProperty(varName)) {
            this[varName] = newValue;
            
            // Handle dependent variables
            if (varName === 'Delta_x') {
                this.Delta_p = this.hbar / newValue;
            } else if (varName === 'M') {
                this.M_visible = newValue;
                this.M_DM = 0.0;
            } else if (varName === 'rho_fluid') {
                this.rho = newValue;
                this.delta_rho = 0.1 * newValue;
            }
            
            this.updateCache();
            return true;
        } else {
            console.warn(`Warning: Variable '${varName}' not found. Creating new variable.`);
            this[varName] = newValue;
            return true;
        }
    }

    /**
     * Add to variable value
     * @param {string} varName - Variable name
     * @param {number} delta - Amount to add
     * @returns {boolean} Success status
     */
    addToVariable(varName, delta) {
        const current = this.getVariable(varName);
        return this.setVariable(varName, current + delta);
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
     * Example computation at t=971 years (for testing)
     * @returns {number} g_Crab at 971 years (m/s²)
     */
    exampleAt971Years() {
        const t_example = 971 * this.year_to_s;
        return this.compute_g_Crab(t_example);
    }

    /**
     * Get equation text (descriptive)
     * @returns {string} Equation description
     */
    getEquationText() {
        return "g_Crab(r, t) = (G × M / r(t)²) × (1 + H(z) × t) × (1 - B / B_crit) × (1 + f_TRZ) + " +
               "(Ug1 + Ug2 + Ug3 + Ug4) + (Λ × c² / 3) + " +
               "(ℏ / √(Δx × Δp)) × ∫(ψ* H ψ dV) × (2π / t_Hubble) + q (v × B) + " +
               "ρ_fluid × V × g + 2A cos(k x) cos(ω t) + (2π / 13.8) A exp(i (k x - ω t)) + " +
               "(M_visible + M_DM) × (δρ/ρ + 3 G M / r³) + a_wind + M_mag\n" +
               "Where r(t) = r0 + v_exp × t; " +
               "a_wind = [P_pulsar / (4π r²) × (1 + v_shock/c)] / ρ × 1e-12; " +
               "M_mag = (q v_shock B) / m_e × 1e-12";
    }

    /**
     * Print parameters to console
     */
    printParameters() {
        console.log("Crab Nebula (M1) Parameters:");
        console.log(`G: ${this.G.toExponential(3)}, M: ${(this.M/this.M_sun).toFixed(1)} M_sun, r0: ${(this.r0/9.461e15).toFixed(2)} ly`);
        console.log(`v_exp: ${(this.v_exp/1e3).toFixed(1)} km/s, P_pulsar: ${this.P_pulsar.toExponential(3)} W`);
        console.log(`z: ${this.z.toFixed(4)}, rho_fluid: ${this.rho_fluid.toExponential(3)} kg/m³`);
        console.log(`v_shock: ${(this.v_shock/1e3).toFixed(1)} km/s, B: ${this.B.toExponential(3)} T`);
        console.log(`Age: ${(this.t/this.year_to_s).toFixed(0)} years, r(t): ${(this.computeRadius(this.t)/9.461e15).toFixed(2)} ly`);
        console.log(`SC correction: ${(1 - this.B/this.B_crit).toFixed(6)}`);
        console.log(`Ug1: ${this.Ug1.toExponential(3)} m/s²`);
    }

    // ========== ENHANCED DYNAMIC CAPABILITIES (25 methods) ==========

    // ----- Variable Management (5 methods) -----

    /**
     * Create or update a variable
     * @param {string} name - Variable name
     * @param {number} value - Variable value
     */
    createVariable(name, value) {
        this[name] = value;
    }

    /**
     * Remove a variable (resets to defaults)
     * @param {string} name - Variable name
     * @returns {boolean} Success status
     */
    removeVariable(name) {
        if (this.hasOwnProperty(name)) {
            delete this[name];
            return true;
        }
        return false;
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
        const vars = [];
        for (const key in this) {
            if (this.hasOwnProperty(key) && typeof this[key] === 'number') {
                vars.push(key);
            }
        }
        return vars;
    }

    /**
     * Get system name
     * @returns {string} System name
     */
    getSystemName() {
        return "CrabNebula_M1";
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

    // ----- Self-Expansion (4 methods - domain-specific for Crab) -----

    /**
     * Expand entire parameter space by a scale factor
     * @param {number} scale_factor - Scaling factor for all scalable parameters
     */
    expandParameterSpace(scale_factor) {
        const scalable = ['M', 'r0', 'v_exp', 'P_pulsar', 'B', 'v_shock', 'rho_fluid'];
        this.scaleVariableGroup(scalable, scale_factor);
    }

    /**
     * Expand nebula mass and initial radius scales
     * @param {number} M_scale - Mass scaling factor
     * @param {number} r0_scale - Initial radius scaling factor
     */
    expandNebulaScale(M_scale, r0_scale) {
        this.setVariable('M', this.M * M_scale);
        this.setVariable('r0', this.r0 * r0_scale);
    }

    /**
     * Expand pulsar scales (UNIQUE to pulsar-driven remnants)
     * @param {number} P_pulsar_scale - Pulsar power scaling
     * @param {number} v_exp_scale - Expansion velocity scaling
     */
    expandPulsarScale(P_pulsar_scale, v_exp_scale) {
        this.setVariable('P_pulsar', this.P_pulsar * P_pulsar_scale);
        this.setVariable('v_exp', this.v_exp * v_exp_scale);
    }

    /**
     * Expand magnetic scales (UNIQUE to magnetized nebulae)
     * @param {number} B_scale - Magnetic field scaling
     * @param {number} v_shock_scale - Shock velocity scaling
     */
    expandMagneticScale(B_scale, v_shock_scale) {
        this.setVariable('B', this.B * B_scale);
        this.setVariable('v_shock', this.v_shock * v_shock_scale);
    }

    // ----- Self-Refinement (3 methods) -----

    /**
     * Auto-refine parameters based on observations
     * @param {Array<Array<number>>} observations - Array of [time, g] pairs
     */
    autoRefineParameters(observations) {
        if (observations.length === 0) return;

        let total_error = 0.0;
        observations.forEach(([t, g_obs]) => {
            const g_model = this.compute_g_Crab(t);
            total_error += Math.abs(g_model - g_obs);
        });

        const avg_error = total_error / observations.length;
        if (avg_error > 1e-3) {
            const adjustment = 1.0 - (avg_error / (avg_error + 1.0)) * 0.1;
            this.setVariable('M', this.M * adjustment);
        }
    }

    /**
     * Calibrate to observational data
     * @param {Array<Array<number>>} obs_data - Array of [time, g] pairs
     */
    calibrateToObservations(obs_data) {
        this.autoRefineParameters(obs_data);
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
        const dt = (t_end - t_start) / (steps - 1);

        for (let i = 0; i < steps; i++) {
            const t = t_start + i * dt;
            const g = this.compute_g_Crab(t);
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
                // Skip constants
                if (['G', 'c', 'hbar', 'pi', 'M_sun', 'Mpc_to_m', 'year_to_s', 'm_e'].includes(varName)) {
                    variant[varName] = base;
                } else {
                    const variation = base * (1.0 + (Math.random() * 2 - 1) * variation_percent / 100.0);
                    variant[varName] = variation;
                }
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
        const mutable_vars = ['M', 'r0', 'v_exp', 'P_pulsar', 'B', 'v_shock', 'rho_fluid'];
        mutable_vars.forEach(varName => {
            const current = this.getVariable(varName);
            const mutated = current * (1.0 + (Math.random() * 2 - 1) * mutation_rate);
            this.setVariable(varName, mutated);
        });
    }

    /**
     * Evolve system over generations using fitness function
     * @param {number} generations - Number of generations
     * @param {Function} fitness - Fitness function (takes NebulaCrab instance)
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
        NebulaCrab.savedStates[label] = state;
    }

    /**
     * Restore state from a label
     * @param {string} label - State label
     * @returns {boolean} Success status
     */
    restoreState(label) {
        const state = NebulaCrab.savedStates[label];
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
        return Object.keys(NebulaCrab.savedStates);
    }

    /**
     * Export current state as string
     * @returns {string} State export string
     */
    exportState() {
        let output = "Crab Nebula (M1) State Export\n";
        output += "================================\n";

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
        const baseline = this.compute_g_Crab(t);

        const test_vars = ['M', 'r0', 'v_exp', 'P_pulsar', 'B', 'v_shock', 'rho_fluid'];
        test_vars.forEach(varName => {
            const original = this.getVariable(varName);
            this.setVariable(varName, original * (1.0 + perturbation));
            const perturbed = this.compute_g_Crab(t);
            sensitivities[varName] = Math.abs(perturbed - baseline) / (baseline + 1e-100);
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
        const yr_to_s = this.year_to_s;
        const r_t = this.computeRadius(t);

        let report = "========== CRAB NEBULA (M1) UQFF REPORT ==========\n";
        report += `Time: ${(t / yr_to_s).toFixed(1)} years\n`;
        report += `System: ${this.getSystemName()}\n\n`;
        
        report += "Key Parameters:\n";
        report += `  Nebula Mass M: ${(this.M / this.M_sun).toFixed(1)} M_sun\n`;
        report += `  Initial Radius r0: ${(this.r0 / 9.461e15).toFixed(2)} ly\n`;
        report += `  Current Radius r(t): ${(r_t / 9.461e15).toFixed(2)} ly\n`;
        report += `  Expansion Velocity: ${(this.v_exp / 1e3).toFixed(1)} km/s\n`;
        report += `  Pulsar Power: ${this.P_pulsar.toExponential(6)} W\n`;
        report += `  Nebula Density: ${this.rho_fluid.toExponential(6)} kg/m^3\n`;
        report += `  Magnetic Field: ${this.B.toExponential(6)} T\n`;
        report += `  Shock Velocity: ${(this.v_shock / 1e3).toFixed(1)} km/s\n\n`;
        
        const g = this.compute_g_Crab(t);
        report += `Computed g_UQFF: ${g.toExponential(6)} m/s^2\n`;
        report += "======================================================\n";
        return report;
    }

    /**
     * Validate parameter consistency
     * @returns {boolean} Validation status
     */
    validateConsistency() {
        let valid = true;
        if (this.M <= 0 || this.r0 <= 0) valid = false;
        if (this.v_exp < 0 || this.P_pulsar < 0) valid = false;
        if (this.B < 0 || this.v_shock < 0) valid = false;
        if (this.rho_fluid < 0 || this.B_crit <= 0) valid = false;
        return valid;
    }

    /**
     * Auto-correct anomalous parameters
     * @returns {boolean} Whether corrections were made
     */
    autoCorrectAnomalies() {
        let corrected = false;

        if (this.M <= 0) { this.M = 4.6 * this.M_sun; corrected = true; }
        if (this.r0 <= 0) { this.r0 = 5.2e16; corrected = true; }
        if (this.v_exp < 0) { this.v_exp = 1.5e6; corrected = true; }
        if (this.P_pulsar < 0) { this.P_pulsar = 5e31; corrected = true; }
        if (this.B < 0) { this.B = 1e-8; corrected = true; }
        if (this.v_shock < 0) { this.v_shock = 1.5e6; corrected = true; }
        if (this.rho_fluid < 0) { this.rho_fluid = 1e-21; corrected = true; }
        if (this.B_crit <= 0) { this.B_crit = 1e11; corrected = true; }
        if (this.M_visible < 0) { this.M_visible = this.M; corrected = true; }
        
        if (corrected) this.updateCache();
        return corrected;
    }
}

// ========== ENHANCED EXAMPLE FUNCTION ==========

function enhanced_crab_example() {
    console.log("=========================================================");
    console.log("ENHANCED CRAB NEBULA (M1) DEMONSTRATION");
    console.log("Pulsar-Driven Supernova Remnant with Expanding Radius");
    console.log("=========================================================\n");

    const crab = new NebulaCrab();
    const yr_to_s = crab.year_to_s;

    // Step 1: Initial state and validation
    console.log("Step 1: Initial State and Validation");
    console.log(`System: ${crab.getSystemName()}`);
    console.log(`Validation: ${crab.validateConsistency() ? 'PASS' : 'FAIL'}`);
    console.log(`Auto-corrected: ${crab.autoCorrectAnomalies() ? 'Yes' : 'No'}\n`);

    // Step 2: Time evolution showing r(t) expansion
    console.log("Step 2: Time Evolution (Radius Expansion)");
    const t_yr_array = [0, 500, 971, 1500, 2000];
    t_yr_array.forEach(t_yr => {
        const t = t_yr * yr_to_s;
        const r = crab.computeRadius(t);
        const g = crab.compute_g_Crab(t);
        console.log(`  t = ${t_yr.toFixed(0)} yr: r(t) = ${(r/9.461e15).toFixed(2)} ly, g = ${g.toExponential(6)} m/s²`);
    });
    console.log();

    // Step 3: Pulsar power scaling (UNIQUE to Crab)
    console.log("Step 3: Pulsar Power Scaling (P_pulsar sweeps) - PULSAR FEATURE");
    crab.saveState('original');
    const P_factors = [0.5, 1.0, 2.0];
    P_factors.forEach(factor => {
        crab.restoreState('original');
        crab.expandPulsarScale(factor, 1.0);
        const t = 971 * yr_to_s;
        const g = crab.compute_g_Crab(t);
        console.log(`  P_pulsar × ${factor}: g = ${g.toExponential(6)} m/s²`);
    });
    crab.restoreState('original');
    console.log();

    // Step 4: Magnetic field scaling (UNIQUE to Crab)
    console.log("Step 4: Magnetic Field Scaling (B sweeps) - MAGNETIC FEATURE");
    const B_factors = [0.5, 1.0, 2.0];
    B_factors.forEach(factor => {
        crab.restoreState('original');
        crab.expandMagneticScale(factor, 1.0);
        const t = 971 * yr_to_s;
        const g = crab.compute_g_Crab(t);
        console.log(`  B × ${factor}: g = ${g.toExponential(6)} m/s²`);
    });
    crab.restoreState('original');
    console.log();

    // Step 5: Generate report
    console.log("Step 5: Full System Report at 971 years");
    const t_report = 971 * yr_to_s;
    const report = crab.generateReport(t_report);
    console.log(report);

    console.log("=========================================================");
    console.log("ENHANCED DEMONSTRATION COMPLETE");
    console.log("=========================================================\n");
}

// Run inline test if executed directly
if (typeof require !== 'undefined' && require.main === module) {
    console.log("Running inline test...\n");
    enhanced_crab_example();
}

// Export for use as module
export default NebulaCrab;
