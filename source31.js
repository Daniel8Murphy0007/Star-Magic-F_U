/**
 * ================================================================================================
 * Module: NebulaM16.js
 *
 * Description: JavaScript ES6 Module for M16 (Eagle Nebula) UQFF Class
 *              This is the twenty-third module in a series of 500+ code files for the Universal
 *              Quantum Field Framework (UQFF) simulations, focusing on M16 Eagle Nebula evolution
 *              and gravity equations derived from Hubble datasets, high-energy lab simulations,
 *              and UQFF refinements (dated May 09, 2025, updated October 09, 2025).
 *
 * Purpose: Encapsulates the Master Universal Gravity Equation (MUGE) for M16 Eagle Nebula.
 *          Includes ALL terms: base gravity with SC correction, UQFF Ug components, Lambda,
 *          quantum uncertainty, EM/Lorentz, fluid dynamics, resonant waves, DM with density
 *          perturbations, star formation M_sf(t), and radiation erosion E_rad(t).
 *
 * Key Features:
 *   - Default values from UQFF document: M = 1200 M☉, r = 3.31×10¹⁷ m (~35 ly),
 *     SFR = 1 M☉/yr, tau_erode = 3×10⁶ yr, v_gas = 1×10⁵ m/s, z = 0.0015.
 *   - Star Formation: M_sf(t) = (SFR × t_yr) / M0 boosts mass over time
 *   - Radiation Erosion: E_rad(t) = E_0 × (1 - exp(-t/τ)) reduces mass via photoevaporation
 *   - Combined mass factor: M(t) = M × (1 + M_sf(t)) × (1 - E_rad(t))
 *   - Dense gas dynamics: rho_fluid = 1×10⁻²⁰ kg/m³, v_gas = 1×10⁵ m/s
 *   - Units handled: All SI (kg, m, s, T)
 *   - Setter methods: setVariable, addToVariable, subtractFromVariable
 *   - Computes g_M16(r, t) with every term explicitly included (8 major terms)
 *   - Enhanced dynamic capabilities: 38 methods total (13 core + 25 enhanced)
 *
 * Conversion: Ported from C++ source31.cpp (M16UQFFModule.h) to JavaScript ES6
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript
 * Date: November 03, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

class NebulaM16 {
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
        this.q = 1.602e-19;                       // Proton charge (C)
        this.pi = Math.PI;                        // π
        this.t_Hubble = 13.8e9 * 3.156e7;         // Hubble time (s) - 13.8 Gyr
        this.year_to_s = 3.156e7;                 // Seconds per year
        
        // M16 nebula parameters
        this.M_sun = 1.989e30;                    // Solar mass (kg)
        this.M = 1200 * this.M_sun;               // Total initial mass (kg)
        this.M0 = this.M;                         // Initial mass for SFR
        this.SFR = 1 * this.M_sun;                // Star formation rate (M☉/yr)
        this.M_visible = this.M;                  // Visible mass (gas + stars)
        this.M_DM = 0.0;                          // No significant DM
        this.r = 3.31e17;                         // Half span (~35 ly) (m)
        
        // Hubble/cosmology
        this.H0 = 70.0;                           // Hubble constant (km/s/Mpc)
        this.Mpc_to_m = 3.086e22;                 // Megaparsec to meters
        this.z = 0.0015;                          // Redshift
        this.Omega_m = 0.3;                       // Matter density parameter
        this.Omega_Lambda = 0.7;                  // Dark energy density parameter
        this.t = 5e6 * this.year_to_s;            // Default time = 5 Myr (s)
        
        // Gas dynamics
        this.rho_fluid = 1e-20;                   // Dense gas (kg/m³)
        this.V = 1e3;                             // Volume scale (m³)
        this.v_gas = 1e5;                         // Gas velocity (m/s)
        this.delta_rho = 0.1 * this.rho_fluid;    // Density perturbation (kg/m³)
        this.rho = this.rho_fluid;                // Mean density (kg/m³)
        
        // EM/magnetic/superconductivity
        this.B = 1e-5;                            // Nebula magnetic field (T)
        this.B_crit = 1e11;                       // Critical field (10¹⁵ G ≈ 1e11 T)
        
        // Quantum terms
        this.Delta_x = 1e-10;                     // Position uncertainty (m)
        this.Delta_p = this.hbar / this.Delta_x;  // Momentum uncertainty (kg·m/s)
        this.integral_psi = 1.0;                  // Normalized wavefunction integral
        
        // Resonant/oscillatory terms
        this.A = 1e-10;                           // Oscillatory amplitude (m/s²)
        this.k = 1e20;                            // Wave number (m⁻¹)
        this.omega = 1e15;                        // Angular frequency (rad/s)
        this.x = 0.0;                             // Position (m) - central
        
        // Star formation and erosion
        this.tau_erode_yr = 3e6;                  // Erosion timescale (yr)
        this.E_0 = 0.3;                           // Fractional erosion max
        
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
        this.Ug1 = (this.G * this.M) / (this.r * this.r);
        this.Ug4 = this.Ug1 * this.f_sc;
        this.Delta_p = this.hbar / this.Delta_x;
    }

    // ========== CORE PHYSICS METHODS (13 methods) ==========

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
     * Ug1 = G M / r², Ug4 = Ug1 × f_sc, others ≈ 0
     * @returns {number} Total Ug (m/s²)
     */
    computeUgSum() {
        this.Ug1 = (this.G * this.M) / (this.r * this.r);
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
     * @returns {number} DM contribution (force-like, needs normalization by mass)
     */
    computeDMTerm() {
        const pert = this.delta_rho / this.rho;
        const curv = (3 * this.G * this.M) / (this.r * this.r * this.r);
        return (this.M_visible + this.M_DM) * (pert + curv);
    }

    /**
     * Star formation factor: (SFR × t_yr) / M0
     * @param {number} t - Time (s)
     * @returns {number} M_sf factor (unitless)
     */
    computeMsfFactor(t) {
        const t_yr = t / this.year_to_s;
        return (this.SFR / this.M0) * t_yr;
    }

    /**
     * Radiation erosion factor: E_0 × (1 - exp(-t/τ))
     * @param {number} t - Time (s)
     * @returns {number} E_rad factor (unitless)
     */
    computeE_rad(t) {
        const tau_s = this.tau_erode_yr * this.year_to_s;
        return this.E_0 * (1.0 - Math.exp(-t / tau_s));
    }

    /**
     * Main MUGE computation - includes ALL 8 major terms with M(t)
     * @param {number} t - Time (s)
     * @returns {number} Total gravitational acceleration g_M16 (m/s²)
     */
    compute_g_M16(t) {
        if (t < 0) {
            console.error("Error: Time t must be non-negative.");
            return 0.0;
        }

        this.t = t;  // Update time
        const Hz = this.computeHz();
        const expansion = 1.0 + Hz * t;
        const sc_correction = 1.0 - (this.B / this.B_crit);
        const tr_factor = 1.0 + this.f_TRZ;
        const msf_factor = this.computeMsfFactor(t);
        const e_rad = this.computeE_rad(t);
        const m_factor = (1.0 + msf_factor) * (1.0 - e_rad);

        // Term 1: Base gravity with expansion, SC, TR, M_sf, E_rad
        const g_base = ((this.G * this.M * m_factor) / (this.r * this.r)) * expansion * sc_correction * tr_factor;

        // Term 2: UQFF Ug sum (Ug1 + Ug2 + Ug3 + Ug4)
        const ug_sum = this.computeUgSum();

        // Term 3: Cosmological constant
        const lambda_term = (this.Lambda * this.c * this.c) / 3.0;

        // Term 4: Quantum uncertainty
        const quantum_term = this.computeQuantumTerm(this.t_Hubble);

        // Term 5: EM Lorentz (q v × B, with UA/SCm correction)
        const em_base = (this.q * this.v_gas * this.B) / 1.673e-27;  // / proton mass for acceleration
        const corr_UA = 1.0 + (this.rho_vac_UA / this.rho_vac_SCm);  // = 1 + 10 = 11
        const em_term = em_base * corr_UA * this.scale_macro;

        // Term 6: Fluid dynamics
        const fluid_term = this.computeFluidTerm(g_base);

        // Term 7: Resonant oscillatory
        const resonant_term = this.computeResonantTerm(t);

        // Term 8: DM with density perturbations (normalized by M for acceleration)
        const dm_force_like = this.computeDMTerm();
        const dm_term = dm_force_like / this.M;  // Convert force-like to acceleration

        // Total g_M16: sum all 8 terms (erosion already in m_factor)
        return g_base + ug_sum + lambda_term + quantum_term + em_term + fluid_term + resonant_term + dm_term;
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
                this.M0 = newValue;
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
     * Example computation at t=5 Myr (for testing)
     * @returns {number} g_M16 at 5 Myr (m/s²)
     */
    exampleAt5Myr() {
        const t_example = 5e6 * this.year_to_s;
        return this.compute_g_M16(t_example);
    }

    /**
     * Get equation text (descriptive)
     * @returns {string} Equation description
     */
    getEquationText() {
        return "g_M16(r, t) = (G × M(t) / r²) × (1 + H(z) × t) × (1 - B / B_crit) × (1 + f_TRZ) + " +
               "(Ug1 + Ug2 + Ug3 + Ug4) + (Λ × c² / 3) + " +
               "(ℏ / √(Δx × Δp)) × ∫(ψ* H ψ dV) × (2π / t_Hubble) + q (v × B) + " +
               "ρ_fluid × V × g + 2A cos(k x) cos(ω t) + (2π / 13.8) A exp(i (k x - ω t)) + " +
               "(M_visible + M_DM) × (δρ/ρ + 3 G M / r³)\n" +
               "Where M(t) = M × (1 + M_sf(t)) × (1 - E_rad(t)); M_sf(t) = (SFR × t_yr) / M0; " +
               "E_rad(t) = E_0 × (1 - exp(-t / τ))";
    }

    /**
     * Print parameters to console
     */
    printParameters() {
        console.log("Eagle Nebula (M16) Parameters:");
        console.log(`G: ${this.G.toExponential(3)}, M: ${(this.M/this.M_sun).toFixed(1)} M_sun, r: ${(this.r/9.461e15).toFixed(1)} ly`);
        console.log(`SFR: ${(this.SFR/this.M_sun).toFixed(1)} M_sun/yr, tau_erode: ${(this.tau_erode_yr/1e6).toFixed(1)} Myr`);
        console.log(`z: ${this.z.toFixed(4)}, rho_fluid: ${this.rho_fluid.toExponential(3)} kg/m³`);
        console.log(`v_gas: ${(this.v_gas/1e3).toFixed(1)} km/s, B: ${this.B.toExponential(3)} T`);
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
        return "EagleNebula_M16";
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

    // ----- Self-Expansion (4 methods - domain-specific for M16) -----

    /**
     * Expand entire parameter space by a scale factor
     * @param {number} scale_factor - Scaling factor for all scalable parameters
     */
    expandParameterSpace(scale_factor) {
        const scalable = ['M', 'r', 'SFR', 'tau_erode_yr', 'rho_fluid', 'v_gas', 'B'];
        this.scaleVariableGroup(scalable, scale_factor);
    }

    /**
     * Expand nebula mass and radius scales
     * @param {number} M_scale - Mass scaling factor
     * @param {number} r_scale - Radius scaling factor
     */
    expandNebulaScale(M_scale, r_scale) {
        this.setVariable('M', this.M * M_scale);
        this.setVariable('r', this.r * r_scale);
    }

    /**
     * Expand star formation scales (UNIQUE to star-forming nebulae)
     * @param {number} SFR_scale - Star formation rate scaling
     * @param {number} tau_erode_scale - Erosion timescale scaling
     */
    expandStarFormationScale(SFR_scale, tau_erode_scale) {
        this.setVariable('SFR', this.SFR * SFR_scale);
        this.setVariable('tau_erode_yr', this.tau_erode_yr * tau_erode_scale);
    }

    /**
     * Expand gas dynamics scales (UNIQUE to nebular gas)
     * @param {number} rho_fluid_scale - Gas density scaling
     * @param {number} v_gas_scale - Gas velocity scaling
     */
    expandGasScale(rho_fluid_scale, v_gas_scale) {
        this.setVariable('rho_fluid', this.rho_fluid * rho_fluid_scale);
        this.setVariable('v_gas', this.v_gas * v_gas_scale);
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
            const g_model = this.compute_g_M16(t);
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
            const g = this.compute_g_M16(t);
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
                if (['G', 'c', 'hbar', 'pi', 'M_sun', 'Mpc_to_m', 'year_to_s'].includes(varName)) {
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
        const mutable_vars = ['M', 'r', 'SFR', 'tau_erode_yr', 'rho_fluid', 'v_gas', 'B'];
        mutable_vars.forEach(varName => {
            const current = this.getVariable(varName);
            const mutated = current * (1.0 + (Math.random() * 2 - 1) * mutation_rate);
            this.setVariable(varName, mutated);
        });
    }

    /**
     * Evolve system over generations using fitness function
     * @param {number} generations - Number of generations
     * @param {Function} fitness - Fitness function (takes NebulaM16 instance)
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
        NebulaM16.savedStates[label] = state;
    }

    /**
     * Restore state from a label
     * @param {string} label - State label
     * @returns {boolean} Success status
     */
    restoreState(label) {
        const state = NebulaM16.savedStates[label];
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
        return Object.keys(NebulaM16.savedStates);
    }

    /**
     * Export current state as string
     * @returns {string} State export string
     */
    exportState() {
        let output = "Eagle Nebula (M16) State Export\n";
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
        const baseline = this.compute_g_M16(t);

        const test_vars = ['M', 'r', 'SFR', 'tau_erode_yr', 'rho_fluid', 'v_gas', 'B', 'E_0'];
        test_vars.forEach(varName => {
            const original = this.getVariable(varName);
            this.setVariable(varName, original * (1.0 + perturbation));
            const perturbed = this.compute_g_M16(t);
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
        const Myr_to_s = 1e6 * this.year_to_s;

        let report = "========== EAGLE NEBULA (M16) UQFF REPORT ==========\n";
        report += `Time: ${(t / Myr_to_s).toFixed(6)} Myr\n`;
        report += `System: ${this.getSystemName()}\n\n`;
        
        report += "Key Parameters:\n";
        report += `  Nebula Mass M: ${(this.M / this.M_sun).toFixed(1)} M_sun\n`;
        report += `  Nebula Radius r: ${(this.r / 9.461e15).toFixed(1)} ly\n`;
        report += `  Star Formation Rate: ${(this.SFR / this.M_sun).toFixed(1)} M_sun/yr\n`;
        report += `  Erosion Timescale: ${(this.tau_erode_yr / 1e6).toFixed(1)} Myr\n`;
        report += `  Gas Density: ${this.rho_fluid.toExponential(6)} kg/m^3\n`;
        report += `  Gas Velocity: ${(this.v_gas / 1e3).toFixed(1)} km/s\n`;
        report += `  Magnetic Field: ${this.B.toExponential(6)} T\n\n`;
        
        const msf_factor = this.computeMsfFactor(t);
        const e_rad = this.computeE_rad(t);
        report += `Star Formation Factor: ${msf_factor.toFixed(6)}\n`;
        report += `Erosion Factor: ${e_rad.toFixed(6)}\n\n`;
        
        const g = this.compute_g_M16(t);
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
        if (this.M <= 0 || this.r <= 0) valid = false;
        if (this.SFR < 0 || this.tau_erode_yr <= 0) valid = false;
        if (this.rho_fluid < 0 || this.v_gas < 0) valid = false;
        if (this.B < 0 || this.B_crit <= 0) valid = false;
        return valid;
    }

    /**
     * Auto-correct anomalous parameters
     * @returns {boolean} Whether corrections were made
     */
    autoCorrectAnomalies() {
        let corrected = false;

        if (this.M <= 0) { this.M = 1200 * this.M_sun; corrected = true; }
        if (this.r <= 0) { this.r = 3.31e17; corrected = true; }
        if (this.SFR < 0) { this.SFR = 1 * this.M_sun; corrected = true; }
        if (this.tau_erode_yr <= 0) { this.tau_erode_yr = 3e6; corrected = true; }
        if (this.rho_fluid < 0) { this.rho_fluid = 1e-20; corrected = true; }
        if (this.v_gas < 0) { this.v_gas = 1e5; corrected = true; }
        if (this.B < 0) { this.B = 1e-5; corrected = true; }
        if (this.B_crit <= 0) { this.B_crit = 1e11; corrected = true; }
        if (this.M_visible < 0) { this.M_visible = this.M; corrected = true; }
        if (this.M0 <= 0) { this.M0 = this.M; corrected = true; }
        
        if (corrected) this.updateCache();
        return corrected;
    }
}

// ========== ENHANCED EXAMPLE FUNCTION ==========

function enhanced_m16_example() {
    console.log("=========================================================");
    console.log("ENHANCED EAGLE NEBULA (M16) DEMONSTRATION");
    console.log("Pillars of Creation - Star Formation & Erosion");
    console.log("=========================================================\n");

    const m16 = new NebulaM16();
    const Myr_to_s = 1e6 * m16.year_to_s;

    // Step 1: Initial state and validation
    console.log("Step 1: Initial State and Validation");
    console.log(`System: ${m16.getSystemName()}`);
    console.log(`Validation: ${m16.validateConsistency() ? 'PASS' : 'FAIL'}`);
    console.log(`Auto-corrected: ${m16.autoCorrectAnomalies() ? 'Yes' : 'No'}\n`);

    // Step 2: Time evolution showing M(t) with SFR and erosion
    console.log("Step 2: Time Evolution (Star Formation vs Erosion)");
    const t_Myr_array = [0.0, 1.0, 3.0, 5.0, 10.0];
    t_Myr_array.forEach(t_Myr => {
        const t = t_Myr * Myr_to_s;
        const msf = m16.computeMsfFactor(t);
        const e_rad = m16.computeE_rad(t);
        const g = m16.compute_g_M16(t);
        console.log(`  t = ${t_Myr.toFixed(1)} Myr: M_sf = ${msf.toFixed(3)}, E_rad = ${e_rad.toFixed(3)}, g = ${g.toExponential(6)} m/s²`);
    });
    console.log();

    // Step 3: Star formation rate scaling (UNIQUE to M16)
    console.log("Step 3: Star Formation Rate Scaling (SFR sweeps) - SF FEATURE");
    m16.saveState('original');
    const SFR_factors = [0.5, 1.0, 2.0];
    SFR_factors.forEach(factor => {
        m16.restoreState('original');
        m16.expandStarFormationScale(factor, 1.0);
        const t = 5 * Myr_to_s;
        const msf = m16.computeMsfFactor(t);
        const g = m16.compute_g_M16(t);
        console.log(`  SFR × ${factor}: M_sf(5 Myr) = ${msf.toFixed(3)}, g = ${g.toExponential(6)} m/s²`);
    });
    m16.restoreState('original');
    console.log();

    // Step 4: Erosion timescale scaling (UNIQUE to M16)
    console.log("Step 4: Erosion Timescale Scaling (tau_erode sweeps) - EROSION FEATURE");
    const tau_factors = [0.5, 1.0, 2.0];
    tau_factors.forEach(factor => {
        m16.restoreState('original');
        m16.expandStarFormationScale(1.0, factor);
        const t = 5 * Myr_to_s;
        const e_rad = m16.computeE_rad(t);
        const g = m16.compute_g_M16(t);
        console.log(`  tau × ${factor}: E_rad(5 Myr) = ${e_rad.toFixed(3)}, g = ${g.toExponential(6)} m/s²`);
    });
    m16.restoreState('original');
    console.log();

    // Step 5: Generate report
    console.log("Step 5: Full System Report at 5 Myr");
    const t_report = 5 * Myr_to_s;
    const report = m16.generateReport(t_report);
    console.log(report);

    console.log("=========================================================");
    console.log("ENHANCED DEMONSTRATION COMPLETE");
    console.log("=========================================================\n");
}

// Run inline test if executed directly
if (typeof require !== 'undefined' && require.main === module) {
    console.log("Running inline test...\n");
    enhanced_m16_example();
}

// Export for use as module
export default NebulaM16;
