/**
 * ================================================================================================
 * Module: GalaxySombrero.js
 *
 * Description: JavaScript ES6 Module for Sombrero Galaxy (M104) UQFF Class
 *              This is the twenty-first module in a series of 500+ code files for the Universal
 *              Quantum Field Framework (UQFF) simulations, focusing on Sombrero galaxy evolution
 *              and gravity equations derived from Hubble datasets, high-energy lab simulations,
 *              and UQFF refinements (dated May 09, 2025, updated October 08, 2025).
 *
 * Purpose: Encapsulates the Master Universal Gravity Equation (MUGE) for Sombrero M104 evolution.
 *          Includes ALL terms: base gravity with expansion and superconductivity correction,
 *          black hole contribution, UQFF Ug components, Lambda, quantum uncertainty, EM/Lorentz,
 *          fluid dynamics, resonant waves, DM with density perturbations, and dust drag.
 *
 * Key Features:
 *   - Default values from UQFF document: M = 1×10¹¹ M☉, r = 2.36×10²⁰ m (25k ly),
 *     M_BH = 1×10⁹ M☉, z = 0.0063 (Virgo Cluster), v_orbit = 2×10⁵ m/s,
 *     B = 1×10⁻⁵ T, B_crit = 1×10¹¹ T (superconductivity threshold).
 *   - Prominent dust lane: Enhanced rho_dust term
 *   - Superconductivity correction: (1 - B/B_crit) on base gravity
 *   - DM fraction: 20% (bulge-dominated system)
 *   - Units handled: M☉ to kg, ly to m, all SI
 *   - Setter methods: setVariable, addToVariable, subtractFromVariable
 *   - Computes g_Sombrero(r, t) with every term explicitly included (10 major terms)
 *   - Enhanced dynamic capabilities: 38 methods total (13 core + 25 enhanced)
 *
 * Conversion: Ported from C++ source29.cpp (SombreroUQFFModule.h) to JavaScript ES6
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript
 * Date: November 03, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

class GalaxySombrero {
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
        
        // Sombrero galaxy parameters
        const M_sun = 1.989e30;                   // Solar mass (kg)
        this.M_sun = M_sun;
        this.M = 1e11 * M_sun;                    // Total mass (kg) - includes DM
        this.M_visible = 0.8 * this.M;            // Visible mass fraction (bulge/arms dominant)
        this.M_DM = 0.2 * this.M;                 // Dark matter mass (halo, lower fraction)
        this.r = 2.36e20;                         // Radius (m) - half diameter ~25k ly
        this.M_BH = 1e9 * M_sun;                  // Supermassive black hole (kg)
        this.r_BH = 1e15;                         // Core scale (m)
        
        // Hubble/cosmology
        this.H0 = 70.0;                           // Hubble constant (km/s/Mpc)
        this.Mpc_to_m = 3.086e22;                 // Megaparsec to meters
        this.z = 0.0063;                          // Redshift (Virgo Cluster)
        this.Omega_m = 0.3;                       // Matter density parameter
        this.Omega_Lambda = 0.7;                  // Dark energy density parameter
        this.t = 10e9 * 3.156e7;                  // Default time = 10 Gyr (s)
        
        // Dust/fluid dynamics (prominent dust lane)
        this.rho_dust = 1e-20;                    // Dust density (kg/m³)
        this.v_orbit = 2e5;                       // Orbital velocity (m/s)
        this.rho_mass = 1e-21;                    // Mean mass density (kg/m³) - ISM
        this.rho_fluid = 1e-21;                   // Fluid density (kg/m³) - dust lane-like
        this.V = 1e3;                             // Volume scale (m³)
        
        // EM/magnetic/superconductivity
        this.B = 1e-5;                            // Galactic magnetic field (T)
        this.B_crit = 1e15 * 1e-4;                // Critical field (10¹⁵ G ≈ 1e11 T)
        
        // Quantum terms
        this.Delta_x = 1e-10;                     // Position uncertainty (m)
        this.Delta_p = this.hbar / this.Delta_x;  // Momentum uncertainty (kg·m/s)
        this.integral_psi = 1.0;                  // Normalized wavefunction integral
        
        // Resonant/oscillatory terms
        this.A = 1e-10;                           // Oscillatory amplitude (m/s²)
        this.k = 1e20;                            // Wave number (m⁻¹)
        this.omega = 1e15;                        // Angular frequency (rad/s)
        this.x = 0.0;                             // Position (m) - central
        
        // DM perturbations
        this.delta_rho = 0.1 * this.rho_mass;     // Density perturbation (kg/m³)
        this.rho = this.rho_mass;                 // Mean density (kg/m³)
        
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
     * Dust term: ρ_dust × v_orbit² / ρ_mass × scale_macro
     * @returns {number} Dust drag acceleration (m/s²)
     */
    computeDustTerm() {
        const force_dust = this.rho_dust * (this.v_orbit * this.v_orbit);
        return (force_dust / this.rho_mass) * this.scale_macro;
    }

    /**
     * Main MUGE computation - includes ALL 10 major terms with superconductivity
     * @param {number} t - Time (s)
     * @returns {number} Total gravitational acceleration g_Sombrero (m/s²)
     */
    compute_g_Sombrero(t) {
        if (t < 0) {
            console.error("Error: Time t must be non-negative.");
            return 0.0;
        }

        this.t = t;  // Update time
        const Hz = this.computeHz();
        const expansion = 1.0 + Hz * t;
        const sc_correction = 1.0 - (this.B / this.B_crit);  // Superconductivity
        const tr_factor = 1.0 + this.f_TRZ;

        // Term 1: Base gravity with expansion, SC correction, and time-reversal
        const g_base = ((this.G * this.M / (this.r * this.r)) * expansion * sc_correction) * tr_factor;

        // Term 2: Black hole contribution
        const g_BH = (this.G * this.M_BH) / (this.r_BH * this.r_BH);

        // Term 3: UQFF Ug sum (Ug1 + Ug2 + Ug3 + Ug4)
        const ug_sum = this.computeUgSum();

        // Term 4: Cosmological constant
        const lambda_term = (this.Lambda * this.c * this.c) / 3.0;

        // Term 5: Quantum uncertainty
        const quantum_term = this.computeQuantumTerm(this.t_Hubble);

        // Term 6: EM Lorentz (q v × B, with UA/SCm correction)
        const em_base = (this.q * this.v_orbit * this.B) / 1.673e-27;  // / proton mass for acceleration
        const corr_UA = 1.0 + (this.rho_vac_UA / this.rho_vac_SCm);  // = 1 + 10 = 11
        const em_term = em_base * corr_UA * this.scale_macro;

        // Term 7: Fluid dynamics
        const fluid_term = this.computeFluidTerm(g_base);

        // Term 8: Resonant oscillatory
        const resonant_term = this.computeResonantTerm(t);

        // Term 9: DM with density perturbations (normalized by M for acceleration)
        const dm_force_like = this.computeDMTerm();
        const dm_term = dm_force_like / this.M;  // Convert force-like to acceleration

        // Term 10: Dust drag
        const dust_term = this.computeDustTerm();

        // Total g_Sombrero: sum all 10 terms
        return g_base + g_BH + ug_sum + lambda_term + quantum_term + em_term + fluid_term + resonant_term + dm_term + dust_term;
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
                this.M_visible = 0.8 * newValue;
                this.M_DM = 0.2 * newValue;
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
     * Example computation at t=10 Gyr (for testing)
     * @returns {number} g_Sombrero at 10 Gyr (m/s²)
     */
    exampleAt10Gyr() {
        const t_example = 10e9 * 3.156e7;
        return this.compute_g_Sombrero(t_example);
    }

    /**
     * Get equation text (descriptive)
     * @returns {string} Equation description
     */
    getEquationText() {
        return "g_Sombrero(r, t) = (G × M / r²) × (1 + H(z) × t) × (1 - B / B_crit) × (1 + f_TRZ) + " +
               "(G × M_BH / r_BH²) + (Ug1 + Ug2 + Ug3 + Ug4) + (Λ × c² / 3) + " +
               "(ℏ / √(Δx × Δp)) × ∫(ψ* H ψ dV) × (2π / t_Hubble) + q (v × B) + " +
               "ρ_fluid × V × g + 2A cos(k x) cos(ω t) + (2π / 13.8) A exp(i (k x - ω t)) + " +
               "(M_visible + M_DM) × (δρ/ρ + 3 G M / r³) + D_dust";
    }

    /**
     * Print parameters to console
     */
    printParameters() {
        console.log("Sombrero Galaxy (M104) Parameters:");
        console.log(`G: ${this.G.toExponential(3)}, M: ${(this.M/this.M_sun).toExponential(3)} M_sun, r: ${(this.r/9.461e15).toFixed(0)} ly`);
        console.log(`M_visible: ${(this.M_visible/this.M_sun).toExponential(3)} M_sun, M_DM: ${(this.M_DM/this.M_sun).toExponential(3)} M_sun`);
        console.log(`M_BH: ${(this.M_BH/this.M_sun).toExponential(3)} M_sun, r_BH: ${this.r_BH.toExponential(3)} m`);
        console.log(`z: ${this.z.toFixed(4)}, v_orbit: ${(this.v_orbit/1e3).toFixed(1)} km/s`);
        console.log(`B: ${this.B.toExponential(3)} T, B_crit: ${this.B_crit.toExponential(3)} T`);
        console.log(`SC correction: ${(1 - this.B/this.B_crit).toFixed(6)}`);
        console.log(`rho_dust: ${this.rho_dust.toExponential(3)} kg/m^3, Hz: ${this.computeHz().toExponential(3)} s^-1`);
        console.log(`Ug1: ${this.Ug1.toExponential(3)} m/s^2`);
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
        return "SombreroGalaxy_M104";
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

    // ----- Self-Expansion (4 methods - domain-specific for Sombrero galaxy) -----

    /**
     * Expand entire parameter space by a scale factor
     * @param {number} scale_factor - Scaling factor for all scalable parameters
     */
    expandParameterSpace(scale_factor) {
        const scalable = ['M', 'r', 'M_BH', 'r_BH', 'rho_dust', 'rho_fluid', 'v_orbit', 'B'];
        this.scaleVariableGroup(scalable, scale_factor);
    }

    /**
     * Expand galaxy mass and radius scales
     * @param {number} M_scale - Mass scaling factor
     * @param {number} r_scale - Radius scaling factor
     */
    expandGalaxyScale(M_scale, r_scale) {
        this.setVariable('M', this.M * M_scale);
        this.setVariable('r', this.r * r_scale);
    }

    /**
     * Expand black hole scales (UNIQUE to galaxies with SMBH)
     * @param {number} M_BH_scale - SMBH mass scaling
     * @param {number} r_BH_scale - SMBH core radius scaling
     */
    expandBlackHoleScale(M_BH_scale, r_BH_scale) {
        this.setVariable('M_BH', this.M_BH * M_BH_scale);
        this.setVariable('r_BH', this.r_BH * r_BH_scale);
    }

    /**
     * Expand dust lane scales (UNIQUE to Sombrero's prominent dust lane)
     * @param {number} rho_dust_scale - Dust density scaling
     * @param {number} B_scale - Magnetic field scaling
     */
    expandDustLaneScale(rho_dust_scale, B_scale) {
        this.setVariable('rho_dust', this.rho_dust * rho_dust_scale);
        this.setVariable('B', this.B * B_scale);
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
            const g_model = this.compute_g_Sombrero(t);
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
            const g = this.compute_g_Sombrero(t);
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
                if (['G', 'c', 'hbar', 'pi', 'M_sun', 'Mpc_to_m'].includes(varName)) {
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
        const mutable_vars = ['M', 'r', 'M_BH', 'r_BH', 'rho_dust', 'v_orbit', 'B'];
        mutable_vars.forEach(varName => {
            const current = this.getVariable(varName);
            const mutated = current * (1.0 + (Math.random() * 2 - 1) * mutation_rate);
            this.setVariable(varName, mutated);
        });
    }

    /**
     * Evolve system over generations using fitness function
     * @param {number} generations - Number of generations
     * @param {Function} fitness - Fitness function (takes GalaxySombrero instance)
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
        GalaxySombrero.savedStates[label] = state;
    }

    /**
     * Restore state from a label
     * @param {string} label - State label
     * @returns {boolean} Success status
     */
    restoreState(label) {
        const state = GalaxySombrero.savedStates[label];
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
        return Object.keys(GalaxySombrero.savedStates);
    }

    /**
     * Export current state as string
     * @returns {string} State export string
     */
    exportState() {
        let output = "Sombrero Galaxy (M104) State Export\n";
        output += "====================================\n";

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
        const baseline = this.compute_g_Sombrero(t);

        const test_vars = ['M', 'r', 'M_BH', 'r_BH', 'rho_dust', 'v_orbit', 'B', 'B_crit'];
        test_vars.forEach(varName => {
            const original = this.getVariable(varName);
            this.setVariable(varName, original * (1.0 + perturbation));
            const perturbed = this.compute_g_Sombrero(t);
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
        const Gyr_to_s = 1e9 * 3.156e7;
        const ly_to_m = 9.461e15;

        let report = "========== SOMBRERO GALAXY (M104) UQFF REPORT ==========\n";
        report += `Time: ${(t / Gyr_to_s).toFixed(6)} Gyr\n`;
        report += `System: ${this.getSystemName()}\n\n`;
        
        report += "Key Parameters:\n";
        report += `  Total Mass M: ${(this.M / this.M_sun).toExponential(6)} M_sun\n`;
        report += `  Visible Mass: ${(this.M_visible / this.M_sun).toExponential(6)} M_sun\n`;
        report += `  Dark Matter: ${(this.M_DM / this.M_sun).toExponential(6)} M_sun\n`;
        report += `  Radius r: ${(this.r / ly_to_m).toFixed(0)} ly\n`;
        report += `  SMBH Mass: ${(this.M_BH / this.M_sun).toExponential(6)} M_sun\n`;
        report += `  Dust Density: ${this.rho_dust.toExponential(6)} kg/m^3\n`;
        report += `  Magnetic Field: ${this.B.toExponential(6)} T\n`;
        report += `  SC Correction: ${(1.0 - this.B / this.B_crit).toFixed(6)}\n\n`;
        
        const g = this.compute_g_Sombrero(t);
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
        if (this.M_BH < 0) valid = false;
        if (this.rho_dust < 0 || this.B < 0) valid = false;
        if (this.B_crit <= 0) valid = false;
        return valid;
    }

    /**
     * Auto-correct anomalous parameters
     * @returns {boolean} Whether corrections were made
     */
    autoCorrectAnomalies() {
        let corrected = false;
        const M_sun = 1.989e30;

        if (this.M <= 0) { this.M = 1e11 * M_sun; corrected = true; }
        if (this.r <= 0) { this.r = 2.36e20; corrected = true; }
        if (this.M_BH < 0) { this.M_BH = 1e9 * M_sun; corrected = true; }
        if (this.rho_dust < 0) { this.rho_dust = 1e-20; corrected = true; }
        if (this.B < 0) { this.B = 1e-5; corrected = true; }
        if (this.B_crit <= 0) { this.B_crit = 1e15 * 1e-4; corrected = true; }
        if (this.M_visible < 0) { this.M_visible = 0.8 * this.M; corrected = true; }
        if (this.M_DM < 0) { this.M_DM = 0.2 * this.M; corrected = true; }
        
        if (corrected) this.updateCache();
        return corrected;
    }
}

// ========== ENHANCED EXAMPLE FUNCTION ==========

function enhanced_sombrero_example() {
    console.log("=========================================================");
    console.log("ENHANCED SOMBRERO GALAXY (M104) DEMONSTRATION");
    console.log("The Hat - Prominent Dust Lane & Superconductivity");
    console.log("=========================================================\n");

    const sombrero = new GalaxySombrero();
    const M_sun = 1.989e30;
    const Gyr_to_s = 1e9 * 3.156e7;
    const ly_to_m = 9.461e15;

    // Step 1: Initial state and validation
    console.log("Step 1: Initial State and Validation");
    console.log(`System: ${sombrero.getSystemName()}`);
    console.log(`Validation: ${sombrero.validateConsistency() ? 'PASS' : 'FAIL'}`);
    console.log(`Auto-corrected: ${sombrero.autoCorrectAnomalies() ? 'Yes' : 'No'}\n`);

    // Step 2: Time evolution with SC correction
    console.log("Step 2: Time Evolution (Galactic Dynamics with SC)");
    const t_Gyr_array = [0.0, 2.5, 5.0, 10.0, 13.8];
    t_Gyr_array.forEach(t_Gyr => {
        const t = t_Gyr * Gyr_to_s;
        const g = sombrero.compute_g_Sombrero(t);
        const sc_corr = 1.0 - sombrero.B / sombrero.B_crit;
        console.log(`  t = ${t_Gyr.toFixed(1)} Gyr: g = ${g.toExponential(6)} m/s², SC = ${sc_corr.toFixed(6)}`);
    });
    console.log();

    // Step 3: Dust lane density scaling (UNIQUE to Sombrero)
    console.log("Step 3: Dust Lane Density Scaling (rho_dust sweeps) - DUST LANE FEATURE");
    sombrero.saveState('original');
    const rho_dust_factors = [0.5, 1.0, 2.0];
    rho_dust_factors.forEach(factor => {
        sombrero.restoreState('original');
        sombrero.expandDustLaneScale(factor, 1.0);
        const t = 10 * Gyr_to_s;
        const g = sombrero.compute_g_Sombrero(t);
        console.log(`  rho_dust × ${factor}: g(10 Gyr) = ${g.toExponential(6)} m/s²`);
    });
    sombrero.restoreState('original');
    console.log();

    // Step 4: Magnetic field scaling (affects SC correction)
    console.log("Step 4: Magnetic Field Scaling (B sweeps, affects SC)");
    const B_factors = [0.5, 1.0, 2.0];
    B_factors.forEach(factor => {
        sombrero.restoreState('original');
        sombrero.expandDustLaneScale(1.0, factor);
        const t = 10 * Gyr_to_s;
        const g = sombrero.compute_g_Sombrero(t);
        const sc_corr = 1.0 - sombrero.B / sombrero.B_crit;
        console.log(`  B × ${factor}: SC = ${sc_corr.toFixed(6)}, g(10 Gyr) = ${g.toExponential(6)} m/s²`);
    });
    sombrero.restoreState('original');
    console.log();

    // Step 5: Generate report
    console.log("Step 5: Full System Report at 10 Gyr");
    const t_report = 10 * Gyr_to_s;
    const report = sombrero.generateReport(t_report);
    console.log(report);

    console.log("=========================================================");
    console.log("ENHANCED DEMONSTRATION COMPLETE");
    console.log("=========================================================\n");
}

// Run inline test if executed directly
if (typeof require !== 'undefined' && require.main === module) {
    console.log("Running inline test...\n");
    enhanced_sombrero_example();
}

// Export for use as module
export default GalaxySombrero;
