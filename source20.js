/**
 * ================================================================================================
 * Module: source20.js (GalaxyNGC2525)
 *
 * Description: JavaScript Module for Galaxy NGC 2525 Class
 *              This is the tenth module in a series of 500+ code files for the Universal Quantum
 *              Field Framework (UQFF) simulations, focusing on barred spiral galaxy evolution and
 *              gravity equations derived from Hubble datasets, high-energy lab simulations, and
 *              UQFF refinements (dated May 09, 2025, updated for full term inclusion on October 08, 2025).
 *
 * Purpose: Encapsulates the Master Universal Gravity Equation (MUGE) for NGC 2525 evolution.
 *          Includes ALL terms: base gravity (static M), cosmic expansion (H(z)), magnetic correction (static B),
 *          black hole influence, UQFF Ug components with f_TRZ, Lambda, quantum uncertainty, scaled EM with [UA],
 *          fluid dynamics, oscillatory waves, DM/density perturbations, and supernova mass loss - (G * M_SN(t)) / r^2.
 *          Supports dynamic variable updates for all parameters.
 *
 * Integration: Designed for Node.js and browser environments.
 *              Instantiate class: const ngc2525 = new GalaxyNGC2525();
 *              Compute: const g = ngc2525.compute_g_NGC2525(t);
 *
 * Key Features:
 *   - Default values from UQFF document: M ≈ 1.0000225e10 Msun, r = 2.836e20 m, M_BH = 2.25e7 Msun,
 *     r_BH = 1.496e11 m, z ≈ 0.016, Hz ≈ 2.19e-18 s^-1, M_SN0 = 1.4 Msun, tau_SN = 1 yr, B = 1e-5 T.
 *   - Units handled: Msun to kg, ly to m; SN term as mass loss acceleration.
 *   - Setter methods for updates: setVariable(name, value) or addToVariable/subtractFromVariable.
 *   - Computes g_NGC2525(r, t) with every term explicitly included.
 *
 * Author: Converted to JavaScript from C++ by GitHub Copilot
 * Original Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript.
 * Date: November 03, 2025 (JS conversion)
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

class GalaxyNGC2525 {
    // Static property for state storage (replaces anonymous namespace)
    static savedStates = {};

    constructor() {
        // Core parameters (mutable for updates)
        this.G = 0;               // Gravitational constant
        this.M = 0;               // Total galaxy mass (kg)
        this.r = 0;               // Galaxy radius (m)
        this.Hz = 0;              // Hubble parameter at z (s^-1)
        this.B = 0;               // Static magnetic field (T)
        this.B_crit = 0;          // Critical B field (T)
        this.Lambda = 0;          // Cosmological constant
        this.c_light = 0;         // Speed of light
        this.q_charge = 0;        // Charge (proton)
        this.gas_v = 0;           // Gas velocity for EM (m/s)
        this.f_TRZ = 0;           // Time-reversal factor
        this.M_BH = 0;            // Black hole mass (kg)
        this.r_BH = 0;            // Black hole influence radius (m)
        this.M_SN0 = 0;           // Initial SN mass (kg)
        this.tau_SN = 0;          // SN decay timescale (s)
        this.rho_vac_UA = 0;      // UA vacuum density (J/m^3)
        this.rho_vac_SCm = 0;     // SCm vacuum density (J/m^3)
        this.scale_EM = 0;        // EM scaling factor
        this.proton_mass = 0;     // Proton mass for EM acceleration
        this.z_gal = 0;           // Galaxy redshift

        // Additional parameters for full inclusion of terms
        this.hbar = 0;            // Reduced Planck's constant
        this.t_Hubble = 0;        // Hubble time (s)
        this.delta_x = 0;         // Position uncertainty (m)
        this.delta_p = 0;         // Momentum uncertainty (kg m/s)
        this.integral_psi = 0;    // Wavefunction integral approximation
        this.rho_fluid = 0;       // Fluid density (kg/m^3)
        this.A_osc = 0;           // Oscillatory amplitude (m/s^2)
        this.k_osc = 0;           // Wave number (1/m)
        this.omega_osc = 0;       // Angular frequency (rad/s)
        this.x_pos = 0;           // Position for oscillation (m)
        this.t_Hubble_gyr = 0;    // Hubble time in Gyr
        this.M_DM_factor = 0;     // Dark matter mass fraction
        this.delta_rho_over_rho = 0; // Density perturbation fraction

        // Computed caches (updated on demand)
        this.ug1_base = 0;        // Cached Ug1 = G*M/r^2
        this.g_BH = 0;            // Cached BH acceleration

        this.initializeDefaults();
    }

    // Initialization method (called in constructor)
    initializeDefaults() {
        this.G = 6.6743e-11;
        const M_sun = 1.989e30;
        this.M = (1e10 + 2.25e7) * M_sun;
        this.r = 2.836e20;
        this.z_gal = 0.016;
        const Hz_kms = 70 * Math.sqrt(0.3 * Math.pow(1 + this.z_gal, 3) + 0.7);  // km/s/Mpc
        this.Hz = (Hz_kms * 1000 / 3.086e19);  // s^-1
        this.B = 1e-5;
        this.B_crit = 1e11;
        this.Lambda = 1.1e-52;
        this.c_light = 3e8;
        this.q_charge = 1.602e-19;
        this.gas_v = 1e5;
        this.f_TRZ = 0.1;
        this.M_BH = 2.25e7 * M_sun;
        this.r_BH = 1.496e11;
        this.M_SN0 = 1.4 * M_sun;
        this.tau_SN = 1 * 3.156e7;
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
        this.rho_fluid = 1e-21;
        this.A_osc = 1e-10;
        this.k_osc = 1.0 / this.r;
        this.omega_osc = 2 * Math.PI / (this.r / this.c_light);
        this.x_pos = this.r;
        this.M_DM_factor = 0.1;
        this.delta_rho_over_rho = 1e-5;

        this.updateCache();
    }

    // Cache update for efficiency (call after parameter changes)
    updateCache() {
        this.ug1_base = (this.G * this.M) / (this.r * this.r);
        this.g_BH = (this.G * this.M_BH) / (this.r_BH * this.r_BH);
    }

    // Universal setter for any variable (by name, for flexibility)
    setVariable(varName, newValue) {
        const varMap = {
            'G': 'G', 'M': 'M', 'r': 'r', 'Hz': 'Hz', 'B': 'B', 'B_crit': 'B_crit',
            'Lambda': 'Lambda', 'c_light': 'c_light', 'q_charge': 'q_charge', 'gas_v': 'gas_v',
            'f_TRZ': 'f_TRZ', 'M_BH': 'M_BH', 'r_BH': 'r_BH', 'M_SN0': 'M_SN0', 'tau_SN': 'tau_SN',
            'rho_vac_UA': 'rho_vac_UA', 'rho_vac_SCm': 'rho_vac_SCm', 'scale_EM': 'scale_EM',
            'proton_mass': 'proton_mass', 'z_gal': 'z_gal', 'hbar': 'hbar', 't_Hubble': 't_Hubble',
            't_Hubble_gyr': 't_Hubble_gyr', 'delta_x': 'delta_x', 'delta_p': 'delta_p',
            'integral_psi': 'integral_psi', 'rho_fluid': 'rho_fluid', 'A_osc': 'A_osc',
            'k_osc': 'k_osc', 'omega_osc': 'omega_osc', 'x_pos': 'x_pos',
            'M_DM_factor': 'M_DM_factor', 'delta_rho_over_rho': 'delta_rho_over_rho'
        };

        if (varMap[varName]) {
            this[varMap[varName]] = newValue;
            this.updateCache();
            return true;
        } else {
            console.error(`Error: Unknown variable '${varName}'.`);
            return false;
        }
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
        const varMap = {
            'G': this.G, 'M': this.M, 'r': this.r, 'Hz': this.Hz, 'B': this.B, 'B_crit': this.B_crit,
            'Lambda': this.Lambda, 'c_light': this.c_light, 'q_charge': this.q_charge, 'gas_v': this.gas_v,
            'f_TRZ': this.f_TRZ, 'M_BH': this.M_BH, 'r_BH': this.r_BH, 'M_SN0': this.M_SN0, 'tau_SN': this.tau_SN,
            'rho_vac_UA': this.rho_vac_UA, 'rho_vac_SCm': this.rho_vac_SCm, 'scale_EM': this.scale_EM,
            'proton_mass': this.proton_mass, 'z_gal': this.z_gal, 'hbar': this.hbar, 't_Hubble': this.t_Hubble,
            't_Hubble_gyr': this.t_Hubble_gyr, 'delta_x': this.delta_x, 'delta_p': this.delta_p,
            'integral_psi': this.integral_psi, 'rho_fluid': this.rho_fluid, 'A_osc': this.A_osc,
            'k_osc': this.k_osc, 'omega_osc': this.omega_osc, 'x_pos': this.x_pos,
            'M_DM_factor': this.M_DM_factor, 'delta_rho_over_rho': this.delta_rho_over_rho
        };

        if (varMap.hasOwnProperty(varName)) {
            return varMap[varName];
        } else {
            console.error(`Error: Unknown variable '${varName}'.`);
            return 0.0;
        }
    }

    // M_SN(t) computation - supernova mass exponential decay
    M_SN_t(t) {
        return this.M_SN0 * Math.exp(-t / this.tau_SN);
    }

    // Ug terms computation
    compute_Ug() {
        const Ug1 = this.ug1_base;
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
    compute_g_NGC2525(t) {
        if (t < 0) {
            console.error("Error: Time t must be non-negative.");
            return 0.0;
        }

        const MSNt = this.M_SN_t(t);

        // Term 1: Base + Hz + B corrections
        const corr_H = 1 + this.Hz * t;
        const corr_B = 1 - this.B / this.B_crit;
        const term1 = this.ug1_base * corr_H * corr_B;

        // BH term
        const term_BH = this.g_BH;

        // Term 2: UQFF Ug with f_TRZ
        const term2 = this.compute_Ug();

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
        const term_fluid = (this.rho_fluid * V * this.ug1_base) / this.M;

        // Oscillatory terms (real parts)
        const term_osc1 = 2 * this.A_osc * Math.cos(this.k_osc * this.x_pos) * Math.cos(this.omega_osc * t);
        const arg = this.k_osc * this.x_pos - this.omega_osc * t;
        const term_osc2 = (2 * Math.PI / this.t_Hubble_gyr) * this.A_osc * Math.cos(arg);
        const term_osc = term_osc1 + term_osc2;

        // DM and density perturbation term (converted to acceleration)
        const M_dm = this.M * this.M_DM_factor;
        const pert1 = this.delta_rho_over_rho;
        const pert2 = 3 * this.G * this.M / (this.r * this.r * this.r);
        const term_dm_force_like = (this.M + M_dm) * (pert1 + pert2);
        const term_DM = term_dm_force_like / this.M;

        // SN mass loss term (negative acceleration)
        const term_SN = -(this.G * MSNt) / (this.r * this.r);

        // Total g_NGC2525 (all terms summed)
        return term1 + term_BH + term2 + term3 + term4 + term_q + term_fluid + term_osc + term_DM + term_SN;
    }

    // Debug/Output method (for transparency)
    printParameters() {
        console.log("NGC 2525 Parameters:");
        console.log(`G: ${this.G.toFixed(3)}, M: ${this.M.toFixed(3)}, r: ${this.r.toFixed(3)}`);
        console.log(`Hz: ${this.Hz.toFixed(3)}, B: ${this.B.toFixed(3)}, B_crit: ${this.B_crit.toFixed(3)}`);
        console.log(`f_TRZ: ${this.f_TRZ.toFixed(3)}, M_BH: ${this.M_BH.toFixed(3)}, r_BH: ${this.r_BH.toFixed(3)}`);
        console.log(`M_SN0: ${this.M_SN0.toFixed(3)}, tau_SN: ${this.tau_SN.toFixed(3)}`);
        console.log(`rho_fluid: ${this.rho_fluid.toFixed(3)}, M_DM_factor: ${this.M_DM_factor.toFixed(3)}`);
        console.log(`A_osc: ${this.A_osc.toFixed(3)}, delta_rho_over_rho: ${this.delta_rho_over_rho.toFixed(3)}`);
        console.log(`ug1_base: ${this.ug1_base.toFixed(3)}, g_BH: ${this.g_BH.toFixed(3)}`);
    }

    // Example computation at t=7 years (for testing)
    exampleAt7Years() {
        const t_example = 7 * 3.156e7;
        return this.compute_g_NGC2525(t_example);
    }

    // ========== ENHANCED DYNAMIC CAPABILITIES (25 methods) ==========

    // --- Variable Management (5 methods) ---
    createVariable(name, value) {
        return this.setVariable(name, value);
    }

    removeVariable(name) {
        console.error(`Warning: Cannot remove built-in variable '${name}' in GalaxyNGC2525 class.`);
        return false;
    }

    cloneVariable(src, dest) {
        const val = this.getVariable(src);
        return this.setVariable(dest, val);
    }

    listVariables() {
        return ["G", "M", "r", "Hz", "B", "B_crit", "Lambda", "c_light", "q_charge",
                "gas_v", "f_TRZ", "M_BH", "r_BH", "M_SN0", "tau_SN", "rho_vac_UA",
                "rho_vac_SCm", "scale_EM", "proton_mass", "z_gal", "hbar", "t_Hubble",
                "t_Hubble_gyr", "delta_x", "delta_p", "integral_psi", "rho_fluid",
                "A_osc", "k_osc", "omega_osc", "x_pos", "M_DM_factor", "delta_rho_over_rho"];
    }

    getSystemName() {
        return "GalaxyNGC2525";
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
        const expandable = ["M", "r", "B", "rho_fluid", "A_osc", "M_DM_factor"];
        this.scaleVariableGroup(expandable, factor);
    }

    expandGalacticScale(M_scale, r_scale) {
        this.setVariable("M", this.getVariable("M") * M_scale);
        this.setVariable("r", this.getVariable("r") * r_scale);
    }

    expandBlackHoleScale(M_BH_scale, r_BH_scale) {
        this.setVariable("M_BH", this.getVariable("M_BH") * M_BH_scale);
        this.setVariable("r_BH", this.getVariable("r_BH") * r_BH_scale);
    }

    expandSupernovaScale(M_SN0_scale, tau_SN_scale) {
        this.setVariable("M_SN0", this.getVariable("M_SN0") * M_SN0_scale);
        this.setVariable("tau_SN", this.getVariable("tau_SN") * tau_SN_scale);
    }

    // --- Self-Refinement (3 methods) ---
    autoRefineParameters(observations) {
        if (observations.length === 0) return;
        
        let sum_error = 0.0;
        for (const obs of observations) {
            const t = obs[0];
            const g_obs = obs[1];
            const g_calc = this.compute_g_NGC2525(t);
            sum_error += Math.abs(g_calc - g_obs);
        }
        const avg_error = sum_error / observations.length;
        
        if (avg_error > 1e-6) {
            const adj_factor = 1.0 - Math.min(0.1, avg_error / 1e6);
            this.setVariable("M", this.getVariable("M") * adj_factor);
            this.setVariable("f_TRZ", this.getVariable("f_TRZ") * (2.0 - adj_factor));
        }
    }

    calibrateToObservations(times, g_obs) {
        if (times.length !== g_obs.length || times.length === 0) return;
        
        const obs = [];
        for (let i = 0; i < times.length; ++i) {
            obs.push([times[i], g_obs[i]]);
        }
        
        for (let iter = 0; iter < 5; ++iter) {
            this.autoRefineParameters(obs);
        }
    }

    optimizeForMetric(metric, t_start, t_end, steps) {
        let best_score = -1e100;
        const dt = (t_end - t_start) / steps;
        
        for (let i = 0; i <= steps; ++i) {
            const t = t_start + i * dt;
            const g = this.compute_g_NGC2525(t);
            const score = metric(g);
            if (score > best_score) best_score = score;
        }
        return best_score;
    }

    // --- Parameter Exploration (1 method) ---
    generateVariations(count, variation_pct) {
        const variations = [];
        const vars = this.listVariables();
        
        for (let i = 0; i < count; ++i) {
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
    mutateParameters(mutation_rate) {
        const vars = this.listVariables();
        
        for (const v of vars) {
            if (v === "c_light" || v === "G" || v === "hbar") continue;
            const val = this.getVariable(v);
            const delta = val * (Math.random() * 2 - 1) * mutation_rate;
            this.setVariable(v, val + delta);
        }
    }

    evolveSystem(generations, fitness) {
        let best_fitness = fitness(this);
        this.saveState("evolution_best");
        
        for (let gen = 0; gen < generations; ++gen) {
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
        GalaxyNGC2525.savedStates[stateName] = state;
        return true;
    }

    restoreState(stateName) {
        const state = GalaxyNGC2525.savedStates[stateName];
        if (!state) return false;
        
        for (const [key, value] of Object.entries(state)) {
            this.setVariable(key, value);
        }
        return true;
    }

    listSavedStates() {
        return Object.keys(GalaxyNGC2525.savedStates);
    }

    exportState() {
        let output = "GalaxyNGC2525 State Export:\n";
        const vars = this.listVariables();
        for (const v of vars) {
            output += `${v} = ${this.getVariable(v).toExponential(6)}\n`;
        }
        return output;
    }

    // --- System Analysis (4 methods) ---
    sensitivityAnalysis(t, delta_pct) {
        const sensitivities = {};
        const g_base = this.compute_g_NGC2525(t);
        
        const vars = this.listVariables();
        for (const v of vars) {
            if (v === "c_light" || v === "G" || v === "hbar") continue;
            
            const original = this.getVariable(v);
            const delta = original * delta_pct / 100.0;
            
            this.setVariable(v, original + delta);
            const g_plus = this.compute_g_NGC2525(t);
            this.setVariable(v, original);
            
            const sensitivity = (g_base !== 0.0) ? Math.abs((g_plus - g_base) / g_base) : 0.0;
            sensitivities[v] = sensitivity;
        }
        
        return sensitivities;
    }

    generateReport(t) {
        let output = "";
        output += "============================================\n";
        output += "GALAXY NGC 2525 REPORT\n";
        output += "Barred Spiral Galaxy with Supernova\n";
        output += "============================================\n";
        output += `Time: t = ${t.toFixed(6)} s (${(t/3.156e7).toFixed(6)} years)\n\n`;
        
        output += "Physical Parameters:\n";
        const M_sun = 1.989e30;
        output += `  Galaxy Mass M = ${this.M.toFixed(6)} kg (${(this.M/M_sun).toFixed(6)} M_sun)\n`;
        const ly_to_m = 9.461e15;
        output += `  Galaxy Radius r = ${this.r.toFixed(6)} m (${(this.r/ly_to_m).toFixed(6)} ly)\n`;
        output += `  Galaxy Redshift z_gal = ${this.z_gal.toFixed(6)}\n`;
        output += `  Hubble Parameter Hz = ${this.Hz.toFixed(6)} s^-1\n`;
        output += `  Magnetic field B = ${this.B.toFixed(6)} T (B_crit = ${this.B_crit.toFixed(6)} T)\n`;
        output += `  Black Hole M_BH = ${this.M_BH.toFixed(6)} kg (${(this.M_BH/M_sun).toFixed(6)} M_sun)\n`;
        output += `  Black Hole Radius r_BH = ${this.r_BH.toFixed(6)} m (${(this.r_BH/1.496e11).toFixed(6)} AU)\n`;
        output += `  f_TRZ = ${this.f_TRZ.toFixed(6)}\n`;
        output += `  Supernova M_SN0 = ${this.M_SN0.toFixed(6)} kg (${(this.M_SN0/M_sun).toFixed(6)} M_sun)\n`;
        output += `  SN Timescale tau_SN = ${this.tau_SN.toFixed(6)} s (${(this.tau_SN/3.156e7).toFixed(6)} years)\n`;
        output += `  M_SN(t) = ${this.M_SN_t(t).toFixed(6)} kg (${(this.M_SN_t(t)/M_sun).toFixed(6)} M_sun)\n`;
        output += `  Fluid density rho_fluid = ${this.rho_fluid.toFixed(6)} kg/m^3\n`;
        output += `  Gas velocity = ${this.gas_v.toFixed(6)} m/s\n`;
        output += `  DM factor = ${this.M_DM_factor.toFixed(6)}\n\n`;
        
        output += "Computed Acceleration:\n";
        output += `  g_NGC2525(t) = ${this.compute_g_NGC2525(t).toFixed(6)} m/s^2\n\n`;
        
        output += "UQFF Terms:\n";
        const corr_H = 1 + this.Hz * t;
        const corr_B = 1 - this.B / this.B_crit;
        output += `  Base (with Hz, B): ${(this.ug1_base * corr_H * corr_B).toFixed(6)} m/s^2\n`;
        output += `  Black Hole: ${this.g_BH.toFixed(6)} m/s^2\n`;
        output += `  Ug total: ${this.compute_Ug().toFixed(6)} m/s^2\n`;
        output += `  Lambda: ${((this.Lambda * this.c_light * this.c_light) / 3.0).toFixed(6)} m/s^2\n`;
        
        const cross_vB = this.gas_v * this.B;
        const em_base = (this.q_charge * cross_vB) / this.proton_mass;
        const corr_UA = 1 + (this.rho_vac_UA / this.rho_vac_SCm);
        output += `  EM (scaled with UA): ${(em_base * corr_UA * this.scale_EM).toFixed(6)} m/s^2\n`;
        
        const sqrt_unc = Math.sqrt(this.delta_x * this.delta_p);
        output += `  Quantum: ${((this.hbar / sqrt_unc) * this.integral_psi * (2 * Math.PI / this.t_Hubble)).toFixed(6)} m/s^2\n`;
        
        const V = this.compute_V();
        output += `  Fluid: ${((this.rho_fluid * V * this.ug1_base) / this.M).toFixed(6)} m/s^2\n`;
        
        output += `  Oscillatory: (combined real parts)\n`;
        
        const M_dm = this.M * this.M_DM_factor;
        const pert1 = this.delta_rho_over_rho;
        const pert2 = 3 * this.G * this.M / (this.r * this.r * this.r);
        const term_dm_force_like = (this.M + M_dm) * (pert1 + pert2);
        output += `  DM: ${(term_dm_force_like / this.M).toFixed(6)} m/s^2\n`;
        
        const MSNt = this.M_SN_t(t);
        output += `  Supernova Mass Loss: ${(-(this.G * MSNt) / (this.r * this.r)).toFixed(6)} m/s^2\n`;
        
        output += "============================================\n";
        return output;
    }

    validateConsistency() {
        let valid = true;
        
        if (this.M <= 0 || this.r <= 0) { console.error("Error: M and r must be positive."); valid = false; }
        if (this.B < 0 || this.B_crit <= 0) { console.error("Error: B, B_crit must be non-negative/positive."); valid = false; }
        if (this.Hz < 0) { console.error("Error: Hz must be non-negative."); valid = false; }
        if (this.z_gal < 0) { console.warn("Warning: Galaxy redshift z_gal is negative."); }
        if (this.M_BH <= 0 || this.r_BH <= 0) { console.error("Error: M_BH and r_BH must be positive."); valid = false; }
        if (this.M_SN0 < 0 || this.tau_SN <= 0) { console.error("Error: M_SN0 must be non-negative, tau_SN positive."); valid = false; }
        if (this.rho_fluid <= 0) { console.error("Error: Fluid density must be positive."); valid = false; }
        if (this.gas_v < 0) { console.error("Error: Gas velocity must be non-negative."); valid = false; }
        if (this.M_DM_factor < 0 || this.M_DM_factor > 1.0) { console.warn("Warning: DM factor outside [0,1]."); }
        
        return valid;
    }

    autoCorrectAnomalies() {
        let corrected = false;
        
        const M_sun = 1.989e30;
        
        if (this.M <= 0) { this.M = (1e10 + 2.25e7) * M_sun; corrected = true; }
        if (this.r <= 0) { this.r = 2.836e20; corrected = true; }
        if (this.B < 0) { this.B = 1e-5; corrected = true; }
        if (this.B_crit <= 0) { this.B_crit = 1e11; corrected = true; }
        if (this.Hz < 0) { this.Hz = 2.19e-18; corrected = true; }
        if (this.z_gal < 0) { this.z_gal = 0.016; corrected = true; }
        if (this.M_BH <= 0) { this.M_BH = 2.25e7 * M_sun; corrected = true; }
        if (this.r_BH <= 0) { this.r_BH = 1.496e11; corrected = true; }
        if (this.M_SN0 < 0) { this.M_SN0 = 1.4 * M_sun; corrected = true; }
        if (this.tau_SN <= 0) { this.tau_SN = 1 * 3.156e7; corrected = true; }
        if (this.rho_fluid <= 0) { this.rho_fluid = 1e-21; corrected = true; }
        if (this.gas_v < 0) { this.gas_v = 1e5; corrected = true; }
        if (this.M_DM_factor < 0) { this.M_DM_factor = 0.1; corrected = true; }
        if (this.M_DM_factor > 1.0) { this.M_DM_factor = 1.0; corrected = true; }
        
        if (corrected) this.updateCache();
        return corrected;
    }
}

// ========== ENHANCED EXAMPLE FUNCTION ==========
function enhanced_ngc2525_example() {
    console.log("=========================================================");
    console.log("ENHANCED GALAXY NGC 2525 DEMONSTRATION");
    console.log("Barred Spiral Galaxy with Supernova Mass Loss");
    console.log("=========================================================\n");
    
    const ngc2525 = new GalaxyNGC2525();
    
    // Step 1: Initial state and validation
    console.log("Step 1: Initial State and Validation");
    console.log(`System: ${ngc2525.getSystemName()}`);
    console.log(`Validation: ${ngc2525.validateConsistency() ? "PASS" : "FAIL"}`);
    console.log(`Auto-corrected: ${ngc2525.autoCorrectAnomalies() ? "Yes" : "No"}\n`);
    
    // Step 2: Time evolution showing supernova mass loss M_SN(t)
    console.log("Step 2: Time Evolution (Supernova Mass Loss)");
    const t_years_array = [0.0, 1.0, 5.0, 10.0, 20.0];
    const M_sun = 1.989e30;
    for (const t_years of t_years_array) {
        const t = t_years * 3.156e7;
        const MSNt = ngc2525.M_SN_t(t);
        const g = ngc2525.compute_g_NGC2525(t);
        console.log(`  t = ${t_years.toFixed(1)} years: M_SN(t) = ${(MSNt/M_sun).toExponential(6)} M_sun, g = ${g.toExponential(6)} m/s^2`);
    }
    console.log();
    
    // Step 3: Variable listing
    console.log("Step 3: Variable Listing");
    const vars = ngc2525.listVariables();
    console.log(`Total variables: ${vars.length}`);
    console.log(`Sample: ${vars[0]}, ${vars[1]}, ${vars[11]} (M_BH), ${vars[13]} (M_SN0)\n`);
    
    // Step 4: Galactic mass scaling
    console.log("Step 4: Galactic Mass Scaling (M sweeps)");
    ngc2525.saveState("original");
    const M_factors = [0.5, 1.0, 2.0];
    for (const factor of M_factors) {
        ngc2525.restoreState("original");
        ngc2525.expandGalacticScale(factor, 1.0);
        const t = 7 * 3.156e7;
        const g = ngc2525.compute_g_NGC2525(t);
        const M = ngc2525.getVariable("M");
        console.log(`  M × ${factor}: M = ${(M/M_sun).toExponential(3)} M_sun, g(7 years) = ${g.toExponential(6)} m/s^2`);
    }
    ngc2525.restoreState("original");
    console.log();
    
    // Step 5: Black hole scaling (UNIQUE to galaxy with central BH)
    console.log("Step 5: Black Hole Scaling (M_BH sweeps) - CENTRAL BH FEATURE");
    const M_BH_factors = [0.5, 1.0, 2.0];
    for (const factor of M_BH_factors) {
        ngc2525.restoreState("original");
        ngc2525.expandBlackHoleScale(factor, 1.0);
        const t = 7 * 3.156e7;
        const g = ngc2525.compute_g_NGC2525(t);
        const M_BH = ngc2525.getVariable("M_BH");
        console.log(`  M_BH × ${factor}: M_BH = ${(M_BH/M_sun).toExponential(3)} M_sun, g(7 years) = ${g.toExponential(6)} m/s^2`);
    }
    ngc2525.restoreState("original");
    console.log();
    
    // Step 6: Supernova mass loss scaling (UNIQUE to SN galaxy)
    console.log("Step 6: Supernova Initial Mass Scaling (M_SN0 sweeps) - SN FEATURE");
    const M_SN0_factors = [0.5, 1.0, 2.0];
    for (const factor of M_SN0_factors) {
        ngc2525.restoreState("original");
        ngc2525.expandSupernovaScale(factor, 1.0);
        const t = 7 * 3.156e7;
        const MSNt = ngc2525.M_SN_t(t);
        const g = ngc2525.compute_g_NGC2525(t);
        console.log(`  M_SN0 × ${factor}: M_SN(7yr) = ${(MSNt/M_sun).toExponential(3)} M_sun, g = ${g.toExponential(6)} m/s^2`);
    }
    ngc2525.restoreState("original");
    console.log();
    
    // Step 7: Supernova decay timescale sweeps
    console.log("Step 7: Supernova Timescale Sweeps (tau_SN)");
    const tau_SN_factors = [0.5, 1.0, 2.0];
    for (const factor of tau_SN_factors) {
        ngc2525.restoreState("original");
        ngc2525.expandSupernovaScale(1.0, factor);
        const t = 7 * 3.156e7;
        const MSNt = ngc2525.M_SN_t(t);
        console.log(`  tau_SN × ${factor}: M_SN(7yr) = ${(MSNt/M_sun).toExponential(3)} M_sun`);
    }
    ngc2525.restoreState("original");
    console.log();
    
    // Step 8: Black hole radius scaling
    console.log("Step 8: Black Hole Influence Radius Scaling");
    const r_BH_factors = [0.5, 1.0, 2.0];
    for (const factor of r_BH_factors) {
        ngc2525.restoreState("original");
        ngc2525.expandBlackHoleScale(1.0, factor);
        const t = 7 * 3.156e7;
        const g = ngc2525.compute_g_NGC2525(t);
        console.log(`  r_BH × ${factor}: g(7 years) = ${g.toExponential(6)} m/s^2`);
    }
    ngc2525.restoreState("original");
    console.log();
    
    // Step 9: Parameter space expansion
    console.log("Step 9: Parameter Space Expansion (all scalable params)");
    ngc2525.expandParameterSpace(1.2);
    const M_after = ngc2525.getVariable("M");
    console.log(`  After 1.2× expansion: M = ${(M_after/M_sun).toExponential(3)} M_sun`);
    ngc2525.restoreState("original");
    console.log();
    
    // Step 10: Batch operations
    console.log("Step 10: Batch Operations (scale multiple variables)");
    const scale_group = ["M", "r", "M_BH"];
    ngc2525.scaleVariableGroup(scale_group, 1.1);
    console.log("  Scaled {M, r, M_BH} by 1.1×");
    ngc2525.restoreState("original");
    console.log();
    
    // Step 11: State management
    console.log("Step 11: State Management");
    ngc2525.saveState("state_A");
    ngc2525.expandGalacticScale(1.5, 1.2);
    ngc2525.saveState("state_B");
    const states = ngc2525.listSavedStates();
    console.log(`  Saved states: ${states.join(" ")}`);
    ngc2525.restoreState("state_A");
    console.log("  Restored state_A\n");
    
    // Step 12: Generate parameter variations
    console.log("Step 12: Generate Parameter Variations (5% variation)");
    const variations = ngc2525.generateVariations(3, 5.0);
    console.log(`  Generated ${variations.length} variants with 5% random variation`);
    console.log(`  Variant 1 M = ${variations[0].M.toExponential(6)} kg\n`);
    
    // Step 13: Sensitivity analysis
    console.log("Step 13: Sensitivity Analysis at 7 years");
    const t_sens = 7 * 3.156e7;
    const sensitivities = ngc2525.sensitivityAnalysis(t_sens, 1.0);
    console.log("  Top sensitivities (1% perturbation):");
    const sens_vec = Object.entries(sensitivities).sort((a, b) => b[1] - a[1]);
    for (let i = 0; i < 5 && i < sens_vec.length; ++i) {
        console.log(`    ${sens_vec[i][0]}: ${sens_vec[i][1].toExponential(3)}`);
    }
    console.log();
    
    // Step 14: Auto-refinement with synthetic observations
    console.log("Step 14: Auto-Refinement (synthetic observations)");
    const obs = [];
    for (let i = 0; i <= 5; ++i) {
        const t_obs = i * 5 * 3.156e7;
        const g_obs = ngc2525.compute_g_NGC2525(t_obs) * (1.0 + 0.01 * (Math.random() * 100 - 50) / 100.0);
        obs.push([t_obs, g_obs]);
    }
    ngc2525.autoRefineParameters(obs);
    console.log(`  Refined parameters based on ${obs.length} observations\n`);
    
    // Step 15: Optimization for maximum acceleration
    console.log("Step 15: Optimize for Maximum Acceleration");
    ngc2525.restoreState("original");
    const metric = (g) => g;
    const t_opt_start = 0.0;
    const t_opt_end = 20 * 3.156e7;
    const best_g = ngc2525.optimizeForMetric(metric, t_opt_start, t_opt_end, 50);
    console.log(`  Best g over 20 years: ${best_g.toExponential(6)} m/s^2\n`);
    
    // Step 16: Evolutionary system adaptation
    console.log("Step 16: Evolutionary System Adaptation (5 generations)");
    ngc2525.restoreState("original");
    const fitness = (n) => {
        const t = 7 * 3.156e7;
        return n.compute_g_NGC2525(t);
    };
    ngc2525.evolveSystem(5, fitness);
    console.log("  Evolved system over 5 generations (fitness = g at 7 years)\n");
    
    // Step 17: Full system report
    console.log("Step 17: Full System Report at 7 years");
    ngc2525.restoreState("original");
    const t_report = 7 * 3.156e7;
    const report = ngc2525.generateReport(t_report);
    console.log(report + "\n");
    
    // Step 18: Full state export
    console.log("Step 18: Full State Export");
    const exported = ngc2525.exportState();
    console.log("Exported state (first 500 chars):");
    console.log(exported.substring(0, 500) + "...\n");
    
    console.log("=========================================================");
    console.log("ENHANCED DEMONSTRATION COMPLETE");
    console.log("=========================================================");
}

// Module exports
if (typeof module !== 'undefined' && module.exports) {
    module.exports = GalaxyNGC2525;
    module.exports.enhanced_ngc2525_example = enhanced_ngc2525_example;
} else if (typeof window !== 'undefined') {
    window.GalaxyNGC2525 = GalaxyNGC2525;
    window.enhanced_ngc2525_example = enhanced_ngc2525_example;
}

// Default export for ES6 modules
export default GalaxyNGC2525;
export { enhanced_ngc2525_example };
