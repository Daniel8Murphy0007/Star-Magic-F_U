/**
 * ================================================================================================
 * Module: MagnetarSGR0501_4516.js
 *
 * Description: JavaScript Module for SGR 0501+4516 Magnetar Class
 *              This is the first module in a series of 500+ code files for the Universal Quantum
 *              Field Framework (UQFF) simulations, focusing on magnetar evolution and gravity
 *              equations derived from Hubble datasets, high-energy lab simulations, and UQFF
 *              refinements (dated May 08, 2025, updated for full term inclusion on October 08, 2025).
 *
 * Purpose: Encapsulates the Master Universal Gravity Equation (MUGE) for SGR 0501+4516 magnetar
 *          evolution. Now includes ALL terms (no omissions): base gravity, cosmic expansion,
 *          magnetic decay, UQFF Ug components with f_TRZ, Lambda, scaled EM, GW, quantum uncertainty,
 *          fluid dynamics (effective acceleration; note: may overlap with base self-gravity),
 *          oscillatory waves (real part; assumes consistent units), and DM/density perturbations
 *          (converted to acceleration via /M). Supports dynamic variable updates for all parameters.
 *
 * Integration: Designed for Node.js and browser environments.
 *              Node.js: const MagnetarSGR0501_4516 = require('./source14.js');
 *              ES6: import MagnetarSGR0501_4516 from './source14.js';
 *              Instantiate: const mag = new MagnetarSGR0501_4516();
 *              Compute: const g = mag.compute_g_Magnetar(t);
 *
 * Key Features:
 *   - Default values from UQFF document, with approximations for previously "negligible" terms.
 *   - Setter methods for updates: setVariable(name, value) or addToVariable(name, delta)/subtractFromVariable(name, delta).
 *   - Computes g_Magnetar(r, t) with every term explicitly included.
 *   - 25 enhanced dynamic capabilities for parameter exploration, self-expansion, optimization, etc.
 *
 * Author: Converted to JavaScript by GitHub Copilot, based on C++ code by Grok (xAI)
 *         Original UQFF manuscript by Daniel T. Murphy
 * Date: November 03, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

class MagnetarSGR0501_4516 {
    constructor() {
        // Core parameters (mutable for updates)
        this.G = 0;               // Gravitational constant
        this.M = 0;               // Magnetar mass
        this.r = 0;               // Radius
        this.H0 = 0;              // Hubble constant (s^-1)
        this.B0 = 0;              // Initial magnetic field
        this.tau_B = 0;           // B decay timescale (s)
        this.B_crit = 0;          // Critical B field
        this.Lambda = 0;          // Cosmological constant
        this.c_light = 0;         // Speed of light
        this.q_charge = 0;        // Charge (proton)
        this.v_surf = 0;          // Surface velocity
        this.f_TRZ = 0;           // Time-reversal factor
        this.rho_vac_UA = 0;      // UA vacuum density
        this.rho_vac_SCm = 0;     // SCm vacuum density
        this.P_init = 0;          // Initial rotation period (s)
        this.tau_Omega = 0;       // Omega decay timescale (s)
        this.scale_EM = 0;        // EM scaling factor
        this.proton_mass = 0;     // Proton mass for EM acceleration

        // Additional parameters for full inclusion of terms
        this.hbar = 0;            // Reduced Planck's constant
        this.t_Hubble = 0;        // Hubble time (s)
        this.delta_x = 0;         // Position uncertainty (m)
        this.delta_p = 0;         // Momentum uncertainty (kg m/s)
        this.integral_psi = 0;    // Wavefunction integral approximation
        this.rho_fluid = 0;       // Fluid density (kg/m^3)
        this.A_osc = 0;           // Oscillatory amplitude (assumed m/s^2)
        this.k_osc = 0;           // Wave number (1/m)
        this.omega_osc = 0;       // Angular frequency (rad/s)
        this.x_pos = 0;           // Position for oscillation (m)
        this.t_Hubble_gyr = 0;    // Hubble time in Gyr (for oscillatory prefactor)
        this.M_DM_factor = 0;     // Dark matter mass fraction
        this.delta_rho_over_rho = 0; // Density perturbation fraction

        // Computed caches (updated on demand)
        this.ug1_base = 0;        // Cached Ug1 = G*M/r^2

        // Initialize with default UQFF values
        this.initializeDefaults();
    }

    // Initialization method (called in constructor)
    initializeDefaults() {
        this.G = 6.6743e-11;
        this.M = 1.4 * 1.989e30;
        this.r = 20e3;
        this.H0 = 2.184e-18;
        this.B0 = 1e10;
        this.tau_B = 4000 * 3.156e7;
        this.B_crit = 1e11;
        this.Lambda = 1.1e-52;
        this.c_light = 3e8;
        this.q_charge = 1.602e-19;
        this.v_surf = 1e6;
        this.f_TRZ = 0.1;
        this.rho_vac_UA = 7.09e-36;
        this.rho_vac_SCm = 7.09e-37;
        this.P_init = 5.0;
        this.tau_Omega = 10000 * 3.156e7;
        this.scale_EM = 1e-12;
        this.proton_mass = 1.673e-27;

        // Full terms defaults
        this.hbar = 1.0546e-34;
        this.t_Hubble = 13.8e9 * 3.156e7;
        this.delta_x = 1e-10;  // Arbitrary for uncertainty principle
        this.delta_p = this.hbar / this.delta_x;
        this.integral_psi = 1.0;
        this.rho_fluid = 1e17;
        this.A_osc = 1e10;     // Arbitrary amplitude to scale ~1e10 m/s^2
        this.k_osc = 1.0 / this.r;
        this.omega_osc = 2 * Math.PI / this.P_init;
        this.x_pos = this.r;
        this.t_Hubble_gyr = 13.8;
        this.M_DM_factor = 0.1;
        this.delta_rho_over_rho = 1e-5;

        this.updateCache();
    }

    // Cache update for efficiency (call after parameter changes)
    updateCache() {
        this.ug1_base = (this.G * this.M) / (this.r * this.r);
        // Update dependent params if needed, e.g., delta_p = hbar / delta_x; but kept independent
    }

    // Universal setter for any variable (by name, for flexibility)
    setVariable(varName, newValue) {
        const validVars = {
            'G': () => { this.G = newValue; },
            'M': () => { this.M = newValue; },
            'r': () => { this.r = newValue; },
            'H0': () => { this.H0 = newValue; },
            'B0': () => { this.B0 = newValue; },
            'tau_B': () => { this.tau_B = newValue; },
            'B_crit': () => { this.B_crit = newValue; },
            'Lambda': () => { this.Lambda = newValue; },
            'c_light': () => { this.c_light = newValue; },
            'q_charge': () => { this.q_charge = newValue; },
            'v_surf': () => { this.v_surf = newValue; },
            'f_TRZ': () => { this.f_TRZ = newValue; },
            'rho_vac_UA': () => { this.rho_vac_UA = newValue; },
            'rho_vac_SCm': () => { this.rho_vac_SCm = newValue; },
            'P_init': () => { this.P_init = newValue; },
            'tau_Omega': () => { this.tau_Omega = newValue; },
            'scale_EM': () => { this.scale_EM = newValue; },
            'proton_mass': () => { this.proton_mass = newValue; },
            'hbar': () => { this.hbar = newValue; },
            't_Hubble': () => { this.t_Hubble = newValue; },
            'delta_x': () => { this.delta_x = newValue; },
            'delta_p': () => { this.delta_p = newValue; },
            'integral_psi': () => { this.integral_psi = newValue; },
            'rho_fluid': () => { this.rho_fluid = newValue; },
            'A_osc': () => { this.A_osc = newValue; },
            'k_osc': () => { this.k_osc = newValue; },
            'omega_osc': () => { this.omega_osc = newValue; },
            'x_pos': () => { this.x_pos = newValue; },
            't_Hubble_gyr': () => { this.t_Hubble_gyr = newValue; },
            'M_DM_factor': () => { this.M_DM_factor = newValue; },
            'delta_rho_over_rho': () => { this.delta_rho_over_rho = newValue; }
        };

        if (validVars[varName]) {
            validVars[varName]();
            this.updateCache();
            return true;
        } else {
            console.error(`Error: Unknown variable '${varName}'.`);
            return false;
        }
    }

    // Addition method for variables
    addToVariable(varName, delta) {
        const currentValue = this.getVariable(varName);
        if (currentValue === null) {
            return false;
        }
        return this.setVariable(varName, currentValue + delta);
    }

    // Subtraction method for variables
    subtractFromVariable(varName, delta) {
        return this.addToVariable(varName, -delta);
    }

    // Getter for any variable (helper for add/subtract)
    getVariable(varName) {
        const validVars = {
            'G': () => this.G,
            'M': () => this.M,
            'r': () => this.r,
            'H0': () => this.H0,
            'B0': () => this.B0,
            'tau_B': () => this.tau_B,
            'B_crit': () => this.B_crit,
            'Lambda': () => this.Lambda,
            'c_light': () => this.c_light,
            'q_charge': () => this.q_charge,
            'v_surf': () => this.v_surf,
            'f_TRZ': () => this.f_TRZ,
            'rho_vac_UA': () => this.rho_vac_UA,
            'rho_vac_SCm': () => this.rho_vac_SCm,
            'P_init': () => this.P_init,
            'tau_Omega': () => this.tau_Omega,
            'scale_EM': () => this.scale_EM,
            'proton_mass': () => this.proton_mass,
            'hbar': () => this.hbar,
            't_Hubble': () => this.t_Hubble,
            'delta_x': () => this.delta_x,
            'delta_p': () => this.delta_p,
            'integral_psi': () => this.integral_psi,
            'rho_fluid': () => this.rho_fluid,
            'A_osc': () => this.A_osc,
            'k_osc': () => this.k_osc,
            'omega_osc': () => this.omega_osc,
            'x_pos': () => this.x_pos,
            't_Hubble_gyr': () => this.t_Hubble_gyr,
            'M_DM_factor': () => this.M_DM_factor,
            'delta_rho_over_rho': () => this.delta_rho_over_rho
        };

        if (validVars[varName]) {
            return validVars[varName]();
        } else {
            console.error(`Error: Unknown variable '${varName}'.`);
            return null;
        }
    }

    // B(t) computation
    B_t(t) {
        return this.B0 * Math.exp(-t / this.tau_B);
    }

    // Omega(t) computation
    Omega_t(t) {
        return (2 * Math.PI / this.P_init) * Math.exp(-t / this.tau_Omega);
    }

    // dOmega/dt computation
    dOmega_dt(t) {
        const omega0 = 2 * Math.PI / this.P_init;
        return omega0 * (-1.0 / this.tau_Omega) * Math.exp(-t / this.tau_Omega);
    }

    // Ug terms computation
    compute_Ug(Bt) {
        const Ug1 = this.ug1_base;
        const Ug2 = 0.0;
        const Ug3 = 0.0;
        const Ug4 = Ug1 * (1 - Bt / this.B_crit);
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1 + this.f_TRZ);
    }

    // Volume computation for fluid
    compute_V() {
        return (4.0 / 3.0) * Math.PI * this.r * this.r * this.r;
    }

    // Main MUGE computation (now includes ALL terms)
    compute_g_Magnetar(t) {
        if (t < 0) {
            console.error('Error: Time t must be non-negative.');
            return 0.0;
        }

        const Bt = this.B_t(t);
        const dOdt = this.dOmega_dt(t);

        // Term 1: Base + H0 + B corrections
        const corr_H = 1 + this.H0 * t;
        const corr_B = 1 - Bt / this.B_crit;
        const term1 = this.ug1_base * corr_H * corr_B;

        // Term 2: UQFF Ug with f_TRZ
        const term2 = this.compute_Ug(Bt);

        // Term 3: Lambda
        const term3 = (this.Lambda * this.c_light * this.c_light) / 3.0;

        // Term 4: Scaled EM
        const cross_vB = this.v_surf * Bt;  // Magnitude, assuming perpendicular
        const em_base = (this.q_charge * cross_vB) / this.proton_mass;
        const corr_UA = 1 + (this.rho_vac_UA / this.rho_vac_SCm);
        const term4 = (em_base * corr_UA) * this.scale_EM;

        // Term 5: GW
        const gw_prefactor = (this.G * this.M * this.M) / (Math.pow(this.c_light, 4) * this.r);
        const term5 = gw_prefactor * (dOdt * dOdt);

        // Quantum uncertainty term
        const sqrt_unc = Math.sqrt(this.delta_x * this.delta_p);
        const term_q = (this.hbar / sqrt_unc) * this.integral_psi * (2 * Math.PI / this.t_Hubble);

        // Fluid term (effective acceleration; note: may overlap with base self-gravity)
        const V = this.compute_V();
        const term_fluid = (this.rho_fluid * V * this.ug1_base) / this.M;

        // Oscillatory terms (real parts)
        const term_osc1 = 2 * this.A_osc * Math.cos(this.k_osc * this.x_pos) * Math.cos(this.omega_osc * t);
        const term_osc2 = (2 * Math.PI / this.t_Hubble) * this.A_osc * Math.cos(this.k_osc * this.x_pos - this.omega_osc * t);
        const term_osc = term_osc1 + term_osc2;

        // DM and density perturbation term (converted to acceleration)
        const M_dm = this.M * this.M_DM_factor;
        const pert1 = this.delta_rho_over_rho;
        const pert2 = 3 * this.G * this.M / (this.r * this.r * this.r);
        const term_dm_force_like = (this.M + M_dm) * (pert1 + pert2);
        const term_DM = term_dm_force_like / this.M;

        // Total g_Magnetar (all terms summed)
        return term1 + term2 + term3 + term4 + term5 + term_q + term_fluid + term_osc + term_DM;
    }

    // Debug/Output method (for transparency in base program)
    printParameters() {
        console.log('SGR 0501+4516 Parameters:');
        console.log(`G: ${this.G.toExponential(3)}, M: ${this.M.toExponential(3)}, r: ${this.r.toExponential(3)}`);
        console.log(`H0: ${this.H0.toExponential(3)}, B0: ${this.B0.toExponential(3)}, tau_B: ${this.tau_B.toExponential(3)}`);
        console.log(`f_TRZ: ${this.f_TRZ.toFixed(3)}, rho_fluid: ${this.rho_fluid.toExponential(3)}, M_DM_factor: ${this.M_DM_factor.toFixed(3)}`);
        console.log(`A_osc: ${this.A_osc.toExponential(3)}, delta_rho_over_rho: ${this.delta_rho_over_rho.toExponential(3)}`);
        console.log(`ug1_base: ${this.ug1_base.toExponential(3)}`);
    }

    // Example computation at t=5000 years (for testing; now includes all terms)
    exampleAt5000Years() {
        const t_example = 5000 * 3.156e7;
        return this.compute_g_Magnetar(t_example);
    }

    // ========== ENHANCED DYNAMIC CAPABILITIES (25 methods) ==========

    // --- Variable Management (5 methods) ---
    createVariable(name, value) {
        return this.setVariable(name, value);
    }

    removeVariable(name) {
        console.warn(`Warning: Cannot remove built-in variable '${name}' in magnetar class.`);
        return false;
    }

    cloneVariable(src, dest) {
        const val = this.getVariable(src);
        if (val === null) return false;
        return this.setVariable(dest, val);
    }

    listVariables() {
        return ['G', 'M', 'r', 'H0', 'B0', 'tau_B', 'B_crit', 'Lambda', 'c_light', 'q_charge',
                'v_surf', 'f_TRZ', 'rho_vac_UA', 'rho_vac_SCm', 'P_init', 'tau_Omega', 'scale_EM', 'proton_mass',
                'hbar', 't_Hubble', 'delta_x', 'delta_p', 'integral_psi', 'rho_fluid', 'A_osc', 'k_osc',
                'omega_osc', 'x_pos', 't_Hubble_gyr', 'M_DM_factor', 'delta_rho_over_rho'];
    }

    getSystemName() {
        return 'MagnetarSGR0501_4516';
    }

    // --- Batch Operations (2 methods) ---
    transformVariableGroup(names, func) {
        for (const name of names) {
            const val = this.getVariable(name);
            if (val === null) return false;
            if (!this.setVariable(name, func(val))) return false;
        }
        return true;
    }

    scaleVariableGroup(names, factor) {
        return this.transformVariableGroup(names, v => v * factor);
    }

    // --- Self-Expansion (4 methods) ---
    expandParameterSpace(factor) {
        const expandable = ['M', 'r', 'B0', 'tau_B', 'rho_fluid', 'A_osc', 'M_DM_factor'];
        return this.scaleVariableGroup(expandable, factor);
    }

    expandMagneticScale(B_factor, tau_factor) {
        this.setVariable('B0', this.getVariable('B0') * B_factor);
        this.setVariable('B_crit', this.getVariable('B_crit') * B_factor);
        this.setVariable('tau_B', this.getVariable('tau_B') * tau_factor);
    }

    expandDecayScale(omega_factor, tau_omega_factor) {
        this.setVariable('P_init', this.getVariable('P_init') / omega_factor);  // Period inversely related
        this.setVariable('tau_Omega', this.getVariable('tau_Omega') * tau_omega_factor);
    }

    expandFluidDMScale(rho_factor, DM_factor) {
        this.setVariable('rho_fluid', this.getVariable('rho_fluid') * rho_factor);
        this.setVariable('M_DM_factor', this.getVariable('M_DM_factor') * DM_factor);
        this.setVariable('delta_rho_over_rho', this.getVariable('delta_rho_over_rho') * DM_factor);
    }

    // --- Self-Refinement (3 methods) ---
    autoRefineParameters(observations) {
        if (observations.length === 0) return;
        
        let sum_error = 0.0;
        for (const obs of observations) {
            const t = obs[0];
            const g_obs = obs[1];
            const g_calc = this.compute_g_Magnetar(t);
            sum_error += Math.abs(g_calc - g_obs);
        }
        const avg_error = sum_error / observations.length;
        
        // Simple refinement: adjust B0 and tau_B based on error
        if (avg_error > 1e-6) {
            const adj_factor = 1.0 - Math.min(0.1, avg_error / 1e6);
            this.setVariable('B0', this.getVariable('B0') * adj_factor);
            this.setVariable('tau_B', this.getVariable('tau_B') * (2.0 - adj_factor));
        }
    }

    calibrateToObservations(times, g_obs) {
        if (times.length !== g_obs.length || times.length === 0) return;
        
        const obs = [];
        for (let i = 0; i < times.length; i++) {
            obs.push([times[i], g_obs[i]]);
        }
        
        // Iterative refinement (5 passes)
        for (let iter = 0; iter < 5; iter++) {
            this.autoRefineParameters(obs);
        }
    }

    optimizeForMetric(metric, t_start, t_end, steps) {
        let best_score = -1e100;
        const dt = (t_end - t_start) / steps;
        
        for (let i = 0; i <= steps; i++) {
            const t = t_start + i * dt;
            const g = this.compute_g_Magnetar(t);
            const score = metric(g);
            if (score > best_score) best_score = score;
        }
        return best_score;
    }

    // --- Parameter Exploration (1 method) ---
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
    mutateParameters(mutation_rate, rng = null) {
        const vars = this.listVariables();
        
        for (const v of vars) {
            if (v === 'c_light' || v === 'G' || v === 'hbar') continue; // Skip constants
            const val = this.getVariable(v);
            const delta = val * ((Math.random() * 2 - 1) * mutation_rate);
            this.setVariable(v, val + delta);
        }
    }

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
    saveState(stateName) {
        const state = {};
        const vars = this.listVariables();
        for (const v of vars) {
            state[v] = this.getVariable(v);
        }
        MagnetarSGR0501_4516.savedStates[stateName] = state;
        return true;
    }

    restoreState(stateName) {
        const state = MagnetarSGR0501_4516.savedStates[stateName];
        if (!state) return false;
        
        for (const [key, value] of Object.entries(state)) {
            this.setVariable(key, value);
        }
        return true;
    }

    listSavedStates() {
        return Object.keys(MagnetarSGR0501_4516.savedStates);
    }

    exportState() {
        let output = 'MagnetarSGR0501_4516 State Export:\n';
        const vars = this.listVariables();
        for (const v of vars) {
            output += `${v} = ${this.getVariable(v).toExponential(6)}\n`;
        }
        return output;
    }

    // --- System Analysis (4 methods) ---
    sensitivityAnalysis(t, delta_pct) {
        const sensitivities = {};
        const g_base = this.compute_g_Magnetar(t);
        
        const vars = this.listVariables();
        for (const v of vars) {
            if (v === 'c_light' || v === 'G' || v === 'hbar') continue; // Skip constants
            
            const original = this.getVariable(v);
            const delta = original * delta_pct / 100.0;
            
            this.setVariable(v, original + delta);
            const g_plus = this.compute_g_Magnetar(t);
            this.setVariable(v, original);
            
            const sensitivity = (g_base !== 0.0) ? Math.abs((g_plus - g_base) / g_base) : 0.0;
            sensitivities[v] = sensitivity;
        }
        
        return sensitivities;
    }

    generateReport(t) {
        let report = '';
        report += '============================================\n';
        report += 'SGR 0501+4516 MAGNETAR SYSTEM REPORT\n';
        report += '============================================\n';
        report += `Time: t = ${t.toExponential(6)} s (${(t/3.156e7).toFixed(2)} years)\n\n`;
        
        report += 'Physical Parameters:\n';
        report += `  Mass M = ${this.M.toExponential(6)} kg (${(this.M/1.989e30).toFixed(3)} M_sun)\n`;
        report += `  Radius r = ${this.r.toExponential(6)} m (${(this.r/1e3).toFixed(3)} km)\n`;
        report += `  B-field B0 = ${this.B0.toExponential(6)} T (B_crit = ${this.B_crit.toExponential(6)} T)\n`;
        report += `  B(t) = ${this.B_t(t).toExponential(6)} T\n`;
        report += `  Period P_init = ${this.P_init.toFixed(3)} s\n`;
        report += `  Omega(t) = ${this.Omega_t(t).toExponential(6)} rad/s\n`;
        report += `  Fluid density rho = ${this.rho_fluid.toExponential(6)} kg/m^3\n`;
        report += `  DM factor = ${this.M_DM_factor.toFixed(3)}\n\n`;
        
        report += 'Computed Acceleration:\n';
        report += `  g_Magnetar(t) = ${this.compute_g_Magnetar(t).toExponential(6)} m/s^2\n\n`;
        
        report += 'UQFF Terms:\n';
        const Bt = this.B_t(t);
        const corr_H = 1 + this.H0 * t;
        const corr_B = 1 - Bt / this.B_crit;
        report += `  Base (with H0, B): ${(this.ug1_base * corr_H * corr_B).toExponential(6)} m/s^2\n`;
        report += `  Ug total: ${this.compute_Ug(Bt).toExponential(6)} m/s^2\n`;
        report += `  Lambda: ${((this.Lambda * this.c_light * this.c_light) / 3.0).toExponential(6)} m/s^2\n`;
        
        const cross_vB = this.v_surf * Bt;
        const em_base = (this.q_charge * cross_vB) / this.proton_mass;
        const corr_UA = 1 + (this.rho_vac_UA / this.rho_vac_SCm);
        report += `  EM (scaled): ${(em_base * corr_UA * this.scale_EM).toExponential(6)} m/s^2\n`;
        
        const dOdt = this.dOmega_dt(t);
        const gw_prefactor = (this.G * this.M * this.M) / (Math.pow(this.c_light, 4) * this.r);
        report += `  GW: ${(gw_prefactor * dOdt * dOdt).toExponential(6)} m/s^2\n`;
        
        const sqrt_unc = Math.sqrt(this.delta_x * this.delta_p);
        report += `  Quantum: ${((this.hbar / sqrt_unc) * this.integral_psi * (2 * Math.PI / this.t_Hubble)).toExponential(6)} m/s^2\n`;
        
        const V = this.compute_V();
        report += `  Fluid: ${((this.rho_fluid * V * this.ug1_base) / this.M).toExponential(6)} m/s^2\n`;
        
        report += '  Oscillatory: (combined real parts)\n';
        
        const M_dm = this.M * this.M_DM_factor;
        const pert1 = this.delta_rho_over_rho;
        const pert2 = 3 * this.G * this.M / (this.r * this.r * this.r);
        const term_dm_force_like = (this.M + M_dm) * (pert1 + pert2);
        report += `  DM: ${(term_dm_force_like / this.M).toExponential(6)} m/s^2\n`;
        
        report += '============================================\n';
        return report;
    }

    validateConsistency() {
        let valid = true;
        
        if (this.M <= 0 || this.r <= 0) { 
            console.error('Error: M and r must be positive.'); 
            valid = false; 
        }
        if (this.B0 < 0 || this.B_crit <= 0) { 
            console.error('Error: B0, B_crit must be non-negative/positive.'); 
            valid = false; 
        }
        if (this.tau_B <= 0 || this.tau_Omega <= 0) { 
            console.error('Error: Decay timescales must be positive.'); 
            valid = false; 
        }
        if (this.rho_fluid < 0) { 
            console.error('Error: Fluid density must be non-negative.'); 
            valid = false; 
        }
        if (this.M_DM_factor < 0 || this.M_DM_factor > 1.0) { 
            console.warn('Warning: DM factor outside [0,1].'); 
        }
        
        return valid;
    }

    autoCorrectAnomalies() {
        let corrected = false;
        
        if (this.M <= 0) { this.M = 1.4 * 1.989e30; corrected = true; }
        if (this.r <= 0) { this.r = 20e3; corrected = true; }
        if (this.B0 < 0) { this.B0 = 1e10; corrected = true; }
        if (this.B_crit <= 0) { this.B_crit = 1e11; corrected = true; }
        if (this.tau_B <= 0) { this.tau_B = 4000 * 3.156e7; corrected = true; }
        if (this.tau_Omega <= 0) { this.tau_Omega = 10000 * 3.156e7; corrected = true; }
        if (this.rho_fluid < 0) { this.rho_fluid = 1e17; corrected = true; }
        if (this.M_DM_factor < 0) { this.M_DM_factor = 0.1; corrected = true; }
        if (this.M_DM_factor > 1.0) { this.M_DM_factor = 1.0; corrected = true; }
        
        if (corrected) this.updateCache();
        return corrected;
    }
}

// Static property for state storage (class-level)
MagnetarSGR0501_4516.savedStates = {};

// ========== ENHANCED EXAMPLE FUNCTION ==========
function enhancedMagnetarSGR0501Example() {
    console.log('\n========== ENHANCED SGR 0501+4516 MAGNETAR UQFF EXAMPLE ==========\n');
    
    const mag = new MagnetarSGR0501_4516();
    
    // Step 1: Initial state
    console.log('Step 1: Initial Configuration');
    mag.printParameters();
    const t0 = 0.0;
    console.log(`g_Magnetar(t=0) = ${mag.compute_g_Magnetar(t0).toExponential(6)} m/s^2\n`);
    
    // Step 2: Time evolution
    console.log('Step 2: Time Evolution (0, 1000, 3000, 5000 years)');
    for (const t_yr of [0.0, 1000.0, 3000.0, 5000.0]) {
        const t = t_yr * 3.156e7;
        console.log(`  t = ${t_yr} yr: g = ${mag.compute_g_Magnetar(t).toExponential(6)} m/s^2, B(t) = ${mag.B_t(t).toExponential(6)} T, Omega(t) = ${mag.Omega_t(t).toExponential(6)} rad/s`);
    }
    console.log('');
    
    // Step 3: Magnetic field scaling
    console.log('Step 3: Magnetic Field Scaling (B0 x1.5, tau_B x0.8)');
    mag.expandMagneticScale(1.5, 0.8);
    const t_test = 2000 * 3.156e7;
    console.log(`After expansion: B0 = ${mag.getVariable('B0').toExponential(6)} T, tau_B = ${mag.getVariable('tau_B').toExponential(6)} s`);
    console.log(`g_Magnetar(t=2000 yr) = ${mag.compute_g_Magnetar(t_test).toExponential(6)} m/s^2\n`);
    
    // Step 4: Rotation decay scaling
    console.log('Step 4: Rotation Decay Scaling (Omega x1.2, tau_Omega x1.5)');
    mag.expandDecayScale(1.2, 1.5);
    console.log(`After expansion: P_init = ${mag.getVariable('P_init').toFixed(3)} s, tau_Omega = ${mag.getVariable('tau_Omega').toExponential(6)} s`);
    console.log(`g_Magnetar(t=2000 yr) = ${mag.compute_g_Magnetar(t_test).toExponential(6)} m/s^2\n`);
    
    // Step 5: Fluid & DM scaling
    console.log('Step 5: Fluid & DM Scaling (rho x2.0, DM factor x1.3)');
    mag.expandFluidDMScale(2.0, 1.3);
    console.log(`After expansion: rho_fluid = ${mag.getVariable('rho_fluid').toExponential(6)} kg/m^3, M_DM_factor = ${mag.getVariable('M_DM_factor').toFixed(3)}`);
    console.log(`g_Magnetar(t=2000 yr) = ${mag.compute_g_Magnetar(t_test).toExponential(6)} m/s^2\n`);
    
    // Step 6: State save/restore
    console.log('Step 6: State Management');
    mag.saveState('expanded_state');
    mag.setVariable('B0', 5e9);
    console.log(`Modified B0 to 5e9 T: g = ${mag.compute_g_Magnetar(t_test).toExponential(6)} m/s^2`);
    mag.restoreState('expanded_state');
    console.log(`Restored state: B0 = ${mag.getVariable('B0').toExponential(6)} T, g = ${mag.compute_g_Magnetar(t_test).toExponential(6)} m/s^2`);
    console.log(`Saved states: ${mag.listSavedStates().join(', ')}\n`);
    
    // Step 7: Sensitivity analysis
    console.log('Step 7: Sensitivity Analysis at t=2000 yr (top 5 parameters)');
    const sens = mag.sensitivityAnalysis(t_test, 1.0);
    const sens_vec = Object.entries(sens).sort((a, b) => b[1] - a[1]);
    for (let i = 0; i < Math.min(5, sens_vec.length); i++) {
        console.log(`  ${sens_vec[i][0]}: ${sens_vec[i][1].toExponential(6)}`);
    }
    console.log('');
    
    // Step 8: Generate variations
    console.log('Step 8: Generate Parameter Variations (3 variants, 10% variation)');
    const variations = mag.generateVariations(3, 10.0);
    for (let i = 0; i < variations.length; i++) {
        console.log(`  Variant ${i+1}: B0 = ${variations[i]['B0'].toExponential(6)} T, M = ${variations[i]['M'].toExponential(6)} kg`);
    }
    console.log('');
    
    // Step 9: Batch transformation
    console.log('Step 9: Batch Transform (scale mass parameters by 1.1)');
    mag.transformVariableGroup(['M', 'rho_fluid'], v => v * 1.1);
    console.log(`After transform: M = ${mag.getVariable('M').toExponential(6)} kg, rho_fluid = ${mag.getVariable('rho_fluid').toExponential(6)} kg/m^3`);
    console.log(`g_Magnetar(t=2000 yr) = ${mag.compute_g_Magnetar(t_test).toExponential(6)} m/s^2\n`);
    
    // Step 10: Consistency validation
    console.log('Step 10: Consistency Validation');
    const valid = mag.validateConsistency();
    console.log(`System is ${valid ? 'VALID' : 'INVALID'}\n`);
    
    // Step 11: Metric optimization
    console.log('Step 11: Optimize for Maximum g (t=0 to 10000 yr, 100 steps)');
    const max_g = mag.optimizeForMetric(g => g, 0.0, 10000 * 3.156e7, 100);
    console.log(`Maximum g found: ${max_g.toExponential(6)} m/s^2\n`);
    
    // Step 12: Full system report
    console.log('Step 12: Full System Report at t=3000 yr');
    const t_report = 3000 * 3.156e7;
    console.log(mag.generateReport(t_report) + '\n');
    
    // Step 13: B-field sweep
    console.log('Step 13: B-Field Sweep (B0 = 0.5e10, 1.0e10, 1.5e10, 2.0e10 T)');
    mag.saveState('before_sweep');
    for (const B0_val of [0.5e10, 1.0e10, 1.5e10, 2.0e10]) {
        mag.setVariable('B0', B0_val);
        console.log(`  B0 = ${B0_val.toExponential(2)} T: g(t=2000yr) = ${mag.compute_g_Magnetar(t_test).toExponential(6)} m/s^2`);
    }
    mag.restoreState('before_sweep');
    console.log('');
    
    // Step 14: Decay timescale sweep
    console.log('Step 14: Magnetic Decay Timescale Sweep (tau_B = 2000, 4000, 6000, 8000 yr)');
    for (const tau_yr of [2000.0, 4000.0, 6000.0, 8000.0]) {
        mag.setVariable('tau_B', tau_yr * 3.156e7);
        console.log(`  tau_B = ${tau_yr} yr: g(t=2000yr) = ${mag.compute_g_Magnetar(t_test).toExponential(6)} m/s^2`);
    }
    mag.restoreState('before_sweep');
    console.log('');
    
    // Step 15: Rotation period sweep
    console.log('Step 15: Initial Rotation Period Sweep (P_init = 3, 5, 7, 10 s)');
    for (const P_val of [3.0, 5.0, 7.0, 10.0]) {
        mag.setVariable('P_init', P_val);
        console.log(`  P_init = ${P_val} s: g(t=2000yr) = ${mag.compute_g_Magnetar(t_test).toExponential(6)} m/s^2, Omega(t=2000yr) = ${mag.Omega_t(t_test).toExponential(6)} rad/s`);
    }
    mag.restoreState('before_sweep');
    console.log('');
    
    // Step 16: DM factor sweep
    console.log('Step 16: Dark Matter Factor Sweep (M_DM_factor = 0.0, 0.1, 0.2, 0.3)');
    for (const dm of [0.0, 0.1, 0.2, 0.3]) {
        mag.setVariable('M_DM_factor', dm);
        console.log(`  M_DM_factor = ${dm}: g(t=2000yr) = ${mag.compute_g_Magnetar(t_test).toExponential(6)} m/s^2`);
    }
    mag.restoreState('before_sweep');
    console.log('');
    
    // Step 17: Export final state
    console.log('Step 17: Export Final State');
    console.log(mag.exportState() + '\n');
    
    console.log('========== ENHANCED EXAMPLE COMPLETE ==========\n');
}

// Run example if executed directly
if (typeof require !== 'undefined' && require.main === module) {
    enhancedMagnetarSGR0501Example();
}

// Export for Node.js and ES6 modules
if (typeof module !== 'undefined' && module.exports) {
    module.exports = MagnetarSGR0501_4516;
}

export default MagnetarSGR0501_4516;
