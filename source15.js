/**
 * ================================================================================================
 * Module: SMBHSgrAStar.js
 *
 * Description: JavaScript Module for Sagittarius A* (Sgr A*) Supermassive Black Hole Class
 *              This is the third module in a series of 500+ code files for the Universal Quantum
 *              Field Framework (UQFF) simulations, focusing on SMBH evolution and gravity
 *              equations derived from Hubble datasets, high-energy lab simulations, and UQFF
 *              refinements (dated May 09, 2025, updated for full term inclusion on October 08, 2025).
 *
 * Purpose: Encapsulates the Master Universal Gravity Equation (MUGE) for Sgr A* evolution.
 *          Includes ALL terms: base gravity with mass growth M(t), cosmic expansion (H_0), magnetic decay,
 *          UQFF Ug components with f_TRZ, Lambda, quantum uncertainty, EM (with B(t)), fluid dynamics,
 *          oscillatory waves, DM/density perturbations with precession sin(30°), and GW term.
 *          Supports dynamic variable updates for all parameters.
 *
 * Integration: Designed for Node.js and browser environments.
 *              Node.js: const SMBHSgrAStar = require('./source15.js');
 *              ES6: import SMBHSgrAStar from './source15.js';
 *              Instantiate: const sgrA = new SMBHSgrAStar();
 *              Compute: const g = sgrA.compute_g_SgrA(t);
 *
 * Key Features:
 *   - Default values from UQFF document, with approximations for all terms.
 *   - Units handled: B(t) converted to T (1 G = 10^-4 T); energy terms converted to effective acceleration.
 *   - Setter methods for updates: setVariable(name, value) or addToVariable(name, delta)/subtractFromVariable(name, delta).
 *   - Computes g_SgrA(r, t) with every term explicitly included.
 *   - 25 enhanced dynamic capabilities for parameter exploration, self-expansion, optimization, etc.
 *
 * Author: Converted to JavaScript by GitHub Copilot, based on C++ code by Grok (xAI)
 *         Original UQFF manuscript by Daniel T. Murphy
 * Date: November 03, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

class SMBHSgrAStar {
    constructor() {
        // Core parameters (mutable for updates)
        this.G = 0;               // Gravitational constant
        this.M_initial = 0;       // Initial SMBH mass
        this.r = 0;               // Schwarzschild radius
        this.H0 = 0;              // Hubble constant (s^-1)
        this.B0_G = 0;            // Initial magnetic field (G)
        this.tau_B = 0;           // B decay timescale (s)
        this.B_crit = 0;          // Critical B field (T)
        this.Lambda = 0;          // Cosmological constant
        this.c_light = 0;         // Speed of light
        this.q_charge = 0;        // Charge (proton)
        this.v_surf = 0;          // Surface velocity (arbitrary for BH)
        this.f_TRZ = 0;           // Time-reversal factor
        this.M_dot_0 = 0;         // Initial mass accretion rate factor
        this.tau_acc = 0;         // Accretion timescale (s)
        this.spin_factor = 0;     // Spin factor (0.3)
        this.tau_Omega = 0;       // Omega decay timescale (s)

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
        this.precession_angle_deg = 0; // Precession angle (degrees)

        // Computed caches (updated on demand)
        this.ug1_base = 0;        // Cached Ug1 for initial M (will recompute with M(t))

        // Initialize with default UQFF values
        this.initializeDefaults();
    }

    // Initialization method (called in constructor)
    initializeDefaults() {
        this.G = 6.6743e-11;
        this.M_initial = 4.3e6 * 1.989e30;
        this.r = 1.27e10;
        this.H0 = 2.184e-18;
        this.B0_G = 1e4;  // G
        this.tau_B = 1e6 * 3.156e7;
        this.B_crit = 1e11;  // T
        this.Lambda = 1.1e-52;
        this.c_light = 3e8;
        this.q_charge = 1.602e-19;
        this.v_surf = 1e6;  // Arbitrary
        this.f_TRZ = 0.1;
        this.M_dot_0 = 0.01;
        this.tau_acc = 9e9 * 3.156e7;
        this.spin_factor = 0.3;
        this.tau_Omega = 9e9 * 3.156e7;

        // Full terms defaults
        this.hbar = 1.0546e-34;
        this.t_Hubble = 13.8e9 * 3.156e7;
        this.t_Hubble_gyr = 13.8;
        this.delta_x = 1e-10;
        this.delta_p = this.hbar / this.delta_x;
        this.integral_psi = 1.0;
        this.rho_fluid = 1e17;  // Arbitrary for accretion disk
        this.A_osc = 1e6;       // Scaled down for BH
        this.k_osc = 1.0 / this.r;
        this.omega_osc = 2 * Math.PI / (this.r / this.c_light);  // Orbital-like
        this.x_pos = this.r;
        this.M_DM_factor = 0.1;
        this.delta_rho_over_rho = 1e-5;
        this.precession_angle_deg = 30.0;

        this.updateCache();
    }

    // Cache update for efficiency (call after parameter changes)
    updateCache() {
        this.ug1_base = (this.G * this.M_initial) / (this.r * this.r);
    }

    // Universal setter for any variable (by name, for flexibility)
    setVariable(varName, newValue) {
        const validVars = {
            'G': () => { this.G = newValue; },
            'M_initial': () => { this.M_initial = newValue; },
            'r': () => { this.r = newValue; },
            'H0': () => { this.H0 = newValue; },
            'B0_G': () => { this.B0_G = newValue; },
            'tau_B': () => { this.tau_B = newValue; },
            'B_crit': () => { this.B_crit = newValue; },
            'Lambda': () => { this.Lambda = newValue; },
            'c_light': () => { this.c_light = newValue; },
            'q_charge': () => { this.q_charge = newValue; },
            'v_surf': () => { this.v_surf = newValue; },
            'f_TRZ': () => { this.f_TRZ = newValue; },
            'M_dot_0': () => { this.M_dot_0 = newValue; },
            'tau_acc': () => { this.tau_acc = newValue; },
            'spin_factor': () => { this.spin_factor = newValue; },
            'tau_Omega': () => { this.tau_Omega = newValue; },
            'hbar': () => { this.hbar = newValue; },
            't_Hubble': () => { this.t_Hubble = newValue; },
            't_Hubble_gyr': () => { this.t_Hubble_gyr = newValue; },
            'delta_x': () => { this.delta_x = newValue; },
            'delta_p': () => { this.delta_p = newValue; },
            'integral_psi': () => { this.integral_psi = newValue; },
            'rho_fluid': () => { this.rho_fluid = newValue; },
            'A_osc': () => { this.A_osc = newValue; },
            'k_osc': () => { this.k_osc = newValue; },
            'omega_osc': () => { this.omega_osc = newValue; },
            'x_pos': () => { this.x_pos = newValue; },
            'M_DM_factor': () => { this.M_DM_factor = newValue; },
            'delta_rho_over_rho': () => { this.delta_rho_over_rho = newValue; },
            'precession_angle_deg': () => { this.precession_angle_deg = newValue; }
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
            'M_initial': () => this.M_initial,
            'r': () => this.r,
            'H0': () => this.H0,
            'B0_G': () => this.B0_G,
            'tau_B': () => this.tau_B,
            'B_crit': () => this.B_crit,
            'Lambda': () => this.Lambda,
            'c_light': () => this.c_light,
            'q_charge': () => this.q_charge,
            'v_surf': () => this.v_surf,
            'f_TRZ': () => this.f_TRZ,
            'M_dot_0': () => this.M_dot_0,
            'tau_acc': () => this.tau_acc,
            'spin_factor': () => this.spin_factor,
            'tau_Omega': () => this.tau_Omega,
            'hbar': () => this.hbar,
            't_Hubble': () => this.t_Hubble,
            't_Hubble_gyr': () => this.t_Hubble_gyr,
            'delta_x': () => this.delta_x,
            'delta_p': () => this.delta_p,
            'integral_psi': () => this.integral_psi,
            'rho_fluid': () => this.rho_fluid,
            'A_osc': () => this.A_osc,
            'k_osc': () => this.k_osc,
            'omega_osc': () => this.omega_osc,
            'x_pos': () => this.x_pos,
            'M_DM_factor': () => this.M_DM_factor,
            'delta_rho_over_rho': () => this.delta_rho_over_rho,
            'precession_angle_deg': () => this.precession_angle_deg
        };

        if (validVars[varName]) {
            return validVars[varName]();
        } else {
            console.error(`Error: Unknown variable '${varName}'.`);
            return null;
        }
    }

    // M(t) computation
    M_t(t) {
        const M_dot = this.M_dot_0 * Math.exp(-t / this.tau_acc);
        return this.M_initial * (1 + M_dot);
    }

    // B(t) in T
    B_t(t) {
        const B_G = this.B0_G * Math.exp(-t / this.tau_B);
        return B_G * 1e-4;  // G to T
    }

    // Omega(t) computation
    Omega_t(t) {
        const omega0 = this.spin_factor * this.c_light / this.r;
        return omega0 * Math.exp(-t / this.tau_Omega);
    }

    // dOmega/dt computation
    dOmega_dt(t) {
        const omega0 = this.spin_factor * this.c_light / this.r;
        return omega0 * (-1.0 / this.tau_Omega) * Math.exp(-t / this.tau_Omega);
    }

    // Ug terms computation
    compute_Ug(Mt, Bt) {
        const Ug1 = (this.G * Mt) / (this.r * this.r);
        const Ug2 = 0.0;
        const Ug3 = 0.0;
        const corr_B = 1 - Bt / this.B_crit;
        const Ug4 = Ug1 * corr_B;
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1 + this.f_TRZ);
    }

    // Volume computation for fluid
    compute_V() {
        return (4.0 / 3.0) * Math.PI * this.r * this.r * this.r;
    }

    // Main MUGE computation (includes ALL terms)
    compute_g_SgrA(t) {
        if (t < 0) {
            console.error('Error: Time t must be non-negative.');
            return 0.0;
        }

        const Mt = this.M_t(t);
        const Bt = this.B_t(t);
        const dOdt = this.dOmega_dt(t);
        const ug1_t = (this.G * Mt) / (this.r * this.r);

        // Term 1: Base + H0 + B corrections
        const corr_H = 1 + this.H0 * t;
        const corr_B = 1 - Bt / this.B_crit;
        const term1 = ug1_t * corr_H * corr_B;

        // Term 2: UQFF Ug with f_TRZ
        const term2 = this.compute_Ug(Mt, Bt);

        // Term 3: Lambda
        const term3 = (this.Lambda * this.c_light * this.c_light) / 3.0;

        // Term 4: EM (v x B, no scaling or UA here)
        const cross_vB = this.v_surf * Bt;  // Magnitude
        const em_base = this.q_charge * cross_vB / 1.673e-27;  // Acceleration
        const term4 = em_base;

        // Term 5: GW
        const gw_prefactor = (this.G * Mt * Mt) / (Math.pow(this.c_light, 4) * this.r);
        const term5 = gw_prefactor * (dOdt * dOdt);

        // Quantum uncertainty term
        const sqrt_unc = Math.sqrt(this.delta_x * this.delta_p);
        const term_q = (this.hbar / sqrt_unc) * this.integral_psi * (2 * Math.PI / this.t_Hubble);

        // Fluid term (effective acceleration)
        const V = this.compute_V();
        const term_fluid = (this.rho_fluid * V * ug1_t) / Mt;

        // Oscillatory terms (real parts)
        const term_osc1 = 2 * this.A_osc * Math.cos(this.k_osc * this.x_pos) * Math.cos(this.omega_osc * t);
        const arg = this.k_osc * this.x_pos - this.omega_osc * t;
        const term_osc2 = (2 * Math.PI / this.t_Hubble_gyr) * this.A_osc * Math.cos(arg);
        const term_osc = term_osc1 + term_osc2;

        // DM and density perturbation term with precession (converted to acceleration)
        const M_dm = Mt * this.M_DM_factor;
        const sin_prec = Math.sin(this.precession_angle_deg * Math.PI / 180.0);
        const pert1 = this.delta_rho_over_rho;
        const pert2 = 3 * this.G * Mt / (this.r * this.r * this.r);
        const term_dm_force_like = (Mt + M_dm) * (pert1 + pert2 * sin_prec);
        const term_DM = term_dm_force_like / Mt;

        // Total g_SgrA (all terms summed)
        return term1 + term2 + term3 + term4 + term5 + term_q + term_fluid + term_osc + term_DM;
    }

    // Debug/Output method (for transparency in base program)
    printParameters() {
        console.log('Sgr A* Parameters:');
        console.log(`G: ${this.G.toExponential(3)}, M_initial: ${this.M_initial.toExponential(3)}, r: ${this.r.toExponential(3)}`);
        console.log(`H0: ${this.H0.toExponential(3)}, B0_G: ${this.B0_G.toExponential(3)}, tau_B: ${this.tau_B.toExponential(3)}`);
        console.log(`f_TRZ: ${this.f_TRZ.toFixed(3)}, M_dot_0: ${this.M_dot_0.toFixed(3)}, tau_acc: ${this.tau_acc.toExponential(3)}`);
        console.log(`rho_fluid: ${this.rho_fluid.toExponential(3)}, M_DM_factor: ${this.M_DM_factor.toFixed(3)}`);
        console.log(`A_osc: ${this.A_osc.toExponential(3)}, precession_angle_deg: ${this.precession_angle_deg.toFixed(3)}`);
        console.log(`ug1_base: ${this.ug1_base.toExponential(3)}`);
    }

    // Example computation at t=4.5 Gyr (for testing)
    exampleAt4_5Gyr() {
        const t_example = 4.5e9 * 3.156e7;
        return this.compute_g_SgrA(t_example);
    }

    // ========== ENHANCED DYNAMIC CAPABILITIES (25 methods) ==========

    // --- Variable Management (5 methods) ---
    createVariable(name, value) {
        return this.setVariable(name, value);
    }

    removeVariable(name) {
        console.warn(`Warning: Cannot remove built-in variable '${name}' in SMBH class.`);
        return false;
    }

    cloneVariable(src, dest) {
        const val = this.getVariable(src);
        if (val === null) return false;
        return this.setVariable(dest, val);
    }

    listVariables() {
        return ['G', 'M_initial', 'r', 'H0', 'B0_G', 'tau_B', 'B_crit', 'Lambda', 'c_light', 'q_charge',
                'v_surf', 'f_TRZ', 'M_dot_0', 'tau_acc', 'spin_factor', 'tau_Omega',
                'hbar', 't_Hubble', 't_Hubble_gyr', 'delta_x', 'delta_p', 'integral_psi', 'rho_fluid',
                'A_osc', 'k_osc', 'omega_osc', 'x_pos', 'M_DM_factor', 'delta_rho_over_rho', 'precession_angle_deg'];
    }

    getSystemName() {
        return 'SMBHSgrAStar';
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
        const expandable = ['M_initial', 'r', 'B0_G', 'tau_B', 'rho_fluid', 'A_osc', 'M_DM_factor'];
        return this.scaleVariableGroup(expandable, factor);
    }

    expandAccretionScale(M_dot_factor, tau_acc_factor) {
        this.setVariable('M_dot_0', this.getVariable('M_dot_0') * M_dot_factor);
        this.setVariable('tau_acc', this.getVariable('tau_acc') * tau_acc_factor);
    }

    expandMagneticScale(B_factor, tau_B_factor) {
        this.setVariable('B0_G', this.getVariable('B0_G') * B_factor);
        this.setVariable('B_crit', this.getVariable('B_crit') * B_factor);
        this.setVariable('tau_B', this.getVariable('tau_B') * tau_B_factor);
    }

    expandDMPrecessionScale(DM_factor, angle_factor) {
        this.setVariable('M_DM_factor', this.getVariable('M_DM_factor') * DM_factor);
        this.setVariable('delta_rho_over_rho', this.getVariable('delta_rho_over_rho') * DM_factor);
        this.setVariable('precession_angle_deg', this.getVariable('precession_angle_deg') * angle_factor);
    }

    // --- Self-Refinement (3 methods) ---
    autoRefineParameters(observations) {
        if (observations.length === 0) return;
        
        let sum_error = 0.0;
        for (const obs of observations) {
            const t = obs[0];
            const g_obs = obs[1];
            const g_calc = this.compute_g_SgrA(t);
            sum_error += Math.abs(g_calc - g_obs);
        }
        const avg_error = sum_error / observations.length;
        
        // Simple refinement: adjust M_dot_0 and tau_acc based on error
        if (avg_error > 1e-6) {
            const adj_factor = 1.0 - Math.min(0.1, avg_error / 1e6);
            this.setVariable('M_dot_0', this.getVariable('M_dot_0') * adj_factor);
            this.setVariable('tau_acc', this.getVariable('tau_acc') * (2.0 - adj_factor));
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
            const g = this.compute_g_SgrA(t);
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
        SMBHSgrAStar.savedStates[stateName] = state;
        return true;
    }

    restoreState(stateName) {
        const state = SMBHSgrAStar.savedStates[stateName];
        if (!state) return false;
        
        for (const [key, value] of Object.entries(state)) {
            this.setVariable(key, value);
        }
        return true;
    }

    listSavedStates() {
        return Object.keys(SMBHSgrAStar.savedStates);
    }

    exportState() {
        let output = 'SMBHSgrAStar State Export:\n';
        const vars = this.listVariables();
        for (const v of vars) {
            output += `${v} = ${this.getVariable(v).toExponential(6)}\n`;
        }
        return output;
    }

    // --- System Analysis (4 methods) ---
    sensitivityAnalysis(t, delta_pct) {
        const sensitivities = {};
        const g_base = this.compute_g_SgrA(t);
        
        const vars = this.listVariables();
        for (const v of vars) {
            if (v === 'c_light' || v === 'G' || v === 'hbar') continue; // Skip constants
            
            const original = this.getVariable(v);
            const delta = original * delta_pct / 100.0;
            
            this.setVariable(v, original + delta);
            const g_plus = this.compute_g_SgrA(t);
            this.setVariable(v, original);
            
            const sensitivity = (g_base !== 0.0) ? Math.abs((g_plus - g_base) / g_base) : 0.0;
            sensitivities[v] = sensitivity;
        }
        
        return sensitivities;
    }

    generateReport(t) {
        let report = '';
        report += '============================================\n';
        report += 'SGR A* SUPERMASSIVE BLACK HOLE REPORT\n';
        report += '============================================\n';
        report += `Time: t = ${t.toExponential(6)} s (${(t/3.156e7/1e9).toFixed(3)} Gyr)\n\n`;
        
        report += 'Physical Parameters:\n';
        report += `  Initial Mass M_initial = ${this.M_initial.toExponential(6)} kg (${(this.M_initial/1.989e30).toExponential(3)} M_sun)\n`;
        report += `  M(t) = ${this.M_t(t).toExponential(6)} kg (${(this.M_t(t)/1.989e30).toExponential(3)} M_sun)\n`;
        report += `  Schwarzschild radius r = ${this.r.toExponential(6)} m (${(this.r/1e3).toExponential(3)} km)\n`;
        report += `  B-field B0 = ${this.B0_G.toExponential(3)} G (B_crit = ${this.B_crit.toExponential(3)} T)\n`;
        report += `  B(t) = ${this.B_t(t).toExponential(6)} T\n`;
        report += `  Accretion M_dot_0 = ${this.M_dot_0.toFixed(3)}, tau_acc = ${this.tau_acc.toExponential(6)} s\n`;
        report += `  Spin factor = ${this.spin_factor.toFixed(3)}, Omega(t) = ${this.Omega_t(t).toExponential(6)} rad/s\n`;
        report += `  Fluid density rho = ${this.rho_fluid.toExponential(6)} kg/m^3\n`;
        report += `  DM factor = ${this.M_DM_factor.toFixed(3)}, precession angle = ${this.precession_angle_deg.toFixed(3)} deg\n\n`;
        
        report += 'Computed Acceleration:\n';
        report += `  g_SgrA(t) = ${this.compute_g_SgrA(t).toExponential(6)} m/s^2\n\n`;
        
        report += 'UQFF Terms:\n';
        const Mt = this.M_t(t);
        const Bt = this.B_t(t);
        const ug1_t = (this.G * Mt) / (this.r * this.r);
        const corr_H = 1 + this.H0 * t;
        const corr_B = 1 - Bt / this.B_crit;
        report += `  Base (with H0, B, M(t)): ${(ug1_t * corr_H * corr_B).toExponential(6)} m/s^2\n`;
        report += `  Ug total: ${this.compute_Ug(Mt, Bt).toExponential(6)} m/s^2\n`;
        report += `  Lambda: ${((this.Lambda * this.c_light * this.c_light) / 3.0).toExponential(6)} m/s^2\n`;
        
        const cross_vB = this.v_surf * Bt;
        const em_base = this.q_charge * cross_vB / 1.673e-27;
        report += `  EM: ${em_base.toExponential(6)} m/s^2\n`;
        
        const dOdt = this.dOmega_dt(t);
        const gw_prefactor = (this.G * Mt * Mt) / (Math.pow(this.c_light, 4) * this.r);
        report += `  GW: ${(gw_prefactor * dOdt * dOdt).toExponential(6)} m/s^2\n`;
        
        const sqrt_unc = Math.sqrt(this.delta_x * this.delta_p);
        report += `  Quantum: ${((this.hbar / sqrt_unc) * this.integral_psi * (2 * Math.PI / this.t_Hubble)).toExponential(6)} m/s^2\n`;
        
        const V = this.compute_V();
        report += `  Fluid: ${((this.rho_fluid * V * ug1_t) / Mt).toExponential(6)} m/s^2\n`;
        
        report += '  Oscillatory: (combined real parts)\n';
        
        const M_dm = Mt * this.M_DM_factor;
        const sin_prec = Math.sin(this.precession_angle_deg * Math.PI / 180.0);
        const pert1 = this.delta_rho_over_rho;
        const pert2 = 3 * this.G * Mt / (this.r * this.r * this.r);
        const term_dm_force_like = (Mt + M_dm) * (pert1 + pert2 * sin_prec);
        report += `  DM (with precession): ${(term_dm_force_like / Mt).toExponential(6)} m/s^2\n`;
        
        report += '============================================\n';
        return report;
    }

    validateConsistency() {
        let valid = true;
        
        if (this.M_initial <= 0 || this.r <= 0) { 
            console.error('Error: M_initial and r must be positive.'); 
            valid = false; 
        }
        if (this.B0_G < 0 || this.B_crit <= 0) { 
            console.error('Error: B0_G, B_crit must be non-negative/positive.'); 
            valid = false; 
        }
        if (this.tau_B <= 0 || this.tau_acc <= 0 || this.tau_Omega <= 0) { 
            console.error('Error: Decay/accretion timescales must be positive.'); 
            valid = false; 
        }
        if (this.rho_fluid < 0) { 
            console.error('Error: Fluid density must be non-negative.'); 
            valid = false; 
        }
        if (this.M_DM_factor < 0 || this.M_DM_factor > 1.0) { 
            console.warn('Warning: DM factor outside [0,1].'); 
        }
        if (this.spin_factor < 0 || this.spin_factor > 1.0) { 
            console.warn('Warning: Spin factor outside [0,1].'); 
        }
        
        return valid;
    }

    autoCorrectAnomalies() {
        let corrected = false;
        
        if (this.M_initial <= 0) { this.M_initial = 4.3e6 * 1.989e30; corrected = true; }
        if (this.r <= 0) { this.r = 1.27e10; corrected = true; }
        if (this.B0_G < 0) { this.B0_G = 1e4; corrected = true; }
        if (this.B_crit <= 0) { this.B_crit = 1e11; corrected = true; }
        if (this.tau_B <= 0) { this.tau_B = 1e6 * 3.156e7; corrected = true; }
        if (this.tau_acc <= 0) { this.tau_acc = 9e9 * 3.156e7; corrected = true; }
        if (this.tau_Omega <= 0) { this.tau_Omega = 9e9 * 3.156e7; corrected = true; }
        if (this.rho_fluid < 0) { this.rho_fluid = 1e17; corrected = true; }
        if (this.M_DM_factor < 0) { this.M_DM_factor = 0.1; corrected = true; }
        if (this.M_DM_factor > 1.0) { this.M_DM_factor = 1.0; corrected = true; }
        if (this.spin_factor < 0) { this.spin_factor = 0.0; corrected = true; }
        if (this.spin_factor > 1.0) { this.spin_factor = 1.0; corrected = true; }
        
        if (corrected) this.updateCache();
        return corrected;
    }
}

// Static property for state storage (class-level)
SMBHSgrAStar.savedStates = {};

// ========== ENHANCED EXAMPLE FUNCTION ==========
function enhancedSgrAExample() {
    console.log('\n========== ENHANCED SGR A* SUPERMASSIVE BLACK HOLE UQFF EXAMPLE ==========\n');
    
    const sgrA = new SMBHSgrAStar();
    
    // Step 1: Initial state
    console.log('Step 1: Initial Configuration');
    sgrA.printParameters();
    const t0 = 0.0;
    console.log(`g_SgrA(t=0) = ${sgrA.compute_g_SgrA(t0).toExponential(6)} m/s^2\n`);
    
    // Step 2: Time evolution (0, 1, 3, 5, 9 Gyr)
    console.log('Step 2: Time Evolution (0, 1, 3, 5, 9 Gyr)');
    for (const t_gyr of [0.0, 1.0, 3.0, 5.0, 9.0]) {
        const t = t_gyr * 1e9 * 3.156e7;
        console.log(`  t = ${t_gyr} Gyr: g = ${sgrA.compute_g_SgrA(t).toExponential(6)} m/s^2, M(t) = ${(sgrA.M_t(t)/1.989e30).toExponential(3)} M_sun, B(t) = ${sgrA.B_t(t).toExponential(3)} T`);
    }
    console.log('');
    
    // Step 3: Accretion scaling
    console.log('Step 3: Accretion Scaling (M_dot x1.5, tau_acc x0.8)');
    sgrA.expandAccretionScale(1.5, 0.8);
    const t_test = 4.5e9 * 3.156e7;
    console.log(`After expansion: M_dot_0 = ${sgrA.getVariable('M_dot_0').toFixed(4)}, tau_acc = ${sgrA.getVariable('tau_acc').toExponential(6)} s`);
    console.log(`g_SgrA(t=4.5 Gyr) = ${sgrA.compute_g_SgrA(t_test).toExponential(6)} m/s^2\n`);
    
    // Step 4: Magnetic field scaling
    console.log('Step 4: Magnetic Field Scaling (B0 x1.2, tau_B x1.3)');
    sgrA.expandMagneticScale(1.2, 1.3);
    console.log(`After expansion: B0_G = ${sgrA.getVariable('B0_G').toExponential(3)} G, tau_B = ${sgrA.getVariable('tau_B').toExponential(6)} s`);
    console.log(`g_SgrA(t=4.5 Gyr) = ${sgrA.compute_g_SgrA(t_test).toExponential(6)} m/s^2\n`);
    
    // Step 5: DM & precession scaling
    console.log('Step 5: DM & Precession Scaling (DM factor x1.4, angle x1.1)');
    sgrA.expandDMPrecessionScale(1.4, 1.1);
    console.log(`After expansion: M_DM_factor = ${sgrA.getVariable('M_DM_factor').toFixed(3)}, precession_angle_deg = ${sgrA.getVariable('precession_angle_deg').toFixed(2)} deg`);
    console.log(`g_SgrA(t=4.5 Gyr) = ${sgrA.compute_g_SgrA(t_test).toExponential(6)} m/s^2\n`);
    
    // Step 6: State save/restore
    console.log('Step 6: State Management');
    sgrA.saveState('expanded_state');
    sgrA.setVariable('M_initial', 5e6 * 1.989e30);
    console.log(`Modified M_initial to 5e6 M_sun: g = ${sgrA.compute_g_SgrA(t_test).toExponential(6)} m/s^2`);
    sgrA.restoreState('expanded_state');
    console.log(`Restored state: M_initial = ${(sgrA.getVariable('M_initial')/1.989e30).toExponential(3)} M_sun, g = ${sgrA.compute_g_SgrA(t_test).toExponential(6)} m/s^2`);
    console.log(`Saved states: ${sgrA.listSavedStates().join(', ')}\n`);
    
    // Step 7: Sensitivity analysis
    console.log('Step 7: Sensitivity Analysis at t=4.5 Gyr (top 5 parameters)');
    const sens = sgrA.sensitivityAnalysis(t_test, 1.0);
    const sens_vec = Object.entries(sens).sort((a, b) => b[1] - a[1]);
    for (let i = 0; i < Math.min(5, sens_vec.length); i++) {
        console.log(`  ${sens_vec[i][0]}: ${sens_vec[i][1].toExponential(6)}`);
    }
    console.log('');
    
    // Step 8: Generate variations
    console.log('Step 8: Generate Parameter Variations (3 variants, 10% variation)');
    const variations = sgrA.generateVariations(3, 10.0);
    for (let i = 0; i < variations.length; i++) {
        console.log(`  Variant ${i+1}: M_initial = ${(variations[i]['M_initial']/1.989e30).toExponential(3)} M_sun, B0_G = ${variations[i]['B0_G'].toExponential(3)} G`);
    }
    console.log('');
    
    // Step 9: Batch transformation
    console.log('Step 9: Batch Transform (scale accretion parameters by 1.15)');
    sgrA.transformVariableGroup(['M_dot_0', 'rho_fluid'], v => v * 1.15);
    console.log(`After transform: M_dot_0 = ${sgrA.getVariable('M_dot_0').toFixed(4)}, rho_fluid = ${sgrA.getVariable('rho_fluid').toExponential(6)} kg/m^3`);
    console.log(`g_SgrA(t=4.5 Gyr) = ${sgrA.compute_g_SgrA(t_test).toExponential(6)} m/s^2\n`);
    
    // Step 10: Consistency validation
    console.log('Step 10: Consistency Validation');
    const valid = sgrA.validateConsistency();
    console.log(`System is ${valid ? 'VALID' : 'INVALID'}\n`);
    
    // Step 11: Metric optimization
    console.log('Step 11: Optimize for Maximum g (t=0 to 10 Gyr, 100 steps)');
    const max_g = sgrA.optimizeForMetric(g => g, 0.0, 10e9 * 3.156e7, 100);
    console.log(`Maximum g found: ${max_g.toExponential(6)} m/s^2\n`);
    
    // Step 12: Full system report
    console.log('Step 12: Full System Report at t=5 Gyr');
    const t_report = 5e9 * 3.156e7;
    console.log(sgrA.generateReport(t_report) + '\n');
    
    console.log('========== ENHANCED EXAMPLE COMPLETE ==========\n');
}

// Run example if executed directly
if (typeof require !== 'undefined' && require.main === module) {
    enhancedSgrAExample();
}

// Export for Node.js and ES6 modules
if (typeof module !== 'undefined' && module.exports) {
    module.exports = SMBHSgrAStar;
}

export default SMBHSgrAStar;
