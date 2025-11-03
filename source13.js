/**
 * ================================================================================================
 * Module: MagnetarSGR1745_2900.js
 *
 * Description: JavaScript Module for SGR 1745-2900 Magnetar Class
 *              This is the second module in a series of 500+ code files for the Universal Quantum
 *              Field Framework (UQFF) simulations, focusing on magnetar evolution and gravity
 *              equations derived from Chandra X-ray Observatory datasets, high-energy lab simulations,
 *              and UQFF refinements (dated May 11, 2025, updated for full term inclusion on October 08, 2025).
 *
 * Purpose: Encapsulates the Master Universal Gravity Equation (MUGE) for SGR 1745-2900 magnetar
 *          evolution, including black hole proximity (Sgr A*), magnetic energy, and outburst decay.
 *          Includes ALL terms: base gravity, cosmic expansion (H(z)), magnetic decay, BH influence,
 *          UQFF Ug components with f_sc, Lambda, quantum uncertainty, EM, fluid, oscillatory waves,
 *          DM/density perturbations, magnetic energy (effective g), and decay power (cumulative energy effective g).
 *          Supports dynamic variable updates for all parameters.
 *
 * Integration: Designed for inclusion in base program or standalone execution.
 *              Instantiate class: const mag = new MagnetarSGR1745_2900();
 *              Compute: const g = mag.compute_g_Magnetar(t);
 *
 * Key Features:
 *   - Default values from UQFF document, with numerical computations for H(z), M_mag, etc.
 *   - Units handled: Energy terms converted to effective acceleration via / (M * r); decay as cumulative energy / (M * r).
 *   - Setter methods for updates: setVariable(name, value) or addToVariable(name, delta)/subtractFromVariable(name, delta).
 *   - Computes g_Magnetar(r, t) with every term explicitly included.
 *
 * Author: Converted to JavaScript by Copilot, based on Daniel T. Murphy's UQFF manuscript.
 * Date: November 03, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

class MagnetarSGR1745_2900 {
    constructor() {
        // Core parameters (mutable for updates)
        this.G = 0;
        this.M = 0;
        this.r = 0;
        this.Hz = 0;
        this.B0 = 0;
        this.tau_B = 0;
        this.B_crit = 0;
        this.Lambda = 0;
        this.c_light = 0;
        this.q_charge = 0;
        this.v_surf = 0;
        this.f_sc = 0;
        this.rho_vac_UA = 0;
        this.rho_vac_SCm = 0;
        this.P_init = 0;
        this.tau_Omega = 0;
        this.scale_EM = 0;
        this.proton_mass = 0;
        this.M_BH = 0;
        this.r_BH = 0;
        this.mu0 = 0;
        this.L0_W = 0;
        this.tau_decay = 0;

        // Additional parameters for full inclusion of terms
        this.hbar = 0;
        this.t_Hubble = 0;
        this.t_Hubble_gyr = 0;
        this.delta_x = 0;
        this.delta_p = 0;
        this.integral_psi = 0;
        this.rho_fluid = 0;
        this.A_osc = 0;
        this.k_osc = 0;
        this.omega_osc = 0;
        this.x_pos = 0;
        this.M_DM_factor = 0;
        this.delta_rho_over_rho = 0;

        // Computed caches
        this.ug1_base = 0;
        this.B = 0;

        // Extra variables storage
        this.extra_variables = {};

        // Initialize with defaults
        this.initializeDefaults();
    }

    initializeDefaults() {
        this.G = 6.6743e-11;
        this.M = 1.4 * 1.989e30;
        this.r = 1e4;
        this.Hz = 2.269e-18;
        this.B0 = 2e10;
        this.B = this.B0;
        this.tau_B = 4000 * 3.15576e7;
        this.B_crit = 1e11;
        this.Lambda = 1.1e-52;
        this.c_light = 3e8;
        this.q_charge = 1.602e-19;
        this.v_surf = 1e6;
        this.f_sc = 1 - (this.B / this.B_crit);
        this.rho_vac_UA = 7.09e-36;
        this.rho_vac_SCm = 7.09e-37;
        this.P_init = 3.76;
        this.tau_Omega = 10000 * 3.15576e7;
        this.scale_EM = 1e-12;
        this.proton_mass = 1.673e-27;
        this.M_BH = 4e6 * 1.989e30;
        this.r_BH = 2.83e16;
        this.mu0 = 4 * Math.PI * 1e-7;
        this.L0_W = 5e28;
        this.tau_decay = 3.5 * 365.25 * 24 * 3600;

        // Full terms defaults
        this.hbar = 1.0546e-34;
        this.t_Hubble = 13.8e9 * 3.15576e7;
        this.t_Hubble_gyr = 13.8;
        this.delta_x = 1e-10;
        this.delta_p = this.hbar / this.delta_x;
        this.integral_psi = 1.0;
        this.rho_fluid = 1e17;
        this.A_osc = 1e10;
        this.k_osc = 1.0 / this.r;
        this.omega_osc = 2 * Math.PI / this.P_init;
        this.x_pos = this.r;
        this.M_DM_factor = 0.1;
        this.delta_rho_over_rho = 1e-5;

        this.updateCache();
    }

    updateCache() {
        this.ug1_base = (this.G * this.M) / (this.r * this.r);
        this.f_sc = 1 - (this.B / this.B_crit);
    }

    setVariable(varName, newValue) {
        const varMap = {
            'G': () => this.G = newValue,
            'M': () => this.M = newValue,
            'r': () => this.r = newValue,
            'Hz': () => this.Hz = newValue,
            'B0': () => { this.B0 = newValue; this.B = newValue; },
            'tau_B': () => this.tau_B = newValue,
            'B_crit': () => this.B_crit = newValue,
            'Lambda': () => this.Lambda = newValue,
            'c_light': () => this.c_light = newValue,
            'q_charge': () => this.q_charge = newValue,
            'v_surf': () => this.v_surf = newValue,
            'f_sc': () => this.f_sc = newValue,
            'rho_vac_UA': () => this.rho_vac_UA = newValue,
            'rho_vac_SCm': () => this.rho_vac_SCm = newValue,
            'P_init': () => this.P_init = newValue,
            'tau_Omega': () => this.tau_Omega = newValue,
            'scale_EM': () => this.scale_EM = newValue,
            'proton_mass': () => this.proton_mass = newValue,
            'M_BH': () => this.M_BH = newValue,
            'r_BH': () => this.r_BH = newValue,
            'mu0': () => this.mu0 = newValue,
            'L0_W': () => this.L0_W = newValue,
            'tau_decay': () => this.tau_decay = newValue,
            'hbar': () => this.hbar = newValue,
            't_Hubble': () => this.t_Hubble = newValue,
            't_Hubble_gyr': () => this.t_Hubble_gyr = newValue,
            'delta_x': () => this.delta_x = newValue,
            'delta_p': () => this.delta_p = newValue,
            'integral_psi': () => this.integral_psi = newValue,
            'rho_fluid': () => this.rho_fluid = newValue,
            'A_osc': () => this.A_osc = newValue,
            'k_osc': () => this.k_osc = newValue,
            'omega_osc': () => this.omega_osc = newValue,
            'x_pos': () => this.x_pos = newValue,
            'M_DM_factor': () => this.M_DM_factor = newValue,
            'delta_rho_over_rho': () => this.delta_rho_over_rho = newValue
        };

        if (varMap[varName]) {
            varMap[varName]();
            this.updateCache();
            return true;
        } else {
            console.error(`Error: Unknown variable '${varName}'.`);
            return false;
        }
    }

    addToVariable(varName, delta) {
        return this.setVariable(varName, this.getVariable(varName) + delta);
    }

    subtractFromVariable(varName, delta) {
        return this.addToVariable(varName, -delta);
    }

    getVariable(varName) {
        const varMap = {
            'G': this.G, 'M': this.M, 'r': this.r, 'Hz': this.Hz,
            'B0': this.B0, 'tau_B': this.tau_B, 'B_crit': this.B_crit,
            'Lambda': this.Lambda, 'c_light': this.c_light, 'q_charge': this.q_charge,
            'v_surf': this.v_surf, 'f_sc': this.f_sc, 'rho_vac_UA': this.rho_vac_UA,
            'rho_vac_SCm': this.rho_vac_SCm, 'P_init': this.P_init, 'tau_Omega': this.tau_Omega,
            'scale_EM': this.scale_EM, 'proton_mass': this.proton_mass, 'M_BH': this.M_BH,
            'r_BH': this.r_BH, 'mu0': this.mu0, 'L0_W': this.L0_W, 'tau_decay': this.tau_decay,
            'hbar': this.hbar, 't_Hubble': this.t_Hubble, 't_Hubble_gyr': this.t_Hubble_gyr,
            'delta_x': this.delta_x, 'delta_p': this.delta_p, 'integral_psi': this.integral_psi,
            'rho_fluid': this.rho_fluid, 'A_osc': this.A_osc, 'k_osc': this.k_osc,
            'omega_osc': this.omega_osc, 'x_pos': this.x_pos, 'M_DM_factor': this.M_DM_factor,
            'delta_rho_over_rho': this.delta_rho_over_rho
        };

        if (varMap.hasOwnProperty(varName)) {
            return varMap[varName];
        } else {
            console.error(`Error: Unknown variable '${varName}'.`);
            return 0.0;
        }
    }

    B_t(t) {
        return this.B;
    }

    Omega_t(t) {
        return (2 * Math.PI / this.P_init) * Math.exp(-t / this.tau_Omega);
    }

    dOmega_dt(t) {
        const omega0 = 2 * Math.PI / this.P_init;
        return omega0 * (-1.0 / this.tau_Omega) * Math.exp(-t / this.tau_Omega);
    }

    compute_Ug() {
        const Ug1 = this.ug1_base;
        const Ug2 = 0.0;
        const Ug3 = 0.0;
        const Ug4 = Ug1 * this.f_sc;
        return Ug1 + Ug2 + Ug3 + Ug4;
    }

    compute_V() {
        return (4.0 / 3.0) * Math.PI * this.r * this.r * this.r;
    }

    compute_M_mag() {
        const V = this.compute_V();
        return (this.B_t(0) * this.B_t(0) / (2 * this.mu0)) * V;
    }

    compute_cumulative_D(t) {
        const exp_term = Math.exp(-t / this.tau_decay);
        return this.L0_W * this.tau_decay * (1 - exp_term);
    }

    compute_g_Magnetar(t) {
        if (t < 0) {
            console.error("Error: Time t must be non-negative.");
            return 0.0;
        }

        const Bt = this.B_t(t);
        const dOdt = this.dOmega_dt(t);

        // f_sc update
        const current_f_sc = 1 - (Bt / this.B_crit);

        // Term 1: Base + H(z) + B corrections
        const corr_H = 1 + this.Hz * t;
        const corr_B = current_f_sc;
        const term1 = this.ug1_base * corr_H * corr_B;

        // BH term
        const term_BH = (this.G * this.M_BH) / (this.r_BH * this.r_BH);

        // Term 2: UQFF Ug
        const term2 = this.compute_Ug();

        // Term 3: Lambda
        const term3 = (this.Lambda * this.c_light * this.c_light) / 3.0;

        // Term 4: Scaled EM (v x B magnitude)
        const cross_vB = this.v_surf * Bt;
        const em_base = (this.q_charge * cross_vB) / this.proton_mass;
        const term4 = em_base * this.scale_EM;

        // Term 5: GW
        const gw_prefactor = (this.G * this.M * this.M) / (Math.pow(this.c_light, 4) * this.r);
        const term5 = gw_prefactor * (dOdt * dOdt);

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

        // Magnetic energy term (effective g)
        const M_mag = this.compute_M_mag();
        const term_mag = M_mag / (this.M * this.r);

        // Decay term (cumulative energy effective g)
        const cum_D = this.compute_cumulative_D(t);
        const term_decay = cum_D / (this.M * this.r);

        // Total g_Magnetar (all terms summed)
        return term1 + term_BH + term2 + term3 + term4 + term5 + term_q + term_fluid + term_osc + term_DM + term_mag + term_decay;
    }

    printParameters() {
        console.log("SGR 1745-2900 Parameters:");
        console.log(`G: ${this.G.toExponential(3)}, M: ${this.M.toExponential(3)}, r: ${this.r.toExponential(3)}`);
        console.log(`Hz: ${this.Hz.toExponential(3)}, B: ${this.B.toExponential(3)}, M_BH: ${this.M_BH.toExponential(3)}, r_BH: ${this.r_BH.toExponential(3)}`);
        console.log(`L0_W: ${this.L0_W.toExponential(3)}, tau_decay: ${this.tau_decay.toExponential(3)}`);
        console.log(`f_sc: ${this.f_sc.toFixed(3)}, rho_fluid: ${this.rho_fluid.toExponential(3)}, M_DM_factor: ${this.M_DM_factor.toFixed(3)}`);
        console.log(`A_osc: ${this.A_osc.toExponential(3)}, delta_rho_over_rho: ${this.delta_rho_over_rho.toExponential(3)}`);
        const M_mag = this.compute_M_mag();
        console.log(`M_mag (J): ${M_mag.toExponential(3)}, ug1_base: ${this.ug1_base.toExponential(3)}`);
    }

    exampleAtOneYear() {
        const t_example = 1.0 * 365.25 * 24 * 3600;
        return this.compute_g_Magnetar(t_example);
    }

    // ===== ENHANCED DYNAMIC CAPABILITIES =====

    // Variable Management
    createVariable(name, value) {
        this.extra_variables[name] = value;
    }

    removeVariable(name) {
        delete this.extra_variables[name];
    }

    cloneVariable(source, dest) {
        if (this.extra_variables.hasOwnProperty(source)) {
            this.extra_variables[dest] = this.extra_variables[source];
        } else {
            const val = this.getVariable(source);
            this.extra_variables[dest] = val;
        }
    }

    listVariables() {
        const names = [
            'G', 'M', 'r', 'Hz', 'B0', 'tau_B', 'B_crit', 'Lambda', 'c_light',
            'q_charge', 'v_surf', 'f_sc', 'rho_vac_UA', 'rho_vac_SCm', 'P_init',
            'tau_Omega', 'scale_EM', 'proton_mass', 'M_BH', 'r_BH', 'mu0', 'L0_W',
            'tau_decay', 'hbar', 't_Hubble', 't_Hubble_gyr', 'delta_x', 'delta_p',
            'integral_psi', 'rho_fluid', 'A_osc', 'k_osc', 'omega_osc', 'x_pos',
            'M_DM_factor', 'delta_rho_over_rho'
        ];
        for (const key in this.extra_variables) {
            names.push(key);
        }
        return names;
    }

    getSystemName() {
        return "SGR 1745-2900 Magnetar near Sgr A* - Full UQFF & SM Integration";
    }

    // Batch Operations
    transformVariableGroup(names, func) {
        for (const name of names) {
            const current = this.getVariable(name);
            this.setVariable(name, func(current));
        }
    }

    scaleVariableGroup(names, factor) {
        this.transformVariableGroup(names, (v) => v * factor);
    }

    // Self-Expansion (domain-specific for Magnetar)
    expandParameterSpace(scale_factor) {
        this.r *= scale_factor;
        this.M *= scale_factor;
        this.L0_W *= scale_factor;
        this.updateCache();
    }

    expandMagneticScale(scale_factor) {
        this.B0 *= scale_factor;
        this.B *= scale_factor;
        this.tau_B *= scale_factor;
        this.updateCache();
    }

    expandDecayScale(scale_factor) {
        this.L0_W *= scale_factor;
        this.tau_decay *= scale_factor;
    }

    expandBlackHoleScale(scale_factor) {
        this.M_BH *= scale_factor;
        this.r_BH *= scale_factor;
    }

    // Self-Refinement
    autoRefineParameters(tolerance) {
        if (this.M <= 0) this.M = 1.4 * 1.989e30;
        if (this.r <= 0) this.r = 1e4;
        if (this.B0 <= 0) this.B0 = 2e10;
        if (this.B <= 0) this.B = this.B0;
        if (this.L0_W <= 0) this.L0_W = 5e28;
        if (this.tau_decay <= 0) this.tau_decay = 3.5 * 365.25 * 24 * 3600;
        if (this.M_BH <= 0) this.M_BH = 4e6 * 1.989e30;
        if (this.r_BH <= 0) this.r_BH = 2.83e16;
        if (this.rho_fluid <= 0) this.rho_fluid = 1e17;
        this.updateCache();
    }

    calibrateToObservations(obs_data) {
        for (const key in obs_data) {
            this.setVariable(key, obs_data[key]);
        }
        this.updateCache();
    }

    optimizeForMetric(var_name, target_value, iterations) {
        const current = this.getVariable(var_name);
        if (current === 0.0 && var_name !== "unknown") {
            return;
        }

        let best_value = current;
        let best_error = Math.abs(current - target_value);

        for (let i = 0; i < iterations; i++) {
            const test_value = current * (0.9 + Math.random() * 0.2);
            const error = Math.abs(test_value - target_value);
            if (error < best_error) {
                best_error = error;
                best_value = test_value;
            }
        }
        this.setVariable(var_name, best_value);
    }

    // Parameter Exploration
    generateVariations(n_variations) {
        const variations = [];
        for (let i = 0; i < n_variations; i++) {
            const variation = {
                'B0': this.B0 * (0.8 + Math.random() * 0.4),
                'L0_W': this.L0_W * (0.8 + Math.random() * 0.4),
                'tau_decay': this.tau_decay * (0.8 + Math.random() * 0.4),
                'M_BH': this.M_BH * (0.8 + Math.random() * 0.4),
                'r_BH': this.r_BH * (0.8 + Math.random() * 0.4)
            };
            variations.push(variation);
        }
        return variations;
    }

    // Adaptive Evolution
    mutateParameters(mutation_rate) {
        this.B0 *= (1.0 + (Math.random() - 0.5) * 2 * mutation_rate);
        this.B *= (1.0 + (Math.random() - 0.5) * 2 * mutation_rate);
        this.L0_W *= (1.0 + (Math.random() - 0.5) * 2 * mutation_rate);
        this.tau_decay *= (1.0 + (Math.random() - 0.5) * 2 * mutation_rate);
        this.updateCache();
    }

    evolveSystem(generations, fitness_function) {
        let best_fitness = fitness_function();

        const saved_B0 = this.B0;
        const saved_B = this.B;
        const saved_L0W = this.L0_W;
        const saved_tau = this.tau_decay;

        for (let gen = 0; gen < generations; gen++) {
            this.mutateParameters(0.1);
            const fitness = fitness_function();
            if (fitness > best_fitness) {
                best_fitness = fitness;
            } else {
                this.B0 = saved_B0;
                this.B = saved_B;
                this.L0_W = saved_L0W;
                this.tau_decay = saved_tau;
                this.updateCache();
            }
        }
    }

    // State Management
    saveState(label) {
        if (!MagnetarSGR1745_2900.saved_states) {
            MagnetarSGR1745_2900.saved_states = {};
        }
        const state = {};
        const var_names = this.listVariables();
        for (const name of var_names) {
            state[name] = this.getVariable(name);
        }
        MagnetarSGR1745_2900.saved_states[label] = state;
    }

    restoreState(label) {
        if (MagnetarSGR1745_2900.saved_states && MagnetarSGR1745_2900.saved_states[label]) {
            const state = MagnetarSGR1745_2900.saved_states[label];
            for (const key in state) {
                this.setVariable(key, state[key]);
            }
        }
    }

    listSavedStates() {
        if (!MagnetarSGR1745_2900.saved_states) {
            return [];
        }
        return Object.keys(MagnetarSGR1745_2900.saved_states);
    }

    exportState(t) {
        let output = `SGR 1745-2900 State Export at t=${t.toExponential(6)} s (${(t/(365.25*24*3600)).toFixed(3)} yr):\n`;
        output += `M=${this.M.toExponential(6)} kg (${(this.M/1.989e30).toFixed(3)} M_sun)\n`;
        output += `r=${this.r.toExponential(6)} m (${(this.r/1e3).toFixed(3)} km)\n`;
        output += `B0=${this.B0.toExponential(6)} T, B(t)=${this.B_t(t).toExponential(6)} T\n`;
        output += `L0_W=${this.L0_W.toExponential(6)} W, tau_decay=${this.tau_decay.toExponential(6)} s (${(this.tau_decay/(365.25*24*3600)).toFixed(3)} yr)\n`;
        output += `M_BH=${this.M_BH.toExponential(6)} kg (${(this.M_BH/1.989e30).toFixed(3)} M_sun), r_BH=${this.r_BH.toExponential(6)} m\n`;
        output += `rho_fluid=${this.rho_fluid.toExponential(6)} kg/m³, A_osc=${this.A_osc.toExponential(6)} m/s²\n`;
        output += `M_mag=${this.compute_M_mag().toExponential(6)} J, cumulative_D=${this.compute_cumulative_D(t).toExponential(6)} J\n`;
        output += `g_total=${this.compute_g_Magnetar(t).toExponential(6)} m/s²\n`;
        return output;
    }

    // System Analysis
    sensitivityAnalysis(param, t, delta) {
        const result = {};
        const original = this.getVariable(param);
        if (original === 0.0 && param !== "unknown") {
            result['error'] = -1;
            return result;
        }

        const g_original = this.compute_g_Magnetar(t);

        this.setVariable(param, original * (1.0 + delta));
        const g_plus = this.compute_g_Magnetar(t);

        this.setVariable(param, original * (1.0 - delta));
        const g_minus = this.compute_g_Magnetar(t);

        this.setVariable(param, original);

        result[`dg/d${param}`] = (g_plus - g_minus) / (2.0 * delta * original);
        result['g_original'] = g_original;
        result['g_plus'] = g_plus;
        result['g_minus'] = g_minus;

        return result;
    }

    generateReport(t) {
        let report = `===== SGR 1745-2900 MAGNETAR UQFF Module Report (t=${t.toExponential(6)} s) =====\n\n`;
        report += `System: ${this.getSystemName()}\n\n`;

        report += `Core Parameters:\n`;
        report += `  M = ${this.M.toExponential(6)} kg (${(this.M/1.989e30).toFixed(3)} M_sun)\n`;
        report += `  r = ${this.r.toExponential(6)} m (${(this.r/1e3).toFixed(3)} km)\n`;
        report += `  B0 = ${this.B0.toExponential(6)} T, B(t) = ${this.B_t(t).toExponential(6)} T\n`;
        report += `  B_crit = ${this.B_crit.toExponential(6)} T\n`;
        report += `  P_init = ${this.P_init.toFixed(3)} s (pulse period)\n`;
        report += `  L0_W = ${this.L0_W.toExponential(6)} W (initial luminosity)\n`;
        report += `  tau_decay = ${this.tau_decay.toExponential(6)} s (${(this.tau_decay/(365.25*24*3600)).toFixed(3)} yr)\n`;
        report += `  M_BH = ${this.M_BH.toExponential(6)} kg (${(this.M_BH/1.989e30).toFixed(3)} M_sun, Sgr A*)\n`;
        report += `  r_BH = ${this.r_BH.toExponential(6)} m (distance to Sgr A*)\n\n`;

        const Bt = this.B_t(t);
        const dOdt = this.dOmega_dt(t);
        const current_f_sc = 1 - (Bt / this.B_crit);
        const corr_H = 1 + this.Hz * t;
        const corr_B = current_f_sc;
        const term1 = this.ug1_base * corr_H * corr_B;
        const term_BH = (this.G * this.M_BH) / (this.r_BH * this.r_BH);
        const term2 = this.compute_Ug();
        const term3 = (this.Lambda * this.c_light * this.c_light) / 3.0;
        const cross_vB = this.v_surf * Bt;
        const em_base = (this.q_charge * cross_vB) / this.proton_mass;
        const term4 = em_base * this.scale_EM;
        const gw_prefactor = (this.G * this.M * this.M) / (Math.pow(this.c_light, 4) * this.r);
        const term5 = gw_prefactor * (dOdt * dOdt);
        const sqrt_unc = Math.sqrt(this.delta_x * this.delta_p);
        const term_q = (this.hbar / sqrt_unc) * this.integral_psi * (2 * Math.PI / this.t_Hubble);
        const V = this.compute_V();
        const term_fluid = (this.rho_fluid * V * this.ug1_base) / this.M;
        const term_osc1 = 2 * this.A_osc * Math.cos(this.k_osc * this.x_pos) * Math.cos(this.omega_osc * t);
        const arg = this.k_osc * this.x_pos - this.omega_osc * t;
        const term_osc2 = (2 * Math.PI / this.t_Hubble_gyr) * this.A_osc * Math.cos(arg);
        const term_osc = term_osc1 + term_osc2;
        const M_dm = this.M * this.M_DM_factor;
        const pert1 = this.delta_rho_over_rho;
        const pert2 = 3 * this.G * this.M / (this.r * this.r * this.r);
        const term_dm_force_like = (this.M + M_dm) * (pert1 + pert2);
        const term_DM = term_dm_force_like / this.M;
        const M_mag = this.compute_M_mag();
        const term_mag = M_mag / (this.M * this.r);
        const cum_D = this.compute_cumulative_D(t);
        const term_decay = cum_D / (this.M * this.r);
        const g_total = this.compute_g_Magnetar(t);

        report += `Term Breakdown:\n`;
        report += `  Base gravity (H(z), B corrections) = ${term1.toExponential(6)} m/s²\n`;
        report += `  BH term (Sgr A*) = ${term_BH.toExponential(6)} m/s²\n`;
        report += `  Ug_sum (Ug1+Ug2+Ug3+Ug4) = ${term2.toExponential(6)} m/s²\n`;
        report += `  Lambda_term = ${term3.toExponential(6)} m/s²\n`;
        report += `  EM_term (v×B) = ${term4.toExponential(6)} m/s²\n`;
        report += `  GW_term (dΩ/dt)² = ${term5.toExponential(6)} m/s²\n`;
        report += `  Quantum_term = ${term_q.toExponential(6)} m/s²\n`;
        report += `  Fluid_term = ${term_fluid.toExponential(6)} m/s²\n`;
        report += `  Oscillatory_term = ${term_osc.toExponential(6)} m/s²\n`;
        report += `  DM_term = ${term_DM.toExponential(6)} m/s²\n`;
        report += `  Magnetic_energy_term = ${term_mag.toExponential(6)} m/s² [M_mag=${M_mag.toExponential(6)} J]\n`;
        report += `  Decay_energy_term = ${term_decay.toExponential(6)} m/s² [cumulative_D=${cum_D.toExponential(6)} J]\n\n`;
        report += `TOTAL g = ${g_total.toExponential(6)} m/s²\n\n`;

        report += `Physics Notes:\n`;
        report += `- SGR 1745-2900 is a magnetar near Sgr A* supermassive black hole\n`;
        report += `- Ultra-strong magnetic field B~${this.B0.toExponential(3)} T with superconductivity correction f_sc=${current_f_sc.toFixed(6)}\n`;
        report += `- X-ray outburst decay L(t)=L0*exp(-t/tau), cumulative energy affects gravity\n`;
        report += `- Black hole proximity (r_BH=${this.r_BH.toExponential(3)} m) adds tidal term\n`;
        report += `- Full UQFF+SM: gravity, BH, Ug1-4, Lambda, EM, GW, quantum, fluid, oscillatory, DM, magnetic energy, decay\n`;
        report += `- Rotation period P=${this.P_init.toFixed(3)} s, slowing with tau_Omega=${this.tau_Omega.toExponential(3)} s\n`;
        report += `- Dense NS interior: rho_fluid=${this.rho_fluid.toExponential(3)} kg/m³\n\n`;

        return report;
    }

    validateConsistency() {
        let valid = true;
        if (this.M <= 0) valid = false;
        if (this.r <= 0) valid = false;
        if (this.B0 <= 0) valid = false;
        if (this.B <= 0) valid = false;
        if (this.L0_W <= 0) valid = false;
        if (this.tau_decay <= 0) valid = false;
        if (this.M_BH <= 0) valid = false;
        if (this.r_BH <= 0) valid = false;
        if (this.rho_fluid <= 0) valid = false;
        return valid;
    }

    autoCorrectAnomalies() {
        if (this.M <= 0) this.M = 1.4 * 1.989e30;
        if (this.r <= 0) this.r = 1e4;
        if (this.B0 <= 0) this.B0 = 2e10;
        if (this.B <= 0) this.B = this.B0;
        if (this.L0_W <= 0) this.L0_W = 5e28;
        if (this.tau_decay <= 0) this.tau_decay = 3.5 * 365.25 * 24 * 3600;
        if (this.M_BH <= 0) this.M_BH = 4e6 * 1.989e30;
        if (this.r_BH <= 0) this.r_BH = 2.83e16;
        if (this.rho_fluid <= 0) this.rho_fluid = 1e17;
        this.updateCache();
    }
}

// Enhanced example usage demonstration
function enhanced_magnetar_example() {
    const mag = new MagnetarSGR1745_2900();
    const t_1yr = 365.25 * 24 * 3600;

    console.log("===== ENHANCED SGR 1745-2900 MAGNETAR MODULE DEMONSTRATION =====\n");

    // Step 1: Variable management
    console.log("Step 1: Variable Management");
    mag.createVariable("custom_outburst_factor", 1.5);
    mag.cloneVariable("B0", "B0_backup");
    const vars = mag.listVariables();
    console.log(`Total variables: ${vars.length}`);
    console.log(`System: ${mag.getSystemName()}\n`);

    // Step 2: Batch scaling
    console.log("Step 2: Batch Scaling (Magnetic parameters)");
    mag.scaleVariableGroup(["B0", "L0_W", "tau_decay"], 1.1);
    console.log("Scaled B0, L0_W, tau_decay by 1.1\n");

    // Step 3: Self-expansion
    console.log("Step 3: Self-Expansion");
    mag.expandMagneticScale(1.05);
    console.log("Expanded magnetic scale +5%");
    mag.expandDecayScale(1.03);
    console.log("Expanded decay scale +3%");
    mag.expandBlackHoleScale(1.02);
    console.log("Expanded black hole scale +2%\n");

    // Step 4: Self-refinement
    console.log("Step 4: Self-Refinement");
    mag.autoRefineParameters(1e-10);
    console.log("Auto-refined parameters");
    const obs_data = {
        'B0': 2.2e10,
        'L0_W': 5.5e28,
        'tau_decay': 3.6 * 365.25 * 24 * 3600
    };
    mag.calibrateToObservations(obs_data);
    console.log("Calibrated to observations\n");

    // Step 5: Optimize
    console.log("Step 5: Optimize for B0~2e10 T");
    mag.optimizeForMetric("B0", 2e10, 50);
    console.log("Optimization complete\n");

    // Step 6: Generate variations
    console.log("Step 6: Generate 15 Parameter Variations");
    const variations = mag.generateVariations(15);
    console.log(`Generated ${variations.length} variations\n`);

    // Step 7: State management
    console.log("Step 7: State Management");
    mag.saveState("initial");
    mag.scaleVariableGroup(["B0", "L0_W"], 1.2);
    mag.saveState("enhanced_magnetic");
    mag.expandDecayScale(0.8);
    mag.saveState("reduced_decay");
    console.log("Saved 3 states\n");

    // Step 8: Sensitivity analysis
    console.log("Step 8: Sensitivity Analysis (B0 at t=1yr)");
    mag.restoreState("initial");
    const sensitivity = mag.sensitivityAnalysis("B0", t_1yr, 0.1);
    console.log(`dg/dB0 = ${sensitivity['dg/dB0'].toExponential(6)} (m/s²)/T\n`);

    // Step 9: System validation
    console.log("Step 9: System Validation");
    let valid = mag.validateConsistency();
    console.log(`System consistency: ${valid ? "VALID" : "INVALID"}`);
    if (!valid) {
        mag.autoCorrectAnomalies();
        console.log("Auto-corrected anomalies");
    }
    console.log();

    // Step 10: Comprehensive report
    console.log("Step 10: Comprehensive Report (t=1yr)");
    const report = mag.generateReport(t_1yr);
    console.log(report);

    // Step 11: Time evolution
    console.log("Step 11: Time Evolution (0 to 10 years)");
    const times = [0.0, 1*t_1yr, 2*t_1yr, 3.5*t_1yr, 5*t_1yr, 10*t_1yr];
    for (const t of times) {
        const g = mag.compute_g_Magnetar(t);
        console.log(`t=${t.toExponential(6)} s (${(t/(365.25*24*3600)).toFixed(3)} yr): g=${g.toExponential(6)} m/s²`);
    }
    console.log();

    console.log("===== DEMONSTRATION COMPLETE =====");
}

// Export for Node.js or browser
if (typeof module !== 'undefined' && module.exports) {
    module.exports = { MagnetarSGR1745_2900, enhanced_magnetar_example };
}

// Run example if executed directly
if (typeof require !== 'undefined' && require.main === module) {
    enhanced_magnetar_example();
}
