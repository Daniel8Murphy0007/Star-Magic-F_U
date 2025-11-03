/**
 * ================================================================================================
 * Module: source22.js (Bubble Nebula - NGC 7635)
 *
 * Description: JavaScript ES6 Module for Bubble Nebula Emission Nebula Class
 *              Converted from C++ source22.cpp - Module 12 in the UQFF series
 *              Focuses on emission nebula evolution with expansion dynamics
 *
 * Purpose: Implements Master Universal Gravity Equation (MUGE) for Bubble Nebula evolution
 *          Including ALL terms: base gravity (static M), cosmic expansion (H_0), magnetic
 *          correction, nebula expansion E(t), UQFF Ug components, Lambda, quantum uncertainty,
 *          scaled EM, fluid dynamics, oscillatory waves, DM/density perturbations, and stellar
 *          wind feedback (UNIQUE: E(t) expansion factor reduces gravity over time)
 *
 * Key Features:
 *   - M = 46 M☉ (static mass - central star)
 *   - r = 5 ly nebula radius
 *   - E(t) = E_0 × (1 - exp(-t/tau_exp)) - Nebula expansion factor (UNIQUE)
 *   - E_0 = 0.1 initial expansion factor
 *   - tau_exp = 4 Myr expansion timescale
 *   - v_wind = 1800 km/s stellar wind velocity
 *   - Gravity reduction: g_effective = g_base × (1 - E(t))
 *   - 32 physics parameters + 38 methods (13 core + 25 enhanced)
 *
 * Author: Converted by GitHub Copilot from Daniel T. Murphy's UQFF framework
 * Date: November 03, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

class BubbleNebula {
    constructor() {
        this.initializeDefaults();
    }

    initializeDefaults() {
        // Physical constants
        this.G = 6.6743e-11;           // Gravitational constant (m³/kg/s²)
        this.c_light = 3e8;            // Speed of light (m/s)
        this.hbar = 1.0546e-34;        // Reduced Planck's constant (J·s)
        this.q_charge = 1.602e-19;     // Elementary charge (C)
        this.proton_mass = 1.673e-27;  // Proton mass (kg)

        // Bubble Nebula specific parameters
        const M_sun = 1.989e30;        // Solar mass (kg)
        const ly_to_m = 9.461e15;      // Light-year to meters
        this.M = 46.0 * M_sun;         // Static mass: 46 M☉ (central star)
        this.r = 5.0 * ly_to_m;        // Nebula radius: 5 light-years

        // Cosmological parameters
        this.H0 = 2.184e-18;           // Hubble constant (s⁻¹)
        this.Lambda = 1.1e-52;         // Cosmological constant (m⁻²)
        this.t_Hubble = 13.8e9 * 3.156e7; // Hubble time (s)
        this.t_Hubble_gyr = 13.8;      // Hubble time (Gyr)

        // Magnetic field parameters
        this.B = 1e-6;                 // Static magnetic field (T)
        this.B_crit = 1e11;            // Critical magnetic field (T)

        // NEBULA EXPANSION parameters (UNIQUE to Bubble Nebula)
        this.E_0 = 0.1;                // Initial expansion factor (dimensionless, 0-1)
        this.tau_exp = 4e6 * 3.156e7;  // Expansion timescale: 4 Myr (s)

        // Stellar wind parameters
        this.rho_wind = 1e-21;         // Wind density (kg/m³)
        this.v_wind = 1.8e6;           // Wind velocity: 1800 km/s (m/s)

        // Fluid/gas parameters
        this.rho_fluid = 1e-21;        // Fluid density (kg/m³)
        this.gas_v = 1e5;              // Gas velocity for EM (m/s)

        // UQFF specific parameters
        this.f_TRZ = 0.1;              // Time-reversal zone factor
        this.rho_vac_UA = 7.09e-36;    // UA vacuum energy density (J/m³)
        this.rho_vac_SCm = 7.09e-37;   // SCm vacuum energy density (J/m³)
        this.scale_EM = 1e-12;         // EM scaling factor

        // Quantum uncertainty parameters
        this.delta_x = 1e-10;          // Position uncertainty (m)
        this.delta_p = this.hbar / this.delta_x; // Momentum uncertainty (kg·m/s)
        this.integral_psi = 1.0;       // Wavefunction integral approximation

        // Oscillatory wave parameters
        this.A_osc = 1e-10;            // Oscillatory amplitude (m/s²)
        this.k_osc = 1.0 / this.r;     // Wave number (1/m)
        this.omega_osc = 2 * Math.PI / (this.r / this.c_light); // Angular frequency (rad/s)
        this.x_pos = this.r;           // Position for oscillation (m)

        // Dark matter and density perturbation parameters
        this.M_DM_factor = 0.1;        // Dark matter mass fraction
        this.delta_rho_over_rho = 1e-5; // Density perturbation fraction

        // Update cached values
        this.updateCache();
    }

    updateCache() {
        // Cache base Ug1 = G*M/r² for efficiency
        this.ug1_base = (this.G * this.M) / (this.r * this.r);
    }

    // ========== CORE PHYSICS METHODS (13 methods) ==========

    /**
     * Compute time-dependent expansion factor E(t)
     * E(t) = E_0 × (1 - exp(-t/tau_exp))
     * UNIQUE to Bubble Nebula - represents nebula expansion reducing gravitational effect
     * E(t) → E_0 as t → ∞ (asymptotic expansion)
     */
    E_t(t) {
        return this.E_0 * (1 - Math.exp(-t / this.tau_exp));
    }

    /**
     * Compute UQFF Ug terms with expansion correction
     * Ug is multiplied by (1 - E(t)) to reduce gravity as nebula expands
     */
    compute_Ug(Et) {
        const Ug1 = this.ug1_base;
        const Ug2 = 0.0;
        const Ug3 = 0.0;
        const corr_B = 1 - this.B / this.B_crit;
        const Ug4 = Ug1 * corr_B;
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1 + this.f_TRZ) * (1 - Et);
    }

    /**
     * Compute nebula volume for fluid dynamics
     */
    compute_V() {
        return (4.0 / 3.0) * Math.PI * this.r * this.r * this.r;
    }

    /**
     * Main MUGE computation: g_Bubble(t) with ALL 9 terms
     * Includes expansion correction E(t) reducing gravity (UNIQUE)
     */
    compute_g_Bubble(t) {
        if (t < 0) {
            console.error('Error: Time t must be non-negative.');
            return 0.0;
        }

        const Et = this.E_t(t);

        // Term 1: Base gravity with H0, B, and E(t) corrections
        const corr_H = 1 + this.H0 * t;
        const corr_B = 1 - this.B / this.B_crit;
        const corr_E = 1 - Et;  // UNIQUE: Expansion reduces gravity
        const term1 = this.ug1_base * corr_H * corr_B * corr_E;

        // Term 2: UQFF Ug with f_TRZ and E(t) correction
        const term2 = this.compute_Ug(Et);

        // Term 3: Lambda (cosmological constant)
        const term3 = (this.Lambda * this.c_light * this.c_light) / 3.0;

        // Term 4: Scaled EM with UA vacuum correction
        const cross_vB = this.gas_v * this.B;  // v × B magnitude
        const em_base = (this.q_charge * cross_vB) / this.proton_mass;
        const corr_UA = 1 + (this.rho_vac_UA / this.rho_vac_SCm);
        const term4 = (em_base * corr_UA) * this.scale_EM;

        // Term 5: Quantum uncertainty
        const sqrt_unc = Math.sqrt(this.delta_x * this.delta_p);
        const term_q = (this.hbar / sqrt_unc) * this.integral_psi * (2 * Math.PI / this.t_Hubble);

        // Term 6: Fluid dynamics
        const V = this.compute_V();
        const term_fluid = (this.rho_fluid * V * this.ug1_base) / this.M;

        // Term 7: Oscillatory waves (real parts)
        const term_osc1 = 2 * this.A_osc * Math.cos(this.k_osc * this.x_pos) * Math.cos(this.omega_osc * t);
        const arg = this.k_osc * this.x_pos - this.omega_osc * t;
        const term_osc2 = (2 * Math.PI / this.t_Hubble_gyr) * this.A_osc * Math.cos(arg);
        const term_osc = term_osc1 + term_osc2;

        // Term 8: DM and density perturbations
        const M_dm = this.M * this.M_DM_factor;
        const pert1 = this.delta_rho_over_rho;
        const pert2 = 3 * this.G * this.M / (this.r * this.r * this.r);
        const term_dm_force_like = (this.M + M_dm) * (pert1 + pert2);
        const term_DM = term_dm_force_like / this.M;

        // Term 9: Stellar wind feedback (pressure/density for acceleration)
        const wind_pressure = this.rho_wind * this.v_wind * this.v_wind;
        const term_wind = wind_pressure / this.rho_fluid;

        // Total acceleration with all 9 terms
        return term1 + term2 + term3 + term4 + term_q + term_fluid + term_osc + term_DM + term_wind;
    }

    // ========== VARIABLE MANAGEMENT (5 methods) ==========

    setVariable(varName, newValue) {
        const validVars = [
            'G', 'M', 'r', 'H0', 'B', 'B_crit', 'Lambda', 'c_light', 'q_charge',
            'gas_v', 'f_TRZ', 'E_0', 'tau_exp', 'rho_wind', 'v_wind', 'rho_fluid',
            'rho_vac_UA', 'rho_vac_SCm', 'scale_EM', 'proton_mass', 'hbar', 't_Hubble',
            't_Hubble_gyr', 'delta_x', 'delta_p', 'integral_psi', 'A_osc', 'k_osc',
            'omega_osc', 'x_pos', 'M_DM_factor', 'delta_rho_over_rho'
        ];

        if (!validVars.includes(varName)) {
            console.error(`Error: Unknown variable '${varName}'.`);
            return false;
        }

        this[varName] = newValue;
        this.updateCache();
        return true;
    }

    getVariable(varName) {
        if (this.hasOwnProperty(varName)) {
            return this[varName];
        }
        console.error(`Error: Unknown variable '${varName}'.`);
        return 0.0;
    }

    addToVariable(varName, delta) {
        return this.setVariable(varName, this.getVariable(varName) + delta);
    }

    subtractFromVariable(varName, delta) {
        return this.addToVariable(varName, -delta);
    }

    listVariables() {
        return [
            'G', 'M', 'r', 'H0', 'B', 'B_crit', 'Lambda', 'c_light', 'q_charge',
            'gas_v', 'f_TRZ', 'E_0', 'tau_exp', 'rho_wind', 'v_wind', 'rho_fluid',
            'rho_vac_UA', 'rho_vac_SCm', 'scale_EM', 'proton_mass', 'hbar', 't_Hubble',
            't_Hubble_gyr', 'delta_x', 'delta_p', 'integral_psi', 'A_osc', 'k_osc',
            'omega_osc', 'x_pos', 'M_DM_factor', 'delta_rho_over_rho'
        ];
    }

    getSystemName() {
        return "BubbleNebula";
    }

    // ========== BATCH OPERATIONS (2 methods) ==========

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

    // ========== SELF-EXPANSION (4 methods) ==========

    expandParameterSpace(factor) {
        const expandable = ['M', 'r', 'B', 'rho_fluid', 'rho_wind', 'A_osc', 'M_DM_factor'];
        return this.scaleVariableGroup(expandable, factor);
    }

    /**
     * Scale nebula mass and radius
     */
    expandNebulaScale(M_scale, r_scale) {
        this.setVariable('M', this.getVariable('M') * M_scale);
        this.setVariable('r', this.getVariable('r') * r_scale);
    }

    /**
     * Scale expansion parameters (UNIQUE to Bubble Nebula)
     * E_0_scale: scales initial expansion factor
     * tau_exp_scale: scales expansion timescale
     */
    expandExpansionScale(E_0_scale, tau_exp_scale) {
        this.setVariable('E_0', this.getVariable('E_0') * E_0_scale);
        this.setVariable('tau_exp', this.getVariable('tau_exp') * tau_exp_scale);
    }

    /**
     * Scale stellar wind and magnetic field parameters
     */
    expandWindMagneticScale(rho_wind_scale, v_wind_scale, B_scale) {
        this.setVariable('rho_wind', this.getVariable('rho_wind') * rho_wind_scale);
        this.setVariable('v_wind', this.getVariable('v_wind') * v_wind_scale);
        this.setVariable('B', this.getVariable('B') * B_scale);
        this.setVariable('B_crit', this.getVariable('B_crit') * B_scale);
    }

    // ========== SELF-REFINEMENT (3 methods) ==========

    autoRefineParameters(observations) {
        if (observations.length === 0) return;

        let sum_error = 0.0;
        for (const obs of observations) {
            const [t, g_obs] = obs;
            const g_calc = this.compute_g_Bubble(t);
            sum_error += Math.abs(g_calc - g_obs);
        }
        const avg_error = sum_error / observations.length;

        if (avg_error > 1e-6) {
            const adj_factor = 1.0 - Math.min(0.1, avg_error / 1e6);
            this.setVariable('E_0', this.getVariable('E_0') * adj_factor);
            this.setVariable('tau_exp', this.getVariable('tau_exp') * (2.0 - adj_factor));
        }
    }

    calibrateToObservations(times, g_obs) {
        if (times.length !== g_obs.length || times.length === 0) return;

        const obs = times.map((t, i) => [t, g_obs[i]]);
        
        for (let iter = 0; iter < 5; iter++) {
            this.autoRefineParameters(obs);
        }
    }

    optimizeForMetric(metric, t_start, t_end, steps) {
        let best_score = -1e100;
        const dt = (t_end - t_start) / steps;

        for (let i = 0; i <= steps; i++) {
            const t = t_start + i * dt;
            const g = this.compute_g_Bubble(t);
            const score = metric(g);
            if (score > best_score) best_score = score;
        }
        return best_score;
    }

    // ========== PARAMETER EXPLORATION (1 method) ==========

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

    // ========== ADAPTIVE EVOLUTION (2 methods) ==========

    mutateParameters(mutation_rate) {
        const vars = this.listVariables();
        const constants = ['c_light', 'G', 'hbar'];

        for (const v of vars) {
            if (constants.includes(v)) continue;
            const val = this.getVariable(v);
            const delta = val * (Math.random() * 2 - 1) * mutation_rate;
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

    // ========== STATE MANAGEMENT (4 methods) ==========

    saveState(stateName) {
        const state = {};
        const vars = this.listVariables();
        for (const v of vars) {
            state[v] = this.getVariable(v);
        }
        BubbleNebula.savedStates[stateName] = state;
        return true;
    }

    restoreState(stateName) {
        const state = BubbleNebula.savedStates[stateName];
        if (!state) return false;

        for (const [varName, value] of Object.entries(state)) {
            this.setVariable(varName, value);
        }
        return true;
    }

    listSavedStates() {
        return Object.keys(BubbleNebula.savedStates);
    }

    exportState() {
        let output = 'BubbleNebula State Export:\n';
        const vars = this.listVariables();
        for (const v of vars) {
            output += `${v} = ${this.getVariable(v).toExponential(6)}\n`;
        }
        return output;
    }

    // ========== SYSTEM ANALYSIS (4 methods) ==========

    sensitivityAnalysis(t, delta_pct) {
        const sensitivities = {};
        const g_base = this.compute_g_Bubble(t);
        const vars = this.listVariables();
        const constants = ['c_light', 'G', 'hbar'];

        for (const v of vars) {
            if (constants.includes(v)) continue;

            const original = this.getVariable(v);
            const delta = original * delta_pct / 100.0;

            this.setVariable(v, original + delta);
            const g_plus = this.compute_g_Bubble(t);
            this.setVariable(v, original);

            const sensitivity = (g_base !== 0.0) ? Math.abs((g_plus - g_base) / g_base) : 0.0;
            sensitivities[v] = sensitivity;
        }

        return sensitivities;
    }

    generateReport(t) {
        const M_sun = 1.989e30;
        const ly_to_m = 9.461e15;
        const Et = this.E_t(t);
        const g = this.compute_g_Bubble(t);

        let report = '============================================\n';
        report += 'BUBBLE NEBULA REPORT\n';
        report += 'NGC 7635 Emission Nebula\n';
        report += '============================================\n';
        report += `Time: t = ${t.toExponential(6)} s (${(t/3.156e7/1e6).toFixed(3)} Myr)\n\n`;

        report += 'Physical Parameters:\n';
        report += `  Mass M = ${this.M.toExponential(6)} kg (${(this.M/M_sun).toFixed(2)} M_sun)\n`;
        report += `  Radius r = ${this.r.toExponential(6)} m (${(this.r/ly_to_m).toFixed(2)} ly)\n`;
        report += `  Hubble constant H0 = ${this.H0.toExponential(6)} s^-1\n`;
        report += `  Magnetic field B = ${this.B.toExponential(6)} T (B_crit = ${this.B_crit.toExponential(3)} T)\n`;
        report += `  f_TRZ = ${this.f_TRZ.toFixed(6)}\n`;
        report += `  Expansion E_0 = ${this.E_0.toFixed(6)}, tau_exp = ${this.tau_exp.toExponential(6)} s\n`;
        report += `  E(t) = ${Et.toFixed(6)} (nebula expansion factor)\n`;
        report += `  Wind density rho_wind = ${this.rho_wind.toExponential(6)} kg/m^3, v_wind = ${this.v_wind.toExponential(6)} m/s\n`;
        report += `  Fluid density rho_fluid = ${this.rho_fluid.toExponential(6)} kg/m^3\n`;
        report += `  Gas velocity = ${this.gas_v.toExponential(6)} m/s\n`;
        report += `  DM factor = ${this.M_DM_factor.toFixed(6)}\n\n`;

        report += 'Computed Acceleration:\n';
        report += `  g_Bubble(t) = ${g.toExponential(6)} m/s^2\n\n`;

        // UQFF term breakdown
        const corr_H = 1 + this.H0 * t;
        const corr_B = 1 - this.B / this.B_crit;
        const corr_E = 1 - Et;
        report += 'UQFF Terms:\n';
        report += `  Base (with H0, B, E(t)): ${(this.ug1_base * corr_H * corr_B * corr_E).toExponential(6)} m/s^2\n`;
        report += `  Ug total: ${this.compute_Ug(Et).toExponential(6)} m/s^2\n`;
        report += `  Lambda: ${((this.Lambda * this.c_light * this.c_light) / 3.0).toExponential(6)} m/s^2\n`;

        const cross_vB = this.gas_v * this.B;
        const em_base = (this.q_charge * cross_vB) / this.proton_mass;
        const corr_UA = 1 + (this.rho_vac_UA / this.rho_vac_SCm);
        report += `  EM (scaled with UA): ${(em_base * corr_UA * this.scale_EM).toExponential(6)} m/s^2\n`;

        const sqrt_unc = Math.sqrt(this.delta_x * this.delta_p);
        report += `  Quantum: ${((this.hbar / sqrt_unc) * this.integral_psi * (2 * Math.PI / this.t_Hubble)).toExponential(6)} m/s^2\n`;

        const V = this.compute_V();
        report += `  Fluid: ${((this.rho_fluid * V * this.ug1_base) / this.M).toExponential(6)} m/s^2\n`;
        report += `  Oscillatory: (combined real parts)\n`;

        const M_dm = this.M * this.M_DM_factor;
        const pert1 = this.delta_rho_over_rho;
        const pert2 = 3 * this.G * this.M / (this.r * this.r * this.r);
        const term_dm_force_like = (this.M + M_dm) * (pert1 + pert2);
        report += `  DM: ${(term_dm_force_like / this.M).toExponential(6)} m/s^2\n`;

        const wind_pressure = this.rho_wind * this.v_wind * this.v_wind;
        report += `  Stellar Wind Feedback: ${(wind_pressure / this.rho_fluid).toExponential(6)} m/s^2\n`;

        report += `\n  Expansion Correction: (1 - E(t)) = ${corr_E.toFixed(6)} (UNIQUE)\n`;

        report += '============================================\n';
        return report;
    }

    validateConsistency() {
        let valid = true;

        if (this.M <= 0 || this.r <= 0) {
            console.error('Error: M and r must be positive.');
            valid = false;
        }
        if (this.B < 0 || this.B_crit <= 0) {
            console.error('Error: B, B_crit must be non-negative/positive.');
            valid = false;
        }
        if (this.H0 < 0) {
            console.error('Error: H0 must be non-negative.');
            valid = false;
        }
        if (this.tau_exp <= 0) {
            console.error('Error: Expansion timescale must be positive.');
            valid = false;
        }
        if (this.E_0 < 0 || this.E_0 > 1.0) {
            console.warn('Warning: Expansion factor E_0 outside [0,1].');
        }
        if (this.rho_fluid <= 0 || this.rho_wind < 0) {
            console.error('Error: Fluid/wind densities must be positive/non-negative.');
            valid = false;
        }
        if (this.v_wind < 0 || this.gas_v < 0) {
            console.error('Error: Velocities must be non-negative.');
            valid = false;
        }
        if (this.M_DM_factor < 0 || this.M_DM_factor > 1.0) {
            console.warn('Warning: DM factor outside [0,1].');
        }

        return valid;
    }

    autoCorrectAnomalies() {
        let corrected = false;
        const M_sun = 1.989e30;
        const ly_to_m = 9.461e15;

        if (this.M <= 0) { this.M = 46.0 * M_sun; corrected = true; }
        if (this.r <= 0) { this.r = 5.0 * ly_to_m; corrected = true; }
        if (this.B < 0) { this.B = 1e-6; corrected = true; }
        if (this.B_crit <= 0) { this.B_crit = 1e11; corrected = true; }
        if (this.H0 < 0) { this.H0 = 2.184e-18; corrected = true; }
        if (this.tau_exp <= 0) { this.tau_exp = 4e6 * 3.156e7; corrected = true; }
        if (this.E_0 < 0) { this.E_0 = 0.0; corrected = true; }
        if (this.E_0 > 1.0) { this.E_0 = 1.0; corrected = true; }
        if (this.rho_fluid <= 0) { this.rho_fluid = 1e-21; corrected = true; }
        if (this.rho_wind < 0) { this.rho_wind = 1e-21; corrected = true; }
        if (this.v_wind < 0) { this.v_wind = 1.8e6; corrected = true; }
        if (this.gas_v < 0) { this.gas_v = 1e5; corrected = true; }
        if (this.M_DM_factor < 0) { this.M_DM_factor = 0.1; corrected = true; }
        if (this.M_DM_factor > 1.0) { this.M_DM_factor = 1.0; corrected = true; }

        if (corrected) this.updateCache();
        return corrected;
    }

    printParameters() {
        console.log('Bubble Nebula Parameters:');
        console.log(`G: ${this.G}, M: ${this.M}, r: ${this.r}`);
        console.log(`H0: ${this.H0}, B: ${this.B}, B_crit: ${this.B_crit}`);
        console.log(`f_TRZ: ${this.f_TRZ}, E_0: ${this.E_0}, tau_exp: ${this.tau_exp}`);
        console.log(`rho_fluid: ${this.rho_fluid}, rho_wind: ${this.rho_wind}, v_wind: ${this.v_wind}`);
        console.log(`gas_v: ${this.gas_v}, M_DM_factor: ${this.M_DM_factor}`);
        console.log(`A_osc: ${this.A_osc}, delta_rho_over_rho: ${this.delta_rho_over_rho}`);
        console.log(`ug1_base: ${this.ug1_base}`);
    }

    exampleAt2Myr() {
        const t_example = 2e6 * 3.156e7; // 2 Myr in seconds
        return this.compute_g_Bubble(t_example);
    }
}

// Static property for saved states (replaces C++ anonymous namespace)
BubbleNebula.savedStates = {};

// ========== INLINE TEST ==========
if (typeof module !== 'undefined' && require.main === module) {
    console.log('='.repeat(60));
    console.log('BUBBLE NEBULA INLINE TEST');
    console.log('NGC 7635 Emission Nebula with Expansion Dynamics');
    console.log('='.repeat(60));

    const bubble = new BubbleNebula();
    const M_sun = 1.989e30;
    const ly_to_m = 9.461e15;

    console.log(`\nInitial Parameters:`);
    console.log(`  M = ${(bubble.M/M_sun).toFixed(1)} M_sun (static central star mass)`);
    console.log(`  r = ${(bubble.r/ly_to_m).toFixed(1)} ly (nebula radius)`);
    console.log(`  E_0 = ${bubble.E_0.toFixed(3)} (expansion factor)`);
    console.log(`  tau_exp = ${(bubble.tau_exp/3.156e7/1e6).toFixed(1)} Myr`);
    console.log(`  v_wind = ${(bubble.v_wind/1e3).toFixed(0)} km/s`);

    // Test at 2 Myr
    const t_test = 2e6 * 3.156e7;
    const g = bubble.compute_g_Bubble(t_test);
    const Et = bubble.E_t(t_test);

    console.log(`\nTest at t = 2 Myr:`);
    console.log(`  g = ${g.toExponential(6)} m/s^2`);
    console.log(`  E(t) = ${Et.toFixed(6)} (expansion factor)`);
    console.log(`  (1 - E(t)) = ${(1-Et).toFixed(6)} (gravity reduction factor)`);

    // Test expansion unique feature
    console.log(`\nNEBULA EXPANSION Evolution (UNIQUE to Bubble Nebula):`);
    const t_values_Myr = [0.0, 1.0, 2.0, 4.0, 8.0];
    for (const t_Myr of t_values_Myr) {
        const t = t_Myr * 1e6 * 3.156e7;
        const E = bubble.E_t(t);
        const reduction = 1 - E;
        console.log(`  t = ${t_Myr.toFixed(1)} Myr: E(t) = ${E.toFixed(6)}, gravity reduction = ${reduction.toFixed(6)}`);
    }

    // Test gravity reduction
    const g_t0 = bubble.compute_g_Bubble(0);
    const g_t2 = bubble.compute_g_Bubble(2e6 * 3.156e7);
    console.log(`\nGravity Reduction Due to Expansion:`);
    console.log(`  g(t=0) = ${g_t0.toExponential(6)} m/s^2`);
    console.log(`  g(t=2 Myr) = ${g_t2.toExponential(6)} m/s^2`);
    console.log(`  Reduction ratio: ${(g_t2/g_t0).toFixed(6)}`);

    console.log('\n' + '='.repeat(60));
    console.log('Bubble Nebula Inline Test Complete ✓');
    console.log('='.repeat(60));
}

// ES6 Module Export
export default BubbleNebula;
