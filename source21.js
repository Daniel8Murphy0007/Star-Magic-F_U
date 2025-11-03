/**
 * ================================================================================================
 * Module: source21.js (NGC 3603 Extreme Star Cluster)
 *
 * Description: JavaScript ES6 Module for NGC 3603 Young Massive Star Cluster Class
 *              Converted from C++ source21.cpp - Module 11 in the UQFF series
 *              Focuses on extreme star-forming region with cavity pressure dynamics
 *
 * Purpose: Implements Master Universal Gravity Equation (MUGE) for NGC 3603 evolution
 *          Including ALL terms: base gravity with mass growth M(t), cosmic expansion (H_0),
 *          magnetic correction, UQFF Ug components, Lambda, quantum uncertainty, scaled EM,
 *          fluid dynamics, oscillatory waves, DM/density perturbations, stellar wind feedback,
 *          and CAVITY PRESSURE P(t)/rho_fluid (unique feature)
 *
 * Key Features:
 *   - M0 = 400,000 M☉ (extremely massive young cluster)
 *   - r = 9.5 ly cluster radius
 *   - tau_SF = 1 Myr star formation timescale
 *   - P(t) = P0 × exp(-t/tau_exp) - Cavity pressure expansion (UNIQUE)
 *   - rho_wind = 1e-20 kg/m³, v_wind = 2e6 m/s (extreme stellar winds)
 *   - M(t) = M0 × (1 + M_dot_factor × exp(-t/tau_SF)) - Mass growth
 *   - 34 physics parameters + 38 methods (13 core + 25 enhanced)
 *
 * Author: Converted by GitHub Copilot from Daniel T. Murphy's UQFF framework
 * Date: November 03, 2025
 * Copyright: Daniel T. Murphy, daniel.murphy00@gmail.com
 * ================================================================================================
 */

class NGC3603 {
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

        // NGC 3603 specific parameters
        const M_sun = 1.989e30;        // Solar mass (kg)
        const ly_to_m = 9.461e15;      // Light-year to meters
        this.M0 = 400000.0 * M_sun;    // Initial cluster mass: 400,000 M☉
        this.r = 9.5 * ly_to_m;        // Cluster radius: 9.5 light-years

        // Cosmological parameters
        this.H0 = 2.184e-18;           // Hubble constant (s⁻¹)
        this.Lambda = 1.1e-52;         // Cosmological constant (m⁻²)
        this.t_Hubble = 13.8e9 * 3.156e7; // Hubble time (s)
        this.t_Hubble_gyr = 13.8;      // Hubble time (Gyr)

        // Magnetic field parameters
        this.B = 1e-5;                 // Static magnetic field (T)
        this.B_crit = 1e11;            // Critical magnetic field (T)

        // Star formation parameters
        this.M_dot_factor = 1.0;       // Star formation rate factor (dimensionless)
        this.tau_SF = 1e6 * 3.156e7;   // Star formation timescale: 1 Myr (s)

        // Stellar wind parameters (EXTREME winds from massive stars)
        this.rho_wind = 1e-20;         // Wind density (kg/m³)
        this.v_wind = 2e6;             // Wind velocity: 2000 km/s (m/s)

        // Fluid/gas parameters
        this.rho_fluid = 1e-20;        // Fluid density (kg/m³)
        this.gas_v = 1e5;              // Gas velocity for EM (m/s)

        // CAVITY PRESSURE parameters (UNIQUE to NGC 3603)
        this.P0 = 4e-8;                // Initial cavity pressure (Pa)
        this.tau_exp = 1e6 * 3.156e7;  // Cavity expansion timescale: 1 Myr (s)

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
        // Cache base Ug1 for efficiency
        this.ug1_base = (this.G * this.M0) / (this.r * this.r);
    }

    // ========== CORE PHYSICS METHODS (13 methods) ==========

    /**
     * Compute time-dependent mass M(t) with exponential star formation
     * M(t) = M0 × (1 + M_dot_factor × exp(-t/tau_SF))
     */
    M_t(t) {
        const M_dot = this.M_dot_factor * Math.exp(-t / this.tau_SF);
        return this.M0 * (1 + M_dot);
    }

    /**
     * Compute time-dependent cavity pressure P(t) with exponential expansion
     * P(t) = P0 × exp(-t/tau_exp)
     * UNIQUE to NGC 3603 - represents expanding cavity from stellar winds
     */
    P_t(t) {
        return this.P0 * Math.exp(-t / this.tau_exp);
    }

    /**
     * Compute UQFF Ug terms (Ug1 + Ug2 + Ug3 + Ug4) with f_TRZ factor
     */
    compute_Ug(Mt) {
        const Ug1 = (this.G * Mt) / (this.r * this.r);
        const Ug2 = 0.0;
        const Ug3 = 0.0;
        const corr_B = 1 - this.B / this.B_crit;
        const Ug4 = Ug1 * corr_B;
        return (Ug1 + Ug2 + Ug3 + Ug4) * (1 + this.f_TRZ);
    }

    /**
     * Compute cluster volume for fluid dynamics
     */
    compute_V() {
        return (4.0 / 3.0) * Math.PI * this.r * this.r * this.r;
    }

    /**
     * Main MUGE computation: g_NGC3603(t) with ALL 10 terms
     * Includes cavity pressure term P(t)/rho_fluid (UNIQUE)
     */
    compute_g_NGC3603(t) {
        if (t < 0) {
            console.error('Error: Time t must be non-negative.');
            return 0.0;
        }

        const Mt = this.M_t(t);
        const Pt = this.P_t(t);
        const ug1_t = (this.G * Mt) / (this.r * this.r);

        // Term 1: Base gravity with H0 and B corrections
        const corr_H = 1 + this.H0 * t;
        const corr_B = 1 - this.B / this.B_crit;
        const term1 = ug1_t * corr_H * corr_B;

        // Term 2: UQFF Ug with f_TRZ
        const term2 = this.compute_Ug(Mt);

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
        const term_fluid = (this.rho_fluid * V * ug1_t) / Mt;

        // Term 7: Oscillatory waves (real parts)
        const term_osc1 = 2 * this.A_osc * Math.cos(this.k_osc * this.x_pos) * Math.cos(this.omega_osc * t);
        const arg = this.k_osc * this.x_pos - this.omega_osc * t;
        const term_osc2 = (2 * Math.PI / this.t_Hubble_gyr) * this.A_osc * Math.cos(arg);
        const term_osc = term_osc1 + term_osc2;

        // Term 8: DM and density perturbations
        const M_dm = Mt * this.M_DM_factor;
        const pert1 = this.delta_rho_over_rho;
        const pert2 = 3 * this.G * Mt / (this.r * this.r * this.r);
        const term_dm_force_like = (Mt + M_dm) * (pert1 + pert2);
        const term_DM = term_dm_force_like / Mt;

        // Term 9: Stellar wind feedback (pressure/density for acceleration)
        const wind_pressure = this.rho_wind * this.v_wind * this.v_wind;
        const term_wind = wind_pressure / this.rho_fluid;

        // Term 10: CAVITY PRESSURE (UNIQUE to NGC 3603)
        // P(t)/rho_fluid provides acceleration from expanding cavity
        const term_pressure = Pt / this.rho_fluid;

        // Total acceleration with all 10 terms
        return term1 + term2 + term3 + term4 + term_q + term_fluid + term_osc + term_DM + term_wind + term_pressure;
    }

    // ========== VARIABLE MANAGEMENT (5 methods) ==========

    setVariable(varName, newValue) {
        const validVars = [
            'G', 'M0', 'r', 'H0', 'B', 'B_crit', 'Lambda', 'c_light', 'q_charge',
            'gas_v', 'f_TRZ', 'M_dot_factor', 'tau_SF', 'rho_wind', 'v_wind',
            'rho_fluid', 'P0', 'tau_exp', 'rho_vac_UA', 'rho_vac_SCm', 'scale_EM',
            'proton_mass', 'hbar', 't_Hubble', 't_Hubble_gyr', 'delta_x', 'delta_p',
            'integral_psi', 'A_osc', 'k_osc', 'omega_osc', 'x_pos', 'M_DM_factor',
            'delta_rho_over_rho'
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
            'G', 'M0', 'r', 'H0', 'B', 'B_crit', 'Lambda', 'c_light', 'q_charge',
            'gas_v', 'f_TRZ', 'M_dot_factor', 'tau_SF', 'rho_wind', 'v_wind',
            'rho_fluid', 'P0', 'tau_exp', 'rho_vac_UA', 'rho_vac_SCm', 'scale_EM',
            'proton_mass', 'hbar', 't_Hubble', 't_Hubble_gyr', 'delta_x', 'delta_p',
            'integral_psi', 'A_osc', 'k_osc', 'omega_osc', 'x_pos', 'M_DM_factor',
            'delta_rho_over_rho'
        ];
    }

    getSystemName() {
        return "NGC3603";
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
        const expandable = ['M0', 'r', 'B', 'rho_fluid', 'rho_wind', 'A_osc', 'M_DM_factor'];
        return this.scaleVariableGroup(expandable, factor);
    }

    /**
     * Scale star formation parameters
     */
    expandStarFormationScale(M_dot_factor_scale, tau_SF_scale) {
        this.setVariable('M_dot_factor', this.getVariable('M_dot_factor') * M_dot_factor_scale);
        this.setVariable('tau_SF', this.getVariable('tau_SF') * tau_SF_scale);
    }

    /**
     * Scale cavity pressure parameters (UNIQUE to NGC 3603)
     * P0_scale: scales initial cavity pressure
     * tau_exp_scale: scales cavity expansion timescale
     */
    expandCavityPressureScale(P0_scale, tau_exp_scale) {
        this.setVariable('P0', this.getVariable('P0') * P0_scale);
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
            const g_calc = this.compute_g_NGC3603(t);
            sum_error += Math.abs(g_calc - g_obs);
        }
        const avg_error = sum_error / observations.length;

        if (avg_error > 1e-6) {
            const adj_factor = 1.0 - Math.min(0.1, avg_error / 1e6);
            this.setVariable('M_dot_factor', this.getVariable('M_dot_factor') * adj_factor);
            this.setVariable('tau_SF', this.getVariable('tau_SF') * (2.0 - adj_factor));
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
            const g = this.compute_g_NGC3603(t);
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
        NGC3603.savedStates[stateName] = state;
        return true;
    }

    restoreState(stateName) {
        const state = NGC3603.savedStates[stateName];
        if (!state) return false;

        for (const [varName, value] of Object.entries(state)) {
            this.setVariable(varName, value);
        }
        return true;
    }

    listSavedStates() {
        return Object.keys(NGC3603.savedStates);
    }

    exportState() {
        let output = 'NGC3603 State Export:\n';
        const vars = this.listVariables();
        for (const v of vars) {
            output += `${v} = ${this.getVariable(v).toExponential(6)}\n`;
        }
        return output;
    }

    // ========== SYSTEM ANALYSIS (4 methods) ==========

    sensitivityAnalysis(t, delta_pct) {
        const sensitivities = {};
        const g_base = this.compute_g_NGC3603(t);
        const vars = this.listVariables();
        const constants = ['c_light', 'G', 'hbar'];

        for (const v of vars) {
            if (constants.includes(v)) continue;

            const original = this.getVariable(v);
            const delta = original * delta_pct / 100.0;

            this.setVariable(v, original + delta);
            const g_plus = this.compute_g_NGC3603(t);
            this.setVariable(v, original);

            const sensitivity = (g_base !== 0.0) ? Math.abs((g_plus - g_base) / g_base) : 0.0;
            sensitivities[v] = sensitivity;
        }

        return sensitivities;
    }

    generateReport(t) {
        const M_sun = 1.989e30;
        const ly_to_m = 9.461e15;
        const Mt = this.M_t(t);
        const Pt = this.P_t(t);
        const g = this.compute_g_NGC3603(t);

        let report = '============================================\n';
        report += 'NGC 3603 REPORT\n';
        report += 'Extreme Young Massive Star Cluster\n';
        report += '============================================\n';
        report += `Time: t = ${t.toExponential(6)} s (${(t/3.156e7/1e6).toFixed(3)} Myr)\n\n`;

        report += 'Physical Parameters:\n';
        report += `  Initial Mass M0 = ${this.M0.toExponential(6)} kg (${(this.M0/M_sun).toExponential(3)} M_sun)\n`;
        report += `  M(t) = ${Mt.toExponential(6)} kg (${(Mt/M_sun).toExponential(3)} M_sun)\n`;
        report += `  Radius r = ${this.r.toExponential(6)} m (${(this.r/ly_to_m).toFixed(2)} ly)\n`;
        report += `  Hubble constant H0 = ${this.H0.toExponential(6)} s^-1\n`;
        report += `  Magnetic field B = ${this.B.toExponential(6)} T (B_crit = ${this.B_crit.toExponential(3)} T)\n`;
        report += `  f_TRZ = ${this.f_TRZ.toFixed(6)}\n`;
        report += `  Star formation M_dot_factor = ${this.M_dot_factor.toFixed(6)}, tau_SF = ${this.tau_SF.toExponential(6)} s\n`;
        report += `  Cavity pressure P0 = ${this.P0.toExponential(6)} Pa, tau_exp = ${this.tau_exp.toExponential(6)} s\n`;
        report += `  P(t) = ${Pt.toExponential(6)} Pa\n`;
        report += `  Wind density rho_wind = ${this.rho_wind.toExponential(6)} kg/m^3, v_wind = ${this.v_wind.toExponential(6)} m/s\n`;
        report += `  Fluid density rho_fluid = ${this.rho_fluid.toExponential(6)} kg/m^3\n`;
        report += `  Gas velocity = ${this.gas_v.toExponential(6)} m/s\n`;
        report += `  DM factor = ${this.M_DM_factor.toFixed(6)}\n\n`;

        report += 'Computed Acceleration:\n';
        report += `  g_NGC3603(t) = ${g.toExponential(6)} m/s^2\n\n`;

        // UQFF term breakdown
        const ug1_t = (this.G * Mt) / (this.r * this.r);
        const corr_H = 1 + this.H0 * t;
        const corr_B = 1 - this.B / this.B_crit;
        report += 'UQFF Terms:\n';
        report += `  Base (with H0, B, M(t)): ${(ug1_t * corr_H * corr_B).toExponential(6)} m/s^2\n`;
        report += `  Ug total: ${this.compute_Ug(Mt).toExponential(6)} m/s^2\n`;
        report += `  Lambda: ${((this.Lambda * this.c_light * this.c_light) / 3.0).toExponential(6)} m/s^2\n`;

        const cross_vB = this.gas_v * this.B;
        const em_base = (this.q_charge * cross_vB) / this.proton_mass;
        const corr_UA = 1 + (this.rho_vac_UA / this.rho_vac_SCm);
        report += `  EM (scaled with UA): ${(em_base * corr_UA * this.scale_EM).toExponential(6)} m/s^2\n`;

        const sqrt_unc = Math.sqrt(this.delta_x * this.delta_p);
        report += `  Quantum: ${((this.hbar / sqrt_unc) * this.integral_psi * (2 * Math.PI / this.t_Hubble)).toExponential(6)} m/s^2\n`;

        const V = this.compute_V();
        report += `  Fluid: ${((this.rho_fluid * V * ug1_t) / Mt).toExponential(6)} m/s^2\n`;
        report += `  Oscillatory: (combined real parts)\n`;

        const M_dm = Mt * this.M_DM_factor;
        const pert1 = this.delta_rho_over_rho;
        const pert2 = 3 * this.G * Mt / (this.r * this.r * this.r);
        const term_dm_force_like = (Mt + M_dm) * (pert1 + pert2);
        report += `  DM: ${(term_dm_force_like / Mt).toExponential(6)} m/s^2\n`;

        const wind_pressure = this.rho_wind * this.v_wind * this.v_wind;
        report += `  Stellar Wind Feedback: ${(wind_pressure / this.rho_fluid).toExponential(6)} m/s^2\n`;

        report += `  Cavity Pressure: ${(Pt / this.rho_fluid).toExponential(6)} m/s^2 (UNIQUE)\n`;

        report += '============================================\n';
        return report;
    }

    validateConsistency() {
        let valid = true;

        if (this.M0 <= 0 || this.r <= 0) {
            console.error('Error: M0 and r must be positive.');
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
        if (this.tau_SF <= 0 || this.tau_exp <= 0) {
            console.error('Error: Timescales must be positive.');
            valid = false;
        }
        if (this.P0 < 0) {
            console.warn('Warning: Initial pressure P0 is negative.');
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

        if (this.M0 <= 0) { this.M0 = 400000.0 * M_sun; corrected = true; }
        if (this.r <= 0) { this.r = 9.5 * ly_to_m; corrected = true; }
        if (this.B < 0) { this.B = 1e-5; corrected = true; }
        if (this.B_crit <= 0) { this.B_crit = 1e11; corrected = true; }
        if (this.H0 < 0) { this.H0 = 2.184e-18; corrected = true; }
        if (this.tau_SF <= 0) { this.tau_SF = 1e6 * 3.156e7; corrected = true; }
        if (this.tau_exp <= 0) { this.tau_exp = 1e6 * 3.156e7; corrected = true; }
        if (this.P0 < 0) { this.P0 = 4e-8; corrected = true; }
        if (this.rho_fluid <= 0) { this.rho_fluid = 1e-20; corrected = true; }
        if (this.rho_wind < 0) { this.rho_wind = 1e-20; corrected = true; }
        if (this.v_wind < 0) { this.v_wind = 2e6; corrected = true; }
        if (this.gas_v < 0) { this.gas_v = 1e5; corrected = true; }
        if (this.M_DM_factor < 0) { this.M_DM_factor = 0.1; corrected = true; }
        if (this.M_DM_factor > 1.0) { this.M_DM_factor = 1.0; corrected = true; }

        if (corrected) this.updateCache();
        return corrected;
    }

    printParameters() {
        console.log('NGC 3603 Parameters:');
        console.log(`G: ${this.G}, M0: ${this.M0}, r: ${this.r}`);
        console.log(`H0: ${this.H0}, B: ${this.B}, B_crit: ${this.B_crit}`);
        console.log(`f_TRZ: ${this.f_TRZ}, M_dot_factor: ${this.M_dot_factor}, tau_SF: ${this.tau_SF}`);
        console.log(`rho_fluid: ${this.rho_fluid}, rho_wind: ${this.rho_wind}, v_wind: ${this.v_wind}`);
        console.log(`P0: ${this.P0}, tau_exp: ${this.tau_exp}`);
        console.log(`gas_v: ${this.gas_v}, M_DM_factor: ${this.M_DM_factor}`);
        console.log(`A_osc: ${this.A_osc}, delta_rho_over_rho: ${this.delta_rho_over_rho}`);
        console.log(`ug1_base: ${this.ug1_base}`);
    }

    exampleAt500kYears() {
        const t_example = 5e5 * 3.156e7; // 500k years in seconds
        return this.compute_g_NGC3603(t_example);
    }
}

// Static property for saved states (replaces C++ anonymous namespace)
NGC3603.savedStates = {};

// ========== INLINE TEST ==========
if (typeof module !== 'undefined' && require.main === module) {
    console.log('='.repeat(60));
    console.log('NGC 3603 INLINE TEST');
    console.log('Extreme Young Massive Star Cluster with Cavity Pressure');
    console.log('='.repeat(60));

    const ngc = new NGC3603();
    const M_sun = 1.989e30;
    const ly_to_m = 9.461e15;

    console.log(`\nInitial Parameters:`);
    console.log(`  M0 = ${(ngc.M0/M_sun).toFixed(0)} M_sun`);
    console.log(`  r = ${(ngc.r/ly_to_m).toFixed(1)} ly`);
    console.log(`  P0 = ${ngc.P0.toExponential(3)} Pa (cavity pressure)`);
    console.log(`  tau_exp = ${(ngc.tau_exp/3.156e7/1e6).toFixed(1)} Myr`);
    console.log(`  v_wind = ${(ngc.v_wind/1e3).toFixed(0)} km/s (extreme stellar winds)`);

    // Test at 500k years
    const t_test = 5e5 * 3.156e7;
    const g = ngc.compute_g_NGC3603(t_test);
    const Mt = ngc.M_t(t_test);
    const Pt = ngc.P_t(t_test);

    console.log(`\nTest at t = 500k years:`);
    console.log(`  g = ${g.toExponential(6)} m/s^2`);
    console.log(`  M(t) = ${(Mt/M_sun).toFixed(0)} M_sun`);
    console.log(`  P(t) = ${Pt.toExponential(3)} Pa (cavity pressure decay)`);
    console.log(`  P(t)/P0 = ${(Pt/ngc.P0).toFixed(6)} (exponential decay)`);

    // Test cavity pressure unique feature
    console.log(`\nCAVITY PRESSURE Evolution (UNIQUE to NGC 3603):`);
    const t_values_Myr = [0.0, 0.5, 1.0, 2.0, 5.0];
    for (const t_Myr of t_values_Myr) {
        const t = t_Myr * 1e6 * 3.156e7;
        const P = ngc.P_t(t);
        console.log(`  t = ${t_Myr.toFixed(1)} Myr: P(t) = ${P.toExponential(3)} Pa, P(t)/P0 = ${(P/ngc.P0).toFixed(6)}`);
    }

    // Test stellar wind contribution
    const wind_pressure = ngc.rho_wind * ngc.v_wind * ngc.v_wind;
    const wind_accel = wind_pressure / ngc.rho_fluid;
    console.log(`\nStellar Wind Feedback (extreme winds):`);
    console.log(`  Wind pressure = ${wind_pressure.toExponential(3)} Pa`);
    console.log(`  Wind acceleration = ${wind_accel.toExponential(3)} m/s^2`);

    console.log('\n' + '='.repeat(60));
    console.log('NGC 3603 Inline Test Complete ✓');
    console.log('='.repeat(60));
}

// ES6 Module Export
export default NGC3603;
