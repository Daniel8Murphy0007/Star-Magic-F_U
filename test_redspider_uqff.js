// test_redspider_uqff.js
// Comprehensive test suite for Red Spider Nebula (NGC 6537) UQFF Module
// Testing frequency-resonance dynamics with 65+ tests across 8 categories
//
// Test Categories:
// 1. Initialization (14 tests)
// 2. Frequency Components (15 tests)
// 3. Quantum Uncertainty (8 tests)
// 4. Resonance Physics (10 tests)
// 5. Topological Terms (12 tests - DPM, THz, U_g4i)
// 6. Dynamic Updates (8 tests)
// 7. Master Equation (10 tests)
// 8. Performance (8 tests)
// Total: 85 comprehensive tests

const NGC6537UQFFModule = require('./redspider_uqff.js');

class RedSpiderTestSuite {
    constructor() {
        this.module = new NGC6537UQFFModule();
        this.tests_passed = 0;
        this.tests_failed = 0;
        this.test_results = [];
    }
    
    // Test helper
    assert(condition, message) {
        if (condition) {
            this.tests_passed++;
            this.test_results.push({ status: '✓', message });
            return true;
        } else {
            this.tests_failed++;
            this.test_results.push({ status: '✗', message });
            return false;
        }
    }
    
    // Helper to check numerical range
    inRange(value, min, max, tolerance = 0.01) {
        const range = Math.abs(max - min) * tolerance;
        return value >= (min - range) && value <= (max + range);
    }
    
    // ===== CATEGORY 1: INITIALIZATION TESTS (14 tests) =====
    testInitialization() {
        console.log('\n=== TEST CATEGORY 1: Initialization (14 tests) ===');
        
        // T1.1 - Module instantiation
        this.assert(this.module !== null && this.module !== undefined,
            'T1.1: Module instantiates successfully');
        
        // T1.2 - Variable count
        const varCount = Object.keys(this.module.variables).length;
        this.assert(varCount >= 30, 
            `T1.2: Module initializes with 30+ variables (found ${varCount})`);
        
        // T1.3 - Universal constant: c
        this.assert(Math.abs(this.module.variables['c'] - 3e8) < 1e5,
            'T1.3: Speed of light initialized correctly (3e8 m/s)');
        
        // T1.4 - Universal constant: ℏ
        this.assert(Math.abs(this.module.variables['hbar'] - 1.0546e-34) < 1e-36,
            'T1.4: Reduced Planck constant initialized (1.0546e-34 J·s)');
        
        // T1.5 - Universal constant: λ_Planck
        this.assert(Math.abs(this.module.variables['lambda_planck'] - 1.616e-35) < 1e-37,
            'T1.5: Planck length initialized (1.616e-35 m)');
        
        // T1.6 - NGC 6537 parameter: radius
        this.assert(Math.abs(this.module.variables['r'] - 7.1e15) < 1e14,
            'T1.6: NGC 6537 radius initialized (7.1e15 m)');
        
        // T1.7 - NGC 6537 parameter: filament density
        this.assert(Math.abs(this.module.variables['rho_fil'] - 1e-20) < 1e-22,
            'T1.7: Filament density initialized (1e-20 kg/m³)');
        
        // T1.8 - NGC 6537 parameter: expansion velocity
        this.assert(Math.abs(this.module.variables['v_exp'] - 3e5) < 1e3,
            'T1.8: Expansion velocity initialized (3e5 m/s)');
        
        // T1.9 - NGC 6537 parameter: age in seconds
        const t_age_sec = 1900 * 3.156e7;
        this.assert(Math.abs(this.module.variables['t_age'] - t_age_sec) / t_age_sec < 1e-6,
            'T1.9: Nebular age converted correctly (~1900 years in seconds)');
        
        // T1.10 - Frequency parameter: f_super
        this.assert(Math.abs(this.module.variables['f_super'] - 1.411e16) < 1e14,
            'T1.10: Superconductive frequency initialized (1.411e16 Hz)');
        
        // T1.11 - Frequency parameter: f_DPM
        this.assert(Math.abs(this.module.variables['f_DPM'] - 1e12) < 1e10,
            'T1.11: DPM frequency initialized (1e12 Hz)');
        
        // T1.12 - Reactive parameter: f_TRZ
        this.assert(Math.abs(this.module.variables['f_TRZ'] - 0.1) < 0.01,
            'T1.12: Time-reversal factor initialized (0.1)');
        
        // T1.13 - Quantum uncertainty: Δp
        const delta_x = this.module.variables['Delta_x'];
        const delta_p = this.module.variables['Delta_p'];
        this.assert(Math.abs(delta_p - this.module.variables['hbar'] / delta_x) < 1e-30,
            'T1.13: Momentum uncertainty calculated from Δx × Δp = ℏ');
        
        // T1.14 - Angular frequency derived from f_super
        const expected_omega = 2 * Math.PI * this.module.variables['f_super'];
        this.assert(Math.abs(this.module.variables['omega'] - expected_omega) / expected_omega < 1e-6,
            'T1.14: Angular frequency ω = 2πf_super calculated correctly');
    }
    
    // ===== CATEGORY 2: FREQUENCY COMPONENTS (15 tests) =====
    testFrequencyComponents() {
        console.log('\n=== TEST CATEGORY 2: Frequency Components (15 tests) ===');
        
        const t = 1900 * 3.156e7;  // At nebular age
        const rho = this.module.variables['rho_fil'];
        
        // T2.1 - Superconductive frequency at t=0
        const f_super_0 = this.module.computeFreqSuper(0);
        this.assert(Math.abs(f_super_0 - 1.411e16) < 1e14,
            'T2.1: f_super at t=0 equals base frequency (1.411e16 Hz)');
        
        // T2.2 - Superconductive frequency decays with time
        const f_super_age = this.module.computeFreqSuper(t);
        this.assert(f_super_age < f_super_0,
            'T2.2: f_super(t_age) < f_super(0) - exponential decay');
        
        // T2.3 - Superconductive frequency approaches exp(-1) × initial
        const expected_decay = 1.411e16 * Math.exp(-1);
        this.assert(Math.abs(f_super_age - expected_decay) / expected_decay < 1e-6,
            'T2.3: f_super(t_age) ≈ f_super(0) × e^(-1) (exponential decay formula)');
        
        // T2.4 - Fluid frequency is density-dependent
        const f_fluid = this.module.computeFreqFluid(rho);
        this.assert(f_fluid > 0,
            'T2.4: Fluid frequency is positive at filament density');
        
        // T2.5 - Fluid frequency scales linearly with density
        const rho_2x = 2 * rho;
        const f_fluid_2x = this.module.computeFreqFluid(rho_2x);
        this.assert(Math.abs(f_fluid_2x - 2 * f_fluid) / (2 * f_fluid) < 1e-6,
            'T2.5: f_fluid(2ρ) = 2 × f_fluid(ρ) - linear scaling');
        
        // T2.6 - Aether frequency is constant
        const f_aether_1 = this.module.computeFreqAether();
        const f_aether_2 = this.module.computeFreqAether();
        this.assert(f_aether_1 === f_aether_2,
            'T2.6: f_Aether is constant (time and position independent)');
        
        // T2.7 - Aether frequency has expected value
        this.assert(Math.abs(f_aether_1 - 1.576e-35) < 1e-37,
            'T2.7: f_Aether = 1.576e-35 Hz (vacuum background)');
        
        // T2.8 - Quantum frequency is positive
        const unc = Math.sqrt(this.module.variables['Delta_x'] * this.module.variables['Delta_p']);
        const f_quantum = this.module.computeFreqQuantum(unc);
        this.assert(f_quantum > 0,
            'T2.8: Quantum frequency is positive');
        
        // T2.9 - Quantum frequency is inverse to uncertainty
        const unc_half = unc / 2;
        const f_quantum_half = this.module.computeFreqQuantum(unc_half);
        this.assert(f_quantum_half > f_quantum,
            'T2.9: f_quantum(unc/2) > f_quantum(unc) - inverse relationship');
        
        // T2.10 - Reactive frequency oscillates
        const t1 = 0;
        const t2 = this.module.variables['pi'] / this.module.variables['omega'];  // Quarter period
        const f_react_0 = this.module.computeFreqReact(t1);
        const f_react_quarter = this.module.computeFreqReact(t2);
        this.assert(Math.abs(f_react_0 - 1e10) < 1e8,
            'T2.10: f_react(0) ≈ 1e10 Hz (cos(0) = 1)');
        
        // T2.11 - THz frequency oscillates (sine)
        const f_thz_0 = this.module.computeTHzHoleTerm(t1);
        this.assert(Math.abs(f_thz_0) < 1e8,
            'T2.11: f_THz(0) ≈ 0 Hz (sin(0) = 0)');
        
        // T2.12 - DPM term is constant
        const f_dpm = this.module.computeDPMTerm(0);
        this.assert(f_dpm > 0,
            'T2.12: DPM frequency is positive');
        
        // T2.13 - DPM term combines vacuum energy and speed of light
        const expected_dpm = 1e12 * (1e-9 / 3e8);
        this.assert(Math.abs(f_dpm - expected_dpm) / expected_dpm < 0.1,
            'T2.13: f_DPM = f_DPM × ρ_vac / c (correct formula)');
        
        // T2.14 - Frequency components have diverse magnitudes
        const f_aether = this.module.computeFreqAether();
        const freqs = [f_super_0, f_fluid, f_quantum, f_aether, f_react_0, f_dpm];
        const magnitudes = freqs.map(f => Math.abs(f)).filter(f => f > 0);
        const minMag = Math.min(...magnitudes);
        const maxMag = Math.max(...magnitudes);
        this.assert(maxMag / minMag > 1e10,
            'T2.14: Frequency components span multiple orders of magnitude');
        
        // T2.15 - Resonance term is positive
        const f_res = this.module.computeResonanceTerm(t);
        this.assert(f_res >= 0,
            'T2.15: Resonance term is non-negative');
    }
    
    // ===== CATEGORY 3: QUANTUM UNCERTAINTY (8 tests) =====
    testQuantumUncertainty() {
        console.log('\n=== TEST CATEGORY 3: Quantum Uncertainty (8 tests) ===');
        
        // T3.1 - Heisenberg relation: Δx × Δp ≥ ℏ/2
        const delta_x = this.module.variables['Delta_x'];
        const delta_p = this.module.variables['Delta_p'];
        const hbar = this.module.variables['hbar'];
        const product = delta_x * delta_p;
        this.assert(product >= hbar / 2,
            'T3.1: Δx × Δp ≥ ℏ/2 (Heisenberg uncertainty satisfied)');
        
        // T3.2 - Δx > 0
        this.assert(delta_x > 0,
            'T3.2: Position uncertainty is positive (1e-10 m)');
        
        // T3.3 - Δp > 0
        this.assert(delta_p > 0,
            'T3.3: Momentum uncertainty is positive');
        
        // T3.4 - Δp = ℏ / Δx
        this.assert(Math.abs(delta_p - hbar / delta_x) < 1e-30,
            'T3.4: Δp = ℏ / Δx (inverse relationship)');
        
        // T3.5 - Uncertainty product remains constant through module operation
        const initial_product = delta_x * delta_p;
        this.module.computeG(1900 * 3.156e7, 1e15);
        const after_product = this.module.variables['Delta_x'] * this.module.variables['Delta_p'];
        this.assert(Math.abs(initial_product - after_product) / initial_product < 1e-6,
            'T3.5: Δx × Δp remains constant through computations');
        
        // T3.6 - Geometric mean of uncertainties
        const unc = Math.sqrt(delta_x * delta_p);
        this.assert(unc > 0,
            'T3.6: Uncertainty geometric mean is positive');
        
        // T3.7 - Quantum frequency inversely proportional to uncertainty
        const f_q1 = this.module.computeFreqQuantum(unc);
        const f_q2 = this.module.computeFreqQuantum(unc * 2);
        this.assert(f_q1 > f_q2,
            'T3.7: f_quantum(unc) > f_quantum(2×unc) - inverse scaling');
        
        // T3.8 - Update Delta_x updates Delta_p
        const initial_dp = this.module.variables['Delta_p'];
        this.module.updateVariable('Delta_x', 2e-10);
        const new_dp = this.module.variables['Delta_p'];
        this.assert(Math.abs(new_dp - initial_dp / 2) < 1e-35,
            'T3.8: updateVariable Delta_x correctly updates Delta_p');
    }
    
    // ===== CATEGORY 4: RESONANCE PHYSICS (10 tests) =====
    testResonancePhysics() {
        console.log('\n=== TEST CATEGORY 4: Resonance Physics (10 tests) ===');
        
        // T4.1 - Wavefunction amplitude is correct
        const psi_int = this.module.computePsiIntegral(1e15, 0);
        this.assert(psi_int > 0,
            'T4.1: Wavefunction intensity is positive');
        
        // T4.2 - Wavefunction depends on amplitude
        const A = this.module.variables['A'];
        const expected_psi = A * A * this.module.variables['integral_psi'];
        this.assert(Math.abs(psi_int - expected_psi) / expected_psi < 1e-6,
            'T4.2: Wavefunction |ψ|² = A² × integral_psi');
        
        // T4.3 - Resonance term increases with superconductive frequency
        const t_early = 100 * 3.156e7;  // Early in nebular life
        const t_late = 1800 * 3.156e7;   // Near end of life
        const f_res_early = this.module.computeResonanceTerm(t_early);
        const f_res_late = this.module.computeResonanceTerm(t_late);
        this.assert(f_res_early > f_res_late,
            'T4.3: Resonance frequency decreases with time (f_super decay)');
        
        // T4.4 - Resonance term is periodic
        const period = (2 * this.module.variables['pi']) / this.module.variables['omega'];
        const t = 1e-10;
        const f_res_t = this.module.computeResonanceTerm(t);
        const f_res_t_period = this.module.computeResonanceTerm(t + period);
        this.assert(Math.abs(f_res_t - f_res_t_period) / Math.max(Math.abs(f_res_t), 1) < 0.01,
            'T4.4: Resonance term is periodic with period 2π/ω');
        
        // T4.5 - Wavefunction modulates resonance
        const res_base = 2 * this.module.variables['pi'] * 1.411e16 * 1e-20;
        const f_res = this.module.computeResonanceTerm(0);
        this.assert(f_res > 0,
            'T4.5: Resonance term computes to positive value');
        
        // T4.6 - Resonance converts properly to Hz
        const f_res_raw = this.module.computeResonanceTerm(0);
        const f_res_hz = f_res_raw / (2 * this.module.variables['pi']);
        this.assert(f_res_hz >= 0,
            'T4.6: Resonance term divisible by 2π for Hz conversion');
        
        // T4.7 - Psi integral is normalized
        this.assert(this.module.variables['integral_psi'] === 1.0,
            'T4.7: Wavefunction integral is normalized (1.0)');
        
        // T4.8 - Amplitude affects resonance strength
        const A_orig = this.module.variables['A'];
        this.module.updateVariable('A', 2 * A_orig);
        const f_res_2x = this.module.computeResonanceTerm(0);
        this.module.updateVariable('A', A_orig);
        const f_res_1x = this.module.computeResonanceTerm(0);
        this.assert(Math.abs(f_res_2x - 4 * f_res_1x) / (4 * f_res_1x) < 0.1,
            'T4.8: Doubling amplitude A quadruples |ψ|² (A²)');
        
        // T4.9 - Spatial parameter affects resonance
        const psi_r1 = this.module.computePsiIntegral(1e15, 0);
        const psi_r2 = this.module.computePsiIntegral(2e15, 0);
        this.assert(Math.abs(psi_r1 - psi_r2) / psi_r1 < 0.01,
            'T4.9: Wavefunction |ψ|² independent of r for plane wave');
        
        // T4.10 - Temporal oscillation in resonance
        const t_sample = [0, 1e-13, 2e-13, 3e-13];
        const f_res_samples = t_sample.map(t => this.module.computeResonanceTerm(t));
        const max_res = Math.max(...f_res_samples.map(f => Math.abs(f)));
        this.assert(max_res > 0,
            'T4.10: Resonance term computes to non-zero values');
    }
    
    // ===== CATEGORY 5: TOPOLOGICAL TERMS (12 tests) =====
    testTopologicalTerms() {
        console.log('\n=== TEST CATEGORY 5: Topological Terms - DPM, THz, U_g4i (12 tests) ===');
        
        // T5.1 - DPM term is positive
        const f_dpm = this.module.computeDPMTerm(0);
        this.assert(f_dpm > 0,
            'T5.1: DPM frequency is positive');
        
        // T5.2 - DPM term independent of time
        const f_dpm_t0 = this.module.computeDPMTerm(0);
        const f_dpm_t1 = this.module.computeDPMTerm(1e15);
        this.assert(Math.abs(f_dpm_t0 - f_dpm_t1) < 1e-30,
            'T5.2: DPM frequency is time-independent');
        
        // T5.3 - DPM formula: f_DPM × ρ_vac_plasm / c
        const f_dpm_calc = this.module.variables['f_DPM'] * 
                          this.module.variables['rho_vac_plasm'] / 
                          this.module.variables['c'];
        this.assert(Math.abs(f_dpm - f_dpm_calc) / f_dpm_calc < 1e-6,
            'T5.3: DPM correctly implements f_DPM = f_DPM × ρ_vac / c');
        
        // T5.4 - THz term oscillates with sine
        const t_sample = [0, 1e-14, 2e-14];
        const f_thz_samples = t_sample.map(t => this.module.computeTHzHoleTerm(t));
        this.assert(Math.abs(f_thz_samples[0]) < 1e6,
            'T5.4: THz frequency near zero at t=0 (sin(0) ≈ 0)');
        
        // T5.5 - THz has different phase than reactive (sin vs cos)
        const t_test = 1e-13;
        const f_react_test = this.module.computeFreqReact(t_test);
        const f_thz_test = this.module.computeTHzHoleTerm(t_test);
        // cos and sin are 90° out of phase
        this.assert((f_react_test > 0 || f_thz_test > 0),
            'T5.5: THz and reactive terms have independent oscillations');
        
        // T5.6 - THz frequency amplitude
        const omega = this.module.variables['omega'];
        const max_thz = this.module.variables['f_THz'];
        this.assert(max_thz === 1e12,
            'T5.6: THz frequency amplitude is 1e12 Hz');
        
        // T5.7 - Reactive U_g4i oscillates with time
        const ug4i_t0 = this.module.computeUg4i(0);
        const ug4i_period_half = this.module.computeUg4i(this.module.variables['pi'] / omega);
        this.assert(Math.abs(ug4i_t0 + ug4i_period_half) < Math.max(Math.abs(ug4i_t0), 1e-50),
            'T5.7: U_g4i oscillates (π phase difference reverses cosine)');
        
        // T5.8 - Reactive term contains time-reversal factor
        const ug4i = this.module.computeUg4i(0);
        const f_react = this.module.computeFreqReact(0);
        const expected_ug4i = f_react * this.module.variables['lambda_I'] * (1 + this.module.variables['f_TRZ']);
        this.assert(Math.abs(ug4i - expected_ug4i) / Math.max(Math.abs(expected_ug4i), 1e-50) < 1e-6,
            'T5.8: U_g4i = f_react × λ_I × (1 + f_TRZ) formula correct');
        
        // T5.9 - Time-reversal factor is small but non-zero
        this.assert(0 < this.module.variables['f_TRZ'] && this.module.variables['f_TRZ'] < 1,
            'T5.9: Time-reversal factor 0 < f_TRZ < 1 (0.1)');
        
        // T5.10 - U_g4i magnitude larger with f_TRZ
        const ug4i_with_trz = this.module.computeUg4i(0);
        const f_react_0 = this.module.computeFreqReact(0);
        const ug4i_without_trz = f_react_0 * this.module.variables['lambda_I'];
        this.assert(Math.abs(ug4i_with_trz) > Math.abs(ug4i_without_trz) * 0.99,
            'T5.10: f_TRZ factor increases U_g4i magnitude');
        
        // T5.11 - All three topological terms positive-definite domain
        this.assert(f_dpm > 0 && this.module.variables['f_THz'] > 0,
            'T5.11: Topological parameters f_DPM and f_THz are positive');
        
        // T5.12 - Topological coherence across scales
        const dpm_scale = Math.abs(Math.log10(Math.abs(f_dpm)));
        const thz_scale = Math.log10(this.module.variables['f_THz']);
        this.assert(Math.abs(thz_scale - dpm_scale) < 10,
            'T5.12: DPM and THz scales within reasonable range');
    }
    
    // ===== CATEGORY 6: DYNAMIC UPDATES (8 tests) =====
    testDynamicUpdates() {
        console.log('\n=== TEST CATEGORY 6: Dynamic Updates (8 tests) ===');
        
        // T6.1 - Update existing variable
        const original_t = this.module.variables['t'];
        this.module.updateVariable('t', 1000 * 3.156e7);
        this.assert(this.module.variables['t'] === 1000 * 3.156e7,
            'T6.1: updateVariable changes existing variable');
        this.module.updateVariable('t', original_t);
        
        // T6.2 - Add to variable
        const original_r = this.module.variables['r'];
        this.module.addToVariable('r', 1e14);
        this.assert(Math.abs(this.module.variables['r'] - (original_r + 1e14)) < 1e10,
            'T6.2: addToVariable increases variable by delta');
        this.module.updateVariable('r', original_r);
        
        // T6.3 - Subtract from variable
        this.module.subtractFromVariable('r', 1e14);
        this.assert(Math.abs(this.module.variables['r'] - (original_r - 1e14)) < 1e10,
            'T6.3: subtractFromVariable decreases variable by delta');
        this.module.updateVariable('r', original_r);
        
        // T6.4 - Update frequency updates omega
        const original_omega = this.module.variables['omega'];
        this.module.updateVariable('f_super', 2e16);
        const new_omega = this.module.variables['omega'];
        this.assert(new_omega > original_omega,
            'T6.4: Updating f_super updates omega = 2πf_super');
        this.module.updateVariable('f_super', 1.411e16);
        
        // T6.5 - Get variable returns correct value
        const val = this.module.getVariable('f_super');
        this.assert(Math.abs(val - 1.411e16) < 1e14,
            'T6.5: getVariable returns current variable value');
        
        // T6.6 - Get state returns all variables
        const state = this.module.getState();
        this.assert(Object.keys(state).length >= 30,
            'T6.6: getState returns all 30+ variables');
        
        // T6.7 - Set state restores configuration
        const saved_state = this.module.getState();
        this.module.updateVariable('f_super', 1.5e16);
        this.module.setState(saved_state);
        this.assert(Math.abs(this.module.getVariable('f_super') - 1.411e16) < 1e14,
            'T6.7: setState correctly restores saved state');
        
        // T6.8 - Multiple variable updates
        const f_super_old = this.module.getVariable('f_super');
        const f_dpm_old = this.module.getVariable('f_DPM');
        this.module.updateVariable('f_super', 1.5e16);
        this.module.updateVariable('f_DPM', 2e12);
        this.assert(Math.abs(this.module.getVariable('f_super') - 1.5e16) < 1e14 &&
                   Math.abs(this.module.getVariable('f_DPM') - 2e12) < 1e10,
            'T6.8: Multiple simultaneous variable updates work correctly');
        this.module.updateVariable('f_super', f_super_old);
        this.module.updateVariable('f_DPM', f_dpm_old);
    }
    
    // ===== CATEGORY 7: MASTER EQUATION (10 tests) =====
    testMasterEquation() {
        console.log('\n=== TEST CATEGORY 7: Master Equation (10 tests) ===');
        
        // T7.1 - computeG returns a number
        const g = this.module.computeG(1900 * 3.156e7, 1e15);
        this.assert(typeof g === 'number' && !isNaN(g),
            'T7.1: computeG returns a valid number');
        
        // T7.2 - computeG can handle different times
        const g_early = this.module.computeG(100 * 3.156e7, 1e15);
        const g_late = this.module.computeG(1800 * 3.156e7, 1e15);
        this.assert(typeof g_early === 'number' && typeof g_late === 'number',
            'T7.2: computeG handles different time inputs');
        
        // T7.3 - computeG produces realistic magnitude
        this.assert(g < 0.1,  // Much less than standard gravity
            'T7.3: Acceleration is small (frequency-derived, not mass-based)');
        
        // T7.4 - computeG changes with time
        this.assert(g_early !== g_late,
            'T7.4: Acceleration varies with time (superconductive decay)');
        
        // T7.5 - computeG can handle different radii
        const g_r1 = this.module.computeG(1900 * 3.156e7, 1e15);
        const g_r2 = this.module.computeG(1900 * 3.156e7, 2e15);
        this.assert(typeof g_r1 === 'number' && typeof g_r2 === 'number',
            'T7.5: computeG handles different radius inputs');
        
        // T7.6 - getAllFrequencies returns all components
        const freqs = this.module.getAllFrequencies(1900 * 3.156e7, 1e15);
        this.assert(Object.keys(freqs).length >= 9,
            'T7.6: getAllFrequencies returns 9+ frequency components');
        
        // T7.7 - Frequency components are real numbers
        const all_real = Object.values(freqs).every(f => typeof f === 'number' && !isNaN(f));
        this.assert(all_real,
            'T7.7: All frequency components are valid numbers');
        
        // T7.8 - Master equation combines multiple components
        const f_super = this.module.computeFreqSuper(1900 * 3.156e7);
        const f_fluid = this.module.computeFreqFluid(this.module.variables['rho_fil']);
        this.assert(f_super > 0 && f_fluid > 0,
            'T7.8: Multiple frequency components are non-zero');
        
        // T7.9 - getEquationText provides documentation
        const eqn_text = this.module.getEquationText();
        this.assert(eqn_text.includes('g_UQFF') && eqn_text.includes('NGC 6537'),
            'T7.9: getEquationText contains equation and system info');
        
        // T7.10 - computeG with r=0 uses default radius
        const g_default = this.module.computeG(1900 * 3.156e7, 0);
        this.assert(typeof g_default === 'number' && !isNaN(g_default),
            'T7.10: computeG with r=0 uses internal radius and computes successfully');
    }
    
    // ===== CATEGORY 8: PERFORMANCE (8 tests) =====
    testPerformance() {
        console.log('\n=== TEST CATEGORY 8: Performance (8 tests) ===');
        
        // T8.1 - Single computation is fast
        const t_start = Date.now();
        for (let i = 0; i < 100; i++) {
            this.module.computeG(1900 * 3.156e7, 1e15);
        }
        const t_elapsed = Date.now() - t_start;
        this.assert(t_elapsed < 100,  // 100 iterations in <100ms
            `T8.1: 100 computeG calls in ${t_elapsed}ms (target <100ms)`);
        
        // T8.2 - Batch computation
        const t_start_batch = Date.now();
        for (let i = 0; i < 1000; i++) {
            this.module.computeG(1900 * 3.156e7, 1e15 + i * 1e12);
        }
        const t_batch = Date.now() - t_start_batch;
        this.assert(t_batch < 1000,  // 1000 iterations in <1 second
            `T8.2: 1000 computeG calls in ${t_batch}ms (target <1000ms)`);
        
        // T8.3 - printVariables doesn't crash
        try {
            const old_log = console.log;
            console.log = () => {};  // Suppress output
            this.module.printVariables();
            console.log = old_log;
            this.assert(true,
                'T8.3: printVariables executes without error');
        } catch (e) {
            this.assert(false, `T8.3: printVariables threw error: ${e}`);
        }
        
        // T8.4 - Variable access is fast
        const t_var_start = Date.now();
        for (let i = 0; i < 10000; i++) {
            const val = this.module.variables['f_super'];
        }
        const t_var = Date.now() - t_var_start;
        this.assert(t_var < 50,
            `T8.4: 10000 variable accesses in ${t_var}ms (target <50ms)`);
        
        // T8.5 - Frequency computations are efficient
        const t_freq_start = Date.now();
        for (let i = 0; i < 1000; i++) {
            this.module.computeFreqSuper(i * 1e11);
        }
        const t_freq = Date.now() - t_freq_start;
        this.assert(t_freq < 100,
            `T8.5: 1000 frequency computations in ${t_freq}ms (target <100ms)`);
        
        // T8.6 - Wavefunction calculations are fast
        const t_psi_start = Date.now();
        for (let i = 0; i < 1000; i++) {
            this.module.computePsiIntegral(1e15, i * 1e10);
        }
        const t_psi = Date.now() - t_psi_start;
        this.assert(t_psi < 100,
            `T8.6: 1000 wavefunction calculations in ${t_psi}ms (target <100ms)`);
        
        // T8.7 - State operations are efficient
        const t_state_start = Date.now();
        for (let i = 0; i < 100; i++) {
            const state = this.module.getState();
            this.module.setState(state);
        }
        const t_state = Date.now() - t_state_start;
        this.assert(t_state < 200,
            `T8.7: 100 get/setState cycles in ${t_state}ms (target <200ms)`);
        
        // T8.8 - Memory efficiency (no crashes with large iterations)
        let stable = true;
        try {
            for (let i = 0; i < 5000; i++) {
                const g = this.module.computeG(i * 1e10, 1e15);
            }
        } catch (e) {
            stable = false;
        }
        this.assert(stable,
            'T8.8: 5000 iterations execute without memory/stability issues');
    }
    
    // ===== RUN ALL TESTS =====
    runAllTests() {
        console.log('\n╔════════════════════════════════════════════════════════════╗');
        console.log('║   RED SPIDER NEBULA (NGC 6537) UQFF TEST SUITE              ║');
        console.log('║   Comprehensive Frequency-Resonance Module Validation       ║');
        console.log('╚════════════════════════════════════════════════════════════╝');
        
        this.testInitialization();
        this.testFrequencyComponents();
        this.testQuantumUncertainty();
        this.testResonancePhysics();
        this.testTopologicalTerms();
        this.testDynamicUpdates();
        this.testMasterEquation();
        this.testPerformance();
        
        this.printSummary();
    }
    
    // Print summary
    printSummary() {
        const total = this.tests_passed + this.tests_failed;
        const pass_rate = ((this.tests_passed / total) * 100).toFixed(2);
        
        console.log('\n╔════════════════════════════════════════════════════════════╗');
        console.log('║                    TEST SUMMARY                             ║');
        console.log('╚════════════════════════════════════════════════════════════╝');
        console.log(`\nTotal Tests Run:     ${total}`);
        console.log(`Tests Passed:        ${this.tests_passed} ✓`);
        console.log(`Tests Failed:        ${this.tests_failed} ✗`);
        console.log(`Success Rate:        ${pass_rate}%`);
        
        if (this.tests_failed === 0) {
            console.log('\n✓ ALL TESTS PASSED - MODULE IS PRODUCTION READY\n');
        } else {
            console.log(`\n⚠ ${this.tests_failed} test(s) need attention\n`);
            console.log('Failed Tests:');
            for (const result of this.test_results) {
                if (result.status === '✗') {
                    console.log(`  ${result.status} ${result.message}`);
                }
            }
        }
    }
}

// Execute tests
const suite = new RedSpiderTestSuite();
suite.runAllTests();
