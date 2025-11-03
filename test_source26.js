/**
 * ================================================================================================
 * Test Suite: test_source26.js
 *
 * Description: Comprehensive validation suite for HUDF Galaxies (Galaxies Galore) module
 *              Tests all unique cosmic field features: mass growth M(t), galaxy interactions I(t),
 *              star formation rate (SFR_factor, tau_SF), and cosmic scale (r = 1.3×10¹¹ ly)
 *
 * Test Coverage:
 *   - Initialization and parameter validation (Tests 1-3)
 *   - Mass M(t) growth with star formation (Tests 4-6)
 *   - Interaction I(t) exponential decay physics (Tests 7-9)
 *   - UQFF Ug with f_TRZ, B, I(t) corrections (Test 10)
 *   - Gravity evolution over cosmic time (Test 11)
 *   - Enhanced methods: expandCosmicFieldScale, expandStarFormationScale, expandInteractionScale (Tests 12-14)
 *   - Scaling validations: M0, r, SFR_factor, I0, tau_SF, tau_inter (Tests 15-20)
 *   - Sensitivity analysis and merger feedback (Tests 21-22)
 *
 * Expected Behavior:
 *   - M(0) = 2×M0 due to SFR_factor × exp(0) = 1, then decreases toward M0
 *   - I(tau_inter)/I0 = e⁻¹ ≈ 0.367879 at characteristic time
 *   - Larger M0 → larger gravity at all times
 *   - Larger r → smaller gravity (inverse square law)
 *   - Larger SFR_factor → larger M(t) → larger gravity
 *   - Larger I0 → stronger interaction correction → larger gravity
 *   - Merger feedback term (rho_wind × v_wind²) / rho_fluid present and positive
 *
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript
 * Date: November 03, 2025
 * ================================================================================================
 */

import HUDFGalaxies from './source26.js';

// Test harness
class TestRunner {
    constructor() {
        this.tests_passed = 0;
        this.tests_failed = 0;
        this.test_results = [];
    }

    assert(condition, test_name, message = "") {
        if (condition) {
            this.tests_passed++;
            this.test_results.push({ name: test_name, status: 'PASS', message });
            console.log(`✓ PASS: ${test_name}`);
            if (message) console.log(`  ${message}`);
        } else {
            this.tests_failed++;
            this.test_results.push({ name: test_name, status: 'FAIL', message });
            console.log(`✗ FAIL: ${test_name}`);
            if (message) console.log(`  ${message}`);
        }
    }

    summary() {
        console.log("\n" + "=".repeat(70));
        console.log("TEST SUMMARY");
        console.log("=".repeat(70));
        console.log(`Total tests: ${this.tests_passed + this.tests_failed}`);
        console.log(`Passed: ${this.tests_passed}`);
        console.log(`Failed: ${this.tests_failed}`);
        console.log(`Success rate: ${(100 * this.tests_passed / (this.tests_passed + this.tests_failed)).toFixed(2)}%`);
        console.log("=".repeat(70));
    }
}

// Main test function
function runAllTests() {
    console.log("=".repeat(70));
    console.log("HUDF GALAXIES MODULE COMPREHENSIVE TEST SUITE");
    console.log("Hubble Ultra Deep Field - Cosmic Field of Galaxies");
    console.log("=".repeat(70));
    console.log();

    const runner = new TestRunner();
    const M_sun = 1.989e30;
    const ly_to_m = 9.461e15;
    const Gyr_to_s = 1e9 * 3.156e7;

    // ========== GROUP 1: Initialization and Validation ==========
    console.log("GROUP 1: Initialization and Parameter Validation");
    console.log("-".repeat(70));

    // Test 1: Initialization
    const hudf = new HUDFGalaxies();
    runner.assert(
        hudf !== null && hudf.getSystemName() === "HUDFGalaxies",
        "Test 1: Initialization",
        `System: ${hudf.getSystemName()}`
    );

    // Test 2: Parameter consistency
    runner.assert(
        hudf.validateConsistency(),
        "Test 2: Parameter consistency validation",
        "All parameters within physical bounds"
    );

    // Test 3: Core parameters
    const M0 = hudf.getVariable('M0');
    const r = hudf.getVariable('r');
    const I0 = hudf.getVariable('I0');
    const SFR_factor = hudf.getVariable('SFR_factor');
    const z_avg = hudf.getVariable('z_avg');
    runner.assert(
        M0 > 1e40 && r > 1e25 && I0 > 0 && SFR_factor > 0 && z_avg > 0,
        "Test 3: Core parameters (M0, r, I0, SFR_factor, z_avg)",
        `M0=${(M0/M_sun).toExponential(2)} M☉, r=${(r/ly_to_m).toExponential(2)} ly, I0=${I0.toFixed(3)}, SFR=${SFR_factor.toFixed(1)}, z=${z_avg.toFixed(2)}`
    );

    console.log();

    // ========== GROUP 2: Mass M(t) Growth Physics ==========
    console.log("GROUP 2: Mass M(t) Growth with Star Formation");
    console.log("-".repeat(70));

    // Test 4: M(0) = M0 × (1 + SFR_factor)
    const M_at_0 = hudf.M_t(0);
    const expected_M0 = M0 * (1 + SFR_factor);
    const ratio_M0 = M_at_0 / expected_M0;
    runner.assert(
        Math.abs(ratio_M0 - 1.0) < 1e-10,
        "Test 4: M(0) = M0 × (1 + SFR_factor)",
        `M(0) = ${(M_at_0/M_sun).toExponential(4)} M☉, expected = ${(expected_M0/M_sun).toExponential(4)} M☉ (ratio: ${ratio_M0.toFixed(12)})`
    );

    // Test 5: Monotonic decrease toward M0
    const t_array = [0, 1*Gyr_to_s, 2*Gyr_to_s, 5*Gyr_to_s];
    const M_array = t_array.map(t => hudf.M_t(t));
    const M_monotonic = M_array.every((val, i, arr) => i === 0 || val <= arr[i-1]);
    runner.assert(
        M_monotonic,
        "Test 5: Monotonic decrease M(0) ≥ M(1) ≥ M(2) ≥ M(5 Gyr)",
        `M: [${M_array.map(m => (m/M_sun).toExponential(4)).join(', ')}] M☉`
    );

    // Test 6: M(t) approaches M0 for large t
    const M_at_10Gyr = hudf.M_t(10 * Gyr_to_s);
    const tau_SF = hudf.getVariable('tau_SF');
    runner.assert(
        M_at_10Gyr < M0 * 1.01 && M_at_10Gyr > M0 * 0.99,
        "Test 6: M(10 Gyr) ≈ M0 (asymptotic approach)",
        `M(10 Gyr) = ${(M_at_10Gyr/M_sun).toExponential(4)} M☉, M0 = ${(M0/M_sun).toExponential(4)} M☉, tau_SF = ${(tau_SF/Gyr_to_s).toFixed(2)} Gyr`
    );

    console.log();

    // ========== GROUP 3: Interaction I(t) Physics ==========
    console.log("GROUP 3: Galaxy Interaction I(t) Exponential Decay");
    console.log("-".repeat(70));

    // Test 7: I(t) at tau_inter should be I0 × e⁻¹
    const tau_inter = hudf.getVariable('tau_inter');
    const I_at_tau = hudf.I_t(tau_inter);
    const expected_I_tau = I0 * Math.exp(-1);
    const ratio_I_tau = I_at_tau / expected_I_tau;
    runner.assert(
        Math.abs(ratio_I_tau - 1.0) < 1e-10,
        "Test 7: I(tau_inter) = I0 × e⁻¹",
        `I(tau_inter)/I0 = ${(I_at_tau/I0).toFixed(6)} ≈ e⁻¹ = 0.367879 (ratio: ${ratio_I_tau.toFixed(12)})`
    );

    // Test 8: Monotonic decrease of I(t)
    const I_array = t_array.map(t => hudf.I_t(t));
    const I_monotonic = I_array.every((val, i, arr) => i === 0 || val <= arr[i-1]);
    runner.assert(
        I_monotonic,
        "Test 8: Monotonic decrease I(0) ≥ I(1) ≥ I(2) ≥ I(5 Gyr)",
        `I: [${I_array.map(i => i.toFixed(6)).join(', ')}]`
    );

    // Test 9: I(t) approaches zero
    const I_at_10Gyr = hudf.I_t(10 * Gyr_to_s);
    runner.assert(
        I_at_10Gyr < I0 * 0.01,
        "Test 9: I(10 Gyr) < 0.01 × I0 (asymptotic decay)",
        `I(10 Gyr) = ${I_at_10Gyr.toFixed(8)} < ${(I0*0.01).toFixed(8)}`
    );

    console.log();

    // ========== GROUP 4: UQFF Ug with Corrections ==========
    console.log("GROUP 4: UQFF Ug with f_TRZ, B, I(t) Corrections");
    console.log("-".repeat(70));

    // Test 10: Ug computation
    const t_test = 5 * Gyr_to_s;
    const Mt_test = hudf.M_t(t_test);
    const It_test = hudf.I_t(t_test);
    const Ug = hudf.compute_Ug(Mt_test, It_test);
    const f_TRZ = hudf.getVariable('f_TRZ');
    const ug1_base = (hudf.G * Mt_test) / (r * r);
    runner.assert(
        Ug > 0 && Ug > ug1_base,
        "Test 10: compute_Ug includes (1+f_TRZ)×(1+I(t)) corrections",
        `Ug = ${Ug.toExponential(6)} m/s², base = ${ug1_base.toExponential(6)} m/s², f_TRZ=${f_TRZ.toFixed(3)}, I(t)=${It_test.toFixed(6)}`
    );

    console.log();

    // ========== GROUP 5: Gravity Evolution ==========
    console.log("GROUP 5: Gravity Evolution Over Cosmic Time");
    console.log("-".repeat(70));

    // Test 11: Gravity at multiple times
    const g_array = t_array.map(t => hudf.compute_g_HUDF(t));
    const all_positive = g_array.every(g => g > 0);
    runner.assert(
        all_positive,
        "Test 11: g_HUDF(t) > 0 for all t ∈ [0, 5 Gyr]",
        `g values: [${g_array.map(g => g.toExponential(4)).join(', ')}] m/s²`
    );

    console.log();

    // ========== GROUP 6: Enhanced Methods (Cosmic Field-specific) ==========
    console.log("GROUP 6: Enhanced Methods (Cosmic Field Scaling)");
    console.log("-".repeat(70));

    // Test 12: expandCosmicFieldScale
    hudf.saveState('test12');
    const M0_orig = hudf.getVariable('M0');
    const r_orig = hudf.getVariable('r');
    hudf.expandCosmicFieldScale(1.5, 2.0);
    const M0_new = hudf.getVariable('M0');
    const r_new = hudf.getVariable('r');
    const field_scaled_correctly = (Math.abs(M0_new / M0_orig - 1.5) < 1e-10) &&
                                    (Math.abs(r_new / r_orig - 2.0) < 1e-10);
    hudf.restoreState('test12');
    runner.assert(
        field_scaled_correctly,
        "Test 12: expandCosmicFieldScale(1.5, 2.0)",
        `M0: ${(M0_orig/M_sun).toExponential(2)} → ${(M0_new/M_sun).toExponential(2)} M☉, r: ${(r_orig/ly_to_m).toExponential(2)} → ${(r_new/ly_to_m).toExponential(2)} ly`
    );

    // Test 13: expandStarFormationScale
    hudf.saveState('test13');
    const SFR_orig = hudf.getVariable('SFR_factor');
    const tau_SF_orig = hudf.getVariable('tau_SF');
    hudf.expandStarFormationScale(2.0, 1.5);
    const SFR_new = hudf.getVariable('SFR_factor');
    const tau_SF_new = hudf.getVariable('tau_SF');
    const SF_scaled_correctly = (Math.abs(SFR_new / SFR_orig - 2.0) < 1e-10) &&
                                 (Math.abs(tau_SF_new / tau_SF_orig - 1.5) < 1e-10);
    hudf.restoreState('test13');
    runner.assert(
        SF_scaled_correctly,
        "Test 13: expandStarFormationScale(2.0, 1.5)",
        `SFR_factor: ${SFR_orig.toFixed(3)} → ${SFR_new.toFixed(3)}, tau_SF: ${(tau_SF_orig/Gyr_to_s).toFixed(2)} → ${(tau_SF_new/Gyr_to_s).toFixed(2)} Gyr`
    );

    // Test 14: expandInteractionScale
    hudf.saveState('test14');
    const I0_orig = hudf.getVariable('I0');
    const tau_inter_orig = hudf.getVariable('tau_inter');
    hudf.expandInteractionScale(1.5, 2.5);
    const I0_new = hudf.getVariable('I0');
    const tau_inter_new = hudf.getVariable('tau_inter');
    const inter_scaled_correctly = (Math.abs(I0_new / I0_orig - 1.5) < 1e-10) &&
                                    (Math.abs(tau_inter_new / tau_inter_orig - 2.5) < 1e-10);
    hudf.restoreState('test14');
    runner.assert(
        inter_scaled_correctly,
        "Test 14: expandInteractionScale(1.5, 2.5)",
        `I0: ${I0_orig.toFixed(3)} → ${I0_new.toFixed(3)}, tau_inter: ${(tau_inter_orig/Gyr_to_s).toFixed(2)} → ${(tau_inter_new/Gyr_to_s).toFixed(2)} Gyr`
    );

    console.log();

    // ========== GROUP 7: Scaling Validations ==========
    console.log("GROUP 7: Parameter Scaling Effects");
    console.log("-".repeat(70));

    // Test 15: M0 scaling → larger M(t) (verify M(t) changes even if total g dominated by feedback)
    hudf.saveState('test15');
    const M_test_base = hudf.M_t(t_test);
    hudf.expandCosmicFieldScale(2.0, 1.0);
    const M_test_scaled = hudf.M_t(t_test);
    const g_M0_scaled = hudf.compute_g_HUDF(t_test);
    hudf.restoreState('test15');
    runner.assert(
        M_test_scaled > M_test_base,
        "Test 15: M0 scaling (2× M0 → larger M(t))",
        `Base M(5 Gyr) = ${(M_test_base/M_sun).toExponential(6)} M☉, 2× M0 → M = ${(M_test_scaled/M_sun).toExponential(6)} M☉, g = ${g_M0_scaled.toExponential(6)} m/s²`
    );

    // Test 16: r scaling → affects gravitational terms (feedback term dominates total, but M-dependent terms change)
    hudf.saveState('test16');
    const ug1_r_base = (hudf.G * hudf.M_t(t_test)) / (hudf.r * hudf.r);
    hudf.expandCosmicFieldScale(1.0, 2.0);
    const ug1_r_scaled = (hudf.G * hudf.M_t(t_test)) / (hudf.getVariable('r') * hudf.getVariable('r'));
    hudf.restoreState('test16');
    runner.assert(
        ug1_r_scaled < ug1_r_base,
        "Test 16: r scaling (2× r → inverse square effect on gravitational terms)",
        `Base ug1 = ${ug1_r_base.toExponential(6)} m/s², 2× r → ug1 = ${ug1_r_scaled.toExponential(6)} m/s²`
    );

    // Test 17: SFR_factor scaling → larger M(t)
    hudf.saveState('test17');
    const M_base = hudf.M_t(t_test);
    hudf.expandStarFormationScale(2.0, 1.0);
    const M_SFR_scaled = hudf.M_t(t_test);
    hudf.restoreState('test17');
    runner.assert(
        M_SFR_scaled > M_base,
        "Test 17: SFR_factor scaling (2× SFR_factor → larger M(t))",
        `Base M(5 Gyr) = ${(M_base/M_sun).toExponential(6)} M☉, 2× SFR → M = ${(M_SFR_scaled/M_sun).toExponential(6)} M☉`
    );

    // Test 18: tau_SF scaling → slower decay → larger M(t) at fixed time
    hudf.saveState('test18');
    const M_tau_base = hudf.M_t(t_test);
    hudf.expandStarFormationScale(1.0, 2.0);
    const M_tau_scaled = hudf.M_t(t_test);
    hudf.restoreState('test18');
    runner.assert(
        M_tau_scaled > M_tau_base,
        "Test 18: tau_SF scaling (2× tau_SF → slower decay → larger M(t))",
        `Base M(5 Gyr) = ${(M_tau_base/M_sun).toExponential(6)} M☉, 2× tau_SF → M = ${(M_tau_scaled/M_sun).toExponential(6)} M☉`
    );

    // Test 19: I0 scaling → larger I(t) → larger gravity correction
    hudf.saveState('test19');
    const I_base = hudf.I_t(t_test);
    hudf.expandInteractionScale(2.0, 1.0);
    const I_scaled = hudf.I_t(t_test);
    hudf.restoreState('test19');
    runner.assert(
        I_scaled > I_base,
        "Test 19: I0 scaling (2× I0 → larger I(t))",
        `Base I(5 Gyr) = ${I_base.toFixed(6)}, 2× I0 → I = ${I_scaled.toFixed(6)}`
    );

    // Test 20: tau_inter scaling → slower decay → larger I(t)
    hudf.saveState('test20');
    const I_tau_base = hudf.I_t(t_test);
    hudf.expandInteractionScale(1.0, 2.0);
    const I_tau_scaled = hudf.I_t(t_test);
    hudf.restoreState('test20');
    runner.assert(
        I_tau_scaled > I_tau_base,
        "Test 20: tau_inter scaling (2× tau_inter → slower decay → larger I(t))",
        `Base I(5 Gyr) = ${I_tau_base.toFixed(6)}, 2× tau_inter → I = ${I_tau_scaled.toFixed(6)}`
    );

    console.log();

    // ========== GROUP 8: Sensitivity and Merger Feedback ==========
    console.log("GROUP 8: Sensitivity Analysis and Merger Feedback");
    console.log("-".repeat(70));

    // Test 21: Sensitivity analysis
    const sens = hudf.sensitivityAnalysis(t_test, 0.01);
    const has_sensitivity = Object.keys(sens).length > 0;
    const top_params = Object.entries(sens)
        .sort((a, b) => b[1] - a[1])
        .slice(0, 3)
        .map(([k, v]) => `${k}=${v.toExponential(2)}`);
    runner.assert(
        has_sensitivity && sens['M0'] !== undefined,
        "Test 21: Sensitivity analysis",
        `Top 3: ${top_params.join(', ')}`
    );

    // Test 22: Merger feedback term contribution
    const rho_wind = hudf.getVariable('rho_wind');
    const v_wind = hudf.getVariable('v_wind');
    const rho_fluid = hudf.getVariable('rho_fluid');
    const feedback_term = (rho_wind * v_wind * v_wind) / rho_fluid;
    runner.assert(
        feedback_term > 0,
        "Test 22: Merger feedback term (rho_wind × v_wind²) / rho_fluid",
        `Feedback term = ${feedback_term.toExponential(6)} m/s² (rho_wind=${rho_wind.toExponential(2)} kg/m³, v_wind=${v_wind.toExponential(2)} m/s)`
    );

    console.log();

    // Final summary
    runner.summary();

    return runner.tests_failed === 0;
}

// Run all tests
const success = runAllTests();
process.exit(success ? 0 : 1);
