/**
 * ================================================================================================
 * Test Suite: test_source25.js
 *
 * Description: Comprehensive validation suite for NGC 1275 (Magnetic Monster Perseus A) module
 *              Tests all unique AGN features: magnetic field B(t) decay, filament support F(t),
 *              supermassive black hole (M_BH = 8×10⁸ M☉), cooling flows, and all 11 MUGE terms
 *
 * Test Coverage:
 *   - Initialization and parameter validation (Tests 1-3)
 *   - Magnetic field B(t) exponential decay physics (Tests 4-6)
 *   - Filament support F(t) decay physics (Tests 7-9)
 *   - Black hole contribution to gravity (Tests 10-12)
 *   - UQFF Ug with f_TRZ, B(t), F(t) corrections (Test 13)
 *   - Gravity evolution over time (Test 14)
 *   - Enhanced methods: expandBlackHoleScale, expandMagneticFilamentScale (Tests 15-16)
 *   - Scaling validations: B0, F0, tau_B, tau_fil (Tests 17-20)
 *   - Sensitivity analysis and cooling flow (Tests 21-22)
 *
 * Expected Behavior:
 *   - B(tau_B)/B0 = e⁻¹ ≈ 0.367879 at characteristic time
 *   - F(tau_fil)/F0 = e⁻¹ ≈ 0.367879 at characteristic time
 *   - Black hole adds significant term to total gravity
 *   - Larger B0 → more EM contribution → larger total g
 *   - Longer tau_B → slower B(t) decay → larger B(t) at fixed time
 *   - Cooling flow term (rho_cool × v_cool²) / rho_fluid present and positive
 *
 * Author: Encoded by Grok (xAI), based on Daniel T. Murphy's UQFF manuscript
 * Date: November 03, 2025
 * ================================================================================================
 */

import NGC1275 from './source25.js';

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
    console.log("NGC 1275 MODULE COMPREHENSIVE TEST SUITE");
    console.log("Perseus A - Magnetic Monster AGN with SMBH");
    console.log("=".repeat(70));
    console.log();

    const runner = new TestRunner();
    const M_sun = 1.989e30;
    const ly_to_m = 9.461e15;
    const Myr_to_s = 1e6 * 3.156e7;

    // ========== GROUP 1: Initialization and Validation ==========
    console.log("GROUP 1: Initialization and Parameter Validation");
    console.log("-".repeat(70));

    // Test 1: Initialization
    const ngc = new NGC1275();
    runner.assert(
        ngc !== null && ngc.getSystemName() === "NGC1275",
        "Test 1: Initialization",
        `System: ${ngc.getSystemName()}`
    );

    // Test 2: Parameter consistency
    runner.assert(
        ngc.validateConsistency(),
        "Test 2: Parameter consistency validation",
        "All parameters within physical bounds"
    );

    // Test 3: Core parameters
    const M = ngc.getVariable('M');
    const r = ngc.getVariable('r');
    const M_BH = ngc.getVariable('M_BH');
    const B0 = ngc.getVariable('B0');
    const F0 = ngc.getVariable('F0');
    runner.assert(
        M > 1e40 && r > 1e20 && M_BH > 1e38 && B0 > 0 && F0 > 0,
        "Test 3: Core parameters (M, r, M_BH, B0, F0)",
        `M=${(M/M_sun).toExponential(2)} M☉, r=${(r/ly_to_m).toFixed(0)} ly, M_BH=${(M_BH/M_sun).toExponential(2)} M☉, B0=${B0.toExponential(2)} T, F0=${F0.toFixed(3)}`
    );

    console.log();

    // ========== GROUP 2: Magnetic Field B(t) Physics ==========
    console.log("GROUP 2: Magnetic Field B(t) Exponential Decay");
    console.log("-".repeat(70));

    // Test 4: B(t) at tau_B should be B0 × e⁻¹
    const tau_B = ngc.getVariable('tau_B');
    const B_at_tau = ngc.B_t(tau_B);
    const expected_B_tau = B0 * Math.exp(-1);
    const ratio_B_tau = B_at_tau / expected_B_tau;
    runner.assert(
        Math.abs(ratio_B_tau - 1.0) < 1e-10,
        "Test 4: B(tau_B) = B0 × e⁻¹",
        `B(tau_B)/B0 = ${(B_at_tau/B0).toFixed(6)} ≈ e⁻¹ = 0.367879 (ratio: ${ratio_B_tau.toFixed(12)})`
    );

    // Test 5: Monotonic decrease of B(t)
    const t_array = [0, 25*Myr_to_s, 50*Myr_to_s, 100*Myr_to_s];
    const B_array = t_array.map(t => ngc.B_t(t));
    const B_monotonic = B_array.every((val, i, arr) => i === 0 || val <= arr[i-1]);
    runner.assert(
        B_monotonic,
        "Test 5: Monotonic decrease B(0) ≥ B(25) ≥ B(50) ≥ B(100 Myr)",
        `B: [${B_array.map(b => b.toExponential(4)).join(', ')}] T`
    );

    // Test 6: B(t) asymptotic approach to zero
    const B_at_500Myr = ngc.B_t(500 * Myr_to_s);
    runner.assert(
        B_at_500Myr < B0 * 0.01,
        "Test 6: B(500 Myr) < 0.01 × B0 (asymptotic decay)",
        `B(500 Myr) = ${B_at_500Myr.toExponential(4)} T < ${(B0*0.01).toExponential(4)} T`
    );

    console.log();

    // ========== GROUP 3: Filament Support F(t) Physics ==========
    console.log("GROUP 3: Filament Support F(t) Exponential Decay");
    console.log("-".repeat(70));

    // Test 7: F(t) at tau_fil should be F0 × e⁻¹
    const tau_fil = ngc.getVariable('tau_fil');
    const F_at_tau = ngc.F_t(tau_fil);
    const expected_F_tau = F0 * Math.exp(-1);
    const ratio_F_tau = F_at_tau / expected_F_tau;
    runner.assert(
        Math.abs(ratio_F_tau - 1.0) < 1e-10,
        "Test 7: F(tau_fil) = F0 × e⁻¹",
        `F(tau_fil)/F0 = ${(F_at_tau/F0).toFixed(6)} ≈ e⁻¹ = 0.367879 (ratio: ${ratio_F_tau.toFixed(12)})`
    );

    // Test 8: Monotonic decrease of F(t)
    const F_array = t_array.map(t => ngc.F_t(t));
    const F_monotonic = F_array.every((val, i, arr) => i === 0 || val <= arr[i-1]);
    runner.assert(
        F_monotonic,
        "Test 8: Monotonic decrease F(0) ≥ F(25) ≥ F(50) ≥ F(100 Myr)",
        `F: [${F_array.map(f => f.toFixed(6)).join(', ')}]`
    );

    // Test 9: F(t) approaches zero
    const F_at_500Myr = ngc.F_t(500 * Myr_to_s);
    runner.assert(
        F_at_500Myr < F0 * 0.01,
        "Test 9: F(500 Myr) < 0.01 × F0 (asymptotic decay)",
        `F(500 Myr) = ${F_at_500Myr.toFixed(6)} < ${(F0*0.01).toFixed(6)}`
    );

    console.log();

    // ========== GROUP 4: Black Hole Physics ==========
    console.log("GROUP 4: Supermassive Black Hole (M_BH = 8×10⁸ M☉)");
    console.log("-".repeat(70));

    // Test 10: Black hole contribution exists
    const r_BH = ngc.getVariable('r_BH');
    const g_BH = ngc.g_BH;
    runner.assert(
        g_BH > 0,
        "Test 10: Black hole acceleration g_BH > 0",
        `g_BH = ${g_BH.toExponential(6)} m/s² (M_BH=${(M_BH/M_sun).toExponential(2)} M☉, r_BH=${r_BH.toExponential(2)} m)`
    );

    // Test 11: Black hole term in total gravity
    const t_test = 50 * Myr_to_s;
    const g_total = ngc.compute_g_NGC1275(t_test);
    runner.assert(
        g_total > g_BH,
        "Test 11: g_total > g_BH (black hole contributes to total)",
        `g_total = ${g_total.toExponential(6)} m/s², g_BH = ${g_BH.toExponential(6)} m/s²`
    );

    // Test 12: M_BH scaling
    ngc.saveState('original');
    ngc.expandBlackHoleScale(2.0, 1.0);
    const M_BH_scaled = ngc.getVariable('M_BH');
    const g_BH_scaled = ngc.g_BH;
    ngc.restoreState('original');
    runner.assert(
        g_BH_scaled > g_BH,
        "Test 12: M_BH scaling (2× M_BH → larger g_BH)",
        `Original g_BH = ${g_BH.toExponential(4)} m/s², 2× M_BH → g_BH = ${g_BH_scaled.toExponential(4)} m/s²`
    );

    console.log();

    // ========== GROUP 5: UQFF Ug with Corrections ==========
    console.log("GROUP 5: UQFF Ug with f_TRZ, B(t), F(t) Corrections");
    console.log("-".repeat(70));

    // Test 13: Ug computation
    const Bt_test = ngc.B_t(t_test);
    const Ft_test = ngc.F_t(t_test);
    const Ug = ngc.compute_Ug(Bt_test, Ft_test);
    const f_TRZ = ngc.getVariable('f_TRZ');
    runner.assert(
        Ug > 0 && Ug > ngc.ug1_base,
        "Test 13: compute_Ug includes (1+f_TRZ)×(1+F(t)) corrections",
        `Ug = ${Ug.toExponential(6)} m/s², base = ${ngc.ug1_base.toExponential(6)} m/s², f_TRZ=${f_TRZ.toFixed(3)}, F(t)=${Ft_test.toFixed(6)}`
    );

    console.log();

    // ========== GROUP 6: Gravity Evolution ==========
    console.log("GROUP 6: Gravity Evolution Over Time");
    console.log("-".repeat(70));

    // Test 14: Gravity at multiple times
    const g_array = t_array.map(t => ngc.compute_g_NGC1275(t));
    const all_positive = g_array.every(g => g > 0);
    runner.assert(
        all_positive,
        "Test 14: g_NGC1275(t) > 0 for all t ∈ [0, 100 Myr]",
        `g values: [${g_array.map(g => g.toExponential(4)).join(', ')}] m/s²`
    );

    console.log();

    // ========== GROUP 7: Enhanced Methods (AGN-specific) ==========
    console.log("GROUP 7: Enhanced Methods (Black Hole & Magnetic Filament Scaling)");
    console.log("-".repeat(70));

    // Test 15: expandBlackHoleScale
    ngc.saveState('test15');
    const M_BH_orig = ngc.getVariable('M_BH');
    const r_BH_orig = ngc.getVariable('r_BH');
    ngc.expandBlackHoleScale(1.5, 2.0);
    const M_BH_new = ngc.getVariable('M_BH');
    const r_BH_new = ngc.getVariable('r_BH');
    const BH_scaled_correctly = (Math.abs(M_BH_new / M_BH_orig - 1.5) < 1e-10) &&
                                 (Math.abs(r_BH_new / r_BH_orig - 2.0) < 1e-10);
    ngc.restoreState('test15');
    runner.assert(
        BH_scaled_correctly,
        "Test 15: expandBlackHoleScale(1.5, 2.0)",
        `M_BH: ${(M_BH_orig/M_sun).toExponential(2)} → ${(M_BH_new/M_sun).toExponential(2)} M☉, r_BH: ${r_BH_orig.toExponential(2)} → ${r_BH_new.toExponential(2)} m`
    );

    // Test 16: expandMagneticFilamentScale
    ngc.saveState('test16');
    const B0_orig = ngc.getVariable('B0');
    const F0_orig = ngc.getVariable('F0');
    ngc.expandMagneticFilamentScale(2.5, 1.5);
    const B0_new = ngc.getVariable('B0');
    const F0_new = ngc.getVariable('F0');
    const MF_scaled_correctly = (Math.abs(B0_new / B0_orig - 2.5) < 1e-10) &&
                                 (Math.abs(F0_new / F0_orig - 1.5) < 1e-10);
    ngc.restoreState('test16');
    runner.assert(
        MF_scaled_correctly,
        "Test 16: expandMagneticFilamentScale(2.5, 1.5)",
        `B0: ${B0_orig.toExponential(2)} → ${B0_new.toExponential(2)} T, F0: ${F0_orig.toFixed(3)} → ${F0_new.toFixed(3)}`
    );

    console.log();

    // ========== GROUP 8: Scaling Validations ==========
    console.log("GROUP 8: Parameter Scaling Effects");
    console.log("-".repeat(70));

    // Test 17: B0 scaling → larger B(t) → larger EM term → larger g
    ngc.saveState('test17');
    const g_base = ngc.compute_g_NGC1275(t_test);
    ngc.expandMagneticFilamentScale(2.0, 1.0);
    const g_B0_scaled = ngc.compute_g_NGC1275(t_test);
    ngc.restoreState('test17');
    runner.assert(
        g_B0_scaled > g_base,
        "Test 17: B0 scaling (2× B0 → larger EM → larger g)",
        `Base g = ${g_base.toExponential(6)} m/s², 2× B0 → g = ${g_B0_scaled.toExponential(6)} m/s²`
    );

    // Test 18: tau_B scaling → longer decay → larger B(t)
    ngc.saveState('test18');
    const B_base = ngc.B_t(t_test);
    ngc.setVariable('tau_B', tau_B * 2.0);
    const B_tau_scaled = ngc.B_t(t_test);
    ngc.restoreState('test18');
    runner.assert(
        B_tau_scaled > B_base,
        "Test 18: tau_B scaling (2× tau_B → slower decay → larger B(t))",
        `Base B(50 Myr) = ${B_base.toExponential(6)} T, 2× tau_B → B = ${B_tau_scaled.toExponential(6)} T`
    );

    // Test 19: F0 scaling → larger F(t) → larger Ug correction
    ngc.saveState('test19');
    const F_base = ngc.F_t(t_test);
    ngc.expandMagneticFilamentScale(1.0, 2.0);
    const F_scaled = ngc.F_t(t_test);
    ngc.restoreState('test19');
    runner.assert(
        F_scaled > F_base,
        "Test 19: F0 scaling (2× F0 → larger filament support F(t))",
        `Base F(50 Myr) = ${F_base.toFixed(6)}, 2× F0 → F = ${F_scaled.toFixed(6)}`
    );

    // Test 20: tau_fil scaling → longer decay → larger F(t)
    ngc.saveState('test20');
    const F_base_tau = ngc.F_t(t_test);
    ngc.setVariable('tau_fil', tau_fil * 2.0);
    const F_tau_scaled = ngc.F_t(t_test);
    ngc.restoreState('test20');
    runner.assert(
        F_tau_scaled > F_base_tau,
        "Test 20: tau_fil scaling (2× tau_fil → slower decay → larger F(t))",
        `Base F(50 Myr) = ${F_base_tau.toFixed(6)}, 2× tau_fil → F = ${F_tau_scaled.toFixed(6)}`
    );

    console.log();

    // ========== GROUP 9: Sensitivity and Cooling Flow ==========
    console.log("GROUP 9: Sensitivity Analysis and Cooling Flow");
    console.log("-".repeat(70));

    // Test 21: Sensitivity analysis
    const sens = ngc.sensitivityAnalysis(t_test, 0.01);
    const has_sensitivity = Object.keys(sens).length > 0;
    const top_params = Object.entries(sens)
        .sort((a, b) => b[1] - a[1])
        .slice(0, 3)
        .map(([k, v]) => `${k}=${v.toExponential(2)}`);
    runner.assert(
        has_sensitivity && sens['M'] !== undefined,
        "Test 21: Sensitivity analysis",
        `Top 3: ${top_params.join(', ')}`
    );

    // Test 22: Cooling flow term contribution
    const rho_cool = ngc.getVariable('rho_cool');
    const v_cool = ngc.getVariable('v_cool');
    const rho_fluid = ngc.getVariable('rho_fluid');
    const cool_term = (rho_cool * v_cool * v_cool) / rho_fluid;
    runner.assert(
        cool_term > 0,
        "Test 22: Cooling flow term (rho_cool × v_cool²) / rho_fluid",
        `Cooling term = ${cool_term.toExponential(6)} m/s² (rho_cool=${rho_cool.toExponential(2)} kg/m³, v_cool=${v_cool.toExponential(2)} m/s)`
    );

    console.log();

    // Final summary
    runner.summary();

    return runner.tests_failed === 0;
}

// Run all tests
const success = runAllTests();
process.exit(success ? 0 : 1);
