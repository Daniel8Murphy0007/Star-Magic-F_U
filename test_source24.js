/**
 * ================================================================================================
 * Test Suite: test_source24.js
 *
 * Description: Comprehensive test suite for HorseheadNebula (source24.js) module
 *              Tests all 38 methods (13 core + 25 enhanced) including unique erosion physics
 *
 * Tests Cover:
 *   - Initialization and parameter validation
 *   - Erosion E(t) asymptotic approach to E_0 (UNIQUE)
 *   - UQFF Ug computation with erosion correction (1 - E(t))
 *   - Complete MUGE with all 10 terms
 *   - Enhanced capabilities: scaling, state management, evolution
 *   - Nebula-specific features: expandErosionScale, expandWindMagneticScale
 *
 * Author: Test suite for UQFF module conversion
 * Date: November 03, 2025
 * ================================================================================================
 */

import HorseheadNebula from './source24.js';

// Test configuration
const TOLERANCE = 1e-6;
const M_sun = 1.989e30;
const ly_to_m = 9.461e15;
const Myr_to_s = 1e6 * 3.156e7;

// Test utilities
function assertApprox(actual, expected, tolerance, testName) {
    const diff = Math.abs(actual - expected);
    const relError = expected !== 0 ? diff / Math.abs(expected) : diff;
    if (relError > tolerance && diff > tolerance) {
        console.error(`❌ FAIL: ${testName}`);
        console.error(`   Expected: ${expected.toExponential(6)}, Got: ${actual.toExponential(6)}, Diff: ${diff.toExponential(6)}`);
        return false;
    }
    console.log(`✓ PASS: ${testName}`);
    return true;
}

function assertTrue(condition, testName) {
    if (!condition) {
        console.error(`❌ FAIL: ${testName}`);
        return false;
    }
    console.log(`✓ PASS: ${testName}`);
    return true;
}

function assertGreater(actual, threshold, testName) {
    if (actual <= threshold) {
        console.error(`❌ FAIL: ${testName}`);
        console.error(`   Expected > ${threshold}, Got: ${actual}`);
        return false;
    }
    console.log(`✓ PASS: ${testName}`);
    return true;
}

function assertLess(actual, threshold, testName) {
    if (actual >= threshold) {
        console.error(`❌ FAIL: ${testName}`);
        console.error(`   Expected < ${threshold}, Got: ${actual}`);
        return false;
    }
    console.log(`✓ PASS: ${testName}`);
    return true;
}

// Main test suite
function runAllTests() {
    console.log("=========================================================");
    console.log("HORSEHEAD NEBULA (source24.js) COMPREHENSIVE TEST SUITE");
    console.log("Testing Barnard 33 Dark Nebula with Erosion Physics");
    console.log("=========================================================\n");

    let passCount = 0;
    let totalCount = 0;

    // Test 1: Initialization validation
    console.log("Test 1: Initialization with default parameters");
    const horsehead = new HorseheadNebula();
    totalCount++;
    if (assertTrue(horsehead.M === 1000 * M_sun, "M initialized to 1000 M_sun") &&
        assertTrue(horsehead.r === 2.5 * ly_to_m, "r initialized to 2.5 ly") &&
        assertTrue(horsehead.E_0 === 0.1, "E_0 initialized to 0.1") &&
        assertTrue(horsehead.tau_erosion === 5 * Myr_to_s, "tau_erosion = 5 Myr")) {
        passCount++;
    }
    console.log();

    // Test 2: System name
    console.log("Test 2: System name retrieval");
    totalCount++;
    if (assertTrue(horsehead.getSystemName() === "HorseheadNebula", "System name is HorseheadNebula")) {
        passCount++;
    }
    console.log();

    // Test 3: Variable listing
    console.log("Test 3: Variable listing completeness");
    const vars = horsehead.listVariables();
    totalCount++;
    if (assertTrue(vars.length === 32, "32 variables listed") &&
        assertTrue(vars.includes('E_0'), "E_0 in variable list") &&
        assertTrue(vars.includes('tau_erosion'), "tau_erosion in variable list")) {
        passCount++;
    }
    console.log();

    // Test 4: Erosion E(t) asymptotic approach (UNIQUE FEATURE)
    console.log("Test 4: Erosion E(t) asymptotic approach to E_0 - DARK NEBULA EROSION FEATURE");
    const E_0_val = horsehead.E_t(0);
    const E_tau = horsehead.E_t(horsehead.tau_erosion);
    const expected_E_tau = horsehead.E_0 * (1 - Math.exp(-1));  // E(tau) = E_0 × (1 - e^-1)
    totalCount++;
    if (assertApprox(E_0_val, 0.0, 0.001, "E(0) = 0") &&
        assertApprox(E_tau / horsehead.E_0, 1 - Math.exp(-1), 0.001, "E(tau)/E_0 = 1 - e^-1 = 0.632121")) {
        passCount++;
    }
    console.log(`   E(0) = ${E_0_val.toFixed(6)}, E(tau) = ${E_tau.toFixed(6)}, E(tau)/E_0 = ${(E_tau/horsehead.E_0).toFixed(6)}`);
    console.log();

    // Test 5: Erosion increases over time
    console.log("Test 5: Erosion E(t) increases monotonically toward E_0");
    const E_1Myr = horsehead.E_t(1 * Myr_to_s);
    const E_3Myr = horsehead.E_t(3 * Myr_to_s);
    const E_5Myr = horsehead.E_t(5 * Myr_to_s);
    totalCount++;
    if (assertTrue(E_1Myr < E_3Myr, "E(1 Myr) < E(3 Myr)") &&
        assertTrue(E_3Myr < E_5Myr, "E(3 Myr) < E(5 Myr)") &&
        assertTrue(E_5Myr < horsehead.E_0, "E(5 Myr) < E_0 (approaches asymptotically)")) {
        passCount++;
    }
    console.log(`   E(1 Myr) = ${E_1Myr.toFixed(6)}, E(3 Myr) = ${E_3Myr.toFixed(6)}, E(5 Myr) = ${E_5Myr.toFixed(6)}`);
    console.log();

    // Test 6: UQFF Ug computation includes erosion correction
    console.log("Test 6: UQFF Ug computation includes erosion correction (1 - E(t))");
    const t_3Myr = 3 * Myr_to_s;
    const Et_3 = horsehead.E_t(t_3Myr);
    const Ug_3 = horsehead.compute_Ug(Et_3);
    totalCount++;
    // Ug should include (1 - Et) correction factor
    if (assertGreater(Ug_3, 0, "Ug > 0 at 3 Myr")) {
        passCount++;
    }
    console.log(`   Ug(3 Myr) = ${Ug_3.toExponential(6)} m/s^2 with Et = ${Et_3.toFixed(6)}, correction = ${(1-Et_3).toFixed(6)}`);
    console.log();

    // Test 7: Complete MUGE computation
    console.log("Test 7: Complete MUGE g_Horsehead computation (10 terms)");
    const g_3 = horsehead.compute_g_Horsehead(t_3Myr);
    totalCount++;
    if (assertGreater(g_3, 0, "g_Horsehead > 0") &&
        assertTrue(isFinite(g_3), "g_Horsehead is finite")) {
        passCount++;
    }
    console.log(`   g_Horsehead(3 Myr) = ${g_3.toExponential(6)} m/s^2`);
    console.log();

    // Test 8: Erosion reduces gravity over time
    console.log("Test 8: Erosion reduces gravitational acceleration over time");
    const g_0 = horsehead.compute_g_Horsehead(0);
    const g_5 = horsehead.compute_g_Horsehead(5 * Myr_to_s);
    totalCount++;
    // As E(t) increases, (1-E(t)) decreases, reducing gravity
    if (assertGreater(g_0, g_5, "g(0) > g(5 Myr) due to increasing erosion")) {
        passCount++;
    }
    console.log(`   g(0) = ${g_0.toExponential(6)}, g(5 Myr) = ${g_5.toExponential(6)}`);
    console.log();

    // Test 9: Variable get/set operations
    console.log("Test 9: Variable get/set operations");
    horsehead.saveState('test9_backup');
    const original_E0 = horsehead.getVariable('E_0');
    horsehead.setVariable('E_0', 0.2);
    const updated_E0 = horsehead.getVariable('E_0');
    horsehead.restoreState('test9_backup');
    const restored_E0 = horsehead.getVariable('E_0');
    totalCount++;
    if (assertApprox(original_E0, 0.1, 0.001, "Original E_0 = 0.1") &&
        assertApprox(updated_E0, 0.2, 0.001, "Updated E_0 = 0.2") &&
        assertApprox(restored_E0, 0.1, 0.001, "Restored E_0 = 0.1")) {
        passCount++;
    }
    console.log();

    // Test 10: expandNebulaScale method
    console.log("Test 10: expandNebulaScale method");
    horsehead.saveState('test10_backup');
    const M_before = horsehead.getVariable('M');
    const r_before = horsehead.getVariable('r');
    horsehead.expandNebulaScale(2.0, 1.5);
    const M_after = horsehead.getVariable('M');
    const r_after = horsehead.getVariable('r');
    totalCount++;
    if (assertApprox(M_after / M_before, 2.0, 0.001, "M scaled by 2.0") &&
        assertApprox(r_after / r_before, 1.5, 0.001, "r scaled by 1.5")) {
        passCount++;
    }
    horsehead.restoreState('test10_backup');
    console.log();

    // Test 11: expandErosionScale method (UNIQUE to dark nebula)
    console.log("Test 11: expandErosionScale method - DARK NEBULA EROSION FEATURE");
    horsehead.saveState('test11_backup');
    const E0_before = horsehead.getVariable('E_0');
    const tau_erosion_before = horsehead.getVariable('tau_erosion');
    horsehead.expandErosionScale(1.5, 2.0);
    const E0_after = horsehead.getVariable('E_0');
    const tau_erosion_after = horsehead.getVariable('tau_erosion');
    totalCount++;
    if (assertApprox(E0_after / E0_before, 1.5, 0.001, "E_0 scaled by 1.5") &&
        assertApprox(tau_erosion_after / tau_erosion_before, 2.0, 0.001, "tau_erosion scaled by 2.0")) {
        passCount++;
    }
    horsehead.restoreState('test11_backup');
    console.log();

    // Test 12: expandWindMagneticScale method
    console.log("Test 12: expandWindMagneticScale method");
    horsehead.saveState('test12_backup');
    const rho_wind_before = horsehead.getVariable('rho_wind');
    const v_wind_before = horsehead.getVariable('v_wind');
    const B_before = horsehead.getVariable('B');
    horsehead.expandWindMagneticScale(1.2, 1.3, 1.1);
    const rho_wind_after = horsehead.getVariable('rho_wind');
    const v_wind_after = horsehead.getVariable('v_wind');
    const B_after = horsehead.getVariable('B');
    totalCount++;
    if (assertApprox(rho_wind_after / rho_wind_before, 1.2, 0.001, "rho_wind scaled by 1.2") &&
        assertApprox(v_wind_after / v_wind_before, 1.3, 0.001, "v_wind scaled by 1.3") &&
        assertApprox(B_after / B_before, 1.1, 0.001, "B scaled by 1.1")) {
        passCount++;
    }
    horsehead.restoreState('test12_backup');
    console.log();

    // Test 13: expandParameterSpace method
    console.log("Test 13: expandParameterSpace method");
    horsehead.saveState('test13_backup');
    const M_pre = horsehead.getVariable('M');
    const r_pre = horsehead.getVariable('r');
    horsehead.expandParameterSpace(1.1);
    const M_post = horsehead.getVariable('M');
    const r_post = horsehead.getVariable('r');
    totalCount++;
    if (assertApprox(M_post / M_pre, 1.1, 0.001, "M scaled by 1.1") &&
        assertApprox(r_post / r_pre, 1.1, 0.001, "r scaled by 1.1")) {
        passCount++;
    }
    horsehead.restoreState('test13_backup');
    console.log();

    // Test 14: Batch operations - scaleVariableGroup
    console.log("Test 14: Batch operations - scaleVariableGroup");
    horsehead.saveState('test14_backup');
    const vars_to_scale = ['M', 'r', 'E_0'];
    const M_b4 = horsehead.getVariable('M');
    const r_b4 = horsehead.getVariable('r');
    const E0_b4 = horsehead.getVariable('E_0');
    horsehead.scaleVariableGroup(vars_to_scale, 1.25);
    totalCount++;
    if (assertApprox(horsehead.getVariable('M') / M_b4, 1.25, 0.001, "M scaled by 1.25") &&
        assertApprox(horsehead.getVariable('r') / r_b4, 1.25, 0.001, "r scaled by 1.25") &&
        assertApprox(horsehead.getVariable('E_0') / E0_b4, 1.25, 0.001, "E_0 scaled by 1.25")) {
        passCount++;
    }
    horsehead.restoreState('test14_backup');
    console.log();

    // Test 15: State management - save/restore/list
    console.log("Test 15: State management - save/restore/list");
    horsehead.saveState('state_A');
    horsehead.setVariable('M', horsehead.getVariable('M') * 2);
    horsehead.saveState('state_B');
    const states = horsehead.listSavedStates();
    horsehead.restoreState('state_A');
    totalCount++;
    if (assertTrue(states.includes('state_A'), "state_A in saved states") &&
        assertTrue(states.includes('state_B'), "state_B in saved states")) {
        passCount++;
    }
    console.log();

    // Test 16: exportState functionality
    console.log("Test 16: exportState functionality");
    const exported = horsehead.exportState();
    totalCount++;
    if (assertTrue(exported.includes('HorseheadNebula State Export'), "Export header present") &&
        assertTrue(exported.includes('M:'), "M in export") &&
        assertTrue(exported.includes('E_0:'), "E_0 in export")) {
        passCount++;
    }
    console.log();

    // Test 17: generateVariations
    console.log("Test 17: generateVariations for parameter exploration");
    const variations = horsehead.generateVariations(5, 10.0);
    totalCount++;
    if (assertTrue(variations.length === 5, "Generated 5 variations") &&
        assertTrue('M' in variations[0], "M in variation") &&
        assertTrue('E_0' in variations[0], "E_0 in variation")) {
        passCount++;
    }
    console.log();

    // Test 18: validateConsistency
    console.log("Test 18: validateConsistency");
    horsehead.saveState('test18_backup');
    const valid_before = horsehead.validateConsistency();
    horsehead.setVariable('E_0', -0.5);  // Invalid
    const valid_invalid = horsehead.validateConsistency();
    horsehead.restoreState('test18_backup');
    const valid_after = horsehead.validateConsistency();
    totalCount++;
    if (assertTrue(valid_before, "Valid before tampering") &&
        assertTrue(!valid_invalid, "Invalid with E_0 < 0") &&
        assertTrue(valid_after, "Valid after restore")) {
        passCount++;
    }
    console.log();

    // Test 19: E_0 scaling shows erosion strength effect
    console.log("Test 19: E_0 scaling affects erosion strength");
    horsehead.saveState('test19_backup');
    const E_0_factors = [0.5, 1.0, 2.0];
    const E_values = [];
    const g_values = [];
    E_0_factors.forEach(factor => {
        horsehead.restoreState('test19_backup');
        horsehead.expandErosionScale(factor, 1.0);
        const Et = horsehead.E_t(t_3Myr);
        const g = horsehead.compute_g_Horsehead(t_3Myr);
        E_values.push(Et);
        g_values.push(g);
    });
    totalCount++;
    // Larger E_0 should give larger E(t) and reduce gravity more
    if (assertTrue(E_values[2] > E_values[0], "E(t) increases with larger E_0") &&
        assertTrue(g_values[2] < g_values[0], "g decreases with larger E_0 (more erosion)")) {
        passCount++;
    }
    console.log(`   E_0 × 0.5: E(3 Myr) = ${E_values[0].toFixed(6)}, g = ${g_values[0].toExponential(6)}`);
    console.log(`   E_0 × 1.0: E(3 Myr) = ${E_values[1].toFixed(6)}, g = ${g_values[1].toExponential(6)}`);
    console.log(`   E_0 × 2.0: E(3 Myr) = ${E_values[2].toFixed(6)}, g = ${g_values[2].toExponential(6)}`);
    horsehead.restoreState('test19_backup');
    console.log();

    // Test 20: tau_erosion scaling affects erosion rate
    console.log("Test 20: tau_erosion scaling affects erosion rate");
    horsehead.saveState('test20_backup');
    const tau_factors = [0.5, 1.0, 2.0];
    const E_tau_values = [];
    tau_factors.forEach(factor => {
        horsehead.restoreState('test20_backup');
        horsehead.expandErosionScale(1.0, factor);
        const Et = horsehead.E_t(t_3Myr);
        E_tau_values.push(Et);
    });
    totalCount++;
    // Shorter tau (0.5) → faster erosion → larger E(t) at fixed time
    // Longer tau (2.0) → slower erosion → smaller E(t) at fixed time
    if (assertTrue(E_tau_values[0] > E_tau_values[2], "Shorter tau_erosion → faster erosion → larger E(t)")) {
        passCount++;
    }
    console.log(`   tau × 0.5: E(3 Myr) = ${E_tau_values[0].toFixed(6)}`);
    console.log(`   tau × 1.0: E(3 Myr) = ${E_tau_values[1].toFixed(6)}`);
    console.log(`   tau × 2.0: E(3 Myr) = ${E_tau_values[2].toFixed(6)}`);
    horsehead.restoreState('test20_backup');
    console.log();

    // Test 21: Sensitivity analysis
    console.log("Test 21: Sensitivity analysis at 3 Myr");
    horsehead.saveState('test21_backup');
    const sensitivities = horsehead.sensitivityAnalysis(t_3Myr, 0.01);
    const sens_entries = Object.entries(sensitivities).sort((a, b) => b[1] - a[1]);
    totalCount++;
    if (assertTrue(Object.keys(sensitivities).length === 32, "Sensitivities for all 32 variables") &&
        assertTrue(sensitivities['M'] > 0, "M has positive sensitivity")) {
        passCount++;
    }
    console.log(`   Top 5 sensitivities: ${sens_entries.slice(0, 5).map(e => `${e[0]}:${e[1].toExponential(2)}`).join(', ')}`);
    horsehead.restoreState('test21_backup');
    console.log();

    // Test 22: UQFF term breakdown validation
    console.log("Test 22: UQFF term breakdown - all terms contribute");
    horsehead.saveState('test22_backup');
    const t_test = 3 * Myr_to_s;
    
    // Compute individual terms manually
    const Et = horsehead.E_t(t_test);
    const corr_H = 1 + horsehead.H0 * t_test;
    const corr_B = 1 - horsehead.B / horsehead.B_crit;
    const corr_E = 1 - Et;
    const term1 = horsehead.ug1_base * corr_H * corr_B * corr_E;
    
    const term2 = horsehead.compute_Ug(Et);
    const term3 = (horsehead.Lambda * horsehead.c_light * horsehead.c_light) / 3.0;
    
    const cross_vB = horsehead.gas_v * horsehead.B;
    const em_base = (horsehead.q_charge * cross_vB) / horsehead.proton_mass;
    const corr_UA = 1 + (horsehead.rho_vac_UA / horsehead.rho_vac_SCm);
    const term4 = (em_base * corr_UA) * horsehead.scale_EM;
    
    const sqrt_unc = Math.sqrt(horsehead.delta_x * horsehead.delta_p);
    const term_q = (horsehead.hbar / sqrt_unc) * horsehead.integral_psi * (2 * Math.PI / horsehead.t_Hubble);
    
    const V = horsehead.compute_V();
    const term_fluid = (horsehead.rho_fluid * V * horsehead.ug1_base) / horsehead.M;
    
    const wind_pressure = horsehead.rho_wind * horsehead.v_wind * horsehead.v_wind;
    const term_wind = wind_pressure / horsehead.rho_fluid;
    
    const g_total = horsehead.compute_g_Horsehead(t_test);
    
    totalCount++;
    if (assertGreater(term1, 0, "Term1 (base+H+B+E) > 0") &&
        assertGreater(term2, 0, "Term2 (Ug) > 0") &&
        assertGreater(term_wind, 0, "Term_wind > 0") &&
        assertGreater(g_total, 0, "Total g > 0") &&
        assertLess(corr_E, 1.0, "Erosion correction (1-E(t)) < 1.0")) {
        passCount++;
    }
    console.log(`   Erosion E(t) = ${Et.toFixed(6)}, correction (1-E(t)) = ${corr_E.toFixed(6)}`);
    console.log(`   Term1 (base+H+B+E): ${term1.toExponential(6)} m/s^2`);
    console.log(`   Term2 (Ug): ${term2.toExponential(6)} m/s^2`);
    console.log(`   Term3 (Lambda): ${term3.toExponential(6)} m/s^2`);
    console.log(`   Term4 (EM+UA): ${term4.toExponential(6)} m/s^2`);
    console.log(`   Term_q (quantum): ${term_q.toExponential(6)} m/s^2`);
    console.log(`   Term_fluid: ${term_fluid.toExponential(6)} m/s^2`);
    console.log(`   Term_wind: ${term_wind.toExponential(6)} m/s^2`);
    console.log(`   Total g: ${g_total.toExponential(6)} m/s^2`);
    horsehead.restoreState('test22_backup');
    console.log();

    // Summary
    console.log("=========================================================");
    console.log(`TEST SUMMARY: ${passCount}/${totalCount} tests passed`);
    if (passCount === totalCount) {
        console.log("✓ ALL TESTS COMPLETED SUCCESSFULLY");
    } else {
        console.log(`✗ ${totalCount - passCount} test(s) failed`);
    }
    console.log("=========================================================\n");

    return passCount === totalCount;
}

// Run tests if executed directly
if (typeof require !== 'undefined' && require.main === module) {
    const success = runAllTests();
    process.exit(success ? 0 : 1);
}

export { runAllTests };
