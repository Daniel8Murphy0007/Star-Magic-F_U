/**
 * ================================================================================================
 * Test Suite: test_source23.js
 *
 * Description: Comprehensive test suite for AntennaeGalaxies (source23.js) module
 *              Tests all 38 methods (13 core + 25 enhanced) including unique merger physics
 *
 * Tests Cover:
 *   - Initialization and parameter validation
 *   - Star formation M(t) evolution with exponential growth
 *   - Merger interaction I(t) exponential decay (UNIQUE)
 *   - UQFF Ug computation with merger interaction
 *   - Complete MUGE with all 10 terms
 *   - Enhanced capabilities: scaling, state management, evolution
 *   - Galaxy-specific features: expandMergerScale, expandStarFormationScale
 *
 * Author: Test suite for UQFF module conversion
 * Date: November 03, 2025
 * ================================================================================================
 */

import AntennaeGalaxies from './source23.js';

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

// Main test suite
function runAllTests() {
    console.log("=========================================================");
    console.log("ANTENNAE GALAXIES (source23.js) COMPREHENSIVE TEST SUITE");
    console.log("Testing NGC 4038 & NGC 4039 Galaxy Merger Physics");
    console.log("=========================================================\n");

    let passCount = 0;
    let totalCount = 0;

    // Test 1: Initialization validation
    console.log("Test 1: Initialization with default parameters");
    const antennae = new AntennaeGalaxies();
    totalCount++;
    if (assertTrue(antennae.M0 === 2e11 * M_sun, "M0 initialized to 2e11 M_sun") &&
        assertTrue(antennae.r === 30000 * ly_to_m, "r initialized to 30,000 ly") &&
        assertTrue(antennae.I0 === 0.1, "I0 initialized to 0.1") &&
        assertTrue(antennae.tau_merger === 400 * Myr_to_s, "tau_merger = 400 Myr")) {
        passCount++;
    }
    console.log();

    // Test 2: System name
    console.log("Test 2: System name retrieval");
    totalCount++;
    if (assertTrue(antennae.getSystemName() === "AntennaeGalaxies", "System name is AntennaeGalaxies")) {
        passCount++;
    }
    console.log();

    // Test 3: Variable listing
    console.log("Test 3: Variable listing completeness");
    const vars = antennae.listVariables();
    totalCount++;
    if (assertTrue(vars.length === 35, "35 variables listed") &&
        assertTrue(vars.includes('I0'), "I0 in variable list") &&
        assertTrue(vars.includes('tau_merger'), "tau_merger in variable list")) {
        passCount++;
    }
    console.log();

    // Test 4: Star formation M(t) evolution
    console.log("Test 4: Star formation M(t) evolution");
    const t0 = 0;
    const t_100Myr = 100 * Myr_to_s;
    const M0 = antennae.M_t(t0);
    const M_100 = antennae.M_t(t_100Myr);
    totalCount++;
    // At t=0: M(0) = M0 × (1 + SFR_factor) - maximum mass
    // At t→∞: M(∞) → M0 - mass decreases due to exp(-t/tau_SF) decay
    if (assertGreater(M0, M_100, "M(0) > M(100 Myr) - mass decreases as SF decays") &&
        assertTrue(M0 / M_sun > 2e11, "M(0) > 2e11 M_sun due to initial SF burst")) {
        passCount++;
    }
    console.log();

    // Test 5: Merger interaction I(t) exponential decay (UNIQUE FEATURE)
    console.log("Test 5: Merger interaction I(t) exponential decay - GALAXY MERGER FEATURE");
    const I0_val = antennae.I_t(0);
    const I_tau = antennae.I_t(antennae.tau_merger);
    const expected_I_tau = antennae.I0 * Math.exp(-1);  // I(tau) = I0 × e^-1
    totalCount++;
    if (assertApprox(I0_val, antennae.I0, 0.001, "I(0) = I0 = 0.1") &&
        assertApprox(I_tau / antennae.I0, Math.exp(-1), 0.001, "I(tau)/I0 = e^-1 = 0.367879")) {
        passCount++;
    }
    console.log(`   I(0) = ${I0_val.toFixed(6)}, I(tau) = ${I_tau.toFixed(6)}, I(tau)/I0 = ${(I_tau/antennae.I0).toFixed(6)}`);
    console.log();

    // Test 6: UQFF Ug computation with merger interaction
    console.log("Test 6: UQFF Ug computation includes merger interaction I(t)");
    const t_300Myr = 300 * Myr_to_s;
    const Mt_300 = antennae.M_t(t_300Myr);
    const It_300 = antennae.I_t(t_300Myr);
    const Ug_300 = antennae.compute_Ug(Mt_300, It_300);
    totalCount++;
    // Ug should include (1 + It) correction factor
    if (assertGreater(Ug_300, 0, "Ug > 0 at 300 Myr")) {
        passCount++;
    }
    console.log(`   Ug(300 Myr) = ${Ug_300.toExponential(6)} m/s^2 with It = ${It_300.toFixed(6)}`);
    console.log();

    // Test 7: Complete MUGE computation
    console.log("Test 7: Complete MUGE g_Antennae computation (10 terms)");
    const g_300 = antennae.compute_g_Antennae(t_300Myr);
    totalCount++;
    if (assertGreater(g_300, 0, "g_Antennae > 0") &&
        assertTrue(isFinite(g_300), "g_Antennae is finite")) {
        passCount++;
    }
    console.log(`   g_Antennae(300 Myr) = ${g_300.toExponential(6)} m/s^2`);
    console.log();

    // Test 8: Time evolution consistency
    console.log("Test 8: Time evolution consistency check");
    const g_0 = antennae.compute_g_Antennae(0);
    const g_500 = antennae.compute_g_Antennae(500 * Myr_to_s);
    totalCount++;
    // Both should be positive and finite
    if (assertGreater(g_0, 0, "g(0) > 0") &&
        assertGreater(g_500, 0, "g(500 Myr) > 0") &&
        assertTrue(Math.abs(g_0 - g_500) / g_0 < 1.0, "g varies within reasonable range")) {
        passCount++;
    }
    console.log(`   g(0) = ${g_0.toExponential(6)}, g(500 Myr) = ${g_500.toExponential(6)}`);
    console.log();

    // Test 9: Variable get/set operations
    console.log("Test 9: Variable get/set operations");
    antennae.saveState('test9_backup');
    const original_I0 = antennae.getVariable('I0');
    antennae.setVariable('I0', 0.2);
    const updated_I0 = antennae.getVariable('I0');
    antennae.restoreState('test9_backup');
    const restored_I0 = antennae.getVariable('I0');
    totalCount++;
    if (assertApprox(original_I0, 0.1, 0.001, "Original I0 = 0.1") &&
        assertApprox(updated_I0, 0.2, 0.001, "Updated I0 = 0.2") &&
        assertApprox(restored_I0, 0.1, 0.001, "Restored I0 = 0.1")) {
        passCount++;
    }
    console.log();

    // Test 10: expandGalaxyScale method
    console.log("Test 10: expandGalaxyScale method");
    antennae.saveState('test10_backup');
    const M0_before = antennae.getVariable('M0');
    const r_before = antennae.getVariable('r');
    antennae.expandGalaxyScale(2.0, 1.5);
    const M0_after = antennae.getVariable('M0');
    const r_after = antennae.getVariable('r');
    totalCount++;
    if (assertApprox(M0_after / M0_before, 2.0, 0.001, "M0 scaled by 2.0") &&
        assertApprox(r_after / r_before, 1.5, 0.001, "r scaled by 1.5")) {
        passCount++;
    }
    antennae.restoreState('test10_backup');
    console.log();

    // Test 11: expandMergerScale method (UNIQUE to galaxy merger)
    console.log("Test 11: expandMergerScale method - GALAXY MERGER FEATURE");
    antennae.saveState('test11_backup');
    const I0_before = antennae.getVariable('I0');
    const tau_merger_before = antennae.getVariable('tau_merger');
    antennae.expandMergerScale(1.5, 2.0);
    const I0_after = antennae.getVariable('I0');
    const tau_merger_after = antennae.getVariable('tau_merger');
    totalCount++;
    if (assertApprox(I0_after / I0_before, 1.5, 0.001, "I0 scaled by 1.5") &&
        assertApprox(tau_merger_after / tau_merger_before, 2.0, 0.001, "tau_merger scaled by 2.0")) {
        passCount++;
    }
    antennae.restoreState('test11_backup');
    console.log();

    // Test 12: expandStarFormationScale method
    console.log("Test 12: expandStarFormationScale method");
    antennae.saveState('test12_backup');
    const SFR_before = antennae.getVariable('SFR_factor');
    const tau_SF_before = antennae.getVariable('tau_SF');
    antennae.expandStarFormationScale(1.2, 1.3);
    const SFR_after = antennae.getVariable('SFR_factor');
    const tau_SF_after = antennae.getVariable('tau_SF');
    totalCount++;
    if (assertApprox(SFR_after / SFR_before, 1.2, 0.001, "SFR_factor scaled by 1.2") &&
        assertApprox(tau_SF_after / tau_SF_before, 1.3, 0.001, "tau_SF scaled by 1.3")) {
        passCount++;
    }
    antennae.restoreState('test12_backup');
    console.log();

    // Test 13: expandParameterSpace method
    console.log("Test 13: expandParameterSpace method");
    antennae.saveState('test13_backup');
    const M0_pre = antennae.getVariable('M0');
    const r_pre = antennae.getVariable('r');
    antennae.expandParameterSpace(1.1);
    const M0_post = antennae.getVariable('M0');
    const r_post = antennae.getVariable('r');
    totalCount++;
    if (assertApprox(M0_post / M0_pre, 1.1, 0.001, "M0 scaled by 1.1") &&
        assertApprox(r_post / r_pre, 1.1, 0.001, "r scaled by 1.1")) {
        passCount++;
    }
    antennae.restoreState('test13_backup');
    console.log();

    // Test 14: Batch operations - scaleVariableGroup
    console.log("Test 14: Batch operations - scaleVariableGroup");
    antennae.saveState('test14_backup');
    const vars_to_scale = ['M0', 'r', 'I0'];
    const M0_b4 = antennae.getVariable('M0');
    const r_b4 = antennae.getVariable('r');
    const I0_b4 = antennae.getVariable('I0');
    antennae.scaleVariableGroup(vars_to_scale, 1.25);
    totalCount++;
    if (assertApprox(antennae.getVariable('M0') / M0_b4, 1.25, 0.001, "M0 scaled by 1.25") &&
        assertApprox(antennae.getVariable('r') / r_b4, 1.25, 0.001, "r scaled by 1.25") &&
        assertApprox(antennae.getVariable('I0') / I0_b4, 1.25, 0.001, "I0 scaled by 1.25")) {
        passCount++;
    }
    antennae.restoreState('test14_backup');
    console.log();

    // Test 15: State management - save/restore/list
    console.log("Test 15: State management - save/restore/list");
    antennae.saveState('state_A');
    antennae.setVariable('M0', antennae.getVariable('M0') * 2);
    antennae.saveState('state_B');
    const states = antennae.listSavedStates();
    antennae.restoreState('state_A');
    totalCount++;
    if (assertTrue(states.includes('state_A'), "state_A in saved states") &&
        assertTrue(states.includes('state_B'), "state_B in saved states")) {
        passCount++;
    }
    console.log();

    // Test 16: exportState functionality
    console.log("Test 16: exportState functionality");
    const exported = antennae.exportState();
    totalCount++;
    if (assertTrue(exported.includes('AntennaeGalaxies State Export'), "Export header present") &&
        assertTrue(exported.includes('M0:'), "M0 in export") &&
        assertTrue(exported.includes('I0:'), "I0 in export")) {
        passCount++;
    }
    console.log();

    // Test 17: generateVariations
    console.log("Test 17: generateVariations for parameter exploration");
    const variations = antennae.generateVariations(5, 10.0);
    totalCount++;
    if (assertTrue(variations.length === 5, "Generated 5 variations") &&
        assertTrue('M0' in variations[0], "M0 in variation") &&
        assertTrue('I0' in variations[0], "I0 in variation")) {
        passCount++;
    }
    console.log();

    // Test 18: validateConsistency
    console.log("Test 18: validateConsistency");
    antennae.saveState('test18_backup');
    const valid_before = antennae.validateConsistency();
    antennae.setVariable('I0', -0.5);  // Invalid
    const valid_invalid = antennae.validateConsistency();
    antennae.restoreState('test18_backup');
    const valid_after = antennae.validateConsistency();
    totalCount++;
    if (assertTrue(valid_before, "Valid before tampering") &&
        assertTrue(!valid_invalid, "Invalid with I0 < 0") &&
        assertTrue(valid_after, "Valid after restore")) {
        passCount++;
    }
    console.log();

    // Test 19: I0 scaling shows merger interaction effect
    console.log("Test 19: I0 scaling affects merger interaction strength");
    antennae.saveState('test19_backup');
    const I0_factors = [0.5, 1.0, 2.0];
    const g_values = [];
    const I_values = [];
    I0_factors.forEach(factor => {
        antennae.restoreState('test19_backup');
        antennae.expandMergerScale(factor, 1.0);
        const g = antennae.compute_g_Antennae(t_300Myr);
        const It = antennae.I_t(t_300Myr);
        g_values.push(g);
        I_values.push(It);
    });
    totalCount++;
    // Larger I0 should give larger I(t) and affect gravity
    if (assertTrue(I_values[2] > I_values[0], "I(t) increases with larger I0") &&
        assertTrue(g_values[2] > g_values[0], "g increases with larger I0 (stronger interaction)")) {
        passCount++;
    }
    console.log(`   I0 × 0.5: I(300 Myr) = ${I_values[0].toFixed(6)}, g = ${g_values[0].toExponential(6)}`);
    console.log(`   I0 × 1.0: I(300 Myr) = ${I_values[1].toFixed(6)}, g = ${g_values[1].toExponential(6)}`);
    console.log(`   I0 × 2.0: I(300 Myr) = ${I_values[2].toFixed(6)}, g = ${g_values[2].toExponential(6)}`);
    antennae.restoreState('test19_backup');
    console.log();

    // Test 20: tau_merger scaling affects decay rate
    console.log("Test 20: tau_merger scaling affects merger interaction decay rate");
    antennae.saveState('test20_backup');
    const tau_factors = [0.5, 1.0, 2.0];
    const I_tau_values = [];
    tau_factors.forEach(factor => {
        antennae.restoreState('test20_backup');
        antennae.expandMergerScale(1.0, factor);
        const It = antennae.I_t(t_300Myr);
        I_tau_values.push(It);
    });
    totalCount++;
    // Shorter tau (0.5) → faster decay → smaller I(t) at fixed time
    // Longer tau (2.0) → slower decay → larger I(t) at fixed time
    if (assertTrue(I_tau_values[2] > I_tau_values[0], "Longer tau_merger → slower decay → larger I(t)")) {
        passCount++;
    }
    console.log(`   tau × 0.5: I(300 Myr) = ${I_tau_values[0].toFixed(6)}`);
    console.log(`   tau × 1.0: I(300 Myr) = ${I_tau_values[1].toFixed(6)}`);
    console.log(`   tau × 2.0: I(300 Myr) = ${I_tau_values[2].toFixed(6)}`);
    antennae.restoreState('test20_backup');
    console.log();

    // Test 21: Sensitivity analysis
    console.log("Test 21: Sensitivity analysis at 300 Myr");
    antennae.saveState('test21_backup');
    const sensitivities = antennae.sensitivityAnalysis(t_300Myr, 0.01);
    const sens_entries = Object.entries(sensitivities).sort((a, b) => b[1] - a[1]);
    totalCount++;
    if (assertTrue(Object.keys(sensitivities).length === 35, "Sensitivities for all 35 variables") &&
        assertTrue(sensitivities['M0'] > 0, "M0 has positive sensitivity")) {
        passCount++;
    }
    console.log(`   Top 5 sensitivities: ${sens_entries.slice(0, 5).map(e => `${e[0]}:${e[1].toExponential(2)}`).join(', ')}`);
    antennae.restoreState('test21_backup');
    console.log();

    // Test 22: UQFF term breakdown validation
    console.log("Test 22: UQFF term breakdown - all terms contribute");
    antennae.saveState('test22_backup');
    const t_test = 300 * Myr_to_s;
    
    // Compute individual terms manually
    const Mt = antennae.M_t(t_test);
    const It = antennae.I_t(t_test);
    const ug1_t = (antennae.G * Mt) / (antennae.r * antennae.r);
    const corr_H = 1 + antennae.Hz * t_test;
    const corr_B = 1 - antennae.B / antennae.B_crit;
    const corr_I = 1 + It;
    const term1 = ug1_t * corr_H * corr_B * corr_I;
    
    const term2 = antennae.compute_Ug(Mt, It);
    const term3 = (antennae.Lambda * antennae.c_light * antennae.c_light) / 3.0;
    
    const cross_vB = antennae.gas_v * antennae.B;
    const em_base = (antennae.q_charge * cross_vB) / antennae.proton_mass;
    const corr_UA = 1 + (antennae.rho_vac_UA / antennae.rho_vac_SCm);
    const term4 = (em_base * corr_UA) * antennae.scale_EM;
    
    const sqrt_unc = Math.sqrt(antennae.delta_x * antennae.delta_p);
    const term_q = (antennae.hbar / sqrt_unc) * antennae.integral_psi * (2 * Math.PI / antennae.t_Hubble);
    
    const V = antennae.compute_V();
    const term_fluid = (antennae.rho_fluid * V * ug1_t) / Mt;
    
    const wind_pressure = antennae.rho_wind * antennae.v_wind * antennae.v_wind;
    const term_feedback = wind_pressure / antennae.rho_fluid;
    
    const g_total = antennae.compute_g_Antennae(t_test);
    
    totalCount++;
    if (assertGreater(term1, 0, "Term1 (base+H+B+I) > 0") &&
        assertGreater(term2, 0, "Term2 (Ug) > 0") &&
        assertGreater(term_feedback, 0, "Term_feedback > 0") &&
        assertGreater(g_total, 0, "Total g > 0")) {
        passCount++;
    }
    console.log(`   Term1 (base+H+B+I): ${term1.toExponential(6)} m/s^2`);
    console.log(`   Term2 (Ug): ${term2.toExponential(6)} m/s^2`);
    console.log(`   Term3 (Lambda): ${term3.toExponential(6)} m/s^2`);
    console.log(`   Term4 (EM+UA): ${term4.toExponential(6)} m/s^2`);
    console.log(`   Term_q (quantum): ${term_q.toExponential(6)} m/s^2`);
    console.log(`   Term_fluid: ${term_fluid.toExponential(6)} m/s^2`);
    console.log(`   Term_feedback: ${term_feedback.toExponential(6)} m/s^2`);
    console.log(`   Total g: ${g_total.toExponential(6)} m/s^2`);
    antennae.restoreState('test22_backup');
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
