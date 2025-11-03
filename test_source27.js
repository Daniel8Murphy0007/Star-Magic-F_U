/**
 * ================================================================================================
 * Comprehensive Test Suite for GalaxyNGC1792 (NGC 1792 "Stellar Forge" Starburst Galaxy)
 * 
 * Tests all physics components, time evolution, starburst features, and UQFF framework
 * ================================================================================================
 */

import GalaxyNGC1792 from './source27.js';

// Constants
const M_sun = 1.989e30;
const ly_to_m = 9.461e15;
const Myr_to_s = 1e6 * 3.156e7;
const Gyr_to_s = 1e9 * 3.156e7;

// Test utilities
let testCount = 0;
let passCount = 0;

function assert(condition, testName) {
    testCount++;
    if (condition) {
        passCount++;
        console.log(`✓ Test ${testCount}: ${testName}`);
        return true;
    } else {
        console.log(`✗ Test ${testCount}: ${testName}`);
        return false;
    }
}

function assertClose(a, b, tolerance, testName) {
    const relDiff = Math.abs(a - b) / Math.abs(b);
    return assert(relDiff < tolerance, `${testName} (rel_diff: ${relDiff.toExponential(2)})`);
}

console.log("=========================================================");
console.log("NGC 1792 COMPREHENSIVE TEST SUITE");
console.log("=========================================================\n");

const ngc = new GalaxyNGC1792();

// ========== TEST GROUP 1: BASIC INITIALIZATION (Tests 1-3) ==========
console.log("Test Group 1: Basic Initialization");
console.log("-----------------------------------");

// Test 1: M0 and r within expected ranges
const M0_Msun = ngc.M0 / M_sun;
const r_ly = ngc.r / ly_to_m;
assert(M0_Msun > 1e9 && M0_Msun < 1e11 && r_ly > 1e4 && r_ly < 1e6,
    "M0 and r within reasonable galaxy ranges");

// Test 2: Hz > 0 from cosmological model
assert(ngc.Hz > 0 && ngc.Hz < 1e-16, "Hz > 0 and reasonable for z = 0.0095");

// Test 3: Static magnetic field B = 1×10⁻⁵ T
assertClose(ngc.B, 1e-5, 1e-6, "B = 1e-5 T (static starburst field)");

console.log();

// ========== TEST GROUP 2: STARBURST MASS GROWTH M(t) (Tests 4-6) ==========
console.log("Test Group 2: Starburst Mass Growth M(t)");
console.log("------------------------------------------");

// Test 4: M(0) = M0 × (1 + SFR_factor) - initial mass boost
const M_at_0 = ngc.M_t(0);
const expected_M0_boost = ngc.M0 * (1 + ngc.SFR_factor);
assertClose(M_at_0, expected_M0_boost, 1e-9, "M(0) = M0 × (1 + SFR_factor)");

// Test 5: Monotonic decrease check: M(0) ≥ M(25 Myr) ≥ M(50 Myr) ≥ M(100 Myr)
const M_0 = ngc.M_t(0);
const M_25 = ngc.M_t(25 * Myr_to_s);
const M_50 = ngc.M_t(50 * Myr_to_s);
const M_100 = ngc.M_t(100 * Myr_to_s);
assert(M_0 >= M_25 && M_25 >= M_50 && M_50 >= M_100,
    "M(t) monotonically decreases: M(0) ≥ M(25) ≥ M(50) ≥ M(100 Myr)");

// Test 6: M(1 Gyr) ≈ M0 (asymptotic approach, SFR term decayed)
const M_1Gyr = ngc.M_t(1.0 * Gyr_to_s);
const decay_at_1Gyr = Math.exp(-(1.0 * Gyr_to_s) / ngc.tau_SF);
const expected_M_1Gyr = ngc.M0 * (1 + ngc.SFR_factor * decay_at_1Gyr);
assertClose(M_1Gyr, expected_M_1Gyr, 1e-9, "M(1 Gyr) ≈ M0 × (1 + SFR_factor × e^(-10)) ≈ M0");

console.log();

// ========== TEST GROUP 3: UQFF COMPONENTS (Tests 7-9) ==========
console.log("Test Group 3: UQFF Components");
console.log("-------------------------------");

// Test 7: UQFF Ug includes (1 + f_TRZ) correction
const Mt_50 = ngc.M_t(50 * Myr_to_s);
const Ug_50 = ngc.compute_Ug(Mt_50);
const ug1_raw = (ngc.G * Mt_50) / (ngc.r * ngc.r);
const corr_B = 1 - ngc.B / ngc.B_crit;
const ug4_raw = ug1_raw * corr_B;
const expected_Ug = (ug1_raw + ug4_raw) * (1 + ngc.f_TRZ);
assertClose(Ug_50, expected_Ug, 1e-9, "UQFF Ug = (Ug1 + Ug4) × (1 + f_TRZ)");

// Test 8: g_NGC1792(t) > 0 for all times tested
const times_to_test = [0, 25, 50, 100, 200, 500].map(t => t * Myr_to_s);
const all_positive = times_to_test.every(t => ngc.compute_g_NGC1792(t) > 0);
assert(all_positive, "g_NGC1792(t) > 0 for all t ∈ [0, 500 Myr]");

// Test 9: Negative time returns 0
const g_neg = ngc.compute_g_NGC1792(-1);
assert(g_neg === 0, "g_NGC1792(t < 0) returns 0");

console.log();

// ========== TEST GROUP 4: STARBURST SCALING METHODS (Tests 10-12) ==========
console.log("Test Group 4: Starburst Scaling Methods (UNIQUE)");
console.log("--------------------------------------------------");

// Test 10: expandGalaxyScale(M0, r)
ngc.saveState('test10');
const M0_orig = ngc.M0;
const r_orig = ngc.r;
ngc.expandGalaxyScale(1.5, 2.0);
const M0_new = ngc.M0;
const r_new = ngc.r;
const test10_pass = assertClose(M0_new / M0_orig, 1.5, 1e-9, "expandGalaxyScale(1.5, 2.0) scales M0") &&
                    assertClose(r_new / r_orig, 2.0, 1e-9, "expandGalaxyScale(1.5, 2.0) scales r");
ngc.restoreState('test10');

// Test 11: expandStarburstScale(SFR_factor, tau_SF) - UNIQUE to starburst
ngc.saveState('test11');
const SFR_orig = ngc.SFR_factor;
const tau_SF_orig = ngc.tau_SF;
ngc.expandStarburstScale(2.0, 1.5);
const SFR_new = ngc.SFR_factor;
const tau_SF_new = ngc.tau_SF;
const test11_pass = assertClose(SFR_new / SFR_orig, 2.0, 1e-9, "expandStarburstScale(2.0, 1.5) scales SFR_factor") &&
                    assertClose(tau_SF_new / tau_SF_orig, 1.5, 1e-9, "expandStarburstScale(2.0, 1.5) scales tau_SF");
ngc.restoreState('test11');

// Test 12: expandWindMagneticScale(rho_wind, v_wind, B) - UNIQUE to supernova starburst
ngc.saveState('test12');
const rho_wind_orig = ngc.rho_wind;
const v_wind_orig = ngc.v_wind;
const B_orig = ngc.B;
ngc.expandWindMagneticScale(1.5, 2.0, 1.2);
const rho_wind_new = ngc.rho_wind;
const v_wind_new = ngc.v_wind;
const B_new = ngc.B;
const test12_pass = assertClose(rho_wind_new / rho_wind_orig, 1.5, 1e-9, "expandWindMagneticScale(1.5, 2.0, 1.2) scales rho_wind") &&
                    assertClose(v_wind_new / v_wind_orig, 2.0, 1e-9, "expandWindMagneticScale(1.5, 2.0, 1.2) scales v_wind") &&
                    assertClose(B_new / B_orig, 1.2, 1e-9, "expandWindMagneticScale(1.5, 2.0, 1.2) scales B");
ngc.restoreState('test12');

console.log();

// ========== TEST GROUP 5: PARAMETER SCALING EFFECTS (Tests 13-18) ==========
console.log("Test Group 5: Parameter Scaling Effects");
console.log("----------------------------------------");

// Test 13: M0 scaling → affects M(t) directly
ngc.saveState('test13');
const t_test = 50 * Myr_to_s;
const M_test_base = ngc.M_t(t_test);
ngc.setVariable('M0', ngc.M0 * 2);
const M_test_scaled = ngc.M_t(t_test);
const M_ratio = M_test_scaled / M_test_base;
assert(Math.abs(M_ratio - 2.0) < 0.01, `M0 × 2 → M(t) × ~2 (ratio: ${M_ratio.toFixed(6)})`);
ngc.restoreState('test13');

// Test 14: r scaling → affects ug1 gravitational term (inverse square)
ngc.saveState('test14');
const Mt_test = ngc.M_t(t_test);
const ug1_r_base = (ngc.G * Mt_test) / (ngc.r * ngc.r);
ngc.setVariable('r', ngc.r * 2);
const ug1_r_scaled = (ngc.G * Mt_test) / (ngc.r * ngc.r);
const ug1_ratio = ug1_r_base / ug1_r_scaled;
assertClose(ug1_ratio, 4.0, 0.01, `r × 2 → ug1 × 1/4 (ratio: ${ug1_ratio.toFixed(6)})`);
ngc.restoreState('test14');

// Test 15: SFR_factor scaling → affects M(t)
ngc.saveState('test15');
const M_SFR_base = ngc.M_t(t_test);
ngc.setVariable('SFR_factor', ngc.SFR_factor * 2);
const M_SFR_scaled = ngc.M_t(t_test);
assert(M_SFR_scaled > M_SFR_base, "SFR_factor × 2 → M(t) increases");
ngc.restoreState('test15');

// Test 16: tau_SF scaling → affects M(t) decay rate
ngc.saveState('test16');
const M_tau_base = ngc.M_t(t_test);
ngc.setVariable('tau_SF', ngc.tau_SF * 2);  // Longer timescale → slower decay
const M_tau_scaled = ngc.M_t(t_test);
assert(M_tau_scaled > M_tau_base, "tau_SF × 2 → M(t) higher (slower decay)");
ngc.restoreState('test16');

// Test 17: rho_wind scaling → affects feedback term
ngc.saveState('test17');
const g_rho_wind_base = ngc.compute_g_NGC1792(t_test);
ngc.setVariable('rho_wind', ngc.rho_wind * 2);
const g_rho_wind_scaled = ngc.compute_g_NGC1792(t_test);
assert(g_rho_wind_scaled > g_rho_wind_base, "rho_wind × 2 → g increases (feedback increases)");
ngc.restoreState('test17');

// Test 18: v_wind scaling → affects feedback term (v² dependence)
ngc.saveState('test18');
const g_v_wind_base = ngc.compute_g_NGC1792(t_test);
ngc.setVariable('v_wind', ngc.v_wind * 2);
const g_v_wind_scaled = ngc.compute_g_NGC1792(t_test);
const g_v_wind_ratio = g_v_wind_scaled / g_v_wind_base;
// Expect ~4× increase in feedback term contribution, but may be dominated by other terms
assert(g_v_wind_ratio > 1.1, `v_wind × 2 → g increases (v² dependence, ratio: ${g_v_wind_ratio.toFixed(6)})`);
ngc.restoreState('test18');

console.log();

// ========== TEST GROUP 6: PHYSICS VALIDATION (Tests 19-21) ==========
console.log("Test Group 6: Physics Validation");
console.log("---------------------------------");

// Test 19: Sensitivity analysis identifies key parameters
ngc.saveState('test19');
const sens = ngc.sensitivityAnalysis(50 * Myr_to_s, 0.01);
const sorted_sens = Object.entries(sens).sort((a, b) => b[1] - a[1]);
const top_3 = sorted_sens.slice(0, 3).map(([k, v]) => k);
console.log(`  Top 3 sensitive parameters: ${top_3.join(', ')}`);
// Expect v_wind, rho_wind, rho_fluid to dominate (feedback term)
const has_wind_params = top_3.some(p => ['v_wind', 'rho_wind', 'rho_fluid'].includes(p));
assert(has_wind_params, "Sensitivity includes wind/feedback parameters");
ngc.restoreState('test19');

// Test 20: Supernova feedback term magnitude (rho_wind × v_wind² / rho_fluid)
const wind_pressure = ngc.rho_wind * ngc.v_wind * ngc.v_wind;
const feedback_term = wind_pressure / ngc.rho_fluid;
console.log(`  Feedback term: ${feedback_term.toExponential(6)} m/s²`);
assert(feedback_term > 1e10, "Feedback term >> 1e10 m/s² (dominates like HUDF)");

// Test 21: Validate consistency check
assert(ngc.validateConsistency(), "validateConsistency() returns true for valid parameters");

console.log();

// ========== TEST GROUP 7: ENHANCED CAPABILITIES (Test 22) ==========
console.log("Test Group 7: Enhanced Capabilities");
console.log("------------------------------------");

// Test 22: State save/restore functionality
ngc.saveState('enhanced_test');
const original_M0 = ngc.M0;
ngc.setVariable('M0', ngc.M0 * 3);
const modified_M0 = ngc.M0;
ngc.restoreState('enhanced_test');
const restored_M0 = ngc.M0;
assert(Math.abs(restored_M0 - original_M0) < 1e-20 && Math.abs(modified_M0 / original_M0 - 3) < 1e-9,
    "State save/restore preserves parameters");

console.log();

// ========== SUMMARY ==========
console.log("=========================================================");
console.log(`TEST SUMMARY: ${passCount}/${testCount} tests passed (${(100*passCount/testCount).toFixed(1)}%)`);
console.log("=========================================================\n");

if (passCount === testCount) {
    console.log("✓ ALL TESTS PASSED - NGC 1792 module validated!");
} else {
    console.log(`✗ ${testCount - passCount} test(s) failed - review required`);
}

// Additional diagnostic output
console.log("\n--- Diagnostic Information ---");
console.log(`M0: ${(ngc.M0 / M_sun).toExponential(4)} M_sun`);
console.log(`r: ${(ngc.r / ly_to_m).toFixed(0)} ly`);
console.log(`SFR_factor: ${ngc.SFR_factor.toExponential(4)}`);
console.log(`tau_SF: ${(ngc.tau_SF / Myr_to_s).toFixed(0)} Myr`);
console.log(`B: ${ngc.B.toExponential(4)} T`);
console.log(`rho_wind: ${ngc.rho_wind.toExponential(4)} kg/m³`);
console.log(`v_wind: ${ngc.v_wind.toExponential(4)} m/s`);
console.log(`Feedback term: ${feedback_term.toExponential(6)} m/s²`);

// Sample time evolution
console.log("\n--- Sample Time Evolution ---");
[0, 25, 50, 100, 200].forEach(t_Myr => {
    const t = t_Myr * Myr_to_s;
    const Mt = ngc.M_t(t);
    const g = ngc.compute_g_NGC1792(t);
    console.log(`t = ${t_Myr} Myr: M(t) = ${(Mt/M_sun).toExponential(4)} M_sun, g = ${g.toExponential(6)} m/s²`);
});
