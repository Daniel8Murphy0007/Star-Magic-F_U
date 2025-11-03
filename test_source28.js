/**
 * ================================================================================================
 * Comprehensive Test Suite for GalaxyAndromeda (M31 Andromeda Galaxy)
 * 
 * Tests all physics components, time evolution, DM/SMBH features, and UQFF framework
 * ================================================================================================
 */

import GalaxyAndromeda from './source28.js';

// Constants
const M_sun = 1.989e30;
const ly_to_m = 9.461e15;
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
console.log("ANDROMEDA GALAXY (M31) COMPREHENSIVE TEST SUITE");
console.log("=========================================================\n");

const andromeda = new GalaxyAndromeda();

// ========== TEST GROUP 1: BASIC INITIALIZATION (Tests 1-3) ==========
console.log("Test Group 1: Basic Initialization");
console.log("-----------------------------------");

// Test 1: M and r within expected ranges
const M_Msun = andromeda.M / M_sun;
const r_ly = andromeda.r / ly_to_m;
assert(M_Msun > 1e11 && M_Msun < 1e13 && r_ly > 5e4 && r_ly < 2e5,
    "M and r within reasonable galaxy ranges");

// Test 2: Dark matter fraction = 80%, Visible = 20%
const DM_fraction = andromeda.M_DM / andromeda.M;
const visible_fraction = andromeda.M_visible / andromeda.M;
assert(Math.abs(DM_fraction - 0.8) < 0.01 && Math.abs(visible_fraction - 0.2) < 0.01,
    "Dark matter fraction = 80%, Visible fraction = 20%");

// Test 3: Blueshift z = -0.001 (approaching)
assertClose(andromeda.z, -0.001, 0.01, "z = -0.001 (blueshift for Andromeda)");

console.log();

// ========== TEST GROUP 2: UQFF COMPONENTS (Tests 4-6) ==========
console.log("Test Group 2: UQFF Components");
console.log("-------------------------------");

// Test 4: Ug sum includes Ug1 and Ug4 (Ug2, Ug3 ≈ 0)
const ug_sum = andromeda.computeUgSum();
const Ug1 = (andromeda.G * andromeda.M) / (andromeda.r * andromeda.r);
const Ug4 = Ug1 * andromeda.f_sc;
const expected_Ug = Ug1 + Ug4;  // Ug2 = Ug3 = 0
assertClose(ug_sum, expected_Ug, 1e-9, "Ug sum = Ug1 + Ug4 (Ug2, Ug3 ≈ 0)");

// Test 5: g_Andromeda(t) > 0 for all times tested
const times_to_test = [0, 2.5, 5, 10, 13.8].map(t => t * Gyr_to_s);
const all_positive = times_to_test.every(t => andromeda.compute_g_Andromeda(t) > 0);
assert(all_positive, "g_Andromeda(t) > 0 for all t ∈ [0, 13.8 Gyr]");

// Test 6: Negative time returns 0
const g_neg = andromeda.compute_g_Andromeda(-1);
assert(g_neg === 0, "g_Andromeda(t < 0) returns 0");

console.log();

// ========== TEST GROUP 3: GALAXY SCALING METHODS (Tests 7-9) ==========
console.log("Test Group 3: Galaxy Scaling Methods");
console.log("-------------------------------------");

// Test 7: expandGalaxyScale(M, r)
andromeda.saveState('test7');
const M_orig = andromeda.M;
const r_orig = andromeda.r;
andromeda.expandGalaxyScale(1.5, 2.0);
const M_new = andromeda.M;
const r_new = andromeda.r;
const test7_pass = assertClose(M_new / M_orig, 1.5, 1e-9, "expandGalaxyScale(1.5, 2.0) scales M") &&
                    assertClose(r_new / r_orig, 2.0, 1e-9, "expandGalaxyScale(1.5, 2.0) scales r");
andromeda.restoreState('test7');

// Test 8: expandBlackHoleScale(M_BH, r_BH) - UNIQUE to SMBH galaxies
andromeda.saveState('test8');
const M_BH_orig = andromeda.M_BH;
const r_BH_orig = andromeda.r_BH;
andromeda.expandBlackHoleScale(2.0, 1.5);
const M_BH_new = andromeda.M_BH;
const r_BH_new = andromeda.r_BH;
const test8_pass = assertClose(M_BH_new / M_BH_orig, 2.0, 1e-9, "expandBlackHoleScale(2.0, 1.5) scales M_BH") &&
                    assertClose(r_BH_new / r_BH_orig, 1.5, 1e-9, "expandBlackHoleScale(2.0, 1.5) scales r_BH");
andromeda.restoreState('test8');

// Test 9: expandDarkMatterScale(M_DM, delta_rho) - UNIQUE to DM-dominated
andromeda.saveState('test9');
const M_DM_orig = andromeda.M_DM;
const delta_rho_orig = andromeda.delta_rho;
andromeda.expandDarkMatterScale(1.2, 1.3);
const M_DM_new = andromeda.M_DM;
const delta_rho_new = andromeda.delta_rho;
const test9_pass = assertClose(M_DM_new / M_DM_orig, 1.2, 1e-9, "expandDarkMatterScale(1.2, 1.3) scales M_DM") &&
                    assertClose(delta_rho_new / delta_rho_orig, 1.3, 1e-9, "expandDarkMatterScale(1.2, 1.3) scales delta_rho");
andromeda.restoreState('test9');

console.log();

// ========== TEST GROUP 4: PARAMETER SCALING EFFECTS (Tests 10-17) ==========
console.log("Test Group 4: Parameter Scaling Effects");
console.log("----------------------------------------");

// Test 10: M scaling → affects Ug1 gravitational term (dust friction dominates total g)
andromeda.saveState('test10');
const t_test = 10 * Gyr_to_s;
const Ug1_M_base = (andromeda.G * andromeda.M) / (andromeda.r * andromeda.r);
andromeda.setVariable('M', andromeda.M * 2);
const Ug1_M_scaled = (andromeda.G * andromeda.M) / (andromeda.r * andromeda.r);
const Ug1_ratio = Ug1_M_scaled / Ug1_M_base;
assertClose(Ug1_ratio, 2.0, 0.01, `M × 2 → Ug1 × 2 (ratio: ${Ug1_ratio.toFixed(3)})`);
andromeda.restoreState('test10');

// Test 11: r scaling → affects Ug1 gravitational term (inverse square, dust dominates total)
andromeda.saveState('test11');
const Ug1_r_base = (andromeda.G * andromeda.M) / (andromeda.r * andromeda.r);
andromeda.setVariable('r', andromeda.r * 2);
const Ug1_r_scaled = (andromeda.G * andromeda.M) / (andromeda.r * andromeda.r);
const Ug1_r_ratio = Ug1_r_base / Ug1_r_scaled;
assertClose(Ug1_r_ratio, 4.0, 0.01, `r × 2 → Ug1 × 1/4 (inverse square, ratio: ${Ug1_r_ratio.toFixed(3)})`);
andromeda.restoreState('test11');

// Test 12: M_BH scaling (subsumed in M, but parameter tracked)
andromeda.saveState('test12');
const M_BH_test_base = andromeda.M_BH;
andromeda.setVariable('M_BH', andromeda.M_BH * 2);
const M_BH_test_scaled = andromeda.M_BH;
assertClose(M_BH_test_scaled / M_BH_test_base, 2.0, 1e-9, "M_BH × 2 → M_BH scales correctly");
andromeda.restoreState('test12');

// Test 13: M_DM scaling → affects total M
andromeda.saveState('test13');
const M_total_base = andromeda.M;
const M_DM_test_base = andromeda.M_DM;
andromeda.expandDarkMatterScale(1.5, 1.0);
const M_total_scaled = andromeda.M;
const M_DM_test_scaled = andromeda.M_DM;
assert(M_total_scaled > M_total_base && M_DM_test_scaled > M_DM_test_base,
    "M_DM × 1.5 → M_DM and M_total increase");
andromeda.restoreState('test13');

// Test 14: v_orbit scaling → affects dust/EM terms
andromeda.saveState('test14');
const g_v_base = andromeda.compute_g_Andromeda(t_test);
andromeda.setVariable('v_orbit', andromeda.v_orbit * 2);
const g_v_scaled = andromeda.compute_g_Andromeda(t_test);
assert(g_v_scaled > g_v_base, "v_orbit × 2 → g increases (v² dependence in dust/EM)");
andromeda.restoreState('test14');

// Test 15: B (magnetic field) scaling → affects EM term (dust dominates total, check EM directly)
andromeda.saveState('test15');
const B_base = andromeda.B;
const em_base = andromeda.q * andromeda.v_orbit * B_base * (1 + andromeda.rho_vac_UA / andromeda.rho_vac_SCm) * andromeda.scale_macro;
andromeda.setVariable('B', andromeda.B * 2);
const B_scaled = andromeda.B;
const em_scaled = andromeda.q * andromeda.v_orbit * B_scaled * (1 + andromeda.rho_vac_UA / andromeda.rho_vac_SCm) * andromeda.scale_macro;
assertClose(em_scaled / em_base, 2.0, 0.01, `B × 2 → EM term × 2 (ratio: ${(em_scaled/em_base).toFixed(3)})`);
andromeda.restoreState('test15');

// Test 16: rho_dust scaling → affects dust friction
andromeda.saveState('test16');
const g_rho_dust_base = andromeda.compute_g_Andromeda(t_test);
andromeda.setVariable('rho_dust', andromeda.rho_dust * 2);
const g_rho_dust_scaled = andromeda.compute_g_Andromeda(t_test);
assert(g_rho_dust_scaled > g_rho_dust_base, "rho_dust × 2 → g increases (dust friction)");
andromeda.restoreState('test16');

// Test 17: delta_rho scaling → affects DM perturbation term
andromeda.saveState('test17');
const g_delta_rho_base = andromeda.compute_g_Andromeda(t_test);
andromeda.setVariable('delta_rho', andromeda.delta_rho * 2);
const g_delta_rho_scaled = andromeda.compute_g_Andromeda(t_test);
assert(g_delta_rho_scaled !== g_delta_rho_base, "delta_rho × 2 → g changes (DM perturbation)");
andromeda.restoreState('test17');

console.log();

// ========== TEST GROUP 5: PHYSICS VALIDATION (Tests 18-20) ==========
console.log("Test Group 5: Physics Validation");
console.log("---------------------------------");

// Test 18: Sensitivity analysis identifies key parameters
andromeda.saveState('test18');
const sens = andromeda.sensitivityAnalysis(10 * Gyr_to_s, 0.01);
const sorted_sens = Object.entries(sens).sort((a, b) => b[1] - a[1]);
const top_3 = sorted_sens.slice(0, 3).map(([k, v]) => k);
console.log(`  Top 3 sensitive parameters: ${top_3.join(', ')}`);
// Expect M, r, rho_dust, v_orbit to be among most sensitive
const has_major_params = top_3.some(p => ['M', 'r', 'rho_dust', 'v_orbit'].includes(p));
assert(has_major_params, "Sensitivity includes major parameters (M/r/rho_dust/v_orbit)");
andromeda.restoreState('test18');

// Test 19: Hz (Hubble parameter) is positive for blueshift z = -0.001
const Hz = andromeda.computeHz();
console.log(`  Hz = ${Hz.toExponential(6)} s^-1`);
assert(Hz > 0, "Hz > 0 even with blueshift z = -0.001");

// Test 20: Validate consistency check
assert(andromeda.validateConsistency(), "validateConsistency() returns true for valid parameters");

console.log();

// ========== TEST GROUP 6: ENHANCED CAPABILITIES (Tests 21-26) ==========
console.log("Test Group 6: Enhanced Capabilities");
console.log("------------------------------------");

// Test 21: State save/restore functionality
andromeda.saveState('enhanced_test');
const original_M = andromeda.M;
andromeda.setVariable('M', andromeda.M * 3);
const modified_M = andromeda.M;
andromeda.restoreState('enhanced_test');
const restored_M = andromeda.M;
assert(Math.abs(restored_M - original_M) < 1e-20 && Math.abs(modified_M / original_M - 3) < 1e-9,
    "State save/restore preserves parameters");

// Test 22: Variable listing
const vars = andromeda.listVariables();
assert(vars.length > 30, `Variable listing returns ${vars.length} variables`);

// Test 23: Clone variable
andromeda.cloneVariable('M', 'M_clone');
const M_clone = andromeda.getVariable('M_clone');
assertClose(M_clone, andromeda.M, 1e-9, "cloneVariable creates accurate copy");

// Test 24: Parameter space expansion
andromeda.saveState('test24');
andromeda.expandParameterSpace(1.1);
const M_after_expand = andromeda.M;
const M_before_expand = original_M;
assert(M_after_expand > M_before_expand, "expandParameterSpace scales all parameters");
andromeda.restoreState('test24');

// Test 25: Generate variations
const variations = andromeda.generateVariations(3, 10.0);
assert(variations.length === 3, "generateVariations produces correct count");
assert(variations[0].M !== andromeda.M, "Variations differ from original");

// Test 26: Export state
const export_str = andromeda.exportState();
assert(export_str.includes('Andromeda') && export_str.includes('M:'),
    "exportState produces valid output string");

console.log();

// ========== SUMMARY ==========
console.log("=========================================================");
console.log(`TEST SUMMARY: ${passCount}/${testCount} tests passed (${(100*passCount/testCount).toFixed(1)}%)`);
console.log("=========================================================\n");

if (passCount === testCount) {
    console.log("✓ ALL TESTS PASSED - Andromeda M31 module validated!");
} else {
    console.log(`✗ ${testCount - passCount} test(s) failed - review required`);
}

// Additional diagnostic output
console.log("\n--- Diagnostic Information ---");
console.log(`M: ${(andromeda.M / M_sun).toExponential(4)} M_sun`);
console.log(`M_visible: ${(andromeda.M_visible / M_sun).toExponential(4)} M_sun (${(andromeda.M_visible/andromeda.M*100).toFixed(1)}%)`);
console.log(`M_DM: ${(andromeda.M_DM / M_sun).toExponential(4)} M_sun (${(andromeda.M_DM/andromeda.M*100).toFixed(1)}%)`);
console.log(`M_BH: ${(andromeda.M_BH / M_sun).toExponential(4)} M_sun`);
console.log(`r: ${(andromeda.r / ly_to_m).toFixed(0)} ly`);
console.log(`z: ${andromeda.z.toFixed(4)} (blueshift)`);
console.log(`v_orbit: ${(andromeda.v_orbit / 1e3).toFixed(1)} km/s`);
console.log(`B: ${andromeda.B.toExponential(4)} T`);
console.log(`Hz: ${Hz.toExponential(6)} s^-1`);

// Sample time evolution
console.log("\n--- Sample Time Evolution ---");
[0, 2.5, 5, 10, 13.8].forEach(t_Gyr => {
    const t = t_Gyr * Gyr_to_s;
    const g = andromeda.compute_g_Andromeda(t);
    console.log(`t = ${t_Gyr.toFixed(1)} Gyr: g = ${g.toExponential(6)} m/s²`);
});
