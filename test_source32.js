/**
 * Comprehensive Test Suite for NebulaCrab (Crab Nebula M1)
 * Tests all 38 methods (15 core + 23 enhanced) and unique pulsar/magnetic features
 */

import NebulaCrab from './source32.js';

const TESTS_PASSED = [];
const TESTS_FAILED = [];

function assert(condition, testName, message = '') {
    if (condition) {
        TESTS_PASSED.push(testName);
        console.log(`✓ Test ${TESTS_PASSED.length + TESTS_FAILED.length}: ${testName} PASSED`);
    } else {
        TESTS_FAILED.push({ name: testName, message });
        console.log(`✗ Test ${TESTS_PASSED.length + TESTS_FAILED.length}: ${testName} FAILED - ${message}`);
    }
}

console.log("========================================");
console.log("COMPREHENSIVE TEST SUITE: Crab Nebula M1");
console.log("========================================\n");

const crab = new NebulaCrab();
const yr_to_s = crab.year_to_s;
const M_sun = crab.M_sun;

// ========== GROUP 1: INITIALIZATION & BASIC PARAMETERS (8 tests) ==========
console.log("GROUP 1: Initialization & Basic Parameters\n");

// Test 1: Mass within range
assert(crab.M > 4 * M_sun && crab.M < 5 * M_sun, 
    "Crab mass within range",
    `M = ${crab.M/M_sun} M☉`);

// Test 2: Initial radius within range  
assert(crab.r0 > 5e16 && crab.r0 < 6e16,
    "Crab initial radius within range",
    `r0 = ${crab.r0/9.461e15} ly`);

// Test 3: Expansion velocity
assert(Math.abs(crab.v_exp / 1e6 - 1.5) < 0.01,
    "v_exp = 1.5×10⁶ m/s",
    `v_exp = ${crab.v_exp/1e3} km/s`);

// Test 4: Pulsar power
assert(Math.abs(crab.P_pulsar / 5e31 - 1.0) < 0.01,
    "P_pulsar = 5×10³¹ W",
    `P_pulsar = ${crab.P_pulsar.toExponential(3)} W`);

// Test 5: Redshift for nearby remnant
assert(Math.abs(crab.z - 0.0015) < 0.0001,
    "z = 0.0015 (nearby)",
    `z = ${crab.z}`);

// Test 6: Shock velocity
assert(Math.abs(crab.v_shock / 1e6 - 1.5) < 0.01,
    "v_shock = 1.5×10⁶ m/s",
    `v_shock = ${crab.v_shock/1e3} km/s`);

// Test 7: Nebula density
assert(crab.rho_fluid === 1e-21,
    "rho_fluid = 1e-21 kg/m³",
    `rho_fluid = ${crab.rho_fluid}`);

// Test 8: No dark matter in remnant
assert(crab.M_DM === 0 && crab.M_visible === crab.M,
    "No DM, M_visible = M",
    `M_DM = ${crab.M_DM}, M_visible = ${crab.M_visible/M_sun} M☉`);

// ========== GROUP 2: TIME-DEPENDENT RADIUS EVOLUTION (6 tests) ==========
console.log("\nGROUP 2: Time-Dependent Radius Evolution (r(t) = r0 + v_exp × t)\n");

// Test 9: r(t) grows linearly with time
const t1 = 500 * yr_to_s;
const t2 = 1000 * yr_to_s;
const r1 = crab.computeRadius(t1);
const r2 = crab.computeRadius(t2);
const dr_expected = crab.v_exp * (t2 - t1);
const dr_actual = r2 - r1;
assert(Math.abs(dr_actual - dr_expected) < 1e10,
    "r(t) grows linearly with time",
    `Δr expected: ${(dr_expected/9.461e15).toFixed(2)} ly, actual: ${(dr_actual/9.461e15).toFixed(2)} ly`);

// Test 10: r(t) at 971 years matches r0 + v_exp × t
const t_971 = 971 * yr_to_s;
const r_971 = crab.computeRadius(t_971);
const r_expected = crab.r0 + crab.v_exp * t_971;
assert(Math.abs(r_971 - r_expected) < 1e6,
    "r(971 yr) = r0 + v_exp × 971 yr",
    `r = ${(r_971/9.461e15).toFixed(2)} ly`);

// Test 11: r(0) = r0
const r_0 = crab.computeRadius(0);
assert(Math.abs(r_0 - crab.r0) < 1.0,
    "r(0) = r0",
    `r(0) = ${(r_0/9.461e15).toFixed(2)} ly`);

// Test 12: Radius increases with time
const t_early = 100 * yr_to_s;
const t_late = 2000 * yr_to_s;
const r_early = crab.computeRadius(t_early);
const r_late = crab.computeRadius(t_late);
assert(r_late > r_early,
    "Radius increases with time",
    `r(100 yr) = ${(r_early/9.461e15).toFixed(2)} ly, r(2000 yr) = ${(r_late/9.461e15).toFixed(2)} ly`);

// Test 13: v_exp affects radius growth rate
crab.saveState('test13');
const r_base = crab.computeRadius(t_971);
crab.setVariable('v_exp', 2 * crab.v_exp);
const r_double_v = crab.computeRadius(t_971);
crab.restoreState('test13');
const ratio = (r_double_v - crab.r0) / (r_base - crab.r0);
assert(Math.abs(ratio - 2.0) < 0.01,
    "Doubling v_exp doubles radius growth",
    `ratio = ${ratio.toFixed(3)}`);

// Test 14: r0 affects initial radius
crab.saveState('test14');
const r0_original = crab.r0;
crab.setVariable('r0', 2 * r0_original);
const r_double_r0 = crab.computeRadius(t_971);
crab.restoreState('test14');
const r_orig = crab.computeRadius(t_971);
assert(r_double_r0 > r_orig,
    "Larger r0 increases r(t)",
    `r increased from ${(r_orig/9.461e15).toFixed(2)} to ${(r_double_r0/9.461e15).toFixed(2)} ly`);

// ========== GROUP 3: CORE PHYSICS COMPUTATIONS (8 tests) ==========
console.log("\nGROUP 3: Core Physics Computations\n");

// Test 15: Hz near zero for small z
const Hz = crab.computeHz();
assert(Hz > 0 && Hz < 1e-16,
    "Hz ≈ 0 for z = 0.0015",
    `Hz = ${Hz.toExponential(3)} s⁻¹`);

// Test 16: Ug sum positive
const ug_sum = crab.computeUgSum(r_971);
assert(ug_sum > 0,
    "Ug sum > 0",
    `Ug = ${ug_sum.toExponential(3)} m/s²`);

// Test 17: Quantum term small
const quantum = crab.computeQuantumTerm(crab.t_Hubble);
assert(quantum > 0 && quantum < 1e-5,
    "Quantum term small",
    `quantum = ${quantum.toExponential(3)} m/s²`);

// Test 18: Fluid term depends on rho_fluid
const g_base_test = 1e-3;
const fluid1 = crab.computeFluidTerm(g_base_test);
crab.setVariable('rho_fluid', 2e-21);
const fluid2 = crab.computeFluidTerm(g_base_test);
crab.setVariable('rho_fluid', 1e-21);  // Reset
const fluid_ratio = fluid2 / fluid1;
assert(Math.abs(fluid_ratio - 2.0) < 0.01,
    "Fluid term scales with rho_fluid",
    `fluid ratio = ${fluid_ratio.toFixed(3)}`);

// Test 19: Resonant term oscillates
const resonant1 = crab.computeResonantTerm(0);
const resonant2 = crab.computeResonantTerm(Math.PI / crab.omega);
assert(resonant1 !== resonant2,
    "Resonant term oscillates with time",
    `resonant(0) = ${resonant1.toExponential(3)}, resonant(π/ω) = ${resonant2.toExponential(3)}`);

// Test 20: DM term reflects density perturbation
const dm_term = crab.computeDMTerm(r_971) / crab.M;
assert(Math.abs(dm_term - 0.1) < 0.01,
    "DM term ≈ δρ/ρ = 0.1 (density perturbation)",
    `DM term = ${dm_term.toExponential(3)} m/s²`);

// Test 21: Wind term inversely proportional to r²
const r_small = 1e16;
const r_large = 4e16;
const wind_small = crab.computeWindTerm(r_small);
const wind_large = crab.computeWindTerm(r_large);
const wind_ratio = wind_small / wind_large;
const expected_ratio = (r_large / r_small) ** 2;
assert(Math.abs(wind_ratio - expected_ratio) < 0.1,
    "Wind term ∝ 1/r²",
    `wind ratio = ${wind_ratio.toFixed(2)}, expected = ${expected_ratio.toFixed(2)}`);

// Test 22: Magnetic term depends on v_shock and B
const mag_base = crab.computeMagTerm();
crab.saveState('test22');
crab.setVariable('v_shock', 2 * crab.v_shock);
crab.setVariable('B', 2 * crab.B);
const mag_scaled = crab.computeMagTerm();
crab.restoreState('test22');
const mag_ratio = mag_scaled / mag_base;
assert(Math.abs(mag_ratio - 4.0) < 0.1,
    "Magnetic term ∝ v_shock × B",
    `mag ratio = ${mag_ratio.toFixed(2)} (expect 4.0)`);

// ========== GROUP 4: UNIQUE EXPANSION METHODS (6 tests) ==========
console.log("\nGROUP 4: Unique Expansion Methods (Pulsar & Magnetic)\n");

// Test 23: expandPulsarScale - P_pulsar scaling
crab.saveState('test23');
const P_orig = crab.P_pulsar;
crab.expandPulsarScale(3.0, 1.0);
assert(Math.abs(crab.P_pulsar / P_orig - 3.0) < 0.01,
    "expandPulsarScale scales P_pulsar",
    `P_pulsar scaled ${(crab.P_pulsar/P_orig).toFixed(2)}×`);
crab.restoreState('test23');

// Test 24: expandPulsarScale - v_exp scaling
crab.saveState('test24');
const v_exp_orig = crab.v_exp;
crab.expandPulsarScale(1.0, 2.5);
assert(Math.abs(crab.v_exp / v_exp_orig - 2.5) < 0.01,
    "expandPulsarScale scales v_exp",
    `v_exp scaled ${(crab.v_exp/v_exp_orig).toFixed(2)}×`);
crab.restoreState('test24');

// Test 25: expandMagneticScale - B scaling
crab.saveState('test25');
const B_orig = crab.B;
crab.expandMagneticScale(4.0, 1.0);
assert(Math.abs(crab.B / B_orig - 4.0) < 0.01,
    "expandMagneticScale scales B",
    `B scaled ${(crab.B/B_orig).toFixed(2)}×`);
crab.restoreState('test25');

// Test 26: expandMagneticScale - v_shock scaling
crab.saveState('test26');
const v_shock_orig = crab.v_shock;
crab.expandMagneticScale(1.0, 2.0);
assert(Math.abs(crab.v_shock / v_shock_orig - 2.0) < 0.01,
    "expandMagneticScale scales v_shock",
    `v_shock scaled ${(crab.v_shock/v_shock_orig).toFixed(2)}×`);
crab.restoreState('test26');

// Test 27: expandNebulaScale affects M and r0
crab.saveState('test27');
const M_orig = crab.M;
const r0_init = crab.r0;
crab.expandNebulaScale(1.5, 2.0);
const M_ratio = crab.M / M_orig;
const r0_ratio = crab.r0 / r0_init;
assert(Math.abs(M_ratio - 1.5) < 0.01 && Math.abs(r0_ratio - 2.0) < 0.01,
    "expandNebulaScale scales M and r0 independently",
    `M scaled ${M_ratio.toFixed(2)}×, r0 scaled ${r0_ratio.toFixed(2)}×`);
crab.restoreState('test27');

// Test 28: P_pulsar affects wind term at fixed time
crab.saveState('test28');
const wind_base = crab.computeWindTerm(r_971);
crab.expandPulsarScale(2.0, 1.0);
const wind_scaled = crab.computeWindTerm(r_971);
crab.restoreState('test28');
const wind_ratio_test = wind_scaled / wind_base;
assert(Math.abs(wind_ratio_test - 2.0) < 0.01,
    "P_pulsar × 2 → wind × 2",
    `wind ratio = ${wind_ratio_test.toFixed(3)}`);

// ========== GROUP 5: PARAMETER SENSITIVITY (4 tests) ==========
console.log("\nGROUP 5: Parameter Sensitivity\n");

// Test 29: v_exp affects r(t) and g
crab.saveState('test29');
const g_base_vexp = crab.compute_g_Crab(t_971);
crab.setVariable('v_exp', 2 * crab.v_exp);
const g_high_vexp = crab.compute_g_Crab(t_971);
crab.restoreState('test29');
// Higher v_exp → larger r(t) → smaller g_base, but wind term also changes
assert(g_high_vexp !== g_base_vexp,
    "v_exp affects g (via r(t))",
    `g changed from ${g_base_vexp.toExponential(3)} to ${g_high_vexp.toExponential(3)}`);

// Test 30: P_pulsar dominates total g
crab.saveState('test30');
const g_total = crab.compute_g_Crab(t_971);
const wind_contrib = crab.computeWindTerm(r_971);
const wind_fraction = wind_contrib / g_total;
crab.restoreState('test30');
assert(wind_fraction > 0.9,
    "Wind term dominates (>90% of total g)",
    `wind fraction = ${(wind_fraction*100).toFixed(1)}%`);

// Test 31: B affects magnetic term
crab.saveState('test31');
const mag_base_B = crab.computeMagTerm();
crab.setVariable('B', 2 * crab.B);
const mag_high_B = crab.computeMagTerm();
crab.restoreState('test31');
const mag_B_ratio = mag_high_B / mag_base_B;
assert(Math.abs(mag_B_ratio - 2.0) < 0.01,
    "B × 2 → M_mag × 2",
    `mag ratio = ${mag_B_ratio.toFixed(3)}`);

// Test 32: Sensitivity analysis returns valid data
const sensitivity = crab.sensitivityAnalysis(t_971, 0.01);
const has_P = 'P_pulsar' in sensitivity;
const has_v_expansion = 'v_exp' in sensitivity;
assert(has_P && has_v_expansion,
    "sensitivityAnalysis includes P_pulsar and v_exp",
    `Keys: ${Object.keys(sensitivity).length}`);

// ========== GROUP 6: ENHANCED CAPABILITIES (6 tests) ==========
console.log("\nGROUP 6: Enhanced Capabilities\n");

// Test 33: State save/restore
crab.saveState('test33');
const M_before = crab.M;
crab.setVariable('M', 10 * M_sun);
crab.restoreState('test33');
assert(crab.M === M_before,
    "State restore works",
    `M restored to ${crab.M/M_sun} M☉`);

// Test 34: Variable list contains key params
const vars = crab.listVariables();
const has_M = vars.includes('M');
const has_P_pulsar = vars.includes('P_pulsar');
const has_v_expansion_var = vars.includes('v_exp');
assert(has_M && has_P_pulsar && has_v_expansion_var,
    "listVariables includes M, P_pulsar, v_exp",
    `Total variables: ${vars.length}`);

// Test 35: Generate variations
const variations = crab.generateVariations(5, 10.0);
assert(variations.length === 5,
    "generateVariations creates 5 variants",
    `Created ${variations.length} variations`);

// Test 36: Validate consistency
const is_valid = crab.validateConsistency();
assert(is_valid,
    "validateConsistency passes",
    `Valid: ${is_valid}`);

// Test 37: System name
const name = crab.getSystemName();
assert(name === "CrabNebula_M1",
    "getSystemName returns CrabNebula_M1",
    `Name: ${name}`);

// Test 38: Generate report contains key info
const report = crab.generateReport(t_971);
const has_pulsar_in_report = report.includes("Pulsar");
const has_expansion_in_report = report.includes("Expansion");
assert(has_pulsar_in_report && has_expansion_in_report,
    "generateReport includes pulsar and expansion data",
    `Report length: ${report.length} chars`);

// ========== SUMMARY ==========
console.log("\n========================================");
console.log("TEST SUMMARY");
console.log("========================================");
console.log(`Total Tests: ${TESTS_PASSED.length + TESTS_FAILED.length}`);
console.log(`Passed: ${TESTS_PASSED.length}`);
console.log(`Failed: ${TESTS_FAILED.length}`);

if (TESTS_FAILED.length > 0) {
    console.log("\nFailed Tests:");
    TESTS_FAILED.forEach(test => {
        console.log(`  - ${test.name}: ${test.message}`);
    });
}

console.log("\n========================================\n");

// Exit with appropriate code
process.exit(TESTS_FAILED.length > 0 ? 1 : 0);
