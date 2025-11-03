// test_source33.js
// Comprehensive test suite for SGR 1745-2900 Magnetar module (source33.js)
// 38 tests covering initialization, physics, unique magnetar features, and enhanced capabilities

import { MagnetarSGR1745 } from './source33.js';

let passed = 0;
let failed = 0;

function assert(condition, testName) {
    if (condition) {
        console.log(`✓ ${testName}`);
        passed++;
    } else {
        console.log(`✗ ${testName}`);
        failed++;
    }
}

console.log("========== SGR 1745-2900 MAGNETAR COMPREHENSIVE TEST SUITE ==========\n");

const sgr = new MagnetarSGR1745();
const t_test = 1000 * 3.156e7; // 1000 years

// ========== GROUP 1: INITIALIZATION & BASIC PARAMETERS (8 tests) ==========
console.log("GROUP 1: Initialization & Basic Parameters");

// Test 1: Mass initialization
assert(
    Math.abs(sgr.M / sgr.M_sun - 1.4) < 0.01,
    "Test 1: Mass M = 1.4 M_sun"
);

// Test 2: Radius initialization
assert(
    Math.abs(sgr.r - 1e4) < 1 && Math.abs(sgr.r / 1e3 - 10) < 0.1,
    "Test 2: Radius r = 10 km"
);

// Test 3: Magnetic field
assert(
    Math.abs(sgr.B - 2e10) < 1e8 && Math.abs(sgr.B / 1e10 - 2.0) < 0.01,
    "Test 3: Magnetic field B = 2×10^10 T"
);

// Test 4: Spin period
assert(
    Math.abs(sgr.P - 3.76) < 0.01,
    "Test 4: Spin period P = 3.76 s"
);

// Test 5: Redshift (Galactic Center)
assert(
    Math.abs(sgr.z - 0.0) < 0.001,
    "Test 5: Redshift z ≈ 0 (Galactic Center)"
);

// Test 6: Crust density
assert(
    Math.abs(sgr.rho_fluid - 1e17) < 1e15,
    "Test 6: Crust density rho_fluid = 1×10^17 kg/m³"
);

// Test 7: Critical field
assert(
    Math.abs(sgr.B_crit - 1e11) < 1e9,
    "Test 7: Critical field B_crit = 1×10^11 T"
);

// Test 8: No dark matter
assert(
    sgr.M_DM === 0 && sgr.M_visible === sgr.M,
    "Test 8: M_DM = 0, M_visible = M"
);

console.log();

// ========== GROUP 2: CORE PHYSICS COMPUTATIONS (8 tests) ==========
console.log("GROUP 2: Core Physics Computations");

// Test 9: Hubble parameter near zero
const Hz = sgr.computeHz();
assert(
    Hz > 0 && Hz < 1e-17,
    "Test 9: Hubble parameter Hz > 0 (small for z ≈ 0)"
);

// Test 10: Ug sum positive
const ug_sum = sgr.computeUgSum();
assert(
    ug_sum > 1e12 && ug_sum < 1e13,
    "Test 10: Ug sum in reasonable range (e12 m/s²)"
);

// Test 11: Quantum term small
const quantum = sgr.computeQuantumTerm();
assert(
    Math.abs(quantum) < 1e-30,
    "Test 11: Quantum term very small (< 1e-30 m/s²)"
);

// Test 12: EM term significant
const em_term = sgr.computeEMTerm();
assert(
    em_term > 1e11 && em_term < 1e12,
    "Test 12: EM term significant (1e11-1e12 m/s²)"
);

// Test 13: Fluid term micro-level
const g_base_test = 1e12;
const fluid_term = sgr.computeFluidTerm(g_base_test);
assert(
    Math.abs(fluid_term) < 1e-3,
    "Test 13: Fluid term micro-level (< 1e-3 m/s²)"
);

// Test 14: Resonant term oscillatory
const resonant_t0 = sgr.computeResonantTerm(0);
const resonant_t1 = sgr.computeResonantTerm(sgr.P);
assert(
    Math.abs(resonant_t0) < 1e-9 && Math.abs(resonant_t1) < 1e-9,
    "Test 14: Resonant term oscillates (small amplitude)"
);

// Test 15: DM term small (no DM)
const dm_term = sgr.computeDMTerm();
assert(
    dm_term > 0 && dm_term < 1e6,
    "Test 15: DM term small (< 1e6 m/s², perturbation only)"
);

// Test 16: Full g computation
const g_full = sgr.compute_g_SGR1745(t_test);
assert(
    g_full > 1e12 && g_full < 1e13,
    "Test 16: Full g_SGR1745 in NS range (1e12-1e13 m/s²)"
);

console.log();

// ========== GROUP 3: MAGNETAR-SPECIFIC FEATURES (6 tests) ==========
console.log("GROUP 3: Magnetar-Specific Features");

// Test 17: SC correction significant
const sc_corr = sgr.computeSCCorrection();
assert(
    Math.abs(sc_corr - 0.8) < 0.01,
    "Test 17: SC correction = 0.8 (significant for B = 2e10 T)"
);

// Test 18: SC correction responds to B
sgr.B = 5e10; // Half of B_crit
const sc_half = sgr.computeSCCorrection();
sgr.B = 2e10; // Restore
assert(
    Math.abs(sc_half - 0.5) < 0.01,
    "Test 18: SC correction = 0.5 when B = B_crit / 2"
);

// Test 19: Spin velocity at equator
const v_expected = (2 * Math.PI * sgr.r) / sgr.P;
assert(
    Math.abs(sgr.v_spin - v_expected) / v_expected < 0.01,
    "Test 19: v_spin = 2πr / P (equatorial velocity)"
);

// Test 20: EM term scales with B
const em_B1 = sgr.computeEMTerm();
sgr.B = 4e10; // Double B
const em_B2 = sgr.computeEMTerm();
sgr.B = 2e10; // Restore
const em_ratio = em_B2 / em_B1;
assert(
    Math.abs(em_ratio - 2.0) < 0.1,
    "Test 20: EM term doubles when B doubles"
);

// Test 21: EM term scales with v_spin
sgr.P = 1.88; // Half period (double omega)
sgr.omega = 2 * sgr.pi / sgr.P;
sgr.v_spin = (2 * sgr.pi * sgr.r) / sgr.P;
const em_v2 = sgr.computeEMTerm();
sgr.P = 3.76; // Restore
sgr.omega = 2 * sgr.pi / sgr.P;
sgr.v_spin = (2 * sgr.pi * sgr.r) / sgr.P;
const em_v1 = sgr.computeEMTerm();
const v_ratio = em_v2 / em_v1;
assert(
    Math.abs(v_ratio - 2.0) < 0.1,
    "Test 21: EM term doubles when v_spin doubles"
);

// Test 22: Fluid term scales with crust density
const rho_orig = sgr.rho_fluid;
sgr.rho_fluid = 3e17;
const fluid_rho2 = sgr.computeFluidTerm(g_base_test);
sgr.rho_fluid = rho_orig;
const fluid_rho1 = sgr.computeFluidTerm(g_base_test);
const rho_ratio = fluid_rho2 / fluid_rho1;
assert(
    Math.abs(rho_ratio - 3.0) < 0.1,
    "Test 22: Fluid term triples when rho triples"
);

console.log();

// ========== GROUP 4: UNIQUE EXPANSION METHODS (6 tests) ==========
console.log("GROUP 4: Unique Expansion Methods");

// Test 23: expandMagnetarScale mass
const M_initial = sgr.M;
sgr.expandMagnetarScale(1.5, 1.0);
const M_scaled = sgr.M;
sgr.M = M_initial;
sgr.M_visible = M_initial;
assert(
    Math.abs(M_scaled / M_initial - 1.5) < 0.01,
    "Test 23: expandMagnetarScale scales mass by factor"
);

// Test 24: expandMagnetarScale radius
const r_initial = sgr.r;
sgr.expandMagnetarScale(1.0, 2.0);
const r_scaled = sgr.r;
sgr.r = r_initial;
sgr.v_spin = (2 * sgr.pi * sgr.r) / sgr.P;
assert(
    Math.abs(r_scaled / r_initial - 2.0) < 0.01,
    "Test 24: expandMagnetarScale scales radius by factor"
);

// Test 25: expandMagneticFieldScale B-field
const B_initial = sgr.B;
sgr.expandMagneticFieldScale(3.0, 1.0);
const B_scaled = sgr.B;
sgr.B = B_initial;
assert(
    Math.abs(B_scaled / B_initial - 3.0) < 0.01,
    "Test 25: expandMagneticFieldScale scales B-field"
);

// Test 26: expandMagneticFieldScale period
const P_initial = sgr.P;
sgr.expandMagneticFieldScale(1.0, 1.5);
const P_scaled = sgr.P;
sgr.P = P_initial;
sgr.omega = 2 * sgr.pi / sgr.P;
sgr.v_spin = (2 * sgr.pi * sgr.r) / sgr.P;
assert(
    Math.abs(P_scaled / P_initial - 1.5) < 0.01,
    "Test 26: expandMagneticFieldScale scales period"
);

// Test 27: expandCrustScale density
const rho_init = sgr.rho_fluid;
sgr.expandCrustScale(2.5, 1.0);
const rho_scaled = sgr.rho_fluid;
sgr.rho_fluid = rho_init;
sgr.rho = rho_init;
sgr.delta_rho = 0.1 * rho_init;
assert(
    Math.abs(rho_scaled / rho_init - 2.5) < 0.01,
    "Test 27: expandCrustScale scales density"
);

// Test 28: expandCrustScale volume
const V_init = sgr.V;
sgr.expandCrustScale(1.0, 2.0);
const V_scaled = sgr.V;
sgr.V = V_init;
assert(
    Math.abs(V_scaled / V_init - 2.0) < 0.01,
    "Test 28: expandCrustScale scales volume"
);

console.log();

// ========== GROUP 5: PARAMETER SENSITIVITY (4 tests) ==========
console.log("GROUP 5: Parameter Sensitivity");

// Test 29: g increases with mass
const g_M1 = sgr.compute_g_SGR1745(t_test);
sgr.M *= 1.2;
sgr.M_visible = sgr.M;
const g_M2 = sgr.compute_g_SGR1745(t_test);
sgr.M /= 1.2;
sgr.M_visible = sgr.M;
assert(
    g_M2 > g_M1,
    "Test 29: g increases with mass"
);

// Test 30: g decreases with radius
const g_r1 = sgr.compute_g_SGR1745(t_test);
sgr.r *= 1.2;
sgr.v_spin = (2 * sgr.pi * sgr.r) / sgr.P;
const g_r2 = sgr.compute_g_SGR1745(t_test);
sgr.r /= 1.2;
sgr.v_spin = (2 * sgr.pi * sgr.r) / sgr.P;
assert(
    g_r2 < g_r1,
    "Test 30: g decreases with radius"
);

// Test 31: EM term (not total g) increases with B-field
// Note: Total g may not increase due to SC correction (1 - B/B_crit) decreasing
const em_B1_test = sgr.computeEMTerm();
sgr.B *= 1.5;
const em_B2_test = sgr.computeEMTerm();
sgr.B /= 1.5;
assert(
    em_B2_test > em_B1_test,
    "Test 31: EM term increases with B-field"
);

// Test 32: Sensitivity analysis includes key params
const sens = sgr.sensitivityAnalysis(t_test, 0.01);
const has_M = 'M' in sens;
const has_B = 'B' in sens;
const has_P = 'P' in sens;
const has_rho = 'rho_fluid' in sens;
assert(
    has_M && has_B && has_P && has_rho,
    "Test 32: sensitivityAnalysis includes M, B, P, rho_fluid"
);

console.log();

// ========== GROUP 6: ENHANCED CAPABILITIES (6 tests) ==========
console.log("GROUP 6: Enhanced Capabilities");

// Test 33: State save and restore
sgr.saveState("test_state");
const M_before = sgr.M;
sgr.M *= 2;
sgr.restoreState("test_state");
const M_after = sgr.M;
assert(
    Math.abs(M_after - M_before) / M_before < 0.01,
    "Test 33: State save/restore preserves M"
);

// Test 34: List variables includes key params
const vars = sgr.listVariables();
const has_M_var = vars.includes('M');
const has_B_var = vars.includes('B');
const has_P_var = vars.includes('P');
assert(
    has_M_var && has_B_var && has_P_var && vars.length > 40,
    "Test 34: listVariables includes M, B, P (40+ variables)"
);

// Test 35: Generate variations
const variations = sgr.generateVariations(5, 5.0);
assert(
    variations.length === 5,
    "Test 35: generateVariations creates 5 variants"
);

// Test 36: Validate consistency
const is_valid = sgr.validateConsistency();
assert(
    is_valid === true,
    "Test 36: validateConsistency returns true for valid state"
);

// Test 37: Get system name
const name = sgr.getSystemName();
assert(
    name === "SGR1745-2900_Magnetar",
    "Test 37: getSystemName returns correct identifier"
);

// Test 38: Generate report includes key info
const report = sgr.generateReport(t_test);
const has_mass = report.includes("Mass");
const has_field = report.includes("B-field");
const has_em = report.includes("EM");
assert(
    has_mass && has_field && has_em,
    "Test 38: generateReport includes Mass, B-field, EM"
);

console.log();

// ========== SUMMARY ==========
console.log("=".repeat(70));
console.log(`RESULTS: ${passed} passed, ${failed} failed out of 38 tests`);
if (failed === 0) {
    console.log("✓ ALL TESTS PASSED");
} else {
    console.log(`✗ ${failed} TEST(S) FAILED`);
}
console.log("=".repeat(70));

process.exit(failed > 0 ? 1 : 0);
