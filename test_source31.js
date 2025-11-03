/**
 * Comprehensive Test Suite for NebulaM16 (Eagle Nebula M16)
 * Tests all 38 methods (13 core + 25 enhanced) and unique star formation/erosion features
 */

import NebulaM16 from './source31.js';

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
console.log("COMPREHENSIVE TEST SUITE: Eagle Nebula M16");
console.log("========================================\n");

const m16 = new NebulaM16();
const Myr_to_s = 1e6 * m16.year_to_s;
const M_sun = m16.M_sun;

// ========== GROUP 1: INITIALIZATION & BASIC PARAMETERS (8 tests) ==========
console.log("GROUP 1: Initialization & Basic Parameters\n");

// Test 1: Mass within range
assert(m16.M > 1000 * M_sun && m16.M < 1500 * M_sun, 
    "M16 mass within range",
    `M = ${m16.M/M_sun} M☉`);

// Test 2: Radius within range  
assert(m16.r > 3e17 && m16.r < 4e17,
    "M16 radius within range",
    `r = ${m16.r/9.461e15} ly`);

// Test 3: Star formation rate
assert(Math.abs(m16.SFR / M_sun - 1.0) < 0.01,
    "SFR = 1 M☉/yr",
    `SFR = ${m16.SFR/M_sun} M☉/yr`);

// Test 4: Erosion timescale
assert(Math.abs(m16.tau_erode_yr / 1e6 - 3.0) < 0.01,
    "tau_erode = 3 Myr",
    `tau_erode = ${m16.tau_erode_yr/1e6} Myr`);

// Test 5: Redshift for nearby nebula
assert(Math.abs(m16.z - 0.0015) < 0.0001,
    "z = 0.0015 (nearby)",
    `z = ${m16.z}`);

// Test 6: Gas velocity
assert(Math.abs(m16.v_gas / 1e5 - 1.0) < 0.01,
    "v_gas = 100 km/s",
    `v_gas = ${m16.v_gas/1e3} km/s`);

// Test 7: Gas density
assert(m16.rho_fluid === 1e-20,
    "rho_fluid = 1e-20 kg/m³",
    `rho_fluid = ${m16.rho_fluid}`);

// Test 8: No dark matter in nebula
assert(m16.M_DM === 0 && m16.M_visible === m16.M,
    "No DM, M_visible = M",
    `M_DM = ${m16.M_DM}, M_visible = ${m16.M_visible/M_sun} M☉`);

// ========== GROUP 2: TIME-DEPENDENT MASS EVOLUTION (6 tests) ==========
console.log("\nGROUP 2: Time-Dependent Mass Evolution (Star Formation & Erosion)\n");

// Test 9: M_sf grows linearly with time
const t1 = 1 * Myr_to_s;
const t2 = 2 * Myr_to_s;
const msf1 = m16.computeMsfFactor(t1);
const msf2 = m16.computeMsfFactor(t2);
const msf_ratio = msf2 / msf1;
assert(Math.abs(msf_ratio - 2.0) < 0.01,
    "M_sf scales linearly with time",
    `M_sf(2Myr)/M_sf(1Myr) = ${msf_ratio.toFixed(3)}`);

// Test 10: M_sf at 5 Myr matches expectation
const t5 = 5 * Myr_to_s;
const msf5 = m16.computeMsfFactor(t5);
const expected_msf5 = 5e6 / 1200;  // 5 Myr / 1200 M☉
assert(Math.abs(msf5 - expected_msf5) < 1.0,
    "M_sf(5 Myr) = 4166.67",
    `M_sf = ${msf5.toFixed(2)}`);

// Test 11: E_rad approaches E_0 asymptotically
const t_large = 30 * Myr_to_s;  // 30 Myr >> 3 Myr tau
const e_rad_large = m16.computeE_rad(t_large);
assert(Math.abs(e_rad_large - m16.E_0) < 0.01,
    "E_rad approaches E_0 = 0.3 at large t",
    `E_rad(30 Myr) = ${e_rad_large.toFixed(3)}`);

// Test 12: E_rad at 5 Myr matches formula
const e_rad5 = m16.computeE_rad(t5);
const expected_e_rad5 = 0.3 * (1 - Math.exp(-5/3));
assert(Math.abs(e_rad5 - expected_e_rad5) < 0.001,
    "E_rad(5 Myr) = 0.243",
    `E_rad = ${e_rad5.toFixed(6)}`);

// Test 13: E_rad increases with time
const e_rad1 = m16.computeE_rad(t1);
const e_rad2 = m16.computeE_rad(t2);
assert(e_rad2 > e_rad1,
    "E_rad increases with time",
    `E_rad(1Myr) = ${e_rad1.toFixed(3)}, E_rad(2Myr) = ${e_rad2.toFixed(3)}`);

// Test 14: Combined mass factor m_factor
const m_factor = (1 + msf5) * (1 - e_rad5);
assert(m_factor > 1.0,
    "m_factor > 1 (SF dominates over erosion at 5 Myr)",
    `m_factor = ${m_factor.toFixed(3)}`);

// ========== GROUP 3: CORE PHYSICS COMPUTATIONS (8 tests) ==========
console.log("\nGROUP 3: Core Physics Computations\n");

// Test 15: Hz near zero for small z
const Hz = m16.computeHz();
assert(Hz > 0 && Hz < 1e-16,
    "Hz ≈ 0 for z = 0.0015",
    `Hz = ${Hz.toExponential(3)} s⁻¹`);

// Test 16: Ug sum positive
const ug_sum = m16.computeUgSum();
assert(ug_sum > 0,
    "Ug sum > 0",
    `Ug = ${ug_sum.toExponential(3)} m/s²`);

// Test 17: Quantum term small
const quantum = m16.computeQuantumTerm(m16.t_Hubble);
assert(quantum > 0 && quantum < 1e-5,
    "Quantum term small",
    `quantum = ${quantum.toExponential(3)} m/s²`);

// Test 18: Fluid term depends on rho_fluid
const g_base_test = 1e-3;
const fluid1 = m16.computeFluidTerm(g_base_test);
m16.setVariable('rho_fluid', 2e-20);
const fluid2 = m16.computeFluidTerm(g_base_test);
m16.setVariable('rho_fluid', 1e-20);  // Reset
const fluid_ratio = fluid2 / fluid1;
assert(Math.abs(fluid_ratio - 2.0) < 0.01,
    "Fluid term scales with rho_fluid",
    `fluid ratio = ${fluid_ratio.toFixed(3)}`);

// Test 19: Resonant term oscillates
const resonant1 = m16.computeResonantTerm(0);
const resonant2 = m16.computeResonantTerm(Math.PI / m16.omega);
assert(resonant1 !== resonant2,
    "Resonant term oscillates with time",
    `resonant(0) = ${resonant1.toExponential(3)}, resonant(π/ω) = ${resonant2.toExponential(3)}`);

// Test 20: DM term reflects delta_rho/rho (even with M_DM = 0)
const dm_term = m16.computeDMTerm() / m16.M;
assert(Math.abs(dm_term - 0.1) < 0.01,
    "DM term ≈ δρ/ρ = 0.1 (density perturbation)",
    `DM term = ${dm_term.toExponential(3)} m/s²`);

// Test 21: SC correction ratio behavior
const B_orig = m16.B;
const B_crit = m16.B_crit;
const ratio_orig = B_orig / B_crit;
m16.setVariable('B', 2 * B_orig);
const ratio_scaled = m16.B / B_crit;
const ratio_change = ratio_scaled / ratio_orig;
m16.setVariable('B', B_orig);  // Reset
assert(Math.abs(ratio_change - 2.0) < 0.01,
    "B/B_crit ratio doubles when B doubles",
    `ratio change = ${ratio_change.toFixed(3)}`);

// Test 22: g_M16 positive at all times
const times = [0, 1e6 * m16.year_to_s, 5e6 * m16.year_to_s, 10e6 * m16.year_to_s];
const all_positive = times.every(t => m16.compute_g_M16(t) > 0);
assert(all_positive,
    "g_M16(t) > 0 for all times",
    `Tested ${times.length} time points`);

// ========== GROUP 4: UNIQUE EXPANSION METHODS (6 tests) ==========
console.log("\nGROUP 4: Unique Expansion Methods (Star Formation & Gas)\n");

// Test 23: expandStarFormationScale - SFR scaling
m16.saveState('test23');
const SFR_orig = m16.SFR;
m16.expandStarFormationScale(3.0, 1.0);
assert(Math.abs(m16.SFR / SFR_orig - 3.0) < 0.01,
    "expandStarFormationScale scales SFR",
    `SFR scaled ${(m16.SFR/SFR_orig).toFixed(2)}×`);
m16.restoreState('test23');

// Test 24: expandStarFormationScale - tau_erode scaling
m16.saveState('test24');
const tau_orig = m16.tau_erode_yr;
m16.expandStarFormationScale(1.0, 2.5);
assert(Math.abs(m16.tau_erode_yr / tau_orig - 2.5) < 0.01,
    "expandStarFormationScale scales tau_erode",
    `tau_erode scaled ${(m16.tau_erode_yr/tau_orig).toFixed(2)}×`);
m16.restoreState('test24');

// Test 25: expandGasScale - rho_fluid scaling
m16.saveState('test25');
const rho_orig = m16.rho_fluid;
m16.expandGasScale(4.0, 1.0);
assert(Math.abs(m16.rho_fluid / rho_orig - 4.0) < 0.01,
    "expandGasScale scales rho_fluid",
    `rho_fluid scaled ${(m16.rho_fluid/rho_orig).toFixed(2)}×`);
m16.restoreState('test25');

// Test 26: expandGasScale - v_gas scaling
m16.saveState('test26');
const v_orig = m16.v_gas;
m16.expandGasScale(1.0, 2.0);
assert(Math.abs(m16.v_gas / v_orig - 2.0) < 0.01,
    "expandGasScale scales v_gas",
    `v_gas scaled ${(m16.v_gas/v_orig).toFixed(2)}×`);
m16.restoreState('test26');

// Test 27: expandNebulaScale affects M and r
m16.saveState('test27');
const M_orig = m16.M;
const r_orig = m16.r;
m16.expandNebulaScale(1.5, 2.0);
const M_ratio = m16.M / M_orig;
const r_ratio = m16.r / r_orig;
assert(Math.abs(M_ratio - 1.5) < 0.01 && Math.abs(r_ratio - 2.0) < 0.01,
    "expandNebulaScale scales M and r independently",
    `M scaled ${M_ratio.toFixed(2)}×, r scaled ${r_ratio.toFixed(2)}×`);
m16.restoreState('test27');

// Test 28: SFR affects M_sf at fixed time
m16.saveState('test28');
const msf_base = m16.computeMsfFactor(t5);
m16.expandStarFormationScale(2.0, 1.0);
const msf_scaled = m16.computeMsfFactor(t5);
const msf_ratio_test = msf_scaled / msf_base;
m16.restoreState('test28');
assert(Math.abs(msf_ratio_test - 2.0) < 0.01,
    "SFR × 2 → M_sf × 2",
    `M_sf ratio = ${msf_ratio_test.toFixed(3)}`);

// ========== GROUP 5: PARAMETER SENSITIVITY (4 tests) ==========
console.log("\nGROUP 5: Parameter Sensitivity\n");

// Test 29: v_gas affects g (via EM term)
m16.saveState('test29');
const g_base_vgas = m16.compute_g_M16(t5);
m16.setVariable('v_gas', 2e5);  // Double v_gas
const g_high_vgas = m16.compute_g_M16(t5);
m16.restoreState('test29');
assert(g_high_vgas > g_base_vgas,
    "Higher v_gas increases g (EM term)",
    `g increased ${((g_high_vgas/g_base_vgas - 1) * 100).toFixed(1)}%`);

// Test 30: rho_fluid affects fluid term directly
m16.saveState('test30');
const g_base_test30 = 1e-3;  // Fixed g_base for direct fluid term comparison
const fluid_base = m16.computeFluidTerm(g_base_test30);
m16.setVariable('rho_fluid', 2e-20);
const fluid_scaled = m16.computeFluidTerm(g_base_test30);
m16.restoreState('test30');
const fluid_ratio_test30 = fluid_scaled / fluid_base;
assert(Math.abs(fluid_ratio_test30 - 2.0) < 0.01,
    "rho_fluid × 2 → fluid term × 2",
    `fluid ratio = ${fluid_ratio_test30.toFixed(3)}`);

// Test 31: tau_erode affects E_rad
m16.saveState('test31');
const e_rad_base = m16.computeE_rad(t5);
m16.setVariable('tau_erode_yr', 6e6);  // Double tau
const e_rad_long = m16.computeE_rad(t5);
m16.restoreState('test31');
assert(e_rad_long < e_rad_base,
    "Longer tau_erode → slower erosion",
    `E_rad: ${e_rad_base.toFixed(3)} → ${e_rad_long.toFixed(3)}`);

// Test 32: Sensitivity analysis returns valid data
const sensitivity = m16.sensitivityAnalysis(t5, 0.01);
const has_sfr = 'SFR' in sensitivity;
const has_v_gas = 'v_gas' in sensitivity;
assert(has_sfr && has_v_gas,
    "sensitivityAnalysis includes SFR and v_gas",
    `Keys: ${Object.keys(sensitivity).length}`);

// ========== GROUP 6: ENHANCED CAPABILITIES (6 tests) ==========
console.log("\nGROUP 6: Enhanced Capabilities\n");

// Test 33: State save/restore
m16.saveState('test33');
const M_before = m16.M;
m16.setVariable('M', 2000 * M_sun);
m16.restoreState('test33');
assert(m16.M === M_before,
    "State restore works",
    `M restored to ${m16.M/M_sun} M☉`);

// Test 34: Variable list contains key params
const vars = m16.listVariables();
const has_M = vars.includes('M');
const has_SFR = vars.includes('SFR');
const has_tau = vars.includes('tau_erode_yr');
assert(has_M && has_SFR && has_tau,
    "listVariables includes M, SFR, tau_erode_yr",
    `Total variables: ${vars.length}`);

// Test 35: Generate variations
const variations = m16.generateVariations(5, 10.0);
assert(variations.length === 5,
    "generateVariations creates 5 variants",
    `Created ${variations.length} variations`);

// Test 36: Validate consistency
const is_valid = m16.validateConsistency();
assert(is_valid,
    "validateConsistency passes",
    `Valid: ${is_valid}`);

// Test 37: System name
const name = m16.getSystemName();
assert(name === "EagleNebula_M16",
    "getSystemName returns EagleNebula_M16",
    `Name: ${name}`);

// Test 38: Generate report contains key info
const report = m16.generateReport(t5);
const has_msf_in_report = report.includes("Star Formation");
const has_erosion_in_report = report.includes("Erosion");
assert(has_msf_in_report && has_erosion_in_report,
    "generateReport includes SF and erosion data",
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
