/**
 * Inline Test for NebulaCrab (Crab Nebula M1) - Quick Validation
 */

import NebulaCrab from './source32.js';

console.log("=== INLINE TEST: Crab Nebula M1 ===\n");

// Test 1: Initialization
const crab = new NebulaCrab();
console.log("Test 1: Initialization");
console.log(`  M = ${(crab.M / crab.M_sun).toFixed(1)} M☉ (expect 4.6 M☉)`);
console.log(`  r0 = ${(crab.r0 / 9.461e15).toFixed(2)} ly (expect ~5.5 ly)`);
console.log(`  v_exp = ${(crab.v_exp / 1e3).toFixed(1)} km/s (expect 1500 km/s)`);
console.log(`  P_pulsar = ${crab.P_pulsar.toExponential(2)} W (expect 5e31 W)`);
console.log(`  v_shock = ${(crab.v_shock / 1e3).toFixed(1)} km/s (expect 1500 km/s)`);
console.log(`  Age = ${(crab.t / crab.year_to_s).toFixed(0)} years (expect 971 years)\n`);

// Test 2: Radius expansion at 971 years
const t_971yr = 971 * crab.year_to_s;
const r_971 = crab.computeRadius(t_971yr);
console.log("Test 2: Radius Expansion");
console.log(`  r(971 yr) = ${(r_971 / 9.461e15).toFixed(2)} ly (expect r0 + v_exp × 971 yr)\n`);

// Test 3: Wind term dominance at 971 years
const wind = crab.computeWindTerm(r_971);
console.log("Test 3: Pulsar Wind Term");
console.log(`  a_wind(971 yr) = ${wind.toExponential(6)} m/s² (expect dominant term)\n`);

// Test 4: Magnetic term at 971 years
const mag = crab.computeMagTerm();
console.log("Test 4: Magnetic Term");
console.log(`  M_mag = ${mag.toExponential(6)} m/s² (expect positive)\n`);

// Test 5: Total gravity at 971 years
const g = crab.compute_g_Crab(t_971yr);
console.log("Test 5: Total Gravity at 971 years");
console.log(`  g_Crab(971 yr) = ${g.toExponential(6)} m/s² (expect positive, wind dominant)\n`);

// Test 6: expandPulsarScale (UNIQUE to Crab)
const P_before = crab.P_pulsar;
const v_before = crab.v_exp;
crab.expandPulsarScale(2.0, 3.0);
console.log("Test 6: expandPulsarScale(2.0, 3.0)");
console.log(`  P_pulsar: ${(P_before).toExponential(2)} → ${(crab.P_pulsar).toExponential(2)} W (expect 2× scale)`);
console.log(`  v_exp: ${(v_before / 1e3).toFixed(1)} → ${(crab.v_exp / 1e3).toFixed(1)} km/s (expect 3× scale)\n`);

// Test 7: expandMagneticScale (UNIQUE to Crab)
crab.initializeDefaults();  // Reset
const B_before = crab.B;
const v_shock_before = crab.v_shock;
crab.expandMagneticScale(2.0, 3.0);
console.log("Test 7: expandMagneticScale(2.0, 3.0)");
console.log(`  B: ${B_before.toExponential(3)} → ${crab.B.toExponential(3)} T (expect 2× scale)`);
console.log(`  v_shock: ${(v_shock_before / 1e3).toFixed(1)} → ${(crab.v_shock / 1e3).toFixed(1)} km/s (expect 3× scale)\n`);

console.log("=== INLINE TEST COMPLETE ===");
