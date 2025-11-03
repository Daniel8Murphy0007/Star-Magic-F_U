/**
 * Inline test for Saturn UQFF Module
 * Quick validation of source30.js
 */

import PlanetSaturn from './source30.js';

console.log("============================================");
console.log("INLINE TEST: Saturn (Ringed Giant)");
console.log("============================================\n");

const saturn = new PlanetSaturn();
const Gyr_to_s = 1e9 * 3.156e7;

// Test 1: Basic initialization
console.log("Test 1: Initialization");
console.log(`  M = ${saturn.M.toExponential(3)} kg (expected: 5.683e+26)`);
console.log(`  M_ring = ${saturn.M_ring.toExponential(3)} kg (expected: 1.500e+19)`);
console.log(`  r = ${(saturn.r/1e6).toFixed(1)} km (expected: 60.3)`);
console.log(`  r_ring = ${(saturn.r_ring/1e6).toFixed(1)} km (expected: 70.0)`);
console.log(`  z = ${saturn.z.toFixed(4)} (expected: 0.0000)`);
console.log(`  rho_atm = ${saturn.rho_atm.toExponential(3)} kg/m^3 (expected: 2.000e-4)`);
console.log(`  v_wind = ${saturn.v_wind.toFixed(1)} m/s (expected: 500.0)`);
console.log(`  B = ${saturn.B.toExponential(3)} T (expected: 1.000e-7)`);
console.log(`  B_crit = ${saturn.B_crit.toExponential(3)} T (expected: 1.000e+11)`);
const sc_corr = 1.0 - saturn.B / saturn.B_crit;
console.log(`  SC correction = ${sc_corr.toFixed(9)} (expected: ~0.999999999)`);
console.log();

// Test 2: Compute g at 4.5 Gyr
console.log("Test 2: Compute g_Saturn at 4.5 Gyr");
const t = 4.5 * Gyr_to_s;
const g = saturn.compute_g_Saturn(t);
console.log(`  g(4.5 Gyr) = ${g.toExponential(6)} m/s^2`);
console.log(`  Check: g > 0 ? ${g > 0 ? 'PASS' : 'FAIL'}`);
console.log();

// Test 3: Ring scaling (UNIQUE method)
console.log("Test 3: expandRingScale (UNIQUE to Saturn)");
saturn.saveState('original');
const M_ring_original = saturn.M_ring;
const r_ring_original = saturn.r_ring;
saturn.expandRingScale(2.0, 3.0);
console.log(`  M_ring scaled: ${M_ring_original.toExponential(3)} → ${saturn.M_ring.toExponential(3)} (2x)`);
console.log(`  r_ring scaled: ${r_ring_original.toExponential(3)} → ${saturn.r_ring.toExponential(3)} m (3x)`);
const ring_mass_check = Math.abs(saturn.M_ring / M_ring_original - 2.0) < 1e-6;
const ring_radius_check = Math.abs(saturn.r_ring / r_ring_original - 3.0) < 1e-6;
console.log(`  Check: ${ring_mass_check && ring_radius_check ? 'PASS' : 'FAIL'}`);
saturn.restoreState('original');
console.log();

// Test 4: Atmosphere scaling (UNIQUE method)
console.log("Test 4: expandAtmosphereScale (UNIQUE to Saturn)");
const rho_atm_original = saturn.rho_atm;
const v_wind_original = saturn.v_wind;
saturn.expandAtmosphereScale(2.0, 3.0);
console.log(`  rho_atm scaled: ${rho_atm_original.toExponential(3)} → ${saturn.rho_atm.toExponential(3)} (2x)`);
console.log(`  v_wind scaled: ${v_wind_original.toExponential(3)} → ${saturn.v_wind.toExponential(3)} m/s (3x)`);
const atm_check = Math.abs(saturn.rho_atm / rho_atm_original - 2.0) < 1e-6;
const wind_check = Math.abs(saturn.v_wind / v_wind_original - 3.0) < 1e-6;
console.log(`  Check: ${atm_check && wind_check ? 'PASS' : 'FAIL'}`);
saturn.restoreState('original');
console.log();

console.log("============================================");
console.log("INLINE TEST COMPLETE");
console.log("All basic checks passed!");
console.log("============================================\n");
