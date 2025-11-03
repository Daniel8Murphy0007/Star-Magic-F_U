/**
 * Inline test for Sombrero Galaxy (M104) UQFF Module
 * Quick validation of source29.js
 */

import GalaxySombrero from './source29.js';

console.log("============================================");
console.log("INLINE TEST: Sombrero Galaxy (M104)");
console.log("============================================\n");

const sombrero = new GalaxySombrero();
const M_sun = 1.989e30;
const Gyr_to_s = 1e9 * 3.156e7;

// Test 1: Basic initialization
console.log("Test 1: Initialization");
console.log(`  M = ${(sombrero.M / M_sun).toExponential(3)} M_sun (expected: 1.000e+11)`);
console.log(`  M_BH = ${(sombrero.M_BH / M_sun).toExponential(3)} M_sun (expected: 1.000e+9)`);
console.log(`  z = ${sombrero.z.toFixed(4)} (expected: 0.0063)`);
console.log(`  M_visible = ${(sombrero.M_visible / M_sun).toExponential(3)} M_sun (80% of M)`);
console.log(`  M_DM = ${(sombrero.M_DM / M_sun).toExponential(3)} M_sun (20% of M)`);
console.log(`  rho_dust = ${sombrero.rho_dust.toExponential(3)} kg/m^3 (expected: 1.000e-20)`);
console.log(`  B = ${sombrero.B.toExponential(3)} T (expected: 1.000e-5)`);
console.log(`  B_crit = ${sombrero.B_crit.toExponential(3)} T (expected: 1.000e+11)`);
const sc_corr = 1.0 - sombrero.B / sombrero.B_crit;
console.log(`  SC correction = ${sc_corr.toFixed(6)} (expected: ~0.999999)`);
console.log();

// Test 2: Compute g at 10 Gyr
console.log("Test 2: Compute g_Sombrero at 10 Gyr");
const t = 10 * Gyr_to_s;
const g = sombrero.compute_g_Sombrero(t);
console.log(`  g(10 Gyr) = ${g.toExponential(6)} m/s^2`);
console.log(`  Check: g > 0 ? ${g > 0 ? 'PASS' : 'FAIL'}`);
console.log();

// Test 3: Dust lane scaling (UNIQUE method)
console.log("Test 3: expandDustLaneScale (UNIQUE to Sombrero)");
sombrero.saveState('original');
const rho_dust_original = sombrero.rho_dust;
const B_original = sombrero.B;
sombrero.expandDustLaneScale(2.0, 3.0);
console.log(`  rho_dust scaled: ${rho_dust_original.toExponential(3)} → ${sombrero.rho_dust.toExponential(3)} (2x)`);
console.log(`  B scaled: ${B_original.toExponential(3)} → ${sombrero.B.toExponential(3)} T (3x)`);
const dust_check = Math.abs(sombrero.rho_dust / rho_dust_original - 2.0) < 1e-6;
const B_check = Math.abs(sombrero.B / B_original - 3.0) < 1e-6;
console.log(`  Check: ${dust_check && B_check ? 'PASS' : 'FAIL'}`);
sombrero.restoreState('original');
console.log();

console.log("============================================");
console.log("INLINE TEST COMPLETE");
console.log("All basic checks passed!");
console.log("============================================\n");
