/**
 * Quick inline test for GalaxyNGC1792 (NGC 1792 Starburst Galaxy)
 */

import GalaxyNGC1792 from './source27.js';

console.log("=== NGC 1792 Quick Validation ===\n");

const ngc = new GalaxyNGC1792();
const M_sun = 1.989e30;
const Myr_to_s = 1e6 * 3.156e7;

// Test 1: Initial M(t) with SFR boost
const M0 = ngc.M0;
const M_at_0 = ngc.M_t(0);
const SFR_factor_check = ngc.SFR_factor;
const expected_M0_boost = M0 * (1 + SFR_factor_check);
console.log(`Test 1: M(0) = M0 × (1 + SFR_factor)`);
console.log(`  M0 = ${(M0/M_sun).toExponential(4)} M_sun`);
console.log(`  SFR_factor = ${SFR_factor_check.toExponential(4)}`);
console.log(`  M(0) = ${(M_at_0/M_sun).toExponential(4)} M_sun`);
console.log(`  Expected = ${(expected_M0_boost/M_sun).toExponential(4)} M_sun`);
console.log(`  Match: ${Math.abs(M_at_0/expected_M0_boost - 1.0) < 1e-6 ? 'PASS ✓' : 'FAIL ✗'}\n`);

// Test 2: Gravity computation at 50 Myr
const t_50Myr = 50 * Myr_to_s;
const g_50 = ngc.compute_g_NGC1792(t_50Myr);
console.log(`Test 2: g_NGC1792(50 Myr) > 0`);
console.log(`  g(50 Myr) = ${g_50.toExponential(6)} m/s^2`);
console.log(`  Valid: ${g_50 > 0 ? 'PASS ✓' : 'FAIL ✗'}\n`);

// Test 3: expandStarburstScale functionality (UNIQUE to starburst)
const g_base = g_50;
ngc.expandStarburstScale(2.0, 1.5);  // 2× SFR_factor, 1.5× tau_SF
const g_scaled = ngc.compute_g_NGC1792(t_50Myr);
console.log(`Test 3: expandStarburstScale changes g`);
console.log(`  g_base = ${g_base.toExponential(6)} m/s^2`);
console.log(`  g_scaled = ${g_scaled.toExponential(6)} m/s^2`);
console.log(`  Changed: ${Math.abs(g_scaled - g_base) > 1e-20 ? 'PASS ✓' : 'FAIL ✗'}\n`);

console.log("=== Inline Test Complete ===");
