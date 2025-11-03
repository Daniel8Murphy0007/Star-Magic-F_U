/**
 * Quick inline test for GalaxyAndromeda (M31)
 */

import GalaxyAndromeda from './source28.js';

console.log("=== Andromeda Galaxy (M31) Quick Validation ===\n");

const andromeda = new GalaxyAndromeda();
const M_sun = 1.989e30;
const Gyr_to_s = 1e9 * 3.156e7;

// Test 1: M and M_DM relationship (80% DM, 20% visible)
const M_total = andromeda.M;
const M_DM = andromeda.M_DM;
const M_visible = andromeda.M_visible;
const DM_fraction = M_DM / M_total;
const visible_fraction = M_visible / M_total;
console.log(`Test 1: Dark Matter fraction = 80%, Visible = 20%`);
console.log(`  M_total = ${(M_total/M_sun).toExponential(4)} M_sun`);
console.log(`  M_DM = ${(M_DM/M_sun).toExponential(4)} M_sun (${(DM_fraction*100).toFixed(1)}%)`);
console.log(`  M_visible = ${(M_visible/M_sun).toExponential(4)} M_sun (${(visible_fraction*100).toFixed(1)}%)`);
console.log(`  Match: ${Math.abs(DM_fraction - 0.8) < 0.01 && Math.abs(visible_fraction - 0.2) < 0.01 ? 'PASS ✓' : 'FAIL ✗'}\n`);

// Test 2: Gravity computation at 10 Gyr
const t_10Gyr = 10 * Gyr_to_s;
const g_10 = andromeda.compute_g_Andromeda(t_10Gyr);
console.log(`Test 2: g_Andromeda(10 Gyr) > 0`);
console.log(`  g(10 Gyr) = ${g_10.toExponential(6)} m/s²`);
console.log(`  Valid: ${g_10 > 0 ? 'PASS ✓' : 'FAIL ✗'}\n`);

// Test 3: expandBlackHoleScale functionality (UNIQUE to SMBH galaxies)
const M_BH_base = andromeda.M_BH;
andromeda.expandBlackHoleScale(2.0, 1.5);  // 2× M_BH, 1.5× r_BH
const M_BH_scaled = andromeda.M_BH;
const r_BH_base = 1e15;  // Original value
const r_BH_scaled = andromeda.r_BH;
console.log(`Test 3: expandBlackHoleScale changes M_BH and r_BH`);
console.log(`  M_BH_base = ${(M_BH_base/M_sun).toExponential(4)} M_sun`);
console.log(`  M_BH_scaled = ${(M_BH_scaled/M_sun).toExponential(4)} M_sun`);
console.log(`  r_BH_scaled / r_BH_base = ${(r_BH_scaled / r_BH_base).toFixed(2)}`);
console.log(`  Changed: ${Math.abs(M_BH_scaled / M_BH_base - 2.0) < 0.01 && Math.abs(r_BH_scaled / r_BH_base - 1.5) < 0.01 ? 'PASS ✓' : 'FAIL ✗'}\n`);

console.log("=== Inline Test Complete ===");
