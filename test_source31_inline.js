/**
 * Inline Test for NebulaM16 (Eagle Nebula M16) - Quick Validation
 */

import NebulaM16 from './source31.js';

console.log("=== INLINE TEST: Eagle Nebula M16 ===\n");

// Test 1: Initialization
const m16 = new NebulaM16();
console.log("Test 1: Initialization");
console.log(`  M = ${(m16.M / m16.M_sun).toFixed(1)} M☉ (expect 1200 M☉)`);
console.log(`  r = ${(m16.r / 9.461e15).toFixed(1)} ly (expect ~35 ly)`);
console.log(`  SFR = ${(m16.SFR / m16.M_sun).toFixed(1)} M☉/yr (expect 1 M☉/yr)`);
console.log(`  tau_erode = ${(m16.tau_erode_yr / 1e6).toFixed(1)} Myr (expect 3 Myr)`);
console.log(`  v_gas = ${(m16.v_gas / 1e3).toFixed(1)} km/s (expect 100 km/s)`);
console.log(`  z = ${m16.z.toFixed(4)} (expect 0.0015)\n`);

// Test 2: Star formation factor at 5 Myr
const t_5Myr = 5e6 * m16.year_to_s;
const msf = m16.computeMsfFactor(t_5Myr);
console.log("Test 2: Star Formation Factor");
console.log(`  M_sf(5 Myr) = ${msf.toFixed(6)} (expect 5×10⁶/1200 = 4166.67)\n`);

// Test 3: Erosion factor at 5 Myr
const e_rad = m16.computeE_rad(t_5Myr);
console.log("Test 3: Radiation Erosion Factor");
console.log(`  E_rad(5 Myr) = ${e_rad.toFixed(6)} (expect ~0.3 × (1 - exp(-5/3)) = 0.243)\n`);

// Test 4: Total gravity at 5 Myr
const g = m16.compute_g_M16(t_5Myr);
console.log("Test 4: Total Gravity at 5 Myr");
console.log(`  g_M16(5 Myr) = ${g.toExponential(6)} m/s² (expect positive, EM likely dominant)\n`);

// Test 5: expandStarFormationScale (UNIQUE to M16)
const SFR_before = m16.SFR;
const tau_before = m16.tau_erode_yr;
m16.expandStarFormationScale(2.0, 3.0);
console.log("Test 5: expandStarFormationScale(2.0, 3.0)");
console.log(`  SFR: ${(SFR_before / m16.M_sun).toFixed(1)} → ${(m16.SFR / m16.M_sun).toFixed(1)} M☉/yr (expect 1 → 2)`);
console.log(`  tau_erode: ${(tau_before / 1e6).toFixed(1)} → ${(m16.tau_erode_yr / 1e6).toFixed(1)} Myr (expect 3 → 9)\n`);

// Test 6: expandGasScale (UNIQUE to M16)
m16.initializeDefaults();  // Reset
const rho_before = m16.rho_fluid;
const v_before = m16.v_gas;
m16.expandGasScale(2.0, 3.0);
console.log("Test 6: expandGasScale(2.0, 3.0)");
console.log(`  rho_fluid: ${rho_before.toExponential(3)} → ${m16.rho_fluid.toExponential(3)} kg/m³ (expect 2× scale)`);
console.log(`  v_gas: ${(v_before / 1e3).toFixed(1)} → ${(m16.v_gas / 1e3).toFixed(1)} km/s (expect 100 → 300)\n`);

console.log("=== INLINE TEST COMPLETE ===");
