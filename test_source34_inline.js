// test_source34_inline.js
// Quick inline tests for SGR 1745-2900 Frequency/Resonance module (source34.js)

import { MagnetarSGR1745Frequency } from './source34.js';

console.log("========== SGR 1745-2900 FREQUENCY/RESONANCE INLINE TESTS ==========\n");

const sgr = new MagnetarSGR1745Frequency();
const t_test = 1000 * 3.156e7; // 1000 years

// Test 1: Initialization
console.log("Test 1: Initialization");
console.log(`  M = ${(sgr.M / sgr.M_sun).toFixed(2)} M_sun (expected: 1.50)`);
console.log(`  r = ${(sgr.r / 1e3).toFixed(2)} km (expected: 10.00)`);
console.log(`  I = ${sgr.I.toExponential(2)} A (expected: 1.00e+21)`);
console.log(`  f_DPM = ${(sgr.f_DPM / 1e12).toFixed(3)} THz (expected: 1.000)`);
console.log(`  f_THz = ${(sgr.f_THz / 1e12).toFixed(3)} THz (expected: 1.000)`);
const pass1 = Math.abs(sgr.M / sgr.M_sun - 1.5) < 0.01 &&
              Math.abs(sgr.r - 1e4) < 1 &&
              Math.abs(sgr.I - 1e21) < 1e19 &&
              Math.abs(sgr.f_DPM - 1e12) < 1e10;
console.log(`  Result: ${pass1 ? "PASS" : "FAIL"}\n`);

// Test 2: DPM term computation
console.log("Test 2: DPM term (base frequency-driven acceleration)");
const a_DPM = sgr.computeDPMTerm();
console.log(`  a_DPM = ${a_DPM.toExponential(6)} m/s²`);
console.log(`  a_DPM > 0? ${a_DPM > 0 ? "Yes" : "No"}`);
const pass2 = a_DPM > 0 && a_DPM < 1e-15;  // Should be tiny (micro-scale)
console.log(`  Result: ${pass2 ? "PASS" : "FAIL"}\n`);

// Test 3: THz term (should dominate per documentation)
console.log("Test 3: THz term dominance");
const a_THz = sgr.computeTHzTerm();
console.log(`  a_THz = ${a_THz.toExponential(6)} m/s²`);
console.log(`  a_THz > a_DPM? ${a_THz > a_DPM ? "Yes" : "No"}`);
const pass3 = a_THz > a_DPM;  // THz should be larger
console.log(`  Result: ${pass3 ? "PASS" : "FAIL"}\n`);

// Test 4: Full g computation (galactic-scale due to U_g4i)
console.log("Test 4: Full g_SGR1745 computation");
const g_total = sgr.compute_g_SGR1745(t_test);
console.log(`  g_UQFF = ${g_total.toExponential(6)} m/s²`);
console.log(`  Expected order: ~10³² m/s² (galactic-scale, U_g4i dominates)`);
const pass4 = g_total > 1e30 && g_total < 1e34;  // Galactic-scale
console.log(`  Result: ${pass4 ? "PASS" : "FAIL"}\n`);

// Test 5: expandDPMScale
console.log("Test 5: expandDPMScale (2.0x current, 1.5x frequency)");
const I_orig = sgr.I;
const f_DPM_orig = sgr.f_DPM;
sgr.expandDPMScale(2.0, 1.5);
const I_new = sgr.I;
const f_DPM_new = sgr.f_DPM;
console.log(`  Original I = ${I_orig.toExponential(2)} A`);
console.log(`  New I = ${I_new.toExponential(2)} A (expected: ${(I_orig * 2.0).toExponential(2)})`);
console.log(`  Original f_DPM = ${(f_DPM_orig / 1e12).toFixed(3)} THz`);
console.log(`  New f_DPM = ${(f_DPM_new / 1e12).toFixed(3)} THz (expected: ${(f_DPM_orig * 1.5 / 1e12).toFixed(3)})`);
const pass5 = Math.abs(I_new / I_orig - 2.0) < 0.01 && Math.abs(f_DPM_new / f_DPM_orig - 1.5) < 0.01;
console.log(`  Result: ${pass5 ? "PASS" : "FAIL"}\n`);

// Reset
sgr.I = I_orig;
sgr.f_DPM = f_DPM_orig;

// Test 6: expandFrequencyScale
console.log("Test 6: expandFrequencyScale (1.5x THz, 2.0x super)");
const f_THz_orig = sgr.f_THz;
const f_super_orig = sgr.f_super;
sgr.expandFrequencyScale(1.5, 2.0);
const f_THz_new = sgr.f_THz;
const f_super_new = sgr.f_super;
console.log(`  Original f_THz = ${(f_THz_orig / 1e12).toFixed(3)} THz`);
console.log(`  New f_THz = ${(f_THz_new / 1e12).toFixed(3)} THz (expected: ${(f_THz_orig * 1.5 / 1e12).toFixed(3)})`);
console.log(`  f_super scaled by 2.0x`);
const pass6 = Math.abs(f_THz_new / f_THz_orig - 1.5) < 0.01 && Math.abs(f_super_new / f_super_orig - 2.0) < 0.01;
console.log(`  Result: ${pass6 ? "PASS" : "FAIL"}\n`);

// Reset
sgr.f_THz = f_THz_orig;
sgr.f_super = f_super_orig;

// Test 7: expandVacuumScale
console.log("Test 7: expandVacuumScale (1.3x E_vac, 2.0x v_exp)");
const E_vac_orig = sgr.E_vac_neb;
const v_exp_orig = sgr.v_exp;
sgr.expandVacuumScale(1.3, 2.0);
const E_vac_new = sgr.E_vac_neb;
const v_exp_new = sgr.v_exp;
console.log(`  Original E_vac_neb = ${E_vac_orig.toExponential(3)} J/m³`);
console.log(`  New E_vac_neb = ${E_vac_new.toExponential(3)} J/m³`);
console.log(`  Original v_exp = ${(v_exp_orig / 1e3).toFixed(3)} km/s`);
console.log(`  New v_exp = ${(v_exp_new / 1e3).toFixed(3)} km/s (expected: ${(v_exp_orig * 2.0 / 1e3).toFixed(3)})`);
const pass7 = Math.abs(E_vac_new / E_vac_orig - 1.3) < 0.01 && Math.abs(v_exp_new / v_exp_orig - 2.0) < 0.01;
console.log(`  Result: ${pass7 ? "PASS" : "FAIL"}\n`);

// Summary
const all_passed = pass1 && pass2 && pass3 && pass4 && pass5 && pass6 && pass7;
console.log("========================================");
console.log(`All tests passed: ${all_passed ? "YES" : "NO"}`);
console.log(`Pass rate: ${[pass1, pass2, pass3, pass4, pass5, pass6, pass7].filter(p => p).length}/7`);
console.log("========================================\n");
