// test_source33_inline.js
// Quick inline tests for SGR 1745-2900 Magnetar module (source33.js)

import { MagnetarSGR1745 } from './source33.js';

console.log("========== SGR 1745-2900 MAGNETAR INLINE TESTS ==========\n");

const sgr = new MagnetarSGR1745();
const t_test = 1000 * 3.156e7; // 1000 years

// Test 1: Initialization
console.log("Test 1: Initialization");
console.log(`  M = ${(sgr.M / sgr.M_sun).toFixed(2)} M_sun (expected: 1.40)`);
console.log(`  r = ${(sgr.r / 1e3).toFixed(2)} km (expected: 10.00)`);
console.log(`  B = ${(sgr.B / 1e10).toFixed(2)} × 10^10 T (expected: 2.00)`);
console.log(`  P = ${sgr.P.toFixed(3)} s (expected: 3.760)`);
console.log(`  v_spin = ${(sgr.v_spin / 1e3).toFixed(2)} km/s`);
console.log(`  rho_fluid = ${sgr.rho_fluid.toExponential(2)} kg/m³ (expected: 1.00e+17)`);
const pass1 = Math.abs(sgr.M / sgr.M_sun - 1.4) < 0.01 &&
              Math.abs(sgr.r - 1e4) < 1 &&
              Math.abs(sgr.B - 2e10) < 1e8 &&
              Math.abs(sgr.P - 3.76) < 0.01;
console.log(`  Result: ${pass1 ? "PASS" : "FAIL"}\n`);

// Test 2: Superconductivity correction
console.log("Test 2: Superconductivity correction (1 - B/B_crit)");
const sc_corr = sgr.computeSCCorrection();
const expected_sc = 1.0 - (2e10 / 1e11); // = 0.8
console.log(`  SC correction = ${sc_corr.toFixed(4)} (expected: ${expected_sc.toFixed(4)})`);
const pass2 = Math.abs(sc_corr - expected_sc) < 0.001;
console.log(`  Result: ${pass2 ? "PASS" : "FAIL"}\n`);

// Test 3: EM term contribution
console.log("Test 3: EM term (q v_spin B / m_p × amplification × scale)");
const em_term = sgr.computeEMTerm();
const g_base_raw = (sgr.G * sgr.M) / (sgr.r * sgr.r);
console.log(`  EM term = ${em_term.toExponential(6)} m/s²`);
console.log(`  Base gravity (raw) = ${g_base_raw.toExponential(6)} m/s²`);
console.log(`  EM/base ratio = ${(em_term / g_base_raw).toFixed(4)}`);
const pass3 = em_term > 1e10 && em_term < 1e13; // EM should be significant (1e11 range)
console.log(`  EM in reasonable range (1e10-1e13)? ${pass3 ? "Yes" : "No"}`);
console.log(`  Result: ${pass3 ? "PASS" : "FAIL"}\n`);

// Test 4: Full g computation
console.log("Test 4: Full g_SGR1745 computation at t = 1000 years");
const g_total = sgr.compute_g_SGR1745(t_test);
console.log(`  g_UQFF = ${g_total.toExponential(6)} m/s²`);
const pass4 = g_total > 1e12 && g_total < 1e13; // Expected order of magnitude (e12)
console.log(`  In expected range (1e12-1e13)? ${pass4 ? "Yes" : "No"}`);
console.log(`  Result: ${pass4 ? "PASS" : "FAIL"}\n`);

// Test 5: expandMagnetarScale
console.log("Test 5: expandMagnetarScale (2.0x mass, 1.5x radius)");
const M_orig = sgr.M;
const r_orig = sgr.r;
sgr.expandMagnetarScale(2.0, 1.5);
const M_new = sgr.M;
const r_new = sgr.r;
console.log(`  Original M = ${(M_orig / sgr.M_sun).toFixed(2)} M_sun`);
console.log(`  New M = ${(M_new / sgr.M_sun).toFixed(2)} M_sun (expected: ${(M_orig * 2.0 / sgr.M_sun).toFixed(2)})`);
console.log(`  Original r = ${(r_orig / 1e3).toFixed(2)} km`);
console.log(`  New r = ${(r_new / 1e3).toFixed(2)} km (expected: ${(r_orig * 1.5 / 1e3).toFixed(2)})`);
const pass5 = Math.abs(M_new / M_orig - 2.0) < 0.01 && Math.abs(r_new / r_orig - 1.5) < 0.01;
console.log(`  Result: ${pass5 ? "PASS" : "FAIL"}\n`);

// Reset for next test
sgr.M = M_orig;
sgr.r = r_orig;
sgr.M_visible = M_orig;
sgr.v_spin = (2 * sgr.pi * sgr.r) / sgr.P;

// Test 6: expandMagneticFieldScale
console.log("Test 6: expandMagneticFieldScale (3.0x B-field, 2.0x period)");
const B_orig = sgr.B;
const P_orig = sgr.P;
sgr.expandMagneticFieldScale(3.0, 2.0);
const B_new = sgr.B;
const P_new = sgr.P;
console.log(`  Original B = ${(B_orig / 1e10).toFixed(2)} × 10^10 T`);
console.log(`  New B = ${(B_new / 1e10).toFixed(2)} × 10^10 T (expected: ${(B_orig * 3.0 / 1e10).toFixed(2)})`);
console.log(`  Original P = ${P_orig.toFixed(3)} s`);
console.log(`  New P = ${P_new.toFixed(3)} s (expected: ${(P_orig * 2.0).toFixed(3)})`);
const pass6 = Math.abs(B_new / B_orig - 3.0) < 0.01 && Math.abs(P_new / P_orig - 2.0) < 0.01;
console.log(`  Result: ${pass6 ? "PASS" : "FAIL"}\n`);

// Reset for next test
sgr.B = B_orig;
sgr.P = P_orig;
sgr.omega = 2 * sgr.pi / P_orig;
sgr.v_spin = (2 * sgr.pi * sgr.r) / sgr.P;

// Test 7: expandCrustScale
console.log("Test 7: expandCrustScale (2.5x density, 1.5x volume)");
const rho_orig = sgr.rho_fluid;
const V_orig = sgr.V;
sgr.expandCrustScale(2.5, 1.5);
const rho_new = sgr.rho_fluid;
const V_new = sgr.V;
console.log(`  Original rho = ${rho_orig.toExponential(2)} kg/m³`);
console.log(`  New rho = ${rho_new.toExponential(2)} kg/m³ (expected: ${(rho_orig * 2.5).toExponential(2)})`);
console.log(`  Original V = ${V_orig.toExponential(2)} m³`);
console.log(`  New V = ${V_new.toExponential(2)} m³ (expected: ${(V_orig * 1.5).toExponential(2)})`);
const pass7 = Math.abs(rho_new / rho_orig - 2.5) < 0.01 && Math.abs(V_new / V_orig - 1.5) < 0.01;
console.log(`  Result: ${pass7 ? "PASS" : "FAIL"}\n`);

// Summary
const all_passed = pass1 && pass2 && pass3 && pass4 && pass5 && pass6 && pass7;
console.log("========================================");
console.log(`All tests passed: ${all_passed ? "YES" : "NO"}`);
console.log(`Pass rate: ${[pass1, pass2, pass3, pass4, pass5, pass6, pass7].filter(p => p).length}/7`);
console.log("========================================\n");
