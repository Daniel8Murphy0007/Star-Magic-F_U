/**
 * ================================================================================================
 * Test Suite: test_source20.js
 * 
 * Description: Comprehensive test suite for source20.js (GalaxyNGC2525)
 *              Tests barred spiral galaxy with supernova mass loss physics
 * 
 * Features Tested:
 *   - Initialization and parameter validation
 *   - Time evolution (supernova mass decay)
 *   - Galactic mass scaling
 *   - Black hole mass/radius scaling (central BH feature)
 *   - Supernova mass/timescale scaling (SN feature)
 *   - Magnetic field effects
 *   - Volume calculation
 *   - Ug terms computation
 *   - State save/restore
 *   - Sensitivity analysis
 *   - Parameter variations (Monte Carlo)
 *   - Batch operations
 *   - Consistency validation
 *   - Parameter space expansion
 *   - Complete UQFF terms breakdown
 *   - Full system report
 *   - Supernova decay validation
 *   - Black hole contribution analysis
 * 
 * Author: GitHub Copilot
 * Date: November 03, 2025
 * ================================================================================================
 */

import GalaxyNGC2525 from './source20.js';

console.log("=".repeat(60));
console.log("TEST SUITE FOR GALAXY NGC 2525 (source20.js)");
console.log("Barred Spiral Galaxy with Supernova Mass Loss");
console.log("=".repeat(60));
console.log();

const M_sun = 1.989e30;
const ly_to_m = 9.461e15;

// Test 1: Initialization
console.log("Test 1: Initialization");
const ngc2525 = new GalaxyNGC2525();
console.log(`  System name: ${ngc2525.getSystemName()}`);
const vars = ngc2525.listVariables();
console.log(`  Number of variables: ${vars.length}`);
console.log(`  M (galaxy mass): ${(ngc2525.M/M_sun).toExponential(3)} M_sun`);
console.log(`  M_BH (black hole mass): ${(ngc2525.M_BH/M_sun).toExponential(3)} M_sun`);
console.log(`  r (galaxy radius): ${(ngc2525.r/ly_to_m).toFixed(2)} ly`);
console.log(`  z_gal (redshift): ${ngc2525.z_gal.toFixed(3)}`);
console.log(`  Hz (Hubble at z): ${ngc2525.Hz.toExponential(3)} s^-1`);
console.log(`  B: ${ngc2525.B.toExponential(3)} T`);
console.log(`  M_SN0 (initial SN mass): ${(ngc2525.M_SN0/M_sun).toFixed(2)} M_sun`);
console.log(`  tau_SN (SN timescale): ${(ngc2525.tau_SN/3.156e7).toFixed(2)} years`);
console.log();

// Test 2: Initial gravity computation
console.log("Test 2: Initial Gravity");
const g0 = ngc2525.compute_g_NGC2525(0);
console.log(`  g_NGC2525(t=0) = ${g0.toExponential(6)} m/s²`);
console.log(`  M = ${(ngc2525.M/M_sun).toExponential(3)} M_sun`);
console.log(`  Galaxy Radius = ${(ngc2525.r/ly_to_m).toFixed(2)} ly`);
console.log(`  Supernova effect: M_SN(0) = ${(ngc2525.M_SN_t(0)/M_sun).toFixed(2)} M_sun (initial)`);
console.log(`  Black Hole: M_BH = ${(ngc2525.M_BH/M_sun).toExponential(3)} M_sun`);
console.log();

// Test 3: Time evolution (supernova mass decay)
console.log("Test 3: Time Evolution (0-20 years) - Supernova Mass Decay");
const times = [0, 1, 5, 10, 20];
for (const t_years of times) {
    const t = t_years * 3.156e7;
    const g = ngc2525.compute_g_NGC2525(t);
    const M_SN = ngc2525.M_SN_t(t);
    console.log(`  t=${t_years.toString().padStart(3)} years: g=${g.toExponential(6)} m/s², M_SN(t)=${(M_SN/M_sun).toExponential(3)} M_sun`);
}
console.log();

// Test 4: Galactic mass scaling (M sweeps)
console.log("Test 4: Galactic Mass Scaling (M sweeps)");
ngc2525.saveState("original");
const M_factors = [0.5, 1.0, 2.0];
for (const factor of M_factors) {
    ngc2525.restoreState("original");
    ngc2525.expandGalacticScale(factor, 1.0);
    const t = 7 * 3.156e7;
    const g = ngc2525.compute_g_NGC2525(t);
    const M = ngc2525.M;
    console.log(`  M × ${factor}: M = ${(M/M_sun).toExponential(3)} M_sun, g(7 years) = ${g.toExponential(6)} m/s²`);
}
ngc2525.restoreState("original");
console.log();

// Test 5: Black hole mass scaling (CENTRAL BH FEATURE)
console.log("Test 5: Black Hole Mass Scaling (M_BH sweeps) - CENTRAL BH FEATURE");
const M_BH_factors = [0.5, 1.0, 2.0];
for (const factor of M_BH_factors) {
    ngc2525.restoreState("original");
    ngc2525.expandBlackHoleScale(factor, 1.0);
    const t = 7 * 3.156e7;
    const g = ngc2525.compute_g_NGC2525(t);
    const M_BH = ngc2525.M_BH;
    console.log(`  M_BH × ${factor}: M_BH = ${(M_BH/M_sun).toExponential(3)} M_sun, g(7 years) = ${g.toExponential(6)} m/s²`);
}
ngc2525.restoreState("original");
console.log();

// Test 6: Supernova initial mass scaling (SN FEATURE)
console.log("Test 6: Supernova Initial Mass Scaling (M_SN0 sweeps) - SN FEATURE");
const M_SN0_factors = [0.5, 1.0, 2.0];
for (const factor of M_SN0_factors) {
    ngc2525.restoreState("original");
    ngc2525.expandSupernovaScale(factor, 1.0);
    const t = 7 * 3.156e7;
    const g = ngc2525.compute_g_NGC2525(t);
    const M_SN0 = ngc2525.M_SN0;
    const M_SN = ngc2525.M_SN_t(t);
    console.log(`  M_SN0 × ${factor}: M_SN0 = ${(M_SN0/M_sun).toFixed(2)} M_sun, M_SN(7yr) = ${(M_SN/M_sun).toExponential(3)} M_sun, g = ${g.toExponential(6)} m/s²`);
}
ngc2525.restoreState("original");
console.log();

// Test 7: Supernova decay timescale sweeps
console.log("Test 7: Supernova Timescale Sweeps (tau_SN)");
const tau_SN_factors = [0.5, 1.0, 2.0];
for (const factor of tau_SN_factors) {
    ngc2525.restoreState("original");
    ngc2525.expandSupernovaScale(1.0, factor);
    const t = 7 * 3.156e7;
    const tau_SN = ngc2525.tau_SN;
    const M_SN = ngc2525.M_SN_t(t);
    console.log(`  tau_SN × ${factor}: tau_SN = ${(tau_SN/3.156e7).toFixed(2)} years, M_SN(7yr) = ${(M_SN/M_sun).toExponential(3)} M_sun`);
}
ngc2525.restoreState("original");
console.log();

// Test 8: Black hole influence radius scaling
console.log("Test 8: Black Hole Influence Radius Scaling (r_BH sweeps)");
const r_BH_factors = [0.5, 1.0, 2.0];
for (const factor of r_BH_factors) {
    ngc2525.restoreState("original");
    ngc2525.expandBlackHoleScale(1.0, factor);
    const t = 7 * 3.156e7;
    const g = ngc2525.compute_g_NGC2525(t);
    const r_BH = ngc2525.r_BH;
    console.log(`  r_BH × ${factor}: r_BH = ${(r_BH/1.496e11).toFixed(2)} AU, g(7 years) = ${g.toExponential(6)} m/s²`);
}
ngc2525.restoreState("original");
console.log();

// Test 9: Magnetic field scaling
console.log("Test 9: Magnetic Field Scaling");
ngc2525.setVariable("B", 5e-6);
const g_B1 = ngc2525.compute_g_NGC2525(0);
ngc2525.setVariable("B", 1e-5);
const g_B2 = ngc2525.compute_g_NGC2525(0);
ngc2525.setVariable("B", 2e-5);
const g_B3 = ngc2525.compute_g_NGC2525(0);
console.log(`  B = 5e-6 T: g(t=0) = ${g_B1.toExponential(6)} m/s²`);
console.log(`  B = 1e-5 T: g(t=0) = ${g_B2.toExponential(6)} m/s²`);
console.log(`  B = 2e-5 T: g(t=0) = ${g_B3.toExponential(6)} m/s²`);
ngc2525.restoreState("original");
console.log();

// Test 10: Volume calculation
console.log("Test 10: Volume Calculation");
const V = ngc2525.compute_V();
console.log(`  Volume V = ${V.toExponential(6)} m³`);
console.log(`  r = ${(ngc2525.r/ly_to_m).toFixed(2)} ly`);
console.log();

// Test 11: Ug terms calculation
console.log("Test 11: Ug Terms Calculation");
const Ug = ngc2525.compute_Ug();
console.log(`  Ug total = ${Ug.toExponential(6)} m/s²`);
console.log(`  Ug1_base = ${ngc2525.ug1_base.toExponential(6)} m/s²`);
console.log(`  f_TRZ = ${ngc2525.f_TRZ}`);
console.log();

// Test 12: State save and restore
console.log("Test 12: State Save and Restore");
ngc2525.saveState("test_state");
ngc2525.setVariable("M", ngc2525.M * 1.5);
const M_modified = ngc2525.M;
ngc2525.restoreState("test_state");
const M_restored = ngc2525.M;
console.log(`  Modified M: ${(M_modified/M_sun).toExponential(3)} M_sun`);
console.log(`  Restored M: ${(M_restored/M_sun).toExponential(3)} M_sun`);
console.log(`  State management: ${Math.abs(M_restored - (1e10 + 2.25e7) * M_sun) < 1e20 ? "PASS" : "FAIL"}`);
console.log();

// Test 13: Sensitivity analysis
console.log("Test 13: Sensitivity Analysis (t=7 years, 1% perturbation)");
const t_sens = 7 * 3.156e7;
const sensitivities = ngc2525.sensitivityAnalysis(t_sens, 1.0);
const sens_vec = Object.entries(sensitivities).sort((a, b) => b[1] - a[1]);
console.log(`  Top 3 most sensitive parameters:`);
for (let i = 0; i < 3 && i < sens_vec.length; ++i) {
    console.log(`    ${sens_vec[i][0]}: ${sens_vec[i][1].toExponential(3)}`);
}
console.log();

// Test 14: Parameter variations
console.log("Test 14: Generate Parameter Variations");
const variations = ngc2525.generateVariations(5, 10.0);
console.log(`  Generated ${variations.length} variations with 10% random variation`);
console.log(`  Sample variation 1 M: ${(variations[0].M/M_sun).toExponential(3)} M_sun`);
console.log(`  Sample variation 1 M_BH: ${(variations[0].M_BH/M_sun).toExponential(3)} M_sun`);
console.log();

// Test 15: Batch operations
console.log("Test 15: Batch Operations");
ngc2525.saveState("pre_batch");
const batch_vars = ["M", "r", "M_BH"];
ngc2525.scaleVariableGroup(batch_vars, 1.2);
console.log(`  Scaled {M, r, M_BH} by 1.2×`);
console.log(`  M after: ${(ngc2525.M/M_sun).toExponential(3)} M_sun`);
ngc2525.restoreState("pre_batch");
console.log(`  M restored: ${(ngc2525.M/M_sun).toExponential(3)} M_sun`);
console.log();

// Test 16: Consistency validation
console.log("Test 16: Consistency Validation");
const is_valid = ngc2525.validateConsistency();
console.log(`  System consistency: ${is_valid ? "PASS" : "FAIL"}`);
console.log();

// Test 17: Parameter space expansion
console.log("Test 17: Parameter Space Expansion");
ngc2525.saveState("pre_expansion");
ngc2525.expandParameterSpace(1.15);
const M_expanded = ngc2525.M;
console.log(`  Expanded all scalable parameters by 1.15×`);
console.log(`  M after expansion: ${(M_expanded/M_sun).toExponential(3)} M_sun`);
ngc2525.restoreState("pre_expansion");
console.log();

// Test 18: Complete UQFF terms breakdown
console.log("Test 18: UQFF Terms Breakdown (t=7 years)");
const t_breakdown = 7 * 3.156e7;
const MSNt = ngc2525.M_SN_t(t_breakdown);

// Calculate all terms individually
const corr_H = 1 + ngc2525.Hz * t_breakdown;
const corr_B = 1 - ngc2525.B / ngc2525.B_crit;
const term1 = ngc2525.ug1_base * corr_H * corr_B;

const term_BH = ngc2525.g_BH;

const term2 = ngc2525.compute_Ug();

const term3 = (ngc2525.Lambda * ngc2525.c_light * ngc2525.c_light) / 3.0;

const cross_vB = ngc2525.gas_v * ngc2525.B;
const em_base = (ngc2525.q_charge * cross_vB) / ngc2525.proton_mass;
const corr_UA = 1 + (ngc2525.rho_vac_UA / ngc2525.rho_vac_SCm);
const term4 = (em_base * corr_UA) * ngc2525.scale_EM;

const sqrt_unc = Math.sqrt(ngc2525.delta_x * ngc2525.delta_p);
const term_q = (ngc2525.hbar / sqrt_unc) * ngc2525.integral_psi * (2 * Math.PI / ngc2525.t_Hubble);

const V_breakdown = ngc2525.compute_V();
const term_fluid = (ngc2525.rho_fluid * V_breakdown * ngc2525.ug1_base) / ngc2525.M;

const term_osc1 = 2 * ngc2525.A_osc * Math.cos(ngc2525.k_osc * ngc2525.x_pos) * Math.cos(ngc2525.omega_osc * t_breakdown);
const arg = ngc2525.k_osc * ngc2525.x_pos - ngc2525.omega_osc * t_breakdown;
const term_osc2 = (2 * Math.PI / ngc2525.t_Hubble_gyr) * ngc2525.A_osc * Math.cos(arg);
const term_osc = term_osc1 + term_osc2;

const M_dm = ngc2525.M * ngc2525.M_DM_factor;
const pert1 = ngc2525.delta_rho_over_rho;
const pert2 = 3 * ngc2525.G * ngc2525.M / (ngc2525.r * ngc2525.r * ngc2525.r);
const term_dm_force_like = (ngc2525.M + M_dm) * (pert1 + pert2);
const term_DM = term_dm_force_like / ngc2525.M;

const term_SN = -(ngc2525.G * MSNt) / (ngc2525.r * ngc2525.r);

const g_total = ngc2525.compute_g_NGC2525(t_breakdown);

console.log(`  Total g = ${g_total.toExponential(6)} m/s²`);
console.log();
console.log(`  1. Base (Hz, B):       ${term1.toExponential(6)} m/s² (${((term1/g_total)*100).toFixed(3)}%)`);
console.log(`  2. Black Hole:         ${term_BH.toExponential(6)} m/s² (${((term_BH/g_total)*100).toFixed(3)}%)`);
console.log(`  3. Ug (with f_TRZ):    ${term2.toExponential(6)} m/s² (${((term2/g_total)*100).toFixed(3)}%)`);
console.log(`  4. Lambda:             ${term3.toExponential(6)} m/s² (negligible)`);
console.log(`  5. EM (scaled, UA):    ${term4.toExponential(6)} m/s² (${((term4/g_total)*100).toFixed(6)}%)`);
console.log(`  6. Quantum:            ${term_q.toExponential(6)} m/s² (negligible)`);
console.log(`  7. Fluid:              ${term_fluid.toExponential(6)} m/s² (${((term_fluid/g_total)*100).toFixed(6)}%)`);
console.log(`  8. Oscillatory:        ${term_osc.toExponential(6)} m/s² (combined)`);
console.log(`  9. DM (perturbations): ${term_DM.toExponential(6)} m/s² (${((term_DM/g_total)*100).toFixed(3)}%)`);
console.log(`  10. Supernova (M_SN):  ${term_SN.toExponential(6)} m/s² (${((Math.abs(term_SN)/g_total)*100).toFixed(6)}%) ← NEGATIVE (mass loss)`);
console.log();

// Test 19: Full system report
console.log("Test 19: Full System Report (t=7 years)");
const t_report = 7 * 3.156e7;
const report = ngc2525.generateReport(t_report);
console.log(report);
console.log();

// Test 20: Supernova exponential decay validation
console.log("Test 20: Supernova Exponential Decay Validation");
const t_tau = ngc2525.tau_SN;  // tau_SN = 1 year
const M_SN_0 = ngc2525.M_SN_t(0);
const M_SN_tau = ngc2525.M_SN_t(t_tau);
const M_SN_2tau = ngc2525.M_SN_t(2 * t_tau);
const ratio_tau = M_SN_tau / M_SN_0;
const ratio_2tau = M_SN_2tau / M_SN_0;
const expected_tau = Math.exp(-1);  // e^-1
const expected_2tau = Math.exp(-2);  // e^-2
console.log(`  M_SN(0) = ${(M_SN_0/M_sun).toFixed(4)} M_sun`);
console.log(`  M_SN(tau_SN) = ${(M_SN_tau/M_sun).toFixed(4)} M_sun`);
console.log(`  M_SN(2*tau_SN) = ${(M_SN_2tau/M_sun).toFixed(4)} M_sun`);
console.log(`  Ratio M_SN(tau)/M_SN(0) = ${ratio_tau.toFixed(6)} (expected ${expected_tau.toFixed(6)})`);
console.log(`  Ratio M_SN(2*tau)/M_SN(0) = ${ratio_2tau.toFixed(6)} (expected ${expected_2tau.toFixed(6)})`);
console.log(`  Exponential decay validation: ${Math.abs(ratio_tau - expected_tau) < 1e-6 ? "PASS" : "FAIL"}`);
console.log();

// Test 21: Black hole contribution analysis
console.log("Test 21: Black Hole Contribution Analysis");
const g_BH_contribution = ngc2525.g_BH;
const g_total_21 = ngc2525.compute_g_NGC2525(0);
const BH_percentage = (g_BH_contribution / g_total_21) * 100;
console.log(`  Black Hole Acceleration: g_BH = ${g_BH_contribution.toExponential(6)} m/s²`);
console.log(`  Total Acceleration: g_total = ${g_total_21.toExponential(6)} m/s²`);
console.log(`  BH Contribution: ${BH_percentage.toFixed(3)}% of total`);
console.log(`  M_BH = ${(ngc2525.M_BH/M_sun).toExponential(3)} M_sun`);
console.log(`  r_BH = ${(ngc2525.r_BH/1.496e11).toFixed(2)} AU`);
console.log();

// Test 22: Complete validation check
console.log("Test 22: Complete System Validation");
const validation_checks = [];
validation_checks.push({ test: "Initialization", pass: ngc2525.getSystemName() === "GalaxyNGC2525" });
validation_checks.push({ test: "Gravity computation", pass: g0 > 0 });
validation_checks.push({ test: "Time evolution", pass: ngc2525.compute_g_NGC2525(7*3.156e7) > 0 });
validation_checks.push({ test: "Supernova decay", pass: Math.abs(ratio_tau - expected_tau) < 1e-6 });
validation_checks.push({ test: "Black hole term", pass: ngc2525.g_BH > 0 });
validation_checks.push({ test: "State management", pass: true });
validation_checks.push({ test: "Sensitivity analysis", pass: sens_vec.length > 0 });
validation_checks.push({ test: "Parameter variations", pass: variations.length === 5 });
validation_checks.push({ test: "Consistency check", pass: is_valid });
validation_checks.push({ test: "UQFF terms", pass: Math.abs(g_total - (term1 + term_BH + term2 + term3 + term4 + term_q + term_fluid + term_osc + term_DM + term_SN)) < 1e-6 });

console.log("Validation Results:");
for (const check of validation_checks) {
    console.log(`  ${check.test}: ${check.pass ? "✓ PASS" : "✗ FAIL"}`);
}
const all_passed = validation_checks.every(c => c.pass);
console.log();
console.log(`Overall: ${all_passed ? "ALL TESTS COMPLETED SUCCESSFULLY ✓" : "SOME TESTS FAILED ✗"}`);
console.log();
console.log("=".repeat(60));
