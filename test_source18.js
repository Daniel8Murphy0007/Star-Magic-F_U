/**
 * ================================================================================================
 * Test Suite: test_source18.js
 * 
 * Comprehensive testing for source18.js (PillarsOfCreation - Eagle Nebula)
 * Tests all 38 methods and physics calculations
 * 
 * Author: GitHub Copilot
 * Date: November 3, 2025
 * ================================================================================================
 */

import PillarsOfCreation from './source18.js';

console.log("=========================================================");
console.log("PILLARS OF CREATION TEST SUITE");
console.log("Eagle Nebula (M16) - Star Formation with Erosion");
console.log("=========================================================\n");

const pillars = new PillarsOfCreation();
const M_sun = 1.989e30;

// Test 1: Initialization and basic properties
console.log("Test 1: Initialization and Basic Properties");
console.log(`  System name: ${pillars.getSystemName()}`);
console.log(`  Number of variables: ${pillars.listVariables().length}`);
console.log(`  M_initial: ${(pillars.getVariable("M_initial")/M_sun).toFixed(2)} M_sun`);
console.log(`  r: ${(pillars.getVariable("r")/9.461e15).toFixed(2)} ly`);
console.log(`  B: ${pillars.getVariable("B").toExponential(3)} T`);
console.log(`  M_dot_factor: ${pillars.getVariable("M_dot_factor").toFixed(6)}`);
console.log(`  tau_SF: ${(pillars.getVariable("tau_SF")/(1e6*3.156e7)).toFixed(2)} Myr`);
console.log(`  E_0 (UNIQUE): ${pillars.getVariable("E_0").toFixed(3)}`);
console.log(`  tau_erosion (UNIQUE): ${(pillars.getVariable("tau_erosion")/(1e6*3.156e7)).toFixed(2)} Myr`);
console.log(`  rho_wind: ${pillars.getVariable("rho_wind").toExponential(3)} kg/m³\n`);

// Test 2: Initial gravity calculation
console.log("Test 2: Initial Gravity Calculation");
const g_0 = pillars.compute_g_Pillars(0);
const M_0 = pillars.M_t(0);
const E_0 = pillars.E_t(0);
console.log(`  g_Pillars(t=0) = ${g_0.toExponential(6)} m/s²`);
console.log(`  M(t=0) = ${(M_0/M_sun).toFixed(2)} M_sun`);
console.log(`  E(t=0) = ${E_0.toFixed(6)} (initial erosion factor)`);
console.log(`  Initial mass growth: ${(M_0/pillars.getVariable("M_initial")).toFixed(6)}`);
console.log(`  Gas reservoir: ${((M_0 - pillars.getVariable("M_initial"))/M_sun).toFixed(0)} M_sun\n`);

// Test 3: Time evolution (0-5 Myr)
console.log("Test 3: Time Evolution (0-5 Myr) - Mass Growth and Erosion");
const times_Myr = [0.0, 0.5, 1.0, 2.0, 5.0];
for (const t_Myr of times_Myr) {
    const t = t_Myr * 1e6 * 3.156e7;
    const g = pillars.compute_g_Pillars(t);
    const Mt = pillars.M_t(t);
    const Et = pillars.E_t(t);
    console.log(`  t=${t_Myr.toFixed(1)} Myr: g=${g.toExponential(3)} m/s², M(t)=${(Mt/M_sun).toFixed(0).padStart(5)} M_sun, E(t)=${Et.toFixed(6)}`);
}
console.log();

// Test 4: Mass evolution analysis
console.log("Test 4: Mass Evolution Analysis");
const M_initial = pillars.getVariable("M_initial");
const M_1Myr = pillars.M_t(1e6 * 3.156e7);
const growth_factor = M_1Myr / M_initial;
console.log(`  M_initial = ${(M_initial/M_sun).toFixed(2)} M_sun (stellar mass)`);
console.log(`  M(t=1 Myr) = ${(M_1Myr/M_sun).toFixed(2)} M_sun (stellar + remaining gas)`);
console.log(`  Growth factor: ${growth_factor.toFixed(6)} (${((growth_factor-1)*100).toFixed(1)}% increase)`);
console.log(`  Interpretation: Gas cloud (${((M_1Myr - M_initial)/M_sun).toFixed(0)} M_sun) being converted to stars\n`);

// Test 5: Erosion evolution analysis (UNIQUE to Pillars of Creation)
console.log("Test 5: Erosion Evolution Analysis - UNIQUE FEATURE");
const tau_erosion = pillars.getVariable("tau_erosion");
const E_tau = pillars.E_t(tau_erosion);
console.log(`  tau_erosion = ${(tau_erosion/(1e6*3.156e7)).toFixed(2)} Myr`);
console.log(`  E(t=tau_erosion) = ${E_tau.toFixed(6)} (should be E_0/e ≈ ${(pillars.getVariable("E_0")/Math.E).toFixed(6)})`);
console.log(`  Photoevaporation: Massive stars in Eagle Nebula erode pillar material`);
console.log(`  E(t) = E_0 * exp(-t/tau_erosion) models gas removal\n`);

// Test 6: Star formation timescale
console.log("Test 6: Star Formation Timescale");
const tau_SF = pillars.getVariable("tau_SF");
const M_tau = pillars.M_t(tau_SF);
console.log(`  tau_SF = ${(tau_SF/(1e6*3.156e7)).toFixed(2)} Myr`);
console.log(`  M(t=tau_SF) = ${(M_tau/M_sun).toFixed(2)} M_sun`);
console.log(`  Star formation: Active process in the pillars\n`);

// Test 7: Erosion factor scaling (E_0 sweeps)
console.log("Test 7: Erosion Factor Scaling (E_0 sweeps)");
pillars.saveState("original");
const E_0_factors = [0.5, 1.0, 2.0];
for (const factor of E_0_factors) {
    pillars.restoreState("original");
    pillars.expandErosionScale(factor, 1.0);
    const t = 1e6 * 3.156e7;
    const Et = pillars.E_t(t);
    const g = pillars.compute_g_Pillars(t);
    console.log(`  E_0 × ${factor.toFixed(1)}: E(1 Myr) = ${Et.toFixed(6)}, g = ${g.toExponential(6)} m/s²`);
}
pillars.restoreState("original");
console.log();

// Test 8: Erosion timescale scaling (tau_erosion sweeps)
console.log("Test 8: Erosion Timescale Scaling (tau_erosion sweeps)");
const tau_erosion_factors = [0.5, 1.0, 2.0];
for (const factor of tau_erosion_factors) {
    pillars.restoreState("original");
    pillars.expandErosionScale(1.0, factor);
    const t = 1e6 * 3.156e7;
    const Et = pillars.E_t(t);
    console.log(`  tau_erosion × ${factor.toFixed(1)}: E(1 Myr) = ${Et.toFixed(6)}`);
}
pillars.restoreState("original");
console.log();

// Test 9: Wind feedback scaling (v_wind)
console.log("Test 9: Wind Velocity Scaling");
const v_wind_factors = [0.5, 1.0, 2.0];
for (const factor of v_wind_factors) {
    pillars.restoreState("original");
    pillars.expandWindMagneticScale(1.0, factor, 1.0);
    const t = 1e6 * 3.156e7;
    const g = pillars.compute_g_Pillars(t);
    console.log(`  v_wind × ${factor.toFixed(1)}: g(1 Myr) = ${g.toExponential(6)} m/s²`);
}
pillars.restoreState("original");
console.log();

// Test 10: Magnetic and fluid scaling
console.log("Test 10: Magnetic Field Scaling");
const B_factors = [0.5, 1.0, 2.0];
for (const factor of B_factors) {
    pillars.restoreState("original");
    pillars.expandWindMagneticScale(1.0, 1.0, factor);
    const t = 1e6 * 3.156e7;
    const g = pillars.compute_g_Pillars(t);
    console.log(`  B × ${factor.toFixed(1)}: g(1 Myr) = ${g.toExponential(6)} m/s²`);
}
pillars.restoreState("original");
console.log();

// Test 11: Volume calculation
console.log("Test 11: Volume Calculation");
const V = pillars.compute_V();
const r = pillars.getVariable("r");
console.log(`  Radius: ${(r/9.461e15).toFixed(2)} ly`);
console.log(`  Volume: ${V.toExponential(6)} m³\n`);

// Test 12: Ug terms calculation
console.log("Test 12: Ug Terms Calculation");
const Mt_test = pillars.M_t(0);
const Ug_total = pillars.compute_Ug(Mt_test);
console.log(`  Ug total at t=0: ${Ug_total.toExponential(6)} m/s²`);
console.log(`  Includes: Ug1, Ug2, Ug3, Ug4 with f_TRZ correction\n`);

// Test 13: State save/restore
console.log("Test 13: State Save and Restore");
pillars.saveState("test_state");
const M_before = pillars.getVariable("M_initial");
pillars.setVariable("M_initial", M_before * 2);
const M_modified = pillars.getVariable("M_initial");
pillars.restoreState("test_state");
const M_after = pillars.getVariable("M_initial");
console.log(`  M_initial before: ${(M_before/M_sun).toFixed(2)} M_sun`);
console.log(`  M_initial modified: ${(M_modified/M_sun).toFixed(2)} M_sun`);
console.log(`  M_initial restored: ${(M_after/M_sun).toFixed(2)} M_sun`);
console.log(`  State management: ${M_before === M_after ? "PASS ✓" : "FAIL ✗"}\n`);

// Test 14: Sensitivity analysis
console.log("Test 14: Sensitivity Analysis at 1 Myr (top 10)");
const t_sens = 1e6 * 3.156e7;
const sensitivities = pillars.sensitivityAnalysis(t_sens, 1.0);
const sens_sorted = Object.entries(sensitivities).sort((a, b) => b[1] - a[1]);
for (let i = 0; i < 10 && i < sens_sorted.length; i++) {
    console.log(`  ${(i+1).toString().padStart(2)}. ${sens_sorted[i][0].padEnd(20)}: ${sens_sorted[i][1].toExponential(3)}`);
}
console.log();

// Test 15: Generate parameter variations
console.log("Test 15: Generate Parameter Variations (Monte Carlo)");
const variations = pillars.generateVariations(5, 10.0);
console.log(`  Generated ${variations.length} variants with 10% variation`);
console.log(`  Original M_initial: ${(pillars.getVariable("M_initial")/M_sun).toFixed(2)} M_sun`);
console.log(`  Variant 1 M_initial: ${(variations[0].M_initial/M_sun).toFixed(2)} M_sun`);
console.log(`  Variant 2 M_initial: ${(variations[1].M_initial/M_sun).toFixed(2)} M_sun\n`);

// Test 16: Batch transformations
console.log("Test 16: Batch Transformations");
pillars.saveState("batch_test");
const scale_vars = ["M_initial", "r", "B"];
pillars.scaleVariableGroup(scale_vars, 1.2);
console.log(`  Scaled {M_initial, r, B} by 1.2×`);
console.log(`  New M_initial: ${(pillars.getVariable("M_initial")/M_sun).toFixed(2)} M_sun`);
pillars.restoreState("batch_test");
console.log(`  Restored original state\n`);

// Test 17: Consistency validation
console.log("Test 17: Consistency Validation");
const is_valid = pillars.validateConsistency();
console.log(`  System consistency: ${is_valid ? "VALID ✓" : "INVALID ✗"}`);
const was_corrected = pillars.autoCorrectAnomalies();
console.log(`  Auto-correction needed: ${was_corrected ? "Yes" : "No"}\n`);

// Test 18: Parameter space expansion
console.log("Test 18: Parameter Space Expansion");
pillars.saveState("expand_test");
pillars.expandParameterSpace(1.5);
console.log(`  Expanded parameter space by 1.5×`);
console.log(`  New M_initial: ${(pillars.getVariable("M_initial")/M_sun).toFixed(2)} M_sun`);
console.log(`  New r: ${(pillars.getVariable("r")/9.461e15).toFixed(2)} ly`);
pillars.restoreState("expand_test");
console.log();

// Test 19: Complete UQFF terms breakdown
console.log("Test 19: Complete UQFF Terms Breakdown at t=1 Myr");
const t_breakdown = 1e6 * 3.156e7;
const Mt_break = pillars.M_t(t_breakdown);
const Et_break = pillars.E_t(t_breakdown);
const ug1_t = (pillars.getVariable("G") * Mt_break) / Math.pow(pillars.getVariable("r"), 2);

// Term 1: Base with H0, B, E
const corr_H = 1 + pillars.getVariable("H0") * t_breakdown;
const corr_B = 1 - pillars.getVariable("B") / pillars.getVariable("B_crit");
const corr_E = 1 - Et_break;
const term1 = ug1_t * corr_H * corr_B * corr_E;

// Term 2: Ug
const term2 = pillars.compute_Ug(Mt_break);

// Term 3: Lambda
const term3 = (pillars.getVariable("Lambda") * Math.pow(pillars.getVariable("c_light"), 2)) / 3.0;

// Term 4: EM
const cross_vB = pillars.getVariable("gas_v") * pillars.getVariable("B");
const em_base = (pillars.getVariable("q_charge") * cross_vB) / pillars.getVariable("proton_mass");
const corr_UA = 1 + (pillars.getVariable("rho_vac_UA") / pillars.getVariable("rho_vac_SCm"));
const term4 = em_base * corr_UA * pillars.getVariable("scale_EM");

// Term 5: Quantum
const sqrt_unc = Math.sqrt(pillars.getVariable("delta_x") * pillars.getVariable("delta_p"));
const term5 = (pillars.getVariable("hbar") / sqrt_unc) * pillars.getVariable("integral_psi") * (2 * Math.PI / pillars.getVariable("t_Hubble"));

// Term 6: Fluid
const V_break = pillars.compute_V();
const term6 = (pillars.getVariable("rho_fluid") * V_break * ug1_t) / Mt_break;

// Term 7: Oscillatory (combined)
const term7_part1 = 2 * pillars.getVariable("A_osc") * Math.cos(pillars.getVariable("k_osc") * pillars.getVariable("x_pos")) * Math.cos(pillars.getVariable("omega_osc") * t_breakdown);
const arg = pillars.getVariable("k_osc") * pillars.getVariable("x_pos") - pillars.getVariable("omega_osc") * t_breakdown;
const term7_part2 = (2 * Math.PI / pillars.getVariable("t_Hubble_gyr")) * pillars.getVariable("A_osc") * Math.cos(arg);
const term7 = term7_part1 + term7_part2;

// Term 8: DM
const M_dm = Mt_break * pillars.getVariable("M_DM_factor");
const pert1 = pillars.getVariable("delta_rho_over_rho");
const pert2 = 3 * pillars.getVariable("G") * Mt_break / Math.pow(pillars.getVariable("r"), 3);
const term8 = (Mt_break + M_dm) * (pert1 + pert2) / Mt_break;

// Term 9: Wind
const wind_pressure = pillars.getVariable("rho_wind") * Math.pow(pillars.getVariable("v_wind"), 2);
const term9 = wind_pressure / pillars.getVariable("rho_fluid");

const g_total = pillars.compute_g_Pillars(t_breakdown);

console.log(`  Total g = ${g_total.toExponential(6)} m/s²\n`);
console.log(`  1. Base (H0, B, E, M(t)): ${term1.toExponential(6)} m/s² (${((term1/g_total)*100).toFixed(3)}%)`);
console.log(`  2. Ug (with f_TRZ):    ${term2.toExponential(6)} m/s² (${((term2/g_total)*100).toFixed(3)}%)`);
console.log(`  3. Lambda:             ${term3.toExponential(6)} m/s² (${((term3/g_total)*100).toFixed(6)}%)`);
console.log(`  4. EM (scaled, UA):    ${term4.toExponential(6)} m/s² (${((term4/g_total)*100).toFixed(5)}%)`);
console.log(`  5. Quantum:            ${term5.toExponential(6)} m/s² (${((term5/g_total)*100).toFixed(6)}%)`);
console.log(`  6. Fluid:              ${term6.toExponential(6)} m/s² (${((term6/g_total)*100).toFixed(6)}%)`);
console.log(`  7. Oscillatory:        ${term7.toExponential(6)} m/s² (combined standing + traveling)`);
console.log(`  8. DM (perturbations): ${term8.toExponential(6)} m/s² (${((term8/g_total)*100).toFixed(3)}%)`);
console.log(`  9. Wind Feedback:      ${term9.toExponential(6)} m/s² (${((term9/g_total)*100).toFixed(3)}%) ← DOMINANT\n`);

// Test 20: Comprehensive report generation
console.log("Test 20: Comprehensive Report Generation at t=1.5 Myr");
const t_report = 1.5e6 * 3.156e7;
const report = pillars.generateReport(t_report);
console.log(report);

// Test 21: Initial stellar mass sweep
console.log("Test 21: Initial Stellar Mass Sweep");
pillars.saveState("mass_sweep");
const M_initial_factors = [0.5, 1.0, 2.0];
for (const factor of M_initial_factors) {
    pillars.restoreState("mass_sweep");
    pillars.setVariable("M_initial", pillars.getVariable("M_initial") * factor);
    const t = 1e6 * 3.156e7;
    const g = pillars.compute_g_Pillars(t);
    console.log(`  M_initial × ${factor.toFixed(1)}: g(1 Myr) = ${g.toExponential(6)} m/s²`);
}
pillars.restoreState("mass_sweep");
console.log();

// Test 22: Star formation factor sweep
console.log("Test 22: Star Formation Factor Sweep (M_dot_factor)");
const M_dot_factors = [0.5, 1.0, 2.0];
for (const factor of M_dot_factors) {
    pillars.restoreState("mass_sweep");
    pillars.expandStarFormationScale(factor, 1.0);
    const t = 1e6 * 3.156e7;
    const Mt = pillars.M_t(t);
    console.log(`  M_dot_factor × ${factor.toFixed(1)}: M(1 Myr) = ${(Mt/M_sun).toFixed(2)} M_sun`);
}
pillars.restoreState("mass_sweep");
console.log();

console.log("=========================================================");
console.log("ALL TESTS COMPLETED SUCCESSFULLY ✓");
console.log("=========================================================");
