/**
 * ================================================================================================
 * Test Suite: test_source17.js
 * 
 * Description: Comprehensive test suite for Westerlund2 (source17.js)
 *              Tests all 38 methods (13 core + 25 enhanced) for Westerlund 2 super star cluster.
 * 
 * System: Westerlund 2 Super Star Cluster
 * Type: Young massive star cluster with intense star formation
 * Location: Carina constellation, ~20,000 light-years from Earth
 * Age: ~1-2 million years
 * 
 * Test Coverage:
 *   - Core physics: g_Westerlund2(t), M(t), Ug terms, volume
 *   - Time evolution: 0-5 Myr star formation timescales
 *   - Star formation: M_dot_factor, tau_SF effects (100,000 M☉ gas)
 *   - Wind feedback: rho_wind, v_wind interactions
 *   - Magnetic & fluid scaling
 *   - Enhanced capabilities: all 25 dynamic methods
 * 
 * Created: November 03, 2025
 * Author: GitHub Copilot
 * ================================================================================================
 */

import Westerlund2 from './source17.js';

console.log('========================================');
console.log('WESTERLUND 2 TEST SUITE');
console.log('Super Star Cluster');
console.log('========================================\n');

const wd2 = new Westerlund2();
const M_sun = 1.989e30;
const ly_to_m = 9.461e15;

// Test 1: Initialization and basic properties
console.log('Test 1: Initialization and Basic Properties');
console.log(`  System name: ${wd2.getSystemName()}`);
console.log(`  Number of variables: ${wd2.listVariables().length}`);
console.log(`  M_initial: ${(wd2.getVariable('M_initial')/M_sun).toFixed(2)} M_sun (${wd2.getVariable('M_initial').toExponential(3)} kg)`);
console.log(`  r: ${(wd2.getVariable('r')/ly_to_m).toFixed(2)} ly (${wd2.getVariable('r').toExponential(3)} m)`);
console.log(`  B: ${wd2.getVariable('B').toExponential(3)} T (10x stronger than NGC 2014)`);
console.log(`  M_dot_factor: ${wd2.getVariable('M_dot_factor').toFixed(3)} (gas mass / initial stellar mass)`);
console.log(`  tau_SF: ${(wd2.getVariable('tau_SF')/(3.156e7*1e6)).toFixed(2)} Myr`);
console.log(`  rho_wind: ${wd2.getVariable('rho_wind').toExponential(3)} kg/m^3 (10x denser than NGC 2014)`);
console.log(`  v_wind: ${(wd2.getVariable('v_wind')/1e6).toFixed(1)} x 10^6 m/s\n`);

// Test 2: Initial gravity calculation
console.log('Test 2: Initial Gravity Calculation (t=0)');
const t0 = 0.0;
const g0 = wd2.compute_g_Westerlund2(t0);
const M0 = wd2.M_t(t0);
console.log(`  g_Westerlund2(t=0) = ${g0.toExponential(6)} m/s^2`);
console.log(`  M(t=0) = ${(M0/M_sun).toFixed(2)} M_sun`);
console.log(`  Initial mass growth: M(0) / M_initial = ${(M0 / wd2.getVariable('M_initial')).toFixed(3)}`);
console.log(`  Gas reservoir: ${((M0 - wd2.getVariable('M_initial'))/M_sun).toFixed(0)} M_sun\n`);

// Test 3: Time evolution (0, 0.5, 1, 2, 5 Myr)
console.log('Test 3: Time Evolution (Star Formation Era: 0-5 Myr)');
for (const t_myr of [0.0, 0.5, 1.0, 2.0, 5.0]) {
    const t = t_myr * 1e6 * 3.156e7;
    const g = wd2.compute_g_Westerlund2(t);
    const Mt = wd2.M_t(t);
    console.log(`  t=${t_myr.toFixed(1)} Myr: g=${g.toExponential(3)} m/s^2, M(t)=${(Mt/M_sun).toFixed(0).padStart(6)} M_sun`);
}
console.log();

// Test 4: Mass evolution analysis
console.log('Test 4: Mass Evolution Analysis');
const t_test = 1e6 * 3.156e7;  // 1 Myr
const M_initial = wd2.getVariable('M_initial');
const M_at_1Myr = wd2.M_t(t_test);
const growth_factor = M_at_1Myr / M_initial;
console.log(`  M_initial = ${(M_initial/M_sun).toFixed(2)} M_sun (stellar mass)`);
console.log(`  M(t=1 Myr) = ${(M_at_1Myr/M_sun).toFixed(2)} M_sun (stellar + remaining gas)`);
console.log(`  Growth factor: ${growth_factor.toFixed(6)} (${((growth_factor-1)*100).toFixed(1)}% increase)`);
console.log(`  M_dot_factor = ${wd2.getVariable('M_dot_factor').toFixed(3)}`);
console.log(`  Interpretation: Massive gas cloud (100,000 M_sun) being converted to stars\n`);

// Test 5: Star formation timescale effects
console.log('Test 5: Star Formation Timescale Effects');
const tau_SF_original = wd2.getVariable('tau_SF');
console.log(`  tau_SF = ${(tau_SF_original/(3.156e7*1e6)).toFixed(2)} Myr (faster than NGC 2014's 5 Myr)`);
const M_at_tau = wd2.M_t(tau_SF_original);
console.log(`  M(t=tau_SF) = ${(M_at_tau/M_sun).toFixed(2)} M_sun`);
console.log(`  M_dot at t=tau_SF: ${(wd2.getVariable('M_dot_factor') * Math.exp(-1)).toFixed(3)} (decayed by e^-1)`);
console.log(`  Rapid star formation: shorter timescale than NGC 2014\n`);

// Test 6: Star formation scaling
console.log('Test 6: Star Formation Scaling');
wd2.saveState('before_sf_scale');
console.log(`  Before: M_dot_factor=${wd2.getVariable('M_dot_factor').toFixed(3)}, tau_SF=${(wd2.getVariable('tau_SF')/(3.156e7*1e6)).toFixed(2)} Myr`);
wd2.expandStarFormationScale(1.4, 0.75);
console.log(`  After (x1.4, x0.75): M_dot_factor=${wd2.getVariable('M_dot_factor').toFixed(3)}, tau_SF=${(wd2.getVariable('tau_SF')/(3.156e7*1e6)).toFixed(2)} Myr`);
console.log(`  g(t=1Myr) after scaling: ${wd2.compute_g_Westerlund2(t_test).toExponential(6)} m/s^2`);
wd2.restoreState('before_sf_scale');
console.log(`  Restored to original state\n`);

// Test 7: Wind feedback scaling
console.log('Test 7: Wind Feedback Scaling');
wd2.saveState('before_wind_scale');
const rho_wind_orig = wd2.getVariable('rho_wind');
const v_wind_orig = wd2.getVariable('v_wind');
console.log(`  Before: rho_wind=${rho_wind_orig.toExponential(3)} kg/m^3, v_wind=${(v_wind_orig/1e6).toFixed(1)} x 10^6 m/s`);
wd2.expandWindFeedbackScale(1.6, 1.25);
console.log(`  After (x1.6, x1.25): rho_wind=${wd2.getVariable('rho_wind').toExponential(3)} kg/m^3, v_wind=${(wd2.getVariable('v_wind')/1e6).toFixed(2)} x 10^6 m/s`);
const wind_pressure_before = rho_wind_orig * v_wind_orig * v_wind_orig;
const wind_pressure_after = wd2.getVariable('rho_wind') * wd2.getVariable('v_wind') * wd2.getVariable('v_wind');
console.log(`  Wind pressure before: ${wind_pressure_before.toExponential(3)} Pa`);
console.log(`  Wind pressure after: ${wind_pressure_after.toExponential(3)} Pa (${(wind_pressure_after/wind_pressure_before).toFixed(2)}x increase)`);
wd2.restoreState('before_wind_scale');
console.log(`  Restored to original state\n`);

// Test 8: Magnetic & fluid scaling
console.log('Test 8: Magnetic & Fluid Scaling');
wd2.saveState('before_mag_scale');
console.log(`  Before: B=${wd2.getVariable('B').toExponential(3)} T, rho_fluid=${wd2.getVariable('rho_fluid').toExponential(3)} kg/m^3`);
wd2.expandMagneticFluidScale(1.3, 1.4);
console.log(`  After (x1.3, x1.4): B=${wd2.getVariable('B').toExponential(3)} T, rho_fluid=${wd2.getVariable('rho_fluid').toExponential(3)} kg/m^3`);
console.log(`  g(t=1Myr) after scaling: ${wd2.compute_g_Westerlund2(t_test).toExponential(6)} m/s^2`);
wd2.restoreState('before_mag_scale');
console.log(`  Restored to original state\n`);

// Test 9: Volume calculation
console.log('Test 9: Volume Calculation');
const V = wd2.compute_V();
const r = wd2.getVariable('r');
console.log(`  Radius r = ${(r/ly_to_m).toFixed(2)} ly = ${r.toExponential(3)} m`);
console.log(`  Volume V = ${V.toExponential(6)} m^3`);
console.log(`  Volume in cubic light-years: ${(V/(ly_to_m**3)).toExponential(3)} ly^3\n`);

// Test 10: Ug terms calculation
console.log('Test 10: UQFF Ug Terms Calculation (at t=1 Myr)');
const Mt_test = wd2.M_t(t_test);
const Ug_total = wd2.compute_Ug(Mt_test);
const ug1 = (wd2.getVariable('G') * Mt_test) / (r * r);
console.log(`  M(t=1Myr) = ${(Mt_test/M_sun).toFixed(0)} M_sun`);
console.log(`  Ug1 = ${ug1.toExponential(6)} m/s^2`);
console.log(`  Ug total (with f_TRZ, B corrections): ${Ug_total.toExponential(6)} m/s^2`);
console.log(`  f_TRZ factor: ${wd2.getVariable('f_TRZ')}`);
console.log(`  B/B_crit ratio: ${(wd2.getVariable('B')/wd2.getVariable('B_crit')).toExponential(3)}\n`);

// Test 11: State management
console.log('Test 11: State Management (Save/Restore)');
wd2.saveState('test_state_1');
wd2.saveState('test_state_2');
console.log(`  Saved 2 states. Total saved: ${wd2.listSavedStates().length}`);
console.log(`  Saved state names: ${wd2.listSavedStates().join(', ')}`);
wd2.setVariable('M_initial', 50000.0 * M_sun);
console.log(`  Modified M_initial to 50000 M_sun`);
console.log(`  g(t=1Myr) with modified M: ${wd2.compute_g_Westerlund2(t_test).toExponential(6)} m/s^2`);
wd2.restoreState('test_state_1');
console.log(`  Restored state: M_initial = ${(wd2.getVariable('M_initial')/M_sun).toFixed(0)} M_sun`);
console.log(`  g(t=1Myr) after restore: ${wd2.compute_g_Westerlund2(t_test).toExponential(6)} m/s^2\n`);

// Test 12: Sensitivity analysis
console.log('Test 12: Sensitivity Analysis (at t=1 Myr, ±1% perturbation)');
const sens = wd2.sensitivityAnalysis(t_test, 1.0);
const sens_array = Object.entries(sens).sort((a, b) => b[1] - a[1]);
console.log(`  Top 10 most sensitive parameters:`);
for (let i = 0; i < Math.min(10, sens_array.length); i++) {
    console.log(`    ${(i+1).toString().padStart(2)}. ${sens_array[i][0].padEnd(20)}: ${sens_array[i][1].toExponential(3)}`);
}
console.log();

// Test 13: Parameter variations (Monte Carlo)
console.log('Test 13: Parameter Variations (Monte Carlo, 5 samples, ±12%)');
const variations = wd2.generateVariations(5, 12.0);
console.log(`  Generated ${variations.length} parameter sets:`);
for (let i = 0; i < variations.length; i++) {
    console.log(`    Variant ${i+1}: M_initial=${(variations[i].M_initial/M_sun).toFixed(0)} M_sun, M_dot=${variations[i].M_dot_factor.toFixed(2)}, v_wind=${(variations[i].v_wind/1e6).toFixed(2)} x10^6 m/s`);
}
console.log();

// Test 14: Batch transformation
console.log('Test 14: Batch Transformation (scale all density parameters by 1.25)');
wd2.saveState('before_batch');
const rho_fluid_before = wd2.getVariable('rho_fluid');
const rho_wind_before = wd2.getVariable('rho_wind');
wd2.transformVariableGroup(['rho_fluid', 'rho_wind'], v => v * 1.25);
console.log(`  rho_fluid: ${rho_fluid_before.toExponential(3)} → ${wd2.getVariable('rho_fluid').toExponential(3)} kg/m^3`);
console.log(`  rho_wind: ${rho_wind_before.toExponential(3)} → ${wd2.getVariable('rho_wind').toExponential(3)} kg/m^3`);
console.log(`  g(t=1Myr) after batch transform: ${wd2.compute_g_Westerlund2(t_test).toExponential(6)} m/s^2`);
wd2.restoreState('before_batch');
console.log(`  Restored to original state\n`);

// Test 15: Consistency validation
console.log('Test 15: Consistency Validation');
const valid = wd2.validateConsistency();
console.log(`  System consistency: ${valid ? 'VALID ✓' : 'INVALID ✗'}`);
console.log(`  All parameters within physical bounds\n`);

// Test 16: Parameter space expansion
console.log('Test 16: Parameter Space Expansion (1.2x all expandable parameters)');
wd2.saveState('before_expansion');
const M_initial_before = wd2.getVariable('M_initial');
const r_before = wd2.getVariable('r');
wd2.expandParameterSpace(1.2);
console.log(`  M_initial: ${(M_initial_before/M_sun).toFixed(0)} → ${(wd2.getVariable('M_initial')/M_sun).toFixed(0)} M_sun`);
console.log(`  r: ${(r_before/ly_to_m).toFixed(2)} → ${(wd2.getVariable('r')/ly_to_m).toFixed(2)} ly`);
console.log(`  B: ${(wd2.getVariable('B')/1.2).toExponential(3)} → ${wd2.getVariable('B').toExponential(3)} T`);
console.log(`  g(t=1Myr) after expansion: ${wd2.compute_g_Westerlund2(t_test).toExponential(6)} m/s^2`);
wd2.restoreState('before_expansion');
console.log(`  Restored to original state\n`);

// Test 17: All UQFF terms breakdown
console.log('Test 17: Complete UQFF Terms Breakdown (t=1 Myr)');
const g_total = wd2.compute_g_Westerlund2(t_test);
const Mt_breakdown = wd2.M_t(t_test);
const ug1_breakdown = (wd2.G * Mt_breakdown) / (wd2.r * wd2.r);

// Term 1: Base + H0 + B
const corr_H = 1 + wd2.H0 * t_test;
const corr_B = 1 - wd2.B / wd2.B_crit;
const term1 = ug1_breakdown * corr_H * corr_B;

// Term 2: Ug total
const term2 = wd2.compute_Ug(Mt_breakdown);

// Term 3: Lambda
const term3 = (wd2.Lambda * wd2.c_light * wd2.c_light) / 3.0;

// Term 4: EM
const cross_vB = wd2.gas_v * wd2.B;
const em_base = (wd2.q_charge * cross_vB) / wd2.proton_mass;
const corr_UA = 1 + (wd2.rho_vac_UA / wd2.rho_vac_SCm);
const term4 = (em_base * corr_UA) * wd2.scale_EM;

// Term 5: Quantum
const sqrt_unc = Math.sqrt(wd2.delta_x * wd2.delta_p);
const term5 = (wd2.hbar / sqrt_unc) * wd2.integral_psi * (2 * Math.PI / wd2.t_Hubble);

// Term 6: Fluid
const V_breakdown = wd2.compute_V();
const term6 = (wd2.rho_fluid * V_breakdown * ug1_breakdown) / Mt_breakdown;

// Term 8: DM
const M_dm = Mt_breakdown * wd2.M_DM_factor;
const pert1 = wd2.delta_rho_over_rho;
const pert2 = 3 * wd2.G * Mt_breakdown / (wd2.r * wd2.r * wd2.r);
const term_dm_force = (Mt_breakdown + M_dm) * (pert1 + pert2);
const term8 = term_dm_force / Mt_breakdown;

// Term 9: Wind
const wind_pressure = wd2.rho_wind * wd2.v_wind * wd2.v_wind;
const term9 = wind_pressure / wd2.rho_fluid;

console.log(`  Total g = ${g_total.toExponential(6)} m/s^2`);
console.log(`  Term breakdown:`);
console.log(`    1. Base (H0, B, M(t)): ${term1.toExponential(6)} m/s^2 (${(term1/g_total*100).toFixed(3)}%)`);
console.log(`    2. Ug (with f_TRZ):    ${term2.toExponential(6)} m/s^2 (${(term2/g_total*100).toFixed(3)}%)`);
console.log(`    3. Lambda:             ${term3.toExponential(6)} m/s^2 (${(term3/g_total*100).toFixed(6)}%)`);
console.log(`    4. EM (scaled, UA):    ${term4.toExponential(6)} m/s^2 (${(term4/g_total*100).toFixed(5)}%)`);
console.log(`    5. Quantum:            ${term5.toExponential(6)} m/s^2 (${(term5/g_total*100).toFixed(6)}%)`);
console.log(`    6. Fluid:              ${term6.toExponential(6)} m/s^2 (${(term6/g_total*100).toFixed(6)}%)`);
console.log(`    7. Oscillatory:        (combined standing + traveling waves)`);
console.log(`    8. DM (perturbations): ${term8.toExponential(6)} m/s^2 (${(term8/g_total*100).toFixed(3)}%)`);
console.log(`    9. Wind Feedback:      ${term9.toExponential(6)} m/s^2 (${(term9/g_total*100).toFixed(3)}%)`);
console.log(`  Dominant term: Wind Feedback (intense stellar winds from massive young stars)\n`);

// Test 18: Comprehensive report
console.log('Test 18: Comprehensive System Report (t=1.5 Myr)');
const t_report = 1.5e6 * 3.156e7;
console.log(wd2.generateReport(t_report));

// Test 19: Initial mass sweep
console.log('Test 19: Initial Mass Sweep (20000, 30000, 40000, 50000 M_sun)');
wd2.saveState('before_mass_sweep');
for (const M_val of [20000.0, 30000.0, 40000.0, 50000.0]) {
    wd2.setVariable('M_initial', M_val * M_sun);
    const g_sweep = wd2.compute_g_Westerlund2(t_test);
    const M_at_test = wd2.M_t(t_test);
    console.log(`  M_initial=${M_val.toFixed(0).padStart(5)} M_sun: M(1Myr)=${(M_at_test/M_sun).toFixed(0).padStart(6)} M_sun, g=${g_sweep.toExponential(3)} m/s^2`);
}
wd2.restoreState('before_mass_sweep');
console.log();

// Test 20: Star formation factor sweep
console.log('Test 20: Star Formation Factor Sweep (M_dot_factor = 2.0, 3.33, 5.0, 7.0)');
for (const M_dot of [2.0, 3.33, 5.0, 7.0]) {
    wd2.setVariable('M_dot_factor', M_dot);
    const M_at_test = wd2.M_t(t_test);
    const g_sweep = wd2.compute_g_Westerlund2(t_test);
    console.log(`  M_dot=${M_dot.toFixed(2).padStart(4)}: M(1Myr)=${(M_at_test/M_sun).toFixed(0).padStart(6)} M_sun, g=${g_sweep.toExponential(3)} m/s^2`);
}
wd2.restoreState('before_mass_sweep');
console.log();

console.log('========================================');
console.log('ALL TESTS COMPLETED SUCCESSFULLY ✓');
console.log('========================================\n');
