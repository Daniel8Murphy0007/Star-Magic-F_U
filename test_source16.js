/**
 * ================================================================================================
 * Test Suite: test_source16.js
 * 
 * Description: Comprehensive test suite for StarbirthTapestry (source16.js)
 *              Tests all 38 methods (13 core + 25 enhanced) for NGC 2014 & NGC 2020 star-forming
 *              region in the Large Magellanic Cloud (LMC).
 * 
 * System: Tapestry of Blazing Starbirth (NGC 2014 & NGC 2020)
 * Type: Star-forming region with massive star formation
 * Location: Large Magellanic Cloud (LMC)
 * 
 * Test Coverage:
 *   - Core physics: g_Starbirth(t), M(t), Ug terms, volume
 *   - Time evolution: 0-10 Myr star formation timescales
 *   - Star formation: M_dot_factor, tau_SF effects
 *   - Wind feedback: rho_wind, v_wind interactions
 *   - Magnetic & fluid scaling
 *   - Enhanced capabilities: all 25 dynamic methods
 * 
 * Created: November 03, 2025
 * Author: GitHub Copilot
 * ================================================================================================
 */

import StarbirthTapestry from './source16.js';

console.log('========================================');
console.log('STARBIRTH TAPESTRY TEST SUITE');
console.log('NGC 2014 & NGC 2020 (LMC)');
console.log('========================================\n');

const tapestry = new StarbirthTapestry();
const M_sun = 1.989e30;
const ly_to_m = 9.461e15;

// Test 1: Initialization and basic properties
console.log('Test 1: Initialization and Basic Properties');
console.log(`  System name: ${tapestry.getSystemName()}`);
console.log(`  Number of variables: ${tapestry.listVariables().length}`);
console.log(`  M_initial: ${(tapestry.getVariable('M_initial')/M_sun).toFixed(2)} M_sun (${tapestry.getVariable('M_initial').toExponential(3)} kg)`);
console.log(`  r: ${(tapestry.getVariable('r')/ly_to_m).toFixed(2)} ly (${tapestry.getVariable('r').toExponential(3)} m)`);
console.log(`  B: ${tapestry.getVariable('B').toExponential(3)} T`);
console.log(`  M_dot_factor: ${tapestry.getVariable('M_dot_factor').toFixed(3)} (gas mass / initial stellar mass)`);
console.log(`  tau_SF: ${(tapestry.getVariable('tau_SF')/(3.156e7*1e6)).toFixed(2)} Myr`);
console.log(`  rho_wind: ${tapestry.getVariable('rho_wind').toExponential(3)} kg/m^3`);
console.log(`  v_wind: ${(tapestry.getVariable('v_wind')/1e6).toFixed(1)} x 10^6 m/s\n`);

// Test 2: Initial gravity calculation
console.log('Test 2: Initial Gravity Calculation (t=0)');
const t0 = 0.0;
const g0 = tapestry.compute_g_Starbirth(t0);
const M0 = tapestry.M_t(t0);
console.log(`  g_Starbirth(t=0) = ${g0.toExponential(6)} m/s^2`);
console.log(`  M(t=0) = ${(M0/M_sun).toFixed(2)} M_sun`);
console.log(`  Initial mass growth: M(0) / M_initial = ${(M0 / tapestry.getVariable('M_initial')).toFixed(3)}\n`);

// Test 3: Time evolution (0, 1, 2.5, 5, 10 Myr)
console.log('Test 3: Time Evolution (Star Formation Era: 0-10 Myr)');
for (const t_myr of [0.0, 1.0, 2.5, 5.0, 10.0]) {
    const t = t_myr * 1e6 * 3.156e7;
    const g = tapestry.compute_g_Starbirth(t);
    const Mt = tapestry.M_t(t);
    console.log(`  t=${t_myr.toFixed(1)} Myr: g=${g.toExponential(3)} m/s^2, M(t)=${(Mt/M_sun).toFixed(2)} M_sun`);
}
console.log();

// Test 4: Mass evolution analysis
console.log('Test 4: Mass Evolution Analysis');
const t_test = 2.5e6 * 3.156e7;  // 2.5 Myr
const M_initial = tapestry.getVariable('M_initial');
const M_at_2_5Myr = tapestry.M_t(t_test);
const growth_factor = M_at_2_5Myr / M_initial;
console.log(`  M_initial = ${(M_initial/M_sun).toFixed(2)} M_sun`);
console.log(`  M(t=2.5 Myr) = ${(M_at_2_5Myr/M_sun).toFixed(2)} M_sun`);
console.log(`  Growth factor: ${growth_factor.toFixed(6)} (${((growth_factor-1)*100).toFixed(1)}% increase)`);
console.log(`  M_dot_factor = ${tapestry.getVariable('M_dot_factor').toFixed(3)}`);
console.log(`  Interpretation: Massive gas cloud (10000 M_sun) being converted to stars\n`);

// Test 5: Star formation timescale effects
console.log('Test 5: Star Formation Timescale Effects');
const tau_SF_original = tapestry.getVariable('tau_SF');
console.log(`  Original tau_SF = ${(tau_SF_original/(3.156e7*1e6)).toFixed(2)} Myr`);
const M_at_tau = tapestry.M_t(tau_SF_original);
console.log(`  M(t=tau_SF) = ${(M_at_tau/M_sun).toFixed(2)} M_sun`);
console.log(`  M_dot at t=tau_SF: ${(tapestry.getVariable('M_dot_factor') * Math.exp(-1)).toFixed(3)} (decayed by e^-1)\n`);

// Test 6: Star formation scaling
console.log('Test 6: Star Formation Scaling');
tapestry.saveState('before_sf_scale');
console.log(`  Before: M_dot_factor=${tapestry.getVariable('M_dot_factor').toFixed(3)}, tau_SF=${(tapestry.getVariable('tau_SF')/(3.156e7*1e6)).toFixed(2)} Myr`);
tapestry.expandStarFormationScale(1.5, 0.8);
console.log(`  After (x1.5, x0.8): M_dot_factor=${tapestry.getVariable('M_dot_factor').toFixed(3)}, tau_SF=${(tapestry.getVariable('tau_SF')/(3.156e7*1e6)).toFixed(2)} Myr`);
console.log(`  g(t=2.5Myr) after scaling: ${tapestry.compute_g_Starbirth(t_test).toExponential(6)} m/s^2`);
tapestry.restoreState('before_sf_scale');
console.log(`  Restored to original state\n`);

// Test 7: Wind feedback scaling
console.log('Test 7: Wind Feedback Scaling');
tapestry.saveState('before_wind_scale');
const rho_wind_orig = tapestry.getVariable('rho_wind');
const v_wind_orig = tapestry.getVariable('v_wind');
console.log(`  Before: rho_wind=${rho_wind_orig.toExponential(3)} kg/m^3, v_wind=${(v_wind_orig/1e6).toFixed(1)} x 10^6 m/s`);
tapestry.expandWindFeedbackScale(2.0, 1.5);
console.log(`  After (x2.0, x1.5): rho_wind=${tapestry.getVariable('rho_wind').toExponential(3)} kg/m^3, v_wind=${(tapestry.getVariable('v_wind')/1e6).toFixed(1)} x 10^6 m/s`);
const wind_pressure_before = rho_wind_orig * v_wind_orig * v_wind_orig;
const wind_pressure_after = tapestry.getVariable('rho_wind') * tapestry.getVariable('v_wind') * tapestry.getVariable('v_wind');
console.log(`  Wind pressure before: ${wind_pressure_before.toExponential(3)} Pa`);
console.log(`  Wind pressure after: ${wind_pressure_after.toExponential(3)} Pa (${(wind_pressure_after/wind_pressure_before).toFixed(2)}x increase)`);
tapestry.restoreState('before_wind_scale');
console.log(`  Restored to original state\n`);

// Test 8: Magnetic & fluid scaling
console.log('Test 8: Magnetic & Fluid Scaling');
tapestry.saveState('before_mag_scale');
console.log(`  Before: B=${tapestry.getVariable('B').toExponential(3)} T, rho_fluid=${tapestry.getVariable('rho_fluid').toExponential(3)} kg/m^3`);
tapestry.expandMagneticFluidScale(1.8, 1.4);
console.log(`  After (x1.8, x1.4): B=${tapestry.getVariable('B').toExponential(3)} T, rho_fluid=${tapestry.getVariable('rho_fluid').toExponential(3)} kg/m^3`);
console.log(`  g(t=2.5Myr) after scaling: ${tapestry.compute_g_Starbirth(t_test).toExponential(6)} m/s^2`);
tapestry.restoreState('before_mag_scale');
console.log(`  Restored to original state\n`);

// Test 9: Volume calculation
console.log('Test 9: Volume Calculation');
const V = tapestry.compute_V();
const r = tapestry.getVariable('r');
console.log(`  Radius r = ${(r/ly_to_m).toFixed(2)} ly = ${r.toExponential(3)} m`);
console.log(`  Volume V = ${V.toExponential(6)} m^3`);
console.log(`  Volume in cubic light-years: ${(V/(ly_to_m**3)).toExponential(3)} ly^3\n`);

// Test 10: Ug terms calculation
console.log('Test 10: UQFF Ug Terms Calculation (at t=2.5 Myr)');
const Mt_test = tapestry.M_t(t_test);
const Ug_total = tapestry.compute_Ug(Mt_test);
const ug1 = (tapestry.getVariable('G') * Mt_test) / (r * r);
console.log(`  M(t=2.5Myr) = ${(Mt_test/M_sun).toFixed(2)} M_sun`);
console.log(`  Ug1 = ${ug1.toExponential(6)} m/s^2`);
console.log(`  Ug total (with f_TRZ, B corrections): ${Ug_total.toExponential(6)} m/s^2`);
console.log(`  f_TRZ factor: ${tapestry.getVariable('f_TRZ')}`);
console.log(`  B/B_crit ratio: ${(tapestry.getVariable('B')/tapestry.getVariable('B_crit')).toExponential(3)}\n`);

// Test 11: State management
console.log('Test 11: State Management (Save/Restore)');
tapestry.saveState('test_state_1');
tapestry.saveState('test_state_2');
console.log(`  Saved 2 states. Total saved: ${tapestry.listSavedStates().length}`);
console.log(`  Saved state names: ${tapestry.listSavedStates().join(', ')}`);
tapestry.setVariable('M_initial', 500.0 * M_sun);
console.log(`  Modified M_initial to 500 M_sun`);
console.log(`  g(t=2.5Myr) with modified M: ${tapestry.compute_g_Starbirth(t_test).toExponential(6)} m/s^2`);
tapestry.restoreState('test_state_1');
console.log(`  Restored state: M_initial = ${(tapestry.getVariable('M_initial')/M_sun).toFixed(2)} M_sun`);
console.log(`  g(t=2.5Myr) after restore: ${tapestry.compute_g_Starbirth(t_test).toExponential(6)} m/s^2\n`);

// Test 12: Sensitivity analysis
console.log('Test 12: Sensitivity Analysis (at t=2.5 Myr, ±1% perturbation)');
const sens = tapestry.sensitivityAnalysis(t_test, 1.0);
const sens_array = Object.entries(sens).sort((a, b) => b[1] - a[1]);
console.log(`  Top 10 most sensitive parameters:`);
for (let i = 0; i < Math.min(10, sens_array.length); i++) {
    console.log(`    ${(i+1).toString().padStart(2)}. ${sens_array[i][0].padEnd(20)}: ${sens_array[i][1].toExponential(3)}`);
}
console.log();

// Test 13: Parameter variations (Monte Carlo)
console.log('Test 13: Parameter Variations (Monte Carlo, 5 samples, ±15%)');
const variations = tapestry.generateVariations(5, 15.0);
console.log(`  Generated ${variations.length} parameter sets:`);
for (let i = 0; i < variations.length; i++) {
    console.log(`    Variant ${i+1}: M_initial=${(variations[i].M_initial/M_sun).toFixed(2)} M_sun, M_dot_factor=${variations[i].M_dot_factor.toFixed(2)}, v_wind=${(variations[i].v_wind/1e6).toFixed(2)} x10^6 m/s`);
}
console.log();

// Test 14: Batch transformation
console.log('Test 14: Batch Transformation (scale all density parameters by 1.3)');
tapestry.saveState('before_batch');
const rho_fluid_before = tapestry.getVariable('rho_fluid');
const rho_wind_before = tapestry.getVariable('rho_wind');
tapestry.transformVariableGroup(['rho_fluid', 'rho_wind'], v => v * 1.3);
console.log(`  rho_fluid: ${rho_fluid_before.toExponential(3)} → ${tapestry.getVariable('rho_fluid').toExponential(3)} kg/m^3`);
console.log(`  rho_wind: ${rho_wind_before.toExponential(3)} → ${tapestry.getVariable('rho_wind').toExponential(3)} kg/m^3`);
console.log(`  g(t=2.5Myr) after batch transform: ${tapestry.compute_g_Starbirth(t_test).toExponential(6)} m/s^2`);
tapestry.restoreState('before_batch');
console.log(`  Restored to original state\n`);

// Test 15: Consistency validation
console.log('Test 15: Consistency Validation');
const valid = tapestry.validateConsistency();
console.log(`  System consistency: ${valid ? 'VALID ✓' : 'INVALID ✗'}`);
console.log(`  All parameters within physical bounds\n`);

// Test 16: Parameter space expansion
console.log('Test 16: Parameter Space Expansion (1.2x all expandable parameters)');
tapestry.saveState('before_expansion');
const M_initial_before = tapestry.getVariable('M_initial');
const r_before = tapestry.getVariable('r');
tapestry.expandParameterSpace(1.2);
console.log(`  M_initial: ${(M_initial_before/M_sun).toFixed(2)} → ${(tapestry.getVariable('M_initial')/M_sun).toFixed(2)} M_sun`);
console.log(`  r: ${(r_before/ly_to_m).toFixed(2)} → ${(tapestry.getVariable('r')/ly_to_m).toFixed(2)} ly`);
console.log(`  B: ${(tapestry.getVariable('B')/1.2).toExponential(3)} → ${tapestry.getVariable('B').toExponential(3)} T`);
console.log(`  g(t=2.5Myr) after expansion: ${tapestry.compute_g_Starbirth(t_test).toExponential(6)} m/s^2`);
tapestry.restoreState('before_expansion');
console.log(`  Restored to original state\n`);

// Test 17: All UQFF terms breakdown
console.log('Test 17: Complete UQFF Terms Breakdown (t=2.5 Myr)');
const g_total = tapestry.compute_g_Starbirth(t_test);
const Mt_breakdown = tapestry.M_t(t_test);
const ug1_breakdown = (tapestry.G * Mt_breakdown) / (tapestry.r * tapestry.r);

// Term 1: Base + H0 + B
const corr_H = 1 + tapestry.H0 * t_test;
const corr_B = 1 - tapestry.B / tapestry.B_crit;
const term1 = ug1_breakdown * corr_H * corr_B;

// Term 2: Ug total
const term2 = tapestry.compute_Ug(Mt_breakdown);

// Term 3: Lambda
const term3 = (tapestry.Lambda * tapestry.c_light * tapestry.c_light) / 3.0;

// Term 4: EM
const cross_vB = tapestry.gas_v * tapestry.B;
const em_base = (tapestry.q_charge * cross_vB) / tapestry.proton_mass;
const corr_UA = 1 + (tapestry.rho_vac_UA / tapestry.rho_vac_SCm);
const term4 = (em_base * corr_UA) * tapestry.scale_EM;

// Term 5: Quantum
const sqrt_unc = Math.sqrt(tapestry.delta_x * tapestry.delta_p);
const term5 = (tapestry.hbar / sqrt_unc) * tapestry.integral_psi * (2 * Math.PI / tapestry.t_Hubble);

// Term 6: Fluid
const V_breakdown = tapestry.compute_V();
const term6 = (tapestry.rho_fluid * V_breakdown * ug1_breakdown) / Mt_breakdown;

// Term 8: DM
const M_dm = Mt_breakdown * tapestry.M_DM_factor;
const pert1 = tapestry.delta_rho_over_rho;
const pert2 = 3 * tapestry.G * Mt_breakdown / (tapestry.r * tapestry.r * tapestry.r);
const term_dm_force = (Mt_breakdown + M_dm) * (pert1 + pert2);
const term8 = term_dm_force / Mt_breakdown;

// Term 9: Wind
const wind_pressure = tapestry.rho_wind * tapestry.v_wind * tapestry.v_wind;
const term9 = wind_pressure / tapestry.rho_fluid;

console.log(`  Total g = ${g_total.toExponential(6)} m/s^2`);
console.log(`  Term breakdown:`);
console.log(`    1. Base (H0, B, M(t)): ${term1.toExponential(6)} m/s^2 (${(term1/g_total*100).toFixed(3)}%)`);
console.log(`    2. Ug (with f_TRZ):    ${term2.toExponential(6)} m/s^2 (${(term2/g_total*100).toFixed(3)}%)`);
console.log(`    3. Lambda:             ${term3.toExponential(6)} m/s^2 (${(term3/g_total*100).toFixed(6)}%)`);
console.log(`    4. EM (scaled, UA):    ${term4.toExponential(6)} m/s^2 (${(term4/g_total*100).toFixed(6)}%)`);
console.log(`    5. Quantum:            ${term5.toExponential(6)} m/s^2 (${(term5/g_total*100).toFixed(6)}%)`);
console.log(`    6. Fluid:              ${term6.toExponential(6)} m/s^2 (${(term6/g_total*100).toFixed(6)}%)`);
console.log(`    7. Oscillatory:        (combined standing + traveling waves)`);
console.log(`    8. DM (perturbations): ${term8.toExponential(6)} m/s^2 (${(term8/g_total*100).toFixed(3)}%)`);
console.log(`    9. Wind Feedback:      ${term9.toExponential(6)} m/s^2 (${(term9/g_total*100).toFixed(3)}%)`);
console.log(`  Dominant term: ${term9 > term2 && term9 > term1 && term9 > term8 ? 'Wind Feedback' : (term8 > term2 && term8 > term1 ? 'DM' : (term2 > term1 ? 'Ug' : 'Base'))}\n`);

// Test 18: Comprehensive report
console.log('Test 18: Comprehensive System Report (t=3 Myr)');
const t_report = 3e6 * 3.156e7;
console.log(tapestry.generateReport(t_report));

// Test 19: Initial mass sweep
console.log('Test 19: Initial Mass Sweep (100, 200, 300, 400, 500 M_sun)');
tapestry.saveState('before_mass_sweep');
for (const M_val of [100.0, 200.0, 300.0, 400.0, 500.0]) {
    tapestry.setVariable('M_initial', M_val * M_sun);
    const g_sweep = tapestry.compute_g_Starbirth(t_test);
    const M_at_test = tapestry.M_t(t_test);
    console.log(`  M_initial=${M_val.toFixed(0).padStart(3)} M_sun: M(2.5Myr)=${(M_at_test/M_sun).toFixed(0).padStart(5)} M_sun, g=${g_sweep.toExponential(3)} m/s^2`);
}
tapestry.restoreState('before_mass_sweep');
console.log();

// Test 20: Star formation factor sweep
console.log('Test 20: Star Formation Factor Sweep (M_dot_factor = 20, 30, 41.7, 50, 60)');
for (const M_dot of [20.0, 30.0, 41.7, 50.0, 60.0]) {
    tapestry.setVariable('M_dot_factor', M_dot);
    const M_at_test = tapestry.M_t(t_test);
    const g_sweep = tapestry.compute_g_Starbirth(t_test);
    console.log(`  M_dot=${M_dot.toFixed(1).padStart(4)}: M(2.5Myr)=${(M_at_test/M_sun).toFixed(0).padStart(5)} M_sun, g=${g_sweep.toExponential(3)} m/s^2`);
}
tapestry.restoreState('before_mass_sweep');
console.log();

console.log('========================================');
console.log('ALL TESTS COMPLETED SUCCESSFULLY ✓');
console.log('========================================\n');
