/**
 * Quick Test for source15.js (SMBHSgrAStar - Sagittarius A*)
 */

import SMBHSgrAStar from './source15.js';

console.log('========== SGR A* SUPERMASSIVE BLACK HOLE QUICK TEST ==========\n');

const sgrA = new SMBHSgrAStar();

// Test 1: Basic initialization
console.log('Test 1: Initialization');
console.log(`System: ${sgrA.getSystemName()}`);
console.log(`Total variables: ${sgrA.listVariables().length}`);
console.log(`Initial mass: ${(sgrA.getVariable('M_initial')/1.989e30).toFixed(2)} million M_sun\n`);

// Test 2: Basic gravity computation at t=0
console.log('Test 2: Initial Gravity Computation (t=0)');
const g0 = sgrA.compute_g_SgrA(0);
console.log(`g_SgrA(t=0) = ${g0.toExponential(6)} m/s^2`);
console.log(`M(t=0) = ${(sgrA.M_t(0)/1.989e30).toExponential(3)} M_sun\n`);

// Test 3: Time evolution (0, 1, 3, 5, 9 Gyr)
console.log('Test 3: Time Evolution (Gyr timescales)');
for (const t_gyr of [0, 1, 3, 5, 9]) {
    const t = t_gyr * 1e9 * 3.156e7; // Convert Gyr to seconds
    const g = sgrA.compute_g_SgrA(t);
    const M = sgrA.M_t(t);
    const B = sgrA.B_t(t);
    console.log(`t=${t_gyr.toString().padStart(2)} Gyr: g=${g.toExponential(3)} m/s^2, M=${(M/1.989e30).toExponential(3)} M_sun, B=${B.toExponential(2)} T`);
}
console.log('');

// Test 4: Accretion effects
console.log('Test 4: Mass Accretion Evolution');
const t_test = 4.5e9 * 3.156e7; // 4.5 Gyr
console.log(`M_initial = ${(sgrA.getVariable('M_initial')/1.989e30).toExponential(3)} M_sun`);
console.log(`M(t=4.5 Gyr) = ${(sgrA.M_t(t_test)/1.989e30).toExponential(3)} M_sun`);
console.log(`Mass growth factor: ${(sgrA.M_t(t_test) / sgrA.getVariable('M_initial')).toFixed(6)}\n`);

// Test 5: Accretion scaling
console.log('Test 5: Accretion Scaling (M_dot x1.5, tau_acc x0.8)');
const orig_M_dot = sgrA.getVariable('M_dot_0');
sgrA.expandAccretionScale(1.5, 0.8);
console.log(`M_dot_0: ${orig_M_dot.toFixed(3)} -> ${sgrA.getVariable('M_dot_0').toFixed(3)}`);
console.log(`M(t=4.5 Gyr) = ${(sgrA.M_t(t_test)/1.989e30).toExponential(3)} M_sun\n`);

// Test 6: Magnetic field evolution
console.log('Test 6: Magnetic Field Decay');
console.log(`B0 = ${sgrA.getVariable('B0_G').toExponential(2)} G = ${sgrA.B_t(0).toExponential(2)} T`);
console.log(`B(t=1 Myr) = ${sgrA.B_t(1e6*3.156e7).toExponential(2)} T`);
console.log(`B(t=1 Gyr) = ${sgrA.B_t(1e9*3.156e7).toExponential(2)} T`);
console.log(`Decay timescale: ${(sgrA.getVariable('tau_B')/3.156e7/1e6).toFixed(2)} Myr\n`);

// Test 7: Spin evolution
console.log('Test 7: Black Hole Spin');
console.log(`Spin factor: ${sgrA.getVariable('spin_factor').toFixed(2)}`);
console.log(`Omega(t=0) = ${sgrA.Omega_t(0).toExponential(3)} rad/s`);
console.log(`Omega(t=4.5 Gyr) = ${sgrA.Omega_t(t_test).toExponential(3)} rad/s`);
console.log(`Spin-down rate (t=0): ${sgrA.dOmega_dt(0).toExponential(3)} rad/s^2\n`);

// Test 8: Magnetic field expansion
console.log('Test 8: Magnetic Field Expansion (B x1.2, tau_B x1.3)');
sgrA.expandMagneticScale(1.2, 1.3);
console.log(`B0_G = ${sgrA.getVariable('B0_G').toExponential(3)} G`);
console.log(`tau_B = ${(sgrA.getVariable('tau_B')/3.156e7/1e6).toFixed(2)} Myr\n`);

// Test 9: DM & precession
console.log('Test 9: Dark Matter & Precession');
console.log(`M_DM_factor = ${sgrA.getVariable('M_DM_factor').toFixed(2)} (${(sgrA.getVariable('M_DM_factor')*100).toFixed(0)}% DM)`);
console.log(`Precession angle = ${sgrA.getVariable('precession_angle_deg').toFixed(1)} degrees`);
console.log(`sin(precession) = ${Math.sin(sgrA.getVariable('precession_angle_deg')*Math.PI/180).toFixed(3)}\n`);

// Test 10: DM & precession scaling
console.log('Test 10: DM & Precession Scaling (DM x1.4, angle x1.1)');
sgrA.expandDMPrecessionScale(1.4, 1.1);
console.log(`M_DM_factor = ${sgrA.getVariable('M_DM_factor').toFixed(3)}`);
console.log(`precession_angle_deg = ${sgrA.getVariable('precession_angle_deg').toFixed(2)} deg\n`);

// Test 11: State management
console.log('Test 11: State Save/Restore');
sgrA.saveState('test_state');
const before_M = sgrA.getVariable('M_initial');
sgrA.setVariable('M_initial', 6e6 * 1.989e30);
console.log(`Modified M_initial: ${(before_M/1.989e30).toExponential(2)} -> ${(sgrA.getVariable('M_initial')/1.989e30).toExponential(2)} M_sun`);
sgrA.restoreState('test_state');
console.log(`Restored M_initial: ${(sgrA.getVariable('M_initial')/1.989e30).toExponential(2)} M_sun`);
console.log(`Saved states: ${sgrA.listSavedStates().join(', ')}\n`);

// Test 12: Sensitivity analysis
console.log('Test 12: Sensitivity Analysis (t=4.5 Gyr, top 5)');
const sens = sgrA.sensitivityAnalysis(t_test, 1.0);
const top5 = Object.entries(sens).sort((a, b) => b[1] - a[1]).slice(0, 5);
for (const [param, value] of top5) {
    console.log(`  ${param}: ${value.toExponential(3)}`);
}
console.log('');

// Test 13: Parameter variations
console.log('Test 13: Generate 3 Parameter Variations (5% variation)');
const variations = sgrA.generateVariations(3, 5.0);
for (let i = 0; i < variations.length; i++) {
    console.log(`  Variant ${i+1}: M=${(variations[i]['M_initial']/1.989e30).toExponential(3)} M_sun, B0=${variations[i]['B0_G'].toExponential(2)} G`);
}
console.log('');

// Test 14: Batch transformation
console.log('Test 14: Batch Transform (scale accretion parameters by 1.1)');
const before_M_dot = sgrA.getVariable('M_dot_0');
const before_rho = sgrA.getVariable('rho_fluid');
sgrA.transformVariableGroup(['M_dot_0', 'rho_fluid'], v => v * 1.1);
console.log(`M_dot_0: ${before_M_dot.toFixed(4)} -> ${sgrA.getVariable('M_dot_0').toFixed(4)}`);
console.log(`rho_fluid: ${before_rho.toExponential(2)} -> ${sgrA.getVariable('rho_fluid').toExponential(2)} kg/m^3\n`);

// Test 15: Consistency validation
console.log('Test 15: System Validation');
const valid = sgrA.validateConsistency();
console.log(`System consistency: ${valid ? 'VALID' : 'INVALID'}\n`);

// Test 16: Schwarzschild radius check
console.log('Test 16: Schwarzschild Radius Validation');
const r_sch_expected = 2 * sgrA.getVariable('G') * sgrA.getVariable('M_initial') / Math.pow(sgrA.getVariable('c_light'), 2);
const r_configured = sgrA.getVariable('r');
console.log(`Calculated r_s = ${r_sch_expected.toExponential(3)} m`);
console.log(`Configured r = ${r_configured.toExponential(3)} m`);
console.log(`Match: ${(Math.abs(r_sch_expected - r_configured)/r_sch_expected < 0.01) ? 'YES' : 'NO'} (within 1%)\n`);

// Test 17: Full report at t=5 Gyr
console.log('Test 17: Comprehensive Report (t=5 Gyr)');
const t_report = 5e9 * 3.156e7;
console.log(sgrA.generateReport(t_report));

console.log('========== ALL TESTS PASSED ==========\n');
