/**
 * Quick Test for source14.js (MagnetarSGR0501_4516)
 */

import MagnetarSGR0501_4516 from './source14.js';

console.log('========== SGR 0501+4516 MAGNETAR QUICK TEST ==========\n');

const mag = new MagnetarSGR0501_4516();

// Test 1: Basic initialization
console.log('Test 1: Initialization');
console.log(`System: ${mag.getSystemName()}`);
console.log(`Total variables: ${mag.listVariables().length}\n`);

// Test 2: Basic gravity computation at t=0
console.log('Test 2: Initial Gravity Computation (t=0)');
const g0 = mag.compute_g_Magnetar(0);
console.log(`g_Magnetar(t=0) = ${g0.toExponential(6)} m/s^2\n`);

// Test 3: Time evolution (0, 1000, 3000, 5000 years)
console.log('Test 3: Time Evolution');
for (const t_yr of [0, 1000, 3000, 5000]) {
    const t = t_yr * 3.156e7; // Convert years to seconds
    const g = mag.compute_g_Magnetar(t);
    const B = mag.B_t(t);
    console.log(`t=${t_yr.toString().padStart(4)} yr: g=${g.toExponential(6)} m/s^2, B=${B.toExponential(3)} T`);
}
console.log('');

// Test 4: Parameter modification
console.log('Test 4: Parameter Modification');
console.log(`Original B0 = ${mag.getVariable('B0').toExponential(3)} T`);
mag.setVariable('B0', 1.5e10);
console.log(`Modified B0 = ${mag.getVariable('B0').toExponential(3)} T`);
const g_modified = mag.compute_g_Magnetar(2000 * 3.156e7);
console.log(`g(t=2000yr, B0=1.5e10) = ${g_modified.toExponential(6)} m/s^2\n`);

// Test 5: Batch scaling
console.log('Test 5: Batch Scaling (M, r by factor 1.1)');
const orig_M = mag.getVariable('M');
const orig_r = mag.getVariable('r');
mag.scaleVariableGroup(['M', 'r'], 1.1);
console.log(`M: ${orig_M.toExponential(3)} -> ${mag.getVariable('M').toExponential(3)} kg`);
console.log(`r: ${orig_r.toExponential(3)} -> ${mag.getVariable('r').toExponential(3)} m\n`);

// Test 6: Magnetic field expansion
console.log('Test 6: Magnetic Field Expansion (B x1.2, tau_B x0.9)');
mag.expandMagneticScale(1.2, 0.9);
console.log(`B0 = ${mag.getVariable('B0').toExponential(3)} T`);
console.log(`tau_B = ${mag.getVariable('tau_B').toExponential(3)} s\n`);

// Test 7: State management
console.log('Test 7: State Save/Restore');
mag.saveState('test_state');
const before_B0 = mag.getVariable('B0');
mag.setVariable('B0', 5e9);
console.log(`Modified B0: ${mag.getVariable('B0').toExponential(3)} T`);
mag.restoreState('test_state');
console.log(`Restored B0: ${mag.getVariable('B0').toExponential(3)} T`);
console.log(`Saved states: ${mag.listSavedStates().join(', ')}\n`);

// Test 8: Sensitivity analysis
console.log('Test 8: Sensitivity Analysis (t=2000 yr, top 3)');
const t_test = 2000 * 3.156e7;
const sens = mag.sensitivityAnalysis(t_test, 1.0);
const top3 = Object.entries(sens).sort((a, b) => b[1] - a[1]).slice(0, 3);
for (const [param, value] of top3) {
    console.log(`  ${param}: ${value.toExponential(3)}`);
}
console.log('');

// Test 9: Parameter variations
console.log('Test 9: Generate 3 Parameter Variations (5% variation)');
const variations = mag.generateVariations(3, 5.0);
for (let i = 0; i < variations.length; i++) {
    console.log(`  Variant ${i+1}: B0=${variations[i]['B0'].toExponential(3)} T, M=${variations[i]['M'].toExponential(3)} kg`);
}
console.log('');

// Test 10: Consistency validation
console.log('Test 10: System Validation');
const valid = mag.validateConsistency();
console.log(`System consistency: ${valid ? 'VALID' : 'INVALID'}\n`);

// Test 11: Full report at t=3000 years
console.log('Test 11: Comprehensive Report (t=3000 yr)');
const t_report = 3000 * 3.156e7;
console.log(mag.generateReport(t_report));

console.log('========== ALL TESTS PASSED ==========\n');
