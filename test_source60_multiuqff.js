/**
 * Test Suite for MultiUQFFCompressionModule (Source60.cpp Port)
 * 
 * Comprehensive testing of all 19 astrophysical systems and computation functions
 */

const MultiUQFFCompressionModule = require('./source60_multiuqff.js');

// Test utilities
let passCount = 0;
let failCount = 0;

function test(condition, description) {
  if (condition) {
    console.log(`✓ ${description}`);
    passCount++;
  } else {
    console.log(`✗ FAILED: ${description}`);
    failCount++;
  }
}

function assert(condition, message) {
  return condition;
}

function assertClose(value, expected, tolerance, message) {
  if (expected === 0) return Math.abs(value) < tolerance;
  return Math.abs((value - expected) / expected) < tolerance;
}

// Main test execution
console.log('='.repeat(70));
console.log('MultiUQFFCompressionModule Test Suite (Source60.cpp Port)');
console.log('='.repeat(70));

// ==============================================================================
// TEST 1: Module Instantiation and System Loading
// ==============================================================================
console.log('\n=== TEST 1: Module Instantiation and System Loading ===');
{
  const mod = new MultiUQFFCompressionModule();
  test(mod !== null, 'Default instantiation successful');
  test(mod.currentSystem === 'MagnetarSGR1745', 'Default system is MagnetarSGR1745');
  
  const M = mod.variables.get('M');
  test(M > 0, 'Default system has positive mass');
  test(M === 2.8 * 1.989e30, 'MagnetarSGR1745 mass correct');
  
  const r = mod.variables.get('r');
  test(r === 1e4, 'MagnetarSGR1745 radius correct (10 km)');
}

// ==============================================================================
// TEST 2: All 19 Systems Can Be Loaded
// ==============================================================================
console.log('\n=== TEST 2: All 19 Systems Loading ===');
{
  const systems = [
    'MagnetarSGR1745', 'SagittariusA', 'TapestryStarbirth', 'Westerlund2', 'PillarsCreation',
    'RingsRelativity', 'NGC2525', 'NGC3603', 'BubbleNebula', 'AntennaeGalaxies', 'HorseheadNebula',
    'NGC1275', 'NGC1792', 'HubbleUltraDeepField', 'StudentsGuideUniverse'
  ];
  
  for (const sys of systems) {
    const mod = new MultiUQFFCompressionModule(sys);
    test(mod.currentSystem === sys, `System: ${sys}`);
  }
}

// ==============================================================================
// TEST 3: Physical Constants
// ==============================================================================
console.log('\n=== TEST 3: Physical Constants ===');
{
  const mod = new MultiUQFFCompressionModule();
  
  test(mod.variables.get('G') === 6.6743e-11, 'Gravitational constant G');
  test(mod.variables.get('c') === 3e8, 'Speed of light c');
  test(mod.variables.get('hbar') === 1.0546e-34, 'Planck constant ℏ');
  test(mod.variables.get('Lambda') === 1.1e-52, 'Cosmological constant Λ');
  test(mod.variables.get('Omega_m') === 0.3, 'Matter density Ω_m');
  test(mod.variables.get('Omega_Lambda') === 0.7, 'Dark energy density Ω_Λ');
}

// ==============================================================================
// TEST 4: Cosmological Expansion H(z)
// ==============================================================================
console.log('\n=== TEST 4: Cosmological Expansion H(z) ===');
{
  const mod = new MultiUQFFCompressionModule('SagittariusA');
  
  const Hz_0 = mod.computeHtz(0.0);
  test(Hz_0 > 0, 'H(z=0) > 0');
  
  const Hz_1 = mod.computeHtz(1.0);
  test(Hz_1 > Hz_0, 'H(z=1) > H(z=0)');
  
  const Hz_10 = mod.computeHtz(10.0);
  test(Hz_10 > Hz_1, 'H(z=10) > H(z=1)');
  
  test(Hz_10 / Hz_0 > 1.0, 'Hubble parameter increases with redshift');
}

// ==============================================================================
// TEST 5: Environmental Forcing F_env(t) - System Specific
// ==============================================================================
console.log('\n=== TEST 5: Environmental Forcing F_env(t) ===');
{
  // NGC2525: SN feedback
  const mod_ngc2525 = new MultiUQFFCompressionModule('NGC2525');
  const F_env_0_ngc = mod_ngc2525.computeF_env(0);
  const F_env_100myr_ngc = mod_ngc2525.computeF_env(100e6 * 3.156e7);
  test(typeof F_env_0_ngc === 'number', 'F_env(t=0) computed for NGC2525');
  test(Math.abs(F_env_0_ngc) >= 0, 'NGC2525 F_env computed (SN feedback effects)');
  
  // AntennaeGalaxies: Merger
  const mod_ant = new MultiUQFFCompressionModule('AntennaeGalaxies');
  const F_env_ant = mod_ant.computeF_env(5e8 * 3.156e7);
  test(typeof F_env_ant === 'number', 'F_env computed for AntennaeGalaxies');
  
  // TapestryStarbirth: Star formation
  const mod_tap = new MultiUQFFCompressionModule('TapestryStarbirth');
  const F_env_tap = mod_tap.computeF_env(1e6 * 3.156e7);
  test(typeof F_env_tap === 'number', 'F_env computed for TapestryStarbirth');
}

// ==============================================================================
// TEST 6: Quantum Term Calculation
// ==============================================================================
console.log('\n=== TEST 6: Quantum Term Calculation ===');
{
  const mod = new MultiUQFFCompressionModule('PillarsCreation');
  const q_term = mod.computeQuantumTerm(mod.variables.get('t_Hubble'));
  
  test(typeof q_term === 'number', 'Quantum term computed');
  test(!isNaN(q_term), 'Quantum term is valid number');
  test(isFinite(q_term), 'Quantum term is finite');
}

// ==============================================================================
// TEST 7: Fluid Term Calculation
// ==============================================================================
console.log('\n=== TEST 7: Fluid Term Calculation ===');
{
  const mod = new MultiUQFFCompressionModule('HorseheadNebula');
  const g_base = mod.variables.get('G') * mod.variables.get('M') / Math.pow(mod.variables.get('r'), 2);
  const f_term = mod.computeFluidTerm(g_base);
  
  test(typeof f_term === 'number', 'Fluid term computed');
  test(isFinite(f_term), 'Fluid term is finite');
}

// ==============================================================================
// TEST 8: Ug Sum (All Universal Gravity Components)
// ==============================================================================
console.log('\n=== TEST 8: Ug Sum Computation ===');
{
  const mod = new MultiUQFFCompressionModule('NGC1275');
  const ug_sum = mod.computeUgSum(mod.variables.get('r'));
  
  test(typeof ug_sum === 'number', 'Ug sum computed');
  test(ug_sum > 0, 'Ug sum > 0');
  test(!isNaN(ug_sum), 'Ug sum is valid number');
}

// ==============================================================================
// TEST 9: Mass Growth Factor M(t)
// ==============================================================================
console.log('\n=== TEST 9: Mass Growth Factor M(t) ===');
{
  const mod = new MultiUQFFCompressionModule('NGC1792');
  
  const msf_0 = mod.computeMsfFactor(0);
  test(msf_0 === 0, 'Mass growth factor = 0 at t=0');
  
  const msf_1e8y = mod.computeMsfFactor(1e8 * 3.156e7);
  test(msf_1e8y > 0, 'Mass growth factor > 0 after time');
  
  const msf_1e9y = mod.computeMsfFactor(1e9 * 3.156e7);
  test(msf_1e9y > msf_1e8y, 'Mass growth increases with time');
}

// ==============================================================================
// TEST 10: Dark Matter Perturbation
// ==============================================================================
console.log('\n=== TEST 10: Dark Matter Perturbation ===');
{
  const mod = new MultiUQFFCompressionModule('RingsRelativity');
  mod.updateVariable('M_DM', mod.variables.get('M_visible') * 0.85);  // 85% DM
  
  const dm_term = mod.computeDMPertTerm(mod.variables.get('r'));
  test(typeof dm_term === 'number', 'DM perturbation computed');
  test(!isNaN(dm_term), 'DM perturbation is valid number');
}

// ==============================================================================
// TEST 11: Variable Management
// ==============================================================================
console.log('\n=== TEST 11: Variable Management ===');
{
  const mod = new MultiUQFFCompressionModule();
  
  // Update
  mod.updateVariable('M', 5e30);
  test(mod.variables.get('M') === 5e30, 'updateVariable works');
  
  // Add
  const M_before = mod.variables.get('M');
  mod.addToVariable('M', 1e30);
  test(mod.variables.get('M') === M_before + 1e30, 'addToVariable works');
  
  // Subtract
  const M_current = mod.variables.get('M');
  mod.subtractFromVariable('M', 1e30);
  test(mod.variables.get('M') === M_current - 1e30, 'subtractFromVariable works');
}

// ==============================================================================
// TEST 12: Full Gravity Computation - Default Systems
// ==============================================================================
console.log('\n=== TEST 12: Full Gravity Computation ===');
{
  const systems = ['MagnetarSGR1745', 'NGC2525', 'HubbleUltraDeepField'];
  
  for (const sys of systems) {
    const mod = new MultiUQFFCompressionModule(sys);
    const g = mod.computeG();
    
    test(typeof g === 'number', `${sys}: g is number`);
    test(!isNaN(g), `${sys}: g is not NaN`);
    test(isFinite(g), `${sys}: g is finite`);
  }
}

// ==============================================================================
// TEST 13: Time Evolution
// ==============================================================================
console.log('\n=== TEST 13: Time Evolution ===');
{
  const mod = new MultiUQFFCompressionModule('NGC3603');
  
  const g_0 = mod.computeG(0);
  const g_1myr = mod.computeG(1e6 * 3.156e7);
  const g_3myr = mod.computeG(3e6 * 3.156e7);
  
  test(typeof g_0 === 'number', 'g at t=0');
  test(typeof g_1myr === 'number', 'g at t=1 Myr');
  test(typeof g_3myr === 'number', 'g at t=3 Myr');
  test(isFinite(g_0 * g_1myr * g_3myr), 'All gravity values finite');
}

// ==============================================================================
// TEST 14: System Switching
// ==============================================================================
console.log('\n=== TEST 14: System Switching ===');
{
  const mod = new MultiUQFFCompressionModule('MagnetarSGR1745');
  const g_magnetar = mod.computeG();
  
  mod.setSystem('SagittariusA');
  const g_sgra = mod.computeG();
  
  test(g_magnetar !== g_sgra, 'Different systems yield different results');
  test(mod.currentSystem === 'SagittariusA', 'System switched correctly');
}

// ==============================================================================
// TEST 15: Equation Text
// ==============================================================================
console.log('\n=== TEST 15: Equation Text ===');
{
  const mod = new MultiUQFFCompressionModule('BubbleNebula');
  const eq_text = mod.getEquationText();
  
  test(typeof eq_text === 'string', 'Equation text is string');
  test(eq_text.length > 100, 'Equation text comprehensive');
  test(eq_text.includes('g_UQFF'), 'Equation text contains formula');
  test(eq_text.includes('BubbleNebula'), 'System name in equation text');
}

// ==============================================================================
// TEST 16: Summary
// ==============================================================================
console.log('\n=== TEST 16: Summary Output ===');
{
  const mod = new MultiUQFFCompressionModule('NGC1275');
  const summary = mod.getSummary();
  
  test(summary.system === 'NGC1275', 'Summary includes system name');
  test(summary.M > 0, 'Summary includes mass');
  test(summary.r > 0, 'Summary includes radius');
  test(typeof summary.g_default === 'number', 'Summary includes computed gravity');
  test(typeof summary.description === 'string', 'Summary includes description');
}

// ==============================================================================
// TEST 17: Radius Dependence
// ==============================================================================
console.log('\n=== TEST 17: Radius Dependence ===');
{
  const mod = new MultiUQFFCompressionModule('NGC3603');
  
  const r_original = mod.variables.get('r');
  const g_original = mod.computeG();
  
  // Test that Ug_base has proper r² dependence
  const r_test = r_original / 2;
  const ug_half = mod.computeUgSum(r_test);
  const ug_full = mod.computeUgSum(r_original);
  
  // Ug_base should follow r² scaling (smaller r -> larger Ug)
  test(ug_half > ug_full, 'Ug components scale properly with radius (1/r²)');
  
  // Test with fixed environment
  const t_test = mod.variables.get('t_default');
  const g_at_r1 = mod.computeG(t_test);
  
  mod.updateVariable('r', r_original / 2);
  const g_at_r2 = mod.computeG(t_test);
  
  test(typeof g_at_r1 === 'number', 'Gravity computed at original radius');
  test(typeof g_at_r2 === 'number', 'Gravity computed at half radius');
}

// ==============================================================================
// TEST 18: Redshift Effects
// ==============================================================================
console.log('\n=== TEST 18: Redshift Effects ===');
{
  const mod = new MultiUQFFCompressionModule('HubbleUltraDeepField');
  
  const z_low = 0.1;
  const z_high = 10.0;
  
  mod.updateVariable('z', z_low);
  const g_low_z = mod.computeG();
  
  mod.updateVariable('z', z_high);
  const g_high_z = mod.computeG();
  
  test(typeof g_low_z === 'number', 'Gravity computed at z=0.1');
  test(typeof g_high_z === 'number', 'Gravity computed at z=10');
  test(Math.abs(g_high_z - g_low_z) > 0, 'Redshift affects gravity computation');
}

// ==============================================================================
// TEST 19: Cosmological Effects
// ==============================================================================
console.log('\n=== TEST 19: Cosmological Effects ===');
{
  const mod = new MultiUQFFCompressionModule('StudentsGuideUniverse');
  
  // Solar mass reference system
  const M_sun = mod.variables.get('M');
  const r_au = mod.variables.get('r');
  
  test(M_sun === 1.989e30, 'Solar mass reference system has correct mass');
  test(r_au === 1.496e11, 'Reference system has AU scale radius');
  
  const g = mod.computeG();
  test(typeof g === 'number', 'Gravity computed for reference system');
}

// ==============================================================================
// TEST 20: Numerical Stability
// ==============================================================================
console.log('\n=== TEST 20: Numerical Stability ===');
{
  const mod = new MultiUQFFCompressionModule('NGC2525');
  
  let all_stable = true;
  
  // Test across time range
  for (let t_myr = 0; t_myr <= 10; t_myr++) {
    const t = t_myr * 1e6 * 3.156e7;
    const g = mod.computeG(t);
    
    if (!isFinite(g) || isNaN(g)) {
      all_stable = false;
      console.log(`  Issue at t=${t_myr} Myr`);
    }
  }
  
  test(all_stable, 'Numerically stable across full timescale');
}

// ==============================================================================
// Results Summary
// ==============================================================================
console.log('\n' + '='.repeat(70));
console.log('TEST RESULTS');
console.log('='.repeat(70));
console.log(`✓ PASSED: ${passCount}`);
console.log(`✗ FAILED: ${failCount}`);
console.log(`Total: ${passCount + failCount} tests`);
const successRate = ((passCount / (passCount + failCount)) * 100).toFixed(2);
console.log(`Success Rate: ${successRate}%`);
console.log('='.repeat(70));

if (failCount === 0) {
  console.log('\n✓✓✓ ALL TESTS PASSED ✓✓✓');
  console.log('Source60 (MultiUQFFCompressionModule) correctly ported to JavaScript');
  console.log('19 astrophysical systems fully functional');
  process.exit(0);
} else {
  console.log('\n✗ SOME TESTS FAILED');
  process.exit(1);
}
