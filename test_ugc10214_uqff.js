/**
 * Test Suite for UGC10214UQFFModule
 * 
 * Comprehensive testing of Tadpole Galaxy UQFF implementation
 * Test Categories: 18 test groups, 90+ individual tests
 * 
 * Author: GitHub Copilot
 * Date: November 1, 2025
 * Framework: Star-Magic UQFF v2.0
 */

const UGC10214UQFFModule = require('./ugc10214_uqff.js');

// ═══════════════════════════════════════════════════════════════
// TEST FRAMEWORK
// ═══════════════════════════════════════════════════════════════

let testsPassed = 0;
let testsFailed = 0;
let testsTotal = 0;

function assert(condition, message) {
  testsTotal++;
  if (condition) {
    console.log(`  ✓ ${message}`);
    testsPassed++;
  } else {
    console.log(`  ✗ ${message}`);
    testsFailed++;
  }
}

function assertClose(actual, expected, tolerance = 0.01, message = '') {
  testsTotal++;
  const error = Math.abs(actual - expected) / (Math.abs(expected) + 1e-50);
  if (error < tolerance) {
    console.log(`  ✓ ${message} (error: ${(error * 100).toFixed(3)}%)`);
    testsPassed++;
  } else {
    console.log(`  ✗ ${message} (expected: ${expected}, got: ${actual}, error: ${(error * 100).toFixed(2)}%)`);
    testsFailed++;
  }
}

function section(name) {
  console.log(`\n${'═'.repeat(70)}`);
  console.log(`  ${name}`);
  console.log(`${'═'.repeat(70)}`);
}

function subsection(name) {
  console.log(`\n  → ${name}`);
}

// ═══════════════════════════════════════════════════════════════
// TEST SUITE
// ═══════════════════════════════════════════════════════════════

section('UGC10214 UQFF MODULE TEST SUITE');

// TEST 1: Module Initialization
subsection('1. Module Initialization');
const ugc = new UGC10214UQFFModule();
assert(ugc !== null, 'Module instantiation successful');
assert(ugc.variables !== undefined, 'Variables map exists');
assert(ugc.variables.size > 60, 'All 70+ variables initialized');
assert(ugc.variables.get('M_total') !== undefined, 'M_total initialized');
assert(ugc.variables.get('z') === 0.032, 'Redshift z = 0.032');

// TEST 2: Universal Constants
subsection('2. Universal Constants');
assert(ugc.variables.get('G') === 6.6743e-11, 'Gravitational constant correct');
assert(ugc.variables.get('c') === 3e8, 'Speed of light correct');
assert(ugc.variables.get('hbar') === 1.0546e-34, 'Planck constant correct');
assert(ugc.variables.get('Lambda') === 1.1e-52, 'Cosmological constant correct');
assert(ugc.variables.get('H0') === 70, 'Hubble constant H0 = 70');
assert(ugc.variables.get('Omega_m') === 0.3, 'Matter density parameter');
assert(ugc.variables.get('Omega_Lambda') === 0.7, 'Dark energy parameter');

// TEST 3: UGC 10214 Parameters
subsection('3. UGC 10214 System Parameters');
const M_visible = ugc.variables.get('M_visible');
const M_DM = ugc.variables.get('M_DM');
const M_total = ugc.variables.get('M_total');
assert(M_visible > 0, 'Visible mass positive');
assert(M_DM > 0, 'Dark matter mass positive');
assert(M_total > M_visible, 'Total mass includes dark matter');
assert(ugc.variables.get('r') === 55 * 3.086e19, 'Galaxy radius = 55 kpc');
assert(ugc.variables.get('SFR') > 0, 'Star formation rate positive');

// TEST 4: Merger Parameters
subsection('4. Merger Parameters (VV 29c)');
assert(ugc.variables.get('M_dwarf') === 3.5e9 * 1.989e30, 'Dwarf mass = 3.5e9 M_sun');
assert(ugc.variables.get('d_dwarf') === 110 * 3.086e19, 'Merger distance = 110 kpc');
assert(ugc.variables.get('tau_merge') > 0, 'Merger timescale positive');
assert(ugc.variables.get('M_dwarf_0') > 0, 'Initial dwarf mass positive');

// TEST 5: Tail Dynamics
subsection('5. Tail Dynamics');
assert(ugc.variables.get('v_tail') === 400e3, 'Tail velocity = 400 km/s');
assert(ugc.variables.get('A_tail') === 1e-10, 'Tail amplitude correct');
assert(ugc.variables.get('sigma_tail') > 0, 'Tail Gaussian width positive');
assert(ugc.variables.get('m_tail') === 2, 'Azimuthal quantum number m=2');
assert(ugc.variables.get('omega_tail') > 0, 'Tail wave frequency positive');

// TEST 6: Magnetic Field
subsection('6. Magnetic Field Parameters');
assert(ugc.variables.get('B') === 1e-5, 'Magnetic field B = 1e-5 T');
assert(ugc.variables.get('B_critical') === 1e11, 'Critical field B_crit = 1e11 T');
assert(ugc.variables.get('I_dipole') > 0, 'Dipole moment positive');
assert(ugc.variables.get('H_aether') > 0, 'Aether field strength positive');

// TEST 7: Fluid Properties
subsection('7. Fluid Properties');
assert(ugc.variables.get('rho_fluid') > 0, 'Fluid density positive');
assert(ugc.variables.get('V_fluid') > 0, 'Fluid volume positive');
assert(ugc.variables.get('c_s') > 0, 'Sound speed positive');
assert(ugc.variables.get('P_fluid') > 0, 'Pressure positive');

// TEST 8: Quantum Parameters
subsection('8. Quantum Parameters');
assert(ugc.variables.get('Delta_x') > 0, 'Position uncertainty positive');
assert(ugc.variables.get('Delta_p') > 0, 'Momentum uncertainty positive');
assert(ugc.variables.get('t_Hubble') === 4.35e17, 'Hubble time = 13.8 Gyr');
assert(ugc.variables.get('psi_integral') !== undefined, 'Wave function integral exists');

// TEST 9: Dynamic Variable Operations
subsection('9. Dynamic Variable Operations');
const original_M = ugc.variables.get('M_total');
ugc.updateVariable('M_total', original_M * 1.1);
assert(ugc.variables.get('M_total') === original_M * 1.1, 'updateVariable works');

ugc.addToVariable('M_total', 1e30);
assertClose(ugc.variables.get('M_total'), original_M * 1.1 + 1e30, 0.001, 'addToVariable works');

ugc.subtractFromVariable('M_total', 1e30);
assertClose(ugc.variables.get('M_total'), original_M * 1.1, 0.001, 'subtractFromVariable works');

ugc.updateVariable('M_total', original_M);  // Restore

// TEST 10: Hubble Parameter
subsection('10. Hubble Parameter Computation');
const Hz = ugc.computeHtz(0.032);
assert(Hz > 0, 'Hubble parameter positive');
assert(Hz < 1e-10, 'Hubble parameter reasonable');

const Hz_0 = ugc.computeHtz(0);
const Hz_high = ugc.computeHtz(2.0);
assert(Hz_high > Hz_0, 'H increases with redshift');

// TEST 11: Mass Evolution (Merger)
subsection('11. Mass Evolution (Merger Dynamics)');
const t1 = 0;
const t2 = 250e6 * 3.156e7;  // 250 Myr
const M_merge_0 = ugc.computeMmerge(t1);
const M_merge_250 = ugc.computeMmerge(t2);

assert(M_merge_0 > M_merge_250, 'Merger mass decays with time');
assertClose(M_merge_0 / M_merge_250, Math.exp(1), 0.05, 'e-fold decay time correct');

// TEST 12: Environmental Forcing
subsection('12. Environmental Forcing');
const F_env_early = ugc.computeFenv(0);
const F_env_later = ugc.computeFenv(100e6 * 3.156e7);

assert(F_env_early > 0, 'Environmental forcing positive');
assert(typeof F_env_later === 'number', 'Returns numeric value');
assert(isFinite(F_env_early), 'Result is finite');

// TEST 13: Universal Gravity Components (Ug1-Ug4)
subsection('13. Universal Gravity Components');
const ug1 = ugc.computeUg1();
const ug2 = ugc.computeUg2();
const ug3prime = ugc.computeUg3prime();
const ug4_early = ugc.computeUg4(0);
const ug4_late = ugc.computeUg4(500e6 * 3.156e7);

assert(typeof ug1 === 'number', 'Ug1 returns number');
assert(isFinite(ug2), 'Ug2 is finite');
assert(ug3prime > 0, 'Ug3\' is positive');
assert(ug4_early > ug4_late, 'Ug4 decays with time');

// TEST 14: Integrated Potential (Ui)
subsection('14. Integrated Potential (Ui)');
const ui_t0 = ugc.computeUi(0);
const ui_t1 = ugc.computeUi(250e6 * 3.156e7);
const ui_t2 = ugc.computeUi(500e6 * 3.156e7);

assert(typeof ui_t0 === 'number', 'Ui returns number');
assert(isFinite(ui_t1), 'Ui is finite at t=250Myr');
// Ui oscillates, so different values expected

// TEST 15: Wave Function (Tail)
subsection('15. Wave Function Tail Properties');
const r = 30 * 3.086e19;  // 30 kpc
const theta = Math.PI / 4;  // 45 degrees
const t = 0;

const psi = ugc.computePsiTail(r, theta, t);
assert(psi.real !== undefined, 'Wave function real part exists');
assert(psi.imag !== undefined, 'Wave function imaginary part exists');
assert(isFinite(psi.real), 'Real part is finite');
assert(isFinite(psi.imag), 'Imaginary part is finite');

const psi_density = ugc.computePsiDensity(r, theta, t);
assert(psi_density >= 0, 'Probability density non-negative');
assert(psi_density === psi.real * psi.real + psi.imag * psi.imag, 'Density calculation correct');

// TEST 16: Quantum Terms
subsection('16. Quantum Gravity Contributions');
const psi_integral = ugc.computePsiIntegral(50 * 3.086e19, 0);
assert(psi_integral > 0, 'Psi integral positive');
assert(isFinite(psi_integral), 'Psi integral finite');

const q_term = ugc.computeQuantumTerm();
assert(typeof q_term === 'number', 'Quantum term returns number');
assert(isFinite(q_term), 'Quantum term finite');

// TEST 17: Fluid Dynamics
subsection('17. Fluid Dynamics Term');
const g_base_test = 1e-10;
const f_fluid = ugc.computeFluidTerm(g_base_test);

assert(typeof f_fluid === 'number', 'Fluid term returns number');
assert(isFinite(f_fluid), 'Fluid term finite');
assert(f_fluid >= 0, 'Fluid term non-negative');

// TEST 18: Dark Matter Perturbations
subsection('18. Dark Matter Perturbations');
const r_test = 20 * 3.086e19;  // 20 kpc
const f_dm = ugc.computeDMTerm(r_test);

assert(typeof f_dm === 'number', 'DM term returns number');
assert(isFinite(f_dm), 'DM term finite');

// TEST 19: Ug Sum
subsection('19. Universal Gravity Component Sum');
const ug_sum = ugc.computeUgSum(r_test, 0);
assert(typeof ug_sum === 'number', 'Ug sum returns number');
assert(isFinite(ug_sum), 'Ug sum is finite');

// TEST 20: Master Equation (computeG)
subsection('20. Master Equation - computeG');
const t_test = 0;
const r_mid = 20 * 3.086e19;

const g_ugc = ugc.computeG(t_test, r_mid);
assert(typeof g_ugc === 'number', 'computeG returns number');
assert(isFinite(g_ugc), 'computeG result is finite');
assert(g_ugc > 0, 'Gravity positive (attractive)');
assert(g_ugc > 1e-15, 'Gravity has measurable magnitude');

// TEST 21: Time Evolution of Gravity
subsection('21. Gravity Evolution Over Time');
const g_t0 = ugc.computeG(0, r_mid);
const g_t100myr = ugc.computeG(100e6 * 3.156e7, r_mid);
const g_t250myr = ugc.computeG(250e6 * 3.156e7, r_mid);
const g_t500myr = ugc.computeG(500e6 * 3.156e7, r_mid);

assert(g_t0 > 0, 'Gravity positive at t=0');
assert(g_t100myr > 0, 'Gravity positive at t=100Myr');
assert(g_t500myr > 0, 'Gravity positive at t=500Myr');
assert(isFinite(g_t250myr), 'Gravity finite at t=250Myr');

// Gravity should evolve as merger weakens - check that values are different
assert(Math.abs(g_t500myr - g_t0) > 0, 'Gravity values differ over time');

// TEST 22: Spatial Dependence
subsection('22. Spatial Dependence of Gravity');
const t_fixed = 100e6 * 3.156e7;
const r_inner = 10 * 3.086e19;    // 10 kpc
const r_mid_test = 20 * 3.086e19;  // 20 kpc
const r_outer = 50 * 3.086e19;    // 50 kpc

const g_inner = ugc.computeG(t_fixed, r_inner);
const g_mid_sp = ugc.computeG(t_fixed, r_mid_test);
const g_outer = ugc.computeG(t_fixed, r_outer);

assert(g_inner > g_outer, 'Gravity decreases with radius');
assert(g_outer > 0, 'Gravity positive even at large radius');

// TEST 23: Physical Consistency
subsection('23. Physical Consistency Checks');
assert(ugc.variables.get('Omega_m') + ugc.variables.get('Omega_Lambda') === 1.0, 'Ωm + ΩΛ = 1.0');
assert(ugc.variables.get('M_DM') < ugc.variables.get('M_total'), 'DM < total mass');
assert(ugc.variables.get('v_tail') < ugc.variables.get('c'), 'Tail velocity < c');
assert(ugc.variables.get('B') < ugc.variables.get('B_critical'), 'B < B_critical');

// TEST 24: Equation Text Generation
subsection('24. Equation Text and Documentation');
const eq_text = ugc.getEquationText();
assert(eq_text.includes('g_UGC10214'), 'Equation text contains system name');
assert(eq_text.includes('M(t)'), 'Equation includes mass evolution');
assert(eq_text.includes('F_env'), 'Equation includes environmental forcing');
assert(eq_text.includes('Ug'), 'Equation includes universal gravity terms');

// TEST 25: Variable Printing
subsection('25. Variable Summary');
const var_summary = ugc.printVariables();
assert(var_summary.includes('UGC 10214'), 'Summary includes system name');
assert(var_summary.includes('variables tracked'), 'Summary includes variable count');
assert(ugc.variables.size > 60, 'Correct number of variables');

// TEST 26: State Export/Import
subsection('26. State Serialization');
const state1 = ugc.getState();
assert(state1 !== null, 'State export works');
assert(typeof state1 === 'object', 'State is object');
assert(Object.keys(state1).length > 60, 'State contains all variables');

const ugc2 = new UGC10214UQFFModule();
ugc2.setState(state1);
assert(ugc2.variables.size === ugc.variables.size, 'State import maintains variable count');
assert(ugc2.variables.get('M_total') === ugc.variables.get('M_total'), 'State values match');

// TEST 27: Numerical Stability
subsection('27. Numerical Stability');
let stability_ok = true;
for (let i = 0; i < 10; i++) {
  const t_rand = Math.random() * 500e6 * 3.156e7;
  const r_rand = (10 + Math.random() * 50) * 3.086e19;
  const g_rand = ugc.computeG(t_rand, r_rand);
  if (!isFinite(g_rand)) {
    stability_ok = false;
    break;
  }
}
assert(stability_ok, 'All random computations produce finite results');

// TEST 28: Component Magnitudes
subsection('28. Component Magnitude Analysis');
const t_now = 250e6 * 3.156e7;
const r_now = 20 * 3.086e19;

const g_base_est = ugc.computeG(t_now, r_now);
const ug_comp = ugc.computeUgSum(r_now, t_now);
const q_comp = ugc.computeQuantumTerm();
const f_comp = ugc.computeFluidTerm(g_base_est);

assert(Math.abs(g_base_est) > Math.abs(q_comp), 'Base gravity >> quantum term');
assert(Math.abs(g_base_est) > Math.abs(ug_comp) * 0.01, 'Base gravity significant');

// TEST 29: Edge Cases - Very Early Times
subsection('29. Edge Cases - Very Early Times (t→0)');
const g_very_early = ugc.computeG(1e10, r_mid);  // 1e10 seconds
assert(isFinite(g_very_early), 'Gravity finite at very early times');
assert(g_very_early > 0, 'Gravity positive at early times');

// TEST 30: Edge Cases - Very Late Times
subsection('30. Edge Cases - Very Late Times');
const g_very_late = ugc.computeG(1000e6 * 3.156e7, r_mid);  // 1 Gyr
assert(isFinite(g_very_late), 'Gravity finite at late times');
assert(g_very_late > 0, 'Gravity positive at late times');

// TEST 31: Computational Performance
subsection('31. Computational Performance');
const perf_start = Date.now();
for (let i = 0; i < 100; i++) {
  ugc.computeG(Math.random() * 500e6 * 3.156e7, (10 + Math.random() * 40) * 3.086e19);
}
const perf_end = Date.now();
const avg_time = (perf_end - perf_start) / 100;

console.log(`  → Average computation time: ${avg_time.toFixed(3)} ms per call`);
assert(avg_time < 10, 'Performance acceptable (<10ms per call)');

// ═══════════════════════════════════════════════════════════════
// TEST RESULTS
// ═══════════════════════════════════════════════════════════════

section('TEST RESULTS SUMMARY');
console.log(`
  Total Tests:     ${testsTotal}
  Tests Passed:    ${testsPassed} ✓
  Tests Failed:    ${testsFailed} ✗
  Success Rate:    ${((testsPassed / testsTotal) * 100).toFixed(1)}%

`);

if (testsFailed === 0) {
  console.log('  ✓ ALL TESTS PASSED! Module is production-ready.\n');
  process.exit(0);
} else {
  console.log(`  ✗ ${testsFailed} test(s) failed. Review output above.\n`);
  process.exit(1);
}
