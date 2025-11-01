/**
 * Test Suite for NGC 2264 UQFF Module
 * 
 * Validates comprehensive star-forming nebula physics with wind, erosion, pillar dynamics,
 * protostar formation, and quantum wave structure.
 * 
 * Tests: 16 comprehensive categories covering all physics components
 * Expected: All tests PASSED
 */

const NGC2264UQFFModule = require('./ngc2264_uqff.js');

function assert(condition, message) {
  if (!condition) {
    console.error(`❌ FAILED: ${message}`);
    return false;
  }
  console.log(`✓ ${message}`);
  return true;
}

function assertClose(actual, expected, tolerance, message) {
  if (Math.abs(actual - expected) <= tolerance) {
    console.log(`✓ ${message}`);
    return true;
  } else {
    console.error(`❌ FAILED: ${message} (expected ${expected}, got ${actual})`);
    return false;
  }
}

let passed = 0;
let failed = 0;

function test(result, message) {
  if (result) {
    passed++;
  } else {
    failed++;
  }
}

// ==============================================================================
// TEST 1: Module Instantiation and Default Initialization
// ==============================================================================
console.log("\n=== TEST 1: Module Instantiation ===");
{
  const module = new NGC2264UQFFModule();
  
  test(assert(module !== null, "Module instantiation successful"), "instantiation");
  test(assert(module.G > 0, "Gravitational constant initialized"), "G constant");
  test(assert(module.c === 3e8, "Speed of light initialized correctly"), "c constant");
  test(assert(module.M_visible > 0, "Visible mass initialized"), "M_visible");
  test(assert(module.M_DM > 0, "Dark matter mass initialized"), "M_DM");
  test(assert(module.M === module.M_visible + module.M_DM, "Total mass correct"), "M total");
  test(assert(module.r > 0, "Radius initialized"), "radius");
  test(assert(module.v_wind === 20e3, "Stellar wind velocity set to 20 km/s"), "v_wind");
  test(assert(module.SFR > 0, "Star formation rate initialized"), "SFR");
}

// ==============================================================================
// TEST 2: NGC 2264 Physical Parameters Validation
// ==============================================================================
console.log("\n=== TEST 2: NGC 2264 Parameter Validation ===");
{
  const module = new NGC2264UQFFModule();
  
  test(assert((module.M / module.M_sun) === 100, "Total mass = 100 M☉"), "mass");
  test(assert((module.M_visible / module.M_sun) === 80, "Visible mass = 80 M☉"), "M_visible");
  test(assert((module.M_DM / module.M_sun) === 20, "Dark matter = 20 M☉"), "M_DM");
  test(assertClose(module.r, 3.31e16, 1e14, "Radius ≈ 3.31e16 m (~3.5 ly)"), "radius");
  test(assert((module.SFR / (module.M_sun / module.year_to_s)) === 0.01, "SFR = 0.01 M☉/yr"), "SFR");
  test(assert(module.z === 0.0008, "Redshift z = 0.0008"), "redshift");
}

// ==============================================================================
// TEST 3: Stellar Wind and Environmental Parameters
// ==============================================================================
console.log("\n=== TEST 3: Stellar Wind & Environmental Parameters ===");
{
  const module = new NGC2264UQFFModule();
  
  test(assert(module.v_wind === 20e3, "Stellar wind velocity = 20 km/s"), "v_wind");
  test(assert(module.v_r === 1e3, "Radial velocity = 1 km/s"), "v_r");
  test(assert(module.rho_fluid === 1e-20, "Dust/gas density = 1e-20 kg/m³"), "rho_fluid");
  test(assert(module.B === 1e-5, "Magnetic field = 1e-5 T"), "B field");
  test(assert(module.omega_spin === 1e-5, "Protostar spin = 1e-5 rad/s"), "omega_spin");
}

// ==============================================================================
// TEST 4: Pillar Wave Parameters
// ==============================================================================
console.log("\n=== TEST 4: Pillar Wave Parameters ===");
{
  const module = new NGC2264UQFFModule();
  
  test(assert(module.omega === 1e-14, "Pillar frequency = 1e-14 rad/s"), "omega");
  test(assert(module.sigma === 1e15, "Gaussian scale σ = 1e15 m"), "sigma");
  test(assert(module.A === 1e-10, "Pillar amplitude A = 1e-10"), "amplitude");
}

// ==============================================================================
// TEST 5: Variable Update and Cascade Effects
// ==============================================================================
console.log("\n=== TEST 5: Variable Update and Cascade ===");
{
  const module = new NGC2264UQFFModule();
  
  module.updateVariable('M_visible', 1e31);
  test(assert(module.M_visible === 1e31, "updateVariable: M_visible set"), "update M_visible");
  test(assert(module.M === 1e31 + module.M_DM, "updateVariable: M cascaded"), "cascade M");
  
  module.updateVariable('Delta_x', 1e-11);
  test(assert(module.Delta_p === module.hbar / 1e-11, "updateVariable: Delta_p recalculated"), "cascade Delta_p");
}

// ==============================================================================
// TEST 6: AddToVariable and SubtractFromVariable
// ==============================================================================
console.log("\n=== TEST 6: Variable Arithmetic ===");
{
  const module = new NGC2264UQFFModule();
  const initial_M = module.M;
  
  module.addToVariable('M', 1e30);
  test(assert(module.M === initial_M + 1e30, "addToVariable: M incremented"), "add to M");
  
  module.subtractFromVariable('M', 5e29);
  test(assert(module.M === initial_M + 1e30 - 5e29, "subtractFromVariable: M decremented"), "subtract from M");
}

// ==============================================================================
// TEST 7: Cosmological and Environmental Computations
// ==============================================================================
console.log("\n=== TEST 7: Cosmological Computations ===");
{
  const module = new NGC2264UQFFModule();
  
  const Hz = module.computeHtz(module.z);
  test(assert(Hz > 0, "Hubble parameter H(z) > 0"), "Hz positive");
  
  const t = 3e6 * module.year_to_s;
  const msf = module.computeMsfFactor(t);
  test(assert(msf > 0, "Mass growth factor > 0"), "msf positive");
  
  const rt = module.computeRt(t);
  test(assertClose(rt, module.r + module.v_r * t, 1e10, "Radius evolution r(t)"), "rt");
  
  const f_env = module.computeFenv(t);
  test(assert(f_env > 0, "Environmental forcing > 0"), "Fenv positive");
}

// ==============================================================================
// TEST 8: Ug Gravity Components
// ==============================================================================
console.log("\n=== TEST 8: Ug Gravity Components ===");
{
  const module = new NGC2264UQFFModule();
  
  const Ug1 = module.computeUg1(module.t);
  test(assert(typeof Ug1 === 'number', "Ug1 (dipole) computed"), "Ug1");
  
  const Ug2 = module.computeUg2(module.t);
  test(assert(Ug2 > 0, "Ug2 (superconductor) > 0"), "Ug2");
  
  const Ug3 = module.computeUg3prime(module.t);
  test(assert(Ug3 > 0, "Ug3' (stellar wind) > 0"), "Ug3");
  
  const Ug4 = module.computeUg4(module.t);
  test(assert(Ug4 > 0, "Ug4 (reaction) > 0"), "Ug4");
}

// ==============================================================================
// TEST 9: Aether Vacuum Integration (Ui)
// ==============================================================================
console.log("\n=== TEST 9: Aether Vacuum Integration ===");
{
  const module = new NGC2264UQFFModule();
  
  const Ui = module.computeUi(module.t);
  test(assert(typeof Ui === 'number', "Ui (aether vacuum) computed"), "Ui");
  test(assert(!isNaN(Ui), "Ui is valid number"), "Ui valid");
}

// ==============================================================================
// TEST 10: Pillar Wave Structure (Psi Integral)
// ==============================================================================
console.log("\n=== TEST 10: Pillar Wave Structure ===");
{
  const module = new NGC2264UQFFModule();
  
  const r_test = module.r / 10;  // Test at smaller radius
  const t_test = 1e6 * module.year_to_s;
  
  const psi = module.computePsiIntegral(r_test, t_test);
  test(assert(psi > 0, "Pillar wave intensity |ψ|² > 0"), "psi positive");
  test(assert(psi <= module.A * module.A, "Pillar intensity bounded by A²"), "psi bounded");
  
  const psi_at_center = module.computePsiIntegral(0, t_test);
  test(assert(psi_at_center > 0, "Pillar intensity at center > 0"), "psi center");
}

// ==============================================================================
// TEST 11: Quantum Gravity Term
// ==============================================================================
console.log("\n=== TEST 11: Quantum Gravity Term ===");
{
  const module = new NGC2264UQFFModule();
  const t = 1e6 * module.year_to_s;
  
  const q_term = module.computeQuantumTerm(module.t_Hubble, module.r);
  test(assert(typeof q_term === 'number', "Quantum term computed"), "quantum");
  test(assert(!isNaN(q_term), "Quantum term is valid"), "quantum valid");
}

// ==============================================================================
// TEST 12: Fluid and Dark Matter Terms
// ==============================================================================
console.log("\n=== TEST 12: Fluid and Dark Matter Terms ===");
{
  const module = new NGC2264UQFFModule();
  
  const f_term = module.computeFluidTerm(1e-5);
  test(assert(typeof f_term === 'number', "Fluid term computed"), "fluid");
  
  const dm_term = module.computeDMTerm(module.r);
  test(assert(dm_term > 0, "Dark matter term > 0"), "DM term");
}

// ==============================================================================
// TEST 13: Ug Sum (All Universal Gravity Components)
// ==============================================================================
console.log("\n=== TEST 13: Ug Sum Computation ===");
{
  const module = new NGC2264UQFFModule();
  
  const ug_sum = module.computeUgSum(module.r);
  test(assert(typeof ug_sum === 'number', "Ug sum computed"), "ug_sum");
  test(assert(ug_sum > 0, "Ug sum > 0"), "ug_sum positive");
}

// ==============================================================================
// TEST 14: Full Gravity Field Computation (computeG)
// ==============================================================================
console.log("\n=== TEST 14: Full Gravity Field (computeG) ===");
{
  const module = new NGC2264UQFFModule();
  const t = 3e6 * module.year_to_s;  // 3 Myr
  
  const g = module.computeG(t);
  test(assert(g > 0 && isFinite(g), "g_NGC2264 > 0 and finite"), "g positive finite");
  test(assert(g > 1e20, "g_NGC2264 > 1e20 m/s² (strong field from 100 M☉)"), "g magnitude");
  
  const g_at_different_r = module.computeG(t, 1e16);
  test(assert(g_at_different_r > g, "g at smaller radius > g at standard radius"), "g alt radius");
}

// ==============================================================================
// TEST 15: Time Evolution of Gravity
// ==============================================================================
console.log("\n=== TEST 15: Time Evolution ===");
{
  const module = new NGC2264UQFFModule();
  
  const times = [0, 1e6 * module.year_to_s, 3e6 * module.year_to_s, 1e7 * module.year_to_s];
  let all_valid = true;
  
  for (const t of times) {
    const g = module.computeG(t);
    if (!isFinite(g)) {
      all_valid = false;
    }
  }
  
  test(assert(all_valid, "Gravity valid at all time points"), "time evolution");
}

// ==============================================================================
// TEST 16: Documentation and Introspection
// ==============================================================================
console.log("\n=== TEST 16: Documentation and Introspection ===");
{
  const module = new NGC2264UQFFModule();
  
  const summary = module.getSummary();
  test(assert(summary.name === 'NGC 2264 Cone Nebula UQFF Module', "getSummary name"), "summary name");
  test(assert(summary.current_system === 'NGC2264', "getSummary system"), "summary system");
  test(assert(Array.isArray(summary.physicsComponents) && summary.physicsComponents.length > 10, 
    "getSummary lists physics components"), "summary components");
  
  const equation = module.getEquationText();
  test(assert(equation.includes('g_NGC2264') && equation.includes('pillar'), "getEquationText complete"), "equation text");
  
  module.printVariables();  // Ensure no crash
  test(true, "printVariables executes without error");
}

// ==============================================================================
// SUMMARY
// ==============================================================================
console.log(`\n${'='.repeat(70)}`);
console.log(`TEST RESULTS: ${passed} PASSED, ${failed} FAILED`);
console.log(`Total: ${passed + failed} tests`);
console.log(`Success Rate: ${((passed / (passed + failed)) * 100).toFixed(2)}%`);
console.log(`${'='.repeat(70)}\n`);

if (failed === 0) {
  console.log(`✓✓✓ ALL TESTS PASSED ✓✓✓`);
  console.log(`NGC 2264 UQFF Module (System #49) Ready for Integration`);
  process.exit(0);
} else {
  console.log(`❌ SOME TESTS FAILED`);
  process.exit(1);
}
