/**
 * Test Suite for UQFF Compressed Resonance Module
 * 
 * Validates multi-system support, dual-mode operation, recursion safety, and physics correctness
 * Tests: 14 comprehensive categories
 * 
 * Expected: All tests PASSED
 */

const UQFFCompressedResonanceModule = require('./uqff_compressed_resonance.js');

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
  const module = new UQFFCompressedResonanceModule();
  
  test(assert(module !== null, "Module instantiation successful"), "instantiation");
  test(assert(module.G > 0, "Gravitational constant initialized"), "G constant");
  test(assert(module.c === 3e8, "Speed of light initialized correctly"), "c constant");
  test(assert(module.current_system === "Guide", "Default system is Guide"), "default system");
  test(assert(module.mode === "compressed", "Default mode is compressed"), "default mode");
  test(assert(module.M > 0, "Default mass initialized"), "default mass");
  test(assert(module.r > 0, "Default radius initialized"), "default radius");
}

// ==============================================================================
// TEST 2: System Loading (All 8 Supported Systems)
// ==============================================================================
console.log("\n=== TEST 2: System Loading (8 Systems) ===");
{
  const module = new UQFFCompressedResonanceModule();
  const systems = ['YoungStars', 'Eagle', 'BigBang', 'M51', 'NGC1316', 'V838Mon', 'NGC1300', 'Guide'];
  
  for (const sys of systems) {
    module.setSystem(sys);
    test(assert(module.current_system === sys, `System '${sys}' loaded successfully`), `load ${sys}`);
    test(assert(module.M > 0, `${sys}: Mass > 0`), `${sys} mass positive`);
    test(assert(module.r > 0, `${sys}: Radius > 0`), `${sys} radius positive`);
  }
}

// ==============================================================================
// TEST 3: System-Specific Parameter Validation
// ==============================================================================
console.log("\n=== TEST 3: System Parameter Validation ===");
{
  const module = new UQFFCompressedResonanceModule();
  
  // YoungStars: Pre-main sequence
  module.setSystem('YoungStars');
  test(assert(module.M / module.M_sun > 500 && module.M / module.M_sun < 2000, 
    "YoungStars: ~1000 M☉"), "YoungStars mass");
  test(assert(module.SFR > 0, "YoungStars: SFR > 0"), "YoungStars SFR");
  
  // BigBang: Very large scale
  module.setSystem('BigBang');
  test(assert(module.z > 1000, "BigBang: z > 1000 (z~1100)"), "BigBang redshift");
  test(assert(module.r > 1e20, "BigBang: r > 1e20 m (vast cosmic scale)"), "BigBang radius");
  
  // NGC1300: Barred spiral
  module.setSystem('NGC1300');
  test(assert(module.M / module.M_sun > 1e10, "NGC1300: > 10^10 M☉"), "NGC1300 mass");
  test(assert(module.SFR > 0, "NGC1300: SFR > 0 (active star formation)"), "NGC1300 SFR");
}

// ==============================================================================
// TEST 4: Mode Switching (Compressed vs Resonance)
// ==============================================================================
console.log("\n=== TEST 4: Mode Switching ===");
{
  const module = new UQFFCompressedResonanceModule();
  module.setSystem('Guide');
  const t = 1e9 * module.year_to_s;
  
  // Compressed mode
  module.setMode('compressed');
  test(assert(module.mode === 'compressed', "Mode set to compressed"), "compressed mode");
  const g_compressed = module.computeG(t);
  test(assert(g_compressed > 0, "g_UQFF > 0 in compressed mode"), "g positive compressed");
  
  // Resonance mode (with resonance terms)
  module.setMode('resonance');
  test(assert(module.mode === 'resonance', "Mode set to resonance"), "resonance mode");
  const g_resonance = module.computeG(t);
  test(assert(g_resonance > 0 || g_resonance === g_compressed, 
    "g_UQFF computed without crashing in resonance mode"), "resonance no crash");
}

// ==============================================================================
// TEST 5: Recursion Safety (CRITICAL: Resonance does not recurse)
// ==============================================================================
console.log("\n=== TEST 5: Recursion Safety (CRITICAL) ===");
{
  const module = new UQFFCompressedResonanceModule();
  
  // Test resonance mode with multiple time steps
  module.setSystem('Guide');
  module.setMode('resonance');
  
  let crashed = false;
  const times = [1e6, 1e7, 1e8, 1e9, 1e10];
  
  for (const t of times) {
    try {
      const g = module.computeG(t);
      test(assert(typeof g === 'number' && !isNaN(g), 
        `Resonance mode: computeG(${t.toExponential(1)}) returned valid number`), 
        `no recursion at t=${t.toExponential(1)}`);
    } catch (e) {
      crashed = true;
      test(false, `Resonance mode crashed at t=${t}: ${e.message}`);
    }
  }
  
  test(assert(!crashed, "No stack overflow or recursion in resonance mode"), "recursion safety");
}

// ==============================================================================
// TEST 6: Variable Update and Cascade Effects
// ==============================================================================
console.log("\n=== TEST 6: Variable Update and Cascade ===");
{
  const module = new UQFFCompressedResonanceModule();
  const original_m_visible = module.M_visible;
  
  module.updateVariable('M', 1e42);
  test(assert(module.M === 1e42, "updateVariable: M set correctly"), "update M");
  test(assert(module.M_visible === 0.7 * 1e42, 
    "updateVariable: M_visible cascaded (0.7 * M)"), "cascade M_visible");
  test(assert(module.M_DM === 0.3 * 1e42, 
    "updateVariable: M_DM cascaded (0.3 * M)"), "cascade M_DM");
  
  const original_Delta_p = module.Delta_p;
  module.updateVariable('Delta_x', 1e-11);
  test(assert(module.Delta_p === module.hbar / 1e-11, 
    "updateVariable: Delta_p recalculated from Delta_x"), "cascade Delta_p");
}

// ==============================================================================
// TEST 7: AddToVariable and SubtractFromVariable
// ==============================================================================
console.log("\n=== TEST 7: Variable Arithmetic ===");
{
  const module = new UQFFCompressedResonanceModule();
  const initial_M = module.M;
  
  module.addToVariable('M', 1e40);
  test(assert(module.M === initial_M + 1e40, "addToVariable: M incremented"), "add to M");
  
  module.subtractFromVariable('M', 5e39);
  test(assert(module.M === initial_M + 1e40 - 5e39, 
    "subtractFromVariable: M decremented"), "subtract from M");
}

// ==============================================================================
// TEST 8: Computation Functions - All Physics Components
// ==============================================================================
console.log("\n=== TEST 8: Physics Component Computation ===");
{
  const module = new UQFFCompressedResonanceModule();
  module.setSystem('M51');
  const t = 1e9 * module.year_to_s;
  
  const Hz = module.computeHtz(module.z);
  test(assert(Hz > 0, "Hubble parameter H(z) > 0"), "Hz positive");
  
  const f_env = module.computeFenv(t);
  test(assert(f_env === 0.1, "Environmental forcing constant 0.1"), "Fenv");
  
  const ug_sum = module.computeUgSum();
  test(assert(ug_sum === 1e-10, "Ug sum placeholder 1e-10"), "UgSum");
  
  const psi = module.computePsiTotal(t);
  test(assert(typeof psi === 'number', "Total perturbation psi computed"), "psi");
  
  const q_term = module.computeQuantumTerm(module.t_Hubble);
  test(assert(typeof q_term === 'number', "Quantum term computed"), "quantum");
  
  const f_term = module.computeFluidTerm(1e-5);
  test(assert(f_term > 0, "Fluid dynamics term > 0"), "fluid term");
  
  const dm_term = module.computeDMTerm();
  test(assert(typeof dm_term === 'number', "Dark matter term computed"), "DM term");
}

// ==============================================================================
// TEST 9: Resonance Term (Safe, No Recursion)
// ==============================================================================
console.log("\n=== TEST 9: Resonance Term Safety ===");
{
  const module = new UQFFCompressedResonanceModule();
  module.setSystem('Guide');
  const t = 1e8 * module.year_to_s;
  const g_base = 1e-2;
  
  // Compressed mode: resonance term should be 0
  module.setMode('compressed');
  const res_compressed = module.computeResonanceTerm(t, g_base);
  test(assert(res_compressed === 0.0, 
    "Compressed mode: resonance term = 0"), "resonance compressed");
  
  // Resonance mode: resonance term should be non-zero (or at least computed)
  module.setMode('resonance');
  const res_resonance = module.computeResonanceTerm(t, g_base);
  test(assert(typeof res_resonance === 'number', 
    "Resonance mode: resonance term computed safely"), "resonance computed");
}

// ==============================================================================
// TEST 10: V838Mon Special Case (Light Echo Intensity)
// ==============================================================================
console.log("\n=== TEST 10: V838Mon Light Echo Special Case ===");
{
  const module = new UQFFCompressedResonanceModule();
  module.setSystem('V838Mon');
  const t = 365.25 * 24 * 3600;  // 1 year
  
  const intensity = module.computeG(t);
  test(assert(intensity > 0, "V838Mon: Light echo intensity > 0"), "I_echo positive");
  test(assert(intensity < 1e10, "V838Mon: Intensity physically reasonable"), "I_echo reasonable");
}

// ==============================================================================
// TEST 11: BigBang Special Case (Cosmic Expansion)
// ==============================================================================
console.log("\n=== TEST 11: BigBang Cosmological Expansion ===");
{
  const module = new UQFFCompressedResonanceModule();
  module.setSystem('BigBang');
  const t = 1000 * module.year_to_s;  // 1000 years after BB (still expanding)
  
  const r_before = module.r;
  const g = module.computeG(t);
  const r_after = module.r;
  
  test(assert(r_after === module.c * t, 
    "BigBang: Radius updated to c·t (light travel distance)"), "BigBang r=ct");
  test(assert(typeof g === 'number', "BigBang: g_UQFF computed"), "BigBang g");
}

// ==============================================================================
// TEST 12: ComputeG with Optional Radius
// ==============================================================================
console.log("\n=== TEST 12: ComputeG with Optional Radius ===");
{
  const module = new UQFFCompressedResonanceModule();
  module.setSystem('Guide');
  const t = 1e8 * module.year_to_s;
  const original_r = module.r;
  
  // Call with r_in
  const g1 = module.computeG(t, 1e12);
  test(assert(module.r === 1e12, "ComputeG with r_in: radius updated"), "r_in update");
  
  // Call with r_in = 0 (uses current r)
  const g2 = module.computeG(t, 0);
  test(assert(module.r === 1e12, "ComputeG with r_in=0: radius unchanged"), "r_in zero");
}

// ==============================================================================
// TEST 13: getSummary and getEquationText
// ==============================================================================
console.log("\n=== TEST 13: Summary and Documentation ===");
{
  const module = new UQFFCompressedResonanceModule();
  module.setSystem('NGC1300');
  module.setMode('resonance');
  
  const summary = module.getSummary();
  test(assert(summary.name === 'UQFF Compressed Resonance Module', 
    "Summary: name field correct"), "summary name");
  test(assert(summary.current_system === 'NGC1300', 
    "Summary: current system NGC1300"), "summary system");
  test(assert(summary.mode === 'resonance', 
    "Summary: mode resonance"), "summary mode");
  test(assert(Array.isArray(summary.systemsAvailable) && summary.systemsAvailable.length === 8,
    "Summary: 8 systems available"), "summary systems");
  
  const equation = module.getEquationText();
  test(assert(equation.includes('g_UQFF'), "Equation text includes g_UQFF formula"), "equation formula");
  test(assert(equation.includes('resonance'), "Equation text mentions resonance mode"), "equation resonance");
  test(assert(equation.includes('NGC1300'), "Equation text shows current system"), "equation system");
}

// ==============================================================================
// TEST 14: All Supported Systems - computeG Performance
// ==============================================================================
console.log("\n=== TEST 14: All Systems Compute Performance ===");
{
  const module = new UQFFCompressedResonanceModule();
  const systems = ['YoungStars', 'Eagle', 'BigBang', 'M51', 'NGC1316', 'V838Mon', 'NGC1300', 'Guide'];
  const t = 1e9 * module.year_to_s;
  
  for (const sys of systems) {
    module.setSystem(sys);
    try {
      const g = module.computeG(t);
      const is_valid = typeof g === 'number' && !isNaN(g) && isFinite(g);
      test(assert(is_valid, 
        `${sys}: computeG returns valid finite number (${g.toExponential(3)})`), 
        `${sys} compute`);
    } catch (e) {
      test(false, `${sys}: computeG crashed: ${e.message}`);
    }
  }
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
  process.exit(0);
} else {
  console.log(`❌ SOME TESTS FAILED`);
  process.exit(1);
}
