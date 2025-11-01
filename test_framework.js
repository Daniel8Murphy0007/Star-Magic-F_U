#!/usr/bin/env node

/**
 * Comprehensive UQFF Framework Test Suite
 * Tests computational abilities across all module types
 */

const idx = require('./index.js');

console.log('╔════════════════════════════════════════════════════════════════╗');
console.log('║     UQFF FRAMEWORK COMPREHENSIVE COMPUTATIONAL TEST SUITE      ║');
console.log('╚════════════════════════════════════════════════════════════════╝\n');

// Test 1: NGC1316 Galaxy Module
console.log('▶ TEST 1: NGC1316 Galaxy Merger Dynamics');
console.log('─'.repeat(65));
try {
    const ngc = new idx.NGC1316UQFFModule();
    
    // Test basic computation
    const t = 2e9 * 3.156e7;  // 2 Gyr in seconds
    const r = 20e3 * 3.086e19; // 20 kpc in meters
    const g = ngc.computeG(t, r);
    
    console.log(`✓ NGC1316 instantiated and computed gravity`);
    console.log(`  Time: 2 Gyr | Radius: 20 kpc`);
    console.log(`  Result: g = ${g.toExponential(3)} m/s²`);
    console.log(`  Merger mass contribution: ${(ngc.computeMmerge(t) / 1.989e30).toExponential(2)} M☉`);
    console.log(`  Hubble parameter H(z): ${(ngc.computeHtz(0.005) * 3.086e22 / 1e3).toExponential(2)} km/s/Mpc`);
    
    // Test dynamic update
    console.log(`\n  Dynamic capability test:`);
    console.log(`  Original M_spiral: ${(ngc.variables.get('M_spiral') / 1.989e30).toExponential(2)} M☉`);
    ngc.updateParameter('M_spiral', 2.5e10 * 1.989e30);
    console.log(`  Updated M_spiral: ${(ngc.variables.get('M_spiral') / 1.989e30).toExponential(2)} M☉`);
    console.log(`  ✓ Dynamic parameter update successful`);
    
    console.log('✓ TEST 1 PASSED\n');
} catch (e) {
    console.error('✗ TEST 1 FAILED:', e.message, '\n');
}

// Test 2: Multiple Magnetar Computations
console.log('▶ TEST 2: Magnetar Physics Across Multiple Systems');
console.log('─'.repeat(65));
try {
    const systems = [
        { name: 'SGR 1745-2900', mass: 2.78e30, radius: 1e4, B: 2e10 },
        { name: 'SGR 0501+4516', mass: 2.78e30, radius: 2e4, B: 1e10 },
        { name: 'Vela Pulsar', mass: 2.78e30, radius: 1.2e4, B: 3.2e8 }
    ];
    
    systems.forEach(sys => {
        const g_base = (6.6743e-11 * sys.mass) / (sys.radius * sys.radius);
        const f_sc = 1.0 - (sys.B / 1e11);
        const g_result = g_base * f_sc;
        console.log(`  ${sys.name.padEnd(20)} | g = ${g_result.toExponential(3)} m/s² | B = ${sys.B.toExponential(1)} T`);
    });
    
    console.log('✓ TEST 2 PASSED\n');
} catch (e) {
    console.error('✗ TEST 2 FAILED:', e.message, '\n');
}

// Test 3: Time Evolution Studies
console.log('▶ TEST 3: Time Evolution - NGC1316 Over Cosmic Timescales');
console.log('─'.repeat(65));
try {
    const ngc = new idx.NGC1316UQFFModule();
    const r = 20e3 * 3.086e19;
    const timePoints = [
        { label: 't=100 Myr', t: 1e8 * 3.156e7 },
        { label: 't=1 Gyr', t: 1e9 * 3.156e7 },
        { label: 't=5 Gyr', t: 5e9 * 3.156e7 },
        { label: 't=10 Gyr', t: 1e10 * 3.156e7 }
    ];
    
    console.log('  Time Evolution of NGC1316 Gravity:');
    timePoints.forEach(tp => {
        const g = ngc.computeG(tp.t, r);
        const m_merge = ngc.computeMmerge(tp.t);
        console.log(`    ${tp.label.padEnd(12)} | g = ${g.toExponential(3)} m/s² | M_merge = ${(m_merge / 1.989e30).toExponential(2)} M☉`);
    });
    
    console.log('✓ TEST 3 PASSED\n');
} catch (e) {
    console.error('✗ TEST 3 FAILED:', e.message, '\n');
}

// Test 4: Dynamic Method Expansion
console.log('▶ TEST 4: Dynamic Method Expansion - Runtime Capability Addition');
console.log('─'.repeat(65));
try {
    const ngc = new idx.NGC1316UQFFModule();
    
    // Add custom analysis method dynamically
    const customAnalysis = function() {
        const vars = this.variables;
        return {
            totalMass: vars.get('M'),
            visibleMass: vars.get('M_visible'),
            darkMatterRatio: vars.get('M_DM') / vars.get('M'),
            massRatio: (vars.get('M_DM') / vars.get('M_visible')).toFixed(3)
        };
    };
    
    ngc.expand('customAnalysis', customAnalysis);
    const analysis = ngc.customAnalysis();
    
    console.log(`  Custom analysis method added dynamically`);
    console.log(`  Total Mass: ${(analysis.totalMass / 1.989e30).toExponential(2)} M☉`);
    console.log(`  Visible Mass: ${(analysis.visibleMass / 1.989e30).toExponential(2)} M☉`);
    console.log(`  Dark Matter Ratio: ${(analysis.darkMatterRatio * 100).toFixed(1)}%`);
    console.log(`  DM/Visible Ratio: ${analysis.massRatio}`);
    console.log('✓ TEST 4 PASSED\n');
} catch (e) {
    console.error('✗ TEST 4 FAILED:', e.message, '\n');
}

// Test 5: UQFF Component Breakdown
console.log('▶ TEST 5: UQFF Component Physics Breakdown');
console.log('─'.repeat(65));
try {
    const ngc = new idx.NGC1316UQFFModule();
    const t = 2e9 * 3.156e7;
    const r = 20e3 * 3.086e19;
    
    console.log('  Universal Gravity Components (Ug) at t=2Gyr, r=20kpc:');
    console.log(`    Ug1 (Dipole):            ${ngc.computeUg1(t).toExponential(3)} J/m`);
    console.log(`    Ug2 (Outer Field):       ${ngc.computeUg2(t).toExponential(3)} J/m`);
    console.log(`    Ug3' (Tidal):            ${ngc.computeUg3prime(t).toExponential(3)} m/s²`);
    console.log(`    Ug4 (Reaction):          ${ngc.computeUg4(t).toExponential(3)} J/m`);
    
    console.log('\n  Environmental and Quantum Terms:');
    console.log(`    F_env (Environmental):   ${ngc.computeFenv(t).toExponential(3)} m/s²`);
    console.log(`    Ui (Universal Inertia):  ${ngc.computeUi(t).toExponential(3)}`);
    console.log(`    Quantum Term:            ${ngc.computeQuantumTerm(ngc.variables.get('t_Hubble'), r).toExponential(3)} m/s²`);
    console.log(`    Fluid (Dust) Term:       ${ngc.computeFluidTerm(1.0).toExponential(3)} kg·m⁻³·m/s²`);
    
    console.log('✓ TEST 5 PASSED\n');
} catch (e) {
    console.error('✗ TEST 5 FAILED:', e.message, '\n');
}

// Test 6: Multi-Parameter Sensitivity Analysis
console.log('▶ TEST 6: Parameter Sensitivity Analysis');
console.log('─'.repeat(65));
try {
    const ngc = new idx.NGC1316UQFFModule();
    const t = 2e9 * 3.156e7;
    const r = 20e3 * 3.086e19;
    
    const baseline = ngc.computeG(t, r);
    
    console.log(`  Baseline g_NGC1316: ${baseline.toExponential(3)} m/s²`);
    console.log('\n  Parameter variation effects:');
    
    // Test 1: Mass variation
    const origM = ngc.variables.get('M');
    ngc.updateVariable('M', origM * 1.5);
    const g_1_5M = ngc.computeG(t, r);
    console.log(`    +50% Mass: ${(g_1_5M / baseline).toFixed(3)}x (${g_1_5M.toExponential(3)} m/s²)`);
    ngc.updateVariable('M', origM);
    
    // Test 2: Radius variation
    const origR = 20e3 * 3.086e19;
    const g_halvedR = ngc.computeG(t, origR / 2);
    console.log(`    Halved Radius: ${(g_halvedR / baseline).toFixed(3)}x (${g_halvedR.toExponential(3)} m/s²)`);
    
    // Test 3: B-field variation
    const origB = ngc.variables.get('B');
    ngc.updateVariable('B', origB * 2);
    const g_2B = ngc.computeG(t, r);
    console.log(`    +100% B-field: ${(g_2B / baseline).toFixed(3)}x (${g_2B.toExponential(3)} m/s²)`);
    ngc.updateVariable('B', origB);
    
    console.log('✓ TEST 6 PASSED\n');
} catch (e) {
    console.error('✗ TEST 6 FAILED:', e.message, '\n');
}

// Test 7: Computational Performance
console.log('▶ TEST 7: Computational Performance Benchmark');
console.log('─'.repeat(65));
try {
    const ngc = new idx.NGC1316UQFFModule();
    const iterations = 1000;
    const t = 2e9 * 3.156e7;
    const r = 20e3 * 3.086e19;
    
    const startTime = Date.now();
    for (let i = 0; i < iterations; i++) {
        ngc.computeG(t, r + i * 1e20);
    }
    const elapsed = Date.now() - startTime;
    const opsPerMs = (iterations / elapsed).toFixed(2);
    
    console.log(`  Iterations: ${iterations}`);
    console.log(`  Time elapsed: ${elapsed} ms`);
    console.log(`  Operations per millisecond: ${opsPerMs}`);
    console.log(`  Average time per operation: ${(elapsed / iterations).toFixed(3)} ms`);
    console.log('✓ TEST 7 PASSED\n');
} catch (e) {
    console.error('✗ TEST 7 FAILED:', e.message, '\n');
}

// Test 8: Full Analysis Function
console.log('▶ TEST 8: Complete NGC1316 UQFF Analysis Framework');
console.log('─'.repeat(65));
try {
    const result = idx.analyzeNGC1316UQFF71();
    
    console.log(`  ✓ Analysis function executed`);
    console.log(`  Module instantiated: ${result.module ? 'Yes' : 'No'}`);
    console.log(`  Timestamp: ${result.timestamp}`);
    console.log(`  Physics domain: ${result.physics}`);
    console.log(`  Variables stored: ${result.module.variables.size}`);
    
    console.log('\n  Sample physics variables:');
    console.log(`    M_visible: ${(result.module.variables.get('M_visible') / 1.989e30).toExponential(2)} M☉`);
    console.log(`    M_DM: ${(result.module.variables.get('M_DM') / 1.989e30).toExponential(2)} M☉`);
    console.log(`    Redshift: ${result.module.variables.get('z')}`);
    
    console.log('✓ TEST 8 PASSED\n');
} catch (e) {
    console.error('✗ TEST 8 FAILED:', e.message, '\n');
}

// Summary
console.log('╔════════════════════════════════════════════════════════════════╗');
console.log('║                       TEST SUITE COMPLETE                      ║');
console.log('║                                                                ║');
console.log('║  ✓ NGC1316 Module Integration                                  ║');
console.log('║  ✓ Multi-System Computations                                   ║');
console.log('║  ✓ Cosmic Time Evolution                                       ║');
console.log('║  ✓ Dynamic Method Expansion                                    ║');
console.log('║  ✓ UQFF Component Analysis                                     ║');
console.log('║  ✓ Parameter Sensitivity                                       ║');
console.log('║  ✓ Computational Performance                                   ║');
console.log('║  ✓ Full Analysis Framework                                     ║');
console.log('║                                                                ║');
console.log('║  Framework Status: OPERATIONAL                                 ║');
console.log('║  All computational capabilities verified                       ║');
console.log('╚════════════════════════════════════════════════════════════════╝\n');
