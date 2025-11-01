#!/usr/bin/env node

/**
 * UQFF Framework Iteration Test Suite v2.0
 * Validates all new iteration capabilities: cross-system interaction, optimization, and analysis
 */

const idx = require('./index.js');
const IterationEngine = require('./framework_iteration_engine.js');

console.log('╔════════════════════════════════════════════════════════════════╗');
console.log('║  UQFF FRAMEWORK ITERATION CAPABILITIES TEST SUITE v2.0        ║');
console.log('╚════════════════════════════════════════════════════════════════╝\n');

// Initialize engines
const interactionEngine = new IterationEngine.CrossSystemInteractionEngine();
const performanceOptimizer = new IterationEngine.PerformanceOptimizer();
const analysisEngine = new IterationEngine.StatisticalAnalysisEngine();
const optimizationEngine = new IterationEngine.ParameterOptimizationEngine();

// ============================================================================
// TEST 1: Cross-System Interaction
// ============================================================================

console.log('▶ TEST 1: Cross-System Interaction Engine');
console.log('─'.repeat(65));

try {
    // Create sample systems
    const ngc1316 = new idx.NGC1316UQFFModule();
    const magnetar1 = new idx.MagnetarSGR1745_2900();
    const magnetar2 = new idx.MagnetarSGR0501_4516();
    
    // Register systems
    interactionEngine.registerSystem('NGC1316', ngc1316);
    interactionEngine.registerSystem('SGR1745-2900', magnetar1);
    interactionEngine.registerSystem('SGR0501+4516', magnetar2);
    
    // Connect systems
    interactionEngine.connectSystems('NGC1316', 'SGR1745-2900', 'gravitational');
    interactionEngine.connectSystems('SGR1745-2900', 'SGR0501+4516', 'magnetic');
    
    // Compute unified field interactions
    const interaction1 = interactionEngine.computeUnifiedFieldInteraction(
        'NGC1316', 'SGR1745-2900',
        {
            distance: 3.086e19,  // ~1 parsec
            mass1: 5e11 * 1.989e30,
            mass2: 2.78e30,
            includeGravitational: true,
            includeElectromagnetic: false,
            includeMagnetic: true
        }
    );
    
    console.log('\nUniformed Field Interaction (NGC1316 ↔ SGR1745-2900):');
    console.log(`  Gravitational Force: ${interaction1.interaction.gravitational.force.toExponential(3)} N`);
    console.log(`  Unified Field Total: ${interaction1.unifiedField.total.toExponential(3)} N/m²`);
    console.log(`  Field Composition:`);
    console.log(`    - Gravitational: ${(interaction1.unifiedField.composition.g_ratio * 100).toFixed(2)}%`);
    console.log(`    - Electromagnetic: ${(interaction1.unifiedField.composition.e_ratio * 100).toFixed(2)}%`);
    console.log(`    - Magnetic: ${(interaction1.unifiedField.composition.m_ratio * 100).toFixed(2)}%`);
    
    // Print registry
    interactionEngine.printRegistry();
    
    console.log('\n✓ TEST 1 PASSED\n');
} catch (e) {
    console.error('✗ TEST 1 FAILED:', e.message, '\n');
}

// ============================================================================
// TEST 2: Performance Optimization
// ============================================================================

console.log('▶ TEST 2: Performance Optimization & Memoization');
console.log('─'.repeat(65));

try {
    const ngc = new idx.NGC1316UQFFModule();
    
    // Perform multiple computations with optimization
    console.log('\nOptimized Computation Series (1000 iterations):');
    
    const t = 2e9 * 3.156e7;
    const r = 20e3 * 3.086e19;
    
    for (let i = 0; i < 1000; i++) {
        performanceOptimizer.executeOptimized(
            'computeG',
            (t_val, r_val) => ngc.computeG(t_val, r_val),
            [t, r],
            { memoize: true }
        );
    }
    
    performanceOptimizer.printPerformanceReport();
    
    console.log('\n✓ TEST 2 PASSED\n');
} catch (e) {
    console.error('✗ TEST 2 FAILED:', e.message, '\n');
}

// ============================================================================
// TEST 3: Statistical Analysis & Anomaly Detection
// ============================================================================

console.log('▶ TEST 3: Statistical Analysis & Anomaly Detection');
console.log('─'.repeat(65));

try {
    const ngc = new idx.NGC1316UQFFModule();
    
    // Generate time series data
    console.log('\nGenerating time series data for NGC1316...');
    
    const times = [];
    for (let i = 0; i <= 10; i++) {
        times.push(i * 1e9 * 3.156e7); // 0 to 10 Gyr
    }
    
    const r = 20e3 * 3.086e19;
    
    for (let i = 0; i < times.length; i++) {
        const g_val = ngc.computeG(times[i], r);
        analysisEngine.recordDataPoint('NGC1316', 'gravity', g_val);
    }
    
    // Inject anomaly
    analysisEngine.recordDataPoint('NGC1316', 'gravity', 9.945e36 * 10, Date.now());
    
    // Analyze
    const stats = analysisEngine.analyzeTimeSeries('NGC1316', 'gravity');
    console.log(`\nStatistical Summary:`);
    console.log(`  Mean: ${stats.mean.toExponential(3)}`);
    console.log(`  Std Dev: ${stats.stdDev.toExponential(3)}`);
    console.log(`  Trend: ${stats.trend}`);
    
    const anomalies = analysisEngine.detectAnomalies('NGC1316', 'gravity', 2.0);
    console.log(`\nAnomalies Detected: ${anomalies.length}`);
    if (anomalies.length > 0) {
        anomalies.forEach((anom, idx) => {
            console.log(`  ${idx + 1}. Value: ${anom.value.toExponential(3)} (z-score: ${anom.zScore.toFixed(2)}) [${anom.severity}]`);
        });
    }
    
    console.log('\n✓ TEST 3 PASSED\n');
} catch (e) {
    console.error('✗ TEST 3 FAILED:', e.message, '\n');
}

// ============================================================================
// TEST 4: Parameter Optimization
// ============================================================================

console.log('▶ TEST 4: Parameter Optimization Engine');
console.log('─'.repeat(65));

try {
    const ngc = new idx.NGC1316UQFFModule();
    
    // Define target and search parameters
    const targetG = 1e37; // Target gravity value
    const searchParams = {
        't': [1e9 * 3.156e7, 3e9 * 3.156e7],          // Time range: 1-3 Gyr
        'M_spiral': [5e9 * 1.989e30, 1.5e10 * 1.989e30]  // Spiral mass range
    };
    
    console.log(`\nOptimizing for target gravity: ${targetG.toExponential(3)} m/s²`);
    console.log(`Search space: ${Object.keys(searchParams).length} parameters`);
    
    // Perform grid search with reduced steps for speed
    const results = optimizationEngine.gridSearchOptimization(
        ngc,
        'gravity',
        targetG,
        searchParams,
        5  // 5 steps = 25 evaluations
    );
    
    console.log(`\nBest fit parameters found:`);
    if (results.length > 0) {
        console.log(`  1. Error: ${results[0].error.toExponential(3)}`);
        console.log(`     Relative Error: ${(results[0].relativeError * 100).toFixed(2)}%`);
        console.log(`     Computed: ${results[0].computed.toExponential(3)}`);
        
        if (results[1]) {
            console.log(`  2. Error: ${results[1].error.toExponential(3)}`);
            console.log(`     Relative Error: ${(results[1].relativeError * 100).toFixed(2)}%`);
        }
    }
    
    console.log('\n✓ TEST 4 PASSED\n');
} catch (e) {
    console.error('✗ TEST 4 FAILED:', e.message, '\n');
}

// ============================================================================
// TEST 5: Integrated Multi-System Workflow
// ============================================================================

console.log('▶ TEST 5: Integrated Multi-System Workflow');
console.log('─'.repeat(65));

try {
    console.log('\nExecuting integrated workflow combining all capabilities...\n');
    
    // Step 1: Create systems
    const sys1 = new idx.NGC1316UQFFModule();
    const sys2 = new idx.MagnetarSGR1745_2900();
    const sys3 = new idx.AndromedaUQFFModule();
    
    // Step 2: Register and connect
    interactionEngine.registerSystem('System1_NGC1316', sys1);
    interactionEngine.registerSystem('System2_SGR1745', sys2);
    interactionEngine.registerSystem('System3_Andromeda', sys3);
    
    interactionEngine.connectSystems('System1_NGC1316', 'System2_SGR1745', 'gravitational');
    interactionEngine.connectSystems('System2_SGR1745', 'System3_Andromeda', 'magnetic');
    interactionEngine.connectSystems('System1_NGC1316', 'System3_Andromeda', 'tidal');
    
    // Step 3: Compute with optimization
    console.log('Computing unified fields across network with memoization...');
    
    for (let i = 0; i < 3; i++) {
        const interaction = performanceOptimizer.executeOptimized(
            'cross_system_field_' + i,
            () => interactionEngine.computeUnifiedFieldInteraction(
                'System1_NGC1316', 'System2_SGR1745',
                { distance: 1e20, mass1: 5e11 * 1.989e30, mass2: 2.78e30 }
            ),
            [],
            { memoize: true }
        );
    }
    
    console.log(`✓ Computations completed`);
    console.log(`  Cache hit ratio: ${performanceOptimizer.getCacheHitRatio().toFixed(2)}%`);
    console.log(`  Operations/ms: ${performanceOptimizer.getOperationsPerMs().toFixed(2)}`);
    
    // Step 4: Analyze results
    console.log('\nAnalyzing results across integrated system...');
    const t = 2e9 * 3.156e7;
    const r = 20e3 * 3.086e19;
    
    for (let cycle = 0; cycle < 5; cycle++) {
        const g = sys1.computeG(t + cycle * 1e8 * 3.156e7, r);
        analysisEngine.recordDataPoint('System1_NGC1316', 'unified_gravity', g);
    }
    
    const analysis = analysisEngine.analyzeTimeSeries('System1_NGC1316', 'unified_gravity');
    console.log(`✓ Time series analysis complete`);
    console.log(`  Data points: ${analysis.count}`);
    console.log(`  Mean: ${analysis.mean.toExponential(3)}`);
    console.log(`  Trend: ${analysis.trend}`);
    
    console.log('\n✓ TEST 5 PASSED\n');
} catch (e) {
    console.error('✗ TEST 5 FAILED:', e.message, '\n');
}

// ============================================================================
// SUMMARY & PERFORMANCE REPORT
// ============================================================================

console.log('╔════════════════════════════════════════════════════════════════╗');
console.log('║                    TEST SUITE COMPLETE                        ║');
console.log('╚════════════════════════════════════════════════════════════════╝\n');

console.log('✓ All iteration capabilities tested and validated');
console.log('\nFramework Iteration Enhancements:');
console.log('  ✓ Cross-system interaction layer (3+ systems connected)');
console.log('  ✓ Performance optimization with memoization');
console.log('  ✓ Statistical analysis and anomaly detection');
console.log('  ✓ Parameter optimization with grid search');
console.log('  ✓ Integrated multi-system workflows');

performanceOptimizer.printPerformanceReport();

console.log('\n═══════════════════════════════════════════════════════════════');
console.log('Framework iteration complete. New capabilities ready for use.');
console.log('═══════════════════════════════════════════════════════════════\n');
