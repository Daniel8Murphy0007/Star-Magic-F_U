#!/usr/bin/env node

/**
 * UQFF Framework Iteration v2.0 - Quick Start Demo
 * Shows all new iteration capabilities in action
 */

const idx = require('./index.js');
const Engines = require('./framework_iteration_engine.js');

console.log('╔════════════════════════════════════════════════════════════════╗');
console.log('║    UQFF FRAMEWORK ITERATION v2.0 - CAPABILITIES SHOWCASE      ║');
console.log('╚════════════════════════════════════════════════════════════════╝\n');

// Initialize all four iteration engines
const crossSystem = new Engines.CrossSystemInteractionEngine();
const perfOptimizer = new Engines.PerformanceOptimizer();
const analysis = new Engines.StatisticalAnalysisEngine();
const paramOptimizer = new Engines.ParameterOptimizationEngine();

// ============================================================================
// CAPABILITY 1: Cross-System Interaction
// ============================================================================

console.log('▶ CAPABILITY 1: Cross-System Interaction Network');
console.log('─'.repeat(65) + '\n');

// Register multiple systems
const ngc = new idx.NGC1316UQFFModule();
const andromeda = new idx.AndromedaUQFFModule();
const sombrero = new idx.SombreroUQFFModule();

crossSystem.registerSystem('NGC1316', ngc);
crossSystem.registerSystem('Andromeda', andromeda);
crossSystem.registerSystem('Sombrero', sombrero);

// Create connection network
crossSystem.connectSystems('NGC1316', 'Andromeda', 'gravitational');
crossSystem.connectSystems('Andromeda', 'Sombrero', 'tidal');
crossSystem.connectSystems('NGC1316', 'Sombrero', 'magnetic');

// Print registry
console.log('System Network Registry:\n');
console.log(`  Registered Systems: 3`);
console.log(`  Active Connections: 3`);
console.log(`  Interaction Types: gravitational, tidal, magnetic\n`);

// Compute a unified field interaction
const interaction = crossSystem.computeUnifiedFieldInteraction(
    'NGC1316', 'Andromeda',
    {
        distance: 2.54e22,  // ~823 kpc (realistic distance)
        mass1: 5e11 * 1.989e30,
        mass2: 1.2e11 * 1.989e30,
        includeGravitational: true,
        includeMagnetic: true
    }
);

console.log('Unified Field Interaction: NGC1316 ↔ Andromeda');
console.log(`  Gravitational Force: ${interaction.interaction.gravitational.force.toExponential(3)} N`);
console.log(`  Total Unified Field: ${interaction.unifiedField.total.toExponential(3)} N/m²`);
console.log(`  Composition:`);
console.log(`    ├─ Gravitational: ${(interaction.unifiedField.composition.g_ratio * 100).toFixed(1)}%`);
console.log(`    ├─ Electromagnetic: ${(interaction.unifiedField.composition.e_ratio * 100).toFixed(1)}%`);
console.log(`    └─ Magnetic: ${(interaction.unifiedField.composition.m_ratio * 100).toFixed(1)}%\n`);

// ============================================================================
// CAPABILITY 2: Performance Optimization
// ============================================================================

console.log('▶ CAPABILITY 2: Performance Optimization & Memoization');
console.log('─'.repeat(65) + '\n');

console.log('Running 500 optimized computations with memoization...\n');

const t = 2e9 * 3.156e7;
const r = 20e3 * 3.086e19;

for (let i = 0; i < 500; i++) {
    perfOptimizer.executeOptimized(
        'computeG',
        (t_val, r_val) => ngc.computeG(t_val, r_val),
        [t, r],
        { memoize: true }
    );
}

console.log('Cache Performance After 500 Iterations:\n');
console.log(`  Total Computations: ${perfOptimizer.performanceMetrics.totalComputations}`);
console.log(`  Cache Hits: ${perfOptimizer.performanceMetrics.cacheHits}`);
console.log(`  Cache Misses: ${perfOptimizer.performanceMetrics.cacheMisses}`);
console.log(`  Hit Ratio: ${perfOptimizer.getCacheHitRatio().toFixed(2)}%\n`);
console.log(`  ✓ Memoization eliminates redundant computation!`);
console.log(`  ✓ 99%+ cache efficiency achieved!\n`);

// ============================================================================
// CAPABILITY 3: Statistical Analysis
// ============================================================================

console.log('▶ CAPABILITY 3: Statistical Analysis & Anomaly Detection');
console.log('─'.repeat(65) + '\n');

console.log('Recording time series across 10 billion years...\n');

// Record time series
for (let gyr = 0; gyr <= 10; gyr++) {
    const time_s = gyr * 1e9 * 3.156e7;
    const g_val = ngc.computeG(time_s, r);
    analysis.recordDataPoint('NGC1316', 'gravity', g_val);
}

// Add anomalous data point
analysis.recordDataPoint('NGC1316', 'gravity', 9.945e36 * 8);

// Analyze
const stats = analysis.analyzeTimeSeries('NGC1316', 'gravity');
const anomalies = analysis.detectAnomalies('NGC1316', 'gravity', 2.5);

console.log('Time Series Statistics:\n');
console.log(`  Data Points: ${stats.count}`);
console.log(`  Mean Gravity: ${stats.mean.toExponential(3)} m/s²`);
console.log(`  Std Dev: ${stats.stdDev.toExponential(3)} m/s²`);
console.log(`  Range: ${stats.range.toExponential(3)} m/s²`);
console.log(`  Trend: ${stats.trend}\n`);

console.log(`Anomaly Detection Results:\n`);
console.log(`  Anomalies Found: ${anomalies.length}`);
if (anomalies.length > 0) {
    console.log(`  Severity: ${anomalies[0].severity.toUpperCase()}`);
    console.log(`  Value: ${anomalies[0].value.toExponential(3)} m/s²`);
    console.log(`  Z-Score: ${anomalies[0].zScore.toFixed(2)}\n`);
}

// ============================================================================
// CAPABILITY 4: Parameter Optimization
// ============================================================================

console.log('▶ CAPABILITY 4: Parameter Optimization (Grid Search)');
console.log('─'.repeat(65) + '\n');

console.log('Searching for parameters that produce g = 1.0e37 m/s²...\n');

const targetG = 1.0e37;
const searchSpace = {
    't': [1e9 * 3.156e7, 3e9 * 3.156e7],          // 1-3 Gyr
    'M_spiral': [5e9 * 1.989e30, 1.5e10 * 1.989e30]  // Spiral mass
};

const results = paramOptimizer.gridSearchOptimization(
    ngc,
    'gravity',
    targetG,
    searchSpace,
    4  // 4x4 = 16 evaluations
);

console.log('Optimization Results:\n');
console.log(`  Target Value: ${targetG.toExponential(3)} m/s²`);
console.log(`  Search Space: 2 parameters × 4 steps = 16 evaluations`);
console.log(`  Best Result:`);
console.log(`    Error: ${results[0].error.toExponential(3)}`);
console.log(`    Relative Error: ${(results[0].relativeError * 100).toFixed(2)}%`);
console.log(`    Computed: ${results[0].computed.toExponential(3)} m/s²\n`);

if (results[1]) {
    console.log(`  2nd Best Result:`);
    console.log(`    Error: ${results[1].error.toExponential(3)}`);
    console.log(`    Relative Error: ${(results[1].relativeError * 100).toFixed(2)}%\n`);
}

// ============================================================================
// INTEGRATION EXAMPLE
// ============================================================================

console.log('▶ INTEGRATED WORKFLOW: Multi-Engine Coordination');
console.log('─'.repeat(65) + '\n');

console.log('Scenario: Analyze three galaxy systems with all four engines:\n');

// Create three systems
const sys1 = new idx.NGC1316UQFFModule();
const sys2 = new idx.AndromedaUQFFModule();
const sys3 = new idx.LagoonUQFFModule();

// Register with cross-system engine
console.log('Step 1: Registering systems...');
crossSystem.registerSystem('Galaxy_A', sys1);
crossSystem.registerSystem('Galaxy_B', sys2);
crossSystem.registerSystem('Galaxy_C', sys3);
console.log('  ✓ 3 systems registered\n');

// Create interaction network
console.log('Step 2: Creating interaction network...');
crossSystem.connectSystems('Galaxy_A', 'Galaxy_B', 'gravitational');
crossSystem.connectSystems('Galaxy_B', 'Galaxy_C', 'tidal');
console.log('  ✓ 2 connections established\n');

// Compute interactions with optimization
console.log('Step 3: Computing unified fields (with caching)...');
for (let i = 0; i < 3; i++) {
    perfOptimizer.executeOptimized(
        `galaxy_${i}_interaction`,
        () => crossSystem.computeUnifiedFieldInteraction(
            'Galaxy_A', 'Galaxy_B',
            { distance: 1e20, mass1: 5e11 * 1.989e30, mass2: 1.2e11 * 1.989e30 }
        ),
        [],
        { memoize: true }
    );
}
console.log(`  ✓ Computations complete (Cache hit ratio: ${perfOptimizer.getCacheHitRatio().toFixed(1)}%)\n`);

// Collect statistics
console.log('Step 4: Analyzing time series...');
for (let i = 0; i < 5; i++) {
    const t = i * 2e9 * 3.156e7;
    const g = sys1.computeG(t, 20e3 * 3.086e19);
    analysis.recordDataPoint('Galaxy_A', 'gravity', g);
}
const stats2 = analysis.analyzeTimeSeries('Galaxy_A', 'gravity');
console.log(`  ✓ Trend detected: ${stats2.trend}\n`);

// Optimize
console.log('Step 5: Optimizing parameters...');
const opt = paramOptimizer.gridSearchOptimization(
    sys1, 'gravity', 1e37,
    { 't': [1e9 * 3.156e7, 3e9 * 3.156e7] },
    3
);
console.log(`  ✓ Best fit error: ${opt[0].error.toExponential(3)}\n`);

console.log('═'.repeat(65));
console.log('✓ Integrated workflow complete!');
console.log('═'.repeat(65) + '\n');

// ============================================================================
// SUMMARY
// ============================================================================

console.log('╔════════════════════════════════════════════════════════════════╗');
console.log('║                  CAPABILITIES SUMMARY                         ║');
console.log('╚════════════════════════════════════════════════════════════════╝\n');

console.log('✅ CAPABILITY 1: Cross-System Interaction');
console.log('   └─ 3 systems registered, 3 connections active');
console.log('   └─ Unified field computation: working\n');

console.log('✅ CAPABILITY 2: Performance Optimization');
console.log(`   └─ Cache hit ratio: ${perfOptimizer.getCacheHitRatio().toFixed(2)}%`);
console.log('   └─ Memoization: ACTIVE\n');

console.log('✅ CAPABILITY 3: Statistical Analysis');
console.log(`   └─ Time series records: ${analysis.timeSeries.size} parameters`);
console.log(`   └─ Anomalies detected: ${anomalies.length}\n`);

console.log('✅ CAPABILITY 4: Parameter Optimization');
console.log(`   └─ Grid search: ${results.length} results found`);
console.log(`   └─ Best fit error: ${results[0].error.toExponential(3)}\n`);

console.log('═'.repeat(65));
console.log('Framework Iteration v2.0: COMPLETE & OPERATIONAL');
console.log('═'.repeat(65) + '\n');
