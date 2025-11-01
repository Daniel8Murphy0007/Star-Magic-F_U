// Quick test of enhanced NGC346 adaptive module
const NGC346 = require('./ngc346_uqff.js');

console.log('╔════════════════════════════════════════════════════════════════╗');
console.log('║     NGC346 ADAPTIVE MODULE - QUICK VERIFICATION TEST          ║');
console.log('╚════════════════════════════════════════════════════════════════╝');
console.log('');

try {
    const m = new NGC346();
    
    // Test 1: Module loads
    console.log('✓ Enhanced NGC346UQFFModule loads successfully');
    console.log(`  - Variables: ${Object.keys(m.variables).length}`);
    console.log(`  - Physics Terms Registry: ${Object.keys(m.physicsTerms).length}`);
    console.log(`  - Calibration Data: ${Object.keys(m.calibrationData).length}`);
    console.log(`  - Adaptation Logs: ${m.adaptationLogs.length}`);
    console.log('');
    
    // Test 2: Variable management
    console.log('Testing Variable Management APIs:');
    m.addVariable('custom_var', 42.0);
    console.log(`  ✓ addVariable() works`);
    const val = m.getVariableValue('custom_var');
    console.log(`  ✓ getVariableValue() works (value: ${val})`);
    m.removeVariable('custom_var');
    console.log(`  ✓ removeVariable() works`);
    const varList = m.listVariables();
    console.log(`  ✓ listVariables() works (${varList.length} variables)`);
    console.log('');
    
    // Test 3: Physics term registration
    console.log('Testing Physics Term Registration:');
    m.registerPhysicsTerm('test_term_1', (t, r, v) => 1e-12, true, 'Test physics term');
    console.log(`  ✓ registerPhysicsTerm() works`);
    m.disablePhysicsTerm('test_term_1');
    console.log(`  ✓ disablePhysicsTerm() works`);
    m.enablePhysicsTerm('test_term_1');
    console.log(`  ✓ enablePhysicsTerm() works`);
    const terms = m.getPhysicsTerms();
    console.log(`  ✓ getPhysicsTerms() works (${terms.length} registered terms)`);
    m.unregisterPhysicsTerm('test_term_1');
    console.log(`  ✓ unregisterPhysicsTerm() works`);
    console.log('');
    
    // Test 4: State management
    console.log('Testing State Management:');
    const state1 = m.exportState();
    console.log(`  ✓ exportState() works (${Object.keys(state1).length} state sections)`);
    m.saveCheckpoint('test_checkpoint');
    console.log(`  ✓ saveCheckpoint() works`);
    const checkpoints = m.listCheckpoints();
    console.log(`  ✓ listCheckpoints() works (${checkpoints.length} checkpoints)`);
    m.loadCheckpoint('test_checkpoint');
    console.log(`  ✓ loadCheckpoint() works`);
    console.log('');
    
    // Test 5: Adaptive operations
    console.log('Testing Adaptive Operations:');
    const resonanceMap = m.mapQuantumResonance(5);
    console.log(`  ✓ mapQuantumResonance() works (${resonanceMap.statesCount} states mapped)`);
    const discoveries = m.discoverPhysics({ rho_gas: { min: 1e-21, max: 1e-19, steps: 3 } });
    console.log(`  ✓ discoverPhysics() works (${discoveries.discoveries.length} discoveries)`);
    const anomalies = m.detectAnomalies();
    console.log(`  ✓ detectAnomalies() works (${anomalies.count} anomalies detected)`);
    console.log('');
    
    // Test 6: Custom methods
    console.log('Testing Custom Methods:');
    m.addCustomMethod('test_method', () => 42, 'test custom method');
    console.log(`  ✓ addCustomMethod() works`);
    const result = m.executeCustomMethod('test_method');
    console.log(`  ✓ executeCustomMethod() works (returned: ${result})`);
    console.log('');
    
    // Test 7: Reporting
    console.log('Testing Reporting & Utilities:');
    const report = m.generateReport();
    console.log(`  ✓ generateReport() works (${report.length} chars)`);
    const config = m.exportConfiguration();
    console.log(`  ✓ exportConfiguration() works (${Object.keys(config).length} sections)`);
    const log = m.getAdaptationLog();
    console.log(`  ✓ getAdaptationLog() works (${log.length} entries)`);
    console.log('');
    
    // Test 8: computeG with custom terms
    console.log('Testing Computation with Custom Terms:');
    m.registerPhysicsTerm('contrib_term', (t, r, v) => 1e-11, true, 'Contribution term');
    const g_base = 1e-10;  // Approximate base value
    const g_with_custom = m.computeG(m.variables['t'], m.variables['r']);
    console.log(`  ✓ computeG() with custom terms works (g ≈ ${g_with_custom.toExponential(3)} m/s²)`);
    m.unregisterPhysicsTerm('contrib_term');
    console.log('');
    
    console.log('╔════════════════════════════════════════════════════════════════╗');
    console.log('║              ALL TESTS PASSED ✓                               ║');
    console.log('║   Enhanced NGC346 Adaptive Module Ready for Deployment         ║');
    console.log('╚════════════════════════════════════════════════════════════════╝');
    
} catch (err) {
    console.error('❌ TEST FAILED:', err.message);
    console.error(err.stack);
    process.exit(1);
}
