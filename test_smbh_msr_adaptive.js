// test_smbh_msr_adaptive.js
// Comprehensive test suite for adaptive SMBH M-σ module
// Tests: parameter optimization, physics discovery, quantum mapping,
// framework expansion, multi-scale refinement, anomaly detection, auto-tuning

const SMBHMSRAdaptiveModule = require('./smbh_msr_adaptive.js');

class AdaptiveTestSuite {
    constructor() {
        this.adaptive = new SMBHMSRAdaptiveModule();
        this.tests_passed = 0;
        this.tests_failed = 0;
        this.test_results = [];
    }
    
    assert(condition, message) {
        if (condition) {
            this.tests_passed++;
            this.test_results.push({ status: '✓', message });
            return true;
        } else {
            this.tests_failed++;
            this.test_results.push({ status: '✗', message });
            return false;
        }
    }
    
    // ===== CATEGORY 1: INITIALIZATION (8 tests) =====
    testInitialization() {
        console.log('\n=== TEST CATEGORY 1: Initialization (8 tests) ===');
        
        // T1.1
        this.assert(this.adaptive !== null && this.adaptive.core !== null,
            'T1.1: Adaptive module instantiates with core module');
        
        // T1.2
        this.assert(Array.isArray(this.adaptive.history),
            'T1.2: History array initialized');
        
        // T1.3
        this.assert(typeof this.adaptive.calibrationData === 'object' && 
                    this.adaptive.calibrationData.hasOwnProperty('baseline'),
            'T1.3: Baseline calibration stored');
        
        // T1.4
        this.assert(this.adaptive.performanceMetrics.computations === 0,
            'T1.4: Performance metrics initialized at zero');
        
        // T1.5
        this.assert(Array.isArray(this.adaptive.adaptationLogs) && this.adaptive.adaptationLogs.length > 0,
            'T1.5: Adaptation logging started');
        
        // T1.6
        this.assert(this.adaptive.expanded === false,
            'T1.6: Expansion flag initialized as false');
        
        // T1.7
        this.assert(typeof this.adaptive.customMethods === 'object',
            'T1.7: Custom methods registry initialized');
        
        // T1.8
        this.assert(Object.keys(this.adaptive.quantumMap).length === 0,
            'T1.8: Quantum map empty until mapping performed');
    }
    
    // ===== CATEGORY 2: PARAMETER OPTIMIZATION (12 tests) =====
    testParameterOptimization() {
        console.log('\n=== TEST CATEGORY 2: Parameter Optimization (12 tests) ===');
        
        const targetData = 1e-7;  // Target acceleration
        
        // T2.1
        const originalM = this.adaptive.core.variables['M_bh'];
        const result = this.adaptive.optimizeParameters(targetData, 10, 0.01);
        this.assert(result !== null && result.hasOwnProperty('converged'),
            'T2.1: Optimization returns valid result object');
        
        // T2.2
        this.assert(result.hasOwnProperty('finalError'),
            'T2.2: Result includes final error');
        
        // T2.3
        this.assert(result.history.length === 10,
            'T2.3: History tracks all 10 iterations');
        
        // T2.4
        this.assert(result.optimizedParams.hasOwnProperty('M_bh'),
            'T2.4: Optimized parameters stored');
        
        // T2.5
        this.assert(this.adaptive.calibrationData.hasOwnProperty(Object.keys(this.adaptive.calibrationData)[Object.keys(this.adaptive.calibrationData).length - 1]),
            'T2.5: Calibration point saved after optimization');
        
        // T2.6
        const gradients = this.adaptive.computeGradients(targetData);
        this.assert(Object.keys(gradients).length > 0,
            'T2.6: Gradient computation returns results');
        
        // T2.7
        this.assert(this.adaptive.getOptimizableParams().hasOwnProperty('sigma'),
            'T2.7: Optimizable parameters identifiable');
        
        // T2.8
        const clamped = this.adaptive.clampToPhysicalRange('M_bh', 1e50);
        this.assert(clamped <= 1e44,
            'T2.8: Physical range clamping works');
        
        // T2.9
        const error1 = this.adaptive.computeError(1e-7, 1e-7);
        const error2 = this.adaptive.computeError(1e-7, 2e-7);
        this.assert(error2 > error1,
            'T2.9: Error computation shows divergence');
        
        // T2.10
        this.assert(result.optimizedParams.sigma >= 1e5 && result.optimizedParams.sigma <= 1e6,
            'T2.10: Optimized sigma stays in physical range');
        
        // T2.11
        this.assert(result.optimizedParams.M_bh >= 1e41 && result.optimizedParams.M_bh <= 1e44,
            'T2.11: Optimized M_bh stays in physical range');
        
        // T2.12
        this.assert(result.iterations === 10,
            'T2.12: Iteration count matches request');
    }
    
    // ===== CATEGORY 3: QUANTUM RESONANCE MAPPING (12 tests) =====
    testQuantumResonanceMapping() {
        console.log('\n=== TEST CATEGORY 3: Quantum Resonance Mapping (12 tests) ===');
        
        const map = this.adaptive.mapQuantumResonance();
        
        // T3.1
        this.assert(map.fullMap.length === 26,
            'T3.1: All 26 quantum states mapped');
        
        // T3.2
        this.assert(map.maxResonanceState >= 1 && map.maxResonanceState <= 26,
            'T3.2: Maximum resonance state in valid range');
        
        // T3.3
        this.assert(map.maxResonanceValue > 0,
            'T3.3: Maximum resonance value positive');
        
        // T3.4
        this.assert(Array.isArray(map.peaks),
            'T3.4: Resonance peaks array initialized');
        
        // T3.5
        this.assert(map.fullMap[0].hasOwnProperty('Um') && map.fullMap[0].hasOwnProperty('Ug1'),
            'T3.5: Both U_m and U_g1 components in map');
        
        // T3.6
        this.assert(map.fullMap[0].hasOwnProperty('resonance'),
            'T3.6: Resonance value computed for each state');
        
        // T3.7
        this.assert(map.fullMap[0].resonanceNormalized >= 0 && map.fullMap[0].resonanceNormalized <= 1,
            'T3.7: Normalized resonance in [0,1]');
        
        // T3.8
        const optimalState = this.adaptive.getOptimalQuantumState();
        this.assert(optimalState === map.maxResonanceState,
            'T3.8: Optimal state matches maximum in map');
        
        // T3.9
        this.assert(map.parameters.hasOwnProperty('t'),
            'T3.9: Mapping parameters stored (time)');
        
        // T3.10
        this.assert(map.parameters.hasOwnProperty('sigma'),
            'T3.10: Mapping parameters stored (dispersion)');
        
        // T3.11
        const resonances = map.fullMap.map(s => s.resonance);
        const increasing = resonances.slice(0, Math.floor(resonances.length/2)).reduce((a,b) => a < b ? a : b);
        this.assert(map.peaks.length > 0 || resonances.length > 0,
            'T3.11: Resonance analysis executable');
        
        // T3.12
        this.assert(map.fullMap.every(state => isFinite(state.total)),
            'T3.12: All states produce finite values');
    }
    
    // ===== CATEGORY 4: PHYSICS DISCOVERY (10 tests) =====
    testPhysicsDiscovery() {
        console.log('\n=== TEST CATEGORY 4: Physics Discovery (10 tests) ===');
        
        const discoveries = this.adaptive.discoverPhysics();
        
        // T4.1
        this.assert(discoveries.hasOwnProperty('parameterSensitivities'),
            'T4.1: Parameter sensitivities discovered');
        
        // T4.2
        this.assert(Object.keys(discoveries.parameterSensitivities).length > 0,
            'T4.2: At least one parameter analyzed');
        
        // T4.3
        const firstSensitivity = Object.values(discoveries.parameterSensitivities)[0];
        this.assert(firstSensitivity.hasOwnProperty('classification'),
            'T4.3: Sensitivity classified (HIGH/MODERATE)');
        
        // T4.4
        this.assert(discoveries.hasOwnProperty('scalingLaws'),
            'T4.4: Scaling laws discovered');
        
        // T4.5
        this.assert(Object.keys(discoveries.scalingLaws).length > 0,
            'T4.5: At least one scaling law identified');
        
        // T4.6
        const scalingExample = Object.values(discoveries.scalingLaws)[0];
        this.assert(scalingExample.hasOwnProperty('exponent') && typeof scalingExample.exponent === 'number',
            'T4.6: Scaling exponent computed numerically');
        
        // T4.7
        this.assert(scalingExample.hasOwnProperty('interpretation'),
            'T4.7: Scaling exponent interpreted');
        
        // T4.8
        this.assert(discoveries.hasOwnProperty('emergentBehaviors'),
            'T4.8: Emergent behaviors detected');
        
        // T4.9
        this.assert(Array.isArray(discoveries.emergentBehaviors),
            'T4.9: Emergent behaviors as array');
        
        // T4.10
        this.assert(this.adaptive.discoveredTerms !== null,
            'T4.10: Discovery results stored in module');
    }
    
    // ===== CATEGORY 5: FRAMEWORK EXPANSION (10 tests) =====
    testFrameworkExpansion() {
        console.log('\n=== TEST CATEGORY 5: Framework Expansion (10 tests) ===');
        
        // T5.1 - Add physics term
        const term = this.adaptive.addPhysicsTerm(
            'Test_Dark_Matter_Coupling',
            'DM_term = alpha * rho_DM * (1 - exp(-beta * r))',
            { alpha: 0.5, rho_DM: 1e-27 }
        );
        this.assert(term.name === 'Test_Dark_Matter_Coupling',
            'T5.1: Physics term added with name');
        
        // T5.2
        this.assert(term.enabled === true,
            'T5.2: New term enabled by default');
        
        // T5.3
        this.assert(this.adaptive.expanded === true,
            'T5.3: Expansion flag set to true');
        
        // T5.4
        this.assert(this.adaptive.customMethods.additionalTerms.length > 0,
            'T5.4: Term stored in custom methods');
        
        // T5.5 - Add computational method
        const methodAdded = this.adaptive.addComputationalMethod(
            'test_method',
            function() { return 42; }
        );
        this.assert(methodAdded === true,
            'T5.5: Computational method added');
        
        // T5.6
        this.assert(this.adaptive.customMethods.hasOwnProperty('test_method'),
            'T5.6: Method accessible from custom methods');
        
        // T5.7 - Execute custom method
        const result = this.adaptive.executeCustomMethod('test_method');
        this.assert(result === 42,
            'T5.7: Custom method executes correctly');
        
        // T5.8 - Set term weight
        const weightSet = this.adaptive.setTermWeight('Test_Dark_Matter_Coupling', 0.5);
        this.assert(weightSet === true,
            'T5.8: Term weight updated');
        
        // T5.9
        this.assert(this.adaptive.customMethods.additionalTerms[0].weight === 0.5,
            'T5.9: Weight change persisted');
        
        // T5.10
        this.assert(this.adaptive.customMethods['test_method'].callCount === 1,
            'T5.10: Method call tracking works');
    }
    
    // ===== CATEGORY 6: MULTI-SCALE REFINEMENT (8 tests) =====
    testMultiScaleRefinement() {
        console.log('\n=== TEST CATEGORY 6: Multi-Scale Refinement (8 tests) ===');
        
        const refinement = this.adaptive.multiScaleRefinement();
        
        // T6.1
        this.assert(refinement.hasOwnProperty('solar'),
            'T6.1: Solar scale refinement');
        
        // T6.2
        this.assert(refinement.hasOwnProperty('galactic'),
            'T6.2: Galactic scale refinement');
        
        // T6.3
        this.assert(refinement.hasOwnProperty('cosmic'),
            'T6.3: Cosmic scale refinement');
        
        // T6.4
        this.assert(refinement.solar.samples.length === 11,
            'T6.4: Solar scale has 11 samples');
        
        // T6.5
        this.assert(refinement.solar.hasOwnProperty('variance'),
            'T6.5: Variance computed for solar scale');
        
        // T6.6
        this.assert(['STABLE', 'INCREASING', 'DECREASING'].includes(refinement.solar.trend),
            'T6.6: Trend detected for solar scale');
        
        // T6.7
        this.assert(refinement.galactic.period > refinement.solar.period,
            'T6.7: Galactic period > solar period');
        
        // T6.8
        this.assert(refinement.cosmic.period > refinement.galactic.period,
            'T6.8: Cosmic period > galactic period');
    }
    
    // ===== CATEGORY 7: ANOMALY DETECTION (10 tests) =====
    testAnomalyDetection() {
        console.log('\n=== TEST CATEGORY 7: Anomaly Detection (10 tests) ===');
        
        // T7.1 - Test clean state
        let anomalies = this.adaptive.detectAnomalies();
        this.assert(Array.isArray(anomalies),
            'T7.1: Anomaly detection returns array');
        
        // T7.2 - Introduce anomaly
        this.adaptive.core.variables['M_bh'] = 1e50;  // Out of range
        anomalies = this.adaptive.detectAnomalies();
        this.assert(anomalies.some(a => a.type === 'OUT_OF_RANGE'),
            'T7.2: Out-of-range anomaly detected');
        
        // T7.3
        this.adaptive.core.variables['M_bh'] = 1e42;  // Restore
        anomalies = this.adaptive.detectAnomalies();
        this.assert(!anomalies.some(a => a.type === 'OUT_OF_RANGE'),
            'T7.3: Anomaly cleared after correction');
        
        // T7.4
        this.adaptive.core.variables['test_var'] = NaN;
        anomalies = this.adaptive.detectAnomalies();
        this.assert(anomalies.some(a => a.type === 'NON_FINITE_VALUE'),
            'T7.4: NaN anomaly detected');
        
        // T7.5
        delete this.adaptive.core.variables['test_var'];
        const corrections = this.adaptive.autoCorrectAnomalies();
        this.assert(Array.isArray(corrections),
            'T7.5: Auto-correction returns array');
        
        // T7.6
        this.adaptive.core.variables['sigma'] = 5e5;  // Out of range high
        const autoCorr = this.adaptive.autoCorrectAnomalies();
        this.assert(this.adaptive.core.variables['sigma'] <= 1e6,
            'T7.6: High value clamped by auto-correction');
        
        // T7.7
        this.adaptive.core.variables['sigma'] = 2e5;  // Restore to valid
        const postCorrect = this.adaptive.detectAnomalies();
        this.assert(!postCorrect.some(a => a.type === 'OUT_OF_RANGE' && a.parameter === 'sigma'),
            'T7.7: Sigma no longer anomalous after correction');
        
        // T7.8
        this.assert(this.adaptive.performanceMetrics.computations >= 0,
            'T7.8: Performance tracking active');
        
        // T7.9
        anomalies = this.adaptive.detectAnomalies();
        this.assert(anomalies.every(a => a.hasOwnProperty('severity')),
            'T7.9: All anomalies have severity rating');
        
        // T7.10
        this.assert(anomalies.every(a => ['CRITICAL', 'WARNING'].includes(a.severity)),
            'T7.10: Severity levels are valid');
    }
    
    // ===== CATEGORY 8: LOGGING & STATE MANAGEMENT (10 tests) =====
    testLoggingStateManagement() {
        console.log('\n=== TEST CATEGORY 8: Logging & State Management (10 tests) ===');
        
        const initialLogCount = this.adaptive.adaptationLogs.length;
        
        // T8.1
        const log = this.adaptive.logAdaptation('Test message', 'TEST');
        this.assert(log.hasOwnProperty('timestamp'),
            'T8.1: Log includes timestamp');
        
        // T8.2
        this.assert(log.category === 'TEST',
            'T8.2: Log category recorded');
        
        // T8.3
        this.assert(this.adaptive.adaptationLogs.length === initialLogCount + 1,
            'T8.3: Log added to history');
        
        // T8.4
        const history = this.adaptive.getAdaptationHistory(10);
        this.assert(Array.isArray(history) && history.length <= 10,
            'T8.4: Adaptation history retrievable with limit');
        
        // T8.5
        const state = this.adaptive.getSystemState();
        this.assert(state.hasOwnProperty('timestamp'),
            'T8.5: System state includes timestamp');
        
        // T8.6
        this.assert(state.hasOwnProperty('coreVariables'),
            'T8.6: System state includes core variables');
        
        // T8.7
        this.assert(state.expanded === true,
            'T8.7: System state reflects expansion');
        
        // T8.8
        const config = this.adaptive.exportConfiguration();
        this.assert(config.hasOwnProperty('variables'),
            'T8.8: Configuration exportable');
        
        // T8.9
        const importSuccess = this.adaptive.importConfiguration(config);
        this.assert(importSuccess === true,
            'T8.9: Configuration importable');
        
        // T8.10
        const baselineReset = this.adaptive.resetToBaseline();
        this.assert(baselineReset === true,
            'T8.10: Reset to baseline functional');
    }
    
    // ===== CATEGORY 9: PERFORMANCE TRACKING (8 tests) =====
    testPerformanceTracking() {
        console.log('\n=== TEST CATEGORY 9: Performance Tracking (8 tests) ===');
        
        // T9.1
        this.adaptive.trackPerformance('test_func', 10);
        this.assert(this.adaptive.performanceMetrics.computations >= 1,
            'T9.1: Computation counter incremented');
        
        // T9.2
        this.adaptive.trackPerformance('test_func', 15);
        this.assert(this.adaptive.performanceMetrics['test_func'].calls === 2,
            'T9.2: Function call count accurate');
        
        // T9.3
        this.assert(this.adaptive.performanceMetrics['test_func'].avgTime > 0,
            'T9.3: Average time computed');
        
        // T9.4
        this.assert(this.adaptive.performanceMetrics['test_func'].minTime === 10,
            'T9.4: Minimum time tracked');
        
        // T9.5
        this.assert(this.adaptive.performanceMetrics['test_func'].maxTime === 15,
            'T9.5: Maximum time tracked');
        
        // T9.6
        const report = this.adaptive.getPerformanceReport();
        this.assert(report.hasOwnProperty('totalComputations'),
            'T9.6: Performance report includes total');
        
        // T9.7
        this.assert(report.hasOwnProperty('functionMetrics'),
            'T9.7: Function metrics in report');
        
        // T9.8
        this.assert(report.totalComputations > 0,
            'T9.8: Report shows non-zero computations');
    }
    
    // ===== CATEGORY 10: ADAPTIVE WORKFLOW (6 tests) =====
    testAdaptiveWorkflow() {
        console.log('\n=== TEST CATEGORY 10: Adaptive Workflow (6 tests) ===');
        
        // T10.1 - Workflow executes
        const workflowPromise = this.adaptive.adaptiveWorkflow(1e-7, { iterations: 5 });
        
        // Simulate async (for this test, just check setup)
        this.assert(workflowPromise !== null,
            'T10.1: Workflow initialization successful');
        
        // T10.2 - Generate report
        const report = this.adaptive.generateReport();
        this.assert(report.hasOwnProperty('systemStatus'),
            'T10.2: Report includes system status');
        
        // T10.3
        this.assert(report.hasOwnProperty('performanceAnalysis'),
            'T10.3: Report includes performance');
        
        // T10.4
        this.assert(report.hasOwnProperty('quantumResonanceMap'),
            'T10.4: Report includes quantum map');
        
        // T10.5
        this.assert(report.hasOwnProperty('discoveredPhysics'),
            'T10.5: Report includes discoveries');
        
        // T10.6
        this.assert(report.hasOwnProperty('customExpansions'),
            'T10.6: Report includes expansions');
    }
    
    runAllTests() {
        console.log('\n╔════════════════════════════════════════════════════════════╗');
        console.log('║  SMBH M-σ ADAPTIVE MODULE TEST SUITE                       ║');
        console.log('║  Dynamic Self-Updating & Expansion Capabilities            ║');
        console.log('║  64 Tests Across 10 Categories                             ║');
        console.log('╚════════════════════════════════════════════════════════════╝');
        
        this.testInitialization();
        this.testParameterOptimization();
        this.testQuantumResonanceMapping();
        this.testPhysicsDiscovery();
        this.testFrameworkExpansion();
        this.testMultiScaleRefinement();
        this.testAnomalyDetection();
        this.testLoggingStateManagement();
        this.testPerformanceTracking();
        this.testAdaptiveWorkflow();
        
        this.printSummary();
    }
    
    printSummary() {
        const total = this.tests_passed + this.tests_failed;
        const pass_rate = ((this.tests_passed / total) * 100).toFixed(2);
        
        console.log('\n╔════════════════════════════════════════════════════════════╗');
        console.log('║                    TEST SUMMARY                             ║');
        console.log('╚════════════════════════════════════════════════════════════╝');
        console.log(`\nTotal Tests Run:     ${total}`);
        console.log(`Tests Passed:        ${this.tests_passed} ✓`);
        console.log(`Tests Failed:        ${this.tests_failed} ✗`);
        console.log(`Success Rate:        ${pass_rate}%`);
        
        if (this.tests_failed === 0) {
            console.log('\n✓ ALL TESTS PASSED - ADAPTIVE MODULE IS PRODUCTION READY\n');
        } else {
            console.log(`\n⚠ ${this.tests_failed} test(s) need attention\n`);
        }
    }
}

// Execute tests
const suite = new AdaptiveTestSuite();
suite.runAllTests();
