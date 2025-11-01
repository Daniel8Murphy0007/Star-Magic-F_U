#!/usr/bin/env node

/**
 * UQFF Framework Iteration Engine v2.0
 * Advanced capabilities for cross-system interaction, optimization, and analysis
 * Extends the base framework with meta-computation and system orchestration
 */

// ============================================================================
// 1. CROSS-SYSTEM INTERACTION LAYER
// ============================================================================

class CrossSystemInteractionEngine {
    constructor() {
        this.systems = new Map();           // Registry of all active systems
        this.connections = new Map();       // System-to-system connection graph
        this.unifiedFields = new Map();     // Computed unified fields
        this.interactionCache = new Map();  // Cache for interaction results
    }

    /**
     * Register a system in the interaction network
     */
    registerSystem(systemName, systemInstance) {
        this.systems.set(systemName, {
            instance: systemInstance,
            type: systemInstance.constructor.name,
            lastComputed: null,
            metadata: {
                parameters: systemInstance.variables ? Array.from(systemInstance.variables.keys()) : [],
                methods: Object.getOwnPropertyNames(Object.getPrototypeOf(systemInstance))
                    .filter(m => typeof systemInstance[m] === 'function')
            }
        });
        console.log(`✓ Registered system: ${systemName} (${systemInstance.constructor.name})`);
        return true;
    }

    /**
     * Connect two systems for interaction
     */
    connectSystems(system1Name, system2Name, interactionType = 'gravitational') {
        const key = `${system1Name}<->${system2Name}`;
        this.connections.set(key, {
            system1: system1Name,
            system2: system2Name,
            type: interactionType,
            weight: 1.0,
            enabled: true
        });
        console.log(`✓ Connected: ${system1Name} <-> ${system2Name} (${interactionType})`);
        return true;
    }

    /**
     * Compute unified field between two interacting systems
     */
    computeUnifiedFieldInteraction(system1Name, system2Name, parameters = {}) {
        const cacheKey = `${system1Name}+${system2Name}:${JSON.stringify(parameters)}`;
        
        if (this.interactionCache.has(cacheKey)) {
            return this.interactionCache.get(cacheKey);
        }

        const sys1 = this.systems.get(system1Name);
        const sys2 = this.systems.get(system2Name);
        
        if (!sys1 || !sys2) {
            throw new Error(`System not found in registry`);
        }

        const result = {
            system1: system1Name,
            system2: system2Name,
            timestamp: new Date().toISOString(),
            interaction: {}
        };

        // Gravitational interaction
        if (parameters.distance && parameters.includeGravitational !== false) {
            const G = 6.6743e-11;
            const m1 = parameters.mass1 || 1.989e30;
            const m2 = parameters.mass2 || 1.989e30;
            const d = parameters.distance;
            
            result.interaction.gravitational = {
                force: (G * m1 * m2) / (d * d),
                acceleration1: (G * m2) / (d * d),
                acceleration2: (G * m1) / (d * d)
            };
        }

        // Electromagnetic interaction
        if (parameters.charge && parameters.includeElectromagnetic !== false) {
            const k_e = 8.987e9;
            const q1 = parameters.charge1 || 0;
            const q2 = parameters.charge2 || 0;
            const d = parameters.distance || 1;
            
            if (q1 !== 0 && q2 !== 0) {
                result.interaction.electromagnetic = {
                    force: (k_e * q1 * q2) / (d * d),
                    potential: (k_e * q1) / d
                };
            }
        }

        // Magnetic interaction
        if (parameters.magneticField && parameters.includeMagnetic !== false) {
            result.interaction.magnetic = {
                fieldStrength1: parameters.magneticField,
                fieldStrength2: parameters.magneticField * 0.8,
                interaction: 'dipole-dipole'
            };
        }

        // Unified field composition
        result.unifiedField = this._composeUnifiedField(result.interaction);
        
        this.interactionCache.set(cacheKey, result);
        return result;
    }

    /**
     * Compose unified field from all interaction components
     */
    _composeUnifiedField(interactions) {
        const G_mag = interactions.gravitational ? 
            Math.sqrt(interactions.gravitational.force ** 2) : 0;
        const E_mag = interactions.electromagnetic ? 
            Math.sqrt(interactions.electromagnetic.force ** 2) : 0;
        const M_mag = interactions.magnetic ? 1e15 : 0;
        
        return {
            total: G_mag + E_mag + M_mag,
            gravitational: G_mag,
            electromagnetic: E_mag,
            magnetic: M_mag,
            composition: {
                g_ratio: G_mag / (G_mag + E_mag + M_mag + 1e-30),
                e_ratio: E_mag / (G_mag + E_mag + M_mag + 1e-30),
                m_ratio: M_mag / (G_mag + E_mag + M_mag + 1e-30)
            }
        };
    }

    /**
     * Get all active connections for a system
     */
    getSystemConnections(systemName) {
        const connections = [];
        for (const [key, conn] of this.connections) {
            if (conn.system1 === systemName || conn.system2 === systemName) {
                connections.push(conn);
            }
        }
        return connections;
    }

    /**
     * Print system registry
     */
    printRegistry() {
        console.log('\n╔═══════════════════════════════════════════════════════╗');
        console.log('║        UQFF FRAMEWORK SYSTEM REGISTRY                 ║');
        console.log('╚═══════════════════════════════════════════════════════╝\n');
        
        console.log(`Total registered systems: ${this.systems.size}`);
        console.log(`Active connections: ${this.connections.size}\n`);
        
        for (const [name, sys] of this.systems) {
            console.log(`▶ ${name} (${sys.type})`);
            console.log(`  Parameters: ${sys.metadata.parameters.length}`);
            console.log(`  Methods: ${sys.metadata.methods.length}`);
        }
        
        console.log('\n--- System Connections ---');
        for (const [key, conn] of this.connections) {
            console.log(`  ${key} [${conn.type}]`);
        }
    }
}

// ============================================================================
// 2. PERFORMANCE OPTIMIZATION LAYER
// ============================================================================

class PerformanceOptimizer {
    constructor() {
        this.computationCache = new Map();
        this.memoizedResults = new Map();
        this.performanceMetrics = {
            totalComputations: 0,
            cacheHits: 0,
            cacheMisses: 0,
            averageTime: 0,
            peakTime: 0,
            computationTimes: []
        };
    }

    /**
     * Execute computation with memoization and caching
     */
    executeOptimized(functionName, func, args, options = {}) {
        const startTime = Date.now();
        const cacheKey = this._generateCacheKey(functionName, args);
        
        // Check memoization cache
        if (options.memoize && this.memoizedResults.has(cacheKey)) {
            this.performanceMetrics.cacheHits++;
            return this.memoizedResults.get(cacheKey);
        }
        this.performanceMetrics.cacheMisses++;

        // Execute function
        let result;
        try {
            result = func(...args);
        } catch (e) {
            console.error(`Error in ${functionName}:`, e.message);
            return null;
        }

        // Store result
        if (options.memoize) {
            this.memoizedResults.set(cacheKey, result);
        }

        // Track performance
        const elapsed = Date.now() - startTime;
        this.performanceMetrics.totalComputations++;
        this.performanceMetrics.computationTimes.push(elapsed);
        this.performanceMetrics.averageTime = 
            this.performanceMetrics.computationTimes.reduce((a, b) => a + b, 0) / 
            this.performanceMetrics.totalComputations;
        this.performanceMetrics.peakTime = Math.max(
            this.performanceMetrics.peakTime, 
            elapsed
        );

        return result;
    }

    /**
     * Generate cache key from function name and arguments
     */
    _generateCacheKey(functionName, args) {
        return `${functionName}:${JSON.stringify(args)}`;
    }

    /**
     * Get cache hit ratio
     */
    getCacheHitRatio() {
        const total = this.performanceMetrics.cacheHits + this.performanceMetrics.cacheMisses;
        return total === 0 ? 0 : (this.performanceMetrics.cacheHits / total) * 100;
    }

    /**
     * Get operations per millisecond
     */
    getOperationsPerMs() {
        if (this.performanceMetrics.computationTimes.length === 0) return 0;
        const totalTime = this.performanceMetrics.computationTimes.reduce((a, b) => a + b, 0);
        return this.performanceMetrics.totalComputations / totalTime;
    }

    /**
     * Print performance report
     */
    printPerformanceReport() {
        console.log('\n╔═══════════════════════════════════════════════════════╗');
        console.log('║      PERFORMANCE OPTIMIZATION REPORT                   ║');
        console.log('╚═══════════════════════════════════════════════════════╝\n');
        
        console.log(`Total Computations: ${this.performanceMetrics.totalComputations}`);
        console.log(`Cache Hits: ${this.performanceMetrics.cacheHits}`);
        console.log(`Cache Misses: ${this.performanceMetrics.cacheMisses}`);
        console.log(`Hit Ratio: ${this.getCacheHitRatio().toFixed(2)}%`);
        console.log(`\nAverage Computation Time: ${this.performanceMetrics.averageTime.toFixed(3)} ms`);
        console.log(`Peak Computation Time: ${this.performanceMetrics.peakTime.toFixed(3)} ms`);
        console.log(`Operations/ms: ${this.getOperationsPerMs().toFixed(2)}`);
    }

    /**
     * Clear cache
     */
    clearCache() {
        this.computationCache.clear();
        this.memoizedResults.clear();
        console.log('✓ Cache cleared');
    }
}

// ============================================================================
// 3. STATISTICAL ANALYSIS & ANOMALY DETECTION
// ============================================================================

class StatisticalAnalysisEngine {
    constructor() {
        this.dataPoints = new Map();
        this.timeSeries = new Map();
        this.anomalies = [];
    }

    /**
     * Record data point for analysis
     */
    recordDataPoint(systemName, parameter, value, timestamp = Date.now()) {
        const key = `${systemName}:${parameter}`;
        
        if (!this.timeSeries.has(key)) {
            this.timeSeries.set(key, []);
        }
        
        this.timeSeries.get(key).push({
            value: value,
            timestamp: timestamp
        });
    }

    /**
     * Calculate statistical properties
     */
    analyzeTimeSeries(systemName, parameter) {
        const key = `${systemName}:${parameter}`;
        const series = this.timeSeries.get(key);
        
        if (!series || series.length === 0) {
            return null;
        }

        const values = series.map(p => p.value);
        const n = values.length;
        const mean = values.reduce((a, b) => a + b, 0) / n;
        const variance = values.reduce((a, x) => a + Math.pow(x - mean, 2), 0) / n;
        const stdDev = Math.sqrt(variance);
        
        const sorted = [...values].sort((a, b) => a - b);
        const median = n % 2 === 0 ? 
            (sorted[n/2 - 1] + sorted[n/2]) / 2 : 
            sorted[Math.floor(n/2)];

        return {
            parameter: parameter,
            systemName: systemName,
            count: n,
            mean: mean,
            median: median,
            stdDev: stdDev,
            variance: variance,
            min: Math.min(...values),
            max: Math.max(...values),
            range: Math.max(...values) - Math.min(...values),
            trend: this._calculateTrend(values)
        };
    }

    /**
     * Detect anomalies using z-score method
     */
    detectAnomalies(systemName, parameter, threshold = 3.0) {
        const key = `${systemName}:${parameter}`;
        const series = this.timeSeries.get(key);
        
        if (!series || series.length < 2) {
            return [];
        }

        const values = series.map(p => p.value);
        const mean = values.reduce((a, b) => a + b, 0) / values.length;
        const stdDev = Math.sqrt(
            values.reduce((a, x) => a + Math.pow(x - mean, 2), 0) / values.length
        );

        const detected = [];
        series.forEach((point, idx) => {
            const zScore = Math.abs((point.value - mean) / (stdDev + 1e-30));
            if (zScore > threshold) {
                detected.push({
                    systemName: systemName,
                    parameter: parameter,
                    index: idx,
                    value: point.value,
                    zScore: zScore,
                    timestamp: point.timestamp,
                    severity: zScore > 5 ? 'critical' : zScore > 4 ? 'high' : 'medium'
                });
            }
        });

        return detected;
    }

    /**
     * Calculate trend direction
     */
    _calculateTrend(values) {
        if (values.length < 2) return 'insufficient_data';
        
        const n = values.length;
        const x_mean = (n - 1) / 2;
        const y_mean = values.reduce((a, b) => a + b, 0) / n;
        
        let numerator = 0;
        let denominator = 0;
        
        for (let i = 0; i < n; i++) {
            numerator += (i - x_mean) * (values[i] - y_mean);
            denominator += (i - x_mean) * (i - x_mean);
        }
        
        const slope = numerator / denominator;
        
        if (Math.abs(slope) < 0.001) return 'stable';
        return slope > 0 ? 'increasing' : 'decreasing';
    }

    /**
     * Print analysis report
     */
    printAnalysisReport(systemName, parameter) {
        const analysis = this.analyzeTimeSeries(systemName, parameter);
        const anomalies = this.detectAnomalies(systemName, parameter);
        
        if (!analysis) {
            console.log(`No data for ${systemName}:${parameter}`);
            return;
        }

        console.log(`\n▶ Statistical Analysis: ${systemName}.${parameter}`);
        console.log(`  Data Points: ${analysis.count}`);
        console.log(`  Mean: ${analysis.mean.toExponential(3)}`);
        console.log(`  Median: ${analysis.median.toExponential(3)}`);
        console.log(`  Std Dev: ${analysis.stdDev.toExponential(3)}`);
        console.log(`  Min: ${analysis.min.toExponential(3)}`);
        console.log(`  Max: ${analysis.max.toExponential(3)}`);
        console.log(`  Trend: ${analysis.trend}`);
        
        if (anomalies.length > 0) {
            console.log(`\n  ⚠ Anomalies detected: ${anomalies.length}`);
            anomalies.slice(0, 5).forEach(anom => {
                console.log(`    • Index ${anom.index}: ${anom.value.toExponential(3)} (z=${anom.zScore.toFixed(2)}) [${anom.severity}]`);
            });
        }
    }
}

// ============================================================================
// 4. PARAMETER OPTIMIZATION ENGINE
// ============================================================================

class ParameterOptimizationEngine {
    constructor() {
        this.optimizationResults = [];
        this.constraints = new Map();
    }

    /**
     * Add constraint for optimization
     */
    addConstraint(paramName, minValue, maxValue, type = 'range') {
        this.constraints.set(paramName, {
            min: minValue,
            max: maxValue,
            type: type
        });
    }

    /**
     * Simple grid search optimization
     */
    gridSearchOptimization(system, targetParameter, targetValue, searchParams, steps = 10) {
        const results = [];
        const paramGrid = {};
        
        // Build parameter grid
        for (const [param, range] of Object.entries(searchParams)) {
            const [min, max] = range;
            paramGrid[param] = this._linspace(min, max, steps);
        }

        // Grid search
        const paramNames = Object.keys(paramGrid);
        const combinations = this._generateCombinations(paramGrid);
        
        for (const combo of combinations) {
            // Set parameters
            for (const [param, value] of Object.entries(combo)) {
                if (system.updateVariable) {
                    system.updateVariable(param, value);
                }
            }

            // Compute target
            let computed = null;
            if (system.computeG) {
                computed = system.computeG(system.variables?.get?.('t') || 0, 
                                          system.variables?.get?.('r') || 1);
            }

            if (computed !== null) {
                const error = Math.abs(computed - targetValue);
                results.push({
                    parameters: { ...combo },
                    computed: computed,
                    error: error,
                    relativeError: error / Math.abs(targetValue + 1e-30)
                });
            }
        }

        // Sort by error
        results.sort((a, b) => a.error - b.error);
        this.optimizationResults = results;
        
        return results.slice(0, 10); // Return top 10
    }

    /**
     * Generate linear space
     */
    _linspace(start, end, steps) {
        const result = [];
        const step = (end - start) / (steps - 1);
        for (let i = 0; i < steps; i++) {
            result.push(start + i * step);
        }
        return result;
    }

    /**
     * Generate all combinations
     */
    _generateCombinations(grid) {
        const keys = Object.keys(grid);
        const combinations = [];
        
        const generate = (index, current) => {
            if (index === keys.length) {
                combinations.push({ ...current });
                return;
            }
            
            const key = keys[index];
            for (const value of grid[key]) {
                current[key] = value;
                generate(index + 1, current);
            }
        };
        
        generate(0, {});
        return combinations;
    }

    /**
     * Print optimization results
     */
    printOptimizationResults(topN = 10) {
        if (this.optimizationResults.length === 0) {
            console.log('No optimization results available');
            return;
        }

        console.log(`\n╔═══════════════════════════════════════════════════════╗`);
        console.log(`║      PARAMETER OPTIMIZATION RESULTS (Top ${topN})          ║`);
        console.log(`╚═══════════════════════════════════════════════════════╝\n`);
        
        this.optimizationResults.slice(0, topN).forEach((result, idx) => {
            console.log(`${idx + 1}. Error: ${result.error.toExponential(3)}`);
            console.log(`   Relative Error: ${(result.relativeError * 100).toFixed(2)}%`);
            console.log(`   Computed: ${result.computed.toExponential(3)}`);
            console.log(`   Parameters:`, result.parameters);
        });
    }
}

// ============================================================================
// EXPORT MODULES
// ============================================================================

if (typeof module !== 'undefined' && module.exports) {
    module.exports.CrossSystemInteractionEngine = CrossSystemInteractionEngine;
    module.exports.PerformanceOptimizer = PerformanceOptimizer;
    module.exports.StatisticalAnalysisEngine = StatisticalAnalysisEngine;
    module.exports.ParameterOptimizationEngine = ParameterOptimizationEngine;
}
