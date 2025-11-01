# Source81 (NGC 346) Dynamic Capability Enhancement - Detailed Patch Plan

**Status**: Option B - Core Module Modification (Non-Invasive Enhancement)  
**Target File**: `ngc346_uqff.js` (currently 448 lines)  
**Objective**: Add runtime variable maps, physics-term registration, state serialization, and adaptive APIs to achieve feature parity with Source82 adaptive capabilities.

---

## Executive Summary

Source81 currently has **parameter-level runtime updates** (e.g., `updateVariable(name, value)`), but lacks:
1. **Generic variable map** for arbitrary named variables
2. **Physics term registration system** (runtime add/remove of computation terms)
3. **State serialization/restore APIs** (snapshots and import/export)
4. **Adaptive workflow orchestration** (auto-calibration, discovery, anomaly detection)
5. **Comprehensive testing** for runtime extensibility

This patch plan introduces non-invasive enhancements to `ngc346_uqff.js` to enable:
- ✅ Dynamic variable storage and retrieval
- ✅ Runtime addition of custom physics terms
- ✅ State snapshots and restoration
- ✅ Adaptive optimization and calibration loops
- ✅ Anomaly detection and auto-correction
- ✅ Full backward compatibility with existing tests

**Estimated Code Changes**: ~500–700 lines added to `ngc346_uqff.js`  
**Estimated Test Coverage**: ~80 new adaptive tests (similar pattern to S82 adaptive tests)  
**Risk Level**: **Low** (new methods, no existing method modification; existing tests unaffected)

---

## Architecture Overview

### Current State (NGC346UQFFModule - Basic)

```
Constructor
├── Initialize 57 static variables
├── Set NGC 346 parameters
└── Ready for `computeG(t, r)` calls

Public Methods (Basic):
├── computeG(t, r) → gravity acceleration
├── getEquationText() → documentation
├── printVariables() → debug output
└── updateVariable(name, value) → parameter update

Private Computation Methods (14):
├── computeHtz(z), computeFenv(t), computeMsfFactor(t)
├── computeUg1/2/3/4, computeUi, computeUm
├── computeEcore, computeTempCore
└── ... 6 other specialized methods

Limitations:
❌ No generic variable registration beyond hardcoded set
❌ No runtime term addition (physics extensibility)
❌ No state export/import (snapshots)
❌ No adaptive loops (optimization, discovery)
❌ No anomaly detection framework
```

### Target State (NGC346UQFFModule - Enhanced)

```
Constructor
├── Initialize 57 static variables (unchanged)
├── Initialize adaptive system (NEW)
│   ├── physicsTerms = {} (registry of custom terms)
│   ├── calibrationData = {} (snapshots & metadata)
│   ├── discoveredTerms = {} (discovered relationships)
│   ├── adaptationLogs = [] (audit trail)
│   └── customMethods = {} (runtime added functions)
└── Ready for both basic and adaptive operations

Public Methods (Enhanced):
├── === BASIC (unchanged) ===
│   ├── computeG(t, r) → gravity acceleration
│   ├── getEquationText() → documentation
│   ├── printVariables() → debug output
│   └── updateVariable(name, value) → parameter update
│
├── === VARIABLE MANAGEMENT (NEW) ===
│   ├── addVariable(name, value) → register new variable
│   ├── getVariable(name) → retrieve with type checking
│   ├── removeVariable(name) → unregister variable
│   ├── getVariables() → return all variables snapshot
│   └── listVariables() → list all registered names
│
├── === PHYSICS TERM REGISTRATION (NEW) ===
│   ├── registerPhysicsTerm(name, computeFn, enabled=true)
│   ├── unregisterPhysicsTerm(name) → remove custom term
│   ├── enablePhysicsTerm(name) → toggle on
│   ├── disablePhysicsTerm(name) → toggle off
│   ├── getPhysicsTerms() → list registered terms
│   └── computeCustomTerms(t, r) → sum all custom terms
│
├── === STATE MANAGEMENT (NEW) ===
│   ├── getState() → export { variables, terms, calibration }
│   ├── setState(state) → import snapshot
│   ├── saveCheckpoint(label) → named snapshot
│   ├── loadCheckpoint(label) → restore named snapshot
│   └── listCheckpoints() → enumerate saved states
│
├── === ADAPTIVE OPERATIONS (NEW) ===
│   ├── optimizeParameters(targetData, iterations, learningRate)
│   ├── mapQuantumResonance(numStates) → resonance landscape
│   ├── discoverPhysics(scanRange) → detect new relationships
│   ├── detectAnomalies(threshold) → identify outliers
│   ├── autoCorrectAnomalies() → apply corrections
│   ├── addCustomMethod(name, fn) → register arbitrary computation
│   ├── executeCustomMethod(name, ...args) → invoke stored method
│   └── adaptiveWorkflow(targetData, config) → full orchestration
│
├── === REPORTING (NEW) ===
│   ├── generateReport() → comprehensive state report
│   ├── exportConfiguration() → JSON-serializable config
│   └── getAdaptationLog() → audit trail
│
└── === UTILITIES (NEW) ===
    ├── clampToPhysicalRange(param, value)
    ├── validateVariable(name, value)
    └── logAdaptation(msg, category)

Private Computation Methods (unchanged):
├── computeHtz(z), computeFenv(t), computeMsfFactor(t)
├── computeUg1/2/3/4, computeUi, computeUm
├── computeEcore, computeTempCore
└── ... (14 methods, no modification)

Key Additions:
✅ Adaptive system initialization
✅ Physics term registry with enable/disable
✅ Calibration data storage with timestamps
✅ Discovered terms and relationships
✅ Quantum resonance map
✅ Anomaly detection metrics
✅ Custom method storage
✅ Comprehensive logging
```

---

## Detailed Patch Sections

### 1. Constructor Initialization (NEW Sections)

**Location**: End of `constructor()` (after all variable initializations, before closing brace)

**Change Type**: Addition (non-breaking)

#### Before:
```javascript
// constructor() last lines (existing)
this.variables['delta_rho_over_rho'] = 1e-5;         // Density perturbation ratio
}
```

#### After:
```javascript
// constructor() last lines (existing)
this.variables['delta_rho_over_rho'] = 1e-5;         // Density perturbation ratio

// ===== ADAPTIVE SYSTEM INITIALIZATION (NEW) =====
this.physicsTerms = {};           // Registry: { termName: { fn, enabled, description } }
this.calibrationData = {};        // Snapshots: { label: { timestamp, variables, metadata } }
this.discoveredTerms = {};        // New physics: { termName: { relationship, sensitivity } }
this.adaptationLogs = [];         // Audit trail: [{ timestamp, category, message }]
this.customMethods = {};          // User functions: { methodName: fn }
this.quantumMap = {};             // Resonance: { stateIndex: { energy, amplitude, ...} }
this.performanceMetrics = {       // Stats
    computations: 0,
    totalTime: 0,
    averageTime: 0
};

this.initializeAdaptiveSystem();  // Initialize calibration baseline
}

initializeAdaptiveSystem() {
    // Store baseline calibration at construction time
    this.calibrationData['_baseline'] = {
        timestamp: Date.now(),
        variables: JSON.parse(JSON.stringify(this.variables)),
        config: 'initialization'
    };
    
    this.logAdaptation('Adaptive system initialized (57 variables, 14 compute methods)', 'INIT');
}
```

**Impact**: ✅ Fully backward-compatible. Existing tests unaffected.

---

### 2. Variable Management APIs (NEW Methods)

**Location**: After `printVariables()` (line ~350)

**Change Type**: Addition (3 new methods)

#### PATCH CONTENT:

```javascript
    // ===== ADAPTIVE: VARIABLE MANAGEMENT (NEW) =====
    
    /**
     * Register a new variable in the system
     * @param {string} name - Variable identifier
     * @param {number} value - Initial value
     * @returns {boolean} - Success flag
     */
    addVariable(name, value) {
        if (this.variables.hasOwnProperty(name)) {
            console.warn(`Variable '${name}' already exists. Use updateVariable() to modify.`);
            return false;
        }
        this.variables[name] = value;
        this.logAdaptation(`Variable added: ${name} = ${value}`, 'VAR_ADD');
        return true;
    }
    
    /**
     * Retrieve a variable with existence check
     * @param {string} name - Variable identifier
     * @returns {number|null} - Value or null if not found
     */
    getVariable(name) {
        return this.variables.hasOwnProperty(name) ? this.variables[name] : null;
    }
    
    /**
     * Remove a variable from the system
     * @param {string} name - Variable identifier
     * @returns {boolean} - Success flag
     */
    removeVariable(name) {
        if (!this.variables.hasOwnProperty(name)) {
            console.warn(`Variable '${name}' not found.`);
            return false;
        }
        delete this.variables[name];
        this.logAdaptation(`Variable removed: ${name}`, 'VAR_REMOVE');
        return true;
    }
    
    /**
     * Get snapshot of all variables (for comparison/export)
     * @returns {Object} - Copy of variables map
     */
    getVariables() {
        return JSON.parse(JSON.stringify(this.variables));
    }
    
    /**
     * List all registered variable names
     * @returns {Array<string>} - Sorted list of variable names
     */
    listVariables() {
        return Object.keys(this.variables).sort();
    }
```

**Impact**: ✅ New public APIs. No breaking changes. Enables dynamic variable access.

---

### 3. Physics Term Registration (NEW Methods)

**Location**: After variable management methods

**Change Type**: Addition (6 new methods)

#### PATCH CONTENT:

```javascript
    // ===== ADAPTIVE: PHYSICS TERM REGISTRATION (NEW) =====
    
    /**
     * Register a custom physics term that is automatically included in computations
     * @param {string} name - Term identifier (e.g., 'custom_dark_energy_correction')
     * @param {Function} computeFn - Async computation function(t, r, variables) => number
     * @param {boolean} enabled - Initially enabled? (default: true)
     * @param {string} description - Human-readable description
     * @returns {boolean} - Success flag
     */
    registerPhysicsTerm(name, computeFn, enabled = true, description = '') {
        if (this.physicsTerms.hasOwnProperty(name)) {
            console.warn(`Physics term '${name}' already registered.`);
            return false;
        }
        
        if (typeof computeFn !== 'function') {
            console.error(`Physics term '${name}': computeFn must be a function.`);
            return false;
        }
        
        this.physicsTerms[name] = {
            fn: computeFn,
            enabled: enabled,
            description: description,
            registeredAt: Date.now()
        };
        
        this.logAdaptation(
            `Physics term registered: ${name} (${enabled ? 'enabled' : 'disabled'}) - ${description}`,
            'TERM_REG'
        );
        return true;
    }
    
    /**
     * Unregister a physics term (remove from registry)
     * @param {string} name - Term identifier
     * @returns {boolean} - Success flag
     */
    unregisterPhysicsTerm(name) {
        if (!this.physicsTerms.hasOwnProperty(name)) {
            console.warn(`Physics term '${name}' not found.`);
            return false;
        }
        delete this.physicsTerms[name];
        this.logAdaptation(`Physics term unregistered: ${name}`, 'TERM_UNREG');
        return true;
    }
    
    /**
     * Enable a registered physics term
     * @param {string} name - Term identifier
     * @returns {boolean} - Success flag
     */
    enablePhysicsTerm(name) {
        if (!this.physicsTerms.hasOwnProperty(name)) {
            console.warn(`Physics term '${name}' not found.`);
            return false;
        }
        this.physicsTerms[name].enabled = true;
        this.logAdaptation(`Physics term enabled: ${name}`, 'TERM_ENABLE');
        return true;
    }
    
    /**
     * Disable a registered physics term (without removing from registry)
     * @param {string} name - Term identifier
     * @returns {boolean} - Success flag
     */
    disablePhysicsTerm(name) {
        if (!this.physicsTerms.hasOwnProperty(name)) {
            console.warn(`Physics term '${name}' not found.`);
            return false;
        }
        this.physicsTerms[name].enabled = false;
        this.logAdaptation(`Physics term disabled: ${name}`, 'TERM_DISABLE');
        return true;
    }
    
    /**
     * Get list of all registered physics terms
     * @returns {Array<Object>} - [ { name, enabled, description } ]
     */
    getPhysicsTerms() {
        return Object.entries(this.physicsTerms).map(([name, term]) => ({
            name: name,
            enabled: term.enabled,
            description: term.description,
            registeredAt: term.registeredAt
        }));
    }
    
    /**
     * Compute sum of all enabled custom physics terms
     * Called from within computeG() to augment result
     * @param {number} t - Time parameter
     * @param {number} r - Radius parameter
     * @returns {number} - Sum of all enabled term contributions
     */
    computeCustomTerms(t, r) {
        let sum = 0.0;
        for (const [name, term] of Object.entries(this.physicsTerms)) {
            if (term.enabled) {
                try {
                    const value = term.fn(t, r, this.variables);
                    sum += (typeof value === 'number' ? value : 0.0);
                } catch (err) {
                    console.error(`Error computing physics term '${name}':`, err);
                    this.logAdaptation(`Physics term error: ${name}: ${err.message}`, 'TERM_ERROR');
                }
            }
        }
        return sum;
    }
```

**Impact**: ✅ New public APIs for runtime physics extensibility. Fully backward-compatible.

---

### 4. State Management APIs (NEW Methods)

**Location**: After physics term registration methods

**Change Type**: Addition (5 new methods)

#### PATCH CONTENT:

```javascript
    // ===== ADAPTIVE: STATE MANAGEMENT (NEW) =====
    
    /**
     * Export complete state (variables, physics terms, calibration metadata)
     * @returns {Object} - Serializable state snapshot
     */
    getState() {
        return {
            variables: JSON.parse(JSON.stringify(this.variables)),
            physicsTerms: Object.entries(this.physicsTerms).map(([name, term]) => ({
                name: name,
                enabled: term.enabled,
                description: term.description,
                registeredAt: term.registeredAt
                // fn is not serialized (functions don't JSON.stringify)
            })),
            calibrationKeys: Object.keys(this.calibrationData),
            discoveredTerms: this.discoveredTerms,
            timestamp: Date.now()
        };
    }
    
    /**
     * Import/restore state from a snapshot
     * Note: Physics term functions must be re-registered separately
     * @param {Object} state - State object from getState()
     * @returns {boolean} - Success flag
     */
    setState(state) {
        try {
            if (!state.variables || typeof state.variables !== 'object') {
                console.error('Invalid state object: missing or invalid variables.');
                return false;
            }
            
            // Restore variables
            this.variables = JSON.parse(JSON.stringify(state.variables));
            
            // Restore discovered terms metadata
            if (state.discoveredTerms && typeof state.discoveredTerms === 'object') {
                this.discoveredTerms = state.discoveredTerms;
            }
            
            // Note: physicsTerms functions cannot be restored from JSON; user must re-register
            if (state.physicsTerms && Array.isArray(state.physicsTerms)) {
                console.info(`State restored. Note: ${state.physicsTerms.length} physics term(s) need re-registration.`);
            }
            
            this.logAdaptation('State restored from snapshot', 'STATE_RESTORE');
            return true;
        } catch (err) {
            console.error('Error restoring state:', err);
            this.logAdaptation(`State restore error: ${err.message}`, 'STATE_ERROR');
            return false;
        }
    }
    
    /**
     * Save named checkpoint (snapshot with label)
     * @param {string} label - Checkpoint identifier
     * @returns {boolean} - Success flag
     */
    saveCheckpoint(label) {
        if (!label || typeof label !== 'string') {
            console.error('Checkpoint label must be a non-empty string.');
            return false;
        }
        
        this.calibrationData[`checkpoint_${label}_${Date.now()}`] = {
            label: label,
            timestamp: Date.now(),
            variables: JSON.parse(JSON.stringify(this.variables)),
            physicsTermsCount: Object.keys(this.physicsTerms).length,
            metadata: 'user checkpoint'
        };
        
        this.logAdaptation(`Checkpoint saved: ${label}`, 'CHECKPOINT_SAVE');
        return true;
    }
    
    /**
     * Load a named checkpoint (by pattern matching)
     * @param {string} label - Checkpoint label (pattern match on label field)
     * @returns {boolean} - Success flag
     */
    loadCheckpoint(label) {
        // Find checkpoint by label
        let foundCheckpoint = null;
        for (const [key, checkpoint] of Object.entries(this.calibrationData)) {
            if (checkpoint.label === label) {
                foundCheckpoint = checkpoint;
                break;
            }
        }
        
        if (!foundCheckpoint) {
            console.warn(`Checkpoint '${label}' not found.`);
            return false;
        }
        
        // Restore variables from checkpoint
        this.variables = JSON.parse(JSON.stringify(foundCheckpoint.variables));
        this.logAdaptation(`Checkpoint loaded: ${label}`, 'CHECKPOINT_LOAD');
        return true;
    }
    
    /**
     * List all saved checkpoints
     * @returns {Array<Object>} - [ { label, timestamp, physicsTermsCount } ]
     */
    listCheckpoints() {
        return Object.values(this.calibrationData)
            .filter(cp => cp.label !== undefined)
            .map(cp => ({
                label: cp.label,
                timestamp: cp.timestamp,
                physicsTermsCount: cp.physicsTermsCount || 0
            }))
            .sort((a, b) => b.timestamp - a.timestamp);
    }
```

**Impact**: ✅ New state serialization APIs. Enables checkpointing and state comparison. No breaking changes.

---

### 5. Adaptive Operations (NEW Methods)

**Location**: After state management methods

**Change Type**: Addition (8 new methods)

#### PATCH CONTENT:

```javascript
    // ===== ADAPTIVE: OPTIMIZATION & DISCOVERY (NEW) =====
    
    /**
     * Auto-optimize tunable parameters toward target observational data
     * Uses gradient descent with adaptive learning rate
     * @param {number|Object} targetData - Target value or { value, weight }
     * @param {number} iterations - Optimization iterations (default: 100)
     * @param {number} learningRate - Initial learning rate (default: 0.01)
     * @returns {Object} - { converged, finalError, iterations, optimizedParams }
     */
    optimizeParameters(targetData, iterations = 100, learningRate = 0.01) {
        this.logAdaptation(`Starting parameter optimization (${iterations} iterations)`, 'OPTIMIZE');
        
        const optimizationHistory = [];
        let bestError = Infinity;
        let bestState = {};
        
        const tunableParams = ['M_visible', 'SFR', 'rho_gas', 'B', 'v_rad'];
        
        for (let iter = 0; iter < iterations; iter++) {
            // Compute current output
            const t = this.variables['t'];
            const r = this.variables['r'];
            const currentOutput = this.computeG(t, r) + this.computeCustomTerms(t, r);
            
            // Compute error vs. target
            const targetVal = typeof targetData === 'number' ? targetData : targetData.value;
            const error = Math.abs((currentOutput - targetVal) / Math.max(Math.abs(targetVal), 1e-20));
            optimizationHistory.push({ iteration: iter, error: error });
            
            if (error < bestError) {
                bestError = error;
                bestState = JSON.parse(JSON.stringify(this.variables));
            }
            
            // Compute simple finite-difference gradients for each tunable param
            const epsilon = 1e-10;
            const gradients = {};
            
            for (const param of tunableParams) {
                if (!this.variables.hasOwnProperty(param)) continue;
                
                const original = this.variables[param];
                
                // Forward
                this.variables[param] = original + epsilon;
                const fPlus = this.computeG(t, r) + this.computeCustomTerms(t, r);
                
                // Backward
                this.variables[param] = original - epsilon;
                const fMinus = this.computeG(t, r) + this.computeCustomTerms(t, r);
                
                this.variables[param] = original;
                
                gradients[param] = (fPlus - fMinus) / (2 * epsilon);
            }
            
            // Gradient descent step
            for (const param of tunableParams) {
                if (gradients.hasOwnProperty(param)) {
                    this.variables[param] -= learningRate * gradients[param];
                    this.variables[param] = this.clampToPhysicalRange(param, this.variables[param]);
                }
            }
            
            // Decay learning rate
            if (iter % 10 === 0 && iter > 0) {
                learningRate *= 0.95;
            }
        }
        
        // Restore best state
        this.variables = bestState;
        
        const result = {
            converged: bestError < 0.1,
            finalError: bestError,
            iterations: iterations,
            history: optimizationHistory,
            optimizedParams: bestState
        };
        
        this.calibrationData[`optimized_${Date.now()}`] = result;
        this.logAdaptation(
            `Optimization complete: error=${bestError.toExponential(2)}, converged=${result.converged}`,
            'OPTIMIZE'
        );
        
        return result;
    }
    
    /**
     * Map quantum resonance across energy levels (26 levels per UQFF)
     * Returns resonance landscape for discovering quantum transitions
     * @param {number} numStates - Number of energy levels to scan (default: 26)
     * @returns {Object} - { states: [...], resonancePeaks: [...] }
     */
    mapQuantumResonance(numStates = 26) {
        this.logAdaptation(`Mapping quantum resonance (${numStates} states)`, 'QUANTUM_MAP');
        
        const states = [];
        const resonancePeaks = [];
        
        const E0 = 1e-40;  // Base energy from UQFF framework
        
        for (let n = 1; n <= numStates; n++) {
            const En = E0 * Math.pow(10, n);  // E_n = E_0 * 10^n
            const amplitude = Math.exp(-n / 8);  // Gaussian decay
            const frequency = n * 1e8;  // Hz
            
            const state = {
                n: n,
                energy: En,
                amplitude: amplitude,
                frequency: frequency,
                resonance: amplitude * Math.cos(frequency * this.variables['t'])
            };
            
            states.push(state);
            
            if (state.resonance > 0.5) {
                resonancePeaks.push({ n: n, resonance: state.resonance });
            }
        }
        
        this.quantumMap = states;
        
        const result = {
            statesCount: numStates,
            states: states,
            resonancePeaks: resonancePeaks,
            timestamp: Date.now()
        };
        
        this.logAdaptation(
            `Quantum map: ${numStates} states mapped, ${resonancePeaks.length} peaks detected`,
            'QUANTUM_MAP'
        );
        
        return result;
    }
    
    /**
     * Discover new physics relationships by scanning parameter space
     * Identifies trends and proposes new physics terms
     * @param {Object} scanRange - { paramName: { min, max, steps } }
     * @returns {Object} - { discoveries: [...], sensitivities: {...} }
     */
    discoverPhysics(scanRange = {}) {
        this.logAdaptation('Starting physics discovery scan', 'DISCOVER');
        
        // Default scan if not provided
        if (Object.keys(scanRange).length === 0) {
            scanRange = {
                'rho_gas': { min: 1e-21, max: 1e-19, steps: 5 },
                'B': { min: 1e-6, max: 1e-4, steps: 5 },
                'v_rad': { min: -20e3, max: -1e3, steps: 5 }
            };
        }
        
        const discoveries = [];
        const sensitivities = {};
        
        // Grid scan
        const paramNames = Object.keys(scanRange);
        for (const paramName of paramNames) {
            const range = scanRange[paramName];
            sensitivities[paramName] = [];
            
            const originalValue = this.variables[paramName];
            let minVal = Infinity, maxVal = -Infinity;
            
            for (let i = 0; i < range.steps; i++) {
                const val = range.min + (range.max - range.min) * (i / range.steps);
                this.variables[paramName] = val;
                
                const result = this.computeG(this.variables['t'], this.variables['r']);
                sensitivities[paramName].push({ paramValue: val, output: result });
                
                minVal = Math.min(minVal, result);
                maxVal = Math.max(maxVal, result);
            }
            
            this.variables[paramName] = originalValue;  // Restore
            
            // Report sensitivity
            const sensitivity = (maxVal - minVal) / Math.max(Math.abs(minVal), 1e-20);
            discoveries.push({
                parameter: paramName,
                sensitivity: sensitivity,
                rangeMin: range.min,
                rangeMax: range.max,
                outputMin: minVal,
                outputMax: maxVal
            });
        }
        
        // Store discovered relationships
        for (const disc of discoveries) {
            this.discoveredTerms[`scaling_${disc.parameter}`] = {
                relationship: `g ∝ ${disc.parameter}^(sensitivity=${disc.sensitivity.toExponential(2)})`,
                sensitivity: disc.sensitivity,
                timestamp: Date.now()
            };
        }
        
        const result = {
            discoveries: discoveries,
            sensitivities: sensitivities,
            timestamp: Date.now()
        };
        
        this.logAdaptation(
            `Physics discovery: ${discoveries.length} parameter relationships mapped`,
            'DISCOVER'
        );
        
        return result;
    }
    
    /**
     * Detect anomalies in system state vs. expected range
     * @param {number} threshold - Anomaly threshold (default: 2.0 sigma)
     * @returns {Object} - { anomalies: [...], count: number }
     */
    detectAnomalies(threshold = 2.0) {
        this.logAdaptation(`Anomaly detection scan (threshold: ${threshold} sigma)`, 'ANOMALY_DETECT');
        
        const anomalies = [];
        const baselineVars = this.calibrationData['_baseline']?.variables || this.variables;
        
        for (const [key, baselineVal] of Object.entries(baselineVars)) {
            const currentVal = this.variables[key];
            
            if (typeof currentVal === 'number' && typeof baselineVal === 'number' && baselineVal !== 0) {
                const relDeviation = Math.abs((currentVal - baselineVal) / Math.max(Math.abs(baselineVal), 1e-20));
                
                if (relDeviation > threshold) {
                    anomalies.push({
                        variable: key,
                        baselineValue: baselineVal,
                        currentValue: currentVal,
                        relativeDeviation: relDeviation,
                        severity: relDeviation > threshold * 2 ? 'HIGH' : 'MEDIUM'
                    });
                }
            }
        }
        
        const result = {
            anomalies: anomalies,
            count: anomalies.length,
            threshold: threshold,
            timestamp: Date.now()
        };
        
        this.logAdaptation(`Anomaly detection: ${anomalies.length} anomalies found`, 'ANOMALY_DETECT');
        
        return result;
    }
    
    /**
     * Auto-correct detected anomalies by restoring toward baseline or applying constraints
     * @returns {Object} - { corrected: [...], totalCorrections: number }
     */
    autoCorrectAnomalies() {
        this.logAdaptation('Auto-correction of anomalies initiated', 'ANOMALY_CORRECT');
        
        const anomalies = this.detectAnomalies();
        const corrected = [];
        
        const baselineVars = this.calibrationData['_baseline']?.variables || {};
        
        for (const anom of anomalies.anomalies) {
            if (baselineVars.hasOwnProperty(anom.variable)) {
                // Blend current value toward baseline
                const blendFactor = 0.3;
                const correctedVal = 
                    this.variables[anom.variable] * (1 - blendFactor) + 
                    baselineVars[anom.variable] * blendFactor;
                
                this.variables[anom.variable] = correctedVal;
                corrected.push({
                    variable: anom.variable,
                    originalValue: anom.currentValue,
                    correctedValue: correctedVal,
                    blendFactor: blendFactor
                });
            }
        }
        
        const result = {
            corrected: corrected,
            totalCorrections: corrected.length,
            timestamp: Date.now()
        };
        
        this.logAdaptation(
            `Auto-correction applied: ${corrected.length} variables corrected`,
            'ANOMALY_CORRECT'
        );
        
        return result;
    }
```

**Impact**: ✅ New optimization, discovery, and anomaly detection APIs. Fully backward-compatible.

---

### 6. Custom Method Registration (NEW Methods)

**Location**: After adaptive operations methods

**Change Type**: Addition (2 new methods)

#### PATCH CONTENT:

```javascript
    // ===== ADAPTIVE: CUSTOM METHOD REGISTRATION (NEW) =====
    
    /**
     * Register a custom computation method (arbitrary function)
     * Enables runtime extension of computational capabilities
     * @param {string} name - Method identifier
     * @param {Function} fn - Function to store
     * @param {string} description - Human-readable description
     * @returns {boolean} - Success flag
     */
    addCustomMethod(name, fn, description = '') {
        if (typeof fn !== 'function') {
            console.error(`Custom method '${name}': argument must be a function.`);
            return false;
        }
        
        this.customMethods[name] = {
            fn: fn,
            description: description,
            registeredAt: Date.now()
        };
        
        this.logAdaptation(
            `Custom method added: ${name} - ${description}`,
            'METHOD_ADD'
        );
        return true;
    }
    
    /**
     * Execute a registered custom method
     * @param {string} name - Method identifier
     * @param {...any} args - Arguments to pass to the method
     * @returns {any} - Return value from custom method
     */
    executeCustomMethod(name, ...args) {
        if (!this.customMethods.hasOwnProperty(name)) {
            console.error(`Custom method '${name}' not found.`);
            this.logAdaptation(`Custom method execution failed: ${name} not found`, 'METHOD_ERROR');
            return null;
        }
        
        try {
            const result = this.customMethods[name].fn.apply(this, args);
            this.logAdaptation(`Custom method executed: ${name}`, 'METHOD_EXEC');
            return result;
        } catch (err) {
            console.error(`Error executing custom method '${name}':`, err);
            this.logAdaptation(`Custom method error: ${name}: ${err.message}`, 'METHOD_ERROR');
            return null;
        }
    }
```

**Impact**: ✅ New method registration APIs. Enables runtime extension. Fully backward-compatible.

---

### 7. Reporting & Utilities (NEW Methods)

**Location**: After custom methods

**Change Type**: Addition (5 new methods)

#### PATCH CONTENT:

```javascript
    // ===== ADAPTIVE: REPORTING & UTILITIES (NEW) =====
    
    /**
     * Generate comprehensive adaptive state report
     * @returns {string} - Formatted report text
     */
    generateReport() {
        const lines = [];
        lines.push('╔════════════════════════════════════════════════════════════════╗');
        lines.push('║       NGC346 ADAPTIVE UQFF MODULE - STATE REPORT               ║');
        lines.push('╚════════════════════════════════════════════════════════════════╝');
        lines.push('');
        
        lines.push('Variables Registered:', Object.keys(this.variables).length);
        lines.push('Physics Terms Registered:', Object.keys(this.physicsTerms).length);
        lines.push('Physics Terms Enabled:', Object.values(this.physicsTerms).filter(t => t.enabled).length);
        lines.push('Checkpoints Saved:', this.listCheckpoints().length);
        lines.push('Adaptation Log Entries:', this.adaptationLogs.length);
        lines.push('Custom Methods:', Object.keys(this.customMethods).length);
        lines.push('');
        
        lines.push('--- Physics Terms Status ---');
        for (const term of this.getPhysicsTerms()) {
            lines.push(`  ${term.enabled ? '✓' : '✗'} ${term.name}: ${term.description}`);
        }
        
        lines.push('');
        lines.push('--- Recent Adaptations ---');
        const recentLogs = this.adaptationLogs.slice(-5);
        for (const log of recentLogs) {
            lines.push(`  [${log.category}] ${log.message}`);
        }
        
        lines.push('');
        lines.push('--- Performance Metrics ---');
        lines.push(`  Computations: ${this.performanceMetrics.computations}`);
        lines.push(`  Total Time: ${(this.performanceMetrics.totalTime / 1000).toFixed(3)} s`);
        lines.push(`  Avg Time/Compute: ${(this.performanceMetrics.averageTime).toFixed(6)} ms`);
        
        return lines.join('\n');
    }
    
    /**
     * Export configuration as JSON-serializable object
     * @returns {Object} - Exportable configuration
     */
    exportConfiguration() {
        return {
            module: 'NGC346UQFFModule',
            timestamp: Date.now(),
            variables: this.getVariables(),
            physicsTermsRegistered: this.getPhysicsTerms(),
            checkpoints: this.listCheckpoints(),
            discoveredTerms: this.discoveredTerms,
            performanceMetrics: this.performanceMetrics,
            adaptationCount: this.adaptationLogs.length
        };
    }
    
    /**
     * Get audit trail of all adaptations
     * @returns {Array<Object>} - Adaptation log entries
     */
    getAdaptationLog() {
        return JSON.parse(JSON.stringify(this.adaptationLogs));
    }
    
    /**
     * Clamp parameter to physically reasonable range
     * (Helper for optimization)
     * @param {string} paramName - Parameter identifier
     * @param {number} value - Proposed value
     * @returns {number} - Clamped value
     */
    clampToPhysicalRange(paramName, value) {
        const ranges = {
            'M_visible': { min: 1e29, max: 1e35 },     // 0.005 to 50,000 M_sun
            'SFR': { min: 1e20, max: 1e25 },           // 0.001 to 100 M_sun/yr
            'rho_gas': { min: 1e-22, max: 1e-18 },     // kg/m³
            'B': { min: 1e-6, max: 1e-2 },             // Tesla
            'v_rad': { min: -100e3, max: 100e3 },      // m/s
            'r': { min: 1e15, max: 1e18 },             // m
            'z': { min: 0, max: 10 },                  // Redshift
            'T': { min: 1e4, max: 1e8 }                // K
        };
        
        if (ranges.hasOwnProperty(paramName)) {
            const range = ranges[paramName];
            return Math.max(range.min, Math.min(range.max, value));
        }
        
        return value;  // No range defined, return as-is
    }
    
    /**
     * Internal logging utility
     * @param {string} message - Log message
     * @param {string} category - Log category (INIT, OPTIMIZE, DISCOVER, etc.)
     */
    logAdaptation(message, category = 'INFO') {
        this.adaptationLogs.push({
            timestamp: Date.now(),
            category: category,
            message: message
        });
        
        // Keep log size bounded (~1000 entries)
        if (this.adaptationLogs.length > 1000) {
            this.adaptationLogs = this.adaptationLogs.slice(-500);
        }
    }
```

**Impact**: ✅ New reporting and utility methods. Fully backward-compatible.

---

### 8. Modify `computeG()` to Include Custom Terms (EXISTING METHOD MODIFICATION)

**Location**: Existing `computeG(t, r)` method (line ~265–335)

**Change Type**: Minimal modification (add 1 line to return statement)

#### BEFORE:
```javascript
    computeG(t, r) {
        this.variables['t'] = t;
        // ... 70 lines of existing computation ...
        const result = g_base + ug_sum + lambda_term + ui_term + quantum_term + fluid_term + dm_term;
        return result;
    }
```

#### AFTER:
```javascript
    computeG(t, r) {
        this.variables['t'] = t;
        // ... 70 lines of existing computation (unchanged) ...
        const result = g_base + ug_sum + lambda_term + ui_term + quantum_term + fluid_term + dm_term;
        
        // MODIFICATION: Add contribution from custom registered physics terms (NEW)
        const customTermsSum = this.computeCustomTerms(t, r);
        
        return result + customTermsSum;  // MODIFIED: was: return result;
    }
```

**Impact**: ✅ Minimal change. Custom terms contribution is additive and non-breaking. Existing behavior unchanged if no custom terms registered (custom sum = 0).

---

### 9. Adaptive Workflow Orchestration (NEW Method)

**Location**: After utilities section

**Change Type**: Addition (1 comprehensive orchestration method)

#### PATCH CONTENT:

```javascript
    /**
     * Full adaptive workflow orchestration
     * Combines optimization, discovery, anomaly detection, and correction
     * @param {number|Object} targetData - Target observational data
     * @param {Object} config - Workflow configuration
     *   - optimize: { enabled, iterations, learningRate }
     *   - discover: { enabled, scanRange }
     *   - detectAnomalies: { enabled, threshold }
     *   - correctAnomalies: { enabled }
     *   - mapQuantumResonance: { enabled, numStates }
     * @returns {Object} - Complete workflow result
     */
    adaptiveWorkflow(targetData, config = {}) {
        const startTime = Date.now();
        const workflow = {};
        
        // Default configuration
        const defaultConfig = {
            optimize: { enabled: true, iterations: 50, learningRate: 0.01 },
            discover: { enabled: true, scanRange: {} },
            detectAnomalies: { enabled: true, threshold: 2.0 },
            correctAnomalies: { enabled: true },
            mapQuantumResonance: { enabled: true, numStates: 26 }
        };
        
        const cfg = { ...defaultConfig, ...config };
        
        this.logAdaptation('Adaptive workflow started', 'WORKFLOW_START');
        
        // Step 1: Optimize parameters
        if (cfg.optimize.enabled) {
            workflow.optimization = this.optimizeParameters(
                targetData,
                cfg.optimize.iterations,
                cfg.optimize.learningRate
            );
        }
        
        // Step 2: Discover physics
        if (cfg.discover.enabled) {
            workflow.discovery = this.discoverPhysics(cfg.discover.scanRange);
        }
        
        // Step 3: Map quantum resonance
        if (cfg.mapQuantumResonance.enabled) {
            workflow.quantumResonance = this.mapQuantumResonance(cfg.mapQuantumResonance.numStates);
        }
        
        // Step 4: Detect anomalies
        if (cfg.detectAnomalies.enabled) {
            workflow.anomalies = this.detectAnomalies(cfg.detectAnomalies.threshold);
        }
        
        // Step 5: Auto-correct anomalies
        if (cfg.correctAnomalies.enabled && workflow.anomalies && workflow.anomalies.count > 0) {
            workflow.corrections = this.autoCorrectAnomalies();
        }
        
        const totalTime = Date.now() - startTime;
        
        workflow.summary = {
            completed: true,
            totalTime: totalTime,
            stepsExecuted: [
                cfg.optimize.enabled && 'optimize',
                cfg.discover.enabled && 'discover',
                cfg.mapQuantumResonance.enabled && 'quantumResonance',
                cfg.detectAnomalies.enabled && 'anomalyDetection',
                cfg.correctAnomalies.enabled && 'anomalyCorrection'
            ].filter(Boolean)
        };
        
        this.logAdaptation(
            `Adaptive workflow completed in ${totalTime}ms with ${workflow.summary.stepsExecuted.length} steps`,
            'WORKFLOW_COMPLETE'
        );
        
        return workflow;
    }
}

module.exports = NGC346UQFFModule;
```

**Impact**: ✅ Adds single high-level orchestration method. Fully backward-compatible.

---

## Summary of Changes

| Category | Type | Quantity | Lines Added | Backward Compatible? |
|----------|------|----------|-------------|----------------------|
| **Constructor Enhancement** | Addition | 1 method | ~15 | ✅ YES |
| **Variable Management** | Addition | 5 methods | ~80 | ✅ YES |
| **Physics Term Registry** | Addition | 6 methods | ~150 | ✅ YES |
| **State Management** | Addition | 5 methods | ~120 | ✅ YES |
| **Optimization & Discovery** | Addition | 8 methods | ~350 | ✅ YES |
| **Custom Methods** | Addition | 2 methods | ~40 | ✅ YES |
| **Reporting & Utilities** | Addition | 5 methods | ~100 | ✅ YES |
| **computeG() Modification** | Minimal | 1 line | ~2 | ✅ YES (custom terms additive) |
| **Workflow Orchestration** | Addition | 1 method | ~80 | ✅ YES |
| **TOTAL** | | **33 new methods** | **~937 lines** | ✅ **YES** |

**Final File Size**: `ngc346_uqff.js` grows from 448 → ~1,385 lines (308% increase)

---

## Testing Strategy

### Test Categories (Estimated 80+ Tests)

1. **Variable Management** (8 tests)
   - Add/get/remove variables
   - List operations
   - Error handling

2. **Physics Term Registration** (12 tests)
   - Register/unregister terms
   - Enable/disable operations
   - Compute custom terms sum
   - Error cases

3. **State Management** (10 tests)
   - Export/import state
   - Checkpoint save/load
   - State restoration with validation

4. **Parameter Optimization** (15 tests)
   - Convergence behavior
   - Learning rate effects
   - Gradient computation
   - Boundary conditions

5. **Physics Discovery** (10 tests)
   - Parameter sensitivity scanning
   - Relationship discovery
   - Edge cases

6. **Quantum Mapping** (8 tests)
   - Energy level mapping
   - Resonance peak detection
   - State count validation

7. **Anomaly Detection & Correction** (10 tests)
   - Threshold sensitivity
   - Anomaly identification
   - Correction application
   - Baseline handling

8. **Custom Methods** (5 tests)
   - Registration/execution
   - Error handling
   - Return value validation

9. **Adaptive Workflow** (8 tests)
   - Full workflow execution
   - Step enablement/disablement
   - Result aggregation
   - Timing validation

### Example Test File Structure

```javascript
// test_ngc346_adaptive.js
// ~900 lines, 80 comprehensive tests

class NGC346AdaptiveTestSuite {
    constructor() {
        this.module = new NGC346UQFFModule();
        this.testsPassed = 0;
        this.testsFailed = 0;
    }
    
    // Test categories as shown above
    testVariableManagement() { ... }
    testPhysicsTermRegistry() { ... }
    testStateManagement() { ... }
    testParameterOptimization() { ... }
    // ... etc
    
    run() { ... }
}
```

---

## Integration Plan

### Phase 1: Core Implementation (4 hours estimated)
1. Add enhanced constructor with adaptive system init
2. Implement variable management methods (5)
3. Implement physics term registry (6)
4. Implement state management (5)

### Phase 2: Advanced Features (6 hours estimated)
5. Implement optimization & discovery (8)
6. Implement custom methods (2)
7. Implement reporting & utilities (5)
8. Add workflow orchestration (1)

### Phase 3: Integration & Testing (4 hours estimated)
9. Modify `computeG()` to call `computeCustomTerms()`
10. Create comprehensive test suite (~80 tests)
11. Run full test suite → verify 100% pass rate
12. Update index.js to export enhanced NGC346UQFFModule

### Phase 4: Documentation & Verification (2 hours estimated)
13. Create NGC346 adaptive integration report
14. Run deployment verification script
15. Create deployment log

**Total Estimated Effort**: 16 hours  
**Risk Level**: Low (non-breaking, fully backward-compatible)

---

## Risks & Mitigation

| Risk | Probability | Mitigation |
|------|-------------|-----------|
| Custom term functions cause performance degradation | Low | Optimize `computeCustomTerms()` loop; profile before/after |
| State serialization misses edge cases | Low | Comprehensive testing of edge cases (null, undefined, large numbers) |
| Physics discovery scan is slow | Medium | Implement early stopping; add timeout capability |
| Gradient descent doesn't converge | Low | Provide convergence monitoring; adaptive step size |

---

## Expected Outcomes

After implementing this patch:

✅ **NGC346UQFFModule feature parity with SMBHMSRAdaptiveModule**
- Same core APIs for adaptive operations
- Same state serialization patterns
- Same testing methodology

✅ **Production-Ready Dynamic Capability**
- Runtime parameter optimization
- Physics term registration and management
- State snapshots and restoration
- Anomaly detection and correction

✅ **100% Test Coverage**
- All existing tests continue to pass (backward compat)
- 80+ new adaptive tests added
- Comprehensive integration verification

✅ **Framework Consistency**
- S77-S82 now all support adaptive capabilities
- Unified API across all systems
- Extensible for future systems (S83+)

---

## Next Steps

**Approval Required For**:
1. ✅ Proceed with implementation of all patches as documented
2. ✅ Add comprehensive test suite (80+ tests)
3. ✅ Update index.js exports
4. ✅ Create final deployment documentation

**Ready to Proceed?** → Confirm "start implementation" or request modifications.

