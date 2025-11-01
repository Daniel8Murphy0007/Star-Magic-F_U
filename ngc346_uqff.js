// ngc346_uqff.js
// NGC 346 Young Stellar Cluster Nebula - Protostar Formation & Cluster Entanglement UQFF Module
// Port of Source81.cpp - Models young stellar cluster evolution, protostar collapse dynamics, 
// cluster entanglement via Ug forces, blueshifted quantum waves, and multi-scale physics
//
// System: NGC 346 (Small Magellanic Cloud, ~60 kpc away)
// Mass: 1200 M☉ (1000 visible + 200 dark matter)
// Radius: 5 pc = 1.543e17 m
// Star Formation Rate: 0.1 M☉/year
// Observation Timescale: 10 Myr
// Key Physics: Protostar formation (Ug3 collapse dominant), cluster entanglement, blueshifted waves
// 57 dynamic variables, 15 methods (14 private computation + master equation + utilities)

class NGC346UQFFModule {
    constructor() {
        this.variables = {};
        
        // ===== UNIVERSAL CONSTANTS (11) =====
        this.variables['G'] = 6.6743e-11;                    // m³ kg⁻¹ s⁻²
        this.variables['c'] = 3e8;                           // m/s
        this.variables['hbar'] = 1.0546e-34;                 // J·s
        this.variables['Lambda'] = 1.1e-52;                  // m⁻²
        this.variables['q'] = 1.602e-19;                     // C (elementary charge)
        this.variables['pi'] = 3.141592653589793;
        this.variables['t_Hubble'] = 13.8e9 * 3.156e7;       // s (~4.35e17 s)
        this.variables['year_to_s'] = 3.156e7;               // s/yr
        this.variables['H0'] = 70.0;                         // km/s/Mpc
        this.variables['Mpc_to_m'] = 3.086e22;               // m/Mpc
        this.variables['mu_0'] = 4 * this.variables['pi'] * 1e-7; // H/m
        
        // ===== NGC 346 SYSTEM PARAMETERS (11) =====
        const M_sun = 1.989e30;                              // kg
        const pc = 3.086e16;                                 // m
        
        this.variables['M_visible'] = 1000 * M_sun;          // kg
        this.variables['M_DM'] = 200 * M_sun;                // kg
        this.variables['M'] = this.variables['M_visible'] + this.variables['M_DM']; // Total
        this.variables['M0'] = this.variables['M'];           // Reference mass
        this.variables['SFR'] = 0.1 * M_sun / this.variables['year_to_s']; // kg/s
        this.variables['r'] = 5 * pc;                        // m
        this.variables['z'] = 0.0006;                        // Redshift (SMC)
        this.variables['rho_gas'] = 1e-20;                   // kg/m³
        this.variables['v_rad'] = -10e3;                     // m/s (blueshift)
        this.variables['t'] = 1e7 * this.variables['year_to_s']; // Default t=10 Myr (s)
        this.variables['V'] = 1e49;                          // m³ (volume)
        
        // ===== MAGNETIC & ELECTROMAGNETIC (4) =====
        this.variables['B'] = 1e-5;                          // Tesla
        this.variables['B_crit'] = 1e11;                     // Tesla (MHD instability threshold)
        this.variables['H_aether'] = 1e-6;                   // A/m (aether field strength)
        
        // ===== QUANTUM WAVE PARAMETERS (7) =====
        this.variables['Delta_x'] = 1e-10;                   // m
        this.variables['Delta_p'] = this.variables['hbar'] / this.variables['Delta_x']; // kg·m/s
        this.variables['integral_psi'] = 1.0;                // Normalized wavefunction
        this.variables['A'] = 1e-10;                         // Wave amplitude
        this.variables['k'] = 1e20;                          // Wavenumber (m⁻¹)
        this.variables['omega'] = 1e-14;                     // Angular frequency (rad/s)
        this.variables['sigma'] = 1e16;                      // Gaussian width (m)
        
        // ===== GRAVITATIONAL FORCE COMPONENTS (8) =====
        this.variables['Ug1'] = 0.0;                         // Dipole
        this.variables['Ug2'] = 0.0;                         // Superconductor
        this.variables['Ug3'] = 0.0;                         // Collapse (dominant)
        this.variables['Ug4'] = 0.0;                         // Reaction
        this.variables['Ui'] = 0.0;                          // Universal inertia
        this.variables['Um'] = 0.0;                          // Universal magnetism
        this.variables['rho_vac_UA'] = 7.09e-36;             // J/m³
        this.variables['lambda_I'] = 1.0;                    // Integral scale
        
        // ===== OSCILLATORY & SCALE PARAMETERS (12) =====
        this.variables['x'] = 0.0;                           // Position
        this.variables['v'] = Math.abs(this.variables['v_rad']); // Velocity magnitude (m/s)
        this.variables['omega_i'] = 1e-8;                    // Inertia frequency (rad/s)
        this.variables['t_n'] = 0.0;                         // Time phase
        this.variables['F_RZ'] = 0.01;                       // Time-reversal factor
        this.variables['k_4'] = 1.0;                         // Reaction amplitude coefficient
        this.variables['k_SF'] = 1e-10;                      // Star formation coupling
        this.variables['scale_macro'] = 1e-12;               // Macro scale
        this.variables['f_TRZ'] = 0.1;                       // Time-reversal factor (large)
        this.variables['f_sc'] = 1.0;                        // Scale factor
        this.variables['v_r'] = 1e3;                         // Radial velocity (m/s)
        this.variables['rho'] = this.variables['rho_gas'];   // Density reference
        
        // ===== COSMOLOGICAL (3) =====
        this.variables['Omega_m'] = 0.3;                     // Matter density parameter
        this.variables['Omega_Lambda'] = 0.7;                // Dark energy density parameter
        this.variables['delta_rho_over_rho'] = 1e-5;         // Density perturbation ratio
        
        // ===== ADAPTIVE SYSTEM INITIALIZATION (NEW) =====
        this.physicsTerms = {};           // Registry: { termName: { fn, enabled, description, registeredAt } }
        this.calibrationData = {};        // Snapshots: { label: { timestamp, variables, metadata } }
        this.discoveredTerms = {};        // New physics: { termName: { relationship, sensitivity } }
        this.adaptationLogs = [];         // Audit trail: [{ timestamp, category, message }]
        this.customMethods = {};          // User functions: { methodName: { fn, description, registeredAt } }
        this.quantumMap = {};             // Resonance: { stateIndex: { energy, amplitude, ...} }
        this.performanceMetrics = {       // Statistics
            computations: 0,
            totalTime: 0,
            averageTime: 0
        };
        
        this.initializeAdaptiveSystem();  // Initialize calibration baseline
    }
    
    // ===== ADAPTIVE SYSTEM INITIALIZATION (NEW) =====
    
    initializeAdaptiveSystem() {
        // Store baseline calibration at construction time
        this.calibrationData['_baseline'] = {
            timestamp: Date.now(),
            variables: JSON.parse(JSON.stringify(this.variables)),
            config: 'initialization'
        };
        
        this.logAdaptation('Adaptive system initialized (57 variables, 14 compute methods)', 'INIT');
    }
    
    // ===== PRIVATE COMPUTATION METHODS (14) =====
    
    // Hubble parameter at redshift z
    computeHtz(z_val) {
        const Hz_kms = this.variables['H0'] * Math.sqrt(
            this.variables['Omega_m'] * Math.pow(1.0 + z_val, 3) + 
            this.variables['Omega_Lambda']
        );
        return (Hz_kms * 1e3) / this.variables['Mpc_to_m'];
    }
    
    // Mass evolution from star formation
    computeMsfFactor(t) {
        return (this.variables['SFR'] * t) / this.variables['M0'];
    }
    
    // Radius evolution (expansion/contraction)
    computeRt(t) {
        return this.variables['r'] + this.variables['v_r'] * t;
    }
    
    // Environmental force (collapse + star formation)
    computeFenv(t) {
        const F_collapse = this.variables['rho_gas'] * Math.pow(this.variables['v_rad'], 2);
        const F_SF = this.variables['k_SF'] * this.variables['SFR'] / 1.989e30;
        return F_collapse + F_SF;
    }
    
    // Ug1: Dipole component
    computeUg1(t) {
        return 1e-10 * Math.cos(this.variables['omega'] * t);
    }
    
    // Ug2: Superconductor/magnetic energy density
    computeUg2(t) {
        const B_super = this.variables['mu_0'] * this.variables['H_aether'];
        return (B_super * B_super) / (2 * this.variables['mu_0']);
    }
    
    // Ug3: COLLAPSE MECHANISM - Magnetic strings disk (DOMINANT TERM)
    // This drives protostar formation
    computeUg3(t) {
        const rho_vac = this.variables['rho_vac_UA'];
        return this.variables['G'] * this.variables['M'] / 
               (this.variables['r'] * this.variables['r']) * 
               (this.variables['rho_gas'] / rho_vac);
    }
    
    // Ug4: Reaction/expansion component
    computeUg4(t) {
        const E_react = 1e40 * Math.exp(-0.0005 * t);
        return this.variables['k_4'] * E_react;
    }
    
    // Ui: Universal inertia
    computeUi(t) {
        return this.variables['lambda_I'] * 
               (this.variables['rho_vac_UA'] / 1e-9) * 
               this.variables['omega_i'] * 
               Math.cos(this.variables['pi'] * this.variables['t_n']);
    }
    
    // Um: Universal magnetism (Lorentz force)
    computeUm(t) {
        return this.variables['q'] * this.variables['v_rad'] * this.variables['B'];
    }
    
    // Psi integral: Quantum wavefunction intensity |ψ|²
    computePsiIntegral(r, t) {
        const A = this.variables['A'];
        const sigma = this.variables['sigma'];
        const omega = this.variables['omega'];
        
        // |ψ|² = A² exp(-r²/σ²)
        // (Gaussian envelope, time oscillation encoded separately)
        const envelope = Math.exp(-(r * r) / (sigma * sigma));
        return A * A * envelope;
    }
    
    // Quantum term: Bridges quantum mechanics with cosmology
    computeQuantumTerm(t_Hubble_val, r) {
        const unc = Math.sqrt(this.variables['Delta_x'] * this.variables['Delta_p']);
        const psi_int = this.computePsiIntegral(r, this.variables['t']);
        return (this.variables['hbar'] / unc) * 
               this.variables['integral_psi'] * 
               (2 * this.variables['pi'] / t_Hubble_val) * 
               psi_int;
    }
    
    // Fluid term: Hydrostatic pressure correction
    computeFluidTerm(g_base) {
        return this.variables['rho_gas'] * this.variables['V'] * g_base;
    }
    
    // DM term: Dark matter perturbation and curvature
    computeDMTerm(r) {
        const pert = this.variables['delta_rho_over_rho'];
        const curv = 3 * this.variables['G'] * this.variables['M'] / (r * r * r);
        return (this.variables['M_visible'] + this.variables['M_DM']) * (pert + curv);
    }
    
    // Total gravitational component (cluster entanglement sum)
    computeUgSum(r) {
        const Ug_base = (this.variables['G'] * this.variables['M']) / (r * r);
        this.variables['Ug1'] = this.computeUg1(this.variables['t']);
        this.variables['Ug2'] = this.computeUg2(this.variables['t']);
        this.variables['Ug3'] = this.computeUg3(this.variables['t']);
        this.variables['Ug4'] = this.computeUg4(this.variables['t']);
        const um = this.computeUm(this.variables['t']);
        return Ug_base + this.variables['Ug1'] + this.variables['Ug2'] + 
               this.variables['Ug3'] + this.variables['Ug4'] + um;
    }
    
    // ===== HELPER METHODS =====
    
    // Core energy from collapse
    computeEcore(rho) {
        const ug3 = this.computeUg3(this.variables['t']);
        const ui = this.computeUi(this.variables['t']);
        return ug3 + ui * rho;
    }
    
    // Core temperature
    computeTempCore(ug3) {
        const rho_vac = this.variables['rho_vac_UA'];
        return 1.424e7 * (ug3 * rho_vac);  // Scaled to K
    }
    
    // ===== PUBLIC INTERFACE METHODS =====
    
    // Dynamic variable operations
    updateVariable(name, value) {
        if (name in this.variables) {
            this.variables[name] = value;
        } else {
            this.variables[name] = value;
        }
        
        // Update dependent variables
        if (name === 'Delta_x') {
            this.variables['Delta_p'] = this.variables['hbar'] / value;
        } else if (name === 'M_visible' || name === 'M_DM') {
            this.variables['M'] = this.variables['M_visible'] + this.variables['M_DM'];
            this.variables['M0'] = this.variables['M'];
        } else if (name === 'rho_gas') {
            this.variables['rho'] = value;
        }
    }
    
    addToVariable(name, delta) {
        if (name in this.variables) {
            this.variables[name] += delta;
        } else {
            this.variables[name] = delta;
        }
    }
    
    subtractFromVariable(name, delta) {
        this.addToVariable(name, -delta);
    }
    
    getVariable(name) {
        return this.variables[name];
    }
    
    getState() {
        return Object.assign({}, this.variables);
    }
    
    setState(state_object) {
        this.variables = Object.assign({}, state_object);
    }
    
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
    getVariableValue(name) {
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
    
    // ===== ADAPTIVE: PHYSICS TERM REGISTRATION (NEW) =====
    
    /**
     * Register a custom physics term that is automatically included in computations
     * @param {string} name - Term identifier (e.g., 'custom_dark_energy_correction')
     * @param {Function} computeFn - Computation function(t, r, variables) => number
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
    
    // ===== ADAPTIVE: STATE MANAGEMENT (NEW) =====
    
    /**
     * Export complete state (variables, physics terms, calibration metadata)
     * @returns {Object} - Serializable state snapshot
     */
    exportState() {
        return {
            variables: JSON.parse(JSON.stringify(this.variables)),
            physicsTerms: Object.entries(this.physicsTerms).map(([name, term]) => ({
                name: name,
                enabled: term.enabled,
                description: term.description,
                registeredAt: term.registeredAt
            })),
            calibrationKeys: Object.keys(this.calibrationData),
            discoveredTerms: this.discoveredTerms,
            timestamp: Date.now()
        };
    }
    
    /**
     * Import/restore state from a snapshot
     * Note: Physics term functions must be re-registered separately
     * @param {Object} state - State object from exportState()
     * @returns {boolean} - Success flag
     */
    importState(state) {
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
    
    /**
     * Full adaptive workflow orchestration
     * Combines optimization, discovery, anomaly detection, and correction
     * @param {number|Object} targetData - Target observational data
     * @param {Object} config - Workflow configuration
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
        
        lines.push(`Variables Registered: ${Object.keys(this.variables).length}`);
        lines.push(`Physics Terms Registered: ${Object.keys(this.physicsTerms).length}`);
        lines.push(`Physics Terms Enabled: ${Object.values(this.physicsTerms).filter(t => t.enabled).length}`);
        lines.push(`Checkpoints Saved: ${this.listCheckpoints().length}`);
        lines.push(`Adaptation Log Entries: ${this.adaptationLogs.length}`);
        lines.push(`Custom Methods: ${Object.keys(this.customMethods).length}`);
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
    
    // ===== MASTER UQFF EQUATION FOR NGC 346 =====
    // g_NGC346(r, t) = Base gravity with expansions/corrections + all multi-physics terms
    computeG(t, r) {
        this.variables['t'] = t;
        
        if (r === 0) r = this.variables['r'];  // Use default radius if not specified
        
        // Mass evolution from star formation
        const msf_factor = this.computeMsfFactor(t);
        const m_factor = 1.0 + msf_factor;
        
        // Cosmological expansion
        const Hz = this.computeHtz(this.variables['z']);
        const expansion = 1.0 + Hz * t;
        
        // Magnetic suppression (usually ~1)
        const sc_correction = 1.0 - (this.variables['B'] / this.variables['B_crit']);
        
        // Environmental feedback
        const f_env = this.computeFenv(t);
        
        // Time-reversal factor
        const tr_factor = 1.0 + this.variables['f_TRZ'];
        
        // Base gravity with all corrections
        const g_base = (this.variables['G'] * this.variables['M'] * m_factor / (r * r)) * 
                      expansion * sc_correction * (1.0 + f_env) * tr_factor;
        
        // Gravitational components (subtract base to avoid double-count)
        const ug_sum = this.computeUgSum(r) - g_base;
        
        // Cosmological constant term
        const lambda_term = this.variables['Lambda'] * (this.variables['c'] * this.variables['c']) / 3.0;
        
        // Universal inertia
        const ui_term = this.computeUi(t);
        
        // Quantum term
        const quantum_term = this.computeQuantumTerm(this.variables['t_Hubble'], r);
        
        // Fluid term
        const fluid_term = this.computeFluidTerm(g_base);
        
        // Dark matter term
        const dm_term = this.computeDMTerm(r);
        
        // Total acceleration
        const result = g_base + ug_sum + lambda_term + ui_term + quantum_term + fluid_term + dm_term;
        
        // MODIFICATION (Phase 2C): Add contribution from custom registered physics terms
        const customTermsSum = this.computeCustomTerms(t, r);
        
        return result + customTermsSum;
    }
    
    // ===== OUTPUT & DOCUMENTATION =====
    
    getAllFrequencies(t, r) {
        this.variables['t'] = t;
        
        return {
            'Ug1': this.computeUg1(t),
            'Ug2': this.computeUg2(t),
            'Ug3': this.computeUg3(t),
            'Ug4': this.computeUg4(t),
            'Ui': this.computeUi(t),
            'Um': this.computeUm(t),
            'F_env': this.computeFenv(t),
            'Ecore': this.computeEcore(this.variables['rho_gas']),
            'Tcore': this.computeTempCore(this.computeUg3(t))
        };
    }
    
    getEquationText() {
        return `
NGC 346 UQFF EQUATION (Source81)
═════════════════════════════════════════════════════════════

g_NGC346(r, t) = (G·M(t)/r²)·(1 + H(t,z))·(1 - B/B_crit)·(1 + F_env(t))·(1 + f_TRZ)
                + [Ug1 + Ug2 + Ug3 + Ug4 + Um] - g_base
                + (Λ·c²/3) + Ui + Q_term + F_fluid + U_DM

COMPONENT DESCRIPTIONS:

Base Gravity:
  g_base = (G·M(t)/r²)·(1 + M_SF(t)/M₀)
  M(t) evolves via star formation: M_SF(t) = SFR·t

Corrections:
  (1 + H(t,z)) - Cosmological expansion (Hubble parameter)
  (1 - B/B_crit) - Magnetic field suppression (~1.0 for weak fields)
  (1 + F_env(t)) - Environmental feedback (collapse + star formation)
  (1 + f_TRZ) - Time-reversal factor (~1.1)

Gravitational Components (Cluster Entanglement):
  Ug1 = 10⁻¹⁰·cos(ω·t) - Dipole oscillation
  Ug2 = B_super²/(2μ₀) - Superconductor energy (~10⁻¹⁸ m/s²)
  Ug3 = (G·M/r²)·(ρ_gas/ρ_vac,UA) - COLLAPSE DOMINANT (~10¹¹ m/s²!)
  Ug4 = k₄·E_react(t) - Reaction/expansion (exp decay)
  Um = q·v_rad·B - Lorentz magnetism (~10⁻¹⁸ m/s²)

Multi-Physics Terms:
  Ui = λ_I·(ρ_vac/ρ_plasm)·ω_i·cos(π·t_n) - Inertial coupling
  Q_term = (ℏ/√(Δx·Δp))·ψ_integral·(2π/t_Hubble)·|ψ|² - Quantum correction
  F_fluid = ρ_gas·V·g_base - Hydrostatic pressure
  U_DM = (M_vis + M_DM)·(Δρ/ρ + 3GM/r³) - Dark matter perturbation

Cosmological:
  Λ·c²/3 - Cosmological constant term

KEY PHYSICS PHENOMENA:
═════════════════════════════════════════════════════════════

1. PROTOSTAR COLLAPSE: Ug3 term dominates (~1e11 m/s²)
   - Drives gravitational instability and star formation
   - Gas density enhancement via ρ_gas/ρ_vac ratio

2. CLUSTER ENTANGLEMENT: Multi-component Ug coupling
   - Sum of Ug1-Ug4 encodes system-wide interactions
   - Non-local effects via Ui (inertial coupling)

3. BLUESHIFTED QUANTUM WAVES: v_rad = -10 km/s (approaching)
   - Wavefunction: ψ(r,t) = A·exp(-r²/σ²)·exp(i(kx - ωt))
   - σ = 1e16 m (large-scale coherence)
   - ω = 1e-14 rad/s (slow modulation, period ~600 Myr)

4. MULTI-TIMESCALE EVOLUTION:
   - Wave timescale: 2π/ω ≈ 600 Myr (very slow)
   - Collapse timescale: ~0.1-1 Myr (protostar formation)
   - Cluster lifetime: ~100 Myr (dispersal)
   - SFR timescale: M/SFR ≈ 12 Myr (mass consumption)

5. REAL NGC 346 PROPERTIES:
   - Location: Small Magellanic Cloud (~60 kpc)
   - Age: ~2-3 Myr (young cluster)
   - SFR: 0.1 M☉/yr (active star formation)
   - ~350 embedded protostars across ~4 pc
   - Redshift z = 0.0006 (distance ~8 Mly)
        `;
    }
    
    printVariables() {
        console.log('NGC 346 UQFF Module Variables:');
        console.log('═'.repeat(60));
        const sorted = Object.keys(this.variables).sort();
        for (const name of sorted) {
            const value = this.variables[name];
            const sci = typeof value === 'number' && Math.abs(value) > 1e5 ? 
                       value.toExponential(6) : value.toFixed(6);
            console.log(`  ${name.padEnd(20)} = ${sci}`);
        }
        console.log('═'.repeat(60));
    }
    
    // Get multi-timescale evolution (similar to coalescence tracking in S80)
    getCollapseEvolution(num_points) {
        const t_max = this.variables['t'];
        const evolution = [];
        
        for (let i = 0; i <= num_points; i++) {
            const t_frac = i / num_points;
            const t_sample = t_frac * t_max;
            const r_sample = this.variables['r'] + this.variables['v_r'] * t_sample;
            
            // Compute state at this point
            const ug1 = this.computeUg1(t_sample);
            const ug3 = this.computeUg3(t_sample);
            const ug4 = this.computeUg4(t_sample);
            const f_env = this.computeFenv(t_sample);
            const acceleration = this.computeG(t_sample, r_sample);
            
            evolution.push({
                time_fraction: t_frac,
                time_seconds: t_sample,
                time_to_end: t_max - t_sample,
                radius_m: r_sample,
                Ug1: ug1,
                Ug3_collapse: ug3,
                Ug4_reaction: ug4,
                F_env: f_env,
                acceleration: acceleration
            });
        }
        
        return evolution;
    }
}

module.exports = NGC346UQFFModule;
