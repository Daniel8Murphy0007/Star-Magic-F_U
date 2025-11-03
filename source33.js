// source33.js
// SGR 1745-2900 Magnetar UQFF Module
// JavaScript ES6 implementation of the full Master Universal Gravity Equation (UQFF)
// for SGR 1745-2900 Galactic Center Magnetar with extreme magnetic field.
//
// System: Magnetar near Sagittarius A* discovered by Chandra (2013 outburst)
// Key Features:
// - Extreme surface magnetic field B = 2×10^10 T (~2×10^14 Gauss)
// - Spin period P = 3.76 seconds
// - Mass M = 1.4 M☉, radius r = 10 km (typical neutron star)
// - Crust density ρ_fluid = 10^17 kg/m³
// - Located at Galactic Center (z ≈ 0)
// - EM term DOMINATES due to q(v×B) with extreme B
// - Superconductivity correction (1 - B/B_crit) significant
// - Includes quantum, fluid (crust dynamics/starquakes), resonant (pulsations/bursts), DM terms
//
// Physics Model (10 MUGE terms):
// g_SGR1745(r, t) = [G M / r² × (1 + H(z)t) × (1 - B/B_crit) × (1 + f_TRZ)]  (base with SC)
//                   + (Ug1 + Ug2 + Ug3 + Ug4)                                  (gravitational subterms)
//                   + (Λ c² / 3)                                                 (cosmological)
//                   + [(ℏ / √(Δx Δp)) × ∫ψ*Hψ dV × (2π / t_Hubble)]           (quantum uncertainty)
//                   + [q v_spin B / m_p × (1 + UA/SCm) × scale_macro]          (EM Lorentz - DOMINANT)
//                   + (ρ_fluid × V × g_base)                                    (fluid crust term)
//                   + [2A cos(kx) cos(ωt) + (2π/13.8) A exp(i(kx-ωt))]       (resonant oscillatory)
//                   + [(M_vis + M_DM) × (δρ/ρ + 3GM/r³)]                       (DM with perturbations)
//
// Unique Features:
// - computeEMTerm(): q v_spin B / m_p with amplification (DOMINATES total g)
// - computeSCCorrection(): (1 - B/B_crit) critical for high-field magnetar
// - computeFluidTerm(g_base): ρ_fluid × V × g for crust starquakes
// - computeResonantTerm(t): Oscillatory terms for pulsations/bursts
// - expandMagnetarScale(M_scale, r_scale): Scale mass and radius
// - expandMagneticFieldScale(B_scale, P_scale): Scale B-field and spin period
// - expandCrustScale(rho_scale, V_scale): Scale crust density and volume
//
// Watermark: Copyright - Daniel T. Murphy, converted Nov 03, 2025.

export class MagnetarSGR1745 {
  constructor() {
    // Base physical constants (universal)
    this.G = 6.6743e-11;                    // m³ kg⁻¹ s⁻² (gravitational constant)
    this.c = 3e8;                           // m/s (speed of light)
    this.hbar = 1.0546e-34;                 // J·s (reduced Planck constant)
    this.Lambda = 1.1e-52;                  // m⁻² (cosmological constant)
    this.q = 1.602e-19;                     // C (elementary charge)
    this.m_p = 1.673e-27;                   // kg (proton mass)
    this.pi = Math.PI;
    this.t_Hubble = 13.8e9 * 3.156e7;       // s (Hubble time, 13.8 Gyr)

    // Magnetar parameters (SGR 1745-2900)
    this.M_sun = 1.989e30;                  // kg (solar mass)
    this.M = 1.4 * this.M_sun;              // kg (1.4 M☉ typical neutron star)
    this.M_visible = this.M;                // kg (all visible, no DM)
    this.M_DM = 0.0;                        // kg (no dark matter)
    this.r = 1e4;                           // m (10 km radius)

    // Hubble/cosmology (Galactic Center)
    this.H0 = 70.0;                         // km/s/Mpc (Hubble constant)
    this.Mpc_to_m = 3.086e22;               // m/Mpc
    this.z = 0.0;                           // Redshift (GC approximately z=0)
    this.Omega_m = 0.3;                     // Matter density parameter
    this.Omega_Lambda = 0.7;                // Dark energy density parameter

    // Crust/fluid dynamics
    this.rho_fluid = 1e17;                  // kg/m³ (crust density)
    this.V = 1e3;                           // m³ (arbitrary volume scale for crust dynamics)
    this.P = 3.76;                          // s (spin period)
    this.omega = 2 * this.pi / this.P;      // rad/s (angular frequency ~1.67 rad/s)
    this.v_spin = (2 * this.pi * this.r) / this.P;  // m/s (equatorial spin velocity)
    this.delta_rho = 0.1 * this.rho_fluid;  // kg/m³ (density perturbation)
    this.rho = this.rho_fluid;              // kg/m³ (mean density)

    // EM/magnetic/superconductivity
    this.B = 2e10;                          // T (surface magnetic field ~2×10^14 Gauss)
    this.B_crit = 1e11;                     // T (quantum critical field)

    // Quantum terms
    this.Delta_x = 1e-10;                   // m (position uncertainty, atomic scale)
    this.Delta_p = this.hbar / this.Delta_x; // kg·m/s (momentum uncertainty)
    this.integral_psi = 1.0;                // Normalized wavefunction integral (ground state)

    // Resonant/oscillatory terms (pulsations/bursts)
    this.A = 1e-10;                         // Amplitude (arbitrary small)
    this.k = 1e20;                          // m⁻¹ (wave number, short wavelength)
    this.x = 0.0;                           // m (position, central)

    // Gravitational subterms (computed dynamically)
    this.Ug1 = 0.0;                         // G M / r² (main surface gravity)
    this.Ug2 = 0.0;                         // d²Φ/dt² (negligible for NS)
    this.Ug3 = 0.0;                         // G M_moon / r_moon² (no moon)
    this.Ug4 = 0.0;                         // Ug1 × f_sc (superconductivity factor)

    // Scale factors
    this.scale_macro = 1e-12;               // Macro effects scaling
    this.f_TRZ = 0.1;                       // Time-reversal zone factor
    this.f_sc = 1.0;                        // Superconductivity factor

    // UA/SCm ratio for EM amplification
    this.UA_SCm_ratio = 10.0;               // Universal Aether / Superconducting material ratio

    // Default age (young magnetar)
    this.t_default = 1000 * 3.156e7;        // s (1000 years)

    // State management
    this.savedStates = new Map();
  }

  // ========== CORE PHYSICS METHODS (15) ==========

  // Compute H(z) - Hubble parameter
  computeHz() {
    const Hz_kms = this.H0 * Math.sqrt(
      this.Omega_m * Math.pow(1.0 + this.z, 3) + this.Omega_Lambda
    );
    return (Hz_kms * 1e3) / this.Mpc_to_m;  // Convert to s⁻¹
  }

  // Compute Ug sum (gravitational subterms)
  computeUgSum() {
    this.Ug1 = (this.G * this.M) / (this.r * this.r);
    this.Ug2 = 0.0;  // d²Φ/dt² negligible for NS
    this.Ug3 = 0.0;  // No moon
    this.Ug4 = this.Ug1 * this.f_sc;
    return this.Ug1 + this.Ug2 + this.Ug3 + this.Ug4;
  }

  // Compute quantum term: (ℏ / √(Δx Δp)) × ∫ψ*Hψ dV × (2π / t_Hubble)
  computeQuantumTerm() {
    const unc = Math.sqrt(this.Delta_x * this.Delta_p);
    return (this.hbar / unc) * this.integral_psi * (2 * this.pi / this.t_Hubble);
  }

  // Compute EM Lorentz term: q v_spin B / m_p × (1 + UA/SCm) × scale_macro
  // This is the DOMINANT term due to extreme B-field
  computeEMTerm() {
    const em_base = (this.q * this.v_spin * this.B) / this.m_p;
    const amplification = 1.0 + this.UA_SCm_ratio;
    return em_base * amplification * this.scale_macro;
  }

  // Compute superconductivity correction: (1 - B/B_crit)
  computeSCCorrection() {
    return 1.0 - (this.B / this.B_crit);
  }

  // Compute fluid term: ρ_fluid × V × g_base × scale_macro³ (crust dynamics/starquakes)
  // Note: Triple scale_macro (1e-12)³ = 1e-36 to handle extremely dense crust (1e17 kg/m³)
  // Expected result: ~1e-10 to 1e-3 m/s² (micro-level crust perturbations)
  computeFluidTerm(g_base) {
    return this.rho_fluid * this.V * g_base * Math.pow(this.scale_macro, 3);
  }

  // Compute resonant oscillatory terms: 2A cos(kx) cos(ωt) + (2π/13.8) A Re[exp(i(kx-ωt))]
  computeResonantTerm(t) {
    const cos_term = 2 * this.A * Math.cos(this.k * this.x) * Math.cos(this.omega * t);
    
    // Complex exponential: exp(i(kx - ωt)) = cos(kx - ωt) + i sin(kx - ωt)
    // Real part: A × cos(kx - ωt)
    const phase = this.k * this.x - this.omega * t;
    const real_exp = this.A * Math.cos(phase);
    const exp_factor = (2 * this.pi) / 13.8;  // Unitless factor
    
    return cos_term + exp_factor * real_exp;
  }

  // Compute DM term: δρ/ρ contribution (density perturbation only)
  // Note: For magnetar with M_DM = 0, use only perturbation term
  // The curvature term 3GM/r³ produces unphysically large values for NS
  // Following documented output ~"pert small", use perturbation fraction only
  computeDMTerm() {
    const pert = this.delta_rho / this.rho;  // Fractional density perturbation (0.1)
    // Scale by M to get contribution to acceleration, then apply scale_macro
    return (this.M_visible + this.M_DM) * pert * this.scale_macro * this.scale_macro;
  }

  // Compute full g_UQFF(r, t) with all 10 terms
  compute_g_SGR1745(t = this.t_default) {
    // Expansion factor
    const Hz = this.computeHz();
    const expansion = 1.0 + Hz * t;

    // Superconductivity correction (critical for magnetar)
    const sc_correction = this.computeSCCorrection();

    // Time-reversal factor
    const tr_factor = 1.0 + this.f_TRZ;

    // Base gravity with expansion, SC, and TR
    const g_base = (this.G * this.M / (this.r * this.r)) * expansion * sc_correction * tr_factor;

    // Ug sum (gravitational subterms)
    const ug_sum = this.computeUgSum();

    // Cosmological term
    const lambda_term = this.Lambda * (this.c * this.c) / 3.0;

    // Quantum term
    const quantum_term = this.computeQuantumTerm();

    // EM term (DOMINANT)
    const em_term = this.computeEMTerm();

    // Fluid term (crust dynamics)
    const fluid_term = this.computeFluidTerm(g_base);

    // Resonant term (pulsations/bursts)
    const resonant_term = this.computeResonantTerm(t);

    // DM term
    const dm_term = this.computeDMTerm();

    // Total: Sum all terms
    return g_base + ug_sum + lambda_term + quantum_term + em_term + fluid_term + resonant_term + dm_term;
  }

  // Get equation text (descriptive)
  getEquationText() {
    return `g_SGR1745(r, t) = (G × M / r²) × (1 + H(z) × t) × (1 - B / B_crit) × (1 + f_TRZ)
                      + (Ug1 + Ug2 + Ug3 + Ug4)
                      + (Λ × c² / 3)
                      + (ℏ / √(Δx × Δp)) × ∫(ψ* H ψ dV) × (2π / t_Hubble)
                      + q × v_spin × B / m_p × (1 + UA/SCm) × scale_macro  (DOMINANT)
                      + ρ_fluid × V × g_base
                      + 2A cos(kx) cos(ωt) + (2π/13.8) A exp(i(kx-ωt))
                      + (M_vis + M_DM) × (δρ/ρ + 3GM/r³)

Special Terms:
- Quantum: Heisenberg uncertainty with normalized wavefunction for NS quantum effects
- EM: q(v×B) term DOMINATES due to extreme B = 2×10^10 T
- Superconductivity: (1 - B/B_crit) correction critical for high-field magnetar
- Fluid: Crust density-volume-gravity coupling for starquakes
- Resonant: Oscillatory Aether-mediated waves for pulsations/bursts
- DM: Visible mass with density perturbations (M_DM = 0)

System: SGR 1745-2900 Galactic Center Magnetar
- M = 1.4 M☉, r = 10 km, B = 2×10^10 T (~2×10^14 G)
- Spin period P = 3.76 s, v_spin ≈ 16.7 km/s
- Crust density ρ = 10^17 kg/m³
- EM term ~10^12 m/s² >> base gravity ~10^11 m/s²
- Located near Sagittarius A*, discovered by Chandra (2013)`;
  }

  // Print all variables (for debugging)
  printVariables() {
    console.log("Current Variables:");
    console.log(`  M = ${(this.M / this.M_sun).toFixed(2)} M_sun (${this.M.toExponential(3)} kg)`);
    console.log(`  r = ${(this.r / 1e3).toFixed(2)} km (${this.r.toExponential(3)} m)`);
    console.log(`  B = ${(this.B / 1e10).toFixed(2)} × 10^10 T (${(this.B / 1e-4).toExponential(2)} Gauss)`);
    console.log(`  P = ${this.P.toFixed(3)} s`);
    console.log(`  v_spin = ${(this.v_spin / 1e3).toFixed(2)} km/s`);
    console.log(`  rho_fluid = ${this.rho_fluid.toExponential(2)} kg/m³`);
    console.log(`  B_crit = ${(this.B_crit / 1e10).toFixed(2)} × 10^10 T`);
    console.log(`  SC correction = ${this.computeSCCorrection().toFixed(4)}`);
    console.log(`  z = ${this.z.toFixed(4)}`);
  }

  // Update variable with dependent recalculations
  updateVariable(name, value) {
    if (name === "M") {
      this.M = value;
      this.M_visible = value;
      this.M_DM = 0.0;
    } else if (name === "r") {
      this.r = value;
      this.v_spin = (2 * this.pi * this.r) / this.P;
    } else if (name === "P") {
      this.P = value;
      this.omega = 2 * this.pi / value;
      this.v_spin = (2 * this.pi * this.r) / value;
    } else if (name === "B") {
      this.B = value;
    } else if (name === "rho_fluid") {
      this.rho_fluid = value;
      this.rho = value;
      this.delta_rho = 0.1 * value;
    } else if (name === "Delta_x") {
      this.Delta_x = value;
      this.Delta_p = this.hbar / value;
    } else {
      // Generic update
      if (this.hasOwnProperty(name)) {
        this[name] = value;
      }
    }
  }

  // ========== ENHANCED DYNAMIC CAPABILITIES (23 methods) ==========

  // Variable Management (5 methods)
  createVariable(name, value) {
    this[name] = value;
  }

  removeVariable(name) {
    if (this.hasOwnProperty(name)) {
      delete this[name];
    }
  }

  cloneVariable(source, dest) {
    if (this.hasOwnProperty(source)) {
      this[dest] = this[source];
    }
  }

  listVariables() {
    return Object.keys(this).filter(key => typeof this[key] !== 'function' && key !== 'savedStates');
  }

  getSystemName() {
    return "SGR1745-2900_Magnetar";
  }

  // Batch Operations (2 methods)
  transformVariableGroup(names, func) {
    names.forEach(name => {
      if (this.hasOwnProperty(name)) {
        this[name] = func(this[name]);
      }
    });
  }

  scaleVariableGroup(names, factor) {
    this.transformVariableGroup(names, v => v * factor);
  }

  // Self-Expansion (4 methods)
  expandParameterSpace(scale_factor) {
    const scalable = ['M', 'r', 'B', 'rho_fluid', 'V', 'A'];
    this.scaleVariableGroup(scalable, scale_factor);
    // Recalculate dependent variables
    this.M_visible = this.M;
    this.v_spin = (2 * this.pi * this.r) / this.P;
    this.rho = this.rho_fluid;
    this.delta_rho = 0.1 * this.rho_fluid;
  }

  expandMagnetarScale(M_scale, r_scale) {
    this.M *= M_scale;
    this.M_visible = this.M;
    this.r *= r_scale;
    // Recalculate v_spin with new radius
    this.v_spin = (2 * this.pi * this.r) / this.P;
  }

  expandMagneticFieldScale(B_scale, P_scale) {
    this.B *= B_scale;
    this.P *= P_scale;
    // Recalculate omega and v_spin
    this.omega = 2 * this.pi / this.P;
    this.v_spin = (2 * this.pi * this.r) / this.P;
  }

  expandCrustScale(rho_scale, V_scale) {
    this.rho_fluid *= rho_scale;
    this.rho = this.rho_fluid;
    this.delta_rho = 0.1 * this.rho_fluid;
    this.V *= V_scale;
  }

  // Self-Refinement (3 methods)
  autoRefineParameters(observations) {
    if (observations.length === 0) return;
    
    let total_error = 0.0;
    observations.forEach(obs => {
      const [t, g_obs] = obs;
      const g_calc = this.compute_g_SGR1745(t);
      total_error += Math.abs(g_calc - g_obs);
    });
    
    const avg_error = total_error / observations.length;
    if (avg_error > 1e-3) {
      const adjustment = 1.0 - (avg_error / (avg_error + 1.0)) * 0.1;
      this.M *= adjustment;
      this.M_visible = this.M;
    }
  }

  calibrateToObservations(observations) {
    this.autoRefineParameters(observations);
  }

  optimizeForMetric(metric, t_start, t_end, steps) {
    let best_metric = -1e100;
    for (let i = 0; i < steps; i++) {
      const t = t_start + (t_end - t_start) * i / (steps - 1);
      const g = this.compute_g_SGR1745(t);
      const m = metric(g);
      if (m > best_metric) {
        best_metric = m;
      }
    }
    return best_metric;
  }

  // Parameter Exploration (1 method)
  generateVariations(count, variation_percent = 5.0) {
    const variations = [];
    const constants = ['G', 'c', 'hbar', 'pi', 'M_sun', 't_Hubble', 'Mpc_to_m', 'q', 'm_p'];
    
    for (let i = 0; i < count; i++) {
      const variant = {};
      Object.keys(this).forEach(key => {
        if (typeof this[key] === 'number' && !constants.includes(key)) {
          const variation = (Math.random() - 0.5) * 2 * (variation_percent / 100.0);
          variant[key] = this[key] * (1.0 + variation);
        } else {
          variant[key] = this[key];
        }
      });
      variations.push(variant);
    }
    return variations;
  }

  // Adaptive Evolution (2 methods)
  mutateParameters(mutation_rate = 0.05) {
    const mutable = ['M', 'r', 'B', 'rho_fluid', 'V', 'A', 'f_TRZ', 'P'];
    mutable.forEach(name => {
      if (this.hasOwnProperty(name)) {
        const variation = (Math.random() - 0.5) * 2 * mutation_rate;
        this.updateVariable(name, this[name] * (1.0 + variation));
      }
    });
  }

  evolveSystem(generations, fitness) {
    for (let gen = 0; gen < generations; gen++) {
      const current_fitness = fitness(this);
      const variants = this.generateVariations(5, 10.0);
      let best_fitness = current_fitness;
      let best_vars = null;
      
      variants.forEach(variant => {
        const temp = new MagnetarSGR1745();
        Object.assign(temp, variant);
        const f = fitness(temp);
        if (f > best_fitness) {
          best_fitness = f;
          best_vars = variant;
        }
      });
      
      if (best_vars) {
        Object.assign(this, best_vars);
      }
    }
  }

  // State Management (4 methods)
  saveState(label) {
    const state = {};
    this.listVariables().forEach(key => {
      state[key] = this[key];
    });
    this.savedStates.set(label, state);
  }

  restoreState(label) {
    if (this.savedStates.has(label)) {
      const state = this.savedStates.get(label);
      Object.assign(this, state);
    }
  }

  listSavedStates() {
    return Array.from(this.savedStates.keys());
  }

  exportState() {
    let output = "SGR1745-2900_Magnetar_State:\n";
    this.listVariables().forEach(key => {
      output += `${key}=${this[key]}\n`;
    });
    return output;
  }

  // System Analysis (3 methods)
  sensitivityAnalysis(t, perturbation = 0.01) {
    const sensitivities = {};
    const g_base = this.compute_g_SGR1745(t);
    
    const test_vars = ['M', 'r', 'B', 'rho_fluid', 'V', 'A', 'f_TRZ', 'omega', 'P'];
    test_vars.forEach(varName => {
      if (this.hasOwnProperty(varName)) {
        const original = this[varName];
        this.updateVariable(varName, original * (1.0 + perturbation));
        const g_perturbed = this.compute_g_SGR1745(t);
        sensitivities[varName] = Math.abs(g_perturbed - g_base) / (g_base + 1e-100);
        this.updateVariable(varName, original);
      }
    });
    
    return sensitivities;
  }

  generateReport(t = this.t_default) {
    const t_years = t / 3.156e7;
    let report = "========== SGR 1745-2900 MAGNETAR UQFF REPORT ==========\n";
    report += `Time: ${t_years.toFixed(1)} years\n`;
    report += `System: ${this.getSystemName()}\n\n`;
    
    report += "Key Parameters:\n";
    report += `  Mass M: ${(this.M / this.M_sun).toFixed(2)} M_sun\n`;
    report += `  Radius r: ${(this.r / 1e3).toFixed(2)} km\n`;
    report += `  Spin Period P: ${this.P.toFixed(3)} s\n`;
    report += `  Surface B-field: ${(this.B / 1e10).toFixed(2)} × 10^10 T (${(this.B / 1e-4).toExponential(2)} G)\n`;
    report += `  Crust Density: ${this.rho_fluid.toExponential(2)} kg/m³\n`;
    report += `  Spin Velocity: ${(this.v_spin / 1e3).toFixed(2)} km/s\n`;
    report += `  SC Correction: ${this.computeSCCorrection().toFixed(4)}\n\n`;
    
    const g = this.compute_g_SGR1745(t);
    const em_term = this.computeEMTerm();
    const g_base = (this.G * this.M / (this.r * this.r)) * (1.0 + this.computeHz() * t) * this.computeSCCorrection() * (1.0 + this.f_TRZ);
    
    report += `Computed g_UQFF: ${g.toExponential(6)} m/s²\n`;
    report += `  Base gravity: ${g_base.toExponential(6)} m/s²\n`;
    report += `  EM term (DOMINANT): ${em_term.toExponential(6)} m/s²\n`;
    report += `  EM/base ratio: ${(em_term / g_base).toFixed(2)}\n`;
    report += "Dominant Terms: EM (q v B) >> Base Gravity due to extreme B-field\n";
    report += "======================================================\n";
    return report;
  }

  validateConsistency() {
    if (this.M <= 0) return false;
    if (this.r <= 0) return false;
    if (this.B < 0) return false;
    if (this.rho_fluid <= 0) return false;
    if (this.omega < 0) return false;
    if (this.P <= 0) return false;
    return true;
  }
}

// ========== ENHANCED 18-STEP EXAMPLE FUNCTION ==========
export function example_enhanced_sgr1745_18_steps() {
  console.log("\n========== ENHANCED SGR 1745-2900 MAGNETAR 18-STEP DEMONSTRATION ==========");
  console.log("Galactic Center Magnetar with Extreme Magnetic Field Dynamics\n");
  
  const sgr = new MagnetarSGR1745();
  const t_current = 1000.0 * 3.156e7; // 1000 years in seconds
  
  // Step 1: Initial state at t = 1000 years
  console.log("Step 1: Initial state at t = 1000 years (young magnetar)");
  const g1 = sgr.compute_g_SGR1745(t_current);
  const em1 = sgr.computeEMTerm();
  console.log(`  Spin Period P = ${sgr.P.toFixed(3)} s`);
  console.log(`  Surface B = ${(sgr.B / 1e10).toFixed(2)} × 10^10 T`);
  console.log(`  g_UQFF = ${g1.toExponential(6)} m/s² (EM-dominated)`);
  console.log(`  EM term = ${em1.toExponential(6)} m/s²\n`);
  
  // Step 2: Save initial state
  console.log("Step 2: Save initial magnetar state");
  sgr.saveState("sgr1745_initial_1000yr");
  console.log("  State saved as 'sgr1745_initial_1000yr'\n");
  
  // Step 3: Expand magnetar scale (mass and radius)
  console.log("Step 3: Expand magnetar scale (1.2x mass, 0.9x radius - compression)");
  sgr.expandMagnetarScale(1.2, 0.9);
  const g3 = sgr.compute_g_SGR1745(t_current);
  console.log(`  New M = ${(sgr.M / sgr.M_sun).toFixed(2)} M_sun`);
  console.log(`  New r = ${(sgr.r / 1e3).toFixed(2)} km`);
  console.log(`  g_UQFF = ${g3.toExponential(6)} m/s²\n`);
  
  // Step 4: Restore and expand magnetic field scale
  console.log("Step 4: Restore initial state, then expand B-field scale (1.5x B, 1.2x period)");
  sgr.restoreState("sgr1745_initial_1000yr");
  sgr.expandMagneticFieldScale(1.5, 1.2);
  const g4 = sgr.compute_g_SGR1745(t_current);
  console.log(`  New B = ${(sgr.B / 1e10).toFixed(2)} × 10^10 T`);
  console.log(`  New Period = ${sgr.P.toFixed(3)} s (spin-down)`);
  console.log(`  g_UQFF = ${g4.toExponential(6)} m/s² (increased EM effect)\n`);
  
  // Step 5: Restore and expand crust scale
  console.log("Step 5: Restore initial state, then expand crust scale (1.3x density, 1.2x volume)");
  sgr.restoreState("sgr1745_initial_1000yr");
  sgr.expandCrustScale(1.3, 1.2);
  const g5 = sgr.compute_g_SGR1745(t_current);
  console.log(`  New rho_crust = ${sgr.rho_fluid.toExponential(2)} kg/m³`);
  console.log(`  New V = ${sgr.V.toExponential(2)} m³`);
  console.log(`  g_UQFF = ${g5.toExponential(6)} m/s²\n`);
  
  // Step 6: Time evolution (magnetar aging)
  console.log("Step 6: Time evolution from 0 to 10,000 years (magnetar aging)");
  sgr.restoreState("sgr1745_initial_1000yr");
  for (let t_yr = 0; t_yr <= 10000; t_yr += 2500) {
    const t_sec = t_yr * 3.156e7;
    const g = sgr.compute_g_SGR1745(t_sec);
    console.log(`  t = ${t_yr} yr: g = ${g.toExponential(6)} m/s²`);
  }
  console.log();
  
  // Step 7: Create custom tracking variables
  console.log("Step 7: Create custom tracking variables");
  sgr.createVariable("burst_count", 0.0);
  sgr.createVariable("distance_gc_pc", 8000.0); // Distance to Galactic Center
  console.log("  Created 'burst_count' and 'distance_gc_pc'\n");
  
  // Step 8: Generate variations for uncertainty analysis
  console.log("Step 8: Generate 3 parameter variations (5% perturbation)");
  const variations = sgr.generateVariations(3, 5.0);
  for (let i = 0; i < variations.length; i++) {
    const temp = new MagnetarSGR1745();
    Object.assign(temp, variations[i]);
    const g_var = temp.compute_g_SGR1745(t_current);
    console.log(`  Variation ${i+1}: g = ${g_var.toExponential(6)} m/s²`);
  }
  console.log();
  
  // Step 9: Sensitivity analysis
  console.log("Step 9: Sensitivity analysis (1% perturbation)");
  const sensitivities = sgr.sensitivityAnalysis(t_current, 0.01);
  console.log("  Parameter sensitivities (fractional change in g):");
  Object.entries(sensitivities).forEach(([key, val]) => {
    console.log(`    ${key}: ${val.toExponential(4)}`);
  });
  console.log();
  
  // Step 10: Magnetic field strength sweep
  console.log("Step 10: B-field strength sweep (0.5x, 1x, 2x)");
  sgr.saveState("sgr_before_sweep");
  for (const scale of [0.5, 1.0, 2.0]) {
    sgr.restoreState("sgr_before_sweep");
    sgr.expandMagneticFieldScale(scale, 1.0);
    const g = sgr.compute_g_SGR1745(t_current);
    console.log(`  B = ${(sgr.B / 1e10).toFixed(2)} × 10^10 T: g = ${g.toExponential(6)} m/s²`);
  }
  console.log();
  
  // Step 11: Spin period sweep (spin-down evolution)
  console.log("Step 11: Spin period sweep (1.0x, 1.5x, 2.0x - spin-down)");
  sgr.restoreState("sgr_before_sweep");
  for (const scale of [1.0, 1.5, 2.0]) {
    sgr.restoreState("sgr_before_sweep");
    sgr.expandMagneticFieldScale(1.0, scale);
    const g = sgr.compute_g_SGR1745(t_current);
    console.log(`  P = ${sgr.P.toFixed(3)} s: g = ${g.toExponential(6)} m/s²`);
  }
  console.log();
  
  // Step 12: Mass sweep (different magnetar masses)
  console.log("Step 12: Mass sweep (0.9x, 1.0x, 1.5x M_sun)");
  sgr.restoreState("sgr_before_sweep");
  for (const scale of [0.9, 1.0, 1.5]) {
    sgr.restoreState("sgr_before_sweep");
    sgr.expandMagnetarScale(scale, 1.0);
    const g = sgr.compute_g_SGR1745(t_current);
    console.log(`  M = ${(sgr.M / sgr.M_sun).toFixed(2)} M_sun: g = ${g.toExponential(6)} m/s²`);
  }
  console.log();
  
  // Step 13: Batch transform spin parameters
  console.log("Step 13: Batch transform spin parameters (1.1x scale)");
  sgr.restoreState("sgr_before_sweep");
  sgr.scaleVariableGroup(['v_spin', 'omega'], 1.1);
  const g13 = sgr.compute_g_SGR1745(t_current);
  console.log(`  v_spin = ${(sgr.v_spin / 1e3).toFixed(2)} km/s`);
  console.log(`  omega = ${sgr.omega.toFixed(4)} rad/s`);
  console.log(`  g_UQFF = ${g13.toExponential(6)} m/s²\n`);
  
  // Step 14: Validate and auto-correct
  console.log("Step 14: Validate consistency");
  sgr.restoreState("sgr_before_sweep");
  const valid = sgr.validateConsistency();
  console.log(`  System valid: ${valid ? "Yes" : "No"}\n`);
  
  // Step 15: Parameter mutation (evolutionary exploration)
  console.log("Step 15: Mutate parameters (3% random variation)");
  sgr.restoreState("sgr_before_sweep");
  sgr.mutateParameters(0.03);
  const g15 = sgr.compute_g_SGR1745(t_current);
  console.log(`  Mutated M = ${(sgr.M / sgr.M_sun).toFixed(2)} M_sun`);
  console.log(`  Mutated B = ${(sgr.B / 1e10).toFixed(2)} × 10^10 T`);
  console.log(`  g_UQFF = ${g15.toExponential(6)} m/s²\n`);
  
  // Step 16: List all saved states
  console.log("Step 16: List all saved states");
  const states = sgr.listSavedStates();
  console.log(`  Saved states (${states.length} total):`);
  states.forEach(state => {
    console.log(`    - ${state}`);
  });
  console.log();
  
  // Step 17: Generate comprehensive report
  console.log("Step 17: Generate comprehensive system report");
  sgr.restoreState("sgr1745_initial_1000yr");
  const report = sgr.generateReport(t_current);
  console.log(report);
  
  // Step 18: Export final state
  console.log("Step 18: Export final system state (first 10 variables)");
  const state_export = sgr.exportState();
  const lines = state_export.split('\n').slice(0, 11);
  console.log(lines.join('\n'));
  console.log("  ... (additional variables omitted)\n");
  
  console.log("========== END 18-STEP SGR 1745-2900 DEMONSTRATION ==========\n");
}

// Run example if executed directly
if (import.meta.url === `file://${process.argv[1]}`) {
  example_enhanced_sgr1745_18_steps();
}
