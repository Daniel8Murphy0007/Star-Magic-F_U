/**
 * UQFF Compressed Resonance Module (JavaScript Port)
 * 
 * Modular implementation of compressed and resonance UQFF equations for multi-system evolution.
 * Supports both compressed gravity mode and resonance oscillatory mode.
 * 
 * Supports multiple astrophysical systems:
 * - YoungStars: Pre-main sequence star outflows
 * - Eagle: Eagle Nebula star formation
 * - BigBang: Cosmic microwave background era (z~1100)
 * - M51: Whirlpool Galaxy
 * - NGC1316: Centaurus A (merger remnant)
 * - V838Mon: Luminous red nova light echo
 * - NGC1300: Barred spiral galaxy
 * - Guide: General educational system
 * 
 * Original C++ Implementation: Source74.cpp
 * JavaScript Port: November 1, 2025
 * 
 * Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.
 */

class UQFFCompressedResonanceModule {
  /**
   * Constructor: Initialize with universal constants and default parameters
   */
  constructor() {
    // Universal constants
    this.G = 6.6743e-11;                      // m³ kg⁻¹ s⁻²
    this.c = 3e8;                             // m/s
    this.hbar = 1.0546e-34;                   // J·s
    this.Lambda = 1.1e-52;                    // m⁻²
    this.q = 1.602e-19;                       // C
    this.pi = Math.PI;                        // π

    // Cosmological parameters
    this.t_Hubble = 13.8e9 * 3.156e7;         // s
    this.year_to_s = 3.156e7;                 // s/yr
    this.H0 = 70.0;                           // km/s/Mpc
    this.Mpc_to_m = 3.086e22;                 // m/Mpc
    this.Omega_m = 0.3;                       // Matter density
    this.Omega_Lambda = 0.7;                  // Dark energy density

    // Solar & distance units
    this.M_sun = 1.989e30;                    // kg
    this.kpc = 3.086e19;                      // m

    // General defaults (overridden by setSystem)
    this.M = 1e41;                            // kg
    this.M0 = this.M;
    this.SFR = 6e19;                          // kg/s (~2 M☉/yr)
    this.r = 1e20;                            // m
    this.z = 0.005;                           // Redshift
    this.M_visible = 0.7 * this.M;
    this.M_DM = 0.3 * this.M;
    this.t = 1e9 * this.year_to_s;            // s
    this.rho_fluid = 1e-21;                   // kg/m³
    this.V = 1e50;                            // m³
    this.B = 1e-5;                            // T
    this.B_crit = 1e11;                       // T
    this.Delta_x = 1e-10;                     // m
    this.Delta_p = this.hbar / this.Delta_x;  // kg·m/s
    this.integral_psi = 1.0;                  // Normalized

    // Resonance & wave parameters
    this.A = 1e-10;                           // Amplitude
    this.k = 1e20;                            // m⁻¹ (wavenumber)
    this.omega = 1e15;                        // rad/s (frequency)
    this.x = 0.0;                             // Position
    this.v = 1e3;                             // m/s (velocity)

    // Ug components
    this.Ug1 = 0.0;
    this.Ug2 = 0.0;
    this.Ug3 = 0.0;
    this.Ug4 = 0.0;

    // Scales and corrections
    this.scale_macro = 1e-12;
    this.f_TRZ = 0.1;
    this.f_sc = 1.0;
    this.delta_rho = 1e-5 * this.rho_fluid;
    this.rho = this.rho_fluid;
    this.F_wind = 0.0;

    // State
    this.current_system = "Guide";
    this.mode = "compressed";  // "compressed" or "resonance"
  }

  /**
   * Set system and load system-specific parameters (DeepSearch auto-load simulation)
   * 
   * @param {string} sys_name - System identifier (YoungStars, Eagle, BigBang, M51, NGC1316, V838Mon, NGC1300, Guide)
   */
  setSystem(sys_name) {
    this.current_system = sys_name;

    const yr_s = this.year_to_s;
    const kpc = this.kpc;
    const Msun = this.M_sun;

    if (sys_name === "YoungStars") {
      // Pre-main sequence outflows
      this.M = 1000 * Msun;
      this.r = 3e17;
      this.SFR = 0.1 * Msun / yr_s;
      this.rho_fluid = 1e-20;
      this.B = 1e-6;
      this.z = 0.0006;
    } else if (sys_name === "Eagle") {
      // Eagle Nebula star formation
      this.M = 1e4 * Msun;
      this.r = 2e17;
      this.SFR = 0.5 * Msun / yr_s;
      this.rho_fluid = 1e-21;
      this.B = 3e-5;
      this.z = 0.002;
    } else if (sys_name === "BigBang") {
      // Cosmic microwave background era (z~1100)
      this.rho_fluid = 8e-27;
      this.r = 1e26;
      this.z = 1100;
      this.SFR = 0;
      this.M = 1e53;  // Observable universe ~10^53 kg
      this.B = 1e-10;
      this.t = 13.8e9 * yr_s;
    } else if (sys_name === "M51") {
      // Whirlpool Galaxy
      this.M = 1.6e11 * Msun;
      this.r = 23e3 * kpc;
      this.SFR = 2 * Msun / yr_s;
      this.rho_fluid = 1e-21;
      this.B = 1e-5;
      this.z = 0.005;
    } else if (sys_name === "NGC1316") {
      // Centaurus A (merger remnant)
      this.M = 5e11 * Msun;
      this.r = 23e3 * kpc;
      this.SFR = 0.1 * Msun / yr_s;
      this.rho_fluid = 1e-22;
      this.B = 1e-5;
      this.z = 0.006;
    } else if (sys_name === "V838Mon") {
      // Luminous red nova
      this.M = 8 * Msun;
      this.r = 2e13;
      this.SFR = 0;
      this.rho_fluid = 1e-22;
      this.B = 1e-6;
      this.z = 0.005;
    } else if (sys_name === "NGC1300") {
      // Barred spiral galaxy
      this.M = 1e11 * Msun;
      this.r = 12e3 * kpc;
      this.SFR = 1 * Msun / yr_s;
      this.rho_fluid = 1e-21;
      this.B = 1e-5;
      this.z = 0.005;
    } else {
      // Default: Guide/educational
      this.M = Msun;
      this.r = 1e11;
      this.SFR = 1e-10 * Msun / yr_s;
      this.rho_fluid = 1e-20;
      this.B = 1e-5;
      this.z = 0;
    }

    // Update dependents
    this.M_visible = 0.7 * this.M;
    this.M_DM = 0.3 * this.M;
    this.M0 = this.M;
    this.Delta_p = this.hbar / this.Delta_x;
  }

  /**
   * Set computation mode: "compressed" or "resonance"
   * 
   * @param {string} m - Mode identifier
   */
  setMode(m) {
    if (m === "compressed" || m === "resonance") {
      this.mode = m;
    } else {
      console.warn(`Unknown mode "${m}". Keeping current mode "${this.mode}".`);
    }
  }

  /**
   * Update a variable
   * 
   * @param {string} name - Variable name
   * @param {number} value - New value
   */
  updateVariable(name, value) {
    this[name] = value;
    
    // Cascade updates for dependent parameters
    if (name === "M") {
      this.M_visible = 0.7 * value;
      this.M_DM = 0.3 * value;
      this.M0 = value;
    } else if (name === "Delta_x") {
      this.Delta_p = this.hbar / value;
    }
  }

  /**
   * Add delta to a variable
   * 
   * @param {string} name - Variable name
   * @param {number} delta - Value to add
   */
  addToVariable(name, delta) {
    this.updateVariable(name, this[name] + delta);
  }

  /**
   * Subtract delta from a variable
   * 
   * @param {string} name - Variable name
   * @param {number} delta - Value to subtract
   */
  subtractFromVariable(name, delta) {
    this.updateVariable(name, this[name] - delta);
  }

  /**
   * Compute Hubble parameter H(z)
   * 
   * @param {number} z_val - Redshift
   * @returns {number} H(z) in s⁻¹
   */
  computeHtz(z_val) {
    const Hz_kms = this.H0 * Math.sqrt(
      this.Omega_m * Math.pow(1 + z_val, 3) + this.Omega_Lambda
    );
    return (Hz_kms * 1e3) / this.Mpc_to_m;
  }

  /**
   * Compute environmental forcing (simplified)
   * Placeholder for F_bar + F_SF + F_wave, or can be overridden
   * 
   * @param {number} t - Time in seconds
   * @returns {number} Environmental force factor
   */
  computeFenv(t) {
    // Simplified: constant factor
    // Full implementations would include bar, star formation, and wave terms
    return 0.1;
  }

  /**
   * Compute sum of Ug gravity components (placeholder)
   * Full implementations would include Ug1-Ug4 with physical models
   * 
   * @returns {number} Sum of Ug components
   */
  computeUgSum() {
    // Placeholder: return constant
    // Real implementation would compute Ug1 (dipole), Ug2 (superconductor), 
    // Ug3 (external), Ug4 (reaction) and sum them
    return 1e-10;
  }

  /**
   * Compute total wave perturbation psi
   * Combines charged particle motion and density wave terms
   * 
   * @param {number} t - Time in seconds
   * @returns {number} Total perturbation
   */
  computePsiTotal(t) {
    // Lorentz term: q·v·B
    const lorentz = this.q * this.v * this.B;
    
    // Density wave: 2·A·cos(k·x + ω·t)
    const wave = 2 * this.A * Math.cos(this.k * this.x + this.omega * t);
    
    return lorentz + wave;
  }

  /**
   * Compute resonance term (SAFE: does not call computeG to avoid recursion)
   * Oscillatory coupling with pre-computed g_base
   * 
   * @param {number} t - Time in seconds
   * @param {number} g_base - Pre-computed base gravitational acceleration
   * @returns {number} Resonance contribution to gravity
   */
  computeResonanceTerm(t, g_base) {
    if (this.mode !== "resonance") {
      return 0.0;
    }

    // Resonant amplitude: A·exp(i(k·x - ω·t))
    // Real part: A·cos(k·x - ω·t)
    const real_part = this.A * Math.cos(this.k * this.x - this.omega * t);
    
    // Resonance coupling: (2π / t_Hubble_Gyr) × real_part × g_base
    // t_Hubble_Gyr = 13.8 Gyr
    const coupling = (2 * this.pi / 13.8) * real_part * g_base;
    
    return coupling;
  }

  /**
   * Compute quantum gravity term (Heisenberg uncertainty + wave evolution)
   * 
   * @param {number} t_Hubble_val - Hubble time in seconds
   * @returns {number} Quantum correction
   */
  computeQuantumTerm(t_Hubble_val) {
    const unc = Math.sqrt(this.Delta_x * this.Delta_p);
    const psi = this.computePsiTotal(this.t);
    
    return (this.hbar / unc) * this.integral_psi * (2 * this.pi / t_Hubble_val) * psi;
  }

  /**
   * Compute fluid dynamics contribution
   * 
   * @param {number} g_base - Base gravitational acceleration
   * @returns {number} Fluid term
   */
  computeFluidTerm(g_base) {
    return this.rho_fluid * this.V * g_base;
  }

  /**
   * Compute dark matter response term
   * Includes density perturbations and gravitational curvature
   * 
   * @returns {number} Dark matter contribution
   */
  computeDMTerm() {
    const pert = this.delta_rho / this.rho;
    const curv = 3 * this.G * this.M / Math.pow(this.r, 3);
    
    return (this.M_visible + this.M_DM) * (pert + curv);
  }

  /**
   * Compute star formation mass growth factor
   * 
   * @param {number} t - Time in seconds
   * @returns {number} Fractional mass increase (SFR·t / M0)
   */
  computeMsfFactor(t) {
    return this.SFR * t / this.M0;
  }

  /**
   * Compute full gravitational/intensity field g_UQFF(r, t)
   * Core function: integrates all physics components with optional recursion-safe resonance
   * 
   * Special case: V838Mon system returns I_echo (light echo intensity) instead
   * Special case: BigBang updates r = c·t (cosmological expansion)
   * 
   * @param {number} t - Time in seconds
   * @param {number} r_in - Optional radius (if 0, uses default this.r)
   * @returns {number} g_UQFF(r, t) [m/s²] or I_echo [W/m²] for V838Mon
   */
  computeG(t, r_in = 0.0) {
    // Update radius if provided
    if (r_in > 0) {
      this.r = r_in;
    }
    this.t = t;

    // Special case: BigBang uses r = c·t (light travel distance)
    if (this.current_system === "BigBang") {
      this.r = this.c * t;
    }

    // Special case: V838Mon returns I_echo (light echo intensity)
    if (this.current_system === "V838Mon") {
      // Dust density modulated by gravity
      const rho_d = this.rho_fluid * Math.exp(-1.0 * (this.G * this.M / (this.r * this.r)));
      
      // I_echo ≈ (L_outburst / 4π·r²) · σ_scatter · ρ_d
      // Approximate: 600,000 L☉ outburst
      const L_outburst = 600000 * 3.826e26;  // W
      const sigma_scatter = 1e-12;            // m²
      
      return (L_outburst / (4 * this.pi * this.r * this.r)) * sigma_scatter * rho_d;
    }

    // Standard UQFF computation:

    // Cosmological expansion factor
    const Hz = this.computeHtz(this.z);
    const expansion = 1.0 + Hz * t;

    // Magnetic suppression factor
    const sc = 1.0 - this.B / this.B_crit;

    // Star formation mass growth
    const msf = this.computeMsfFactor(t);
    const m_fact = 1.0 + msf;

    // Environmental forcing
    const f_env = this.computeFenv(t);

    // Base gravity: (G·M(t)/r²) × expansion × mag_suppression × env_forcing
    const g_base = (this.G * this.M * m_fact / (this.r * this.r)) * 
                   expansion * sc * (1.0 + f_env);

    // Sum of Ug components
    const ug_sum = this.computeUgSum();

    // Cosmological constant term
    const lambda_t = this.Lambda * this.c * this.c / 3.0;

    // Quantum gravity term
    const q_term = this.computeQuantumTerm(this.t_Hubble);

    // Fluid dynamics term
    const f_term = this.computeFluidTerm(g_base);

    // Dark matter term
    const dm_term = this.computeDMTerm();

    // RECURSION FIX: Compute resonance term with g_base (does not call computeG)
    const res_term = this.computeResonanceTerm(t, g_base);

    // Total: all components summed
    return g_base + ug_sum + lambda_t + q_term + f_term + dm_term + res_term;
  }

  /**
   * Compute g_UQFF at standard system radius
   * 
   * @param {number} t - Time in seconds
   * @returns {number} g_UQFF at standard radius
   */
  computeGStandard(t) {
    return this.computeG(t, this.r);
  }

  /**
   * Get complete equation description with mode-specific details
   * 
   * @returns {string} Formatted equation text
   */
  getEquationText() {
    let eq = `g_UQFF(r,t) = [G·M(t)/r²] · (1 + H(t,z)) · (1 - B/B_crit) · (1 + F_env)
                + ΣUg_i + Λ·c²/3
                + (ℏ/√(Δx·Δp)) · ∫|ψ|² dV · (2π/t_Hubble)
                + ρ_fluid · V · g + (M_visible + M_DM) · (Δρ/ρ + 3GM/r³)`;

    if (this.mode === "resonance") {
      eq += `
                + [resonance mode: A·cos(k·x - ω·t) + (2π/13.8)·Re[A·exp(i(k·x - ω·t))]·g_base]`;
    }

    eq += `

Where:
  M(t) = M₀(1 + SFR·t/M₀)  [Star formation mass growth]
  H(t,z) = H₀√(Ω_m(1+z)³ + Ω_Λ)  [Hubble parameter]
  F_env ≈ 0.1  [Environmental forcing (bar, SF, waves)]
  
System: ${this.current_system}
Mode: ${this.mode}
Physics: Compressed unified field with optional resonance oscillations
Learning: Diverse scales (young stars to BigBang) refine UQFF accuracy
Advancing: Resonance coupling explains outflow acceleration and cosmic expansion patterns
    `;

    return eq;
  }

  /**
   * Get physics summary for the object
   * 
   * @returns {Object} Summary of key properties
   */
  getSummary() {
    return {
      name: 'UQFF Compressed Resonance Module',
      physicsType: 'Multi-system compressed and resonance gravity',
      current_system: this.current_system,
      mode: this.mode,
      mass_Msun: this.M / this.M_sun,
      radius_m: this.r,
      starFormationRate_Msun_yr: this.SFR / (this.M_sun / this.year_to_s),
      redshift: this.z,
      magnetic_field_T: this.B,
      systemsAvailable: ['YoungStars', 'Eagle', 'BigBang', 'M51', 'NGC1316', 'V838Mon', 'NGC1300', 'Guide'],
      modesAvailable: ['compressed', 'resonance'],
      uqffComponentsActive: ['g_base', 'Ug_sum', 'Lambda', 'quantum', 'fluid', 'dark_matter', 'resonance_if_enabled'],
      equation: 'g_UQFF = [G·M(t)/r²]·expansion·mag_suppression·env_forcing + all_additional_terms',
      typicalOutput_ms2: this.computeGStandard(this.t),
      description: 'Versatile multi-system UQFF module with compressed and resonance modes for astrophysical modeling'
    };
  }

  /**
   * Print formatted variable summary
   */
  printVariables() {
    console.log(`\n=== UQFF Compressed Resonance Module ===`);
    console.log(`System: ${this.current_system}, Mode: ${this.mode}\n`);
    console.log(`Constants:
  G=${this.G.toExponential(6)}, c=${this.c.toExponential(6)}, ℏ=${this.hbar.toExponential(6)}
  t_Hubble=${(this.t_Hubble / this.year_to_s / 1e9).toFixed(1)} Gyr
  
Cosmology:
  H₀=${this.H0} km/s/Mpc, Ω_m=${this.Omega_m}, Ω_Λ=${this.Omega_Lambda}, z=${this.z}
  
System Parameters:
  M=${(this.M / this.M_sun).toExponential(6)} M☉
  r=${this.r.toExponential(6)} m
  SFR=${(this.SFR / (this.M_sun / this.year_to_s)).toExponential(6)} M☉/yr
  
Physics:
  B=${this.B.toExponential(6)} T, B_crit=${this.B_crit.toExponential(6)} T
  ρ_fluid=${this.rho_fluid.toExponential(6)} kg/m³
  f_TRZ=${this.f_TRZ}
  
Resonance (if enabled):
  A=${this.A.toExponential(6)}, k=${this.k.toExponential(6)} m⁻¹, ω=${this.omega.toExponential(6)} rad/s
`);
  }
}

// Module export for integration
module.exports = UQFFCompressedResonanceModule;
