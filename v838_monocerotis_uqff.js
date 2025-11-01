/**
 * V838 Monocerotis UQFF Module (JavaScript Port)
 * 
 * Modular implementation of the Master Universal Gravity Equation (MUGE & UQFF Integration)
 * for V838 Monocerotis Light Echo Evolution.
 * 
 * This module models light echo intensity evolution incorporating:
 * - Outburst luminosity (600,000 L☉)
 * - Dust scattering with gravitational modulation (Ug1)
 * - Time-reversal effects (TRZ factor)
 * - Aetheric vacuum energy corrections (UA/SCm)
 * 
 * Physics Domain: Luminous Red Nova with dramatic light echo phenomenon
 * Astrophysical Context: V838 Mon located at ~6.1 kpc in constellation Monoceros
 * 
 * Original C++ Implementation: Source72.cpp
 * JavaScript Port: November 1, 2025
 * 
 * Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.
 */

class V838MonocerotisUQFFModule {
  /**
   * Constructor: Initialize with V838 Monocerotis defaults
   */
  constructor() {
    // Universal constants
    this.c = 3e8;                           // m/s (speed of light)
    this.G = 6.6743e-11;                    // m³ kg⁻¹ s⁻² (gravitational constant)
    this.hbar = 1.0546e-34;                 // J s (reduced Planck constant)
    this.pi = Math.PI;                      // π

    // Solar reference
    const M_sun_val = 1.989e30;             // kg
    const L_sun_val = 3.826e26;             // W

    // V838 Mon astrophysical parameters
    this.M_s = 8 * M_sun_val;               // kg (stellar mass: 8 M☉)
    this.L_outburst = 600000 * L_sun_val;   // W (~2.3e38 W, ~600,000 L☉)
    this.distance = 6.1 * 3.086e22;         // m (6.1 kpc converted to meters)
    
    // Dust parameters
    this.rho_0 = 1e-22;                     // kg/m³ (dust density at reference)
    this.sigma_scatter = 1e-12;             // m² (dust grain scattering cross-section)
    
    // UQFF coupling parameters
    this.k1 = 1.0;                          // Ug1 scaling factor
    this.mu_s = 1.0;                        // Superconductive permeability
    this.alpha = 0.0005;                    // Exponential decay rate
    this.beta = 1.0;                        // Dust modulation strength
    
    // Temporal phase
    this.t_n = 0.0;                         // Phase component
    
    // Time-reversal zone factor
    this.f_TRZ = 0.1;                       // Time-reversal correction
    
    // Vacuum energy densities (UQFF aetheric components)
    this.rho_vac_UA = 7.09e-36;             // J/m³ (Universal Aether)
    this.rho_vac_SCm = 7.09e-37;            // J/m³ (Superconductive material)
    
    // Default time
    this.t = 3 * 3.156e7;                   // Default t = 3 years in seconds
    
    // Periodic deflection
    this.delta_def = 0.01 * Math.sin(0.001 * 1e7);  // Periodic component at t=0
    
    // Scaling factor
    this.scale_macro = 1e-12;               // Macroscopic scale correction
  }

  /**
   * Update a variable dynamically
   * @param {string} name - Variable name
   * @param {number} value - New value
   */
  updateVariable(name, value) {
    if (this.hasOwnProperty(name)) {
      this[name] = value;
      // Update delta_def when t changes
      if (name === 't') {
        this.delta_def = 0.01 * Math.sin(0.001 * value);
      }
    } else {
      console.warn(`Variable '${name}' not found. Adding as new property.`);
      this[name] = value;
    }
  }

  /**
   * Add delta to a variable
   * @param {string} name - Variable name
   * @param {number} delta - Value to add
   */
  addToVariable(name, delta) {
    if (this.hasOwnProperty(name)) {
      this[name] += delta;
    } else {
      this[name] = delta;
    }
  }

  /**
   * Subtract delta from a variable
   * @param {string} name - Variable name
   * @param {number} delta - Value to subtract
   */
  subtractFromVariable(name, delta) {
    this.addToVariable(name, -delta);
  }

  /**
   * Compute Ug1: Universal Gravity component 1
   * Represents gravitational gradient with exponential decay and phase modulation
   * 
   * @param {number} t - Time in seconds
   * @param {number} r - Radial distance in meters
   * @returns {number} Ug1 value (N/m²)
   */
  computeUg1(t, r) {
    // Simplified gravitational gradient: ∇(M_s/r) ≈ M_s/r³
    const grad_term = this.M_s / (r * r * r);
    
    // Exponential decay envelope
    const exp_decay = Math.exp(-this.alpha * t);
    
    // Phase cosine modulation
    const cos_phase = Math.cos(this.pi * this.t_n);
    
    // Periodic deflection (updates with t)
    const delta = 0.01 * Math.sin(0.001 * t);
    
    // Full Ug1 computation
    return this.k1 * this.mu_s * grad_term * exp_decay * cos_phase * (1 + delta);
  }

  /**
   * Compute dust density modulated by Ug1
   * 
   * @param {number} r - Radial distance in meters
   * @param {number} t - Time in seconds
   * @returns {number} Dust density ρ_dust (kg/m³)
   */
  computeRhodust(r, t) {
    const ug1 = this.computeUg1(t, r);
    return this.rho_0 * Math.exp(-this.beta * ug1);
  }

  /**
   * Compute base light echo intensity (without modulation)
   * Simple inverse-square law from outburst luminosity
   * 
   * @param {number} r - Radial distance in meters
   * @returns {number} Base intensity I_base (W/m²)
   */
  computeIechoBase(r) {
    return this.L_outburst / (4 * this.pi * r * r);
  }

  /**
   * Compute time-reversal zone correction factor
   * TRZ accounts for temporal symmetry breaking in light echo region
   * 
   * @returns {number} TRZ correction factor (dimensionless)
   */
  computeTRZCorrection() {
    return 1.0 + this.f_TRZ;
  }

  /**
   * Compute Aether/Superconductive Material correction
   * UA/SCm ratio captures aetheric propagation effects
   * 
   * @returns {number} UA/SCm correction factor (dimensionless)
   */
  computeUAscCorrection() {
    return 1.0 + (this.rho_vac_UA / this.rho_vac_SCm);
  }

  /**
   * Compute full light echo intensity I_echo(r, t)
   * Complete UQFF equation with all physical effects:
   * I_echo = I_base × σ_scatter × ρ_dust × TRZ_correction × UA_correction
   * 
   * @param {number} t - Time in seconds since outburst
   * @param {number} r - Radial distance in meters (typically r = c·t for light echo)
   * @returns {number} Light echo intensity (W/m²)
   */
  computeIecho(t, r) {
    // Update temporal phase
    this.t = t;
    this.delta_def = 0.01 * Math.sin(0.001 * t);

    // Compute all physical components
    const rho_dust = this.computeRhodust(r, t);
    const i_base = this.computeIechoBase(r);
    const trz = this.computeTRZCorrection();
    const ua_sc = this.computeUAscCorrection();

    // Full equation with all modulations
    return i_base * this.sigma_scatter * rho_dust * trz * ua_sc;
  }

  /**
   * Compute standard light echo at characteristic radius r = c·t
   * 
   * @param {number} t - Time in seconds
   * @returns {number} Light echo intensity at r_echo = c·t (W/m²)
   */
  computeIechoStandard(t) {
    const r_echo = this.c * t;
    return this.computeIecho(t, r_echo);
  }

  /**
   * Get complete equation description with all components
   * 
   * @returns {string} Formatted equation text with physical interpretation
   */
  getEquationText() {
    return `I_echo(r, t) = [L_outburst / (4π r²)] × σ_scatter × ρ_0 × exp(-β [k₁ μ_s (M_s/r³) e^(-αt) cos(πt_n) (1 + δ_def)]) × (1 + f_TRZ) × (1 + ρ_vac,UA/ρ_vac,SCm)

Where:
  r_echo(t) = c·t (light echo radius at time t)
  δ_def = 0.01 sin(0.001·t) (periodic temporal modulation)
  ∇(M_s/r) ≈ M_s/r³ (simplified gravitational gradient)
  L_outburst ≈ 2.3×10³⁸ W (V838 Mon ~600,000 L☉)
  ρ_0 = 1×10⁻²² kg/m³ (dust density)
  f_TRZ = 0.1 (time-reversal zone factor)
  ρ_vac,UA = 7.09×10⁻³⁶ J/m³ (Universal Aether)
  ρ_vac,SCm = 7.09×10⁻³⁷ J/m³ (Superconductive material)

Physics Interpretation:
  • Attractive force (Ug1) modulates dust density and scattering efficiency
  • Time-reversal symmetry (TRZ) breaks causality locally, enhancing propagation
  • Aetheric correction (UA) increases effective wave propagation beyond standard electromagnetism
  • Dust density decreases exponentially with Ug1 strength (negative feedback)
  • Result: I_echo ~1×10⁻²⁰ W/m² at t=3 yr, r=9×10¹⁵ m (dust scattering dominant)

Astronomical Context:
  V838 Monocerotis: M_s=8 M☉, d=6.1 kpc, outburst observed 2002
  Hubble ACS data validation: 2004 light echo measurements
  Unique phenomenon: Dramatic shell expansion with complex dust interaction
  UQFF insight: Gravity (Ug1) + Aether (UA) explain observed echo dimming patterns`;
  }

  /**
   * Get all current variables and their values
   * 
   * @returns {Object} Object containing all module variables
   */
  getVariables() {
    return {
      c: this.c,
      G: this.G,
      hbar: this.hbar,
      pi: this.pi,
      M_s: this.M_s,
      L_outburst: this.L_outburst,
      distance: this.distance,
      rho_0: this.rho_0,
      sigma_scatter: this.sigma_scatter,
      k1: this.k1,
      mu_s: this.mu_s,
      alpha: this.alpha,
      beta: this.beta,
      t_n: this.t_n,
      f_TRZ: this.f_TRZ,
      rho_vac_UA: this.rho_vac_UA,
      rho_vac_SCm: this.rho_vac_SCm,
      t: this.t,
      delta_def: this.delta_def,
      scale_macro: this.scale_macro
    };
  }

  /**
   * Print formatted variable summary
   */
  printVariables() {
    console.log('\n=== V838 Monocerotis UQFF Module Variables ===\n');
    const vars = this.getVariables();
    Object.entries(vars).forEach(([name, value]) => {
      console.log(`${name.padEnd(20)} = ${value.toExponential(6)}`);
    });
    console.log('\n');
  }

  /**
   * Get physics summary for the object
   * 
   * @returns {Object} Summary of key physical properties
   */
  getSummary() {
    return {
      name: 'V838 Monocerotis UQFF Module',
      physicsType: 'Luminous Red Nova Light Echo Evolution',
      stellarMass_Msun: this.M_s / 1.989e30,
      outburstLuminosity_Lsun: this.L_outburst / 3.826e26,
      distance_kpc: this.distance / 3.086e22,
      dustDensity: this.rho_0,
      uqffComponentsActive: ['Ug1', 'TRZ', 'UA', 'SCm'],
      equation: 'I_echo(r,t) = I_base × σ_scatter × ρ_dust(Ug1) × TRZ × UA/SCm',
      typicalOutput_Wm2: this.computeIechoStandard(this.t),
      description: 'Models light echo intensity evolution with gravitational dust modulation and aetheric effects'
    };
  }
}

// Module export for integration
module.exports = V838MonocerotisUQFFModule;
