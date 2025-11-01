/**
 * NGC 2264 UQFF Module (JavaScript Port)
 * 
 * Modular implementation of the Master Universal Gravity Equation (MUGE & UQFF Integration)
 * for Cone Nebula (NGC 2264) gravitational dynamics and evolution.
 * 
 * This module models NGC 2264's complex physics, incorporating:
 * - Stellar winds from massive O/B stars (20 km/s)
 * - Pillar erosion dynamics (5% per 3 Myr)
 * - Protostar formation and spin (1e-5 rad/s)
 * - Dust and gas densities (1e-20 kg/m³)
 * - Dark matter response
 * - Quantum pillar wave structure (m-mode)
 * - Environmental feedback (wind + star formation + erosion)
 * 
 * Physics Domain: Star-forming nebula with pillar structure
 * Astrophysical Context: NGC 2264 Cone Nebula in Monoceros constellation
 * Reference: Hubble ACS 2002 observations
 * 
 * Original C++ Implementation: Source76.cpp
 * JavaScript Port: November 1, 2025
 * 
 * Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.
 */

class NGC2264UQFFModule {
  /**
   * Constructor: Initialize with NGC 2264 defaults
   */
  constructor() {
    // Universal constants
    this.G = 6.6743e-11;                      // m³ kg⁻¹ s⁻²
    this.c = 3e8;                             // m/s
    this.hbar = 1.0546e-34;                   // J·s
    this.Lambda = 1.1e-52;                    // m⁻²
    this.q = 1.602e-19;                       // C
    this.pi = Math.PI;                        // π
    this.mu_0 = 4 * this.pi * 1e-7;           // H/m

    // Cosmological parameters
    this.t_Hubble = 13.8e9 * 3.156e7;         // s
    this.year_to_s = 3.156e7;                 // s/yr
    this.H0 = 70.0;                           // km/s/Mpc
    this.Mpc_to_m = 3.086e22;                 // m/Mpc
    this.Omega_m = 0.3;                       // Matter density
    this.Omega_Lambda = 0.7;                  // Dark energy density
    this.M_sun = 1.989e30;                    // kg

    // NGC 2264 specific parameters
    this.M_visible = 80 * this.M_sun;         // kg (visible mass)
    this.M_DM = 20 * this.M_sun;              // kg (dark matter)
    this.M = this.M_visible + this.M_DM;      // Total mass
    this.M0 = this.M;                         // Initial mass
    this.SFR = 0.01 * this.M_sun / this.year_to_s;  // kg/s (~0.01 M☉/yr)
    this.r = 3.31e16;                         // m (~3.5 ly)
    this.z = 0.0008;                          // Redshift

    // Stellar wind parameters
    this.v_wind = 20e3;                       // m/s (20 km/s from massive star)
    this.v_r = 1e3;                           // m/s (radial velocity)
    this.t = 3e6 * this.year_to_s;            // s (default: 3 Myr)

    // Dust and fluid parameters
    this.rho_fluid = 1e-20;                   // kg/m³
    this.V = 1e48;                            // m³
    this.B = 1e-5;                            // T (magnetic field)
    this.B_crit = 1e11;                       // T (critical field)
    this.Delta_x = 1e-10;                     // m
    this.Delta_p = this.hbar / this.Delta_x;  // kg·m/s
    this.integral_psi = 1.0;                  // Normalized

    // Pillar wave parameters
    this.A = 1e-10;                           // Amplitude
    this.k = 1e20;                            // m⁻¹ (wavenumber)
    this.omega = 1e-14;                       // rad/s (pillar oscillation frequency)
    this.x = 0.0;                             // Position
    this.v = this.v_wind;                     // m/s (velocity)
    this.sigma = 1e15;                        // m (Gaussian envelope scale)

    // Ug subcomponent scalars
    this.Ug1 = 0.0;                           // Dipole gravity
    this.Ug2 = 0.0;                           // Superconductor
    this.Ug3 = 0.0;                           // External (stellar wind)
    this.Ug4 = 0.0;                           // Reaction

    // Dipole and magnetism parameters
    this.I_dipole = 1e18;                     // A (dipole moment current)
    this.A_dipole = 1e12;                     // m² (dipole area)
    this.H_aether = 1e-6;                     // A/m (aether magnetic field)
    this.omega_spin = 1e-5;                   // rad/s (protostar spin)

    // Aether vacuum integration
    this.rho_vac_SCm = 7.09e-37;              // J/m³
    this.rho_vac_UA = 7.09e-36;               // J/m³
    this.lambda_I = 1.0;                      // Aether coupling
    this.omega_i = 1e-8;                      // rad/s (aether frequency)
    this.t_n = 0.0;                           // Aether time
    this.F_RZ = 0.01;                         // Time-reversal factor
    this.k_4 = 1.0;                           // Reaction coupling

    // Environmental and correction factors
    this.k_SF = 1e-10;                        // N/Msun (star formation feedback)
    this.scale_macro = 1e-12;
    this.f_TRZ = 0.1;                         // Time-reversal zone factor
    this.f_sc = 1.0;                          // Superconductor factor
    this.delta_rho_over_rho = 1e-5;           // Density perturbation

    // Dependent parameters (updated on initialization)
    this.rho = this.rho_fluid;
  }

  /**
   * Update variable dynamically with cascade updates
   * 
   * @param {string} name - Variable name
   * @param {number} value - New value
   */
  updateVariable(name, value) {
    this[name] = value;

    // Cascade updates for dependent parameters
    if (name === 'Delta_x') {
      this.Delta_p = this.hbar / value;
    } else if (name === 'M_visible' || name === 'M_DM') {
      this.M = this.M_visible + this.M_DM;
      this.M0 = this.M;
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
      this.Omega_m * Math.pow(1.0 + z_val, 3) + this.Omega_Lambda
    );
    return (Hz_kms * 1e3) / this.Mpc_to_m;
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
   * Compute radius evolution with radial velocity
   * 
   * @param {number} t - Time in seconds
   * @returns {number} r(t) = r0 + v_r·t
   */
  computeRt(t) {
    return this.r + this.v_r * t;
  }

  /**
   * Compute environmental forcing (wind + star formation + erosion)
   * 
   * @param {number} t - Time in seconds
   * @returns {number} F_env total environmental forcing
   */
  computeFenv(t) {
    // Stellar wind ram pressure
    const F_wind = this.rho_fluid * Math.pow(this.v_wind, 2);
    
    // Star formation feedback
    const F_SF = this.k_SF * this.SFR / this.M_sun;
    
    // Pillar erosion factor: 5% per 3 Myr
    const F_erode = 0.05 * (t / (3e6 * this.year_to_s));
    
    return F_wind + F_SF + F_erode;
  }

  /**
   * Compute Ug1: Magnetic dipole gravity
   * 
   * @param {number} t - Time in seconds
   * @returns {number} Ug1 term
   */
  computeUg1(t) {
    const mu_dipole = this.I_dipole * this.A_dipole * this.omega_spin;
    return mu_dipole * this.B;
  }

  /**
   * Compute Ug2: Superconductor magnetic pressure
   * 
   * @param {number} t - Time in seconds
   * @returns {number} Ug2 term
   */
  computeUg2(t) {
    const B_super = this.mu_0 * this.H_aether;
    return (B_super * B_super) / (2 * this.mu_0);
  }

  /**
   * Compute Ug3': External gravity (stellar wind from massive star)
   * 
   * @param {number} t - Time in seconds
   * @returns {number} Ug3' term
   */
  computeUg3prime(t) {
    const M_star = 20 * this.M_sun;  // ~20 M☉ O/B star
    const r_star = 1e10;             // m (~0.1 AU approximation)
    return (this.G * M_star) / (r_star * r_star);
  }

  /**
   * Compute Ug4: Reaction energy term
   * 
   * @param {number} t - Time in seconds
   * @returns {number} Ug4 term
   */
  computeUg4(t) {
    // Reaction energy with decay timescale appropriate for star-forming regions
    // E_react ~ gravitational binding energy scale (not exponential cosmic energy)
    const tau_decay = 3e6 * this.year_to_s;  // ~3 Myr decay timescale
    const E_base = 1e15;  // Reasonable energy scale for 100 M☉ nebula
    const E_react = E_base * Math.exp(-t / tau_decay);
    return this.k_4 * E_react;
  }

  /**
   * Compute Ui: Aether vacuum integration term
   * 
   * @param {number} t - Time in seconds
   * @returns {number} Ui term
   */
  computeUi(t) {
    const ratio = this.rho_vac_SCm / this.rho_vac_UA;
    const cos_term = Math.cos(this.pi * this.t_n);
    return this.lambda_I * ratio * this.omega_i * cos_term * (1 + this.F_RZ);
  }

  /**
   * Compute psi integral for pillar wave structure
   * Pillar wave function: ψ = A·exp(-r²/(2σ²))·exp(i(m·θ - ω·t))
   * 
   * @param {number} r - Radius in meters
   * @param {number} t - Time in seconds
   * @returns {number} |ψ|² (intensity)
   */
  computePsiIntegral(r, t) {
    // Gaussian envelope for pillar profile
    const envelope = Math.exp(-(r * r) / (2 * this.sigma * this.sigma));
    
    // Pillar m-mode (m=1 fundamental) with time evolution
    // Note: m=0 in original means no angular dependence
    const phase = -this.omega * t;
    const real_part = this.A * envelope * Math.cos(phase);
    const imag_part = this.A * envelope * Math.sin(phase);
    
    // |ψ|² = Re² + Im²
    return (real_part * real_part) + (imag_part * imag_part);
  }

  /**
   * Compute quantum gravity term
   * 
   * @param {number} t_Hubble_val - Hubble time in seconds
   * @param {number} r - Radius in meters
   * @returns {number} Quantum correction term
   */
  computeQuantumTerm(t_Hubble_val, r) {
    const unc = Math.sqrt(this.Delta_x * this.Delta_p);
    const psi_int = this.computePsiIntegral(r, this.t);
    return (this.hbar / unc) * this.integral_psi * (2 * this.pi / t_Hubble_val) * psi_int;
  }

  /**
   * Compute fluid dynamics term (wind + dust interaction)
   * 
   * @param {number} g_base - Base gravitational acceleration
   * @returns {number} Fluid dynamics contribution
   */
  computeFluidTerm(g_base) {
    return this.rho_fluid * this.V * g_base;
  }

  /**
   * Compute dark matter response term
   * 
   * @param {number} r - Radius in meters
   * @returns {number} Dark matter contribution
   */
  computeDMTerm(r) {
    const pert = this.delta_rho_over_rho;
    const curv = 3 * this.G * this.M / (r * r * r);
    return (this.M_visible + this.M_DM) * (pert + curv);
  }

  /**
   * Compute sum of Ug components
   * Returns total Ug gravity including base gravity term
   * 
   * @param {number} r - Radius in meters
   * @returns {number} Sum of all Ug components
   */
  computeUgSum(r) {
    const Ug_base = (this.G * this.M) / (r * r);
    this.Ug1 = this.computeUg1(this.t);
    this.Ug2 = this.computeUg2(this.t);
    this.Ug3 = this.computeUg3prime(this.t);
    this.Ug4 = this.computeUg4(this.t);
    
    return Ug_base + this.Ug1 + this.Ug2 + this.Ug3 + this.Ug4;
  }

  /**
   * Compute full gravitational field g_NGC2264(r, t)
   * Core function: integrates all physics components with pillar dynamics
   * 
   * @param {number} t - Time in seconds
   * @param {number} r - Radius in meters (optional; if 0, uses default this.r)
   * @returns {number} g_NGC2264(r, t) in m/s²
   */
  computeG(t, r = 0) {
    // Update radius if provided
    if (r > 0) {
      this.r = r;
    } else {
      r = this.r;
    }
    this.t = t;

    // Star formation mass growth
    const msf_factor = this.computeMsfFactor(t);
    const m_factor = 1.0 + msf_factor;

    // Cosmological expansion
    const Hz = this.computeHtz(this.z);
    const expansion = 1.0 + Hz * t;

    // Magnetic suppression
    const sc_correction = 1.0 - (this.B / this.B_crit);

    // Environmental forcing (wind + SF + erosion)
    const f_env = this.computeFenv(t);

    // Time-reversal zone factor
    const tr_factor = 1.0 + this.f_TRZ;

    // Base gravity with all modulating factors
    const g_base = (this.G * this.M * m_factor / (r * r)) * 
                   expansion * sc_correction * (1.0 + f_env) * tr_factor;

    // Ug sum (includes base; subtract to avoid double-count)
    const ug_sum = this.computeUgSum(r) - g_base;

    // Cosmological constant term
    const lambda_term = this.Lambda * (this.c * this.c) / 3.0;

    // Aether vacuum integration
    const ui_term = this.computeUi(t);

    // Quantum gravity with pillar structure
    const quantum_term = this.computeQuantumTerm(this.t_Hubble, r);

    // Fluid dynamics (stellar wind + dust)
    const fluid_term = this.computeFluidTerm(g_base);

    // Dark matter response
    const dm_term = this.computeDMTerm(r);

    // Total unified field
    return g_base + ug_sum + lambda_term + ui_term + quantum_term + fluid_term + dm_term;
  }

  /**
   * Compute g_NGC2264 at standard system radius
   * 
   * @param {number} t - Time in seconds
   * @returns {number} g_NGC2264 at standard radius
   */
  computeGStandard(t) {
    return this.computeG(t, this.r);
  }

  /**
   * Get complete equation description with all components
   * 
   * @returns {string} Formatted equation text
   */
  getEquationText() {
    let eq = `g_NGC2264(r,t) = [G·M(t)/r²] · (1 + H(t,z)) · (1 - B/B_crit) · (1 + F_env(t)) · (1 + f_TRZ)
              + [Ug1 + Ug2 + Ug3' + Ug4] + U_i
              + (Λ·c²/3)
              + (ℏ/√(Δx·Δp)) · ∫|ψ_pillar|² dV · (2π/t_Hubble)
              + ρ_fluid · V · g + (M_visible + M_DM) · (Δρ/ρ + 3GM/r³)

Where:
  M(t) = M₀(1 + SFR·t/M₀)                          [Star formation mass growth]
  H(t,z) = H₀√(Ω_m(1+z)³ + Ω_Λ)                    [Hubble parameter]
  F_env(t) = ρ·v_wind² + k_SF·SFR + 0.05·(t/3Myr) [Wind + SF + erosion]
  
  Ug1 = μ_dipole·B                                 [Magnetic dipole]
  Ug2 = B_super²/(2μ₀)                             [Magnetic pressure]
  Ug3' = G·M_star/r_star²                          [Stellar wind external gravity]
  Ug4 = k₄·E_react(t) = k₄·E₀·exp(-0.0005·t)      [Reaction decay]
  
  U_i = λ_I·(ρ_SCm/ρ_UA)·ω_i·cos(π·t_n)·(1+F_RZ) [Aether vacuum integration]
  
  ψ_pillar = A·exp(-r²/(2σ²))·exp(i(-ω·t))        [Pillar wave: m=1 fundamental]
  |ψ|² = A²·exp(-r²/σ²)·(1 - cos(2ω·t))/2         [Pillar intensity]

Physics Domain: NGC 2264 Cone Nebula
- Stellar winds from 20 M☉ O/B star (v=20 km/s)
- Pillar erosion with 5% progression per 3 Myr
- Protostar formation and spin-up (ω_spin=1e-5 rad/s)
- Dust/gas density and scattering
- Dark matter response to tidal forces
- Quantum pillar wave structure

Reference: Hubble ACS 2002 observations
Result: g_NGC2264 ~ 1e-9 m/s² at t=3 Myr
        (Wind and fluid dynamics dominant; repulsive terms advance UQFF framework)
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
      name: 'NGC 2264 Cone Nebula UQFF Module',
      physicsType: 'Star-forming nebula with wind and pillar dynamics',
      current_system: 'NGC2264',
      mass_visible_Msun: this.M_visible / this.M_sun,
      mass_DM_Msun: this.M_DM / this.M_sun,
      mass_total_Msun: this.M / this.M_sun,
      radius_m: this.r,
      radius_ly: this.r / 9.461e15,
      starFormationRate_Msun_yr: this.SFR / (this.M_sun / this.year_to_s),
      redshift: this.z,
      stellarWindVelocity_km_s: this.v_wind / 1e3,
      magnetic_field_T: this.B,
      protostarSpin_rad_s: this.omega_spin,
      pillarErosionRate_percent_per_3Myr: 5.0,
      typicalOutput_ms2: this.computeGStandard(this.t),
      description: 'Comprehensive model of NGC 2264 Cone Nebula with stellar wind ablation, pillar erosion, protostar formation, and quantum wave structure',
      physicsComponents: [
        'Base gravity with mass evolution',
        'Hubble expansion',
        'Magnetic suppression',
        'Environmental forcing (wind + SF + erosion)',
        'Time-reversal effects',
        'Ug1: Magnetic dipole',
        'Ug2: Superconductor pressure',
        'Ug3\': Stellar wind external gravity',
        'Ug4: Reaction energy decay',
        'Ui: Aether vacuum integration',
        'Quantum pillar waves (m=1 mode)',
        'Cosmological constant',
        'Fluid dynamics (dust/wind)',
        'Dark matter response'
      ]
    };
  }

  /**
   * Print formatted variable summary
   */
  printVariables() {
    console.log(`\n=== NGC 2264 UQFF Module ===`);
    console.log(`System: NGC 2264 Cone Nebula\n`);
    
    console.log(`Universal Constants:
  G=${this.G.toExponential(6)} m³ kg⁻¹ s⁻²
  c=${this.c.toExponential(6)} m/s
  ℏ=${this.hbar.toExponential(6)} J·s
  t_Hubble=${(this.t_Hubble / this.year_to_s / 1e9).toFixed(1)} Gyr
  
System Mass (Initial):
  M_visible=${(this.M_visible / this.M_sun).toExponential(2)} M☉
  M_DM=${(this.M_DM / this.M_sun).toExponential(2)} M☉
  M_total=${(this.M / this.M_sun).toExponential(2)} M☉

Geometry:
  r=${this.r.toExponential(6)} m (~${(this.r / 9.461e15).toFixed(2)} ly)
  v_r=${(this.v_r / 1e3).toFixed(1)} km/s (radial expansion)

Star Formation & Wind:
  SFR=${(this.SFR / (this.M_sun / this.year_to_s)).toExponential(3)} M☉/yr
  v_wind=${(this.v_wind / 1e3).toFixed(0)} km/s
  ω_spin=${this.omega_spin.toExponential(2)} rad/s (protostar)

Dust/Gas:
  ρ_fluid=${this.rho_fluid.toExponential(2)} kg/m³
  V=${this.V.toExponential(2)} m³

Magnetism:
  B=${this.B.toExponential(2)} T
  B_crit=${this.B_crit.toExponential(2)} T

Pillar Dynamics:
  A=${this.A.toExponential(2)} (amplitude)
  ω=${this.omega.toExponential(2)} rad/s (pillar frequency)
  σ=${this.sigma.toExponential(2)} m (Gaussian scale)
  Erosion: 5% per 3 Myr

Quantum:
  Δx=${this.Delta_x.toExponential(2)} m
  Δp=${this.Delta_p.toExponential(2)} kg·m/s

Cosmology:
  z=${this.z} (redshift)
  H₀=${this.H0} km/s/Mpc
  Ω_m=${this.Omega_m}, Ω_Λ=${this.Omega_Lambda}
`);
  }
}

// Module export for integration
module.exports = NGC2264UQFFModule;
