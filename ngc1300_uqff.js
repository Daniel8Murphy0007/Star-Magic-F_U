/**
 * NGC 1300 Barred Spiral Galaxy UQFF Module (JavaScript Port)
 * 
 * Modular implementation of the Master Universal Gravity Equation (MUGE & UQFF Integration)
 * for NGC 1300 barred spiral galaxy dynamics.
 * 
 * This module models gravitational dynamics incorporating:
 * - Bar-driven gas funneling (active nucleus dynamics)
 * - Spiral arm density waves (m=2 mode, 200 km/s gas velocity)
 * - Star formation (SFR=1 M☉/yr)
 * - Dark matter component (3×10¹⁰ M☉)
 * - Cosmological expansion (Hubble correction)
 * - Advanced UQFF gravity components (Ug1-Ug4)
 * 
 * Physics Domain: Barred Spiral Galaxy Gravitational Dynamics
 * Astrophysical Context: NGC 1300 (M=1×10¹¹ M☉, z=0.005, SFR=1 M☉/yr)
 * Distance: 11.79 Mpc (from Hubble data)
 * 
 * Original C++ Implementation: Source73.cpp
 * JavaScript Port: November 1, 2025
 * 
 * Copyright - Daniel T. Murphy, analyzed Oct 10, 2025.
 */

class NGC1300UQFFModule {
  /**
   * Constructor: Initialize with NGC 1300 barred spiral galaxy defaults
   */
  constructor() {
    // Universal constants
    this.G = 6.6743e-11;                      // m³ kg⁻¹ s⁻² (gravitational constant)
    this.c = 3e8;                             // m/s (speed of light)
    this.hbar = 1.0546e-34;                   // J s (reduced Planck constant)
    this.Lambda = 1.1e-52;                    // m⁻² (cosmological constant)
    this.q = 1.602e-19;                       // C (elementary charge)
    this.pi = Math.PI;                        // π

    // Cosmological parameters
    this.t_Hubble = 13.8e9 * 3.156e7;         // s (Hubble time, ~13.8 Gyr)
    this.year_to_s = 3.156e7;                 // s/yr
    this.H0 = 70.0;                           // km/s/Mpc (Hubble constant)
    this.Mpc_to_m = 3.086e22;                 // m/Mpc
    this.Omega_m = 0.3;                       // Matter density parameter
    this.Omega_Lambda = 0.7;                  // Dark energy density parameter

    // Solar/distance conversions
    const M_sun_val = 1.989e30;               // kg
    const kpc_val = 3.086e19;                 // m

    // NGC 1300 astrophysical parameters
    this.M_visible = 7e10 * M_sun_val;        // kg (visible mass)
    this.M_DM = 3e10 * M_sun_val;             // kg (dark matter: 3×10¹⁰ M☉)
    this.M = this.M_visible + this.M_DM;      // kg (total initial mass: 1×10¹¹ M☉)
    this.M0 = this.M;                         // Reference mass
    this.SFR = 1 * M_sun_val / this.year_to_s; // kg/s (star formation rate: 1 M☉/yr)
    this.r = 11.79e3 * kpc_val;               // m (galaxy radius: 11.79 kpc)
    this.z = 0.005;                           // Redshift
    this.v_arm = 200e3;                       // m/s (spiral arm gas velocity: 200 km/s)
    this.t = 1e9 * this.year_to_s;            // Default t = 1 Gyr in seconds

    // Fluid/magnetic dynamics
    this.rho_fluid = 1e-21;                   // kg/m³ (fluid density)
    this.V = 1e50;                            // m³ (volume element)
    this.B = 1e-5;                            // T (magnetic field)
    this.B_crit = 1e11;                       // T (critical magnetic field)
    this.Delta_x = 1e-10;                     // m (position uncertainty)
    this.Delta_p = this.hbar / this.Delta_x;  // kg·m/s (momentum uncertainty, Heisenberg)

    // Spiral arm parameters (m=2 mode density waves)
    this.A = 1e-10;                           // Amplitude
    this.k = 1e20;                            // m⁻¹ (wavenumber)
    this.omega = 1e-15;                       // rad/s (density wave frequency)
    this.x = 0.0;                             // Position (used in wave calculation)
    this.v = this.v_arm;                      // m/s (velocity)
    this.sigma = 1e3 * kpc_val;               // m (Gaussian width: 1 kpc)
    this.integral_psi = 1.0;                  // Normalized wave integral

    // Ug subcomponents initialization
    this.Ug1 = 0.0;                           // Dipole term
    this.Ug2 = 0.0;                           // Superconductor term
    this.Ug3 = 0.0;                           // External (bar) term
    this.Ug4 = 0.0;                           // Reaction term
    this.Ui = 0.0;                            // Internal universal term

    // Magnetic/vacuum parameters
    this.mu_0 = 4 * this.pi * 1e-7;           // H/m (vacuum permeability)
    this.rho_vac_SCm = 7.09e-37;              // J/m³ (superconductive material vacuum energy)
    this.rho_vac_UA = 7.09e-36;               // J/m³ (universal aether vacuum energy)
    this.lambda_I = 1.0;                      // Internal scaling factor
    this.omega_i = 1e-8;                      // rad/s (internal oscillation frequency)
    this.t_n = 0.0;                           // Phase
    this.F_RZ = 0.01;                         // Reaction zone factor
    this.k_4 = 1.0;                           // Ug4 scaling

    // Environmental/bar forcing
    this.k_SF = 1e-10;                        // Star formation force scale (N/M☉ → m/s²)
    this.omega_spin = 1e-4;                   // rad/s (bar rotation rate)
    this.I_dipole = 1e20;                     // A (dipole current)
    this.A_dipole = 1e15;                     // m² (dipole area)
    this.H_aether = 1e-6;                     // A/m (aether magnetic field)
    this.delta_rho_over_rho = 1e-5;           // Density perturbation

    // Scales and corrections
    this.scale_macro = 1e-12;                 // Macroscopic scale factor
    this.f_TRZ = 0.1;                         // Time-reversal zone factor
    this.f_sc = 1.0;                          // Superconductor correction
    this.v_r = 1e3;                           // m/s (radial expansion velocity)
    this.rho = this.rho_fluid;                // Current density
  }

  /**
   * Update a variable dynamically with cascading updates for dependent parameters
   * @param {string} name - Variable name
   * @param {number} value - New value
   */
  updateVariable(name, value) {
    if (this.hasOwnProperty(name)) {
      this[name] = value;
      
      // Cascade updates for dependent parameters
      if (name === 'Delta_x') {
        this.Delta_p = this.hbar / value;
      } else if (name === 'M') {
        this.M_visible = 0.7 * value;
        this.M_DM = 0.3 * value;
        this.M0 = value;
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
   * Compute Hubble parameter H(t, z)
   * Describes cosmological expansion rate
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
   * Compute mass factor from star formation
   * M(t) = M₀(1 + SFR·t/M₀)
   * 
   * @param {number} t - Time in seconds
   * @returns {number} Fractional mass increase
   */
  computeMsfFactor(t) {
    return this.SFR * t / this.M0;
  }

  /**
   * Compute galaxy radius evolution with expansion
   * r(t) = r₀ + v_r·t
   * 
   * @param {number} t - Time in seconds
   * @returns {number} Galaxy radius at time t (m)
   */
  computeRt(t) {
    return this.r + this.v_r * t;
  }

  /**
   * Compute environmental forcing (bar + star formation + density waves)
   * F_env = F_bar + F_SF + F_wave
   * 
   * @param {number} t - Time in seconds
   * @returns {number} Environmental force factor
   */
  computeFenv(t) {
    // Bar funneling force (10% of central gravity)
    const F_bar = 0.1 * (this.G * this.M) / (this.r * this.r);
    
    // Star formation feedback force
    const F_SF = this.k_SF * this.SFR / 1.989e30;
    
    // Density wave pressure (from spiral arms)
    const F_wave = this.rho_fluid * this.v_arm * this.v_arm;
    
    return F_bar + F_SF + F_wave;
  }

  /**
   * Compute Ug1: Dipole gravity component
   * Related to stellar/galactic magnetic dipole moment
   * 
   * @param {number} t - Time in seconds
   * @returns {number} Ug1 component (m/s²)
   */
  computeUg1(t) {
    const mu_dipole = this.I_dipole * this.A_dipole * this.omega_spin;
    return mu_dipole * this.B;
  }

  /**
   * Compute Ug2: Superconductor gravity component
   * Magnetic energy density from aetheric superconductor
   * 
   * @param {number} t - Time in seconds
   * @returns {number} Ug2 component (m/s²)
   */
  computeUg2(t) {
    const B_super = this.mu_0 * this.H_aether;
    return (B_super * B_super) / (2 * this.mu_0);
  }

  /**
   * Compute Ug3: External gravity (bar-driven)
   * Models non-axisymmetric perturbation from galactic bar
   * 
   * @param {number} t - Time in seconds
   * @returns {number} Ug3 component (m/s²)
   */
  computeUg3prime(t) {
    const M_bar = 0.2 * this.M;          // Bar mass: 20% of galaxy
    const r_bar = 0.3 * this.r;          // Bar radius: 30% of galaxy radius
    return (this.G * M_bar) / (r_bar * r_bar);
  }

  /**
   * Compute Ug4: Reaction gravity component
   * Energy release from nuclear reactions (decaying with time)
   * 
   * @param {number} t - Time in seconds
   * @returns {number} Ug4 component (m/s²)
   */
  computeUg4(t) {
    const E_react = 1e46 * Math.exp(-0.0005 * t);
    return this.k_4 * E_react;
  }

  /**
   * Compute Ui: Internal universal gravity component
   * Oscillatory aetheric contribution with time-reversal correction
   * 
   * @param {number} t - Time in seconds
   * @returns {number} Ui component (m/s²)
   */
  computeUi(t) {
    return this.lambda_I * (this.rho_vac_SCm / this.rho_vac_UA) * 
           this.omega_i * Math.cos(this.pi * this.t_n) * (1 + this.F_RZ);
  }

  /**
   * Compute spiral arm density wave perturbation
   * ψ_spiral = A exp(-r²/2σ²) exp(i(m·φ - ω·t))
   * Returns |ψ|² for probability density
   * 
   * @param {number} r - Galactic radius (m)
   * @param {number} t - Time (s)
   * @returns {number} |ψ|² (dimensionless perturbation strength)
   */
  computePsiIntegral(r, t) {
    const m = 2.0;  // m=2 mode (typical for barred spirals)
    
    // Gaussian envelope: exp(-r²/2σ²)
    const gaussian = Math.exp(-(r * r) / (2 * this.sigma * this.sigma));
    
    // Spiral phase: exp(i(m·φ - ω·t))
    // Real part: cos(m·φ - ω·t), Imaginary part: sin(m·φ - ω·t)
    // For φ=0: real = cos(-ω·t) = cos(ω·t)
    const phase = Math.cos(this.omega * t);
    
    // Full wave: A·gaussian·exp(i·phase)
    // |ψ|² = A²·gaussian²·1 (since |exp(iθ)| = 1)
    const psi_squared = (this.A * this.A) * (gaussian * gaussian) * (phase * phase);
    
    return psi_squared;
  }

  /**
   * Compute quantum gravity term
   * Incorporates Heisenberg uncertainty and wave function evolution
   * 
   * @param {number} t_Hubble_val - Hubble time (s)
   * @param {number} r - Galactic radius (m)
   * @returns {number} Quantum correction (m/s²)
   */
  computeQuantumTerm(t_Hubble_val, r) {
    const unc = Math.sqrt(this.Delta_x * this.Delta_p);
    const psi_int = this.computePsiIntegral(r, this.t);
    
    return (this.hbar / unc) * this.integral_psi * 
           (2 * this.pi / t_Hubble_val) * psi_int;
  }

  /**
   * Compute fluid dynamics contribution
   * Pressure forces from gas in galaxy
   * 
   * @param {number} g_base - Base gravitational acceleration (m/s²)
   * @returns {number} Fluid term (m/s²)
   */
  computeFluidTerm(g_base) {
    return this.rho_fluid * this.V * g_base;
  }

  /**
   * Compute dark matter contribution
   * Includes density perturbations and gravitational response
   * 
   * @param {number} r - Galactic radius (m)
   * @returns {number} Dark matter term (m/s²)
   */
  computeDMTerm(r) {
    const pert = this.delta_rho_over_rho;
    const curv = 3 * this.G * this.M / (r * r * r);
    
    return (this.M_visible + this.M_DM) * (pert + curv);
  }

  /**
   * Compute sum of Ug components
   * Updates individual Ug values and returns total
   * 
   * @param {number} r - Galactic radius (m)
   * @returns {number} Sum of Ug1 + Ug2 + Ug3 + Ug4 (m/s²)
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
   * Compute full gravitational acceleration for NGC 1300
   * Complete UQFF equation with all terms:
   * g = g_base + Ug_sum + λ_term + Ui + quantum + fluid + DM
   * 
   * @param {number} t - Time in seconds
   * @param {number} r - Galactic radius in meters
   * @returns {number} Total gravitational acceleration (m/s²)
   */
  computeG(t, r) {
    this.t = t;
    
    // Star formation mass growth factor
    const msf_factor = this.computeMsfFactor(t);
    const m_factor = 1.0 + msf_factor;
    
    // Cosmological expansion factor
    const Hz = this.computeHtz(this.z);
    const expansion = 1.0 + Hz * t;
    
    // Magnetic field suppression factor
    const sc_correction = 1.0 - (this.B / this.B_crit);
    
    // Environmental forcing
    const f_env = this.computeFenv(t);
    
    // Time-reversal correction
    const tr_factor = 1.0 + this.f_TRZ;
    
    // Base gravitational acceleration with all multiplicative corrections
    const g_base = (this.G * this.M * m_factor / (r * r)) * 
                   expansion * sc_correction * (1.0 + f_env) * tr_factor;
    
    // Ug sum (subtract base to avoid double-counting)
    const ug_sum = this.computeUgSum(r) - g_base;
    
    // Cosmological constant term
    const lambda_term = this.Lambda * (this.c * this.c) / 3.0;
    
    // Internal universal term
    const ui_term = this.computeUi(t);
    
    // Quantum gravity correction
    const quantum_term = this.computeQuantumTerm(this.t_Hubble, r);
    
    // Fluid dynamics
    const fluid_term = this.computeFluidTerm(g_base);
    
    // Dark matter response
    const dm_term = this.computeDMTerm(r);
    
    // Total: all components summed
    return g_base + ug_sum + lambda_term + ui_term + quantum_term + fluid_term + dm_term;
  }

  /**
   * Compute gravitational acceleration at standard radius (11.79 kpc)
   * 
   * @param {number} t - Time in seconds
   * @returns {number} g_NGC1300(t) at standard radius (m/s²)
   */
  computeGStandard(t) {
    return this.computeG(t, this.r);
  }

  /**
   * Get complete equation description with all components
   * 
   * @returns {string} Formatted equation text with physical interpretation
   */
  getEquationText() {
    return `g_NGC1300(r, t) = [G·M(t)/r²] · (1 + H(t,z)) · (1 - B/B_crit) · (1 + F_env(t)) · (1 + f_TRZ) 
                 + (Ug1 + Ug2 + Ug3' + Ug4) + Ui + (Λ·c²/3) 
                 + (ℏ/√(ΔxΔp)) · ∫|ψ_total|² dV · (2π/t_Hubble) 
                 + ρ_fluid · V · g + (M_visible + M_DM) · (Δρ/ρ + 3GM/r³)

Where:
  M(t) = M₀(1 + SFR·t/M₀) [Star formation growth]
  r(t) = r₀ + v_r·t [Galactic expansion]
  H(t,z) = H₀·√(Ω_m(1+z)³ + Ω_Λ) [Hubble parameter]
  
  F_env(t) = F_bar + F_SF + F_wave [Environmental forcing]
    F_bar = 0.1·G·M/r² [Bar-driven funneling]
    F_wave = ρ_fluid·v_arm² [Density wave pressure]
  
  Ug1 = μ_dipole·B [Dipole gravity]
  Ug2 = B_super²/(2μ₀) [Superconductor gravity]
  Ug3' = G·M_bar/r_bar² [External bar gravity, M_bar=0.2M, r_bar=0.3r]
  Ug4 = k₄·E_react(t)·exp(-0.0005·t) [Reaction gravity, decaying]
  Ui = λ_I·(ρ_vac,SCm/ρ_vac,UA)·ω_i·cos(π·t_n)·(1+F_RZ) [Internal universal term]
  
  ψ_spiral = A·exp(-r²/2σ²)·exp(i(m·φ-ω·t)) [m=2 density wave mode]
  
Physics Interpretation:
  • Attractive forces (g_base, Ug1, Ug3') drive central concentration
  • Repulsive forces (Ug2, Λ) oppose collapse, enable disk stability
  • Environmental terms (F_bar, F_wave) drive gas dynamics and star formation
  • Quantum term reflects wave nature of galactic matter distribution
  • Dark matter and fluid terms represent extended halo + ISM coupling
  • Time-reversal zone (f_TRZ) breaks causality locally, enhancing dynamics
  • Star formation (SFR) grows central mass over ~1 Gyr timescale
  
Astronomical Context:
  NGC 1300: M = 1×10¹¹ M☉, r = 11.79 kpc, z = 0.005, SFR = 1 M☉/yr
  Visible: 7×10¹⁰ M☉, Dark Matter: 3×10¹⁰ M☉
  Bar: 20% of mass, 30% of radius
  Spiral arms: m=2 mode, velocity 200 km/s
  Result: g_NGC1300 ~ 2×10⁻¹⁰ m/s² at t=1 Gyr (DM/fluid dominated)
  
UQFF Insights:
  Bar and spiral arms are external symmetry-breaking perturbations
  Dark matter halo (Ug3, DM term) stabilizes against collapse
  Cosmological expansion (H term) provides universal tension
  Reaction gravity (Ug4) decays on Gyr timescale
  Framework predicts enhanced star formation in bar region`;
  }

  /**
   * Get all current variables and their values
   * 
   * @returns {Object} Object containing all module variables
   */
  getVariables() {
    return {
      G: this.G,
      c: this.c,
      hbar: this.hbar,
      Lambda: this.Lambda,
      q: this.q,
      pi: this.pi,
      t_Hubble: this.t_Hubble,
      year_to_s: this.year_to_s,
      H0: this.H0,
      Mpc_to_m: this.Mpc_to_m,
      Omega_m: this.Omega_m,
      Omega_Lambda: this.Omega_Lambda,
      M_visible: this.M_visible,
      M_DM: this.M_DM,
      M: this.M,
      M0: this.M0,
      SFR: this.SFR,
      r: this.r,
      z: this.z,
      v_arm: this.v_arm,
      t: this.t,
      rho_fluid: this.rho_fluid,
      V: this.V,
      B: this.B,
      B_crit: this.B_crit,
      Delta_x: this.Delta_x,
      Delta_p: this.Delta_p,
      A: this.A,
      k: this.k,
      omega: this.omega,
      x: this.x,
      v: this.v,
      sigma: this.sigma,
      integral_psi: this.integral_psi,
      Ug1: this.Ug1,
      Ug2: this.Ug2,
      Ug3: this.Ug3,
      Ug4: this.Ug4,
      Ui: this.Ui,
      mu_0: this.mu_0,
      rho_vac_SCm: this.rho_vac_SCm,
      rho_vac_UA: this.rho_vac_UA,
      lambda_I: this.lambda_I,
      omega_i: this.omega_i,
      t_n: this.t_n,
      F_RZ: this.F_RZ,
      k_4: this.k_4,
      k_SF: this.k_SF,
      omega_spin: this.omega_spin,
      I_dipole: this.I_dipole,
      A_dipole: this.A_dipole,
      H_aether: this.H_aether,
      delta_rho_over_rho: this.delta_rho_over_rho,
      scale_macro: this.scale_macro,
      f_TRZ: this.f_TRZ,
      f_sc: this.f_sc,
      v_r: this.v_r,
      rho: this.rho
    };
  }

  /**
   * Print formatted variable summary
   */
  printVariables() {
    console.log('\n=== NGC 1300 Barred Spiral Galaxy UQFF Module Variables ===\n');
    const vars = this.getVariables();
    Object.entries(vars).forEach(([name, value]) => {
      console.log(`${name.padEnd(25)} = ${value.toExponential(6)}`);
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
      name: 'NGC 1300 Barred Spiral Galaxy UQFF Module',
      physicsType: 'Barred Spiral Galaxy Gravitational Dynamics',
      totalMass_Msun: this.M / 1.989e30,
      visibleMass_Msun: this.M_visible / 1.989e30,
      darkMatterMass_Msun: this.M_DM / 1.989e30,
      starFormationRate_Msun_yr: this.SFR / (1.989e30 / this.year_to_s),
      galaxyRadius_kpc: this.r / 3.086e19 / 1e3,
      distance_Mpc: 11.79,
      redshift: this.z,
      spiralArmVelocity_km_s: this.v_arm / 1e3,
      barMass_fraction: 0.2,
      barRadius_fraction: 0.3,
      uqffComponentsActive: ['Ug1', 'Ug2', 'Ug3 (bar)', 'Ug4', 'Ui', 'Quantum', 'Fluid', 'DM'],
      densityWaveMode: 'm = 2 spiral arms',
      equation: 'g_NGC1300(r,t) = g_base(expansion, SFR, bar) + Ug_sum + Λ·c²/3 + Ui + quantum + fluid + DM',
      typicalOutput_ms2: this.computeGStandard(this.t),
      description: 'Models barred spiral galaxy with bar-driven gas funneling, density waves, star formation, and dark matter'
    };
  }
}

// Module export for integration
module.exports = NGC1300UQFFModule;
