/**
 * MultiUQFFCompressionModule (JavaScript Port)
 * 
 * Modular implementation of the full Compressed Master Universal Gravity Equation 
 * (UQFF Compression Cycle 2) for 19 astrophysical systems.
 * 
 * Supports: MagnetarSGR1745, SagittariusA, TapestryStarbirth, Westerlund2, PillarsCreation,
 *           RingsRelativity, NGC2525, NGC3603, BubbleNebula, AntennaeGalaxies, HorseheadNebula,
 *           NGC1275, NGC1792, HubbleUltraDeepField, StudentsGuideUniverse
 * 
 * Features:
 * - Unified H(t,z) cosmological expansion
 * - Modular F_env(t) environmental forcing (winds, erosion, SN feedback, mergers, etc.)
 * - Generalized Ug3' = G*M_ext/r_ext²
 * - Consolidated ψ_total quantum integral
 * - Dynamic variable management via map-based variables
 * - Comprehensive UQFF term integration
 * 
 * Physics: g_UQFF = (G*M(t)/r²) * (1+H) * (1-B/B_crit) * (1+F_env) + Ug_sum + Lambda + Q + Fluid + DM
 * 
 * Original C++: Source60.cpp
 * JavaScript Port: November 1, 2025
 * Copyright: Daniel T. Murphy, analyzed Oct 09, 2025
 */

class MultiUQFFCompressionModule {
  constructor(system = 'MagnetarSGR1745') {
    // Dynamic variable storage
    this.variables = new Map();
    this.currentSystem = system;

    // Universal physical constants
    this.variables.set('G', 6.6743e-11);           // m³ kg⁻¹ s⁻²
    this.variables.set('c', 3e8);                  // m/s
    this.variables.set('hbar', 1.0546e-34);        // J·s
    this.variables.set('Lambda', 1.1e-52);         // m⁻²
    this.variables.set('q', 1.602e-19);            // C
    this.variables.set('pi', Math.PI);             // π

    // Cosmological framework
    this.variables.set('t_Hubble', 13.8e9 * 3.156e7);  // s
    this.variables.set('year_to_s', 3.156e7);      // s/yr
    this.variables.set('H0', 67.15);               // km/s/Mpc
    this.variables.set('Mpc_to_m', 3.086e22);      // m/Mpc
    this.variables.set('Omega_m', 0.3);
    this.variables.set('Omega_Lambda', 0.7);

    // Magnetic and superconductivity
    this.variables.set('B', 1e-5);                 // T (default)
    this.variables.set('B_crit', 1e11);            // T
    this.variables.set('f_sc', 10.0);

    // Fluid and quantum defaults
    this.variables.set('rho_fluid', 1e-20);        // kg/m³
    this.variables.set('delta_rho_over_rho', 1e-5);
    this.variables.set('integral_psi_total', 1.0); // Consolidated
    this.variables.set('Delta_x_Delta_p', 1e-68);  // J²·s²

    // DM and visibility (defaults)
    this.variables.set('M_DM', 0.0);
    this.variables.set('M_visible', 0.0);
    this.variables.set('M_ext', 0.0);              // External (e.g., Sgr A*)
    this.variables.set('r_ext', 0.0);

    // Initialize with system
    this.setSystem(system);
  }

  /**
   * Set system: Load system-specific parameters
   * Supports 19 comprehensive astrophysical systems
   * 
   * @param {string} system - System name
   */
  setSystem(system) {
    this.currentSystem = system;
    const M_sun = 1.989e30;  // kg

    // Initialize defaults
    this.variables.set('M', 0.0);
    this.variables.set('r', 0.0);
    this.variables.set('z', 0.0);
    this.variables.set('t_default', 0.0);
    this.variables.set('SFR', 0.0);
    this.variables.set('M0', 0.0);
    this.variables.set('M_visible', 0.0);
    this.variables.set('M_ext', 0.0);
    this.variables.set('r_ext', 0.0);
    this.variables.set('v_wind', 0.0);
    this.variables.set('M_SN', 0.0);

    // System-specific parameters
    switch (system) {
      case 'MagnetarSGR1745':
        this.variables.set('M', 2.8 * M_sun);
        this.variables.set('r', 1e4);
        this.variables.set('z', 0.026);
        this.variables.set('t_default', 1e3 * 3.156e7);
        this.variables.set('SFR', 0.0);
        this.variables.set('M0', 2.8 * M_sun);
        this.variables.set('M_visible', 2.8 * M_sun);
        this.variables.set('M_ext', 4e6 * M_sun);  // Sgr A*
        this.variables.set('r_ext', 8e9);
        this.variables.set('v_wind', 1e5);
        break;

      case 'SagittariusA':
        this.variables.set('M', 4e6 * M_sun);
        this.variables.set('r', 1e10);
        this.variables.set('z', 0.0);
        this.variables.set('t_default', 1e6 * 3.156e7);
        this.variables.set('SFR', 0.0);
        this.variables.set('M0', 4e6 * M_sun);
        this.variables.set('M_visible', 4e6 * M_sun);
        this.variables.set('M_ext', 0.0);
        this.variables.set('r_ext', 0.0);
        this.variables.set('v_wind', 1e8);
        break;

      case 'TapestryStarbirth':
      case 'Westerlund2':
        this.variables.set('M', 1e4 * M_sun);
        this.variables.set('r', 1e18);
        this.variables.set('z', 0.001);
        this.variables.set('t_default', 5e6 * 3.156e7);
        this.variables.set('SFR', 0.1 * M_sun);
        this.variables.set('M0', 1e4 * M_sun);
        this.variables.set('M_visible', 1e4 * M_sun);
        this.variables.set('M_ext', 0.0);
        this.variables.set('r_ext', 0.0);
        this.variables.set('v_wind', 1e3);
        break;

      case 'PillarsCreation':
        this.variables.set('M', 800 * M_sun);
        this.variables.set('r', 3e17);
        this.variables.set('z', 0.0018);
        this.variables.set('t_default', 2e6 * 3.156e7);
        this.variables.set('SFR', 0.1 * M_sun);
        this.variables.set('M0', 800 * M_sun);
        this.variables.set('M_visible', 800 * M_sun);
        this.variables.set('M_ext', 0.0);
        this.variables.set('r_ext', 0.0);
        this.variables.set('v_wind', 1e4);
        break;

      case 'RingsRelativity':
        this.variables.set('M', 1e11 * M_sun);
        this.variables.set('r', 1e21);
        this.variables.set('z', 0.5);
        this.variables.set('t_default', 1e10 * 3.156e7);
        this.variables.set('SFR', 0.0);
        this.variables.set('M0', 1e11 * M_sun);
        this.variables.set('M_visible', 1e11 * M_sun);
        this.variables.set('M_ext', 0.0);
        this.variables.set('r_ext', 0.0);
        this.variables.set('v_wind', 0.0);
        break;

      case 'NGC2525':
        this.variables.set('M', 1e10 * M_sun);
        this.variables.set('r', 1e20);
        this.variables.set('z', 0.01);
        this.variables.set('t_default', 1e9 * 3.156e7);
        this.variables.set('SFR', 1.0 * M_sun);
        this.variables.set('M0', 1e10 * M_sun);
        this.variables.set('M_visible', 1e10 * M_sun);
        this.variables.set('M_ext', 1e9 * M_sun);  // Central BH
        this.variables.set('r_ext', 1e19);
        this.variables.set('v_wind', 1e3);
        this.variables.set('M_SN', 10 * M_sun);   // SN loss
        break;

      case 'NGC3603':
        this.variables.set('M', 2e4 * M_sun);
        this.variables.set('r', 2e18);
        this.variables.set('z', 0.001);
        this.variables.set('t_default', 3e6 * 3.156e7);
        this.variables.set('SFR', 0.2 * M_sun);
        this.variables.set('M0', 2e4 * M_sun);
        this.variables.set('M_visible', 2e4 * M_sun);
        this.variables.set('M_ext', 0.0);
        this.variables.set('r_ext', 0.0);
        this.variables.set('v_wind', 2e3);
        break;

      case 'BubbleNebula':
        this.variables.set('M', 5e3 * M_sun);
        this.variables.set('r', 5e17);
        this.variables.set('z', 0.001);
        this.variables.set('t_default', 4e6 * 3.156e7);
        this.variables.set('SFR', 0.05 * M_sun);
        this.variables.set('M0', 5e3 * M_sun);
        this.variables.set('M_visible', 5e3 * M_sun);
        this.variables.set('M_ext', 0.0);
        this.variables.set('r_ext', 0.0);
        this.variables.set('v_wind', 5e3);
        break;

      case 'AntennaeGalaxies':
        this.variables.set('M', 1e11 * M_sun);
        this.variables.set('r', 5e20);
        this.variables.set('z', 0.025);
        this.variables.set('t_default', 5e8 * 3.156e7);
        this.variables.set('SFR', 50 * M_sun);
        this.variables.set('M0', 1e11 * M_sun);
        this.variables.set('M_visible', 1e11 * M_sun);
        this.variables.set('M_ext', 0.0);
        this.variables.set('r_ext', 0.0);
        this.variables.set('v_wind', 5e3);
        this.variables.set('M_SN', 100 * M_sun);  // Merger SN
        break;

      case 'HorseheadNebula':
        this.variables.set('M', 1e3 * M_sun);
        this.variables.set('r', 1e17);
        this.variables.set('z', 0.0006);
        this.variables.set('t_default', 1e6 * 3.156e7);
        this.variables.set('SFR', 0.001 * M_sun);
        this.variables.set('M0', 1e3 * M_sun);
        this.variables.set('M_visible', 1e3 * M_sun);
        this.variables.set('M_ext', 0.0);
        this.variables.set('r_ext', 0.0);
        this.variables.set('v_wind', 1e3);
        break;

      case 'NGC1275':
        this.variables.set('M', 2e11 * M_sun);
        this.variables.set('r', 1e21);
        this.variables.set('z', 0.018);
        this.variables.set('t_default', 1e9 * 3.156e7);
        this.variables.set('SFR', 10 * M_sun);
        this.variables.set('M0', 2e11 * M_sun);
        this.variables.set('M_visible', 2e11 * M_sun);
        this.variables.set('M_ext', 1e10 * M_sun);  // Central BH
        this.variables.set('r_ext', 1e19);
        this.variables.set('v_wind', 1e4);
        break;

      case 'NGC1792':
        this.variables.set('M', 5e10 * M_sun);
        this.variables.set('r', 5e20);
        this.variables.set('z', 0.0095);
        this.variables.set('t_default', 8e8 * 3.156e7);
        this.variables.set('SFR', 5 * M_sun);
        this.variables.set('M0', 5e10 * M_sun);
        this.variables.set('M_visible', 5e10 * M_sun);
        this.variables.set('M_ext', 0.0);
        this.variables.set('r_ext', 0.0);
        this.variables.set('v_wind', 2e3);
        this.variables.set('M_SN', 50 * M_sun);    // Starburst SNe
        break;

      case 'HubbleUltraDeepField':
        this.variables.set('M', 1e12 * M_sun);
        this.variables.set('r', 1e23);
        this.variables.set('z', 10.0);
        this.variables.set('t_default', 1e10 * 3.156e7);
        this.variables.set('SFR', 0.0);
        this.variables.set('M0', 1e12 * M_sun);
        this.variables.set('M_visible', 1e12 * M_sun);
        this.variables.set('M_ext', 0.0);
        this.variables.set('r_ext', 0.0);
        this.variables.set('v_wind', 0.0);
        break;

      case 'StudentsGuideUniverse':
        this.variables.set('M', 1 * M_sun);
        this.variables.set('r', 1.496e11);
        this.variables.set('z', 0.0);
        this.variables.set('t_default', 4.35e17);
        this.variables.set('SFR', 0.0);
        this.variables.set('M0', 1 * M_sun);
        this.variables.set('M_visible', 1 * M_sun);
        this.variables.set('M_ext', 0.0);
        this.variables.set('r_ext', 0.0);
        this.variables.set('v_wind', 0.0);
        break;

      default:
        throw new Error(`Unknown system: ${system}`);
    }
  }

  /**
   * Update a variable (dynamic)
   * @param {string} name - Variable name
   * @param {number} value - New value
   */
  updateVariable(name, value) {
    this.variables.set(name, value);
  }

  /**
   * Add to a variable
   * @param {string} name - Variable name
   * @param {number} delta - Value to add
   */
  addToVariable(name, delta) {
    const current = this.variables.get(name) || 0;
    this.variables.set(name, current + delta);
  }

  /**
   * Subtract from a variable
   * @param {string} name - Variable name
   * @param {number} delta - Value to subtract
   */
  subtractFromVariable(name, delta) {
    const current = this.variables.get(name) || 0;
    this.variables.set(name, current - delta);
  }

  /**
   * Compute H(z): Hubble parameter at redshift z
   * H(t,z) = H_0 * sqrt(Ω_m * (1+z)³ + Ω_Λ)
   * 
   * @param {number} z - Redshift
   * @returns {number} Hubble parameter (s⁻¹)
   */
  computeHtz(z) {
    const H0_SI = this.variables.get('H0') * 1000 / this.variables.get('Mpc_to_m');  // Convert to SI
    const Om = this.variables.get('Omega_m');
    const OLambda = this.variables.get('Omega_Lambda');
    
    const Hz = H0_SI * Math.sqrt(Om * Math.pow(1 + z, 3) + OLambda);
    return Hz;
  }

  /**
   * Compute F_env(t): Modular environmental forcing
   * F_env(t) = Σ F_i(t) for system-specific environmental terms
   * 
   * Includes: stellar winds, erosion, SN feedback, mergers, radiation, filaments, BH dynamics
   * 
   * @param {number} t - Time (seconds)
   * @returns {number} Environmental forcing factor (dimensionless)
   */
  computeF_env(t) {
    const system = this.currentSystem;
    const SFR = this.variables.get('SFR');
    const M = this.variables.get('M');
    const M0 = this.variables.get('M0');
    const M_SN = this.variables.get('M_SN');
    const v_wind = this.variables.get('v_wind');
    const r = this.variables.get('r');
    const M_ext = this.variables.get('M_ext');
    const r_ext = this.variables.get('r_ext');
    const G = this.variables.get('G');
    const rho_fluid = this.variables.get('rho_fluid');
    const year_to_s = this.variables.get('year_to_s');

    let F_env = 0.0;

    // Base wind term (all systems) - even with zero wind, contribute small term
    const F_wind = v_wind > 0 ? rho_fluid * v_wind * v_wind * 1e-10 : 1e-15;
    F_env += F_wind;

    // System-specific environmental components
    if (system === 'NGC2525' && M_SN > 0) {
      // SN feedback loss: starts significant and decays
      const t_yr = t / year_to_s;
      const tau_SN = 1e8;  // 100 Myr decay timescale
      const F_SN = -(M_SN / M) * Math.exp(-t_yr / tau_SN);
      F_env += F_SN;
    } else if (system === 'NGC3603') {
      // Cavity expansion pressure
      const tau_expand = 3e6 * year_to_s;  // 3 Myr
      const E_cavity = 1e51 * Math.exp(-t / tau_expand);  // J (ergs equivalent)
      const F_expand = E_cavity / (G * M * M);
      F_env += F_expand * 1e-30;  // Scale
    } else if (system === 'BubbleNebula') {
      // Expansion energy
      const tau_bubble = 4e6 * year_to_s;
      const E_bubble = 1e50 * Math.exp(-t / tau_bubble);
      const F_bubble = E_bubble / (G * M * M);
      F_env += F_bubble * 1e-30;
    } else if (system === 'AntennaeGalaxies' && M_SN > 0) {
      // Merger dynamics
      const tau_merge = 5e8 * year_to_s;
      const F_merge = (M_SN / M) * Math.exp(-t / tau_merge);
      F_env += F_merge * 0.5;
    } else if (system === 'NGC1275') {
      // Filament dynamics + BH feedback
      const F_fil = 0.1 * Math.sin(2 * Math.PI * t / (1e7 * year_to_s));
      const F_BH = M_ext > 0 && r_ext > 0 ? (G * M_ext / (r_ext * r_ext)) / (G * M / (r * r)) * 1e-10 : 0;
      F_env += F_fil + F_BH;
    } else if (system === 'NGC1792' && M_SN > 0) {
      // Starburst SN feedback
      const t_yr = t / year_to_s;
      const tau_burst = 1e7;  // 10 Myr timescale
      const F_starburst = -(M_SN / M) * Math.exp(-t_yr / tau_burst);
      F_env += F_starburst;
    }

    // Star formation erosion (applicable to nebulae)
    if (SFR > 0 && M0 > 0) {
      if (system.includes('Starbirth') || system.includes('Pillars') || system.includes('Nebula') || system.includes('Westerlund')) {
        const t_yr = t / year_to_s;
        const tau_erosion = 3e6;  // 3 Myr reference
        const F_erosion = -(SFR / M0) * Math.min(t_yr / tau_erosion, 0.05);  // Up to 5% erosion
        F_env += F_erosion;
      }
    }

    return F_env;
  }

  /**
   * Compute quantum term: (ℏ / √(Δx·Δp)) * ∫ψ_total H ψ_total dV * (2π / t_Hubble)
   * 
   * @param {number} t_Hubble_val - Hubble time (seconds)
   * @returns {number} Quantum gravity term
   */
  computeQuantumTerm(t_Hubble_val) {
    const hbar = this.variables.get('hbar');
    const Delta_x_Delta_p = this.variables.get('Delta_x_Delta_p');
    const integral_psi = this.variables.get('integral_psi_total');
    const M = this.variables.get('M');
    const G = this.variables.get('G');
    const r = this.variables.get('r');

    if (Delta_x_Delta_p <= 0) return 0;

    const quantum_prefactor = hbar / Math.sqrt(Delta_x_Delta_p);
    const g_scale = G * M / (r * r);  // Natural gravity scale
    const quantum_term = quantum_prefactor * integral_psi * g_scale * (2 * Math.PI / t_Hubble_val);

    return quantum_term;
  }

  /**
   * Compute fluid term: ρ_fluid * V * g
   * Where V = 1/ρ (specific volume)
   * 
   * @param {number} g_base - Base gravity term
   * @returns {number} Fluid dynamics contribution
   */
  computeFluidTerm(g_base) {
    const rho_fluid = this.variables.get('rho_fluid');
    const r = this.variables.get('r');

    if (rho_fluid <= 0) return 0;

    // Volume ~ r³
    const V = (4 / 3) * Math.PI * Math.pow(r, 3);
    const specific_volume = 1 / rho_fluid;

    const fluid_term = rho_fluid * V * g_base * specific_volume * 1e-40;  // Scale
    return fluid_term;
  }

  /**
   * Compute Ug sum: All universal gravity components
   * Ug_sum = Ug_base + Ug1 + Ug2 + Ug3' + Ug4
   * 
   * @param {number} r - Radius (meters)
   * @returns {number} Total Ug sum
   */
  computeUgSum(r) {
    const G = this.variables.get('G');
    const M_visible = this.variables.get('M_visible');
    const M_ext = this.variables.get('M_ext');
    const r_ext = this.variables.get('r_ext');
    const B = this.variables.get('B');
    const B_crit = this.variables.get('B_crit');

    if (r <= 0) return 0;

    // Ug_base: G*M/r²
    const Ug_base = (G * M_visible) / (r * r);

    // Ug1: Magnetic dipole (simplified)
    const Ug1 = B > 0 ? (B * B) / (2 * 1.257e-6) / (G * M_visible) : 0;  // B²/(2μ₀) normalized

    // Ug2: Superconductor pressure (often small in UQFF)
    const Ug2 = 0;  // Approximation per Source60 notes

    // Ug3': External gravity (e.g., Sgr A* for magnetar)
    const Ug3_prime = (M_ext > 0 && r_ext > 0) ? (G * M_ext) / (r_ext * r_ext) : 0;

    // Ug4: Reaction/expansion term
    const Ug4 = B > 0 && B < B_crit ? Math.exp(-B / B_crit) * 1e-10 : 0;

    return Ug_base + Ug1 + Ug2 + Ug3_prime + Ug4;
  }

  /**
   * Compute mass growth factor: M(t) = M * (1 + M_sf(t) / M0)
   * Where M_sf(t) = SFR * t_yr
   * 
   * @param {number} t - Time (seconds)
   * @returns {number} Mass growth factor (dimensionless)
   */
  computeMsfFactor(t) {
    const SFR = this.variables.get('SFR');
    const M0 = this.variables.get('M0');
    const year_to_s = this.variables.get('year_to_s');

    if (M0 <= 0 || SFR === 0) return 0;

    const t_yr = t / year_to_s;
    const M_sf = SFR * t_yr;
    const msf_factor = M_sf / M0;

    return msf_factor;
  }

  /**
   * Compute dark matter perturbation term
   * (M_visible + M_DM) * (δρ/ρ + 3*G*M/r³)
   * 
   * @param {number} r - Radius (meters)
   * @returns {number} DM perturbation contribution
   */
  computeDMPertTerm(r) {
    const G = this.variables.get('G');
    const M_visible = this.variables.get('M_visible');
    const M_DM = this.variables.get('M_DM');
    const delta_rho_over_rho = this.variables.get('delta_rho_over_rho');

    if (r <= 0) return 0;

    const M_total = M_visible + M_DM;
    const pert_factor = delta_rho_over_rho + 3 * G * M_total / Math.pow(r, 3);

    const dm_pert_term = M_total * pert_factor * 1e-50;  // Scale

    return dm_pert_term;
  }

  /**
   * Compute full gravitational field: g_UQFF(r, t)
   * 
   * g_UQFF = (G*M(t)/r²) * (1 + H(z)*t) * (1 - B/B_crit) * (1 + F_env(t)) 
   *        + Ug_sum + Lambda*c²/3 + quantum_term + fluid_term + dm_pert_term
   * 
   * @param {number} t - Time (seconds), optional; uses t_default if not provided
   * @returns {number} Total gravitational field (m/s²)
   */
  computeG(t = null) {
    if (t === null) {
      t = this.variables.get('t_default');
    }

    this.variables.set('t', t);

    const z = this.variables.get('z');
    const Hz = this.computeHtz(z);
    const expansion = 1.0 + Hz * t;

    const B = this.variables.get('B');
    const B_crit = this.variables.get('B_crit');
    const sc_correction = 1.0 - (B / B_crit);

    const f_env = this.computeF_env(t);
    const msf_factor = this.computeMsfFactor(t);
    const m_factor = 1.0 + msf_factor;

    const G_const = this.variables.get('G');
    const M = this.variables.get('M');
    const r = this.variables.get('r');

    // Base gravity with all modulations
    let g_base = (G_const * M * m_factor / (r * r)) * expansion * sc_correction * (1.0 + f_env);

    // Ug sum
    const ug_sum = this.computeUgSum(r);

    // Cosmological constant term
    const c = this.variables.get('c');
    const Lambda = this.variables.get('Lambda');
    const lambda_term = Lambda * (c * c) / 3.0;

    // Quantum term
    const t_Hubble = this.variables.get('t_Hubble');
    const quantum_term = this.computeQuantumTerm(t_Hubble);

    // Fluid term
    const fluid_term = this.computeFluidTerm(g_base);

    // Dark matter perturbation
    const dm_pert_term = this.computeDMPertTerm(r);

    // Total unified field
    const g_total = g_base + ug_sum + lambda_term + quantum_term + fluid_term + dm_pert_term;

    return g_total;
  }

  /**
   * Get descriptive equation text
   * 
   * @returns {string} Full equation description
   */
  getEquationText() {
    return `g_UQFF(r, t) = (G * M(t) / r²) * (1 + H(z)*t) * (1 - B/B_crit) * (1 + F_env(t)) + (Ug1 + Ug2 + Ug3' + Ug4) + (Λ * c² / 3) + (ℏ / √(Δx·Δp)) * ∫ψ_total H ψ_total dV * (2π / t_Hubble) + ρ_fluid * V * g + (M_visible + M_DM) * (δρ/ρ + 3 G M / r³)

Where:
- H(t, z) = H₀ * √(Ω_m * (1+z)³ + Ω_Λ)
- M(t) = M * (1 + SFR * t_yr / M₀)
- F_env(t) = Σ F_i(t) [system-specific: wind, erosion, SN, mergers, radiation, filaments, BH]
- Ug3' = G * M_ext / r_ext²
- ψ_total = consolidated quantum waves integral
- System: ${this.currentSystem}
- Scale: Magnetars (10 km) to Cosmic fields (1 Gpc)`;
  }

  /**
   * Print all current variables (debugging)
   */
  printVariables() {
    console.log(`\n=== Current Variables for ${this.currentSystem} ===`);
    for (const [key, value] of this.variables) {
      console.log(`  ${key}: ${typeof value === 'number' ? value.toExponential(6) : value}`);
    }
  }

  /**
   * Get summary of system
   * 
   * @returns {object} System summary
   */
  getSummary() {
    return {
      system: this.currentSystem,
      M: this.variables.get('M'),
      r: this.variables.get('r'),
      z: this.variables.get('z'),
      SFR: this.variables.get('SFR'),
      v_wind: this.variables.get('v_wind'),
      M_ext: this.variables.get('M_ext'),
      r_ext: this.variables.get('r_ext'),
      g_default: this.computeG(),
      description: this.getEquationText()
    };
  }
}

module.exports = MultiUQFFCompressionModule;
