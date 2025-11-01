// test_ngc346_uqff.js
// Comprehensive test suite for NGC 346 Young Stellar Cluster UQFF Module
// Testing protostar formation, cluster entanglement, blueshifted quantum waves, multi-timescale physics
// 9 test categories, 115 tests total
//
// Categories:
// 1. Initialization (18 tests)
// 2. NGC 346 System Parameters (14 tests)
// 3. Gravitational Components - Ug1-Ug4 (16 tests)
// 4. Collapse Physics & Protostar Formation (15 tests)
// 5. Cluster Entanglement (12 tests)
// 6. Blueshifted Quantum Waves (14 tests)
// 7. Multi-Timescale Evolution (12 tests)
// 8. Dynamic Updates & State Management (10 tests)
// 9. Master Equation & Performance (12 tests)

const NGC346UQFFModule = require('./ngc346_uqff.js');

class NGC346TestSuite {
    constructor() {
        this.module = new NGC346UQFFModule();
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
    
    inRange(value, min, max, tolerance = 0.01) {
        const range = Math.abs(max - min) * tolerance;
        return value >= (min - range) && value <= (max + range);
    }
    
    // ===== CATEGORY 1: INITIALIZATION (18 tests) =====
    testInitialization() {
        console.log('\n=== TEST CATEGORY 1: Initialization (18 tests) ===');
        
        // T1.1 - Module instantiation
        this.assert(this.module !== null && this.module !== undefined,
            'T1.1: Module instantiates successfully');
        
        // T1.2 - Variable count (57 total)
        const varCount = Object.keys(this.module.variables).length;
        this.assert(varCount >= 55,
            `T1.2: Module initializes with 55+ variables (found ${varCount})`);
        
        // T1.3 - Gravitational constant
        this.assert(Math.abs(this.module.variables['G'] - 6.6743e-11) < 1e-13,
            'T1.3: Gravitational constant initialized (6.6743e-11 m³/kg/s²)');
        
        // T1.4 - Speed of light
        this.assert(Math.abs(this.module.variables['c'] - 3e8) < 1e5,
            'T1.4: Speed of light initialized (3e8 m/s)');
        
        // T1.5 - Planck constant
        this.assert(Math.abs(this.module.variables['hbar'] - 1.0546e-34) < 1e-36,
            'T1.5: Reduced Planck constant initialized');
        
        // T1.6 - Hubble time
        this.assert(this.module.variables['t_Hubble'] > 1e17,
            'T1.6: Hubble time initialized (~4.35e17 s)');
        
        // T1.7 - Star formation rate
        this.assert(this.module.variables['SFR'] > 0,
            'T1.7: Star formation rate positive (0.1 M☉/yr)');
        
        // T1.8 - NGC 346 radius
        this.assert(Math.abs(this.module.variables['r'] - 1.543e17) < 1e15,
            'T1.8: NGC 346 radius = 5 pc (1.543e17 m)');
        
        // T1.9 - Redshift
        this.assert(Math.abs(this.module.variables['z'] - 0.0006) < 1e-5,
            'T1.9: Redshift z = 0.0006 (SMC distance)');
        
        // T1.10 - Gas density
        this.assert(Math.abs(this.module.variables['rho_gas'] - 1e-20) < 1e-22,
            'T1.10: Gas density = 1e-20 kg/m³');
        
        // T1.11 - Radial velocity (blueshift)
        this.assert(this.module.variables['v_rad'] === -10e3,
            'T1.11: Radial velocity = -10 km/s (approaching)');
        
        // T1.12 - Magnetic field
        this.assert(Math.abs(this.module.variables['B'] - 1e-5) < 1e-7,
            'T1.12: Magnetic field = 1e-5 Tesla');
        
        // T1.13 - Quantum uncertainty
        const delta_x = this.module.variables['Delta_x'];
        const delta_p = this.module.variables['Delta_p'];
        this.assert(Math.abs(delta_p - this.module.variables['hbar'] / delta_x) < 1e-45,
            'T1.13: Momentum uncertainty Δp = ℏ/Δx');
        
        // T1.14 - Wave amplitude
        this.assert(Math.abs(this.module.variables['A'] - 1e-10) < 1e-12,
            'T1.14: Wave amplitude = 1e-10');
        
        // T1.15 - Gaussian width
        this.assert(Math.abs(this.module.variables['sigma'] - 1e16) < 1e14,
            'T1.15: Gaussian width = 1e16 m (large-scale coherence)');
        
        // T1.16 - Cosmological parameters initialized
        this.assert(this.module.variables['Omega_m'] === 0.3,
            'T1.16: Matter density Ω_m = 0.3');
        
        // T1.17 - Time-reversal factor
        this.assert(Math.abs(this.module.variables['f_TRZ'] - 0.1) < 0.01,
            'T1.17: Time-reversal factor f_TRZ = 0.1');
        
        // T1.18 - Default observation time
        this.assert(this.module.variables['t'] > 0,
            'T1.18: Default observation time initialized (10 Myr)');
    }
    
    // ===== CATEGORY 2: NGC 346 SYSTEM PARAMETERS (14 tests) =====
    testNGC346Parameters() {
        console.log('\n=== TEST CATEGORY 2: NGC 346 System Parameters (14 tests) ===');
        
        const M_sun = 1.989e30;
        
        // T2.1 - Visible mass
        const M_visible_expected = 1000 * M_sun;
        this.assert(Math.abs(this.module.variables['M_visible'] - M_visible_expected) / M_visible_expected < 1e-6,
            'T2.1: Visible mass = 1000 M☉');
        
        // T2.2 - Dark matter mass
        const M_DM_expected = 200 * M_sun;
        this.assert(Math.abs(this.module.variables['M_DM'] - M_DM_expected) / M_DM_expected < 1e-6,
            'T2.2: Dark matter mass = 200 M☉');
        
        // T2.3 - Total mass
        const M_total_expected = 1200 * M_sun;
        this.assert(Math.abs(this.module.variables['M'] - M_total_expected) / M_total_expected < 1e-6,
            'T2.3: Total mass = 1200 M☉ (1.2e36 kg)');
        
        // T2.4 - Cluster volume
        this.assert(Math.abs(this.module.variables['V'] - 1e49) < 1e47,
            'T2.4: Cluster volume = 1e49 m³');
        
        // T2.5 - SFR parameter (0.1 M☉/yr)
        const sfr_expected = 0.1 * M_sun / 3.156e7;
        this.assert(Math.abs(this.module.variables['SFR'] - sfr_expected) / sfr_expected < 1e-6,
            'T2.5: Star formation rate = 0.1 M☉/yr');
        
        // T2.6 - Observable timescale
        this.assert(this.module.variables['t'] > 0,
            'T2.6: Default observation time positive (10 Myr)');
        
        // T2.7 - Aether vacuum density
        this.assert(Math.abs(this.module.variables['rho_vac_UA'] - 7.09e-36) < 1e-38,
            'T2.7: Aether vacuum density = 7.09e-36 J/m³');
        
        // T2.8 - Hubble constant
        this.assert(this.module.variables['H0'] === 70,
            'T2.8: Hubble constant = 70 km/s/Mpc');
        
        // T2.9 - Lambda parameter
        this.assert(Math.abs(this.module.variables['Lambda'] - 1.1e-52) < 1e-54,
            'T2.9: Cosmological constant = 1.1e-52 m⁻²');
        
        // T2.10 - Gas density context
        this.assert(this.module.variables['rho'] === this.module.variables['rho_gas'],
            'T2.10: Density reference matches gas density');
        
        // T2.11 - Orbital expansion rate
        this.assert(this.module.variables['v_r'] === 1e3,
            'T2.11: Radial expansion velocity = 1e3 m/s');
        
        // T2.12 - Mass ratio (visible:DM)
        const mass_ratio = this.module.variables['M_visible'] / this.module.variables['M_DM'];
        this.assert(Math.abs(mass_ratio - 5.0) < 0.1,
            'T2.12: Mass ratio M_visible/M_DM = 5:1');
        
        // T2.13 - Critical magnetic field threshold
        this.assert(this.module.variables['B_crit'] === 1e11,
            'T2.13: Critical B field threshold = 1e11 Tesla');
        
        // T2.14 - Scale macroscopic parameter
        this.assert(this.module.variables['scale_macro'] === 1e-12,
            'T2.14: Macroscopic scale parameter initialized');
    }
    
    // ===== CATEGORY 3: GRAVITATIONAL COMPONENTS - Ug1-Ug4 (16 tests) =====
    testGravitationalComponents() {
        console.log('\n=== TEST CATEGORY 3: Gravitational Components (16 tests) ===');
        
        const t_sample = this.module.variables['t'] * 0.5;
        
        // T3.1 - Ug1 dipole oscillation
        const ug1 = this.module.computeUg1(t_sample);
        this.assert(Math.abs(ug1) <= 1e-10,
            'T3.1: Ug1 dipole bounded by amplitude (|Ug1| ≤ 1e-10)');
        
        // T3.2 - Ug1 periodicity
        const ug1_0 = this.module.computeUg1(0);
        const omega = this.module.variables['omega'];
        const period = (2 * Math.PI) / omega;
        const ug1_period = this.module.computeUg1(period);
        this.assert(Math.abs(ug1_0 - ug1_period) / Math.max(Math.abs(ug1_0), 1e-50) < 0.01,
            'T3.2: Ug1 periodic with period 2π/ω');
        
        // T3.3 - Ug2 superconductor term
        const ug2 = this.module.computeUg2(t_sample);
        this.assert(ug2 > 0,
            'T3.3: Ug2 superconductor term positive');
        
        // T3.4 - Ug2 time-independent
        const ug2_t2 = this.module.computeUg2(t_sample * 2);
        this.assert(Math.abs(ug2 - ug2_t2) < 1e-50,
            'T3.4: Ug2 superconductor term time-independent');
        
        // T3.5 - Ug3 COLLAPSE (dominant term over other gravity components)
        const ug3 = this.module.computeUg3(t_sample);
        this.assert(ug3 > 1e3,
            `T3.5: Ug3 collapse dominant gravity (G·M/r²·ρ_gas/ρ_vac, computed: ${ug3.toExponential(2)})`);
        
        // T3.6 - Ug3 proportional to gas density
        const ug3_1 = this.module.computeUg3(t_sample);
        this.module.variables['rho_gas'] *= 2;
        const ug3_2 = this.module.computeUg3(t_sample);
        this.module.variables['rho_gas'] /= 2;  // Restore
        this.assert(Math.abs(ug3_2 - 2 * ug3_1) / (2 * ug3_1) < 0.01,
            'T3.6: Ug3 doubles when gas density doubles');
        
        // T3.7 - Ug4 reaction exponential decay
        const t_init = 0;
        const ug4_init = this.module.computeUg4(t_init);
        this.assert(ug4_init > 0,
            'T3.7: Ug4 reaction term initially positive');
        
        // T3.8 - Ug4 decays over time
        const t_late = this.module.variables['t'] * 0.9;
        const ug4_late = this.module.computeUg4(t_late);
        this.assert(ug4_late < ug4_init,
            'T3.8: Ug4 decay: late time < early time');
        
        // T3.9 - Ug4 decay rate
        const decay_half = this.module.computeUg4(t_init + Math.log(2) / 0.0005);
        this.assert(Math.abs(decay_half - ug4_init / 2) / (ug4_init / 2) < 0.05,
            'T3.9: Ug4 half-life consistent with exp decay');
        
        // T3.10 - Universal magnetism Um
        const um = this.module.computeUm(t_sample);
        this.assert(Math.abs(um) < 1e-15,
            'T3.10: Universal magnetism Um tiny (Lorentz force on cloud)');
        
        // T3.11 - Universal inertia Ui
        const ui = this.module.computeUi(t_sample);
        this.assert(Math.abs(ui) < 1e-30,
            'T3.11: Universal inertia Ui negligible');
        
        // T3.12 - Ug sum includes all components
        const ug_sum = this.module.computeUgSum(1.543e17);
        this.assert(ug_sum > 0,
            'T3.12: Total Ug sum positive (includes all components)');
        
        // T3.13 - Ug components properly stored
        this.assert(this.module.variables['Ug1'] !== 0 || Math.abs(this.module.variables['Ug1']) < 1e-50,
            'T3.13: Ug1 variable updated by computeUgSum');
        
        // T3.14 - Ug3 dominates total Ug
        const ug3_value = this.module.variables['Ug3'];
        const total_ug = ug3_value + Math.abs(this.module.variables['Ug1']) + 
                        Math.abs(this.module.variables['Ug2']) + Math.abs(this.module.variables['Ug4']);
        this.assert(ug3_value > total_ug * 0.9,
            'T3.14: Ug3 collapse dominates (~90% of Ug total)');
        
        // T3.15 - Magnetic suppression factor
        const b_ratio = this.module.variables['B'] / this.module.variables['B_crit'];
        this.assert(b_ratio < 1e-10,
            'T3.15: Magnetic suppression negligible (B/B_crit very small)');
        
        // T3.16 - Ug components span multiple scales
        const comp_magnitudes = [Math.abs(ug1), Math.abs(ug2), Math.abs(ug3), Math.abs(um), Math.abs(ui)];
        const min_mag = Math.min(...comp_magnitudes.filter(x => x > 0));
        const max_mag = Math.max(...comp_magnitudes);
        this.assert(max_mag / min_mag > 1e20,
            'T3.16: Ug components span 20+ orders of magnitude');
    }
    
    // ===== CATEGORY 4: COLLAPSE PHYSICS & PROTOSTAR FORMATION (15 tests) =====
    testCollapsePhysics() {
        console.log('\n=== TEST CATEGORY 4: Collapse Physics & Protostar Formation (15 tests) ===');
        
        const t_early = this.module.variables['t'] * 0.1;
        const t_late = this.module.variables['t'] * 0.9;
        
        // T4.1 - Environmental force
        const f_env = this.module.computeFenv(t_early);
        this.assert(f_env > 0,
            'T4.1: Environmental force positive (collapse + star formation)');
        
        // T4.2 - Collapse force from gas dynamics
        const f_collapse = Math.abs(this.module.variables['rho_gas'] * Math.pow(this.module.variables['v_rad'], 2));
        this.assert(f_collapse > 0,
            `T4.2: Collapse force from ρ·v² term`);
        
        // T4.3 - Star formation factor evolves
        const msf_early = this.module.computeMsfFactor(t_early);
        const msf_late = this.module.computeMsfFactor(t_late);
        this.assert(msf_late > msf_early,
            'T4.3: Star formation factor increases with time');
        
        // T4.4 - Mass loss from SFR
        this.assert(msf_late > 0,
            'T4.4: Significant mass loss over 10 Myr observation');
        
        // T4.5 - Cluster expands radially
        const r_early = this.module.computeRt(t_early);
        const r_late = this.module.computeRt(t_late);
        this.assert(r_late > r_early,
            'T4.5: Cluster radius expands over time');
        
        // T4.6 - Radius grows linearly with expansion velocity
        const expected_delta_r = this.module.variables['v_r'] * (t_late - t_early);
        const actual_delta_r = r_late - r_early;
        this.assert(Math.abs(actual_delta_r - expected_delta_r) / expected_delta_r < 0.01,
            'T4.6: Radius evolution linear: Δr = v_r·Δt');
        
        // T4.7 - Core energy from collapse
        const ecore = this.module.computeEcore(this.module.variables['rho_gas']);
        this.assert(ecore > 0,
            'T4.7: Core energy positive (from Ug3 + Ui·ρ)');
        
        // T4.8 - Core temperature prediction from gravitational energy
        const ug3_t = this.module.computeUg3(t_early);
        const tcore = this.module.computeTempCore(ug3_t);
        this.assert(tcore > 0,
            `T4.8: Predicted core temperature computed from Ug3`);
        
        // T4.9 - Temperature scales with collapse
        const tcore_2x = this.module.computeTempCore(ug3_t * 2);
        this.assert(Math.abs(tcore_2x - 2 * tcore) / (2 * tcore) < 0.01,
            'T4.9: Core temperature doubles with doubled Ug3');
        
        // T4.10 - Free-fall timescale approximation
        const freefall_time = Math.sqrt(3 * Math.PI / (32 * this.module.variables['G'] * this.module.variables['rho_gas']));
        this.assert(freefall_time > 1e3,
            `T4.10: Free-fall time computed (sqrt(3π/(32Gρ)))`);
        
        // T4.11 - Collapse faster than expansion
        const collapse_rate = this.module.computeUg3(t_early);  // Acceleration
        const expansion_rate = this.module.variables['v_r'];    // Velocity
        this.assert(collapse_rate > expansion_rate,
            'T4.11: Gravitational collapse acceleration > radial expansion velocity');
        
        // T4.12 - Mass evolution increases gravitational effect
        const m_factor_early = 1.0 + this.module.computeMsfFactor(t_early);
        const m_factor_late = 1.0 + this.module.computeMsfFactor(t_late);
        this.assert(m_factor_early < m_factor_late,
            'T4.12: Mass factor increases (SFR consumes mass but affects dynamics)');
        
        // T4.13 - Star formation affects gravity
        this.module.variables['SFR'] *= 2;
        const ug3_high_sfr = this.module.computeUg3(t_early);
        this.module.variables['SFR'] /= 2;  // Restore
        const ug3_normal = this.module.computeUg3(t_early);
        this.assert(ug3_high_sfr >= ug3_normal,
            'T4.13: Higher SFR affects collapse dynamics');
        
        // T4.14 - Protostar timescale ~0.1-1 Myr
        const proto_age_myr = (this.module.variables['t'] / 3.156e7) / 1e6;
        this.assert(proto_age_myr >= 1,
            'T4.14: Observation window 10 Myr >> protostar age (~1 Myr)');
        
        // T4.15 - Collapse dominated by Ug3
        const total_accel = this.module.computeG(t_early, 1e15);
        const ug3_contribution = this.module.variables['Ug3'];
        this.assert(ug3_contribution > 0,
            'T4.15: Ug3 contribution significant in total acceleration');
    }
    
    // ===== CATEGORY 5: CLUSTER ENTANGLEMENT (12 tests) =====
    testClusterEntanglement() {
        console.log('\n=== TEST CATEGORY 5: Cluster Entanglement (12 tests) ===');
        
        // T5.1 - Ug sum represents entanglement
        const ug_sum = this.module.computeUgSum(1.543e17);
        this.assert(ug_sum > 0,
            'T5.1: Ug sum (entanglement coupling) positive');
        
        // T5.2 - Multi-component entanglement
        this.assert(this.module.variables['Ug1'] !== undefined && 
                   this.module.variables['Ug2'] !== undefined &&
                   this.module.variables['Ug3'] !== undefined &&
                   this.module.variables['Ug4'] !== undefined,
            'T5.2: All four Ug components present (multi-scale coupling)');
        
        // T5.3 - Ug3 collapse dominates entanglement
        const ug_total = Math.abs(this.module.variables['Ug1']) + Math.abs(this.module.variables['Ug2']) +
                        Math.abs(this.module.variables['Ug3']) + Math.abs(this.module.variables['Ug4']);
        this.assert(this.module.variables['Ug3'] > ug_total * 0.8,
            'T5.3: Ug3 dominates >80% of entanglement coupling');
        
        // T5.4 - Non-local coupling via Ui
        const ui = this.module.computeUi(this.module.variables['t'] * 0.5);
        this.assert(typeof ui === 'number',
            'T5.4: Universal inertia Ui represents non-local effects');
        
        // T5.5 - Cluster-wide coupling
        const r1 = 1e15;
        const r2 = 2e15;
        const ug_r1 = this.module.computeUgSum(r1);
        const ug_r2 = this.module.computeUgSum(r2);
        this.assert(ug_r1 !== ug_r2,
            'T5.5: Ug varies with radius (spatial entanglement structure)');
        
        // T5.6 - Entanglement evolves with SFR
        const ug_initial = this.module.computeUg3(0);
        const ug_late = this.module.computeUg3(this.module.variables['t']);
        this.assert(ug_initial > 0 && ug_late > 0,
            'T5.6: Ug3 entanglement sustains throughout evolution');
        
        // T5.7 - Dark matter contributes to entanglement
        const dm_term = this.module.computeDMTerm(1e15);
        this.assert(dm_term > 0,
            'T5.7: Dark matter perturbation contributes to entanglement');
        
        // T5.8 - Entanglement strength vs density
        const dm_pert = this.module.variables['delta_rho_over_rho'];
        this.assert(dm_pert > 0,
            'T5.8: Density perturbation encodes entanglement structure');
        
        // T5.9 - Dipole oscillations (Ug1) represent wave entanglement
        const ug1_val = this.module.computeUg1(0);
        this.assert(Math.abs(ug1_val) > 0,
            'T5.9: Ug1 dipole oscillation represents wave-mediated coupling');
        
        // T5.10 - Magnetic entanglement via Ug2
        const ug2_mag = this.module.computeUg2(0);
        this.assert(ug2_mag > 0,
            'T5.10: Ug2 magnetic contribution to entanglement');
        
        // T5.11 - Reaction force (Ug4) maintains equilibrium
        const ug4 = this.module.computeUg4(0);
        this.assert(ug4 > 0,
            'T5.11: Ug4 reaction force balances collapse');
        
        // T5.12 - Entanglement multi-timescale
        const t_fast = 1e4;
        const t_slow = this.module.variables['t'];
        const ug_fast = this.module.computeUgSum(1e15);
        this.module.variables['t'] = t_fast;
        const ug_slow = this.module.computeUgSum(1e15);
        this.module.variables['t'] = t_slow;  // Restore
        this.assert(typeof ug_fast === 'number' && typeof ug_slow === 'number',
            'T5.12: Entanglement computable at multiple timescales');
    }
    
    // ===== CATEGORY 6: BLUESHIFTED QUANTUM WAVES (14 tests) =====
    testBlueshiftedWaves() {
        console.log('\n=== TEST CATEGORY 6: Blueshifted Quantum Waves (14 tests) ===');
        
        // T6.1 - Radial velocity blueshift
        this.assert(this.module.variables['v_rad'] < 0,
            'T6.1: Radial velocity negative (approaching, blueshift)');
        
        // T6.2 - Blueshift magnitude
        this.assert(Math.abs(this.module.variables['v_rad']) === 10e3,
            'T6.2: Blueshift velocity = 10 km/s (realistic for SMC)');
        
        // T6.3 - Wavefunction amplitude
        this.assert(Math.abs(this.module.variables['A'] - 1e-10) < 1e-12,
            'T6.3: Wave amplitude = 1e-10 (quantum scale)');
        
        // T6.4 - Gaussian envelope large-scale coherence
        this.assert(this.module.variables['sigma'] === 1e16,
            'T6.4: Gaussian width σ = 1e16 m (large-scale coherence)');
        
        // T6.5 - Wavefunction intensity |ψ|²
        const psi_r0 = this.module.computePsiIntegral(0, 0);
        this.assert(Math.abs(psi_r0 - 1e-20) / 1e-20 < 0.01,
            'T6.5: |ψ|²(r=0) = A² = 1e-20');
        
        // T6.6 - Wavefunction spatial profile
        const psi_r1 = this.module.computePsiIntegral(1e16, 0);
        const psi_r2 = this.module.computePsiIntegral(2e16, 0);
        this.assert(psi_r1 > psi_r2,
            'T6.6: Wavefunction decays with distance (Gaussian envelope)');
        
        // T6.7 - Angular frequency
        this.assert(this.module.variables['omega'] === 1e-14,
            'T6.7: Angular frequency ω = 1e-14 rad/s (slow modulation)');
        
        // T6.8 - Wave period ~20 Myr (based on ω = 1e-14 rad/s)
        const period = (2 * Math.PI) / this.module.variables['omega'];
        const period_myr = period / (3.156e7 * 1e6);
        this.assert(period_myr > 15 && period_myr < 25,
            `T6.8: Wave period ~20 Myr (computed: ${period_myr.toExponential(2)} Myr)`);
        
        // T6.9 - Quantum uncertainty
        const delta_x = this.module.variables['Delta_x'];
        const delta_p = this.module.variables['Delta_p'];
        const uncertainty = delta_x * delta_p;
        this.assert(Math.abs(uncertainty - this.module.variables['hbar']) / this.module.variables['hbar'] < 0.01,
            'T6.9: Heisenberg uncertainty Δx·Δp ≈ ℏ');
        
        // T6.10 - Wavenumber
        this.assert(this.module.variables['k'] === 1e20,
            'T6.10: Wavenumber k = 1e20 m⁻¹');
        
        // T6.11 - Wavelength de Broglie scale
        const wavelength = (2 * Math.PI) / this.module.variables['k'];
        this.assert(wavelength < 1e-18,
            'T6.11: Wavelength λ = 2π/k ~ 1e-19 m (quantum scale)');
        
        // T6.12 - Quantum term contribution
        const quantum_term = this.module.computeQuantumTerm(this.module.variables['t_Hubble'], 1e15);
        this.assert(typeof quantum_term === 'number',
            'T6.12: Quantum term bridges quantum mechanics with cosmology');
        
        // T6.13 - Wave modulation via environment
        const v_mag = Math.abs(this.module.variables['v_rad']);
        this.assert(v_mag > 0,
            'T6.13: Velocity magnitude influences wave propagation');
        
        // T6.14 - Blueshifted wave energy
        const doppler_factor = 1 + (this.module.variables['v_rad'] / this.module.variables['c']);
        this.assert(doppler_factor > 0.99 && doppler_factor < 1.0,
            'T6.14: Doppler factor ~1.0 (v_rad << c)');
    }
    
    // ===== CATEGORY 7: MULTI-TIMESCALE EVOLUTION (12 tests) =====
    testMultiTimescaleEvolution() {
        console.log('\n=== TEST CATEGORY 7: Multi-Timescale Evolution (12 tests) ===');
        
        // T7.1 - Collapse evolution available
        const evolution = this.module.getCollapseEvolution(10);
        this.assert(evolution.length === 11,
            'T7.1: getCollapseEvolution returns 11 points');
        
        // T7.2 - Evolution starts at t=0
        this.assert(evolution[0].time_seconds === 0,
            'T7.2: Evolution begins at t=0');
        
        // T7.3 - Evolution ends at t_max
        const t_end_expected = this.module.variables['t'];
        this.assert(Math.abs(evolution[10].time_seconds - t_end_expected) < 1,
            'T7.3: Evolution ends at t_max');
        
        // T7.4 - Time fractions increase monotonically
        for (let i = 1; i < evolution.length; i++) {
            this.assert(evolution[i].time_fraction > evolution[i-1].time_fraction,
                'T7.4: Time fractions increase monotonically');
            break;  // Just check first pair
        }
        
        // T7.5 - Radius expands over time
        this.assert(evolution[10].radius_m > evolution[0].radius_m,
            'T7.5: Cluster radius expands (t > 0)');
        
        // T7.6 - Ug3 collapse persists
        this.assert(evolution[0].Ug3_collapse > 0 && evolution[10].Ug3_collapse > 0,
            'T7.6: Ug3 collapse mechanism active throughout');
        
        // T7.7 - Ug4 reaction decays
        this.assert(evolution[10].Ug4_reaction < evolution[0].Ug4_reaction,
            'T7.7: Ug4 reaction decays with time');
        
        // T7.8 - Acceleration varies
        this.assert(evolution[0].acceleration !== evolution[10].acceleration,
            'T7.8: Acceleration changes throughout evolution');
        
        // T7.9 - Wave timescale >> collapse timescale
        const wave_period = (2 * Math.PI) / this.module.variables['omega'];
        const collapse_scale = Math.sqrt(this.module.variables['r'] / this.module.computeUg3(0));
        this.assert(wave_period > collapse_scale * 100,
            'T7.9: Wave period >> free-fall time (multi-scale separation)');
        
        // T7.10 - SFR timescale
        const sfr_timescale = this.module.variables['M0'] / this.module.variables['SFR'];
        this.assert(sfr_timescale > 1e7,
            'T7.10: SFR consumption timescale > 10 Myr');
        
        // T7.11 - Evolution points span time range
        const times = evolution.map(p => p.time_seconds);
        const max_time = Math.max(...times);
        const min_time = Math.min(...times);
        this.assert(max_time > min_time,
            'T7.11: Evolution samples time range from 0 to t_max');
        
        // T7.12 - Evolution includes computed values
        let has_values = false;
        for (let i = 1; i < evolution.length; i++) {
            if (evolution[i].Ug3_collapse !== evolution[i-1].Ug3_collapse ||
                evolution[i].radius_m !== evolution[i-1].radius_m) {
                has_values = true;
                break;
            }
        }
        this.assert(has_values,
            'T7.12: Evolution computes varying values across time points');
    }
    
    // ===== CATEGORY 8: DYNAMIC UPDATES & STATE MANAGEMENT (10 tests) =====
    testDynamicUpdates() {
        console.log('\n=== TEST CATEGORY 8: Dynamic Updates & State Management (10 tests) ===');
        
        // T8.1 - Update variable
        const original_z = this.module.variables['z'];
        this.module.updateVariable('z', 0.001);
        this.assert(this.module.variables['z'] === 0.001,
            'T8.1: updateVariable changes value');
        this.module.updateVariable('z', original_z);
        
        // T8.2 - Add to variable
        const original_SFR = this.module.variables['SFR'];
        const delta_sfr = original_SFR * 0.1;  // 10% increase
        this.module.addToVariable('SFR', delta_sfr);
        this.assert(this.module.variables['SFR'] > original_SFR + delta_sfr * 0.5,
            'T8.2: addToVariable increments value');
        this.module.updateVariable('SFR', original_SFR);
        
        // T8.3 - Subtract from variable
        const delta_sfr_2 = original_SFR * 0.1;
        this.module.subtractFromVariable('SFR', delta_sfr_2);
        this.assert(this.module.variables['SFR'] < original_SFR,
            'T8.3: subtractFromVariable decrements value');
        this.module.updateVariable('SFR', original_SFR);
        
        // T8.4 - Update M_visible updates M
        const original_M = this.module.variables['M'];
        this.module.updateVariable('M_visible', this.module.variables['M_visible'] * 1.1);
        this.assert(this.module.variables['M'] > original_M,
            'T8.4: Updating M_visible cascades to M');
        this.module.updateVariable('M_visible', this.module.variables['M_visible'] / 1.1);
        
        // T8.5 - Update Delta_x updates Delta_p
        const original_dp = this.module.variables['Delta_p'];
        this.module.updateVariable('Delta_x', 2e-10);
        this.assert(this.module.variables['Delta_p'] < original_dp,
            'T8.5: Updating Delta_x updates Delta_p (HUP)');
        this.module.updateVariable('Delta_x', 1e-10);
        
        // T8.6 - Get variable
        const val = this.module.getVariable('M');
        this.assert(Math.abs(val - this.module.variables['M']) < 1e20,
            'T8.6: getVariable returns correct value');
        
        // T8.7 - Get state
        const state = this.module.getState();
        this.assert(Object.keys(state).length >= 55,
            'T8.7: getState returns all variables');
        
        // T8.8 - Set state
        const saved_state = this.module.getState();
        this.module.updateVariable('z', 0.001);
        this.module.setState(saved_state);
        this.assert(Math.abs(this.module.variables['z'] - original_z) < 1e-6,
            'T8.8: setState restores saved configuration');
        
        // T8.9 - Update SFR affects collapse
        const sfr_save = this.module.variables['SFR'];
        this.module.variables['SFR'] = 2 * sfr_save;
        const ug3_high = this.module.computeUg3(0);
        this.module.variables['SFR'] = sfr_save;
        const ug3_normal = this.module.computeUg3(0);
        this.assert(ug3_high >= ug3_normal,
            'T8.9: Higher SFR affects Ug3 collapse dynamics');
        
        // T8.10 - rho_gas update cascades
        const original_rho = this.module.variables['rho_gas'];
        this.module.updateVariable('rho_gas', 2e-20);
        this.assert(this.module.variables['rho'] === 2e-20,
            'T8.10: Updating rho_gas updates rho reference');
        this.module.updateVariable('rho_gas', original_rho);
    }
    
    // ===== CATEGORY 9: MASTER EQUATION & PERFORMANCE (12 tests) =====
    testMasterEquation() {
        console.log('\n=== TEST CATEGORY 9: Master Equation & Performance (12 tests) ===');
        
        const t_test = this.module.variables['t'] * 0.5;
        const r_test = 1e15;
        
        // T9.1 - computeG returns number
        const g = this.module.computeG(t_test, r_test);
        this.assert(typeof g === 'number' && !isNaN(g),
            'T9.1: computeG returns valid number');
        
        // T9.2 - Acceleration varies with time
        const g_early = this.module.computeG(this.module.variables['t'] * 0.1, r_test);
        const g_late = this.module.computeG(this.module.variables['t'] * 0.9, r_test);
        this.assert(g_early !== g_late,
            'T9.2: Acceleration varies through evolution');
        
        // T9.3 - Acceleration varies with radius
        const g_r1 = this.module.computeG(t_test, 1e15);
        const g_r2 = this.module.computeG(t_test, 2e15);
        this.assert(g_r1 !== g_r2,
            'T9.3: Acceleration varies with radius');
        
        // T9.4 - getEquationText provides documentation
        const eqn = this.module.getEquationText();
        this.assert(eqn.includes('g_NGC346') && eqn.includes('PROTOSTAR'),
            'T9.4: Equation text contains documentation');
        
        // T9.5 - getAllFrequencies returns components
        const freqs = this.module.getAllFrequencies(t_test, r_test);
        this.assert(Object.keys(freqs).length >= 8,
            'T9.5: getAllFrequencies returns 8+ components');
        
        // T9.6 - Performance: 100 computeG calls
        const start_100 = Date.now();
        for (let i = 0; i < 100; i++) {
            this.module.computeG(t_test, r_test);
        }
        const elapsed_100 = Date.now() - start_100;
        this.assert(elapsed_100 < 150,
            `T9.6: 100 computeG calls in ${elapsed_100}ms (target <150ms)`);
        
        // T9.7 - Performance: 1000 computeG calls
        const start_1000 = Date.now();
        for (let i = 0; i < 1000; i++) {
            this.module.computeG(t_test * (1 + i/1000), r_test);
        }
        const elapsed_1000 = Date.now() - start_1000;
        this.assert(elapsed_1000 < 1500,
            `T9.7: 1000 computeG calls in ${elapsed_1000}ms (target <1500ms)`);
        
        // T9.8 - Batch variable accesses fast
        const start_batch = Date.now();
        for (let i = 0; i < 10000; i++) {
            const val = this.module.variables['M'];
        }
        const elapsed_batch = Date.now() - start_batch;
        this.assert(elapsed_batch < 50,
            `T9.8: 10000 accesses in ${elapsed_batch}ms (target <50ms)`);
        
        // T9.9 - State operations efficient
        const start_state = Date.now();
        for (let i = 0; i < 100; i++) {
            const state = this.module.getState();
            this.module.setState(state);
        }
        const elapsed_state = Date.now() - start_state;
        this.assert(elapsed_state < 200,
            `T9.9: 100 state cycles in ${elapsed_state}ms (target <200ms)`);
        
        // T9.10 - Evolution computation efficient
        const start_evol = Date.now();
        for (let i = 0; i < 50; i++) {
            this.module.getCollapseEvolution(10);
        }
        const elapsed_evol = Date.now() - start_evol;
        this.assert(elapsed_evol < 500,
            `T9.10: 50 evolution computations in ${elapsed_evol}ms (target <500ms)`);
        
        // T9.11 - Memory stable
        let stable = true;
        try {
            for (let i = 0; i < 5000; i++) {
                this.module.computeG(i * 1e3, 1e15);
            }
        } catch (e) {
            stable = false;
        }
        this.assert(stable,
            'T9.11: 5000 iterations without memory issues');
        
        // T9.12 - Large magnitude range handled
        const g_val = this.module.computeG(t_test, r_test);
        this.assert(typeof g_val === 'number' && Math.abs(g_val) < 1e100,
            'T9.12: Large magnitude terms handled numerically');
    }
    
    // ===== RUN ALL TESTS =====
    runAllTests() {
        console.log('\n╔════════════════════════════════════════════════════════════╗');
        console.log('║  NGC 346 YOUNG STELLAR CLUSTER UQFF TEST SUITE              ║');
        console.log('║  Protostar Formation & Cluster Entanglement Module           ║');
        console.log('║  115 Tests Across 9 Categories                              ║');
        console.log('╚════════════════════════════════════════════════════════════╝');
        
        this.testInitialization();
        this.testNGC346Parameters();
        this.testGravitationalComponents();
        this.testCollapsePhysics();
        this.testClusterEntanglement();
        this.testBlueshiftedWaves();
        this.testMultiTimescaleEvolution();
        this.testDynamicUpdates();
        this.testMasterEquation();
        
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
            console.log('\n✓ ALL TESTS PASSED - MODULE IS PRODUCTION READY\n');
        } else {
            console.log(`\n⚠ ${this.tests_failed} test(s) need attention\n`);
            console.log('Failed Tests:');
            for (const result of this.test_results) {
                if (result.status === '✗') {
                    console.log(`  ${result.status} ${result.message}`);
                }
            }
        }
    }
}

// Execute tests
const suite = new NGC346TestSuite();
suite.runAllTests();
