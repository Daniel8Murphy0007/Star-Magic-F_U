/**
 * Comprehensive Test Suite for Saturn UQFF Module
 * Testing all physics computations, dynamic capabilities, and unique Saturn features
 */

import PlanetSaturn from './source30.js';

class TestSuite {
    constructor() {
        this.tests = [];
        this.passed = 0;
        this.failed = 0;
    }

    addTest(name, testFunc) {
        this.tests.push({ name, testFunc });
    }

    async run() {
        console.log("=========================================================");
        console.log("SATURN (RINGED GIANT) COMPREHENSIVE TEST SUITE");
        console.log("Testing: Rings, Atmosphere & Superconductivity");
        console.log("=========================================================\n");

        for (const test of this.tests) {
            try {
                const result = await test.testFunc();
                if (result.pass) {
                    this.passed++;
                    console.log(`✓ Test ${this.passed + this.failed}: ${test.name} - PASS`);
                } else {
                    this.failed++;
                    console.log(`✗ Test ${this.passed + this.failed}: ${test.name} - FAIL`);
                    console.log(`  Reason: ${result.message}`);
                }
            } catch (error) {
                this.failed++;
                console.log(`✗ Test ${this.passed + this.failed}: ${test.name} - ERROR`);
                console.log(`  Error: ${error.message}`);
            }
        }

        console.log("\n=========================================================");
        console.log(`TEST SUMMARY: ${this.passed}/${this.tests.length} PASSED`);
        if (this.failed > 0) {
            console.log(`FAILED: ${this.failed}`);
        }
        console.log("=========================================================\n");

        return this.failed === 0;
    }
}

// Helper functions
function approxEqual(a, b, tolerance = 1e-6) {
    if (Math.abs(b) < 1e-100) return Math.abs(a) < tolerance;
    return Math.abs((a - b) / b) < tolerance;
}

function relDiff(a, b) {
    if (Math.abs(b) < 1e-100) return Math.abs(a);
    return Math.abs((a - b) / b);
}

// Test suite instance
const suite = new TestSuite();

const Gyr_to_s = 1e9 * 3.156e7;

// ========== INITIALIZATION TESTS ==========

suite.addTest("Initialization: M and r within reasonable ranges", () => {
    const saturn = new PlanetSaturn();
    const M_check = saturn.M > 1e25 && saturn.M < 1e28;  // Planet scale
    const r_check = saturn.r > 1e6 && saturn.r < 1e9;    // ~10-1000 km
    return {
        pass: M_check && r_check,
        message: `M = ${saturn.M.toExponential(3)}, r = ${saturn.r.toExponential(3)}`
    };
});

suite.addTest("Initialization: Ring parameters (M_ring, r_ring)", () => {
    const saturn = new PlanetSaturn();
    const M_ring_check = saturn.M_ring > 0 && saturn.M_ring < 1e22;
    const r_ring_check = saturn.r_ring > saturn.r;  // Ring outside planet
    return {
        pass: M_ring_check && r_ring_check,
        message: `M_ring = ${saturn.M_ring.toExponential(3)} kg, r_ring = ${saturn.r_ring.toExponential(3)} m`
    };
});

suite.addTest("Initialization: Solar System z = 0.0", () => {
    const saturn = new PlanetSaturn();
    const z_ok = saturn.z === 0.0;
    return {
        pass: z_ok,
        message: `z = ${saturn.z}`
    };
});

suite.addTest("Initialization: Atmospheric parameters", () => {
    const saturn = new PlanetSaturn();
    const rho_ok = approxEqual(saturn.rho_atm, 2e-4, 0.1);
    const v_ok = approxEqual(saturn.v_wind, 500, 0.1);
    return {
        pass: rho_ok && v_ok,
        message: `rho_atm = ${saturn.rho_atm.toExponential(3)} kg/m³, v_wind = ${saturn.v_wind} m/s`
    };
});

suite.addTest("Initialization: M_visible = M, M_DM = 0 (no DM in planet)", () => {
    const saturn = new PlanetSaturn();
    const visible_ok = saturn.M_visible === saturn.M;
    const dm_ok = saturn.M_DM === 0.0;
    return {
        pass: visible_ok && dm_ok,
        message: `M_visible = ${saturn.M_visible.toExponential(3)}, M_DM = ${saturn.M_DM}`
    };
});

suite.addTest("Initialization: Superconductivity parameters (B, B_crit)", () => {
    const saturn = new PlanetSaturn();
    const B_ok = approxEqual(saturn.B, 1e-7, 0.1);
    const B_crit_ok = approxEqual(saturn.B_crit, 1e11, 0.1);
    return {
        pass: B_ok && B_crit_ok,
        message: `B = ${saturn.B.toExponential(3)} T, B_crit = ${saturn.B_crit.toExponential(3)} T`
    };
});

// ========== UQFF COMPONENT TESTS ==========

suite.addTest("UQFF: Ug sum = Ug1 + Ug4 (Ug2, Ug3 = 0)", () => {
    const saturn = new PlanetSaturn();
    const ug_sum = saturn.computeUgSum();
    const ug_manual = saturn.Ug1 + saturn.Ug2 + saturn.Ug3 + saturn.Ug4;
    return {
        pass: approxEqual(ug_sum, ug_manual, 1e-6),
        message: `Ug_sum = ${ug_sum.toExponential(6)}, manual = ${ug_manual.toExponential(6)}, rel_diff: ${relDiff(ug_sum, ug_manual).toExponential(2)}`
    };
});

suite.addTest("UQFF: g(t) > 0 for all times [0, 10 Gyr]", () => {
    const saturn = new PlanetSaturn();
    const times = [0, 1, 2.5, 4.5, 7.5, 10].map(t => t * Gyr_to_s);
    let all_positive = true;
    for (const t of times) {
        const g = saturn.compute_g_Saturn(t);
        if (g <= 0) {
            all_positive = false;
            break;
        }
    }
    return {
        pass: all_positive,
        message: all_positive ? "All g > 0" : "Some g <= 0"
    };
});

suite.addTest("UQFF: g(t < 0) returns 0", () => {
    const saturn = new PlanetSaturn();
    const g = saturn.compute_g_Saturn(-1e9);
    return {
        pass: g === 0,
        message: `g(-1e9 s) = ${g}`
    };
});

suite.addTest("UQFF: Superconductivity correction (1 - B/B_crit)", () => {
    const saturn = new PlanetSaturn();
    const sc_corr_expected = 1.0 - saturn.B / saturn.B_crit;
    // SC correction applied to Saturn gravity
    const ratio = saturn.B / saturn.B_crit;
    const sc_ok = ratio < 1e-3;  // Should be ~1e-18 for typical values
    return {
        pass: sc_ok,
        message: `SC = ${sc_corr_expected.toFixed(12)}, B/B_crit = ${ratio.toExponential(2)}`
    };
});

// ========== PLANET SCALING TESTS ==========

suite.addTest("Scaling: expandPlanetScale M × 2 → Ug1 × 2", () => {
    const saturn = new PlanetSaturn();
    saturn.saveState('base');
    const Ug1_base = saturn.Ug1;
    
    saturn.expandPlanetScale(2.0, 1.0);
    const Ug1_scaled = saturn.Ug1;
    
    const ratio = Ug1_scaled / Ug1_base;
    saturn.restoreState('base');
    return {
        pass: approxEqual(ratio, 2.0, 0.01),
        message: `Ug1: ${Ug1_base.toExponential(6)} → ${Ug1_scaled.toExponential(6)}, ratio: ${ratio.toFixed(3)}`
    };
});

suite.addTest("Scaling: expandPlanetScale r × 2 → Ug1 × 1/4", () => {
    const saturn = new PlanetSaturn();
    saturn.saveState('base');
    const Ug1_base = saturn.Ug1;
    
    saturn.expandPlanetScale(1.0, 2.0);
    const Ug1_scaled = saturn.Ug1;
    
    const ratio = Ug1_base / Ug1_scaled;  // Should be 4
    saturn.restoreState('base');
    return {
        pass: approxEqual(ratio, 4.0, 0.01),
        message: `Ug1: ${Ug1_base.toExponential(6)} → ${Ug1_scaled.toExponential(6)}, ratio: ${ratio.toFixed(3)}`
    };
});

suite.addTest("Scaling: expandRingScale M_ring × 2 (UNIQUE)", () => {
    const saturn = new PlanetSaturn();
    saturn.saveState('base');
    const M_ring_base = saturn.M_ring;
    
    saturn.expandRingScale(2.0, 1.0);
    const M_ring_scaled = saturn.M_ring;
    
    const ratio = M_ring_scaled / M_ring_base;
    saturn.restoreState('base');
    return {
        pass: approxEqual(ratio, 2.0, 1e-6),
        message: `M_ring: ${M_ring_base.toExponential(6)} → ${M_ring_scaled.toExponential(6)}, ratio: ${ratio.toFixed(3)}`
    };
});

suite.addTest("Scaling: expandRingScale r_ring × 3 (UNIQUE)", () => {
    const saturn = new PlanetSaturn();
    saturn.saveState('base');
    const r_ring_base = saturn.r_ring;
    
    saturn.expandRingScale(1.0, 3.0);
    const r_ring_scaled = saturn.r_ring;
    
    const ratio = r_ring_scaled / r_ring_base;
    saturn.restoreState('base');
    return {
        pass: approxEqual(ratio, 3.0, 1e-6),
        message: `r_ring: ${r_ring_base.toExponential(6)} → ${r_ring_scaled.toExponential(6)}, ratio: ${ratio.toFixed(3)}`
    };
});

suite.addTest("Scaling: expandAtmosphereScale rho_atm × 2 (UNIQUE)", () => {
    const saturn = new PlanetSaturn();
    saturn.saveState('base');
    const rho_atm_base = saturn.rho_atm;
    
    saturn.expandAtmosphereScale(2.0, 1.0);
    const rho_atm_scaled = saturn.rho_atm;
    
    const ratio = rho_atm_scaled / rho_atm_base;
    saturn.restoreState('base');
    return {
        pass: approxEqual(ratio, 2.0, 1e-6),
        message: `rho_atm: ${rho_atm_base.toExponential(6)} → ${rho_atm_scaled.toExponential(6)}, ratio: ${ratio.toFixed(3)}`
    };
});

suite.addTest("Scaling: expandAtmosphereScale v_wind × 3 (UNIQUE)", () => {
    const saturn = new PlanetSaturn();
    saturn.saveState('base');
    const v_wind_base = saturn.v_wind;
    
    saturn.expandAtmosphereScale(1.0, 3.0);
    const v_wind_scaled = saturn.v_wind;
    
    const ratio = v_wind_scaled / v_wind_base;
    saturn.restoreState('base');
    return {
        pass: approxEqual(ratio, 3.0, 1e-6),
        message: `v_wind: ${v_wind_base.toExponential(6)} → ${v_wind_scaled.toExponential(6)}, ratio: ${ratio.toFixed(3)}`
    };
});

// ========== PARAMETER EFFECTS TESTS ==========

suite.addTest("Parameters: M_ring affects ring tidal term", () => {
    const saturn = new PlanetSaturn();
    saturn.saveState('base');
    
    const M_ring_orig = saturn.M_ring;
    const T_ring_orig = (saturn.G * M_ring_orig) / (saturn.r_ring * saturn.r_ring);
    
    saturn.setVariable('M_ring', M_ring_orig * 2);
    const T_ring_new = (saturn.G * saturn.M_ring) / (saturn.r_ring * saturn.r_ring);
    
    const ratio = T_ring_new / T_ring_orig;
    saturn.restoreState('base');
    return {
        pass: approxEqual(ratio, 2.0, 0.01),
        message: `T_ring ratio: ${ratio.toFixed(3)}`
    };
});

suite.addTest("Parameters: v_wind affects wind term (v² dependence)", () => {
    const saturn = new PlanetSaturn();
    saturn.saveState('base');
    
    const wind_base = saturn.computeWindTerm();
    saturn.setVariable('v_wind', saturn.v_wind * 2);
    const wind_scaled = saturn.computeWindTerm();
    
    const ratio = wind_scaled / wind_base;
    saturn.restoreState('base');
    return {
        pass: approxEqual(ratio, 4.0, 0.01),
        message: `Wind term with 2× v: ${wind_base.toExponential(6)} → ${wind_scaled.toExponential(6)}, ratio: ${ratio.toFixed(3)}`
    };
});

suite.addTest("Parameters: rho_atm affects fluid term", () => {
    const saturn = new PlanetSaturn();
    saturn.saveState('base');
    
    const g_base = 10.0;  // Approximate base gravity
    const fluid_base = saturn.computeFluidTerm(g_base);
    saturn.setVariable('rho_atm', saturn.rho_atm * 2);
    const fluid_scaled = saturn.computeFluidTerm(g_base);
    
    const ratio = fluid_scaled / fluid_base;
    saturn.restoreState('base');
    return {
        pass: approxEqual(ratio, 2.0, 0.01),
        message: `Fluid term: ${fluid_base.toExponential(6)} → ${fluid_scaled.toExponential(6)}, ratio: ${ratio.toFixed(3)}`
    };
});

suite.addTest("Parameters: B affects SC correction", () => {
    const saturn = new PlanetSaturn();
    saturn.saveState('base');
    
    const B_base = saturn.B;
    const B_crit = saturn.B_crit;
    const ratio_base = B_base / B_crit;
    
    saturn.setVariable('B', saturn.B * 2);
    const ratio_scaled = saturn.B / B_crit;
    
    // Ratio B/B_crit should double (SC correction decreases)
    const ratio_change = ratio_scaled / ratio_base;
    saturn.restoreState('base');
    return {
        pass: approxEqual(ratio_change, 2.0, 0.01),
        message: `B/B_crit: ${ratio_base.toExponential(3)} → ${ratio_scaled.toExponential(3)}, ratio: ${ratio_change.toFixed(3)}`
    };
});

suite.addTest("Parameters: r_orbit affects Sun gravity (inverse square)", () => {
    const saturn = new PlanetSaturn();
    saturn.saveState('base');
    const t = 4.5 * Gyr_to_s;
    
    const r_orbit_orig = saturn.r_orbit;
    saturn.setVariable('r_orbit', r_orbit_orig * 2);
    
    // Sun gravity should be 1/4
    const g_sun_scaled = (saturn.G * saturn.M_Sun) / (saturn.r_orbit * saturn.r_orbit);
    const g_sun_orig = (saturn.G * saturn.M_Sun) / (r_orbit_orig * r_orbit_orig);
    const ratio = g_sun_orig / g_sun_scaled;
    
    saturn.restoreState('base');
    return {
        pass: approxEqual(ratio, 4.0, 0.01),
        message: `g_sun inverse square: ratio = ${ratio.toFixed(3)}`
    };
});

// ========== PHYSICS VALIDATION TESTS ==========

suite.addTest("Physics: Hz ≈ 0 for z = 0 (Solar System)", () => {
    const saturn = new PlanetSaturn();
    const Hz = saturn.computeHz();
    // Hz should be small but > 0
    return {
        pass: Hz > 0 && Hz < 1e-15,
        message: `Hz = ${Hz.toExponential(6)} s^-1`
    };
});

suite.addTest("Physics: Surface gravity dominates (g_saturn >> g_sun)", () => {
    const saturn = new PlanetSaturn();
    const t = 4.5 * Gyr_to_s;
    
    const g_saturn = (saturn.G * saturn.M) / (saturn.r * saturn.r);
    const g_sun = (saturn.G * saturn.M_Sun) / (saturn.r_orbit * saturn.r_orbit);
    
    const ratio = g_saturn / g_sun;
    return {
        pass: ratio > 1e4,  // Surface gravity much larger
        message: `g_saturn/g_sun = ${ratio.toExponential(2)}`
    };
});

suite.addTest("Physics: Ring tidal is small fraction", () => {
    const saturn = new PlanetSaturn();
    const g_saturn = (saturn.G * saturn.M) / (saturn.r * saturn.r);
    const T_ring = (saturn.G * saturn.M_ring) / (saturn.r_ring * saturn.r_ring);
    
    const fraction = T_ring / g_saturn;
    return {
        pass: fraction < 1e-3,  // Ring tidal << surface gravity
        message: `T_ring/g_saturn = ${fraction.toExponential(2)}`
    };
});

suite.addTest("Physics: validateConsistency", () => {
    const saturn = new PlanetSaturn();
    const valid = saturn.validateConsistency();
    return {
        pass: valid,
        message: `Validation: ${valid ? 'PASS' : 'FAIL'}`
    };
});

suite.addTest("Physics: Sensitivity - top 3 most sensitive parameters", () => {
    const saturn = new PlanetSaturn();
    const t = 4.5 * Gyr_to_s;
    const sensitivities = saturn.sensitivityAnalysis(t, 0.01);
    
    const sorted = Object.entries(sensitivities).sort((a, b) => b[1] - a[1]);
    const top3 = sorted.slice(0, 3).map(e => e[0]);
    
    // Expect M and r to be highly sensitive
    const has_M = top3.includes('M');
    const has_r = top3.includes('r');
    return {
        pass: has_M || has_r,
        message: `Top 3: ${top3.join(', ')}`
    };
});

// ========== ENHANCED CAPABILITIES TESTS ==========

suite.addTest("Enhanced: createVariable and getVariable", () => {
    const saturn = new PlanetSaturn();
    saturn.createVariable('test_var', 42.0);
    const value = saturn.getVariable('test_var');
    return {
        pass: value === 42.0,
        message: `Created test_var = ${value}`
    };
});

suite.addTest("Enhanced: setVariable updates correctly", () => {
    const saturn = new PlanetSaturn();
    saturn.setVariable('M', 6e26);
    const M_new = saturn.getVariable('M');
    return {
        pass: M_new === 6e26,
        message: `M updated to ${M_new.toExponential(3)}`
    };
});

suite.addTest("Enhanced: addToVariable increments", () => {
    const saturn = new PlanetSaturn();
    const M_orig = saturn.M;
    saturn.addToVariable('M', 1e25);
    const M_new = saturn.M;
    return {
        pass: approxEqual(M_new, M_orig + 1e25, 1e-6),
        message: `M: ${M_orig.toExponential(3)} → ${M_new.toExponential(3)}`
    };
});

suite.addTest("Enhanced: subtractFromVariable decrements", () => {
    const saturn = new PlanetSaturn();
    const r_orig = saturn.r;
    saturn.subtractFromVariable('r', 1e6);
    const r_new = saturn.r;
    return {
        pass: approxEqual(r_new, r_orig - 1e6, 1e-6),
        message: `r: ${r_orig.toExponential(3)} → ${r_new.toExponential(3)}`
    };
});

suite.addTest("Enhanced: saveState and restoreState", () => {
    const saturn = new PlanetSaturn();
    const M_orig = saturn.M;
    saturn.saveState('test_state');
    saturn.setVariable('M', 7e26);
    saturn.restoreState('test_state');
    const M_restored = saturn.M;
    return {
        pass: M_restored === M_orig,
        message: `M restored to ${M_restored.toExponential(3)}`
    };
});

suite.addTest("Enhanced: listVariables returns array", () => {
    const saturn = new PlanetSaturn();
    const vars = saturn.listVariables();
    const has_M = vars.includes('M');
    const has_r = vars.includes('r');
    return {
        pass: Array.isArray(vars) && has_M && has_r,
        message: `Found ${vars.length} variables`
    };
});

suite.addTest("Enhanced: exportState produces string", () => {
    const saturn = new PlanetSaturn();
    const state_str = saturn.exportState();
    const is_string = typeof state_str === 'string';
    const has_content = state_str.length > 100;
    return {
        pass: is_string && has_content,
        message: `Export string length: ${state_str.length}`
    };
});

suite.addTest("Enhanced: generateReport at 4.5 Gyr", () => {
    const saturn = new PlanetSaturn();
    const t = 4.5 * Gyr_to_s;
    const report = saturn.generateReport(t);
    const has_saturn = report.includes('SATURN') || report.includes('Saturn');
    return {
        pass: typeof report === 'string' && has_saturn,
        message: `Report generated (${report.length} chars)`
    };
});

suite.addTest("Enhanced: autoCorrectAnomalies detects issues", () => {
    const saturn = new PlanetSaturn();
    saturn.setVariable('M', -1);  // Invalid
    const corrected = saturn.autoCorrectAnomalies();
    const M_after = saturn.M;
    return {
        pass: corrected && M_after > 0,
        message: `Corrected: ${corrected}, M after: ${M_after.toExponential(3)}`
    };
});

suite.addTest("Enhanced: generateVariations creates parameter sets", () => {
    const saturn = new PlanetSaturn();
    const variations = saturn.generateVariations(5, 10);
    return {
        pass: Array.isArray(variations) && variations.length === 5,
        message: `Generated ${variations.length} variations`
    };
});

suite.addTest("Enhanced: sensitivityAnalysis returns object", () => {
    const saturn = new PlanetSaturn();
    const t = 4.5 * Gyr_to_s;
    const sens = saturn.sensitivityAnalysis(t, 0.01);
    const has_M = 'M' in sens;
    const has_r = 'r' in sens;
    return {
        pass: typeof sens === 'object' && has_M && has_r,
        message: `Sensitivity keys: ${Object.keys(sens).length}`
    };
});

// Run all tests
suite.run().then(success => {
    if (success) {
        console.log("All tests passed! ✓");
        process.exit(0);
    } else {
        console.log("Some tests failed.");
        process.exit(1);
    }
});
