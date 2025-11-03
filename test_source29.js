/**
 * Comprehensive Test Suite for Sombrero Galaxy (M104) UQFF Module
 * Testing all physics computations, dynamic capabilities, and unique Sombrero features
 */

import GalaxySombrero from './source29.js';

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
        console.log("SOMBRERO GALAXY (M104) COMPREHENSIVE TEST SUITE");
        console.log("Testing: The Hat - Dust Lane & Superconductivity");
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

const M_sun = 1.989e30;
const Gyr_to_s = 1e9 * 3.156e7;
const ly_to_m = 9.461e15;

// ========== INITIALIZATION TESTS ==========

suite.addTest("Initialization: M and r within reasonable ranges", () => {
    const sombrero = new GalaxySombrero();
    const M_check = sombrero.M > 1e40 && sombrero.M < 1e43;  // 1e10-1e12 M_sun
    const r_check = sombrero.r > 1e19 && sombrero.r < 1e22;  // ~1-100 kly
    return {
        pass: M_check && r_check,
        message: `M = ${sombrero.M.toExponential(3)}, r = ${sombrero.r.toExponential(3)}`
    };
});

suite.addTest("Initialization: Mass distribution (80% visible, 20% DM)", () => {
    const sombrero = new GalaxySombrero();
    const visible_frac = sombrero.M_visible / sombrero.M;
    const dm_frac = sombrero.M_DM / sombrero.M;
    const visible_ok = approxEqual(visible_frac, 0.8, 0.01);
    const dm_ok = approxEqual(dm_frac, 0.2, 0.01);
    return {
        pass: visible_ok && dm_ok,
        message: `Visible: ${visible_frac.toFixed(3)}, DM: ${dm_frac.toFixed(3)}`
    };
});

suite.addTest("Initialization: Redshift z = 0.0063 (Virgo Cluster)", () => {
    const sombrero = new GalaxySombrero();
    const z_ok = approxEqual(sombrero.z, 0.0063, 0.001);
    return {
        pass: z_ok,
        message: `z = ${sombrero.z}, rel_diff: ${relDiff(sombrero.z, 0.0063).toExponential(2)}`
    };
});

suite.addTest("Initialization: Massive SMBH (M_BH = 1e9 M_sun)", () => {
    const sombrero = new GalaxySombrero();
    const M_BH_check = approxEqual(sombrero.M_BH, 1e9 * M_sun, 0.1);
    return {
        pass: M_BH_check,
        message: `M_BH = ${(sombrero.M_BH / M_sun).toExponential(3)} M_sun`
    };
});

suite.addTest("Initialization: Prominent dust lane (rho_dust = 1e-20)", () => {
    const sombrero = new GalaxySombrero();
    const dust_ok = approxEqual(sombrero.rho_dust, 1e-20, 0.1);
    return {
        pass: dust_ok,
        message: `rho_dust = ${sombrero.rho_dust.toExponential(3)} kg/m^3`
    };
});

suite.addTest("Initialization: Superconductivity parameters (B, B_crit)", () => {
    const sombrero = new GalaxySombrero();
    const B_ok = approxEqual(sombrero.B, 1e-5, 0.1);
    const B_crit_ok = approxEqual(sombrero.B_crit, 1e11, 0.1);
    return {
        pass: B_ok && B_crit_ok,
        message: `B = ${sombrero.B.toExponential(3)} T, B_crit = ${sombrero.B_crit.toExponential(3)} T`
    };
});

// ========== UQFF COMPONENT TESTS ==========

suite.addTest("UQFF: Ug sum = Ug1 + Ug4 (Ug2, Ug3 ≈ 0)", () => {
    const sombrero = new GalaxySombrero();
    const ug_sum = sombrero.computeUgSum();
    const ug_manual = sombrero.Ug1 + sombrero.Ug2 + sombrero.Ug3 + sombrero.Ug4;
    return {
        pass: approxEqual(ug_sum, ug_manual, 1e-6),
        message: `Ug_sum = ${ug_sum.toExponential(6)}, manual = ${ug_manual.toExponential(6)}, rel_diff: ${relDiff(ug_sum, ug_manual).toExponential(2)}`
    };
});

suite.addTest("UQFF: g(t) > 0 for all times [0, 13.8 Gyr]", () => {
    const sombrero = new GalaxySombrero();
    const times = [0, 2.5, 5, 7.5, 10, 13.8].map(t => t * Gyr_to_s);
    let all_positive = true;
    for (const t of times) {
        const g = sombrero.compute_g_Sombrero(t);
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
    const sombrero = new GalaxySombrero();
    const g = sombrero.compute_g_Sombrero(-1e9);
    return {
        pass: g === 0,
        message: `g(-1e9 s) = ${g}`
    };
});

suite.addTest("UQFF: Superconductivity correction (1 - B/B_crit)", () => {
    const sombrero = new GalaxySombrero();
    const sc_corr_expected = 1.0 - sombrero.B / sombrero.B_crit;
    // SC correction applied to base gravity, verify it's close to 1 for small B/B_crit
    const ratio = sombrero.B / sombrero.B_crit;
    const sc_ok = ratio < 0.01;  // Should be ~1e-6 for typical values
    return {
        pass: sc_ok,
        message: `SC = ${sc_corr_expected.toFixed(9)}, B/B_crit = ${ratio.toExponential(2)}`
    };
});

// ========== GALAXY SCALING TESTS ==========

suite.addTest("Scaling: expandGalaxyScale M × 2 → Ug1 × 2", () => {
    const sombrero = new GalaxySombrero();
    sombrero.saveState('base');
    const Ug1_base = sombrero.Ug1;
    
    sombrero.expandGalaxyScale(2.0, 1.0);
    const Ug1_scaled = sombrero.Ug1;
    
    const ratio = Ug1_scaled / Ug1_base;
    sombrero.restoreState('base');
    return {
        pass: approxEqual(ratio, 2.0, 0.01),
        message: `Ug1: ${Ug1_base.toExponential(6)} → ${Ug1_scaled.toExponential(6)}, ratio: ${ratio.toFixed(3)}`
    };
});

suite.addTest("Scaling: expandGalaxyScale r × 2 → Ug1 × 1/4", () => {
    const sombrero = new GalaxySombrero();
    sombrero.saveState('base');
    const Ug1_base = sombrero.Ug1;
    
    sombrero.expandGalaxyScale(1.0, 2.0);
    const Ug1_scaled = sombrero.Ug1;
    
    const ratio = Ug1_base / Ug1_scaled;  // Should be 4
    sombrero.restoreState('base');
    return {
        pass: approxEqual(ratio, 4.0, 0.01),
        message: `Ug1: ${Ug1_base.toExponential(6)} → ${Ug1_scaled.toExponential(6)}, ratio: ${ratio.toFixed(3)}`
    };
});

suite.addTest("Scaling: expandBlackHoleScale M_BH × 2", () => {
    const sombrero = new GalaxySombrero();
    sombrero.saveState('base');
    const M_BH_base = sombrero.M_BH;
    
    sombrero.expandBlackHoleScale(2.0, 1.0);
    const M_BH_scaled = sombrero.M_BH;
    
    const ratio = M_BH_scaled / M_BH_base;
    sombrero.restoreState('base');
    return {
        pass: approxEqual(ratio, 2.0, 1e-6),
        message: `M_BH: ${M_BH_base.toExponential(6)} → ${M_BH_scaled.toExponential(6)}, ratio: ${ratio.toFixed(3)}`
    };
});

suite.addTest("Scaling: expandBlackHoleScale r_BH × 3", () => {
    const sombrero = new GalaxySombrero();
    sombrero.saveState('base');
    const r_BH_base = sombrero.r_BH;
    
    sombrero.expandBlackHoleScale(1.0, 3.0);
    const r_BH_scaled = sombrero.r_BH;
    
    const ratio = r_BH_scaled / r_BH_base;
    sombrero.restoreState('base');
    return {
        pass: approxEqual(ratio, 3.0, 1e-6),
        message: `r_BH: ${r_BH_base.toExponential(6)} → ${r_BH_scaled.toExponential(6)}, ratio: ${ratio.toFixed(3)}`
    };
});

suite.addTest("Scaling: expandDustLaneScale rho_dust × 2 (UNIQUE)", () => {
    const sombrero = new GalaxySombrero();
    sombrero.saveState('base');
    const rho_dust_base = sombrero.rho_dust;
    
    sombrero.expandDustLaneScale(2.0, 1.0);
    const rho_dust_scaled = sombrero.rho_dust;
    
    const ratio = rho_dust_scaled / rho_dust_base;
    sombrero.restoreState('base');
    return {
        pass: approxEqual(ratio, 2.0, 1e-6),
        message: `rho_dust: ${rho_dust_base.toExponential(6)} → ${rho_dust_scaled.toExponential(6)}, ratio: ${ratio.toFixed(3)}`
    };
});

suite.addTest("Scaling: expandDustLaneScale B × 3 (UNIQUE)", () => {
    const sombrero = new GalaxySombrero();
    sombrero.saveState('base');
    const B_base = sombrero.B;
    
    sombrero.expandDustLaneScale(1.0, 3.0);
    const B_scaled = sombrero.B;
    
    const ratio = B_scaled / B_base;
    sombrero.restoreState('base');
    return {
        pass: approxEqual(ratio, 3.0, 1e-6),
        message: `B: ${B_base.toExponential(6)} → ${B_scaled.toExponential(6)}, ratio: ${ratio.toFixed(3)}`
    };
});

// ========== PARAMETER EFFECTS TESTS ==========

suite.addTest("Parameters: M_BH affects black hole contribution", () => {
    const sombrero = new GalaxySombrero();
    sombrero.saveState('base');
    const t = 10 * Gyr_to_s;
    
    const M_BH_orig = sombrero.M_BH;
    sombrero.setVariable('M_BH', M_BH_orig * 2);
    const g_BH_new = (sombrero.G * sombrero.M_BH) / (sombrero.r_BH * sombrero.r_BH);
    
    sombrero.restoreState('base');
    const g_BH_orig = (sombrero.G * sombrero.M_BH) / (sombrero.r_BH * sombrero.r_BH);
    
    const ratio = g_BH_new / g_BH_orig;
    return {
        pass: approxEqual(ratio, 2.0, 0.01),
        message: `g_BH ratio: ${ratio.toFixed(3)}`
    };
});

suite.addTest("Parameters: rho_dust affects dust term directly", () => {
    const sombrero = new GalaxySombrero();
    sombrero.saveState('base');
    
    const dust_base = sombrero.computeDustTerm();
    sombrero.setVariable('rho_dust', sombrero.rho_dust * 2);
    const dust_scaled = sombrero.computeDustTerm();
    
    const ratio = dust_scaled / dust_base;
    sombrero.restoreState('base');
    return {
        pass: approxEqual(ratio, 2.0, 0.01),
        message: `Dust term: ${dust_base.toExponential(6)} → ${dust_scaled.toExponential(6)}, ratio: ${ratio.toFixed(3)}`
    };
});

suite.addTest("Parameters: v_orbit affects dust term (v² dependence)", () => {
    const sombrero = new GalaxySombrero();
    sombrero.saveState('base');
    
    const dust_base = sombrero.computeDustTerm();
    sombrero.setVariable('v_orbit', sombrero.v_orbit * 2);
    const dust_scaled = sombrero.computeDustTerm();
    
    const ratio = dust_scaled / dust_base;
    sombrero.restoreState('base');
    return {
        pass: approxEqual(ratio, 4.0, 0.01),
        message: `Dust term with 2× v: ${dust_base.toExponential(6)} → ${dust_scaled.toExponential(6)}, ratio: ${ratio.toFixed(3)}`
    };
});

suite.addTest("Parameters: B affects SC correction (1 - B/B_crit)", () => {
    const sombrero = new GalaxySombrero();
    sombrero.saveState('base');
    
    const sc_base = 1.0 - sombrero.B / sombrero.B_crit;
    sombrero.setVariable('B', sombrero.B * 2);
    const sc_scaled = 1.0 - sombrero.B / sombrero.B_crit;
    
    // SC correction should decrease (approach 0 as B → B_crit)
    const decreases = sc_scaled < sc_base;
    sombrero.restoreState('base');
    return {
        pass: decreases,
        message: `SC: ${sc_base.toFixed(9)} → ${sc_scaled.toFixed(9)}, decreases: ${decreases}`
    };
});

suite.addTest("Parameters: B_crit affects SC correction inversely", () => {
    const sombrero = new GalaxySombrero();
    sombrero.saveState('base');
    
    const sc_base = 1.0 - sombrero.B / sombrero.B_crit;
    sombrero.setVariable('B_crit', sombrero.B_crit * 2);
    const sc_scaled = 1.0 - sombrero.B / sombrero.B_crit;
    
    // SC correction should increase (closer to 1) when B_crit increases
    const increases = sc_scaled > sc_base;
    sombrero.restoreState('base');
    return {
        pass: increases,
        message: `SC: ${sc_base.toFixed(9)} → ${sc_scaled.toFixed(9)}, increases: ${increases}`
    };
});

// ========== PHYSICS VALIDATION TESTS ==========

suite.addTest("Physics: Hz > 0 for z = 0.0063", () => {
    const sombrero = new GalaxySombrero();
    const Hz = sombrero.computeHz();
    return {
        pass: Hz > 0 && Hz < 1e-15,  // Reasonable range for H0
        message: `Hz = ${Hz.toExponential(6)} s^-1`
    };
});

suite.addTest("Physics: Dust term dominates (like source28)", () => {
    const sombrero = new GalaxySombrero();
    const t = 10 * Gyr_to_s;
    const g_total = sombrero.compute_g_Sombrero(t);
    const dust_term = sombrero.computeDustTerm();
    
    // Dust should be significant portion (similar to source28)
    const dust_fraction = Math.abs(dust_term / g_total);
    return {
        pass: dust_fraction > 0.5,  // Dust likely dominates
        message: `g_total: ${g_total.toExponential(6)}, dust: ${dust_term.toExponential(6)}, fraction: ${dust_fraction.toFixed(3)}`
    };
});

suite.addTest("Physics: Consistency check validateConsistency", () => {
    const sombrero = new GalaxySombrero();
    const valid = sombrero.validateConsistency();
    return {
        pass: valid,
        message: `Validation: ${valid ? 'PASS' : 'FAIL'}`
    };
});

suite.addTest("Physics: Sensitivity - top 3 most sensitive parameters", () => {
    const sombrero = new GalaxySombrero();
    const t = 10 * Gyr_to_s;
    const sensitivities = sombrero.sensitivityAnalysis(t, 0.01);
    
    const sorted = Object.entries(sensitivities).sort((a, b) => b[1] - a[1]);
    const top3 = sorted.slice(0, 3).map(e => e[0]);
    
    // Expect dust-related parameters to be highly sensitive
    const has_dust = top3.includes('rho_dust') || top3.includes('v_orbit');
    return {
        pass: has_dust,
        message: `Top 3: ${top3.join(', ')}`
    };
});

// ========== ENHANCED CAPABILITIES TESTS ==========

suite.addTest("Enhanced: createVariable and getVariable", () => {
    const sombrero = new GalaxySombrero();
    sombrero.createVariable('test_var', 42.0);
    const value = sombrero.getVariable('test_var');
    return {
        pass: value === 42.0,
        message: `Created test_var = ${value}`
    };
});

suite.addTest("Enhanced: setVariable updates correctly", () => {
    const sombrero = new GalaxySombrero();
    sombrero.setVariable('M', 2e41);
    const M_new = sombrero.getVariable('M');
    return {
        pass: M_new === 2e41,
        message: `M updated to ${M_new.toExponential(3)}`
    };
});

suite.addTest("Enhanced: addToVariable increments", () => {
    const sombrero = new GalaxySombrero();
    const M_orig = sombrero.M;
    sombrero.addToVariable('M', 1e40);
    const M_new = sombrero.M;
    return {
        pass: approxEqual(M_new, M_orig + 1e40, 1e-6),
        message: `M: ${M_orig.toExponential(3)} → ${M_new.toExponential(3)}`
    };
});

suite.addTest("Enhanced: subtractFromVariable decrements", () => {
    const sombrero = new GalaxySombrero();
    const r_orig = sombrero.r;
    sombrero.subtractFromVariable('r', 1e19);
    const r_new = sombrero.r;
    return {
        pass: approxEqual(r_new, r_orig - 1e19, 1e-6),
        message: `r: ${r_orig.toExponential(3)} → ${r_new.toExponential(3)}`
    };
});

suite.addTest("Enhanced: saveState and restoreState", () => {
    const sombrero = new GalaxySombrero();
    const M_orig = sombrero.M;
    sombrero.saveState('test_state');
    sombrero.setVariable('M', 9e41);
    sombrero.restoreState('test_state');
    const M_restored = sombrero.M;
    return {
        pass: M_restored === M_orig,
        message: `M restored to ${M_restored.toExponential(3)}`
    };
});

suite.addTest("Enhanced: listVariables returns array", () => {
    const sombrero = new GalaxySombrero();
    const vars = sombrero.listVariables();
    const has_M = vars.includes('M');
    const has_r = vars.includes('r');
    return {
        pass: Array.isArray(vars) && has_M && has_r,
        message: `Found ${vars.length} variables`
    };
});

suite.addTest("Enhanced: exportState produces string", () => {
    const sombrero = new GalaxySombrero();
    const state_str = sombrero.exportState();
    const is_string = typeof state_str === 'string';
    const has_content = state_str.length > 100;
    return {
        pass: is_string && has_content,
        message: `Export string length: ${state_str.length}`
    };
});

suite.addTest("Enhanced: generateReport at 10 Gyr", () => {
    const sombrero = new GalaxySombrero();
    const t = 10 * Gyr_to_s;
    const report = sombrero.generateReport(t);
    const has_M104 = report.includes('M104') || report.includes('SOMBRERO');
    return {
        pass: typeof report === 'string' && has_M104,
        message: `Report generated (${report.length} chars)`
    };
});

suite.addTest("Enhanced: autoCorrectAnomalies detects issues", () => {
    const sombrero = new GalaxySombrero();
    sombrero.setVariable('M', -1);  // Invalid
    const corrected = sombrero.autoCorrectAnomalies();
    const M_after = sombrero.M;
    return {
        pass: corrected && M_after > 0,
        message: `Corrected: ${corrected}, M after: ${M_after.toExponential(3)}`
    };
});

suite.addTest("Enhanced: generateVariations creates parameter sets", () => {
    const sombrero = new GalaxySombrero();
    const variations = sombrero.generateVariations(5, 10);
    return {
        pass: Array.isArray(variations) && variations.length === 5,
        message: `Generated ${variations.length} variations`
    };
});

suite.addTest("Enhanced: sensitivityAnalysis returns object", () => {
    const sombrero = new GalaxySombrero();
    const t = 10 * Gyr_to_s;
    const sens = sombrero.sensitivityAnalysis(t, 0.01);
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
