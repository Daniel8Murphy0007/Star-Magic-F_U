/**
 * Simple inline test for HUDF Galaxies module
 */

import HUDFGalaxies from './source26.js';

console.log("Running inline test for HUDF Galaxies...\n");

const hudf = new HUDFGalaxies();
const M_sun = 1.989e30;
const Gyr_to_s = 1e9 * 3.156e7;

console.log("Step 1: Initial State");
console.log(`System: ${hudf.getSystemName()}`);
console.log(`Validation: ${hudf.validateConsistency() ? 'PASS' : 'FAIL'}\n`);

console.log("Step 2: Time Evolution (M(t) and I(t))");
const t_test = [0, 1 * Gyr_to_s, 5 * Gyr_to_s];
t_test.forEach((t, idx) => {
    const t_Gyr = [0, 1, 5][idx];
    const Mt = hudf.M_t(t);
    const It = hudf.I_t(t);
    const g = hudf.compute_g_HUDF(t);
    console.log(`  t = ${t_Gyr} Gyr: M(t) = ${(Mt/M_sun).toExponential(4)} M_sun, I(t) = ${It.toFixed(6)}, g = ${g.toExponential(4)} m/s^2`);
});

console.log("\n✓ Inline test complete - HUDF Galaxies module operational\n");
