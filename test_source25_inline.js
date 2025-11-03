/**
 * Simple inline test for NGC 1275 module
 */

import NGC1275 from './source25.js';

console.log("Running inline test for NGC 1275...\n");

const ngc = new NGC1275();
const M_sun = 1.989e30;
const Myr_to_s = 1e6 * 3.156e7;

console.log("Step 1: Initial State");
console.log(`System: ${ngc.getSystemName()}`);
console.log(`Validation: ${ngc.validateConsistency() ? 'PASS' : 'FAIL'}\n`);

console.log("Step 2: Time Evolution (B(t) and F(t))");
const t_test = [0, 50 * Myr_to_s, 100 * Myr_to_s];
t_test.forEach((t, idx) => {
    const t_Myr = [0, 50, 100][idx];
    const Bt = ngc.B_t(t);
    const Ft = ngc.F_t(t);
    const g = ngc.compute_g_NGC1275(t);
    console.log(`  t = ${t_Myr} Myr: B(t) = ${Bt.toExponential(4)} T, F(t) = ${Ft.toFixed(6)}, g = ${g.toExponential(4)} m/s^2`);
});

console.log("\n✓ Inline test complete - NGC 1275 module operational\n");
