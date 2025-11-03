// debug_source33.js
// Debug individual terms in SGR 1745-2900 computation

import { MagnetarSGR1745 } from './source33.js';

const sgr = new MagnetarSGR1745();
const t = 1000 * 3.156e7; // 1000 years

console.log("========== TERM-BY-TERM ANALYSIS ==========\n");

// Base components
const Hz = sgr.computeHz();
const expansion = 1.0 + Hz * t;
const sc_correction = sgr.computeSCCorrection();
const tr_factor = 1.0 + sgr.f_TRZ;

console.log(`Hz = ${Hz.toExponential(6)} s⁻¹`);
console.log(`expansion = ${expansion.toExponential(6)}`);
console.log(`sc_correction = ${sc_correction.toFixed(6)}`);
console.log(`tr_factor = ${tr_factor.toFixed(6)}\n`);

// Base gravity
const g_base = (sgr.G * sgr.M / (sgr.r * sgr.r)) * expansion * sc_correction * tr_factor;
console.log(`g_base = ${g_base.toExponential(6)} m/s²\n`);

// Individual terms
const ug_sum = sgr.computeUgSum();
console.log(`Ug sum = ${ug_sum.toExponential(6)} m/s²`);
console.log(`  Ug1 = ${sgr.Ug1.toExponential(6)}`);
console.log(`  Ug4 = ${sgr.Ug4.toExponential(6)}\n`);

const lambda_term = sgr.Lambda * (sgr.c * sgr.c) / 3.0;
console.log(`Lambda term = ${lambda_term.toExponential(6)} m/s²\n`);

const quantum_term = sgr.computeQuantumTerm();
console.log(`Quantum term = ${quantum_term.toExponential(6)} m/s²\n`);

const em_term = sgr.computeEMTerm();
console.log(`EM term = ${em_term.toExponential(6)} m/s²\n`);

const fluid_term = sgr.computeFluidTerm(g_base);
console.log(`Fluid term = ${fluid_term.toExponential(6)} m/s²\n`);

const resonant_term = sgr.computeResonantTerm(t);
console.log(`Resonant term = ${resonant_term.toExponential(6)} m/s²\n`);

const dm_term = sgr.computeDMTerm();
console.log(`DM term = ${dm_term.toExponential(6)} m/s²\n`);

// Total
const g_total = g_base + ug_sum + lambda_term + quantum_term + em_term + fluid_term + resonant_term + dm_term;
console.log(`========== TOTAL ==========`);
console.log(`g_UQFF = ${g_total.toExponential(6)} m/s²\n`);

console.log(`Expected from C++: ~1.2e12 m/s²`);
console.log(`Actual: ${g_total.toExponential(2)} m/s²\n`);

// Find dominant term
const terms = [
  { name: "g_base", value: g_base },
  { name: "ug_sum", value: ug_sum },
  { name: "lambda", value: lambda_term },
  { name: "quantum", value: quantum_term },
  { name: "EM", value: em_term },
  { name: "fluid", value: fluid_term },
  { name: "resonant", value: resonant_term },
  { name: "DM", value: dm_term }
];

terms.sort((a, b) => Math.abs(b.value) - Math.abs(a.value));

console.log("Terms sorted by magnitude:");
terms.forEach((t, i) => {
  console.log(`${i+1}. ${t.name}: ${t.value.toExponential(4)} m/s²`);
});
