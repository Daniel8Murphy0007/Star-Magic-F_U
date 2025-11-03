// debug_source34.js
// Debug individual terms

import { MagnetarSGR1745Frequency } from './source34.js';

const sgr = new MagnetarSGR1745Frequency();
const t = 1000 * 3.156e7;

console.log("========== TERM-BY-TERM ANALYSIS ==========\n");

const a_DPM = sgr.computeDPMTerm();
console.log(`a_DPM = ${a_DPM.toExponential(6)} m/s²`);

const a_THz = sgr.computeTHzTerm();
console.log(`a_THz = ${a_THz.toExponential(6)} m/s²`);

const a_vac_diff = sgr.computeVacDiffTerm();
console.log(`a_vac_diff = ${a_vac_diff.toExponential(6)} m/s²`);

const a_super = sgr.computeSuperFreqTerm();
console.log(`a_super = ${a_super.toExponential(6)} m/s²`);

const a_aether_res = sgr.computeAetherResTerm();
console.log(`a_aether_res = ${a_aether_res.toExponential(6)} m/s²`);

const a_u_g4i = sgr.computeU_g4iTerm();
console.log(`a_u_g4i = ${a_u_g4i.toExponential(6)} m/s² <<<< CHECK THIS`);

const a_quantum = sgr.computeQuantumFreqTerm();
console.log(`a_quantum = ${a_quantum.toExponential(6)} m/s²`);

const a_aether_freq = sgr.computeAetherFreqTerm();
console.log(`a_aether_freq = ${a_aether_freq.toExponential(6)} m/s²`);

const a_fluid = sgr.computeFluidFreqTerm();
console.log(`a_fluid = ${a_fluid.toExponential(6)} m/s²`);

const a_osc = sgr.computeOscTerm();
console.log(`a_osc = ${a_osc.toExponential(6)} m/s²`);

const a_exp = sgr.computeExpFreqTerm();
console.log(`a_exp = ${a_exp.toExponential(6)} m/s²`);

console.log("\n========== TOTAL ==========");
const g_total = sgr.compute_g_SGR1745(t);
console.log(`g_UQFF = ${g_total.toExponential(6)} m/s²\n`);

const terms = [
  { name: "DPM", value: a_DPM },
  { name: "THz", value: a_THz },
  { name: "vac_diff", value: a_vac_diff },
  { name: "super", value: a_super },
  { name: "aether_res", value: a_aether_res },
  { name: "u_g4i", value: a_u_g4i },
  { name: "quantum", value: a_quantum },
  { name: "aether_freq", value: a_aether_freq },
  { name: "fluid", value: a_fluid },
  { name: "osc", value: a_osc },
  { name: "exp", value: a_exp }
];

terms.sort((a, b) => Math.abs(b.value) - Math.abs(a.value));

console.log("Terms sorted by magnitude:");
terms.forEach((t, i) => {
  console.log(`${i+1}. ${t.name}: ${t.value.toExponential(4)} m/s²`);
});
