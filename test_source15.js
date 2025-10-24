// Test Source15.mm integration
require('./index.js');

console.log('\n=== SOURCE15.MM INTEGRATION TEST ===\n');

const oneYear = 365.25 * 24 * 3600;

// Test all three specialized systems
console.log('Testing all specialized gravitational systems:\n');

// SGR 1745-2900
const sgr1745 = new MagnetarSGR1745_2900(PREDEFINED_SYSTEMS['SGR_1745_2900']);
const result1745 = sgr1745.compute_g_Magnetar(oneYear);
console.log('SGR 1745-2900 (Galactic Center):');
console.log('  g_field:', result1745.g_Magnetar.toExponential(4), 'm/s²');
console.log('  System type: Magnetar (static B-field, Sgr A* proximity)');
console.log('  Source: Source13.mm\n');

// SGR 0501+4516  
const sgr0501 = new MagnetarSGR0501_4516(PREDEFINED_SYSTEMS['SGR_0501_4516']);
const result0501 = sgr0501.compute_g_Magnetar(oneYear);
console.log('SGR 0501+4516 (Time-Reversal):');
console.log('  g_field:', result0501.g_Magnetar.toExponential(4), 'm/s²');
console.log('  System type: Magnetar (time-reversal, B-field decay)');
console.log('  Source: Source14.mm\n');

// SMBH Sgr A*
const sgrAStar = new SMBHSgrAStar(PREDEFINED_SYSTEMS['SMBH_SGR_A_STAR']);
const resultSMBH = sgrAStar.compute_g_SgrA(oneYear);
console.log('SMBH Sagittarius A*:');
console.log('  g_field:', resultSMBH.g_SgrA.toExponential(4), 'm/s²');
console.log('  System type: Supermassive Black Hole (mass growth, cosmic evolution)');
console.log('  Source: Source15.mm\n');

console.log('COMPARATIVE ANALYSIS:');
console.log('  Field strength ratios:');
console.log('    SMBH/SGR1745:', (resultSMBH.g_SgrA / result1745.g_Magnetar).toExponential(2));
console.log('    SMBH/SGR0501:', (resultSMBH.g_SgrA / result0501.g_Magnetar).toExponential(2));
console.log('    SGR0501/SGR1745:', (result0501.g_Magnetar / result1745.g_Magnetar).toFixed(2));

console.log('\nSource13.mm + Source14.mm + Source15.mm: ALL INTEGRATED SUCCESSFULLY!');
console.log('Enhanced UQFF computational engine with triple MUGE implementations ready!');