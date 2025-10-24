
// === Source105.mm: GalacticBlackHoleModule Integration ===
class GalacticBlackHoleModule {
    constructor() {
        this.variables = new Map();
        this.variables.set('M_sun', 1.989e30);
        this.variables.set('M_bh', 8.15e36);
        this.variables.set('beta_1', 0.6);
        this.variables.set('U_g1', 1.39e26);
        this.variables.set('Omega_g', 7.3e-16);
        this.variables.set('d_g', 2.55e20);
        this.variables.set('epsilon_sw', 0.001);
        this.variables.set('rho_vac_sw', 8e-21);
        this.variables.set('U_UA', 1.0);
        this.variables.set('t_n', 0.0);
        this.variables.set('pi', Math.PI);
        this.variables.set('k_4', 1.0);
        this.variables.set('rho_vac_SCm', 7.09e-37);
        this.variables.set('alpha', 0.001);
        this.variables.set('f_feedback', 0.1);
    }

    updateVariable(name, value) {
        this.variables.set(name, value);
    }

    computeM_bh() { return this.variables.get('M_bh'); }
    computeM_bhInMsun() { return this.computeM_bh() / this.variables.get('M_sun'); }
    computeMbhOverDg() { return this.computeM_bh() / this.variables.get('d_g'); }

    computeU_b1() {
        const beta_1 = this.variables.get('beta_1');
        const U_g1 = this.variables.get('U_g1');
        const Omega_g = this.variables.get('Omega_g');
        const mbh_over_dg = this.computeMbhOverDg();
        const swirl_factor = 1.0 + this.variables.get('epsilon_sw') * this.variables.get('rho_vac_sw');
        const U_UA = this.variables.get('U_UA');
        const cos_term = Math.cos(this.variables.get('pi') * this.variables.get('t_n'));
        return -beta_1 * U_g1 * Omega_g * mbh_over_dg * swirl_factor * U_UA * cos_term;
    }

    computeU_g4() {
        const k_4 = this.variables.get('k_4');
        const rho_vac_SCm = this.variables.get('rho_vac_SCm');
        const mbh_over_dg = this.computeMbhOverDg();
        const exp_term = Math.exp(-this.variables.get('alpha') * this.variables.get('t_n'));
        const cos_term = Math.cos(this.variables.get('pi') * this.variables.get('t_n'));
        const feedback_factor = 1.0 + this.variables.get('f_feedback');
        return k_4 * (rho_vac_SCm * this.computeM_bh() / this.variables.get('d_g')) * exp_term * cos_term * feedback_factor;
    }

    getEquationText() {
        return 'U_bi = -beta_i U_gi Omega_g (M_bh / d_g) (1 + epsilon_sw rho_vac,sw) U_UA cos(pi t_n)\nU_g4 = k_4 (rho_vac,[SCm] M_bh / d_g) exp(-alpha t) cos(pi t_n) (1 + f_feedback)\nWhere M_bh = 8.15e36 kg (4.1e6 M_sun) for Sgr A*';
    }
}

function analyzeGalacticBlackHoleUQFF105(timePoints = [0, 86400, 31536000]) {
    console.log('\n ANALYZING Galactic Black Hole UQFF (Source105.mm)');
    console.log('===================================================');
    const gbhm = new GalacticBlackHoleModule();
    console.log('M_bh = 8.15e36 kg; M_bh/d_g = 3.20e16 kg/m');
    console.log(gbhm.getEquationText());
    
    timePoints.forEach((t_seconds, index) => {
        const t_days = t_seconds / 86400;
        console.log('--- Time Point ' + (index + 1) + ': t = ' + t_days.toFixed(1) + ' days ---');
        gbhm.updateVariable('t_n', t_seconds);
        const m_bh = gbhm.computeM_bh();
        const u_b1 = gbhm.computeU_b1();
        const u_g4 = gbhm.computeU_g4();
        console.log('   M_bh = ' + m_bh.toExponential(3) + ' kg');
        console.log('   U_b1 = ' + u_b1.toExponential(3) + ' J/m^3');
        console.log('   U_g4 = ' + u_g4.toExponential(3) + ' J/m^3');
    });
    return gbhm;
}

// Export Source105.mm classes and functions
if (typeof module !== 'undefined' && module.exports) {
    module.exports.GalacticBlackHoleModule = GalacticBlackHoleModule;
    module.exports.analyzeGalacticBlackHoleUQFF105 = analyzeGalacticBlackHoleUQFF105;
}
