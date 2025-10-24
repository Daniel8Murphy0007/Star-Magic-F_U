// CelestialBody.h
#pragma once

#include <string>
#include <vector>

struct CelestialBody {
    std::string name;
    double Ms;          // Mass (kg)
    double Rs;          // Radius (m)
    double Rb;          // Bubble radius (e.g., heliosphere or magnetosphere, m)
    double Ts_surface;  // Surface temperature (K)
    double omega_s;     // Rotation rate (rad/s)
    double Bs_avg;      // Average surface magnetic field (T)
    double SCm_density; // SCm density (kg/m^3)
    double QUA;         // Trapped Universal Aether charge (C)
    double Pcore;       // Planetary core penetration factor
    double PSCm;        // SCm penetration factor
    double omega_c;     // Cycle frequency (rad/s)
};

// Function declarations
double compute_Ug1(const CelestialBody& body, double r, double t, double tn, double alpha, double delta_def, double k1);
double compute_Ug2(const CelestialBody& body, double r, double t, double tn, double k2, double QA, double delta_sw, double v_sw, double HSCm, double rho_A, double kappa);
double compute_Ug3(const CelestialBody& body, double r, double t, double tn, double theta, double rho_A, double kappa, double k3);
double compute_Um(const CelestialBody& body, double t, double tn, double rj, double gamma, double rho_A, double kappa, double num_strings, double phi_hat = 1.0);
void output_json_params(const CelestialBody& body);
std::vector<CelestialBody> load_bodies(const std::string& filename);

// CelestialBody.cpp
#include "CelestialBody.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>

extern const double PI;
extern const double G;
extern double v_SCm;
extern double rho_A;
extern double rho_sw;
extern double QA;
extern double Qs;
extern double kappa;
extern double alpha;
extern double gamma;
extern double delta_sw;
extern double epsilon_sw;
extern double delta_def;
extern double HSCm;
extern double UUA;
extern double eta;
extern double k1, k2, k3, k4;
extern double beta_i;
extern double rho_v;
extern double C_concentration;
extern double f_feedback;
extern double num_strings;
extern double Ts00;
extern std::vector<std::vector<double>> g_mu_nu;
extern double Omega_g;
extern double Mbh;
extern double dg;

// Helper functions
double step_function(double r, double Rb) {
    return (r > Rb) ? 1.0 : 0.0;
}

double compute_Ereact(double t, double rho_SCm, double v_SCm, double rho_A, double kappa) {
    if (rho_A == 0.0) throw std::runtime_error("Division by zero in rho_A");
    return (rho_SCm * v_SCm * v_SCm / rho_A) * std::exp(-kappa * t);
}

double compute_mu_s(double t, double Bs, double omega_c, double Rs, double SCm_contrib = 1e3) {
    double Bs_t = Bs + 0.4 * std::sin(omega_c * t) + SCm_contrib;
    return Bs_t * std::pow(Rs, 3);
}

double compute_grad_Ms_r(double Ms, double Rs) {
    if (Rs == 0.0) throw std::runtime_error("Division by zero in Rs");
    return G * Ms / (Rs * Rs);
}

double compute_Bj(double t, double omega_c, double SCm_contrib = 1e3) {
    return 1e-3 + 0.4 * std::sin(omega_c * t) + SCm_contrib;
}

double compute_omega_s_t(double t, double omega_s, double omega_c) {
    return omega_s - 0.4e-6 * std::sin(omega_c * t);
}

double compute_mu_j(double t, double omega_c, double Rs, double SCm_contrib = 1e3) {
    double Bj = compute_Bj(t, omega_c, SCm_contrib);
    return Bj * std::pow(Rs, 3);
}

// Main functions
double compute_Ug1(const CelestialBody& body, double r, double t, double tn, double alpha, double delta_def, double k1) {
    if (r <= 0.0) throw std::runtime_error("Invalid r value");
    double mu_s = compute_mu_s(t, body.Bs_avg, body.omega_c, body.Rs);
    double grad_Ms_r = compute_grad_Ms_r(body.Ms, body.Rs);
    double defect = 1.0 + delta_def * std::sin(0.001 * t);
    return k1 * mu_s * grad_Ms_r * std::exp(-alpha * t) * std::cos(PI * tn) * defect;
}

double compute_Ug2(const CelestialBody& body, double r, double t, double tn, double k2, double QA, double delta_sw, double v_sw, double HSCm, double rho_A, double kappa) {
    if (r == 0.0) throw std::runtime_error("Division by zero in r");
    double Ereact = compute_Ereact(t, body.SCm_density, v_SCm, rho_A, kappa);
    double S = step_function(r, body.Rb);
    double wind_mod = 1.0 + delta_sw * v_sw;
    return k2 * (QA + body.QUA) * body.Ms / (r * r) * S * wind_mod * HSCm * Ereact;
}

double compute_Ug3(const CelestialBody& body, double r, double t, double tn, double theta, double rho_A, double kappa, double k3) {
    double Ereact = compute_Ereact(t, body.SCm_density, v_SCm, rho_A, kappa);
    double omega_s_t = compute_omega_s_t(t, body.omega_s, body.omega_c);
    double Bj = compute_Bj(t, body.omega_c);
    return k3 * Bj * std::cos(omega_s_t * t * PI) * body.Pcore * Ereact;
}

double compute_Um(const CelestialBody& body, double t, double tn, double rj, double gamma, double rho_A, double kappa, double num_strings, double phi_hat) {
    if (rj == 0.0) throw std::runtime_error("Division by zero in rj");
    double Ereact = compute_Ereact(t, body.SCm_density, v_SCm, rho_A, kappa);
    double mu_j = compute_mu_j(t, body.omega_c, body.Rs);
    double decay = 1.0 - std::exp(-gamma * t * std::cos(PI * tn));
    double single = mu_j / rj * decay * phi_hat;
    return single * num_strings * body.PSCm * Ereact;
}

void output_json_params(const CelestialBody& body) {
    std::cout << "{" << std::endl;
    std::cout << "  \"name\": \"" << body.name << "\"," << std::endl;
    std::cout << "  \"SCm_density\": " << body.SCm_density << "," << std::endl;
    std::cout << "  \"UA\": " << body.QUA << "," << std::endl;
    std::cout << "  \"Qs\": " << Qs << std::endl;
    std::cout << "}" << std::endl;
}

std::vector<CelestialBody> load_bodies(const std::string& filename) {
    std::vector<CelestialBody> bodies;
    std::ifstream in(filename);
    if (!in.is_open()) {
        throw std::runtime_error("Failed to open bodies file: " + filename);
    }
    std::string line;
    while (std::getline(in, line)) {
        std::stringstream ss(line);
        CelestialBody body;
        std::string token;
        std::getline(ss, body.name, ',');
        std::getline(ss, token, ','); body.Ms = std::stod(token);
        std::getline(ss, token, ','); body.Rs = std::stod(token);
        std::getline(ss, token, ','); body.Rb = std::stod(token);
        std::getline(ss, token, ','); body.Ts_surface = std::stod(token);
        std::getline(ss, token, ','); body.omega_s = std::stod(token);
        std::getline(ss, token, ','); body.Bs_avg = std::stod(token);
        std::getline(ss, token, ','); body.SCm_density = std::stod(token);
        std::getline(ss, token, ','); body.QUA = std::stod(token);
        std::getline(ss, token, ','); body.Pcore = std::stod(token);
        std::getline(ss, token, ','); body.PSCm = std::stod(token);
        std::getline(ss, token, ','); body.omega_c = std::stod(token);
        bodies.push_back(body);
    }
    return bodies;
}

// MUGE.h
#pragma once

#include <string>

struct ResonanceParams {
    double fDPM = 1e12;
    double fTHz = 1e12;
    double Evac_neb = 7.09e-36;
    double Evac_ISM = 7.09e-37;
    double Delta_Evac = 6.381e-36;
    double Fsuper = 6.287e-19;
    double UA_SCM = 10;
    double omega_i = 1e-8;
    double k4_res = 1.0;
    double freact = 1e10;
    double fquantum = 1.445e-17;
    double fAether = 1.576e-35;
    double fosc = 4.57e14;
    double fTRZ = 0.1;
    double c_res = 3e8;
};

struct MUGESystem {
    std::string name;
    double I;
    double A;
    double omega1;
    double omega2;
    double Vsys;
    double vexp;
    double t;
    double z;
    double ffluid;
    double M;
    double r;
    double B;
    double Bcrit;
    double rho_fluid;
    double g_local;
    double M_DM;
    double delta_rho_rho;
};

// Compressed MUGE term functions
double compute_compressed_base(const MUGESystem& sys);
double compute_compressed_expansion(const MUGESystem& sys, double H0 = 2.269e-18);
double compute_compressed_super_adj(const MUGESystem& sys);
double compute_compressed_env();
double compute_compressed_Ug_sum();
double compute_compressed_cosm(double Lambda = 1.1e-52);
double compute_compressed_quantum(double hbar = 1.0546e-34, double Delta_x_p = 1e-68, double integral_psi = 2.176e-18, double tHubble = 4.35e17);
double compute_compressed_fluid(const MUGESystem& sys);
double compute_compressed_perturbation(const MUGESystem& sys);

// Resonance MUGE term functions
double compute_aDPM(const MUGESystem& sys, const ResonanceParams& res);
double compute_aTHz(double aDPM, const MUGESystem& sys, const ResonanceParams& res);
double compute_avac_diff(double aDPM, const MUGESystem& sys, const ResonanceParams& res);
double compute_asuper_freq(double aDPM, const ResonanceParams& res);
double compute_aaether_res(double aDPM, const ResonanceParams& res);
double compute_Ug4i(double aDPM, const MUGESystem& sys, const ResonanceParams& res);
double compute_aquantum_freq(double aDPM, const ResonanceParams& res);
double compute_aAether_freq(double aDPM, const ResonanceParams& res);
double compute_afluid_freq(const MUGESystem& sys, const ResonanceParams& res);
double compute_Osc_term();
double compute_aexp_freq(double aDPM, const MUGESystem& sys, const ResonanceParams& res, double H_z = 2.270e-18);
double compute_fTRZ(const ResonanceParams& res);
double compute_a_wormhole(double r, double b = 1.0, double f_worm = 1.0, double Evac_neb = 7.09e-36);

// Full MUGE functions
double compute_compressed_MUGE(const MUGESystem& sys);
double compute_resonance_MUGE(const MUGESystem& sys, const ResonanceParams& res);

std::vector<MUGESystem> load_muge_systems(const std::string& filename);

// MUGE.cpp
#include "MUGE.h"
#include <stdexcept>
#include <fstream>
#include <sstream>

extern const double PI;
extern const double c;
extern const double G;

double compute_compressed_base(const MUGESystem& sys) {
    if (sys.r == 0.0) throw std::runtime_error("Division by zero in r");
    return G * sys.M / (sys.r * sys.r);
}

double compute_compressed_expansion(const MUGESystem& sys, double H0) {
    double H_tz = H0 * sys.t;
    return 1 + H_tz;
}

double compute_compressed_super_adj(const MUGESystem& sys) {
    if (sys.Bcrit == 0.0) throw std::runtime_error("Division by zero in Bcrit");
    return 1 - sys.B / sys.Bcrit;
}

double compute_compressed_env() {
    return 1.0;
}

double compute_compressed_Ug_sum() {
    return 0.0;
}

double compute_compressed_cosm(double Lambda) {
    return Lambda * c * c / 3.0;
}

double compute_compressed_quantum(double hbar, double Delta_x_p, double integral_psi, double tHubble) {
    if (Delta_x_p == 0.0) throw std::runtime_error("Division by zero in Delta_x_p");
    return (hbar / Delta_x_p) * integral_psi * (2 * PI / tHubble);
}

double compute_compressed_fluid(const MUGESystem& sys) {
    return sys.rho_fluid * sys.Vsys * sys.g_local;
}

double compute_compressed_perturbation(const MUGESystem& sys) {
    if (sys.r == 0.0) throw std::runtime_error("Division by zero in r^3");
    return (sys.M + sys.M_DM) * (sys.delta_rho_rho + 3 * G * sys.M / (sys.r * sys.r * sys.r));
}

double compute_compressed_MUGE(const MUGESystem& sys) {
    double base = compute_compressed_base(sys);
    double expansion = compute_compressed_expansion(sys);
    double super_adj = compute_compressed_super_adj(sys);
    double env = compute_compressed_env();
    double adjusted_base = base * expansion * super_adj * env;

    double Ug_sum = compute_compressed_Ug_sum();

    double cosm = compute_compressed_cosm();

    double quantum = compute_compressed_quantum();

    double fluid = compute_compressed_fluid(sys);

    double perturbation = compute_compressed_perturbation(sys);

    return adjusted_base + Ug_sum + cosm + quantum + fluid + perturbation;
}

double compute_aDPM(const MUGESystem& sys, const ResonanceParams& res) {
    double FDPM = sys.I * sys.A * (sys.omega1 - sys.omega2);
    return FDPM * res.fDPM * res.Evac_neb * res.c_res * sys.Vsys;
}

double compute_aTHz(double aDPM, const MUGESystem& sys, const ResonanceParams& res) {
    return res.fTHz * res.Evac_neb * sys.vexp * aDPM / res.Evac_ISM / res.c_res;
}

double compute_avac_diff(double aDPM, const MUGESystem& sys, const ResonanceParams& res) {
    return res.Delta_Evac * sys.vexp * sys.vexp * aDPM / res.Evac_neb / (res.c_res * res.c_res);
}

double compute_asuper_freq(double aDPM, const ResonanceParams& res) {
    return res.Fsuper * res.fTHz * aDPM / res.Evac_neb / res.c_res;
}

double compute_aaether_res(double aDPM, const ResonanceParams& res) {
    return res.UA_SCM * res.omega_i * res.fTHz * aDPM * (1 + res.fTRZ);
}

double compute_Ug4i(double aDPM, const MUGESystem& sys, const ResonanceParams& res) {
    double Ereact = 1046 * std::exp(-0.0005 * sys.t);
    return res.k4_res * Ereact * res.freact * aDPM / res.Evac_neb * res.c_res;
}

double compute_aquantum_freq(double aDPM, const ResonanceParams& res) {
    return res.fquantum * res.Evac_neb * aDPM / res.Evac_ISM / res.c_res;
}

double compute_aAether_freq(double aDPM, const ResonanceParams& res) {
    return res.fAether * res.Evac_neb * aDPM / res.Evac_ISM / res.c_res;
}

double compute_afluid_freq(const MUGESystem& sys, const ResonanceParams& res) {
    return sys.ffluid * res.Evac_neb * sys.Vsys / res.Evac_ISM / res.c_res;
}

double compute_Osc_term() {
    return 0.0;
}

double compute_aexp_freq(double aDPM, const MUGESystem& sys, const ResonanceParams& res, double H_z) {
    double fexp = 2 * PI * H_z * sys.t;
    return fexp * res.Evac_neb * aDPM / res.Evac_ISM / res.c_res;
}

double compute_fTRZ(const ResonanceParams& res) {
    return res.fTRZ;
}

double compute_a_wormhole(double r, double b, double f_worm, double Evac_neb) {
    return f_worm * Evac_neb * (1.0 / (b * b + r * r));
}

double compute_resonance_MUGE(const MUGESystem& sys, const ResonanceParams& res) {
    double aDPM = compute_aDPM(sys, res);
    double aTHz = compute_aTHz(aDPM, sys, res);
    double avac_diff = compute_avac_diff(aDPM, sys, res);
    double asuper_freq = compute_asuper_freq(aDPM, res);
    double aaether_res = compute_aaether_res(aDPM, res);
    double Ug4i = compute_Ug4i(aDPM, sys, res);
    double aquantum_freq = compute_aquantum_freq(aDPM, res);
    double aAether_freq = compute_aAether_freq(aDPM, res);
    double afluid_freq = compute_afluid_freq(sys, res);
    double Osc_term = compute_Osc_term();
    double aexp_freq = compute_aexp_freq(aDPM, sys, res);
    double fTRZ = compute_fTRZ(res);
    double a_worm = compute_a_wormhole(sys.r);

    return aDPM + aTHz + avac_diff + asuper_freq + aaether_res + Ug4i + aquantum_freq + aAether_freq + afluid_freq + Osc_term + aexp_freq + fTRZ + a_worm;
}

std::vector<MUGESystem> load_muge_systems(const std::string& filename) {
    std::vector<MUGESystem> systems;
    std::ifstream in(filename);
    if (!in.is_open()) {
        throw std::runtime_error("Failed to open MUGE systems file: " + filename);
    }
    std::string line;
    while (std::getline(in, line)) {
        std::stringstream ss(line);
        MUGESystem sys;
        std::string token;
        std::getline(ss, sys.name, ',');
        std::getline(ss, token, ','); sys.I = std::stod(token);
        std::getline(ss, token, ','); sys.A = std::stod(token);
        std::getline(ss, token, ','); sys.omega1 = std::stod(token);
        std::getline(ss, token, ','); sys.omega2 = std::stod(token);
        std::getline(ss, token, ','); sys.Vsys = std::stod(token);
        std::getline(ss, token, ','); sys.vexp = std::stod(token);
        std::getline(ss, token, ','); sys.t = std::stod(token);
        std::getline(ss, token, ','); sys.z = std::stod(token);
        std::getline(ss, token, ','); sys.ffluid = std::stod(token);
        std::getline(ss, token, ','); sys.M = std::stod(token);
        std::getline(ss, token, ','); sys.r = std::stod(token);
        std::getline(ss, token, ','); sys.B = std::stod(token);
        std::getline(ss, token, ','); sys.Bcrit = std::stod(token);
        std::getline(ss, token, ','); sys.rho_fluid = std::stod(token);
        std::getline(ss, token, ','); sys.g_local = std::stod(token);
        std::getline(ss, token, ','); sys.M_DM = std::stod(token);
        std::getline(ss, token, ','); sys.delta_rho_rho = std::stod(token);
        systems.push_back(sys);
    }
    return systems;
}

// FluidSolver.h
#pragma once

#include <vector>

class FluidSolver {
public:
    std::vector<double> u, v, u_prev, v_prev, dens, dens_prev;

    FluidSolver();
    void add_source(std::vector<double>& x, std::vector<double>& s);
    void diffuse(int b, std::vector<double>& x, std::vector<double>& x0, double diff);
    void advect(int b, std::vector<double>& d, std::vector<double>& d0);
    void project(std::vector<double>& u, std::vector<double>& v, std::vector<double>& p, std::vector<double>& div);
    void set_bnd(int b, std::vector<double>& x);
    void step(double uqff_g = 0.0);
    void add_jet_force(double force);
    void print_velocity_field();
};

extern const int N;
extern const double dt_ns;
extern const double visc;
extern const double force_jet;

// FluidSolver.cpp
#include "FluidSolver.h"
#include <cmath>
#include <iostream>

#define IX(i, j) ((i) + (N + 2) * (j))

const int N = 32;
const double dt_ns = 0.1;
const double visc = 0.0001;
const double force_jet = 10.0;

FluidSolver::FluidSolver() {
    int size = (N + 2) * (N + 2);
    u.resize(size, 0.0);
    v.resize(size, 0.0);
    u_prev.resize(size, 0.0);
    v_prev.resize(size, 0.0);
    dens.resize(size, 0.0);
    dens_prev.resize(size, 0.0);
}

void FluidSolver::add_source(std::vector<double>& x, std::vector<double>& s) {
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] += dt_ns * s[i];
    }
}

void FluidSolver::diffuse(int b, std::vector<double>& x, std::vector<double>& x0, double diff) {
    double a = dt_ns * diff * N * N;
    for (int k = 0; k < 20; ++k) {
        for (int i = 1; i <= N; ++i) {
            for (int j = 1; j <= N; ++j) {
                x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] +
                    x[IX(i, j - 1)] + x[IX(i, j + 1)])) / (1 + 4 * a);
            }
        }
        set_bnd(b, x);
    }
}

void FluidSolver::advect(int b, std::vector<double>& d, std::vector<double>& d0) {
    int i0, j0, i1, j1;
    double x, y, s0, t0, s1, t1;
    for (int i = 1; i <= N; ++i) {
        for (int j = 1; j <= N; ++j) {
            x = i - dt_ns * N * u[IX(i, j)];
            y = j - dt_ns * N * v[IX(i, j)];
            if (x < 0.5) x = 0.5; if (x > N + 0.5) x = N + 0.5;
            if (y < 0.5) y = 0.5; if (y > N + 0.5) y = N + 0.5;
            i0 = (int)x; i1 = i0 + 1;
            j0 = (int)y; j1 = j0 + 1;
            s1 = x - i0; s0 = 1 - s1;
            t1 = y - j0; t0 = 1 - t1;
            d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
                s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
        }
    }
    set_bnd(b, d);
}

void FluidSolver::project(std::vector<double>& u, std::vector<double>& v, std::vector<double>& p, std::vector<double>& div) {
    double h = 1.0 / N;
    for (int i = 1; i <= N; ++i) {
        for (int j = 1; j <= N; ++j) {
            div[IX(i, j)] = -0.5 * h * (u[IX(i + 1, j)] - u[IX(i - 1, j)] + v[IX(i, j + 1)] - v[IX(i, j - 1)]);
            p[IX(i, j)] = 0;
        }
    }
    set_bnd(0, div); set_bnd(0, p);
    for (int k = 0; k < 20; ++k) {
        for (int i = 1; i <= N; ++i) {
            for (int j = 1; j <= N; ++j) {
                p[IX(i, j)] = (div[IX(i, j)] + p[IX(i - 1, j)] + p[IX(i + 1, j)] +
                    p[IX(i, j - 1)] + p[IX(i, j + 1)]) / 4;
            }
        }
        set_bnd(0, p);
    }
    for (int i = 1; i <= N; ++i) {
        for (int j = 1; j <= N; ++j) {
            u[IX(i, j)] -= 0.5 * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) / h;
            v[IX(i, j)] -= 0.5 * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) / h;
        }
    }
    set_bnd(1, u); set_bnd(2, v);
}

void FluidSolver::set_bnd(int b, std::vector<double>& x) {
    for (int i = 1; i <= N; ++i) {
        x[IX(0, i)] = (b == 1) ? -x[IX(1, i)] : x[IX(1, i)];
        x[IX(N + 1, i)] = (b == 1) ? -x[IX(N, i)] : x[IX(N, i)];
        x[IX(i, 0)] = (b == 2) ? -x[IX(i, 1)] : x[IX(i, 1)];
        x[IX(i, N + 1)] = (b == 2) ? -x[IX(i, N)] : x[IX(i, N)];
    }
    x[IX(0, 0)] = 0.5 * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, N + 1)] = 0.5 * (x[IX(1, N + 1)] + x[IX(0, N)]);
    x[IX(N + 1, 0)] = 0.5 * (x[IX(N, 0)] + x[IX(N + 1, 1)]);
    x[IX(N + 1, N + 1)] = 0.5 * (x[IX(N, N + 1)] + x[IX(N + 1, N)]);
}

void FluidSolver::step(double uqff_g) {
    for (int i = 1; i <= N; ++i) {
        for (int j = 1; j <= N; ++j) {
            v[IX(i, j)] += dt_ns * uqff_g;
        }
    }

    diffuse(1, u_prev, u, visc);
    diffuse(2, v_prev, v, visc);
    project(u_prev, v_prev, u, v);
    advect(1, u, u_prev);
    advect(2, v, v_prev);
    project(u, v, u_prev, v_prev);
}

void FluidSolver::add_jet_force(double force) {
    for (int i = N / 4; i <= 3 * N / 4; ++i) {
        v[IX(i, N / 2)] += force;
    }
}

void FluidSolver::print_velocity_field() {
    std::cout << "Velocity field (magnitude):" << std::endl;
    for (int j = N; j >= 1; --j) {
        for (int i = 1; i <= N; ++i) {
            double mag = std::sqrt(u[IX(i, j)] * u[IX(i, j)] + v[IX(i, j)] * v[IX(i, j)]);
            char sym = (mag > 1.0) ? '#' : (mag > 0.5) ? '+' : (mag > 0.1) ? '.' : ' ';
            std::cout << sym;
        }
        std::cout << std::endl;
    }
}

// UnitTests.h
#pragma once

void run_unit_tests();

// UnitTests.cpp
#include "UnitTests.h"
#include "CelestialBody.h"
#include "MUGE.h"
#include <cassert>
#include <cmath>

void test_compute_compressed_base() {
    MUGESystem test_sys;
    test_sys.M = 1.989e30;
    test_sys.r = 1.496e11;
    double expected = G * test_sys.M / (test_sys.r * test_sys.r);
    double result = compute_compressed_base(test_sys);
    assert(std::abs(result - expected) < 1e-6);
}

void test_compute_compressed_expansion() {
    MUGESystem test_sys;
    test_sys.t = 0.0;
    double expected = 1.0;
    double result = compute_compressed_expansion(test_sys);
    assert(std::abs(result - expected) < 1e-6);
}

void test_compute_compressed_super_adj() {
    MUGESystem test_sys;
    test_sys.B = 1e10;
    test_sys.Bcrit = 1e11;
    double expected = 0.9;
    double result = compute_compressed_super_adj(test_sys);
    assert(std::abs(result - expected) < 1e-6);
}

void test_compute_compressed_fluid() {
    MUGESystem test_sys;
    test_sys.rho_fluid = 1e-15;
    test_sys.Vsys = 4.189e12;
    test_sys.g_local = 10.0;
    double expected = 4.189e-2;
    double result = compute_compressed_fluid(test_sys);
    assert(std::abs(result - expected) < 1e-6);
}

void test_compute_compressed_env() {
    double expected = 1.0;
    double result = compute_compressed_env();
    assert(std::abs(result - expected) < 1e-6);
}

void test_compute_compressed_Ug_sum() {
    double expected = 0.0;
    double result = compute_compressed_Ug_sum();
    assert(std::abs(result - expected) < 1e-6);
}

void test_compute_compressed_cosm() {
    double expected = 1.1e-52 * c * c / 3.0;
    double result = compute_compressed_cosm();
    assert(std::abs(result - expected) / expected < 1e-6);
}

void test_compute_compressed_quantum() {
    double expected = (1.0546e-34 / 1e-68) * 2.176e-18 * (2 * PI / 4.35e17);
    double result = compute_compressed_quantum();
    assert(std::abs(result - expected) / expected < 1e-6);
}

void test_compute_compressed_perturbation() {
    MUGESystem test_sys;
    test_sys.M = 2.984e30;
    test_sys.r = 1e4;
    test_sys.M_DM = 0.0;
    test_sys.delta_rho_rho = 1e-5;
    double expected = test_sys.M * (1e-5 + 3 * G * test_sys.M / (1e4 * 1e4 * 1e4));
    double result = compute_compressed_perturbation(test_sys);
    assert(std::abs(result - expected) / expected < 1e-6);
}

void test_compute_aDPM() {
    MUGESystem test_sys;
    ResonanceParams res;
    test_sys.I = 1e21;
    test_sys.A = 3.142e8;
    test_sys.omega1 = 1e-3;
    test_sys.omega2 = -1e-3;
    test_sys.Vsys = 4.189e12;
    double FDPM = test_sys.I * test_sys.A * (test_sys.omega1 - test_sys.omega2);
    double expected = FDPM * res.fDPM * res.Evac_neb * res.c_res * test_sys.Vsys;
    double result = compute_aDPM(test_sys, res);
    assert(std::abs((result - expected) / expected) < 1e-6);
}

void test_compute_aTHz() {
    ResonanceParams res;
    MUGESystem test_sys;
    double aDPM = 3.545e-42;
    test_sys.vexp = 1e3;
    double expected = 1.182e-33;
    double result = compute_aTHz(aDPM, test_sys, res);
    assert(std::abs((result - expected) / expected) < 1e-6);
}

void test_compute_avac_diff() {
    ResonanceParams res;
    MUGESystem test_sys;
    double aDPM = 3.545e-42;
    test_sys.vexp = 1e3;
    double expected = 3.545e-53;
    double result = compute_avac_diff(aDPM, test_sys, res);
    assert(std::abs((result - expected) / expected) < 1e-6);
}

void test_compute_asuper_freq() {
    ResonanceParams res;
    double aDPM = 3.545e-42;
    double expected = 1.048e-21;
    double result = compute_asuper_freq(aDPM, res);
    assert(std::abs((result - expected) / expected) < 1e-6);
}

void test_compute_aaether_res() {
    ResonanceParams res;
    double aDPM = 3.545e-42;
    double expected = 3.900e-38;
    double result = compute_aaether_res(aDPM, res);
    assert(std::abs((result - expected) / expected) < 1e-6);
}

void test_compute_Ug4i() {
    ResonanceParams res;
    MUGESystem test_sys;
    double aDPM = 3.545e-42;
    test_sys.t = 3.799e10;
    double expected = 0.0;
    double result = compute_Ug4i(aDPM, test_sys, res);
    assert(std::abs(result - expected) < 1e-6);
}

void test_compute_aquantum_freq() {
    ResonanceParams res;
    double aDPM = 3.545e-42;
    double expected = 1.708e-66;
    double result = compute_aquantum_freq(aDPM, res);
    assert(std::abs((result - expected) / expected) < 1e-6);
}

void test_compute_aAether_freq() {
    ResonanceParams res;
    double aDPM = 3.545e-42;
    double expected = 1.863e-84;
    double result = compute_aAether_freq(aDPM, res);
    assert(std::abs((result - expected) / expected) < 1e-6);
}

void test_compute_afluid_freq() {
    ResonanceParams res;
    MUGESystem test_sys;
    test_sys.ffluid = 1.269e-14;
    test_sys.Vsys = 4.189e12;
    double expected = 1.773e-9;
    double result = compute_afluid_freq(test_sys, res);
    assert(std::abs((result - expected) / expected) < 1e-6);
}

void test_compute_Osc_term() {
    double expected = 0.0;
    double result = compute_Osc_term();
    assert(std::abs(result - expected) < 1e-6);
}

void test_compute_aexp_freq() {
    ResonanceParams res;
    MUGESystem test_sys;
    double aDPM = 3.545e-42;
    test_sys.t = 3.799e10;
    double expected = 1.623e-57;
    double result = compute_aexp_freq(aDPM, test_sys, res);
    assert(std::abs((result - expected) / expected) < 1e-6);
}

void test_compute_fTRZ() {
    ResonanceParams res;
    double expected = 0.1;
    double result = compute_fTRZ(res);
    assert(std::abs(result - expected) < 1e-6);
}

void test_compute_compressed_MUGE() {
    MUGESystem test_sys;
    test_sys.name = "Test";
    test_sys.M = 2.984e30;
    test_sys.r = 1e4;
    test_sys.B = 1e10;
    test_sys.Bcrit = 1e11;
    test_sys.rho_fluid = 1e-15;
    test_sys.Vsys = 4.189e12;
    test_sys.g_local = 10.0;
    test_sys.M_DM = 0.0;
    test_sys.delta_rho_rho = 1e-5;
    test_sys.t = 3.799e10;
    test_sys.z = 0.0009;
    test_sys.ffluid = 1.269e-14;
    double expected = 1.782e39;  // Approximate from attachment
    double result = compute_compressed_MUGE(test_sys);
    assert(std::abs((result - expected) / expected) < 1e-3);
}

void test_compute_resonance_MUGE() {
    ResonanceParams res;
    MUGESystem test_sys;
    test_sys.name = "Test";
    test_sys.I = 1e21;
    test_sys.A = 3.142e8;
    test_sys.omega1 = 1e-3;
    test_sys.omega2 = -1e-3;
    test_sys.Vsys = 4.189e12;
    test_sys.vexp = 1e3;
    test_sys.t = 3.799e10;
    test_sys.ffluid = 1.269e-14;
    test_sys.r = 1e4;
    double expected = 1.773e-9;  // Approximate from attachment
    double result = compute_resonance_MUGE(test_sys, res);
    assert(std::abs((result - expected) / expected) < 1e-3);
}

void test_compute_a_wormhole() {
    double r = 1e4;
    double b = 1.0;
    double expected = 7.09e-36 / (1.0 + r * r);
    double result = compute_a_wormhole(r);
    assert(std::abs((result - expected) / expected) < 1e-6);
}

void run_unit_tests() {
    try {
        test_compute_compressed_base();
        test_compute_compressed_expansion();
        test_compute_compressed_super_adj();
        test_compute_compressed_fluid();
        test_compute_compressed_env();
        test_compute_compressed_Ug_sum();
        test_compute_compressed_cosm();
        test_compute_compressed_quantum();
        test_compute_compressed_perturbation();
        test_compute_compressed_MUGE();
        test_compute_aDPM();
        test_compute_aTHz();
        test_compute_avac_diff();
        test_compute_asuper_freq();
        test_compute_aaether_res();
        test_compute_Ug4i();
        test_compute_aquantum_freq();
        test_compute_aAether_freq();
        test_compute_afluid_freq();
        test_compute_Osc_term();
        test_compute_aexp_freq();
        test_compute_fTRZ();
        test_compute_resonance_MUGE();
        test_compute_a_wormhole();
        std::cout << "All unit tests passed!" << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "Unit test failed: " << e.what() << std::endl;
    }
}

// main.cpp
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include "CelestialBody.h"
#include "MUGE.h"
#include "FluidSolver.h"
#include "UnitTests.h"

const double PI = 3.141592653589793;
const double c = 3.0e8;
const double G = 6.67430e-11;

double Omega_g = 7.3e-16;
double Mbh = 8.15e36;
double dg = 2.55e20;

double v_SCm = 0.99 * c;
double rho_A = 1e-23;
double rho_sw = 8e-21;
double v_sw = 5e5;
double QA = 1e-10;
double Qs = 0.0;
double kappa = 0.0005;
double alpha = 0.001;
double gamma = 0.00005;
double delta_sw = 0.01;
double epsilon_sw = 0.001;
double delta_def = 0.01;
double HSCm = 1.0;
double UUA = 1.0;
double eta = 1e-22;
double k1 = 1.5, k2 = 1.2, k3 = 1.8, k4 = 2.0;
double beta_i = 0.6;
double rho_v = 6e-27;
double C_concentration = 1.0;
double f_feedback = 0.1;
const double num_strings = 1e9;
double Ts00 = 1.27e3 + 1.11e7;
std::vector<std::vector<double>> g_mu_nu = {
    {1.0, 0.0, 0.0, 0.0},
    {0.0, -1.0, 0.0, 0.0},
    {0.0, 0.0, -1.0, 0.0},
    {0.0, 0.0, 0.0, -1.0}
};

double compute_Ug4(double t, double tn, double rho_v, double C_concentration, double Mbh, double dg, double alpha, double f_feedback, double k4) {
    double decay = std::exp(-alpha * t);
    double cycle = std::cos(PI * tn);
    return k4 * rho_v * C_concentration * Mbh / dg * decay * cycle * (1 + f_feedback);
}

double compute_Ubi(double Ugi, double beta_i, double Omega_g, double Mbh, double dg, double epsilon_sw, double rho_sw, double UUA, double tn) {
    double wind_mod = 1.0 + epsilon_sw * rho_sw;
    return -beta_i * Ugi * Omega_g * Mbh / dg * wind_mod * UUA * std::cos(PI * tn);
}

std::vector<std::vector<double>> compute_A_mu_nu(double tn, double eta, double Ts00) {
    std::vector<std::vector<double>> A = g_mu_nu;
    double mod = eta * Ts00 * std::cos(PI * tn);
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            A[i][j] += mod;
        }
    }
    return A;
}

double compute_FU(const CelestialBody& body, double r, double t, double tn, double theta) {
    try {
        double Ug1 = compute_Ug1(body, r, t, tn, alpha, delta_def, k1);
        double Ug2 = compute_Ug2(body, r, t, tn, k2, QA, delta_sw, v_sw, HSCm, rho_A, kappa);
        double Ug3 = compute_Ug3(body, r, t, tn, theta, rho_A, kappa, k3);
        double Ug4 = compute_Ug4(t, tn, rho_v, C_concentration, Mbh, dg, alpha, f_feedback, k4);
        double sum_Ugi = Ug1 + Ug2 + Ug3 + Ug4;

        double Ubi1 = compute_Ubi(Ug1, beta_i, Omega_g, Mbh, dg, epsilon_sw, rho_sw, UUA, tn);
        double Ubi2 = compute_Ubi(Ug2, beta_i, Omega_g, Mbh, dg, epsilon_sw, rho_sw, UUA, tn);
        double Ubi3 = compute_Ubi(Ug3, beta_i, Omega_g, Mbh, dg, epsilon_sw, rho_sw, UUA, tn);
        double Ubi4 = compute_Ubi(Ug4, beta_i, Omega_g, Mbh, dg, epsilon_sw, rho_sw, UUA, tn);
        double sum_Ubi = Ubi1 + Ubi2 + Ubi3 + Ubi4;

        double Um = compute_Um(body, t, tn, body.Rb, gamma, rho_A, kappa, num_strings);

        auto A = compute_A_mu_nu(tn, eta, Ts00);
        double A_scalar = A[0][0] + A[1][1] + A[2][2] + A[3][3];

        return sum_Ugi + sum_Ubi + Um + A_scalar;
    }
    catch (const std::exception& e) {
        std::cerr << "Error in compute_FU: " << e.what() << std::endl;
        return 0.0;
    }
}

void simulate_quasar_jet(double initial_velocity) {
    try {
        FluidSolver solver;
        solver.add_jet_force(initial_velocity / 10.0);

        ResonanceParams res;
        MUGESystem sagA;  // Placeholder, define as needed
        sagA.name = "Sagittarius A*";
        sagA.I = 1e23;
        sagA.A = 2.813e30;
        sagA.omega1 = 1e-5;
        sagA.omega2 = -1e-5;
        sagA.Vsys = 3.552e45;
        sagA.vexp = 5e6;
        sagA.t = 3.786e14;
        sagA.ffluid = 3.465e-8;
        sagA.r = 1e12;
        double uqff_g = compute_resonance_MUGE(sagA, res);

        std::cout << "Simulating quasar jet with Navier-Stokes (10 steps) using UQFF g=" << uqff_g << "..." << std::endl;
        for (int step = 0; step < 10; ++step) {
            solver.step(uqff_g / 1e30);
        }
        solver.print_velocity_field();
    }
    catch (const std::exception& e) {
        std::cerr << "Error in simulate_quasar_jet: " << e.what() << std::endl;
    }
}

void write_velocity_to_csv(const FluidSolver& solver, const std::string& filename) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        throw std::runtime_error("Failed to open CSV file: " + filename);
    }
    for (int j = 1; j <= N; ++j) {
        for (int i = 1; i <= N; ++i) {
            double mag = std::sqrt(solver.u[IX(i, j)] * solver.u[IX(i, j)] + solver.v[IX(i, j)] * solver.v[IX(i, j)]);
            out << mag << ",";
        }
        out << std::endl;
    }
}

void print_summary_stats(const std::vector<double>& values, const std::string& name) {
    if (values.empty()) return;
    double min = *std::min_element(values.begin(), values.end());
    double max = *std::max_element(values.begin(), values.end());
    double sum = std::accumulate(values.begin(), values.end(), 0.0);
    double mean = sum / values.size();
    std::cout << name << " summary - Min: " << min << ", Max: " << max << ", Mean: " << mean << std::endl;
}

int main(int argc, char** argv) {
    std::string input_file_bodies, input_file_muge, output_file;
    for (int i = 1; i < argc; i += 2) {
        std::string arg = argv[i];
        if (arg == "--input-bodies" && i + 1 < argc) {
            input_file_bodies = argv[i + 1];
        }
        else if (arg == "--input-muge" && i + 1 < argc) {
            input_file_muge = argv[i + 1];
        }
        else if (arg == "--output" && i + 1 < argc) {
            output_file = argv[i + 1];
        }
    }

    std::vector<double> fu_values, compressed_values, resonance_values;

    try {
        std::vector<CelestialBody> bodies = input_file_bodies.empty() ? std::vector<CelestialBody>() : load_bodies(input_file_bodies);
        if (bodies.empty()) {
            CelestialBody sun = { "Sun", 1.989e30, 6.96e8, 1.496e13, 5778.0, 2.5e-6, 1e-4, 1e15, 1e-11, 1.0, 1.0, 2 * PI / (11.0 * 365.25 * 24 * 3600) };
            CelestialBody earth = { "Earth", 5.972e24, 6.371e6, 1e7, 288.0, 7.292e-5, 3e-5, 1e12, 1e-12, 1e-3, 1e-3, 2 * PI / (1.0 * 365.25 * 24 * 3600) };
            CelestialBody jupiter = { "Jupiter", 1.898e27, 6.9911e7, 1e8, 165.0, 1.76e-4, 4e-4, 1e13, 1e-11, 1e-3, 1e-3, 2 * PI / (11.86 * 365.25 * 24 * 3600) };
            CelestialBody neptune = { "Neptune", 1.024e26, 2.4622e7, 5e7, 72.0, 1.08e-4, 1e-4, 1e11, 1e-13, 1e-3, 1e-3, 2 * PI / (164.8 * 365.25 * 24 * 3600) };
            bodies = { sun, earth, jupiter, neptune };
        }

        double r = 1e13;
        double t = 0.0;
        double tn = t;
        double theta = 0.0;

        for (const auto& body : bodies) {
            r = body.Rb;
            double FU = compute_FU(body, r, t, tn, theta);
            fu_values.push_back(FU);
            std::cout << "Unified Field Strength (FU) for " << body.name << " at t=" << t << ", r=" << r << ": " << FU << " (normalized units)" << std::endl;

            double Ug1 = compute_Ug1(body, r, t, tn, alpha, delta_def, k1);
            std::cout << "Ug1: " << Ug1 << std::endl;
            double Ug2 = compute_Ug2(body, r, t, tn, k2, QA, delta_sw, v_sw, HSCm, rho_A, kappa);
            std::cout << "Ug2: " << Ug2 << std::endl;
            double Ug3 = compute_Ug3(body, r, t, tn, theta, rho_A, kappa, k3);
            std::cout << "Ug3: " << Ug3 << std::endl;
            double Ug4 = compute_Ug4(t, tn, rho_v, C_concentration, Mbh, dg, alpha, f_feedback, k4);
            std::cout << "Ug4: " << Ug4 << std::endl;
            double Um = compute_Um(body, t, tn, body.Rb, gamma, rho_A, kappa, num_strings);
            std::cout << "Um: " << Um << std::endl;

            auto A = compute_A_mu_nu(tn, eta, Ts00);
            std::cout << "A_mu_nu trace: " << A[0][0] + A[1][1] + A[2][2] + A[3][3] << std::endl;

            std::cout << "JSON parameters for " << body.name << ":" << std::endl;
            output_json_params(body);
            std::cout << std::endl;
        }

        print_summary_stats(fu_values, "FU");

        simulate_quasar_jet(v_SCm);

        if (!output_file.empty()) {
            FluidSolver solver;  // Assume solver is populated
            write_velocity_to_csv(solver, output_file + ".csv");
        }

        std::vector<MUGESystem> muge_systems = input_file_muge.empty() ? std::vector<MUGESystem>() : load_muge_systems(input_file_muge);
        if (muge_systems.empty()) {
            MUGESystem sgr1745 = { "Magnetar SGR 1745-2900", 1e21, 3.142e8, 1e-3, -1e-3, 4.189e12, 1e3, 3.799e10, 0.0009, 1.269e-14, 2.984e30, 1e4, 1e10, 1e11, 1e-15, 10.0, 0.0, 1e-5 };
            MUGESystem sagA = { "Sagittarius A*", 1e23, 2.813e30, 1e-5, -1e-5, 3.552e45, 5e6, 3.786e14, 0.0009, 3.465e-8, 8.155e36, 1e12, 1e-5, 1e-4, 1e-20, 1e-5, 1e37, 1e-3 };
            MUGESystem tapestry = { "Tapestry of Blazing Starbirth", 1e22, 1e35, 1e-4, -1e-4, 1e53, 1e4, 3.156e13, 0.0, 1e-12, 1.989e35, 3.086e17, 1e-4, 1e-3, 1e-21, 1e-8, 1e35, 1e-4 };
            MUGESystem westerlund = tapestry; westerlund.name = "Westerlund 2";
            MUGESystem pillars = { "Pillars of Creation", 1e21, 2.813e32, 1e-3, -1e-3, 3.552e48, 2e3, 3.156e13, 0.0, 8.457e-14, 1.989e32, 9.46e15, 1e-4, 1e-3, 1e-21, 1e-8, 0.0, 1e-5 };
            MUGESystem rings = { "Rings of Relativity", 1e22, 1e35, 1e-4, -1e-4, 1e54, 1e5, 3.156e14, 0.01, 1e-9, 1.989e36, 3.086e17, 1e-5, 1e-4, 1e-20, 1e-5, 1e36, 1e-3 };
            MUGESystem student_guide = { "Student’s Guide to the Universe", 1e24, 1e52, 1e-6, -1e-6, 1e80, 3e8, 4.35e17, 0.0, 1e-18, 1e53, 1e26, 1e-10, 1e-9, 1e-30, 1e-10, 1e53, 1e-6 };
            muge_systems = { sgr1745, sagA, tapestry, westerlund, pillars, rings, student_guide };
        }

        for (const auto& sys : muge_systems) {
            double compressed_g = compute_compressed_MUGE(sys);
            compressed_values.push_back(compressed_g);
            double resonance_g = compute_resonance_MUGE(sys, ResonanceParams{});
            resonance_values.push_back(resonance_g);
            std::cout << "Compressed MUGE g for " << sys.name << ": " << compressed_g << " m/s2" << std::endl;
            std::cout << "Resonance MUGE g for " << sys.name << ": " << resonance_g << " m/s2" << std::endl;
        }

        print_summary_stats(compressed_values, "Compressed MUGE");
        print_summary_stats(resonance_values, "Resonance MUGE");

        run_unit_tests();
    }
    catch (const std::exception& e) {
        std::cerr << "Main error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

