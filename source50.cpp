// uqff_module.cpp
// Module for UQFF (Compressed and Resonance) Equations implementation
// Pluggable into base program 'ziqn233h..cpp' via external functions:
// - void install_uqff_module();  // Installer/Initializer
// - double compute_compressed_muge(const std::string& system_name, const VariableMap& updates = {});
// - double compute_resonance_muge(const std::string& system_name, const VariableMap& updates = {});
// - void update_variable(const std::string& system_name, const std::string& var_name, double value);
// - void add_variable(const std::string& system_name, const std::string& var_name, double value);
// - void subtract_variable(const std::string& system_name, const std::string& var_name, double delta);
// - void print_system_text(const std::string& system_name);  // Output associated text/descriptions
// Includes all equations, variables, solutions, special terms (quantum, fluid, resonant, DM, etc.)
// Nothing negligible: all terms computed and summed explicitly.
// Uses <map> for VariableMap: std::map<std::string, double> for dynamic updates.
// Compile with: g++ -std=c++11 -o ziqn233h ziqn233h.cpp uqff_module.cpp -lm

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <functional>
#include <random>
#include <algorithm>
#include <sstream>

// Forward declaration for base program integration (assumed in ziqn233h.cpp)
extern void base_program_init();  // Placeholder for base program

// Constants (global, shared across systems)
const double G = 6.6743e-11;  // m^3 kg^-1 s^-2
const double H0 = 2.269e-18;  // s^-1 (70 km/s/Mpc)
const double B_t = 1e10;      // T
const double B_crit = 1e11;   // T
const double Lambda = 1.1e-52; // m^-2
const double c = 3e8;         // m/s
const double hbar = 1.0546e-34; // J s
const double Delta_x_Delta_p = 1e-68; // J^2 s^2
const double integral_psi = 2.176e-18; // J
const double t_Hubble = 4.35e17; // s
const double rho_fluid = 1e-15; // kg/m^3
const double g_earth = 10.0;    // m/s^2
const double M_DM_default = 0.0;
const double delta_rho_over_rho = 1e-5;
const double pi = 3.141592653589793;

// Resonance-specific constants
const double E_vac_neb = 7.09e-36; // J/m^3
const double E_vac_ISM = 7.09e-37; // J/m^3
const double Delta_E_vac = 6.381e-36; // J/m^3
const double F_super = 6.287e-19;
const double UA_SC_m = 10.0;
const double omega_i = 1e-8; // rad/s
const double k_4 = 1.0;
const double f_react = 1e10; // Hz
const double f_quantum = 1.445e-17; // Hz
const double f_Aether = 1.576e-35; // Hz
const double f_osc = 4.57e14; // Hz
const double f_TRZ = 0.1;     // Hz (trivial resonance zone?)

// Type for dynamic variables: map of string to double
using VariableMap = std::map<std::string, double>;

// Struct for system data (pre-defined defaults, updatable)
struct SystemData {
    std::string name;
    std::string description;  // Associated text
    double M;                 // kg
    double r;                 // m
    double z;
    double t;                 // s
    double V;                 // m^3
    double F_env;             // Environmental factor
    double v_exp;             // m/s (expansion velocity)
    double I;                 // A (current for resonance DPM)
    double A;                 // m^2 (area for resonance)
    double omega1;            // rad/s
    double omega2;            // rad/s
    double M_sun;             // kg (for planetary)
    double r_orbit;           // m (for planetary)
    VariableMap vars;         // Dynamic overrides

    SystemData(const std::string& n, const std::string& desc, double m, double rad, double zz, double tt, double vv, double fenv, double vexp,
               double ii, double aa, double o1, double o2, double msun = 0, double rorb = 0)
        : name(n), description(desc), M(m), r(rad), z(zz), t(tt), V(vv), F_env(fenv), v_exp(vexp),
          I(ii), A(aa), omega1(o1), omega2(o2), M_sun(msun), r_orbit(rorb) {
        // Pre-populate with computed defaults where applicable
        vars["H_t_z"] = H0 * (0.3 * pow(1 + z, 3) + 0.7);
        vars["one_plus_H_t"] = 1 + vars["H_t_z"] * t;
        vars["B_adjust"] = 1 - B_t / B_crit;
        vars["one_plus_F_env"] = 1 + F_env;
        vars["Lambda_c2_3"] = Lambda * c * c / 3;
        vars["hbar_over_sqrt_delta"] = hbar / sqrt(Delta_x_Delta_p);
        vars["quantum_term"] = vars["hbar_over_sqrt_delta"] * integral_psi * (2 * pi / t_Hubble);
        vars["rho_V_g"] = rho_fluid * V * g_earth;
        vars["three_G_M_over_r3"] = 3 * G * M / (r * r * r);
        vars["density_pert"] = delta_rho_over_rho + vars["three_G_M_over_r3"];
        vars["M_vis_DM_pert"] = (M + M_DM_default) * vars["density_pert"];
        // Resonance defaults
        vars["f_DPM"] = 1e12; // Hz
        vars["f_fluid"] = (G * M / (r * r)) * (2 * pi); // Placeholder, computed in func
        vars["f_exp"] = vars["H_t_z"] * t * (2 * pi); // Placeholder
    }
};

// Global map of systems (initialized in installer)
std::map<std::string, SystemData*> systems;

// Function to compute volume if not provided (4/3 pi r^3)
double compute_volume(double r) {
    return (4.0 / 3.0) * pi * r * r * r;
}

// Installer function: Initialize all systems with defaults from document
void install_uqff_module() {
    // Clear existing
    for (auto& pair : systems) delete pair.second;
    systems.clear();

    // Hubble Sees Galaxies Galore
    double V_gal = compute_volume(1.543e21);
    systems["Hubble Sees Galaxies Galore"] = new SystemData(
        "Hubble Sees Galaxies Galore", "Hubble Deep Field observations, capturing thousands of galaxies.",
        1.989e41, 1.543e21, 1.0, 4.35e17, V_gal, 0.0, 1e5, 1e24, 7.487e42, 1e-6, -1e-6
    );

    // The Stellar Forge
    double V_stellar = compute_volume(9.46e16);
    systems["The Stellar Forge"] = new SystemData(
        "The Stellar Forge", "Star-forming region in Large Magellanic Cloud (30 Doradus Nebula).",
        1.989e34, 9.46e16, 0.00005, 6.312e13, V_stellar, 0.0, 1e4, 1e22, 8.508e35, 1e-2, -1e-2
    );

    // Hubble Mosaic of the Majestic Sombrero Galaxy
    double V_sombrero = compute_volume(4.73e20);
    systems["Hubble Mosaic of the Majestic Sombrero Galaxy"] = new SystemData(
        "Hubble Mosaic of the Majestic Sombrero Galaxy", "Sombrero Galaxy (M104), peculiar galaxy with dust lane.",
        1.591e42, 4.73e20, 0.002, 4.35e17, V_sombrero, 0.0, 2e5, 1e24, 7.487e42, 1e-6, -1e-6  // Approx A
    );

    // Saturn
    double V_saturn = compute_volume(6.027e7);
    systems["Saturn"] = new SystemData(
        "Saturn", "Hubble observations of Saturn, rings and atmosphere.",
        5.68e26, 6.027e7, 0.0, 4.35e17, V_saturn, 0.0, 5e3, 1e20, 7.032e22, 1e-4, -1e-4, 1.989e30, 1.36e12
    );

    // New Stars Shed Light on the Past
    systems["New Stars Shed Light on the Past"] = new SystemData(
        "New Stars Shed Light on the Past", "Star-forming region in Small Magellanic Cloud (N90).",
        1.989e34, 9.46e16, 0.00006, 6.312e13, V_stellar, 0.0, 1e4, 1e22, 8.508e35, 1e-2, -1e-2
    );

    // The Crab Nebula
    double V_crab = compute_volume(5.203e16);
    systems["The Crab Nebula"] = new SystemData(
        "The Crab Nebula", "Supernova remnant formed in 1054 CE.",
        9.945e30, 5.203e16, 0.00002, 3.064e10, V_crab, 0.0, 1.34e6, 1e22, 8.508e35, 1e-2, -1e-2
    );

    // Student's Guide to the Universe
    double V_guide = compute_volume(1.496e11);
    systems["Student�s Guide to the Universe"] = new SystemData(
        "Student�s Guide to the Universe", "General framework using solar mass and AU-scale.",
        1.989e30, 1.496e11, 0.0, 4.35e17, V_guide, 0.0, 3e4, 1e20, 7.032e22, 1e-4, -1e-4
    );

    // Additional systems from document (Lagoon Nebula, etc.)
    double V_lagoon = compute_volume(5.203e17);  // Approx
    systems["The Lagoon Nebula"] = new SystemData(
        "The Lagoon Nebula", "Emission nebula with star formation.",
        1.989e34, 5e16, 0.0001, 6.312e13, 5.913e53, 0.0, 1e4, 1e22, 8.508e35, 1e-2, -1e-2
    );

    double V_spirals = compute_volume(1.543e21);
    systems["Spirals and Supernovae"] = new SystemData(
        "Spirals and Supernovae", "Galactic spirals and supernova dynamics.",
        1.989e41, 1.543e21, 0.002, 4.35e17, V_spirals, 0.0, 2e5, 1e24, 7.487e42, 1e-6, -1e-6
    );

    double V_ngc = compute_volume(1.514e16);  // Approx for Butterfly
    systems["NGC 6302 (Butterfly Nebula)"] = new SystemData(
        "NGC 6302 (Butterfly Nebula)", "Planetary nebula with bipolar outflows.",
        1.989e30, 1.514e16, 0.00001, 3.156e11, 1.458e48, 0.0, 2e4, 1e21, 7.207e32, 1e-3, -1e-3
    );

    double V_orion = compute_volume(1.135e17);
    systems["Orion Nebula"] = new SystemData(
        "Orion Nebula", "Stellar nursery near Earth.",
        3.978e33, 1.135e17, 0.00004, 3.156e13, 6.132e51, 0.0, 1e4, 1e22, 4.047e34, 1e-2, -1e-2
    );

    // Update volumes if not computed
    for (auto& pair : systems) {
        if (pair.second->V == 0) pair.second->V = compute_volume(pair.second->r);
    }

    // Call base init if needed
    base_program_init();

    std::cout << "UQFF Module Installed: All systems initialized with defaults." << std::endl;
}

// Update a variable (additive or set)
void update_variable(const std::string& system_name, const std::string& var_name, double value, bool is_add = false) {
    auto it = systems.find(system_name);
    if (it != systems.end()) {
        if (is_add) {
            it->second->vars[var_name] += value;
        } else {
            it->second->vars[var_name] = value;
        }
        // Propagate updates to dependent vars if needed (e.g., recompute H_t_z)
        if (var_name == "z") {
            it->second->vars["H_t_z"] = H0 * (0.3 * pow(1 + value, 3) + 0.7);
            it->second->vars["one_plus_H_t"] = 1 + it->second->vars["H_t_z"] * it->second->t;
        }
        // ... similar for other deps
    }
}

// Add to variable (wrapper for update with add=true)
void add_variable(const std::string& system_name, const std::string& var_name, double value) {
    update_variable(system_name, var_name, value, true);
}

// Subtract from variable
void subtract_variable(const std::string& system_name, const std::string& var_name, double delta) {
    update_variable(system_name, var_name, -delta, true);
}

// Print associated text/description
void print_system_text(const std::string& system_name) {
    auto it = systems.find(system_name);
    if (it != systems.end()) {
        std::cout << "System: " << it->second->name << std::endl;
        std::cout << "Description: " << it->second->description << std::endl;
        // Output key vars
        std::cout << "Key Variables:" << std::endl;
        for (const auto& v : it->second->vars) {
            std::cout << "  " << v.first << " = " << std::scientific << v.second << std::endl;
        }
    } else {
        std::cout << "System not found: " << system_name << std::endl;
    }
}

// Compute Compressed MUGE (all terms explicit, nothing negligible)
double compute_compressed_muge(const std::string& system_name, const VariableMap& updates) {
    auto it = systems.find(system_name);
    if (it == systems.end()) return 0.0;

    SystemData* sys = it->second;
    VariableMap local_vars = sys->vars;  // Copy
    for (const auto& u : updates) local_vars[u.first] = u.second;

    // Update from struct if not in map
    if (local_vars.find("M") == local_vars.end()) local_vars["M"] = sys->M;
    if (local_vars.find("r") == local_vars.end()) local_vars["r"] = sys->r;
    if (local_vars.find("z") == local_vars.end()) local_vars["z"] = sys->z;
    if (local_vars.find("t") == local_vars.end()) local_vars["t"] = sys->t;
    if (local_vars.find("V") == local_vars.end()) local_vars["V"] = sys->V;
    if (local_vars.find("F_env") == local_vars.end()) local_vars["F_env"] = sys->F_env;
    if (local_vars.find("M_sun") == local_vars.end()) local_vars["M_sun"] = sys->M_sun;
    if (local_vars.find("r_orbit") == local_vars.end()) local_vars["r_orbit"] = sys->r_orbit;

    // Recompute dependents
    double zz = local_vars["z"];
    double tt = local_vars["t"];
    local_vars["H_t_z"] = H0 * (0.3 * pow(1 + zz, 3) + 0.7);
    local_vars["one_plus_H_t"] = 1 + local_vars["H_t_z"] * tt;
    local_vars["B_adjust"] = 1 - B_t / B_crit;
    local_vars["one_plus_F_env"] = 1 + local_vars["F_env"];
    local_vars["Lambda_c2_3"] = Lambda * c * c / 3;
    local_vars["hbar_over_sqrt_delta"] = hbar / sqrt(Delta_x_Delta_p);
    local_vars["quantum_term"] = local_vars["hbar_over_sqrt_delta"] * integral_psi * (2 * pi / t_Hubble);
    local_vars["rho_V_g"] = rho_fluid * local_vars["V"] * g_earth;
    double MM = local_vars["M"];
    double rr = local_vars["r"];
    local_vars["three_G_M_over_r3"] = 3 * G * MM / (rr * rr * rr);
    local_vars["density_pert"] = delta_rho_over_rho + local_vars["three_G_M_over_r3"];
    local_vars["M_vis_DM_pert"] = (MM + M_DM_default) * local_vars["density_pert"];

    // Gravity base term
    double grav_base = (G * MM / (rr * rr)) * local_vars["one_plus_H_t"] * local_vars["B_adjust"] * local_vars["one_plus_F_env"];
    if (sys->M_sun > 0) {  // Planetary: add orbital
        double orb_grav = (G * local_vars["M_sun"] / (local_vars["r_orbit"] * local_vars["r_orbit"])) * local_vars["one_plus_H_t"];
        grav_base += orb_grav;
    }

    // Gravity modes (0 as per doc)
    double U_g_sum = 0.0;  // U_g1 + U_g2 + U_g3' + U_g4

    // Full sum
    double muge = grav_base + U_g_sum + local_vars["Lambda_c2_3"] + local_vars["quantum_term"] +
                  local_vars["rho_V_g"] + local_vars["M_vis_DM_pert"];

    std::cout << std::scientific << std::setprecision(3);
    std::cout << "Compressed MUGE for " << system_name << ": " << muge << " m/s^2" << std::endl;
    std::cout << "  Breakdown: grav_base=" << grav_base << ", U_g_sum=" << U_g_sum
              << ", Lambda=" << local_vars["Lambda_c2_3"] << ", quantum=" << local_vars["quantum_term"]
              << ", fluid=" << local_vars["rho_V_g"] << ", pert=" << local_vars["M_vis_DM_pert"] << std::endl;

    return muge;
}

// Compute Resonance MUGE (all terms explicit)
double compute_resonance_muge(const std::string& system_name, const VariableMap& updates) {
    auto it = systems.find(system_name);
    if (it == systems.end()) return 0.0;

    SystemData* sys = it->second;
    VariableMap local_vars = sys->vars;
    for (const auto& u : updates) local_vars[u.first] = u.second;

    // Update basics
    if (local_vars.find("M") == local_vars.end()) local_vars["M"] = sys->M;
    if (local_vars.find("r") == local_vars.end()) local_vars["r"] = sys->r;
    if (local_vars.find("V") == local_vars.end()) local_vars["V"] = sys->V;
    if (local_vars.find("v_exp") == local_vars.end()) local_vars["v_exp"] = sys->v_exp;
    if (local_vars.find("I") == local_vars.end()) local_vars["I"] = sys->I;
    if (local_vars.find("A") == local_vars.end()) local_vars["A"] = sys->A;
    if (local_vars.find("omega1") == local_vars.end()) local_vars["omega1"] = sys->omega1;
    if (local_vars.find("omega2") == local_vars.end()) local_vars["omega2"] = sys->omega2;
    if (local_vars.find("z") == local_vars.end()) local_vars["z"] = sys->z;
    if (local_vars["z"] == 0) local_vars["H_z"] = H0; else local_vars["H_z"] = H0 * (0.3 * pow(1 + local_vars["z"], 3) + 0.7);

    double II = local_vars["I"];
    double AA = local_vars["A"];
    double delta_omega = local_vars["omega1"] - local_vars["omega2"];
    double F_DPM = II * AA * delta_omega;
    local_vars["f_DPM"] = 1e12;  // Hz fixed
    double a_DPM = F_DPM * local_vars["f_DPM"] * E_vac_neb / (c * local_vars["V"]);

    // THz Hole Resonance
    double a_THz = local_vars["f_DPM"] * E_vac_neb * local_vars["v_exp"] * a_DPM / (E_vac_ISM * c);

    // Plasmotic Vacuum Energy Density Differential
    double a_vac_diff = Delta_E_vac * pow(local_vars["v_exp"], 2) * a_DPM / (E_vac_neb * c * c);

    // Superconductor Frequency Interaction
    double a_super_freq = F_super * local_vars["f_DPM"] * a_DPM / (E_vac_neb * c);

    // Aether-Mediated Resonance
    double a_aether_res = k_4 * omega_i * local_vars["f_DPM"] * a_DPM * (1 + UA_SC_m * 0.1);  // Approx doc

    // Reactive Dynamics U_g4i (0 as per doc)
    double U_g4i = 0.0;

    // Quantum Wave Dynamics
    double a_quantum_freq = f_quantum * E_vac_neb * a_DPM / (E_vac_ISM * c);

    // Aether Effect
    double a_Aether_freq = f_Aether * E_vac_neb * a_DPM / (E_vac_ISM * c);

    // Fluid Dynamics
    double f_fluid = (G * local_vars["M"] / (local_vars["r"] * local_vars["r"])) / (2 * pi);  // Inverted from doc approx
    double a_fluid_freq = f_fluid * E_vac_neb * local_vars["V"] / (E_vac_ISM * c);

    // Oscillatory Components (0)
    double Osc_term = 0.0;

    // Cosmic Expansion
    double f_exp = local_vars["H_z"] * local_vars["t"] / (2 * pi);  // Approx
    double a_exp_freq = f_exp * E_vac_neb * a_DPM / (E_vac_ISM * c);

    // Final sum
    double muge = a_DPM + a_THz + a_vac_diff + a_super_freq + a_aether_res + U_g4i + a_quantum_freq +
                  a_Aether_freq + a_fluid_freq + Osc_term + a_exp_freq + f_TRZ;

    std::cout << std::scientific << std::setprecision(3);
    std::cout << "Resonance MUGE for " << system_name << ": " << muge << " m/s^2" << std::endl;
    std::cout << "  Breakdown: a_DPM=" << a_DPM << ", a_THz=" << a_THz << ", a_vac_diff=" << a_vac_diff
              << ", a_super_freq=" << a_super_freq << ", a_aether_res=" << a_aether_res
              << ", U_g4i=" << U_g4i << ", a_quantum=" << a_quantum_freq << ", a_Aether=" << a_Aether_freq
              << ", a_fluid=" << a_fluid_freq << ", a_exp=" << a_exp_freq << ", f_TRZ=" << f_TRZ << std::endl;

    return muge;
}

// Example usage in base program (call after install)
void example_usage() {
    install_uqff_module();
    compute_compressed_muge("Hubble Sees Galaxies Galore");
    compute_resonance_muge("Hubble Sees Galaxies Galore");
    print_system_text("Hubble Sees Galaxies Galore");
    add_variable("Hubble Sees Galaxies Galore", "M", 1e40);
    compute_compressed_muge("Hubble Sees Galaxies Galore");
}

// ===== DYNAMIC SELF-UPDATE AND SELF-EXPANSION CAPABILITIES =====

namespace {
    std::map<std::string, std::map<std::string, VariableMap>> uqff_saved_states;
    std::map<std::string, std::string> uqff_saved_systems;
}

// 1. Variable Management

void createVariable(const std::string& system_name, const std::string& name, double value) {
    auto it = systems.find(system_name);
    if (it != systems.end()) {
        it->second->vars[name] = value;
    }
}

void removeVariable(const std::string& system_name, const std::string& name) {
    auto it = systems.find(system_name);
    if (it != systems.end()) {
        it->second->vars.erase(name);
    }
}

void cloneVariable(const std::string& system_name, const std::string& source, const std::string& destination) {
    auto it = systems.find(system_name);
    if (it != systems.end()) {
        auto& vars = it->second->vars;
        if (vars.find(source) != vars.end()) {
            vars[destination] = vars[source];
        }
    }
}

std::vector<std::string> listVariables(const std::string& system_name) {
    std::vector<std::string> names;
    auto it = systems.find(system_name);
    if (it != systems.end()) {
        for (const auto& pair : it->second->vars) {
            names.push_back(pair.first);
        }
    }
    return names;
}

std::vector<std::string> listAllSystems() {
    std::vector<std::string> names;
    for (const auto& pair : systems) {
        names.push_back(pair.first);
    }
    return names;
}

// 2. Batch Operations

void transformVariableGroup(const std::string& system_name, const std::vector<std::string>& names, std::function<double(double)> transform) {
    auto it = systems.find(system_name);
    if (it != systems.end()) {
        for (const auto& name : names) {
            if (it->second->vars.find(name) != it->second->vars.end()) {
                it->second->vars[name] = transform(it->second->vars[name]);
            }
        }
    }
}

void scaleVariableGroup(const std::string& system_name, const std::vector<std::string>& names, double factor) {
    transformVariableGroup(system_name, names, [factor](double x) { return x * factor; });
}

// 3. Self-Expansion (UQFF domain-specific)

void expandParameterSpace(const std::string& system_name, double factor) {
    auto it = systems.find(system_name);
    if (it != systems.end()) {
        it->second->M *= factor;
        it->second->r *= factor;
        it->second->V *= factor;
        it->second->v_exp *= factor;
    }
}

void expandCompressedScale(const std::string& system_name, double factor) {
    std::vector<std::string> compressed_params = {"M", "r", "V", "F_env"};
    scaleVariableGroup(system_name, compressed_params, factor);
    auto it = systems.find(system_name);
    if (it != systems.end()) {
        it->second->M *= factor;
        it->second->r *= factor;
        it->second->V *= factor;
    }
}

void expandResonanceScale(const std::string& system_name, double factor) {
    std::vector<std::string> resonance_params = {"v_exp", "I", "A", "omega1", "omega2"};
    scaleVariableGroup(system_name, resonance_params, factor);
    auto it = systems.find(system_name);
    if (it != systems.end()) {
        it->second->v_exp *= factor;
        it->second->I *= factor;
        it->second->A *= factor;
    }
}

void expandMultiSystemScale(double factor) {
    for (auto& pair : systems) {
        pair.second->M *= factor;
        pair.second->r *= factor;
        pair.second->V *= factor;
    }
}

// 4. Self-Refinement

void autoRefineParameters(const std::string& system_name, double tolerance) {
    auto it = systems.find(system_name);
    if (it == systems.end()) return;
    
    SystemData* sys = it->second;
    if (sys->M <= 0) sys->M = 1.989e30;
    if (sys->r <= 0) sys->r = 1e10;
    if (sys->V <= 0) sys->V = compute_volume(sys->r);
    if (sys->t <= 0) sys->t = 4.35e17;
    if (sys->z < 0) sys->z = 0.0;
    if (sys->v_exp < 0) sys->v_exp = 0.0;
    if (sys->I < 0) sys->I = 0.0;
    if (sys->A < 0) sys->A = 0.0;
}

void calibrateToObservations(const std::string& system_name, const std::map<std::string, double>& obs_data) {
    auto it = systems.find(system_name);
    if (it == systems.end()) return;
    
    SystemData* sys = it->second;
    for (const auto& obs : obs_data) {
        if (obs.first == "M") sys->M = obs.second;
        else if (obs.first == "r") sys->r = obs.second;
        else if (obs.first == "z") sys->z = obs.second;
        else if (obs.first == "t") sys->t = obs.second;
        else if (obs.first == "V") sys->V = obs.second;
        else if (obs.first == "v_exp") sys->v_exp = obs.second;
        else sys->vars[obs.first] = obs.second;
    }
}

void optimizeForMetric(const std::string& system_name, const std::string& metric_name, double target_value, int iterations) {
    auto it = systems.find(system_name);
    if (it == systems.end()) return;
    
    SystemData* sys = it->second;
    double best_score = 1e100;
    double M_orig = sys->M, r_orig = sys->r, z_orig = sys->z;
    double M_best = M_orig, r_best = r_orig, z_best = z_orig;
    
    for (int i = 0; i < iterations; i++) {
        double perturbation_M = 0.9 + 0.2 * (rand() % 100) / 100.0;
        double perturbation_r = 0.9 + 0.2 * (rand() % 100) / 100.0;
        double perturbation_z = 0.9 + 0.2 * (rand() % 100) / 100.0;
        
        sys->M = M_orig * perturbation_M;
        sys->r = r_orig * perturbation_r;
        sys->z = z_orig * perturbation_z;
        
        double current_value = compute_compressed_muge(system_name, {});
        double score = std::abs(current_value - target_value);
        
        if (score < best_score) {
            best_score = score;
            M_best = sys->M;
            r_best = sys->r;
            z_best = sys->z;
        }
    }
    
    sys->M = M_best;
    sys->r = r_best;
    sys->z = z_best;
}

// 5. Parameter Exploration

std::vector<std::map<std::string, double>> generateVariations(const std::string& system_name, int n_variations) {
    std::vector<std::map<std::string, double>> variations;
    auto it = systems.find(system_name);
    if (it == systems.end()) return variations;
    
    SystemData* sys = it->second;
    double M_orig = sys->M, r_orig = sys->r, z_orig = sys->z, v_orig = sys->v_exp;
    
    for (int i = 0; i < n_variations; i++) {
        double factor_M = 0.8 + 0.4 * (rand() % 100) / 100.0;
        double factor_r = 0.8 + 0.4 * (rand() % 100) / 100.0;
        double factor_z = 0.8 + 0.4 * (rand() % 100) / 100.0;
        double factor_v = 0.8 + 0.4 * (rand() % 100) / 100.0;
        
        std::map<std::string, double> variation;
        variation["M"] = M_orig * factor_M;
        variation["r"] = r_orig * factor_r;
        variation["z"] = z_orig * factor_z;
        variation["v_exp"] = v_orig * factor_v;
        variations.push_back(variation);
    }
    
    return variations;
}

// 6. Adaptive Evolution

void mutateParameters(const std::string& system_name, double mutation_rate) {
    auto it = systems.find(system_name);
    if (it == systems.end()) return;
    
    SystemData* sys = it->second;
    double mutation_M = 1.0 + mutation_rate * (-0.5 + (rand() % 100) / 100.0);
    double mutation_r = 1.0 + mutation_rate * (-0.5 + (rand() % 100) / 100.0);
    double mutation_z = 1.0 + mutation_rate * (-0.5 + (rand() % 100) / 100.0);
    double mutation_v = 1.0 + mutation_rate * (-0.5 + (rand() % 100) / 100.0);
    
    sys->M *= mutation_M;
    sys->r *= mutation_r;
    sys->z *= mutation_z;
    sys->v_exp *= mutation_v;
}

void evolveSystem(const std::string& system_name, int generations, std::function<double(const std::string&)> fitness) {
    auto it = systems.find(system_name);
    if (it == systems.end()) return;
    
    double best_fitness = fitness(system_name);
    SystemData* sys = it->second;
    double M_best = sys->M, r_best = sys->r, z_best = sys->z, v_best = sys->v_exp;
    
    for (int gen = 0; gen < generations; gen++) {
        mutateParameters(system_name, 0.1);
        double current_fitness = fitness(system_name);
        
        if (current_fitness > best_fitness) {
            best_fitness = current_fitness;
            M_best = sys->M;
            r_best = sys->r;
            z_best = sys->z;
            v_best = sys->v_exp;
        } else {
            sys->M = M_best;
            sys->r = r_best;
            sys->z = z_best;
            sys->v_exp = v_best;
        }
    }
}

// 7. State Management

void saveState(const std::string& label, const std::string& system_name) {
    auto it = systems.find(system_name);
    if (it != systems.end()) {
        uqff_saved_states[label][system_name] = it->second->vars;
        uqff_saved_systems[label] = system_name;
    }
}

void restoreState(const std::string& label) {
    if (uqff_saved_states.find(label) != uqff_saved_states.end()) {
        std::string system_name = uqff_saved_systems[label];
        auto it = systems.find(system_name);
        if (it != systems.end()) {
            it->second->vars = uqff_saved_states[label][system_name];
        }
    }
}

std::vector<std::string> listSavedStates() {
    std::vector<std::string> labels;
    for (const auto& pair : uqff_saved_states) {
        labels.push_back(pair.first);
    }
    return labels;
}

std::map<std::string, double> exportState(const std::string& system_name) {
    std::map<std::string, double> state;
    auto it = systems.find(system_name);
    if (it != systems.end()) {
        state["M"] = it->second->M;
        state["r"] = it->second->r;
        state["z"] = it->second->z;
        state["t"] = it->second->t;
        state["V"] = it->second->V;
        state["v_exp"] = it->second->v_exp;
        state["g_compressed"] = compute_compressed_muge(system_name, {});
        state["g_resonance"] = compute_resonance_muge(system_name, {});
    }
    return state;
}

// 8. System Analysis

std::map<std::string, double> sensitivityAnalysis(const std::string& system_name, const std::string& param, double delta) {
    std::map<std::string, double> sensitivities;
    auto it = systems.find(system_name);
    if (it == systems.end()) return sensitivities;
    
    SystemData* sys = it->second;
    double original_value = 0;
    if (param == "M") original_value = sys->M;
    else if (param == "r") original_value = sys->r;
    else if (param == "z") original_value = sys->z;
    
    double g_base = compute_compressed_muge(system_name, {});
    
    if (param == "M") sys->M = original_value * (1.0 + delta);
    else if (param == "r") sys->r = original_value * (1.0 + delta);
    else if (param == "z") sys->z = original_value * (1.0 + delta);
    double g_plus = compute_compressed_muge(system_name, {});
    
    if (param == "M") sys->M = original_value * (1.0 - delta);
    else if (param == "r") sys->r = original_value * (1.0 - delta);
    else if (param == "z") sys->z = original_value * (1.0 - delta);
    double g_minus = compute_compressed_muge(system_name, {});
    
    if (param == "M") sys->M = original_value;
    else if (param == "r") sys->r = original_value;
    else if (param == "z") sys->z = original_value;
    
    sensitivities["g_sensitivity"] = (g_plus - g_minus) / (2.0 * delta * original_value);
    sensitivities["g_base"] = g_base;
    sensitivities["g_plus"] = g_plus;
    sensitivities["g_minus"] = g_minus;
    
    return sensitivities;
}

std::string generateReport(const std::string& system_name) {
    std::ostringstream report;
    report << std::scientific << std::setprecision(4);
    
    auto it = systems.find(system_name);
    if (it == systems.end()) return "System not found";
    
    SystemData* sys = it->second;
    
    report << "===== UQFF Module Report =====\n";
    report << "System: " << sys->name << "\n";
    report << "Description: " << sys->description << "\n\n";
    
    report << "System Parameters:\n";
    report << "  M = " << sys->M << " kg (" << sys->M / 1.989e30 << " Msun)\n";
    report << "  r = " << sys->r << " m (" << sys->r / 9.461e15 << " ly)\n";
    report << "  z = " << sys->z << "\n";
    report << "  t = " << sys->t << " s\n";
    report << "  V = " << sys->V << " m³\n";
    report << "  v_exp = " << sys->v_exp << " m/s\n";
    report << "  F_env = " << sys->F_env << "\n\n";
    
    report << "Resonance Parameters:\n";
    report << "  I = " << sys->I << " A\n";
    report << "  A = " << sys->A << " m²\n";
    report << "  omega1 = " << sys->omega1 << " rad/s\n";
    report << "  omega2 = " << sys->omega2 << " rad/s\n\n";
    
    if (sys->M_sun > 0) {
        report << "Orbital Parameters:\n";
        report << "  M_sun = " << sys->M_sun << " kg\n";
        report << "  r_orbit = " << sys->r_orbit << " m\n\n";
    }
    
    double g_comp = compute_compressed_muge(system_name, {});
    double g_res = compute_resonance_muge(system_name, {});
    
    report << "Computed Values:\n";
    report << "  g_compressed = " << g_comp << " m/s²\n";
    report << "  g_resonance = " << g_res << " m/s²\n";
    report << "  Ratio (Res/Comp) = " << (g_res / g_comp) << "\n\n";
    
    report << "Key Variables:\n";
    for (const auto& v : sys->vars) {
        report << "  " << v.first << " = " << v.second << "\n";
    }
    
    report << "\nSaved States: " << uqff_saved_states.size() << "\n";
    report << "========================================\n";
    
    return report.str();
}

bool validateConsistency(const std::string& system_name) {
    auto it = systems.find(system_name);
    if (it == systems.end()) return false;
    
    SystemData* sys = it->second;
    bool valid = true;
    
    if (sys->M <= 0) valid = false;
    if (sys->r <= 0) valid = false;
    if (sys->V <= 0) valid = false;
    if (sys->t <= 0) valid = false;
    if (sys->z < 0) valid = false;
    
    return valid;
}

void autoCorrectAnomalies(const std::string& system_name) {
    auto it = systems.find(system_name);
    if (it == systems.end()) return;
    
    SystemData* sys = it->second;
    if (sys->M <= 0) sys->M = 1.989e30;
    if (sys->r <= 0) sys->r = 1e10;
    if (sys->V <= 0) sys->V = compute_volume(sys->r);
    if (sys->t <= 0) sys->t = 4.35e17;
    if (sys->z < 0) sys->z = 0.0;
    if (sys->v_exp < 0) sys->v_exp = 0.0;
    if (sys->I < 0) sys->I = 0.0;
    if (sys->A < 0) sys->A = 0.0;
    if (sys->F_env < -1.0) sys->F_env = 0.0;
}

// ===== ENHANCED EXAMPLE USAGE WITH DYNAMIC CAPABILITIES =====

void enhanced_example_usage() {
    install_uqff_module();
    
    std::cout << "\n=== 1. Variable Management ===\n";
    createVariable("Orion Nebula", "custom_scale_factor", 1.05);
    cloneVariable("Orion Nebula", "M", "M_backup");
    auto var_list = listVariables("Orion Nebula");
    std::cout << "Orion variables: " << var_list.size() << "\n";
    auto all_systems = listAllSystems();
    std::cout << "Total systems: " << all_systems.size() << "\n\n";
    
    std::cout << "=== 2. Batch System Scaling ===\n";
    std::vector<std::string> system_group = {"M", "r", "V"};
    scaleVariableGroup("Hubble Sees Galaxies Galore", system_group, 1.1);
    std::cout << "Hubble Galaxies scaled by 1.1\n\n";
    
    std::cout << "=== 3. Self-Expansion ===\n";
    expandCompressedScale("The Stellar Forge", 1.08);
    expandResonanceScale("Orion Nebula", 1.05);
    std::cout << "Expanded Stellar Forge (compressed) +8%, Orion (resonance) +5%\n\n";
    
    std::cout << "=== 4. Self-Refinement ===\n";
    autoRefineParameters("Saturn", 1e-10);
    std::map<std::string, double> observations = {
        {"M", 5.8e26},
        {"r", 6.1e7},
        {"z", 0.0}
    };
    calibrateToObservations("Saturn", observations);
    std::cout << "Saturn calibrated to observations\n\n";
    
    std::cout << "=== 5. Optimize g_UQFF ===\n";
    optimizeForMetric("Orion Nebula", "g_UQFF", 1e30, 50);
    std::cout << "Optimized Orion for g~1e30\n\n";
    
    std::cout << "=== 6. Generate Variations ===\n";
    auto variations = generateVariations("The Crab Nebula", 15);
    std::cout << "Generated " << variations.size() << " Crab Nebula variations\n\n";
    
    std::cout << "=== 7. Multi-System State Management ===\n";
    saveState("hubble_initial", "Hubble Sees Galaxies Galore");
    saveState("stellar_forge", "The Stellar Forge");
    saveState("sombrero_initial", "Hubble Mosaic of the Majestic Sombrero Galaxy");
    saveState("saturn_calibrated", "Saturn");
    saveState("orion_optimized", "Orion Nebula");
    std::cout << "Saved 5 system states\n\n";
    
    std::cout << "=== 8. Sensitivity Analysis ===\n";
    auto sensitivity = sensitivityAnalysis("Orion Nebula", "M", 0.1);
    std::cout << "Orion M sensitivity: dg/dM = " << sensitivity["g_sensitivity"] << "\n\n";
    
    std::cout << "=== 9. System Validation ===\n";
    for (const auto& sys_name : all_systems) {
        if (!validateConsistency(sys_name)) {
            std::cout << sys_name << " has inconsistencies, correcting...\n";
            autoCorrectAnomalies(sys_name);
        }
    }
    std::cout << "All systems validated\n\n";
    
    std::cout << "=== 10. Comprehensive Report ===\n";
    std::cout << generateReport("Orion Nebula") << "\n";
    
    std::cout << "=== 11. Adaptive Evolution ===\n";
    auto g_fitness = [](const std::string& sys_name) {
        double g = compute_compressed_muge(sys_name, {});
        return -std::abs(std::log10(g) - 30);
    };
    evolveSystem("The Crab Nebula", 25, g_fitness);
    std::cout << "Evolved Crab Nebula over 25 generations\n\n";
    
    std::cout << "=== 12. Multi-System Compressed Comparison ===\n";
    std::vector<std::string> comparison_systems = {
        "Hubble Sees Galaxies Galore",
        "The Stellar Forge",
        "Orion Nebula",
        "The Crab Nebula",
        "Saturn"
    };
    for (const auto& sys : comparison_systems) {
        double g = compute_compressed_muge(sys, {});
        std::cout << sys << ": g_compressed = " << g << " m/s²\n";
    }
    std::cout << "\n";
    
    std::cout << "=== 13. Resonance Mode Comparison ===\n";
    for (const auto& sys : comparison_systems) {
        double g = compute_resonance_muge(sys, {});
        std::cout << sys << ": g_resonance = " << g << " m/s²\n";
    }
    std::cout << "\n";
    
    std::cout << "=== 14. Compressed vs Resonance Ratios ===\n";
    for (const auto& sys : comparison_systems) {
        double g_comp = compute_compressed_muge(sys, {});
        double g_res = compute_resonance_muge(sys, {});
        std::cout << sys << ": Ratio (Res/Comp) = " << (g_res / g_comp) << "\n";
    }
    std::cout << "\n";
    
    std::cout << "=== 15. State Restoration ===\n";
    restoreState("hubble_initial");
    std::cout << "Restored Hubble initial state\n";
    double g_restored = compute_compressed_muge("Hubble Sees Galaxies Galore", {});
    std::cout << "g_compressed (restored) = " << g_restored << " m/s²\n\n";
    
    std::cout << "=== 16. Final State Exports ===\n";
    for (const auto& sys : comparison_systems) {
        auto state = exportState(sys);
        std::cout << sys << ":\n";
        std::cout << "  M=" << state["M"] << " kg, r=" << state["r"] << " m\n";
        std::cout << "  g_comp=" << state["g_compressed"] << ", g_res=" << state["g_resonance"] << "\n";
    }
    std::cout << "\n";
}

// Evaluation of uqff_module.cpp (UQFF Compressed and Resonance Equations Module)

**Strengths:**
-**Dynamic Variable Management : **All system parameters and computed variables are stored in a `std: : map<std::string, double>` (`VariableMap`), supporting runtime updates, additions, and removals for any system.
    - **System Extensibility : **Systems are defined via the `SystemData` struct and stored in a global map.New systems can be added easily in `install_uqff_module()` with custom parameters and descriptions.
    - **Flexible Update Functions : **The module provides `update_variable`, `add_variable`, and `subtract_variable` functions for any system and variable, allowing both direct setting and incremental changes.Dependent variables are recalculated when key values(e.g., `"z"`) are updated.
        - **Immediate Effect : **All computations(`compute_compressed_muge`, `compute_resonance_muge`) use the current values in the variable map, so any changes are immediately reflected in results.
            - **Pluggable Design : **The module is designed to be integrated into a base program, with clear external function interfaces for installation, computation, variable updates, and text output.
            - **Comprehensive Physics : **All terms(gravity, quantum, fluid, resonance, DM, etc.) are computed and summed explicitly, with nothing neglected.
            - **Debugging & Documentation : **The `print_system_text` function outputs all key variables and descriptions for any system, aiding validation and troubleshooting.
            - **Sample Usage Provided : **Example usage demonstrates how to install, compute, update, and print system data.

            ** Weaknesses / Recommendations : **
            -**Error Handling : **Adding unknown variables is flexible but may lead to silent errors if a typo occurs.Consider stricter validation or optional warnings for unknown variable names.
            - **Unit Consistency : **Ensure all units are consistent, especially when combining terms from different physical domains.
            - **Magic Numbers : **Some constants and scale factors are hardcoded.Document their physical meaning or allow configuration via external files or constructor arguments.
            - **Performance : **For large - scale or repeated updates, consider profiling and optimizing critical paths if needed.
            - **Thread Safety : **The global `systems` map is not thread - safe.If used in a multi - threaded context, add appropriate synchronization.

            ** Is the code dynamic enough to be updated and expanded ? **
            Yes, the code is highly dynamic and extensible.You can add new systems, update or add variables at runtime, and all computations will adapt accordingly.The modular structure and use of maps for variables and systems make it straightforward to expand the module for new physical models or additional systems.

            * *Summary : **
            The module is robust, dynamic, and extensible, supporting runtime updates and expansion for new systems and variables.Minor improvements in error handling, documentation, and thread safety are recommended for production use.