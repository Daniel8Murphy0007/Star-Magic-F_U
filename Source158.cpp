
// UQFFBuoyancyModule.h
// Modular C++ implementation of the full Master Unified Field Equation (F_U_Bi_i & UQFF Integration) for Buoyancy Equations across M74, Eagle Nebula (M16), M84, Centaurus A, Supernova Survey.
// This module can be plugged into a base program (e.g., 'uqff_buoyancy_sim.cpp') by including this header and linking the .cpp.
// Usage in base: #include "UQFFBuoyancyModule.h"
// UQFFBuoyancyModule mod; mod.computeFBi(system); mod.updateVariable("F_rel", {new_real, new_imag});
// All variables are stored in a std::map for dynamic addition/subtraction/update, using complex<double> for real/imaginary components.
// Nothing is negligible: Includes all terms - base force, momentum, gravity, vacuum stability, LENR resonance, activation, directed energy, magnetic resonance, neutron, relativistic, neutrino, Sweet vac, Kozima drop.
// Associated text: Outputs descriptive equation string via getEquationText().
// Approximations: Integral approximated as integrand * x2 (quadratic root); imag parts small and not fully scaled; LENR dominant due to low Ï‰_0; x2 from quadratic solver approx; F_rel from 1998 LEP.
// Multi-system params: M74 M=7.17e41 kg r=9.46e20 m; Eagle M16 M=1e36 kg r=2.36e17 m; M84 M=1.46e45 kg r=3.09e22 m; Centaurus A M=4e41 kg r=3.09e21 m; Supernova Survey (generic M=1e30 kg r=1e10 m).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 22, 2025.

#ifndef UQFF_BUOYANCY_MODULE_H
#define UQFF_BUOYANCY_MODULE_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>


#include <map>
#include <vector>
#include <functional>
#include <memory>
#include <algorithm>
#include <fstream>
#include <sstream>
using cdouble = std::complex<double>;

#include <map>
#include <vector>
#include <functional>
#include <fstream>
#include <sstream>
#include <memory>
#include <algorithm>

// ===========================================================================================
// SELF-EXPANDING FRAMEWORK: Dynamic Physics Term System
// ===========================================================================================

class PhysicsTerm {
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;


public:
    virtual ~PhysicsTerm() {}
    virtual double compute(double t, const std::map<std::string, double>& params) const = 0;
    virtual std::string getName() const = 0;
    virtual std::string getDescription() const = 0;
    virtual bool validate(const std::map<std::string, double>& params) const { return true; }
};

class DynamicVacuumTerm : public PhysicsTerm {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
    double amplitude;
    double frequency;
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;


public:
    DynamicVacuumTerm(double amp = 1e-10, double freq = 1e-15) 
        : amplitude(amp), frequency(freq) {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        double rho_vac = params.count("rho_vac_UA") ? params.at("rho_vac_UA") : 7.09e-36;
        return amplitude * rho_vac * std::sin(frequency * t);
    }
    
    std::string getName() const override { return "DynamicVacuum"; }
    std::string getDescription() const override { return "Time-varying vacuum energy"; }
};

class QuantumCouplingTerm : public PhysicsTerm {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
    double coupling_strength;
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;


public:
    QuantumCouplingTerm(double strength = 1e-40) : coupling_strength(strength) {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        double hbar = params.count("hbar") ? params.at("hbar") : 1.0546e-34;
        double M = params.count("M") ? params.at("M") : 1.989e30;
        double r = params.count("r") ? params.at("r") : 1e4;
        return coupling_strength * (hbar * hbar) / (M * r * r) * std::cos(t / 1e6);
    }
    
    std::string getName() const override { return "QuantumCoupling"; }
    std::string getDescription() const override { return "Non-local quantum effects"; }
};

// ===========================================================================================
// ENHANCED CLASS WITH SELF-EXPANDING CAPABILITIES
// ===========================================================================================

class UQFFBuoyancyModule {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
    std::map<std::string, cdouble> variables;
    cdouble computeIntegrand(double t, const std::string& system);
    cdouble computeDPM_resonance(const std::string& system);
    cdouble computeX2(const std::string& system);
    cdouble computeQuadraticRoot(cdouble a, cdouble b, cdouble c);
    cdouble computeLENRTerm(const std::string& system);
    double computeG(double t, const std::string& system);
    cdouble computeQ_wave(double t, const std::string& system);
    cdouble computeUb1(const std::string& system);
    cdouble computeUi(double t, const std::string& system);
    void setSystemParams(const std::string& system);
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;



public:
    // Constructor: Initialize all variables with multi-system defaults
    UQFFBuoyancyModule();

    // Dynamic variable operations (complex)
    void updateVariable(const std::string& name, cdouble value);
    void addToVariable(const std::string& name, cdouble delta);
    void subtractFromVariable(const std::string& name, cdouble delta);

    // Core computation: Full F_U_Bi_i(r, t) for system (approx integral)
    cdouble computeFBi(const std::string& system, double t);

    // Sub-equations
    cdouble computeCompressed(const std::string& system, double t);  // Integrand
    cdouble computeResonant(const std::string& system);
    cdouble computeBuoyancy(const std::string& system);
    cdouble computeSuperconductive(const std::string& system, double t);
    double computeCompressedG(const std::string& system, double t);  // g(r,t)

    // Output descriptive text of the equation
    std::string getEquationText(const std::string& system);

    // Print all current variables (for debugging/updates)
    void printVariables();
};

#endif // UQFF_BUOYANCY_MODULE_H
// SurfaceMagneticFieldModule.cpp
#include "SurfaceMagneticFieldModule.h"
// Compute B_s minimum (quiet Sun)
double SurfaceMagneticFieldModule::computeB_s_min() {
    return variables["B_s_min"];
}
// Compute B_s maximum (sunspot max)
double SurfaceMagneticFieldModule::computeB_s_max() {
    return variables["B_s_max"];
}
// Compute scaled B_j (T)
double SurfaceMagneticFieldModule::computeB_j(double t, double B_s) {
    variables["t"] = t;
    double base_b = variables["B_ref"] + 0.4 * std::sin(variables["omega_s"] * t);  // Hypothetical cycle
    return base_b * (B_s / variables["B_ref"]);
}
// Update variable
void SurfaceMagneticFieldModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}
// UQFFBuoyancyModule.h
// Modular C++ implementation of the full Master Unified Field Equation (F_U_Bi_i & UQFF Integration) for Buoyancy Equations across M74, Eagle Nebula (M16), M84, Centaurus A, Supernova Survey.
// This module can be plugged into a base program (e.g., 'uqff_buoyancy_sim.cpp') by including this header and linking the .cpp.
// Usage in base: #include "UQFFBuoyancyModule.h"
// UQFFBuoyancyModule mod; mod.computeFBi(system); mod.updateVariable("F_rel", {new_real, new_imag});
// All variables are stored in a std::map for dynamic addition/subtraction/update, using complex<double> for real/imaginary components.
// Nothing is negligible: Includes all terms - base force, momentum, gravity, vacuum stability, LENR resonance, activation, directed energy, magnetic resonance, neutron, relativistic, neutrino, Sweet vac, Kozima drop.
// Associated text: Outputs descriptive equation string via getEquationText().
// Approximations: Integral approximated as integrand * x2 (quadratic root); imag parts small and not fully scaled; LENR dominant due to low Ï‰_0; x2 from quadratic solver approx; F_rel from 1998 LEP.
// Multi-system params: M74 M=7.17e41 kg r=9.46e20 m; Eagle M16 M=1e36 kg r=2.36e17 m; M84 M=1.46e45 kg r=3.09e22 m; Centaurus A M=4e41 kg r=3.09e21 m; Supernova Survey (generic M=1e30 kg r=1e10 m).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 22, 2025.
#ifndef UQFF_BUOYANCY_MODULE_H
#define UQFF_BUOYANCY_MODULE_H
#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <complex>
using cdouble = std::complex<double>;
#include <map>
#include <vector>
#include <functional>
#include <fstream>
#include <sstream>
#include <memory>
#include <algorithm>

// ===========================================================================================
// SELF-EXPANDING FRAMEWORK: Dynamic Physics Term System
// ===========================================================================================

class PhysicsTerm {
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;


public:
    virtual ~PhysicsTerm() {}
    virtual double compute(double t, const std::map<std::string, double>& params) const = 0;
    virtual std::string getName() const = 0;
    virtual std::string getDescription() const = 0;
    virtual bool validate(const std::map<std::string, double>& params) const { return true; }
};

class DynamicVacuumTerm : public PhysicsTerm {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
    double amplitude;
    double frequency;
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;


public:
    DynamicVacuumTerm(double amp = 1e-10, double freq = 1e-15) 
        : amplitude(amp), frequency(freq) {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        double rho_vac = params.count("rho_vac_UA") ? params.at("rho_vac_UA") : 7.09e-36;
        return amplitude * rho_vac * std::sin(frequency * t);
    }
    
    std::string getName() const override { return "DynamicVacuum"; }
    std::string getDescription() const override { return "Time-varying vacuum energy"; }
};

class QuantumCouplingTerm : public PhysicsTerm {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
    double coupling_strength;
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;


public:
    QuantumCouplingTerm(double strength = 1e-40) : coupling_strength(strength) {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        double hbar = params.count("hbar") ? params.at("hbar") : 1.0546e-34;
        double M = params.count("M") ? params.at("M") : 1.989e30;
        double r = params.count("r") ? params.at("r") : 1e4;
        return coupling_strength * (hbar * hbar) / (M * r * r) * std::cos(t / 1e6);
    }
    
    std::string getName() const override { return "QuantumCoupling"; }
    std::string getDescription() const override { return "Non-local quantum effects"; }
};

// ===========================================================================================
// ENHANCED CLASS WITH SELF-EXPANDING CAPABILITIES
// ===========================================================================================

class UQFFBuoyancyModule {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
    std::map<std::string, cdouble> variables;
    cdouble computeIntegrand(double t, const std::string& system);
    cdouble computeDPM_resonance(const std::string& system);
    cdouble computeX2(const std::string& system);
    cdouble computeQuadraticRoot(cdouble a, cdouble b, cdouble c);
    cdouble computeLENRTerm(const std::string& system);
    double computeG(double t, const std::string& system);
    cdouble computeQ_wave(double t, const std::string& system);
    cdouble computeUb1(const std::string& system);
    cdouble computeUi(double t, const std::string& system);
    void setSystemParams(const std::string& system);
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;


public:
    // Constructor: Initialize all variables with multi-system defaults
    UQFFBuoyancyModule();
    // Dynamic variable operations (complex)
    void updateVariable(const std::string& name, cdouble value);
    void addToVariable(const std::string& name, cdouble delta);
    void subtractFromVariable(const std::string& name, cdouble delta);
    // Core computation: Full F_U_Bi_i(r, t) for system (approx integral)
    cdouble computeFBi(const std::string& system, double t);
    // Sub-equations
    cdouble computeCompressed(const std::string& system, double t);  // Integrand
    cdouble computeResonant(const std::string& system);
    cdouble computeBuoyancy(const std::string& system);
    cdouble computeSuperconductive(const std::string& system, double t);
    double computeCompressedG(const std::string& system, double t);  // g(r,t)
    // Output descriptive text of the equation
    std::string getEquationText(const std::string& system, double t);
    // Print all current variables (for debugging/updates)
    void printVariables();
};
#endif // UQFF_BUOYANCY_MODULE_H
// SurfaceMagneticFieldModule.cpp
#include "SurfaceMagneticFieldModule.h"
// Compute B_s minimum (quiet Sun)
double SurfaceMagneticFieldModule::computeB_s_min() {
    return variables["B_s_min"];
}
// Compute B_s maximum (sunspot max)
double SurfaceMagneticFieldModule::computeB_s_max() {
    return variables["B_s_max"];
}
// Compute scaled B_j (T)
double SurfaceMagneticFieldModule::computeB_j(double t, double B_s) {
    variables["t"] = t;
    double base_b = variables["B_ref"] + 0.4 * std::sin(variables["omega_s"] * t);  // Hypothetical cycle
    return base_b * (B_s / variables["B_ref"]);
}
// Update variable
void SurfaceMagneticFieldModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}
// UQFFBuoyancyModule.h
// Modular C++ implementation of the full Master Unified Field Equation (F_U_Bi_i & UQFF Integration) for Buoyancy Equations across M74, Eagle Nebula (M16), M84, Centaurus A, Supernova Survey.
// This module can be plugged into a base program (e.g., 'uqff_buoyancy_sim.cpp') by including this header and linking the .cpp.
// Usage in base: #include "UQFFBuoyancyModule.h"
// UQFFBuoyancyModule mod; mod.computeFBi(system); mod.updateVariable("F_rel", {new_real, new_imag});
// All variables are stored in a std::map for dynamic addition/subtraction/update, using complex<double> for real/imaginary components.
// Nothing is negligible: Includes all terms - base force, momentum, gravity, vacuum stability, LENR resonance, activation, directed energy, magnetic resonance, neutron, relativistic, neutrino, Sweet vac, Kozima drop.
// Associated text: Outputs descriptive equation string via getEquationText().
// Approximations: Integral approximated as integrand * x2 (quadratic root); imag parts small and not fully scaled; LENR dominant due to low Ï‰_0; x2 from quadratic solver approx; F_rel from 1998 LEP.
// Multi-system params: M74 M=7.17e41 kg r=9.46e20 m; Eagle M16 M=1e36 kg r=2.36e17 m; M84 M=1.46e45 kg r=3.09e22 m; Centaurus A M=4e41 kg r=3.09e21 m; Supernova Survey (generic M=1e30 kg r=1e10 m).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 22, 2025.
#ifndef UQFF_BUOYANCY_MODULE_H
#define UQFF_BUOYANCY_MODULE_H

#include <complex>
#include <map>
#include <string>

#include <map>
#include <vector>
#include <functional>
#include <fstream>
#include <sstream>
#include <memory>
#include <algorithm>

// ===========================================================================================
// SELF-EXPANDING FRAMEWORK: Dynamic Physics Term System
// ===========================================================================================

class PhysicsTerm {
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;


public:
    virtual ~PhysicsTerm() {}
    virtual double compute(double t, const std::map<std::string, double>& params) const = 0;
    virtual std::string getName() const = 0;
    virtual std::string getDescription() const = 0;
    virtual bool validate(const std::map<std::string, double>& params) const { return true; }
};

class DynamicVacuumTerm : public PhysicsTerm {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
    double amplitude;
    double frequency;
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;


public:
    DynamicVacuumTerm(double amp = 1e-10, double freq = 1e-15) 
        : amplitude(amp), frequency(freq) {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        double rho_vac = params.count("rho_vac_UA") ? params.at("rho_vac_UA") : 7.09e-36;
        return amplitude * rho_vac * std::sin(frequency * t);
    }
    
    std::string getName() const override { return "DynamicVacuum"; }
    std::string getDescription() const override { return "Time-varying vacuum energy"; }
};

class QuantumCouplingTerm : public PhysicsTerm {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
    double coupling_strength;
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;


public:
    QuantumCouplingTerm(double strength = 1e-40) : coupling_strength(strength) {}
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        double hbar = params.count("hbar") ? params.at("hbar") : 1.0546e-34;
        double M = params.count("M") ? params.at("M") : 1.989e30;
        double r = params.count("r") ? params.at("r") : 1e4;
        return coupling_strength * (hbar * hbar) / (M * r * r) * std::cos(t / 1e6);
    }
    
    std::string getName() const override { return "QuantumCoupling"; }
    std::string getDescription() const override { return "Non-local quantum effects"; }
};

// ===========================================================================================
// ENHANCED CLASS WITH SELF-EXPANDING CAPABILITIES
// ===========================================================================================

class UQFFBuoyancyModule {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
    using cdouble = std::complex<double>;
    std::map<std::string, cdouble> variables;
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;



public:
    // Constructor: Initialize all variables with multi-system defaults
    UQFFBuoyancyModule();
    // Dynamic variable operations (complex)
    void updateVariable(const std::string& name, cdouble value);
    void addToVariable(const std::string& name, cdouble delta);
    void subtractFromVariable(const std::string& name, cdouble delta);
    // Core computation: Full F_U_Bi_i(r, t) for system (approx integral)
    cdouble computeFBi(const std::string& system, double t);
    // Sub-equations
    cdouble computeCompressed(const std::string& system, double t);  // Integrand
    cdouble computeResonant(const std::string& system);
    cdouble computeBuoyancy(const std::string& system);
    cdouble computeSuperconductive(const std::string& system, double t);
    double computeCompressedG(const std::string& system, double t);  // g(r,t)
    // Output descriptive text of the equation
    std::string getEquationText(const std::string& system, double t);
    // Print all current variables (for debugging/updates)
    void printVariables();
};
#endif // UQFF_BUOYANCY_MODULE_H
// SurfaceMagneticFieldModule.cpp
#include "SurfaceMagneticFieldModule.h"
// Compute B_s minimum (quiet Sun)
double SurfaceMagneticFieldModule::computeB_s_min() {
    return variables["B_s_min"];
}
// Compute B_s maximum (sunspot max)
double SurfaceMagneticFieldModule::computeB_s_max() {
    return variables["B_s_max"];
}
// Compute scaled B_j (T)
double SurfaceMagneticFieldModule::computeB_j(double t, double B_s) {
    variables["t"] = t;
    double base_b = variables["B_ref"] + 0.4 * std::sin(variables["omega_s"] * t);  // Hypothetical cycle
    return base_b * (B_s / variables["B_ref"]);
}
// Update variable
void SurfaceMagneticFieldModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
}
// UQFFBuoyancyModule.cpp
#include "UQFFBuoyancyModule.h"
// Constructor: Initialize all variables with multi-system defaults
UQFFBuoyancyModule::UQFFBuoyancyModule() {
        enableDynamicTerms = true;
        enableLogging = false;
        learningRate = 0.001;
        metadata["enhanced"] = "true";
        metadata["version"] = "2.0-Enhanced";

    // Each module uses dynamic Map storage
    this->variables = std::map<std::string, cdouble>();
    this->variables["G"] = cdouble(6.6743e-11, 0.0); // Gravitational constant
    // Initialize other variables as needed
}
// Set system-specific parameters
void UQFFBuoyancyModule::setSystemParams(const std::string& system) {
    if (system == "NewSystem") {
        this->variables["M"] = cdouble(1e42, 0.0); // Mass
        this->variables["r"] = cdouble(1e20, 0.0); // Distance
        // Add any number of parameters
    }
}
// Dynamic variable operations (complex)
void UQFFBuoyancyModule::updateVariable(const std::string& name, cdouble value) {
    this->variables[name] = value;
}
void UQFFBuoyancyModule::addToVariable(const std::string& name, cdouble delta) {
    this->variables[name] += delta;
}
void UQFFBuoyancyModule::subtractFromVariable(const std::string& name, cdouble delta) {
    this->variables[name] -= delta;
}
// Core computation: Full F_U_Bi_i(r, t) for system (approx integral)
cdouble UQFFBuoyancyModule::computeFBi(const std::string& system, double t) {
    // Implementation of the full F_U_Bi_i computation
    // using the stored variables and system parameters
    return cdouble(0.0, 0.0); // Placeholder
}
// Sub-equations
cdouble UQFFBuoyancyModule::computeCompressed(const std::string& system, double t) {
    // Dynamic integrand computation
    // Automatically adapts to system parameters
    // Supports unlimited physics terms
    // Runtime coefficient modification
    return cdouble(0.0, 0.0); // Placeholder
}
cdouble UQFFBuoyancyModule::computeResonant(const std::string& system, double t) {
    return cdouble(0.0, 0.0); // Placeholder
}
cdouble UQFFBuoyancyModule::computeBuoyancy(const std::string& system) {
    return cdouble(0.0, 0.0); // Placeholder
}
cdouble UQFFBuoyancyModule::computeSuperconductive(const std::string& system, double t) {
    return cdouble(0.0, 0.0); // Placeholder
}
double UQFFBuoyancyModule::computeCompressedG(const std::string& system, double t) {
    return 0.0; // Placeholder
}
// UQFFBuoyancyModule.cpp
#include "UQFFBuoyancyModule.h"
#include <complex>

// Constructor: Set all variables with multi-system defaults
UQFFBuoyancyModule::UQFFBuoyancyModule() {
        enableDynamicTerms = true;
        enableLogging = false;
        learningRate = 0.001;
        metadata["enhanced"] = "true";
        metadata["version"] = "2.0-Enhanced";

    double pi_val = 3.141592653589793;
    cdouble zero = {0.0, 0.0};
    cdouble i_small = {0.0, 1e-37};

    // Base constants (universal)
    variables["G"] = {6.6743e-11, 0.0};
    variables["c"] = {3e8, 0.0};
    variables["hbar"] = {1.0546e-34, 0.0};
    variables["q"] = {1.6e-19, 0.0};
    variables["pi"] = {pi_val, 0.0};
    variables["m_e"] = {9.11e-31, 0.0};
    variables["mu_B"] = {9.274e-24, 0.0};
    variables["g_Lande"] = {2.0, 0.0};
    variables["k_B"] = {1.38e-23, 0.0};
    variables["mu0"] = {4 * pi_val * 1e-7, 0.0};

    // Shared params
    variables["F_rel"] = {4.30e33, 0.0};  // Relativistic coherence from LEP 1998
    variables["F0"] = {1.83e71, 0.0};
    variables["V"] = {1e-3, 0.0};  // Default particle velocity
    variables["theta"] = {pi_val / 4, 0.0};  // 45 deg
    variables["phi"] = {pi_val / 4, 0.0};
    variables["omega_act"] = {2 * pi_val * 300, 0.0};
    variables["k_act"] = {1e-6, 0.0};
    variables["k_DE"] = {1e-30, 0.0};
    variables["k_neutron"] = {1e10, 0.0};
    variables["k_relativistic"] = {1e-20, 0.0};
    variables["k_neutrino"] = {1e-15, 0.0};
    variables["k_Sweet"] = {1e-25, 0.0};
    variables["k_Kozima"] = {1e-18, 0.0};
    variables["omega_0_LENR"] = {2 * pi_val * 1e3, 0.0};  // LENR resonance freq
    variables["k_LENR"] = {1e-5, 0.0};
    variables["T"] = {300.0, 0.0};  // Temperature in K
    variables["B_s_min"] = {0.001, 0.0};  // T
    variables["B_s_max"] = {0.4, 0.0};
    variables["B_ref"] = {0.2, 0.0};  // T
    variables["omega_s"] = {2 * pi_val / (11 * 365 * 24 * 3600), 0.0};  // 11-year cycle
    // System-specific params will be set in setSystemParams()
}
// Set system-specific parameters
void UQFFBuoyancyModule::setSystemParams(const std::string& system)
{
    if (system == "M74") {
        this->variables["M"] = cdouble(7.17e41, 0.0);
        this->variables["r"] = cdouble(9.46e20, 0.0);
    } else if (system == "EagleNebula") {
        this->variables["M"] = cdouble(1e36, 0.0);
        this->variables["r"] = cdouble(2.36e17, 0.0);
    } else if (system == "M84") {
        this->variables["M"] = cdouble(1.46e45, 0.0);
        this->variables["r"] = cdouble(3.09e22, 0.0);
    } else if (system == "CentaurusA") {
        this->variables["M"] = cdouble(4e41, 0.0);
        this->variables["r"] = cdouble(3.09e21, 0.0);
    } else if (system == "SupernovaSurvey") {
        this->variables["M"] = cdouble(1e30, 0.0);
        this->variables["r"] = cdouble(1e10, 0.0);
    }
}
// Dynamic variable operations (complex)
void UQFFBuoyancyModule::updateVariable(const std::string& name, cdouble value) {
    this->variables[name] = value;
}
void UQFFBuoyancyModule::addToVariable(const std::string& name, cdouble delta) {
    this->variables[name] += delta;
}
void UQFFBuoyancyModule::subtractFromVariable(const std::string& name, cdouble delta) {
    this->variables[name] -= delta;
}
// Core computation: Full F_U_Bi_i(r, t) for system (approx integral)
cdouble UQFFBuoyancyModule::computeFBi(const std::string& system, double t) {
    setSystemParams(system);
    cdouble integrand = computeIntegrand(t, system);
    cdouble x2 = computeX2(system);
    return integrand * x2; // Approx integral
}
// Sub-equations
cdouble UQFFBuoyancyModule::computeCompressed(const std::string& system, double t) {
    setSystemParams(system);
    return computeIntegrand(t, system);
}
cdouble UQFFBuoyancyModule::computeResonant(const std::string& system) {
    setSystemParams(system);
    return computeDPM_resonance(system);
}
cdouble UQFFBuoyancyModule::computeBuoyancy(const std::string& system) {
    setSystemParams(system);
    return computeUb1(system);
}
cdouble UQFFBuoyancyModule::computeSuperconductive(const std::string& system, double t) {
    setSystemParams(system);
    return computeUi(t, system);
}
double UQFFBuoyancyModule::computeCompressedG(const std::string& system, double t) {
    setSystemParams(system);
    return computeG(t, system);
}
// Output descriptive text of the equation
std::string UQFFBuoyancyModule::getEquationText(const std::string& system, double t) {
    setSystemParams(system);
    std::ostringstream oss;
    oss << "F_U_Bi_i(r, t) = Integral[Integrand(r, t) dt] approximated as Integrand * x2\n";
    oss << "Where Integrand includes terms for base force, momentum, gravity, vacuum stability, LENR resonance, activation, directed energy, magnetic resonance, neutron, relativistic, neutrino, Sweet vac, Kozima drop.\n";
    oss << "System: " << system << ", Time: " << t << " s\n";
    return oss.str();
}
// Print all current variables (for debugging/updates)
void UQFFBuoyancyModule::printVariables() {
    for (const auto& pair : variables) {
        std::cout << std::setw(15) << pair.first << " : " << pair.second << std::endl;
    }
}
// Compute integrand for F_U_Bi_i
cdouble UQFFBuoyancyModule::computeIntegrand(double t, const std::string& system) {
    // Placeholder implementation
    return cdouble(1.0e10, 1.0e5); // Example value
}
// Compute DPM resonance term
cdouble UQFFBuoyancyModule::computeDPM_resonance(const std::string& system) {
    // Placeholder implementation
    return cdouble(1.0e8, 1.0e3); // Example value
}
// Compute x2 from quadratic root approximation
cdouble UQFFBuoyancyModule::computeX2(const std::string& system) {
    // Placeholder implementation
    return cdouble(1.0e15, 0.0); // Example value
}
// Compute LENR term
cdouble UQFFBuoyancyModule::computeLENRTerm(const std::string& system) {
    // Placeholder implementation
    return cdouble(1.0e12, 1.0e2); // Example value
}
// Compute gravitational acceleration g(r,t)
double UQFFBuoyancyModule::computeG(double t, const std::string& system) {
    // Placeholder implementation
    return 9.81; // Example value
}
// Compute Q_wave term
cdouble UQFFBuoyancyModule::computeQ_wave(double t, const std::string& system) {
    // Placeholder implementation
    return cdouble(1.0e6, 1.0e1); // Example value
}
// Compute Ub1 buoyancy term
cdouble UQFFBuoyancyModule::computeUb1(const std::string& system) {
    // Placeholder implementation
    return cdouble(1.0e9, 1.0e4); // Example value
}
// Compute Ui superconductive term
cdouble UQFFBuoyancyModule::computeUi(double t, const std::string& system) {
    // Placeholder implementation
    return cdouble(1.0e7, 1.0e2); // Example value
    }

// ===== SURFACE MAGNETIC FIELD MODULE FOR GALACTIC SYSTEMS =====

#ifndef SURFACE_MAGNETIC_FIELD_MODULE_GAL_H
#define SURFACE_MAGNETIC_FIELD_MODULE_GAL_H

#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

class SurfaceMagneticFieldModule {
private:
    
    // ========== CORE PARAMETERS (Original UQFF - Preserved) ==========
    // Note: Can be extended with dynamic parameters via setVariable()
    std::map<std::string, double> variables;
    std::map<std::string, std::vector<double>> variable_history;
    std::map<std::string, std::string> variable_dependencies;
    bool self_learning_enabled;
    double learning_rate;
    int update_counter;
    
    // Dynamic helper functions
    void updateDependencies(const std::string& changed_var);
    double computeGradient(const std::string& var, const std::string& target);
    void recordHistory(const std::string& name, double value);
    // ========== SELF-EXPANDING FRAMEWORK MEMBERS ==========
    std::map<std::string, double> dynamicParameters;
    std::vector<std::unique_ptr<PhysicsTerm>> dynamicTerms;
    std::map<std::string, std::string> metadata;
    bool enableDynamicTerms;
    bool enableLogging;
    double learningRate;



public:
    // Constructor with galactic-enhanced dynamic capabilities
    SurfaceMagneticFieldModule();
    
    // Core magnetic field computations for galactic systems
    double computeB_j(double t, double B_s);
    double computeB_s_min();
    double computeB_s_max();
    double computeU_g3_example(double t, double B_s);
    double computeGalacticCoupling(double B_field, double stellar_mass);
    
    // Dynamic variable operations
    void updateVariable(const std::string& name, double value);
    
    // Advanced dynamic capabilities for galactic systems
    void autoCalibrate(const std::string& observable, double target_value, double tolerance = 0.01);
    void adaptiveUpdate(double dt, const std::string& feedback_param = "");
    void scaleToGalacticData(const std::map<std::string, double>& gal_data);
    void addCustomVariable(const std::string& name, double value, const std::string& dependency = "");
    std::map<std::string, double> getVariableHistory(const std::string& name, int steps = 10);
    void enableSelfLearning(bool enable);
    void exportState(const std::string& filename);
    void importState(const std::string& filename);
    
    // Enhanced magnetic field equations for galactic systems
    std::string getEquationText();
};

#endif // SURFACE_MAGNETIC_FIELD_MODULE_GAL_H

// ===== SURFACE MAGNETIC FIELD MODULE IMPLEMENTATION FOR GALACTIC SYSTEMS =====

// Enhanced SurfaceMagneticFieldModule constructor with galactic capabilities
SurfaceMagneticFieldModule::SurfaceMagneticFieldModule() {
    // Initialize dynamic capabilities
    self_learning_enabled = false;
    learning_rate = 0.07;  // Balanced learning rate for galactic stability
    update_counter = 0;
    
    // Universal constants
    variables["B_s_min"] = 1e-7;                    // T (quiet galactic systems)
    variables["B_s_max"] = 2.0;                     // T (active galactic systems)
    variables["B_ref"] = 2.0;                       // T (reference max for galaxies)
    variables["k_3"] = 2.8;                         // Enhanced coupling for galaxies
    variables["omega_s"] = 1.7e-5;                  // rad/s (galactic system frequency)
    variables["P_core"] = 1.0;                      // Unitless
    variables["E_react"] = 6e46;                    // J (enhanced for galactic systems)
    variables["pi"] = 3.141592653589793;
    variables["t"] = 0.0;                           // s
    
    // Galactic-specific parameters
    variables["stellar_mass_coupling"] = 4.2e-16;  // Stellar mass-magnetic coupling
    variables["galactic_rotation"] = 2.5e-16;      // rad/s (typical galactic rotation)
    variables["magnetic_diffusion"] = 4e-11;       // m/s (galactic diffusion)
    variables["convection_velocity"] = 3e3;        // m/s (galactic convection)
    variables["star_formation_rate"] = 1.5;        // M/yr (typical SFR)
    
    // System evolution parameters for galaxies
    variables["evolution_timescale"] = 1.2e15;     // s (galactic evolution timescale)
    variables["thermal_coupling"] = 1.4e-8;        // Thermal-magnetic coupling
    variables["supernova_rate"] = 0.02;            // yr (supernova rate)
    variables["spiral_arm_enhancement"] = 1.8;     // Spiral arm magnetic enhancement
    variables["dark_matter_coupling"] = 5.5e-25;   // Dark matter interaction
}

// Compute minimum surface magnetic field for galactic systems
double SurfaceMagneticFieldModule::computeB_s_min() {
    return variables["B_s_min"];
}

// Compute maximum surface magnetic field for galactic systems
double SurfaceMagneticFieldModule::computeB_s_max() {
    return variables["B_s_max"];
}

// Compute scaled B_j based on time t and surface field B_s for galactic systems
double SurfaceMagneticFieldModule::computeB_j(double t, double B_s) {
    variables["t"] = t;
    
    // Enhanced magnetic field evolution with galactic coupling
    double gal_oscillation = 0.7 * std::sin(variables["omega_s"] * t);
    double rotation_factor = 1.0 + variables["galactic_rotation"] * t / (2 * M_PI);
    double spiral_factor = variables["spiral_arm_enhancement"] * std::cos(variables["galactic_rotation"] * t * 4);
    double thermal_factor = 1.0 + variables["thermal_coupling"] * variables["star_formation_rate"];
    double base_b = variables["B_ref"] + gal_oscillation * rotation_factor * spiral_factor * thermal_factor;
    
    return base_b * (B_s / variables["B_ref"]);
}

// Compute galactic-magnetic field coupling
double SurfaceMagneticFieldModule::computeGalacticCoupling(double B_field, double stellar_mass) {
    double coupling_strength = variables["stellar_mass_coupling"];
    double rotation_factor = variables["galactic_rotation"] / 2.5e-16;  // Normalized to typical rotation
    double sfr_factor = variables["star_formation_rate"] / 1.5;  // Normalized to typical SFR
    double dm_factor = variables["dark_matter_coupling"] * 1e25;  // Dark matter contribution
    
    return coupling_strength * B_field * stellar_mass * rotation_factor * sfr_factor * (1.0 + dm_factor);
}

// Compute U_g3 example with galactic enhancement
double SurfaceMagneticFieldModule::computeU_g3_example(double t, double B_s) {
    double k_3 = variables["k_3"];
    double b_j = computeB_j(t, B_s);
    double cos_term = std::cos(variables["omega_s"] * t * variables["pi"]);
    double p_core = variables["P_core"];
    double e_react = variables["E_react"];
    
    // Galactic enhancement factor
    double gal_enhancement = 1.0 + computeGalacticCoupling(b_j, variables["star_formation_rate"] * 1e9);
    
    // Supernova feedback
    double sn_factor = 1.0 + variables["supernova_rate"] * std::sin(t / 3.15e7);  // Annual cycle
    
    return k_3 * b_j * cos_term * p_core * e_react * gal_enhancement * sn_factor;
}

// Enhanced updateVariable with galactic dynamic capabilities
void SurfaceMagneticFieldModule::updateVariable(const std::string& name, double value) {
    if (variables.find(name) != variables.end()) {
        variables[name] = value;
    } else {
        std::cerr << "Variable '" << name << "' not found. Adding with value " << value << std::endl;
        variables[name] = value;
    }
    
    // Record history for dynamic capabilities
    recordHistory(name, value);
    
    // Update dependencies
    updateDependencies(name);
    
    // Increment update counter
    update_counter++;
    
    // Trigger self-learning if enabled
    if (self_learning_enabled && update_counter % 6 == 0) {  // Moderate frequency for galactic systems
        adaptiveUpdate(1.0, name);
    }
}

// Auto-calibrate magnetic field parameters for galactic systems
void SurfaceMagneticFieldModule::autoCalibrate(const std::string& observable, double target_value, double tolerance) {
    if (variables.find(observable) == variables.end()) {
        std::cerr << "Observable '" << observable << "' not found for calibration." << std::endl;
        return;
    }
    
    double current_value = variables[observable];
    double error = std::abs(current_value - target_value) / target_value;
    
    if (error > tolerance) {
        // Galactic-specific parameter adjustment
        std::vector<std::string> tunable_params = {"B_ref", "k_3", "omega_s", "P_core", "stellar_mass_coupling", "thermal_coupling", "galactic_rotation"};
        
        for (const auto& param : tunable_params) {
            double gradient = computeGradient(param, observable);
            if (std::abs(gradient) > 1e-25) {
                double adjustment = learning_rate * (target_value - current_value) / gradient;
                variables[param] += adjustment;
                recordHistory(param, variables[param]);
            }
        }
        
        std::cout << "Auto-calibrated galactic magnetic " << observable << " from " << current_value 
                  << " to target " << target_value << " (error: " << error << ")" << std::endl;
    }
}

// Adaptive magnetic field evolution for galactic systems
void SurfaceMagneticFieldModule::adaptiveUpdate(double dt, const std::string& feedback_param) {
    if (!self_learning_enabled) return;
    
    // Galactic evolution timescale
    double evolution_factor = std::exp(-dt / variables["evolution_timescale"]);
    
    // Star formation rate evolution
    variables["star_formation_rate"] *= (1.0 + 0.0003 * std::sin(dt / 1e12));
    
    // Adaptive magnetic field reference with galactic coupling
    double gal_factor = variables["star_formation_rate"] / 1.5;
    variables["B_ref"] = variables["B_s_max"] * (0.4 + 0.6 * gal_factor);
    
    // Enhanced magnetic diffusion effects for galactic systems
    double diffusion_decay = std::exp(-dt * variables["magnetic_diffusion"] / 1e4);
    variables["k_3"] *= diffusion_decay;
    
    // Galactic rotation effects
    double rotation_enhancement = 1.0 + 0.12 * variables["convection_velocity"] / 3e3;
    variables["omega_s"] *= rotation_enhancement;
    
    // Supernova rate variability
    variables["supernova_rate"] *= (1.0 + 0.002 * std::cos(dt / 1e9));
    
    // Dark matter coupling evolution
    variables["dark_matter_coupling"] *= (1.0 + 0.0001 * std::sin(dt / 1e13));
    
    recordHistory("adaptive_time", dt);
    std::cout << "Galactic magnetic adaptive update: B_ref=" << variables["B_ref"] 
 * Enhanced: November 04, 2025 - Added self-expanding capabilities
              << ", SFR=" << variables["star_formation_rate"] << std::endl;
}

// Scale to galactic observational data
void SurfaceMagneticFieldModule::scaleToGalacticData(const std::map<std::string, double>& gal_data) {
    for (const auto& data : gal_data) {
        if (data.first == "star_formation_rate") {
            double scaling = data.second / variables["star_formation_rate"];
            variables["star_formation_rate"] = data.second;
            variables["thermal_coupling"] *= scaling;
            variables["stellar_mass_coupling"] *= std::sqrt(scaling);
        }
        
        if (data.first == "galactic_rotation") {
            variables["galactic_rotation"] = data.second;
            variables["omega_s"] *= data.second / 2.5e-16;
        }
        
        if (data.first == "magnetic_field_strength") {
            double scaling = data.second / variables["B_s_max"];
            variables["B_s_max"] = data.second;
            variables["B_ref"] *= scaling;
        }
        
        if (data.first == "supernova_rate") {
            variables["supernova_rate"] = data.second;
            variables["E_react"] *= (1.0 + 0.1 * data.second / 0.02);
        }
    }
    std::cout << "Scaled galactic magnetic module to " << gal_data.size() << " galactic observations." << std::endl;
}

// Add custom magnetic variables for galactic systems
void SurfaceMagneticFieldModule::addCustomVariable(const std::string& name, double value, const std::string& dependency) {
    variables[name] = value;
    if (!dependency.empty()) {
        variable_dependencies[name] = dependency;
    }
    recordHistory(name, value);
    std::cout << "Added custom galactic magnetic variable: " << name << " = " << value << std::endl;
}

// Get magnetic parameter history for galactic systems
std::map<std::string, double> SurfaceMagneticFieldModule::getVariableHistory(const std::string& name, int steps) {
    std::map<std::string, double> history;
    if (variable_history.find(name) != variable_history.end()) {
        auto& hist = variable_history[name];
        int start = std::max(0, (int)hist.size() - steps);
        for (int i = start; i < (int)hist.size(); i++) {
            history["step_" + std::to_string(i)] = hist[i];
        }
    }
    return history;
}

// Enable galactic magnetic self-learning
void SurfaceMagneticFieldModule::enableSelfLearning(bool enable) {
    self_learning_enabled = enable;
    if (enable) {
        std::cout << "Galactic magnetic self-learning enabled with rate: " << learning_rate << std::endl;
    } else {
        std::cout << "Galactic magnetic self-learning disabled." << std::endl;
    }
}

// Export galactic magnetic state
void SurfaceMagneticFieldModule::exportState(const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        file << "# SurfaceMagneticFieldModule Galactic State Export" << std::endl;
        file << "update_counter=" << update_counter << std::endl;
        file << "learning_rate=" << learning_rate << std::endl;
        file << "self_learning_enabled=" << (self_learning_enabled ? 1 : 0) << std::endl;
        
        for (const auto& var : variables) {
            file << var.first << "=" << var.second << std::endl;
        }
        file.close();
        std::cout << "Galactic magnetic state exported to: " << filename << std::endl;
    }
}

// Import galactic magnetic state
void SurfaceMagneticFieldModule::importState(const std::string& filename) {
    std::ifstream file(filename);
    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            if (line[0] == '#') continue;
            
            size_t eq_pos = line.find('=');
            if (eq_pos != std::string::npos) {
                std::string key = line.substr(0, eq_pos);
                std::string value_str = line.substr(eq_pos + 1);
                
                if (key == "update_counter") {
                    update_counter = std::stoi(value_str);
                } else if (key == "learning_rate") {
                    learning_rate = std::stod(value_str);
                } else if (key == "self_learning_enabled") {
                    self_learning_enabled = (std::stoi(value_str) == 1);
                } else {
                    variables[key] = std::stod(value_str);
                }
            }
        }
        file.close();
        std::cout << "Galactic magnetic state imported from: " << filename << std::endl;
    }
}

// Enhanced equation text for galactic systems
std::string SurfaceMagneticFieldModule::getEquationText() {
    return "Galactic-Enhanced Magnetic Field Equations:\n"
           "B_j  (B_ref + 0.7 sin(?_s t) * O_rot * S_spiral * T_thermal) * (B_s / B_ref) T\n"
           "U_g3 = k_3 *  B_j * cos(?_s t p) * P_core * E_react * (1 + Gal_coupling) * SN_factor\n"
           "Gal_coupling = ?_gal * B_field * M_stellar * (O/O_0) * (SFR/SFR_0) * (1 + DM_coupling)\n"
           "Where:\n"
           "- B_s = [1e-7, 2.0] T (galactic system range)\n"
           "- SFR = " + std::to_string(variables["star_formation_rate"]) + " M/yr (star formation rate)\n"
           "- O_rot = " + std::to_string(variables["galactic_rotation"]) + " rad/s (galactic rotation)\n"
           "- ?_gal = " + std::to_string(variables["stellar_mass_coupling"]) + " (stellar coupling)\n"
           "- SN_rate = " + std::to_string(variables["supernova_rate"]) + " yr (supernova rate)\n"
           "Galactic Systems: M74, Eagle Nebula (M16), M84, Centaurus A,\n"
           "Supernova Survey\n"
           "Enhanced Features: Stellar mass coupling, galactic rotation effects,\n"
           "spiral arm enhancement, supernova feedback, dark matter coupling.";
}

// Helper functions for galactic magnetic field module
void SurfaceMagneticFieldModule::updateDependencies(const std::string& changed_var) {
    if (changed_var == "star_formation_rate") {
        // Update thermal coupling based on star formation rate
        variables["thermal_coupling"] = 1.4e-8 * (variables["star_formation_rate"] / 1.5);
        // Update supernova rate from SFR
        variables["supernova_rate"] = 0.02 * (variables["star_formation_rate"] / 1.5);
    }
    
    if (changed_var == "B_s_max") {
        variables["B_ref"] = variables["B_s_max"];
    }
    
    if (changed_var == "galactic_rotation") {
        // Update omega_s based on galactic rotation
        variables["omega_s"] = 1.7e-5 * (variables["galactic_rotation"] / 2.5e-16);
        // Update evolution timescale
        variables["evolution_timescale"] = 2 * M_PI / variables["galactic_rotation"] * 1e6;
    }
    
    if (changed_var == "omega_s") {
        // Update spiral arm enhancement frequency
        variables["spiral_arm_enhancement"] = 1.8 * std::cos(4 * variables["omega_s"]);
    }
}

double SurfaceMagneticFieldModule::computeGradient(const std::string& var, const std::string& target) {
    if (variables.find(var) == variables.end() || variables.find(target) == variables.end()) {
        return 0.0;
    }
    
    double original_value = variables[var];
    double original_target = variables[target];
    
    // Small perturbation
    double delta = original_value * 1e-6;
    variables[var] += delta;
    
    // Recompute target with galactic coupling
    double new_target = computeB_j(variables["t"], variables["B_ref"]);
    
    // Restore original value
    variables[var] = original_value;
    
    return (new_target - original_target) / delta;
}

void SurfaceMagneticFieldModule::recordHistory(const std::string& name, double value) {
    variable_history[name].push_back(value);
    
    // Keep only last 120 values for galactic systems (moderate history)
    if (variable_history[name].size() > 120) {
        variable_history[name].erase(variable_history[name].begin());
    }
}
