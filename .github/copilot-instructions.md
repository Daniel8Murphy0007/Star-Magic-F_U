
# Star-Magic: Unified Quantum Field Force (UQFF) – AI Agent Instructions

This project is a theoretical physics and computational modeling platform for the Unified Quantum Field Force (UQFF) and Aetheric Propulsion, integrating quantum, relativistic, and astrophysical phenomena. It is not a standard software project—AI agents must follow these project-specific conventions and workflows:

## 1. Project Structure & Key Files
- **`index.js`**: Main computational engine. Implements 26-layer UQFF gravity, magnetism, buoyancy, and aether models. Integrates mathematical frameworks from `MAIN_1.cpp` and references C++ modules (e.g., `Source13.cpp`, `Source14.cpp`, etc.).
- **`Star-Magic.md`**: Theoretical documentation. Contains the full UQFF mathematical framework, force definitions, and scientific context. Always reference this for scientific constants, notation, and theory.
- **`README.md`**: High-level overview and project intent.
- **Module JS files** (e.g., `ngc2264_uqff.js`, `smbhbinary_uqff.js`): Each models a specific astrophysical system, using UQFF patterns and formulas. See `FOLLOW_UP.md` for known formula issues.
- **`MAIN_1.cpp`**: Source of core mathematical frameworks and comments. Many constants and equations in JS are ported from here.
- **`SETUP.md`/`WHAT_WAS_SETUP.md`**: Setup, usage, and module architecture documentation. Review for build/run/test instructions and project status.

## 2. Development & Execution Workflows
- **Install dependencies**: `npm install`
- **Run main engine**: `npm start` or `node index.js`
- **Run tests**: `npm test` (integration), `npm verify` (verification), `npm demo` (quick demo)
- **Module development**: Add/modify modules in JS, following the structure of existing system modules. Use scientific constants and patterns from `index.js` and `MAIN_1.cpp`.
- **Debugging**: Use in-code `console.log`/`console.error` for variable inspection. Many modules and the main engine print detailed intermediate results for debugging.

## 3. Project-Specific Patterns & Conventions
- **26-layer polynomial structure**: All force calculations (gravity, magnetism, buoyancy) are modeled as sums over 26 quantum layers. See `index.js` and `MAIN_1.cpp` for implementation.
- **Scientific notation**: Use Greek letters, subscripts, and SI units in all calculations and comments. E.g., `ΔUg_i`, `ρ_vac`, `ω_s`, `M_s = 1.989e30 kg`.
- **Module interface**: Each system module exports a class with methods for core calculations (e.g., `computeStarFormation`, `computeCoalescence`). Support dynamic parameter updates and method expansion (see `updateParameter`, `expand`).
- **Cross-language integration**: Many JS modules are direct ports or wrappers for C++ code in `MAIN_1.cpp` and `Source*.cpp`. Comments often reference the original C++ source for traceability.
- **Known formula issues**: See `FOLLOW_UP.md` for modules needing formula refinement (e.g., Jeans mass, Peters-Mathews, M-σ relation). Always check this file before modifying scientific calculations.

## 4. Testing, Validation, and Documentation
- **No automated tests yet**: Testing is manual via `npm test`, `npm verify`, and direct script execution. See `SETUP.md` for details.
- **Validation**: Compare results with astronomical data and published literature. See `FOLLOW_UP.md` and `WHAT_WAS_SETUP.md` for validation priorities.
- **Documentation**: All scientific and mathematical context is in `Star-Magic.md`. Update this file with any new theoretical developments or changes to the unified field equation.

## 5. External Dependencies & Data
- **@anthropic-ai/sdk**: Installed for future AI integration (not yet used in main code).
- **Node.js built-ins**: Some modules use `fs` for file I/O.
- **No web or API endpoints yet**: Planned for future expansion (see `SETUP.md`).

## 6. AI Agent Guidance
- **Always prioritize scientific rigor and traceability to theory**. All code must support, validate, or extend the UQFF framework as described in `Star-Magic.md` and `MAIN_1.cpp`.
- **Do not introduce generic software patterns** unless justified by the project’s scientific context.
- **When in doubt, reference the original C++ or markdown documentation** for constants, formulas, and intent.
- **Document all changes and new modules with clear scientific rationale**.

---
For any unclear or incomplete section, consult the project owner or request clarification. This project is a living research platform—AI agents should iterate and refine instructions as the theory and codebase evolve.