# AeroFlex Program Analysis and Plan for Experimental Wing Implementation

## Program Structure and Workflow (based on `examples/example1/example1.m`)

The program appears to be written in MATLAB. The general workflow for setting up and running a simulation is as follows:

1.  **Setup Global Parameters (`softPARAMS`):**
    *   Define system units (e.g., `isIS`).
    *   Define boundary conditions (e.g., `isPINNED`, `vecFREEDEG`).
    *   Include/exclude gravity (`isGRAV`, `g`).
    *   Set parameters for equilibrium calculation (e.g., `isITER`, `numITER`).
    *   Select aerodynamic model (`modAED`).
    *   Control structural Jacobian updates (`updateStrJac`).
    *   Enable/disable 3D plotting during equilibrium calculation (`plota3d`).

2.  **Structure Initialization (`load_structure` function):**
    *   This function creates an `Airplane` object.
    *   It takes the number of elements (`numele`) and a damping coefficient (`damping`) as input.
    *   It calls `create_flexible_member` to define the properties of the wing (or other flexible components).
    *   `create_flexible_member` defines:
        *   `Length`: Span of the beam.
        *   `KG`: Sectional rigidity matrix. This is a diagonal matrix `diag([EA, GJ, EI_flap, EI_chordwise])`.
        *   `CG`: Sectional damping matrix. In `example1`, this is `damp_ratio * KG`.
        *   `pos_cg`: Position of the section's center of gravity relative to its elastic axis `[x_offset, y_offset, z_offset]`.
        *   `Inertia`: Mass moment of inertia matrix per unit length `diag([Ixx, Iyy, Izz])` for the section, where Ixx is torsional inertia about elastic axis, Iyy is about chordwise axis, Izz is about flapwise axis.
        *   `mcs`: Mass per unit length (kg/m).
        *   `rot0`: A struct defining initial `dihedral`, `sweep`, and `twist` angles (in radians) of the beam.
        *   `aeroparams`: Aerodynamic parameters (empty in `example1`).
        *   `geometry`: Basic geometric parameters like chord, thickness (not fully clear from `example1`).
        *   `num_elements`: Number of discrete elements along the beam.
    *   `create_flexible_member` then calls `create_uniform_structure` to generate the member object based on these properties. For non-uniform structures, this function would need to be adapted or replaced.
    *   The `Airplane` object can consist of multiple members, a fuselage, and engines/stores, e.g., `ap = Airplane({flexible_member_list}, fuselage_obj, {engine_list});`.

3.  **Modal Analysis (Optional):**
    *   Uses `StructuralModes(ap)` to compute natural frequencies and mode shapes.
    *   Allows plotting of mode shapes.

4.  **Equilibrium Calculation:**
    *   Uses `trimairplane(ap, V, altitude, Vwind, throttle, deltaflap)` to find the static equilibrium state of the aircraft under given flight conditions and gravity.
    *   Outputs `rb_eq` (rigid body equilibrium state) and `strain_eq` (structural equilibrium strains).
    *   Allows plotting of the deformed structure.

5.  **Nonlinear Simulation:**
    *   Uses `simulate(ap, tspan, strain0, beta0, k0, Vwind, engine_input_func, aero_surface_input_func, solver_type)` to run a time-domain simulation.
    *   `tspan`: Simulation time interval `[t_start, t_end]`.
    *   `strain0`: Initial structural strains.
    *   `beta0`: Initial rigid body velocities `[u, v, w, p, q, r]`.
    *   `k0`: Initial rigid body position/orientation (not fully clear what `k0` represents, `example1` uses `[0;0;0;0]`).
    *   `Vwind`: Wind velocity.
    *   The simulation outputs time histories of strains, strain rates, Lagrange multipliers, rigid body velocities, and kinetic/potential energies.
    *   Results can be post-processed (e.g., calculate tip displacement) and animated (`airplanemovie`).

## Plan for Implementing the Experimental Wing

1.  **Create a New Example Directory:**
    *   Create `examples/experimental_wing/`.
    *   Copy `examples/example1/example1.m` to `examples/experimental_wing/experimental_wing_sim.m` as a template.
    *   (Optional) Copy `README.md` from `example1` and adapt it later.

2.  **Modify `experimental_wing_sim.m`:**

    *   **Update `load_structure` and `create_flexible_member` functions:**
        *   **Wing Properties:**
            *   `Length = 0.4508` (m)
            *   `mcs = 0.2351` (kg/m)
            *   `pos_cg`: Elastic axis is at 50% chord, CG at 49% chord. Chord *c* = 0.0508 m.
                *   x-offset = (0.49 - 0.50) * 0.0508 = -0.000508 m.
                *   Assuming y-offset = 0, z-offset = 0. So, `pos_cg = [-0.000508, 0, 0]`.
            *   `KG` matrix:
                *   *EI<sub>1</sub>* (Flap bending) = 0.4186 Nm<sup>2</sup>. Let this be `K33`.
                *   *EI<sub>2</sub>* (Chordwise bending) = 0.1844 × 10<sup>2</sup> = 18.44 Nm<sup>2</sup>. Let this be `K44`.
                *   *GJ* (Torsional rigidity) = 0.9539 Nm<sup>2</sup>. Let this be `K22`.
                *   Axial rigidity *EA* is not given. I will need to ask the user for this value or make an assumption (e.g., use a large value if it's not expected to be significant, or a value scaled from bending/torsional). **For now, I'll use the `K11` value from `example1` (1e10) as a placeholder.**
                *   So, `KG = diag([EA_placeholder, 0.9539, 0.4186, 18.44])`.
            *   `Inertia` matrix (per unit length):
                *   Moment of inertia (50% chord) = 0.2056 × 10<sup>-4</sup> kgm. This is likely the pitch moment of inertia per unit length *around the elastic axis*. For a wing section, the elastic axis is often near the quarter-chord or mid-chord. If this is \(I_{xx}\) (torsional inertia about spanwise axis through EA), then `I11 = 0.2056e-4`. The other two principal moments of inertia per unit length (for bending) are not directly given for the wing section itself. These would be mass distribution properties. **I may need to ask the user or derive them if possible, or assume they are small/negligible compared to the slender body if it's attached.**
                *   The `example1.m` uses `Inertia = diag([I11 I22 I33]);` where `I11` is about the beam's x-axis (torsion), `I22` about y-axis (chordwise bending inertia term), and `I33` about z-axis (flapwise bending inertia term).
                *   The given "Moment of inertia (50% chord)" value is \(0.2056 \times 10^{-4} \text{ kgm}\). This unit is ambiguous. If it's \( \text{kg} \cdot \text{m}^2/\text{m} \), then it's a distributed inertia. I'll assume this is \(I_{xx}\) (torsional inertia per unit length). The other components \(I_{yy}\) and \(I_{zz}\) for the wing section are not provided. These are mass distribution properties. **I will assume `I22` and `I33` for the wing section are zero for now, or ask the user, as these relate to how the mass `mcs` is distributed around the elastic axis.**
            *   `CG` (Damping matrix):
                *   Modal damping ratios are given: *ξ<sub>w</sub>* (flap) = 0.02, *ξ<sub>v</sub>* (chordwise) = 0.025, *ξ<sub>φ</sub>* (torsion) = 0.031.
                *   The current example uses `CG = damp_ratio * KG`. This implies the same damping ratio for all modes if KG is diagonal. This is a simplification.
                *   To use specific modal damping, the program might need a different way to input damping (e.g. directly providing a C matrix or using modal damping transformation). I will first try using an average of the given damping ratios for `damp_ratio`, e.g., `(0.02 + 0.025 + 0.031)/3`. Or, use the flap damping `0.02` as a starting point if a single value is easiest. **For now, I'll use `damp_ratio = 0.02` (flap damping).**
            *   `rot0`: Assume `dihedral = 0, sweep = 0, twist = 0` unless specified otherwise. The table doesn't provide these.
        *   **Slender Body (Tip Mass/Store):**
            *   Radius *R* = 0.4762 × 10<sup>-2</sup> m
            *   Chord length *c<sub>SB</sub>* = 0.1406 m
            *   Mass *M* = 0.0417 kg
            *   Moment of inertia *I<sub>x</sub>* = 0.9753 × 10<sup>-4</sup> kgm<sup>2</sup>
            *   Moment of inertia *I<sub>y</sub>* = 0.3783 × 10<sup>-5</sup> kgm<sup>2</sup>
            *   Moment of inertia *I<sub>z</sub>* = 0.9753 × 10<sup>-4</sup> kgm<sup>2</sup>
            *   This needs to be added to the `Airplane` object. `example1` shows `motor1 = []`. This might be the place to add a point mass. I'll need to search the codebase or documentation for how to define and attach a point mass/store with its own inertial properties at the tip of the wing.
            *   **Action: Search codebase for "point mass", "tip mass", "store", or how `motor` objects are defined and if they can represent a generic mass with inertia.**

    *   **Global Parameters (`softPARAMS`):**
        *   Verify `isIS = 1` (International System).
        *   `isPINNED = 1` (assuming cantilevered wing, pinned at root).
        *   `isGRAV = 1`, `g = 9.81`.
        *   Other `softPARAMS` can be kept as in `example1` initially.

    *   **Simulation Parameters:**
        *   Set `numele` (e.g., 10 or 20 elements).
        *   Adjust simulation time `tSIM` as needed.
        *   Initial conditions (`strain0`, `beta0`, `k0`) likely zero for a test starting from undeformed state.

3.  **Testing and Verification:**
    *   Run the simulation.
    *   Compare results (e.g., natural frequencies, equilibrium deflection, dynamic response) with any available experimental data or expected behavior if possible. The provided table is "Experimental wing model data," suggesting there might be results to compare against.

## Missing Information / Assumptions:

*   **Axial Rigidity (EA) for the wing:** Not provided. Assuming `1e10` (from `example1`) as a placeholder.
*   **Mass moments of inertia *I<sub>yy</sub>* and *I<sub>zz</sub>* per unit length for the wing section:** Not provided. Assuming zero or small for now. The provided "Moment of inertia (50% chord)" is assumed to be \(I_{xx}\) (torsional inertia) per unit length.
*   **Damping Implementation:** Using a single `damp_ratio` based on flap damping (0.02) as an initial approach due to `CG = damp_ratio * KG` structure in `example1`. Implementing distinct modal damping ratios might require deeper changes.
*   **How to model the "Slender body":** Need to find out how to add a tip mass/store with its own M and I properties.
*   **Initial orientation `rot0` (dihedral, sweep, twist):** Assumed to be zero.
*   **The nature of "Moment of inertia (50% chord) = 0.2056e-4 kgm":** Interpreted as \(I_{xx}\) (torsional inertia per unit length, so \(0.2056 \times 10^{-4} \text{ kgm}^2/\text{m}\)). If it's total inertia for the wing, the input to `create_uniform_structure` needs adjustment. The unit "kgm" is unusual for MOI. Usually kg.m^2. Given "mass per unit length" is also given, "MOI (50% chord)" is likely also a "per unit length" property.

Next, I will search the codebase to understand how to implement the "Slender body" as a tip mass. 