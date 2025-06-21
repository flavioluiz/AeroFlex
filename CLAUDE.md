# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

AeroFlex is an academic MATLAB toolbox for studying the flight dynamics of highly flexible airplanes. It models structural dynamics using strain-based geometrically non-linear beams and aerodynamics using strip theory with quasi-steady, quasi-steady with apparent mass, and full unsteady approaches.

## Architecture

### Core Components

- **`main/`** - Core library files containing all primary functionality
  - `Airplane.m` - Main class that integrates all airplane components (flexible members, engines, fuselage)
  - `Element.m` - Individual structural element representation
  - `Node.m` - Structural node definition
  - `RigidFuselage.m` - Rigid fuselage modeling
  - `Engine.m` - Engine/propulsion modeling
  - `StructuralModes.m` - Modal analysis functionality

### Key Methods in Airplane Class

- `update()` - Updates strain information
- `updateStrJac()` - Updates structural Jacobians
- `linearize()` - Linear stability analysis
- `trimairplane()` - Longitudinal equilibrium finding
- `trimairplanefull()` - Full 3D equilibrium finding
- `simulate()` - Numerical integration of equations of motion
- `plotairplane3D()` - 3D visualization
- `airplanemovie()` - Animation generation

### Analysis Functions

- `flutter_speed.m` - Flutter speed calculation
- `dinamicaflex.m` - Flexible dynamics simulation
- `aeroforceandmoment.m` - Aerodynamic force/moment calculation
- `atmosfera.m` - Atmospheric modeling

## Development Workflow

### Running Examples

Examples are located in `examples/` folder. Each example is self-contained:

```matlab
cd examples/example1
example1  % Run the example
```

### Global Parameters

All examples require setting global `softPARAMS` structure with:
- `isIS` - International System flag
- `isPINNED` - Rigid body DOF pinning
- `vecFREEDEG` - Free degrees of freedom [u v w p q r]
- `isGRAV` - Gravity inclusion flag
- `modAED` - Aerodynamic model selection

### Path Setup

Add main folder to MATLAB path:
```matlab
addpath('../../main');  % From example folders
addpath('./main');      % From root
```

## Example Types

### Core Examples
1. **Structural Dynamics** (`example1/`) - Modal analysis, equilibrium, structural simulation
2. **Aeroelasticity** (`example2/`) - Flutter analysis, linear/nonlinear comparison
3. **Control Reversal** (`example3/`) - Aileron reversal studies
4. **Flight Dynamics** (`example4/`) - Full rigid-body + flexible dynamics
5. **Stability Analysis** (`example5/`) - Work in progress

### Recently Added Test Cases (for Leonardo's PhD Research)
6. **Static AOA Analysis** (`example_static_aoa/`) - Static equilibrium analysis at specified angle of attack
   - Finds equilibrium with aerodynamics + gravity at given AoA
   - Extracts and plots spanwise vertical deflection and torsion
   - Includes lift coefficient vs angle of attack data (`data_cl_alpha.mat`)
   - Uses quasi-steady aerodynamic model with variable cl_alpha along span

7. **Experimental Wing** (`experimental_wing/`) - Experimental wing configuration with tip mass
   - Modal analysis of experimental wing structure with realistic properties
   - Length: 0.4508m, Chord: 0.0508m, Mass: 0.2351 kg/m
   - Tip mass: 0.0417 kg with realistic moment of inertia values
   - Structural dynamics with tip mass effects and gravity equilibrium

8. **PAZY Wing** (`pazy_wing/`) - Based on PAZY wing benchmark case
   - High-fidelity flexible wing model for research validation
   - Modal analysis and equilibrium computation
   - Standard benchmark case used in aeroelastic research community
   - Contains validation data and visualization results

## Key Analysis Capabilities

- Static aeroelastic analysis (divergence, control reversal)
- Dynamic aeroelastic analysis (flutter)
- Nonlinear geometric effects
- Flexible aircraft flight dynamics
- Modal analysis and visualization
- Time-domain simulation
- Stability analysis

## File Structure Notes

- All `.m` files use MATLAB syntax
- Examples include visualization (`.jpg`, `.gif` files)
- Each example folder contains its own README with specific instructions
- Test cases in `tests/` folder provide validation examples
- Recently added examples include experimental data files (`.mat` format)
- Documentation files: `AeroFlex_Program_Analysis.md`, `experimental_wing_data.md`

## Recent Additions for PhD Research

The following files were added to support Leonardo's PhD work:

### New Examples
- `examples/example_static_aoa/` - Static angle of attack analysis with aerodynamic data
  - `example_static_aoa_sim.m` - Main simulation script with variable AoA capability
  - `data_cl_alpha.mat` - Spanwise lift coefficient data for realistic aerodynamic modeling
- `examples/experimental_wing/` - Experimental wing with tip mass modeling
  - `experimental_wing_sim.m` - Complete structural analysis with realistic wing properties
  - Includes modal analysis, equilibrium computation, and time-domain simulation
- `examples/pazy_wing/` - PAZY wing benchmark implementation
  - Standard research benchmark for aeroelastic validation
  - Contains modal analysis and structural dynamics examples

### Documentation
- `AeroFlex_Program_Analysis.md` - Comprehensive program analysis documentation
- `experimental_wing_data.md` - Experimental wing data and configuration details

### Key Features of New Examples
- **Realistic Material Properties**: Based on experimental data and literature values
- **Advanced Aerodynamic Modeling**: Variable lift coefficient distribution along span
- **Tip Mass Effects**: Proper modeling of concentrated masses at wing tips
- **Validation Cases**: Standard benchmarks for research verification
- **Enhanced Visualization**: Improved plotting of deflections, torsion, and modal shapes

These additions focus on advanced aeroelastic analysis, experimental validation, and benchmark test cases commonly used in flexible aircraft research.