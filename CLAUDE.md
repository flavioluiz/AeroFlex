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

1. **Structural Dynamics** (`example1/`) - Modal analysis, equilibrium, structural simulation
2. **Aeroelasticity** (`example2/`) - Flutter analysis, linear/nonlinear comparison
3. **Control Reversal** (`example3/`) - Aileron reversal studies
4. **Flight Dynamics** (`example4/`) - Full rigid-body + flexible dynamics
5. **Stability Analysis** (`example5/`) - Work in progress

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