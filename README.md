# Post-Processing Laminar Oil Transport in a Pipeline

MATLAB-based post-processing tools for Computational Fluid Dynamics (CFD) analysis of laminar oil flow in a cylindrical pipeline.

## Overview

This repository contains MATLAB scripts for analyzing CFD simulation results of laminar oil transport through a pipeline. The tools process simulation data to extract key flow parameters, velocity profiles, pressure distributions, and wall shear stress along the pipeline.

## Features

- **Flow Characterization**
  - Reynolds number calculation
  - Entrance length estimation using multiple formulas
  - Fully developed flow verification at multiple axial locations

- **Data Processing**
  - Import and process CFD simulation files (.phi format)
  - Extract velocity, pressure, and derivative fields
  - Handle both Cartesian and Body-Fitted Coordinate (BFC) grids

- **Flow Analysis**
  - Axial velocity profile extraction
  - Wall shear stress computation
  - Pressure distribution analysis
  - Flow development visualization

## Files

- **CFD1_Processing.m**: Main post-processing script for analyzing pipeline flow data
  - Calculates Reynolds number and entrance length
  - Processes velocity and pressure fields
  - Computes wall shear stress distribution
  
- **XYZ_reduced_19_20.m**: Utility script for importing CFD data files
  - Imports .phi files from CFD simulations
  - Handles grid coordinates and cell data
  - Supports both Cartesian and BFC grids
  - Allows selective variable processing

## Physical Parameters

Default simulation parameters (can be modified in the scripts):
- **Pipe Diameter (D)**: 0.15 m
- **Bulk Velocity (Wb)**: 0.45 m/s
- **Oil Density (ρ)**: 910 kg/m³
- **Kinematic Viscosity (ν)**: 3.5×10⁻⁴ m²/s
- **Pipe Length (L)**: 3.0 m

## Usage

### CFD Data Processing

1. Update the file paths in `XYZ_reduced_19_20.m`:
   ```matlab
   directoryLoad = 'path/to/your/simulation.phi';
   directorySave = 'path/to/save/output.mat';
   ```

2. Specify which variables to process:
   ```matlab
   variables = {'P1', 'V1', 'W1', 'DWDY'};
   ```

3. Run the import script to convert CFD data:
   ```matlab
   XYZ_reduced_19_20
   ```

### Flow Analysis

1. Ensure the processed data (`main.mat`) is available

2. Run the main processing script:
   ```matlab
   CFD1_Processing
   ```

3. The script will output:
   - Reynolds number
   - Entrance length estimates
   - Velocity profiles at specified locations
   - Wall shear stress distribution

## Output Variables

After processing, the following variables are available in the workspace:

### Grid Information
- `NX`, `NY`, `NZ`: Number of cells in each direction
- `X_C`, `Y_C`, `Z_C`: Cell center coordinates
- `X_E`, `Y_N`, `Z_H`: Cell face coordinates

### Flow Fields
- `P1`: Pressure field (NX × NY × NZ)
- `U1`, `V1`, `W1`: Velocity components (NX × NY × NZ)
- `DWDY`: Axial velocity gradient in radial direction

### Derived Quantities
- `tau_wall`: Wall shear stress distribution
- `Re`: Reynolds number
- `Le1`, `Le2`: Entrance length estimates

## Requirements

- MATLAB (tested on R2019a or later)
- Sufficient memory for 3D field data processing

## Applications

This toolset is designed for:
- Academic research in fluid mechanics
- CFD simulation validation
- Pipeline flow analysis
- Educational purposes in computational fluid dynamics

## License

This project is licensed under the GPL-3.0 License - see the [LICENSE](LICENSE) file for details.

## Author

Ziad Kandil

## Acknowledgments

Developed as part of MSc Mathematical Engineering coursework in Computational Fluid Dynamics.
