# SSA_python
Python based optimised implementation for ray tracing and SSA (Single Scattering Albedo) calculations
To run: python/pypy SSA_retrive_parallel.py <params>.txt <ssa_output> <output_SSA_hdr>

## Overview
This codebase provides optimized implementations for multi-dimensional grid interpolation and ray tracing using Numba for Just-In-Time (JIT) compilation. The functions are designed to improve computational performance in applications that require fast interpolation and spectral data processing.

## Features
Efficient multi-dimensional interpolation with the multidimgrid function.
Ray interpolation with the rayterpolate function.
Bilinear interpolation through the bilint_new function.
Uses Numba to accelerate numerical computations.

## Installation
Requirements
Python 3.6 or higher
NumPy
Numba

## Install Dependencies
pip install -r requirements.txt

## Run
python SSA_retrive_parallel.py <params>.txt <ssa_op_file> <ssa_hdr_op_file>

## Usage

### 1. multidimgrid(ndim_interp, unit_grid, axes, vector)
#### Parameters:

ndim_interp: The number of dimensions for interpolation (int).
unit_grid: A 2D array representing the grid values to interpolate.
axes: A 2D array defining the axes for interpolation.
vector: A vector to be interpolated.
Returns: The interpolated value (float).

### 2. bilint_new(BAM, nbands, ddcols, icolx, dxf, dyf, incr, numu, mxband, mxwaves, mxcols, mxline, mxsamp, mxcolssm, out)
#### Parameters:

BAM: Input data array.
nbands: Number of bands (int).
ddcols: Number of columns (int).
icolx: Column index (int).
dxf: X-coordinate (float).
dyf: Y-coordinate (float).
incr: Increment (int).
numu: Number of units (int).
mxband, mxwaves, mxcols, mxline, mxsamp, mxcolssm: Various maximum dimensions for the inputs.
out: Output array.
Returns: The output array after bilinear interpolation.

### 3. rayterpolate(trdrsp, disortout, mlst, nlst, endboundsubband, pressure, xpres, nsurfpres, alb, nalb, mxcolssm, mxline, mxband, mxsamp, mxwaves, max_val, ssa_out)
#### Parameters:

trdrsp: Input spectral data (1D array).
disortout: Output data for interpolation (2D array).
mlst, nlst: Indices for the input data (arrays).
endboundsubband: Number of subbands (int).
pressure: Pressure value (float).
xpres: Array of x-values (1D array).
nsurfpres: Number of surface pressures (int).
alb: Albedo values (1D array).
nalb: Number of albedos (int).
mxcolssm, mxline, mxband, mxsamp, mxwaves: Maximum dimensions.
max_val: Maximum value for calculations.
ssa_out: Output array for results.
Returns: The updated ssa_out array containing results after ray interpolation.

## Performance Considerations
The code is optimized for performance using Numba's JIT compilation. Ensure that Numba is correctly installed and configured to fully utilize these optimizations.


## Contact
For questions or support, please contact:

Your Name: [aniketdeenanath.singh@stonybrook.edu]
