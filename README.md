# LILY

A custom Computational Fluid Dynamics (CFD) solver based on the Finite Volume Method (FVM), designed for high-performance parallel computing (OpenMP).

## Features

* **Advanced Algorithms:** Implements robust pressure-velocity coupling algorithms including PISO, PIMPLE, and SIMPLE for both steady-state and highly transient flows.
* **High-Resolution Schemes:** Utilizes the QUICK scheme for accurate convection modeling with built-in boundary fallbacks.
* **Multiphase & Non-Newtonian Capabilities:** Supports Volume of Fluid (VOF) phase tracking and Carreau viscosity models for complex physics like blood flow and Rayleigh-Taylor instabilities.
* **Parallelized Computation:** Fully optimized with OpenMP to leverage multi-core processors for fast matrix assembly and solving.
* **Standardized Output:** Generates VTK (`.vtk`) files for seamless 2D flow field visualization.

---

## Prerequisites

To compile and visualize the results from this solver, you will need:
1. A Fortran compiler that supports OpenMP (e.g., `gfortran`, Intel `ifort`, or `ifx`).
2. [ParaView](https://www.paraview.org/) (for post-processing and visualization).
3. MATLAB and Python to run plot files.

---

## Compilation & Execution

### 1. Compiling the Code
Open your terminal or command prompt and compile the `.f90` source file. Be sure to enable OpenMP and optimization flags for maximum performance.

**Using GNU Fortran (`gfortran`):**
```bash
gfortran -fopenmp -O3 solver.f90 -o lily
