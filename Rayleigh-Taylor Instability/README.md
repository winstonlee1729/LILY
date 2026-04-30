# Solve
This case solves the 2D Rayleigh-Taylor Instability (a heavy fluid falling into a light fluid under gravity).

# Spatial Discretization Scheme
1. Convection: QUICK Scheme (First-order upwind near the boundaries).
2. Diffusion: Central Differencing Scheme

# Temporal Discretization Scheme
This solver uses an explicit Euler time integration scheme.

# Pressure-Velocity Coupling Scheme
The pressure-velocity coupling is handled via a Fractional Step method. The resulting pressure Poisson equation is solved using the Red-Black Successive Over-Relaxation (SOR) iterative method.

# Additional Features
* **Multiphase Interface Tracking:** The solver tracks the interface between the two fluids by advecting the density (`rho`) field explicitly.
* **Parallelization:** Matrix assembly and equation solving are highly parallelized using OpenMP (`!$omp parallel do`) for multi-core performance.
* **Standardized Output:** Generates `.vtk` files directly, allowing you to seamlessly visualize the density, pressure, and velocity vector fields in ParaView.