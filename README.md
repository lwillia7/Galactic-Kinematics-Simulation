# Galactic-Kinematics-Simulation
Python-based N-body simulation modeling the Milky Way–Andromeda collision.

This repository contains a suite of Python scripts designed to model the gravitational interaction and eventual merger of the Milky Way and Andromeda galaxies.

##Technical Implementation
* **Numerical Integration:** Utilizes a **4th-order Runge-Kutta (RK4)** scheme to solve coupled ordinary differential equations (ODEs) for position and velocity state vectors.
* **Dynamical Friction:** Implements Chandrasekhar's dynamical friction formula to account for momentum exchange within dark matter halos.
* **Reference Frames:** Includes scripts optimized for different barycentric perspectives (MW-centered, Andromeda-centered, and Center-of-Mass).

##Features
* **Vectorized Computations:** Uses NumPy for high-performance linear algebra operations.
* **Animated Visualizations:** Generated using `Matplotlib.animation` to track trajectories over a 3.2-billion-year timescale.
* **Adaptive Simulation:** Configurable time-steps (`dt`) to balance computational speed with integration accuracy.
