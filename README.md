# Diffusion2D

## Overview
Diffusion2D is a numerical solver for 2D diffusion equations, implemented using central differencing for spatial discretization and the forward Euler method for time integration.

## Equations
The solver handles the following equations:

1. **Laplace Equation:**
   \[
   \nabla^2 T = 0
   \]
   where \( \nabla^2 \) denotes the Laplacian operator.

2. **Poisson Equation:**
   \[
   \nabla^2 T + f = 0
   \]
   Here, \( f \) represents the source term.

3. **Unsteady Diffusion Equation:**
   \[
   \nabla^2 T + f = \frac{\partial T}{\partial t}
   \]
   This equation describes the time-dependent behavior of the diffusion process, where \( \frac{\partial T}{\partial t} \) is the time derivative of \( T \).

## Methods
The numerical methods used are:
- **Spatial Discretization:** Central differencing
- **Time Integration:** Forward Euler method

---

Feel free to explore and contribute to the development of this solver for various 2D diffusion problems.
