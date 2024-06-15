# Diffusion2D

## Overview
Diffusion2D is a numerical solver for 2D diffusion equations, implemented using central differencing for spatial discretization and the forward Euler method for time integration.

## Equations
The solver handles the following equations:

1. **Laplace Equation:**
   ```math
   \nabla^2 T = 0
   ```

2. **Poisson Equation:**
   ```math
   \nabla^2 T + f = 0
   ```
   Here, `f` represents the source term.

3. **Unsteady Diffusion Equation:**
   ```math
   \nabla^2 T + f = \frac{\partial T}{\partial t}
   ```

## Methods
The numerical methods used are:
- **Spatial Discretization:** Central differencing
- **Time Integration:** Forward differencing

---
