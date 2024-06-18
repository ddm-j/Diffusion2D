# Diffusion2D

## Overview
Diffusion2D is a numerical solver for 2D diffusion equations, implemented using central differencing for spatial discretization and the forward Euler method for time integration.
**Features**
- Steady state solver (uses `Eigen`'s conjugate gradient linear solver)
- Transient solver (forward differencing)
- UDF support (functions of time) using LuaJIT
- Paraview VTR (rectilinear) mesh output
- Python bindings

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

## Example Solutions
### Steady State Laplace
<img src="images/laplace.png" width="500">

### Transient Diffusion
Using the Python bindings, arbitrary initial conditions can be specified to the solver.

<img src="images/text_animation.gif" width="500">

### Transient Simulation with UDFs
UDFs can be specified using valid lua code, either by direct text (string) input, or by specifying a lua script. The lua code is jit compiled (LuaJIT) for **_near C_** performance without need to recompile.

<img src="images/UDF_animation.gif" width="500">
