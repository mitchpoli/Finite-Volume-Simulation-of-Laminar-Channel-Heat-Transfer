# Finite Volume Simulation of Laminar Channel Heat Transfer  
**ME-474 – Numerical Flow Simulation (EPFL)**  

## Overview

This project implements a 2D finite volume solver for the steady convection–diffusion equation in a laminar Poiseuille channel flow.

The objective is to compute the temperature field **T(x,y)** and analyze:

- Thermal entrance length  
- Local and mean temperatures  
- Local Nusselt number  
- Grid convergence  
- Effect of convection schemes (UD vs QUICK)  
- Energy balance consistency  

The solver is fully developed in MATLAB using a structured mesh.

## Physical Model

- Geometry: Plane channel (L = 10 m, H = 1 m)
- Laminar Poiseuille flow
- Passive scalar temperature
- Governing equation:

\[
\nabla \cdot (\rho \mathbf{u} T) = \nabla \cdot (\Gamma \nabla T)
\]

- Finite Volume Method (FVM)
- Structured uniform grid

## Numerical Methods

### Diffusion
- Central Differencing (2nd order)

### Convection
- Upwind Differencing (UD)
- QUICK scheme (bounded deferred correction)

### Linear Solvers
- Direct solver (MATLAB `\`)
- Iterative SOR method

## Boundary Conditions

- Inlet: Dirichlet (fixed temperature)
- Outlet: Neumann (∂T/∂x = 0)
- Walls:
  - Isothermal (Dirichlet)
  - Constant heat flux (Neumann)
  - Localized heating segment

## Key Results

- Thermal entrance length ≈ 3.0–3.3 m (grid dependent)
- Fully developed Nusselt number:
  - Isothermal walls → Nu ≈ 7.54 (theoretical match on refined grid)
  - Constant heat flux → Nu ≈ 8.24
- QUICK reduces numerical diffusion compared to UD
- Grid refinement in y-direction is critical for accurate wall heat flux
- Global energy balance error < 2%

## Numerical Insights

- Coarse meshes underpredict Nusselt number due to poor wall gradient resolution.
- High Péclet numbers require finer grids for stability and accuracy.
- Bounded QUICK improves accuracy without sacrificing stability.
- Refining in the wall-normal direction is more important than refining in x.

## Tools

- MATLAB
- Sparse matrix assembly
- Finite Volume Method
- Structured grid implementation

## Report

A complete report is present in the Finite_Volume_Simulation_of_Laminar_Channel_Heat_Transfer.pdf file.
