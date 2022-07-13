# SUNSET code

The **S**calable, **U**nstructured **N**ode-**SET** code, for mesh-free DNS.

Developed by Dr Jack King, University of Manchester.

## Overview and features

- Numerically solves the compressible Navier-Stokes equations.
- Domain discretised with unstructured node-set.
- Spatial discretisation between 4th and 10th order using **LABFM**.
- Temporal discretisation 3rd order explicit Runge-Kutta.
- Characteristic based boundary conditions
   + Walls (can be curved)
   + Inflow (subsonic, non-reflecting or hard)
   + Outflow (subsonic, non-reflecting)
   + Symmetry
   + Periodic
- Parallelised with OpenMP + MPI.
- Two-dimensional domain discretisation, tested up to 1024 cores.



