# Numerical Laplacians and Perfect Discretization

This repository contains the computational work from my undergraduate physics senior thesis on numerical discretizations of the Laplace operator.

I independently implemented and compared standard 5-point and 9-point numerical Laplacians in MATLAB. I also attempted to implement a perfect Laplacian based on the work of Simon Hauswirth.

## Project Overview

For a scalar function $u(x,y)$, the two-dimensional Laplace operator is

```math
\Delta u
=
\frac{\partial^2 u}{\partial x^2}
+
\frac{\partial^2 u}{\partial y^2}.
```

The project investigated three discrete representations of this operator:

- a 5-point finite-difference Laplacian;
- a 9-point finite-difference Laplacian;
- an attempted implementation of a perfect Laplacian following Simon Hauswirth's treatment of fixed-point lattice operators.

The 5-point and 9-point implementations, numerical experiments, and comparison code were my own work.

## Perfect Laplacian and Attribution

The mathematical construction of the perfect Laplacian was taken from the work of **Simon Hauswirth**.

I attempted to translate Hauswirth's formulation into MATLAB as part of my senior thesis. That implementation was my own code, but I did not reproduce Hauswirth's method successfully enough for it to be considered a faithful reference implementation.

Accordingly, this repository should distinguish between:

- my independent implementations of the standard numerical Laplacians; and
- my attempted implementation of Hauswirth's perfect-Laplacian construction.

The latter is best viewed as exploratory undergraduate work undertaken while learning the theory.

## Computational Work

MATLAB was used throughout the project.

My work included:

- implementing the 5-point Laplacian;
- implementing the 9-point Laplacian;
- constructing numerical experiments to compare the discretizations;
- analyzing the resulting behavior;
- attempting to implement Hauswirth's perfect-Laplacian formulation.

## Limitations

The standard finite-difference implementations were completed independently.

The perfect-Laplacian portion should be treated with caution. Although I wrote the MATLAB implementation myself, the underlying construction was Hauswirth's, and my implementation did not reproduce that construction as accurately as intended.

## Mathematical and Computational Topics

- Numerical analysis
- Partial differential equations
- Finite-difference methods
- Laplace operator
- Perfect discretization
- Scientific computing
- MATLAB

## Future Improvements

- Revisit Hauswirth's derivation in detail.
- Reimplement the perfect Laplacian more faithfully.
- Document the correspondence between the code and Hauswirth's equations.
- Add reproducible numerical experiments and figures.
- Quantify error for the 5-point and 9-point schemes.
- Compare those results with a corrected perfect-Laplacian implementation.
- Reproduce selected computations in Python.

## Author

**Latimer Galvan Harris-Ward**

M.S. Mathematics  
B.S. Applied Mathematics  
B.S. Physics, Concentration in Mathematical Physics
