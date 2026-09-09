# Numerical Laplacians and Perfect Discretization

This repository contains the computational work from my undergraduate physics senior thesis on numerical discretizations of the Laplace operator.

The project compares standard 5-point and 9-point numerical Laplacians with a perfect Laplacian derived from the Klein–Gordon equation. The computations were implemented in MATLAB.

## Project Overview

Numerical solutions of partial differential equations require continuous differential operators to be replaced by discrete approximations that can be evaluated computationally.

For a scalar function $u(x,y)$, the two-dimensional Laplace operator is

```math
\Delta u
=
\frac{\partial^2 u}{\partial x^2}
+
\frac{\partial^2 u}{\partial y^2}.
```

This project investigates different discrete representations of this operator and compares conventional finite-difference Laplacians with a perfect-discretization approach.

The methods studied were:

- 5-point numerical Laplacian
- 9-point numerical Laplacian
- Perfect Laplacian derived from the Klein–Gordon equation

## Motivation

Discretization replaces a continuous mathematical model with a finite computational representation. Different discretizations can preserve the behavior of the underlying continuous system with different levels of accuracy and computational complexity.

The purpose of this project was to study how several representations of the Laplacian behave computationally and to explore the idea of a perfect discretization.

## Computational Work

MATLAB was used to implement the numerical operators and perform the computations used in the thesis.

The project involved:

- construction of discrete Laplacian operators;
- numerical comparison of multiple discretization schemes;
- implementation of the methods in MATLAB;
- analysis of the resulting numerical behavior.

## Mathematical and Computational Topics

- Numerical analysis
- Partial differential equations
- Finite-difference methods
- Laplace operator
- Klein–Gordon equation
- Perfect discretization
- Scientific computing
- MATLAB

## Repository Contents

This repository contains the MATLAB files used for the numerical computations presented in the senior thesis.

## Future Improvements

The original code was written as part of an undergraduate research project. Planned improvements to the repository include:

- reorganizing the MATLAB code into clearly documented components;
- adding instructions for reproducing the numerical experiments;
- adding figures illustrating the comparisons;
- documenting the mathematical derivation used in the project;
- adding comments and function documentation to the MATLAB source;
- reproducing selected computations in Python.

## Author

**Latimer Galvan Harris-Ward**

M.S. Mathematics  
B.S. Applied Mathematics  
B.S. Physics, Concentration in Mathematical Physics
