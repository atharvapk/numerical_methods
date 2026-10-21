# Numerical Methods

C++ implementations and numerical investigations of iterative methods for solving linear systems.

## Iterative Methods for Linear Systems

This project began with a numerical analysis exam problem involving the Jacobi and Gauss-Seidel methods. I extended the original calculations into a computational investigation of their convergence behavior.

The investigation covers:

* Implementing Jacobi and Gauss-Seidel solvers in C++.
* Comparing convergence rates and residual norms.
* Constructing iteration matrices and examining their spectral radii.
* Using power iteration and the Rayleigh quotient to estimate dominant eigenvalues.
* Investigating the limitations of ordinary power iteration when the dominant eigenvalues form a complex-conjugate pair.

For the example system, Jacobi converged in 13 iterations and Gauss-Seidel in 9, using the same initial guess and stopping tolerance.

## Full Technical Report

**[Read the report (PDF)](Iterative_Methods_Report.pdf)**

The report contains the mathematical derivations, numerical results, convergence plots, and selected C++ implementations.

## Implementation

The repository contains the C++ source code used for the numerical experiments.
