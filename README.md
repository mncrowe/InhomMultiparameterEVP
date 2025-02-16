# Inhomogeneous Multiparameter Eigenvalue Problems

Consider the weird eigenvalue problem
```math
\left[ A - \sum_{i = 1}^n K_i B_i \right] a = c_0 + \sum_{i = 1}^n K_i c_i, \quad \textrm{where} \quad d_j^T a = 0 \quad \textrm{for} \quad j \in \{1, 2, \dots, n\},
```
where $A$ and the $B_i$ are $N \times N$ matrices, $a$, the $c_i$ and the $d_j$ are vectors in $R^N$ and the $K_i$ are eigenvalues. Solutions to this system consist of an n-tuple of eigenvalues, $(K_1, K_2, \dots, K_n)$, and an eigenvector, $a$.

This repository contains MATLAB scripts examining this system. We build up from simple cases and present both a linear algebra method and a root finding method for solving the problem. The linear algebra method (conversion to a linear system is of size $[(n+1)N]^n \times [(n+1)N]^n$) scales badly with $n$ and $N$ so is not feasible for large problems. The root finding method scales well and can be accelerated by analytically determining the Jacobian required for gradient-based methods.

This problem arises from the [N-layer quasi-geostrophic (QG) vortex problem](https://doi.org/10.1017/jfm.2024.619) and the numerical methods described here are implemented in [QGDipoles.jl](https://github.com/mncrowe/QGDipoles.jl) and [QGDipoles.m](https://github.com/mncrowe/QGDipoles.m).
