# NURBS.jl

A comprehensive Julia implementation of Non-Uniform Rational B-Splines based on
*The NURBS Book* by Les Piegl and Wayne Tiller (Springer, 2nd Edition, 1997).

## Features

- Complete coverage of algorithms A2.1 through A10.4
- Zero non-standard dependencies (only `LinearAlgebra`)
- Fully 1-based Julia-idiomatic indexing
- Parametric types supporting `Float32`, `Float64`, and `BigFloat`
- Comprehensive docstrings referencing book algorithms, equations, and pages

## Getting Started

```julia
using NURBS

# Define a cubic B-spline curve
knots = KnotVector([0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.0])
P = [[0.0, 0.0], [1.0, 2.0], [3.0, 3.0], [5.0, 2.0], [6.0, 0.0], [7.0, -1.0]]
crv = BSplineCurve(3, knots, P)

# Evaluate
pt = curve_point(crv, 1.5)
```

## Reference

> L. Piegl and W. Tiller, *The NURBS Book*, 2nd ed., Monographs in Visual
> Communication. Berlin, Heidelberg: Springer, 1997.
"""
