# NURBSBOOK.jl

[![Build Status](https://github.com/adam-t-burke/NURBSBOOK.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/adam-t-burke/NURBSBOOK.jl/actions/workflows/CI.yml)
[![Dev Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://adam-t-burke.github.io/NURBSBOOK.jl/dev/)
[![Stable Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://adam-t-burke.github.io/NURBSBOOK.jl/stable/)

A comprehensive Julia implementation of Non-Uniform Rational B-Splines (NURBS) based on
*The NURBS Book* by Les Piegl and Wayne Tiller (Springer, 2nd Edition, 1997).

## Features

- **Complete coverage** of algorithms A2.1 through A10.4 from *The NURBS Book*
- **Zero non-standard dependencies** — only `LinearAlgebra` from Julia's standard library
- **Parametric types** supporting `Float32`, `Float64`, and `BigFloat`
- **Comprehensive docstrings** with references to book algorithms, equations, and pages
- **Documenter.jl-ready** documentation

## Installation

```julia
using Pkg
Pkg.add("NURBSBOOK")
```

## Quick Start

```julia
using NURBSBOOK

# Define a cubic B-spline curve
knots = KnotVector([0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 5.0, 5.0, 5.0])
ctrl_pts = [
    [0.0, 0.0], [1.0, 2.0], [3.0, 3.0], [5.0, 2.0],
    [6.0, 0.0], [7.0, -1.0], [8.0, 0.0], [9.0, 1.0],
]
curve = BSplineCurve(3, knots, ctrl_pts)

# Evaluate a point on the curve
pt = curve_point(curve, 2.5)
```

## Algorithms Implemented

| Chapter | Topic | Algorithms |
|---------|-------|------------|
| 2 | B-spline Basis Functions | A2.1–A2.5 |
| 3 | B-spline Curves & Surfaces | A3.1–A3.8 |
| 4 | NURBS Curves & Surfaces | A4.1–A4.5 |
| 5 | Fundamental Geometric Algorithms | A5.1–A5.11 |
| 6 | Advanced Geometric Algorithms | A6.1–A6.5 |
| 7 | Conics & Circles | A7.1–A7.8 |
| 8 | Common Surfaces | A8.1–A8.4 |
| 9 | Curve & Surface Fitting | A9.1–A9.8 |
| 10 | Advanced Surface Construction | A10.1–A10.4 |
| 11 | Shape Modification Tools | — |

## Reference

> L. Piegl and W. Tiller, *The NURBS Book*, 2nd ed., Monographs in Visual Communication.
> Berlin, Heidelberg: Springer, 1997.

## License

MIT
