"""
    NURBS

A comprehensive Julia implementation of Non-Uniform Rational B-Splines based on
*The NURBS Book* by Les Piegl and Wayne Tiller (Springer, 2nd Edition, 1997).

All algorithms are faithfully translated from the book's pseudocode into idiomatic
Julia, with full docstrings referencing algorithm numbers, equations, and page numbers.

Only standard-library dependencies are used (`LinearAlgebra`).
"""
module NURBS

using LinearAlgebra

# Utility helpers (tolerances, binomial coefficients)
include("utils.jl")

# Core types: KnotVector, BSplineCurve, NURBSCurve, BSplineSurface, NURBSSurface
include("types.jl")

# Ch 2: B-spline basis functions (A2.1–A2.5)
include("basis.jl")

# Ch 3: B-spline curves (A3.1–A3.4)
include("curves.jl")

# Ch 3: B-spline surfaces (A3.5–A3.8)
include("surfaces.jl")

# Ch 4: NURBS (rational) curves & surfaces (A4.1–A4.5)
include("rational.jl")

# Ch 5: Knot insertion, refinement, removal (A5.1–A5.4)
include("knots.jl")

# Ch 5: Degree elevation & reduction (A5.5–A5.11)
include("degree.jl")

# Ch 6: Advanced geometric algorithms (A6.1–A6.5)
include("advanced.jl")

# Ch 7: Conics & circles (A7.1–A7.8)
include("conics.jl")

# Ch 8: Common surface constructions (A8.1–A8.4)
include("surfaces_common.jl")

# Ch 9: Curve & surface fitting (A9.1–A9.8)
include("fitting.jl")

# Ch 10: Advanced surface construction (A10.1–A10.4)
include("surfaces_advanced.jl")

# Ch 11: Shape modification tools
include("modification.jl")

# ---- Exports ------------------------------------------------------------- #

# Types
export KnotVector, BSplineCurve, NURBSCurve, BSplineSurface, NURBSSurface
export AbstractSplineCurve, AbstractSplineSurface

# Basis functions (Ch 2)
export find_span, basis_functions, basis_function_derivatives
export one_basis_function, one_basis_function_derivatives
export all_basis_functions

# Curves (Ch 3)
export curve_point, curve_derivatives, curve_deriv_control_points

# Surfaces (Ch 3)
export surface_point, surface_derivatives, surface_deriv_control_points

# Rational (Ch 4) — uses same names via dispatch

# Knot operations (Ch 5)
export insert_knot, refine_knots, remove_knot

# Degree operations (Ch 5)
export degree_elevate, degree_reduce

# Advanced algorithms (Ch 6)
export point_inversion_curve, point_inversion_surface
export curve_reverse, surface_reverse

# Conics & circles (Ch 7)
export make_circle, make_arc, make_conic_arc, make_one_arc

# Common surfaces (Ch 8)
export make_bilinear_surface, make_cylinder, make_ruled_surface
export make_revolved_surface

# Fitting (Ch 9)
export global_curve_interpolation, global_curve_approximation
export global_surface_interpolation
export local_curve_interpolation

# Advanced surfaces (Ch 10)
export make_swung_surface, make_skinned_surface
export make_swept_surface, make_gordon_surface, make_coons_surface

# Shape modification (Ch 11)
export move_control_point, modify_weight

end # module
