# ---------------------------------------------------------------------------- #
#  Chapter 5: Degree Elevation & Reduction — A5.9, A5.11                       #
#                                                                              #
#  Uses Bezier decomposition for robustness:                                   #
#    1. Decompose B-spline into Bezier segments via knot insertion             #
#    2. Degree-elevate each Bezier segment                                     #
#    3. Reassemble into a single B-spline curve                                #
# ---------------------------------------------------------------------------- #

"""
    _bezier_degree_elevate(P::Vector{Vector{T}}, t::Int) -> Vector{Vector{T}}

Degree-elevate a single Bezier curve (given by control points `P` of degree
`length(P)-1`) by `t` degrees.

Uses the identity (Eq. 5.36 from *The NURBS Book*):
```math
Q_i = \\sum_{j=\\max(0,i-t)}^{\\min(p,i)} \\frac{\\binom{p}{j}\\binom{t}{i-j}}{\\binom{p+t}{i}} P_j
```
"""
function _bezier_degree_elevate(P::Vector{Vector{T}}, t::Int) where {T}
    p = length(P) - 1
    ph = p + t
    Q = [zeros(T, length(P[1])) for _ in 1:(ph + 1)]
    for i in 0:ph
        for j in max(0, i - t):min(p, i)
            coeff = T(binomial_coefficient(p, j)) * T(binomial_coefficient(t, i - j)) /
                    T(binomial_coefficient(ph, i))
            Q[i + 1] .+= coeff .* P[j + 1]
        end
    end
    return Q
end

"""
    degree_elevate(crv::BSplineCurve{T}, t::Int=1) -> BSplineCurve{T}

Elevate the degree of a B-spline curve by ``t`` without changing its shape.

The algorithm works in three steps:

1. **Decompose** the B-spline into Bézier segments via knot insertion
   (every interior knot is raised to multiplicity ``p``).
2. **Degree-elevate** each Bézier segment using (Eq. 5.36):
```math
Q_i = \\sum_{j=\\max(0,\\,i-t)}^{\\min(p,\\,i)}
  \\frac{\\binom{p}{j}\\,\\binom{t}{i-j}}{\\binom{p+t}{i}}\\, P_j
  \\qquad i = 0, \\ldots, p+t
```
3. **Reassemble** the elevated Bézier segments into a single
   ``(p+t)``-degree B-spline by removing shared endpoints.

Equivalent to **Algorithm A5.9** from Piegl & Tiller, *The NURBS Book*,
2nd ed., p. 206.

# Arguments
- `crv::BSplineCurve{T}`: input curve of degree ``p``
- `t::Int`: degree increment (default 1)

# Returns
- `BSplineCurve{T}`: new curve of degree ``p + t``

See also: [`degree_reduce`](@ref)
"""
function degree_elevate(crv::BSplineCurve{T}, t::Int=1) where {T}
    p = crv.degree
    U = crv.knots
    P = crv.controlpoints
    n = length(P)
    dim = length(P[1])
    t <= 0 && return crv

    ph = p + t

    # Step 1: Find distinct interior knots and their multiplicities
    interior_knots = T[]
    interior_mults = Int[]
    i = p + 2  # first interior knot position
    while i <= n  # last interior knot = U[n]
        val = U[i]
        mult = 1
        while i + mult <= n && abs(U[i + mult] - val) <= NURBS_EPSILON
            mult += 1
        end
        push!(interior_knots, val)
        push!(interior_mults, mult)
        i += mult
    end

    # Step 2: Insert knots to decompose into Bezier segments
    # Each interior knot needs multiplicity p
    knots_to_insert = T[]
    for (kval, kmult) in zip(interior_knots, interior_mults)
        for _ in 1:(p - kmult)
            push!(knots_to_insert, kval)
        end
    end

    if isempty(knots_to_insert)
        bezier_crv = crv
    else
        sort!(knots_to_insert)
        bezier_crv = refine_knots(crv, knots_to_insert)
    end

    Pb = bezier_crv.controlpoints
    Ub = bezier_crv.knots
    n_bezier = (length(Pb) - 1) ÷ p  # number of Bezier segments
    n_bezier >= 1 || return crv

    # Step 3: Degree-elevate each Bezier segment
    elevated_segs = Vector{Vector{Vector{T}}}(undef, n_bezier)
    for seg in 1:n_bezier
        start_idx = (seg - 1) * p + 1
        seg_pts = Pb[start_idx:(start_idx + p)]
        elevated_segs[seg] = _bezier_degree_elevate(seg_pts, t)
    end

    # Step 4: Reassemble — consecutive Bezier segments share an endpoint
    new_P = Vector{Vector{T}}()
    for seg in 1:n_bezier
        first_j = seg == 1 ? 1 : 2
        for j in first_j:(ph + 1)
            push!(new_P, elevated_segs[seg][j])
        end
    end

    # Step 5: Build knot vector — use Bezier-level multiplicities (ph at each
    # interior breakpoint), matching the control point count exactly.
    new_knots = T[]
    for _ in 1:(ph + 1)
        push!(new_knots, U[1])
    end
    for kval in interior_knots
        for _ in 1:ph
            push!(new_knots, kval)
        end
    end
    for _ in 1:(ph + 1)
        push!(new_knots, U[end])
    end

    return BSplineCurve(ph, KnotVector(new_knots), new_P)
end

"""
    degree_elevate(crv::NURBSCurve{T}, t::Int=1) -> NURBSCurve{T}

Elevate the degree of a NURBS curve by ``t``.

The NURBS curve is lifted to homogeneous coordinates
``P_i^w = (w_i P_i,\\, w_i)``, B-spline degree elevation is applied
(Eq. 5.36), and the result is projected back.

See also: [`degree_elevate(::BSplineCurve, ::Int)`](@ref)
"""
function degree_elevate(crv::NURBSCurve{T}, t::Int=1) where {T}
    Pw = _to_homogeneous(crv.controlpoints, crv.weights)
    bcrv = BSplineCurve(crv.degree, crv.knots, Pw)
    newcrv = degree_elevate(bcrv, t)
    new_pts, new_w = _from_homogeneous(newcrv.controlpoints)
    return NURBSCurve(newcrv.degree, newcrv.knots, new_pts, new_w)
end

"""
    degree_reduce(crv::BSplineCurve{T}; tol::Real=1e-6) -> Tuple{BSplineCurve{T}, T}

Approximate degree reduction of a B-spline curve by 1.

Degree reduction seeks a curve ``\\hat{C}(u)`` of degree ``p - 1`` that
approximates the original degree-``p`` curve ``C(u)``:

```math
\\max_{u} \\lVert C(u) - \\hat{C}(u) \\rVert \\le \\varepsilon
```

This implementation uses Greville abscissa sampling to place interior
control points:

```math
\\bar{u}_i = \\frac{1}{p-1} \\sum_{j=i+1}^{i+p-1} u_j
```

and reports the maximum deviation ``\\varepsilon`` over a dense parameter
sample.

Based on the degree reduction framework in Piegl & Tiller, *The NURBS Book*,
2nd ed., Section 5.6, p. 220.

# Arguments
- `crv::BSplineCurve{T}`: input curve of degree ``p \\ge 2``
- `tol::Real`: tolerance for error estimate

# Returns
- `(reduced_curve, max_error)`: the ``(p-1)``-degree curve and max deviation

See also: [`degree_elevate`](@ref)
"""
function degree_reduce(crv::BSplineCurve{T}; tol::Real=1e-6) where {T}
    p = crv.degree
    p >= 2 || throw(ArgumentError("Cannot reduce degree below 1"))

    U = crv.knots
    Pw = crv.controlpoints
    n = length(Pw)
    dim = length(Pw[1])
    ph = p - 1

    # Build reduced knot vector: decrease interior multiplicities by 1
    unique_knots = T[]
    mults = Int[]
    idx = 1
    m = length(U)
    while idx <= m
        val = U[idx]
        c = 0
        while idx + c <= m && abs(U[idx + c] - val) <= NURBS_EPSILON
            c += 1
        end
        push!(unique_knots, val)
        push!(mults, c)
        idx += c
    end

    new_knots = T[]
    for i in eachindex(unique_knots)
        mult = (i == 1 || i == length(unique_knots)) ? ph + 1 : max(mults[i] - 1, 1)
        for _ in 1:mult
            push!(new_knots, unique_knots[i])
        end
    end

    new_n = length(new_knots) - ph - 1
    new_pts = [zeros(T, dim) for _ in 1:new_n]
    new_pts[1] = copy(Pw[1])
    new_pts[new_n] = copy(Pw[n])

    # Approximate interior control points via Greville abscissae
    for i in 2:(new_n - 1)
        t_param = sum(new_knots[(i + 1):(i + ph)]) / ph
        new_pts[i] = curve_point(crv, t_param)
    end

    newcrv = BSplineCurve(ph, KnotVector(new_knots), new_pts)

    max_err = zero(T)
    for t_param in range(U[1], U[end]; length=100)
        err = norm(curve_point(crv, t_param) .- curve_point(newcrv, t_param))
        max_err = max(max_err, err)
    end

    return (newcrv, max_err)
end
