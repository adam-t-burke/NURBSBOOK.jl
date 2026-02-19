# ---------------------------------------------------------------------------- #
#  Chapter 7: Conics and Circles — A7.1–A7.8                                   #
# ---------------------------------------------------------------------------- #

"""
    make_arc(O::Vector{T}, X::Vector{T}, Y::Vector{T}, r::T,
             ths::T, the::T) -> NURBSCurve{T}

Construct a circular arc as a NURBS curve.

Points on the arc are:

```math
P(\\theta) = O + r\\cos\\theta\\, \\mathbf{X} + r\\sin\\theta\\, \\mathbf{Y},
\\qquad \\theta_s \\le \\theta \\le \\theta_e
```

The arc is represented as one or more rational quadratic Bézier segments
(degree 2 NURBS). Each segment subtends at most ``\\pi/2`` radians, with
interior weights ``w = \\cos(\\Delta\\theta / 2)`` (Eq. 7.12). The number
of arcs is:
- 1 if ``\\theta_e - \\theta_s \\le \\pi/2``
- 2 if ``\\le \\pi``
- 3 if ``\\le 3\\pi/2``
- 4 otherwise

Implements **Algorithm A7.1** from Piegl & Tiller, *The NURBS Book*, 2nd ed., p. 308.

# Arguments
- `O`: center
- `X`, `Y`: orthonormal axes defining the plane of the arc
- `r`: radius
- `ths`: start angle (radians)
- `the`: end angle (radians)

# Returns
- `NURBSCurve{T}`: the circular arc (degree 2)
"""
function make_arc(O::Vector{T}, X::Vector{T}, Y::Vector{T}, r::T,
                  ths::T, the::T) where {T}
    theta = the - ths
    if theta < 0
        theta += 2 * T(π)
    end

    narcs = theta <= T(π)/2 ? 1 : theta <= T(π) ? 2 : theta <= 3*T(π)/2 ? 3 : 4
    dtheta = theta / narcs
    w1 = cos(dtheta / 2)

    n_pts = 2 * narcs + 1
    Pws = [zeros(T, length(O)) for _ in 1:n_pts]
    wts = ones(T, n_pts)

    P0 = O .+ r * cos(ths) .* X .+ r * sin(ths) .* Y
    Pws[1] = copy(P0)

    angle = ths
    idx = 1
    for _ in 1:narcs
        angle += dtheta
        P2 = O .+ r * cos(angle) .* X .+ r * sin(angle) .* Y
        mid_angle = angle - dtheta / 2
        P1 = O .+ (r / w1) .* (cos(mid_angle) .* X .+ sin(mid_angle) .* Y)

        Pws[idx + 1] = copy(P1)
        wts[idx + 1] = w1
        Pws[idx + 2] = copy(P2)

        idx += 2
    end

    # Build knot vector
    knot_vals = T[0]
    for i in 1:narcs
        push!(knot_vals, T(i) / narcs)
    end
    knots_full = T[]
    for _ in 1:3; push!(knots_full, knot_vals[1]); end
    for i in 2:narcs
        for _ in 1:2; push!(knots_full, knot_vals[i]); end
    end
    for _ in 1:3; push!(knots_full, knot_vals[end]); end

    return NURBSCurve(2, KnotVector(knots_full), Pws, wts)
end

"""
    make_circle(O::Vector{T}, X::Vector{T}, Y::Vector{T}, r::T) -> NURBSCurve{T}

Construct a full circle as a NURBS curve.

The circle is represented as four rational quadratic arcs, each subtending
``\\pi/2`` radians:

```math
\\{P(\\theta) : 0 \\le \\theta < 2\\pi\\}, \\qquad
P(\\theta) = O + r\\cos\\theta\\, \\mathbf{X} + r\\sin\\theta\\, \\mathbf{Y}
```

This is a convenience wrapper around [`make_arc`](@ref) spanning ``[0, 2\\pi)``.

Piegl & Tiller, *The NURBS Book*, 2nd ed., Section 7.5, p. 308.
"""
function make_circle(O::Vector{T}, X::Vector{T}, Y::Vector{T}, r::T) where {T}
    return make_arc(O, X, Y, r, zero(T), T(2π) - T(1e-10))
end

"""
    make_conic_arc(P0::Vector{T}, T0::Vector{T}, P2::Vector{T}, T2::Vector{T},
                   w::T) -> NURBSCurve{T}

Construct a conic arc as a rational quadratic Bézier curve.

The resulting curve is (Eq. 7.2):

```math
C(u) = \\frac{(1-u)^2 P_0 + 2u(1-u)w_1 P_1 + u^2 P_2}
             {(1-u)^2 + 2u(1-u)w_1 + u^2}
```

where ``P_1`` is the intersection of the tangent lines at ``P_0`` and
``P_2``, and ``w_1`` is the shoulder weight. The conic type is determined by:

- ``0 < w_1 < 1``: **ellipse**
- ``w_1 = 1``: **parabola**
- ``w_1 > 1``: **hyperbola**

Piegl & Tiller, *The NURBS Book*, 2nd ed., Section 7.3, p. 296.

# Arguments
- `P0`: start point
- `T0`: tangent direction at start
- `P2`: end point
- `T2`: tangent direction at end
- `w`: shoulder weight
"""
function make_conic_arc(P0::Vector{T}, T0::Vector{T}, P2::Vector{T}, T2::Vector{T},
                        w::T) where {T}
    # Intersect tangent lines to find P1
    denom = T0[1] * T2[2] - T0[2] * T2[1]
    if abs(denom) > NURBS_EPSILON
        s = ((P2[1] - P0[1]) * T2[2] - (P2[2] - P0[2]) * T2[1]) / denom
        P1 = P0 .+ s .* T0
    else
        P1 = (P0 .+ P2) ./ 2
    end
    knots = KnotVector(T[0, 0, 0, 1, 1, 1])
    return NURBSCurve(2, knots, [copy(P0), copy(P1), copy(P2)], T[1, w, 1])
end

"""
    make_one_arc(P0::Vector{T}, T0::Vector{T}, P2::Vector{T}, T2::Vector{T},
                 P::Vector{T}) -> NURBSCurve{T}

Construct a rational quadratic Bézier arc through ``P_0``, shoulder point
``P``, and ``P_2`` with prescribed tangent directions.

The control point ``P_1`` is the intersection of the tangent lines from
``P_0`` and ``P_2``. The weight ``w_1`` is computed so that the curve passes
through the shoulder point ``P`` at ``u = 1/2`` (Eq. 7.7):

```math
C\\!\\left(\\tfrac{1}{2}\\right) =
  \\frac{\\tfrac{1}{4} P_0 + \\tfrac{1}{2} w_1 P_1 + \\tfrac{1}{4} P_2}
       {\\tfrac{1}{2} + \\tfrac{1}{2} w_1} = P
```

Piegl & Tiller, *The NURBS Book*, 2nd ed., Section 7.3, p. 295.
"""
function make_one_arc(P0::Vector{T}, T0::Vector{T}, P2::Vector{T}, T2::Vector{T},
                      P::Vector{T}) where {T}
    # Find P1 as intersection of tangent lines
    denom = T0[1] * T2[2] - T0[2] * T2[1]
    if abs(denom) > NURBS_EPSILON
        s = ((P2[1] - P0[1]) * T2[2] - (P2[2] - P0[2]) * T2[1]) / denom
        P1 = P0 .+ s .* T0
    else
        P1 = (P0 .+ P2) ./ 2
    end

    # Compute weight from shoulder point
    mid = T(0.25) .* P0 .+ T(0.5) .* P1 .+ T(0.25) .* P2
    diff = P .- mid
    denom2 = dot(diff, P1 .- P)
    if abs(denom2) > NURBS_EPSILON
        w1 = dot(P .- T(0.5) .* (P0 .+ P2), P1 .- P) / denom2
    else
        w1 = one(T)
    end
    w1 = max(w1, T(NURBS_EPSILON))

    knots = KnotVector(T[0, 0, 0, 1, 1, 1])
    return NURBSCurve(2, knots, [copy(P0), copy(P1), copy(P2)], T[1, w1, 1])
end
