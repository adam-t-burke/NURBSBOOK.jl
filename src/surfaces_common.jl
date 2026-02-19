# ---------------------------------------------------------------------------- #
#  Chapter 8: Construction of Common Surfaces — A8.1–A8.4                      #
# ---------------------------------------------------------------------------- #

"""
    make_bilinear_surface(P00::Vector{T}, P10::Vector{T},
                          P01::Vector{T}, P11::Vector{T}) -> BSplineSurface{T}

Construct a bilinear surface from four corner points.

The bilinear surface is the simplest tensor-product surface (Eq. 8.1):

```math
S(u,v) = (1-u)(1-v)\\, P_{00} + u(1-v)\\, P_{10}
       + (1-u)v\\, P_{01} + uv\\, P_{11}
```

This is a degree ``(1, 1)`` B-spline surface with knot vectors
``\\{0, 0, 1, 1\\}`` in both directions.

Piegl & Tiller, *The NURBS Book*, 2nd ed., Algorithm A8.1, p. 334.
"""
function make_bilinear_surface(P00::Vector{T}, P10::Vector{T},
                               P01::Vector{T}, P11::Vector{T}) where {T}
    uk = KnotVector(T[0, 0, 1, 1])
    vk = KnotVector(T[0, 0, 1, 1])
    P = Matrix{Vector{T}}(undef, 2, 2)
    P[1, 1] = copy(P00); P[2, 1] = copy(P10)
    P[1, 2] = copy(P01); P[2, 2] = copy(P11)
    return BSplineSurface(1, 1, uk, vk, P)
end

"""
    make_ruled_surface(crv1::BSplineCurve{T},
                       crv2::BSplineCurve{T}) -> BSplineSurface{T}

Construct a ruled surface linearly interpolating between two B-spline curves.

A ruled surface through ``C_1(u)`` and ``C_2(u)`` is (Eq. 8.7):

```math
S(u,v) = (1 - v)\\, C_1(u) + v\\, C_2(u)
```

The two curves are made compatible (same degree and knot vector) via degree
elevation and knot refinement before the control net is constructed.

Piegl & Tiller, *The NURBS Book*, 2nd ed., Algorithm A8.3, p. 337.
"""
function make_ruled_surface(crv1::BSplineCurve{T}, crv2::BSplineCurve{T}) where {T}
    c1, c2 = crv1, crv2

    # Match degrees
    if c1.degree < c2.degree
        c1 = degree_elevate(c1, c2.degree - c1.degree)
    elseif c2.degree < c1.degree
        c2 = degree_elevate(c2, c1.degree - c2.degree)
    end

    # Unify knot vectors
    if length(c1.controlpoints) != length(c2.controlpoints)
        extras1 = T[]
        for i in 1:length(c2.knots)
            if !any(abs(c1.knots[j] - c2.knots[i]) <= NURBS_EPSILON for j in 1:length(c1.knots))
                push!(extras1, c2.knots[i])
            end
        end
        extras2 = T[]
        for i in 1:length(c1.knots)
            if !any(abs(c2.knots[j] - c1.knots[i]) <= NURBS_EPSILON for j in 1:length(c2.knots))
                push!(extras2, c1.knots[i])
            end
        end
        isempty(extras1) || (c1 = refine_knots(c1, sort(extras1)))
        isempty(extras2) || (c2 = refine_knots(c2, sort(extras2)))
    end

    n = length(c1.controlpoints)
    vk = KnotVector(T[0, 0, 1, 1])
    P = Matrix{Vector{T}}(undef, n, 2)
    for i in 1:n
        P[i, 1] = copy(c1.controlpoints[i])
        P[i, 2] = copy(c2.controlpoints[i])
    end
    return BSplineSurface(c1.degree, 1, c1.knots, vk, P)
end

"""
    make_cylinder(crv::BSplineCurve{T}, dir::Vector{T},
                  len::T) -> BSplineSurface{T}

Construct a cylinder by extruding a B-spline curve along a direction.

The extruded (translational) surface is (Eq. 8.3):

```math
S(u,v) = C(u) + v\\, L\\, \\hat{T}
```

where ``\\hat{T}`` is the extrusion direction and ``L`` is the length. This
is implemented as a ruled surface between ``C(u)`` and
``C(u) + L\\hat{T}``.

Piegl & Tiller, *The NURBS Book*, 2nd ed., Algorithm A8.2, p. 336.
"""
function make_cylinder(crv::BSplineCurve{T}, dir::Vector{T}, len::T) where {T}
    offset = len .* dir
    crv2 = BSplineCurve(crv.degree, crv.knots,
                        [pt .+ offset for pt in crv.controlpoints])
    return make_ruled_surface(crv, crv2)
end

"""
    make_revolved_surface(crv::NURBSCurve{T}, S::Vector{T}, axis::Vector{T},
                          theta::T) -> NURBSSurface{T}

Construct a surface of revolution by rotating a NURBS profile about an axis.

Each control point ``P_i`` of the generatrix is projected onto the axis,
yielding center ``O_i`` and radius ``r_i = \\lVert P_i - O_i \\rVert``.
A circular arc of angle ``\\theta`` is constructed for each point, and the
surface is assembled as:

```math
S(u, \\theta) = \\sum_{i=0}^{n} \\sum_{j=0}^{m}
  N_{i,p}(u)\\, N_{j,2}(\\theta)\\, w_{i,j}\\, P_{i,j}
```

where the ``v``-direction is degree 2 with weights derived from
``w_j = \\cos(\\Delta\\theta / 2)`` for interior arc control points.

Piegl & Tiller, *The NURBS Book*, 2nd ed., Algorithm A8.1, p. 340.

# Arguments
- `crv`: generatrix (profile) NURBS curve
- `S`: a point on the axis of revolution
- `axis`: unit direction vector of the axis
- `theta`: revolution angle in radians (up to ``2\\pi``)
"""
function make_revolved_surface(crv::NURBSCurve{T}, S::Vector{T}, axis::Vector{T},
                               theta::T) where {T}
    narcs = theta <= T(π)/2 ? 1 : theta <= T(π) ? 2 : theta <= 3*T(π)/2 ? 3 : 4
    dtheta = theta / narcs
    wm = cos(dtheta / 2)
    n_prof = length(crv.controlpoints)
    n_v = 2 * narcs + 1

    # V knot vector
    vk_arr = T[]
    for _ in 1:3; push!(vk_arr, zero(T)); end
    for i in 1:(narcs - 1)
        val = T(i) / narcs
        push!(vk_arr, val); push!(vk_arr, val)
    end
    for _ in 1:3; push!(vk_arr, one(T)); end

    Pij = Matrix{Vector{T}}(undef, n_prof, n_v)
    wij = ones(T, n_prof, n_v)

    for i in 1:n_prof
        Pi = crv.controlpoints[i]
        wi = crv.weights[i]
        O_proj = S .+ dot(Pi .- S, axis) .* axis
        X_dir = Pi .- O_proj
        r_i = norm(X_dir)

        Pij[i, 1] = copy(Pi)
        wij[i, 1] = wi

        if r_i < NURBS_EPSILON
            for j in 1:n_v
                Pij[i, j] = copy(Pi)
                wij[i, j] = wi
            end
        else
            X_hat = X_dir ./ r_i
            Y_hat = cross(axis, X_hat)
            idx = 1
            for arc in 1:narcs
                mid_angle = (arc - T(0.5)) * dtheta
                end_angle = arc * dtheta
                P2 = O_proj .+ r_i .* (cos(end_angle) .* X_hat .+ sin(end_angle) .* Y_hat)
                P1 = O_proj .+ (r_i / wm) .* (cos(mid_angle) .* X_hat .+ sin(mid_angle) .* Y_hat)

                Pij[i, idx + 1] = P1
                wij[i, idx + 1] = wm * wi
                Pij[i, idx + 2] = P2
                wij[i, idx + 2] = wi
                idx += 2
            end
        end
    end

    return NURBSSurface(crv.degree, 2, crv.knots, KnotVector(vk_arr), Pij, wij)
end
