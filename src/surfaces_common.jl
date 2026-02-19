# ---------------------------------------------------------------------------- #
#  Chapter 8: Construction of Common Surfaces — A8.1–A8.4                      #
# ---------------------------------------------------------------------------- #

"""
    make_bilinear_surface(P00::Vector{T}, P10::Vector{T},
                          P01::Vector{T}, P11::Vector{T}) -> BSplineSurface{T}

Construct a bilinear surface from four corner points.

Piegl & Tiller, *The NURBS Book*, 2nd ed., Algorithm A8.1, p. 334.

``S(u,v) = (1-u)(1-v)P_{00} + u(1-v)P_{10} + (1-u)v P_{01} + uv P_{11}``
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

Piegl & Tiller, *The NURBS Book*, 2nd ed., Algorithm A8.3, p. 337.

The curves must have the same degree and knot vector (after compatibility
processing, which this function performs automatically).
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

Piegl & Tiller, *The NURBS Book*, 2nd ed., Algorithm A8.1, p. 340.

# Arguments
- `crv`: generatrix curve
- `S`: point on axis
- `axis`: unit direction of axis
- `theta`: revolution angle (radians, up to 2π)
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
