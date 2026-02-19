# ---------------------------------------------------------------------------- #
#  Chapter 10: Advanced Surface Construction — A10.1–A10.4                     #
# ---------------------------------------------------------------------------- #

"""
    make_swung_surface(profile::NURBSCurve{T},
                       trajectory::NURBSCurve{T}) -> NURBSSurface{T}

Construct a swung surface from a profile (xz-plane) and trajectory (xy-plane).

A swung surface is a generalization of a surface of revolution where the
profile curve is scaled by the trajectory. The control points and weights
are (Eq. 10.7):

```math
P_{i,j} = (P_i^x \\cdot T_j^x,\\; P_i^x \\cdot T_j^y,\\; P_i^z),
\\qquad w_{i,j} = w_i^{\\text{prof}} \\cdot w_j^{\\text{traj}}
```

where ``P_i = (P_i^x, *, P_i^z)`` are profile control points and
``T_j = (T_j^x, T_j^y, *)`` are trajectory control points.

Piegl & Tiller, *The NURBS Book*, 2nd ed., Algorithm A10.1, p. 458.
"""
function make_swung_surface(profile::NURBSCurve{T}, trajectory::NURBSCurve{T}) where {T}
    n_prof = length(profile.controlpoints)
    n_traj = length(trajectory.controlpoints)
    dim = length(profile.controlpoints[1])

    Pij = Matrix{Vector{T}}(undef, n_prof, n_traj)
    wij = Matrix{T}(undef, n_prof, n_traj)

    for i in 1:n_prof, j in 1:n_traj
        Pi = profile.controlpoints[i]
        Tj = trajectory.controlpoints[j]
        Pij[i, j] = dim >= 3 ?
            [Pi[1] * Tj[1], Pi[1] * Tj[2], Pi[3]] :
            [Pi[1] * Tj[1], Pi[1] * Tj[2]]
        wij[i, j] = profile.weights[i] * trajectory.weights[j]
    end

    return NURBSSurface(profile.degree, trajectory.degree,
                        profile.knots, trajectory.knots, Pij, wij)
end

"""
    make_skinned_surface(curves::Vector{BSplineCurve{T}},
                         q::Int=3) -> BSplineSurface{T}

Construct a skinned surface through cross-section curves.

Given ``K`` compatible cross-section curves ``C_k(u)``, ``k = 0, \\ldots, K-1``,
the skinned surface interpolates them all (Eq. 10.19):

```math
S(u, v_k) = C_k(u), \\qquad k = 0, \\ldots, K - 1
```

The ``v``-direction interpolation is performed independently for each
``u``-direction control point index via global curve interpolation.
Cross-section curves are first made compatible (same degree and knot vector).

Piegl & Tiller, *The NURBS Book*, 2nd ed., Algorithm A10.2, p. 472.

# Arguments
- `curves`: B-spline cross-sections (at least 2)
- `q::Int`: v-direction degree (default 3)
"""
function make_skinned_surface(curves::Vector{BSplineCurve{T}}, q::Int=3) where {T}
    K = length(curves)
    K >= 2 || throw(ArgumentError("Need at least 2 curves"))

    # Elevate to common degree
    max_deg = maximum(c.degree for c in curves)
    elevated = [c.degree < max_deg ? degree_elevate(c, max_deg - c.degree) : c
                for c in curves]

    # Unify knot vectors
    all_knots = Set{T}()
    for c in elevated, k in c.knots
        push!(all_knots, k)
    end

    unified = BSplineCurve{T}[]
    for c in elevated
        extras = T[k for k in all_knots
                   if !any(abs(c.knots[j] - k) <= NURBS_EPSILON for j in 1:length(c.knots))]
        push!(unified, isempty(extras) ? c : refine_knots(c, sort(extras)))
    end

    n = length(unified[1].controlpoints)
    dim = length(unified[1].controlpoints[1])

    q_actual = min(q, K - 1)

    # Interpolate in v for each control point index
    v_params = [T(i - 1) / T(K - 1) for i in 1:K]
    v_knots = _compute_knot_vector(v_params, q_actual, K)

    final_P = Matrix{Vector{T}}(undef, n, K)
    for i in 1:n
        pts = [unified[j].controlpoints[i] for j in 1:K]
        interp = global_curve_interpolation(pts, q_actual; method=:chord)
        for j in 1:K
            final_P[i, j] = interp.controlpoints[j]
        end
    end

    return BSplineSurface(max_deg, q_actual, unified[1].knots, v_knots, final_P)
end

"""
    make_swept_surface(section::BSplineCurve{T},
                       trajectory::BSplineCurve{T}) -> BSplineSurface{T}

Construct a swept surface by translating a cross-section along a trajectory.

The swept surface is defined by (Eq. 10.26):

```math
S(u,v) = C_{\\text{section}}(u) + T(v) - T(0)
```

where ``T(v)`` is the trajectory curve. The implementation samples the
trajectory at multiple ``v``-values, translates the cross-section to each
position, and skins the resulting family of curves.

Piegl & Tiller, *The NURBS Book*, 2nd ed., Algorithm A10.3, p. 485.
"""
function make_swept_surface(section::BSplineCurve{T},
                            trajectory::BSplineCurve{T}) where {T}
    n_samples = max(length(trajectory.controlpoints), 5)
    v_params = range(trajectory.knots[1], trajectory.knots[end]; length=n_samples)

    sections = BSplineCurve{T}[]
    origin0 = section.controlpoints[1]
    for v in v_params
        origin = curve_point(trajectory, v)
        new_pts = [pt .+ origin .- origin0 for pt in section.controlpoints]
        push!(sections, BSplineCurve(section.degree, section.knots, new_pts))
    end

    return make_skinned_surface(sections, trajectory.degree)
end

"""
    make_gordon_surface(u_curves::Vector{BSplineCurve{T}},
                        v_curves::Vector{BSplineCurve{T}},
                        Q::Matrix{Vector{T}},
                        p::Int, q::Int) -> BSplineSurface{T}

Construct a Gordon surface from a curve network.

Given a network of ``u``-curves and ``v``-curves that intersect at points
``Q_{k,l}``, the Gordon surface is (Eq. 10.34):

```math
S(u,v) = S_u(u,v) + S_v(u,v) - S_{tp}(u,v)
```

where ``S_u`` is the skin through the ``u``-curves, ``S_v`` is the skin
through the ``v``-curves, and ``S_{tp}`` is the tensor-product surface
interpolating the intersection points.

Piegl & Tiller, *The NURBS Book*, 2nd ed., Algorithm A10.4, p. 496.
"""
function make_gordon_surface(u_curves::Vector{BSplineCurve{T}},
                             v_curves::Vector{BSplineCurve{T}},
                             Q::Matrix{Vector{T}}, p::Int, q::Int) where {T}
    Su = make_skinned_surface(u_curves, q)
    Sv = make_skinned_surface(v_curves, p)
    Stp = global_surface_interpolation(Q, p, q)

    nu = min(size(Su.controlpoints, 1), size(Sv.controlpoints, 1), size(Stp.controlpoints, 1))
    nv = min(size(Su.controlpoints, 2), size(Sv.controlpoints, 2), size(Stp.controlpoints, 2))

    P = [Su.controlpoints[i, j] .+ Sv.controlpoints[i, j] .- Stp.controlpoints[i, j]
         for i in 1:nu, j in 1:nv]

    return BSplineSurface(p, q, Su.uknots, Su.vknots, P)
end

"""
    make_coons_surface(c1::BSplineCurve{T}, c2::BSplineCurve{T},
                       c3::BSplineCurve{T}, c4::BSplineCurve{T}) -> BSplineSurface{T}

Construct a Coons patch from four boundary curves.

The Coons surface combines a ruled surface through the ``u``-boundaries, a
ruled surface through the ``v``-boundaries, and a bilinear correction
(Eq. 10.38):

```math
S(u,v) = S_u(u,v) + S_v(u,v) - S_{bl}(u,v)
```

where ``S_u`` is the ruled surface between ``c_3`` (``u=0``) and ``c_4``
(``u=1``), ``S_v`` is the ruled surface between ``c_1`` (``v=0``) and
``c_2`` (``v=1``), and ``S_{bl}`` is the bilinear surface through the
four corners.

Piegl & Tiller, *The NURBS Book*, 2nd ed., Section 10.6, p. 509.

Boundaries: `c1` (``v=0``), `c2` (``v=1``), `c3` (``u=0``), `c4` (``u=1``).
"""
function make_coons_surface(c1::BSplineCurve{T}, c2::BSplineCurve{T},
                            c3::BSplineCurve{T}, c4::BSplineCurve{T}) where {T}
    P00 = curve_point(c1, c1.knots[1])
    P10 = curve_point(c1, c1.knots[end])
    P01 = curve_point(c2, c2.knots[1])
    P11 = curve_point(c2, c2.knots[end])

    Srv = make_ruled_surface(c1, c2)
    Sru = make_ruled_surface(c3, c4)
    Sbl = make_bilinear_surface(P00, P10, P01, P11)

    nu = min(size(Srv.controlpoints, 1), size(Sru.controlpoints, 1))
    nv = min(size(Srv.controlpoints, 2), size(Sru.controlpoints, 2))

    P = Matrix{Vector{T}}(undef, nu, nv)
    for i in 1:nu, j in 1:nv
        pu = i <= size(Sru.controlpoints, 1) && j <= size(Sru.controlpoints, 2) ?
             Sru.controlpoints[i, j] : zeros(T, length(P00))
        pv = i <= size(Srv.controlpoints, 1) && j <= size(Srv.controlpoints, 2) ?
             Srv.controlpoints[i, j] : zeros(T, length(P00))
        pb = i <= size(Sbl.controlpoints, 1) && j <= size(Sbl.controlpoints, 2) ?
             Sbl.controlpoints[i, j] : zeros(T, length(P00))
        P[i, j] = pu .+ pv .- pb
    end

    return BSplineSurface(Srv.udegree, Srv.vdegree, Srv.uknots, Srv.vknots, P)
end
