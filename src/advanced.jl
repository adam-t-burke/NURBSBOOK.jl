# ---------------------------------------------------------------------------- #
#  Chapter 6: Advanced Geometric Algorithms — A6.1–A6.5                        #
# ---------------------------------------------------------------------------- #

"""
    point_inversion_curve(crv::AbstractSplineCurve{T}, P::AbstractVector{<:Real};
                          u0::Real=NaN, tol::Real=1e-8, maxiter::Int=50) -> T

Newton iteration to find parameter ``u`` such that ``C(u) \\approx P``.

Given a target point ``P``, this solves the minimization problem

```math
\\min_u \\lVert C(u) - P \\rVert^2
```

via Newton's method. Setting ``f(u) = C'(u) \\cdot (C(u) - P)``, the update
is (Eq. 6.3–6.4):

```math
u_{i+1} = u_i - \\frac{C'(u_i) \\cdot (C(u_i) - P)}
  {C''(u_i) \\cdot (C(u_i) - P) + \\lVert C'(u_i) \\rVert^2}
```

Convergence is checked on both ``\\lVert C(u) - P \\rVert`` and
``|\\Delta u|``.

Piegl & Tiller, *The NURBS Book*, 2nd ed., Section 6.1, p. 230.

# Arguments
- `crv`: spline curve
- `P`: target point
- `u0`: initial guess (default: domain midpoint)
- `tol`: convergence tolerance
- `maxiter`: max iterations

# Returns
- Parameter value ``u``
"""
function point_inversion_curve(crv::AbstractSplineCurve{T}, P::AbstractVector{<:Real};
                               u0::Real=NaN, tol::Real=1e-8, maxiter::Int=50) where {T}
    U = crv.knots
    a, b = U[1], U[end]
    u = isnan(u0) ? (a + b) / 2 : T(u0)
    u = clamp(u, a, b)

    for _ in 1:maxiter
        CK = curve_derivatives(crv, u, 2)
        diff = CK[1] .- P
        norm(diff) <= tol && return u

        denom = dot(CK[2], CK[2]) + dot(diff, CK[3])
        abs(denom) < NURBS_EPSILON && break

        du = dot(CK[2], diff) / denom
        u_new = clamp(u - du, a, b)
        abs(du) <= tol && return u_new
        u = u_new
    end
    return u
end

"""
    point_inversion_surface(surf::AbstractSplineSurface{T}, P::AbstractVector{<:Real};
                            u0::Real=NaN, v0::Real=NaN,
                            tol::Real=1e-8, maxiter::Int=50) -> Tuple{T, T}

Newton iteration to find ``(u, v)`` such that ``S(u,v) \\approx P``.

Solves ``\\min_{u,v} \\lVert S(u,v) - P \\rVert^2`` via the 2D Newton update:

```math
\\begin{pmatrix} \\Delta u \\\\ \\Delta v \\end{pmatrix}
= -J^{-1}
\\begin{pmatrix} S_u \\cdot (S - P) \\\\ S_v \\cdot (S - P) \\end{pmatrix}
```

where ``J`` is the ``2 \\times 2`` Hessian of the squared-distance objective:

```math
J = \\begin{pmatrix}
  S_u \\cdot S_u + (S - P) \\cdot S_{uu} & S_u \\cdot S_v + (S - P) \\cdot S_{uv} \\\\
  S_u \\cdot S_v + (S - P) \\cdot S_{uv} & S_v \\cdot S_v + (S - P) \\cdot S_{vv}
\\end{pmatrix}
```

Piegl & Tiller, *The NURBS Book*, 2nd ed., Section 6.1, p. 232.

# Arguments
- `surf`: spline surface
- `P`: target point
- `u0`, `v0`: initial guess (default: domain midpoints)
- `tol`: convergence tolerance
- `maxiter`: max iterations

# Returns
- ``(u, v)`` parameter pair
"""
function point_inversion_surface(surf::AbstractSplineSurface{T}, P::AbstractVector{<:Real};
                                 u0::Real=NaN, v0::Real=NaN,
                                 tol::Real=1e-8, maxiter::Int=50) where {T}
    Uk, Vk = surf.uknots, surf.vknots
    ua, ub = Uk[1], Uk[end]
    va, vb = Vk[1], Vk[end]

    u = isnan(u0) ? (ua + ub) / 2 : T(u0)
    v = isnan(v0) ? (va + vb) / 2 : T(v0)
    u, v = clamp(u, ua, ub), clamp(v, va, vb)

    for _ in 1:maxiter
        SKL = surface_derivatives(surf, u, v, 2)
        diff = SKL[1, 1] .- P
        norm(diff) <= tol && return (u, v)

        Su, Sv = SKL[2, 1], SKL[1, 2]
        Suu, Suv, Svv = SKL[3, 1], SKL[2, 2], SKL[1, 3]

        J11 = dot(Su, Su) + dot(diff, Suu)
        J12 = dot(Su, Sv) + dot(diff, Suv)
        J22 = dot(Sv, Sv) + dot(diff, Svv)
        det_J = J11 * J22 - J12 * J12
        abs(det_J) < NURBS_EPSILON && break

        fu, fv = dot(Su, diff), dot(Sv, diff)
        du = (J22 * fu - J12 * fv) / det_J
        dv = (J11 * fv - J12 * fu) / det_J

        u_new = clamp(u - du, ua, ub)
        v_new = clamp(v - dv, va, vb)
        (abs(du) + abs(dv)) <= tol && return (u_new, v_new)
        u, v = u_new, v_new
    end
    return (u, v)
end

"""
    curve_reverse(crv::BSplineCurve{T}) -> BSplineCurve{T}

Reverse the parameterization of a B-spline curve.

Constructs ``\\bar{C}(u) = C(a + b - u)`` where ``[a, b]`` is the parameter
domain. The new knots and control points are (Eq. 6.48):

```math
\\bar{u}_i = a + b - u_{m-i}, \\qquad \\bar{P}_i = P_{n-i}
```

Piegl & Tiller, *The NURBS Book*, 2nd ed., Section 6.5, p. 263.
"""
function curve_reverse(crv::BSplineCurve{T}) where {T}
    U = crv.knots
    a, b = U[1], U[end]
    new_knots = [a + b - U[length(U) + 1 - i] for i in 1:length(U)]
    return BSplineCurve(crv.degree, KnotVector(new_knots),
                        [copy(pt) for pt in reverse(crv.controlpoints)])
end

"""
    curve_reverse(crv::NURBSCurve{T}) -> NURBSCurve{T}

Reverse the parameterization of a NURBS curve.

Applies the reversal formulae ``\\bar{u}_i = a + b - u_{m-i}``,
``\\bar{P}_i = P_{n-i}``, ``\\bar{w}_i = w_{n-i}`` (Eq. 6.48).
"""
function curve_reverse(crv::NURBSCurve{T}) where {T}
    U = crv.knots
    a, b = U[1], U[end]
    new_knots = [a + b - U[length(U) + 1 - i] for i in 1:length(U)]
    return NURBSCurve(crv.degree, KnotVector(new_knots),
                      [copy(pt) for pt in reverse(crv.controlpoints)],
                      reverse(copy(crv.weights)))
end

"""
    surface_reverse(surf::BSplineSurface{T}, dir::Symbol) -> BSplineSurface{T}

Reverse parameterization of a B-spline surface in direction `:u` or `:v`.

For reversal in the ``u``-direction (Eq. 6.49):

```math
\\bar{u}_i = a + b - u_{r-i}, \\qquad \\bar{P}_{i,j} = P_{n-i,\\,j}
```

and analogously for ``v``.

Piegl & Tiller, *The NURBS Book*, 2nd ed., Section 6.5, p. 264.
"""
function surface_reverse(surf::BSplineSurface{T}, dir::Symbol) where {T}
    P = surf.controlpoints
    nu, nv = size(P)

    if dir == :u
        a, b = surf.uknots[1], surf.uknots[end]
        new_uk = KnotVector([a + b - surf.uknots[length(surf.uknots) + 1 - i]
                             for i in 1:length(surf.uknots)])
        new_P = [P[nu + 1 - i, j] for i in 1:nu, j in 1:nv]
        return BSplineSurface(surf.udegree, surf.vdegree, new_uk, surf.vknots, new_P)
    elseif dir == :v
        a, b = surf.vknots[1], surf.vknots[end]
        new_vk = KnotVector([a + b - surf.vknots[length(surf.vknots) + 1 - i]
                             for i in 1:length(surf.vknots)])
        new_P = [P[i, nv + 1 - j] for i in 1:nu, j in 1:nv]
        return BSplineSurface(surf.udegree, surf.vdegree, surf.uknots, new_vk, new_P)
    else
        throw(ArgumentError("direction must be :u or :v, got :$dir"))
    end
end
