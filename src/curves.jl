# ---------------------------------------------------------------------------- #
#  Chapter 3: B-spline Curves — Algorithms A3.1–A3.4                           #
# ---------------------------------------------------------------------------- #

"""
    curve_point(crv::BSplineCurve{T}, u::Real) -> Vector{T}

Compute a point on a B-spline curve at parameter `u`.

Implements Algorithm A3.1 from Piegl & Tiller, *The NURBS Book*, 2nd ed., p. 82.

```math
C(u) = \\sum_{i=1}^{n} N_{i,p}(u)\\, P_i
```

# Arguments
- `crv::BSplineCurve{T}`: B-spline curve
- `u::Real`: parameter value

# Returns
- `Vector{T}`: the point ``C(u)``

See also: [`curve_derivatives`](@ref)
"""
function curve_point(crv::BSplineCurve{T}, u::Real) where {T}
    p = crv.degree
    U = crv.knots
    P = crv.controlpoints
    n = length(P)

    span = find_span(n, p, u, U)
    N = basis_functions(span, u, p, U)

    C = zeros(T, length(P[1]))
    for j in 1:(p + 1)
        C .+= N[j] .* P[span - p + j - 1]
    end
    return C
end

"""
    curve_derivatives(crv::BSplineCurve{T}, u::Real, d::Int) -> Vector{Vector{T}}

Compute derivatives of a B-spline curve up to order `d` at parameter `u`.

Implements Algorithm A3.2 from Piegl & Tiller, *The NURBS Book*, 2nd ed., p. 93.

# Arguments
- `crv::BSplineCurve{T}`: B-spline curve
- `u::Real`: parameter value
- `d::Int`: maximum derivative order

# Returns
- `Vector{Vector{T}}` of length `d + 1`. Entry `k` is the `(k-1)`-th derivative.
  Entry 1 is the point, entry 2 is the first derivative, etc.

See also: [`curve_point`](@ref), [`curve_deriv_control_points`](@ref)
"""
function curve_derivatives(crv::BSplineCurve{T}, u::Real, d::Int) where {T}
    p = crv.degree
    U = crv.knots
    P = crv.controlpoints
    n = length(P)
    dim = length(P[1])
    du = min(d, p)

    CK = [zeros(T, dim) for _ in 1:(d + 1)]

    span = find_span(n, p, u, U)
    nders = basis_function_derivatives(span, u, p, du, U)

    for k in 1:(du + 1)
        for j in 1:(p + 1)
            CK[k] .+= nders[k, j] .* P[span - p + j - 1]
        end
    end

    return CK
end

"""
    curve_deriv_control_points(n::Int, p::Int, U::KnotVector{T},
                               P::Vector{Vector{T}}, d::Int,
                               r1::Int, r2::Int) -> Vector{Vector{Vector{T}}}

Compute control points of all derivative curves up to order `d`, restricted
to control point indices `r1` through `r2` (1-based).

Implements Algorithm A3.3 from Piegl & Tiller, *The NURBS Book*, 2nd ed., p. 98.

# Arguments
- `n::Int`: number of control points
- `p::Int`: degree
- `U::KnotVector{T}`: knot vector
- `P::Vector{Vector{T}}`: control points
- `d::Int`: maximum derivative order
- `r1::Int`, `r2::Int`: range of control point indices (1-based)

# Returns
- Nested vector `PK` where `PK[k][i]` is the `i`-th control point
  of the `(k-1)`-th derivative curve. `k` ranges from 1 to `d + 1`.

See also: [`curve_derivatives`](@ref)
"""
function curve_deriv_control_points(n::Int, p::Int, U::KnotVector{T},
                                    P::Vector{Vector{T}}, d::Int,
                                    r1::Int, r2::Int) where {T}
    r = r2 - r1
    PK = Vector{Vector{Vector{T}}}(undef, d + 1)

    # Zeroth derivative = original control points
    PK[1] = [copy(P[r1 + i - 1]) for i in 1:(r + 1)]

    for k in 1:d
        tmp = r - k + 1
        pkvec = Vector{Vector{T}}(undef, tmp)
        for i in 1:tmp
            denom = U[r1 + i - 1 + p + 1] - U[r1 + i - 1 + k]
            pkvec[i] = T(p - k + 1) / denom .* (PK[k][i + 1] .- PK[k][i])
        end
        PK[k + 1] = pkvec
    end

    return PK
end
