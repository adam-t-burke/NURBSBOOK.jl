# ---------------------------------------------------------------------------- #
#  Chapter 3: B-spline Surfaces — Algorithms A3.5–A3.8                         #
# ---------------------------------------------------------------------------- #

"""
    surface_point(surf::BSplineSurface{T}, u::Real, v::Real) -> Vector{T}

Compute a point on a B-spline surface at parameters `(u, v)`.

Implements Algorithm A3.5 from Piegl & Tiller, *The NURBS Book*, 2nd ed., p. 103.

# Arguments
- `surf::BSplineSurface{T}`: B-spline surface
- `u::Real`, `v::Real`: parameter values

# Returns
- `Vector{T}`: the point ``S(u,v)``

See also: [`surface_derivatives`](@ref)
"""
function surface_point(surf::BSplineSurface{T}, u::Real, v::Real) where {T}
    p = surf.udegree
    q = surf.vdegree
    Uk = surf.uknots
    Vk = surf.vknots
    P = surf.controlpoints
    nu = size(P, 1)
    nv = size(P, 2)
    dim = length(P[1, 1])

    uspan = find_span(nu, p, u, Uk)
    Nu = basis_functions(uspan, u, p, Uk)
    vspan = find_span(nv, q, v, Vk)
    Nv = basis_functions(vspan, v, q, Vk)

    S = zeros(T, dim)
    for l in 1:(q + 1)
        temp = zeros(T, dim)
        for k in 1:(p + 1)
            temp .+= Nu[k] .* P[uspan - p + k - 1, vspan - q + l - 1]
        end
        S .+= Nv[l] .* temp
    end
    return S
end

"""
    surface_derivatives(surf::BSplineSurface{T}, u::Real, v::Real,
                        d::Int) -> Matrix{Vector{T}}

Compute all partial derivatives of a B-spline surface up to total order `d`.

Implements Algorithm A3.6 from Piegl & Tiller, *The NURBS Book*, 2nd ed., p. 111.

# Arguments
- `surf::BSplineSurface{T}`: B-spline surface
- `u::Real`, `v::Real`: parameter values
- `d::Int`: maximum total derivative order

# Returns
- `Matrix{Vector{T}}` of size `(d+1) × (d+1)`. Entry `[k, l]` is
  ``\\frac{\\partial^{k+l-2}}{\\partial u^{k-1} \\partial v^{l-1}} S(u,v)``.

See also: [`surface_point`](@ref)
"""
function surface_derivatives(surf::BSplineSurface{T}, u::Real, v::Real,
                             d::Int) where {T}
    p = surf.udegree
    q = surf.vdegree
    Uk = surf.uknots
    Vk = surf.vknots
    P = surf.controlpoints
    nu = size(P, 1)
    nv = size(P, 2)
    dim = length(P[1, 1])

    du = min(d, p)
    dv = min(d, q)

    SKL = [zeros(T, dim) for _ in 1:(d + 1), _ in 1:(d + 1)]

    uspan = find_span(nu, p, u, Uk)
    Nu = basis_function_derivatives(uspan, u, p, du, Uk)
    vspan = find_span(nv, q, v, Vk)
    Nv = basis_function_derivatives(vspan, v, q, dv, Vk)

    for k in 1:(du + 1)
        temp = [zeros(T, dim) for _ in 1:(q + 1)]
        for s in 1:(q + 1)
            for r in 1:(p + 1)
                temp[s] .+= Nu[k, r] .* P[uspan - p + r - 1, vspan - q + s - 1]
            end
        end
        dd = min(d - k + 1, dv)
        for l in 1:(dd + 1)
            for s in 1:(q + 1)
                SKL[k, l] .+= Nv[l, s] .* temp[s]
            end
        end
    end

    return SKL
end

"""
    surface_deriv_control_points(surf::BSplineSurface{T}, d::Int,
                                 r1::Int, r2::Int, s1::Int, s2::Int)

Compute control points of partial derivative surfaces up to order `d`.

Implements Algorithm A3.7 from Piegl & Tiller, *The NURBS Book*, 2nd ed., p. 114.

# Arguments
- `surf::BSplineSurface{T}`: B-spline surface
- `d::Int`: max derivative order
- `r1, r2`: u-direction control point range (1-based)
- `s1, s2`: v-direction control point range (1-based)

# Returns
- Nested structure `PKL[k][l]` as a matrix of control point vectors.
"""
function surface_deriv_control_points(surf::BSplineSurface{T}, d::Int,
                                      r1::Int, r2::Int, s1::Int, s2::Int) where {T}
    p = surf.udegree
    q = surf.vdegree
    Uk = surf.uknots
    Vk = surf.vknots
    P = surf.controlpoints

    r = r2 - r1
    s = s2 - s1

    PKL = [[Matrix{Vector{T}}(undef, 0, 0) for _ in 1:(d + 1)] for _ in 1:(d + 1)]

    # PKL[1][1] = original control points in the range
    PKL[1][1] = [copy(P[r1 + i - 1, s1 + j - 1]) for i in 1:(r + 1), j in 1:(s + 1)]

    # Differentiate in u
    for k in 1:min(d, p)
        rows = r - k + 1
        cols = s + 1
        tmp = Matrix{Vector{T}}(undef, rows, cols)
        for j in 1:cols
            for i in 1:rows
                denom = Uk[r1 + i - 1 + p + 1] - Uk[r1 + i - 1 + k]
                tmp[i, j] = T(p - k + 1) / denom .* (PKL[k][1][i + 1, j] .- PKL[k][1][i, j])
            end
        end
        PKL[k + 1][1] = tmp
    end

    # Differentiate in v
    for k in 1:min(d, p)
        dd = min(d - k, q)
        for l in 1:dd
            prev = PKL[k + 1][l]
            rows = size(prev, 1)
            cols = size(prev, 2) - 1
            tmp = Matrix{Vector{T}}(undef, rows, cols)
            for i in 1:rows
                for j in 1:cols
                    denom = Vk[s1 + j - 1 + q + 1] - Vk[s1 + j - 1 + l]
                    tmp[i, j] = T(q - l + 1) / denom .* (prev[i, j + 1] .- prev[i, j])
                end
            end
            PKL[k + 1][l + 1] = tmp
        end
    end

    # Also differentiate PKL[1][·] in v
    for l in 1:min(d, q)
        prev = PKL[1][l]
        rows = size(prev, 1)
        cols = size(prev, 2) - 1
        tmp = Matrix{Vector{T}}(undef, rows, cols)
        for i in 1:rows
            for j in 1:cols
                denom = Vk[s1 + j - 1 + q + 1] - Vk[s1 + j - 1 + l]
                tmp[i, j] = T(q - l + 1) / denom .* (prev[i, j + 1] .- prev[i, j])
            end
        end
        PKL[1][l + 1] = tmp
    end

    return PKL
end
