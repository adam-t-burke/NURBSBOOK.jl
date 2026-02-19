# ---------------------------------------------------------------------------- #
#  Chapter 4: Rational B-spline (NURBS) Curves & Surfaces — A4.1–A4.5         #
# ---------------------------------------------------------------------------- #

# ---- Helper: project from homogeneous to Cartesian coordinates ------------ #

function _to_homogeneous(P::Vector{Vector{T}}, w::Vector{T}) where {T}
    return [vcat(w[i] .* P[i], w[i]) for i in eachindex(P)]
end

function _from_homogeneous(Pw::Vector{Vector{T}}) where {T}
    dim = length(Pw[1]) - 1
    n = length(Pw)
    pts = Vector{Vector{T}}(undef, n)
    wts = Vector{T}(undef, n)
    for i in 1:n
        wts[i] = Pw[i][end]
        pts[i] = Pw[i][1:dim] ./ wts[i]
    end
    return pts, wts
end

# ---- NURBS Curve ---------------------------------------------------------- #

"""
    curve_point(crv::NURBSCurve{T}, u::Real) -> Vector{T}

Compute a point on a NURBS curve at parameter ``u``.

A NURBS curve of degree ``p`` is defined by (Eq. 4.2):

```math
C(u) = \\frac{\\sum_{i=0}^{n} N_{i,p}(u)\\, w_i\\, P_i}
             {\\sum_{i=0}^{n} N_{i,p}(u)\\, w_i}
     = \\sum_{i=0}^{n} R_{i,p}(u)\\, P_i
```

where the rational basis functions are (Eq. 4.3):

```math
R_{i,p}(u) = \\frac{N_{i,p}(u)\\, w_i}{\\sum_{j=0}^{n} N_{j,p}(u)\\, w_j}
```

This function works in homogeneous (weighted) coordinates
``P_i^w = (w_i P_i,\\, w_i)`` and projects back to Cartesian space.

Implements **Algorithm A4.1** from Piegl & Tiller, *The NURBS Book*, 2nd ed., p. 124.

# Arguments
- `crv::NURBSCurve{T}`: NURBS curve
- `u::Real`: parameter value

# Returns
- `Vector{T}`: the point ``C(u)``

See also: [`curve_derivatives(::NURBSCurve, ::Real, ::Int)`](@ref)
"""
function curve_point(crv::NURBSCurve{T}, u::Real) where {T}
    p = crv.degree
    U = crv.knots
    P = crv.controlpoints
    w = crv.weights
    n = length(P)
    dim = length(P[1])

    span = find_span(n, p, u, U)
    N = basis_functions(span, u, p, U)

    Sw = zeros(T, dim + 1)
    for j in 1:(p + 1)
        idx = span - p + j - 1
        wi = w[idx]
        for d in 1:dim
            Sw[d] += N[j] * wi * P[idx][d]
        end
        Sw[dim + 1] += N[j] * wi
    end

    return Sw[1:dim] ./ Sw[dim + 1]
end

"""
    curve_derivatives(crv::NURBSCurve{T}, u::Real, d::Int) -> Vector{Vector{T}}

Compute derivatives of a NURBS curve up to order ``d`` at parameter ``u``.

Let ``C^w(u) = \\sum N_{i,p}(u)\\, P_i^w`` denote the curve in homogeneous
space, with components ``A(u)`` (spatial, dimension ``d``) and ``w(u)``
(weight). The ``k``-th derivative of the rational curve is obtained via
(Eq. 4.8):

```math
C^{(k)}(u) = \\frac{A^{(k)}(u) - \\sum_{i=1}^{k} \\binom{k}{i}\\,
  w^{(i)}(u)\\, C^{(k-i)}(u)}{w(u)}
```

This recurrence "projects" the homogeneous derivatives to Cartesian space,
starting from ``C^{(0)}(u) = A(u)/w(u)``.

Implements **Algorithm A4.2** from Piegl & Tiller, *The NURBS Book*, 2nd ed., p. 127.

# Arguments
- `crv::NURBSCurve{T}`: NURBS curve
- `u::Real`: parameter value
- `d::Int`: maximum derivative order

# Returns
- `Vector{Vector{T}}` of length ``d + 1``

See also: [`curve_point(::NURBSCurve, ::Real)`](@ref)
"""
function curve_derivatives(crv::NURBSCurve{T}, u::Real, d::Int) where {T}
    p = crv.degree
    U = crv.knots
    P = crv.controlpoints
    w = crv.weights
    n = length(P)
    dim = length(P[1])

    # Work in homogeneous coordinates
    Pw = _to_homogeneous(P, w)
    bcrv = BSplineCurve(p, U, Pw)
    CKw = curve_derivatives(bcrv, u, d)

    # Rational projection (Eq. 4.8)
    CK = [zeros(T, dim) for _ in 1:(d + 1)]
    for k in 1:(d + 1)
        v = CKw[k][1:dim]
        for i in 2:k
            v .-= binomial_coefficient(k - 1, i - 1) .* CKw[i][dim + 1] .* CK[k - i + 1]
        end
        CK[k] = v ./ CKw[1][dim + 1]
    end

    return CK
end

# ---- NURBS Surface -------------------------------------------------------- #

"""
    surface_point(surf::NURBSSurface{T}, u::Real, v::Real) -> Vector{T}

Compute a point on a NURBS surface at parameters ``(u, v)``.

A NURBS surface of degree ``p`` in ``u`` and ``q`` in ``v`` is defined by
(Eq. 4.15):

```math
S(u,v) = \\frac{\\sum_{i=0}^{n} \\sum_{j=0}^{m}
  N_{i,p}(u)\\, N_{j,q}(v)\\, w_{i,j}\\, P_{i,j}}
  {\\sum_{i=0}^{n} \\sum_{j=0}^{m}
  N_{i,p}(u)\\, N_{j,q}(v)\\, w_{i,j}}
```

Evaluation is performed in homogeneous coordinates and projected to
Cartesian space.

Implements **Algorithm A4.3** from Piegl & Tiller, *The NURBS Book*, 2nd ed., p. 134.

# Arguments
- `surf::NURBSSurface{T}`: NURBS surface
- `u::Real`, `v::Real`: parameter values

# Returns
- `Vector{T}`: the point ``S(u,v)``
"""
function surface_point(surf::NURBSSurface{T}, u::Real, v::Real) where {T}
    p = surf.udegree
    q = surf.vdegree
    Uk = surf.uknots
    Vk = surf.vknots
    P = surf.controlpoints
    w = surf.weights
    nu = size(P, 1)
    nv = size(P, 2)
    dim = length(P[1, 1])

    uspan = find_span(nu, p, u, Uk)
    Nu = basis_functions(uspan, u, p, Uk)
    vspan = find_span(nv, q, v, Vk)
    Nv = basis_functions(vspan, v, q, Vk)

    Sw = zeros(T, dim + 1)
    for l in 1:(q + 1)
        temp = zeros(T, dim + 1)
        vi = vspan - q + l - 1
        for k in 1:(p + 1)
            ui = uspan - p + k - 1
            wi = w[ui, vi]
            for dd in 1:dim
                temp[dd] += Nu[k] * wi * P[ui, vi][dd]
            end
            temp[dim + 1] += Nu[k] * wi
        end
        Sw .+= Nv[l] .* temp
    end
    return Sw[1:dim] ./ Sw[dim + 1]
end

"""
    surface_derivatives(surf::NURBSSurface{T}, u::Real, v::Real,
                        d::Int) -> Matrix{Vector{T}}

Compute partial derivatives of a NURBS surface up to total order ``d``.

The ``(k,l)``-th partial derivative of a rational surface is obtained by
projecting from homogeneous coordinates via (Eq. 4.20):

```math
S^{(k,l)} = \\frac{1}{w}\\biggl[
  A^{(k,l)}
  - \\sum_{j=1}^{l} \\binom{l}{j}\\, w^{(0,j)}\\, S^{(k,l-j)}
  - \\sum_{i=1}^{k} \\binom{k}{i}\\!\\left(
      w^{(i,0)}\\, S^{(k-i,l)}
    + \\sum_{j=1}^{l} \\binom{l}{j}\\, w^{(i,j)}\\, S^{(k-i,l-j)}
  \\right)
\\biggr]
```

where ``A^{(k,l)}`` and ``w^{(k,l)}`` are the spatial and weight components
of the ``(k,l)``-th derivative of the surface in homogeneous space.

Implements **Algorithm A4.4** from Piegl & Tiller, *The NURBS Book*, 2nd ed., p. 137.

# Arguments
- `surf::NURBSSurface{T}`: NURBS surface
- `u::Real`, `v::Real`: parameter values
- `d::Int`: maximum total derivative order

# Returns
- `Matrix{Vector{T}}` of size ``(d+1) \\times (d+1)``
"""
function surface_derivatives(surf::NURBSSurface{T}, u::Real, v::Real,
                             d::Int) where {T}
    p = surf.udegree
    q = surf.vdegree
    P = surf.controlpoints
    w = surf.weights
    nu = size(P, 1)
    nv = size(P, 2)
    dim = length(P[1, 1])

    Pw = [vcat(w[i, j] .* P[i, j], w[i, j]) for i in 1:nu, j in 1:nv]
    bsurf = BSplineSurface(p, q, surf.uknots, surf.vknots, Pw)
    SKLw = surface_derivatives(bsurf, u, v, d)

    SKL = [zeros(T, dim) for _ in 1:(d + 1), _ in 1:(d + 1)]
    for k in 1:(d + 1)
        for l in 1:(d - k + 2)
            vd = SKLw[k, l][1:dim]
            for j in 2:l
                vd .-= binomial_coefficient(l - 1, j - 1) .* SKLw[1, j][dim + 1] .* SKL[k, l - j + 1]
            end
            for i in 2:k
                vd .-= binomial_coefficient(k - 1, i - 1) .* SKLw[i, 1][dim + 1] .* SKL[k - i + 1, l]
                v2 = zeros(T, dim)
                for j in 2:l
                    v2 .+= binomial_coefficient(l - 1, j - 1) .* SKLw[i, j][dim + 1] .* SKL[k - i + 1, l - j + 1]
                end
                vd .-= binomial_coefficient(k - 1, i - 1) .* v2
            end
            SKL[k, l] = vd ./ SKLw[1, 1][dim + 1]
        end
    end

    return SKL
end
