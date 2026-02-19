# ---------------------------------------------------------------------------- #
#                            Core type definitions                             #
# ---------------------------------------------------------------------------- #

using LinearAlgebra

# ---- Abstract hierarchy ---------------------------------------------------- #

"""
    AbstractSplineCurve{T}

Abstract supertype for all spline curve types parameterized by numeric type `T`.
"""
abstract type AbstractSplineCurve{T} end

"""
    AbstractSplineSurface{T}

Abstract supertype for all spline surface types parameterized by numeric type `T`.
"""
abstract type AbstractSplineSurface{T} end

# ---- KnotVector ------------------------------------------------------------ #

"""
    KnotVector{T<:Real} <: AbstractVector{T}

A clamped (nonperiodic / open) knot vector ``U = \\{u_0, \\ldots, u_m\\}`` satisfying:

- Non-decreasing: ``u_i \\le u_{i+1}``
- At least two distinct values

The knot vector is stored as a plain `Vector{T}` and exposed with Julia's
`AbstractVector` interface for convenient indexing.

Piegl & Tiller, *The NURBS Book*, 2nd ed., Eq. 2.13, p. 66.

# Examples
```julia
U = KnotVector([0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0])
```
"""
struct KnotVector{T<:Real} <: AbstractVector{T}
    knots::Vector{T}

    function KnotVector{T}(knots::Vector{T}) where {T<:Real}
        length(knots) >= 2 || throw(ArgumentError("Knot vector must have at least 2 entries"))
        for i in 1:(length(knots) - 1)
            knots[i] <= knots[i + 1] || throw(ArgumentError(
                "Knot vector must be non-decreasing; got U[$i]=$(knots[i]) > U[$(i+1)]=$(knots[i+1])"))
        end
        return new{T}(knots)
    end
end

KnotVector(knots::Vector{T}) where {T<:Real} = KnotVector{T}(knots)
KnotVector(knots::AbstractVector{T}) where {T<:Real} = KnotVector{T}(collect(knots))

Base.size(U::KnotVector) = size(U.knots)
Base.getindex(U::KnotVector, i::Int) = U.knots[i]
Base.IndexStyle(::Type{<:KnotVector}) = IndexLinear()
Base.similar(U::KnotVector{T}, ::Type{S}, dims::Dims) where {T,S} = similar(U.knots, S, dims)

# ---- BSplineCurve --------------------------------------------------------- #

"""
    BSplineCurve{T<:Real} <: AbstractSplineCurve{T}

A non-rational B-spline curve of degree `p` defined by

```math
C(u) = \\sum_{i=0}^{n} N_{i,p}(u)\\, P_i
```

where ``\\{P_i\\}`` are the control points and ``N_{i,p}`` are B-spline basis
functions on the knot vector ``U``.

Piegl & Tiller, *The NURBS Book*, 2nd ed., Eq. 3.1, p. 82.

# Fields
- `degree::Int` — polynomial degree `p`
- `knots::KnotVector{T}` — knot vector of length `m + 1 = n + p + 2`
- `controlpoints::Vector{Vector{T}}` — `n + 1` control points in ``\\mathbb{R}^d``
"""
struct BSplineCurve{T<:Real} <: AbstractSplineCurve{T}
    degree::Int
    knots::KnotVector{T}
    controlpoints::Vector{Vector{T}}

    function BSplineCurve{T}(degree::Int, knots::KnotVector{T},
                             controlpoints::Vector{Vector{T}}) where {T<:Real}
        p = degree
        n = length(controlpoints) - 1
        m = length(knots) - 1
        p >= 0 || throw(ArgumentError("Degree must be non-negative, got $p"))
        m == n + p + 1 || throw(ArgumentError(
            "Knot vector length ($(m+1)) must equal n + p + 2 = $(n + p + 2)"))
        return new{T}(degree, knots, controlpoints)
    end
end

function BSplineCurve(degree::Int, knots::KnotVector{T},
                      controlpoints::Vector{Vector{T}}) where {T<:Real}
    return BSplineCurve{T}(degree, knots, controlpoints)
end

# ---- NURBSCurve ----------------------------------------------------------- #

"""
    NURBSCurve{T<:Real} <: AbstractSplineCurve{T}

A Non-Uniform Rational B-Spline (NURBS) curve of degree `p` defined by

```math
C(u) = \\frac{\\sum_{i=0}^{n} N_{i,p}(u)\\, w_i\\, P_i}
             {\\sum_{i=0}^{n} N_{i,p}(u)\\, w_i}
```

Piegl & Tiller, *The NURBS Book*, 2nd ed., Eq. 4.2, p. 118.

# Fields
- `degree::Int` — polynomial degree `p`
- `knots::KnotVector{T}` — knot vector of length `n + p + 2`
- `controlpoints::Vector{Vector{T}}` — `n + 1` control points in ``\\mathbb{R}^d``
- `weights::Vector{T}` — `n + 1` positive weights
"""
struct NURBSCurve{T<:Real} <: AbstractSplineCurve{T}
    degree::Int
    knots::KnotVector{T}
    controlpoints::Vector{Vector{T}}
    weights::Vector{T}

    function NURBSCurve{T}(degree::Int, knots::KnotVector{T},
                           controlpoints::Vector{Vector{T}},
                           weights::Vector{T}) where {T<:Real}
        p = degree
        n = length(controlpoints) - 1
        m = length(knots) - 1
        p >= 0 || throw(ArgumentError("Degree must be non-negative, got $p"))
        m == n + p + 1 || throw(ArgumentError(
            "Knot vector length ($(m+1)) must equal n + p + 2 = $(n + p + 2)"))
        length(weights) == n + 1 || throw(ArgumentError(
            "Number of weights ($(length(weights))) must equal number of control points ($(n+1))"))
        all(w -> w > 0, weights) || throw(ArgumentError("All weights must be positive"))
        return new{T}(degree, knots, controlpoints, weights)
    end
end

function NURBSCurve(degree::Int, knots::KnotVector{T},
                    controlpoints::Vector{Vector{T}},
                    weights::Vector{T}) where {T<:Real}
    return NURBSCurve{T}(degree, knots, controlpoints, weights)
end

# ---- BSplineSurface ------------------------------------------------------- #

"""
    BSplineSurface{T<:Real} <: AbstractSplineSurface{T}

A non-rational B-spline surface defined by

```math
S(u,v) = \\sum_{i=0}^{n} \\sum_{j=0}^{m} N_{i,p}(u)\\, N_{j,q}(v)\\, P_{i,j}
```

Piegl & Tiller, *The NURBS Book*, 2nd ed., Eq. 3.11, p. 100.

# Fields
- `udegree::Int` — degree in u-direction (`p`)
- `vdegree::Int` — degree in v-direction (`q`)
- `uknots::KnotVector{T}` — knot vector in u (length `n + p + 2`)
- `vknots::KnotVector{T}` — knot vector in v (length `m + q + 2`)
- `controlpoints::Matrix{Vector{T}}` — `(n+1) × (m+1)` control net
"""
struct BSplineSurface{T<:Real} <: AbstractSplineSurface{T}
    udegree::Int
    vdegree::Int
    uknots::KnotVector{T}
    vknots::KnotVector{T}
    controlpoints::Matrix{Vector{T}}

    function BSplineSurface{T}(udeg::Int, vdeg::Int,
                               uknots::KnotVector{T}, vknots::KnotVector{T},
                               controlpoints::Matrix{Vector{T}}) where {T<:Real}
        n = size(controlpoints, 1) - 1
        m = size(controlpoints, 2) - 1
        r = length(uknots) - 1
        s = length(vknots) - 1
        udeg >= 0 || throw(ArgumentError("u-degree must be non-negative"))
        vdeg >= 0 || throw(ArgumentError("v-degree must be non-negative"))
        r == n + udeg + 1 || throw(ArgumentError(
            "u-knot vector length ($(r+1)) must equal n + p + 2 = $(n + udeg + 2)"))
        s == m + vdeg + 1 || throw(ArgumentError(
            "v-knot vector length ($(s+1)) must equal m + q + 2 = $(m + vdeg + 2)"))
        return new{T}(udeg, vdeg, uknots, vknots, controlpoints)
    end
end

function BSplineSurface(udeg::Int, vdeg::Int,
                        uknots::KnotVector{T}, vknots::KnotVector{T},
                        controlpoints::Matrix{Vector{T}}) where {T<:Real}
    return BSplineSurface{T}(udeg, vdeg, uknots, vknots, controlpoints)
end

# ---- NURBSSurface --------------------------------------------------------- #

"""
    NURBSSurface{T<:Real} <: AbstractSplineSurface{T}

A NURBS surface defined by

```math
S(u,v) = \\frac{\\sum_{i=0}^{n}\\sum_{j=0}^{m} N_{i,p}(u)\\,N_{j,q}(v)\\,w_{i,j}\\,P_{i,j}}
              {\\sum_{i=0}^{n}\\sum_{j=0}^{m} N_{i,p}(u)\\,N_{j,q}(v)\\,w_{i,j}}
```

Piegl & Tiller, *The NURBS Book*, 2nd ed., Eq. 4.15, p. 134.

# Fields
- `udegree::Int`, `vdegree::Int` — degrees in u and v
- `uknots::KnotVector{T}`, `vknots::KnotVector{T}` — knot vectors
- `controlpoints::Matrix{Vector{T}}` — `(n+1) × (m+1)` control net
- `weights::Matrix{T}` — `(n+1) × (m+1)` positive weights
"""
struct NURBSSurface{T<:Real} <: AbstractSplineSurface{T}
    udegree::Int
    vdegree::Int
    uknots::KnotVector{T}
    vknots::KnotVector{T}
    controlpoints::Matrix{Vector{T}}
    weights::Matrix{T}

    function NURBSSurface{T}(udeg::Int, vdeg::Int,
                             uknots::KnotVector{T}, vknots::KnotVector{T},
                             controlpoints::Matrix{Vector{T}},
                             weights::Matrix{T}) where {T<:Real}
        n = size(controlpoints, 1) - 1
        m = size(controlpoints, 2) - 1
        r = length(uknots) - 1
        s = length(vknots) - 1
        udeg >= 0 || throw(ArgumentError("u-degree must be non-negative"))
        vdeg >= 0 || throw(ArgumentError("v-degree must be non-negative"))
        r == n + udeg + 1 || throw(ArgumentError(
            "u-knot vector length ($(r+1)) must equal n + p + 2 = $(n + udeg + 2)"))
        s == m + vdeg + 1 || throw(ArgumentError(
            "v-knot vector length ($(s+1)) must equal m + q + 2 = $(m + vdeg + 2)"))
        size(weights) == size(controlpoints) || throw(ArgumentError(
            "Weight matrix size must match control point matrix size"))
        all(w -> w > 0, weights) || throw(ArgumentError("All weights must be positive"))
        return new{T}(udeg, vdeg, uknots, vknots, controlpoints, weights)
    end
end

function NURBSSurface(udeg::Int, vdeg::Int,
                      uknots::KnotVector{T}, vknots::KnotVector{T},
                      controlpoints::Matrix{Vector{T}},
                      weights::Matrix{T}) where {T<:Real}
    return NURBSSurface{T}(udeg, vdeg, uknots, vknots, controlpoints, weights)
end
