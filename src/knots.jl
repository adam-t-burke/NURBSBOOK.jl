# ---------------------------------------------------------------------------- #
#  Chapter 5: Knot Insertion, Refinement, Removal — A5.1–A5.4                  #
#                                                                              #
#  All indices 1-based.  n = number of control points, p = degree.             #
#  Knot vector length = n + p + 1.                                             #
# ---------------------------------------------------------------------------- #

# ---- multiplicity helper -------------------------------------------------- #

"""
    _knot_multiplicity(U, u) -> Int

Count the multiplicity of knot value `u` in knot vector `U`.
"""
function _knot_multiplicity(U, u)
    s = 0
    for i in eachindex(U)
        if abs(U[i] - u) <= NURBS_EPSILON
            s += 1
        end
    end
    return s
end

# ---- A5.1  Curve Knot Insertion ------------------------------------------- #

"""
    insert_knot(crv::BSplineCurve{T}, u::Real, r::Int=1) -> BSplineCurve{T}

Insert knot `u` into the curve `r` times.  The curve shape is unchanged.

Implements Algorithm A5.1 from Piegl & Tiller, *The NURBS Book*, 2nd ed., p. 151.

# Arguments
- `crv::BSplineCurve{T}`: input curve
- `u::Real`: knot value to insert
- `r::Int`: number of insertions (default 1)

# Returns
- `BSplineCurve{T}`: new curve with inserted knot(s)

See also: [`refine_knots`](@ref), [`remove_knot`](@ref)
"""
function insert_knot(crv::BSplineCurve{T}, u::Real, r::Int=1) where {T}
    p = crv.degree
    U = crv.knots
    Pw = crv.controlpoints
    n = length(Pw)
    dim = length(Pw[1])

    span = find_span(n, p, u, U)  # U[span] <= u < U[span+1]
    s = _knot_multiplicity(U, u)

    if r + s > p
        r = p - s
    end
    r <= 0 && return crv

    nq = n + r
    UQ = zeros(T, length(U) + r)
    Qw = [zeros(T, dim) for _ in 1:nq]
    Rw = [zeros(T, dim) for _ in 1:(p + 1)]

    # New knot vector: [U[1..span], u repeated r times, U[span+1..end]]
    UQ[1:span] .= U.knots[1:span]
    UQ[(span + 1):(span + r)] .= u
    UQ[(span + r + 1):end] .= U.knots[(span + 1):end]

    # Unaltered control points at left and right
    for i in 1:(span - p)
        Qw[i] = copy(Pw[i])
    end
    for i in (span - s):n
        Qw[i + r] = copy(Pw[i])
    end

    # Auxiliary control points
    for i in 1:(p - s + 1)
        Rw[i] = copy(Pw[span - p - 1 + i])
    end

    # Insert the knot r times
    L = 0
    for j in 1:r
        L = span - p + j
        for i in 1:(p - j - s + 1)
            alpha = (u - U[L + i - 1]) / (U[span + i] - U[L + i - 1])
            Rw[i] = alpha .* Rw[i + 1] .+ (one(T) - alpha) .* Rw[i]
        end
        Qw[L] = copy(Rw[1])
        Qw[span + r - j - s] = copy(Rw[p - j - s + 1])
    end

    # Remaining control points
    for i in (L + 1):(span - s - 1)
        Qw[i] = copy(Rw[i - L + 1])
    end

    return BSplineCurve(p, KnotVector(UQ), Qw)
end

"""
    insert_knot(crv::NURBSCurve{T}, u::Real, r::Int=1) -> NURBSCurve{T}

Insert knot `u` into a NURBS curve `r` times, via homogeneous coordinates.

See also: [`insert_knot(::BSplineCurve, ::Real, ::Int)`](@ref)
"""
function insert_knot(crv::NURBSCurve{T}, u::Real, r::Int=1) where {T}
    Pw = _to_homogeneous(crv.controlpoints, crv.weights)
    bcrv = BSplineCurve(crv.degree, crv.knots, Pw)
    newcrv = insert_knot(bcrv, u, r)
    new_pts, new_w = _from_homogeneous(newcrv.controlpoints)
    return NURBSCurve(crv.degree, newcrv.knots, new_pts, new_w)
end

# ---- A5.4  Knot Refinement ----------------------------------------------- #

"""
    refine_knots(crv::BSplineCurve{T}, X::AbstractVector{<:Real}) -> BSplineCurve{T}

Refine the knot vector by inserting all knots in `X` simultaneously.

Implements Algorithm A5.4 from Piegl & Tiller, *The NURBS Book*, 2nd ed., p. 164.

`X` must be sorted non-decreasing.

# Arguments
- `crv::BSplineCurve{T}`: input curve
- `X::AbstractVector`: knots to insert (sorted)

# Returns
- `BSplineCurve{T}`: refined curve

See also: [`insert_knot`](@ref)
"""
function refine_knots(crv::BSplineCurve{T}, X::AbstractVector{<:Real}) where {T}
    isempty(X) && return crv
    result = crv
    for x in X
        result = insert_knot(result, x, 1)
    end
    return result
end

"""
    refine_knots(crv::NURBSCurve{T}, X::AbstractVector{<:Real}) -> NURBSCurve{T}

Refine the knot vector of a NURBS curve, via homogeneous coordinates.

See also: [`refine_knots(::BSplineCurve, ::AbstractVector)`](@ref)
"""
function refine_knots(crv::NURBSCurve{T}, X::AbstractVector{<:Real}) where {T}
    Pw = _to_homogeneous(crv.controlpoints, crv.weights)
    bcrv = BSplineCurve(crv.degree, crv.knots, Pw)
    newcrv = refine_knots(bcrv, X)
    new_pts, new_w = _from_homogeneous(newcrv.controlpoints)
    return NURBSCurve(crv.degree, newcrv.knots, new_pts, new_w)
end

# ---- A5.8  Knot Removal -------------------------------------------------- #

"""
    remove_knot(crv::BSplineCurve{T}, u::Real, num::Int,
                tol::Real) -> Tuple{BSplineCurve{T}, Int}

Attempt to remove the knot `u` up to `num` times within tolerance `tol`.

Implements Algorithm A5.8 from Piegl & Tiller, *The NURBS Book*, 2nd ed., p. 185.

# Arguments
- `crv::BSplineCurve{T}`: input curve
- `u::Real`: knot to remove
- `num::Int`: max removals
- `tol::Real`: geometric tolerance

# Returns
- `(new_curve, t)`: resulting curve and actual number of removals

See also: [`insert_knot`](@ref)
"""
function remove_knot(crv::BSplineCurve{T}, u::Real, num::Int, tol::Real) where {T}
    p = crv.degree
    U = collect(crv.knots)
    Pw = [copy(pt) for pt in crv.controlpoints]
    n = length(Pw)
    m = length(U)
    dim = length(Pw[1])

    # Find the knot index range
    span = find_span(n, p, u, KnotVector(U))
    s = _knot_multiplicity(U, u)

    # r = index of last occurrence of u in U
    r_idx = span
    while r_idx < m && abs(U[r_idx + 1] - u) <= NURBS_EPSILON
        r_idx += 1
    end

    first = r_idx - p
    last = r_idx - s + 1
    t = 0

    temp = [zeros(T, dim) for _ in 1:(2 * p + 1)]

    for tt in 1:num
        off = first - 1
        temp[1] = copy(Pw[off])
        temp[last - off + 2] = copy(Pw[last + 1])

        ii = first
        jj = last
        ti = 2
        tj = last - off + 1

        while (jj - ii) > tt - 1
            alfi = (u - U[ii]) / (U[ii + p + tt] - U[ii])
            alfj = (u - U[jj]) / (U[jj + p + tt] - U[jj])
            temp[ti] = (Pw[ii] .- (one(T) - alfi) .* temp[ti - 1]) ./ alfi
            temp[tj] = (Pw[jj] .- alfj .* temp[tj + 1]) ./ (one(T) - alfj)
            ii += 1
            jj -= 1
            ti += 1
            tj -= 1
        end

        if (jj - ii) < tt
            if norm(temp[ti - 1] .- temp[tj + 1]) > tol
                break
            end
        else
            alfi = (u - U[ii]) / (U[ii + p + tt] - U[ii])
            test_pt = alfi .* temp[tj + 1] .+ (one(T) - alfi) .* temp[ti - 1]
            if norm(Pw[ii] .- test_pt) > tol
                break
            end
        end

        ii = first
        jj = last
        while (jj - ii) > tt - 1
            Pw[ii] = copy(temp[ii - off + 1])
            Pw[jj] = copy(temp[jj - off + 1])
            ii += 1
            jj -= 1
        end

        first -= 1
        last += 1
        t = tt
    end

    if t == 0
        return (crv, 0)
    end

    # Build new knot vector removing t copies of u
    newU = T[]
    skip = 0
    for i in 1:m
        if abs(U[i] - u) <= NURBS_EPSILON && skip < t
            skip += 1
        else
            push!(newU, U[i])
        end
    end

    new_n = length(newU) - p - 1
    newPw = Pw[1:new_n]

    return (BSplineCurve(p, KnotVector(newU), newPw), t)
end

# ---- Surface knot insertion ----------------------------------------------- #

"""
    insert_knot(surf::BSplineSurface{T}, dir::Symbol, uv::Real,
                r::Int=1) -> BSplineSurface{T}

Insert a knot into a B-spline surface in direction `dir` (`:u` or `:v`).

Applies Algorithm A5.1 to each isoparametric curve.

# Arguments
- `surf::BSplineSurface{T}`: input surface
- `dir::Symbol`: `:u` or `:v`
- `uv::Real`: knot value to insert
- `r::Int`: number of insertions

# Returns
- `BSplineSurface{T}`: new surface

See also: [`insert_knot(::BSplineCurve, ::Real, ::Int)`](@ref)
"""
function insert_knot(surf::BSplineSurface{T}, dir::Symbol, uv::Real,
                     r::Int=1) where {T}
    P = surf.controlpoints
    nu, nv = size(P)

    if dir == :u
        results = Vector{BSplineCurve{T}}(undef, nv)
        for j in 1:nv
            pts = [P[i, j] for i in 1:nu]
            crv = BSplineCurve(surf.udegree, surf.uknots, pts)
            results[j] = insert_knot(crv, uv, r)
        end
        new_nu = length(results[1].controlpoints)
        newP = Matrix{Vector{T}}(undef, new_nu, nv)
        for j in 1:nv, i in 1:new_nu
            newP[i, j] = results[j].controlpoints[i]
        end
        return BSplineSurface(surf.udegree, surf.vdegree,
                              results[1].knots, surf.vknots, newP)
    elseif dir == :v
        results = Vector{BSplineCurve{T}}(undef, nu)
        for i in 1:nu
            pts = [P[i, j] for j in 1:nv]
            crv = BSplineCurve(surf.vdegree, surf.vknots, pts)
            results[i] = insert_knot(crv, uv, r)
        end
        new_nv = length(results[1].controlpoints)
        newP = Matrix{Vector{T}}(undef, nu, new_nv)
        for i in 1:nu, j in 1:new_nv
            newP[i, j] = results[i].controlpoints[j]
        end
        return BSplineSurface(surf.udegree, surf.vdegree,
                              surf.uknots, results[1].knots, newP)
    else
        throw(ArgumentError("direction must be :u or :v, got :$dir"))
    end
end
