# ---------------------------------------------------------------------------- #
#  Chapter 11: Shape Modification Tools                                        #
# ---------------------------------------------------------------------------- #

"""
    move_control_point(crv::BSplineCurve{T}, index::Int, delta::Vector{T}) -> BSplineCurve{T}

Create a new curve with control point `index` displaced by `delta`.

Piegl & Tiller, *The NURBS Book*, 2nd ed., Section 11.2, p. 518.

# Arguments
- `crv::BSplineCurve{T}`: input curve
- `index::Int`: 1-based index of the control point to move
- `delta::Vector{T}`: displacement vector

# Returns
- `BSplineCurve{T}`: modified curve

See also: [`modify_weight`](@ref)
"""
function move_control_point(crv::BSplineCurve{T}, index::Int, delta::Vector{T}) where {T}
    1 <= index <= length(crv.controlpoints) || throw(BoundsError(crv.controlpoints, index))
    new_pts = [copy(pt) for pt in crv.controlpoints]
    new_pts[index] = new_pts[index] .+ delta
    return BSplineCurve(crv.degree, crv.knots, new_pts)
end

"""
    move_control_point(crv::NURBSCurve{T}, index::Int, delta::Vector{T}) -> NURBSCurve{T}

Create a new NURBS curve with control point `index` displaced by `delta`.

See also: [`move_control_point(::BSplineCurve, ::Int, ::Vector)`](@ref)
"""
function move_control_point(crv::NURBSCurve{T}, index::Int, delta::Vector{T}) where {T}
    1 <= index <= length(crv.controlpoints) || throw(BoundsError(crv.controlpoints, index))
    new_pts = [copy(pt) for pt in crv.controlpoints]
    new_pts[index] = new_pts[index] .+ delta
    return NURBSCurve(crv.degree, crv.knots, new_pts, copy(crv.weights))
end

"""
    modify_weight(crv::NURBSCurve{T}, index::Int, new_weight::T) -> NURBSCurve{T}

Create a new NURBS curve with the weight at `index` changed to `new_weight`.

Piegl & Tiller, *The NURBS Book*, 2nd ed., Section 11.3, p. 520.

Increasing a weight pulls the curve toward the corresponding control point;
decreasing it pushes the curve away.

# Arguments
- `crv::NURBSCurve{T}`: input curve
- `index::Int`: 1-based index of the weight to modify
- `new_weight::T`: new weight value (must be positive)

# Returns
- `NURBSCurve{T}`: modified curve

See also: [`move_control_point`](@ref)
"""
function modify_weight(crv::NURBSCurve{T}, index::Int, new_weight::T) where {T}
    1 <= index <= length(crv.weights) || throw(BoundsError(crv.weights, index))
    new_weight > 0 || throw(ArgumentError("Weight must be positive, got $new_weight"))
    new_w = copy(crv.weights)
    new_w[index] = new_weight
    return NURBSCurve(crv.degree, crv.knots, [copy(pt) for pt in crv.controlpoints], new_w)
end

"""
    modify_weight(surf::NURBSSurface{T}, i::Int, j::Int, new_weight::T) -> NURBSSurface{T}

Create a new NURBS surface with the weight at `(i, j)` changed.

Piegl & Tiller, *The NURBS Book*, 2nd ed., Section 11.3.3, p. 533.

# Arguments
- `surf::NURBSSurface{T}`: input surface
- `i, j`: 1-based indices of the weight
- `new_weight::T`: new weight value (must be positive)

# Returns
- `NURBSSurface{T}`: modified surface
"""
function modify_weight(surf::NURBSSurface{T}, i::Int, j::Int, new_weight::T) where {T}
    new_weight > 0 || throw(ArgumentError("Weight must be positive, got $new_weight"))
    new_w = copy(surf.weights)
    new_w[i, j] = new_weight
    new_P = [copy(surf.controlpoints[ii, jj]) for ii in 1:size(surf.controlpoints, 1),
                                                    jj in 1:size(surf.controlpoints, 2)]
    return NURBSSurface(surf.udegree, surf.vdegree, surf.uknots, surf.vknots, new_P, new_w)
end

"""
    move_control_point(surf::BSplineSurface{T}, i::Int, j::Int,
                       delta::Vector{T}) -> BSplineSurface{T}

Create a new B-spline surface with control point `(i, j)` displaced by `delta`.

Piegl & Tiller, *The NURBS Book*, 2nd ed., Section 11.2, p. 518.

# Arguments
- `surf::BSplineSurface{T}`: input surface
- `i, j`: 1-based indices of the control point
- `delta::Vector{T}`: displacement vector

# Returns
- `BSplineSurface{T}`: modified surface
"""
function move_control_point(surf::BSplineSurface{T}, i::Int, j::Int,
                            delta::Vector{T}) where {T}
    new_P = [copy(surf.controlpoints[ii, jj]) for ii in 1:size(surf.controlpoints, 1),
                                                    jj in 1:size(surf.controlpoints, 2)]
    new_P[i, j] = new_P[i, j] .+ delta
    return BSplineSurface(surf.udegree, surf.vdegree, surf.uknots, surf.vknots, new_P)
end
