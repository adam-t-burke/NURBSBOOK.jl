# ---------------------------------------------------------------------------- #
#  Chapter 11: Shape Modification Tools                                        #
# ---------------------------------------------------------------------------- #

"""
    move_control_point(crv::BSplineCurve{T}, index::Int, delta::Vector{T}) -> BSplineCurve{T}

Create a new curve with control point ``P_i`` displaced by ``\\Delta P_i``.

The modified curve is (Section 11.2):

```math
\\tilde{C}(u) = C(u) + N_{i,p}(u)\\, \\Delta P_i
```

Since basis function ``N_{i,p}`` has local support on
``[u_i,\\, u_{i+p+1}]``, the curve shape changes only in this interval.

Piegl & Tiller, *The NURBS Book*, 2nd ed., Section 11.2, p. 518.

# Arguments
- `crv::BSplineCurve{T}`: input curve
- `index::Int`: 1-based index of the control point to move
- `delta::Vector{T}`: displacement vector ``\\Delta P_i``

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

Create a new NURBS curve with control point ``P_i`` displaced by
``\\Delta P_i``. The weights remain unchanged.

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

Create a new NURBS curve with the weight ``w_i`` changed to ``w_i^*``.

The effect of changing a single weight is (Eq. 11.1):

```math
C(u;\\, w_i^*) = C(u;\\, w_i) +
  \\frac{R_{i,p}(u)\\,(w_i^* - w_i)}
       {w_i^*\\, R_{i,p}(u) + w_i\\,(1 - R_{i,p}(u))}
  \\bigl(P_i - C(u;\\, w_i)\\bigr)
```

Increasing ``w_i`` pulls the curve toward ``P_i``; decreasing pushes it away.
The deformation is local, confined to the support ``[u_i, u_{i+p+1}]``.

Piegl & Tiller, *The NURBS Book*, 2nd ed., Section 11.3, p. 520.

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

Create a new NURBS surface with the weight ``w_{i,j}`` changed to
``w_{i,j}^*``.

Analogous to the curve case, modifying a surface weight produces a local
deformation toward (or away from) the control point ``P_{i,j}``:

```math
S(u,v;\\, w_{i,j}^*) = S(u,v;\\, w_{i,j}) +
  \\frac{R_{i,j}(u,v)\\,(w_{i,j}^* - w_{i,j})}
       {w_{i,j}^*\\, R_{i,j}(u,v) + w_{i,j}\\,(1 - R_{i,j}(u,v))}
  \\bigl(P_{i,j} - S(u,v;\\, w_{i,j})\\bigr)
```

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

Create a new B-spline surface with control point ``P_{i,j}`` displaced by
``\\Delta P_{i,j}``.

The modified surface is:

```math
\\tilde{S}(u,v) = S(u,v) + N_{i,p}(u)\\, N_{j,q}(v)\\, \\Delta P_{i,j}
```

The deformation is local, confined to the tensor-product support region
``[u_i, u_{i+p+1}] \\times [v_j, v_{j+q+1}]``.

Piegl & Tiller, *The NURBS Book*, 2nd ed., Section 11.2, p. 518.

# Arguments
- `surf::BSplineSurface{T}`: input surface
- `i, j`: 1-based indices of the control point
- `delta::Vector{T}`: displacement vector ``\\Delta P_{i,j}``

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
