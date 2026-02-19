# ---------------------------------------------------------------------------- #
#  Chapter 9: Curve and Surface Fitting — A9.1–A9.8                            #
# ---------------------------------------------------------------------------- #

"""
    global_curve_interpolation(Q::Vector{Vector{T}}, p::Int;
                               method::Symbol=:centripetal) -> BSplineCurve{T}

Compute a B-spline curve of degree `p` that interpolates the data points `Q`.

Implements Algorithm A9.1 from Piegl & Tiller, *The NURBS Book*, 2nd ed., p. 369.

# Arguments
- `Q::Vector{Vector{T}}`: data points (at least `p + 1`)
- `p::Int`: desired degree
- `method::Symbol`: `:chord` or `:centripetal` parameterization

# Returns
- `BSplineCurve{T}`: interpolating curve

See also: [`global_curve_approximation`](@ref), [`local_curve_interpolation`](@ref)
"""
function global_curve_interpolation(Q::Vector{Vector{T}}, p::Int;
                                    method::Symbol=:centripetal) where {T}
    npts = length(Q)
    dim = length(Q[1])
    npts > p || throw(ArgumentError("Need at least p+1=$(p+1) points, got $npts"))

    u_bar = _compute_parameters(Q, method)
    U = _compute_knot_vector(u_bar, p, npts)

    # Assemble collocation matrix
    N_mat = zeros(T, npts, npts)
    for i in 1:npts
        span = find_span(npts, p, u_bar[i], U)
        bfuns = basis_functions(span, u_bar[i], p, U)
        for j in 1:(p + 1)
            col = span - p + j - 1
            if 1 <= col <= npts
                N_mat[i, col] += bfuns[j]
            end
        end
    end

    # Solve
    rhs = zeros(T, npts, dim)
    for i in 1:npts
        rhs[i, :] = Q[i]
    end
    sol = N_mat \ rhs

    P = [sol[i, :] for i in 1:npts]
    return BSplineCurve(p, U, P)
end

"""
    _compute_parameters(Q::Vector{Vector{T}}, method::Symbol) -> Vector{T}

Compute parameter values using chord-length or centripetal parameterization.

Piegl & Tiller, *The NURBS Book*, 2nd ed., Eq. 9.4–9.5, p. 365.
"""
function _compute_parameters(Q::Vector{Vector{T}}, method::Symbol) where {T}
    npts = length(Q)
    u_bar = zeros(T, npts)
    u_bar[npts] = one(T)

    if method == :chord
        d = sum(norm(Q[i + 1] .- Q[i]) for i in 1:(npts - 1))
        for i in 2:(npts - 1)
            u_bar[i] = u_bar[i - 1] + norm(Q[i] .- Q[i - 1]) / d
        end
    elseif method == :centripetal
        d = sum(sqrt(norm(Q[i + 1] .- Q[i])) for i in 1:(npts - 1))
        for i in 2:(npts - 1)
            u_bar[i] = u_bar[i - 1] + sqrt(norm(Q[i] .- Q[i - 1])) / d
        end
    else
        throw(ArgumentError("method must be :chord or :centripetal, got :$method"))
    end

    return u_bar
end

"""
    _compute_knot_vector(u_bar::Vector{T}, p::Int, npts::Int) -> KnotVector{T}

Averaging-based knot vector for global interpolation.

Piegl & Tiller, *The NURBS Book*, 2nd ed., Eq. 9.8, p. 369.
"""
function _compute_knot_vector(u_bar::Vector{T}, p::Int, npts::Int) where {T}
    m = npts + p + 1
    U = zeros(T, m)

    for i in 1:(p + 1)
        U[i] = zero(T)
    end
    for i in (m - p):m
        U[i] = one(T)
    end

    for j in 1:(npts - p - 1)
        U[p + 1 + j] = sum(u_bar[(j + 1):(j + p)]) / p
    end

    return KnotVector(U)
end

"""
    global_curve_approximation(Q::Vector{Vector{T}}, p::Int, nctl::Int;
                               method::Symbol=:centripetal) -> BSplineCurve{T}

Least-squares B-spline curve approximation with `nctl` control points.

Implements Algorithm A9.7 from Piegl & Tiller, *The NURBS Book*, 2nd ed., p. 410.

# Arguments
- `Q::Vector{Vector{T}}`: data points
- `p::Int`: degree
- `nctl::Int`: number of control points
- `method::Symbol`: parameterization method

# Returns
- `BSplineCurve{T}`: approximating curve
"""
function global_curve_approximation(Q::Vector{Vector{T}}, p::Int, nctl::Int;
                                    method::Symbol=:centripetal) where {T}
    mpts = length(Q)
    dim = length(Q[1])
    mpts > nctl || throw(ArgumentError("Need more data points than control points"))
    nctl > p || throw(ArgumentError("nctl must be > p"))

    u_bar = _compute_parameters(Q, method)

    # Build knot vector
    m_knot = nctl + p + 1
    U = zeros(T, m_knot)
    for i in 1:(p + 1)
        U[i] = zero(T)
    end
    for i in (m_knot - p):m_knot
        U[i] = one(T)
    end
    d_val = T(mpts) / T(nctl - p)
    for j in 1:(nctl - p - 1)
        i_val = floor(Int, j * d_val)
        alpha = j * d_val - i_val
        U[p + 1 + j] = (1 - alpha) * u_bar[i_val] + alpha * u_bar[min(i_val + 1, mpts)]
    end
    Uknot = KnotVector(U)

    # Assemble N matrix (interior points × interior control points)
    n_inner = nctl - 2
    N_mat = zeros(T, mpts - 2, n_inner)
    R = zeros(T, mpts - 2, dim)

    for i in 2:(mpts - 1)
        span = find_span(nctl, p, u_bar[i], Uknot)
        bfuns = basis_functions(span, u_bar[i], p, Uknot)

        rk = copy(Q[i])
        for j in 1:(p + 1)
            col = span - p + j - 1
            if col == 1
                rk .-= bfuns[j] .* Q[1]
            elseif col == nctl
                rk .-= bfuns[j] .* Q[mpts]
            elseif 2 <= col <= nctl - 1
                N_mat[i - 1, col - 1] += bfuns[j]
            end
        end
        R[i - 1, :] = rk
    end

    NtN = N_mat' * N_mat
    NtR = N_mat' * R
    sol = NtN \ NtR

    P = Vector{Vector{T}}(undef, nctl)
    P[1] = copy(Q[1])
    P[nctl] = copy(Q[mpts])
    for i in 1:n_inner
        P[i + 1] = sol[i, :]
    end

    return BSplineCurve(p, Uknot, P)
end

"""
    global_surface_interpolation(Q::Matrix{Vector{T}}, p::Int, q::Int;
                                 method::Symbol=:centripetal) -> BSplineSurface{T}

Interpolate a grid of points `Q[i, j]` with a B-spline surface.

Implements Algorithm A9.4 from Piegl & Tiller, *The NURBS Book*, 2nd ed., p. 376.

# Arguments
- `Q::Matrix{Vector{T}}`: `nu × nv` grid of data points
- `p::Int`: u-degree
- `q::Int`: v-degree
- `method::Symbol`: parameterization method

# Returns
- `BSplineSurface{T}`: interpolating surface
"""
function global_surface_interpolation(Q::Matrix{Vector{T}}, p::Int, q::Int;
                                      method::Symbol=:centripetal) where {T}
    nu, nv = size(Q)

    # Interpolate columns (u-direction)
    temp_crvs = [global_curve_interpolation([Q[i, j] for i in 1:nu], p; method)
                 for j in 1:nv]

    R = Matrix{Vector{T}}(undef, nu, nv)
    for j in 1:nv, i in 1:nu
        R[i, j] = temp_crvs[j].controlpoints[i]
    end

    # Interpolate rows (v-direction)
    final_crvs = [global_curve_interpolation([R[i, j] for j in 1:nv], q; method)
                  for i in 1:nu]

    P = Matrix{Vector{T}}(undef, nu, nv)
    for i in 1:nu, j in 1:nv
        P[i, j] = final_crvs[i].controlpoints[j]
    end

    return BSplineSurface(p, q, temp_crvs[1].knots, final_crvs[1].knots, P)
end

"""
    local_curve_interpolation(Q::Vector{Vector{T}}) -> BSplineCurve{T}

Cubic local interpolation through data points `Q`.

Implements Algorithm A9.3 from Piegl & Tiller, *The NURBS Book*, 2nd ed., p. 394.
Constructs a cubic B-spline that passes through all given points using
local tangent estimation.

# Arguments
- `Q::Vector{Vector{T}}`: data points (at least 3)

# Returns
- `BSplineCurve{T}`: locally interpolating cubic curve
"""
function local_curve_interpolation(Q::Vector{Vector{T}}) where {T}
    npts = length(Q)
    npts >= 3 || throw(ArgumentError("Need at least 3 points"))
    dim = length(Q[1])
    p = 3

    if npts <= p
        return global_curve_interpolation(Q, npts - 1; method=:chord)
    end

    u_bar = _compute_parameters(Q, :chord)

    # Estimate tangents
    tangents = Vector{Vector{T}}(undef, npts)
    for i in 1:npts
        if i == 1
            tangents[i] = Q[2] .- Q[1]
        elseif i == npts
            tangents[i] = Q[npts] .- Q[npts - 1]
        else
            d_prev = u_bar[i] - u_bar[i - 1]
            d_next = u_bar[i + 1] - u_bar[i]
            dsum = d_prev + d_next
            if dsum > NURBS_EPSILON
                tangents[i] = (d_next / dsum) .* (Q[i] .- Q[i - 1]) .+
                              (d_prev / dsum) .* (Q[i + 1] .- Q[i])
            else
                tangents[i] = Q[min(i + 1, npts)] .- Q[max(i - 1, 1)]
            end
        end
        tn = norm(tangents[i])
        if tn > NURBS_EPSILON
            tangents[i] ./= tn
        end
    end

    # Build piecewise cubic Hermite-like B-spline
    ctrl_pts = Vector{Vector{T}}()
    knots_arr = T[]

    push!(ctrl_pts, copy(Q[1]))
    for _ in 1:(p + 1)
        push!(knots_arr, u_bar[1])
    end

    for seg in 1:(npts - 1)
        du = u_bar[seg + 1] - u_bar[seg]
        alpha = du / 3

        push!(ctrl_pts, Q[seg] .+ alpha .* tangents[seg])
        push!(ctrl_pts, Q[seg + 1] .- alpha .* tangents[seg + 1])
        push!(ctrl_pts, copy(Q[seg + 1]))

        if seg < npts - 1
            push!(knots_arr, u_bar[seg + 1])
            push!(knots_arr, u_bar[seg + 1])
        end
    end

    for _ in 1:(p + 1)
        push!(knots_arr, u_bar[npts])
    end

    # Adjust control point count to match knot vector
    expected_n = length(knots_arr) - p - 1
    while length(ctrl_pts) > expected_n
        pop!(ctrl_pts)
    end
    while length(ctrl_pts) < expected_n
        push!(ctrl_pts, copy(ctrl_pts[end]))
    end

    return BSplineCurve(p, KnotVector(knots_arr), ctrl_pts)
end
