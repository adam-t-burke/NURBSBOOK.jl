# ---------------------------------------------------------------------------- #
#  Chapter 2: B-spline Basis Functions — Algorithms A2.1–A2.5                  #
#                                                                              #
#  All indices are 1-based following Julia convention.                          #
#  n = number of control points, p = degree                                    #
#  Knot vector U has length n + p + 1                                          #
#  Clamped knot vectors: U[1] = … = U[p+1], U[n+1] = … = U[n+p+1]            #
# ---------------------------------------------------------------------------- #

"""
    find_span(n::Int, p::Int, u::Real, U::KnotVector) -> Int

Return the knot span index ``i`` such that ``U[i] \\le u < U[i+1]``.

For the special case ``u = U[n+1]`` (the right end of the parameter domain)
the function returns ``n`` so that the last non-degenerate span is used.

Implements Algorithm A2.1 from Piegl & Tiller, *The NURBS Book*, 2nd ed., p. 68.

# Arguments
- `n::Int`: number of control points
- `p::Int`: degree
- `u::Real`: parameter value
- `U::KnotVector`: knot vector of length `n + p + 1`

# Returns
- `Int`: knot span index (1-based)

See also: [`basis_functions`](@ref), [`basis_function_derivatives`](@ref)
"""
function find_span(n::Int, p::Int, u::Real, U::KnotVector)
    # Special case: u at right end of domain
    if u >= U[n + 1]
        return n
    end
    # Special case: u at left end of domain
    if u <= U[p + 1]
        return p + 1
    end
    # Binary search
    low = p + 1
    high = n + 1
    mid = (low + high) ÷ 2
    while u < U[mid] || u >= U[mid + 1]
        if u < U[mid]
            high = mid
        else
            low = mid
        end
        mid = (low + high) ÷ 2
    end
    return mid
end

"""
    basis_functions(span::Int, u::Real, p::Int, U::KnotVector{T}) -> Vector{T}

Compute the ``p + 1`` nonvanishing B-spline basis functions at parameter `u`.

Returns a vector `N` of length ``p+1`` where `N[j]` is the value of the
basis function whose support includes the knot span `span`, ordered from
left to right: the function associated with control point `span - p + j - 1`.

Implements Algorithm A2.2 from Piegl & Tiller, *The NURBS Book*, 2nd ed., p. 70.

# Arguments
- `span::Int`: knot span index (1-based), from [`find_span`](@ref)
- `u::Real`: parameter value
- `p::Int`: degree
- `U::KnotVector{T}`: knot vector

# Returns
- `Vector{T}` of length `p + 1`

See also: [`find_span`](@ref), [`basis_function_derivatives`](@ref)
"""
function basis_functions(span::Int, u::Real, p::Int, U::KnotVector{T}) where {T}
    N = zeros(T, p + 1)
    left = zeros(T, p + 1)
    right = zeros(T, p + 1)
    N[1] = one(T)
    for j in 1:p
        left[j] = u - U[span + 1 - j]
        right[j] = U[span + j] - u
        saved = zero(T)
        for r in 1:j
            temp = N[r] / (right[r] + left[j - r + 1])
            N[r] = saved + right[r] * temp
            saved = left[j - r + 1] * temp
        end
        N[j + 1] = saved
    end
    return N
end

"""
    basis_function_derivatives(span::Int, u::Real, p::Int, nders::Int,
                               U::KnotVector{T}) -> Matrix{T}

Compute the nonzero basis functions and their derivatives up to order `nders`.

Implements Algorithm A2.3 from Piegl & Tiller, *The NURBS Book*, 2nd ed., p. 72.

# Arguments
- `span::Int`: knot span index (1-based)
- `u::Real`: parameter value
- `p::Int`: degree
- `nders::Int`: maximum derivative order (clamped to `p`)
- `U::KnotVector{T}`: knot vector

# Returns
- `Matrix{T}` of size `(nders+1) × (p+1)` where `ders[k, j]` is the
  `(k-1)`-th derivative of the `j`-th nonzero basis function.
  Row 1 = function values, row 2 = first derivatives, etc.

See also: [`find_span`](@ref), [`basis_functions`](@ref)
"""
function basis_function_derivatives(span::Int, u::Real, p::Int, nders::Int,
                                    U::KnotVector{T}) where {T}
    nd = min(nders, p)
    ndu = zeros(T, p + 1, p + 1)
    a = zeros(T, 2, p + 1)
    ders = zeros(T, nd + 1, p + 1)
    left = zeros(T, p + 1)
    right = zeros(T, p + 1)

    ndu[1, 1] = one(T)
    for j in 1:p
        left[j] = u - U[span + 1 - j]
        right[j] = U[span + j] - u
        saved = zero(T)
        for r in 1:j
            # Lower triangle (knot differences)
            ndu[j + 1, r] = right[r] + left[j - r + 1]
            temp = ndu[r, j] / ndu[j + 1, r]
            # Upper triangle (basis functions)
            ndu[r, j + 1] = saved + right[r] * temp
            saved = left[j - r + 1] * temp
        end
        ndu[j + 1, j + 1] = saved
    end

    # Load basis function values
    for j in 1:(p + 1)
        ders[1, j] = ndu[j, p + 1]
    end

    # Compute derivatives
    for r in 1:(p + 1)
        s1 = 1
        s2 = 2
        a[1, 1] = one(T)
        for k in 1:nd
            d = zero(T)
            rk = r - k
            pk = p - k
            if r > k
                a[s2, 1] = a[s1, 1] / ndu[pk + 2, rk]
                d = a[s2, 1] * ndu[rk, pk + 1]
            end
            j1 = rk >= 0 ? 1 : 1 - rk
            j2 = (r - 2) <= pk ? k - 1 : p - r + 1
            for j in j1:j2
                a[s2, j + 1] = (a[s1, j + 1] - a[s1, j]) / ndu[pk + 2, rk + j]
                d += a[s2, j + 1] * ndu[rk + j, pk + 1]
            end
            if r <= pk + 1
                a[s2, k + 1] = -a[s1, k] / ndu[pk + 2, r]
                d += a[s2, k + 1] * ndu[r, pk + 1]
            end
            ders[k + 1, r] = d
            s1, s2 = s2, s1
        end
    end

    # Multiply by correct factors
    fac = p
    for k in 1:nd
        for j in 1:(p + 1)
            ders[k + 1, j] *= fac
        end
        fac *= (p - k)
    end

    return ders
end

"""
    one_basis_function(p::Int, U::KnotVector{T}, i::Int, u::Real) -> T

Compute a single basis function ``N_{i,p}(u)`` where `i` is the 1-based
index of the basis function (i.e., associated with control point `i`).

Implements Algorithm A2.4 from Piegl & Tiller, *The NURBS Book*, 2nd ed., p. 74.

# Arguments
- `p::Int`: degree
- `U::KnotVector{T}`: knot vector of length `m`
- `i::Int`: 1-based basis function index
- `u::Real`: parameter value

# Returns
- `T`: the value ``N_{i,p}(u)``

See also: [`basis_functions`](@ref), [`one_basis_function_derivatives`](@ref)
"""
function one_basis_function(p::Int, U::KnotVector{T}, i::Int, u::Real) where {T}
    m = length(U)
    # Special cases for endpoints
    if (i == 1 && u == U[1]) || (i == m - p - 1 && u == U[m])
        return one(T)
    end
    # Local property
    if u < U[i] || u >= U[i + p + 1]
        return zero(T)
    end

    N = zeros(T, p + 1)
    # Initialize zeroth-degree functions
    for j in 1:(p + 1)
        if u >= U[i + j - 1] && u < U[i + j]
            N[j] = one(T)
        end
    end

    # Compute triangular table
    for k in 1:p
        if N[1] == zero(T)
            saved = zero(T)
        else
            saved = ((u - U[i]) * N[1]) / (U[i + k] - U[i])
        end
        for j in 1:(p - k + 1)
            Uleft = U[i + j]
            Uright = U[i + j + k]
            if N[j + 1] == zero(T)
                N[j] = saved
                saved = zero(T)
            else
                temp = N[j + 1] / (Uright - Uleft)
                N[j] = saved + (Uright - u) * temp
                saved = (u - Uleft) * temp
            end
        end
    end

    return N[1]
end

"""
    one_basis_function_derivatives(p::Int, U::KnotVector{T}, i::Int,
                                   u::Real, nders::Int) -> Vector{T}

Compute derivatives of a single basis function ``N_{i,p}(u)`` up to order `nders`.

Implements Algorithm A2.5 from Piegl & Tiller, *The NURBS Book*, 2nd ed., p. 76.

# Arguments
- `p::Int`: degree
- `U::KnotVector{T}`: knot vector
- `i::Int`: 1-based basis function index
- `u::Real`: parameter value
- `nders::Int`: maximum derivative order

# Returns
- `Vector{T}` of length `nders + 1` where entry `k` is the `(k-1)`-th derivative.

See also: [`one_basis_function`](@ref), [`basis_function_derivatives`](@ref)
"""
function one_basis_function_derivatives(p::Int, U::KnotVector{T}, i::Int,
                                        u::Real, nders::Int) where {T}
    ders = zeros(T, nders + 1)

    # Local property
    if u < U[i] || u >= U[i + p + 1]
        return ders
    end

    # Build triangular table of basis function values N[j, k]
    # N[j, k] stores the (k-1)-degree function value for the j-th local function
    N = zeros(T, p + 1, p + 1)
    for j in 1:(p + 1)
        if u >= U[i + j - 1] && u < U[i + j]
            N[j, 1] = one(T)
        end
    end

    # Compute full triangular table
    for k in 1:p
        if N[1, k] == zero(T)
            saved = zero(T)
        else
            saved = ((u - U[i]) * N[1, k]) / (U[i + k] - U[i])
        end
        for j in 1:(p - k + 1)
            Uleft = U[i + j]
            Uright = U[i + j + k]
            if N[j + 1, k] == zero(T)
                N[j, k + 1] = saved
                saved = zero(T)
            else
                temp = N[j + 1, k] / (Uright - Uleft)
                N[j, k + 1] = saved + (Uright - u) * temp
                saved = (u - Uleft) * temp
            end
        end
    end

    # Function value
    ders[1] = N[1, p + 1]

    # Compute derivatives using Eq. (2.9)
    for k in 1:min(nders, p)
        ND = zeros(T, k + 1)
        for j in 1:(k + 1)
            ND[j] = N[j, p - k + 1]
        end
        for jj in 1:k
            if ND[1] == zero(T)
                saved = zero(T)
            else
                saved = ND[1] / (U[i + p - k + jj] - U[i])
            end
            for j in 1:(k - jj + 1)
                Uleft = U[i + j]
                Uright = U[i + j + p - k + jj]
                if ND[j + 1] == zero(T)
                    ND[j] = (p - k + jj) * saved
                    saved = zero(T)
                else
                    temp = ND[j + 1] / (Uright - Uleft)
                    ND[j] = (p - k + jj) * (saved - temp)
                    saved = temp
                end
            end
        end
        ders[k + 1] = ND[1]
    end

    return ders
end
