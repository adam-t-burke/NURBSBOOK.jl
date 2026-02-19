# ---------------------------------------------------------------------------- #
#                            Utility routines                                  #
# ---------------------------------------------------------------------------- #

const NURBS_EPSILON = 1.0e-10

"""
    binomial_coefficient(n::Int, k::Int) -> Int

Compute the binomial coefficient using the multiplicative formula:

```math
\\binom{n}{k} = \\prod_{i=0}^{k-1} \\frac{n - i}{i + 1}
```

Returns 0 when ``k < 0`` or ``k > n``.

Used throughout derivative computations (Eq. 4.8, 4.20, 5.36, etc.) in
Piegl & Tiller, *The NURBS Book*, 2nd ed.
"""
function binomial_coefficient(n::Int, k::Int)
    (k < 0 || k > n) && return 0
    k = min(k, n - k)
    result = 1
    for i in 0:(k - 1)
        result = result * (n - i) ÷ (i + 1)
    end
    return result
end

"""
    allclose(a, b; atol=NURBS_EPSILON) -> Bool

Element-wise approximate equality check for vectors or scalars.
"""
function allclose(a::AbstractVector, b::AbstractVector; atol=NURBS_EPSILON)
    length(a) == length(b) || return false
    return all(abs(a[i] - b[i]) <= atol for i in eachindex(a))
end

function allclose(a::Real, b::Real; atol=NURBS_EPSILON)
    return abs(a - b) <= atol
end
