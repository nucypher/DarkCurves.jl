function gcd_eqs(a::T, b::T) where T
    u = a
    v = b
    x1 = one(T)
    y1 = zero(T)
    x2 = zero(T)
    y2 = one(T)

    coeffs = [(x1, y1, a)]

    while !iszero(u)
        q, r = divrem(v, u)
        x = x2 - q * x1
        y = y2 - q * y1

        push!(coeffs, (x, y, r))

        v = u
        u = r
        x2 = x1
        x1 = x
        y2 = y1
        y1 = y
    end

    d = v
    x = x2
    y = y2

    coeffs
end


function decomposition_coeffs(n::T, l::T) where T
    coeffs = gcd_eqs(n, l)
    sqrt_n = ceil(T, sqrt(n))

    idx = 0
    for (i, triple) in enumerate(coeffs)
        _, _, r = triple
        if r < sqrt_n
            idx = i - 1
            break
        end
    end

    _, ti0, ri0 = coeffs[idx]
    _, ti1, ri1 = coeffs[idx+1]
    _, ti2, ri2 = coeffs[idx+2]

    a1, b1 = ri1, ti1
    if ri0 * ri0 + ti0 * ti0 <= ri2 * ri2 + ti2 * ti2
        a2, b2 = ri0, ti0
    else
        a2, b2 = ri2, ti2
    end

    a1, b1, a2, b2
end


struct DecompositionCoeffs{T <: Unsigned}
    a1 :: T
    b1 :: T
    a2 :: T
    b2 :: T
    n :: T

    b1_signbit :: Bool
    b2_signbit :: Bool

    function DecompositionCoeffs{T}(a1::V, b1::V, a2::V, b2::V, n::V) where {T, V}
        new{T}(a1, abs(b1), a2, abs(b2), n, signbit(b1), signbit(b2))
    end
end


@generated function get_decomposition_coeffs(
        ::Type{T}, ::Type{C}) where {T, C <: EndomorphismType4}

    l = DarkCurves.curve_endomorphism_lambda(C)
    n = curve_order(C)
    a1, b1, a2, b2 = decomposition_coeffs(n, l)
    dc_typed = DecompositionCoeffs{T}(a1, b1, a2, b2, n)

    :( $dc_typed )
end


apply_signbit(x, s) = s ? -x : x


"""
Returns `k1, k2` such that `k = k1 + k2 * l mod n`.
"""
function balanced_decomposition(::Type{C}, k::T) where {C <: EndomorphismType4, T}
    dc = get_decomposition_coeffs(T, C)

    c1 = convert(T, widemul(dc.b2, k) ÷ dc.n)
    c2 = convert(T, widemul(dc.b1, k) ÷ dc.n)

    p1 = c2 * dc.a2
    p2 = c1 * dc.a1
    k1 = k - apply_signbit(p1, dc.b1_signbit) + apply_signbit(p2, dc.b2_signbit)

    p1 = c2 * dc.b2
    p2 = c1 * dc.b1
    p_signbit = p2 > p1

    k2 = apply_signbit(p1 - p2, p_signbit)

    k2_signbit = xor(xor(dc.b1_signbit, dc.b2_signbit), p_signbit)

    k1, k2, k2_signbit
end
