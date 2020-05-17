function repeated_double(p, m)
    for i in 1:m
        p = double(p)
    end
    p
end


function mul_double_and_add(p::P, y::V) where {P <: EllipticCurvePoint, V <: Union{Integer, BigInt}}
    if iszero(y)
        zero(P)
    elseif isone(y)
        p
    else
        acc = zero(P)
        while true
            if isodd(y)
                acc += p
            end
            if iszero(y)
                break
            else
                p = double(p)
                y >>= 1
            end
        end
        acc
    end
end


function mul_windowed(p::P, y::V, w::Int=4) where {P <: EllipticCurvePoint, V <: Union{Integer, BigInt}}

    if iszero(y)
        return zero(P)
    elseif isone(y)
        return p
    end

    precomp = Array{P}(undef, 1 << w)
    precomp[1] = p
    for i in 2:(1 << w)
        precomp[i] = p + precomp[i-1]
    end

    nb = num_bits(y)

    acc = zero(P)
    for i in (nb ÷ w):-1:0
        acc = repeated_double(acc, w)
        d = y >> (i * w)
        if !iszero(d)
            acc += precomp[d]
        end
        y = y - (d << (i * w))
    end
    acc
end


function mul_sliding_window(p::P, y::V, w::Int=4) where {P <: EllipticCurvePoint, V <: Union{Integer, BigInt}}

    if iszero(y)
        return zero(P)
    elseif isone(y)
        return p
    end

    precomp = Array{P}(undef, 1 << (w - 1))
    precomp[1] = repeated_double(p, w - 1)
    for i in 2:(1 << (w - 1))
        precomp[i] = p + precomp[i-1]
    end

    nb = num_bits(y)

    acc = zero(P)
    while true
        new_nb = num_bits(y)
        acc = repeated_double(acc, new_nb < w ? nb - w : nb - new_nb)

        if new_nb < w
            return acc + p * y
        end

        t = y >> (new_nb - w)
        acc += precomp[t - (1 << (w - 1) - 1)]

        nb = new_nb
        y = y - (t << (new_nb - w))
    end
    acc
end


function get_wnaf(val::T, w::Int) where T
    mask = (one(T) << w) - one(T)
    res = Int[]

    tz = trailing_zeros(val)
    val >>= tz
    for j in 1:tz
        push!(res, 0)
    end

    while !iszero(val)
        tz = trailing_zeros(val)
        val >>= tz
        for j in 1:tz-1
            push!(res, 0)
        end

        di = convert(Int, val & mask)
        if di >= 1 << (w - 1)
            val += ((1 << w) - di)
        else
            val -= di
        end

        push!(res, di)
    end
    res
end


function mul_wnaf(p::P, y::V, w::Int=4) where {P <: EllipticCurvePoint, V <: Union{Integer, BigInt}}

    if iszero(y)
        return zero(P)
    elseif isone(y)
        return p
    end

    l = 1 << (w - 1)
    precomp = Array{P}(undef, l) # corresponds to di = [1, 3, 5, ..., 2^(w-1)-1, -2^(w-1)-1, ..., -3, -1]

    dp = double(p)
    precomp[1] = p
    precomp[end] = -p

    for i in 2:(l>>1)
        precomp[i] = precomp[i-1] + dp
        precomp[end-i+1] = -precomp[i]
    end

    ds = get_wnaf(y, w)

    acc = zero(P)
    for idx in length(ds):-1:1
        acc = double(acc)
        if !iszero(ds[idx])
            acc += precomp[(ds[idx] >> 1) + 1]
        end
    end

    acc
end


function mul_endomorphism_wnaf(
        p::P, y::V, w::Int=4
        ) where {P <: EllipticCurvePoint{C, T}, V <: Union{Integer, BigInt}} where {C, T}

    w1 = w
    w2 = w

    k1, k2, k2_signbit = balanced_decomposition(C, y)

    if iszero(k1)
        return apply_signbit(mul_wnaf(curve_endomorphism_type_4(p), k2), k2_signbit)
    elseif iszero(k2)
        return mul_wnaf(p, k1)
    end

    p2 = apply_signbit(endomorphism(p), k2_signbit)

    # corresponds to di = [1, 3, 5, ..., 2^(w-1)-1, -2^(w-1)-1, ..., -3, -1]
    l1 = 1 << (w1 - 1)
    precomp1 = Array{P}(undef, l1)
    l2 = 1 << (w2 - 1)
    precomp2 = Array{P}(undef, l2)

    dp1 = double(p)
    precomp1[1] = p
    precomp1[end] = -p

    for i in 2:(l1>>1)
        precomp1[i] = precomp1[i-1] + dp1
        precomp1[end-i+1] = -precomp1[i]
    end

    dp2 = double(p2)
    precomp2[1] = p2
    precomp2[end] = -p2

    for i in 2:(l2>>1)
        precomp2[i] = precomp2[i-1] + dp2
        precomp2[end-i+1] = -precomp2[i]
    end

    ds1 = get_wnaf(k1, w1)
    ds2 = get_wnaf(k2, w2)

    if length(ds1) > length(ds2)
        ds2 = vcat(ds2, zeros(Int, length(ds1) - length(ds2)))
    elseif length(ds2) > length(ds1)
        ds1 = vcat(ds1, zeros(Int, length(ds2) - length(ds1)))
    end

    acc = zero(P)
    for idx in length(ds1):-1:1
        acc = double(acc)
        if !iszero(ds1[idx])
            acc += precomp1[(ds1[idx] >> 1) + 1]
        end
        if !iszero(ds2[idx])
            acc += precomp2[(ds2[idx] >> 1) + 1]
        end
    end

    acc
end


@inline function Base.:*(
        p::P, y::Union{ModUInt, MgModUInt}
        ) where {P <: EllipticCurvePoint{C, T}} where {C <: EndomorphismType4, T}
    p * value(y)
end


@inline function Base.:*(
        p::P, y::Union{Integer, BigInt}
        ) where {P <: EllipticCurvePoint{C, T}} where {C <: EndomorphismType4, T}
    mul_endomorphism_wnaf(p, y)
end


@inline function Base.:*(p::P, y::Union{ModUInt, MgModUInt}) where {P <: EllipticCurvePoint}
    p * value(y)
end


@inline function Base.:*(p::P, y::Union{Integer, BigInt}) where {P <: EllipticCurvePoint}
    #record_curve_muls!()
    mul_wnaf(p, y)
end
