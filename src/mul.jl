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


function Base.:*(p::P, y::Union{ModUInt, MgModUInt}) where {P <: EllipticCurvePoint}
    p * value(y)
end


function Base.:*(p::P, y::V) where {P <: EllipticCurvePoint, V <: Union{Integer, BigInt}}
    mul_wnaf(p, y)
end
