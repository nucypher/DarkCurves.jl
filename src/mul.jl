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


struct NAFIterator{T}
    val :: T
    window :: Int
end


Base.IteratorSize(::Type{NAFIterator{T}}) where T = Base.SizeUnknown()
Base.IteratorEltype(::Type{NAFIterator{T}}) where T = Base.HasEltype()
Base.eltype(::Type{NAFIterator{T}}) where T = Tuple{Int, Int}


struct NAFIteratorState{T}
    val :: T
    idx :: Int
end


function next(iter::NAFIterator{T}, state::NAFIteratorState{T}) where T
    val = state.val
    idx = state.idx
    w = iter.window

    tz = trailing_zeros(val)
    val >>= tz
    idx += tz

    di = convert(Int, val & ((one(T) << w) - one(T)))
    if di >= 1 << (w - 1)
        val += ((1 << w) - di)
    else
        val -= di
    end

    new_state = NAFIteratorState{T}(val, idx)
    (idx, di), new_state
end


function Base.iterate(iter::NAFIterator{T}, state::Union{Nothing, NAFIteratorState{T}}=nothing) where T
    if state === nothing
        state = NAFIteratorState{T}(iter.val, 0)
    end

    if iszero(state.val)
        nothing
    else
        next(iter, state)
    end
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

    ds = collect(NAFIterator(y, w))

    acc = zero(P)
    for j in length(ds):-1:1
        i, di = ds[j]
        acc = acc + precomp[(di >> 1) + 1]
        acc = repeated_double(acc, j == 1 ? i : i - ds[j-1][1])
    end

    acc
end


function Base.:*(p::P, y::Union{ModUInt, MgModUInt}) where {P <: EllipticCurvePoint}
    p * value(y)
end


function Base.:*(p::P, y::V) where {P <: EllipticCurvePoint, V <: Union{Integer, BigInt}}
    mul_wnaf(p, y)
end
