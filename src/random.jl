struct CurvePointSampler{P} <: Random.Sampler{P}
    base_powers :: Array{P, 1}

    function CurvePointSampler{P}() where P <: EllipticCurvePoint{C, T} where {C, T}
        base_powers = Array{P}(undef, num_bits(curve_order(C)))
        base_powers[1] = one(P)
        for i in 2:length(base_powers)
            base_powers[i] = double(base_powers[i - 1])
        end

        new{P}(base_powers)
    end
end


const _SAMPLERS = Dict{Type, CurvePointSampler}()


function Random.Sampler(
        RNG::Type{<:AbstractRNG}, ::Type{P}, n::Union{Val{1}, Val{Inf}}) where P <: EllipticCurvePoint
    if haskey(_SAMPLERS, P)
        _SAMPLERS[P]
    else
        sampler = CurvePointSampler{P}()
        _SAMPLERS[P] = CurvePointSampler{P}()
        sampler
    end
end


function Base.rand(rng::AbstractRNG, sampler::CurvePointSampler{P}) where P <: EllipticCurvePoint{C, T} where {C, T}
    order = value(curve_order(C, T))
    pwr = rand(rng, one(order):order)

    res = zero(P)
    for i in 1:num_bits(curve_order(C))
        if isodd(pwr)
            res += sampler.base_powers[i]
        end
        pwr >>= 1
    end

    res
end
