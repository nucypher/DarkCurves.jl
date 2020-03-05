abstract type EllipticCurve end


abstract type EllipticCurvePoint{C, T} end


curve_modulus(::Type{<:EllipticCurve}) = error("not implemented")
curve_order(::Type{<:EllipticCurve}) = error("not implemented")
curve_base(::Type{<:EllipticCurve}) = error("not implemented")
curve_coeff_a(::Type{<:EllipticCurve}) = error("not implemented")
curve_coeff_b(::Type{<:EllipticCurve}) = error("not implemented")

curve_endomorphism_lambda(::Type{<:EllipticCurve}) = error("not implemented")
curve_endomorphism_beta(::Type{<:EllipticCurve}) = error("not implemented")


@generated function curve_base(::Type{C}, ::Type{T}) where {C <: EllipticCurve, T}
    res = convert.(T, curve_base(C))
    :( $res )
end


@generated function curve_coeff_a(::Type{C}, ::Type{T}) where {C <: EllipticCurve, T}
    res = convert(T, curve_coeff_a(C))
    :( $res )
end


@generated function curve_coeff_b(::Type{C}, ::Type{T}) where {C <: EllipticCurve, T}
    res = convert(T, curve_coeff_b(C))
    :( $res )
end


function curve_scalar_type(
        ::Type{C}, ::Type{B}, ::Type{T}) where {C <: EllipticCurve, B <: AbstractModUInt, T <: Unsigned}
    @assert bitsizeof(T) >= num_bits(curve_modulus(C))
    B{T, convert(T, curve_modulus(C))}
end


# Generic functions


double(x) = x + x


triple(x) = double(x) + x


square(x) = x * x


Base.Broadcast.broadcastable(x::EllipticCurvePoint) = (x,)


function Base.:-(p::P, q::P) where P <: EllipticCurvePoint
    p + (-q)
end
