"""
Abstract type for an elliptic curve "bundle" that includes a point type,
a scalar type, and operations on them.
"""
abstract type EllipticCurve end


"""
A subtype of [`EllipticCurve`](@ref) for curves defined by an equation `y^2 = x^3 + a * x + b`.
"""
abstract type WeierstrassCurve <: EllipticCurve end


"""
A subtype of [`WeierstrassCurve`](@ref) for curves of the form `y^2 = x^3 + b`
with the modulus `p = 1 mod 3`.
These curves support an endomorphism type 4.
"""
abstract type EndomorphismType4Curve <: WeierstrassCurve end


"""
Abstract type of a curve point.
"""
abstract type EllipticCurvePoint{C <: EllipticCurve} end


"""
Abstract type for built-in curve point types parameterized by coordinate types.
"""
abstract type StandardEllipticCurvePoint{C} <: EllipticCurvePoint{C} end


# Main public interface


"""
    curve_point_type(::Type{<:EllipticCurve})

Returns the type of a curve point (a subtype of [`EllipticCurvePoint`](@ref))
for an object of type [`EllipticCurve`](@ref).

Objects of this type support `zero()`, `iszero()`, `one()`, `rand()`,
`+`, `-`, `==`, `*` (by a scalar type).

Base point is returned by `one()`, and infinity point by `zero()`.
"""
curve_point_type(::Type{<:EllipticCurve}) = error("not implemented")


"""
    curve_scalar_type(::Type{<:EllipticCurve})

Returns the type of an associated scalar for an object of type [`EllipticCurve`](@ref).
All operations with objects of this type are done modulo curve order.

Objects of this type support `zero()`, `iszero()`, `one()`, `rand()`,
`+`, `-`, `==`, `*` (by a scalar type or a point type),
`iseven()`, `isodd()`, `>>`, `inv()`, `divrem()`, `trailing_zeros()`, `DarkIntegers.num_bits()`.
"""
curve_scalar_type(::Type{<:EllipticCurve}) = error("not implemented")


# Helper functions to construct default point and scalar types out of DarkInteger types

@generated function curve_scalar_type(
        ::Type{C}, ::Type{B}, ::Type{T}
        ) where {C <: EllipticCurve, B <: AbstractModUInt, T <: Unsigned}
    @assert bitsizeof(T) >= num_bits(curve_order(C))
    m = curve_order(C, T)
    tp = B{T, m}
    :( $tp )
end

@generated function curve_point_type(
        ::Type{C}, ::Type{P}, ::Type{B}, ::Type{T}
        ) where {C <: EllipticCurve, P <: StandardEllipticCurvePoint, B <: AbstractModUInt, T <: Unsigned}
    @assert bitsizeof(T) >= num_bits(curve_coordinate_modulus(C))
    m = curve_coordinate_modulus(C, T)
    coord_tp = B{T, m}
    point_tp = P{C, coord_tp}
    :( $point_tp )
end


# Generic methods, for use in generated functions and debugging.
# Can be slow, since they return BigInt numbers.

curve_order(::Type{<:EllipticCurve}) = error("not implemented")
curve_base(::Type{<:EllipticCurve}) = error("not implemented")
curve_coordinate_modulus(::Type{<:EllipticCurve}) = error("not implemented")


curve_weierstrass_coeff_a(::Type{<:WeierstrassCurve}) = error("not implemented")
curve_weierstrass_coeff_b(::Type{<:WeierstrassCurve}) = error("not implemented")


# Parameters for curves having endomorphism type 4
curve_endomorphism_type_4(::Type{<:EllipticCurvePoint{C}}) where C <: EndomorphismType4Curve =
    error("not implemented")
curve_endomorphism_type_4_lambda(::Type{<:EndomorphismType4Curve}) = error("not implemented")
curve_endomorphism_type_4_beta(::Type{<:EndomorphismType4Curve}) = error("not implemented")


@generated function curve_base(::Type{C}, ::Type{T}) where {C <: EllipticCurve, T}
    res = convert.(T, curve_base(C))
    :( $res )
end


@generated function curve_order(::Type{C}, ::Type{T}) where {C <: EllipticCurve, T}
    res = convert(T, curve_order(C))
    :( $res )
end


@generated function curve_coordinate_modulus(::Type{C}, ::Type{T}) where {C <: EllipticCurve, T}
    res = convert(T, curve_coordinate_modulus(C))
    :( $res )
end


@generated function curve_weierstrass_coeff_a(::Type{C}, ::Type{T}) where {C <: WeierstrassCurve, T}
    res = convert(T, curve_weierstrass_coeff_a(C))
    :( $res )
end


@generated function curve_weierstrass_coeff_b(::Type{C}, ::Type{T}) where {C <: WeierstrassCurve, T}
    res = convert(T, curve_weierstrass_coeff_b(C))
    :( $res )
end


@generated function curve_endomorphism_type_4_lambda(::Type{C}, ::Type{T}) where {C <: EndomorphismType4Curve, T}
    l = convert(T, curve_endomorphism_type_4_lambda(C))
    :( $l )
end


@generated function curve_endomorphism_type_4_beta(::Type{C}, ::Type{T}) where {C <: EndomorphismType4Curve, T}
    b = convert(T, curve_endomorphism_type_4_beta(C))
    :( $b )
end


# Generic functions


double(x) = x + x


triple(x) = double(x) + x


square(x) = x * x


Base.Broadcast.broadcastable(x::EllipticCurvePoint) = (x,)


function Base.:-(p::P, q::P) where P <: EllipticCurvePoint
    p + (-q)
end
