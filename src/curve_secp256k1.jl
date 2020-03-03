struct Curve_secp256k1 <: EllipticCurve end

curve_modulus(::Type{Curve_secp256k1}) = big(2)^256 - 2^32 - 2^9 - 2^8 - 2^7 - 2^6 - 2^4 - 1
curve_order(::Type{Curve_secp256k1}) = as_builtin(MLUInt{4, UInt64}(reverse((
    0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFE, 0xBAAEDCE6AF48A03B, 0xBFD25E8CD0364141))))
curve_base(::Type{Curve_secp256k1}) = (
        as_builtin(MLUInt{4, UInt64}(reverse((
                0x79BE667EF9DCBBAC, 0x55A06295CE870B07, 0x029BFCDB2DCE28D9, 0x59F2815B16F81798)))),
        as_builtin(MLUInt{4, UInt64}(reverse((
                0x483ADA7726A3C465, 0x5DA4FBFC0E1108A8, 0xFD17B448A6855419, 0x9C47D08FFB10D4B8)))))
curve_coeff_a(::Type{Curve_secp256k1}) = 0
curve_coeff_b(::Type{Curve_secp256k1}) = 7

curve_endomorphism_lambda(::Type{Curve_secp256k1}) =
    37718080363155996902926221483475020450927657555482586988616620542887997980018
curve_endomorphism_beta(::Type{Curve_secp256k1}) =
    55594575648329892869085402983802832744385952214688224221778511981742606582254
