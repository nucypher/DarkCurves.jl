module DarkCurves

using DarkIntegers

include("types.jl")
export curve_order
export curve_base
export curve_modulus
export curve_coeff_a
export curve_coeff_b
export curve_scalar_type
export EllipticCurve
export EllipticCurvePoint

include("point_affine.jl")
export AffinePoint

include("point_jacobian.jl")
export JacobianPoint

include("point_chudnovsky.jl")
export ChudnovskyPoint

include("mul.jl")

include("lincomb.jl")
export lincomb

include("curve_secp256k1.jl")
export Curve_secp256k1

end
