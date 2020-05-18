module DarkCurves

using DarkIntegers
using Random

include("types.jl")
export curve_point_type
export curve_scalar_type
export EllipticCurve
export EllipticCurvePoint

include("point_affine.jl")
include("point_jacobian.jl")
include("point_chudnovsky.jl")

include("mul.jl")

include("balanced_decomposition.jl")

include("batch_operations.jl")
export lin_comb
export batch_mul

include("random.jl")

include("curve_secp256k1.jl")
export Curve_secp256k1

end
