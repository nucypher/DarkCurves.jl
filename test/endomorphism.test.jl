using DarkCurves: curve_endomorphism_type_4_lambda, curve_endomorphism_type_4


@testgroup "Endomorphisms" begin


curve_types = [Curve_secp256k1] => ["SECP256k1"]
point_types =
    [DarkCurves.AffinePoint, DarkCurves.JacobianPoint, DarkCurves.ChudnovskyPoint] =>
    ["affine", "Jacobian", "Chudnovsky"]


@testcase "Endomorphism type 4" for curve_type in curve_types, point_type in point_types

    stp = curve_scalar_type(curve_type)
    ptp = curve_point_type(curve_type, point_type)

    p = one(ptp) * convert(stp, 123)
    l = curve_endomorphism_type_4_lambda(curve_type, MLUInt{4, UInt64})
    @test curve_endomorphism_type_4(p) == p * l

    @test curve_endomorphism_type_4(zero(ptp)) == zero(ptp)
end


end
