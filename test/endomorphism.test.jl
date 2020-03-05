using DarkCurves: curve_endomorphism_lambda, endomorphism


@testgroup "Endomorphisms" begin


curve_types = [Curve_secp256k1] => ["SECP256k1"]
point_types = [AffinePoint, JacobianPoint, ChudnovskyPoint] => ["affine", "Jacobian", "Chudnovsky"]


@testcase "Endomorphism function" for curve_type in curve_types, point_type in point_types
    stp = curve_scalar_type(curve_type, MgModUInt, MLUInt{4, UInt64})
    ref_ptp = AffinePoint{curve_type, stp}
    ptp = point_type{curve_type, stp}
    order = curve_order(curve_type)

    p = one(ptp) * 123
    l = curve_endomorphism_lambda(curve_type, MLUInt{4, UInt64})
    @test endomorphism(p) == p * l

    @test endomorphism(zero(ptp)) == zero(ptp)
end


end
