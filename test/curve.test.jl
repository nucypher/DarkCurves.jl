@testgroup "Curve" begin


curve_types = [Curve_secp256k1] => ["SECP256k1"]
point_types =
    [DarkCurves.AffinePoint, DarkCurves.JacobianPoint, DarkCurves.ChudnovskyPoint] =>
    ["affine", "Jacobian", "Chudnovsky"]


@testcase "Basic operations" for curve_type in curve_types, point_type in point_types

    stp = curve_scalar_type(curve_type, ModUInt, MLUInt{4, UInt64})
    ptp = curve_point_type(curve_type, point_type, MgModUInt, MLUInt{4, UInt64})
    ref_ptp = curve_point_type(curve_type, DarkCurves.AffinePoint, MgModUInt, MLUInt{4, UInt64})

    b = one(ptp)
    b_ref = one(ref_ptp)

    # Trivial tests for point_type == AffinePoint
    @test double(b) == convert(ptp, double(b_ref))
    @test b + b == convert(ptp, b_ref + b_ref)
    @test double(b) + b == convert(ptp, double(b_ref) + b_ref)

    @test ref_mul(b + b, 23) == (b + b) * convert(stp, 23)

    @test iszero(b * zero(stp))
    @test b * one(stp) == b
    @test iszero(b * (-one(stp)) + b)

    @test double(double(b)) - b == double(b) + b
end


(@testcase tags=[:performance] "Addition performance" for
        curve_type in curve_types,
        point_type in point_types

    ptp = curve_point_type(curve_type, point_type)

    b1 = one(ptp)
    b2 = double(b1)
    b4 = double(b2)

    trial = @benchmark $b2 + $b4
    @test_result benchmark_result(trial)
end)


mul_funcs = (
    [DarkCurves.mul_double_and_add, DarkCurves.mul_windowed,
    DarkCurves.mul_sliding_window, DarkCurves.mul_wnaf, DarkCurves.mul_endomorphism_wnaf]
    => ["double-and-add", "windowed", "sliding window", "wNAF", "endomorphism + wNAF"])


(@testcase "Multiplication" for
        curve_type in curve_types,
        point_type in point_types,
        func in mul_funcs

    ptp = curve_point_type(curve_type, point_type)

    b1 = one(ptp)
    p = b1 + b1

    x_bi = 123
    x = convert(MLUInt{2, UInt128}, x_bi)

    @test func(p, x) == ref_mul(p, x)
end)


(@testcase tags=[:performance] "Multiplication performance" for
        curve_type in curve_types,
        point_type in point_types,
        func in mul_funcs

    stp = curve_scalar_type(curve_type)
    ptp = curve_point_type(curve_type, point_type)

    b1 = one(ptp)
    b2 = double(b1)

    x_bi = 115047236638587805833081834189719086745649315857841928574581145752217906325686
    x = convert(stp, x_bi)

    @test func(b2, x) == DarkCurves.mul_double_and_add(b2, x)

    trial = @benchmark $func($b2, $x)
    @test_result benchmark_result(trial)
end)


end
