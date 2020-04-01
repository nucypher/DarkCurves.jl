@testgroup "Curve" begin


curve_types = [Curve_secp256k1] => ["SECP256k1"]
point_types = [AffinePoint, JacobianPoint, ChudnovskyPoint] => ["affine", "Jacobian", "Chudnovsky"]


@testcase "Basic operations" for curve_type in curve_types, point_type in point_types

    stp = curve_scalar_type(curve_type, MgModUInt, MLUInt{4, UInt64})
    ref_ptp = AffinePoint{curve_type, stp}
    ptp = point_type{curve_type, stp}
    order = curve_order(curve_type)

    b = one(ptp)
    b_ref = one(ref_ptp)

    # Trivial tests for point_type == AffinePoint
    @test double(b) == convert(ptp, double(b_ref))
    @test b + b == convert(ptp, b_ref + b_ref)
    @test double(b) + b == convert(ptp, double(b_ref) + b_ref)

    @test ref_mul(b + b, 23) == (b + b) * 23
    @test iszero(b * (order - 1) + b)
    @test iszero(b * order)
    @test b * (order + 1) == b

    @test double(double(b)) - b == double(b) + b
end


(@testcase tags=[:performance] "Addition performance" for
        curve_type in curve_types,
        point_type in point_types

    stp = curve_scalar_type(curve_type, MgModUInt, MLUInt{4, UInt64})
    ptp = point_type{curve_type, stp}

    b1 = one(ptp)
    b2 = double(b1)
    b4 = double(b2)

    trial = @benchmark $b2 + $b4
    @test_result benchmark_result(trial)
end)


mul_funcs = (
    [DarkCurves.mul_double_and_add, DarkCurves.mul_windowed, DarkCurves.mul_sliding_window, DarkCurves.mul_wnaf]
    => ["double-and-add", "windowed", "sliding window", "wNAF"])


(@testcase "Multiplication" for
        curve_type in curve_types,
        point_type in point_types,
        func in mul_funcs

    stp = curve_scalar_type(curve_type, MgModUInt, MLUInt{4, UInt64})
    ptp = point_type{curve_type, stp}

    b1 = one(ptp)
    p = b1 + b1

    @test func(p, 123) == ref_mul(p, 123)
end)


(@testcase tags=[:performance] "Scalar multiplication performance" for
        curve_type in curve_types

    stp = curve_scalar_type(curve_type, MgModUInt, MLUInt{2, UInt128})

    x_bi = 115047236638587805833081834189719086745649315857841928574581145752217906325686
    x = convert(stp, x_bi)

    trial = @benchmark $x * $x
    @test_result benchmark_result(trial)
end)


(@testcase tags=[:performance] "Multiplication performance" for
        curve_type in curve_types,
        point_type in point_types,
        func in mul_funcs

    stp = curve_scalar_type(curve_type, MgModUInt, MLUInt{2, UInt128})
    ptp = point_type{curve_type, stp}

    b1 = one(ptp)
    b2 = double(b1)

    x_bi = 115047236638587805833081834189719086745649315857841928574581145752217906325686
    x = convert(MLUInt{2, UInt128}, x_bi)

    trial = @benchmark $func($b2, $x)
    @test_result benchmark_result(trial)
end)


end
