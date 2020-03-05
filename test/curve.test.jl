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


(@testcase tags=[:performance] "Multiplication" for
        curve_type in curve_types,
        point_type in point_types,
        func in mul_funcs

    stp = curve_scalar_type(curve_type, MgModUInt, MLUInt{4, UInt64})
    ptp = point_type{curve_type, stp}

    b1 = one(ptp)
    p = b1 + b1

    @test func(p, 123) == ref_mul(p, 123)
end)


(@testcase tags=[:performance] "Multiplication performance" for
        curve_type in curve_types,
        point_type in point_types,
        func in mul_funcs

    stp = curve_scalar_type(curve_type, MgModUInt, MLUInt{4, UInt64})
    ptp = point_type{curve_type, stp}

    b1 = one(ptp)
    b2 = double(b1)

    x_bi = 115047236638587805833081834189719086745649315857841928574581145752217906325686
    x = convert(MLUInt{4, UInt64}, x_bi)

    trial = @benchmark $func($b2, $x)
    @test_result benchmark_result(trial)
end)


function prepare_lincomb_dataset(curve_type, point_type, len)
    stp = curve_scalar_type(curve_type, MgModUInt, MLUInt{4, UInt64})
    ptp = point_type{curve_type, stp}

    rng = MersenneTwister(123)
    point_vals = rand(rng, big(1):curve_order(curve_type), len)
    coeffs_bi = rand(rng, big(0):curve_modulus(curve_type)-1, len)

    b1 = one(ptp)
    points = b1 .* point_vals

    coeffs = convert.(MLUInt{4, UInt64}, coeffs_bi)

    points, coeffs
end


lincomb_funcs = (
    [DarkCurves.lincomb_windowed]
    => ["windowed"])
fast_point_types = [JacobianPoint, ChudnovskyPoint] => ["Jacobian", "Chudnovsky"]


(@testcase "Linear combination" for
        curve_type in curve_types,
        point_type in point_types,
        func in lincomb_funcs

    points, coeffs = prepare_lincomb_dataset(curve_type, point_type, 16)
    ref = sum(points .* coeffs)
    res = func(points, coeffs)
    @test ref == res
end)


(@testcase tags=[:performance] "Linear combination baseline performance" for
        curve_type in curve_types,
        point_type in fast_point_types

    points, coeffs = prepare_lincomb_dataset(curve_type, point_type, 1024)
    trial = @benchmark sum($points .* $coeffs)
    @test_result benchmark_result(trial)
end)


(@testcase tags=[:performance] "Linear combination performance" for
        curve_type in curve_types,
        point_type in fast_point_types,
        func in lincomb_funcs,
        width in [3, 4, 5]

    points, coeffs = prepare_lincomb_dataset(curve_type, point_type, 1024)
    trial = @benchmark $func($points, $coeffs, $width)
    @test_result benchmark_result(trial)
end)

end
