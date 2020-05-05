@testgroup "Batch operations" begin


curve_types = [Curve_secp256k1] => ["SECP256k1"]
fast_point_types = [JacobianPoint, ChudnovskyPoint] => ["Jacobian", "Chudnovsky"]

lin_comb_funcs = (
    [DarkCurves.lin_comb_windowed]
    => ["windowed"])


(@testcase "Linear combination" for
        curve_type in curve_types,
        point_type in fast_point_types,
        func in lin_comb_funcs

    points, coeffs = prepare_lin_comb_dataset(curve_type, point_type, 16)
    ref = sum(points .* coeffs)
    res = func(points, coeffs)
    @test ref == res
end)


(@testcase tags=[:performance] "Linear combination baseline performance" for
        curve_type in curve_types,
        point_type in fast_point_types

    points, coeffs = prepare_lin_comb_dataset(curve_type, point_type, 1024)
    trial = @benchmark sum($points .* $coeffs)
    @test_result benchmark_result(trial)
end)


(@testcase tags=[:performance] "Linear combination performance" for
        curve_type in curve_types,
        point_type in fast_point_types,
        func in lin_comb_funcs,
        width in [3, 4, 5]

    points, coeffs = prepare_lin_comb_dataset(curve_type, point_type, 1024)
    trial = @benchmark $func($points, $coeffs, $width)
    @test_result benchmark_result(trial)
end)


batch_mul_funcs = (
    [DarkCurves.batch_mul_endomorhism_wnaf, DarkCurves.batch_mul_addition_chain, DarkCurves.batch_mul_wnaf]
    => ["endomorphism+wNAF", "addition chain", "wNAF"])


(@testcase "Batch multiplication" for
        curve_type in curve_types,
        point_type in fast_point_types,
        func in batch_mul_funcs

    points, coeffs = prepare_lin_comb_dataset(curve_type, point_type, 16)
    ref = points .* coeffs[1]
    res = func(points, coeffs[1])
    @test ref == res
end)


batch_reference(points, coeff) = points .* coeff


(@testcase tags=[:performance] "Batch multiplication baseline performance" for
        curve_type in curve_types,
        point_type in fast_point_types

    points, coeffs = prepare_lin_comb_dataset(curve_type, point_type, 1024)
    coeff = coeffs[1]

    trial = @benchmark $batch_reference($points, $coeff)
    @test_result benchmark_result(trial)
end)


(@testcase tags=[:performance] "Batch multiplication performance" for
        curve_type in curve_types,
        point_type in fast_point_types,
        func in batch_mul_funcs,
        width in [3, 4, 5]

    points, coeffs = prepare_lin_comb_dataset(curve_type, point_type, 1024)
    coeff = coeffs[1]

    trial = @benchmark $func($points, $coeff, $width)
    @test_result benchmark_result(trial)
end)


end
