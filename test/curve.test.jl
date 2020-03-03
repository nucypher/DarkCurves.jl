using DarkIntegers
using BenchmarkTools


function benchmark_result(trial)
    time_str = BenchmarkTools.prettytime(minimum(trial.times))

    if trial.allocs > 0
        mem_str = BenchmarkTools.prettymemory(trial.memory)
        alloc_str = ", $mem_str ($(trial.allocs) allocs)"
    else
        alloc_str = ""
    end

    time_str * alloc_str
end


function ref_mul(p, y)
    res = p
    for i in 2:y
        res += p
    end
    res
end


double(x) = x + x


@testgroup "Curve" begin


point_types = [AffinePoint, JacobianPoint, ChudnovskyPoint] => ["affine", "Jacobian", "Chudnovsky"]


@testcase "SECP256k1" for point_type in point_types

    stp = curve_scalar_type(Curve_secp256k1, MgModUInt, MLUInt{4, UInt64})
    ref_ptp = AffinePoint{Curve_secp256k1, stp}
    ptp = point_type{Curve_secp256k1, stp}
    order = curve_order(Curve_secp256k1)

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
end


@testcase tags=[:performance] "SECP256k1, addition performance" for point_type in point_types

    stp = curve_scalar_type(Curve_secp256k1, MgModUInt, MLUInt{4, UInt64})
    ptp = point_type{Curve_secp256k1, stp}

    b1 = one(ptp)
    b2 = double(b1)
    b4 = double(b2)

    trial = @benchmark $b2 + $b4
    @test_result benchmark_result(trial)
end


mul_funcs = (
    [DarkCurves.mul_double_and_add, DarkCurves.mul_windowed, DarkCurves.mul_sliding_window, DarkCurves.mul_wnaf]
    => ["double-and-add", "windowed", "sliding window", "wNAF"])


@testcase tags=[:performance] "SECP256k1, multiplication performance" for point_type in point_types, func in mul_funcs

    stp = curve_scalar_type(Curve_secp256k1, MgModUInt, MLUInt{4, UInt64})
    ptp = point_type{Curve_secp256k1, stp}

    b1 = one(ptp)
    b2 = double(b1)

    x_bi = 115047236638587805833081834189719086745649315857841928574581145752217906325686
    x = convert(MLUInt{4, UInt64}, x_bi)

    trial = @benchmark $func($b2, $x)
    @test_result benchmark_result(trial)
end


function prepare_lincomb_dataset(point_type, len)
    stp = curve_scalar_type(Curve_secp256k1, MgModUInt, MLUInt{4, UInt64})
    ptp = point_type{Curve_secp256k1, stp}

    rng = MersenneTwister(123)
    point_vals = rand(rng, big(1):curve_order(Curve_secp256k1), len)
    coeffs_bi = rand(rng, big(0):curve_modulus(Curve_secp256k1)-1, len)

    b1 = one(ptp)
    points = b1 .* point_vals

    coeffs = convert.(MLUInt{4, UInt64}, coeffs_bi)

    points, coeffs
end


lincomb_funcs = (
    [DarkCurves.lincomb_windowed]
    => ["windowed"])
fast_point_types = [JacobianPoint, ChudnovskyPoint] => ["Jacobian", "Chudnovsky"]


@testcase "SECP256k1, linear combination" for point_type in point_types, func in lincomb_funcs
    points, coeffs = prepare_lincomb_dataset(point_type, 16)
    ref = sum(points .* coeffs)
    res = func(points, coeffs)
    @test ref == res
end


(@testcase tags=[:performance] "SECP256k1, linear combination baseline performance" for
        point_type in fast_point_types

    points, coeffs = prepare_lincomb_dataset(point_type, 1024)
    trial = @benchmark sum($points .* $coeffs)
    @test_result benchmark_result(trial)
end)


(@testcase tags=[:performance] "SECP256k1, linear combination performance" for
        point_type in fast_point_types,
        func in lincomb_funcs,
        width in [3, 4, 5]

    points, coeffs = prepare_lincomb_dataset(point_type, 1024)
    trial = @benchmark $func($points, $coeffs, $width)
    @test_result benchmark_result(trial)
end)

end
