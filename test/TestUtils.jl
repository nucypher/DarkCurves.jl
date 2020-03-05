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
