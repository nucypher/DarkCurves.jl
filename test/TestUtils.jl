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


function prepare_lin_comb_dataset(curve_type, point_type, len)
    stp = curve_scalar_type(curve_type)
    ptp = curve_point_type(curve_type, point_type)

    rng = MersenneTwister(123)
    points = rand(rng, ptp, len)
    coeffs = rand(rng, stp, len)

    points, coeffs
end
