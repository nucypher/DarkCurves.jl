"""
Linear combination using Shamir's trick and a fixed window.
"""
function lin_comb_windowed(
        points::Array{P, 1}, coeffs::Array{T, 1}, w::Int=4) where {P, T}

    @assert length(points) == length(coeffs)

    precomp = Array{P}(undef, length(points), 1 << w)
    for j in 1:length(points)
        precomp[j, 1] = points[j]
    end
    for i in 2:(1 << w)
        for j in 1:length(points)
            precomp[j, i] = points[j] + precomp[j, i-1]
        end
    end

    nb = maximum(num_bits.(coeffs))
    ys = copy(coeffs)

    acc = zero(P)
    for i in (nb ÷ w):-1:0
        acc = repeated_double(acc, w)
        for j in 1:length(points)
            d = ys[j] >> (i * w)
            di = d % Int
            if !iszero(d)
                acc += precomp[j, di]
            end
            ys[j] -= d << (i * w)
        end
    end
    acc
end


# For each subset in `subsets` (provided as a list of indices into `numbers`),
# compute the sum of that subset of `numbers`. More efficient than the naive method.
function multisubset(numbers::Array{P, 1}, subsets_list::Array{Set{Int}, 1}) where P
    numbers = copy(numbers)
    subsets = Dict((i-1) => copy(subset) for (i, subset) in enumerate(subsets_list))
    output = zeros(P, length(subsets_list))

    for roundcount in 1:9999999
        # Compute counts of every pair of indices in the subset list
        pair_count = Dict{Tuple{Int, Int}, Int}()
        for (index, subset) in subsets
            for x in subset
                for y in subset
                    if y > x
                        pair_count[(x, y)] = get(pair_count, (x, y), 0) + 1
                    end
                end
            end
        end

        # Determine pairs with highest count. The cutoff parameter [:len(numbers)]
        # determines a tradeoff between group operation count and other forms of overhead
        cutoff = min(length(pair_count), length(numbers) * trunc(Int, log(length(numbers))))
        pairs_by_count = collect(keys(pair_count))
        sort!(pairs_by_count; alg=PartialQuickSort(cutoff), by=el->pair_count[el], rev=true)
        pairs_by_count = pairs_by_count[1:cutoff]

        # Exit condition: all subsets have size 1, no pairs
        if length(pairs_by_count) == 0
            for (key, subset) in subsets
                for index in subset
                    output[key+1] += numbers[index+1]
                end
            end
            return output
        end

        # In each of the highest-count pairs, take the sum of the numbers at those indices,
        # and add the result as a new value, and modify `subsets` to include the new value
        # wherever possible
        used = Set{Int}()
        for (maxx, maxy) in pairs_by_count
            if maxx in used || maxy in used
                continue
            end
            push!(used, maxx)
            push!(used, maxy)
            push!(numbers, numbers[maxx+1] + numbers[maxy+1])
            for (key, subset) in subsets
                if maxx in subset && maxy in subset
                    delete!(subset, maxx)
                    delete!(subset, maxy)
                    if isempty(subset)
                        output[key+1] = numbers[end]
                        delete!(subsets, key)
                    else
                        push!(subset, length(numbers)-1)
                    end
                end
            end
        end
    end
end


function lin_comb_vitalik(numbers::Array{P, 1}, factors::Array{T, 1}) where {P, T}
    # Maximum bit length of a number; how many subsets we need to make
    maxbitlen = maximum(num_bits.(factors)) # TODO: -1? otherwise the last subset is empty
    # Compute the subsets: the ith subset contains the numbers whose corresponding factor
    # has a 1 at the ith bit
    subsets = [Set(i-1 for i in 1:length(numbers) if isodd(factors[i] >> j)) for j in 0:maxbitlen]
    subset_sums = multisubset(numbers, subsets)
    # For example, suppose a value V has factor 6 (011 in increasing-order binary). Subset 0
    # will not have V, subset 1 will, and subset 2 will. So if we multiply the output of adding
    # subset 0 with twice the output of adding subset 1, with four times the output of adding
    # subset 2, then V will be represented 0 + 2 + 4 = 6 times. This reasoning applies for every
    # value. So `subset_0_sum + 2 * subset_1_sum + 4 * subset_2_sum` gives us the result we want.
    # Here, we compute this as `((subset_2_sum * 2) + subset_1_sum) * 2 + subset_0_sum` for
    # efficiency: an extra `maxbitlen * 2` group operations.
    o = zero(P)
    for i in length(subsets):-1:1
        o = double(o) + subset_sums[i]
    end
    o
end


"""
Calculates a linear combination of curve points
given an array of points and an array of coefficients.
"""
function lin_comb(
        points::Array{P, 1}, coeffs::Array{T, 1}, w::Int=4) where {P <: EllipticCurvePoint, T <: Integer}
    #record_curve_muls!(length(points))
    lin_comb_windowed(points, coeffs, w)
end


@Base.propagate_inbounds function batch_mul_addition_chain(x::Array{P, 1}, n::T, log_m::Int=4) where {P, T}

    mod_m = (1 << log_m) - 1

    N = n
    Z = copy(x)
    Ys = zeros(P, length(x), 1 << (log_m - 1))

    while true

        # At every step here
        #     x * n == Z * N + sum(Ys .* collect(1:2:2^log_m))

        k = (N & mod_m) % Int
        N >>= log_m

        if k == 0

            for i in 1:length(Z)
                Z[i] = repeated_double(Z[i], log_m)
            end
            #Z .= repeated_double.(Z, log_m)
            continue
        else
            p = trailing_zeros(k)
            q = k >> p

            for i in 1:length(Z)
                Z[i] = repeated_double(Z[i], p)
            end
            #Z .= repeated_double.(Z, p)

            for i in 1:length(Z)
                Ys[i, (q >> 1) + 1] += Z[i]
            end
            #Ys[(q >> 1) + 1,:] .+= Z
        end

        if N > 0
            for i in 1:length(Z)
                Z[i] = repeated_double(Z[i], log_m - p)
            end
            #Z .= repeated_double.(Z, log_m - p)
        else
            break
        end
    end

    for i in (1<<(log_m-1))-1:-1:1
        for j in 1:length(Z)
            Ys[j,i] += Ys[j,i+1]
        end
        #Ys[i,:] .+= Ys[i+1,:]
    end

    for j in 1:length(Z)
        Z[j] = Ys[j,2]
    end

    for k in 3:(1<<(log_m-1))
        for j in 1:length(Z)
            Z[j] += Ys[j,k]
        end
    end

    for j in 1:length(Z)
        Z[j] = double(Z[j]) + Ys[j,1]
    end

    Z
    #Ys[1,:] .+ double.(sum(Ys[2:end,:], dims=1)[:])
end


function batch_mul_wnaf(points::Array{P, 1}, y::T, w::Int=4) where {P <: EllipticCurvePoint, T <: Union{Integer, BigInt}}

    if iszero(y)
        return zeros(P, length(points))
    elseif isone(y)
        return points
    end

    l = 1 << (w - 1)
    precomp = Array{P}(undef, length(points), l) # corresponds to di = [1, 3, 5, ..., 2^(w-1)-1, -2^(w-1)-1, ..., -3, -1]

    dpoints = double.(points)

    for j in 1:length(points)
        precomp[j,1] = points[j]
        precomp[j,end] = -points[j]
    end

    for i in 2:(l>>1)
        for j in 1:length(points)
            precomp[j,i] = precomp[j,i-1] + dpoints[j]
            precomp[j,end-i+1] = -precomp[j,i]
        end
    end

    ds = get_wnaf(y, w)

    acc = zeros(P, length(points))
    for idx in length(ds):-1:1
        for j in 1:length(points)
            acc[j] = double(acc[j])
        end

        if !iszero(ds[idx])
            for j in 1:length(points)
                acc[j] += precomp[j,(ds[idx] >> 1) + 1]
            end
        end

    end

    acc
end


function batch_mul_endomorphism_wnaf(
        points::Array{P, 1}, coeff::T, w::Int=4,
        ) where {P <: EllipticCurvePoint{C}, T <: Integer} where C

    w1 = w
    w2 = w

    k1, k2, k2_signbit = balanced_decomposition(C, coeff)

    if iszero(k1)
        return apply_signbit.(endomorphism.(points) .* k2, k2_signbit)
    elseif iszero(k2)
        return points .* k1
    end

    points2 = Array{P}(undef, length(points))
    for j in 1:length(points)
        points2[j] = apply_signbit(curve_endomorphism_type_4(points[j]), k2_signbit)
    end


    # corresponds to di = [1, 3, 5, ..., 2^(w-1)-1, -2^(w-1)-1, ..., -3, -1]
    l1 = 1 << (w1 - 1)
    precomp1 = Array{P}(undef, length(points), l1)
    l2 = 1 << (w2 - 1)
    precomp2 = Array{P}(undef, length(points), l2)

    dpoints = Array{P}(undef, length(points))
    for j in 1:length(points)
        dpoints[j] = double(points[j])
        precomp1[j,1] = points[j]
        precomp1[j,end] = -points[j]
    end

    for i in 2:(l1>>1)
        for j in 1:length(points)
            precomp1[j,i] = precomp1[j,i-1] + dpoints[j]
            precomp1[j,end-i+1] = -precomp1[j,i]
        end
    end

    dpoints2 = Array{P}(undef, length(points))
    for j in 1:length(points)
        dpoints2[j] = double(points2[j])
        precomp2[j,1] = points2[j]
        precomp2[j,end] = -points2[j]
    end

    for i in 2:(l2>>1)
        for j in 1:length(points)
            precomp2[j,i] = precomp2[j,i-1] + dpoints2[j]
            precomp2[j,end-i+1] = -precomp2[j,i]
        end
    end

    ds1 = get_wnaf(k1, w1)
    ds2 = get_wnaf(k2, w2)

    if length(ds1) > length(ds2)
        ds2 = vcat(ds2, zeros(Int, length(ds1) - length(ds2)))
    elseif length(ds2) > length(ds1)
        ds1 = vcat(ds1, zeros(Int, length(ds2) - length(ds1)))
    end

    acc = zeros(P, length(points))
    for idx in length(ds1):-1:1
        for j in 1:length(points)
            acc[j] = double(acc[j])
        end
        if !iszero(ds1[idx])
            for j in 1:length(points)
                acc[j] += precomp1[j,(ds1[idx] >> 1) + 1]
            end
        end
        if !iszero(ds2[idx])
            for j in 1:length(points)
                acc[j] += precomp2[j,(ds2[idx] >> 1) + 1]
            end
        end
    end

    acc
end


"""
Returns `points .* coeff`.
"""
function batch_mul(
        points::Array{P, 1}, coeff::T, w::Int=4,
        ) where {P <: EllipticCurvePoint{C}, T <: Integer} where {C <: EndomorphismType4Curve}
    batch_mul_endomorphism_wnaf(points, coeff, w)
end


function batch_mul(
        points::Array{P, 1}, coeff::T, w::Int=4,
        ) where {P <: EllipticCurvePoint, T <: Integer} where C
    batch_mul_wnaf(points, coeff, w)
end
