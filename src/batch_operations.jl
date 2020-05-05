"""
Linear combination using Shamir's trick and a fixed window.
"""
function lin_comb_windowed(
        points::Array{P, 1}, coeffs::Array{T, 1}, w::Int=4) where {P <: EllipticCurvePoint, T <: Integer}

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
            if !iszero(d)
                acc += precomp[j, d]
            end
            ys[j] -= d << (i * w)
        end
    end
    acc
end


"""
Calculates a linear combination of curve points
given an array of points and an array of coefficients.
"""
function lin_comb(
        points::Array{P, 1}, coeffs::Array{T, 1}, w::Int=4) where {P <: EllipticCurvePoint, T <: Integer}
    lin_comb_windowed(points, coeffs, w)
end


function batch_mul_wnaf(
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


function batch_mul_endomorhism_wnaf(
        points::Array{P, 1}, coeff::T, w1::Int=4, w2::Int=4,
        ) where {P <: EllipticCurvePoint{C, V}, T <: Integer} where {C, V}

    k1, k2, k2_signbit = balanced_decomposition(C, coeff)

    if iszero(k1)
        return apply_signbit.(endomorphism.(points) .* k2, k2_signbit)
    elseif iszero(k2)
        return points .* k1
    end

    points2 = Array{P}(undef, length(points))
    for j in 1:length(points)
        points2[j] = apply_signbit(endomorphism(points[j]), k2_signbit)
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
        ) where {P <: EllipticCurvePoint{C, V}, T <: Integer} where {C <: EndomorphismType4, V}
    batch_mul_endomorphism_wnaf(points, coeff, w, w)
end


function batch_mul(
        points::Array{P, 1}, coeff::T, w::Int=4,
        ) where {P <: EllipticCurvePoint{C, V}, T <: Integer} where {C, V}
    batch_mul_wnaf(points, coeff, w)
end
