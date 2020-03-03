"""
Linear combination using Shamir's trick and a fixed window.
"""
function lincomb_windowed(
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
function lincomb(
        points::Array{P, 1}, coeffs::Array{T, 1}, w::Int=4) where {P <: EllipticCurvePoint, T <: Integer}
    lincomb_windowed(points, coeffs, w)
end
