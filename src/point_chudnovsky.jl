struct ChudnovskyPoint{C <: EllipticCurve, T <: Number} <: EllipticCurvePoint{C, T}
    x :: T
    y :: T
    z :: T
    z2 :: T
    z3 :: T
    inf :: Bool

    ChudnovskyPoint{C, T}(x::T, y::T, z::T, z2::T, z3::T) where {C, T} = new{C, T}(x, y, z, z2, z3, false)
    ChudnovskyPoint{C, T}(x::T, y::T) where {C, T} = new{C, T}(x, y, one(T), one(T), one(T), false)
    ChudnovskyPoint{C, T}() where {C, T} = new{C, T}(one(T), one(T), zero(T), zero(T), zero(T), true)
end


Base.convert(::Type{ChudnovskyPoint{C, T}}, p::AffinePoint{C, T}) where {C, T} =
    ChudnovskyPoint{C, T}(p.x, p.y)


function Base.:(==)(p::ChudnovskyPoint{C, T}, q::ChudnovskyPoint{C, T}) where {C, T}
    if iszero(p)
        return iszero(q)
    elseif iszero(q)
        return false
    end

    if p.x * q.z2 != q.x * p.z2
        return false
    end

    p.y * q.z3 == q.y * p.z3
end


Base.zero(::Type{ChudnovskyPoint{C, T}}) where {C, T} = ChudnovskyPoint{C, T}()


Base.iszero(p::ChudnovskyPoint{C, T}) where {C, T} = p.inf


function Base.one(::Type{ChudnovskyPoint{C, T}}) where {C, T}
    bx, by = curve_base(C, T)
    ChudnovskyPoint{C, T}(bx, by)
end


function Base.:-(p::ChudnovskyPoint{C, T}) where {C, T}
    if iszero(p)
        p
    else
        ChudnovskyPoint{C, T}(p.x, -p.y, p.z, p.z2, p.z3)
    end
end


function double(p::ChudnovskyPoint{C, T}) where {C, T}

    if iszero(p)
        return p
    end

    a = curve_coeff_a(C)

    S = double(double(p.x * square(p.y)))
    M = triple(square(p.x))
    if !iszero(a)
        M += convert(T, a) * square(p.z2)
    end
    Xp = square(M) - double(S)
    Yp = M * (S - Xp) - double(double(double(square(square(p.y)))))
    Zp = double(p.y * p.z)
    Zp2 = square(Zp)
    Zp3 = Zp2 * Zp

    ChudnovskyPoint{C, T}(Xp, Yp, Zp, Zp2, Zp3)
end


function Base.:+(p::ChudnovskyPoint{C, T}, q::ChudnovskyPoint{C, T}) where {C, T}
    if iszero(p)
        return q
    end

    if iszero(q)
        return p
    end

    X1, Y1, Z1 = p.x, p.y, p.z
    X2, Y2, Z2 = q.x, q.y, q.z

    U1 = X1 * q.z2
    U2 = X2 * p.z2
    S1 = Y1 * q.z3
    S2 = Y2 * p.z3
    if U1 == U2
        if S1 != S2
            return ChudnovskyPoint{C, T}()
        else
            return double(p)
        end
    end

    H = U2 - U1
    I = square(double(H))
    J = H * I
    r = double(S2 - S1)
    V = U1 * I
    X3 = square(r) - J - double(V)
    Y3 = r * (V - X3) - double(S1 * J)
    Z3 = (square(Z1 + Z2) - p.z2 - q.z2) * H

    Z3_2 = square(Z3)
    Z3_3 = Z3_2 * Z3
    ChudnovskyPoint{C, T}(X3, Y3, Z3, Z3_2, Z3_3)
end
