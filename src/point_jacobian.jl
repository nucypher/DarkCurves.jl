struct JacobianPoint{C <: EllipticCurve, T <: Number} <: StandardEllipticCurvePoint{C}
    x :: T
    y :: T
    z :: T
    inf :: Bool

    JacobianPoint{C, T}(x::T, y::T, z::T) where {C, T} = new{C, T}(x, y, z, false)
    JacobianPoint{C, T}(x::T, y::T) where {C, T} = new{C, T}(x, y, one(T), false)
    JacobianPoint{C, T}() where {C, T} = new{C, T}(one(T), one(T), zero(T), true)
end


function Base.:(==)(p::JacobianPoint{C, T}, q::JacobianPoint{C, T}) where {C, T}
    if iszero(p)
        return iszero(q)
    elseif iszero(q)
        return false
    end

    p_z_squared = p.z^2
    q_z_squared = q.z^2
    if p.x * q_z_squared != q.x * p_z_squared
        return false
    end

    p_z_cubed = p_z_squared * p.z
    q_z_cubed = q_z_squared * q.z
    p.y * q_z_cubed == q.y * p_z_cubed
end


Base.convert(::Type{JacobianPoint{C, T}}, p::AffinePoint{C, T}) where {C, T} =
    JacobianPoint{C, T}(p.x, p.y)


Base.zero(::Type{JacobianPoint{C, T}}) where {C, T} = JacobianPoint{C, T}()
Base.zero(::JacobianPoint{C, T}) where {C, T} = zero(JacobianPoint{C, T})


Base.iszero(p::JacobianPoint{C, T}) where {C, T} = p.inf


function Base.one(::Type{JacobianPoint{C, T}}) where {C, T}
    bx, by = curve_base(C, T)
    JacobianPoint{C, T}(bx, by)
end


function Base.:-(p::JacobianPoint{C, T}) where {C, T}
    if iszero(p)
        p
    else
        JacobianPoint{C, T}(p.x, -p.y, p.z)
    end
end


function double(p::JacobianPoint{C, T}) where {C, T}

    if iszero(p)
        return p
    end

    a = curve_weierstrass_coeff_a(C)

    XX = square(p.x)
    YY = square(p.y)
    YYYY = square(YY)
    ZZ = square(p.z)
    S = double(square(p.x + YY) - XX - YYYY)
    M = triple(XX)
    if !iszero(a)
        M += convert(T, a) * square(ZZ)
    end
    T_ = square(M) - double(S)
    X3 = T_
    Y3 = M * (S - T_) - double(double(double(YYYY)))
    Z3 = square(p.y + p.z) - YY - ZZ

    JacobianPoint{C, T}(X3, Y3, Z3)
end


function Base.:+(p::JacobianPoint{C, T}, q::JacobianPoint{C, T}) where {C, T}
    if iszero(p)
        return q
    end

    if iszero(q)
        return p
    end

    X1, Y1, Z1 = p.x, p.y, p.z
    X2, Y2, Z2 = q.x, q.y, q.z

    Z1Z1 = square(Z1)
    Z2Z2 = square(Z2)
    U1 = X1 * Z2Z2
    U2 = X2 * Z1Z1
    S1 = Y1 * Z2 * Z2Z2
    S2 = Y2 * Z1 * Z1Z1

    if U1 == U2
        if S1 != S2
            return JacobianPoint{C, T}()
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
    Z3 = (square(Z1 + Z2) - Z1Z1 - Z2Z2) * H

    if iszero(Z3)
        JacobianPoint{C, T}()
    else
        JacobianPoint{C, T}(X3, Y3, Z3)
    end
end


curve_endomorphism_type_4(p::JacobianPoint{C, T}) where {C <: EndomorphismType4Curve, T} =
    iszero(p) ? p : JacobianPoint{C, T}(curve_endomorphism_type_4_beta(C, T) * p.x, p.y, p.z)
