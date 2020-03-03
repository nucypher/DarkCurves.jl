struct AffinePoint{C <: EllipticCurve, T <: Number} <: EllipticCurvePoint{C, T}
    x :: T
    y :: T
    inf :: Bool

    AffinePoint{C, T}(x::T, y::T) where {C, T} = new{C, T}(x, y, false)
    AffinePoint{C, T}() where {C, T} = new{C, T}(zero(T), zero(T), true)
end


Base.zero(::Type{AffinePoint{C, T}}) where {C, T} = AffinePoint{C, T}()


Base.iszero(p::AffinePoint{C, T}) where {C, T} = p.inf


function Base.one(::Type{AffinePoint{C, T}}) where {C, T}
    bx, by = curve_base(C, T)
    AffinePoint{C, T}(bx, by)
end


function Base.:-(p::AffinePoint{C, T}) where {C, T}
    if iszero(p)
        p
    else
        AffinePoint{C, T}(p.x, -p.y)
    end
end


function Base.:+(p::AffinePoint{C, T}, q::AffinePoint{C, T}) where {C, T}
    if iszero(p)
        return q
    end

    if iszero(q)
        return p
    end

    if p.x == q.x && (p.y != q.y || p.y == 0)
        return zero(AffinePoint{C, T})
    elseif p.x != q.x
        l = (q.y - p.y) * inv(q.x - p.x)
    else
        t = p.x^2
        l = (triple(t) + curve_coeff_a(C, T)) * inv(p.y + p.y)
    end
    x = l^2 - p.x - q.x
    y = l * (p.x - x) - p.y
    AffinePoint{C, T}(x, y)
end
