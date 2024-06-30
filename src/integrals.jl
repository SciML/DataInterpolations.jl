function integral(A::AbstractInterpolation, t::Number)
    ((t < A.t[1] || t > A.t[end]) && !A.extrapolate) && throw(ExtrapolationError())
    integral(A, A.t[1], t)
end

function integral(A::AbstractInterpolation, t1::Number, t2::Number)
    ((t1 < A.t[1] || t1 > A.t[end]) && !A.extrapolate) && throw(ExtrapolationError())
    ((t2 < A.t[1] || t2 > A.t[end]) && !A.extrapolate) && throw(ExtrapolationError())
    # the index less than or equal to t1
    idx1 = max(1, min(searchsortedlast(A.t, t1), length(A.t) - 1))
    # the index less than t2
    idx2 = max(1, min(searchsortedlast(A.t, t2), length(A.t) - 1))
    if A.t[idx2] == t2
        idx2 -= 1
    end
    total = zero(eltype(A.u))
    for idx in idx1:idx2
        lt1 = idx == idx1 ? t1 : A.t[idx]
        lt2 = idx == idx2 ? t2 : A.t[idx + 1]
        total += _integral(A, idx, lt2) - _integral(A, idx, lt1)
    end
    total
end

function _integral(A::LinearInterpolation{<:AbstractVector{<:Number}},
        idx::Number,
        t::Number)
    t1 = A.t[idx]
    t2 = A.t[idx + 1]
    u1 = A.u[idx]
    u2 = A.u[idx + 1]
    t^2 * (u1 - u2) / (2 * t1 - 2 * t2) + t * (t1 * u2 - t2 * u1) / (t1 - t2)
end

function _integral(A::ConstantInterpolation{<:AbstractVector}, idx::Number, t::Number)
    if A.dir === :left
        # :left means that value to the left is used for interpolation
        return A.u[idx] * t
    else
        # :right means that value to the right is used for interpolation
        return A.u[idx + 1] * t
    end
end

function _integral(A::QuadraticInterpolation{<:AbstractVector{<:Number}},
        idx::Number,
        t::Number)
    A.mode == :Backward && idx > 1 && (idx -= 1)
    idx = min(length(A.t) - 2, idx)
    t1 = A.t[idx]
    t2 = A.t[idx + 1]
    t3 = A.t[idx + 2]
    u1 = A.u[idx]
    u2 = A.u[idx + 1]
    u3 = A.u[idx + 2]
    (t^3 * (-t1 * u2 + t1 * u3 + t2 * u1 - t2 * u3 - t3 * u1 + t3 * u2) /
     (3 * t1^2 * t2 - 3 * t1^2 * t3 - 3 * t1 * t2^2 + 3 * t1 * t3^2 + 3 * t2^2 * t3 -
      3 * t2 * t3^2) +
     t^2 * (t1^2 * u2 - t1^2 * u3 - t2^2 * u1 + t2^2 * u3 + t3^2 * u1 - t3^2 * u2) /
     (2 * t1^2 * t2 - 2 * t1^2 * t3 - 2 * t1 * t2^2 + 2 * t1 * t3^2 + 2 * t2^2 * t3 -
      2 * t2 * t3^2) +
     t *
     (t1^2 * t2 * u3 - t1^2 * t3 * u2 - t1 * t2^2 * u3 + t1 * t3^2 * u2 + t2^2 * t3 * u1 -
      t2 * t3^2 * u1) /
     (t1^2 * t2 - t1^2 * t3 - t1 * t2^2 + t1 * t3^2 + t2^2 * t3 - t2 * t3^2))
end

function _integral(A::QuadraticSpline{<:AbstractVector{<:Number}}, idx::Number, t::Number)
    t1 = A.t[idx]
    t2 = A.t[idx + 1]
    u1 = A.u[idx]
    z1 = A.z[idx]
    z2 = A.z[idx + 1]
    t^3 * (z1 - z2) / (6 * t1 - 6 * t2) + t^2 * (t1 * z2 - t2 * z1) / (2 * t1 - 2 * t2) +
    t * (-t1^2 * z1 - t1^2 * z2 + 2 * t1 * t2 * z1 + 2 * t1 * u1 - 2 * t2 * u1) /
    (2 * t1 - 2 * t2)
end

function _integral(A::CubicSpline{<:AbstractVector{<:Number}}, idx::Number, t::Number)
    t1 = A.t[idx]
    t2 = A.t[idx + 1]
    u1 = A.u[idx]
    u2 = A.u[idx + 1]
    z1 = A.z[idx]
    z2 = A.z[idx + 1]
    h2 = A.h[idx + 1]
    (t^4 * (-z1 + z2) / (24 * h2) + t^3 * (-t1 * z2 + t2 * z1) / (6 * h2) +
     t^2 * (h2^2 * z1 - h2^2 * z2 + 3 * t1^2 * z2 - 3 * t2^2 * z1 - 6 * u1 + 6 * u2) /
     (12 * h2) +
     t *
     (h2^2 * t1 * z2 - h2^2 * t2 * z1 - t1^3 * z2 - 6 * t1 * u2 + t2^3 * z1 + 6 * t2 * u1) /
     (6 * h2))
end

function _integral(A::AkimaInterpolation{<:AbstractVector{<:Number}},
        idx::Number,
        t::Number)
    t1 = A.t[idx]
    A.u[idx] * (t - t1) + A.b[idx] * ((t - t1)^2 / 2) + A.c[idx] * ((t - t1)^3 / 3) +
    A.d[idx] * ((t - t1)^4 / 4)
end

integral(A::LagrangeInterpolation, t1::Number, t2::Number) = throw(IntegralNotFoundError())
integral(A::LagrangeInterpolation, t::Number) = throw(IntegralNotFoundError())

function integral(A::BSplineInterpolation{<:AbstractVector{<:Number}},
        t1::Number,
        t2::Number)
    throw(IntegralNotFoundError())
end
function integral(A::BSplineInterpolation{<:AbstractVector{<:Number}}, t::Number)
    throw(IntegralNotFoundError())
end

function integral(A::BSplineApprox{<:AbstractVector{<:Number}}, t1::Number, t2::Number)
    throw(IntegralNotFoundError())
end
function integral(A::BSplineApprox{<:AbstractVector{<:Number}}, t::Number)
    throw(IntegralNotFoundError())
end

# Quintic Hermite Interpolation
function _integral(
        A::QuinticHermiteInterpolation{<:AbstractVector{<:Number}}, idx::Number, t::Number)
    Δt₀ = t - A.t[idx]
    Δt₁ = t - A.t[idx + 1]
    out = Δt₀ * (A.u[idx] + A.du[idx] * Δt₀ / 2 + A.ddu[idx] * Δt₀^2 / 6)
    p = A.p.c₁[idx] + A.p.c₂[idx] * Δt₁ + A.p.c₃[idx] * Δt₁^2
    dp = A.p.c₂[idx] + 2A.p.c₃[idx] * Δt₁
    ddp = 2A.p.c₃[idx]
    out += Δt₀^4 / 4 * (p - Δt₀ / 5 * dp + Δt₀^2 / 30 * ddp)
    out
end
