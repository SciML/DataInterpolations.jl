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
    Δt = t - A.t[idx]
    Δt * (A.u[idx] + A.p.slope[idx] * Δt / 2)
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
    t₀ = A.t[idx]
    t₁ = A.t[idx + 1]
    t₂ = A.t[idx + 2]

    t_sq = (t^2) / 3
    Iu₀ = A.p.l₀[idx] * t * (t_sq - t * (t₁ + t₂) / 2 + t₁ * t₂)
    Iu₁ = A.p.l₁[idx] * t * (t_sq - t * (t₀ + t₂) / 2 + t₀ * t₂)
    Iu₂ = A.p.l₂[idx] * t * (t_sq - t * (t₀ + t₁) / 2 + t₀ * t₁)
    return Iu₀ + Iu₁ + Iu₂
end

function _integral(A::QuadraticSpline{<:AbstractVector{<:Number}}, idx::Number, t::Number)
    Cᵢ = A.u[idx]
    Δt = t - A.t[idx]
    return A.z[idx] * Δt^2 / 2 + A.p.σ[idx] * Δt^3 / 3 + Cᵢ * Δt
end

function _integral(A::CubicSpline{<:AbstractVector{<:Number}}, idx::Number, t::Number)
    Δt₁sq = (t - A.t[idx])^2 / 2
    Δt₂sq = (A.t[idx + 1] - t)^2 / 2
    II = (-A.z[idx] * Δt₂sq^2 + A.z[idx + 1] * Δt₁sq^2) / (6A.h[idx + 1])
    IC = A.p.c₁[idx] * Δt₁sq
    ID = -A.p.c₂[idx] * Δt₂sq
    II + IC + ID
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

# Cubic Hermite Spline
function _integral(
        A::CubicHermiteSpline{<:AbstractVector{<:Number}}, idx::Number, t::Number)
    Δt₀ = t - A.t[idx]
    Δt₁ = t - A.t[idx + 1]
    out = Δt₀ * (A.u[idx] + Δt₀ * A.du[idx] / 2)
    p = A.p.c₁[idx] + Δt₁ * A.p.c₂[idx]
    dp = A.p.c₂[idx]
    out += Δt₀^3 / 3 * (p - dp * Δt₀ / 4)
    out
end

# Quintic Hermite Spline
function _integral(
        A::QuinticHermiteSpline{<:AbstractVector{<:Number}}, idx::Number, t::Number)
    Δt₀ = t - A.t[idx]
    Δt₁ = t - A.t[idx + 1]
    out = Δt₀ * (A.u[idx] + A.du[idx] * Δt₀ / 2 + A.ddu[idx] * Δt₀^2 / 6)
    p = A.p.c₁[idx] + A.p.c₂[idx] * Δt₁ + A.p.c₃[idx] * Δt₁^2
    dp = A.p.c₂[idx] + 2A.p.c₃[idx] * Δt₁
    ddp = 2A.p.c₃[idx]
    out += Δt₀^4 / 4 * (p - Δt₀ / 5 * dp + Δt₀^2 / 30 * ddp)
    out
end
