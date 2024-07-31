function integral(A::AbstractInterpolation, t::Number)
    ((t < A.t[1] || t > A.t[end]) && !A.extrapolate) && throw(ExtrapolationError())
    integral(A, A.t[1], t)
end

function integral(A::AbstractInterpolation, t1::Number, t2::Number)
    ((t1 < A.t[1] || t1 > A.t[end]) && !A.extrapolate) && throw(ExtrapolationError())
    ((t2 < A.t[1] || t2 > A.t[end]) && !A.extrapolate) && throw(ExtrapolationError())
    !hasfield(typeof(A), :I) && throw(IntegralNotFoundError())
    # the index less than or equal to t1
    idx1 = get_idx(A, t1, 0)
    # the index less than t2
    idx2 = get_idx(A, t2, 0; idx_shift = -1, side = :first)

    if A.cache_parameters
        total = A.I[idx2] - A.I[idx1]
        return if t1 == t2
            zero(total)
        else
            total += _integral(A, idx1, A.t[idx1])
            total -= _integral(A, idx1, t1)
            total += _integral(A, idx2, t2)
            total -= _integral(A, idx2, A.t[idx2])
            total
        end
    else
        total = zero(eltype(A.u))
        for idx in idx1:idx2
            lt1 = idx == idx1 ? t1 : A.t[idx]
            lt2 = idx == idx2 ? t2 : A.t[idx + 1]
            total += _integral(A, idx, lt2) - _integral(A, idx, lt1)
        end
        total
    end
end

function _integral(A::LinearInterpolation{<:AbstractVector{<:Number}},
        idx::Number,
        t::Number)
    Δt = t - A.t[idx]
    slope = get_parameters(A, idx)
    Δt * (A.u[idx] + slope * Δt / 2)
end

function _integral(
        A::ConstantInterpolation{<:AbstractVector{<:Number}}, idx::Number, t::Number)
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
    l₀, l₁, l₂ = get_parameters(A, idx)
    Iu₀ = l₀ * t * (t_sq - t * (t₁ + t₂) / 2 + t₁ * t₂)
    Iu₁ = l₁ * t * (t_sq - t * (t₀ + t₂) / 2 + t₀ * t₂)
    Iu₂ = l₂ * t * (t_sq - t * (t₀ + t₁) / 2 + t₀ * t₁)
    return Iu₀ + Iu₁ + Iu₂
end

function _integral(A::QuadraticSpline{<:AbstractVector{<:Number}}, idx::Number, t::Number)
    Cᵢ = A.u[idx]
    Δt = t - A.t[idx]
    σ = get_parameters(A, idx)
    return A.z[idx] * Δt^2 / 2 + σ * Δt^3 / 3 + Cᵢ * Δt
end

function _integral(A::CubicSpline{<:AbstractVector{<:Number}}, idx::Number, t::Number)
    Δt₁sq = (t - A.t[idx])^2 / 2
    Δt₂sq = (A.t[idx + 1] - t)^2 / 2
    II = (-A.z[idx] * Δt₂sq^2 + A.z[idx + 1] * Δt₁sq^2) / (6A.h[idx + 1])
    c₁, c₂ = get_parameters(A, idx)
    IC = c₁ * Δt₁sq
    ID = -c₂ * Δt₂sq
    II + IC + ID
end

function _integral(A::AkimaInterpolation{<:AbstractVector{<:Number}},
        idx::Number,
        t::Number)
    t1 = A.t[idx]
    A.u[idx] * (t - t1) + A.b[idx] * ((t - t1)^2 / 2) + A.c[idx] * ((t - t1)^3 / 3) +
    A.d[idx] * ((t - t1)^4 / 4)
end

_integral(A::LagrangeInterpolation, idx::Number, t::Number) = throw(IntegralNotFoundError())
_integral(A::BSplineInterpolation, idx::Number, t::Number) = throw(IntegralNotFoundError())
_integral(A::BSplineApprox, idx::Number, t::Number) = throw(IntegralNotFoundError())

# Cubic Hermite Spline
function _integral(
        A::CubicHermiteSpline{<:AbstractVector{<:Number}}, idx::Number, t::Number)
    Δt₀ = t - A.t[idx]
    Δt₁ = t - A.t[idx + 1]
    out = Δt₀ * (A.u[idx] + Δt₀ * A.du[idx] / 2)
    c₁, c₂ = get_parameters(A, idx)
    p = c₁ + Δt₁ * c₂
    dp = c₂
    out += Δt₀^3 / 3 * (p - dp * Δt₀ / 4)
    out
end

# Quintic Hermite Spline
function _integral(
        A::QuinticHermiteSpline{<:AbstractVector{<:Number}}, idx::Number, t::Number)
    Δt₀ = t - A.t[idx]
    Δt₁ = t - A.t[idx + 1]
    out = Δt₀ * (A.u[idx] + A.du[idx] * Δt₀ / 2 + A.ddu[idx] * Δt₀^2 / 6)
    c₁, c₂, c₃ = get_parameters(A, idx)
    p = c₁ + c₂ * Δt₁ + c₃ * Δt₁^2
    dp = c₂ + 2c₃ * Δt₁
    ddp = 2c₃
    out += Δt₀^4 / 4 * (p - Δt₀ / 5 * dp + Δt₀^2 / 30 * ddp)
    out
end
