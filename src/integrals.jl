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
            total += _integral(A, idx1, t1, A.t[idx1 + 1])
            total += _integral(A, idx2, A.t[idx2], t2)
            total
        end
    else
        total = zero(eltype(A.u))
        for idx in idx1:idx2
            lt1 = idx == idx1 ? t1 : A.t[idx]
            lt2 = idx == idx2 ? t2 : A.t[idx + 1]
            total += _integral(A, idx, lt1, lt2)
        end
        total
    end
end

function _integral(A::LinearInterpolation{<:AbstractVector{<:Number}},
        idx::Number, t1::Number, t2::Number)
    slope = get_parameters(A, idx)
    u_mean = A.u[idx] + slope * ((t1 + t2) / 2 - A.t[idx])
    u_mean * (t2 - t1)
end

function _integral(
        A::ConstantInterpolation{<:AbstractVector{<:Number}}, idx::Number, t1::Number, t2::Number)
    Δt = t2 - t1
    if A.dir === :left
        # :left means that value to the left is used for interpolation
        return A.u[idx] * Δt
    else
        # :right means that value to the right is used for interpolation
        return A.u[idx + 1] * Δt
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

function _integral(
        A::QuadraticSpline{<:AbstractVector{<:Number}}, idx::Number, t1::Number, t2::Number)
    α, β = get_parameters(A, idx)
    uᵢ = A.u[idx]
    tᵢ = A.t[idx]
    t1_rel = t1 - tᵢ
    t2_rel = t2 - tᵢ
    Δt = t2 - t1
    Δt * (α * (t2_rel^2 + t1_rel * t2_rel + t1_rel^2) / 3 + β * (t2_rel + t1_rel) / 2 + uᵢ)
end

function _integral(
        A::CubicSpline{<:AbstractVector{<:Number}}, idx::Number, t1::Number, t2::Number)
    tᵢ = A.t[idx]
    tᵢ₊₁ = A.t[idx + 1]
    c₁, c₂ = get_parameters(A, idx)
    integrate_cubic_polynomial(t1, t2, tᵢ, 0, c₁, 0, A.z[idx + 1] / (6A.h[idx + 1])) +
    integrate_cubic_polynomial(t1, t2, tᵢ₊₁, 0, -c₂, 0, -A.z[idx] / (6A.h[idx + 1]))
end

function _integral(A::AkimaInterpolation{<:AbstractVector{<:Number}},
        idx::Number, t1::Number, t2::Number)
    integrate_cubic_polynomial(t1, t2, A.t[idx], A.u[idx], A.b[idx], A.c[idx], A.d[idx])
end

_integral(A::LagrangeInterpolation, idx::Number, t::Number) = throw(IntegralNotFoundError())
_integral(A::BSplineInterpolation, idx::Number, t::Number) = throw(IntegralNotFoundError())
_integral(A::BSplineApprox, idx::Number, t::Number) = throw(IntegralNotFoundError())

# Cubic Hermite Spline
function _integral(
        A::CubicHermiteSpline{<:AbstractVector{<:Number}}, idx::Number, t1::Number, t2::Number)
    c₁, c₂ = get_parameters(A, idx)
    tᵢ = A.t[idx]
    tᵢ₊₁ = A.t[idx + 1]
    c = c₁ - c₂ * (tᵢ₊₁ - tᵢ)
    integrate_cubic_polynomial(t1, t2, tᵢ, A.u[idx], A.du[idx], c, c₂)
end

# Quintic Hermite Spline
function _integral(
        A::QuinticHermiteSpline{<:AbstractVector{<:Number}}, idx::Number, t1::Number, t2::Number)
    tᵢ = A.t[idx]
    tᵢ₊₁ = A.t[idx + 1]
    Δt = tᵢ₊₁ - tᵢ
    c₁, c₂, c₃ = get_parameters(A, idx)
    integrate_quintic_polynomial(t1, t2, tᵢ, A.u[idx], A.du[idx], A.ddu[idx] / 2,
        c₁ + Δt * (-c₂ + c₃ * Δt), c₂ - 2c₃ * Δt, c₃)
end
