function integral(A::AbstractInterpolation, t::Number)
    ((t < A.t[1] || t > A.t[end]) && !A.extrapolate) && throw(ExtrapolationError())
    integral(A, A.t[1], t)
end

function integral(A::AbstractInterpolation, t1::Number, t2::Number)
    ((t1 < A.t[1] || t1 > A.t[end]) && !A.extrapolate) && throw(ExtrapolationError())
    ((t2 < A.t[1] || t2 > A.t[end]) && !A.extrapolate) && throw(ExtrapolationError())
    !hasfield(typeof(A), :I) && throw(IntegralNotFoundError())
    (t2 < t1) && return -integral(A, t2, t1)
    # the index less than or equal to t1
    idx1 = get_idx(A, t1, 0)
    # the index less than t2
    idx2 = get_idx(A, t2, 0; idx_shift = -1, side = :first)

    if A.cache_parameters
        total = A.I[max(1, idx2 - 1)] - A.I[idx1]
        return if t1 == t2
            zero(total)
        else
            if idx1 == idx2
                total += _integral(A, idx1, t1, t2)
            else
                total += _integral(A, idx1, t1, A.t[idx1 + 1])
                total += _integral(A, idx2, A.t[idx2], t2)
            end
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

function _integral(A::SmoothedConstantInterpolation{<:AbstractVector},
        idx::Number, t1::Number, t2::Number)
    d_lower, d_upper, c_lower, c_upper = get_parameters(A, idx)

    bound_lower = A.t[idx] + d_lower
    bound_upper = A.t[idx + 1] - d_upper

    out = A.u[idx] * (t2 - t1)

    # Fix extrapolation behavior as constant for now
    if t1 <= first(A.t)
        t1 = first(A.t)
    elseif t2 >= last(A.t)
        t2 = last(A.t)
    end

    if t1 < bound_lower
        t2_ = min(t2, bound_lower)
        out -= c_lower * d_lower *
               (((t2_ - A.t[idx]) / d_lower - 1)^3 - ((t1 - A.t[idx]) / d_lower - 1)^3) / 3
    end

    if t2 > bound_upper
        t1_ = max(t1, bound_upper)
        out += c_upper * d_upper *
               ((1 - (A.t[idx + 1] - t2) / d_upper)^3 -
                (1 - (A.t[idx + 1] - t1_) / d_upper)^3) / 3
    end

    out
end

function _integral(A::QuadraticInterpolation{<:AbstractVector{<:Number}},
        idx::Number, t1::Number, t2::Number)
    α, β = get_parameters(A, idx)
    uᵢ = A.u[idx]
    tᵢ = A.t[idx]
    t1_rel = t1 - tᵢ
    t2_rel = t2 - tᵢ
    Δt = t2 - t1
    Δt * (α * (t2_rel^2 + t1_rel * t2_rel + t1_rel^2) / 3 + β * (t2_rel + t1_rel) / 2 + uᵢ)
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

function _integral(A::LagrangeInterpolation, idx::Number, t1::Number, t2::Number)
    throw(IntegralNotFoundError())
end
function _integral(A::BSplineInterpolation, idx::Number, t1::Number, t2::Number)
    throw(IntegralNotFoundError())
end
function _integral(A::BSplineApprox, idx::Number, t1::Number, t2::Number)
    throw(IntegralNotFoundError())
end

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
