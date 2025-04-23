function integral(A::AbstractInterpolation, t::Number)
    integral(A, first(A.t), t)
end

function integral(A::AbstractInterpolation, t1::Number, t2::Number)
    !hasfield(typeof(A), :I) && throw(IntegralNotFoundError())

    if t1 == t2
        # If the integration interval is trivial then the result is 0
        return zero(eltype(A.I))
    elseif t1 > t2
        # Make sure that t1 < t2
        return -integral(A, t2, t1)
    end

    # the index less than or equal to t1
    idx1 = get_idx(A, t1, 0)
    # the index less than t2
    idx2 = get_idx(A, t2, 0; idx_shift = -1, side = :first)

    total = zero(eltype(A.I))

    # Lower potentially incomplete interval
    if t1 < first(A.t)
        if t2 < first(A.t)
            # If interval is entirely below data
            return _extrapolate_integral_left(A, t1) - _extrapolate_integral_left(A, t2)
        end

        idx1 -= 1 # Make sure lowest complete interval is included
        total += _extrapolate_integral_left(A, t1)
    else
        total += _integral(A, idx1, t1, A.t[idx1 + 1])
    end

    # Upper potentially incomplete interval
    if t2 > last(A.t)
        if t1 > last(A.t)
            # If interval is entirely above data
            return _extrapolate_integral_right(A, t2) - _extrapolate_integral_right(A, t1)
        end

        idx2 += 1 # Make sure highest complete interval is included
        total += _extrapolate_integral_right(A, t2)
    else
        total += _integral(A, idx2, A.t[idx2], t2)
    end

    if idx1 == idx2
        return _integral(A, idx1, t1, t2)
    end

    # Complete intervals
    if A.cache_parameters
        if idx2 > 1
            total += A.I[idx2 - 1]
        end
        if idx1 > 0
            total -= A.I[idx1]
        end
    else
        for idx in (idx1 + 1):(idx2 - 1)
            total += _integral(A, idx, A.t[idx], A.t[idx + 1])
        end
    end

    return total
end

function _extrapolate_integral_left(A, t)
    (; extrapolation_left) = A
    if extrapolation_left == ExtrapolationType.None
        throw(LeftExtrapolationError())
    elseif extrapolation_left == ExtrapolationType.Constant
        first(A.u) * (first(A.t) - t)
    elseif extrapolation_left == ExtrapolationType.Linear
        slope = derivative(A, first(A.t))
        Δt = first(A.t) - t
        (first(A.u) - slope * Δt / 2) * Δt
    elseif extrapolation_left == ExtrapolationType.Extension
        _integral(A, 1, t, first(A.t))
    elseif extrapolation_left == ExtrapolationType.Periodic
        t_, n = transformation_periodic(A, t)
        out = -integral(A, t_)
        if !iszero(n)
            out -= n * integral(A, first(A.t), last(A.t))
        end
        out
    else
        # extrapolation_left == ExtrapolationType.Reflective
        t_, n = transformation_reflective(A, t)
        out = if isodd(n)
            -integral(A, t_, last(A.t))
        else
            -integral(A, t_)
        end
        if !iszero(n)
            out -= n * integral(A, first(A.t), last(A.t))
        end
        out
    end
end

function _extrapolate_integral_right(A, t)
    (; extrapolation_right) = A
    if extrapolation_right == ExtrapolationType.None
        throw(RightExtrapolationError())
    elseif extrapolation_right == ExtrapolationType.Constant
        last(A.u) * (t - last(A.t))
    elseif extrapolation_right == ExtrapolationType.Linear
        slope = derivative(A, last(A.t))
        Δt = t - last(A.t)
        (last(A.u) + slope * Δt / 2) * Δt
    elseif extrapolation_right == ExtrapolationType.Extension
        _integral(A, length(A.t) - 1, last(A.t), t)
    elseif extrapolation_right == ExtrapolationType.Periodic
        t_, n = transformation_periodic(A, t)
        out = integral(A, first(A.t), t_)
        if !iszero(n)
            out += n * integral(A, first(A.t), last(A.t))
        end
        out
    else
        # extrapolation_right == ExtrapolationType.Reflective
        t_, n = transformation_reflective(A, t)
        out = if iseven(n)
            integral(A, t_, last(A.t))
        else
            integral(A, t_)
        end
        if !iszero(n)
            out += n * integral(A, first(A.t), last(A.t))
        end
        out
    end
end

function _extrapolate_integral_right(A::SmoothedConstantInterpolation, t)
    (; extrapolation_right) = A
    if extrapolation_right == ExtrapolationType.None
        throw(RightExtrapolationError())
    elseif A.extrapolation_right in (
        ExtrapolationType.Constant, ExtrapolationType.Extension)
        d = min(A.t[end] - A.t[end - 1], 2A.d_max) / 2
        c = (A.u[end] - A.u[end - 1]) / 2

        Δt_transition = min(t - A.t[end], d)
        Δt_constant = max(0, t - A.t[end] - d)
        out = Δt_transition * A.u[end - 1] -
              c *
              (((Δt_transition / d)^3) / (3 / d) - ((Δt_transition^2) / d) -
               Δt_transition) +
              Δt_constant * A.u[end]

    elseif extrapolation_right == ExtrapolationType.Linear
        slope = derivative(A, last(A.t))
        Δt = t - last(A.t)
        (last(A.u) + slope * Δt / 2) * Δt
        _extrapolate_other(A, t, A.extrapolation_right)
    elseif extrapolation_right == ExtrapolationType.Periodic
        t_, n = transformation_periodic(A, t)
        out = integral(A, first(A.t), t_)
        if !iszero(n)
            out += n * integral(A, first(A.t), last(A.t))
        end
        out
    else
        # extrapolation_right == ExtrapolationType.Reflective
        t_, n = transformation_reflective(A, t)
        out = if iseven(n)
            integral(A, t_, last(A.t))
        else
            integral(A, t_)
        end
        if !iszero(n)
            out += n * integral(A, first(A.t), last(A.t))
        end
        out
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
