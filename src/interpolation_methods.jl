function _interpolate(A, t)
    if t < first(A.t)
        _extrapolate_left(A, t)
    elseif t > last(A.t)
        _extrapolate_right(A, t)
    else
        _interpolate(A, t, A.iguesser)
    end
end

function _extrapolate_left(A, t)
    (; extrapolation_left) = A
    if extrapolation_left == ExtrapolationType.None
        throw(LeftExtrapolationError())
    elseif extrapolation_left == ExtrapolationType.Constant
        slope = _derivative(A, first(A.t), 1)
        first(A.u) + zero(slope * t)
    elseif extrapolation_left == ExtrapolationType.Linear
        slope = _derivative(A, first(A.t), 1)
        first(A.u) + slope * (t - first(A.t))
    else
        _extrapolate_other(A, t, extrapolation_left)
    end
end

function _extrapolate_right(A, t)
    (; extrapolation_right) = A
    if extrapolation_right == ExtrapolationType.None
        throw(RightExtrapolationError())
    elseif extrapolation_right == ExtrapolationType.Constant
        slope = _derivative(A, last(A.t), length(A.t))
        last(A.u) + zero(slope * t)
    elseif extrapolation_right == ExtrapolationType.Linear
        slope = _derivative(A, last(A.t), length(A.t))
        last(A.u) + slope * (t - last(A.t))
    else
        _extrapolate_other(A, t, extrapolation_right)
    end
end

function _extrapolate_other(A, t, extrapolation)
    if extrapolation == ExtrapolationType.Extension
        _interpolate(A, t, A.iguesser)
    elseif extrapolation == ExtrapolationType.Periodic
        t_, _ = transformation_periodic(A, t)
        _interpolate(A, t_, A.iguesser)
    elseif extrapolation == ExtrapolationType.Reflective
        t_, _ = transformation_reflective(A, t)
        _interpolate(A, t_, A.iguesser)
    else
        throw(ExtrapolationNotImplementedError())
    end
end

function _extrapolate_left(A::ConstantInterpolation, t)
    (; extrapolation_left) = A
    if extrapolation_left == ExtrapolationType.None
        throw(LeftExtrapolationError())
    elseif extrapolation_left in (ExtrapolationType.Constant, ExtrapolationType.Linear)
        first(A.u)
    else
        _extrapolate_other(A, t, extrapolation_left)
    end
end

function _extrapolate_right(A::ConstantInterpolation, t)
    (; extrapolation_right) = A
    if extrapolation_right == ExtrapolationType.None
        throw(RightExtrapolationError())
    elseif extrapolation_right in (ExtrapolationType.Constant, ExtrapolationType.Linear)
        last(A.u)
    else
        _extrapolate_other(A, t, extrapolation_right)
    end
end

function _extrapolate_right(A::SmoothedConstantInterpolation, t)
    if A.extrapolation_right == ExtrapolationType.None
        throw(RightExtrapolationError())
    elseif A.extrapolation_right in (
        ExtrapolationType.Constant, ExtrapolationType.Extension)
        d = min(A.t[end] - A.t[end - 1], 2A.d_max) / 2
        if A.t[end] + d < t
            A.u[end]
        else
            c = (A.u[end] - A.u[end - 1]) / 2
            A.u[end - 1] - c * (((t - A.t[end]) / d)^2 - 2 * ((t - A.t[end]) / d) - 1)
        end

    else
        _extrapolate_other(A, t, A.extrapolation_right)
    end
end

# Linear Interpolation
function _interpolate(A::LinearInterpolation{<:AbstractVector}, t::Number, iguess)
    if isnan(t)
        # For correct derivative with NaN
        idx = firstindex(A.u)
        t1 = t2 = oneunit(eltype(A.t))
        u1 = u2 = oneunit(eltype(A.u))
        slope = t / t * get_parameters(A, idx)
    else
        idx = get_idx(A, t, iguess)
        t1, t2 = A.t[idx], A.t[idx + 1]
        u1, u2 = A.u[idx], A.u[idx + 1]
        slope = get_parameters(A, idx)
    end

    Δt = t - t1
    Δu = slope * Δt
    val = u1
    Δu_nan = any(isnan.(Δu))
    if t == t2 && Δu_nan
        val = u2
    elseif !(iszero(Δt) && Δu_nan)
        val += Δu
    end
    val = oftype(Δu, val)

    val
end

function _interpolate(A::LinearInterpolation{<:AbstractArray}, t::Number, iguess)
    idx = get_idx(A, t, iguess)
    Δt = t - A.t[idx]
    slope = get_parameters(A, idx)
    ax = axes(A.u)[1:(end - 1)]
    return A.u[ax..., idx] + slope * Δt
end

# Quadratic Interpolation
function _interpolate(A::QuadraticInterpolation, t::Number, iguess)
    idx = get_idx(A, t, iguess)
    Δt = t - A.t[idx]
    α, β = get_parameters(A, idx)
    out = A.u isa AbstractMatrix ? A.u[:, idx] : A.u[idx]
    out += @. Δt * (α * Δt + β)
    out
end

# Lagrange Interpolation
function _interpolate(A::LagrangeInterpolation{<:AbstractVector}, t::Number, iguess)
    idx = get_idx(A, t, iguess)
    findRequiredIdxs!(A, t, idx)
    if A.t[A.idxs[1]] == t
        return A.u[A.idxs[1]]
    end
    N = zero(A.u[1])
    D = zero(A.t[1])
    tmp = N
    for i in 1:length(A.idxs)
        if isnan(A.bcache[A.idxs[i]])
            mult = one(A.t[1])
            for j in 1:(i - 1)
                mult *= (A.t[A.idxs[i]] - A.t[A.idxs[j]])
            end
            for j in (i + 1):length(A.idxs)
                mult *= (A.t[A.idxs[i]] - A.t[A.idxs[j]])
            end
            A.bcache[A.idxs[i]] = mult
        else
            mult = A.bcache[A.idxs[i]]
        end
        tmp = inv((t - A.t[A.idxs[i]]) * mult)
        D += tmp
        N += (tmp * A.u[A.idxs[i]])
    end
    N / D
end

function _interpolate(A::LagrangeInterpolation{<:AbstractMatrix}, t::Number, iguess)
    idx = get_idx(A, t, iguess)
    findRequiredIdxs!(A, t, idx)
    if A.t[A.idxs[1]] == t
        return A.u[:, A.idxs[1]]
    end
    N = zero(A.u[:, 1])
    D = zero(A.t[1])
    tmp = D
    for i in 1:length(A.idxs)
        if isnan(A.bcache[A.idxs[i]])
            mult = one(A.t[1])
            for j in 1:(i - 1)
                mult *= (A.t[A.idxs[i]] - A.t[A.idxs[j]])
            end
            for j in (i + 1):length(A.idxs)
                mult *= (A.t[A.idxs[i]] - A.t[A.idxs[j]])
            end
            A.bcache[A.idxs[i]] = mult
        else
            mult = A.bcache[A.idxs[i]]
        end
        tmp = inv((t - A.t[A.idxs[i]]) * mult)
        D += tmp
        @. N += (tmp * A.u[:, A.idxs[i]])
    end
    N / D
end

function _interpolate(A::AkimaInterpolation{<:AbstractVector}, t::Number, iguess)
    idx = get_idx(A, t, iguess)
    wj = t - A.t[idx]
    @evalpoly wj A.u[idx] A.b[idx] A.c[idx] A.d[idx]
end

# Constant Interpolation
function _interpolate(A::ConstantInterpolation{<:AbstractVector}, t::Number, iguess)
    if A.dir === :left
        # :left means that value to the left is used for interpolation
        idx = get_idx(A, t, iguess; lb = 1, ub_shift = 0)
    else
        # :right means that value to the right is used for interpolation
        idx = get_idx(A, t, iguess; side = :first, lb = 1, ub_shift = 0)
    end
    A.u[idx]
end

function _interpolate(A::ConstantInterpolation{<:AbstractMatrix}, t::Number, iguess)
    if A.dir === :left
        # :left means that value to the left is used for interpolation
        idx = get_idx(A, t, iguess; lb = 1, ub_shift = 0)
    else
        # :right means that value to the right is used for interpolation
        idx = get_idx(A, t, iguess; side = :first, lb = 1, ub_shift = 0)
    end
    A.u[:, idx]
end

# Smoothed constant Interpolation
function _interpolate(A::SmoothedConstantInterpolation{<:AbstractVector}, t::Number, iguess)
    idx = get_idx(A, t, iguess)
    d_lower, d_upper, c_lower, c_upper = get_parameters(A, idx)

    out = A.u[idx]

    if (t - A.t[idx]) < d_lower
        out -= c_lower * ((t - A.t[idx]) / d_lower - 1)^2
    elseif (A.t[idx + 1] - t) < d_upper
        out += c_upper * (1 - (A.t[idx + 1] - t) / d_upper)^2
    end

    out
end

# QuadraticSpline Interpolation
function _interpolate(A::QuadraticSpline{<:AbstractVector}, t::Number, iguess)
    idx = get_idx(A, t, iguess)
    α, β = get_parameters(A, idx)
    uᵢ = A.u[idx]
    Δt_scaled = (t - A.t[idx]) / (A.t[idx + 1] - A.t[idx])
    Δt_scaled * (α * Δt_scaled + β) + uᵢ
end

# CubicSpline Interpolation
function _interpolate(A::CubicSpline{<:AbstractVector}, t::Number, iguess)
    idx = get_idx(A, t, iguess)
    Δt₁ = t - A.t[idx]
    Δt₂ = A.t[idx + 1] - t
    I = (A.z[idx] * Δt₂^3 + A.z[idx + 1] * Δt₁^3) / (6A.h[idx + 1])
    c₁, c₂ = get_parameters(A, idx)
    C = c₁ * Δt₁
    D = c₂ * Δt₂
    I + C + D
end

function _interpolate(A::CubicSpline{<:AbstractArray}, t::Number, iguess)
    idx = get_idx(A, t, iguess)
    Δt₁ = t - A.t[idx]
    Δt₂ = A.t[idx + 1] - t
    ax = axes(A.z)[1:(end - 1)]
    I = (A.z[ax..., idx] * Δt₂^3 + A.z[ax..., idx + 1] * Δt₁^3) / (6A.h[idx + 1])
    c₁, c₂ = get_parameters(A, idx)
    C = c₁ * Δt₁
    D = c₂ * Δt₂
    I + C + D
end

# BSpline Curve Interpolation
function _interpolate(A::BSplineInterpolation{<:AbstractVector{<:Number}},
        t::Number,
        iguess)
    t < A.t[1] && return A.u[1]
    t > A.t[end] && return A.u[end]
    # change t into param [0 1]
    idx = get_idx(A, t, iguess)
    t = A.p[idx] + (t - A.t[idx]) / (A.t[idx + 1] - A.t[idx]) * (A.p[idx + 1] - A.p[idx])
    n = length(A.t)
    sc = t isa ForwardDiff.Dual ? zeros(eltype(t), n) : A.sc
    nonzero_coefficient_idxs = spline_coefficients!(sc, A.d, A.k, t)
    ucum = zero(eltype(A.u))
    for i in nonzero_coefficient_idxs
        ucum += sc[i] * A.c[i]
    end
    ucum
end

function _interpolate(A::BSplineInterpolation{<:AbstractArray{<:Number}},
        t::Number,
        iguess)
    ax_u = axes(A.u)[1:(end - 1)]
    t < A.t[1] && return A.u[ax_u..., 1]
    t > A.t[end] && return A.u[ax_u..., end]
    # change t into param [0 1]
    idx = get_idx(A, t, iguess)
    t = A.p[idx] + (t - A.t[idx]) / (A.t[idx + 1] - A.t[idx]) * (A.p[idx + 1] - A.p[idx])
    n = length(A.t)
    sc = t isa ForwardDiff.Dual ? zeros(eltype(t), n) : A.sc
    nonzero_coefficient_idxs = spline_coefficients!(sc, A.d, A.k, t)
    ucum = zeros(eltype(A.u), size(A.u)[1:(end - 1)]...)
    for i in nonzero_coefficient_idxs
        ucum = ucum + (sc[i] * A.c[ax_u..., i])
    end
    ucum
end

# BSpline Curve Approx
function _interpolate(A::BSplineApprox{<:AbstractVector{<:Number}}, t::Number, iguess)
    t < A.t[1] && return A.u[1]
    t > A.t[end] && return A.u[end]
    # change t into param [0 1]
    idx = get_idx(A, t, iguess)
    t = A.p[idx] + (t - A.t[idx]) / (A.t[idx + 1] - A.t[idx]) * (A.p[idx + 1] - A.p[idx])
    sc = t isa ForwardDiff.Dual ? zeros(eltype(t), A.h) : A.sc
    nonzero_coefficient_idxs = spline_coefficients!(sc, A.d, A.k, t)
    ucum = zero(eltype(A.u))
    for i in nonzero_coefficient_idxs
        ucum += sc[i] * A.c[i]
    end
    ucum
end

function _interpolate(
        A::BSplineApprox{<:AbstractArray{<:Number}}, t::Number, iguess)
    ax_u = axes(A.u)[1:(end - 1)]
    t < A.t[1] && return A.u[ax_u..., 1]
    t > A.t[end] && return A.u[ax_u..., end]
    # change t into param [0 1]
    idx = get_idx(A, t, iguess)
    t = A.p[idx] + (t - A.t[idx]) / (A.t[idx + 1] - A.t[idx]) * (A.p[idx + 1] - A.p[idx])
    sc = t isa ForwardDiff.Dual ? zeros(eltype(t), A.h) : A.sc
    nonzero_coefficient_idxs = spline_coefficients!(sc, A.d, A.k, t)
    ucum = zeros(eltype(A.u), size(A.u)[1:(end - 1)]...)
    for i in nonzero_coefficient_idxs
        ucum = ucum + (sc[i] * A.c[ax_u..., i])
    end
    ucum
end

# Cubic Hermite Spline
function _interpolate(
        A::CubicHermiteSpline{<:AbstractVector}, t::Number, iguess)
    idx = get_idx(A, t, iguess)
    Δt₀ = t - A.t[idx]
    Δt₁ = t - A.t[idx + 1]
    out = A.u[idx] + Δt₀ * A.du[idx]
    c₁, c₂ = get_parameters(A, idx)
    out += Δt₀^2 * (c₁ + Δt₁ * c₂)
    out
end

# Quintic Hermite Spline
function _interpolate(
        A::QuinticHermiteSpline{<:AbstractVector{<:Number}}, t::Number, iguess)
    idx = get_idx(A, t, iguess)
    Δt₀ = t - A.t[idx]
    Δt₁ = t - A.t[idx + 1]
    out = A.u[idx] + Δt₀ * (A.du[idx] + A.ddu[idx] * Δt₀ / 2)
    c₁, c₂, c₃ = get_parameters(A, idx)
    out += Δt₀^3 * (c₁ + Δt₁ * (c₂ + c₃ * Δt₁))
    out
end

function _interpolate(A::SmoothArcLengthInterpolation, t::Number, iguess)
    (; out, in_place) = A
    out = in_place ? out : similar(out, typeof(t))
    idx = get_idx(A, t, iguess)
    Δt_circ_seg = A.Δt_circle_segment[idx]
    Δt_line_seg = A.Δt_line_segment[idx]
    short_side_left = A.short_side_left[idx]
    Δt = t - A.t[idx]

    in_circle_arc = if short_side_left
        Δt < Δt_circ_seg
    else
        Δt > Δt_line_seg
    end

    if in_circle_arc
        t_circle_seg = short_side_left ? Δt : Δt - Δt_line_seg
        Rⱼ = A.radius[idx]
        S, C = sincos(t_circle_seg / Rⱼ)
        c = view(A.center, :, idx)
        v₁ = view(A.dir_1, :, idx)
        v₂ = view(A.dir_2, :, idx)
        @. out = c + C * v₁ + S * v₂
    else
        if short_side_left
            u₁ = view(A.u, :, idx + 1)
            d₁ = view(A.d, :, idx + 1)
            t_line_seg = A.t[idx + 1] - t
            @. out = u₁ - t_line_seg * d₁
        else
            u₀ = view(A.u, :, idx)
            d₀ = view(A.d, :, idx)
            @. out = u₀ + Δt * d₀
        end
    end

    out
end
