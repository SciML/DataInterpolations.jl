function derivative(A, t, order = 1)
    (order ∉ (1, 2)) && throw(DerivativeNotFoundError())
    if t < first(A.t)
        _extrapolate_derivative_left(A, t, order)
    elseif t > last(A.t)
        _extrapolate_derivative_right(A, t, order)
    else
        iguess = A.iguesser
        (order == 1) ? _derivative(A, t, iguess) :
        ForwardDiff.derivative(t -> begin
                _derivative(A, t, iguess)
            end, t)
    end
end

function _extrapolate_derivative_left(A, t, order)
    (; extrapolation_left) = A
    if extrapolation_left == ExtrapolationType.None
        throw(LeftExtrapolationError())
    elseif extrapolation_left == ExtrapolationType.Constant
        zero(first(A.u) / one(A.t[1]))
    elseif extrapolation_left == ExtrapolationType.Linear
        (order == 1) ? _derivative(A, first(A.t), 1) : zero(first(A.u) / one(A.t[1]))
    elseif extrapolation_left == ExtrapolationType.Extension
        (order == 1) ? _derivative(A, t, length(A.t)) :
        ForwardDiff.derivative(t -> begin
                _derivative(A, t, length(A.t))
            end, t)
    elseif extrapolation_left == ExtrapolationType.Periodic
        t_, _ = transformation_periodic(A, t)
        (order == 1) ? _derivative(A, t_, A.iguesser) :
        ForwardDiff.derivative(t -> begin
                _derivative(A, t, A.iguesser)
            end, t_)
    else
        # extrapolation_left == ExtrapolationType.Reflective
        t_, n = transformation_reflective(A, t)
        sign = isodd(n) ? -1 : 1
        (order == 1) ? sign * _derivative(A, t_, A.iguesser) :
        ForwardDiff.derivative(t -> begin
                sign * _derivative(A, t, A.iguesser)
            end, t_)
    end
end

function _extrapolate_derivative_right(A, t, order)
    (; extrapolation_right) = A
    if extrapolation_right == ExtrapolationType.None
        throw(RightExtrapolationError())
    elseif extrapolation_right == ExtrapolationType.Constant
        zero(first(A.u) / one(A.t[1]))
    elseif extrapolation_right == ExtrapolationType.Linear
        (order == 1) ? _derivative(A, last(A.t), length(A.t)) :
        zero(first(A.u) / one(A.t[1]))
    elseif extrapolation_right == ExtrapolationType.Extension
        (order == 1) ? _derivative(A, t, length(A.t)) :
        ForwardDiff.derivative(t -> begin
                _derivative(A, t, length(A.t))
            end, t)
    elseif extrapolation_right == ExtrapolationType.Periodic
        t_, _ = transformation_periodic(A, t)
        (order == 1) ? _derivative(A, t_, A.iguesser) :
        ForwardDiff.derivative(t -> begin
                _derivative(A, t, A.iguesser)
            end, t_)
    else
        # extrapolation_right == ExtrapolationType.Reflective
        t_, n = transformation_reflective(A, t)
        sign = iseven(n) ? -1 : 1
        (order == 1) ? sign * _derivative(A, t_, A.iguesser) :
        ForwardDiff.derivative(t -> begin
                sign * _derivative(A, t, A.iguesser)
            end, t_)
    end
end

function _derivative(A::LinearInterpolation, t::Number, iguess)
    idx = get_idx(A, t, iguess; idx_shift = -1, ub_shift = -1, side = :first)
    slope = get_parameters(A, idx)
    slope
end

function _derivative(A::QuadraticInterpolation, t::Number, iguess)
    idx = get_idx(A, t, iguess)
    Δt = t - A.t[idx]
    α, β = get_parameters(A, idx)
    return 2α * Δt + β
end

function _derivative(A::LagrangeInterpolation{<:AbstractVector}, t::Number)
    der = zero(A.u[1])
    for j in eachindex(A.t)
        tmp = zero(A.t[1])
        if isnan(A.bcache[j])
            mult = one(A.t[1])
            for i in 1:(j - 1)
                mult *= (A.t[j] - A.t[i])
            end
            for i in (j + 1):length(A.t)
                mult *= (A.t[j] - A.t[i])
            end
            A.bcache[j] = mult
        else
            mult = A.bcache[j]
        end
        for l in eachindex(A.t)
            if l != j
                k = one(A.t[1])
                for m in eachindex(A.t)
                    if m != j && m != l
                        k *= (t - A.t[m])
                    end
                end
                k *= inv(mult)
                tmp += k
            end
        end
        der += A.u[j] * tmp
    end
    der
end

function _derivative(A::LagrangeInterpolation{<:AbstractMatrix}, t::Number)
    der = zero(A.u[:, 1])
    for j in eachindex(A.t)
        tmp = zero(A.t[1])
        if isnan(A.bcache[j])
            mult = one(A.t[1])
            for i in 1:(j - 1)
                mult *= (A.t[j] - A.t[i])
            end
            for i in (j + 1):length(A.t)
                mult *= (A.t[j] - A.t[i])
            end
            A.bcache[j] = mult
        else
            mult = A.bcache[j]
        end
        for l in eachindex(A.t)
            if l != j
                k = one(A.t[1])
                for m in eachindex(A.t)
                    if m != j && m != l
                        k *= (t - A.t[m])
                    end
                end
                k *= inv(mult)
                tmp += k
            end
        end
        der += A.u[:, j] * tmp
    end
    der
end

function _derivative(A::LagrangeInterpolation{<:AbstractVector}, t::Number, idx)
    _derivative(A, t)
end
function _derivative(A::LagrangeInterpolation{<:AbstractMatrix}, t::Number, idx)
    _derivative(A, t)
end

function _derivative(A::AkimaInterpolation{<:AbstractVector}, t::Number, iguess)
    idx = get_idx(A, t, iguess; idx_shift = -1, side = :first)
    j = min(idx, length(A.c))  # for smooth derivative at A.t[end]
    wj = t - A.t[idx]
    @evalpoly wj A.b[idx] 2A.c[j] 3A.d[j]
end

function _derivative(A::ConstantInterpolation, t::Number, iguess)
    return zero(first(A.u))
end

function _derivative(A::ConstantInterpolation{<:AbstractVector}, t::Number, iguess)
    return isempty(searchsorted(A.t, t)) ? zero(A.u[1]) : typed_nan(A.u)
end

function _derivative(A::ConstantInterpolation{<:AbstractMatrix}, t::Number, iguess)
    return isempty(searchsorted(A.t, t)) ? zero(A.u[:, 1]) : typed_nan(A.u) .* A.u[:, 1]
end

# QuadraticSpline Interpolation
function _derivative(A::QuadraticSpline{<:AbstractVector}, t::Number, iguess)
    idx = get_idx(A, t, iguess)
    α, β = get_parameters(A, idx)
    Δt = t - A.t[idx]
    Δt_full = A.t[idx + 1] - A.t[idx]
    2α * Δt / Δt_full^2 + β / Δt_full
end

# CubicSpline Interpolation
function _derivative(A::CubicSpline{<:AbstractVector}, t::Number, iguess)
    idx = get_idx(A, t, iguess)
    Δt₁ = t - A.t[idx]
    Δt₂ = A.t[idx + 1] - t
    dI = (-A.z[idx] * Δt₂^2 + A.z[idx + 1] * Δt₁^2) / (2A.h[idx + 1])
    c₁, c₂ = get_parameters(A, idx)
    dC = c₁
    dD = -c₂
    dI + dC + dD
end

function _derivative(A::BSplineInterpolation{<:AbstractVector{<:Number}}, t::Number, iguess)
    # change t into param [0 1]
    t < A.t[1] && return zero(A.u[1])
    t > A.t[end] && return zero(A.u[end])
    idx = get_idx(A, t, iguess)
    n = length(A.t)
    scale = (A.p[idx + 1] - A.p[idx]) / (A.t[idx + 1] - A.t[idx])
    t_ = A.p[idx] + (t - A.t[idx]) * scale
    sc = t isa ForwardDiff.Dual ? zeros(eltype(t), n) : A.sc
    spline_coefficients!(sc, A.d - 1, A.k, t_)
    ducum = zero(eltype(A.u))
    if t == A.t[1]
        ducum = (A.c[2] - A.c[1]) / (A.k[A.d + 2])
    else
        for i in 1:(n - 1)
            ducum += sc[i + 1] * (A.c[i + 1] - A.c[i]) / (A.k[i + A.d + 1] - A.k[i + 1])
        end
    end
    ducum * A.d * scale
end

function _derivative(
        A::BSplineInterpolation{<:AbstractArray{<:Number}}, t::Number, iguess)
    # change t into param [0 1]
    ax_u = axes(A.u)[1:(end - 1)]
    t < A.t[1] && return zeros(size(A.u)[1:(end - 1)]...)
    t > A.t[end] && return zeros(size(A.u)[1:(end - 1)]...)
    idx = get_idx(A, t, iguess)
    n = length(A.t)
    scale = (A.p[idx + 1] - A.p[idx]) / (A.t[idx + 1] - A.t[idx])
    t_ = A.p[idx] + (t - A.t[idx]) * scale
    sc = t isa ForwardDiff.Dual ? zeros(eltype(t), n) : A.sc
    spline_coefficients!(sc, A.d - 1, A.k, t_)
    ducum = zeros(size(A.u)[1:(end - 1)]...)
    if t == A.t[1]
        ducum = (A.c[ax_u..., 2] - A.c[ax_u..., 1]) / (A.k[A.d + 2])
    else
        for i in 1:(n - 1)
            ducum = ducum +
                    sc[i + 1] * (A.c[ax_u..., i + 1] - A.c[ax_u..., i]) /
                    (A.k[i + A.d + 1] - A.k[i + 1])
        end
    end
    ducum * A.d * scale
end
# BSpline Curve Approx
function _derivative(A::BSplineApprox{<:AbstractVector{<:Number}}, t::Number, iguess)
    # change t into param [0 1]
    t < A.t[1] && return zero(A.u[1])
    t > A.t[end] && return zero(A.u[end])
    idx = get_idx(A, t, iguess)
    scale = (A.p[idx + 1] - A.p[idx]) / (A.t[idx + 1] - A.t[idx])
    t_ = A.p[idx] + (t - A.t[idx]) * scale
    sc = t isa ForwardDiff.Dual ? zeros(eltype(t), A.h) : A.sc
    spline_coefficients!(sc, A.d - 1, A.k, t_)
    ducum = zero(eltype(A.u))
    if t == A.t[1]
        ducum = (A.c[2] - A.c[1]) / (A.k[A.d + 2])
    else
        for i in 1:(A.h - 1)
            ducum += sc[i + 1] * (A.c[i + 1] - A.c[i]) / (A.k[i + A.d + 1] - A.k[i + 1])
        end
    end
    ducum * A.d * scale
end

function _derivative(
        A::BSplineApprox{<:AbstractArray{<:Number}}, t::Number, iguess)
    # change t into param [0 1]
    ax_u = axes(A.u)[1:(end - 1)]
    t < A.t[1] && return zeros(size(A.u)[1:(end - 1)]...)
    t > A.t[end] && return zeros(size(A.u)[1:(end - 1)]...)
    idx = get_idx(A, t, iguess)
    scale = (A.p[idx + 1] - A.p[idx]) / (A.t[idx + 1] - A.t[idx])
    t_ = A.p[idx] + (t - A.t[idx]) * scale
    sc = t isa ForwardDiff.Dual ? zeros(eltype(t), A.h) : A.sc
    spline_coefficients!(sc, A.d - 1, A.k, t_)
    ducum = zeros(size(A.u)[1:(end - 1)]...)
    if t == A.t[1]
        ducum = (A.c[ax_u..., 2] - A.c[ax_u..., 1]) / (A.k[A.d + 2])
    else
        for i in 1:(A.h - 1)
            ducum = ducum +
                    sc[i + 1] * (A.c[ax_u..., i + 1] - A.c[ax_u..., i]) /
                    (A.k[i + A.d + 1] - A.k[i + 1])
        end
    end
    ducum * A.d * scale
end
# Cubic Hermite Spline
function _derivative(
        A::CubicHermiteSpline{<:AbstractVector{<:Number}}, t::Number, iguess)
    idx = get_idx(A, t, iguess)
    Δt₀ = t - A.t[idx]
    Δt₁ = t - A.t[idx + 1]
    out = A.du[idx]
    c₁, c₂ = get_parameters(A, idx)
    out += Δt₀ * (Δt₀ * c₂ + 2(c₁ + Δt₁ * c₂))
    out
end

# Quintic Hermite Spline
function _derivative(
        A::QuinticHermiteSpline{<:AbstractVector{<:Number}}, t::Number, iguess)
    idx = get_idx(A, t, iguess)
    Δt₀ = t - A.t[idx]
    Δt₁ = t - A.t[idx + 1]
    out = A.du[idx] + A.ddu[idx] * Δt₀
    c₁, c₂, c₃ = get_parameters(A, idx)
    out += Δt₀^2 *
           (3c₁ + (3Δt₁ + Δt₀) * c₂ + (3Δt₁^2 + Δt₀ * 2Δt₁) * c₃)
    out
end

function _derivative(A::SmoothArcLengthInterpolation, t::Number, iguess)
    (; out) = A
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
        S, C = sincos(t_circle_seg / A.radius[idx])
        v₁ = view(A.dir_1, :, idx)
        v₂ = view(A.dir_2, :, idx)
        Rⱼ = A.radius[idx]
        @. out = (-S * v₁ + C * v₂) / Rⱼ
    else
        if short_side_left
            d₁ = view(A.d, :, idx + 1)
            @. out = d₁
        else
            d₀ = view(A.d, :, idx)
            @. out = d₀
        end
    end

    out
end
