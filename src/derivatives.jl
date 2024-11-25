function derivative(A, t, order = 1)
    ((t < A.t[1] || t > A.t[end]) && !A.extrapolate) && throw(ExtrapolationError())
    iguess = A.iguesser

    return if order == 1
        _derivative(A, t, iguess)
    elseif order == 2
        ForwardDiff.derivative(t -> begin
                _derivative(A, t, iguess)
            end, t)
    else
        throw(DerivativeNotFoundError())
    end
end

function _derivative(A::LinearInterpolation, t::Number, iguess)
    idx = get_idx(A, t, iguess; idx_shift = -1, ub_shift = -1, side = :first)
    slope = get_parameters(A, idx)
    slope
end

function _derivative(A::SmoothedConstantInterpolation{<:AbstractVector}, t::Number, iguess)
    idx = get_idx(A, t, iguess)
    d_lower, d_upper, c_lower, c_upper = get_parameters(A, idx)

    if (t - A.t[idx]) < d_lower
        -2c_lower * ((t - A.t[idx]) / d_lower - 1) / d_lower
    elseif (A.t[idx + 1] - t) < d_upper
        2c_upper * (1 - (A.t[idx + 1] - t) / d_upper) / d_upper
    else
        zero(c_upper / oneunit(t))
    end
end

function _derivative(A::QuadraticInterpolation, t::Number, iguess)
    idx = get_idx(A, t, iguess)
    Δt = t - A.t[idx]
    α, β = get_parameters(A, idx)
    return 2α * Δt + β
end

function _derivative(A::LagrangeInterpolation{<:AbstractVector}, t::Number)
    ((t < A.t[1] || t > A.t[end]) && !A.extrapolate) && throw(ExtrapolationError())
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
    ((t < A.t[1] || t > A.t[end]) && !A.extrapolate) && throw(ExtrapolationError())
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
    ((t < A.t[1] || t > A.t[end]) && !A.extrapolate) && throw(ExtrapolationError())
    return isempty(searchsorted(A.t, t)) ? zero(A.u[1]) : eltype(A.u)(NaN)
end

function _derivative(A::ConstantInterpolation{<:AbstractMatrix}, t::Number, iguess)
    ((t < A.t[1] || t > A.t[end]) && !A.extrapolate) && throw(ExtrapolationError())
    return isempty(searchsorted(A.t, t)) ? zero(A.u[:, 1]) : eltype(A.u)(NaN) .* A.u[:, 1]
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
        A::BSplineInterpolation{<:AbstractArray{<:Number, N}}, t::Number, iguess) where {N}
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
        A::BSplineApprox{<:AbstractArray{<:Number, N}}, t::Number, iguess) where {N}
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
