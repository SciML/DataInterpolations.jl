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

function _derivative(A::QuadraticInterpolation, t::Number, iguess)
    i₀, i₁, i₂ = _quad_interp_indices(A, t, iguess)
    l₀, l₁, l₂ = get_parameters(A, i₀)
    du₀ = l₀ * (2t - A.t[i₁] - A.t[i₂])
    du₁ = l₁ * (2t - A.t[i₀] - A.t[i₂])
    du₂ = l₂ * (2t - A.t[i₀] - A.t[i₁])
    return @views @. du₀ + du₁ + du₂
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

function _derivative(A::LagrangeInterpolation{<:AbstractArray}, t::Number)
    ((t < A.t[1] || t > A.t[end]) && !A.extrapolate) && throw(ExtrapolationError())
    ax = axes(A.u)[1:(end - 1)]
    der = zero(A.u[ax..., 1])
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
        der += A.u[ax..., j] * tmp
    end
    der
end

function _derivative(A::LagrangeInterpolation{<:AbstractVector}, t::Number, idx)
    _derivative(A, t)
end
function _derivative(A::LagrangeInterpolation{<:AbstractArray}, t::Number, idx)
    _derivative(A, t)
end

function _derivative(A::AkimaInterpolation{<:AbstractVector}, t::Number, iguess)
    idx = get_idx(A, t, iguess; idx_shift = -1, side = :first)
    j = min(idx, length(A.c))  # for smooth derivative at A.t[end]
    wj = t - A.t[idx]
    @evalpoly wj A.b[idx] 2A.c[j] 3A.d[j]
end

function _derivative(A::AkimaInterpolation{<:AbstractArray}, t::Number, iguess)
    idx = get_idx(A, t, iguess; idx_shift = -1, side = :first)
    j = min(idx, length(A.c))  # for smooth derivative at A.t[end]
    wj = t - A.t[idx]
    ax = axes(A.u)[1:(end - 1)]
    @. @evalpoly wj A.b[ax..., idx] 2A.c[ax..., j] 3A.d[ax..., j]
end

function _derivative(A::ConstantInterpolation, t::Number, iguess)
    return zero(first(A.u))
end

function _derivative(A::ConstantInterpolation{<:AbstractVector}, t::Number, iguess)
    ((t < A.t[1] || t > A.t[end]) && !A.extrapolate) && throw(ExtrapolationError())
    return isempty(searchsorted(A.t, t)) ? zero(A.u[1]) : eltype(A.u)(NaN)
end

function _derivative(A::ConstantInterpolation{<:AbstractArray}, t::Number, iguess)
    ((t < A.t[1] || t > A.t[end]) && !A.extrapolate) && throw(ExtrapolationError())
    ax = axes(A.u)[1:(end - 1)]
    return isempty(searchsorted(A.t, t)) ? zero(A.u[ax..., 1]) :
           eltype(A.u)(NaN) .* A.u[ax..., 1]
end

# QuadraticSpline Interpolation
function _derivative(A::QuadraticSpline, t::Number, iguess)
    idx = get_idx(A, t, iguess; lb = 2, ub_shift = 0, side = :first)
    σ = get_parameters(A, idx - 1)
    A.z[idx - 1] + 2σ * (t - A.t[idx - 1])
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

function _derivative(A::CubicSpline{<:AbstractArray}, t::Number, iguess)
    idx = get_idx(A, t, iguess)
    Δt₁ = t - A.t[idx]
    Δt₂ = A.t[idx + 1] - t
    ax = axes(A.u)[1:(end - 1)]
    dI = (-A.z[ax..., idx] * Δt₂^2 + A.z[ax..., idx + 1] * Δt₁^2) / (2A.h[idx + 1])
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

function _derivative(
        A::CubicHermiteSpline{<:AbstractArray}, t::Number, iguess)
    idx = get_idx(A, t, iguess)
    Δt₀ = t - A.t[idx]
    Δt₁ = t - A.t[idx + 1]
    ax = axes(A.u)[1:(end - 1)]
    out = A.du[ax..., idx]
    c₁, c₂ = get_parameters(A, idx)
    out .+= Δt₀ .* (Δt₀ .* c₂ .+ 2(c₁ .+ Δt₁ .* c₂))
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

function _derivative(
        A::QuinticHermiteSpline{<:AbstractArray}, t::Number, iguess)
    idx = get_idx(A, t, iguess)
    Δt₀ = t - A.t[idx]
    Δt₁ = t - A.t[idx + 1]
    ax = axes(A.u)[1:(end - 1)]
    out = A.du[ax..., idx] + A.ddu[ax..., idx] * Δt₀
    c₁, c₂, c₃ = get_parameters(A, idx)
    out .+= Δt₀^2 .*
            (3c₁ .+ (3Δt₁ .+ Δt₀) .* c₂ + (3Δt₁^2 .+ Δt₀ .* 2Δt₁) .* c₃)
    out
end
