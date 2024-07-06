function derivative(A, t, order = 1)
    ((t < A.t[1] || t > A.t[end]) && !A.extrapolate) && throw(ExtrapolationError())
    iguess = A.idx_prev[]

    return if order == 1
        val, idx = _derivative(A, t, iguess)
        A.idx_prev[] = idx
        val
    elseif order == 2
        ForwardDiff.derivative(t -> begin
                val, idx = _derivative(A, t, iguess)
                A.idx_prev[] = idx
                val
            end, t)
    else
        throw(DerivativeNotFoundError())
    end
end

function _derivative(A::LinearInterpolation, t::Number, iguess)
    idx = get_idx(A.t, t, iguess; idx_shift = -1, ub_shift = -2, side = :first)
    A.p.slope[idx], idx
end

function _derivative(A::QuadraticInterpolation, t::Number, iguess)
    i₀, i₁, i₂ = _quad_interp_indices(A, t, iguess)
    du₀ = A.p.l₀[i₀] * (2t - A.t[i₁] - A.t[i₂])
    du₁ = A.p.l₁[i₀] * (2t - A.t[i₀] - A.t[i₂])
    du₂ = A.p.l₂[i₀] * (2t - A.t[i₀] - A.t[i₁])
    return @views @. du₀ + du₁ + du₂, i₀
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
    _derivative(A, t), idx
end
function _derivative(A::LagrangeInterpolation{<:AbstractMatrix}, t::Number, idx)
    _derivative(A, t), idx
end

function _derivative(A::AkimaInterpolation{<:AbstractVector}, t::Number, iguess)
    idx = get_idx(A.t, t, iguess; idx_shift = -1, side = :first)
    j = min(idx, length(A.c))  # for smooth derivative at A.t[end]
    wj = t - A.t[idx]
    (@evalpoly wj A.b[idx] 2A.c[j] 3A.d[j]), idx
end

function _derivative(A::ConstantInterpolation, t::Number, iguess)
    return zero(first(A.u)), iguess
end

function _derivative(A::ConstantInterpolation{<:AbstractVector}, t::Number)
    ((t < A.t[1] || t > A.t[end]) && !A.extrapolate) && throw(ExtrapolationError())
    return isempty(searchsorted(A.t, t)) ? zero(A.u[1]) : eltype(A.u)(NaN)
end

function _derivative(A::ConstantInterpolation{<:AbstractMatrix}, t::Number)
    ((t < A.t[1] || t > A.t[end]) && !A.extrapolate) && throw(ExtrapolationError())
    return isempty(searchsorted(A.t, t)) ? zero(A.u[:, 1]) : eltype(A.u)(NaN) .* A.u[:, 1]
end

# QuadraticSpline Interpolation
function _derivative(A::QuadraticSpline{<:AbstractVector}, t::Number, iguess)
    idx = get_idx(A.t, t, iguess; lb = 2, ub_shift = 0, side = :first)
    σ = A.p.σ[idx - 1]
    A.z[idx - 1] + 2σ * (t - A.t[idx - 1]), idx
end

# CubicSpline Interpolation
function _derivative(A::CubicSpline{<:AbstractVector}, t::Number, iguess)
    idx = get_idx(A.t, t, iguess)
    Δt₁ = t - A.t[idx]
    Δt₂ = A.t[idx + 1] - t
    dI = (-A.z[idx] * Δt₂^2 + A.z[idx + 1] * Δt₁^2) / (2A.h[idx + 1])
    dC = A.p.c₁[idx]
    dD = -A.p.c₂[idx]
    dI + dC + dD, idx
end

function _derivative(A::BSplineInterpolation{<:AbstractVector{<:Number}}, t::Number, iguess)
    # change t into param [0 1]
    t < A.t[1] && return zero(A.u[1]), 1
    t > A.t[end] && return zero(A.u[end]), lastindex(t)
    idx = get_idx(A.t, t, iguess)
    n = length(A.t)
    scale = (A.p[idx + 1] - A.p[idx]) / (A.t[idx + 1] - A.t[idx])
    t_ = A.p[idx] + (t - A.t[idx]) * scale
    N = t isa ForwardDiff.Dual ? zeros(eltype(t), n) : A.N
    spline_coefficients!(N, A.d - 1, A.k, t_)
    ducum = zero(eltype(A.u))
    if t == A.t[1]
        ducum = (A.c[2] - A.c[1]) / (A.k[A.d + 2])
    else
        for i in 1:(n - 1)
            ducum += N[i + 1] * (A.c[i + 1] - A.c[i]) / (A.k[i + A.d + 1] - A.k[i + 1])
        end
    end
    ducum * A.d * scale, idx
end

# BSpline Curve Approx
function _derivative(A::BSplineApprox{<:AbstractVector{<:Number}}, t::Number, iguess)
    # change t into param [0 1]
    t < A.t[1] && return zero(A.u[1]), 1
    t > A.t[end] && return zero(A.u[end]), lastindex(t)
    idx = get_idx(A.t, t, iguess)
    scale = (A.p[idx + 1] - A.p[idx]) / (A.t[idx + 1] - A.t[idx])
    t_ = A.p[idx] + (t - A.t[idx]) * scale
    N = t isa ForwardDiff.Dual ? zeros(eltype(t), A.h) : A.N
    spline_coefficients!(N, A.d - 1, A.k, t_)
    ducum = zero(eltype(A.u))
    if t == A.t[1]
        ducum = (A.c[2] - A.c[1]) / (A.k[A.d + 2])
    else
        for i in 1:(A.h - 1)
            ducum += N[i + 1] * (A.c[i + 1] - A.c[i]) / (A.k[i + A.d + 1] - A.k[i + 1])
        end
    end
    ducum * A.d * scale, idx
end

# Cubic Hermite Spline
function _derivative(
        A::CubicHermiteSpline{<:AbstractVector{<:Number}}, t::Number, iguess)
    idx = get_idx(A.t, t, iguess)
    Δt₀ = t - A.t[idx]
    Δt₁ = t - A.t[idx + 1]
    out = A.du[idx]
    out += Δt₀ * (Δt₀ * A.p.c₂[idx] + 2(A.p.c₁[idx] + Δt₁ * A.p.c₂[idx]))
    out, idx
end

# Quintic Hermite Spline
function _derivative(
        A::QuinticHermiteSpline{<:AbstractVector{<:Number}}, t::Number, iguess)
    idx = get_idx(A.t, t, iguess)
    Δt₀ = t - A.t[idx]
    Δt₁ = t - A.t[idx + 1]
    out = A.du[idx] + A.ddu[idx] * Δt₀
    out += Δt₀^2 *
           (3A.p.c₁[idx] + (3Δt₁ + Δt₀) * A.p.c₂[idx] + (3Δt₁^2 + Δt₀ * 2Δt₁) * A.p.c₃[idx])
    out, idx
end
