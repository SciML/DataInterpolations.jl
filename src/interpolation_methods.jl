function _interpolate(A, t)
    ((t < A.t[1] || t > A.t[end]) && !A.extrapolate) &&
        throw(ExtrapolationError())
    val, idx_prev = _interpolate(A, t, A.idx_prev[])
    A.idx_prev[] = idx_prev
    return val
end

# Linear Interpolation
function _interpolate(A::LinearInterpolation{<:AbstractVector}, t::Number, iguess)
    if isnan(t)
        # For correct derivative with NaN
        idx = firstindex(A.u) - 1
        t1 = t2 = one(eltype(A.t))
        u1 = u2 = one(eltype(A.u))
        slope = t * one(eltype(A.p.slope))
    else
        idx = get_idx(A.t, t, iguess)
        t1, t2 = A.t[idx], A.t[idx + 1]
        u1, u2 = A.u[idx], A.u[idx + 1]
        slope = A.p.slope[idx]
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

    val, idx
end

function _interpolate(A::LinearInterpolation{<:AbstractMatrix}, t::Number, iguess)
    idx = get_idx(A.t, t, iguess)
    Δt = t - A.t[idx]
    return A.u[:, idx] + A.p.slope[idx] * Δt, idx
end

# Quadratic Interpolation
_quad_interp_indices(A, t) = _quad_interp_indices(A, t, firstindex(A.t) - 1)
function _quad_interp_indices(A::QuadraticInterpolation, t::Number, iguess)
    idx = get_idx(A.t, t, iguess; idx_shift = A.mode == :Backward ? -1 : 0, ub_shift = -2)
    idx, idx + 1, idx + 2
end

function _interpolate(A::QuadraticInterpolation, t::Number, iguess)
    i₀, i₁, i₂ = _quad_interp_indices(A, t, iguess)
    u₀ = A.p.l₀[i₀] * (t - A.t[i₁]) * (t - A.t[i₂])
    u₁ = A.p.l₁[i₀] * (t - A.t[i₀]) * (t - A.t[i₂])
    u₂ = A.p.l₂[i₀] * (t - A.t[i₀]) * (t - A.t[i₁])
    return u₀ + u₁ + u₂, i₀
end

# Lagrange Interpolation
function _interpolate(A::LagrangeInterpolation{<:AbstractVector}, t::Number, iguess)
    idx = get_idx(A.t, t, iguess)
    findRequiredIdxs!(A, t, idx)
    if A.t[A.idxs[1]] == t
        return A.u[A.idxs[1]], idx
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
    N / D, idx
end

function _interpolate(A::LagrangeInterpolation{<:AbstractMatrix}, t::Number, iguess)
    idx = get_idx(A.t, t, iguess)
    findRequiredIdxs!(A, t, idx)
    if A.t[A.idxs[1]] == t
        return A.u[:, A.idxs[1]], idx
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
    N / D, idx
end

function _interpolate(A::AkimaInterpolation{<:AbstractVector}, t::Number, iguess)
    idx = get_idx(A.t, t, iguess)
    wj = t - A.t[idx]
    (@evalpoly wj A.u[idx] A.b[idx] A.c[idx] A.d[idx]), idx
end

# ConstantInterpolation Interpolation
function _interpolate(A::ConstantInterpolation{<:AbstractVector}, t::Number, iguess)
    if A.dir === :left
        # :left means that value to the left is used for interpolation
        idx = get_idx(A.t, t, iguess; lb = 1, ub_shift = 0)
    else
        # :right means that value to the right is used for interpolation
        idx = get_idx(A.t, t, iguess; side = :first, lb = 1, ub_shift = 0)
    end
    A.u[idx], idx
end

function _interpolate(A::ConstantInterpolation{<:AbstractMatrix}, t::Number, iguess)
    if A.dir === :left
        # :left means that value to the left is used for interpolation
        idx = get_idx(A.t, t, iguess; lb = 1, ub_shift = 0)
    else
        # :right means that value to the right is used for interpolation
        idx = get_idx(A.t, t, iguess; side = :first, lb = 1, ub_shift = 0)
    end
    A.u[:, idx], idx
end

# QuadraticSpline Interpolation
function _interpolate(A::QuadraticSpline{<:AbstractVector}, t::Number, iguess)
    idx = get_idx(A.t, t, iguess)
    Cᵢ = A.u[idx]
    Δt = t - A.t[idx]
    return A.p.z[idx] * Δt + A.p.σ[idx] * Δt^2 + Cᵢ, idx
end

# CubicSpline Interpolation
function _interpolate(A::CubicSpline{<:AbstractVector}, t::Number, iguess)
    idx = get_idx(A.t, t, iguess)
    Δt₁ = t - A.t[idx]
    Δt₂ = A.t[idx + 1] - t
    I = (A.z[idx] * Δt₂^3 + A.z[idx + 1] * Δt₁^3) / (6A.h[idx + 1])
    C = A.p.c₁[idx] * Δt₁
    D = A.p.c₂[idx] * Δt₂
    I + C + D, idx
end

# BSpline Curve Interpolation
function _interpolate(A::BSplineInterpolation{<:AbstractVector{<:Number}},
        t::Number,
        iguess)
    t < A.t[1] && return A.u[1], 1
    t > A.t[end] && return A.u[end], lastindex(t)
    # change t into param [0 1]
    idx = get_idx(A.t, t, iguess)
    t = A.p[idx] + (t - A.t[idx]) / (A.t[idx + 1] - A.t[idx]) * (A.p[idx + 1] - A.p[idx])
    n = length(A.t)
    N = t isa ForwardDiff.Dual ? zeros(eltype(t), n) : A.N
    nonzero_coefficient_idxs = spline_coefficients!(N, A.d, A.k, t)
    ucum = zero(eltype(A.u))
    for i in nonzero_coefficient_idxs
        ucum += N[i] * A.c[i]
    end
    ucum, idx
end

# BSpline Curve Approx
function _interpolate(A::BSplineApprox{<:AbstractVector{<:Number}}, t::Number, iguess)
    t < A.t[1] && return A.u[1], 1
    t > A.t[end] && return A.u[end], lastindex(t)
    # change t into param [0 1]
    idx = get_idx(A.t, t, iguess)
    t = A.p[idx] + (t - A.t[idx]) / (A.t[idx + 1] - A.t[idx]) * (A.p[idx + 1] - A.p[idx])
    N = t isa ForwardDiff.Dual ? zeros(eltype(t), A.h) : A.N
    nonzero_coefficient_idxs = spline_coefficients!(N, A.d, A.k, t)
    ucum = zero(eltype(A.u))
    for i in nonzero_coefficient_idxs
        ucum += N[i] * A.c[i]
    end
    ucum, idx
end

# Cubic Hermite Spline
function _interpolate(
        A::CubicHermiteSpline{<:AbstractVector{<:Number}}, t::Number, iguess)
    idx = get_idx(A.t, t, iguess)
    Δt₀ = t - A.t[idx]
    Δt₁ = t - A.t[idx + 1]
    out = A.u[idx] + Δt₀ * A.du[idx]
    out += Δt₀^2 * (A.p.c₁[idx] + Δt₁ * A.p.c₂[idx])
    out, idx
end

# Quintic Hermite Spline
function _interpolate(
        A::QuinticHermiteSpline{<:AbstractVector{<:Number}}, t::Number, iguess)
    idx = get_idx(A.t, t, iguess)
    Δt₀ = t - A.t[idx]
    Δt₁ = t - A.t[idx + 1]
    out = A.u[idx] + Δt₀ * (A.du[idx] + A.ddu[idx] * Δt₀ / 2)
    out += Δt₀^3 * (A.p.c₁[idx] + Δt₁ * (A.p.c₂[idx] + A.p.c₃[idx] * Δt₁))
    out, idx
end
