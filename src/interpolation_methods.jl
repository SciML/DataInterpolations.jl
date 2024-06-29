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
function _interpolate(A::LagrangeInterpolation{<:AbstractVector}, t::Number)
    ((t < A.t[1] || t > A.t[end]) && !A.extrapolate) && throw(ExtrapolationError())
    idxs = findRequiredIdxs(A, t)
    if A.t[idxs[1]] == t
        return A.u[idxs[1]]
    end
    N = zero(A.u[1])
    D = zero(A.t[1])
    tmp = N
    for i in 1:length(idxs)
        if isnan(A.bcache[idxs[i]])
            mult = one(A.t[1])
            for j in 1:(i - 1)
                mult *= (A.t[idxs[i]] - A.t[idxs[j]])
            end
            for j in (i + 1):length(idxs)
                mult *= (A.t[idxs[i]] - A.t[idxs[j]])
            end
            A.bcache[idxs[i]] = mult
        else
            mult = A.bcache[idxs[i]]
        end
        tmp = inv((t - A.t[idxs[i]]) * mult)
        D += tmp
        N += (tmp * A.u[idxs[i]])
    end
    N / D
end

function _interpolate(A::LagrangeInterpolation{<:AbstractMatrix}, t::Number)
    ((t < A.t[1] || t > A.t[end]) && !A.extrapolate) && throw(ExtrapolationError())
    idxs = findRequiredIdxs(A, t)
    if A.t[idxs[1]] == t
        return A.u[:, idxs[1]]
    end
    N = zero(A.u[:, 1])
    D = zero(A.t[1])
    tmp = D
    for i in 1:length(idxs)
        if isnan(A.bcache[idxs[i]])
            mult = one(A.t[1])
            for j in 1:(i - 1)
                mult *= (A.t[idxs[i]] - A.t[idxs[j]])
            end
            for j in (i + 1):length(idxs)
                mult *= (A.t[idxs[i]] - A.t[idxs[j]])
            end
            A.bcache[idxs[i]] = mult
        else
            mult = A.bcache[idxs[i]]
        end
        tmp = inv((t - A.t[idxs[i]]) * mult)
        D += tmp
        @. N += (tmp * A.u[:, idxs[i]])
    end
    N / D
end

function _interpolate(A::LagrangeInterpolation{<:AbstractVector}, t::Number, idx)
    _interpolate(A, t), idx
end

function _interpolate(A::LagrangeInterpolation{<:AbstractMatrix}, t::Number, idx)
    _interpolate(A, t), idx
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
    idx = get_idx(A.t, t, iguess; lb = 2, ub_shift = 0, side = :first)
    Cᵢ = A.u[idx - 1]
    σ = 1 // 2 * (A.z[idx] - A.z[idx - 1]) / (A.t[idx] - A.t[idx - 1])
    return A.z[idx - 1] * (t - A.t[idx - 1]) + σ * (t - A.t[idx - 1])^2 + Cᵢ, idx
end

# CubicSpline Interpolation
function _interpolate(A::CubicSpline{<:AbstractVector}, t::Number, iguess)
    idx = get_idx(A.t, t, iguess)
    I = A.z[idx] * (A.t[idx + 1] - t)^3 / (6A.h[idx + 1]) +
        A.z[idx + 1] * (t - A.t[idx])^3 / (6A.h[idx + 1])
    C = (A.u[idx + 1] / A.h[idx + 1] - A.z[idx + 1] * A.h[idx + 1] / 6) * (t - A.t[idx])
    D = (A.u[idx] / A.h[idx + 1] - A.z[idx] * A.h[idx + 1] / 6) * (A.t[idx + 1] - t)
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
    N = spline_coefficients(n, A.d, A.k, t)
    ucum = zero(eltype(A.u))
    for i in 1:n
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
    N = spline_coefficients(A.h, A.d, A.k, t)
    ucum = zero(eltype(A.u))
    for i in 1:(A.h)
        ucum += N[i] * A.c[i]
    end
    ucum, idx
end
