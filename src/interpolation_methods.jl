function _interpolate(A, t)
    ((t < A.t[1] || t > A.t[end]) && !A.extrapolate) &&
        throw(ExtrapolationError())
    return _interpolate(A, t, A.iguesser)
end

# Linear Interpolation
function _interpolate(A::LinearInterpolation{<:AbstractVector}, t::Number, iguess)
    if isnan(t)
        # For correct derivative with NaN
        idx = firstindex(A.u)
        t1 = t2 = oneunit(eltype(A.t))
        u1 = u2 = oneunit(eltype(A.u))
        slope = t/t * get_parameters(A, idx)
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

# ConstantInterpolation Interpolation
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

function _interpolate(A::CubicSpline{<:AbstractArray{T, N}}, t::Number, iguess) where {T, N}
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

function _interpolate(A::BSplineInterpolation{<:AbstractArray{T, N}},
        t::Number,
        iguess) where {T <: Number, N}
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
        A::BSplineApprox{<:AbstractArray{T, N}}, t::Number, iguess) where {T <: Number, N}
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
        A::CubicHermiteSpline{<:AbstractVector{<:Number}}, t::Number, iguess)
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
