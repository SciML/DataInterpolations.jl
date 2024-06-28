function derivative(A, t, order = 1)
    order > 2 && throw(DerivativeNotFoundError())
    ((t < A.t[1] || t > A.t[end]) && !A.extrapolate) && throw(ExtrapolationError())
    order == 1 && return _derivative(A, t, firstindex(A.t) - 1)[1]
    return ForwardDiff.derivative(t -> _derivative(A, t, firstindex(A.t) - 1)[1], t)
end

function _derivative(A::LinearInterpolation{<:AbstractVector}, t::Number, iguess)
    idx = searchsortedfirstcorrelated(A.t, t, iguess)
    idx > length(A.t) ? idx -= 1 : nothing
    idx -= 1
    idx == 0 ? idx += 1 : nothing
    A.p.slope[idx], idx
end

function _derivative(A::LinearInterpolation{<:AbstractMatrix}, t::Number, iguess)
    idx = searchsortedfirstcorrelated(A.t, t, iguess)
    idx > length(A.t) ? idx -= 1 : nothing
    idx -= 1
    idx == 0 ? idx += 1 : nothing
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

_derivative(A::LagrangeInterpolation{<:AbstractVector}, t::Number, i) = _derivative(A, t), i
_derivative(A::LagrangeInterpolation{<:AbstractMatrix}, t::Number, i) = _derivative(A, t), i

function _derivative(A::AkimaInterpolation{<:AbstractVector}, t::Number, iguess)
    i = searchsortedfirstcorrelated(A.t, t, iguess)
    i > length(A.t) ? i -= 1 : nothing
    i -= 1
    i == 0 ? i += 1 : nothing
    j = min(i, length(A.c))  # for smooth derivative at A.t[end]
    wj = t - A.t[i]
    (@evalpoly wj A.b[i] 2A.c[j] 3A.d[j]), i
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
    idx = searchsortedfirstcorrelated(A.t, t, iguess)
    idx > length(A.t) ? idx -= 1 : nothing
    idx == 1 ? idx += 1 : nothing
    σ = 1 // 2 * (A.z[idx] - A.z[idx - 1]) / (A.t[idx] - A.t[idx - 1])
    A.z[idx - 1] + 2σ * (t - A.t[idx - 1]), idx
end

# CubicSpline Interpolation
function _derivative(A::CubicSpline{<:AbstractVector}, t::Number, iguess)
    i = searchsortedfirstcorrelated(A.t, t, iguess)
    i > length(A.t) ? i -= 1 : nothing
    i -= 1
    i == 0 ? i += 1 : nothing
    dI = -3A.z[i] * (A.t[i + 1] - t)^2 / (6A.h[i + 1]) +
         3A.z[i + 1] * (t - A.t[i])^2 / (6A.h[i + 1])
    dC = A.u[i + 1] / A.h[i + 1] - A.z[i + 1] * A.h[i + 1] / 6
    dD = -(A.u[i] / A.h[i + 1] - A.z[i] * A.h[i + 1] / 6)
    dI + dC + dD, i
end

function _derivative(A::BSplineInterpolation{<:AbstractVector{<:Number}}, t::Number, iguess)
    # change t into param [0 1]
    t < A.t[1] && return zero(A.u[1]), 1
    t > A.t[end] && return zero(A.u[end]), lastindex(t)
    idx = searchsortedlastcorrelated(A.t, t, iguess)
    idx == length(A.t) ? idx -= 1 : nothing
    n = length(A.t)
    scale = (A.p[idx + 1] - A.p[idx]) / (A.t[idx + 1] - A.t[idx])
    t_ = A.p[idx] + (t - A.t[idx]) * scale
    N = DataInterpolations.spline_coefficients(n, A.d - 1, A.k, t_)
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
    idx = searchsortedlastcorrelated(A.t, t, iguess)
    idx == length(A.t) ? idx -= 1 : nothing
    scale = (A.p[idx + 1] - A.p[idx]) / (A.t[idx + 1] - A.t[idx])
    t_ = A.p[idx] + (t - A.t[idx]) * scale
    N = spline_coefficients(A.h, A.d - 1, A.k, t_)
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
