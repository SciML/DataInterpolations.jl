"""
    integral(A::AbstractInterpolation, t::Number)
    integral(A::AbstractInterpolation, t1::Number, t2::Number)

Evaluate the integral of the interpolation `A` over `[first(A.t), t]`, or over `[t1, t2]`
when both bounds are given. The integral is computed analytically per interval, so it is
exact for the interpolation. `integral(A, t2, t1) == -integral(A, t1, t2)`.

Parts of the integration interval outside the range of `A.t` are integrated over the
extrapolation, as determined by the `extrapolation_left` and `extrapolation_right`
settings of `A`.

Interpolation types with no analytical antiderivative — `LagrangeInterpolation` and
`Curvefit` — throw an `IntegralNotFoundError`; use a numerical integrator such as
[Integrals.jl](https://docs.sciml.ai/Integrals/stable/) for those.

# Arguments

  - `A`: the interpolation object.
  - `t`: the upper bound, with `first(A.t)` used as the lower bound.
  - `t1`, `t2`: the lower and upper bounds.

# Examples

```jldoctest
using DataInterpolations

u = [1.0, 2.0, 3.0]
t = [1.0, 2.0, 3.0]
A = LinearInterpolation(u, t)

DataInterpolations.integral(A, 1.0, 3.0)

# output

4.0
```
"""
function integral(A::AbstractInterpolation, t::Number)
    return integral(A, first(A.t), t)
end

# For array valued `u` each entry of `A.I` is itself an array, which has no
# `zero`, so build the accumulator from the shape of a data point instead.
_integral_zero(A) = _integral_zero(A.u, eltype(A.I))
_integral_zero(::Any, ::Type{T}) where {T} = zero(T)
function _integral_zero(u::AbstractArray{<:Number}, ::Type{T}) where {T <: AbstractArray}
    return zeros(eltype(T), size(_first(u)))
end
function _integral_zero(
        u::AbstractVector{<:AbstractVector}, ::Type{T}
    ) where {T <: AbstractArray}
    return zeros(eltype(T), size(_first(u)))
end

# `total` is allocated by `_integral_zero` for this call, so array valued
# integrals accumulate into it rather than allocating once per interval.
_add_integral(total::Number, Δ) = total + Δ
_add_integral(total::AbstractArray, Δ) = (total .+= Δ)

_sub_integral(total::Number, Δ) = total - Δ
_sub_integral(total::AbstractArray, Δ) = (total .-= Δ)

function integral(A::AbstractInterpolation, t1::Number, t2::Number)
    !hasfield(typeof(A), :I) && throw(IntegralNotFoundError())

    if t1 == t2
        # If the integration interval is trivial then the result is 0
        return _integral_zero(A)
    elseif t1 > t2
        # Make sure that t1 < t2
        return -integral(A, t2, t1)
    end

    # the index less than or equal to t1
    idx1 = get_idx(A, t1, 0)
    # the index less than t2
    idx2 = get_idx(A, t2, 0; idx_shift = -1, side = :first)

    total = _integral_zero(A)

    # Lower potentially incomplete interval
    if t1 < first(A.t)
        if t2 < first(A.t)
            # If interval is entirely below data
            return _extrapolate_integral_left(A, t1) - _extrapolate_integral_left(A, t2)
        end

        idx1 -= 1 # Make sure lowest complete interval is included
        total = _add_integral(total, _extrapolate_integral_left(A, t1))
    else
        total = _add_integral(total, _integral(A, idx1, t1, A.t[idx1 + 1]))
    end

    # Upper potentially incomplete interval
    if t2 > last(A.t)
        if t1 > last(A.t)
            # If interval is entirely above data
            return _extrapolate_integral_right(A, t2) - _extrapolate_integral_right(A, t1)
        end

        idx2 += 1 # Make sure highest complete interval is included
        total = _add_integral(total, _extrapolate_integral_right(A, t2))
    else
        total = _add_integral(total, _integral(A, idx2, A.t[idx2], t2))
    end

    if idx1 == idx2
        return _integral(A, idx1, t1, t2)
    end

    # Complete intervals
    if A.cache_parameters
        if idx2 > 1
            total = _add_integral(total, A.I[idx2 - 1])
        end
        if idx1 > 0
            total = _sub_integral(total, A.I[idx1])
        end
    else
        for idx in (idx1 + 1):(idx2 - 1)
            total = _add_integral(total, _integral(A, idx, A.t[idx], A.t[idx + 1]))
        end
    end

    return total
end

function _extrapolate_integral_left(A, t)
    (; extrapolation_left) = A
    return if extrapolation_left == ExtrapolationType.None
        throw(LeftExtrapolationError())
    elseif extrapolation_left == ExtrapolationType.Constant
        _first(A.u) * (first(A.t) - t)
    elseif extrapolation_left == ExtrapolationType.Linear
        slope = derivative(A, first(A.t))
        Δt = first(A.t) - t
        u₁ = _first(A.u)
        @. (u₁ - slope * Δt / 2) * Δt
    elseif extrapolation_left == ExtrapolationType.Extension
        _integral(A, 1, t, first(A.t))
    else
        _extrapolate_integral_other_left(A, t, extrapolation_left)
    end
end

function _extrapolate_integral_other_left(A, t, extrapolation_left)
    return if extrapolation_left == ExtrapolationType.Periodic
        t_, n = transformation_periodic(A, t)
        out = -integral(A, t_)
        if !iszero(n)
            out -= n * integral(A, first(A.t), last(A.t))
        end
        out
    elseif extrapolation_left == ExtrapolationType.Reflective
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
    else
        throw(ExtrapolationNotImplementedError())
    end
end

function _extrapolate_integral_right(A, t)
    (; extrapolation_right) = A
    return if extrapolation_right == ExtrapolationType.None
        throw(RightExtrapolationError())
    elseif extrapolation_right == ExtrapolationType.Constant
        _last(A.u) * (t - last(A.t))
    elseif extrapolation_right == ExtrapolationType.Linear
        slope = derivative(A, last(A.t))
        Δt = t - last(A.t)
        uₙ = _last(A.u)
        @. (uₙ + slope * Δt / 2) * Δt
    elseif extrapolation_right == ExtrapolationType.Extension
        _integral(A, length(A.t) - 1, last(A.t), t)
    else
        _extrapolate_integral_other_right(A, t, extrapolation_right)
    end
end

function _extrapolate_integral_other_right(A, t, extrapolation_right)
    return if extrapolation_right == ExtrapolationType.Periodic
        t_, n = transformation_periodic(A, t)
        out = integral(A, first(A.t), t_)
        if !iszero(n)
            out += n * integral(A, first(A.t), last(A.t))
        end
        out
    elseif extrapolation_right == ExtrapolationType.Reflective
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
    else
        throw(ExtrapolationNotImplementedError())
    end
end

# `ConstantInterpolation` extrapolates with the boundary value for `Constant`, `Linear`
# and `Extension` alike (see `_extrapolate_left`/`_extrapolate_right`), so neither the
# generic slope-based `Linear` branch nor the boundary-interval `Extension` branch (which
# picks a `dir`-dependent index) describes what is actually being integrated.
function _extrapolate_integral_left(A::ConstantInterpolation, t)
    (; extrapolation_left) = A
    return if extrapolation_left == ExtrapolationType.None
        throw(LeftExtrapolationError())
    elseif extrapolation_left in (
            ExtrapolationType.Constant, ExtrapolationType.Linear,
            ExtrapolationType.Extension,
        )
        _first(A.u) * (first(A.t) - t)
    else
        _extrapolate_integral_other_left(A, t, extrapolation_left)
    end
end

function _extrapolate_integral_right(A::ConstantInterpolation, t)
    (; extrapolation_right) = A
    return if extrapolation_right == ExtrapolationType.None
        throw(RightExtrapolationError())
    elseif extrapolation_right in (
            ExtrapolationType.Constant, ExtrapolationType.Linear,
            ExtrapolationType.Extension,
        )
        _last(A.u) * (t - last(A.t))
    else
        _extrapolate_integral_other_right(A, t, extrapolation_right)
    end
end

function _extrapolate_integral_right(A::SmoothedConstantInterpolation, t)
    (; extrapolation_right) = A
    return if extrapolation_right == ExtrapolationType.None
        throw(RightExtrapolationError())
    elseif A.extrapolation_right in (
            ExtrapolationType.Constant, ExtrapolationType.Extension,
        )
        n = length(A.t)
        d = min(A.t[end] - A.t[end - 1], 2A.d_max) / 2
        Δt_constant = max(0, t - A.t[end] - d)
        out = Δt_constant * _u_view(A.u, n)

        if !iszero(d)
            c = (_u_view(A.u, n) - _u_view(A.u, n - 1)) / 2
            Δt_transition = min(t - A.t[end], d)
            out = out + Δt_transition * _u_view(A.u, n - 1) -
                c *
                (
                ((Δt_transition / d)^3) / (3 / d) - ((Δt_transition^2) / d) -
                    Δt_transition
            )
        end
        out
    elseif extrapolation_right == ExtrapolationType.Linear
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

"""
    _integral(A::AbstractInterpolation, idx::Integer, t1::Number, t2::Number)

Developer hook used by [`integral`](@ref) to integrate one interval of an interpolation.
Implement this method for a new interpolation subtype that supports analytic integration.
`idx` identifies the interval `[A.t[idx], A.t[idx + 1]]`; the method must return the
integral over the requested subinterval `[t1, t2]`, preserving the shape of one data point.

The subtype must also provide an `I` field containing cached cumulative integrals when
`cache_parameters` is true. If analytic integration is not supported, omit `I` and let the
generic public method report `IntegralNotFoundError`.

User code should call [`integral`](@ref) rather than this developer hook.

# Arguments

- `A`: an [`AbstractInterpolation`](@ref) subtype.
- `idx`: the interval index, with bounds `A.t[idx]` and `A.t[idx + 1]`.
- `t1`, `t2`: the requested interval bounds.

# Examples

```julia
DataInterpolations._integral(A::MyInterpolation, idx::Integer, t1::Number, t2::Number) =
    2 * (t2 - t1)
```
"""
function _integral(
        A::LinearInterpolation{
            <:Union{AbstractArray{<:Number}, AbstractVector{<:AbstractVector}},
        },
        idx::Number, t1::Number, t2::Number
    )
    slope = get_parameters(A, idx)
    uᵢ = _u_view(A.u, idx)
    Δt_mean = (t1 + t2) / 2 - A.t[idx]
    Δt = t2 - t1
    return @. (uᵢ + slope * Δt_mean) * Δt
end

function _integral(
        A::ConstantInterpolation{
            <:Union{AbstractArray{<:Number}, AbstractVector{<:AbstractVector}},
        },
        idx::Number, t1::Number, t2::Number
    )
    # :left/:right means that the value to the left/right is used for interpolation
    uᵢ = _u_view(A.u, A.dir === :left ? idx : idx + 1)
    return uᵢ * (t2 - t1)
end

function _integral(
        A::SmoothedConstantInterpolation,
        idx::Number, t1::Number, t2::Number
    )
    d_lower, d_upper, c_lower, c_upper = get_parameters(A, idx)

    bound_lower = A.t[idx] + d_lower
    bound_upper = A.t[idx + 1] - d_upper

    out = _u_view(A.u, idx) * (t2 - t1)

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
            (
            (1 - (A.t[idx + 1] - t2) / d_upper)^3 -
                (1 - (A.t[idx + 1] - t1_) / d_upper)^3
        ) / 3
    end

    return out
end

function _integral(
        A::QuadraticInterpolation{
            <:Union{AbstractArray{<:Number}, AbstractVector{<:AbstractVector}},
        },
        idx::Number, t1::Number, t2::Number
    )
    α, β = get_parameters(A, idx)
    uᵢ = _u_view(A.u, idx)
    tᵢ = A.t[idx]
    t1_rel = t1 - tᵢ
    t2_rel = t2 - tᵢ
    Δt = t2 - t1
    return @. Δt *
        (α * (t2_rel^2 + t1_rel * t2_rel + t1_rel^2) / 3 + β * (t2_rel + t1_rel) / 2 + uᵢ)
end

function _integral(
        A::QuadraticSpline, idx::Number, t1::Number, t2::Number
    )
    α, β = get_parameters(A, idx)
    uᵢ = _u_view(A.u, idx)
    tᵢ = A.t[idx]
    t1_rel = t1 - tᵢ
    t2_rel = t2 - tᵢ
    Δt = t2 - t1
    return @. Δt *
        (α * (t2_rel^2 + t1_rel * t2_rel + t1_rel^2) / 3 + β * (t2_rel + t1_rel) / 2 + uᵢ)
end

function _integral(
        A::CubicSpline{<:AbstractVector}, idx::Number, t1::Number, t2::Number
    )
    tᵢ = A.t[idx]
    tᵢ₊₁ = A.t[idx + 1]
    c₁, c₂ = get_parameters(A, idx)
    zero_ = zero(c₁)
    return integrate_cubic_polynomial(
        t1, t2, tᵢ, zero_, c₁, zero_, A.z[idx + 1] / (6A.h[idx + 1])
    ) +
        integrate_cubic_polynomial(
        t1, t2, tᵢ₊₁, zero_, -c₂, zero_, -A.z[idx] / (6A.h[idx + 1])
    )
end

function _integral(
        A::CubicSpline{<:AbstractArray}, idx::Number, t1::Number, t2::Number
    )
    tᵢ = A.t[idx]
    tᵢ₊₁ = A.t[idx + 1]
    c₁, c₂ = get_parameters(A, idx)
    zᵢ = _u_view(A.z, idx)
    zᵢ₊₁ = _u_view(A.z, idx + 1)
    zero_ = zero(c₁)
    return integrate_cubic_polynomial(
        t1, t2, tᵢ, zero_, c₁, zero_, zᵢ₊₁ / (6A.h[idx + 1])
    ) +
        integrate_cubic_polynomial(
        t1, t2, tᵢ₊₁, zero_, -c₂, zero_, -zᵢ / (6A.h[idx + 1])
    )
end

function _integral(
        A::AkimaInterpolation{<:AbstractVector{<:Number}},
        idx::Number, t1::Number, t2::Number
    )
    return integrate_cubic_polynomial(t1, t2, A.t[idx], A.u[idx], A.b[idx], A.c[idx], A.d[idx])
end

function _integral(
        A::AkimaInterpolation{<:AbstractArray},
        idx::Number, t1::Number, t2::Number
    )
    uᵢ = _u_view(A.u, idx)
    bᵢ = _u_view(A.b, idx)
    cᵢ = _u_view(A.c, idx)
    dᵢ = _u_view(A.d, idx)
    return integrate_cubic_polynomial(t1, t2, A.t[idx], uᵢ, bᵢ, cᵢ, dᵢ)
end

function _integral(A::LagrangeInterpolation, idx::Number, t1::Number, t2::Number)
    throw(IntegralNotFoundError())
end
# Evaluate the antiderivative of a B-spline at a point.
# The antiderivative of a degree-d B-spline is a degree-(d+1) B-spline with
# extended knot vector and coefficients derived from the original.
function _bspline_antiderivative_val(c, d, k, t_eval)
    nc = length(c)
    dp1 = d + 1
    # Antiderivative coefficients: C[1] = 0, C[i+1] = C[i] + c[i]*(k[i+d+1]-k[i])/(d+1)
    T = promote_type(eltype(c), eltype(k), typeof(t_eval))
    C = zeros(T, nc + 1)
    for i in 1:nc
        C[i + 1] = C[i] + c[i] * (k[i + dp1] - k[i]) / dp1
    end
    # Extended knot vector: prepend k[1], append k[end]
    nk = length(k)
    k_ext = zeros(eltype(k), nk + 2)
    k_ext[1] = k[1]
    for i in 1:nk
        k_ext[i + 1] = k[i]
    end
    k_ext[end] = k[end]
    # Evaluate degree-(d+1) B-spline with coefficients C on k_ext
    sc = zeros(T, nc + 1)
    nonzero = spline_coefficients!(sc, dp1, k_ext, t_eval)
    result = zero(T)
    for i in nonzero
        result += sc[i] * C[i]
    end
    return result
end

function _integral(
        A::BSplineInterpolation{<:AbstractVector{<:Number}}, idx::Number, t1::Number, t2::Number
    )
    return _bspline_antiderivative_val(A.c, A.d, A.k, t2) -
        _bspline_antiderivative_val(A.c, A.d, A.k, t1)
end
function _integral(A::BSplineApprox{<:AbstractVector{<:Number}}, idx::Number, t1::Number, t2::Number)
    return _bspline_antiderivative_val(A.c, A.d, A.k, t2) -
        _bspline_antiderivative_val(A.c, A.d, A.k, t1)
end

# Same antiderivative recursion as the scalar `_bspline_antiderivative_val` above, but
# accumulating array/vector-of-vector-valued control points instead of `Number`s.
function _bspline_antiderivative_val(c::AbstractArray, d, k, t_eval)
    nc = size(c)[end]
    dp1 = d + 1
    c₁ = _u_view(c, 1)
    C = [zero(c₁) for _ in 1:(nc + 1)]
    for i in 1:nc
        C[i + 1] = C[i] + _u_view(c, i) * (k[i + dp1] - k[i]) / dp1
    end
    nk = length(k)
    k_ext = zeros(eltype(k), nk + 2)
    k_ext[1] = k[1]
    for i in 1:nk
        k_ext[i + 1] = k[i]
    end
    k_ext[end] = k[end]
    sc = zeros(eltype(k), nc + 1)
    nonzero = spline_coefficients!(sc, dp1, k_ext, t_eval)
    result = zero(c₁)
    for i in nonzero
        result = result + sc[i] * C[i]
    end
    return result
end

function _integral(
        A::BSplineInterpolation{
            <:Union{AbstractArray{<:Number}, AbstractVector{<:AbstractVector}},
        },
        idx::Number, t1::Number, t2::Number
    )
    return _bspline_antiderivative_val(A.c, A.d, A.k, t2) -
        _bspline_antiderivative_val(A.c, A.d, A.k, t1)
end
function _integral(
        A::BSplineApprox{
            <:Union{AbstractArray{<:Number}, AbstractVector{<:AbstractVector}},
        },
        idx::Number, t1::Number, t2::Number
    )
    return _bspline_antiderivative_val(A.c, A.d, A.k, t2) -
        _bspline_antiderivative_val(A.c, A.d, A.k, t1)
end

# Override integral to bypass the hasfield(:I) check in the generic method.
# The antiderivative is computed on the fly, so no cached I field is needed.
const _BSplineTypes = Union{
    BSplineInterpolation{<:Union{AbstractArray{<:Number}, AbstractVector{<:AbstractVector}}},
    BSplineApprox{<:Union{AbstractArray{<:Number}, AbstractVector{<:AbstractVector}}},
}

function integral(A::_BSplineTypes, t::Number)
    return integral(A, first(A.t), t)
end

function integral(A::_BSplineTypes, t1::Number, t2::Number)
    t1 == t2 && return zero(_first(A.u))
    t1 > t2 && return -integral(A, t2, t1)

    total = zero(_first(A.u))
    lo = t1
    hi = t2

    if lo < first(A.t)
        if hi <= first(A.t)
            return _extrapolate_integral_left(A, lo) - _extrapolate_integral_left(A, hi)
        end
        total = total + _extrapolate_integral_left(A, lo)
        lo = first(A.t)
    end

    if hi > last(A.t)
        if lo >= last(A.t)
            return _extrapolate_integral_right(A, hi) - _extrapolate_integral_right(A, lo)
        end
        total = total + _extrapolate_integral_right(A, hi)
        hi = last(A.t)
    end

    total = total + _bspline_antiderivative_val(A.c, A.d, A.k, hi) -
        _bspline_antiderivative_val(A.c, A.d, A.k, lo)

    return total
end

# Cubic Hermite Spline
function _integral(
        A::CubicHermiteSpline{<:AbstractVector{<:Number}}, idx::Number, t1::Number, t2::Number
    )
    c₁, c₂ = get_parameters(A, idx)
    tᵢ = A.t[idx]
    tᵢ₊₁ = A.t[idx + 1]
    c = c₁ - c₂ * (tᵢ₊₁ - tᵢ)
    return integrate_cubic_polynomial(t1, t2, tᵢ, A.u[idx], A.du[idx], c, c₂)
end

function _integral(
        A::CubicHermiteSpline{<:AbstractArray}, idx::Number, t1::Number, t2::Number
    )
    c₁, c₂ = get_parameters(A, idx)
    tᵢ = A.t[idx]
    tᵢ₊₁ = A.t[idx + 1]
    c = c₁ - c₂ * (tᵢ₊₁ - tᵢ)
    uᵢ = _u_view(A.u, idx)
    duᵢ = _u_view(A.du, idx)
    return integrate_cubic_polynomial(t1, t2, tᵢ, uᵢ, duᵢ, c, c₂)
end

# Quintic Hermite Spline
function _integral(
        A::QuinticHermiteSpline{<:AbstractVector{<:Number}}, idx::Number, t1::Number, t2::Number
    )
    tᵢ = A.t[idx]
    tᵢ₊₁ = A.t[idx + 1]
    Δt = tᵢ₊₁ - tᵢ
    c₁, c₂, c₃ = get_parameters(A, idx)
    return integrate_quintic_polynomial(
        t1, t2, tᵢ, A.u[idx], A.du[idx], A.ddu[idx] / 2,
        c₁ + Δt * (-c₂ + c₃ * Δt), c₂ - 2c₃ * Δt, c₃
    )
end

function _integral(
        A::QuinticHermiteSpline{<:AbstractArray}, idx::Number, t1::Number, t2::Number
    )
    tᵢ = A.t[idx]
    tᵢ₊₁ = A.t[idx + 1]
    Δt = tᵢ₊₁ - tᵢ
    c₁, c₂, c₃ = get_parameters(A, idx)
    uᵢ = _u_view(A.u, idx)
    duᵢ = _u_view(A.du, idx)
    dduᵢ = _u_view(A.ddu, idx)
    return integrate_quintic_polynomial(
        t1, t2, tᵢ, uᵢ, duᵢ, dduᵢ / 2,
        c₁ + Δt * (-c₂ + c₃ * Δt), c₂ - 2c₃ * Δt, c₃
    )
end

# Each `A.t` interval of a `SmoothArcLengthInterpolation` is a circle segment followed by
# a line segment (or vice versa, per `short_side_left[idx]`); both are analytically
# integrable in the arc-length parameter. `_interpolate`/`_derivative` for this type pick
# a branch via a single threshold `τ_break` with no clamp on the other side (so extrapolation
# beyond `[0, total]` continues whichever branch governs at the boundary, e.g. the circle
# keeps rotating rather than the piece vanishing) — mirror that here: split `[τ1, τ2]` at
# `τ_break` only, without an extra clamp to `[0, total]`, so `Extension` extrapolation
# (which reuses this same `_integral` outside the segment's natural range) matches.
function _integral(A::SmoothArcLengthInterpolation, idx::Number, t1::Number, t2::Number)
    Δt_circle_seg = A.Δt_circle_segment[idx]
    Δt_line_seg = A.Δt_line_segment[idx]
    short_side_left = A.short_side_left[idx]
    t0 = A.t[idx]
    τ1 = t1 - t0
    τ2 = t2 - t0

    Rⱼ = A.radius[idx]
    c = view(A.center, :, idx)
    v₁ = view(A.dir_1, :, idx)
    v₂ = view(A.dir_2, :, idx)

    τ_break = short_side_left ? Δt_circle_seg : Δt_line_seg

    if short_side_left
        circle_offset = zero(τ1)
        line_offset = Δt_circle_seg + Δt_line_seg
        u_line = view(A.u, :, idx + 1)
        d_line = view(A.d, :, idx + 1)
        # circle governs τ < τ_break, line governs τ >= τ_break
        circ_a, circ_b = τ1, min(τ2, τ_break)
        line_a, line_b = max(τ1, τ_break), τ2
    else
        circle_offset = Δt_line_seg
        line_offset = zero(τ1)
        u_line = view(A.u, :, idx)
        d_line = view(A.d, :, idx)
        # line governs τ < τ_break, circle governs τ >= τ_break
        line_a, line_b = τ1, min(τ2, τ_break)
        circ_a, circ_b = max(τ1, τ_break), τ2
    end

    # Clamp each piece's own bounds so it contributes zero when it doesn't overlap [τ1, τ2].
    line_b = max(line_a, line_b)
    circ_b = max(circ_a, circ_b)

    # Line piece: out(τ) = u_line + (τ - line_offset) * d_line
    Δt_line = line_b - line_a
    Δt_mean_line = (line_a + line_b) / 2 - line_offset
    line_contribution = @. (u_line + d_line * Δt_mean_line) * Δt_line

    # Circle piece: out(τ) = c + cos((τ - circle_offset) / Rⱼ) * v₁ + sin((τ - circle_offset) / Rⱼ) * v₂
    τ_a = circ_a - circle_offset
    τ_b = circ_b - circle_offset
    Sa, Ca = sincos(τ_a / Rⱼ)
    Sb, Cb = sincos(τ_b / Rⱼ)
    circle_contribution = @. c * (τ_b - τ_a) + Rⱼ * (Sb - Sa) * v₁ - Rⱼ * (Cb - Ca) * v₂

    return line_contribution .+ circle_contribution
end
