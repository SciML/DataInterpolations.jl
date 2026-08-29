"""
    derivative(A::AbstractInterpolation, t::Number, order = 1)

Evaluate the `order`-th derivative of the interpolation `A` at `t`.

First order derivatives are computed analytically; second order derivatives are computed
from the analytical first order derivative with ForwardDiff.jl. The derivative is taken as
a left derivative, so at a data point where the derivative jumps the value from the
interval to the left of the point is returned.

Outside of the range of `A.t` the derivative of the extrapolation is returned, as
determined by the `extrapolation_left` and `extrapolation_right` settings of `A`.

# Arguments

  - `A`: the interpolation object.
  - `t`: the point at which to evaluate the derivative.
  - `order`: the order of the derivative, either `1` or `2`. Defaults to `1`.

# Examples

```jldoctest
using DataInterpolations

u = [1.0, 4.0, 9.0]
t = [1.0, 2.0, 3.0]
A = QuadraticInterpolation(u, t)

DataInterpolations.derivative(A, 2.5)

# output

5.0
```
"""
function derivative(A, t, order = 1)
    (order ∉ (1, 2)) && throw(DerivativeNotFoundError())
    if t < first(A.t)
        _extrapolate_derivative_left(A, t, order)
    elseif t > last(A.t)
        _extrapolate_derivative_right(A, t, order)
    else
        iguess = A.iguesser
        if order == 1
            return _derivative(A, t, iguess)
        end
        return ForwardDiff.derivative(
            t -> begin
                -_derivative(A, -t, iguess)
            end, -t
        ) # take derivative backwards in t to make it a left rather than right derivative
    end
end

# Zero of the derivative's value type, shaped like one data point.
function _zero_derivative(A)
    u₁ = _first(A.u)
    t₁ = one(A.t[1])
    return @. zero(u₁ / t₁)
end

function _extrapolate_derivative_left(A, t, order)
    (; extrapolation_left) = A
    return if extrapolation_left == ExtrapolationType.None
        throw(LeftExtrapolationError())
    elseif extrapolation_left == ExtrapolationType.Constant
        _zero_derivative(A)
    elseif extrapolation_left == ExtrapolationType.Linear
        (order == 1) ? _derivative(A, first(A.t), 1) : _zero_derivative(A)
    elseif extrapolation_left == ExtrapolationType.Extension
        (order == 1) ? _derivative(A, t, length(A.t)) :
            ForwardDiff.derivative(
                t -> begin
                    _derivative(A, t, length(A.t))
                end, t
            )
    elseif extrapolation_left == ExtrapolationType.Periodic
        t_, _ = transformation_periodic(A, t)
        (order == 1) ? _derivative(A, t_, A.iguesser) :
            ForwardDiff.derivative(
                t -> begin
                    _derivative(A, t, A.iguesser)
                end, t_
            )
    else
        # extrapolation_left == ExtrapolationType.Reflective
        t_, n = transformation_reflective(A, t)
        sign = isodd(n) ? -1 : 1
        (order == 1) ? sign * _derivative(A, t_, A.iguesser) :
            ForwardDiff.derivative(
                t -> begin
                    sign * _derivative(A, t, A.iguesser)
                end, t_
            )
    end
end

function _extrapolate_derivative_right(A, t, order)
    (; extrapolation_right) = A
    return if extrapolation_right == ExtrapolationType.None
        throw(RightExtrapolationError())
    elseif extrapolation_right == ExtrapolationType.Constant
        _zero_derivative(A)
    elseif extrapolation_right == ExtrapolationType.Linear
        (order == 1) ? _derivative(A, last(A.t), length(A.t)) : _zero_derivative(A)
    elseif extrapolation_right == ExtrapolationType.Extension
        (order == 1) ? _derivative(A, t, length(A.t)) :
            ForwardDiff.derivative(
                t -> begin
                    _derivative(A, t, length(A.t))
                end, t
            )
    elseif extrapolation_right == ExtrapolationType.Periodic
        t_, _ = transformation_periodic(A, t)
        (order == 1) ? _derivative(A, t_, A.iguesser) :
            ForwardDiff.derivative(
                t -> begin
                    _derivative(A, t, A.iguesser)
                end, t_
            )
    else
        # extrapolation_right == ExtrapolationType.Reflective
        t_, n = transformation_reflective(A, t)
        sign = iseven(n) ? -1 : 1
        (order == 1) ? sign * _derivative(A, t_, A.iguesser) :
            ForwardDiff.derivative(
                t -> begin
                    sign * _derivative(A, t, A.iguesser)
                end, t_
            )
    end
end

function _extrapolate_derivative_right(A::SmoothedConstantInterpolation, t, order)
    return if A.extrapolation_right == ExtrapolationType.None
        throw(RightExtrapolationError())
    elseif A.extrapolation_right in (
            ExtrapolationType.Constant, ExtrapolationType.Extension,
        )
        d = min(A.t[end] - A.t[end - 1], 2A.d_max) / 2
        if A.t[end] + d < t
            zero(_first(A.u))
        else
            n = length(A.t)
            c = (_u_view(A.u, n) - _u_view(A.u, n - 1)) / 2
            -c * (2((t - A.t[end]) / d) - 2) / d
        end

    else
        _extrapolate_other(A, t, A.extrapolation_right)
    end
end

"""
    _derivative(A::AbstractInterpolation, t::Number, iguess)

Developer hook used by [`derivative`](@ref) for an in-range scalar evaluation. Implement
this method for a new [`AbstractInterpolation`](@ref) subtype. The returned value must
have the same shape as one interpolated data point, and the value at a knot must be the
left derivative because the public derivative interface is left-continuous at knots.

`iguess` is the reusable search hint stored by the interpolation. Implementations should
use it when locating the interval and must not replace it with a new mutable search state
for every call.

User code should call [`derivative`](@ref) rather than this developer hook.

# Arguments

- `A`: an [`AbstractInterpolation`](@ref) subtype.
- `t`: an in-range scalar evaluation point.
- `iguess`: the reusable interval-search state from `A.iguesser`.

# Examples

```julia
DataInterpolations._derivative(A::MyInterpolation, t::Number, iguess) = 2.0
```
"""
function _derivative(A::LinearInterpolation, t::Number, iguess)
    idx = get_idx(A, t, iguess; idx_shift = -1, ub_shift = -1, side = :first)
    slope = get_parameters(A, idx)
    return slope
end

function _derivative(A::SmoothedConstantInterpolation, t::Number, iguess)
    idx = get_idx(A, t, iguess)
    d_lower, d_upper, c_lower, c_upper = get_parameters(A, idx)

    return if (t - A.t[idx]) < d_lower
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
    return @. 2α * Δt + β
end

# Differentiation-matrix formula at node `A.t[k]`, generic over `values`: with `A.u` it
# gives the first derivative; with the vector of first derivatives it gives the second
# (D² = D·D, the standard spectral-differentiation trick).
function _lagrange_node_derivative(A::LagrangeInterpolation{<:AbstractVector}, k, values)
    der = zero(values[k])
    invw_k = inv(A.p.w[k])
    for j in eachindex(A.t)
        j == k && continue
        coef = (A.p.w[j] * invw_k) / (A.t[k] - A.t[j])
        der += coef * (values[j] - values[k])
    end
    return der
end

function _lagrange_node_derivative(
        A::LagrangeInterpolation{<:AbstractMatrix}, k, values::AbstractMatrix
    )
    der = zero(values[:, k])
    invw_k = inv(A.p.w[k])
    @views for j in eachindex(A.t)
        j == k && continue
        coef = (A.p.w[j] * invw_k) / (A.t[k] - A.t[j])
        der = der + coef * (values[:, j] - values[:, k])
    end
    return der
end

function _lagrange_second_derivative_at_node(A::LagrangeInterpolation{<:AbstractVector}, k)
    du = [_lagrange_node_derivative(A, j, A.u) for j in eachindex(A.t)]
    return _lagrange_node_derivative(A, k, du)
end

function _lagrange_second_derivative_at_node(A::LagrangeInterpolation{<:AbstractMatrix}, k)
    du = stack(_lagrange_node_derivative(A, j, A.u) for j in eachindex(A.t))
    return _lagrange_node_derivative(A, k, du)
end

# p(t) = N(t)/D(t), so p'(t) = (N'D - ND')/D^2; off-node only, since N and D individually
# diverge at a node (handled above via `_lagrange_node_derivative`).
function _derivative(A::LagrangeInterpolation{<:AbstractVector}, t::Number)
    idx = _searchsortedfirst(A.t, t)
    !isnothing(idx) && return _lagrange_node_derivative(A, idx, A.u)
    N = zero(A.p.wu[1])
    D = zero(A.t[1])
    N′ = zero(A.p.wu[1])
    D′ = zero(A.t[1])
    for i in eachindex(A.t)
        invti = inv(t - A.t[i])
        wi_inv = A.p.w[i] * invti
        wui_inv = A.p.wu[i] * invti
        D += wi_inv
        N += wui_inv
        D′ -= wi_inv * invti
        N′ -= wui_inv * invti
    end
    return (N′ * D - N * D′) / D^2
end

function _derivative(A::LagrangeInterpolation{<:AbstractMatrix}, t::Number)
    idx = _searchsortedfirst(A.t, t)
    !isnothing(idx) && return _lagrange_node_derivative(A, idx, A.u)
    N = zero(A.p.wu[:, 1])
    D = zero(A.t[1])
    N′ = zero(A.p.wu[:, 1])
    D′ = zero(A.t[1])
    @views for i in eachindex(A.t)
        invti = inv(t - A.t[i])
        wi_inv = A.p.w[i] * invti
        D += wi_inv
        D′ -= wi_inv * invti
        N = N + A.p.wu[:, i] * invti
        N′ = N′ - A.p.wu[:, i] * invti^2
    end
    return @. (N′ * D - N * D′) / D^2
end

function _derivative(A::LagrangeInterpolation{<:AbstractVector}, t::Number, idx)
    return _derivative(A, t)
end
function _derivative(A::LagrangeInterpolation{<:AbstractMatrix}, t::Number, idx)
    return _derivative(A, t)
end

# Autodiffing `_derivative` again (the generic `order == 2` path) hits a 0/0 at nodes,
# since ForwardDiff's `Dual` equality compares partials too so the node check misses
# under nested AD. Use the exact second differentiation matrix there instead; off-node,
# the generic nested-autodiff path already works.
function derivative(A::LagrangeInterpolation, t::Number, order::Int = 1)
    (order ∉ (1, 2)) && throw(DerivativeNotFoundError())
    if t < first(A.t)
        return _extrapolate_derivative_left(A, t, order)
    elseif t > last(A.t)
        return _extrapolate_derivative_right(A, t, order)
    end
    iguess = A.iguesser
    order == 1 && return _derivative(A, t, iguess)
    idx = _searchsortedfirst(A.t, t)
    !isnothing(idx) && return _lagrange_second_derivative_at_node(A, idx)
    return ForwardDiff.derivative(τ -> -_derivative(A, -τ, iguess), -t)
end

function _derivative(A::AkimaInterpolation{<:AbstractVector{<:Number}}, t::Number, iguess)
    idx = get_idx(A, t, iguess; idx_shift = -1, side = :first)
    j = min(idx, length(A.c))  # for smooth derivative at A.t[end]
    wj = t - A.t[idx]
    return @evalpoly wj A.b[idx] 2A.c[j] 3A.d[j]
end

function _derivative(A::AkimaInterpolation{<:AbstractArray}, t::Number, iguess)
    idx = get_idx(A, t, iguess; idx_shift = -1, side = :first)
    j = min(idx, length(A.t) - 1)  # for smooth derivative at A.t[end]
    wj = t - A.t[idx]
    bᵢ = _u_view(A.b, idx)
    cⱼ = _u_view(A.c, j)
    dⱼ = _u_view(A.d, j)
    return @. @evalpoly wj bᵢ 2cⱼ 3dⱼ
end

function _derivative(A::ConstantInterpolation, t::Number, iguess)
    return zero(_first(A.u))
end

function _derivative(A::ConstantInterpolation{<:AbstractVector}, t::Number, iguess)
    return isempty(searchsorted(A.t, t)) ? zero(A.u[1]) : typed_nan(A.u)
end

function _derivative(A::ConstantInterpolation{<:AbstractMatrix}, t::Number, iguess)
    u₁ = _first(A.u)
    return isempty(searchsorted(A.t, t)) ? zero(u₁) : fill(typed_nan(A.u), size(u₁))
end

function _derivative(
        A::ConstantInterpolation{<:AbstractVector{<:AbstractVector}}, t::Number, iguess
    )
    u₁ = _first(A.u)
    return isempty(searchsorted(A.t, t)) ? zero(u₁) : fill(typed_nan(u₁), size(u₁))
end

# QuadraticSpline Interpolation
function _derivative(A::QuadraticSpline, t::Number, iguess)
    idx = get_idx(A, t, iguess)
    α, β = get_parameters(A, idx)
    Δt = t - A.t[idx]
    Δt_full = A.t[idx + 1] - A.t[idx]
    return 2α * Δt / Δt_full^2 + β / Δt_full
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
    return dI + dC + dD
end

function _derivative(A::CubicSpline{<:AbstractArray}, t::Number, iguess)
    idx = get_idx(A, t, iguess)
    Δt₁ = t - A.t[idx]
    Δt₂ = A.t[idx + 1] - t
    ax = axes(A.z)[1:(end - 1)]
    dI = (-A.z[ax..., idx] * Δt₂^2 + A.z[ax..., idx + 1] * Δt₁^2) / (2A.h[idx + 1])
    c₁, c₂ = get_parameters(A, idx)
    dC = c₁
    dD = -c₂
    return dI + dC + dD
end

function _derivative(A::BSplineInterpolation{<:AbstractVector{<:Number}}, t::Number, iguess)
    if A.d == 0
        return isempty(searchsorted(A.t, t)) ? zero(A.u[1]) : typed_nan(A.u)
    end
    n = length(A.t)
    # Stack-allocated basis window (see `_interpolate`): must be reentrant, #532.
    vals, offset, m = bspline_nonzero_coefficients(A.d - 1, A.k, t, n)
    ducum = zero(eltype(A.u))
    if t == A.t[1]
        denom = A.k[A.d + 2] - A.k[2]
        ducum = denom != 0 ? (A.c[2] - A.c[1]) / denom : zero(eltype(A.u))
    else
        @inbounds for i in 1:(n - 1)
            denom = A.k[i + A.d + 1] - A.k[i + 1]
            l = i + 1 - offset
            if denom != 0 && 1 <= l <= m
                ducum += vals[l] * (A.c[i + 1] - A.c[i]) / denom
            end
        end
    end
    return ducum * A.d
end

function _derivative(
        A::BSplineInterpolation{<:AbstractVector{<:AbstractVector}}, t::Number, iguess
    )
    if A.d == 0
        return isempty(searchsorted(A.t, t)) ? zero(A.u[1]) :
            fill(typed_nan(A.u[1]), size(A.u[1]))
    end
    n = length(A.t)
    # Stack-allocated basis window (see `_interpolate`): must be reentrant, #532.
    vals, offset, m = bspline_nonzero_coefficients(A.d - 1, A.k, t, n)
    ducum = zero(A.u[1])
    if t == A.t[1]
        denom = A.k[A.d + 2] - A.k[2]
        ducum = denom != 0 ? (A.c[2] - A.c[1]) / denom : zero(A.u[1])
    else
        @inbounds for i in 1:(n - 1)
            denom = A.k[i + A.d + 1] - A.k[i + 1]
            l = i + 1 - offset
            if denom != 0 && 1 <= l <= m
                ducum = ducum + vals[l] * (A.c[i + 1] - A.c[i]) / denom
            end
        end
    end
    return ducum * A.d
end

function _derivative(
        A::BSplineInterpolation{<:AbstractArray{<:Number}}, t::Number, iguess
    )
    if A.d == 0
        return isempty(searchsorted(A.t, t)) ? zero(A.u[:, 1]) :
            typed_nan(A.u) .* A.u[:, 1]
    end
    ax_u = axes(A.u)[1:(end - 1)]
    n = length(A.t)
    # Stack-allocated basis window (see `_interpolate`): must be reentrant, #532.
    vals, offset, m = bspline_nonzero_coefficients(A.d - 1, A.k, t, n)
    ducum = zeros(size(A.u)[1:(end - 1)]...)
    if t == A.t[1]
        denom = A.k[A.d + 2] - A.k[2]
        if denom != 0
            ducum = (A.c[ax_u..., 2] - A.c[ax_u..., 1]) / denom
        end
    else
        @inbounds for i in 1:(n - 1)
            denom = A.k[i + A.d + 1] - A.k[i + 1]
            l = i + 1 - offset
            if denom != 0 && 1 <= l <= m
                ducum = ducum +
                    vals[l] * (A.c[ax_u..., i + 1] - A.c[ax_u..., i]) / denom
            end
        end
    end
    return ducum * A.d
end
# BSpline Curve Approx
function _derivative(A::BSplineApprox{<:AbstractVector{<:Number}}, t::Number, iguess)
    if A.d == 0
        return isempty(searchsorted(A.t, t)) ? zero(A.u[1]) : typed_nan(A.u)
    end
    # Stack-allocated basis window (see `_interpolate`): must be reentrant, #532.
    vals, offset, m = bspline_nonzero_coefficients(A.d - 1, A.k, t, A.h)
    ducum = zero(eltype(A.u))
    if t == A.t[1]
        denom = A.k[A.d + 2] - A.k[2]
        ducum = denom != 0 ? (A.c[2] - A.c[1]) / denom : zero(eltype(A.u))
    else
        @inbounds for i in 1:(A.h - 1)
            denom = A.k[i + A.d + 1] - A.k[i + 1]
            l = i + 1 - offset
            if denom != 0 && 1 <= l <= m
                ducum += vals[l] * (A.c[i + 1] - A.c[i]) / denom
            end
        end
    end
    return ducum * A.d
end

function _derivative(
        A::BSplineApprox{<:AbstractVector{<:AbstractVector}}, t::Number, iguess
    )
    if A.d == 0
        return isempty(searchsorted(A.t, t)) ? zero(A.u[1]) :
            fill(typed_nan(A.u[1]), size(A.u[1]))
    end
    # Stack-allocated basis window (see `_interpolate`): must be reentrant, #532.
    vals, offset, m = bspline_nonzero_coefficients(A.d - 1, A.k, t, A.h)
    ducum = zero(A.u[1])
    if t == A.t[1]
        denom = A.k[A.d + 2] - A.k[2]
        ducum = denom != 0 ? (A.c[2] - A.c[1]) / denom : zero(A.u[1])
    else
        @inbounds for i in 1:(A.h - 1)
            denom = A.k[i + A.d + 1] - A.k[i + 1]
            l = i + 1 - offset
            if denom != 0 && 1 <= l <= m
                ducum = ducum + vals[l] * (A.c[i + 1] - A.c[i]) / denom
            end
        end
    end
    return ducum * A.d
end

function _derivative(
        A::BSplineApprox{<:AbstractArray{<:Number}}, t::Number, iguess
    )
    if A.d == 0
        return isempty(searchsorted(A.t, t)) ? zero(A.u[:, 1]) :
            typed_nan(A.u) .* A.u[:, 1]
    end
    ax_u = axes(A.u)[1:(end - 1)]
    # Stack-allocated basis window (see `_interpolate`): must be reentrant, #532.
    vals, offset, m = bspline_nonzero_coefficients(A.d - 1, A.k, t, A.h)
    ducum = zeros(size(A.u)[1:(end - 1)]...)
    if t == A.t[1]
        denom = A.k[A.d + 2] - A.k[2]
        if denom != 0
            ducum = (A.c[ax_u..., 2] - A.c[ax_u..., 1]) / denom
        end
    else
        @inbounds for i in 1:(A.h - 1)
            denom = A.k[i + A.d + 1] - A.k[i + 1]
            l = i + 1 - offset
            if denom != 0 && 1 <= l <= m
                ducum += vals[l] * (A.c[ax_u..., i + 1] - A.c[ax_u..., i]) / denom
            end
        end
    end
    return ducum * A.d
end
# Cubic Hermite Spline
function _derivative(
        A::CubicHermiteSpline{<:AbstractVector{<:Number}}, t::Number, iguess
    )
    idx = get_idx(A, t, iguess)
    Δt₀ = t - A.t[idx]
    Δt₁ = t - A.t[idx + 1]
    out = A.du[idx]
    c₁, c₂ = get_parameters(A, idx)
    out += Δt₀ * (Δt₀ * c₂ + 2(c₁ + Δt₁ * c₂))
    return out
end

function _derivative(
        A::CubicHermiteSpline{<:AbstractArray}, t::Number, iguess
    )
    idx = get_idx(A, t, iguess)
    Δt₀ = t - A.t[idx]
    Δt₁ = t - A.t[idx + 1]
    out = _u_view(A.du, idx)
    c₁, c₂ = get_parameters(A, idx)
    out = out + Δt₀ * (Δt₀ * c₂ + 2(c₁ + Δt₁ * c₂))
    return out
end

# Quintic Hermite Spline
function _derivative(
        A::QuinticHermiteSpline{<:AbstractVector{<:Number}}, t::Number, iguess
    )
    idx = get_idx(A, t, iguess)
    Δt₀ = t - A.t[idx]
    Δt₁ = t - A.t[idx + 1]
    out = A.du[idx] + A.ddu[idx] * Δt₀
    c₁, c₂, c₃ = get_parameters(A, idx)
    out += Δt₀^2 *
        (3c₁ + (3Δt₁ + Δt₀) * c₂ + (3Δt₁^2 + Δt₀ * 2Δt₁) * c₃)
    return out
end

function _derivative(
        A::QuinticHermiteSpline{<:AbstractArray}, t::Number, iguess
    )
    idx = get_idx(A, t, iguess)
    Δt₀ = t - A.t[idx]
    Δt₁ = t - A.t[idx + 1]
    out = _u_view(A.du, idx) + _u_view(A.ddu, idx) * Δt₀
    c₁, c₂, c₃ = get_parameters(A, idx)
    out = out + Δt₀^2 *
        (3c₁ + (3Δt₁ + Δt₀) * c₂ + (3Δt₁^2 + Δt₀ * 2Δt₁) * c₃)
    return out
end

function _derivative(A::SmoothArcLengthInterpolation, t::Number, iguess)
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

    # See `_interpolate` above: return per branch rather than mutate a shared buffer, so
    # this stays differentiable with reverse-mode AD. `+ zero(Δt)` on the two constant
    # branches (which don't otherwise involve `t`) keeps the element type consistent with
    # `typeof(t)` across all three branches, matching what the old `Vector{typeof(t)}`
    # buffer + broadcast-assignment did for e.g. ForwardDiff `Dual` inputs.
    return if in_circle_arc
        t_circle_seg = short_side_left ? Δt : Δt - Δt_line_seg
        Rⱼ = A.radius[idx]
        S, C = sincos(t_circle_seg / Rⱼ)
        v₁ = view(A.dir_1, :, idx)
        v₂ = view(A.dir_2, :, idx)
        @. (-S * v₁ + C * v₂) / Rⱼ
    elseif short_side_left
        d₁ = view(A.d, :, idx + 1)
        @. d₁ + zero(Δt)
    else
        d₀ = view(A.d, :, idx)
        @. d₀ + zero(Δt)
    end
end
