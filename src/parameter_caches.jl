function _map_cache(prototype::StaticArray{Tuple{N}, <:Any, 1}, f, ::Val{M}) where {N, M}
    vals = ntuple(f, Val(M))
    return similar_type(prototype, eltype(typeof(vals)), Size(M))(vals)
end

function _map_cache(::AbstractArray, f, n::Integer)
    return map(f, 1:n)
end

_map_segments(t::StaticArray{Tuple{N}, <:Any, 1}, f) where {N} = _map_cache(t, f, Val(N - 1))
_map_segments(t::AbstractArray, f) = _map_cache(t, f, length(t) - 1)
_map_knots(t::StaticArray{Tuple{N}, <:Any, 1}, f) where {N} = _map_cache(t, f, Val(N))
_map_knots(t::AbstractArray, f) = _map_cache(t, f, length(t))

function _empty_cache(prototype::StaticArray, ::Type{T}) where {T}
    return similar_type(prototype, T, Size(0))()
end

function _empty_cache(::AbstractArray, ::Type{T}) where {T}
    return T[]
end

function _unzip_parameters(parameters, ::Val{K}) where {K}
    return ntuple(j -> map(p -> p[j], parameters), Val(K))
end

struct LinearParameterCache{pType}
    slope::pType
end

function LinearParameterCache(u, t, cache_parameters)
    return if cache_parameters
        slope = _map_segments(t, i -> linear_interpolation_parameters(u, t, i))
        LinearParameterCache(slope)
    else
        slope = linear_interpolation_parameters(u, t, 1)
        LinearParameterCache(_empty_cache(t, typeof(slope)))
    end
end

# Prevent e.g. Inf - Inf = NaN
function safe_diff(b, a::T) where {T}
    return isequal(b, a) ? zero(T) : b - a
end

function linear_interpolation_parameters(u::AbstractVector, t, idx)
    Δu = safe_diff(u[idx + 1], u[idx])
    Δt = t[idx + 1] - t[idx]
    slope = Δu / Δt
    slope = iszero(Δt) ? zero(slope) : slope
    return slope
end

function linear_interpolation_parameters(u::AbstractArray, t, idx)
    u₀ = _u_view(u, idx)
    u₁ = _u_view(u, idx + 1)
    Δt = t[idx + 1] - t[idx]
    return if iszero(Δt)
        @. zero(safe_diff(u₁, u₀) / oneunit(Δt))
    else
        @. safe_diff(u₁, u₀) / Δt
    end
end

struct SmoothedConstantParameterCache{dType, cType}
    d::dType
    c::cType
end

function SmoothedConstantParameterCache(
        u, t, cache_parameters, d_max, extrapolation_left, extrapolation_right
    )
    return if cache_parameters
        parameters = _map_knots(
            t,
            i -> smoothed_constant_interpolation_parameters(
                u, t, d_max, i, extrapolation_left, extrapolation_right
            ),
        )
        d, c = _unzip_parameters(parameters, Val(2))
        SmoothedConstantParameterCache(d, c)
    else
        d, c = smoothed_constant_interpolation_parameters(
            u, t, d_max, 1, extrapolation_left, extrapolation_right
        )
        SmoothedConstantParameterCache(_empty_cache(t, typeof(d)), _empty_cache(t, typeof(c)))
    end
end

function smoothed_constant_interpolation_parameters(
        u, t, d_max, idx, extrapolation_left, extrapolation_right
    )
    n = length(t)
    return if isone(idx) || (idx == n)
        # If extrapolation is periodic, make the transition differentiable
        if extrapolation_left == extrapolation_right == ExtrapolationType.Periodic
            min(t[end] - t[end - 1], t[2] - t[1], 2d_max) / 2,
                (_u_view(u, 1) - _u_view(u, n - 1)) / 2
        elseif (idx == n) && (
                extrapolation_right in (
                    ExtrapolationType.Constant, ExtrapolationType.Extension,
                )
            )
            min(t[end] - t[end - 1], 2d_max) / 2, (_u_view(u, n) - _u_view(u, n - 1)) / 2
        else
            d = isone(idx) ? min(t[2] - t[1], 2d_max) / 2 :
                min(t[end] - t[end - 1], 2d_max) / 2
            d, zero(_u_view(u, 1) / 2)
        end
    else
        min(t[idx] - t[idx - 1], t[idx + 1] - t[idx], 2d_max) / 2,
            (_u_view(u, idx) - _u_view(u, idx - 1)) / 2
    end
end

struct QuadraticParameterCache{αType, βType}
    α::αType
    β::βType
end

function QuadraticParameterCache(u, t, cache_parameters, mode)
    return if cache_parameters
        parameters = _map_segments(
            t, i -> quadratic_interpolation_parameters(u, t, i, mode)
        )
        α, β = _unzip_parameters(parameters, Val(2))
        QuadraticParameterCache(α, β)
    else
        α, β = quadratic_interpolation_parameters(u, t, 1, mode)
        QuadraticParameterCache(_empty_cache(t, typeof(α)), _empty_cache(t, typeof(β)))
    end
end

function quadratic_interpolation_parameters(u, t, idx, mode)
    # Adjust mode at boundaries
    if idx == 1
        mode = :Forward
    elseif idx == length(t) - 1
        mode = :Backward
    end

    t₀ = t[idx]
    u₀ = _u_view(u, idx)

    t₁ = t[idx + 1]
    u₁ = _u_view(u, idx + 1)

    t₂, u₂ = if mode == :Backward
        t[idx - 1], _u_view(u, idx - 1)
    else
        t[idx + 2], _u_view(u, idx + 2)
    end

    Δt₁ = t₁ - t₀
    Δt₂ = t₂ - t₀
    Δt = t₂ - t₁
    s₁ = @. (u₁ - u₀) / Δt₁
    s₂ = @. (u₂ - u₀) / Δt₂
    α = @. (s₂ - s₁) / Δt
    β = @. s₁ - α * Δt₁

    return α, β
end

struct QuadraticSplineParameterCache{αType, βType}
    α::αType
    β::βType
end

function QuadraticSplineParameterCache(u, t, k, c, cache_parameters)
    return if cache_parameters
        parameters = _map_segments(
            t, i -> quadratic_spline_parameters(u, t, k, c, i)
        )
        α, β = _unzip_parameters(parameters, Val(2))
        QuadraticSplineParameterCache(α, β)
    else
        α, β = quadratic_spline_parameters(u, t, k, c, 1)
        QuadraticSplineParameterCache(_empty_cache(t, typeof(α)), _empty_cache(t, typeof(β)))
    end
end

function quadratic_spline_parameters(u, t, k, c, idx)
    uᵢ₊ = if length(t) == 2
        # For 2 data points the knot vector has boundary multiplicity 2 < degree + 1,
        # so the B-spline basis cannot be evaluated at interior points; the spline
        # degenerates to the linear interpolant (α = 0, β = u₂ - u₁).
        (_u_view(u, 1) + _u_view(u, 2)) / 2
    else
        tᵢ₊ = (t[idx] + t[idx + 1]) / 2
        # Value of the spline at the segment midpoint via the (degree 2) B-spline basis.
        # The three nonzero basis values are evaluated into local variables rather than a
        # shared scratch buffer so that evaluation is reentrant / thread-safe (#532).
        # `tᵢ₊` is always interior, so only the non-boundary branch of the Cox-de Boor
        # recursion is needed (cf. `spline_coefficients!`).
        i = findfirst(x -> x > tᵢ₊, k)::Int - 1
        w₁ = (k[i + 1] - tᵢ₊) / (k[i + 1] - k[i])
        w₂ = (tᵢ₊ - k[i]) / (k[i + 1] - k[i])
        N₁ = (k[i + 1] - tᵢ₊) / (k[i + 1] - k[i - 1]) * w₁
        N₂ = (tᵢ₊ - k[i - 1]) / (k[i + 1] - k[i - 1]) * w₁ +
            (k[i + 2] - tᵢ₊) / (k[i + 2] - k[i]) * w₂
        N₃ = (tᵢ₊ - k[i]) / (k[i + 2] - k[i]) * w₂
        # Seed the accumulator with `zero(_u_view(u, 1))` (as the buffer-based version did)
        # so `uᵢ₊` keeps the element type of `u` for vector/array-valued data.
        uᵢ₊ = zero(_u_view(u, 1))
        uᵢ₊ += N₁ * _u_view(c, i - 2)
        uᵢ₊ += N₂ * _u_view(c, i - 1)
        uᵢ₊ += N₃ * _u_view(c, i)
        uᵢ₊
    end
    α = 2 * (_u_view(u, idx + 1) + _u_view(u, idx)) - 4uᵢ₊
    β = 4 * (uᵢ₊ - _u_view(u, idx)) - (_u_view(u, idx + 1) - _u_view(u, idx))
    return α, β
end

struct CubicSplineParameterCache{c1Type, c2Type}
    c₁::c1Type
    c₂::c2Type
end

function CubicSplineParameterCache(u, h, z, cache_parameters)
    return if cache_parameters
        parameters = if u isa AbstractVector
            _map_segments(u, i -> cubic_spline_parameters(u, h, z, i))
        else
            _map_cache(h, i -> cubic_spline_parameters(u, h, z, i), size(u)[end] - 1)
        end
        c₁, c₂ = _unzip_parameters(parameters, Val(2))
        CubicSplineParameterCache(c₁, c₂)
    else
        c₁, c₂ = cubic_spline_parameters(u, h, z, 1)
        proto = u isa AbstractVector ? u : h
        CubicSplineParameterCache(
            _empty_cache(proto, typeof(c₁)), _empty_cache(proto, typeof(c₂))
        )
    end
end

function cubic_spline_parameters(u::AbstractVector, h, z, idx)
    c₁ = (u[idx + 1] / h[idx + 1] - z[idx + 1] * h[idx + 1] / 6)
    c₂ = (u[idx] / h[idx + 1] - z[idx] * h[idx + 1] / 6)
    return c₁, c₂
end

function cubic_spline_parameters(u::AbstractArray, h, z, idx)
    c₁ = (_u_view(u, idx + 1) / h[idx + 1] - _u_view(z, idx + 1) * h[idx + 1] / 6)
    c₂ = (_u_view(u, idx) / h[idx + 1] - _u_view(z, idx) * h[idx + 1] / 6)
    return c₁, c₂
end

struct CubicHermiteParameterCache{c1Type, c2Type}
    c₁::c1Type
    c₂::c2Type
end

function CubicHermiteParameterCache(du, u, t, cache_parameters)
    return if cache_parameters
        parameters = _map_segments(
            t, i -> cubic_hermite_spline_parameters(du, u, t, i)
        )
        c₁, c₂ = _unzip_parameters(parameters, Val(2))
        CubicHermiteParameterCache(c₁, c₂)
    else
        c₁, c₂ = cubic_hermite_spline_parameters(du, u, t, 1)
        CubicHermiteParameterCache(
            _empty_cache(t, typeof(c₁)), _empty_cache(t, typeof(c₂))
        )
    end
end

function cubic_hermite_spline_parameters(du, u, t, idx)
    Δt = t[idx + 1] - t[idx]
    u₀ = _u_view(u, idx)
    u₁ = _u_view(u, idx + 1)
    du₀ = _u_view(du, idx)
    du₁ = _u_view(du, idx + 1)
    c₁ = (u₁ - u₀ - du₀ * Δt) / Δt^2
    c₂ = (du₁ - du₀ - 2c₁ * Δt) / Δt^2
    return c₁, c₂
end

struct QuinticHermiteParameterCache{c1Type, c2Type, c3Type}
    c₁::c1Type
    c₂::c2Type
    c₃::c3Type
end

function QuinticHermiteParameterCache(ddu, du, u, t, cache_parameters)
    return if cache_parameters
        parameters = _map_segments(
            t, i -> quintic_hermite_spline_parameters(ddu, du, u, t, i)
        )
        c₁, c₂, c₃ = _unzip_parameters(parameters, Val(3))
        QuinticHermiteParameterCache(c₁, c₂, c₃)
    else
        c₁, c₂, c₃ = quintic_hermite_spline_parameters(ddu, du, u, t, 1)
        QuinticHermiteParameterCache(
            _empty_cache(t, typeof(c₁)),
            _empty_cache(t, typeof(c₂)),
            _empty_cache(t, typeof(c₃)),
        )
    end
end

function quintic_hermite_spline_parameters(ddu, du, u, t, idx)
    Δt = t[idx + 1] - t[idx]
    u₀ = _u_view(u, idx)
    u₁ = _u_view(u, idx + 1)
    du₀ = _u_view(du, idx)
    du₁ = _u_view(du, idx + 1)
    ddu₀ = _u_view(ddu, idx)
    ddu₁ = _u_view(ddu, idx + 1)
    c₁ = (u₁ - u₀ - du₀ * Δt - ddu₀ * Δt^2 / 2) / Δt^3
    c₂ = (3u₀ - 3u₁ + 2(du₀ + du₁ / 2)Δt + ddu₀ * Δt^2 / 2) / Δt^4
    c₃ = (6u₁ - 6u₀ - 3(du₀ + du₁)Δt + (ddu₁ - ddu₀)Δt^2 / 2) / Δt^5
    return c₁, c₂, c₃
end

struct LagrangeParameterCache{wType, wuType}
    w::wType
    wu::wuType
end

# Barycentric weights `w[i] = 1 / ∏_{j≠i} (t[i] - t[j])`, computed once so that
# `_interpolate`/`_derivative` can evaluate in O(n) instead of O(n²) per call.
function LagrangeParameterCache(u::AbstractVector, t::AbstractVector)
    n = length(t)
    w = zeros(eltype(t), n)
    for i in 1:n
        mult = one(eltype(t))
        for j in 1:n
            i != j && (mult *= (t[i] - t[j]))
        end
        w[i] = inv(mult)
    end
    wu = w .* u
    return LagrangeParameterCache{typeof(w), typeof(wu)}(w, wu)
end

function LagrangeParameterCache(u::AbstractMatrix, t::AbstractVector)
    n = length(t)
    w = zeros(eltype(t), n)
    for i in 1:n
        mult = one(eltype(t))
        for j in 1:n
            i != j && (mult *= (t[i] - t[j]))
        end
        w[i] = inv(mult)
    end
    wu = u .* w'
    return LagrangeParameterCache{typeof(w), typeof(wu)}(w, wu)
end
