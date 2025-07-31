struct LinearParameterCache{pType}
    slope::pType
end

function LinearParameterCache(u, t, cache_parameters)
    if cache_parameters
        slope = linear_interpolation_parameters.(Ref(u), Ref(t), 1:(length(t) - 1))
        LinearParameterCache(slope)
    else
        # Compute parameters once to infer types
        slope = linear_interpolation_parameters(u, t, 1)
        LinearParameterCache(typeof(slope)[])
    end
end

# Prevent e.g. Inf - Inf = NaN
function safe_diff(b, a::T) where {T}
    isequal(b, a) ? zero(T) : b - a
end

function linear_interpolation_parameters(u::AbstractArray{T, N}, t, idx) where {T, N}
    Δu = if N > 1
        ax = axes(u)
        safe_diff.(
            u[ax[1:(end - 1)]..., (idx + 1)], u[ax[1:(end - 1)]..., idx])
    else
        safe_diff(u[idx + 1], u[idx])
    end
    Δt = t[idx + 1] - t[idx]
    slope = Δu / Δt
    slope = iszero(Δt) ? zero(slope) : slope
    return slope
end

struct SmoothedLinearParameterCache{uType, tType, λType}
    Δu::uType
    Δt::tType
    ΔΔu::uType
    ΔΔt::tType
    u_tilde::uType
    t_tilde::tType
    slope::uType
    # Whether ΔΔt is sufficiently close to 0
    degenerate_ΔΔt::Vector{Bool}
    λ::λType
end

function get_spline_ends(u, Δu, λ)
    u_tilde = zeros(2 * length(u))
    u_tilde[1] = u[1]
    u_tilde[2:2:(end - 1)] = u[1:(end - 1)] .+ (λ / 2) .* Δu[2:(end - 1)]
    u_tilde[3:2:end] = u[2:end] .- (λ / 2) .* Δu[2:(end - 1)]
    u_tilde[end] = u[end]
    return u_tilde
end

function SmoothedLinearParameterCache(u::AbstractVector, t, λ, cache_parameters)
    @assert cache_parameters "Parameter caching is mandatory for SmoothedLinearInterpolation"

    Δu = diff(u)
    Δt = diff(t)
    pushfirst!(Δt, last(Δt))
    push!(Δt, first(Δt))
    pushfirst!(Δu, last(Δu))
    push!(Δu, first(Δu))
    ΔΔu = diff(Δu)
    ΔΔt = diff(Δt)
    u_tilde = get_spline_ends(u, Δu, λ)
    t_tilde = get_spline_ends(t, Δt, λ)
    slope = Δu ./ Δt
    degenerate_ΔΔt = collect(isapprox.(ΔΔt, 0, atol = 1e-5))

    return SmoothedLinearParameterCache(
        Δu, Δt, ΔΔu, ΔΔt, u_tilde, t_tilde, slope, degenerate_ΔΔt, λ)
end

struct SmoothedConstantParameterCache{dType, cType}
    d::dType
    c::cType
end

function SmoothedConstantParameterCache(
        u, t, cache_parameters, d_max, extrapolation_left, extrapolation_right)
    if cache_parameters
        parameters = smoothed_constant_interpolation_parameters.(
            Ref(u), Ref(t), d_max, eachindex(t), extrapolation_left, extrapolation_right)
        d, c = collect.(eachrow(stack(collect.(parameters))))
        SmoothedConstantParameterCache(d, c)
    else
        SmoothedConstantParameterCache(eltype(t)[], eltype(u)[])
    end
end

function smoothed_constant_interpolation_parameters(
        u, t, d_max, idx, extrapolation_left, extrapolation_right)
    if isone(idx) || (idx == length(t))
        # If extrapolation is periodic, make the transition differentiable
        if extrapolation_left == extrapolation_right == ExtrapolationType.Periodic
            min(t[end] - t[end - 1], t[2] - t[1], 2d_max) / 2, (u[1] - u[end - 1]) / 2
        elseif (idx == length(t)) && (extrapolation_right in (
            ExtrapolationType.Constant, ExtrapolationType.Extension))
            min(t[end] - t[end - 1], 2d_max) / 2, (u[end] - u[end - 1]) / 2
        else
            d = isone(idx) ? min(t[2] - t[1], 2d_max) / 2 :
                min(t[end] - t[end - 1], 2d_max) / 2
            d, zero(one(eltype(u)) / 2)
        end
    else
        min(t[idx] - t[idx - 1], t[idx + 1] - t[idx], 2d_max) / 2, (u[idx] - u[idx - 1]) / 2
    end
end

struct QuadraticParameterCache{pType}
    α::pType
    β::pType
end

function QuadraticParameterCache(u, t, cache_parameters, mode)
    if cache_parameters
        parameters = quadratic_interpolation_parameters.(
            Ref(u), Ref(t), 1:(length(t) - 1), mode)
        α, β = collect.(eachrow(stack(collect.(parameters))))
        QuadraticParameterCache(α, β)
    else
        # Compute parameters once to infer types
        α, β = quadratic_interpolation_parameters(u, t, 1, mode)
        pType = typeof(α)
        QuadraticParameterCache(pType[], pType[])
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
    u₀ = u isa AbstractMatrix ? view(u, :, idx) : u[idx]

    t₁ = t[idx + 1]
    u₁ = u isa AbstractMatrix ? view(u, :, idx + 1) : u[idx + 1]

    t₂, u₂ = if mode == :Backward
        t[idx - 1], u isa AbstractMatrix ? view(u, :, idx - 1) : u[idx - 1]
    else
        t[idx + 2], u isa AbstractMatrix ? view(u, :, idx + 2) : u[idx + 2]
    end

    Δt₁ = t₁ - t₀
    Δt₂ = t₂ - t₀
    Δt = t₂ - t₁
    s₁ = (u₁ - u₀) / Δt₁
    s₂ = (u₂ - u₀) / Δt₂
    α = (s₂ - s₁) / Δt
    β = s₁ - α * Δt₁

    α, β
end

struct QuadraticSplineParameterCache{pType}
    α::pType
    β::pType
end

function QuadraticSplineParameterCache(u, t, k, c, sc, cache_parameters)
    if cache_parameters
        parameters = quadratic_spline_parameters.(
            Ref(u), Ref(t), Ref(k), Ref(c), Ref(sc), 1:(length(t) - 1))
        α, β = collect.(eachrow(stack(collect.(parameters))))
        QuadraticSplineParameterCache(α, β)
    else
        # Compute parameters once to infer types
        α, β = quadratic_spline_parameters(u, t, k, c, sc, 1)
        QuadraticSplineParameterCache(typeof(α)[], typeof(β)[])
    end
end

function quadratic_spline_parameters(u, t, k, c, sc, idx)
    tᵢ₊ = (t[idx] + t[idx + 1]) / 2
    nonzero_coefficient_idxs = spline_coefficients!(sc, 2, k, tᵢ₊)
    uᵢ₊ = zero(first(u))
    for j in nonzero_coefficient_idxs
        uᵢ₊ += sc[j] * c[j]
    end
    α = 2 * (u[idx + 1] + u[idx]) - 4uᵢ₊
    β = 4 * (uᵢ₊ - u[idx]) - (u[idx + 1] - u[idx])
    return α, β
end

struct CubicSplineParameterCache{pType}
    c₁::pType
    c₂::pType
end

function CubicSplineParameterCache(u, h, z, cache_parameters)
    if cache_parameters
        parameters = cubic_spline_parameters.(
            Ref(u), Ref(h), Ref(z), 1:(size(u)[end] - 1))
        c₁, c₂ = collect.(eachrow(stack(collect.(parameters))))
        CubicSplineParameterCache(c₁, c₂)
    else
        # Compute parameters once to infer types
        c₁, c₂ = cubic_spline_parameters(u, h, z, 1)
        pType = typeof(c₁)
        CubicSplineParameterCache(pType[], pType[])
    end
end

function cubic_spline_parameters(u::AbstractVector, h, z, idx)
    c₁ = (u[idx + 1] / h[idx + 1] - z[idx + 1] * h[idx + 1] / 6)
    c₂ = (u[idx] / h[idx + 1] - z[idx] * h[idx + 1] / 6)
    return c₁, c₂
end

function cubic_spline_parameters(u::AbstractArray, h, z, idx)
    ax = axes(u)[1:(end - 1)]
    c₁ = (u[ax..., idx + 1] / h[idx + 1] - z[ax..., idx + 1] * h[idx + 1] / 6)
    c₂ = (u[ax..., idx] / h[idx + 1] - z[ax..., idx] * h[idx + 1] / 6)
    return c₁, c₂
end

struct CubicHermiteParameterCache{pType}
    c₁::pType
    c₂::pType
end

function CubicHermiteParameterCache(du, u, t, cache_parameters)
    if cache_parameters
        parameters = cubic_hermite_spline_parameters.(
            Ref(du), Ref(u), Ref(t), 1:(length(t) - 1))
        c₁, c₂ = collect.(eachrow(stack(collect.(parameters))))
        CubicHermiteParameterCache(c₁, c₂)
    else
        # Compute parameters once to infer types
        c₁, c₂ = cubic_hermite_spline_parameters(du, u, t, 1)
        pType = typeof(c₁)
        CubicHermiteParameterCache(pType[], pType[])
    end
end

function cubic_hermite_spline_parameters(du, u, t, idx)
    Δt = t[idx + 1] - t[idx]
    u₀ = u[idx]
    u₁ = u[idx + 1]
    du₀ = du[idx]
    du₁ = du[idx + 1]
    c₁ = (u₁ - u₀ - du₀ * Δt) / Δt^2
    c₂ = (du₁ - du₀ - 2c₁ * Δt) / Δt^2
    return c₁, c₂
end

struct QuinticHermiteParameterCache{pType}
    c₁::pType
    c₂::pType
    c₃::pType
end

function QuinticHermiteParameterCache(ddu, du, u, t, cache_parameters)
    if cache_parameters
        parameters = quintic_hermite_spline_parameters.(
            Ref(ddu), Ref(du), Ref(u), Ref(t), 1:(length(t) - 1))
        c₁, c₂, c₃ = collect.(eachrow(stack(collect.(parameters))))
        QuinticHermiteParameterCache(c₁, c₂, c₃)
    else
        # Compute parameters once to infer types
        c₁, c₂, c₃ = quintic_hermite_spline_parameters(ddu, du, u, t, 1)
        pType = typeof(c₁)
        QuinticHermiteParameterCache(pType[], pType[], pType[])
    end
end

function quintic_hermite_spline_parameters(ddu, du, u, t, idx)
    Δt = t[idx + 1] - t[idx]
    u₀ = u[idx]
    u₁ = u[idx + 1]
    du₀ = du[idx]
    du₁ = du[idx + 1]
    ddu₀ = ddu[idx]
    ddu₁ = ddu[idx + 1]
    c₁ = (u₁ - u₀ - du₀ * Δt - ddu₀ * Δt^2 / 2) / Δt^3
    c₂ = (3u₀ - 3u₁ + 2(du₀ + du₁ / 2)Δt + ddu₀ * Δt^2 / 2) / Δt^4
    c₃ = (6u₁ - 6u₀ - 3(du₀ + du₁)Δt + (ddu₁ - ddu₀)Δt^2 / 2) / Δt^5
    return c₁, c₂, c₃
end
