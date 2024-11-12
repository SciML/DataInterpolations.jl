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
    b == a ? zero(T) : b - a
end

function linear_interpolation_parameters(u::AbstractArray{T}, t, idx) where {T}
    Δu = if u isa AbstractMatrix
        [safe_diff(u[j, idx + 1], u[j, idx]) for j in 1:size(u)[1]]
    else
        safe_diff(u[idx + 1], u[idx])
    end
    Δt = t[idx + 1] - t[idx]
    slope = Δu / Δt
    slope = iszero(Δt) ? zero(slope) : slope
    return slope
end

struct QuadraticParameterCache{pType}
    l₀::pType
    l₁::pType
    l₂::pType
end

function QuadraticParameterCache(u, t, cache_parameters)
    if cache_parameters
        parameters = quadratic_interpolation_parameters.(
            Ref(u), Ref(t), 1:(length(t) - 2))
        l₀, l₁, l₂ = collect.(eachrow(stack(collect.(parameters))))
        QuadraticParameterCache(l₀, l₁, l₂)
    else
        # Compute parameters once to infer types
        l₀, l₁, l₂ = quadratic_interpolation_parameters(u, t, 1)
        pType = typeof(l₀)
        QuadraticParameterCache(pType[], pType[], pType[])
    end
end

function quadratic_interpolation_parameters(u, t, idx)
    if u isa AbstractMatrix
        u₀ = u[:, idx]
        u₁ = u[:, idx + 1]
        u₂ = u[:, idx + 2]
    else
        u₀ = u[idx]
        u₁ = u[idx + 1]
        u₂ = u[idx + 2]
    end
    t₀ = t[idx]
    t₁ = t[idx + 1]
    t₂ = t[idx + 2]
    Δt₀ = t₁ - t₀
    Δt₁ = t₂ - t₁
    Δt₂ = t₂ - t₀
    l₀ = u₀ / (Δt₀ * Δt₂)
    l₁ = -u₁ / (Δt₀ * Δt₁)
    l₂ = u₂ / (Δt₂ * Δt₁)
    return l₀, l₁, l₂
end

struct QuadraticSplineParameterCache{pType}
    α::pType
    β::pType
end

function QuadraticSplineParameterCache(u, t, k, c, N, cache_parameters)
    if cache_parameters
        parameters = quadratic_spline_parameters.(
            Ref(u), Ref(t), Ref(k), Ref(c), Ref(N), 1:(length(t) - 1))
        α, β = collect.(eachrow(stack(collect.(parameters))))
        QuadraticSplineParameterCache(α, β)
    else
        # Compute parameters once to infer types
        α, β = quadratic_spline_parameters(u, t, k, c, N, 1)
        QuadraticSplineParameterCache(typeof(α)[], typeof(β)[])
    end
end

function quadratic_spline_parameters(u, t, k, c, N, idx)
    tᵢ₊ = (t[idx] + t[idx + 1]) / 2
    nonzero_coefficient_idxs = spline_coefficients!(N, 2, k, tᵢ₊)
    uᵢ₊ = zero(first(u))
    for j in nonzero_coefficient_idxs
        uᵢ₊ += N[j] * c[j]
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

function cubic_spline_parameters(u, h, z, idx)
    c₁ = (u[idx + 1] / h[idx + 1] - z[idx + 1] * h[idx + 1] / 6)
    c₂ = (u[idx] / h[idx + 1] - z[idx] * h[idx + 1] / 6)
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
