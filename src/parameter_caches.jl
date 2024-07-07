struct LinearParameterCache{pType}
    slope::pType
end

function LinearParameterCache(u, t)
    slope_prototype = linear_interpolation_parameters(u, t, 1)
    idxs = 1:(length(t) - 1)
    slope = [linear_interpolation_parameters(u, t, idx) for idx in idxs]
    return LinearParameterCache(slope)
end

function linear_interpolation_parameters(
        u::AbstractVector, t::AbstractVector, idx::Integer)::Number
    Δu = u isa AbstractMatrix ? u[:, idx + 1] - u[:, idx] : u[idx + 1] - u[idx]
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

function QuadraticParameterCache(u, t)
    parameters = quadratic_interpolation_parameters.(
        Ref(u), Ref(t), 1:(length(t) - 2))
    l₀, l₁, l₂ = collect.(eachrow(hcat(collect.(parameters)...)))
    return QuadraticParameterCache(l₀, l₁, l₂)
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
    σ::pType
end

function QuadraticSplineParameterCache(z, t)
    σ = quadratic_spline_parameters.(Ref(z), Ref(t), 1:(length(t) - 1))
    return QuadraticSplineParameterCache(σ)
end

function quadratic_spline_parameters(z, t, idx)
    σ = 1 // 2 * (z[idx + 1] - z[idx]) / (t[idx + 1] - t[idx])
    return σ
end

struct CubicSplineParameterCache{pType}
    c₁::pType
    c₂::pType
end

function CubicSplineParameterCache(u, h, z)
    parameters = cubic_spline_parameters.(
        Ref(u), Ref(h), Ref(z), 1:(size(u)[end] - 1))
    c₁, c₂ = collect.(eachrow(hcat(collect.(parameters)...)))
    return CubicSplineParameterCache(c₁, c₂)
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

function CubicHermiteParameterCache(du, u, t)
    parameters = cubic_hermite_spline_parameters.(
        Ref(du), Ref(u), Ref(t), 1:(length(t) - 1))
    c₁, c₂ = collect.(eachrow(hcat(collect.(parameters)...)))
    return CubicHermiteParameterCache(c₁, c₂)
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

function QuinticHermiteParameterCache(ddu, du, u, t)
    parameters = quintic_hermite_spline_parameters.(
        Ref(ddu), Ref(du), Ref(u), Ref(t), 1:(length(t) - 1))
    c₁, c₂, c₃ = collect.(eachrow(hcat(collect.(parameters)...)))
    return QuinticHermiteParameterCache(c₁, c₂, c₃)
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
