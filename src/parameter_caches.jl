struct LinearParameterCache{pType}
    slope::pType
end

function LinearParameterCache(u, t)
    slope = LinearInterpolationParameters.(Ref(u), Ref(t), 1:(length(t) - 1))
    return LinearParameterCache(slope)
end

function LinearInterpolationParameters(u, t, idx)
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
    parameters = QuadraticInterpolationParameters.(
        Ref(u), Ref(t), 1:(length(t) - 2))
    l₀, l₁, l₂ = collect.(eachrow(hcat(collect.(parameters)...)))
    return QuadraticParameterCache(l₀, l₁, l₂)
end

function QuadraticInterpolationParameters(u, t, idx)
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
    σ = QuadraticSplineInterpolationParameters.(Ref(z), Ref(t), 1:(length(t) - 1))
    return QuadraticSplineParameterCache(σ)
end

function QuadraticSplineInterpolationParameters(z, t, idx)
    σ = 1 // 2 * (z[idx + 1] - z[idx]) / (t[idx + 1] - t[idx])
    return σ
end

struct CubicSplineParameterCache{pType}
    c₁::pType
    c₂::pType
end

function CubicSplineParameterCache(u, h, z)
    parameters = CubicSplineInterpolationParameters.(
        Ref(u), Ref(h), Ref(z), 1:(size(u)[end] - 1))
    c₁, c₂ = collect.(eachrow(hcat(collect.(parameters)...)))
    return CubicSplineParameterCache(c₁, c₂)
end

function CubicSplineInterpolationParameters(u, h, z, idx)
    c₁ = (u[idx + 1] / h[idx + 1] - z[idx + 1] * h[idx + 1] / 6)
    c₂ = (u[idx] / h[idx + 1] - z[idx] * h[idx + 1] / 6)
    return c₁, c₂
end
