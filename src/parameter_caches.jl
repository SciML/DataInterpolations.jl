struct LinearParameterCache{pType}
    slope::pType
end

function LinearParameterCache(u, t)
    slope = LinearInterpolationParameters.(Ref(u), Ref(t), Base.OneTo(length(t) - 1))
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
        Ref(u), Ref(t), Base.OneTo(length(t) - 2))
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
