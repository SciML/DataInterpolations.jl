struct CubicHermiteParameterCache{pType}
    c₁::pType
    c₂::pType
end

function CubicHermiteParameterCache(du, u, t)
    parameters = CubicHermiteInterpolationParameters.(
        Ref(du), Ref(u), Ref(t), 1:(length(t) - 1))
    c₁, c₂ = collect.(eachrow(hcat(collect.(parameters)...)))
    return CubicHermiteParameterCache(c₁, c₂)
end

function CubicHermiteInterpolationParameters(du, u, t, idx)
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
    parameters = QuinticHermiteInterpolationParameters.(
        Ref(ddu), Ref(du), Ref(u), Ref(t), 1:(length(t) - 1))
    c₁, c₂, c₃ = collect.(eachrow(hcat(collect.(parameters)...)))
    return QuinticHermiteParameterCache(c₁, c₂, c₃)
end

function QuinticHermiteInterpolationParameters(ddu, du, u, t, idx)
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
