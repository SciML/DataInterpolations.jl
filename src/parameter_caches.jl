struct QuinticHermiteParameterCache{pType}
    c₁::pType
    c₂::pType
    c₃::pType
end

function QuinticHermiteParameterCache(u, t, du, ddu)
    parameters = QuinticHermiteInterpolationParameters.(
        Ref(u), Ref(t), Ref(du), Ref(ddu), 1:(length(t) - 1))
    c₁, c₂, c₃ = collect.(eachrow(hcat(collect.(parameters)...)))
    return QuinticHermiteParameterCache(c₁, c₂, c₃)
end

function QuinticHermiteInterpolationParameters(u, t, du, ddu, idx)
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
