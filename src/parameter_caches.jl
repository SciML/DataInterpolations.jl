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

function quadratic_interpolation_parameters(u::AbstractArray{T, N}, t, idx) where {T, N}
    if N > 1
        ax = axes(u)
        u₀ = u[ax[1:(end - 1)]..., idx]
        u₁ = u[ax[1:(end - 1)]..., idx + 1]
        u₂ = u[ax[1:(end - 1)]..., idx + 2]
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

struct AkimaParameterCache{pType}
    b::pType
    c::pType
    d::pType
end

function AkimaParameterCache(u, t)
    b, c, d = akima_interpolation_parameters(u, t)
    AkimaParameterCache(b, c, d)
end

function akima_interpolation_parameters(u::AbstractVector, t)
    n = length(t)
    dt = diff(t)
    m = Array{eltype(u)}(undef, n + 3)
    m[3:(end - 2)] = diff(u) ./ dt
    m[2] = 2m[3] - m[4]
    m[1] = 2m[2] - m[3]
    m[end - 1] = 2m[end - 2] - m[end - 3]
    m[end] = 2m[end - 1] - m[end - 2]
    b = 0.5 .* (m[4:end] .+ m[1:(end - 3)])
    dm = abs.(diff(m))
    f1 = dm[3:(n + 2)]
    f2 = dm[1:n]
    f12 = f1 + f2
    ind = findall(f12 .> 1e-9 * maximum(f12))
    b[ind] = (f1[ind] .* m[ind .+ 1] .+
              f2[ind] .* m[ind .+ 2]) ./ f12[ind]
    c = (3.0 .* m[3:(end - 2)] .- 2.0 .* b[1:(end - 1)] .- b[2:end]) ./ dt
    d = (b[1:(end - 1)] .+ b[2:end] .- 2.0 .* m[3:(end - 2)]) ./ dt .^ 2
    return b, c, d
end

function akima_interpolation_parameters(u::AbstractArray, t)
    n = length(t)
    dt = diff(t)
    ax = axes(u)[1:(end - 1)]
    su = size(u)
    m = zeros(eltype(u), su[1:(end - 1)]..., n + 3)
    m[ax..., 3:(end - 2)] .= mapslices(
        x -> x ./ dt, diff(u, dims = length(su)); dims = length(su))
    m[ax..., 2] .= 2m[ax..., 3] .- m[ax..., 4]
    m[ax..., 1] .= 2m[ax..., 2] .- m[3]
    m[ax..., end - 1] .= 2m[ax..., end - 2] - m[ax..., end - 3]
    m[ax..., end] .= 2m[ax..., end - 1] .- m[ax..., end - 2]
    b = 0.5 .* (m[ax..., 4:end] .+ m[ax..., 1:(end - 3)])
    dm = abs.(diff(m, dims = length(su)))
    f1 = dm[ax..., 3:(n + 2)]
    f2 = dm[ax..., 1:n]
    f12 = f1 .+ f2
    ind = findall(f12 .> 1e-9 * maximum(f12))
    indi = map(i -> i.I, ind)
    b[ind] .= (f1[ind] .*
               m[CartesianIndex.(map(i -> (i[1:(end - 1)]..., i[end] + 1), indi))] .+
               f2[ind] .*
               m[CartesianIndex.(map(i -> (i[1:(end - 1)]..., i[end] + 2), indi))]) ./
              f12[ind]
    c = mapslices(x -> x ./ dt,
        (3.0 .* m[ax..., 3:(end - 2)] .- 2.0 .* b[ax..., 1:(end - 1)] .- b[ax..., 2:end]);
        dims = length(su))
    d = mapslices(x -> x ./ dt .^ 2,
        (b[ax..., 1:(end - 1)] .+ b[ax..., 2:end] .- 2.0 .* m[ax..., 3:(end - 2)]);
        dims = length(su))
    return b, c, d
end

struct QuadraticSplineParameterCache{pType}
    σ::pType
end

function QuadraticSplineParameterCache(z, t, cache_parameters)
    if cache_parameters
        σ = quadratic_spline_parameters.(Ref(z), Ref(t), 1:(length(t) - 1))
        QuadraticSplineParameterCache(σ)
    else
        # Compute parameters once to infer types
        σ = quadratic_spline_parameters(z, t, 1)
        QuadraticSplineParameterCache(typeof(σ)[])
    end
end

function quadratic_spline_parameters(z::AbstractVector, t, idx)
    σ = 1 // 2 * (z[idx + 1] - z[idx]) / (t[idx + 1] - t[idx])
    return σ
end

function quadratic_spline_parameters(z::AbstractArray, t, idx)
    ax = axes(z)[1:(end - 1)]
    σ = 1 // 2 * (z[ax..., idx + 1] - z[ax..., idx]) / (t[idx + 1] - t[idx])
    return σ
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

function cubic_hermite_spline_parameters(du::AbstractArray, u, t, idx)
    ax = axes(u)[1:(end - 1)]
    Δt = t[idx + 1] - t[idx]
    u₀ = u[ax..., idx]
    u₁ = u[ax..., idx + 1]
    du₀ = du[ax..., idx]
    du₁ = du[ax..., idx + 1]
    c₁ = (u₁ - u₀ - du₀ * Δt) / Δt^2
    c₂ = (du₁ - du₀ - 2c₁ * Δt) / Δt^2
    return c₁, c₂
end

function cubic_hermite_spline_parameters(du::AbstractVector, u, t, idx)
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

function quintic_hermite_spline_parameters(ddu::AbstractVector, du, u, t, idx)
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

function quintic_hermite_spline_parameters(ddu::AbstractArray, du, u, t, idx)
    ax = axes(ddu)[1:(end - 1)]
    Δt = t[idx + 1] - t[idx]
    u₀ = u[ax..., idx]
    u₁ = u[ax..., idx + 1]
    du₀ = du[ax..., idx]
    du₁ = du[ax..., idx + 1]
    ddu₀ = ddu[ax..., idx]
    ddu₁ = ddu[ax..., idx + 1]
    c₁ = (u₁ - u₀ - du₀ * Δt - ddu₀ * Δt^2 / 2) / Δt^3
    c₂ = (3u₀ - 3u₁ + 2(du₀ + du₁ / 2)Δt + ddu₀ * Δt^2 / 2) / Δt^4
    c₃ = (6u₁ - 6u₀ - 3(du₀ + du₁)Δt + (ddu₁ - ddu₀)Δt^2 / 2) / Δt^5
    return c₁, c₂, c₃
end
