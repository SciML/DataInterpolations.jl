struct LinearParameterCache{pType}
    slope::pType
end

function LinearInterpolationParameters(u, t, idx)
    Δu = u[idx + 1] - u[idx]
    Δt = t[idx + 1] - t[idx]
    slope = Δu / Δt
    slope = iszero(Δt) ? zero(slope) : slope
    return slope
end

function LinearParameterCache(u, t)
    slope = LinearInterpolationParameters.(Ref(u), Ref(t), Base.OneTo(size(t)[1] - 1))
    return LinearParameterCache(slope)
end