module DataInterpolationsChainRulesCoreExt

using DataInterpolations: _interpolate, derivative, AbstractInterpolation,
    LinearInterpolation, QuadraticInterpolation,
    LagrangeInterpolation, AkimaInterpolation,
    BSplineInterpolation, BSplineApprox, get_idx, get_parameters,
    munge_data
using ChainRulesCore

function ChainRulesCore.rrule(::typeof(munge_data), u, t)
    u_out, t_out = munge_data(u, t)

    # For now modifications by munge_data not supported
    @assert (u == u_out && t == t_out)

    munge_data_pullback = Δ -> (NoTangent(), Δ[1], Δ[2])
    return (u_out, t_out), munge_data_pullback
end

function ChainRulesCore.rrule(
        ::Type{LinearInterpolation}, u, t, I, p, extrapolate, cache_parameters
    )
    A = LinearInterpolation(u, t, I, p, extrapolate, cache_parameters)
    function LinearInterpolation_pullback(ΔA)
        df = NoTangent()
        du = ΔA.u
        dt = NoTangent()
        dI = NoTangent()
        dp = NoTangent()
        dextrapolate = NoTangent()
        dcache_parameters = NoTangent()
        return df, du, dt, dI, dp, dextrapolate, dcache_parameters
    end

    return A, LinearInterpolation_pullback
end

function ChainRulesCore.rrule(
        ::Type{QuadraticInterpolation}, u, t, I, p, mode, extrapolate, cache_parameters
    )
    A = QuadraticInterpolation(u, t, I, p, mode, extrapolate, cache_parameters)
    function LinearInterpolation_pullback(ΔA)
        df = NoTangent()
        du = ΔA.u
        dt = NoTangent()
        dI = NoTangent()
        dp = NoTangent()
        dmode = NoTangent()
        dextrapolate = NoTangent()
        dcache_parameters = NoTangent()
        return df, du, dt, dI, dp, dmode, dextrapolate, dcache_parameters
    end

    return A, LinearInterpolation_pullback
end

function u_tangent(A::LinearInterpolation, t, Δ)
    out = zero.(A.u)
    idx = get_idx(A, t, A.iguesser)
    t_factor = (t - A.t[idx]) / (A.t[idx + 1] - A.t[idx])
    if out isa AbstractVector{<:Number}
        out[idx] = Δ * (one(eltype(out)) - t_factor)
        out[idx + 1] = Δ * t_factor
    elseif out isa AbstractMatrix
        out[:, idx] = Δ * (one(eltype(out)) - t_factor)
        out[:, idx + 1] = Δ * t_factor
    else
        @. out[idx] = Δ * (true - t_factor)
        @. out[idx + 1] = Δ * t_factor
    end
    return out
end

function _quad_interp_indices(A::QuadraticInterpolation, t::Number, iguess)
    idx = get_idx(A, t, iguess; idx_shift = A.mode == :Backward ? -1 : 0, ub_shift = -2)
    return idx, idx + 1, idx + 2
end

function u_tangent(A::QuadraticInterpolation, t, Δ)
    out = zero.(A.u)
    i₀, i₁, i₂ = _quad_interp_indices(A, t, A.iguesser)
    t₀ = A.t[i₀]
    t₁ = A.t[i₁]
    t₂ = A.t[i₂]
    Δt₀ = t₁ - t₀
    Δt₁ = t₂ - t₁
    Δt₂ = t₂ - t₀
    Δt_rel₀ = t - A.t[i₀]
    Δt_rel₁ = t - A.t[i₁]
    Δt_rel₂ = t - A.t[i₂]
    if out isa AbstractVector{<:Number}
        out[i₀] = Δ * Δt_rel₁ * Δt_rel₂ / (Δt₀ * Δt₂)
        out[i₁] = -Δ * Δt_rel₀ * Δt_rel₂ / (Δt₀ * Δt₁)
        out[i₂] = Δ * Δt_rel₀ * Δt_rel₁ / (Δt₂ * Δt₁)
    elseif out isa AbstractMatrix
        out[:, i₀] = Δ * Δt_rel₁ * Δt_rel₂ / (Δt₀ * Δt₂)
        out[:, i₁] = -Δ * Δt_rel₀ * Δt_rel₂ / (Δt₀ * Δt₁)
        out[:, i₂] = Δ * Δt_rel₀ * Δt_rel₁ / (Δt₂ * Δt₁)
    else
        @. out[i₀] = Δ * Δt_rel₁ * Δt_rel₂ / (Δt₀ * Δt₂)
        @. out[i₁] = -Δ * Δt_rel₀ * Δt_rel₂ / (Δt₀ * Δt₁)
        @. out[i₂] = Δ * Δt_rel₀ * Δt_rel₁ / (Δt₂ * Δt₁)
    end
    return out
end

function u_tangent(A, t, Δ)
    return NoTangent()
end

function ChainRulesCore.rrule(
        ::typeof(_interpolate),
        A::Union{
            LinearInterpolation,
            QuadraticInterpolation,
            LagrangeInterpolation,
            AkimaInterpolation,
            BSplineInterpolation,
            BSplineApprox,
        },
        t::Number
    )
    deriv = derivative(A, t)
    function interpolate_pullback(Δ)
        return (NoTangent(), Tangent{typeof(A)}(; u = u_tangent(A, t, Δ)), sum(deriv .* Δ))
    end
    return _interpolate(A, t), interpolate_pullback
end

function ChainRulesCore.frule(
        (_, _, Δt), ::typeof(_interpolate), A::AbstractInterpolation,
        t::Number
    )
    return _interpolate(A, t), derivative(A, t) * Δt
end

function ChainRulesCore.frule((_, Δt), A::AbstractInterpolation, t::Number)
    return A(t), derivative(A, t) * Δt
end

end # module
