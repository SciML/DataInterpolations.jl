module DataInterpolationsChainRulesCoreExt
if isdefined(Base, :get_extension)
    using DataInterpolations: _interpolate, derivative, AbstractInterpolation,
                              LinearInterpolation, QuadraticInterpolation,
                              LagrangeInterpolation, AkimaInterpolation,
                              BSplineInterpolation, BSplineApprox, get_idx, get_parameters,
                              _quad_interp_indices
    using ChainRulesCore
else
    using ..DataInterpolations: _interpolate, derivative, AbstractInterpolation,
                                LinearInterpolation, QuadraticInterpolation,
                                LagrangeInterpolation, AkimaInterpolation,
                                BSplineInterpolation, BSplineApprox, get_parameters,
                                _quad_interp_indices
    using ..ChainRulesCore
end

function ChainRulesCore.rrule(
        ::Type{LinearInterpolation}, u, t, I, p, extrapolate, cache_parameters)
    A = LinearInterpolation(u, t, I, p, extrapolate, cache_parameters)
    function LinearInterpolation_pullback(ΔA)
        df = NoTangent()
        du = ΔA.u
        dt = NoTangent()
        dI = NoTangent()
        dp = NoTangent()
        dextrapolate = NoTangent()
        dcache_parameters = NoTangent()
        df, du, dt, dI, dp, dextrapolate, dcache_parameters
    end

    A, LinearInterpolation_pullback
end

function ChainRulesCore.rrule(
        ::Type{QuadraticInterpolation}, u, t, I, p, mode, extrapolate, cache_parameters)
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
        df, du, dt, dI, dp, dmode, dextrapolate, dcache_parameters
    end

    A, LinearInterpolation_pullback
end

function u_tangent(A::LinearInterpolation, t, Δ)
    out = zero(A.u)
    idx = get_idx(A, t, A.iguesser)
    t_factor = (t - A.t[idx]) / (A.t[idx + 1] - A.t[idx])
    out[idx] = Δ * (one(eltype(out)) - t_factor)
    out[idx + 1] = Δ * t_factor
    out
end

function u_tangent(A::QuadraticInterpolation, t, Δ)
    out = zero(A.u)
    i₀, i₁, i₂ = _quad_interp_indices(A, t, A.iguesser)
    t₀ = A.t[i₀]
    t₁ = A.t[i₁]
    t₂ = A.t[i₂]
    Δt₀ = t₁ - t₀
    Δt₁ = t₂ - t₁
    Δt₂ = t₂ - t₀
    out[i₀] = Δ * (t - A.t[i₁]) * (t - A.t[i₂]) / (Δt₀ * Δt₂)
    out[i₁] = -Δ * (t - A.t[i₀]) * (t - A.t[i₂]) / (Δt₀ * Δt₁)
    out[i₂] = Δ * (t - A.t[i₀]) * (t - A.t[i₁]) / (Δt₂ * Δt₁)
    out
end

function u_tangent(A, t, Δ)
    NoTangent()
end

function ChainRulesCore.rrule(::typeof(_interpolate),
        A::Union{
            LinearInterpolation,
            QuadraticInterpolation,
            LagrangeInterpolation,
            AkimaInterpolation,
            BSplineInterpolation,
            BSplineApprox
        },
        t::Number)
    deriv = derivative(A, t)
    function interpolate_pullback(Δ)
        (NoTangent(), Tangent{typeof(A)}(; u = u_tangent(A, t, Δ)), deriv * Δ)
    end
    return _interpolate(A, t), interpolate_pullback
end

function ChainRulesCore.frule((_, _, Δt), ::typeof(_interpolate), A::AbstractInterpolation,
        t::Number)
    return _interpolate(A, t), derivative(A, t) * Δt
end

end # module
