module DataInterpolationsChainRulesCoreExt

if isdefined(Base, :get_extension)
    using DataInterpolations: _interpolate, derivative, AbstractInterpolation,
        LagrangeInterpolation, AkimaInterpolation, BSplineInterpolation, BSplineApprox
    using ChainRulesCore
else
    using ..DataInterpolations: _interpolate, derivative, AbstractInterpolation,
        LagrangeInterpolation, AkimaInterpolation, BSplineInterpolation, BSplineApprox
    using ..ChainRulesCore
end

function ChainRulesCore.rrule(::typeof(_interpolate),
    A::Union{
        LagrangeInterpolation,
        AkimaInterpolation,
        BSplineInterpolation,
        BSplineApprox,
    },
    t::Number)
    deriv = derivative(A, t)
    interpolate_pullback(Δ) = (NoTangent(), NoTangent(), deriv * Δ)
    return _interpolate(A, t), interpolate_pullback
end

function ChainRulesCore.frule((_, _, Δt), ::typeof(_interpolate), A::AbstractInterpolation,
    t::Number)
    return _interpolate(A, t), derivative(A, t) * Δt
end

end # module
