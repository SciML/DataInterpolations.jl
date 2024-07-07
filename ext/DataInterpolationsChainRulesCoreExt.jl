module DataInterpolationsChainRulesCoreExt

if isdefined(Base, :get_extension)
    using DataInterpolations: _interpolate, derivative, AbstractInterpolation,
                              LagrangeInterpolation, AkimaInterpolation,
                              BSplineInterpolation, BSplineApprox, LinearInterpolation,
                              linear_interpolation_parameters
    using ChainRulesCore
else
    using ..DataInterpolations: _interpolate, derivative, AbstractInterpolation,
                                LagrangeInterpolation, AkimaInterpolation,
                                BSplineInterpolation, BSplineApprox, LinearInterpolation,
                                linear_interpolation_parameters
    using ..ChainRulesCore
end

function ChainRulesCore.rrule(::typeof(linear_interpolation_parameters), u, t, idx)
    slope = linear_interpolation_parameters(u, t, idx)

    function linear_interpolation_parameters_pullback(Δslope)
        du = @thunk(Δslope*...) # TODO: How to handle sparsity?
        dt = @thunk(Δslope*...) # TODO: How to handle sparsity?
        df = NoTangent()
        didx = NoTangent()

        return (df, du, dt, didx)
    end

    return slope, linear_interpolation_parameters_pullback
end

function _tangent_u(A::LinearInterpolation, t)
    ... # TODO: How to handle sparsity?
end

function _tangent_t(A::LinearInterpolation, t)
    ... # TODO: How to handle sparsity?
end

function ChainRulesCore.rrule(::typeof(_interpolate), A::AType, t, iguess) where {AType}
    u = _interpolate(A, t, iguess)[1]
    function _interpolate_pullback(Δ)
        df = NoTangent()
        dA = Tangent{AType}(; u = _tangent_u(A, t), t = _tangent_t(A, t))
        dt = @thunk(derivative(A, t)*Δ)
        diguess = NoTangent()
        return df, dA, dt, diguess
    end
    u, _interpolate_pullback
end

function ChainRulesCore.frule((_, _, Δt), ::typeof(_interpolate), A::AbstractInterpolation,
        t::Number)
    return _interpolate(A, t), derivative(A, t) * Δt
end

end # module
