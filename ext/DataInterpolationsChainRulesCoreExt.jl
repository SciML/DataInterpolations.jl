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

function ChainRulesCore.rrule(::typeof(linear_interpolation_parameters), u::AbstractVector, t::AbstractVector, idx::Integer)
    slope = linear_interpolation_parameters(u, t, idx)
    Δt = t[idx + 1] - t[idx]
    Δu = u[idx + 1] - t[idx]
    # TODO: use sparse arrays
    du = zero(u)
    dt = zero(t)

    function linear_interpolation_parameters_pullback(Δ)
        df = NoTangent()
        du .= zero(eltype(u))
        du[idx] = - Δ / Δt
        du[idx + 1] = Δslope / Δt
        dt .= zero(eltype(t))
        dt[idx] = Δ * Δu / Δt^2
        dt[idx + 1] = - Δ * Δu / Δt^2
        didx = NoTangent()
        return (df, du, dt, didx)
    end

    return slope, linear_interpolation_parameters_pullback
end

function _tangent_u!(Δu::AbstractVector, A::LinearInterpolation)
    Δu .= zero(eltype(A.u))
    Δu
end

function _tangent_t!(Δt::AbstractVector, A::LinearInterpolation)
    idx = A.idx_prev[]
    Δt .= zero(eltype(Δt))
    Δt[idx] = one(eltype(Δt))
    Δt
end

function ChainRulesCore.rrule(::typeof(_interpolate), A::AType, t, iguess) where {AType}
    u = _interpolate(A, t)
    # TODO: use sparse arrays
    Δu = zero(A.u)
    Δt = zero(A.t)

    function _interpolate_pullback(Δ)
        df = NoTangent()
        dA = Tangent{AType}(; u = _tangent_u!(Δu, A), t = _tangent_t!(Δt, A))
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
