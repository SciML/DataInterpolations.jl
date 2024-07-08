module DataInterpolationsChainRulesCoreExt

if isdefined(Base, :get_extension)
    using DataInterpolations: _interpolate, derivative, AbstractInterpolation,
                              LagrangeInterpolation, AkimaInterpolation,
                              BSplineInterpolation, BSplineApprox, LinearInterpolation,
                              linear_interpolation_parameters, get_idx
    using ChainRulesCore
else
    using ..DataInterpolations: _interpolate, derivative, AbstractInterpolation,
                                LagrangeInterpolation, AkimaInterpolation,
                                BSplineInterpolation, BSplineApprox, LinearInterpolation,
                                linear_interpolation_parameters, get_idx
    using ..ChainRulesCore
end

## Linear interpolation

function ChainRulesCore.rrule(::typeof(linear_interpolation_parameters),
        u::AbstractVector, t::AbstractVector, idx::Integer)
    slope = linear_interpolation_parameters(u, t, idx)
    # TODO: use sparse arrays
    du = zero(u)

    function linear_interpolation_parameters_pullback(Δ)
        df = NoTangent()
        du .= zero(eltype(u))
        Δt = t[idx + 1] - t[idx]
        du[idx] = -Δ / Δt
        du[idx + 1] = Δ / Δt
        dt = NoTangent()
        didx = NoTangent()
        return (df, du, dt, didx)
    end

    return slope, linear_interpolation_parameters_pullback
end

function _tangent_u!(Δu::AbstractVector, A::LinearInterpolation, Δ)
    idx = A.idx_prev[]
    Δu[idx] = Δ
    Δu
end

function _tangent_p!(parameter_tangents::NamedTuple, A::LinearInterpolation, t, Δ)::Nothing
    (; slope) = parameter_tangents
    idx = A.idx_prev[]
    slope[idx] = (t - A.t[idx]) * Δ
    return nothing
end

function allocate_parameter_tangents(A::LinearInterpolation)
    return (; slope = zero(A.p.slope))
end

## Quadratic Spline 



## generic

function ChainRulesCore.rrule(::typeof(_interpolate), A::AType, t) where {AType}
    u = _interpolate(A, t)
    # TODO: use sparse arrays
    Δu = zero(A.u)
    parameter_tangents = allocate_parameter_tangents(A)

    function _interpolate_pullback(Δ)
        A.idx_prev[] = get_idx(A.t, t, A.idx_prev[])
        df = NoTangent()
        _tangent_p!(parameter_tangents, A, t, Δ)
        dA = Tangent{AType}(; u = _tangent_u!(Δu, A, Δ), t = NoTangent(),
            p = Tangent{typeof(A.p)}(; parameter_tangents...))
        dt = @thunk(derivative(A, t)*Δ)
        return df, dA, dt
    end

    u, _interpolate_pullback
end

function ChainRulesCore.frule((_, _, Δt), ::typeof(_interpolate), A::AbstractInterpolation,
        t::Number)
    return _interpolate(A, t), derivative(A, t) * Δt
end

end # module
