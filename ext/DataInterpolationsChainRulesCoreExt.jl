module DataInterpolationsChainRulesCoreExt

if isdefined(Base, :get_extension)
    using DataInterpolations: _interpolate, derivative, AbstractInterpolation, get_idx,
                              interpolation_parameters,
                              LagrangeInterpolation, AkimaInterpolation,
                              BSplineInterpolation, BSplineApprox, LinearInterpolation,
                              QuadraticSpline
    using ChainRulesCore
else
    using ..DataInterpolations: _interpolate, derivative, AbstractInterpolation, get_idx,
                                interpolation_parameters,
                                LagrangeInterpolation, AkimaInterpolation,
                                BSplineInterpolation, BSplineApprox, LinearInterpolation,
                                QuadraticSpline
    using ..ChainRulesCore
end

## Linear interpolation

function ChainRulesCore.rrule(::typeof(interpolation_parameters),
        ::Val{:LinearInterpolation},
        u::AbstractVector, t::AbstractVector, idx::Integer)
    slope = interpolation_parameters(Val(:LinearInterpolation), u, t, idx)
    # TODO: use sparse arrays
    du = zero(u)
    Δt = t[idx + 1] - t[idx]

    function interpolation_parameters_pullback(Δ)
        df = NoTangent()
        dmethod = NoTangent()
        du[idx] = -Δ / Δt
        du[idx + 1] = Δ / Δt
        dt = NoTangent()
        didx = NoTangent()
        return (df, dmethod, du, dt, didx)
    end

    return slope, interpolation_parameters_pullback
end

function allocate_direct_field_tangents(A::LinearInterpolation)
    (; u = zero(A.u))
end

function allocate_parameter_tangents(A::LinearInterpolation)
    return (; slope = zero(A.p.slope))
end

function _tangent_direct_fields!(
        direct_field_tangents::NamedTuple, A::LinearInterpolation, Δt, Δ)
    (; u) = direct_field_tangents
    idx = A.idx_prev[]
    u[idx] = Δ
end

function _tangent_p!(parameter_tangents::NamedTuple, A::LinearInterpolation, Δt, Δ)
    (; slope) = parameter_tangents
    idx = A.idx_prev[]
    slope[idx] = Δt * Δ
end

## Quadratic Spline 

function ChainRulesCore.rrule(::typeof(interpolation_parameters),
        ::Val{:QuadraticSpline},
        z::AbstractVector, t::AbstractVector, idx::Integer)
    σ = interpolation_parameters(Val(:QuadraticSpline), z, t, idx)
    # TODO: use sparse arrays
    dz = zero(z)
    Δt = t[idx + 1] - t[idx]

    function interpolation_parameters_pullback(Δ)
        df = NoTangent()
        dmethod = NoTangent()
        dz[idx] = -1 // 2 * Δ / Δt
        dz[idx + 1] = 1 // 2 * Δ / Δt
        dt = NoTangent()
        didx = NoTangent()
        return (df, dmethod, dz, dt, didx)
    end

    return σ, interpolation_parameters_pullback
end

function allocate_direct_field_tangents(A::QuadraticSpline)
    (; u = zero(A.u), z = zero(A.z))
end

function allocate_parameter_tangents(A::QuadraticSpline)
    return (; σ = zero(A.p.σ))
end

function _tangent_direct_fields!(
        direct_field_tangents::NamedTuple, A::QuadraticSpline, Δt, Δ)
    (; u, z) = direct_field_tangents
    idx = A.idx_prev[]
    u[idx] = Δ
    z[idx] = Δt * Δ
end

function _tangent_p!(parameter_tangents::NamedTuple, A::QuadraticSpline, Δt, Δ)
    (; σ) = parameter_tangents
    idx = A.idx_prev[]
    σ[idx] = Δ * Δt^2
end

## generic

function ChainRulesCore.rrule(::typeof(_interpolate), A::AType, t) where {AType}
    u = _interpolate(A, t)
    # TODO: use sparse arrays
    direct_field_tangents = allocate_direct_field_tangents(A)
    parameter_tangents = allocate_parameter_tangents(A)

    function _interpolate_pullback(Δ)
        idx = get_idx(A.t, t, A.idx_prev[])
        A.idx_prev[] = idx
        Δt = t - A.t[idx]
        df = NoTangent()
        _tangent_direct_fields!(direct_field_tangents, A, Δt, Δ)
        _tangent_p!(parameter_tangents, A, Δt, Δ)
        dA = Tangent{AType}(; direct_field_tangents...,
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
