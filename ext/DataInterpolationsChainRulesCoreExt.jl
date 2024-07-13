module DataInterpolationsChainRulesCoreExt

if isdefined(Base, :get_extension)
    using DataInterpolations: _interpolate, derivative, AbstractInterpolation, get_idx,
                              cumulative_integral, LinearParameterCache,
                              QuadraticSplineParameterCache,
                              LagrangeInterpolation, AkimaInterpolation,
                              BSplineInterpolation, BSplineApprox, LinearInterpolation,
                              QuadraticSpline
    using ChainRulesCore
    using LinearAlgebra
    using SparseArrays
    using ReadOnlyArrays
else
    using ..DataInterpolations: _interpolate, derivative, AbstractInterpolation, get_idx,
                                cumulative_integral, LinearParameterCache,
                                QuadraticSplineParameterCache,
                                LagrangeInterpolation, AkimaInterpolation,
                                BSplineInterpolation, BSplineApprox, LinearInterpolation,
                                QuadraticSpline
    using ..ChainRulesCore
    using ..LinearAlgebra
    using ..SparseArrays
    using ..ReadOnlyArrays
end

## Linear interpolation

function ChainRulesCore.rrule(
        ::Type{LinearParameterCache}, u::AbstractArray, t::AbstractVector)
    p = LinearParameterCache(u, t)
    du = zeros(eltype(p.slope), length(u))

    function LinearParameterCache_pullback(Δp)
        df = NoTangent()
        du[2:end] += Δp.slope
        du[1:(end - 1)] -= Δp.slope
        dt = NoTangent()
        return (df, du, dt)
    end

    p, LinearParameterCache_pullback
end

function ChainRulesCore.rrule(
        ::Type{LinearInterpolation}, u, t, I, p, extrapolate, safetycopy)
    A = LinearInterpolation(u, t, I, p, extrapolate, safetycopy)

    function LinearInterpolation_pullback(ΔA)
        df = NoTangent()
        du = ΔA.u
        dt = NoTangent()
        dI = NoTangent()
        dp = ΔA.p
        dextrapolate = NoTangent()
        dsafetycopy = NoTangent()
        return df, du, dt, dI, dp, dextrapolate, dsafetycopy
    end

    A, LinearInterpolation_pullback
end

function allocate_direct_field_tangents(A::LinearInterpolation)
    idx = A.idx_prev[]
    u = SparseVector(length(A.u), [idx], zeros(1))
    (; u)
end

function allocate_parameter_tangents(A::LinearInterpolation)
    idx = A.idx_prev[]
    slope = SparseVector(length(A.p.slope), [idx], zeros(1))
    return (; slope)
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

function ChainRulesCore.rrule(::Type{QuadraticSplineParameterCache}, u, t)
    p = QuadraticSplineParameterCache(u, t)
    n = length(u)

    Δt = diff(t)
    diagonal_main = [zero(eltype(Δt)), 2 ./ Δt...]
    diagonal_down = -diagonal_main[2:end]
    diagonal_up = zero(diagonal_down)
    ∂d_∂u = Tridiagonal(diagonal_down, diagonal_main, diagonal_up)

    ∂σ_∂z = spzeros(n, n - 1)
    for i in 1:(n - 1)
        ∂σ_∂z[i, i] = -0.5 / Δt[i]
        ∂σ_∂z[i + 1, i] = 0.5 / Δt[i]
    end

    function QuadraticSplineParameterCache_pullback(Δp)
        df = NoTangent()
        temp1 = Δp.z + ∂σ_∂z * Δp.σ
        temp2 = p.tA' \ temp1
        du = ∂d_∂u' * temp2
        dt = NoTangent()
        return (df, du, dt)
    end

    p, QuadraticSplineParameterCache_pullback
end

function ChainRulesCore.rrule(::Type{QuadraticSpline}, u, t, I, p, extrapolate, safetycopy)
    A = QuadraticSpline(u, t, I, p, extrapolate, safetycopy)

    function QuadraticSpline_pullback(ΔA)
        df = NoTangent()
        du = ΔA.u
        dt = NoTangent()
        dI = NoTangent()
        dp = ΔA.p
        dextrapolate = NoTangent()
        dsafetycopy = NoTangent()
        return df, du, dt, dI, dp, dextrapolate, dsafetycopy
    end

    A, QuadraticSpline_pullback
end

function allocate_direct_field_tangents(A::QuadraticSpline)
    idx = A.idx_prev[]
    u = SparseVector(length(A.u), [idx], zeros(1))
    (; u)
end

function allocate_parameter_tangents(A::QuadraticSpline)
    idx = A.idx_prev[]
    z = SparseVector(length(A.p.z), [idx], zeros(1))
    σ = SparseVector(length(A.p.σ), [idx], zeros(1))
    return (; z, σ)
end

function _tangent_direct_fields!(
        direct_field_tangents::NamedTuple, A::QuadraticSpline, Δt, Δ)
    (; u) = direct_field_tangents
    idx = A.idx_prev[]
    u[idx] = Δ
end

function _tangent_p!(parameter_tangents::NamedTuple, A::QuadraticSpline, Δt, Δ)
    (; z, σ) = parameter_tangents
    idx = A.idx_prev[]
    z[idx] = Δ * Δt
    σ[idx] = Δ * Δt^2
end

## generic

function ChainRulesCore.rrule(A::AType, t::Number) where {AType <: AbstractInterpolation}
    u = A(t)
    idx = get_idx(A.t, t, A.idx_prev[])
    direct_field_tangents = allocate_direct_field_tangents(A)
    parameter_tangents = allocate_parameter_tangents(A)

    function _interpolate_pullback(Δ)
        A.idx_prev[] = idx
        Δt = t - A.t[idx]
        _tangent_direct_fields!(direct_field_tangents, A, Δt, Δ)
        _tangent_p!(parameter_tangents, A, Δt, Δ)
        dA = Tangent{AType}(; direct_field_tangents...,
            p = Tangent{typeof(A.p)}(; parameter_tangents...))
        dt = @thunk(derivative(A, t)*Δ)
        return dA, dt
    end

    u, _interpolate_pullback
end

function ChainRulesCore.frule((_, _, Δt), ::typeof(_interpolate), A::AbstractInterpolation,
        t::Number)
    return _interpolate(A, t), derivative(A, t) * Δt
end

function ChainRulesCore.rrule(::Type{ReadOnlyArray}, parent)
    read_only_array = ReadOnlyArray(parent)
    ReadOnlyArray_pullback(Δ) = NoTangent(), Δ
    read_only_array, ReadOnlyArray_pullback
end

function ChainRulesCore.rrule(::typeof(cumulative_integral), A, u)
    I = cumulative_integral(A, u)
    cumulative_integral_pullback(Δ) = NoTangent(), NoTangent()
    I, cumulative_integral_pullback
end

end # module
