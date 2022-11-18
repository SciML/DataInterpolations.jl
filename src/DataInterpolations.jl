module DataInterpolations

### Interface Functionality

abstract type AbstractInterpolation{FT,T} <: AbstractVector{T} end

Base.size(A::AbstractInterpolation) = size(A.u)
Base.size(A::AbstractInterpolation{true}) = length(A.u) .+ size(A.t)
Base.getindex(A::AbstractInterpolation,i) = A.u[i]
Base.getindex(A::AbstractInterpolation{true},i) = i<=length(A.u) ? A.u[i] : A.t[i-length(A.u)]
Base.setindex!(A::AbstractInterpolation,x,i) = A.u[i] = x
Base.setindex!(A::AbstractInterpolation{true},x,i) =
                                   i <= length(A.u) ? (A.u[i] = x) : (A.t[i-length(A.u)] = x)

using ChainRulesCore, LinearAlgebra, RecursiveArrayTools, RecipesBase, Reexport
@reexport using Optim

include("interpolation_caches.jl")
include("interpolation_utils.jl")
include("interpolation_methods.jl")
include("plot_rec.jl")
include("derivatives.jl")
include("integrals.jl")
include("online.jl")

function ChainRulesCore.rrule(::typeof(_interpolate), A::Union{LagrangeInterpolation,AkimaInterpolation,
                                                               BSplineInterpolation,BSplineApprox}, t::Number)
    interpolate_pullback(Δ) = (NoTangent(), NoTangent(), derivative(A, t) * Δ)
    return _interpolate(A, t), interpolate_pullback
end

ChainRulesCore.frule((_, _, Δt), ::typeof(_interpolate), A::AbstractInterpolation, t::Number) = _interpolate(A, t), derivative(A, t) * Δt

(interp::AbstractInterpolation)(t::Number) = _interpolate(interp, t)

using Symbolics: Num, unwrap, SymbolicUtils
(interp::AbstractInterpolation)(t::Num) = SymbolicUtils.term(interp, unwrap(t))
SymbolicUtils.promote_symtype(t::AbstractInterpolation, _...) = Real


export LinearInterpolation, QuadraticInterpolation, LagrangeInterpolation, AkimaInterpolation,
       ConstantInterpolation, QuadraticSpline, CubicSpline, BSplineInterpolation, BSplineApprox, Curvefit

# Deprecated April 2020
export ZeroSpline

end # module
