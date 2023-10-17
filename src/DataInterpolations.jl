module DataInterpolations

### Interface Functionality

abstract type AbstractInterpolation{FT, T} <: AbstractVector{T} end

Base.size(A::AbstractInterpolation) = size(A.u)
Base.size(A::AbstractInterpolation{true}) = length(A.u) .+ size(A.t)
Base.getindex(A::AbstractInterpolation, i) = A.u[i]
function Base.getindex(A::AbstractInterpolation{true}, i)
    i <= length(A.u) ? A.u[i] : A.t[i - length(A.u)]
end
Base.setindex!(A::AbstractInterpolation, x, i) = A.u[i] = x
function Base.setindex!(A::AbstractInterpolation{true}, x, i)
    i <= length(A.u) ? (A.u[i] = x) : (A.t[i - length(A.u)] = x)
end

using LinearAlgebra, RecursiveArrayTools, RecipesBase

include("interpolation_caches.jl")
include("interpolation_utils.jl")
include("interpolation_methods.jl")
include("plot_rec.jl")
include("derivatives.jl")
include("integrals.jl")
include("online.jl")

(interp::AbstractInterpolation)(t::Number) = _interpolate(interp, t)
(interp::AbstractInterpolation)(t::Number, i::Integer) = _interpolate(interp, t, i)
(interp::AbstractInterpolation)(t::AbstractVector) = interp(similar(t, eltype(interp)), t)
function (interp::AbstractInterpolation)(u::AbstractVector, t::AbstractVector)
    iguess = firstindex(interp.t)
    @inbounds for i in eachindex(u, t)
        u[i], iguess = interp(t[i], iguess)
    end
    u
end

const EXTRAPOLATION_ERROR = "Cannot extrapolate as `extrapolate` keyword passed was `false`"
struct ExtrapolationError <: Exception end
function Base.showerror(io::IO, e::ExtrapolationError)
    print(io, EXTRAPOLATION_ERROR)
end

export LinearInterpolation, QuadraticInterpolation, LagrangeInterpolation,
    AkimaInterpolation, ConstantInterpolation, QuadraticSpline, CubicSpline,
    BSplineInterpolation, BSplineApprox

# added for RegularizationSmooth, JJS 11/27/21
### Regularization data smoothing and interpolation
struct RegularizationSmooth{uType, tType, FT, T, T2} <: AbstractInterpolation{FT, T}
    u::uType
    û::uType
    t::tType
    t̂::tType
    wls::uType
    wr::uType
    d::Int       # derivative degree used to calculate the roughness
    λ::T2        # regularization parameter
    alg::Symbol  # how to determine λ: `:fixed`, `:gcv_svd`, `:gcv_tr`, `L_curve`
    Aitp::AbstractInterpolation{FT, T}
    function RegularizationSmooth{FT}(u, û, t, t̂, wls, wr, d, λ, alg, Aitp) where {FT}
        new{typeof(u), typeof(t), FT, eltype(u), typeof(λ)}(u,
            û,
            t,
            t̂,
            wls,
            wr,
            d,
            λ,
            alg,
            Aitp)
    end
end

export RegularizationSmooth

@static if !isdefined(Base, :get_extension)
    using Requires
end

@static if !isdefined(Base, :get_extension)
    function __init__()
        Requires.@require ChainRulesCore="d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4" begin
            include("../ext/DataInterpolationsChainRulesCoreExt.jl")
        end
        Requires.@require Optim="429524aa-4258-5aef-a3af-852621145aeb" begin
            include("../ext/DataInterpolationsOptimExt.jl")
        end
        Requires.@require RegularizationTools="29dad682-9a27-4bc3-9c72-016788665182" begin
            include("../ext/DataInterpolationsRegularizationToolsExt.jl")
        end
        Requires.@require Symbolics="0c5d862f-8b57-4792-8d23-62f2024744c7" begin
            include("../ext/DataInterpolationsSymbolicsExt.jl")
        end
    end
end

# Define an empty function, so that it can be extended via `DataInterpolationsOptimExt`
Curvefit() = error("CurveFit requires loading Optim, e.g. `using Optim`")

export Curvefit

# Deprecated April 2020
export ZeroSpline

end # module
