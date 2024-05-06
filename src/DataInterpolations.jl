module DataInterpolations

### Interface Functionality

abstract type AbstractInterpolation{FT, T} end

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

using LinearAlgebra, RecipesBase
using PrettyTables
using ForwardDiff
import FindFirstFunctions: searchsortedfirstcorrelated, searchsortedlastcorrelated,
                           bracketstrictlymontonic

include("interpolation_caches.jl")
include("interpolation_utils.jl")
include("interpolation_methods.jl")
include("plot_rec.jl")
include("derivatives.jl")
include("integrals.jl")
include("online.jl")
include("show.jl")

(interp::AbstractInterpolation)(t::Number) = _interpolate(interp, t)
(interp::AbstractInterpolation)(t::Number, i::Integer) = _interpolate(interp, t, i)
function (interp::AbstractInterpolation)(t::AbstractVector)
    interp(similar(t, promote_type(eltype(interp), eltype(t))), t)
end
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

const INTEGRAL_NOT_FOUND_ERROR = "Cannot integrate it analytically. Please use Numerical Integration methods."
struct IntegralNotFoundError <: Exception end
function Base.showerror(io::IO, e::IntegralNotFoundError)
    print(io, INTEGRAL_NOT_FOUND_ERROR)
end

const DERIVATIVE_NOT_FOUND_ERROR = "Derivatives greater than second order is not supported."
struct DerivativeNotFoundError <: Exception end
function Base.showerror(io::IO, e::DerivativeNotFoundError)
    print(io, DERIVATIVE_NOT_FOUND_ERROR)
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
    extrapolate::Bool
    function RegularizationSmooth{FT}(u,
            û,
            t,
            t̂,
            wls,
            wr,
            d,
            λ,
            alg,
            Aitp,
            extrapolate) where {FT}
        new{typeof(u), typeof(t), FT, eltype(u), typeof(λ)}(u,
            û,
            t,
            t̂,
            wls,
            wr,
            d,
            λ,
            alg,
            Aitp,
            extrapolate)
    end
end

export RegularizationSmooth

# CurveFit
struct CurvefitCache{
    uType,
    tType,
    mType,
    p0Type,
    ubType,
    lbType,
    algType,
    pminType,
    FT,
    T
} <: AbstractInterpolation{FT, T}
    u::uType
    t::tType
    m::mType        # model type
    p0::p0Type      # initial params
    ub::ubType      # upper bound of params
    lb::lbType      # lower bound of params
    alg::algType    # alg to optimize cost function
    pmin::pminType  # optimized params
    extrapolate::Bool
    function CurvefitCache{FT}(u, t, m, p0, ub, lb, alg, pmin, extrapolate) where {FT}
        new{typeof(u), typeof(t), typeof(m),
            typeof(p0), typeof(ub), typeof(lb),
            typeof(alg), typeof(pmin), FT, eltype(u)}(u,
            t,
            m,
            p0,
            ub,
            lb,
            alg,
            pmin,
            extrapolate)
    end
end

# Define an empty function, so that it can be extended via `DataInterpolationsOptimExt`
function Curvefit()
    error("CurveFit requires loading Optim and ForwardDiff, e.g. `using Optim, ForwardDiff`")
end

export Curvefit

# Deprecated April 2020
export ZeroSpline

end # module
